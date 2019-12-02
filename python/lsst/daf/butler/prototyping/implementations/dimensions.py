from __future__ import annotations

__all__ = ["DatabaseDimensionTableRecords",
           "SkyPixDimensionTableRecords",
           "DatabaseDimensionTableManager"]

from typing import (
    Optional,
)

import sqlalchemy

from ...core.dimensions import (
    DataCoordinate,
    DimensionElement,
    DimensionRecord,
    DimensionUniverse,
    SkyPixDimension,
)
from ...core.dimensions.schema import OVERLAP_TABLE_NAME_PATTERN, makeOverlapTableSpec
from ...core.schema import TableSpec, FieldSpec
from ...core.utils import NamedKeyDict
from ..interfaces import Database, DimensionTableRecords, DimensionTableManager, StaticTablesContext


class DatabaseDimensionTableRecords(DimensionTableRecords):

    def __init__(self, *, db: Database, element: DimensionElement,
                 table: sqlalchemy.schema.Table,
                 commonSkyPixOverlapTable: Optional[sqlalchemy.schema.Table] = None):
        super().__init__(element=element)
        self._db = db
        self._table = table
        self._commonSkyPixOverlapTable = commonSkyPixOverlapTable

    def insert(self, *records: DimensionRecord):
        # Build lists of dicts to insert first, before any database operations,
        # to minimize the time spent in the transaction.
        elementRows = []
        if self.element.spatial:
            commonSkyPixRows = []
            commonSkyPix = self.element.universe.commonSkyPix
        for record in records:
            elementRows.append(record.toDict())
            if self.element.spatial:
                if record.region is None:
                    # TODO: should we warn about this case?
                    continue
                base = record.dataId.byName()
                for begin, end in commonSkyPix.pixelization.envelope(record.region):
                    for skypix in range(begin, end):
                        row = base.copy()
                        row[commonSkyPix.name] = skypix
                        commonSkyPixRows.append(row)
        with self._db.transaction():
            self._db.insert(self._elementTable, *elementRows)
            if self.element.spatial and commonSkyPixRows:
                self._db.insert(self._commonSkyPixOverlapTable, *commonSkyPixRows)

    def fetch(self, dataId: DataCoordinate) -> DimensionRecord:
        RecordClass = self.element.RecordClass
        # I don't know how expensive it is to construct the query below, and
        # hence how much gain there might be to caching it, so I'm going to
        # wait for it to appear as a hotspot in a profile before trying that.
        nRequired = len(self.element.graph.required)
        if self.element.viewOf is not None:
            whereColumns = [self._table.columns[dimension.name]
                            for dimension in self.element.graph.required]
        else:
            whereColumns = [self._table.columns[fieldName]
                            for fieldName in RecordClass.__slots__[:nRequired]]
        selectColumns = whereColumns + [self._table.columns[name]
                                        for name in RecordClass.__slots__[nRequired:]]
        sql = sqlalchemy.sql.select(
            selectColumns
        ).select_from(self._table).where(
            sqlalchemy.sql.and_(*[column == dataId[dimension.name]
                                for column, dimension in zip(whereColumns, self.element.graph.required)])
        )
        row = self._db.query(sql).fetchone()
        if row is None:
            return None
        return RecordClass(*row)

    def select(self) -> sqlalchemy.sql.FromClause:
        return self._table

    def selectCommonSkyPixOverlap(self) -> Optional[sqlalchemy.sql.FromClause]:
        return self._commonSkyPixOverlapTable


class SkyPixDimensionTableRecords(DimensionTableRecords):
    """A storage implementation specialized for `SkyPixDimension` records.

    `SkyPixDimension` records are never stored in a database, but are instead
    generated-on-the-fly from a `sphgeom.Pixelization` instance.

    Parameters
    ----------
    dimension : `SkyPixDimension`
        The dimension for which this instance will simulate storage.
    """

    def __init__(self, element: SkyPixDimension):
        super().__init__(element=element)

    def select(self) -> None:
        return None

    def selectCommonSkyPixOverlap(self) -> None:
        return None

    def insert(self, *records: DimensionRecord):
        # Docstring inherited from RegistryLayerDimensions.insert.
        raise TypeError(f"Cannot insert into SkyPix dimension {self.element.name}.")

    def fetch(self, dataId: DataCoordinate) -> Optional[DimensionRecord]:
        # Docstring inherited from RegistryLayerDimensions.fetch.
        return self.element.RecordClass(dataId[self.element.name],
                                        self.element.pixelization.pixel(dataId[self.element.name]))


class DatabaseDimensionTableManager(DimensionTableManager):

    _META_TABLE_NAME = "dimension_meta"

    _META_TABLE_SPEC = TableSpec(
        fields=[
            FieldSpec("element_name", dtype=sqlalchemy.String, length=64, primaryKey=True),
        ],
    )

    def __init__(self, db: Database, *, metaTable: sqlalchemy.schema.Table, universe: DimensionUniverse):
        super().__init__(universe=universe)
        self._db = db
        self._
        self._records = NamedKeyDict({})
        self.refresh()

    @classmethod
    def initialize(cls, db: Database, context: StaticTablesContext, *,
                   universe: DimensionUniverse) -> DimensionTableManager:
        metaTable = context.addTable(cls._META_TABLE_NAME, cls._META_TABLE_SPEC)
        return cls(db=db, metaTable=metaTable)

    def refresh(self):
        records = NamedKeyDict({})
        for row in self._db.query(self._metaTable.select()).fetchall():
            element = self.universe[row[self._metaTable.columns.element_name]]
            table = self._db.getExistingTable(element.name, element.makeTableSpec())
            if element.spatial:
                commonSkyPixOverlapTable = self._db.getExistingTable(
                    OVERLAP_TABLE_NAME_PATTERN.format(element.name, self.universe.commonSkyPix.name),
                    makeOverlapTableSpec(element, element.universe.commonSkyPix)
                )
            else:
                commonSkyPixOverlapTable = None
            records[element] = DatabaseDimensionTableRecords(
                db=self._db, element=element, table=table,
                commonSkyPixOverlapTable=commonSkyPixOverlapTable
            )
        self._records = records

    def get(self, element: DimensionElement) -> Optional[DimensionTableRecords]:
        return self._records.get(element)

    def register(self, element: DimensionElement) -> DimensionTableRecords:
        result = self._records.get(element)
        if result is None:
            if isinstance(element, SkyPixDimension):
                result = SkyPixDimensionTableRecords(element)
            elif not element.haElement():
                raise RuntimeError(f"Dimension element subclass {element.__class__} is not supported.")
            else:
                # Create the table itself.  If it already exists but wasn't in
                # the dict because it was added by another client since this
                # one was initialized, that's fine.
                table = self._db.ensureTableExists(element.name, element.makeTableSpec())
                if element.spatial:
                    # Also create a corresponding overlap table with the common
                    # skypix dimension, if this is a spatial element.
                    commonSkyPixOverlapTable = self._db.ensureTableExists(
                        OVERLAP_TABLE_NAME_PATTERN.format(element.name, element.universe.commonSkyPix.name),
                        makeOverlapTableSpec(element, element.universe.commonSkyPix)
                    )
                else:
                    commonSkyPixOverlapTable = None
                # Add a row to the layer_meta table so we can find this table
                # in the future.  Also okay if that already exists, so we use
                # sync.
                self._db.sync(self._metaTable, keys={"element_name": element.name})
                result = DatabaseDimensionTableRecords(
                    db=self._db,
                    element=element,
                    table=table,
                    commonSkyPixOverlapTable=commonSkyPixOverlapTable
                )
            self._records[element] = result
        return result

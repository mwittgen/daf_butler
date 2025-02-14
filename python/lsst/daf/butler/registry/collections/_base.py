# This file is part of daf_butler.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import annotations

__all__ = ()

import itertools
from abc import abstractmethod
from collections import namedtuple
from typing import TYPE_CHECKING, Any, Dict, Generic, Iterable, Iterator, Optional, Tuple, Type, TypeVar

import sqlalchemy

from ...core import DimensionUniverse, Timespan, TimespanDatabaseRepresentation, ddl
from .._collectionType import CollectionType
from .._exceptions import MissingCollectionError
from ..interfaces import ChainedCollectionRecord, CollectionManager, CollectionRecord, RunRecord
from ..wildcards import CollectionSearch

if TYPE_CHECKING:
    from ..interfaces import Database, DimensionRecordStorageManager


def _makeCollectionForeignKey(
    sourceColumnName: str, collectionIdName: str, **kwargs: Any
) -> ddl.ForeignKeySpec:
    """Define foreign key specification that refers to collections table.

    Parameters
    ----------
    sourceColumnName : `str`
        Name of the column in the referring table.
    collectionIdName : `str`
        Name of the column in collections table that identifies it (PK).
    **kwargs
        Additional keyword arguments passed directly to `ddl.ForeignKeySpec`.

    Returns
    -------
    spec : `ddl.ForeignKeySpec`
        Foreign key specification.

    Notes
    -----
    This method assumes fixed name ("collection") of a collections table.
    There is also a general assumption that collection primary key consists
    of a single column.
    """
    return ddl.ForeignKeySpec("collection", source=(sourceColumnName,), target=(collectionIdName,), **kwargs)


CollectionTablesTuple = namedtuple("CollectionTablesTuple", ["collection", "run", "collection_chain"])


def makeRunTableSpec(
    collectionIdName: str, collectionIdType: type, TimespanReprClass: Type[TimespanDatabaseRepresentation]
) -> ddl.TableSpec:
    """Define specification for "run" table.

    Parameters
    ----------
    collectionIdName : `str`
        Name of the column in collections table that identifies it (PK).
    collectionIdType
        Type of the PK column in the collections table, one of the
        `sqlalchemy` types.
    TimespanReprClass : `type` [ `TimespanDatabaseRepresentation` ]
        Subclass of `TimespanDatabaseRepresentation` that encapsulates how
        timespans are stored in this database.


    Returns
    -------
    spec : `ddl.TableSpec`
        Specification for run table.

    Notes
    -----
    Assumption here and in the code below is that the name of the identifying
    column is the same in both collections and run tables. The names of
    non-identifying columns containing run metadata are fixed.
    """
    result = ddl.TableSpec(
        fields=[
            ddl.FieldSpec(collectionIdName, dtype=collectionIdType, primaryKey=True),
            ddl.FieldSpec("host", dtype=sqlalchemy.String, length=128),
        ],
        foreignKeys=[
            _makeCollectionForeignKey(collectionIdName, collectionIdName, onDelete="CASCADE"),
        ],
    )
    for fieldSpec in TimespanReprClass.makeFieldSpecs(nullable=True):
        result.fields.add(fieldSpec)
    return result


def makeCollectionChainTableSpec(collectionIdName: str, collectionIdType: type) -> ddl.TableSpec:
    """Define specification for "collection_chain" table.

    Parameters
    ----------
    collectionIdName : `str`
        Name of the column in collections table that identifies it (PK).
    collectionIdType
        Type of the PK column in the collections table, one of the
        `sqlalchemy` types.

    Returns
    -------
    spec : `ddl.TableSpec`
        Specification for collection chain table.

    Notes
    -----
    Collection chain is simply an ordered one-to-many relation between
    collections. The names of the columns in the table are fixed and
    also hardcoded in the code below.
    """
    return ddl.TableSpec(
        fields=[
            ddl.FieldSpec("parent", dtype=collectionIdType, primaryKey=True),
            ddl.FieldSpec("position", dtype=sqlalchemy.SmallInteger, primaryKey=True),
            ddl.FieldSpec("child", dtype=collectionIdType, nullable=False),
        ],
        foreignKeys=[
            _makeCollectionForeignKey("parent", collectionIdName, onDelete="CASCADE"),
            _makeCollectionForeignKey("child", collectionIdName),
        ],
    )


class DefaultRunRecord(RunRecord):
    """Default `RunRecord` implementation.

    This method assumes the same run table definition as produced by
    `makeRunTableSpec` method. The only non-fixed name in the schema
    is the PK column name, this needs to be passed in a constructor.

    Parameters
    ----------
    db : `Database`
        Registry database.
    key
        Unique collection ID, can be the same as ``name`` if ``name`` is used
        for identification. Usually this is an integer or string, but can be
        other database-specific type.
    name : `str`
        Run collection name.
    table : `sqlalchemy.schema.Table`
        Table for run records.
    idColumnName : `str`
        Name of the identifying column in run table.
    host : `str`, optional
        Name of the host where run was produced.
    timespan : `Timespan`, optional
        Timespan for this run.
    """

    def __init__(
        self,
        db: Database,
        key: Any,
        name: str,
        *,
        table: sqlalchemy.schema.Table,
        idColumnName: str,
        host: Optional[str] = None,
        timespan: Optional[Timespan] = None,
    ):
        super().__init__(key=key, name=name, type=CollectionType.RUN)
        self._db = db
        self._table = table
        self._host = host
        if timespan is None:
            timespan = Timespan(begin=None, end=None)
        self._timespan = timespan
        self._idName = idColumnName

    def update(self, host: Optional[str] = None, timespan: Optional[Timespan] = None) -> None:
        # Docstring inherited from RunRecord.
        if timespan is None:
            timespan = Timespan(begin=None, end=None)
        row = {
            self._idName: self.key,
            "host": host,
        }
        self._db.getTimespanRepresentation().update(timespan, result=row)
        count = self._db.update(self._table, {self._idName: self.key}, row)
        if count != 1:
            raise RuntimeError(f"Run update affected {count} records; expected exactly one.")
        self._host = host
        self._timespan = timespan

    @property
    def host(self) -> Optional[str]:
        # Docstring inherited from RunRecord.
        return self._host

    @property
    def timespan(self) -> Timespan:
        # Docstring inherited from RunRecord.
        return self._timespan


class DefaultChainedCollectionRecord(ChainedCollectionRecord):
    """Default `ChainedCollectionRecord` implementation.

    This method assumes the same chain table definition as produced by
    `makeCollectionChainTableSpec` method. All column names in the table are
    fixed and hard-coded in the methods.

    Parameters
    ----------
    db : `Database`
        Registry database.
    key
        Unique collection ID, can be the same as ``name`` if ``name`` is used
        for identification. Usually this is an integer or string, but can be
        other database-specific type.
    name : `str`
        Collection name.
    table : `sqlalchemy.schema.Table`
        Table for chain relationship records.
    universe : `DimensionUniverse`
        Object managing all known dimensions.
    """

    def __init__(
        self,
        db: Database,
        key: Any,
        name: str,
        *,
        table: sqlalchemy.schema.Table,
        universe: DimensionUniverse,
    ):
        super().__init__(key=key, name=name, universe=universe)
        self._db = db
        self._table = table
        self._universe = universe

    def _update(self, manager: CollectionManager, children: CollectionSearch) -> None:
        # Docstring inherited from ChainedCollectionRecord.
        rows = []
        position = itertools.count()
        for child in children.iter(manager, flattenChains=False):
            rows.append(
                {
                    "parent": self.key,
                    "child": child.key,
                    "position": next(position),
                }
            )
        with self._db.transaction():
            self._db.delete(self._table, ["parent"], {"parent": self.key})
            self._db.insert(self._table, *rows)

    def _load(self, manager: CollectionManager) -> CollectionSearch:
        # Docstring inherited from ChainedCollectionRecord.
        sql = (
            sqlalchemy.sql.select(
                self._table.columns.child,
            )
            .select_from(self._table)
            .where(self._table.columns.parent == self.key)
            .order_by(self._table.columns.position)
        )
        return CollectionSearch.fromExpression(
            [manager[row._mapping[self._table.columns.child]].name for row in self._db.query(sql)]
        )


K = TypeVar("K")


class DefaultCollectionManager(Generic[K], CollectionManager):
    """Default `CollectionManager` implementation.

    This implementation uses record classes defined in this module and is
    based on the same assumptions about schema outlined in the record classes.

    Parameters
    ----------
    db : `Database`
        Interface to the underlying database engine and namespace.
    tables : `CollectionTablesTuple`
        Named tuple of SQLAlchemy table objects.
    collectionIdName : `str`
        Name of the column in collections table that identifies it (PK).
    dimensions : `DimensionRecordStorageManager`
        Manager object for the dimensions in this `Registry`.

    Notes
    -----
    Implementation uses "aggressive" pre-fetching and caching of the records
    in memory. Memory cache is synchronized from database when `refresh`
    method is called.
    """

    def __init__(
        self,
        db: Database,
        tables: CollectionTablesTuple,
        collectionIdName: str,
        *,
        dimensions: DimensionRecordStorageManager,
    ):
        super().__init__()
        self._db = db
        self._tables = tables
        self._collectionIdName = collectionIdName
        self._records: Dict[K, CollectionRecord] = {}  # indexed by record ID
        self._dimensions = dimensions

    def refresh(self) -> None:
        # Docstring inherited from CollectionManager.
        sql = sqlalchemy.sql.select(
            *(list(self._tables.collection.columns) + list(self._tables.run.columns))
        ).select_from(self._tables.collection.join(self._tables.run, isouter=True))
        # Put found records into a temporary instead of updating self._records
        # in place, for exception safety.
        records = []
        chains = []
        TimespanReprClass = self._db.getTimespanRepresentation()
        for row in self._db.query(sql).mappings():
            collection_id = row[self._tables.collection.columns[self._collectionIdName]]
            name = row[self._tables.collection.columns.name]
            type = CollectionType(row["type"])
            record: CollectionRecord
            if type is CollectionType.RUN:
                record = DefaultRunRecord(
                    key=collection_id,
                    name=name,
                    db=self._db,
                    table=self._tables.run,
                    idColumnName=self._collectionIdName,
                    host=row[self._tables.run.columns.host],
                    timespan=TimespanReprClass.extract(row),
                )
            elif type is CollectionType.CHAINED:
                record = DefaultChainedCollectionRecord(
                    db=self._db,
                    key=collection_id,
                    table=self._tables.collection_chain,
                    name=name,
                    universe=self._dimensions.universe,
                )
                chains.append(record)
            else:
                record = CollectionRecord(key=collection_id, name=name, type=type)
            records.append(record)
        self._setRecordCache(records)
        for chain in chains:
            try:
                chain.refresh(self)
            except MissingCollectionError:
                # This indicates a race condition in which some other client
                # created a new collection and added it as a child of this
                # (pre-existing) chain between the time we fetched all
                # collections and the time we queried for parent-child
                # relationships.
                # Because that's some other unrelated client, we shouldn't care
                # about that parent collection anyway, so we just drop it on
                # the floor (a manual refresh can be used to get it back).
                self._removeCachedRecord(chain)

    def register(
        self, name: str, type: CollectionType, doc: Optional[str] = None
    ) -> Tuple[CollectionRecord, bool]:
        # Docstring inherited from CollectionManager.
        registered = False
        record = self._getByName(name)
        if record is None:
            row, inserted_or_updated = self._db.sync(
                self._tables.collection,
                keys={"name": name},
                compared={"type": int(type)},
                extra={"doc": doc},
                returning=[self._collectionIdName],
            )
            assert isinstance(inserted_or_updated, bool)
            registered = inserted_or_updated
            assert row is not None
            collection_id = row[self._collectionIdName]
            if type is CollectionType.RUN:
                TimespanReprClass = self._db.getTimespanRepresentation()
                row, _ = self._db.sync(
                    self._tables.run,
                    keys={self._collectionIdName: collection_id},
                    returning=("host",) + TimespanReprClass.getFieldNames(),
                )
                assert row is not None
                record = DefaultRunRecord(
                    db=self._db,
                    key=collection_id,
                    name=name,
                    table=self._tables.run,
                    idColumnName=self._collectionIdName,
                    host=row["host"],
                    timespan=TimespanReprClass.extract(row),
                )
            elif type is CollectionType.CHAINED:
                record = DefaultChainedCollectionRecord(
                    db=self._db,
                    key=collection_id,
                    name=name,
                    table=self._tables.collection_chain,
                    universe=self._dimensions.universe,
                )
            else:
                record = CollectionRecord(key=collection_id, name=name, type=type)
            self._addCachedRecord(record)
        return record, registered

    def remove(self, name: str) -> None:
        # Docstring inherited from CollectionManager.
        record = self._getByName(name)
        if record is None:
            raise MissingCollectionError(f"No collection with name '{name}' found.")
        # This may raise
        self._db.delete(
            self._tables.collection, [self._collectionIdName], {self._collectionIdName: record.key}
        )
        self._removeCachedRecord(record)

    def find(self, name: str) -> CollectionRecord:
        # Docstring inherited from CollectionManager.
        result = self._getByName(name)
        if result is None:
            raise MissingCollectionError(f"No collection with name '{name}' found.")
        return result

    def __getitem__(self, key: Any) -> CollectionRecord:
        # Docstring inherited from CollectionManager.
        try:
            return self._records[key]
        except KeyError as err:
            raise MissingCollectionError(f"Collection with key '{key}' not found.") from err

    def __iter__(self) -> Iterator[CollectionRecord]:
        yield from self._records.values()

    def getDocumentation(self, key: Any) -> Optional[str]:
        # Docstring inherited from CollectionManager.
        sql = (
            sqlalchemy.sql.select(self._tables.collection.columns.doc)
            .select_from(self._tables.collection)
            .where(self._tables.collection.columns[self._collectionIdName] == key)
        )
        return self._db.query(sql).scalar()

    def setDocumentation(self, key: Any, doc: Optional[str]) -> None:
        # Docstring inherited from CollectionManager.
        self._db.update(self._tables.collection, {self._collectionIdName: "key"}, {"key": key, "doc": doc})

    def _setRecordCache(self, records: Iterable[CollectionRecord]) -> None:
        """Set internal record cache to contain given records,
        old cached records will be removed.
        """
        self._records = {}
        for record in records:
            self._records[record.key] = record

    def _addCachedRecord(self, record: CollectionRecord) -> None:
        """Add single record to cache."""
        self._records[record.key] = record

    def _removeCachedRecord(self, record: CollectionRecord) -> None:
        """Remove single record from cache."""
        del self._records[record.key]

    @abstractmethod
    def _getByName(self, name: str) -> Optional[CollectionRecord]:
        """Find collection record given collection name."""
        raise NotImplementedError()

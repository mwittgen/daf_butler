from __future__ import annotations

__all__ = ("ByDimensionsDatasetRecordStorage",)

import uuid
from typing import TYPE_CHECKING, Any, Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

import sqlalchemy
from lsst.daf.butler import (
    CollectionType,
    DataCoordinate,
    DataCoordinateSet,
    DatasetId,
    DatasetRef,
    DatasetType,
    SimpleQuery,
    Timespan,
    ddl,
)
from lsst.daf.butler.registry import (
    CollectionTypeError,
    ConflictingDefinitionError,
    UnsupportedIdGeneratorError,
)
from lsst.daf.butler.registry.interfaces import DatasetIdGenEnum, DatasetRecordStorage

from ...summaries import GovernorDimensionRestriction
from .tables import makeTagTableSpec

if TYPE_CHECKING:
    from ...interfaces import CollectionManager, CollectionRecord, Database, RunRecord
    from .summaries import CollectionSummaryManager
    from .tables import StaticDatasetTablesTuple


class ByDimensionsDatasetRecordStorage(DatasetRecordStorage):
    """Dataset record storage implementation paired with
    `ByDimensionsDatasetRecordStorageManager`; see that class for more
    information.

    Instances of this class should never be constructed directly; use
    `DatasetRecordStorageManager.register` instead.
    """

    def __init__(
        self,
        *,
        datasetType: DatasetType,
        db: Database,
        dataset_type_id: int,
        collections: CollectionManager,
        static: StaticDatasetTablesTuple,
        summaries: CollectionSummaryManager,
        tags: sqlalchemy.schema.Table,
        calibs: Optional[sqlalchemy.schema.Table],
    ):
        super().__init__(datasetType=datasetType)
        self._dataset_type_id = dataset_type_id
        self._db = db
        self._collections = collections
        self._static = static
        self._summaries = summaries
        self._tags = tags
        self._calibs = calibs
        self._runKeyColumn = collections.getRunForeignKeyName()

    def find(
        self, collection: CollectionRecord, dataId: DataCoordinate, timespan: Optional[Timespan] = None
    ) -> Optional[DatasetRef]:
        # Docstring inherited from DatasetRecordStorage.
        assert dataId.graph == self.datasetType.dimensions
        if collection.type is CollectionType.CALIBRATION and timespan is None:
            raise TypeError(
                f"Cannot search for dataset in CALIBRATION collection {collection.name} "
                f"without an input timespan."
            )
        sql = self.select(
            collection, dataId=dataId, id=SimpleQuery.Select, run=SimpleQuery.Select, timespan=timespan
        )
        results = self._db.query(sql)
        row = results.fetchone()
        if row is None:
            return None
        if collection.type is CollectionType.CALIBRATION:
            # For temporal calibration lookups (only!) our invariants do not
            # guarantee that the number of result rows is <= 1.
            # They would if `select` constrained the given timespan to be
            # _contained_ by the validity range in the self._calibs table,
            # instead of simply _overlapping_ it, because we do guarantee that
            # the validity ranges are disjoint for a particular dataset type,
            # collection, and data ID.  But using an overlap test and a check
            # for multiple result rows here allows us to provide a more useful
            # diagnostic, as well as allowing `select` to support more general
            # queries where multiple results are not an error.
            if results.fetchone() is not None:
                raise RuntimeError(
                    f"Multiple matches found for calibration lookup in {collection.name} for "
                    f"{self.datasetType.name} with {dataId} overlapping {timespan}. "
                )
        return DatasetRef(
            datasetType=self.datasetType,
            dataId=dataId,
            id=row.id,
            run=self._collections[row._mapping[self._runKeyColumn]].name,
        )

    def delete(self, datasets: Iterable[DatasetRef]) -> None:
        # Docstring inherited from DatasetRecordStorage.
        # Only delete from common dataset table; ON DELETE foreign key clauses
        # will handle the rest.
        self._db.delete(
            self._static.dataset,
            ["id"],
            *[{"id": dataset.getCheckedId()} for dataset in datasets],
        )

    def associate(self, collection: CollectionRecord, datasets: Iterable[DatasetRef]) -> None:
        # Docstring inherited from DatasetRecordStorage.
        if collection.type is not CollectionType.TAGGED:
            raise TypeError(
                f"Cannot associate into collection '{collection.name}' "
                f"of type {collection.type.name}; must be TAGGED."
            )
        protoRow = {
            self._collections.getCollectionForeignKeyName(): collection.key,
            "dataset_type_id": self._dataset_type_id,
        }
        rows = []
        governorValues = GovernorDimensionRestriction.makeEmpty(self.datasetType.dimensions.universe)
        for dataset in datasets:
            row = dict(protoRow, dataset_id=dataset.getCheckedId())
            for dimension, value in dataset.dataId.items():
                row[dimension.name] = value
            governorValues.update_extract(dataset.dataId)
            rows.append(row)
        # Update the summary tables for this collection in case this is the
        # first time this dataset type or these governor values will be
        # inserted there.
        self._summaries.update(collection, self.datasetType, self._dataset_type_id, governorValues)
        # Update the tag table itself.
        self._db.replace(self._tags, *rows)

    def disassociate(self, collection: CollectionRecord, datasets: Iterable[DatasetRef]) -> None:
        # Docstring inherited from DatasetRecordStorage.
        if collection.type is not CollectionType.TAGGED:
            raise TypeError(
                f"Cannot disassociate from collection '{collection.name}' "
                f"of type {collection.type.name}; must be TAGGED."
            )
        rows = [
            {
                "dataset_id": dataset.getCheckedId(),
                self._collections.getCollectionForeignKeyName(): collection.key,
            }
            for dataset in datasets
        ]
        self._db.delete(self._tags, ["dataset_id", self._collections.getCollectionForeignKeyName()], *rows)

    def _buildCalibOverlapQuery(
        self, collection: CollectionRecord, dataIds: Optional[DataCoordinateSet], timespan: Timespan
    ) -> SimpleQuery:
        assert self._calibs is not None
        # Start by building a SELECT query for any rows that would overlap
        # this one.
        query = SimpleQuery()
        query.join(self._calibs)
        # Add a WHERE clause matching the dataset type and collection.
        query.where.append(self._calibs.columns.dataset_type_id == self._dataset_type_id)
        query.where.append(
            self._calibs.columns[self._collections.getCollectionForeignKeyName()] == collection.key
        )
        # Add a WHERE clause matching any of the given data IDs.
        if dataIds is not None:
            dataIds.constrain(
                query,
                lambda name: self._calibs.columns[name],  # type: ignore
            )
        # Add WHERE clause for timespan overlaps.
        TimespanReprClass = self._db.getTimespanRepresentation()
        query.where.append(
            TimespanReprClass.fromSelectable(self._calibs).overlaps(TimespanReprClass.fromLiteral(timespan))
        )
        return query

    def certify(
        self, collection: CollectionRecord, datasets: Iterable[DatasetRef], timespan: Timespan
    ) -> None:
        # Docstring inherited from DatasetRecordStorage.
        if self._calibs is None:
            raise CollectionTypeError(
                f"Cannot certify datasets of type {self.datasetType.name}, for which "
                f"DatasetType.isCalibration() is False."
            )
        if collection.type is not CollectionType.CALIBRATION:
            raise CollectionTypeError(
                f"Cannot certify into collection '{collection.name}' "
                f"of type {collection.type.name}; must be CALIBRATION."
            )
        TimespanReprClass = self._db.getTimespanRepresentation()
        protoRow = {
            self._collections.getCollectionForeignKeyName(): collection.key,
            "dataset_type_id": self._dataset_type_id,
        }
        rows = []
        governorValues = GovernorDimensionRestriction.makeEmpty(self.datasetType.dimensions.universe)
        dataIds: Optional[Set[DataCoordinate]] = (
            set() if not TimespanReprClass.hasExclusionConstraint() else None
        )
        for dataset in datasets:
            row = dict(protoRow, dataset_id=dataset.getCheckedId())
            for dimension, value in dataset.dataId.items():
                row[dimension.name] = value
            TimespanReprClass.update(timespan, result=row)
            governorValues.update_extract(dataset.dataId)
            rows.append(row)
            if dataIds is not None:
                dataIds.add(dataset.dataId)
        # Update the summary tables for this collection in case this is the
        # first time this dataset type or these governor values will be
        # inserted there.
        self._summaries.update(collection, self.datasetType, self._dataset_type_id, governorValues)
        # Update the association table itself.
        if TimespanReprClass.hasExclusionConstraint():
            # Rely on database constraint to enforce invariants; we just
            # reraise the exception for consistency across DB engines.
            try:
                self._db.insert(self._calibs, *rows)
            except sqlalchemy.exc.IntegrityError as err:
                raise ConflictingDefinitionError(
                    f"Validity range conflict certifying datasets of type {self.datasetType.name} "
                    f"into {collection.name} for range [{timespan.begin}, {timespan.end})."
                ) from err
        else:
            # Have to implement exclusion constraint ourselves.
            # Start by building a SELECT query for any rows that would overlap
            # this one.
            query = self._buildCalibOverlapQuery(
                collection,
                DataCoordinateSet(dataIds, graph=self.datasetType.dimensions),  # type: ignore
                timespan,
            )
            query.columns.append(sqlalchemy.sql.func.count())
            sql = query.combine()
            # Acquire a table lock to ensure there are no concurrent writes
            # could invalidate our checking before we finish the inserts.  We
            # use a SAVEPOINT in case there is an outer transaction that a
            # failure here should not roll back.
            with self._db.transaction(lock=[self._calibs], savepoint=True):
                # Run the check SELECT query.
                conflicting = self._db.query(sql).scalar()
                if conflicting > 0:
                    raise ConflictingDefinitionError(
                        f"{conflicting} validity range conflicts certifying datasets of type "
                        f"{self.datasetType.name} into {collection.name} for range "
                        f"[{timespan.begin}, {timespan.end})."
                    )
                # Proceed with the insert.
                self._db.insert(self._calibs, *rows)

    def decertify(
        self,
        collection: CollectionRecord,
        timespan: Timespan,
        *,
        dataIds: Optional[Iterable[DataCoordinate]] = None,
    ) -> None:
        # Docstring inherited from DatasetRecordStorage.
        if self._calibs is None:
            raise CollectionTypeError(
                f"Cannot decertify datasets of type {self.datasetType.name}, for which "
                f"DatasetType.isCalibration() is False."
            )
        if collection.type is not CollectionType.CALIBRATION:
            raise CollectionTypeError(
                f"Cannot decertify from collection '{collection.name}' "
                f"of type {collection.type.name}; must be CALIBRATION."
            )
        TimespanReprClass = self._db.getTimespanRepresentation()
        # Construct a SELECT query to find all rows that overlap our inputs.
        dataIdSet: Optional[DataCoordinateSet]
        if dataIds is not None:
            dataIdSet = DataCoordinateSet(set(dataIds), graph=self.datasetType.dimensions)
        else:
            dataIdSet = None
        query = self._buildCalibOverlapQuery(collection, dataIdSet, timespan)
        query.columns.extend(self._calibs.columns)
        sql = query.combine()
        # Set up collections to populate with the rows we'll want to modify.
        # The insert rows will have the same values for collection and
        # dataset type.
        protoInsertRow = {
            self._collections.getCollectionForeignKeyName(): collection.key,
            "dataset_type_id": self._dataset_type_id,
        }
        rowsToDelete = []
        rowsToInsert = []
        # Acquire a table lock to ensure there are no concurrent writes
        # between the SELECT and the DELETE and INSERT queries based on it.
        with self._db.transaction(lock=[self._calibs], savepoint=True):
            for row in self._db.query(sql).mappings():
                rowsToDelete.append({"id": row["id"]})
                # Construct the insert row(s) by copying the prototype row,
                # then adding the dimension column values, then adding what's
                # left of the timespan from that row after we subtract the
                # given timespan.
                newInsertRow = protoInsertRow.copy()
                newInsertRow["dataset_id"] = row["dataset_id"]
                for name in self.datasetType.dimensions.required.names:
                    newInsertRow[name] = row[name]
                rowTimespan = TimespanReprClass.extract(row)
                assert rowTimespan is not None, "Field should have a NOT NULL constraint."
                for diffTimespan in rowTimespan.difference(timespan):
                    rowsToInsert.append(TimespanReprClass.update(diffTimespan, result=newInsertRow.copy()))
            # Run the DELETE and INSERT queries.
            self._db.delete(self._calibs, ["id"], *rowsToDelete)
            self._db.insert(self._calibs, *rowsToInsert)

    def select(
        self,
        *collections: CollectionRecord,
        dataId: SimpleQuery.Select.Or[DataCoordinate] = SimpleQuery.Select,
        id: SimpleQuery.Select.Or[Optional[int]] = SimpleQuery.Select,
        run: SimpleQuery.Select.Or[None] = SimpleQuery.Select,
        timespan: SimpleQuery.Select.Or[Optional[Timespan]] = SimpleQuery.Select,
        ingestDate: SimpleQuery.Select.Or[Optional[Timespan]] = None,
    ) -> sqlalchemy.sql.Selectable:
        # Docstring inherited from DatasetRecordStorage.
        collection_types = {collection.type for collection in collections}
        assert CollectionType.CHAINED not in collection_types, "CHAINED collections must be flattened."
        #
        # There are two kinds of table in play here:
        #
        #  - the static dataset table (with the dataset ID, dataset type ID,
        #    run ID/name, and ingest date);
        #
        #  - the dynamic tags/calibs table (with the dataset ID, dataset type
        #    type ID, collection ID/name, data ID, and possibly validity
        #    range).
        #
        # That means that we might want to return a query against either table
        # or a JOIN of both, depending on which quantities the caller wants.
        # But this method is documented/typed such that ``dataId`` is never
        # `None` - i.e. we always constrain or retreive the data ID.  That
        # means we'll always include the tags/calibs table and join in the
        # static dataset table only if we need things from it that we can't get
        # from the tags/calibs table.
        #
        # Note that it's important that we include a WHERE constraint on both
        # tables for any column (e.g. dataset_type_id) that is in both when
        # it's given explicitly; not doing can prevent the query planner from
        # using very important indexes.  At present, we don't include those
        # redundant columns in the JOIN ON expression, however, because the
        # FOREIGN KEY (and its index) are defined only on dataset_id.
        #
        # We'll start by accumulating kwargs to pass to SimpleQuery.join when
        # we bring in the tags/calibs table.  We get the data ID or constrain
        # it in the tags/calibs table(s), but that's multiple columns, not one,
        # so we need to transform the one Select.Or argument into a dictionary
        # of them.
        kwargs: Dict[str, Any]
        if dataId is SimpleQuery.Select:
            kwargs = {dim.name: SimpleQuery.Select for dim in self.datasetType.dimensions.required}
        else:
            kwargs = dict(dataId.byName())
        # We always constrain (never retrieve) the dataset type in at least the
        # tags/calibs table.
        kwargs["dataset_type_id"] = self._dataset_type_id
        # Join in the tags and/or calibs tables, turning those 'kwargs' entries
        # into WHERE constraints or SELECT columns as appropriate.
        if collection_types != {CollectionType.CALIBRATION}:
            # We'll need a subquery for the tags table if any of the given
            # collections are not a CALIBRATION collection.  This intentionally
            # also fires when the list of collections is empty as a way to
            # create a dummy subquery that we know will fail.
            tags_query = SimpleQuery()
            tags_query.join(self._tags, **kwargs)
            self._finish_single_select(
                tags_query, self._tags, collections, id=id, run=run, ingestDate=ingestDate
            )
        else:
            tags_query = None
        if CollectionType.CALIBRATION in collection_types:
            # If at least one collection is a CALIBRATION collection, we'll
            # need a subquery for the calibs table, and could include the
            # timespan as a result or constraint.
            calibs_query = SimpleQuery()
            assert (
                self._calibs is not None
            ), "DatasetTypes with isCalibration() == False can never be found in a CALIBRATION collection."
            TimespanReprClass = self._db.getTimespanRepresentation()
            # Add the timespan column(s) to the result columns, or constrain
            # the timespan via an overlap condition.
            if timespan is SimpleQuery.Select:
                kwargs.update({k: SimpleQuery.Select for k in TimespanReprClass.getFieldNames()})
            elif timespan is not None:
                calibs_query.where.append(
                    TimespanReprClass.fromSelectable(self._calibs).overlaps(
                        TimespanReprClass.fromLiteral(timespan)
                    )
                )
            calibs_query.join(self._calibs, **kwargs)
            self._finish_single_select(
                calibs_query, self._calibs, collections, id=id, run=run, ingestDate=ingestDate
            )
        else:
            calibs_query = None
        if calibs_query is not None:
            if tags_query is not None:
                if timespan is not None:
                    raise TypeError(
                        "Cannot query for timespan when the collections include both calibration and "
                        "non-calibration collections."
                    )
                return tags_query.combine().union(calibs_query.combine())
            else:
                return calibs_query.combine()
        else:
            assert tags_query is not None, "Earlier logic should guaranteed at least one is not None."
            return tags_query.combine()

    def _finish_single_select(
        self,
        query: SimpleQuery,
        table: sqlalchemy.schema.Table,
        collections: Sequence[CollectionRecord],
        id: SimpleQuery.Select.Or[Optional[int]],
        run: SimpleQuery.Select.Or[None],
        ingestDate: SimpleQuery.Select.Or[Optional[Timespan]],
    ) -> None:
        dataset_id_col = table.columns.dataset_id
        collection_col = table.columns[self._collections.getCollectionForeignKeyName()]
        # We always constrain (never retrieve) the collection(s) in the
        # tags/calibs table.
        if len(collections) == 1:
            query.where.append(collection_col == collections[0].key)
        elif len(collections) == 0:
            # We support the case where there are no collections as a way to
            # generate a valid SQL query that can't yield results.  This should
            # never get executed, but lots of downstream code will still try
            # to access the SQLAlchemy objects representing the columns in the
            # subquery.  That's not ideal, but it'd take a lot of refactoring
            # to fix it (DM-31725).
            query.where.append(sqlalchemy.sql.literal(False))
        else:
            query.where.append(collection_col.in_([collection.key for collection in collections]))
        # We can always get the dataset_id from the tags/calibs table or
        # constrain it there.  Can't use kwargs for that because we need to
        # alias it to 'id'.
        if id is SimpleQuery.Select:
            query.columns.append(dataset_id_col.label("id"))
        elif id is not None:
            query.where.append(dataset_id_col == id)
        # It's possible we now have everything we need, from just the
        # tags/calibs table.  The things we might need to get from the static
        # dataset table are the run key and the ingest date.
        need_static_table = False
        static_kwargs: Dict[str, Any] = {}
        if run is not None:
            assert run is SimpleQuery.Select, "To constrain the run name, pass a RunRecord as a collection."
            if len(collections) == 1 and collections[0].type is CollectionType.RUN:
                # If we are searching exactly one RUN collection, we
                # know that if we find the dataset in that collection,
                # then that's the datasets's run; we don't need to
                # query for it.
                query.columns.append(sqlalchemy.sql.literal(collections[0].key).label(self._runKeyColumn))
            else:
                static_kwargs[self._runKeyColumn] = SimpleQuery.Select
                need_static_table = True
        # Ingest date can only come from the static table.
        if ingestDate is not None:
            need_static_table = True
            if ingestDate is SimpleQuery.Select:
                static_kwargs["ingest_date"] = SimpleQuery.Select
            else:
                assert isinstance(ingestDate, Timespan)
                # Timespan is astropy Time (usually in TAI) and ingest_date is
                # TIMESTAMP, convert values to Python datetime for sqlalchemy.
                if ingestDate.isEmpty():
                    raise RuntimeError("Empty timespan constraint provided for ingest_date.")
                if ingestDate.begin is not None:
                    begin = ingestDate.begin.utc.datetime  # type: ignore
                    query.where.append(self._static.dataset.columns.ingest_date >= begin)
                if ingestDate.end is not None:
                    end = ingestDate.end.utc.datetime  # type: ignore
                    query.where.append(self._static.dataset.columns.ingest_date < end)
        # If we need the static table, join it in via dataset_id and
        # dataset_type_id
        if need_static_table:
            query.join(
                self._static.dataset,
                onclause=(dataset_id_col == self._static.dataset.columns.id),
                **static_kwargs,
            )
            # Also constrain dataset_type_id in static table in case that helps
            # generate a better plan.
            # We could also include this in the JOIN ON clause, but my guess is
            # that that's a good idea IFF it's in the foreign key, and right
            # now it isn't.
            query.where.append(self._static.dataset.columns.dataset_type_id == self._dataset_type_id)

    def getDataId(self, id: DatasetId) -> DataCoordinate:
        """Return DataId for a dataset.

        Parameters
        ----------
        id : `DatasetId`
            Unique dataset identifier.

        Returns
        -------
        dataId : `DataCoordinate`
            DataId for the dataset.
        """
        # This query could return multiple rows (one for each tagged collection
        # the dataset is in, plus one for its run collection), and we don't
        # care which of those we get.
        sql = (
            self._tags.select()
            .where(
                sqlalchemy.sql.and_(
                    self._tags.columns.dataset_id == id,
                    self._tags.columns.dataset_type_id == self._dataset_type_id,
                )
            )
            .limit(1)
        )
        row = self._db.query(sql).mappings().fetchone()
        assert row is not None, "Should be guaranteed by caller and foreign key constraints."
        return DataCoordinate.standardize(
            {dimension.name: row[dimension.name] for dimension in self.datasetType.dimensions.required},
            graph=self.datasetType.dimensions,
        )


class ByDimensionsDatasetRecordStorageInt(ByDimensionsDatasetRecordStorage):
    """Implementation of ByDimensionsDatasetRecordStorage which uses integer
    auto-incremented column for dataset IDs.
    """

    def insert(
        self,
        run: RunRecord,
        dataIds: Iterable[DataCoordinate],
        idMode: DatasetIdGenEnum = DatasetIdGenEnum.UNIQUE,
    ) -> Iterator[DatasetRef]:
        # Docstring inherited from DatasetRecordStorage.

        # We only support UNIQUE mode for integer dataset IDs
        if idMode != DatasetIdGenEnum.UNIQUE:
            raise UnsupportedIdGeneratorError("Only UNIQUE mode can be used with integer dataset IDs.")

        # Transform a possibly-single-pass iterable into a list.
        dataIdList = list(dataIds)
        yield from self._insert(run, dataIdList)

    def import_(
        self,
        run: RunRecord,
        datasets: Iterable[DatasetRef],
        idGenerationMode: DatasetIdGenEnum = DatasetIdGenEnum.UNIQUE,
        reuseIds: bool = False,
    ) -> Iterator[DatasetRef]:
        # Docstring inherited from DatasetRecordStorage.

        # We only support UNIQUE mode for integer dataset IDs
        if idGenerationMode != DatasetIdGenEnum.UNIQUE:
            raise UnsupportedIdGeneratorError("Only UNIQUE mode can be used with integer dataset IDs.")

        # Make a list of dataIds and optionally dataset IDs.
        dataIdList: List[DataCoordinate] = []
        datasetIdList: List[int] = []
        for dataset in datasets:
            dataIdList.append(dataset.dataId)

            # We only accept integer dataset IDs, but also allow None.
            datasetId = dataset.id
            if datasetId is None:
                # if reuseIds is set then all IDs must be known
                if reuseIds:
                    raise TypeError("All dataset IDs must be known if `reuseIds` is set")
            elif isinstance(datasetId, int):
                if reuseIds:
                    datasetIdList.append(datasetId)
            else:
                raise TypeError(f"Unsupported type of dataset ID: {type(datasetId)}")

        yield from self._insert(run, dataIdList, datasetIdList)

    def _insert(
        self, run: RunRecord, dataIdList: List[DataCoordinate], datasetIdList: Optional[List[int]] = None
    ) -> Iterator[DatasetRef]:
        """Common part of implementation of `insert` and `import_` methods."""

        # Remember any governor dimension values we see.
        governorValues = GovernorDimensionRestriction.makeEmpty(self.datasetType.dimensions.universe)
        for dataId in dataIdList:
            governorValues.update_extract(dataId)

        staticRow = {
            "dataset_type_id": self._dataset_type_id,
            self._runKeyColumn: run.key,
        }
        with self._db.transaction():
            # Insert into the static dataset table, generating autoincrement
            # dataset_id values.
            if datasetIdList:
                # reuse existing IDs
                rows = [dict(staticRow, id=datasetId) for datasetId in datasetIdList]
                self._db.insert(self._static.dataset, *rows)
            else:
                # use auto-incremented IDs
                datasetIdList = self._db.insert(
                    self._static.dataset, *([staticRow] * len(dataIdList)), returnIds=True
                )
                assert datasetIdList is not None
            # Update the summary tables for this collection in case this is the
            # first time this dataset type or these governor values will be
            # inserted there.
            self._summaries.update(run, self.datasetType, self._dataset_type_id, governorValues)
            # Combine the generated dataset_id values and data ID fields to
            # form rows to be inserted into the tags table.
            protoTagsRow = {
                "dataset_type_id": self._dataset_type_id,
                self._collections.getCollectionForeignKeyName(): run.key,
            }
            tagsRows = [
                dict(protoTagsRow, dataset_id=dataset_id, **dataId.byName())
                for dataId, dataset_id in zip(dataIdList, datasetIdList)
            ]
            # Insert those rows into the tags table.  This is where we'll
            # get any unique constraint violations.
            self._db.insert(self._tags, *tagsRows)

        for dataId, datasetId in zip(dataIdList, datasetIdList):
            yield DatasetRef(
                datasetType=self.datasetType,
                dataId=dataId,
                id=datasetId,
                run=run.name,
            )


class ByDimensionsDatasetRecordStorageUUID(ByDimensionsDatasetRecordStorage):
    """Implementation of ByDimensionsDatasetRecordStorage which uses UUID for
    dataset IDs.
    """

    NS_UUID = uuid.UUID("840b31d9-05cd-5161-b2c8-00d32b280d0f")
    """Namespace UUID used for UUID5 generation. Do not change. This was
    produced by `uuid.uuid5(uuid.NAMESPACE_DNS, "lsst.org")`.
    """

    def insert(
        self,
        run: RunRecord,
        dataIds: Iterable[DataCoordinate],
        idMode: DatasetIdGenEnum = DatasetIdGenEnum.UNIQUE,
    ) -> Iterator[DatasetRef]:
        # Docstring inherited from DatasetRecordStorage.

        # Remember any governor dimension values we see.
        governorValues = GovernorDimensionRestriction.makeEmpty(self.datasetType.dimensions.universe)

        # Iterate over data IDs, transforming a possibly-single-pass iterable
        # into a list.
        dataIdList = []
        rows = []
        for dataId in dataIds:
            dataIdList.append(dataId)
            rows.append(
                {
                    "id": self._makeDatasetId(run, dataId, idMode),
                    "dataset_type_id": self._dataset_type_id,
                    self._runKeyColumn: run.key,
                }
            )
            governorValues.update_extract(dataId)

        with self._db.transaction():
            # Insert into the static dataset table.
            self._db.insert(self._static.dataset, *rows)
            # Update the summary tables for this collection in case this is the
            # first time this dataset type or these governor values will be
            # inserted there.
            self._summaries.update(run, self.datasetType, self._dataset_type_id, governorValues)
            # Combine the generated dataset_id values and data ID fields to
            # form rows to be inserted into the tags table.
            protoTagsRow = {
                "dataset_type_id": self._dataset_type_id,
                self._collections.getCollectionForeignKeyName(): run.key,
            }
            tagsRows = [
                dict(protoTagsRow, dataset_id=row["id"], **dataId.byName())
                for dataId, row in zip(dataIdList, rows)
            ]
            # Insert those rows into the tags table.
            self._db.insert(self._tags, *tagsRows)

        for dataId, row in zip(dataIdList, rows):
            yield DatasetRef(
                datasetType=self.datasetType,
                dataId=dataId,
                id=row["id"],
                run=run.name,
            )

    def import_(
        self,
        run: RunRecord,
        datasets: Iterable[DatasetRef],
        idGenerationMode: DatasetIdGenEnum = DatasetIdGenEnum.UNIQUE,
        reuseIds: bool = False,
    ) -> Iterator[DatasetRef]:
        # Docstring inherited from DatasetRecordStorage.

        # Remember any governor dimension values we see.
        governorValues = GovernorDimensionRestriction.makeEmpty(self.datasetType.dimensions.universe)

        # Iterate over data IDs, transforming a possibly-single-pass iterable
        # into a list.
        dataIds = {}
        for dataset in datasets:
            # Ignore unknown ID types, normally all IDs have the same type but
            # this code supports mixed types or missing IDs.
            datasetId = dataset.id if isinstance(dataset.id, uuid.UUID) else None
            if datasetId is None:
                datasetId = self._makeDatasetId(run, dataset.dataId, idGenerationMode)
            dataIds[datasetId] = dataset.dataId
            governorValues.update_extract(dataset.dataId)

        with self._db.session() as session:

            # insert all new rows into a temporary table
            tableSpec = makeTagTableSpec(
                self.datasetType, type(self._collections), ddl.GUID, constraints=False
            )
            tmp_tags = session.makeTemporaryTable(tableSpec)

            collFkName = self._collections.getCollectionForeignKeyName()
            protoTagsRow = {
                "dataset_type_id": self._dataset_type_id,
                collFkName: run.key,
            }
            tmpRows = [
                dict(protoTagsRow, dataset_id=dataset_id, **dataId.byName())
                for dataset_id, dataId in dataIds.items()
            ]

            with self._db.transaction():

                # store all incoming data in a temporary table
                self._db.insert(tmp_tags, *tmpRows)

                # There are some checks that we want to make for consistency
                # of the new datasets with existing ones.
                self._validateImport(tmp_tags, run)

                # Before we merge temporary table into dataset/tags we need to
                # drop datasets which are already there (and do not conflict).
                self._db.deleteWhere(
                    tmp_tags,
                    tmp_tags.columns.dataset_id.in_(sqlalchemy.sql.select(self._static.dataset.columns.id)),
                )

                # Copy it into dataset table, need to re-label some columns.
                self._db.insert(
                    self._static.dataset,
                    select=sqlalchemy.sql.select(
                        tmp_tags.columns.dataset_id.label("id"),
                        tmp_tags.columns.dataset_type_id,
                        tmp_tags.columns[collFkName].label(self._runKeyColumn),
                    ),
                )

                # Update the summary tables for this collection in case this
                # is the first time this dataset type or these governor values
                # will be inserted there.
                self._summaries.update(run, self.datasetType, self._dataset_type_id, governorValues)

                # Copy it into tags table.
                self._db.insert(self._tags, select=tmp_tags.select())

                # Return refs in the same order as in the input list.
                for dataset_id, dataId in dataIds.items():
                    yield DatasetRef(
                        datasetType=self.datasetType,
                        id=dataset_id,
                        dataId=dataId,
                        run=run.name,
                    )

    def _validateImport(self, tmp_tags: sqlalchemy.schema.Table, run: RunRecord) -> None:
        """Validate imported refs against existing datasets.

        Parameters
        ----------
        tmp_tags : `sqlalchemy.schema.Table`
            Temporary table with new datasets and the same schema as tags
            table.
        run : `RunRecord`
            The record object describing the `~CollectionType.RUN` collection.

        Raises
        ------
        ConflictingDefinitionError
            Raise if new datasets conflict with existing ones.
        """
        dataset = self._static.dataset
        tags = self._tags
        collFkName = self._collections.getCollectionForeignKeyName()

        # Check that existing datasets have the same dataset type and
        # run.
        query = (
            sqlalchemy.sql.select(
                dataset.columns.id.label("dataset_id"),
                dataset.columns.dataset_type_id.label("dataset_type_id"),
                tmp_tags.columns.dataset_type_id.label("new dataset_type_id"),
                dataset.columns[self._runKeyColumn].label("run"),
                tmp_tags.columns[collFkName].label("new run"),
            )
            .select_from(dataset.join(tmp_tags, dataset.columns.id == tmp_tags.columns.dataset_id))
            .where(
                sqlalchemy.sql.or_(
                    dataset.columns.dataset_type_id != tmp_tags.columns.dataset_type_id,
                    dataset.columns[self._runKeyColumn] != tmp_tags.columns[collFkName],
                )
            )
        )
        result = self._db.query(query)
        if (row := result.first()) is not None:
            # Only include the first one in the exception message
            raise ConflictingDefinitionError(
                f"Existing dataset type or run do not match new dataset: {row._asdict()}"
            )

        # Check that matching dataset in tags table has the same DataId.
        query = (
            sqlalchemy.sql.select(
                tags.columns.dataset_id,
                tags.columns.dataset_type_id.label("type_id"),
                tmp_tags.columns.dataset_type_id.label("new type_id"),
                *[tags.columns[dim] for dim in self.datasetType.dimensions.required.names],
                *[
                    tmp_tags.columns[dim].label(f"new {dim}")
                    for dim in self.datasetType.dimensions.required.names
                ],
            )
            .select_from(tags.join(tmp_tags, tags.columns.dataset_id == tmp_tags.columns.dataset_id))
            .where(
                sqlalchemy.sql.or_(
                    tags.columns.dataset_type_id != tmp_tags.columns.dataset_type_id,
                    *[
                        tags.columns[dim] != tmp_tags.columns[dim]
                        for dim in self.datasetType.dimensions.required.names
                    ],
                )
            )
        )
        result = self._db.query(query)
        if (row := result.first()) is not None:
            # Only include the first one in the exception message
            raise ConflictingDefinitionError(
                f"Existing dataset type or dataId do not match new dataset: {row._asdict()}"
            )

        # Check that matching run+dataId have the same dataset ID.
        query = (
            sqlalchemy.sql.select(
                tags.columns.dataset_type_id.label("dataset_type_id"),
                *[tags.columns[dim] for dim in self.datasetType.dimensions.required.names],
                tags.columns.dataset_id,
                tmp_tags.columns.dataset_id.label("new dataset_id"),
                tags.columns[collFkName],
                tmp_tags.columns[collFkName].label(f"new {collFkName}"),
            )
            .select_from(
                tags.join(
                    tmp_tags,
                    sqlalchemy.sql.and_(
                        tags.columns.dataset_type_id == tmp_tags.columns.dataset_type_id,
                        tags.columns[collFkName] == tmp_tags.columns[collFkName],
                        *[
                            tags.columns[dim] == tmp_tags.columns[dim]
                            for dim in self.datasetType.dimensions.required.names
                        ],
                    ),
                )
            )
            .where(tags.columns.dataset_id != tmp_tags.columns.dataset_id)
        )
        result = self._db.query(query)
        if (row := result.first()) is not None:
            # only include the first one in the exception message
            raise ConflictingDefinitionError(
                f"Existing dataset type and dataId does not match new dataset: {row._asdict()}"
            )

    def _makeDatasetId(
        self, run: RunRecord, dataId: DataCoordinate, idGenerationMode: DatasetIdGenEnum
    ) -> uuid.UUID:
        """Generate dataset ID for a dataset.

        Parameters
        ----------
        run : `RunRecord`
            The record object describing the RUN collection for the dataset.
        dataId : `DataCoordinate`
            Expanded data ID for the dataset.
        idGenerationMode : `DatasetIdGenEnum`
            ID generation option. `~DatasetIdGenEnum.UNIQUE` make a random
            UUID4-type ID. `~DatasetIdGenEnum.DATAID_TYPE` makes a
            deterministic UUID5-type ID based on a dataset type name and
            ``dataId``.  `~DatasetIdGenEnum.DATAID_TYPE_RUN` makes a
            deterministic UUID5-type ID based on a dataset type name, run
            collection name, and ``dataId``.

        Returns
        -------
        datasetId : `uuid.UUID`
            Dataset identifier.
        """
        if idGenerationMode is DatasetIdGenEnum.UNIQUE:
            return uuid.uuid4()
        else:
            # WARNING: If you modify this code make sure that the order of
            # items in the `items` list below never changes.
            items: List[Tuple[str, str]] = []
            if idGenerationMode is DatasetIdGenEnum.DATAID_TYPE:
                items = [
                    ("dataset_type", self.datasetType.name),
                ]
            elif idGenerationMode is DatasetIdGenEnum.DATAID_TYPE_RUN:
                items = [
                    ("dataset_type", self.datasetType.name),
                    ("run", run.name),
                ]
            else:
                raise ValueError(f"Unexpected ID generation mode: {idGenerationMode}")

            for name, value in sorted(dataId.byName().items()):
                items.append((name, str(value)))
            data = ",".join(f"{key}={value}" for key, value in items)
            return uuid.uuid5(self.NS_UUID, data)

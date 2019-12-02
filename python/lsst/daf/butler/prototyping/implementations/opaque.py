from __future__ import annotations

__all__ = ["ByNameOpaqueTableRecords", "ByNameOpaqueTableManager"]

from typing import (
    Any,
    Iterator,
    Optional,
)

import sqlalchemy

from ...core.schema import TableSpec, FieldSpec
from ..interfaces import Database, OpaqueTableManager, OpaqueTableRecords, StaticTablesContext


class ByNameOpaqueTableRecords(OpaqueTableRecords):

    def __init__(self, *, db: Database, name: str, table: sqlalchemy.schema.Table):
        super().__init__(name=name)
        self._db = db
        self._table = table

    def insert(self, *data: dict):
        self._db.insert(self._table, *data)

    def fetch(self, **where: Any) -> Iterator[dict]:
        sql = self._table.select().where(
            sqlalchemy.sql.and_(*[self._table.columns[k] == v for k, v in where.items()])
        )
        yield from self._db.query(sql).fetchall()

    def delete(self, **where: Any):
        self._db.delete(self._table, where.keys(), where)


class ByNameOpaqueTableManager(OpaqueTableManager):

    _META_TABLE_NAME = "opaque_meta"

    _META_TABLE_SPEC = TableSpec(
        fields=[
            FieldSpec("table_name", dtype=sqlalchemy.String, length=128, primaryKey=True),
        ],
    )

    def __init__(self, db: Database, metaTable: sqlalchemy.schema.Table):
        self._db = db
        self._metaTable = metaTable
        self._records = {}

    @classmethod
    def initialize(cls, db: Database, context: StaticTablesContext) -> OpaqueTableManager:
        metaTable = context.addTable(cls._META_TABLE_NAME, cls._META_TABLE_SPEC)
        return cls(db=db, metaTable=metaTable)

    def get(self, name: str) -> Optional[OpaqueTableRecords]:
        return self._records.get(name)

    def register(self, name: str, spec: TableSpec) -> OpaqueTableRecords:
        result = self._records.get(name)
        if result is None:
            # Create the table itself.  If it already exists but wasn't in
            # the dict because it was added by another client since this one
            # was initialized, that's fine.
            table = self._db.ensureTableExists(name, spec)
            # Add a row to the meta table so we can find this table in the
            # future.  Also okay if that already exists, so we use sync.
            self._db.sync(self._metaTable, keys={"table_name": name})
            result = ByNameOpaqueTableRecords(name=name, table=table, db=self._db)
            self._records[name] = result
        return result

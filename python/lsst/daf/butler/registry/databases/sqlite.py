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

__all__ = ["SqliteDatabase"]

import copy
import os
import sqlite3
import urllib.parse
from contextlib import closing
from typing import Any, ContextManager, Iterable, List, Optional

import sqlalchemy
import sqlalchemy.dialects.sqlite
import sqlalchemy.ext.compiler

from ...core import ddl
from ..interfaces import Database, StaticTablesContext


def _onSqlite3Connect(
    dbapiConnection: sqlite3.Connection, connectionRecord: sqlalchemy.pool._ConnectionRecord
) -> None:
    assert isinstance(dbapiConnection, sqlite3.Connection)
    # Prevent pysqlite from emitting BEGIN and COMMIT statements.
    dbapiConnection.isolation_level = None
    # Enable foreign keys
    with closing(dbapiConnection.cursor()) as cursor:
        cursor.execute("PRAGMA foreign_keys=ON;")
        cursor.execute("PRAGMA busy_timeout = 300000;")  # in ms, so 5min (way longer than should be needed)


def _onSqlite3Begin(connection: sqlalchemy.engine.Connection) -> sqlalchemy.engine.Connection:
    assert connection.dialect.name == "sqlite"
    # Replace pysqlite's buggy transaction handling that never BEGINs with our
    # own that does, and tell SQLite to try to acquire a lock as soon as we
    # start a transaction (this should lead to more blocking and fewer
    # deadlocks).
    connection.execute(sqlalchemy.text("BEGIN IMMEDIATE"))
    return connection


class SqliteDatabase(Database):
    """An implementation of the `Database` interface for SQLite3.

    Parameters
    ----------
    connection : `sqlalchemy.engine.Connection`
        An existing connection created by a previous call to `connect`.
    origin : `int`
        An integer ID that should be used as the default for any datasets,
        quanta, or other entities that use a (autoincrement, origin) compound
        primary key.
    namespace : `str`, optional
        The namespace (schema) this database is associated with.  If `None`,
        the default schema for the connection is used (which may be `None`).
    writeable : `bool`, optional
        If `True`, allow write operations on the database, including
        ``CREATE TABLE``.

    Notes
    -----
    The case where ``namespace is not None`` is not yet tested, and may be
    broken; we need an API for attaching to different databases in order to
    write those tests, but haven't yet worked out what is common/different
    across databases well enough to define it.
    """

    def __init__(
        self,
        *,
        engine: sqlalchemy.engine.Engine,
        origin: int,
        namespace: Optional[str] = None,
        writeable: bool = True,
    ):
        super().__init__(origin=origin, engine=engine, namespace=namespace)
        # Get the filename from a call to 'PRAGMA database_list'.
        with engine.connect() as connection:
            with closing(connection.connection.cursor()) as cursor:
                dbList = list(cursor.execute("PRAGMA database_list").fetchall())
        if len(dbList) == 0:
            raise RuntimeError("No database in connection.")
        if namespace is None:
            namespace = "main"
        for _, dbname, filename in dbList:
            if dbname == namespace:
                break
        else:
            raise RuntimeError(f"No '{namespace}' database in connection.")
        if not filename:
            self.filename = None
        else:
            self.filename = filename
        self._writeable = writeable

    @classmethod
    def makeDefaultUri(cls, root: str) -> Optional[str]:
        return "sqlite:///" + os.path.join(root, "gen3.sqlite3")

    @classmethod
    def makeEngine(
        cls, uri: Optional[str] = None, *, filename: Optional[str] = None, writeable: bool = True
    ) -> sqlalchemy.engine.Engine:
        """Create a `sqlalchemy.engine.Engine` from a SQLAlchemy URI or
        filename.

        Parameters
        ----------
        uri : `str`
            A SQLAlchemy URI connection string.
        filename : `str`
            Name of the SQLite database file, or `None` to use an in-memory
            database.  Ignored if ``uri is not None``.
        writeable : `bool`, optional
            If `True`, allow write operations on the database, including
            ``CREATE TABLE``.

        Returns
        -------
        engine : `sqlalchemy.engine.Engine`
            A database engine.
        """
        # In order to be able to tell SQLite that we want a read-only or
        # read-write connection, we need to make the SQLite DBAPI connection
        # with a "URI"-based connection string.  SQLAlchemy claims it can do
        # this
        # (https://docs.sqlalchemy.org/en/13/dialects/sqlite.html#uri-connections),
        # but it doesn't seem to work as advertised.  To work around this, we
        # use the 'creator' argument to sqlalchemy.engine.create_engine, which
        # lets us pass a callable that creates the DBAPI connection.
        if uri is None:
            if filename is None:
                target = ":memory:"
                uri = "sqlite://"
            else:
                target = f"file:{filename}"
                uri = f"sqlite:///{filename}"
        else:
            parsed = urllib.parse.urlparse(uri)
            queries = parsed.query.split("&")
            if "uri=true" in queries:
                # This is a SQLAlchemy URI that is already trying to make a
                # SQLite connection via a SQLite URI, and hence there may
                # be URI components for both SQLite and SQLAlchemy.  We
                # don't need to support that, and it'd be a
                # reimplementation of all of the (broken) logic in
                # SQLAlchemy for doing this, so we just don't.
                raise NotImplementedError("SQLite connection strings with 'uri=true' are not supported.")
            # This is just a SQLAlchemy URI with a non-URI SQLite
            # connection string inside it.  Pull that out so we can use it
            # in the creator call.
            if parsed.path.startswith("/"):
                filename = parsed.path[1:]
                target = f"file:{filename}"
            else:
                filename = None
                target = ":memory:"
        if filename is None:
            if not writeable:
                raise NotImplementedError("Read-only :memory: databases are not supported.")
        else:
            if writeable:
                target += "?mode=rwc&uri=true"
            else:
                target += "?mode=ro&uri=true"

        def creator() -> sqlite3.Connection:
            return sqlite3.connect(target, check_same_thread=False, uri=True)

        engine = sqlalchemy.engine.create_engine(uri, creator=creator)

        sqlalchemy.event.listen(engine, "connect", _onSqlite3Connect)
        sqlalchemy.event.listen(engine, "begin", _onSqlite3Begin)

        return engine

    @classmethod
    def fromEngine(
        cls,
        engine: sqlalchemy.engine.Engine,
        *,
        origin: int,
        namespace: Optional[str] = None,
        writeable: bool = True,
    ) -> Database:
        return cls(engine=engine, origin=origin, writeable=writeable, namespace=namespace)

    def isWriteable(self) -> bool:
        return self._writeable

    def __str__(self) -> str:
        if self.filename:
            return f"SQLite3@{self.filename}"
        else:
            return "SQLite3@:memory:"

    def _lockTables(
        self, connection: sqlalchemy.engine.Connection, tables: Iterable[sqlalchemy.schema.Table] = ()
    ) -> None:
        # Docstring inherited.
        # Our SQLite database always acquires full-database locks at the
        # beginning of a transaction, so there's no need to acquire table-level
        # locks - which is good, because SQLite doesn't have table-level
        # locking.
        pass

    # MyPy claims that the return type here isn't covariant with the return
    # type of the base class method, which is formally correct but irrelevant
    # - the base class return type is _GeneratorContextManager, but only
    # because it's generated by the contextmanager decorator.
    def declareStaticTables(self, *, create: bool) -> ContextManager[StaticTablesContext]:  # type: ignore
        # If the user asked for an in-memory, writeable database, then we may
        # need to re-create schema even if create=False because schema can be
        # lost on re-connect. This is only really relevant for tests, and it's
        # convenient there.
        if self.filename is None and self.isWriteable():
            inspector = sqlalchemy.inspect(self._engine)
            tables = inspector.get_table_names(schema=self.namespace)
            if not tables:
                create = True
        return super().declareStaticTables(create=create)

    def _convertFieldSpec(
        self, table: str, spec: ddl.FieldSpec, metadata: sqlalchemy.MetaData, **kwargs: Any
    ) -> sqlalchemy.schema.Column:
        if spec.autoincrement:
            if not spec.primaryKey:
                raise RuntimeError(
                    f"Autoincrement field {table}.{spec.name} that is not a primary key is not supported."
                )
            if spec.dtype != sqlalchemy.Integer:
                # SQLite's autoincrement is really limited; it only works if
                # the column type is exactly "INTEGER".  But it also doesn't
                # care about the distinctions between different integer types,
                # so it's safe to change it.
                spec = copy.copy(spec)
                spec.dtype = sqlalchemy.Integer
        return super()._convertFieldSpec(table, spec, metadata, **kwargs)

    def _makeColumnConstraints(self, table: str, spec: ddl.FieldSpec) -> List[sqlalchemy.CheckConstraint]:
        # For sqlite we force constraints on all string columns since sqlite
        # ignores everything otherwise and this leads to problems with
        # other databases.

        constraints = []
        if spec.isStringType():
            name = self.shrinkDatabaseEntityName("_".join([table, "len", spec.name]))
            constraints.append(
                sqlalchemy.CheckConstraint(
                    f"length({spec.name})<={spec.length}"
                    # Oracle converts
                    # empty strings to
                    # NULL so check
                    f" AND length({spec.name})>=1",
                    name=name,
                )
            )

        constraints.extend(super()._makeColumnConstraints(table, spec))
        return constraints

    def _convertTableSpec(
        self, name: str, spec: ddl.TableSpec, metadata: sqlalchemy.MetaData, **kwargs: Any
    ) -> sqlalchemy.schema.Table:
        primaryKeyFieldNames = set(field.name for field in spec.fields if field.primaryKey)
        autoincrFieldNames = set(field.name for field in spec.fields if field.autoincrement)
        if len(autoincrFieldNames) > 1:
            raise RuntimeError("At most one autoincrement field per table is allowed.")
        if len(primaryKeyFieldNames) > 1 and len(autoincrFieldNames) > 0:
            # SQLite's default rowid-based autoincrement doesn't work if the
            # field is just one field in a compound primary key.  As a
            # workaround, we create an extra table with just one column that
            # we'll insert into to generate those IDs.  That's only safe if
            # that single-column table's records are already unique with just
            # the autoincrement field, not the rest of the primary key.  In
            # practice, that means the single-column table's records are those
            # for which origin == self.origin.
            (autoincrFieldName,) = autoincrFieldNames
            otherPrimaryKeyFieldNames = primaryKeyFieldNames - autoincrFieldNames
            if otherPrimaryKeyFieldNames != {"origin"}:
                # We need the only other field in the key to be 'origin'.
                raise NotImplementedError(
                    "Compound primary keys with an autoincrement are only supported in SQLite "
                    "if the only non-autoincrement primary key field is 'origin'."
                )
        if not spec.recycleIds:
            kwargs = dict(kwargs, sqlite_autoincrement=True)
        return super()._convertTableSpec(name, spec, metadata, **kwargs)

    def replace(self, table: sqlalchemy.schema.Table, *rows: dict) -> None:
        self.assertTableWriteable(table, f"Cannot replace into read-only table {table}.")
        if not rows:
            return
        query = sqlalchemy.dialects.sqlite.insert(table)
        excluded = query.excluded
        data = {
            column.name: getattr(excluded, column.name)
            for column in table.columns
            if column.name not in table.primary_key
        }
        query = query.on_conflict_do_update(index_elements=table.primary_key, set_=data)
        with self._connection() as connection:
            connection.execute(query, rows)

    def ensure(self, table: sqlalchemy.schema.Table, *rows: dict, primary_key_only: bool = False) -> int:
        self.assertTableWriteable(table, f"Cannot ensure into read-only table {table}.")
        if not rows:
            return 0
        query = sqlalchemy.dialects.sqlite.insert(table)
        if primary_key_only:
            query = query.on_conflict_do_nothing(index_elements=table.primary_key)
        else:
            query = query.on_conflict_do_nothing()
        with self._connection() as connection:
            return connection.execute(query, rows).rowcount

    filename: Optional[str]
    """Name of the file this database is connected to (`str` or `None`).

    Set to `None` for in-memory databases.
    """

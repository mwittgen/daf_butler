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

from .utils import iterable
from sqlalchemy import Column, String, Integer, Boolean, LargeBinary, DateTime,\
    Float, ForeignKey, ForeignKeyConstraint, Table, MetaData
from .config import Config

metadata = None  # Needed to make disabled test_hsc not fail on import

__all__ = ("SchemaConfig", "Schema", "SchemaBuilder")


class SchemaConfig(Config):
    @property
    def tables(self):
        allTables = {}
        allTables.update(self['tables'])
        for dataUnitDescription in self['dataUnits']:
            if 'tables' in dataUnitDescription:
                allTables.update(dataUnitDescription['tables'])
        return allTables


class Schema:
    """The SQL schema for a Butler Registry.

    Parameters
    ----------
    config : `SchemaConfig` or `str`
        Load configuration

    Attributes
    ----------
    metadata : `sqlalchemy.MetaData`
        The sqlalchemy schema description.
    """
    def __init__(self, config):
        if isinstance(config, str):
            config = SchemaConfig(config)
        self.config = config
        self.builder = SchemaBuilder()
        self.buildFromConfig(config)

    def buildFromConfig(self, config):
        for tableName, tableDescription in self.config['tables'].items():
            self.builder.addTable(tableName, tableDescription)
        datasetTable = self.builder.metadata.tables['Dataset']
        self.dataUnits = {}
        for dataUnitDescription in self.config['dataUnits'].values():
            if 'tables' in dataUnitDescription:
                for tableName, tableDescription in dataUnitDescription['tables'].items():
                    self.builder.addTable(tableName, tableDescription)
            if 'link' in dataUnitDescription:
                for dataUnitLinkDescription in dataUnitDescription['link']:
                    linkColumn = self.builder.makeColumn(dataUnitLinkDescription)
                    self.dataUnits[dataUnitLinkDescription['name']] = linkColumn
                    datasetTable.append_column(linkColumn)
        self.metadata = self.builder.metadata


class SchemaBuilder:
    """Builds a Schema step-by-step.

    Attributes
    ----------
    metadata : `sqlalchemy.MetaData`
        The sqlalchemy schema description.
    """
    VALID_COLUMN_TYPES = {'string': String, 'int': Integer, 'float': Float,
                          'bool': Boolean, 'blob': LargeBinary, 'datetime': DateTime}

    def __init__(self):
        self.metadata = MetaData()

    def addTable(self, tableName, tableDescription):
        """Add a table to the schema metadata.

        Parameters
        ----------
        tableName : `str`
            Key of the table.
        tableDescription : `dict`
            Table description.

            Requires:
            - columns, a list of column descriptions
            - foreignKeys, a list of foreign-key constraint descriptions
        """
        table = Table(tableName, self.metadata)
        if "columns" not in tableDescription:
            raise ValueError("No columns in table: {}".format(tableName))
        for columnDescription in tableDescription["columns"]:
            table.append_column(self.makeColumn(columnDescription))
        if "foreignKeys" in tableDescription:
            for constraintDescription in tableDescription["foreignKeys"]:
                table.append_constraint(self.makeForeignKeyConstraint(constraintDescription))

    def makeColumn(self, columnDescription):
        """Make a Column entry for addition to a Table.

        Parameters
        ----------
        columnDescription : `dict`
            Description of the column to be created.
            Should always contain:
            - name, descriptive name
            - type, valid column type
            May contain:
            - nullable, entry can be null
            - primary_key, mark this column as primary key
            - foreign_key, link to other table
            - doc, docstring

        Returns
        -------
        c : `sqlalchemy.Column`
            The created `Column` entry.

        Raises
        ------
        ValueError
            If the column description contains unsupported arguments
        """
        description = columnDescription.copy()
        # required
        columnName = description.pop("name")
        args = (columnName, self.VALID_COLUMN_TYPES[description.pop("type")])
        # foreign_key is special
        if "foreign_key" in description:
            args += (ForeignKey(description.pop("foreign_key")), )
        # additional optional arguments can be passed through directly
        kwargs = {}
        for opt in ("nullable", "primary_key", "doc"):
            if opt in description:
                value = description.pop(opt)
                kwargs[opt] = value
        if description:
            raise ValueError("Unhandled extra kwargs: {} for column: {}".format(description, columnName))
        return Column(*args, **kwargs)

    def makeForeignKeyConstraint(self, constraintDescription):
        """Make a ForeignKeyConstraint for addition to a Table.

        Parameters
        ----------
        constraintDescription : `dict`
            Description of the ForeignKeyConstraint to be created.
            Should always contain:
            - src, list of source column names
            - tgt, list of target column names
        """
        src = tuple(iterable(constraintDescription["src"]))
        tgt = tuple(iterable(constraintDescription["tgt"]))
        return ForeignKeyConstraint(src, tgt)

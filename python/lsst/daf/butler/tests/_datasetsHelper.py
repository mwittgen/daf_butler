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

__all__ = (
    "DatasetTestHelper",
    "DatastoreTestHelper",
    "BadWriteFormatter",
    "BadNoWriteFormatter",
    "MultiDetectorFormatter",
)

import os

from lsst.daf.butler import DatasetRef, DatasetType, StorageClass
from lsst.daf.butler.formatters.yaml import YamlFormatter


class DatasetTestHelper:
    """Helper methods for Datasets"""

    def makeDatasetRef(
        self, datasetTypeName, dimensions, storageClass, dataId, *, id=None, run=None, conform=True
    ):
        """Make a DatasetType and wrap it in a DatasetRef for a test"""
        return self._makeDatasetRef(
            datasetTypeName, dimensions, storageClass, dataId, id=id, run=run, conform=conform
        )

    def _makeDatasetRef(
        self, datasetTypeName, dimensions, storageClass, dataId, *, id=None, run=None, conform=True
    ):
        # helper for makeDatasetRef

        # Pretend we have a parent if this looks like a composite
        compositeName, componentName = DatasetType.splitDatasetTypeName(datasetTypeName)
        parentStorageClass = StorageClass("component") if componentName else None

        datasetType = DatasetType(
            datasetTypeName, dimensions, storageClass, parentStorageClass=parentStorageClass
        )
        if id is None:
            self.id += 1
            id = self.id
        if run is None:
            run = "dummy"
        return DatasetRef(datasetType, dataId, id=id, run=run, conform=conform)


class DatastoreTestHelper:
    """Helper methods for Datastore tests"""

    def setUpDatastoreTests(self, registryClass, configClass):
        """Shared setUp code for all Datastore tests"""
        self.registry = registryClass()

        # Need to keep ID for each datasetRef since we have no butler
        # for these tests
        self.id = 1

        self.config = configClass(self.configFile)

        # Some subclasses override the working root directory
        if self.root is not None:
            self.datastoreType.setConfigRoot(self.root, self.config, self.config.copy())

    def makeDatastore(self, sub=None):
        """Make a new Datastore instance of the appropriate type.

        Parameters
        ----------
        sub : str, optional
            If not None, the returned Datastore will be distinct from any
            Datastore constructed with a different value of ``sub``.  For
            PosixDatastore, for example, the converse is also true, and ``sub``
            is used as a subdirectory to form the new root.

        Returns
        -------
        datastore : `Datastore`
            Datastore constructed by this routine using the supplied
            optional subdirectory if supported.
        """
        config = self.config.copy()
        if sub is not None and self.root is not None:
            self.datastoreType.setConfigRoot(os.path.join(self.root, sub), config, self.config)
        if sub is not None:
            # Ensure that each datastore gets its own registry
            registryClass = type(self.registry)
            registry = registryClass()
        else:
            registry = self.registry
        return self.datastoreType(config=config, bridgeManager=registry.getDatastoreBridgeManager())


class BadWriteFormatter(YamlFormatter):
    """A formatter that never works but does leave a file behind."""

    def _readFile(self, path, pytype=None):
        raise NotImplementedError("This formatter can not read anything")

    def _writeFile(self, inMemoryDataset):
        """Write an empty file and then raise an exception."""
        with open(self.fileDescriptor.location.path, "wb"):
            pass
        raise RuntimeError("Did not succeed in writing file")


class BadNoWriteFormatter(BadWriteFormatter):
    """A formatter that always fails without writing anything."""

    def _writeFile(self, inMemoryDataset):
        raise RuntimeError("Did not writing anything at all")


class MultiDetectorFormatter(YamlFormatter):
    def _writeFile(self, inMemoryDataset):
        raise NotImplementedError("Can not write")

    def _fromBytes(self, serializedDataset, pytype=None):
        data = super()._fromBytes(serializedDataset)
        if self.dataId is None:
            raise RuntimeError("This formatter requires a dataId")
        if "detector" not in self.dataId:
            raise RuntimeError("This formatter requires detector to be present in dataId")
        key = f"detector{self.dataId['detector']}"
        if key in data:
            return pytype(data[key])
        raise RuntimeError(f"Could not find '{key}' in data file")

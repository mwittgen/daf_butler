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

import os
import unittest

import lsst.utils.tests

from lsst.daf.butler import Butler
from lsst.daf.butler.core.datasets import DatasetType, DatasetRef

from datasetsHelper import FitsCatalogDatasetsHelper

"""Tests for Butler.
"""


class ButlerTestCase(lsst.utils.tests.TestCase, FitsCatalogDatasetsHelper):
    """Test for Butler.
    """
    def setUp(self):
        self.testDir = os.path.dirname(__file__)
        self.configFile = os.path.join(self.testDir, "config/basic/butler.yaml")

    def testConstructor(self):
        """Independent test of constructor.
        """
        butler = Butler(self.configFile)
        self.assertIsInstance(butler, Butler)

    @unittest.expectedFailure
    def testBasicPutGet(self):
        butler = Butler(self.configFile)
        # Create and register a DatasetType
        datasetTypeName = "catalog"
        dataUnits = ("camera", "visit")
        storageClass = "SourceCatalog"
        datasetType = DatasetType(datasetTypeName, dataUnits, storageClass)
        butler.registry.registerDatasetType(datasetType)
        # Create and store a dataset
        catalog = self.makeExampleCatalog()
        dataId = {"camera": "DummyCam", "visit": 42}
        ref = butler.put(catalog, datasetTypeName, dataId)
        self.assertIsInstance(ref, DatasetRef)
        catalogOut = catalog.copy()
        self.assertCatalogEqual(catalog, catalogOut)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()

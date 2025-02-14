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

"""Unit tests for `lsst.daf.butler.tests.testRepo`, a module for creating
test repositories or butlers.
"""

import os
import shutil
import unittest

import lsst.daf.butler
from lsst.daf.butler.tests import (
    MetricsExample,
    addDataIdValue,
    addDatasetType,
    expandUniqueId,
    makeTestCollection,
    makeTestRepo,
    registerMetricsExample,
)
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir, safeTestTempDir

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class ButlerUtilsTestSuite(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Repository should be re-created for each test case, but
        # this has a prohibitive run-time cost at present
        cls.root = makeTestTempDir(TESTDIR)

        cls.creatorButler = makeTestRepo(cls.root)
        addDataIdValue(cls.creatorButler, "instrument", "notACam")
        addDataIdValue(cls.creatorButler, "instrument", "dummyCam")
        addDataIdValue(cls.creatorButler, "physical_filter", "k2020", band="k", instrument="notACam")
        addDataIdValue(cls.creatorButler, "physical_filter", "l2019", instrument="dummyCam")
        addDataIdValue(cls.creatorButler, "visit", 101, instrument="notACam", physical_filter="k2020")
        addDataIdValue(cls.creatorButler, "visit", 102, instrument="notACam", physical_filter="k2020")
        addDataIdValue(cls.creatorButler, "detector", 5)
        # Leave skymap/patch/tract undefined so that tests can assume
        # they're missing.

        registerMetricsExample(cls.creatorButler)
        addDatasetType(cls.creatorButler, "DataType1", {"instrument"}, "StructuredDataNoComponents")
        addDatasetType(cls.creatorButler, "DataType2", {"instrument", "visit"}, "StructuredData")

    @classmethod
    def tearDownClass(cls):
        # TODO: use addClassCleanup rather than tearDownClass in Python 3.8
        # to keep the addition and removal together and make it more robust
        removeTestTempDir(cls.root)

    def setUp(self):
        # TestCase.id() is unique for each test method
        self.butler = makeTestCollection(self.creatorButler, uniqueId=self.id())

    def testButlerValid(self):
        self.butler.validateConfiguration()

    def testButlerKwargs(self):
        # outfile has the most obvious effects of any Butler.makeRepo keyword
        with safeTestTempDir(TESTDIR) as temp:
            path = os.path.join(temp, "oddConfig.json")
            makeTestRepo(temp, {}, outfile=path)
            self.assertTrue(os.path.isfile(path))

    def _checkButlerDimension(self, dimensions, query, expected):
        result = list(self.butler.registry.queryDataIds(dimensions, where=query, check=False))
        self.assertEqual(len(result), 1)
        self.assertIn(result[0].byName(), expected)

    def testButlerDimensions(self):
        self._checkButlerDimension(
            {"instrument"}, "instrument='notACam'", [{"instrument": "notACam"}, {"instrument": "dummyCam"}]
        )
        self._checkButlerDimension(
            {"visit", "instrument"},
            "visit=101",
            [{"instrument": "notACam", "visit": 101}, {"instrument": "dummyCam", "visit": 101}],
        )
        self._checkButlerDimension(
            {"visit", "instrument"},
            "visit=102",
            [{"instrument": "notACam", "visit": 102}, {"instrument": "dummyCam", "visit": 102}],
        )
        self._checkButlerDimension(
            {"detector", "instrument"},
            "detector=5",
            [{"instrument": "notACam", "detector": 5}, {"instrument": "dummyCam", "detector": 5}],
        )

    def testAddDataIdValue(self):
        addDataIdValue(self.butler, "visit", 1, instrument="notACam", physical_filter="k2020")
        self._checkButlerDimension(
            {"visit", "instrument"}, "visit=1", [{"instrument": "notACam", "visit": 1}]
        )
        addDataIdValue(self.butler, "visit", 2, instrument="dummyCam", physical_filter="l2019")
        self._checkButlerDimension(
            {"visit", "instrument"}, "visit=2", [{"instrument": "dummyCam", "visit": 2}]
        )

        with self.assertRaises(ValueError):
            addDataIdValue(self.butler, "NotADimension", 42)
        with self.assertRaises(ValueError):
            addDataIdValue(self.butler, "detector", "nonNumeric")
        with self.assertRaises(ValueError):
            addDataIdValue(self.butler, "detector", 101, nonsenseField="string")

        # Keywords imply different instruments
        with self.assertRaises(RuntimeError):
            addDataIdValue(self.butler, "exposure", 101, instrument="dummyCam", physical_filter="k2020")

        # No skymap defined
        with self.assertRaises(RuntimeError):
            addDataIdValue(self.butler, "tract", 42)
        # Didn't create skymap "map" first.
        with self.assertRaises(RuntimeError):
            addDataIdValue(self.butler, "tract", 43, skymap="map")

    def testAddDatasetType(self):
        # 1 for StructuredDataNoComponents, 4 for StructuredData
        self.assertEqual(len(list(self.butler.registry.queryDatasetTypes(components=True))), 5)

        # Testing the DatasetType objects is not practical, because all tests
        # need a DimensionUniverse. So just check that we have the dataset
        # types we expect.
        self.butler.registry.getDatasetType("DataType1")
        self.butler.registry.getDatasetType("DataType2")

        with self.assertRaises(ValueError):
            addDatasetType(self.butler, "DataType3", {"4thDimension"}, "NumpyArray")
        with self.assertRaises(ValueError):
            addDatasetType(self.butler, "DataType3", {"instrument"}, "UnstorableType")

    def testRegisterMetricsExample(self):
        id1 = {"instrument": "notACam"}
        id2 = expandUniqueId(self.butler, {"visit": 101})
        data = MetricsExample(summary={"answer": 42, "question": "unknown"})

        self.butler.put(data, "DataType1", id1)
        self.assertEqual(self.butler.get("DataType1", id1), data)

        self.butler.put(data, "DataType2", id2)
        self.assertEqual(self.butler.get("DataType2", id2), data)
        self.assertEqual(self.butler.get("DataType2.summary", id2), data.summary)

    def testRegisterMetricsExampleChained(self):
        """Regression test for registerMetricsExample having no effect
        on ChainedDatastore.
        """
        temp = makeTestTempDir(TESTDIR)
        try:
            config = lsst.daf.butler.Config()
            config["datastore", "cls"] = "lsst.daf.butler.datastores.chainedDatastore.ChainedDatastore"
            config["datastore", "datastores"] = [
                {
                    "cls": "lsst.daf.butler.datastores.fileDatastore.FileDatastore",
                }
            ]

            repo = lsst.daf.butler.Butler.makeRepo(temp, config=config)
            butler = lsst.daf.butler.Butler(repo, run="chainedExample")
            registerMetricsExample(butler)
            addDatasetType(butler, "DummyType", {}, "StructuredDataNoComponents")

            data = MetricsExample(summary={})
            # Should not raise
            butler.put(data, "DummyType")
        finally:
            shutil.rmtree(temp, ignore_errors=True)

    def testUniqueButler(self):
        dataId = {"instrument": "notACam"}
        self.butler.put(MetricsExample({"answer": 42, "question": "unknown"}), "DataType1", dataId)
        self.assertTrue(self.butler.datasetExists("DataType1", dataId))

        newButler = makeTestCollection(self.creatorButler)
        with self.assertRaises(LookupError):
            newButler.datasetExists("DataType1", dataId)

    def testExpandUniqueId(self):
        self.assertEqual(
            dict(expandUniqueId(self.butler, {"instrument": "notACam"})), {"instrument": "notACam"}
        )
        self.assertIn(
            dict(expandUniqueId(self.butler, {"visit": 101})),
            [{"instrument": "notACam", "visit": 101}, {"instrument": "dummyCam", "visit": 101}],
        )
        self.assertIn(
            dict(expandUniqueId(self.butler, {"detector": 5})),
            [{"instrument": "notACam", "detector": 5}, {"instrument": "dummyCam", "detector": 5}],
        )
        self.assertIn(
            dict(expandUniqueId(self.butler, {"physical_filter": "k2020"})),
            [
                {"instrument": "notACam", "physical_filter": "k2020"},
                {"instrument": "notACam", "physical_filter": "k2020"},
            ],
        )
        with self.assertRaises(ValueError):
            expandUniqueId(self.butler, {"tract": 42})


if __name__ == "__main__":
    unittest.main()

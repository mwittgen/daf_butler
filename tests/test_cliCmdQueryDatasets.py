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

"""Unit tests for daf_butler CLI query-datasets command.
"""

import os
import unittest

from astropy.table import Table as AstropyTable
from lsst.daf.butler import StorageClassFactory, script
from lsst.daf.butler.tests import addDatasetType
from lsst.daf.butler.tests.utils import ButlerTestHelper, MetricTestRepo, makeTestTempDir, removeTestTempDir
from lsst.resources import ResourcePath
from numpy import array

TESTDIR = os.path.abspath(os.path.dirname(__file__))


def expectedFilesystemDatastoreTables(root: ResourcePath):
    return (
        AstropyTable(
            array(
                (
                    (
                        "test_metric_comp.data",
                        "ingest/run",
                        "R",
                        "DummyCamComp",
                        "d-r",
                        "423",
                        root.join(
                            "ingest/run/test_metric_comp.data/"
                            "test_metric_comp_v00000423_fDummyCamComp_data.yaml"
                        ),
                    ),
                    (
                        "test_metric_comp.data",
                        "ingest/run",
                        "R",
                        "DummyCamComp",
                        "d-r",
                        "424",
                        root.join(
                            "ingest/run/test_metric_comp.data/"
                            "test_metric_comp_v00000424_fDummyCamComp_data.yaml"
                        ),
                    ),
                )
            ),
            names=("type", "run", "band", "instrument", "physical_filter", "visit", "URI"),
        ),
        AstropyTable(
            array(
                (
                    (
                        "test_metric_comp.output",
                        "ingest/run",
                        "R",
                        "DummyCamComp",
                        "d-r",
                        "423",
                        root.join(
                            "ingest/run/test_metric_comp.output/"
                            "test_metric_comp_v00000423_fDummyCamComp_output.yaml"
                        ),
                    ),
                    (
                        "test_metric_comp.output",
                        "ingest/run",
                        "R",
                        "DummyCamComp",
                        "d-r",
                        "424",
                        root.join(
                            "ingest/run/test_metric_comp.output/"
                            "test_metric_comp_v00000424_fDummyCamComp_output.yaml"
                        ),
                    ),
                )
            ),
            names=("type", "run", "band", "instrument", "physical_filter", "visit", "URI"),
        ),
        AstropyTable(
            array(
                (
                    (
                        "test_metric_comp.summary",
                        "ingest/run",
                        "R",
                        "DummyCamComp",
                        "d-r",
                        "423",
                        root.join(
                            "ingest/run/test_metric_comp.summary/"
                            "test_metric_comp_v00000423_fDummyCamComp_summary.yaml"
                        ),
                    ),
                    (
                        "test_metric_comp.summary",
                        "ingest/run",
                        "R",
                        "DummyCamComp",
                        "d-r",
                        "424",
                        root.join(
                            "ingest/run/test_metric_comp.summary/"
                            "test_metric_comp_v00000424_fDummyCamComp_summary.yaml"
                        ),
                    ),
                )
            ),
            names=("type", "run", "band", "instrument", "physical_filter", "visit", "URI"),
        ),
    )


class QueryDatasetsTest(unittest.TestCase, ButlerTestHelper):

    configFile = os.path.join(TESTDIR, "config/basic/butler.yaml")
    storageClassFactory = StorageClassFactory()

    @staticmethod
    def _queryDatasets(repo, glob=(), collections=(), where="", find_first=False, show_uri=False):
        return script.QueryDatasets(glob, collections, where, find_first, show_uri, repo=repo).getTables()

    def setUp(self):
        self.testdir = makeTestTempDir(TESTDIR)
        self.repoDir = os.path.join(self.testdir, "repo")

    def tearDown(self):
        removeTestTempDir(self.testdir)

    def testChained(self):
        root = makeTestTempDir(TESTDIR)
        testRepo = MetricTestRepo(root, configFile=os.path.join(TESTDIR, "config/basic/butler-chained.yaml"))

        tables = self._queryDatasets(repo=root, show_uri=True)

        self.assertAstropyTablesEqual(
            tables,
            expectedFilesystemDatastoreTables(testRepo.butler.datastore.datastores[1].root),
            filterColumns=True,
        )

    def testShowURI(self):
        """Test for expected output with show_uri=True."""
        testRepo = MetricTestRepo(self.repoDir, configFile=self.configFile)

        tables = self._queryDatasets(repo=self.repoDir, show_uri=True)

        self.assertAstropyTablesEqual(
            tables, expectedFilesystemDatastoreTables(testRepo.butler.datastore.root), filterColumns=True
        )

    def testNoShowURI(self):
        """Test for expected output without show_uri (default is False)."""
        _ = MetricTestRepo(self.repoDir, configFile=self.configFile)

        tables = self._queryDatasets(repo=self.repoDir)

        expectedTables = (
            AstropyTable(
                array(
                    (
                        ("test_metric_comp", "ingest/run", "R", "DummyCamComp", "d-r", "423"),
                        ("test_metric_comp", "ingest/run", "R", "DummyCamComp", "d-r", "424"),
                    )
                ),
                names=("type", "run", "band", "instrument", "physical_filter", "visit"),
            ),
        )

        self.assertAstropyTablesEqual(tables, expectedTables, filterColumns=True)

    def testWhere(self):
        """Test using the where clause to reduce the number of rows returned by
        queryDatasets.
        """
        _ = MetricTestRepo(self.repoDir, configFile=self.configFile)

        tables = self._queryDatasets(repo=self.repoDir, where="instrument='DummyCamComp' AND visit=423")

        expectedTables = (
            AstropyTable(
                array(("test_metric_comp", "ingest/run", "R", "DummyCamComp", "d-r", "423")),
                names=("type", "run", "band", "instrument", "physical_filter", "visit"),
            ),
        )

        self.assertAstropyTablesEqual(tables, expectedTables, filterColumns=True)

    def testGlobDatasetType(self):
        """Test specifying dataset type."""
        # Create and register an additional DatasetType
        testRepo = MetricTestRepo(self.repoDir, configFile=self.configFile)

        testRepo.butler.registry.insertDimensionData(
            "visit",
            {"instrument": "DummyCamComp", "id": 425, "name": "fourtwentyfive", "physical_filter": "d-r"},
        )

        datasetType = addDatasetType(
            testRepo.butler, "alt_test_metric_comp", ("instrument", "visit"), "StructuredCompositeReadComp"
        )

        testRepo.addDataset(dataId={"instrument": "DummyCamComp", "visit": 425}, datasetType=datasetType)

        # verify the new dataset type increases the number of tables found:
        tables = self._queryDatasets(repo=self.repoDir)

        expectedTables = (
            AstropyTable(
                array(
                    (
                        ("test_metric_comp", "ingest/run", "R", "DummyCamComp", "d-r", "423"),
                        ("test_metric_comp", "ingest/run", "R", "DummyCamComp", "d-r", "424"),
                    )
                ),
                names=("type", "run", "band", "instrument", "physical_filter", "visit"),
            ),
            AstropyTable(
                array((("alt_test_metric_comp", "ingest/run", "R", "DummyCamComp", "d-r", "425"))),
                names=("type", "run", "band", "instrument", "physical_filter", "visit"),
            ),
        )

        self.assertAstropyTablesEqual(tables, expectedTables, filterColumns=True)

    def testFindFirstAndCollections(self):
        """Test the find-first option, and the collections option, since it
        is required for find-first."""

        testRepo = MetricTestRepo(self.repoDir, configFile=self.configFile)

        # Add a new run, and add a dataset to shadow an existing dataset.
        testRepo.addDataset(run="foo", dataId={"instrument": "DummyCamComp", "visit": 424})

        # Verify that without find-first, duplicate datasets are returned
        tables = self._queryDatasets(repo=self.repoDir, collections=["foo", "ingest/run"], show_uri=True)

        expectedTables = (
            AstropyTable(
                array(
                    (
                        (
                            "test_metric_comp.data",
                            "foo",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "424",
                            testRepo.butler.datastore.root.join(
                                "foo/test_metric_comp.data/"
                                "test_metric_comp_v00000424_fDummyCamComp_data.yaml"
                            ),
                        ),
                        (
                            "test_metric_comp.data",
                            "ingest/run",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "423",
                            testRepo.butler.datastore.root.join(
                                "ingest/run/test_metric_comp.data/"
                                "test_metric_comp_v00000423_fDummyCamComp_data.yaml"
                            ),
                        ),
                        (
                            "test_metric_comp.data",
                            "ingest/run",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "424",
                            testRepo.butler.datastore.root.join(
                                "ingest/run/test_metric_comp.data/"
                                "test_metric_comp_v00000424_fDummyCamComp_data.yaml"
                            ),
                        ),
                    )
                ),
                names=("type", "run", "band", "instrument", "physical_filter", "visit", "URI"),
            ),
            AstropyTable(
                array(
                    (
                        (
                            "test_metric_comp.output",
                            "foo",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "424",
                            testRepo.butler.datastore.root.join(
                                "foo/test_metric_comp.output/"
                                "test_metric_comp_v00000424_fDummyCamComp_output.yaml"
                            ),
                        ),
                        (
                            "test_metric_comp.output",
                            "ingest/run",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "423",
                            testRepo.butler.datastore.root.join(
                                "ingest/run/test_metric_comp.output/"
                                "test_metric_comp_v00000423_fDummyCamComp_output.yaml"
                            ),
                        ),
                        (
                            "test_metric_comp.output",
                            "ingest/run",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "424",
                            testRepo.butler.datastore.root.join(
                                "ingest/run/test_metric_comp.output/"
                                "test_metric_comp_v00000424_fDummyCamComp_output.yaml"
                            ),
                        ),
                    )
                ),
                names=("type", "run", "band", "instrument", "physical_filter", "visit", "URI"),
            ),
            AstropyTable(
                array(
                    (
                        (
                            "test_metric_comp.summary",
                            "foo",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "424",
                            testRepo.butler.datastore.root.join(
                                "foo/test_metric_comp.summary/"
                                "test_metric_comp_v00000424_fDummyCamComp_summary.yaml"
                            ),
                        ),
                        (
                            "test_metric_comp.summary",
                            "ingest/run",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "423",
                            testRepo.butler.datastore.root.join(
                                "ingest/run/test_metric_comp.summary/"
                                "test_metric_comp_v00000423_fDummyCamComp_summary.yaml"
                            ),
                        ),
                        (
                            "test_metric_comp.summary",
                            "ingest/run",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "424",
                            testRepo.butler.datastore.root.join(
                                "ingest/run/test_metric_comp.summary/"
                                "test_metric_comp_v00000424_fDummyCamComp_summary.yaml"
                            ),
                        ),
                    )
                ),
                names=("type", "run", "band", "instrument", "physical_filter", "visit", "URI"),
            ),
        )

        self.assertAstropyTablesEqual(tables, expectedTables, filterColumns=True)

        # Verify that with find first the duplicate dataset is eliminated and
        # the more recent dataset is returned.
        tables = self._queryDatasets(
            repo=self.repoDir, collections=["foo", "ingest/run"], show_uri=True, find_first=True
        )

        expectedTables = (
            AstropyTable(
                array(
                    (
                        (
                            "test_metric_comp.data",
                            "foo",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "424",
                            testRepo.butler.datastore.root.join(
                                "foo/test_metric_comp.data/test_metric_comp_v00000424_fDummyCamComp_data.yaml"
                            ),
                        ),
                        (
                            "test_metric_comp.data",
                            "ingest/run",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "423",
                            testRepo.butler.datastore.root.join(
                                "ingest/run/test_metric_comp.data/"
                                "test_metric_comp_v00000423_fDummyCamComp_data.yaml"
                            ),
                        ),
                    )
                ),
                names=("type", "run", "band", "instrument", "physical_filter", "visit", "URI"),
            ),
            AstropyTable(
                array(
                    (
                        (
                            "test_metric_comp.output",
                            "foo",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "424",
                            testRepo.butler.datastore.root.join(
                                "foo/test_metric_comp.output/"
                                "test_metric_comp_v00000424_fDummyCamComp_output.yaml"
                            ),
                        ),
                        (
                            "test_metric_comp.output",
                            "ingest/run",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "423",
                            testRepo.butler.datastore.root.join(
                                "ingest/run/test_metric_comp.output/"
                                "test_metric_comp_v00000423_fDummyCamComp_output.yaml"
                            ),
                        ),
                    )
                ),
                names=("type", "run", "band", "instrument", "physical_filter", "visit", "URI"),
            ),
            AstropyTable(
                array(
                    (
                        (
                            "test_metric_comp.summary",
                            "foo",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "424",
                            testRepo.butler.datastore.root.join(
                                "foo/test_metric_comp.summary/"
                                "test_metric_comp_v00000424_fDummyCamComp_summary.yaml"
                            ),
                        ),
                        (
                            "test_metric_comp.summary",
                            "ingest/run",
                            "R",
                            "DummyCamComp",
                            "d-r",
                            "423",
                            testRepo.butler.datastore.root.join(
                                "ingest/run/test_metric_comp.summary/"
                                "test_metric_comp_v00000423_fDummyCamComp_summary.yaml"
                            ),
                        ),
                    )
                ),
                names=("type", "run", "band", "instrument", "physical_filter", "visit", "URI"),
            ),
        )

        self.assertAstropyTablesEqual(tables, expectedTables, filterColumns=True)


if __name__ == "__main__":
    unittest.main()

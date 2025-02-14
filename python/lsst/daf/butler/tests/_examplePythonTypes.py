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

"""
Python classes that can be used to test datastores without requiring
large external dependencies on python classes such as afw or serialization
formats such as FITS or HDF5.
"""

__all__ = (
    "ListDelegate",
    "MetricsDelegate",
    "MetricsExample",
    "registerMetricsExample",
    "MetricsExampleModel",
)


import copy
from typing import Any, Dict, List, Optional

from lsst.daf.butler import StorageClass, StorageClassDelegate
from pydantic import BaseModel


def registerMetricsExample(butler):
    """Modify a repository to support reading and writing
    `MetricsExample` objects.

    This method allows `MetricsExample` to be used with test repositories
    in any package without needing to provide a custom configuration there.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        The repository that needs to support `MetricsExample`.

    Notes
    -----
    This method enables the following storage classes:

    ``StructuredData``
        A `MetricsExample` whose ``summary``, ``output``, and ``data`` members
        can be retrieved as dataset components.
    ``StructuredDataNoComponents``
        A monolithic write of a `MetricsExample`.
    """
    yamlDict = _addFullStorageClass(
        butler,
        "StructuredDataDictYaml",
        "lsst.daf.butler.formatters.yaml.YamlFormatter",
        pytype=dict,
    )

    yamlList = _addFullStorageClass(
        butler,
        "StructuredDataListYaml",
        "lsst.daf.butler.formatters.yaml.YamlFormatter",
        pytype=list,
        parameters={"slice"},
        delegate="lsst.daf.butler.tests.ListDelegate",
    )

    _addFullStorageClass(
        butler,
        "StructuredDataNoComponents",
        "lsst.daf.butler.formatters.pickle.PickleFormatter",
        pytype=MetricsExample,
        parameters={"slice"},
        delegate="lsst.daf.butler.tests.MetricsDelegate",
    )

    _addFullStorageClass(
        butler,
        "StructuredData",
        "lsst.daf.butler.formatters.yaml.YamlFormatter",
        pytype=MetricsExample,
        components={
            "summary": yamlDict,
            "output": yamlDict,
            "data": yamlList,
        },
        delegate="lsst.daf.butler.tests.MetricsDelegate",
    )


def _addFullStorageClass(butler, name, formatter, *args, **kwargs):
    """Create a storage class-formatter pair in a repository if it does not
    already exist.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        The repository that needs to contain the class.
    name : `str`
        The name to use for the class.
    formatter : `str`
        The formatter to use with the storage class. Ignored if ``butler``
        does not use formatters.
    *args
    **kwargs
        Arguments, other than ``name``, to the `~lsst.daf.butler.StorageClass`
        constructor.

    Returns
    -------
    class : `lsst.daf.butler.StorageClass`
        The newly created storage class, or the class of the same name
        previously found in the repository.
    """
    storageRegistry = butler.datastore.storageClassFactory

    storage = StorageClass(name, *args, **kwargs)
    try:
        storageRegistry.registerStorageClass(storage)
    except ValueError:
        storage = storageRegistry.getStorageClass(name)

    for registry in _getAllFormatterRegistries(butler.datastore):
        registry.registerFormatter(storage, formatter)

    return storage


def _getAllFormatterRegistries(datastore):
    """Return all formatter registries used by a datastore.

    Parameters
    ----------
    datastore : `lsst.daf.butler.Datastore`
        A datastore containing zero or more formatter registries.

    Returns
    -------
    registries : `list` [`lsst.daf.butler.FormatterRegistry`]
        A possibly empty list of all formatter registries used
        by ``datastore``.
    """
    try:
        datastores = datastore.datastores
    except AttributeError:
        datastores = [datastore]

    registries = []
    for datastore in datastores:
        try:
            # Not all datastores have a formatterFactory
            formatterRegistry = datastore.formatterFactory
        except AttributeError:
            pass  # no formatter needed
        else:
            registries.append(formatterRegistry)
    return registries


class MetricsExample:
    """Smorgasboard of information that might be the result of some
    processing.

    Parameters
    ----------
    summary : `dict`
        Simple dictionary mapping key performance metrics to a scalar
        result.
    output : `dict`
        Structured nested data.
    data : `list`, optional
        Arbitrary array data.
    """

    def __init__(self, summary=None, output=None, data=None):
        self.summary = summary
        self.output = output
        self.data = data

    def __eq__(self, other):
        return self.summary == other.summary and self.output == other.output and self.data == other.data

    def __str__(self):
        return str(self.exportAsDict())

    def __repr__(self):
        return f"MetricsExample({self.exportAsDict()})"

    def exportAsDict(self):
        """Convert object contents to a single python dict."""
        exportDict = {"summary": self.summary, "output": self.output}
        if self.data is not None:
            exportDict["data"] = list(self.data)
        else:
            exportDict["data"] = None
        return exportDict

    def _asdict(self):
        """Convert object contents to a single Python dict.

        This interface is used for JSON serialization.

        Returns
        -------
        exportDict : `dict`
            Object contents in the form of a dict with keys corresponding
            to object attributes.
        """
        return self.exportAsDict()

    @classmethod
    def makeFromDict(cls, exportDict):
        """Create a new object from a dict that is compatible with that
        created by `exportAsDict`.

        Parameters
        ----------
        exportDict : `dict`
            `dict` with keys "summary", "output", and (optionally) "data".

        Returns
        -------
        newobject : `MetricsExample`
            New `MetricsExample` object.
        """
        data = None
        if "data" in exportDict:
            data = exportDict["data"]
        return cls(exportDict["summary"], exportDict["output"], data)


class MetricsExampleModel(BaseModel):
    """A variant of `MetricsExample` based on model."""

    summary: Optional[Dict[str, Any]]
    output: Optional[Dict[str, Any]]
    data: Optional[List[Any]]

    @classmethod
    def from_metrics(cls, metrics: MetricsExample) -> "MetricsExampleModel":
        """Create a model based on an example."""
        return cls.parse_obj(metrics.exportAsDict())


class ListDelegate(StorageClassDelegate):
    """Parameter handler for list parameters"""

    def handleParameters(self, inMemoryDataset, parameters=None):
        """Modify the in-memory dataset using the supplied parameters,
        returning a possibly new object.

        Parameters
        ----------
        inMemoryDataset : `object`
            Object to modify based on the parameters.
        parameters : `dict`
            Parameters to apply. Values are specific to the parameter.
            Supported parameters are defined in the associated
            `StorageClass`.  If no relevant parameters are specified the
            inMemoryDataset will be return unchanged.

        Returns
        -------
        inMemoryDataset : `object`
            Updated form of supplied in-memory dataset, after parameters
            have been used.
        """
        inMemoryDataset = copy.deepcopy(inMemoryDataset)
        use = self.storageClass.filterParameters(parameters, subset={"slice"})
        if use:
            inMemoryDataset = inMemoryDataset[use["slice"]]
        return inMemoryDataset


class MetricsDelegate(StorageClassDelegate):
    """Parameter handler for parameters using Metrics"""

    def handleParameters(self, inMemoryDataset, parameters=None):
        """Modify the in-memory dataset using the supplied parameters,
        returning a possibly new object.

        Parameters
        ----------
        inMemoryDataset : `object`
            Object to modify based on the parameters.
        parameters : `dict`
            Parameters to apply. Values are specific to the parameter.
            Supported parameters are defined in the associated
            `StorageClass`.  If no relevant parameters are specified the
            inMemoryDataset will be return unchanged.

        Returns
        -------
        inMemoryDataset : `object`
            Updated form of supplied in-memory dataset, after parameters
            have been used.
        """
        inMemoryDataset = copy.deepcopy(inMemoryDataset)
        use = self.storageClass.filterParameters(parameters, subset={"slice"})
        if use:
            inMemoryDataset.data = inMemoryDataset.data[use["slice"]]
        return inMemoryDataset

    def getComponent(self, composite, componentName: str):
        if componentName == "counter":
            return len(composite.data)
        return super().getComponent(composite, componentName)

    @classmethod
    def selectResponsibleComponent(cls, readComponent: str, fromComponents) -> str:
        forwarderMap = {
            "counter": "data",
        }
        forwarder = forwarderMap.get(readComponent)
        if forwarder is not None and forwarder in fromComponents:
            return forwarder
        raise ValueError(f"Can not calculate read component {readComponent} from {fromComponents}")

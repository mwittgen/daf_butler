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

__all__ = ("ingest_files",)

import logging
from collections import defaultdict
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from astropy.table import Table
from lsst.resources import ResourcePath
from lsst.utils import doImport

from .._butler import Butler
from ..core import DatasetRef, FileDataset
from ..registry import DatasetIdGenEnum

if TYPE_CHECKING:
    from ..core import DatasetType, DimensionUniverse

log = logging.getLogger(__name__)


def ingest_files(
    repo: str,
    dataset_type: str,
    run: str,
    table_file: str,
    data_id: Tuple[str, ...] = (),
    formatter: Optional[str] = None,
    id_generation_mode: str = "UNIQUE",
    prefix: Optional[str] = None,
    transfer: str = "auto",
) -> None:
    """Ingest files from a table.

    Parameters
    ----------
    repo : `str`
        URI string of the Butler repo to use.
    dataset_type : `str`
        The name of the dataset type for the files to be ingested. This
        dataset type must exist.
    run : `str`
        The run in which the files should be ingested.
    table_file : `str`
        Path to a table file to read. This file can be in any format that
        can be read by Astropy so long as Astropy can determine the format
        itself.
    data_id : `tuple` of `str`
        Tuple of strings of the form ``keyword=value`` that can be used
        to define dataId elements that are fixed for all ingested files
        found in the table file. This allows those columns to be missing
        from the table file. Dimensions given here override table columns.
    formatter : `str`, optional
        Fully-qualified python class name for the `Formatter` to use
        to read the ingested files. If `None` the formatter is read from
        datastore configuration based on the dataset type.
    id_generation_mode : `str`, optional
        Mode to use for generating IDs.  Should map to `DatasetGenIdEnum`.
    prefix : `str`, optional
        Prefix to use when resolving relative paths in table files. The default
        is to use the current working directory.
    transfer : `str`, optional
        Transfer mode to use for ingest.
    """

    # Check that the formatter can be imported -- validate this as soon
    # as possible before we read a potentially large table file.
    if formatter:
        doImport(formatter)
    else:
        formatter = None

    # Force empty string prefix (from click) to None for API compatibility.
    if not prefix:
        prefix = None

    # Convert the dataset ID gen mode string to enum.
    id_gen_mode = DatasetIdGenEnum.__members__[id_generation_mode]

    # Create the butler with the relevant run attached.
    butler = Butler(repo, run=run)

    datasetType = butler.registry.getDatasetType(dataset_type)

    # Convert the k=v strings into a dataId dict.
    universe = butler.registry.dimensions
    common_data_id = parse_data_id_tuple(data_id, universe)

    # Read the table assuming that Astropy can work out the format.
    table = Table.read(table_file)

    datasets = extract_datasets_from_table(table, common_data_id, datasetType, formatter, prefix)

    butler.ingest(*datasets, transfer=transfer, run=run, idGenerationMode=id_gen_mode)


def extract_datasets_from_table(
    table: Table,
    common_data_id: Dict,
    datasetType: DatasetType,
    formatter: Optional[str] = None,
    prefix: Optional[str] = None,
) -> List[FileDataset]:
    """Extract datasets from the supplied table.

    Parameters
    ----------
    table : `astropy.table.Table`
        Table containing the datasets. The first column is assumed to be
        the file URI and the remaining columns are dimensions.
    common_data_id : `dict`
        Data ID values that are common to every row in the table. These
        take priority if a dimension in this dataId is also present as
        a column in the table.
    datasetType : `DatasetType`
        The dataset type to be associated with the ingested data.
    formatter : `str`, optional
        Fully-qualified python class name for the `Formatter` to use
        to read the ingested files. If `None` the formatter is read from
        datastore configuration based on the dataset type.
    prefix : `str`
        Prefix to be used for relative paths. Can be `None` for current
        working directory.

    Returns
    -------
    datasets : `list` of `FileDataset`
        The `FileDataset` objects corresponding to the rows in the table.
        The number of elements in this list can be smaller than the number
        of rows in the file because one file can appear in multiple rows
        with different dataIds.
    """
    # The file is the first column and everything else is assumed to
    # be dimensions so we need to know the name of that column.
    file_column = table.colnames[0]

    # Handle multiple dataIds per file by grouping by file.
    refs_by_file = defaultdict(list)
    n_dataset_refs = 0
    for row in table:

        # Convert the row to a dataId, remembering to extract the
        # path column.
        dataId = dict(row)
        path = dataId.pop(file_column)

        # The command line can override a column.
        dataId.update(common_data_id)

        # Create the dataset ref that is to be ingested.
        ref = DatasetRef(datasetType, dataId)  # type: ignore

        # Convert path to absolute (because otherwise system will
        # assume relative to datastore root and that is almost certainly
        # never the right default here).
        path_uri = ResourcePath(path, root=prefix, forceAbsolute=True)

        refs_by_file[path_uri].append(ref)
        n_dataset_refs += 1

    datasets = [
        FileDataset(
            path=file_uri,
            refs=refs,
            formatter=formatter,
        )
        for file_uri, refs in refs_by_file.items()
    ]

    log.info("Ingesting %d dataset ref(s) from %d file(s)", n_dataset_refs, len(datasets))

    return datasets


def parse_data_id_tuple(data_ids: Tuple[str, ...], universe: DimensionUniverse) -> Dict[str, Any]:
    # Convert any additional k=v strings in the dataId tuple to dict
    # form.
    data_id: Dict[str, Any] = {}
    for id_str in data_ids:
        dimension_str, value = id_str.split("=")

        try:
            dimension = universe.getStaticDimensions()[dimension_str]
        except KeyError:
            raise ValueError(f"DataID dimension '{dimension_str}' is not known to this universe.") from None

        # Cast the value to the right python type (since they will be
        # strings at this point).
        value = dimension.primaryKey.getPythonType()(value)

        data_id[dimension_str] = value
    return data_id

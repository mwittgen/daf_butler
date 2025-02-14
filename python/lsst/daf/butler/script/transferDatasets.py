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

__all__ = ("transferDatasets",)

import logging
from typing import Tuple

from .._butler import Butler
from ..registry.queries import DatasetQueryResults

log = logging.getLogger(__name__)


def transferDatasets(
    source: str,
    dest: str,
    dataset_type: Tuple[str, ...],
    collections: Tuple[str, ...],
    where: str,
    find_first: bool,
    transfer: str,
    register_dataset_types: bool,
    transfer_dimensions: bool = True,
) -> int:
    """Transfer datasets from run in source to dest.

    Parameters
    ----------
    source : `str`
        URI string of the source Butler repo.
    dest : `str`
        URI string of the destination Butler repo.
    dataset_type : `tuple` of `str`
        Dataset type names. An empty tuple implies all dataset types.
    collections : `tuple` of `str`
        Names of collection globs to match. An empty tuple implies all
        collections.
    where : `str`
        Query modification string.
    find_first : `bool`
        Whether only the first match should be used.
    transfer : `str`
        Transfer mode to use when placing artifacts in the destination.
    register_dataset_types : `bool`
        Indicate whether missing dataset types should be registered.
    transfer_dimensions : `bool`
        Indicate whether dimensions should be transferred along with
        datasets. It can be more efficient to disable this if it is known
        that all dimensions exist.
    """
    source_butler = Butler(source, writeable=False)
    dest_butler = Butler(dest, writeable=True)

    dataset_type_expr = ... if not dataset_type else dataset_type
    collections_expr = ... if not collections else collections

    source_refs = source_butler.registry.queryDatasets(
        datasetType=dataset_type_expr, collections=collections_expr, where=where, findFirst=find_first
    )

    # Might need expanded results if datastore records have to be derived.
    # Not all registries return the same form for results.
    if isinstance(source_refs, DatasetQueryResults):
        source_refs = source_refs.expanded()

    # Place results in a set to remove duplicates
    source_refs_set = set(source_refs)

    transferred = dest_butler.transfer_from(
        source_butler,
        source_refs_set,
        transfer=transfer,
        register_dataset_types=register_dataset_types,
        transfer_dimensions=transfer_dimensions,
    )
    return len(transferred)

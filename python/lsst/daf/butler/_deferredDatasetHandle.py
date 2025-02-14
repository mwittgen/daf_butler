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
Module containing classes used with deferring dataset loading
"""
from __future__ import annotations

__all__ = ("DeferredDatasetHandle",)

import dataclasses
from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:
    from ._limited_butler import LimitedButler
    from .core import DataCoordinate, DatasetRef


@dataclasses.dataclass(frozen=True)
class DeferredDatasetHandle:
    """Proxy class that provides deferred loading of datasets from a butler."""

    def get(
        self, *, component: Optional[str] = None, parameters: Optional[dict] = None, **kwargs: dict
    ) -> Any:
        """Retrieves the dataset pointed to by this handle

        This handle may be used multiple times, possibly with different
        parameters.

        Parameters
        ----------
        component : `str` or None
            If the deferred object is a component dataset type, this parameter
            may specify the name of the component to use in the get operation.
        parameters : `dict` or None
            The parameters argument will be passed to the butler get method.
            It defaults to None. If the value is not None,  this dict will
            be merged with the parameters dict used to construct the
            `DeferredDatasetHandle` class.
        **kwargs
            This argument is deprecated and only exists to support legacy
            gen2 butler code during migration. It is completely ignored
            and will be removed in the future.

        Returns
        -------
        return : `Object`
            The dataset pointed to by this handle
        """
        if self.parameters is not None:
            mergedParameters = self.parameters.copy()
            if parameters is not None:
                mergedParameters.update(parameters)
        elif parameters is not None:
            mergedParameters = parameters
        else:
            mergedParameters = {}

        ref = self.ref.makeComponentRef(component) if component is not None else self.ref
        return self.butler.getDirect(ref, parameters=mergedParameters)

    @property
    def dataId(self) -> DataCoordinate:
        """The full data ID associated with the dataset
        (`DataCoordinate`).

        Guaranteed to contain records.
        """
        return self.ref.dataId

    butler: LimitedButler
    """The butler that will be used to fetch the dataset (`LimitedButler`).
    """

    ref: DatasetRef
    """Reference to the dataset (`DatasetRef`).
    """

    parameters: Optional[dict]
    """Optional parameters that may be used to specify a subset of the dataset
    to be loaded (`dict` or `None`).
    """

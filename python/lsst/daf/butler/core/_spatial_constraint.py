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

__all__ = (
    "SpatialConstraint",
)

from typing import TYPE_CHECKING

import lsst.sphgeom

if TYPE_CHECKING:
    from .dimensions import SkyPixDimension


class SpatialConstraint:
    """A wrapper around `lsst.sphgeom.Region` that remembers skypix envelope
    ranges obtained from it.

    Parameters
    ----------
    arg : `lsst.sphgeom.Region` or `SpatialConstraint`
        A spatial region or another `SpatialConstraint`.  Regions are copied
        and then never modified.
    """

    def __init__(self, arg: lsst.sphgeom.Region | SpatialConstraint):
        self._ranges: dict[SkyPixDimension, lsst.sphgeom.RangeSet] = {}
        if isinstance(arg, SpatialConstraint):
            self.region: lsst.sphgeom.Region = arg.region
            self._ranges.update(arg._ranges)
        else:
            self.region = arg.clone()

    @classmethod
    def make_empty(cls) -> SpatialConstraint:
        """Make a constraint that includes no points.

        Returns
        -------
        constraint : `SpatialConstraint`
            New constraint object.
        """
        return cls(lsst.sphgeom.Box.empty())

    @classmethod
    def make_full(cls) -> SpatialConstraint:
        """Make a constraint that includes the entire sphere.

        Returns
        -------
        constraint : `SpatialConstraint`
            New constraint object.
        """
        return cls(lsst.sphgeom.Box.full())

    def __bool__(self) -> bool:
        return not self.region.getBoundingBox().isEmpty()

    def ranges(self, dimension: SkyPixDimension) -> lsst.sphgeom.RangeSet:
        """Return pixel ranges that include (but may extend beyond) the region.

        Parameters
        ----------
        dimensions : `SkyPixDimension`
            Dimension that defines the pixelization for the ranges.

        Returns
        -------
        ranges : `lsst.sphgeom.RangeSet`
            Set of pixel ranges that include the canonical area of
            the instance.

        Notes
        -----
        This is computed on first use for each dimension and then cached.
        """
        if (ranges := self._ranges.get(dimension)) is None:
            ranges = dimension.pixelization.envelope(self.region)
            self._ranges[dimension] = ranges
        return ranges

    def union(*args: lsst.sphgeom.Region | SpatialConstraint) -> SpatialConstraint:
        """Compute the union of zero or more timespans or temporal constraints.

        Parameters
        ----------
        *args : `lsst.sphgeom.Region` or `SpatialConstraint`
            Operands.

        Returns
        -------
        union : `SpatialConstraint`
            The union of the operands; empty if there are no operands.
        """
        if not args:
            return SpatialConstraint.make_empty()
        first, *rest = args
        first = SpatialConstraint(first)
        region: lsst.sphgeom.Region = first.region
        # We'll only actually union args that aren't contained by the one we've
        # built so far; this is an important optimization for unioning a bunch
        # of full constraints, as will be common.
        kept_args: list[SpatialConstraint] = [first]
        # We'll propagate cached ranges only for skypix dimension where all
        # kept args have cached ranges.
        skypix_in_kept_args = set(first._ranges.keys())
        for arg in rest:
            arg = SpatialConstraint(arg)
            relationship = region.relate(arg.region)
            if relationship & lsst.sphgeom.CONTAINS:
                # This arg is contained by the args so far, so it does not
                # actually participate.
                continue
            elif relationship & lsst.sphgeom.WITHIN:
                # This arg contains all previous args, so they don't
                # participate.
                region = arg.region
                kept_args = [arg]
                skypix_in_kept_args = set(arg._ranges.keys())
            else:
                region = lsst.sphgeom.UnionRegion(region, arg.region)
                skypix_in_kept_args.intersection_update(arg._ranges.keys())
                kept_args.append(arg)
        result = SpatialConstraint(region)
        if not skypix_in_kept_args:
            return result
        first, *rest = kept_args
        ranges = {skypix: lsst.sphgeom.RangeSet(first._ranges[skypix]) for skypix in skypix_in_kept_args}
        for arg in rest:
            for skypix in skypix_in_kept_args:
                ranges[skypix] |= arg._ranges[skypix]
        result._ranges.update(ranges)
        return result

    def intersection(*args: lsst.sphgeom.Region | SpatialConstraint) -> SpatialConstraint:
        """Compute the intersection of zero or more timespans or temporal
        constraints.

        Parameters
        ----------
        *args : `lsst.sphgeom.Region` or `SpatialConstraint`
            Operands.

        Returns
        -------
        union : `SpatialConstraint`
            The intersection of the operands; full if there are no operands.
        """
        if not args:
            return SpatialConstraint.make_full()
        first, *rest = args
        first = SpatialConstraint(first)
        region = first.region
        # We'll only actually intersect args that aren't both contained and
        # within the one we've built so far; this is an important optimization
        # for intersecting a bunch of full constraints, as will be common.
        kept_args: list[SpatialConstraint] = [first]
        # We'll propagate cached ranges only for skypix dimension where all
        # kept args have cached ranges.
        skypix_in_kept_args = set(first._ranges.keys())
        for arg in rest:
            arg = SpatialConstraint(arg)
            relationship = region.relate(arg.region)
            if relationship & lsst.sphgeom.CONTAINS and relationship & lsst.sphgeom.WITHIN:
                # This arg is equal to the args so far, so it does not need to
                # participate.
                continue
            elif relationship & lsst.sphgeom.DISJOINT:
                # This arg is disjoint with the args so far, so the result
                # is empty.
                return SpatialConstraint.make_empty()
            else:
                region = lsst.sphgeom.IntersectionRegion(region, arg.region)
                skypix_in_kept_args.intersection_update(arg._ranges.keys())
                kept_args.append(arg)
        result = SpatialConstraint(region)
        if not skypix_in_kept_args:
            return result
        first, *rest = kept_args
        ranges = {skypix: lsst.sphgeom.RangeSet(first._ranges[skypix]) for skypix in skypix_in_kept_args}
        for arg in rest:
            for skypix in skypix_in_kept_args:
                ranges[skypix] &= arg._ranges[skypix]
        result._ranges.update(ranges)
        return result

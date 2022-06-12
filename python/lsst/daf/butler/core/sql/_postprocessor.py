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

__all__ = ("Postprocessor",)

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, AbstractSet, Any, Iterable, Optional

from ._column_tags import ColumnTag, DimensionRecordColumnTag, ResultRow, ResultTag

if TYPE_CHECKING:
    from .._spatial_regions import SpatialConstraint
    from ..dimensions import DimensionRecord


class Postprocessor(ABC):
    @property
    @abstractmethod
    def columns_required(self) -> AbstractSet[ResultTag]:
        raise NotImplementedError()

    @property
    @abstractmethod
    def columns_provided(self) -> AbstractSet[ResultTag]:
        raise NotImplementedError()

    @property
    @abstractmethod
    def row_multiplier(self) -> float:
        raise NotImplementedError()

    @abstractmethod
    def apply(self, row: ResultRow) -> Optional[ResultRow]:
        raise NotImplementedError()

    def __eq__(self, rhs: Any) -> bool:
        return (self.columns_required == rhs.columns_required) and (
            self.columns_provided == rhs.columns_provided
        )

    def __hash__(self) -> int:
        return hash((frozenset(self.columns_required), frozenset(self.columns_provided)))

    @staticmethod
    def from_spatial_join(a: ColumnTag, b: ColumnTag) -> Postprocessor:
        return _SpatialJoinPostprocessor(a, b)

    @staticmethod
    def from_spatial_constraint(constraint: SpatialConstraint, tag: ColumnTag) -> Postprocessor:
        return _SpatialConstraintPostprocessor(constraint, tag)

    @staticmethod
    def make_dimension_column_extractor(
        record_cls: type[DimensionRecord], columns: AbstractSet[str]
    ) -> Postprocessor:
        return _ExtractDimensionRecordPostprocessor(record_cls, columns)

    def _tiebreaker_sort_key(self) -> float:
        return self.row_multiplier

    @staticmethod
    def sort_and_check(
        postprocessors: Iterable[Postprocessor],
        columns_provided: AbstractSet[ColumnTag],
    ) -> tuple[list[Postprocessor], AbstractSet[Postprocessor], AbstractSet[ResultTag]]:
        todo = set(postprocessors)
        done: list[Postprocessor] = []
        full_columns_provided: set[ResultTag] = set(columns_provided)
        while todo:
            candidates_to_include: set[Postprocessor] = set()
            columns_to_provide: set[ResultTag] = set()
            for candidate in todo:
                if full_columns_provided.issuperset(candidate.columns_required):
                    candidates_to_include.add(candidate)
                    columns_to_provide.update(candidate.columns_provided)
            if not candidates_to_include:
                missing: set[ResultTag] = set()
                for failed in todo:
                    missing.update(failed.columns_required)
                    missing.difference_update(failed.columns_provided)
                return done, todo, missing
            todo.difference_update(candidates_to_include)
            done.extend(sorted(candidates_to_include, key=Postprocessor._tiebreaker_sort_key))
            full_columns_provided.update(columns_to_provide)
        return done, frozenset(), frozenset()

    @staticmethod
    def sort_and_assert(
        postprocessors: Iterable[Postprocessor],
        columns_provided: AbstractSet[ColumnTag],
    ) -> list[Postprocessor]:
        result, _, missing = Postprocessor.sort_and_check(postprocessors, columns_provided)
        assert not missing, "Preprocessors should have been checked at a lower level."
        return result


class _SpatialJoinPostprocessor(Postprocessor):
    def __init__(self, a: ColumnTag, b: ColumnTag):
        self._a = a
        self._b = b

    __slots__ = ("_a", "_b")

    @property
    def columns_required(self) -> AbstractSet[ColumnTag]:
        return {self._a, self._b}

    @property
    def columns_provided(self) -> AbstractSet[ColumnTag]:
        return frozenset()

    @property
    def row_multiplier(self) -> float:
        return 0.5

    def apply(self, row: ResultRow) -> Optional[ResultRow]:
        return None if row[self._a].isDisjointFrom(row[self._b]) else row


class _SpatialConstraintPostprocessor(Postprocessor):
    def __init__(self, constraint: SpatialConstraint, tag: ColumnTag):
        self._constraint = constraint
        self._tag = tag

    __slots__ = ("_constraint", "_constraint")

    @property
    def columns_required(self) -> AbstractSet[ColumnTag]:
        return {self._tag}

    @property
    def columns_provided(self) -> AbstractSet[ColumnTag]:
        return frozenset()

    @property
    def row_multiplier(self) -> float:
        return 0.5

    def apply(self, row: ResultRow) -> Optional[ResultRow]:
        return None if self._constraint.region.isDisjointWith(row[self._tag]) else row


class _ExtractDimensionRecordPostprocessor(Postprocessor):
    def __init__(self, record_cls: type[DimensionRecord], columns: AbstractSet[str]):
        self._record_cls = record_cls
        self._column_tags = frozenset(DimensionRecordColumnTag.generate(record_cls.definition.name, columns))

    __slots__ = ("_record_cls", "_column_tags")

    @property
    def columns_provided(self) -> AbstractSet[ResultTag]:
        return self._column_tags

    @property
    def columns_required(self) -> AbstractSet[ResultTag]:
        return {self._record_cls}

    @property
    def row_multiplier(self) -> float:
        return 1.0

    def apply(self, row: ResultRow) -> ResultRow:
        record = row[self._record_cls]
        for tag in self._column_tags:
            row[tag] = getattr(record, tag.column)
        return row

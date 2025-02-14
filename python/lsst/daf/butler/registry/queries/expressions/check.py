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
    "CheckVisitor",
    "InspectionVisitor",
    "InspectionSummary",
)

import dataclasses
from typing import TYPE_CHECKING, Any, List, Mapping, Optional, Sequence, Set, Tuple

from ....core import (
    DataCoordinate,
    Dimension,
    DimensionElement,
    DimensionGraph,
    DimensionUniverse,
    GovernorDimension,
    NamedKeyDict,
    NamedValueSet,
)
from ..._exceptions import UserExpressionError
from ...summaries import GovernorDimensionRestriction
from .categorize import ExpressionConstant, categorizeConstant, categorizeElementId
from .normalForm import NormalForm, NormalFormVisitor
from .parser import Node, TreeVisitor

if TYPE_CHECKING:
    import astropy.time


@dataclasses.dataclass
class InspectionSummary:
    """Base class for objects used by `CheckVisitor` and `InspectionVisitor`
    to gather information about a parsed expression.
    """

    def update(self, other: InspectionSummary) -> None:
        """Update ``self`` with all dimensions and columns from ``other``.

        Parameters
        ----------
        other : `InspectionSummary`
            The other summary object.
        """
        self.dimensions.update(other.dimensions)
        for element, columns in other.columns.items():
            self.columns.setdefault(element, set()).update(columns)
        self.hasIngestDate = self.hasIngestDate or other.hasIngestDate

    dimensions: NamedValueSet[Dimension] = dataclasses.field(default_factory=NamedValueSet)
    """Dimensions whose primary keys or dependencies were referenced anywhere
    in this branch (`NamedValueSet` [ `Dimension` ]).
    """

    columns: NamedKeyDict[DimensionElement, Set[str]] = dataclasses.field(default_factory=NamedKeyDict)
    """Dimension element tables whose columns were referenced anywhere in this
    branch (`NamedKeyDict` [ `DimensionElement`, `set` [ `str` ] ]).
    """

    hasIngestDate: bool = False
    """Whether this expression includes the special dataset ingest date
    identifier (`bool`).
    """


@dataclasses.dataclass
class TreeSummary(InspectionSummary):
    """Result object used by `InspectionVisitor` to gather information about
    a parsed expression.

    Notes
    -----
    TreeSummary adds attributes that allow dimension equivalence expressions
    (e.g. "tract=4") to be recognized when they appear in simple contexts
    (surrounded only by ANDs and ORs).  When `InspectionVisitor` is used on its
    own (i.e. when ``check=False`` in the query code), these don't do anything,
    but they don't cost much, either.  They are used by `CheckVisitor` when it
    delegates to `InspectionVisitor` to see what governor dimension values are
    set in a branch of the normal-form expression.
    """

    def merge(self, other: TreeSummary, isEq: bool = False) -> TreeSummary:
        """Merge ``other`` into ``self``, making ``self`` a summary of both
        expression tree branches.

        Parameters
        ----------
        other : `TreeSummary`
            The other summary object.
        isEq : `bool`, optional
            If `True` (`False` is default), these summaries are being combined
            via the equality operator.

        Returns
        -------
        self : `TreeSummary`
            The merged summary (updated in-place).
        """
        self.update(other)
        if isEq and self.isDataIdKeyOnly() and other.isDataIdValueOnly():
            self.dataIdValue = other.dataIdValue
        elif isEq and self.isDataIdValueOnly() and other.isDataIdKeyOnly():
            self.dataIdKey = other.dataIdKey
        else:
            self.dataIdKey = None
            self.dataIdValue = None
        return self

    def isDataIdKeyOnly(self) -> bool:
        """Test whether this branch is _just_ a data ID key identifier."""
        return self.dataIdKey is not None and self.dataIdValue is None

    def isDataIdValueOnly(self) -> bool:
        """Test whether this branch is _just_ a literal value that may be
        used as the value in a data ID key-value pair.
        """
        return self.dataIdKey is None and self.dataIdValue is not None

    dataIdKey: Optional[Dimension] = None
    """A `Dimension` that is (if `dataIdValue` is not `None`) or may be
    (if `dataIdValue` is `None`) fully identified by a literal value in this
    branch.
    """

    dataIdValue: Optional[str] = None
    """A literal value that constrains (if `dataIdKey` is not `None`) or may
    constrain (if `dataIdKey` is `None`) a dimension in this branch.

    This is always a `str` or `None`, but it may need to be coerced to `int`
    to reflect the actual user intent.
    """


class InspectionVisitor(TreeVisitor[TreeSummary]):
    """Implements TreeVisitor to identify dimension elements that need
    to be included in a query, prior to actually constructing a SQLAlchemy
    WHERE clause from it.

    Parameters
    ----------
    universe : `DimensionUniverse`
        All known dimensions.
    bind : `Mapping` [ `str`, `object` ]
        Mapping containing literal values that should be injected into the
        query expression, keyed by the identifiers they replace.
    """

    def __init__(self, universe: DimensionUniverse, bind: Mapping[str, Any]):
        self.universe = universe
        self.bind = bind

    def visitNumericLiteral(self, value: str, node: Node) -> TreeSummary:
        # Docstring inherited from TreeVisitor.visitNumericLiteral
        return TreeSummary(dataIdValue=value)

    def visitStringLiteral(self, value: str, node: Node) -> TreeSummary:
        # Docstring inherited from TreeVisitor.visitStringLiteral
        return TreeSummary(dataIdValue=value)

    def visitTimeLiteral(self, value: astropy.time.Time, node: Node) -> TreeSummary:
        # Docstring inherited from TreeVisitor.visitTimeLiteral
        return TreeSummary()

    def visitIdentifier(self, name: str, node: Node) -> TreeSummary:
        # Docstring inherited from TreeVisitor.visitIdentifier
        if name in self.bind:
            return TreeSummary(dataIdValue=self.bind[name])
        constant = categorizeConstant(name)
        if constant is ExpressionConstant.INGEST_DATE:
            return TreeSummary(hasIngestDate=True)
        elif constant is ExpressionConstant.NULL:
            return TreeSummary()
        assert constant is None, "Enum variant conditionals should be exhaustive."
        element, column = categorizeElementId(self.universe, name)
        if column is None:
            assert isinstance(element, Dimension)
            return TreeSummary(
                dimensions=NamedValueSet(element.graph.dimensions),
                dataIdKey=element,
            )
        else:
            return TreeSummary(
                dimensions=NamedValueSet(element.graph.dimensions), columns=NamedKeyDict({element: {column}})
            )

    def visitUnaryOp(self, operator: str, operand: TreeSummary, node: Node) -> TreeSummary:
        # Docstring inherited from TreeVisitor.visitUnaryOp
        return operand

    def visitBinaryOp(self, operator: str, lhs: TreeSummary, rhs: TreeSummary, node: Node) -> TreeSummary:
        # Docstring inherited from TreeVisitor.visitBinaryOp
        return lhs.merge(rhs, isEq=(operator == "="))

    def visitIsIn(self, lhs: TreeSummary, values: List[TreeSummary], not_in: bool, node: Node) -> TreeSummary:
        # Docstring inherited from TreeVisitor.visitIsIn
        for v in values:
            lhs.merge(v)
        return lhs

    def visitParens(self, expression: TreeSummary, node: Node) -> TreeSummary:
        # Docstring inherited from TreeVisitor.visitParens
        return expression

    def visitTupleNode(self, items: Tuple[TreeSummary, ...], node: Node) -> TreeSummary:
        # Docstring inherited from base class
        result = TreeSummary()
        for i in items:
            result.merge(i)
        return result

    def visitRangeLiteral(self, start: int, stop: int, stride: Optional[int], node: Node) -> TreeSummary:
        # Docstring inherited from TreeVisitor.visitRangeLiteral
        return TreeSummary()

    def visitPointNode(self, ra: TreeSummary, dec: TreeSummary, node: Node) -> TreeSummary:
        # Docstring inherited from base class
        return TreeSummary()


@dataclasses.dataclass
class InnerSummary(InspectionSummary):
    """Result object used by `CheckVisitor` to gather referenced dimensions
    and tables from an inner group of AND'd together expression branches, and
    check them for consistency and completeness.
    """

    governors: NamedKeyDict[GovernorDimension, str] = dataclasses.field(default_factory=NamedKeyDict)
    """Mapping containing the values of all governor dimensions that are
    equated with literal values in this expression branch.
    """

    defaultsNeeded: NamedValueSet[GovernorDimension] = dataclasses.field(default_factory=NamedValueSet)
    """Governor dimensions whose values are needed by the query, not provided
    in the query itself, and present in the default data ID.

    These should be added to the query's data ID when finalizing the WHERE
    clause.
    """


@dataclasses.dataclass
class OuterSummary(InspectionSummary):
    """Result object used by `CheckVisitor` to gather referenced dimensions,
    tables, and governor dimension values from the entire expression.
    """

    governors: GovernorDimensionRestriction = dataclasses.field(
        default_factory=GovernorDimensionRestriction.makeFull
    )
    """Mapping containing all values that appear in this expression for
    governor dimension relevant to the query.

    Governor dimensions that are absent from this dict are not constrained by
    this expression.
    """

    defaultsNeeded: NamedValueSet[GovernorDimension] = dataclasses.field(default_factory=NamedValueSet)
    """Governor dimensions whose values are needed by the query, not provided
    in the query itself, and present in the default data ID.

    These should be added to the query's data ID when finalizing the WHERE
    clause.
    """


class CheckVisitor(NormalFormVisitor[TreeSummary, InnerSummary, OuterSummary]):
    """An implementation of `NormalFormVisitor` that identifies the dimensions
    and tables that need to be included in a query while performing some checks
    for completeness and consistency.

    Parameters
    ----------
    dataId : `DataCoordinate`
        Dimension values that are fully known in advance.
    graph : `DimensionGraph`
        The dimensions the query would include in the absence of this
        expression.
    bind : `Mapping` [ `str`, `object` ]
        Mapping containing literal values that should be injected into the
        query expression, keyed by the identifiers they replace.
    defaults : `DataCoordinate`
        A data ID containing default for governor dimensions.
    """

    def __init__(
        self,
        dataId: DataCoordinate,
        graph: DimensionGraph,
        bind: Mapping[str, Any],
        defaults: DataCoordinate,
    ):
        self.dataId = dataId
        self.graph = graph
        self.defaults = defaults
        self._branchVisitor = InspectionVisitor(dataId.universe, bind)

    def visitBranch(self, node: Node) -> TreeSummary:
        # Docstring inherited from NormalFormVisitor.
        return node.visit(self._branchVisitor)

    def visitInner(self, branches: Sequence[TreeSummary], form: NormalForm) -> InnerSummary:
        # Docstring inherited from NormalFormVisitor.
        # Disjunctive normal form means inner branches are AND'd together...
        assert form is NormalForm.DISJUNCTIVE
        # ...and that means each branch we iterate over together below
        # constrains the others, and they all need to be consistent.  Moreover,
        # because outer branches are OR'd together, we also know that if
        # something is missing from one of these branches (like a governor
        # dimension value like the instrument or skymap needed to interpret a
        # visit or tract number), it really is missing, because there's no way
        # some other inner branch can constraint it.
        #
        # That is, except the data ID the visitor was passed at construction;
        # that's AND'd to the entire expression later, and thus it affects all
        # branches.  To take care of that, we add any governor values it
        # contains to the summary in advance.
        summary = InnerSummary()
        summary.governors.update((k, self.dataId[k]) for k in self.dataId.graph.governors)  # type: ignore
        # Finally, we loop over those branches.
        for branch in branches:
            # Update the sets of dimensions and columns we've seen anywhere in
            # the expression in any context.
            summary.update(branch)
            # Test whether this branch has a form like '<dimension>=<value'
            # (or equivalent; categorizeIdentifier is smart enough to see that
            # e.g. 'detector.id=4' is equivalent to 'detector=4').
            # If so, and it's a governor dimension, remember that we've
            # constrained it on this branch, and make sure it's consistent
            # with any other constraints on any other branches its AND'd with.
            if isinstance(branch.dataIdKey, GovernorDimension) and branch.dataIdValue is not None:
                governor = branch.dataIdKey
                value = summary.governors.setdefault(governor, branch.dataIdValue)
                if value != branch.dataIdValue:
                    # Expression says something like "instrument='HSC' AND
                    # instrument='DECam'", or data ID has one and expression
                    # has the other.
                    if governor in self.dataId:
                        raise UserExpressionError(
                            f"Conflict between expression containing {governor.name}={branch.dataIdValue!r} "
                            f"and data ID with {governor.name}={value!r}."
                        )
                    else:
                        raise UserExpressionError(
                            f"Conflicting literal values for {governor.name} in expression: "
                            f"{value!r} != {branch.dataIdValue!r}."
                        )
        # Now that we know which governor values we've constrained, see if any
        # are missing, i.e. if the expression contains something like "visit=X"
        # without saying what instrument that visit corresponds to.  This rules
        # out a lot of accidents, but it also rules out possibly-legitimate
        # multi-instrument queries like "visit.seeing < 0.7".  But it's not
        # unreasonable to ask the user to be explicit about the instruments
        # they want to consider to work around this restriction, and that's
        # what we do.  Note that if someone does write an expression like
        #
        #  (instrument='HSC' OR instrument='DECam') AND visit.seeing < 0.7
        #
        # then in disjunctive normal form that will become
        #
        #  (instrument='HSC' AND visit.seeing < 0.7)
        #    OR (instrument='DECam' AND visit.seeing < 0.7)
        #
        # i.e. each instrument will get its own outer branch and the logic here
        # still works (that sort of thing is why we convert to normal form,
        # after all).
        governorsNeededInBranch: NamedValueSet[GovernorDimension] = NamedValueSet()
        for dimension in summary.dimensions:
            governorsNeededInBranch.update(dimension.graph.governors)
        if not governorsNeededInBranch.issubset(summary.governors.keys()):
            missing = NamedValueSet(governorsNeededInBranch - summary.governors.keys())
            if missing <= self.defaults.keys():
                summary.defaultsNeeded.update(missing)
            else:
                raise UserExpressionError(
                    f"No value(s) for governor dimensions {missing - self.defaults.keys()} in expression "
                    "that references dependent dimensions. 'Governor' dimensions must always be specified "
                    "completely in either the query expression (via simple 'name=<value>' terms, not 'IN' "
                    "terms) or in a data ID passed to the query method."
                )
        return summary

    def visitOuter(self, branches: Sequence[InnerSummary], form: NormalForm) -> OuterSummary:
        # Docstring inherited from NormalFormVisitor.
        # Disjunctive normal form means outer branches are OR'd together.
        assert form is NormalForm.DISJUNCTIVE
        # Iterate over branches in first pass to gather all dimensions and
        # columns referenced.  This aggregation is for the full query, so we
        # don't care whether things are joined by AND or OR (or + or -, etc).
        summary = OuterSummary()
        if branches:
            # To make an OR of branch constraints start with empty selection.
            summary.governors = GovernorDimensionRestriction.makeEmpty(self.graph.universe)
            for branch in branches:
                summary.update(branch)
                summary.governors = summary.governors.union(branch.governors)
                summary.defaultsNeeded.update(branch.defaultsNeeded)
        # See if we've referenced any dimensions that weren't in the original
        # query graph; if so, we update that to include them.  This is what
        # lets a user say "tract=X" on the command line (well, "skymap=Y AND
        # tract=X" - logic in visitInner checks for that) when running a task
        # like ISR that has nothing to do with skymaps.
        if not summary.dimensions.issubset(self.graph.dimensions):
            self.graph = DimensionGraph(
                self.graph.universe,
                dimensions=(summary.dimensions | self.graph.dimensions),
            )
        for governor, values in summary.governors.items():
            if governor in summary.defaultsNeeded:
                # One branch contained an explicit value for this dimension
                # while another needed to refer to the default data ID.
                # Even if these refer to the same value, that inconsistency
                # probably indicates user error.
                raise UserExpressionError(
                    f"Governor dimension {governor.name} is explicitly "
                    f"constrained to {values} in one or more branches of "
                    "this query where expression, but is left to default "
                    f"to {self.defaults[governor]!r} in another branch.  "
                    "Defaults and explicit constraints cannot be mixed."
                )
        # If any default data ID values were needed, update self.dataId with
        # them, and then update the governor restriction with them.
        if summary.defaultsNeeded:
            defaultsNeededGraph = DimensionGraph(self.graph.universe, summary.defaultsNeeded)
            self.dataId = self.dataId.union(self.defaults.subset(defaultsNeededGraph))
            assert self.dataId.hasRecords(), (
                "Should be a union of two data IDs with records, "
                "in which one only adds governor dimension values."
            )
            summary.governors.intersection_update(
                # We know the value for a governor dimension is always a str,
                # and that's all self.defaults should contain, but MyPy doesn't
                # know that.
                {dimension: self.defaults[dimension] for dimension in summary.defaultsNeeded}  # type: ignore
            )
        return summary

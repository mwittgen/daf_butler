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
    "stripIfNotNone",
    "transactional",
)

import fnmatch
import functools
import logging
import re
from typing import TYPE_CHECKING, Any, Callable, List, Optional, Pattern, TypeVar, Union

from lsst.utils.iteration import ensure_iterable

if TYPE_CHECKING:
    from ..registry.wildcards import Ellipsis, EllipsisType


_LOG = logging.getLogger(__name__)


F = TypeVar("F", bound=Callable)


def transactional(func: F) -> F:
    """Decorate a method and makes it transactional.

    This depends on the class also defining a `transaction` method
    that takes no arguments and acts as a context manager.
    """

    @functools.wraps(func)
    def inner(self: Any, *args: Any, **kwargs: Any) -> Any:
        with self.transaction():
            return func(self, *args, **kwargs)

    return inner  # type: ignore


def stripIfNotNone(s: Optional[str]) -> Optional[str]:
    """Strip leading and trailing whitespace if the given object is not None.

    Parameters
    ----------
    s : `str`, optional
        Input string.

    Returns
    -------
    r : `str` or `None`
        A string with leading and trailing whitespace stripped if `s` is not
        `None`, or `None` if `s` is `None`.
    """
    if s is not None:
        s = s.strip()
    return s


def globToRegex(
    expressions: Union[str, EllipsisType, None, List[str]]
) -> Union[List[Union[str, Pattern]], EllipsisType]:
    """Translate glob-style search terms to regex.

    If a stand-alone '``*``' is found in ``expressions``, or expressions is
    empty or `None`, then the special value ``...`` will be returned,
    indicating that any string will match.

    Parameters
    ----------
    expressions : `str` or `list` [`str`]
        A list of glob-style pattern strings to convert.

    Returns
    -------
    expressions : `list` [`str` or `re.Pattern`] or ``...``
        A list of regex Patterns or simple strings. Returns ``...`` if
        the provided expressions would match everything.
    """
    if expressions is Ellipsis or expressions is None:
        return Ellipsis
    expressions = list(ensure_iterable(expressions))
    if not expressions or "*" in expressions:
        return Ellipsis

    # List of special glob characters supported by fnmatch.
    # See: https://docs.python.org/3/library/fnmatch.html
    # The complication is that "[" on its own is not a glob
    # unless there is a match "]".
    magic = re.compile(r"[\*\?]|\[.*\]|\[!.*\]")

    # Try not to convert simple string to a regex.
    results: List[Union[str, Pattern]] = []
    for e in expressions:
        res: Union[str, Pattern]
        if magic.search(e):
            res = re.compile(fnmatch.translate(e))
        else:
            res = e
        results.append(res)
    return results

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


__all__ = (
    "CollectionTypeCallback",
    "collection_type_option",
    "collections_option",
    "components_option",
    "config_option",
    "config_file_option",
    "confirm_option",
    "dataset_type_option",
    "datasets_option",
    "log_level_option",
    "long_log_option",
    "log_file_option",
    "log_label_option",
    "log_tty_option",
    "options_file_option",
    "processes_option",
    "regex_option",
    "register_dataset_types_option",
    "run_option",
    "transfer_option",
    "verbose_option",
    "where_option",
    "order_by_option",
    "limit_option",
    "offset_option",
)


from functools import partial

import click
from lsst.daf.butler.registry import CollectionType

from ..cliLog import CliLog
from ..utils import MWOptionDecorator, MWPath, split_commas, split_kv, unwrap, yaml_presets


class CollectionTypeCallback:

    collectionTypes = tuple(collectionType.name for collectionType in CollectionType.all())

    @staticmethod
    def makeCollectionTypes(context, param, value):
        if not value:
            # Click seems to demand that the default be an empty tuple, rather
            # than a sentinal like None.  The behavior that we want is that
            # not passing this option at all passes all collection types, while
            # passing it uses only the passed collection types.  That works
            # fine for now, since there's no command-line option to subtract
            # collection types, and hence the only way to get an empty tuple
            # is as the default.
            return tuple(CollectionType.all())

        return tuple(CollectionType.from_name(item) for item in split_commas(context, param, value))


collection_type_option = MWOptionDecorator(
    "--collection-type",
    callback=CollectionTypeCallback.makeCollectionTypes,
    multiple=True,
    help="If provided, only list collections of this type.",
    type=click.Choice(choices=CollectionTypeCallback.collectionTypes, case_sensitive=False),
)


collections_option = MWOptionDecorator(
    "--collections",
    help=unwrap(
        """One or more expressions that fully or partially identify
                                                   the collections to search for datasets. If not provided all
                                                   datasets are returned."""
    ),
    multiple=True,
    callback=split_commas,
)


components_option = MWOptionDecorator(
    "--components/--no-components",
    default=None,
    help=unwrap(
        """For --components, apply all expression patterns to
                                                  component dataset type names as well. For --no-components,
                                                  never apply patterns to components. Default (where neither
                                                  is specified) is to apply patterns to components only if
                                                  their parent datasets were not matched by the expression.
                                                  Fully-specified component datasets (`str` or `DatasetType`
                                                  instances) are always included."""
    ),
)


config_option = MWOptionDecorator(
    "-c",
    "--config",
    callback=split_kv,
    help="Config override, as a key-value pair.",
    metavar="TEXT=TEXT",
    multiple=True,
)


config_file_option = MWOptionDecorator(
    "-C",
    "--config-file",
    help=unwrap(
        """Path to a pex config override to be included after the
                                       Instrument config overrides are applied."""
    ),
)


confirm_option = MWOptionDecorator(
    "--confirm/--no-confirm",
    default=True,
    help="Print expected action and a confirmation prompt before executing. Default is --confirm.",
)


dataset_type_option = MWOptionDecorator(
    "-d", "--dataset-type", callback=split_commas, help="Specific DatasetType(s) to validate.", multiple=True
)


datasets_option = MWOptionDecorator("--datasets")


logLevelChoices = ["CRITICAL", "ERROR", "WARNING", "INFO", "VERBOSE", "DEBUG", "TRACE"]
log_level_option = MWOptionDecorator(
    "--log-level",
    callback=partial(
        split_kv,
        choice=click.Choice(choices=logLevelChoices, case_sensitive=False),
        normalize=True,
        unseparated_okay=True,
        add_to_default=True,
        default_key=None,  # No separator
    ),
    help=f"The logging level. Without an explicit logger name, will only affect the default root loggers "
    f"({', '.join(CliLog.root_loggers())}). To modify the root logger use '.=LEVEL'. "
    f"Supported levels are [{'|'.join(logLevelChoices)}]",
    is_eager=True,
    metavar="LEVEL|COMPONENT=LEVEL",
    multiple=True,
)


long_log_option = MWOptionDecorator(
    "--long-log", help="Make log messages appear in long format.", is_flag=True
)

log_file_option = MWOptionDecorator(
    "--log-file",
    default=None,
    multiple=True,
    callback=split_commas,
    type=MWPath(file_okay=True, dir_okay=False, writable=True),
    help="File(s) to write log messages. If the path ends with '.json' then"
    " JSON log records will be written, else formatted text log records"
    " will be written. This file can exist and records will be appended.",
)

log_label_option = MWOptionDecorator(
    "--log-label",
    default=None,
    multiple=True,
    callback=split_kv,
    type=str,
    help="Keyword=value pairs to add to MDC of log records.",
)

log_tty_option = MWOptionDecorator(
    "--log-tty/--no-log-tty",
    default=True,
    help="Log to terminal (default). If false logging to terminal is disabled.",
)

options_file_option = MWOptionDecorator(
    "--options-file",
    "-@",
    expose_value=False,  # This option should not be forwarded
    help=unwrap(
        """URI to YAML file containing overrides
                                                    of command line options. The YAML should be organized
                                                    as a hierarchy with subcommand names at the top
                                                    level options for that subcommand below."""
    ),
    callback=yaml_presets,
)


processes_option = MWOptionDecorator(
    "-j", "--processes", default=1, help="Number of processes to use.", type=click.IntRange(min=1)
)


regex_option = MWOptionDecorator("--regex")


register_dataset_types_option = MWOptionDecorator(
    "--register-dataset-types",
    help=unwrap(
        """Register DatasetTypes that do not already
                                                              exist in the Registry."""
    ),
    is_flag=True,
)

run_option = MWOptionDecorator("--output-run", help="The name of the run datasets should be output to.")


transfer_option = MWOptionDecorator(
    "-t",
    "--transfer",
    default="auto",  # set to `None` if using `required=True`
    help="The external data transfer mode.",
    type=click.Choice(
        choices=["auto", "link", "symlink", "hardlink", "copy", "move", "relsymlink", "direct"],
        case_sensitive=False,
    ),
)


verbose_option = MWOptionDecorator("-v", "--verbose", help="Increase verbosity.", is_flag=True)


where_option = MWOptionDecorator("--where", help="A string expression similar to a SQL WHERE clause.")


order_by_option = MWOptionDecorator(
    "--order-by",
    help=unwrap(
        """One or more comma-separated names used to order records. Names can be dimension names,
                metadata names optionally prefixed by a dimension name and dot, or
                timestamp_begin/timestamp_end (with optional dimension name). To reverse ordering for a name
                prefix it with a minus sign."""
    ),
    multiple=True,
    callback=split_commas,
)


limit_option = MWOptionDecorator(
    "--limit",
    help=unwrap("Limit the number of records, by default all records are shown."),
    type=int,
    default=0,
)

offset_option = MWOptionDecorator(
    "--offset",
    help=unwrap("Skip initial number of records, only used when --limit is specified."),
    type=int,
    default=0,
)

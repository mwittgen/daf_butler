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

"""Unit tests for the daf_butler CliLog utility. Code is implemented in
daf_butler but some only runs if lsst.log.Log can be imported so these parts of
it can't be tested there because daf_butler does not directly depend on
lsst.log, and only uses it if it has been setup by another package."""

import logging
import os
import re
import subprocess
import tempfile
import unittest
from collections import namedtuple
from functools import partial
from io import StringIO
from logging import DEBUG, INFO, WARNING

import click
from lsst.daf.butler.cli.butler import cli as butlerCli
from lsst.daf.butler.cli.cliLog import CliLog
from lsst.daf.butler.cli.utils import LogCliRunner, clickResultMsg, command_test_env
from lsst.daf.butler.core.logging import ButlerLogRecords
from lsst.utils.logging import TRACE

try:
    import lsst.log as lsstLog

    lsstLog_INFO = lsstLog.INFO
    lsstLog_DEBUG = lsstLog.DEBUG
    lsstLog_WARN = lsstLog.WARN
except ModuleNotFoundError:
    lsstLog = None
    lsstLog_INFO = 0
    lsstLog_DEBUG = 0
    lsstLog_WARN = 0


@click.command()
@click.option("--expected-pyroot-level", type=int)
@click.option("--expected-pylsst-level", type=int)
@click.option("--expected-pybutler-level", type=int)
@click.option("--expected-lsstroot-level", type=int)
@click.option("--expected-lsstbutler-level", type=int)
@click.option("--expected-lsstx-level", type=int)
def command_log_settings_test(
    expected_pyroot_level,
    expected_pylsst_level,
    expected_pybutler_level,
    expected_lsstroot_level,
    expected_lsstbutler_level,
    expected_lsstx_level,
):

    LogLevel = namedtuple("LogLevel", ("expected", "actual", "name"))

    logLevels = [
        LogLevel(expected_pyroot_level, logging.getLogger().level, "pyRoot"),
        LogLevel(expected_pylsst_level, logging.getLogger("lsst").getEffectiveLevel(), "pyLsst"),
        LogLevel(
            expected_pybutler_level, logging.getLogger("lsst.daf.butler").getEffectiveLevel(), "pyButler"
        ),
        LogLevel(expected_lsstx_level, logging.getLogger("lsstx").getEffectiveLevel(), "pyLsstx"),
    ]
    if lsstLog is not None:
        logLevels.extend(
            [
                LogLevel(expected_lsstroot_level, lsstLog.getLogger("lsst").getEffectiveLevel(), "lsstRoot"),
                LogLevel(
                    expected_lsstbutler_level,
                    lsstLog.getLogger("lsst.daf.butler").getEffectiveLevel(),
                    "lsstButler",
                ),
            ]
        )
    for expected, actual, name in logLevels:
        if expected != actual:
            raise (
                click.ClickException(message=f"expected {name} level to be {expected!r}, actual:{actual!r}")
            )


class CliLogTestBase:
    """Tests log initialization, reset, and setting log levels."""

    lsstLogHandlerId = None

    def setUp(self):
        self.runner = LogCliRunner()

    def tearDown(self):
        self.lsstLogHandlerId = None

    class PythonLogger:
        """Keeps track of log level of a component and number of handlers
        attached to it at the time this object was initialized."""

        def __init__(self, component):
            self.logger = logging.getLogger(component)
            self.initialLevel = self.logger.level

    class LsstLogger:
        """Keeps track of log level for a component at the time this object was
        initialized."""

        def __init__(self, component):
            self.logger = lsstLog.getLogger(component) if lsstLog else None
            self.initialLevel = self.logger.getLevel() if lsstLog else None

    def runTest(self, cmd):
        """Test that the log context manager works with the butler cli to
        initialize the logging system according to cli inputs for the duration
        of the command execution and resets the logging system to its previous
        state or expected state when command execution finishes."""
        pyRoot = self.PythonLogger(None)
        pyButler = self.PythonLogger("lsst.daf.butler")
        pyLsstRoot = self.PythonLogger("lsst")
        lsstRoot = self.LsstLogger("")
        lsstButler = self.LsstLogger("lsst.daf.butler")

        with command_test_env(
            self.runner, "lsst.daf.butler.tests.cliLogTestBase", "command-log-settings-test"
        ):
            result = cmd()
        self.assertEqual(result.exit_code, 0, clickResultMsg(result))

        self.assertEqual(pyRoot.logger.level, logging.WARNING)
        self.assertEqual(pyLsstRoot.logger.level, logging.INFO)
        self.assertEqual(pyButler.logger.level, pyButler.initialLevel)
        if lsstLog is not None:
            self.assertEqual(lsstRoot.logger.getLevel(), lsstLog.INFO)
            # lsstLogLevel can either be the initial level, or uninitialized or
            # the defined default value.
            expectedLsstLogLevel = (
                (lsstButler.initialLevel,)
                if lsstButler.initialLevel != -1
                else (-1, CliLog.defaultLsstLogLevel)
            )
            self.assertIn(lsstButler.logger.getLevel(), expectedLsstLogLevel)

    def test_butlerCliLog(self):
        """Test that the log context manager works with the butler cli to
        initialize the logging system according to cli inputs for the duration
        of the command execution and resets the logging system to its previous
        state or expected state when command execution finishes."""

        # Run with two different log level settings.
        log_levels = (
            # --log-level / --log-level / expected pyroot, pylsst, pybutler,
            # lsstroot, lsstbutler, lsstx
            (
                "WARNING",
                "lsst.daf.butler=DEBUG",
                WARNING,
                WARNING,
                DEBUG,
                lsstLog_WARN,
                lsstLog_DEBUG,
                WARNING,
            ),
            ("DEBUG", "lsst.daf.butler=TRACE", WARNING, DEBUG, TRACE, lsstLog_DEBUG, lsstLog_DEBUG, WARNING),
            (".=DEBUG", "lsst.daf.butler=WARNING", DEBUG, INFO, WARNING, lsstLog_INFO, lsstLog_WARN, DEBUG),
            (".=DEBUG", "DEBUG", DEBUG, DEBUG, DEBUG, lsstLog_DEBUG, lsstLog_DEBUG, DEBUG),
            (".=DEBUG", "conda=DEBUG", DEBUG, INFO, INFO, lsstLog_INFO, lsstLog_INFO, DEBUG),
        )

        self._test_levels(log_levels)

        # Check that the environment variable can set additional roots.
        log_levels = (
            # --log-level / --log-level / expected pyroot, pylsst, pybutler,
            # lsstroot, lsstbutler, lsstx
            (
                "WARNING",
                "lsst.daf.butler=DEBUG",
                WARNING,
                WARNING,
                DEBUG,
                lsstLog_WARN,
                lsstLog_DEBUG,
                WARNING,
            ),
            ("DEBUG", "lsst.daf.butler=TRACE", WARNING, DEBUG, TRACE, lsstLog_DEBUG, lsstLog_DEBUG, DEBUG),
            (".=DEBUG", "lsst.daf.butler=WARNING", DEBUG, INFO, WARNING, lsstLog_INFO, lsstLog_WARN, INFO),
            (".=DEBUG", "DEBUG", DEBUG, DEBUG, DEBUG, lsstLog_DEBUG, lsstLog_DEBUG, DEBUG),
            (".=DEBUG", "conda=DEBUG", DEBUG, INFO, INFO, lsstLog_INFO, lsstLog_INFO, INFO),
        )

        with unittest.mock.patch.dict(os.environ, {"DAF_BUTLER_ROOT_LOGGER": "lsstx"}):
            self._test_levels(log_levels)

    def _test_levels(self, log_levels):
        for level1, level2, x_pyroot, x_pylsst, x_pybutler, x_lsstroot, x_lsstbutler, x_lsstx in log_levels:
            with self.subTest("Test different log levels", level1=level1, level2=level2):
                self.runTest(
                    partial(
                        self.runner.invoke,
                        butlerCli,
                        [
                            "--log-level",
                            level1,
                            "--log-level",
                            level2,
                            "command-log-settings-test",
                            "--expected-pyroot-level",
                            x_pyroot,
                            "--expected-pylsst-level",
                            x_pylsst,
                            "--expected-pybutler-level",
                            x_pybutler,
                            "--expected-lsstroot-level",
                            x_lsstroot,
                            "--expected-lsstbutler-level",
                            x_lsstbutler,
                            "--expected-lsstx-level",
                            x_lsstx,
                        ],
                    )
                )

    def test_helpLogReset(self):
        """Verify that when a command does not execute, like when the help menu
        is printed instead, that CliLog is still reset."""

        self.runTest(partial(self.runner.invoke, butlerCli, ["command-log-settings-test", "--help"]))

    def testLongLog(self):
        """Verify the timestamp is in the log messages when the --long-log
        flag is set."""

        # When longlog=True, loglines start with the log level and a
        # timestamp with the following format:
        # "year-month-day T hour-minute-second.millisecond-zoneoffset"
        # For example: "DEBUG 2020-10-28T10:20:31-07:00 ...""
        # The log level name can change, we verify there is an all
        # caps word there but do not verify the word. We do not verify
        # the rest of the log string, assume that if the timestamp is
        # in the string that the rest of the string will appear as
        # expected.
        timestampRegex = re.compile(
            r".*[A-Z]+ [0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}(.[0-9]{3})"
            "?([-,+][01][0-9]:[034][05]|Z) .*"
        )

        # When longlog=False, log lines start with the module name and
        # log level, for example:
        # lsst.daf.butler.core.config DEBUG: ...
        modulesRegex = re.compile(r".* ([a-z]+\.)+[a-z]+ [A-Z]+: .*")

        with self.runner.isolated_filesystem():
            for longlog in (True, False):
                # The pytest log handler interferes with the log configuration
                # settings set up by initLog -- therefore test by using
                # a subprocess.
                if longlog:
                    args = ("butler", "--log-level", "DEBUG", "--long-log", "create", "here")
                else:
                    args = ("butler", "--log-level", "DEBUG", "create", "here")
                result = subprocess.run(args, capture_output=True)
                # There are cases where the newlines are stripped from the log
                # output (like in Jenkins), since we can't depend on newlines
                # in log output they are removed here from test output.
                output = StringIO((result.stderr.decode().replace("\n", " ")))
                startedWithTimestamp = any([timestampRegex.match(line) for line in output.readlines()])
                output.seek(0)
                startedWithModule = any(modulesRegex.match(line) for line in output.readlines())
                if longlog:
                    self.assertTrue(
                        startedWithTimestamp, msg=f"did not find timestamp in: \n{output.getvalue()}"
                    )
                    self.assertFalse(
                        startedWithModule, msg=f"found lines starting with module in: \n{output.getvalue()}"
                    )
                else:
                    self.assertFalse(startedWithTimestamp, msg=f"found timestamp in: \n{output.getvalue()}")
                    self.assertTrue(
                        startedWithModule,
                        msg=f"did not find lines starting with module in: \n{output.getvalue()}",
                    )

    def testFileLogging(self):
        """Test --log-file option."""
        with self.runner.isolated_filesystem():
            for i, suffix in enumerate([".json", ".log"]):
                # Get a temporary file name and immediately close it
                fd = tempfile.NamedTemporaryFile(suffix=suffix)
                filename = fd.name
                fd.close()

                args = (
                    "--log-level",
                    "DEBUG",
                    "--log-file",
                    filename,
                    "--log-label",
                    "k1=v1,k2=v2",
                    "--log-label",
                    "k3=v3",
                    "create",
                    f"here{i}",
                )

                result = self.runner.invoke(butlerCli, args)
                self.assertEqual(result.exit_code, 0, clickResultMsg(result))

                # Record to test. Test one in the middle that we know is
                # a DEBUG message. The first message might come from
                # python itself since warnings are redirected to log
                # messages.
                num = 4

                if suffix == ".json":
                    records = ButlerLogRecords.from_file(filename)
                    self.assertGreater(len(records), num)
                    self.assertEqual(records[num].levelname, "DEBUG", str(records[num]))
                    self.assertEqual(records[0].MDC, dict(K1="v1", K2="v2", K3="v3"))
                else:
                    with open(filename) as fd:
                        records = fd.readlines()
                    self.assertGreater(len(records), num)
                    self.assertIn("DEBUG", records[num], str(records[num]))
                    self.assertNotIn("{", records[num], str(records[num]))

                self.assertGreater(len(records), 5)

    def testLogTty(self):
        """Verify that log output to terminal can be suppressed."""

        with self.runner.isolated_filesystem():
            for log_tty in (True, False):
                # The pytest log handler interferes with the log configuration
                # settings set up by initLog -- therefore test by using
                # a subprocess.
                if log_tty:
                    args = ("butler", "--log-level", "DEBUG", "--log-tty", "create", "here")
                else:
                    args = ("butler", "--log-level", "DEBUG", "--no-log-tty", "create", "here2")
                result = subprocess.run(args, capture_output=True)

                output = result.stderr.decode()
                if log_tty:
                    self.assertIn("DEBUG", output)
                else:
                    self.assertNotIn("DEBUG", output)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python2.7
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
"""
# Owner: UM System Development Team
# Purpose: Generate wallclock/memory usage plots from rose-stem tests
# Non-standard Dependencies: Matplotlib requires Python 2.7
"""
import glob
import os
import re
import subprocess
import sys
import traceback

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


TRUNK = 'fcm:um.xm_tr@head'

DATADIR = {
    'meto': os.path.join(os.environ['HOME'], 'monitoring'),
    }

PLOTDIR = {
    'meto': os.path.join(os.environ['HOME'], 'public_html', 'monitoring'),
    }

BAD_DATA = -9999
DATA_FILE_EXT = '.dat'
PLOT_FILE_EXT = '.png'

WALLCLOCK_REGEXP = re.compile(r'PE\s*0\s*Elapsed Wallclock Time:\s*(\S+)')
MEMORY_REGEXP = {}


def meto_memory(lines):
    """Met Office-specific code to determine maximum memory from PBS epilogue.
    """
    find_total_mem = re.compile(r'Total Mem\s*(\d+)')
    find_um_atmos_exe = re.compile(r'um-atmos.exe')
    check_for_percentage = re.compile('[0-9]+[%]')
    find_mem_n_units = re.compile(r'(?P<num>[0-9]*\.[0-9]*)(?P<unit>[A-Za-z])')
    memory = BAD_DATA
    for line in lines:
        result = find_total_mem.search(line)
        if result:
            return int(result.group(1))
        if find_um_atmos_exe.match(line):
            split_line = line.split()
            if check_for_percentage.search(split_line[6]):
                mem = find_mem_n_units.search(split_line[5])
                memory = float(mem.group('num'))
                if mem.group('unit') == 'G':
                    memory *= 1000000
                elif mem.group('unit') == 'M':
                    memory *= 1000
            else:
                memory = int(line.split()[6])
    return memory

MEMORY_REGEXP["meto"] = meto_memory


def read_file(fname):
    '''Return contents of a file given filename'''
    with open(fname, 'r') as fhandle:
        lines = fhandle.readlines()
    return lines


def write_file(fname, lines):
    '''Write a file given names and contents'''
    with open(fname, 'w') as fhandle:
        for line in lines:
            fhandle.write("{0}\n".format(line))


def run_command(command):
    '''Run a single command'''

    # Turn into a list
    cmd = command.split()

    # Create the Popen object and connect out and err to pipes
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # Do the communiate and wait to get the results of the command
    stdout, stderr = proc.communicate()
    rcode = proc.wait()

    # Reformat stdout
    stdout = ''.join(stdout)
    stdout = stdout.split("\n")

    stderr = ''.join(stderr)
    stderr = stderr.split("\n")

    if rcode is not 0:
        print "Command {0} failed:".format(command)
        print stdout
        print stderr
        sys.exit(rc)

    return stdout


class FileDoesNotExistError(Exception):
    """Exception class when file does not exist."""
    def __init__(self, filename):
        self.fname = filename

    def __repr__(self):
        return "Expected to find file {0} but it does not exist.".format(
            self.fname)

    __str__ = __repr__


class UnableToParseOutputError(Exception):
    """Exception class when expected contents not found in output."""
    def __init__(self, cmd, descr, stdout):
        self.cmd = cmd
        self.descr = descr
        self.stdout = stdout

    def __repr__(self):
        return "Unable to parse output of {0} for {1}:\n{2}".format(
            self.cmd, self.descr, self.stdout)

    __str__ = __repr__


class JobResources(object):
    """Class representing the resources used by a single job."""

    def __init__(self, name, suitedir, site):
        self.name = name
        self.site = site
        self.datafile = os.path.join(DATADIR[site], name + DATA_FILE_EXT)
        self.plotfile = os.path.join(PLOTDIR[site], name + PLOT_FILE_EXT)
        self.stdout = os.path.join(suitedir, 'log', 'job',
                                   '1', name, 'NN', 'job.out')

    def parse_output(self):
        """Parse the job.log file from a UM atmos task."""
        if not os.path.isfile(self.stdout):
            raise FileDoesNotExistError(self.stdout)
        stdout = read_file(self.stdout)
        self.wallclock = BAD_DATA
        self.memory = BAD_DATA
        for line in stdout:
            result = WALLCLOCK_REGEXP.search(line)
            if result:
                self.wallclock = result.group(1)
        try:
            self.memory = MEMORY_REGEXP[self.site](stdout)
        except Exception as err:
            print "Memory REGEXP function failed"
            print traceback.format_exc()
            print err
            self.memory = BAD_DATA

    def read_old_data(self):
        """Read in old data from disk."""
        self.wallclocks = []
        self.memories = []
        self.trunk_revisions = []

        # If the file doesn't exist, assume the job is new.
        if not os.path.isfile(self.datafile):
            return

        for line in read_file(self.datafile):
            self.add_new_data(*line.split())

    def add_new_data(self, trunk_revision, wallclock, memory):
        """Add current data point to arrays."""
        self.wallclocks.append(wallclock)
        self.memories.append(memory)
        self.trunk_revisions.append(trunk_revision)

    def write_new_data(self):
        """Write new data to disk."""
        newlines = []
        for rev, clock, mem in zip(self.trunk_revisions, self.wallclocks,
                                   self.memories):
            line = "{0} {1} {2}".format(rev, clock, mem)
            newlines.append(line)
        write_file(self.datafile, newlines)

    def generate_plot(self):
        """Given stored data create a png plot."""
        # Convert to numerical data arrays
        arr_trunk_revisions = np.array(self.trunk_revisions, dtype=int)
        arr_wallclocks = np.array(self.wallclocks, dtype=float)
        arr_memories = np.array(self.memories, dtype=float)

        # Mask out bad data
        ma_arr_trunk_revisions = np.array(arr_trunk_revisions)
        ma_arr_wallclocks = np.ma.masked_equal(arr_wallclocks, BAD_DATA)
        ma_arr_memories = np.ma.masked_equal(arr_memories, BAD_DATA)

        # Plot data
        fig = plt.figure()

        # Test to see if the memory array contains valid data - if it
        # doesn't, don't plot it.
        if ma_arr_memories[-1] is not np.ma.masked:
            # Two sub-plots for wallclock and memory

            # Top plot
            ax1 = fig.add_subplot(211)
            fig.suptitle(self.name)
            ax1.plot(ma_arr_trunk_revisions, ma_arr_wallclocks, 'bx')
            ax1.set_xlabel('Trunk revision')
            ax1.set_ylabel('Wallclock')

            # Bottom plot
            ax2 = fig.add_subplot(212)
            ax2.plot(ma_arr_trunk_revisions, ma_arr_memories, 'rx')
            ax2.set_xlabel('Trunk revision')
            ax2.set_ylabel('Max Memory')

        else:
            # Single plot for wallclock only (memory not available)
            fig.suptitle(self.name)
            plt.plot(ma_arr_trunk_revisions, ma_arr_wallclocks, 'bx')
            plt.ylabel('Wallclock')
            plt.xlabel("Trunk revision")

        fig.savefig(self.plotfile)


def find_atmos_tasks(suitedir):
    """Given a suite name, find atmos tasks."""
    logdir = os.path.join(suitedir, 'log', 'job',
                          '1')
    tasks = glob.glob("{0}/atmos*".format(logdir))
    return tasks


def ascertain_site():
    """Find SITE."""
    cmd = 'rose config'
    stdout = run_command(cmd)
    for line in stdout:
        result = re.search(r'SITE=(\w+)', line)
        if result:
            return result.group(1)
    raise UnableToParseOutputError('rose config', 'site', stdout)


def ascertain_trunk_revision(suitedir):
    """Get trunk revision. First try the rose-suite-run.conf and look for
    SOURCE_UM_MIRROR revision, else fallback to 'fcm info'."""

    rose_suite_run_conf = os.path.join(suitedir, 'log', 'rose-suite-run.conf')
    if os.path.isfile(rose_suite_run_conf):
        lines = read_file(rose_suite_run_conf)
        for line in lines:
            result = re.search(r'SOURCE_UM_MIRROR=.*@(\d+)', line)
            if result:
                return result.group(1)

    cmd = 'fcm info {0}'.format(TRUNK)
    stdout = run_command(cmd)
    for line in stdout:
        result = re.search(r'Last Changed Rev:\s*(\S+)', line)
        if result:
            return result.group(1)
    raise UnableToParseOutputError('fcm info', 'most recent trunk revision',
                                   stdout)


def main():
    """Main program."""

    if len(sys.argv) < 2:
        sys.exit("Syntax: monitoring.py </path/to/cylc-run/suitename>")

    suitedir = sys.argv[1]
    site = ascertain_site()
    trunk_revision = ascertain_trunk_revision(suitedir)
    tasks = find_atmos_tasks(suitedir)

    for task in tasks:
        name = os.path.basename(task)
        job = JobResources(name, suitedir, site)
        job.parse_output()
        job.read_old_data()
        job.add_new_data(trunk_revision, job.wallclock, job.memory)
        job.write_new_data()
        job.generate_plot()


if __name__ == '__main__':
    main()

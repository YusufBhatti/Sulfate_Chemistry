#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

from __future__ import print_function
from optparse import OptionParser
from run_rose_stem import (source_profile, read_file)
import glob
import os
import re
import subprocess
import sys

HOME = os.environ['HOME']

OUTPUT_DIR = 'automated_testing/'


def parse_options():
    parser = OptionParser(description="Takes a suite name as an argument")
    return parser.parse_args()


def get_command(fname, suite_name):
    lines = read_file(fname)
    for line in lines:
        if re.search(r'--name={0}\W'.format(suite_name), line):
            command = re.findall(r"'([^']*)'", line)
            return command
    return None


def next_version(fname):
    '''Ascertain if the original suite ran with FCM_VERSION=next.'''
    lines = read_file(fname)
    for line in lines:
        if 'Modifying FCM_VERSION' in line:
            return True
    return False

if __name__ == '__main__':
    # Get the suite name from the command line
    opts, args = parse_options()
    if len(args) < 1:
        sys.exit("Syntax: restart_stem.py <suite name>")
    suite_name = args[0]

    # Get a list of the log files in mtime order
    search_dir = os.path.join(HOME, OUTPUT_DIR)
    if not os.path.isdir(search_dir):
        sys.exit("{0} is not a directory".format(search_dir))
    logfiles = filter(os.path.isfile, glob.glob(search_dir + "*.out"))
    logfiles.sort(key=os.path.getmtime)

    # For each file, see if it contains the command we're looking for
    for most_recent_file in reversed(logfiles):
        command = get_command(most_recent_file, suite_name)
        # If the command is found, rerun the suite
        if command:
            print ("Found {0} in log file {1}".format(suite_name,
                                                      most_recent_file))
            env_profile = source_profile()
            # If FCM_VERSION was next, Rose and Cylc should be too
            if next_version(most_recent_file):
                print ("Setting ROSE/CYLC_VERSION=next")
                env_profile['ROSE_VERSION'] = 'next'
                env_profile['CYLC_VERSION'] = 'next'

            # Replace --new with --restart
            command = ['--restart' if x == '--new' else x for x in command]
            print ("Running: {0}".format(" ".join(command)))
            subprocess.check_call(command, env=env_profile)
            sys.exit(0)
    print ("Failed to find suite {0} in {1}".format(suite_name, search_dir))

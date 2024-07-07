#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
"""
# Owner: UM System Development Team
# Purpose: Generate wallclock/memory tables
"""

import glob
import os
import re
import sys

BAD_DATA = -9999
WALLCLOCK_REGEXP = re.compile(r'PE\s*0\s*Elapsed Wallclock Time:\s*(\S+)')


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
                memory = float(line.split()[6])
    return memory


def read_file(fname):
    '''Return contents of a file given filename'''
    with open(fname, 'r') as fhandle:
        lines = fhandle.readlines()
    return lines


def get_wallclock_and_memory(fname):
    '''Returns and memory usage given a UM pe0 filename.'''
    lines = read_file(fname)

    wallclock = BAD_DATA
    memory = BAD_DATA

    # Get wallclock
    for line in lines:
        result = WALLCLOCK_REGEXP.search(line)
        if result:
            wallclock = float(result.group(1))

    # Get memory
    memory = meto_memory(lines)

    return (wallclock, memory)


def find_output_files(suitedir):
    '''Find all UM and coupled pe0 files in a suite.'''
    joblogdir = os.path.join(suitedir, 'log', 'job', '1')

    um_pe0_files = glob.glob("{0}/atmos*/NN/job.out".format(joblogdir))
    coupled_pe0_files = glob.glob("{0}/coupled*/NN/job.out".format(joblogdir))
    return um_pe0_files + coupled_pe0_files


def get_name_from_pe0(fname):
    '''Return job name from path to pe0 file.'''
    elements = fname.split("/")
    for element in elements:
        if 'atmos' in element:
            return element
        elif 'coupled' in element:
            return element
    return fname


def main():
    if len(sys.argv) < 2:
        sys.exit("Syntax: produce_resources.py <suitedir>")
    suitedir = sys.argv[1]
    pe0_files = find_output_files(suitedir)

    if len(pe0_files) == 0:
        sys.exit("No output files found in suite directory {0}".format(
                 suitedir))

    print "|| {0:>72s} || {1:>20s} || {2:>20s} ||".format("'''Name'''",
                                                          "'''Wallclock'''",
                                                          "'''Max Memory'''")
    for output_file in pe0_files:
        (wallclock, memory) = get_wallclock_and_memory(output_file)
        name = get_name_from_pe0(output_file)
        if wallclock == BAD_DATA:
            wallclock = 'Unavailable'
        else:
            wallclock = '{0:>20.5f}'.format(wallclock)
        if memory == BAD_DATA:
            memory = 'Unavailable'
        else:
            memory = '{0:>20.0f}'.format(memory)

        print "|| {0:>72s} || {1:>20s} || {2:>20s} ||".format(
                name, wallclock, memory)


if __name__ == '__main__':
    main()

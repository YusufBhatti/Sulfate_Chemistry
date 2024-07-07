#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

from __future__ import print_function
from optparse import (OptionParser, BadOptionError, AmbiguousOptionError)
import os
import re
import subprocess
import sys

USER = os.environ['USER']
LOCALDATA = os.environ['SCRATCH']
FCM = "fcm"
ROSE = "rose"


class PassThroughOptionParser(OptionParser, object):
    """
    An unknown option pass-through implementation of OptionParser.

    When unknown arguments are encountered, bundle with largs and try again,
    until rargs is depleted.

    sys.exit(status) will still be called if a known argument is passed
    incorrectly (e.g. missing arguments or bad argument types, etc.)
    """
    def _process_args(self, largs, rargs, values):
        while rargs:
            try:
                super(PassThroughOptionParser, self)._process_args(largs,
                                                                   rargs,
                                                                   values)
            except (BadOptionError, AmbiguousOptionError), e:
                largs.append(e.opt_str)


def parse_options():
    parser = PassThroughOptionParser()
    parser.add_option("-N", "--next", dest="version", action="store_true",
                      help="Use Rose/Cylc/FCM next")
    parser.add_option("-n", "--name", dest="name", action="store",
                      help="Name of suite")
    parser.add_option("-s", "--source", dest="source", action="append",
                      help="Source tree to run")
    parser.add_option("-V", "--variable", dest="variable", action="append",
                      help="variables.rc modifications in key=value pairs")

    return parser.parse_args()


def read_file(fname):
    '''Return contents of a file given filename.'''
    with open(fname, 'r') as fh:
        lines = fh.readlines()
    return lines


def write_file(fname, lines):
    '''Write a file given names and contents. The optional newline argument
    adds a newline at the end of each element of the list.'''
    with open(fname, 'w') as fh:
        for line in lines:
            fh.write(line)


def source_profile():
    '''Source /etc/profile and return the environment dictionary.'''
    env_dict = {}
    source_command = (
        subprocess.Popen(["bash", "-c", ". /etc/profile && env"],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE))
    source_command.wait()
    if source_command.returncode == 0:
        for line in source_command.stdout:
            if re.match(r"\s*\w+=", line):
                key, _, value = line.partition("=")
                env_dict[key] = value.rstrip("\n")
    else:
        print ("Error sourcing /etc/profile")
        print (stderr)
        sys.exit(1)
    return env_dict


def modify_file(fname, variable, new_value):
    '''Modify a file to set a variable to a new value.'''
    lines = read_file(fname)
    newlines = []
    for line in lines:
        result = re.search(r'{0}=\s*(\S+)'.format(variable), line)
        if result:
            print (line)
            line = re.sub(result.group(1), new_value, line)
            print (line)
        newlines.append(line)
    write_file(fname, newlines)


if __name__ == '__main__':
    opts, args = parse_options()

    # Recreate the rose stem arguments
    argstring = " ".join(args)
    # The first source will be the working copy checked out later
    if opts.source is None:
        sys.exit("No --source specified.")
    if len(opts.source) > 1:
        for sourcetree in opts.source[1:]:
            argstring += " --source={0}".format(sourcetree)
    if opts.name is None:
        sys.exit("No --name specified.")
    argstring += " --name={0}".format(opts.name)

    # Create a new working copy, deleting any previous one
    workdir = "{0}/rose_{1}".format(LOCALDATA, opts.name)
    if os.path.isdir(workdir):
        subprocess.check_call(["rm", "-rf", workdir])
    subprocess.check_call(["mkdir", "-p", workdir])

    # Check out the first source tree
    url = opts.source[0]
    subprocess.check_call([FCM, "co", url, workdir])

    # Get /etc/profile
    env_profile = source_profile()

    # Set up changes for Rose/Cylc/FCM next
    if opts.version:
        fname = "{0}/rose-stem/rose-suite.conf".format(workdir)
        print ("Modifying FCM_VERSION in {0}".format(fname))
        modify_file(fname, 'FCM_VERSION', "'next'")
        env_profile['ROSE_VERSION'] = 'next'
        env_profile['CYLC_VERSION'] = 'next'

    # Modify host if applicable
    if opts.variable:
        for hostpair in opts.variable:
            elements = hostpair.split("=")
            print ("Modifying {0} to {1}".format(elements[0], elements[1]))
            fname = "{0}/rose-stem/site/meto/variables.rc".format(workdir)
            modify_file(fname, elements[0], "'{0}'".format(elements[1]))

    arglist = argstring.split()
    stemcmd = [ROSE, "stem", "-v", "--new",
               "--source={0}".format(workdir), "--name={0}".format(opts.name),
               "-C", "{0}/rose-stem".format(workdir)] + arglist
    print (stemcmd)
    subprocess.check_call(stemcmd, env=env_profile)

#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
"""Script to ensure the version of Rose and Cylc is compatible with the
   rose-stem suite..

   Owner: UM System Development Team
"""

import os
import sys


VERSION_CHECKS = [("Rose", "ROSE_VERSION", "2017.02.0"),
                  ("Cylc", "CYLC_VERSION", "7.3.0")]


def compare_versions(current, minimum, name):
    """Compare version against minimum."""
    min_elements = minimum.split(".")
    cur_elements = current.split(".")

    if len(min_elements) != len(cur_elements):
        sys.exit("Incompatible lengths of versions: {0} and {1}".format(
                 min_elements, cur_elements))

    for c, m in zip(cur_elements, min_elements):
        c = int(c)
        m = int(m)
        if c > m:
            break
        elif c < m:
            print "Incompatible {0} version detected (current ".format(name) \
                  + "{0}, minimum {1})".format(current, minimum)
            sys.exit("Incompatible {0} version detected (current ".format(name)
                     + "{0}, minimum {1})".format(current, minimum))
    print "{0} version compatible ({1} >= {2})".format(name, current, minimum)


def main():
    """Main program block."""
    for name, env_var, min_version in VERSION_CHECKS:
        if env_var not in os.environ:
            sys.exit("{0} not in environment".format(env_var))
        compare_versions(os.environ[env_var], min_version, name)


if __name__ == '__main__':
    main()

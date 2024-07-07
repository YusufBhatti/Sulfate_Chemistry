#!/usr/bin/env python
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
#
"""site_validator.py

Loop through a dictionary of sites with known working groups
(WORKING_CONFIGS) to validate each one.
"""

from __future__ import print_function
import os
import sys
import subprocess

WORKING_CONFIGS = {'meto': ('developer',
                            'nightly',
                            'weekly',
                            'check_all_groups'),
                   'nci': ('developer',
                           'nightly',
                           'weekly'),
                   'niwa': ('developer',
                            'nightly',
                            'weekly'),
                   'vm': ('all',
                          'developer'),
                   }


def run_command(command, env=None):
    """
    Given a command as a string, run it and return the exit code, standard
    out and standard error. This does not run a command through the
    shell, so only one command can be run.

    A dictionary of enviornment variables can be passed as env.

    """

    # Turn command into a list
    command_list = command.split()

    # Create the Popen object and connect out and err to pipes
    p = subprocess.Popen(command_list, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, env=env)

    # Do the communiate and wait to get the results of the command
    stdout, stderr = p.communicate()
    rc = p.wait()
    return rc, stdout, stderr


def install_suite(site, group, int_testing=False):
    """
    Install a rose stem suite with options
    -S SITE=<site> and --group=<group>.

    The int_testing logical specifies whether INTEGRATION_TESTING should
    be true or false.

    Returns the integer exit code (rc = 0 is assumed to be a successful
    execution) and a string containing the name of the suite.

    """
    install_cmd = "rose stem --group={g} -S SITE='{s}'".format(g=group,
                                                               s=site)
    install_cmd = install_cmd + " --local-install-only --new"
    suitename = "_".join(("validate", site, group))

    parent_suitename = "ROSE_SUITE_NAME"
    source_env = "SOURCE_UM_BASE"
    # If this script is called from within rose stem then
    # we will have to specify the UM source and a nested directory
    # structure for the output.
    # If the rose stem variables were not present, then we assume
    # that this script is being run standalone from a UM working copy.
    if source_env in os.environ and parent_suitename in os.environ:
        install_cmd = install_cmd + \
            " --source={s}".format(s=os.environ[source_env])
        suitename = os.path.join(os.environ[parent_suitename],
                                 suitename)

    if int_testing:
        install_cmd = install_cmd + ' -S INTEGRATION_TESTING=true'

    install_cmd = install_cmd + " --name={s}".format(s=suitename)

    # copy local envionment variables (including ROSE_VERSION,
    # CYLC_VERSION etc.)
    my_env = os.environ.copy()

    rc, stdout, stderr = run_command(install_cmd, env=my_env)

    if rc:
        print("[FAIL] Problem in suite installation.")
        print("[FAIL] SuiteName=", suitename)
        print("[FAIL] INTEGRATION_TESTING =", int_testing)
        print(stderr)

    return rc, suitename


def run_validation(suitename):
    """
    Run cylc validate on the suite with name suitename. Note that
    cylc is left to find the absolute path to suitename and perform
    the validation on the suite.rc file.

    """
    validate_cmd = "cylc validate --strict " + suitename
    rc, stdout, stderr = run_command(validate_cmd)
    if rc:
        print("[FAIL] Suite not valid")
        print(stderr)
    return rc


def housekeep_suite(suitename):
    """
    Tidy up the temporary suite (suitename) with rose suite-clean.
    Return any exit code.

    """
    cleanup = "rose suite-clean --yes " + suitename
    rc, stdout, stderr = run_command(cleanup)
    if rc:
        print("[FAIL] Error in housekeeping suite", suitename)
        print(stdout)
        print(stderr)
    return rc


def validate_rosestem_suite(site, group):
    """
    For the rose stem group <group> at site <site> the suite is installed,
    validated by cylc and housekept. If any of the three tasks fails, this
    function will return the error code at the earliest opportunity.

    """

    for int_testing in (False, True):
        # Run rose stem suite in install-only mode.
        rc, suitename = install_suite(site, group,
                                      int_testing=int_testing)
        if rc:
            return rc

        # Run cylc validate on the suite.rc file.
        rc = run_validation(suitename)
        if rc:
            return rc

        # Clean up the created suite.
        rc = housekeep_suite(suitename)
        if rc:
            return rc

    return rc

if __name__ == '__main__':
    # loop over sites and run validate_rosestem_suite on each group.
    # This will install a rose stem suite for each group, validate
    # the group and then clean up the suite.
    errorflag = False
    for site in WORKING_CONFIGS:
        for group in WORKING_CONFIGS[site]:
            print('[INFO] validating', site, group)
            rc = validate_rosestem_suite(site, group)
            if rc == 0:
                print('[INFO]', site, group, 'validated')
            else:
                print('[FAIL]', site, group, 'did not validate')
                errorflag = True
    if errorflag:
        sys.exit(1)

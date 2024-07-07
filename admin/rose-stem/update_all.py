#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
'''
Run a Rose upgrade macro on all UM and coupled apps. The path to the working
copy base directory and the version to upgrade to should be provided as
arguments.

Syntax: update_all.py <--path=/path/to/working/copy> [--um=version]
                      [--makeum=version] [--createbc=version] [--nproc=nproc]
                      [--force] [--sort] [--delprofs]

Where "version" is the AFTER_TAG of the relevant upgrade macro. Note that 
at least one of the upgrade macros or the --force argument must be specified.

Owner:  UM System Development Team
'''

from optparse import OptionParser
import glob
import os
import re
import subprocess
import sys

ROSE = "rose"
METADIRS = ["um-atmos", "um-fcm-make", "um-createbc", ]


def parse_options(opt=None):
    '''Parse command line options.'''
    parser = OptionParser(
        description='Either the --force argument or one of the upgrade ' +
                    'macro arguments are required.',
        epilog="Use the force to all apps fix you need")
    parser.add_option(
        "-P", "--path", dest="path", action="store",
        help="Path to working copy (defaults to \".\")", default=".")
    parser.add_option(
        "-u", "--um", dest="um", action="store",
        help="Version of UM atmos to upgrade to")
    parser.add_option(
        "-U", "--makeum", dest="makeum", action="store",
        help="Version of FCM_make UM to upgrade to")
    parser.add_option(
        "-S", "--sort", dest="sort", action="store_true",
        help="Option to hash sort and remove duplicate namelists")
    parser.add_option(
        "-D", "--delprofs", dest="profs", action="store_true",
        help="Option to remove unused profile namelists")
    parser.add_option(
        "-C", "--createbc", dest="createbc", action="store",
        help="Version of CreateBC to upgrade to")
    parser.add_option(
        "-F", "--force", dest="force", action="store_true",
        help="Use the force to ensure apps and metadata consistent they are")
    parser.add_option(
        "-j", "--nproc", dest="n_proc", action="store",
        default="4", help="Number of concurrent processes to use")

    (opts, args) = parser.parse_args()

    if (not opts.force and not opts.um and
            not opts.makeum and not opts.createbc):
        parser.print_help()
        sys.exit(1)

    return opts


def run_command(command):
    '''Given a command as a string, run it and return the exit code, standard
    out and standard error. This does not run a command through the
    shell, so only one command can be run.'''

    # Turn command into a list
    command_list = command.split()

    # Create the Popen object and connect out and err to pipes
    p = subprocess.Popen(command_list, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    # Do the communiate and wait to get the results of the command
    stdout, stderr = p.communicate()
    rc = p.wait()
    return rc, stdout, stderr


def read_file(fname):
    '''Return contents of a file given filename'''
    with open(fname, 'r') as fh:
        lines = fh.readlines()
    return lines


def write_file(fname, lines):
    '''Write a file given names and contents'''
    with open(fname, 'w') as fh:
        for line in lines:
            fh.write(line)


def upgrade_apps(apps_in, cmd, version, n_procs, test_transform=False):
    ''''For each valid app, change into that directory and run the provided
         command.
         
        Returns number of successes and number of failures.'''
    # The command should be split for passing to subprocess
    command = cmd.split()

    # 3 lists - one is a copy of the input apps, the others are to store
    # in-progress and finished Popen objects respectively
    apps = list(apps_in)
    running = []
    stopped = []
    num = 0

    # Parallel processing loop (implementing a simple queue-like system based
    # on the above lists)
    while True:
        if (len(running) < n_procs and len(apps) > 0):
            # If unprocessed apps exist and there aren't enough running already
            # take the next app from the list
            app = apps.pop()
            app_name = os.path.basename(app)

            # Check to see if the app needs processing
            if already_at_this_version(app, version):
                print(" * App {0:s} already at version {1:s}".format
                      (app_name, version))
                continue

            # Start running the command - redirect stdout to devnull
            # explicitly, due to the problems with subprocess.PIPE in case
            # of large output.  Note that this launching statement handles
            # the changing of directory to the given app
            with open(os.devnull, "w") as FNULL:
                proc = subprocess.Popen(command, cwd=app, close_fds=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
            # Add the app to the list
            running.append((app_name, proc))

        # Having launched any new apps move on to check if any existing apps
        # have finished running
        if (len(running) > 0):
            for app_name, proc in list(running):
                # Move finished apps from one list to the other
                if proc.poll() is not None:
                    num += 1
                    print(" * Processed ({0}/{1}): {2}"
                          .format(num, len(apps_in), app_name))
                    stopped.append((app_name, proc))
                    running.remove((app_name, proc))
        else:
            # If no more apps are running we should be done
            break

    # Now gather up the return codes for output and report any errors
    success, err = 0, 0
    for app_name, proc in stopped:
        if test_transform:
            stdout = proc.stdout.read()
            result = re.search(
                r'rose.macros.DefaultTransforms:\s*changes:\s*(\d+)',
                stdout)
            if result:
                num_transforms = int(result.group(1))
                if num_transforms > 0:
                    # Failure because of non-zero transforms
                    print("[FAIL] Job: {0:s}\n".format(app_name))
                    print("Non-zero number of transforms made by fixer:")
                    print(stdout)
                    print("Your upgrade macro is likely to be incomplete.")
                    err = err + 1
                elif proc.returncode == 0:
                    # Success if zero transforms and returned ok
                    success = success + 1
                else:
                    # Failure if zero transforms but non-zero return code
                    err = err + 1
                    print("[FAIL] Job: {0:s}\n".format(app_name))
                    print(proc.stderr.read())
        elif proc.returncode == 0:
            success = success + 1
        else:
            err = err + 1
            print("[FAIL] Job: {0:s}\n".format(app_name))
            print(proc.stderr.read())

    return success, err


def get_meta_version(dname):
    fname = os.path.join(dname, 'rose-app.conf')
    lines = read_file(fname)
    meta_ver = 'HEAD'
    for line in lines:
        result = re.search(r'^\s*meta=.*/(\S+)', line)
        if result:
            meta_ver = result.group(1)
            break
    return meta_ver


def replace_meta_ver(apps, ver):
    for app in apps:
        fname = os.path.join(app, 'rose-app.conf')
        newlines = []
        lines = read_file(fname)
        for line in lines:
            if re.search(r'^\s*meta=', line):
                line = re.sub(r'meta=(.*)/.*', r'meta=\1/%s' % (ver), line)
            newlines.append(line)
        write_file(fname, newlines)


def force_metadata_consistency(apps, metadir, n_procs=1):
    macro_cmd = "{0:s} macro --fix -y -M {1:s}".format(ROSE, metadir)
    meta_ver = get_meta_version(apps[0])
    replace_meta_ver(apps, 'HEAD')
    num, err = upgrade_apps(apps, macro_cmd, None, n_procs,
                            test_transform=True)
    replace_meta_ver(apps, meta_ver)
    return num, err


def metadata_check(opts):
    err = 0
    num = 0
    for subdir in METADIRS:
        metapath = os.path.join(opts.path, 'rose-meta', subdir, 'HEAD')
        print("Running 'rose metadata-check' on {0:s}".format(metapath))
        rc, stdout, stderr = run_command(
            "{0:s} metadata-check -C {1:s}".format(ROSE, metapath))
        if rc:
            print("[FAIL] meta-data for {0:s} is invalid".format(subdir))
            err = err + 1
        else:
            num = num + 1
    return num, err


def already_at_this_version(appdir, version):
    """Return true or false depending on whether the app in appdir is already
    at the version specified."""

    if version is None:
        return False
    appfile = os.path.join(appdir, 'rose-app.conf')
    cmd = "{0:s} config --file={1:s} meta".format(ROSE, appfile)
    rc, stdout, stderr = run_command(cmd)
    meta = stdout.rstrip().split("/")
    if meta[1] == version:
        return True
    else:
        return False


if __name__ == '__main__':

    # Check the argument list
    if len(sys.argv) < 2:
        parse_options(['-h'])
    opts = parse_options()

    # Having checked valid options were supplied, Use the force
    # to ensure app/metadata consistency throughout.
    opts.force = True

    # Convert the supplied path to ensure it is absolute
    opts.path = os.path.realpath(opts.path)

    # Convert the n_proc value to an integer
    if opts.n_proc.isdigit():
        opts.n_proc = int(opts.n_proc)
        print ("Using {0} process{1}..."
               .format(opts.n_proc, {1: ""}.get(opts.n_proc, "es")))
    else:
        raise ValueError("Invalid value for n_proc, must be integer")

    # Check that the rose-stem/app subdirectory exists
    appdir = os.path.join(opts.path, 'rose-stem', 'app')
    if not os.path.isdir(appdir):
        sys.exit("Cannot find rose-stem/app subdirectory in {0:s}".format
                 (opts.path))

    # Save the user's current working directory
    cwd = os.getcwd()

    # Generate the meta-data path
    metadir = os.path.join(opts.path, 'rose-meta')

    # Generate list of apps of various sorts
    um_apps = (glob.glob("{0:s}/um*".format(appdir)) +
               glob.glob("{0:s}/coupled*".format(appdir)) +
               glob.glob("{0:s}/scm*".format(appdir)) +
               glob.glob("{0:s}/recon*".format(appdir)))
    createbc_apps = glob.glob("{0:s}/createbc*".format(appdir))
    make_apps = glob.glob("{0:s}/fcm_make*".format(appdir))
    atmos_make_apps = []
    for app in make_apps:
        if re.search(r'drivers', os.path.basename(app)):
            continue
        if re.search(r'ocean', os.path.basename(app)):
            continue
        elif re.search(r'install', os.path.basename(app)):
            continue
        else:
            atmos_make_apps.append(app)

    # Initialise error counter
    total_errors = 0

    if opts.force:
        num, errs = metadata_check(opts)
        print("Successfully tested {0:d} rose-meta.conf files".format(num))
        total_errors = total_errors + errs

    # Upgrade UM/Coupled apps if requested
    if opts.um:
        print ("Upgrading UM atmos apps in {0:s} to version {1:s}".format
               (appdir, opts.um))

        # Generate the command using the working copy's rose-meta directory
        upgrade_cmd = ("{0:s} app-upgrade -a -y -M {1:s} {2:s}".format
                       (ROSE, metadir, opts.um))

        num, err = upgrade_apps(um_apps, upgrade_cmd, opts.um, opts.n_proc)

        print("Successfully upgraded {0:d} apps".format(num))
        if err:
            print("Upgrade of {0:d} apps failed".format(err))
            total_errors = total_errors + err

    # hash sort the stash requests and profile requests within an app.
    if opts.sort:
        print("Hash sorting the stash requests and profiles in {0:s}".format
              (appdir))

        # Generate the command macro
        sort_cmd = ("{0:s} macro -y -M {1:s} ".format(ROSE, metadir) +
                    "stash_indices.TidyStashTransformPruneDuplicated")

        num, err = upgrade_apps(um_apps, sort_cmd, None, opts.n_proc)

        print("Successfully hash sorted {0:d} apps".format(num))
        if err:
            print("Hash Sorting of {0:d} apps failed".format(err))
            total_errors = total_errors + err

    # remove unused profile requests within an app.
    if opts.profs:
        print("Removing unused profiles in {0:s}".format(appdir))

        # Generate the command macro
        sort_cmd = ("{0:s} macro -y -M {1:s} ".format(ROSE, metadir) +
                    "stash_requests.StashProfilesRemoveUnused")

        num, err = upgrade_apps(um_apps, sort_cmd, None, opts.n_proc)

        print("Successfully removed unused profiles in {0:d} apps".format(
            num))
        if err:
            print("Removing unused profiles in {0:d} apps failed".format(err))
            total_errors = total_errors + err

    if opts.force:
        print("Forcing UM atmos apps to be consistent with metadata in " +
              "{0:s}".format(metadir))
        num, err = force_metadata_consistency(um_apps, metadir, opts.n_proc)

        print("Successfully fixed {0:d} apps".format(num))
        if err:
            print("Fixing of {0:d} apps failed".format(err))
            total_errors = total_errors + err

    # Upgrade FCM make apps for UM if requested
    if opts.makeum:
        print("Upgrading FCM_make UM apps in " +
              "{0:s} to version {1:s}".format(appdir, opts.makeum))

        # Generate the command using the working copy's rose-meta directory
        upgrade_cmd = ("{0:s} app-upgrade -a -y -M {1:s} {2:s}".format
                       (ROSE, metadir, opts.makeum))

        num, err = upgrade_apps(atmos_make_apps, upgrade_cmd, opts.makeum,
                                opts.n_proc)

        print("Successfully upgraded {0:d} apps".format(num))
        if err:
            print("Upgrade of {0:d} apps failed".format(err))
            total_errors = total_errors + err

    if opts.force:
        print("Forcing FCM_make UM apps to be consistent with metadata in " +
              "{0:s}".format(metadir))
        num, err = force_metadata_consistency(atmos_make_apps, metadir,
                                              opts.n_proc)

        print("Successfully fixed {0:d} apps".format(num))
        if err:
            print("Fixing of {0:d} apps failed".format(err))
            total_errors = total_errors + err

    # Upgrade CreateBC apps if requested
    if opts.createbc:
        print("Upgrading Createbc apps in {0:s} to version {1:s}".format
              (appdir, opts.createbc))

        # Generate the command using the working copy's rose-meta directory
        upgrade_cmd = ("{0:s} app-upgrade -a -y -M {1:s} {2:s}".format
                       (ROSE, metadir, opts.createbc))

        num, err = upgrade_apps(createbc_apps, upgrade_cmd, opts.createbc,
                                opts.n_proc)

        print("Successfully upgraded {0:d} apps".format(num))
        if err:
            print("Upgrade of {0:d} apps failed".format(err))
            total_errors = total_errors + err

    if opts.force:
        print("Forcing createbc apps to be consistent with metadata in " +
              "{0:s}".format(metadir))
        num, err = force_metadata_consistency(createbc_apps, metadir,
                                              opts.n_proc)

        print("Successfully fixed {0:d} apps".format(num))
        if err:
            print("Fixing of {0:d} apps failed".format(err))
            total_errors = total_errors + err

    # Return user to original working directory
    os.chdir(cwd)

    # Print final output
    print("Run complete.")
    if total_errors > 0:
        print("Errors were detected. Please fix them in your dev branch and " +
              "create a new test\nbranch.")
    else:
        print("No errors detected. Please commit your changes and run " +
              "rose-stem.")

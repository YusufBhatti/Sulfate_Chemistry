#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
"""
This script will parse the rose-ana-comparisons.db produces by a valid
rose-stem suite, writing and executing a sub-script which performs a
KGO update (creating new KGO folders and installing the required files)

"""
import os
import re
import sys
import sqlite3
import argparse
import subprocess
from collections import defaultdict

# Global criteria for updating - any statuses matched by this list will
# be treated as "failed" for the purposes of updating.  When running a
# "new version" update this list will be updated to include the "OK"
# status, causing all tested files to be updated.
UPDATE_CRITERIA = ["WARN", "FAIL"]

NEW_KGO_DIRNAME = "New/Unknown KGO"


def banner(message):
    "Print a simple banner message"
    return "{0}\n* {1} *\n{0}".format("%"*(len(message)+4), message)


def confirm(message, skip=False):
    "Prompt the user for a response"
    if skip:
        print ("Non-interactive mode, skipping question: ")
        print (message)
        return True
    while True:
        check = raw_input(message)
        if check in ["Y", "y", "yes"]:
            print("")
            return True
        elif check in ["N", "n", "no"]:
            print("")
            return False


def write_update_script(kgo_dirs, new_dirname, script):
    "Write the update script to a file"
    total_filesize = 0

    # Start the script with a header
    script.write("#!/bin/bash\n")
    script.write("set -eux\n")

    for kgo_dir in sorted(kgo_dirs.keys()):
        update_dict = kgo_dirs[kgo_dir]
        # Construct the new KGO directory name
        if kgo_dir == NEW_KGO_DIRNAME:
            new_kgo_dir = None
            script.write("#########################\n"
                         "# NEW/UNKNOWN KGO FILES #\n"
                         "#########################\n")
        else:
            new_kgo_dir = os.path.join(os.path.dirname(kgo_dir), new_dirname)

            # Write a sub-header for this KGO directory to help visually split
            # up the file (for when the user checks it by eye)
            script.write("#"*(len(new_kgo_dir)+4)+"\n")
            script.write("# {0} #\n".format(new_kgo_dir))
            script.write("#"*(len(new_kgo_dir)+4)+"\n")

        # Build up several different text sections in the upcoming loop
        # (so that they can be printed in groups)
        keep_description = ["# Files to be kept from previous KGO: "]
        keep_commands = []
        copy_description = ["# Files to be copied from suite: "]
        copy_commands = []
        mkdir_commands = []
        if new_kgo_dir is not None:
            mkdir_commands.append("mkdir -p {0}".format(new_kgo_dir))

        kgo_files = sorted(update_dict.keys())
        for kgo_file in kgo_files:
            source = update_dict[kgo_file]

            # Ensure any subdirectories are created as needed
            subdir = os.path.dirname(kgo_file)
            if subdir != "":
                if new_kgo_dir is not None:
                    full_subdir = os.path.join(new_kgo_dir, subdir)
                    mkdir_commands.append("mkdir -p {0}".format(full_subdir))
                else:
                    mkdir_com = "mkdir -p {0}".format(subdir)
                    if not mkdir_commands:
                        mkdir_commands.append(mkdir_com)
                    elif mkdir_com != mkdir_commands[-1]:
                        mkdir_commands.append(mkdir_com)

            if source is None:
                # Files being kept should be sym-linked back to the previous
                # revision, to save space
                if "OK" in UPDATE_CRITERIA:
                    # However, if this is a version upgrade any tasks
                    # without a source here must have been retired from
                    # the tests
                    continue
                keep_description.append(
                    "#    * {0}".format(kgo_file))

                # Construct the paths
                old_file = os.path.join(kgo_dir, kgo_file)
                new_file = os.path.join(new_kgo_dir, kgo_file)

                keep_commands.append(
                    "ln -s {0} {1}".format(
                        os.path.relpath(old_file,
                                        os.path.dirname(new_file)),
                        new_file))
            else:
                # Files from the suite should be copied from the suite
                copy_description.append(
                    "#    * {0}".format(kgo_file))
                new_file = os.path.join(new_kgo_dir, kgo_file)
                copy_commands.append(
                    "cp {0} {1}".format(source, new_file))
                # In this case more disk-space is needed, so update
                # the running total for reporting later
                total_filesize += os.path.getsize(source)

        # Can now write this tasks output to the file:
        if len(keep_description) > 1:
            script.write("\n".join(keep_description) + "\n\n")
        if len(copy_description) > 1:
            script.write("\n".join(copy_description) + "\n\n")
        script.write("# Creating new KGO directories:\n")
        script.write("\n".join(mkdir_commands) + "\n\n")
        if len(keep_commands) > 0:
            script.write("# Linking back to unchanged KGO files:\n")
            script.write("\n".join(keep_commands) + "\n\n")
        if len(copy_commands) > 0:
            script.write("# Copying new KGO files from suite:\n")
            script.write("\n".join(copy_commands) + "\n\n")

    return total_filesize


def report_space_required(total_filesize, skip=False):
    "Let the user know how much space they need for the update"
    # Report the amount of disk space required
    units = ["TB", "GB", "MB", "kB", "B"]
    while total_filesize > 1024.0 and len(units) > 1:
        total_filesize = total_filesize / 1024.0
        units.pop()
    unit = units[-1]

    print(banner("Update will require: {0:.2f} {1} of disk space"
                 .format(total_filesize, unit)))
    if not confirm("Please confirm this much space is available (y/n)? ",
                   skip=skip):
        sys.exit("Aborting...")


def add_untested_kgo_files(kgo_dirs):
    "Add any files present in the KGO directory but not tested by the suite"
    for kgo_dir, update_dict in kgo_dirs.iteritems():
        for path, _, filenames in os.walk(kgo_dir):
            for filename in filenames:
                kgo_file = (
                    os.path.relpath(os.path.join(path, filename), kgo_dir))
                if kgo_file not in update_dict:
                    kgo_dirs[kgo_dir][kgo_file] = None
    return kgo_dirs


def group_comparisons_by_dir(comparisons, skip=False):
    """
    Group the key details from the list of comparisons by the KGO
    directory they are referencing.

    """
    kgo_dirs = defaultdict(dict)
    for _, kgo_file, suite_file, status, _ in comparisons:

        # If the kgo_file isn't set, move on (in some cases the test is not
        # comparing two files at all
        if kgo_file is None:
            continue

        # If it looks like this KGO file never existed at all, we can install
        # it directly to the location given in the database
        if not os.path.exists(kgo_file) and status.strip() in UPDATE_CRITERIA:
            if confirm("KGO file {0} doesn't appear to exist, should we "
                       "install it from the suite (y/n)? "
                       .format(kgo_file), skip=skip):
                kgo_dirs[NEW_KGO_DIRNAME][kgo_file] = suite_file
            continue

        # Extract the base KGO directory by moving backwards through the path
        # until a directory that matches the KGO directory naming style
        basedir = os.path.dirname(kgo_file)
        while (basedir != "/" and not
               re.match(r".*(vn\d+\.\d+(_t\d+|))$", basedir)):
            basedir = os.path.dirname(basedir)

        # If the above goes wrong it will eventually hit root; if this happens
        # we cannot continue as something is seriously not right
        if basedir == "/":
            msg = ("Problem locating KGO directory - "
                   "is this actually a KGO file?\n  {0}")
            sys.exit(msg.format(kgo_file))

        # Otherwise add the result to the list - for each entry we store the
        # relative path to the KGO file and the full path to the file in the
        # suite which should be used to update it (if it needs updating)
        if status.strip() in UPDATE_CRITERIA:

            # Here we could check if this update has already been lodged,
            # and if it has perhaps we can do some additional checking
            # (for instance, are the changes in answers the same?)
            relative_kgo_path = os.path.relpath(kgo_file, basedir)
            if (relative_kgo_path in kgo_dirs[basedir] and
                    kgo_dirs[basedir][relative_kgo_path] != suite_file):
                # Or not... it isn't clear what could be checked here
                continue

            kgo_dirs[basedir][relative_kgo_path] = suite_file

    return kgo_dirs


def get_all_kgo_comparisons(conn):
    "Retrieve all comparisons related to KGO tasks"
    res = conn.execute("SELECT comp_task, kgo_file, suite_file, "
                       "status, comparison FROM comparisons")
    return res.fetchall()


def check_for_incomplete_tasks(conn, skip=False):
    "Check the database to see if all tasks are complete"
    res = conn.execute("SELECT task_name FROM tasks WHERE completed==?", (1,))
    failed_tasks = res.fetchall()
    if len(failed_tasks) > 0:
        print("WARNING - Some comparisons are either still running or "
              "failed to complete properly:")
        for task in failed_tasks:
            print("  * "+task[0])
        print("Please investigate this manually - no KGO updates will be "
              "made for any of the tasks listed above")
        if not confirm("Continue anyway (y/n)? ", skip=skip):
            sys.exit()


def connect_to_kgo_database(suite_dir):
    "Make a connection to the comparisons database in a given suite"
    db_filename = os.path.join(suite_dir, "log", "rose-ana-comparisons.db")
    if not os.path.exists(db_filename):
        msg = "ERROR: Suite comparison database not found at {0}"
        sys.exit(msg.format(db_filename))
    return sqlite3.connect(db_filename)


def get_suite_dir(user, suite_name):
    "Returns the path to the given suite directory of a given user"
    suite_dir = "~{0}/cylc-run/{1}".format(user, suite_name)
    expansion = os.path.expanduser(suite_dir)
    if expansion == suite_dir or not os.path.exists(expansion):
        msg = "ERROR: Unable to find suite ({0})"
        sys.exit(msg.format(expansion))
    else:
        return expansion


def prompt_user_options():
    "Interactively ask the user for the required update details"

    cached_options = read_cached_options()

    if cached_options is None:
        # If no previous options were loaded above, prompt the user here
        user = raw_input("Suite username      : ")
        suite_name = raw_input("Suite name          : ")

        print ("\nHow should the new KGO directory be named? \n"
               " (1) Standard KGO directory (vnX.X_tYYYY)\n"
               " (2) Some other custom directory name")
        response = ""
        while response not in ["1", "2"]:
            response = raw_input("Use option          : ")
            if response == "1":
                change_num = raw_input("Version number (X.X): ")
                ticket_num = raw_input("Ticket number       : ")
                new_kgo_dir = "vn{0}_t{1}".format(change_num, ticket_num)
                break
            elif response == "2":
                new_kgo_dir = raw_input("New directory       : ")
                break
        print("")
        # Save these options to the cache file
        write_cached_options(user, suite_name, new_kgo_dir)
    else:
        # Otherwise use the cached options found above
        user, suite_name, new_kgo_dir = cached_options

    # Check that the new directory looks like a sensible name - don't refuse
    # to run if not but at least query the user
    if not re.match(r"(vn\d+\.\d+(_t\d+|))", new_kgo_dir):
        print("WARNING: '{0}' doesn't look like a normal version or "
              "KGO directory name".format(new_kgo_dir))
        if not confirm("Continue anyway (y/n)? "):
            sys.exit("Aborting...")

    # Obtain the suite directory from the username and suite name
    suite_dir = get_suite_dir(user, suite_name)

    return suite_dir, new_kgo_dir


def get_cache_filename():
    "Return the filename of the options cache file"
    # Use a cache file to remember the last set of options used
    # It will be removed only after a successful run, so for any
    # failed run it can be used to re-run with the same settings
    return os.path.expanduser("~/.kgo_update_opts")


def remove_cache_file():
    "Deletes the cache file"
    cache_file = get_cache_filename()
    os.remove(cache_file)


def read_cached_options():
    "Read in previously cached options from the cache file"
    cache_file = get_cache_filename()
    # Read in the options from the cached file if found
    if os.path.exists(cache_file):
        print("Found cached options file: ")
        try:
            with open(cache_file, "r") as cache:
                user = cache.readline().strip()
                suite_name = cache.readline().strip()
                new_kgo_dir = cache.readline().strip()
                print("\nSuite username : {0}".format(user))
                print("Suite name     : {0}".format(suite_name))
                print("New directory  : {0}".format(new_kgo_dir))
        except:
            # If the cache file is broken in some way, delete it and
            # continue
            print("ERROR: failed to read cached options file, removing it")
            os.remove(cache_file)
            return None
        finally:
            if confirm("Re-run with these settings (y/n)? "):
                return user, suite_name, new_kgo_dir
    else:
        # If the cache file never existed return None
        return None


def write_cached_options(user, suite_name, new_kgo_dir):
    "Write a set of cached options to the cache file"
    cache_file = get_cache_filename()
    # Save the options to the cache file
    with open(cache_file, "w") as cache:
        cache.write(user+"\n")
        cache.write(suite_name+"\n")
        cache.write(new_kgo_dir+"\n")


def main():
    "Toplevel script function - process arguments and run update"

    print(banner("Starting KGO Update"))

    # Create a quick version of the regular raw description formatter which
    # adds spaces between the option help text
    class BlankLinesHelpFormatter(argparse.HelpFormatter):
        "Formatter which adds blank lines between options"
        def _split_lines(self, text, width):
            return super(
                BlankLinesHelpFormatter, self)._split_lines(text, width) + ['']

    parser = argparse.ArgumentParser(
        usage="%(prog)s [--new-release]",
        description="""
        KGO Update - A semi-interactive tool for updating UM KGO.

        This script will scan the rose_ana comparisons database to
        prepare and run a shell script to update the KGO files for
        any failed tasks (or all tasks if the --new-release option
        is used)
        """,
        formatter_class=BlankLinesHelpFormatter,
        )

    parser.add_argument('--new-release',
                        help="if set, indicates that the task being "
                        "performed is installing a full set of new KGO "
                        "for a major model release instead of simply "
                        "updating failed KGO tasks.",
                        action="store_true",
                        )

    parser.add_argument('--non-interactive',
                        help="if set, no questions will be asked and "
                        "the script will proceed without any checking "
                        "or warnings.  Use with caution.",
                        action="store_true",
                        )

    parser.add_argument('-S', "--suite-name",
                        help="specifies the (cylc-run dir) name of "
                        "the suite containing the KGO database; for "
                        "use with non-interactive mode.  If not given "
                        "non-interactive mode will try to take "
                        "this from the environment ($CYLC_SUITE_NAME)."
                        )

    parser.add_argument('-U', "--suite-user",
                        help="specifies the (cylc-run dir) username "
                        "where the KGO database suite can be found; for "
                        "use with non-interactive mode.  If not given "
                        "non-interactive mode will try to take "
                        "this from the environment ($USER)."
                        )

    parser.add_argument('-N', "--new-kgo-dir",
                        help="specifies the name of the new KGO "
                        "subdirectory; for use with non-interactive "
                        "mode. If not given the name will be set to "
                        "'non_interactive'."
                        )

    args = parser.parse_args()

    if args.new_release:
        # If releasing a new version, extend the criteria for performing
        # a KGO update to include tasks which completed successfully
        UPDATE_CRITERIA.append("OK")

    # Prompt the user for the key options
    if args.non_interactive:

        new_kgo_dir = args.new_kgo_dir
        if new_kgo_dir is None:
            new_kgo_dir = "non_interactive"

        suite_name = args.suite_name
        if suite_name is None:
            suite_name = os.environ.get("CYLC_SUITE_NAME", None)
            if suite_name is None:
                sys.exit("No suite name given")

        suite_user = args.suite_user
        if suite_user is None:
            suite_user = os.environ.get("USER", None)
            if suite_user is None:
                sys.exit("No suite user given")
        suite_dir = get_suite_dir(suite_user, suite_name)
        confirm_skip = True
        interactive = False

    else:

        if (args.suite_name is not None or
                args.new_kgo_dir is not None or
                args.suite_user is not None):
            sys.exit("Running in interactive mode but was passed an "
                     "argument used by non-interactive mode")

        suite_dir, new_kgo_dir = prompt_user_options()
        confirm_skip = False
        interactive = True

    # Make a connection to the database
    conn = connect_to_kgo_database(suite_dir)

    # See if any tasks have not finished yet
    check_for_incomplete_tasks(conn, skip=confirm_skip)

    # Extract the comparisons
    comparisons = get_all_kgo_comparisons(conn)

    # Group these by KGO directory
    kgo_dirs = group_comparisons_by_dir(comparisons, skip=confirm_skip)

    # Add in any files in the previous KGO untested by the suite
    kgo_dirs = add_untested_kgo_files(kgo_dirs)

    # Create a file to hold the update script
    script_path = os.path.expanduser("~/kgo_update_{0}.sh".format(new_kgo_dir))
    if os.path.exists(script_path):
        print("WARNING: Script file for this update already exists at {0} "
              "and will be overwritten".format(script_path))
        if not confirm("Okay to overwrite script file (y/n)? ",
                       skip=confirm_skip):
            sys.exit("Aborting...")

    # Write the script file
    with open(script_path, "w") as script:
        total_filesize = write_update_script(kgo_dirs, new_kgo_dir, script)

    print("Script file written to: {0}\n".format(script_path))

    # Report on the required space
    report_space_required(total_filesize, skip=confirm_skip)

    # Open the file for viewing in user's editor
    if interactive:
        editor = os.environ.get("EDITOR", None)
        if editor is None:
            editor = raw_input("Please provide a suitable text-editor command "
                               "to inspect the script: \n> ")
        print("Opening {0} with {1}\n"
              "Please carefully check the commands in this file and close the "
              "editor when ready to continue...".format(script_path, editor))
        error_stat = subprocess.Popen([editor, script_path]).wait()
        if error_stat != 0:
            sys.exit("Failed to open script file")

    # Final confirmation before running
    print("\nWARNING: Once launched the script will make permanent changes")
    if not confirm("Are you sure you want to continue (y/n)? ",
                   skip=confirm_skip):
        sys.exit("Aborting...")

    print(banner("Running KGO Update commands from {0}".format(script_path)))

    error_stat = subprocess.Popen(["bash", script_path]).wait()
    if error_stat != 0:
        sys.exit("Script did not run successfully, update failed")

    # Print a final message
    print(banner("KGO update completed successfully"))

    print("\nYou should now update the relevant KGO variables to point to "
          "the new KGO folder ({0}), then reload your suite and re-trigger "
          "the failed KGO tasks to double-check that the new KGO is "
          "valid\n".format(new_kgo_dir))

    # Delete the options cache file on successful completion
    if interactive:
        remove_cache_file()

if __name__ == "__main__":
    main()

#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved. 
# For further details please refer to the file COPYRIGHT.txt 
# which you should have received as part of this distribution. 
# *****************************COPYRIGHT******************************* 
"""
Script to assist in rose-stem version number updating
Owner: UM System Development Team
Syntax: release_new_version.py <previous version> <new version> <next version>

Run within the rose-stem subdirectory of a working copy.
Example: The last stable release was UM 10.0 and UM 10.1 is being prepared. In
this case the command is:
$ release_new_version.py 10.0 10.1 10.2
or equivalently:
$ release_new_version.py vn10.0 vn10.1 vn10.2

This script prepares the UM for a new release. It must be run inside the 
rose-stem directory of a working copy (I recommend using a branch rather than 
the trunk) and takes three mandatory arguments:
 * The version number of the previous release.
 * The version number of the new release which is now being created.
 * The version number of the next release. This is solely used to generate the
   filenames of the upgrade macros between the version being created and the 
   next version.

The previous version specified on the command line must be a previous major
release (e.g. vn9.2). The prebuild path will be updated automatically to cope
with interim releases (those with new prebuilds but no new upgrade macro
file and the same metadata).


Method:

 * Update the version in the STASHmaster_A files.
 * Copies the meta-data for um-atmos, um-fcm-make and um-createbc in HEAD
   into a new version directory
 * Updates the site/*/variables.rc file to:
    * Change the KGO version directories to the new version
    * Updates PREBUILD_*_DIR to point to the new version
 * Updates the rose-suite.conf file to:
    * Updates HOST_SOURCE_CASIM to point to the appropriate new keyword
    * Updates HOST_SOURCE_JULES to point to the appropriate new keyword
    * Updates HOST_SOURCE_SOCRATES to point to the appropriate new keyword
    * Update VN to the new version
    * Updates BASE_UM_REV to point to the appropriate new keyword
    * Updates BASE_CASIM_REV to point to the appropriate new keyword
    * Updates BASE_JULES_REV to point to the appropriate new keyword
    * Updates BASE_SOCRATES_REV to point to the appropriate new keyword
 * Updates the version numbers in the UM source code and scripts in:
    * src/control/top_level/um_version_mod.F90
    * bin/um_script_functions
 * Adds final upgrade macros to bring the upgrade pathway to the new version.
    * um-fcm-make:
      * Changes the location fcm-make looks for the config files unless
        the app looks like it's using the rose-stem $HOST_SOURCE_UM_BASE 
        variable or the VM $CONFIG_ROOT_PATH & $CONFIG_REVISION variables
        used for offline mode in standalone suites.
      * Changes the environment variable socrates_rev to the new version number
        unless the app looks like it's using the rose-stem BASE_SOCRATES_REV
        variable.
      * Changes the environment variable jules_rev to the new version number
        unless the app looks like it's using the rose-stem BASE_JULES_REV
        variable.
      * Changes the environment variable jules_rev to the new version number
        unless the app looks like it's using the rose-stem BASE_CASIM_REV
        variable.
      * Changes the environment variable um_rev to the new version number
        unless the app looks like it's using the rose-stem BASE_UM_REV variable.
    * um-atmos:
      * If the environment variable VN is present in the app, update it to the
        new version number.
 * Create empty versionYY_ZZ.py files for the three different upgrade macro 
   types ready for use at the next release.
 * Adds the required "import versionYY_ZZ" to the versions.py files in each of
   the three upgrade macro directories to use the new empty files.
 * Execute the update_all script to upgrade all apps to the new version.
   * This runs metadata-check on the local metadata.
   * Then it executes the "rose app-upgrade" command using the local metadata.
   * It then runs "rose macro --fix" on each app.

The script is then complete. Fcm-make apps for coupled ocean execs are NOT
modified as these work on an independent release cycle using metadata from
another repository.
"""


import glob
import os
import re
import shutil
import subprocess
import sys

FCM="/opt/ukmo/utils/bin/fcm"
ROSE="/opt/ukmo/utils/bin/rose"
UPDATE_ALL="../admin/rose-stem/update_all.py"

def read_file(fname):
    '''Return contents of a file given filename.'''
    with open(fname, 'r') as fh:
        lines = fh.readlines()
    return lines


def write_file(fname, lines, newline=False):
    '''Write a file given names and contents. The optional newline argument
    adds a newline at the end of each element of the list.'''
    with open(fname, 'w') as fh:
        for line in lines:
            if newline:
                fh.write("%s\n"%(line))
            else:
                fh.write(line)


def run_command(command, shell=False):
    '''Given a command as a string, run it and return the exit code, standard
    out and standard error. The optional shell argument allows a shell to
    be spawned to allow multiple commands to be run.'''

    import subprocess

    if shell:
        # Create the Popen object and connect out and err to pipes using
        # the shell=True option.
        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE, shell=True)
    else:
        # Turn command into a list
        command_list = command.split()

        # Create the Popen object and connect out and err to pipes
        p = subprocess.Popen(command_list, stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
                                       
    # Do the communicate and wait to get the results of the command
    stdout, stderr = p.communicate()
    rc = p.wait()

    # Reformat stdout
    stdout = ''.join(stdout)
    stdout = stdout.split("\n")
    
    return rc, stdout, stderr

def update_stashmaster(fname, version):
    '''Update version number contained in STASHmaster file'''
    lines = read_file(fname)
    newlines = []
    for line in lines:
        line = re.sub(r'UM_VERSION=\S*', r'UM_VERSION=%s'%(version), line)
        newlines.append(line)
    write_file(fname, newlines)

def update_suite_variables(fname, version):
    '''Update variables.rc to version. This changes the KGO versions,
    and modifies the PREBUILD_*_DIR paths.'''
    lines = read_file(fname)
    newlines = []
    for line in lines:
        if re.search(r'_KGO\s*=', line):
            line = re.sub(r"=\s*'.*'", r"='vn%s'"%(version), line)
        if re.search(r'PREBUILD.*DIR', line):
            # meto and other 'cylc-run/vnx.y_prebuilds' installations:
            line = re.sub(r'cylc-run/.*?_', r'cylc-run/vn%s_'%(version), line)
            # niwa and other 'UMDIR/vnx.y/' installations:
            line = re.sub(r"UMDIR\+'/vn.*?/", "UMDIR+'/vn%s/"%(version), line)
        newlines.append(line)
    write_file(fname, newlines)

def update_rose_suite_conf(fname, version):
    '''Update rose-suite.conf to version. This updates the HOST_SOURCE_CASIM,
    HOST_SOURCE_JULES and HOST_SOURCE_SOCRATES lines, their non-host counterparts,
    and VN and BASE_*_REV versions.'''
    lines = read_file(fname)
    newlines = []
    for line in lines:
        if re.search(r'VN=', line):
            line = "VN='%s'"%(version)
        if re.search(r'HOST_SOURCE_CASIM=', line):
            line = re.sub(r"HOST_SOURCE_CASIM=.*", r"HOST_SOURCE_CASIM='fcm:casim.xm_tr@um%s'"%(version), line)
        if re.search(r'^SOURCE_CASIM=', line):
            line = re.sub(r"SOURCE_CASIM=.*", r"SOURCE_CASIM='fcm:casim.xm_tr@um%s'"%(version), line)
        if re.search(r'HOST_SOURCE_JULES=', line):
            line = re.sub(r"HOST_SOURCE_JULES=.*", r"HOST_SOURCE_JULES='fcm:jules.xm_tr@um%s'"%(version), line)
        if re.search(r'^SOURCE_JULES=', line):
            line = re.sub(r"SOURCE_JULES=.*", r"SOURCE_JULES='fcm:jules.xm_tr@um%s'"%(version), line)
        if re.search(r'HOST_SOURCE_SOCRATES=', line):
            line = re.sub(r"HOST_SOURCE_SOCRATES=.*", r"HOST_SOURCE_SOCRATES='fcm:socrates.xm_tr@um%s'"%(version), line)
        if re.search(r'^SOURCE_SOCRATES=', line):
            line = re.sub(r"SOURCE_SOCRATES=.*", r"SOURCE_SOCRATES='fcm:socrates.xm_tr@um%s'"%(version), line)
        if re.search(r'BASE_CASIM_REV=', line):
            line = re.sub(r"BASE_CASIM_REV=.*", r"BASE_CASIM_REV='um%s'"%(version), line)
        if re.search(r'BASE_JULES_REV=', line):
            line = re.sub(r"BASE_JULES_REV=.*", r"BASE_JULES_REV='um%s'"%(version), line)
        if re.search(r'BASE_SOCRATES_REV=', line):
            line = re.sub(r"BASE_SOCRATES_REV=.*", r"BASE_SOCRATES_REV='um%s'"%(version), line)
        if re.search(r'BASE_UM_REV=', line):
            line = re.sub(r"BASE_UM_REV=.*", r"BASE_UM_REV='vn%s'"%(version), line)
        newlines.append(line)
    write_file(fname, newlines)

def copy_metadata(version):
    '''Create a vnX.Y metadata. This copies the HEAD metadata to vnX.Y for
    um-atmos, um-fcm-make and um-createbc.'''    
    run_command("%s copy ../rose-meta/um-atmos/HEAD ../rose-meta/um-atmos/vn%s"%(FCM, version))
    run_command("%s copy ../rose-meta/um-fcm-make/HEAD ../rose-meta/um-fcm-make/vn%s"%(FCM, version))
    run_command("%s copy ../rose-meta/um-createbc/HEAD ../rose-meta/um-createbc/vn%s"%(FCM, version))


def update_um_version_mod(fname, version):
    '''Update version number in src/control/top_level/um_version_mod.F90'''
    lines = read_file(fname)
    newlines = []

    major_ver = int(version)
    minor_ver = int((10*(version)) - int(10*major_ver))
    int_version = major_ver*100 + minor_ver

    for line in lines:
        if re.search(r'CHARACTER.*um_version_char\s*=', line):
            line = (
                re.sub(r"(::\s*\w+\s*=)\s*.*", r"\1 '{0}'".format(version), line))
        elif re.search(r'INTEGER.*um_version_int\s*=', line):
            line = (
                re.sub(r"(::\s*\w+\s*=)\s*.*", r"\1 {0}".format(int_version), line))
        newlines.append(line)
    write_file(fname, newlines)

def append_to_file(fname, contents):
    '''Appends contents to file given filename and contents.'''
    lines = read_file(fname)
    lines = lines + contents
    write_file(fname, lines)


def ascertain_most_recent_after_tag(fname):
    '''Returns the most recent AFTER_TAG for a given file. This works by
    creating a dictionary containing all the AFTER_TAGS, then deleting each
    one which is also used in a BEFORE_TAG. The remaining item in the 
    dictionary is therefore the final AFTER_TAG.'''
    index = {}
    lines = read_file(fname)
    num_tags = 0
    for line in lines:
        result = re.search(r'^\s*AFTER_TAG\s*=\s*"(\S+)"', line)
        if result:
            after = result.group(1)
            index[after] = True
            num_tags += 1

    
    for line in lines:
        result = re.search(r'^\s*BEFORE_TAG\s*=\s*"(\S+)"', line)
        if result:
            before = result.group(1)
            for item in index.keys():
                if re.search(r'%s'%(before), item):
                    if before in index:
                        del index[before]

    if num_tags < 2 and 'tXXXX' in index:
        empty = True
    else:
        empty = False

    if len(index.keys()) > 1:
      for after_tag in index.keys():
          if after_tag.endswith('_tXXXX'):
              del index[after_tag]

    for after_tag in index.keys():
        answer = after_tag
        return answer, empty


def generate_upgrade_macro(fname, before_tag, version, ticket, empty):
    '''Generates an upgrade macro for a given version. This upgrade macro
    updates apps to the next major version.'''
    
    before_tag = re.sub(r'_tXXXX', r'', before_tag)

    # Initialise lists
    additional_content = ['\n']
    additional_content2 = ['']

    # um-fcm-make-specific content
    if re.search(r'um-fcm-make', fname):
        # This is a um-fcm-make upgrade macro; update revisions
        additional_content = [
          '        socrates_rev = self.get_setting_value(config, ["env", "socrates_rev"])\n',
          '        if socrates_rev != "$BASE_SOCRATES_REV":\n',
          '            self.change_setting_value(config, ["env", "socrates_rev"],\n',
          '                                      "um%s", forced=True)\n'%(version),
          '        jules_rev = self.get_setting_value(config, ["env", "jules_rev"])\n',
          '        if jules_rev != "$BASE_JULES_REV":\n',
          '            self.change_setting_value(config, ["env", "jules_rev"],\n',
          '                                      "um%s", forced=True)\n'%(version),
          '        casim_rev = self.get_setting_value(config, ["env", "casim_rev"])\n',
          '        if casim_rev != "$BASE_CASIM_REV":\n',
          '            self.change_setting_value(config, ["env", "casim_rev"],\n',
          '                                      "um%s", forced=True)\n'%(version),
          '        um_rev = self.get_setting_value(config, ["env", "um_rev"])\n',
          '        if um_rev != "$BASE_UM_REV":\n',
          '            self.change_setting_value(config, ["env", "um_rev"],\n',
          '                                      "vn%s", forced=True)\n'%(version)
                             ]

    # Content common to both fcm-make metadata types
    if re.search(r'um-fcm-make', fname):
        additional_content2 = [

"        config_revision = self.get_setting_value(config, ['env', 'config_revision'])\n",
"        config_root = self.get_setting_value(config, ['env', 'config_root_path'])\n",
"        if config_root == '$HOST_SOURCE_UM_BASE':\n",
"            pass\n",
"        elif (config_root == '$CONFIG_ROOT_PATH' and\n",
"              config_revision == '$CONFIG_REVISION'):\n",
"            pass\n",
"        else:\n",
"            config_root = 'fcm:um.xm_tr'\n",
"            config_revision = '@vn%s'\n"%(version),
"            self.add_report('env', 'config_root_path', config_root, \n",
"                info='Upgrading fcm_make config version to trunk@vn%s',\n"%(
                                                              version),
"                is_warning=True)\n",
"        self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)\n",
"        self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)\n",
    "\n",
"        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])\n",
"        if prebuild_path == r'$PREBUILD':\n",
"            pass\n",
"        elif re.search(r'/vn\d+\.\d+_', prebuild_path):\n",
"            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn%s_', prebuild_path)\n"%(version),
"        elif re.search(r'/r\d+_', prebuild_path):\n",
"            prebuild_path = re.sub(r'/r\d+_', r'/vn%s_', prebuild_path)\n"%(version),
"        else:\n",
"            prebuild_path = ''\n",
"        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path, forced=True)\n",
    "\n",
                              ]

    # Content specific to um-atmos
    if re.search(r'um-atmos', fname):
        additional_content = [
          '        if self.get_setting_value(config, ["env", "VN"]):\n',
          '            self.change_setting_value(config, ["env", "VN"],\n',
          '                                      "%s", forced=True)\n'%(version),
          '        stashmaster_path = self.get_setting_value(config, ["env", "STASHMASTER"])\n',
          '        if stashmaster_path:\n',
          '           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",\n', 
          '                                     "vn%s/ctldata",\n'%(version), 
          '                                      stashmaster_path)\n',
          '           self.change_setting_value(config, ["env", "STASHMASTER"],\n', 
          '                                     stashmaster_path, forced=True)\n',
                             ]
    
    # Generate the upgrade macro
    tentimes = int(float(version) * 10)
    contents1 = [
      '\nclass vn%s_t%s(rose.upgrade.MacroUpgrade):\n\n'%(tentimes, ticket),
      '    BEFORE_TAG = "%s"\n'%(before_tag),
      '    AFTER_TAG = "vn%s"\n\n'%(version),
      '    def upgrade(self, config, meta_config=None):\n',
      '        """Upgrade configuration to next major version."""\n',
              ] 
              
    contents2 = [
      '        return config, self.reports\n\n'
                ]


    if empty:
        # There is no macro in the file yet; overwrite it
        contents = [
                    'import rose.upgrade\n',
                    'import re\n', 
                    'import sys\n\n',
                    'class UpgradeError(Exception):\n\n',
                    '    """Exception created when an upgrade fails."""\n\n',
                    '    def __init__(self, msg):\n',
                    '        self.msg = msg\n\n',      
                    '    def __repr__(self):\n',
                    '        sys.tracebacklimit = 0\n',
                    '        return self.msg\n\n',
                    '    __str__ = __repr__\n\n\n',
                   ]
        contents = contents + contents1 + additional_content + \
                   additional_content2 + contents2
        write_file(fname, contents)
    else:
        append_to_file(fname, contents1)
        append_to_file(fname, additional_content)
        append_to_file(fname, additional_content2)
        append_to_file(fname, contents2)
    clean_file(fname)
    return


def clean_file(fname):
    '''Clean the file of any remaining tXXXX macros'''
    lines = read_file(fname)
    newlines = []
    in_macro_to_be_deleted = False
    for line in lines:
        if re.search(r'^\s*class', line):
            in_macro_to_be_deleted = False        
        if re.search(r'tXXXX', line):
            in_macro_to_be_deleted = True
        if not in_macro_to_be_deleted:
            newlines.append(line)
    write_file(fname, newlines)


def update_um_script_functions(fname, version):
    '''Update bin/um_script_functions to contain new version number.'''
    lines = read_file(fname)
    newlines = []
    found_vn = 0
    for line in lines:
        if re.search(r'^\s*local\s*_VN=', line):
            line = re.sub(r"=.*", r"=%s"%(version), line)
        newlines.append(line)
    write_file(fname, newlines)


def update_stash_tool(fname, version):
    '''Update admin/stash to contain new version number.'''
    lines = read_file(fname)
    newlines = []
    for line in lines:
        if re.search(r"^DefaultVersion=", line):
            line = "DefaultVersion=%s\n"%(version)
        newlines.append(line)
    write_file(fname, newlines)

    
def execute_upgrade_macro(basedir, version):
    """Run the upgrade_all script on all apps."""

    if not os.path.isfile(UPDATE_ALL):
        print "Error: Unable to find %s"%(UPDATE_ALL)
        sys.exit(1)
    
    command = "%s --path=%s --um=vn%s"%(UPDATE_ALL, basedir, version) + \
              " --makeum=vn%s"%(version) + \
              " --createbc=vn%s"%(version)
    print "Running %s"%(command)
    
    rc, stdout, stderr = run_command(command)
    if rc:
        print "Error: Problem running command.\n"
    return


def add_import(dname, newverstr):
    """Add a line to the versions.py file to import the next versionYY_ZZ.py 
    file."""
    fname = os.path.join(dname, 'versions.py')
    append_to_file(fname, ["from .version%s import *\n"%(newverstr)])
    


def initialise_new_upgrade_macro_file(versionspy, version, descr):
    '''Create a new upgrade macro file. This includes the UpgradeError
    exception and a dummy macro. This file is added to the versions.py file
    using the add_import method above.'''
    
    header = [
      'import rose.upgrade\n',
      'import re\n', 
      'import sys\n\n',
      'class UpgradeError(Exception):\n\n',
      '    """Exception created when an upgrade fails."""\n\n',
      '    def __init__(self, msg):\n',
      '        self.msg = msg\n\n',      
      '    def __repr__(self):\n',
      '        sys.tracebacklimit = 0\n',
      '        return self.msg\n\n',
      '    __str__ = __repr__\n\n\n',
             ]

    vn = int(float(version) * 10)
    contents = [
      '\nclass vn%s_tXXXX(rose.upgrade.MacroUpgrade):\n'%(vn),
      '\n    """Upgrade macro for ticket #XXXX by <author>."""\n\n'
        '    BEFORE_TAG = "vn%s"\n'%(version),
        '    AFTER_TAG = "vn%s_tXXXX"\n\n'%(version),
        '    def upgrade(self, config, meta_config=None):\n',
        '        """Upgrade a %s app configuration."""\n'%(descr),
        '        # Input your macro commands here\n',
        '        return config, self.reports\n',
        '\n'
               ]

    final_contents = header + contents
    write_file(versionspy, final_contents)
    run_command("fcm add %s"%(versionspy))

def ascertain_rose_suite_conf_settings(fname):
    '''Print the current values of INTEGRATION_TESTING et al to screen.'''
    VARIABLES = [ 'HOUSEKEEPING', 'COMPARE_OUTPUT', 'CENTRAL_INSTALL', 
                  'DESKTOP_MODE', 'INSTALL_META_VERSION', 'INTEGRATION_TESTING',
                  'PREBUILDS', 'SITE', 'VN' ]

    lines = read_file(fname)
    for line in lines:
        elements = line.split('=')
        if len(elements) < 2:
            continue
        for var in VARIABLES:
            if elements[0] == var:
                print "%s: %s"%(elements[0], elements[1].rstrip('\n'))

if __name__ == '__main__':
    
    # Check that $PWD is rose-stem directory
    cwd = os.getcwd()
    if not re.search(r'rose-stem$', cwd):
        sys.exit("Please run in the rose-stem subdirectory")


    # Get version number information from command line
    if len(sys.argv) < 3:
        sys.exit("Syntax: update_version.py <previous version> <this version> <next version>\ne.g. 'update_version.py 9.2 10.0 10.1'")

    for i in range(1, 4):
      if 'vn' in sys.argv[i]:
          sys.argv[i] = re.sub(r'vn', r'', sys.argv[i])

    prev_version = float(sys.argv[1])
    this_version = float(sys.argv[2])
    next_version = float(sys.argv[3])

    verstr = "%s_%s"%(int(prev_version*10), int(this_version*10))
    newverstr = "%s_%s"%(int(this_version*10), int(next_version*10))

    print "Previous version: ",prev_version
    print "New version: ",this_version
    print "Next version: ",next_version

    upgrade_macro_file = "../rose-meta/um-atmos/version%s.py"%(verstr)
    
    if os.path.isfile(upgrade_macro_file):
        print "Verified %s exists"%(upgrade_macro_file)
    else: 
        sys.exit("Unable to confirm correct previous or new version number")

    upgrade_macro_file = "../rose-meta/um-atmos/version%s.py"%(newverstr)

    if os.path.isfile(upgrade_macro_file):
        sys.exit("Unable to confirm new or next version number")
    else:
        print "Verified %s does not already exist"%(upgrade_macro_file)
    
    # Update atmos STASHmaster
    print "Updating STASHmaster"
    update_stashmaster("../rose-meta/um-atmos/HEAD/etc/stash/STASHmaster/STASHmaster_A", this_version)
    
    # Copy meta-data for UM atmos, UM make and UM CreateBC
    print "Copying metadata"
    copy_metadata(this_version)
    
    for fname in glob.glob("site/*/variables.rc"):
        print "Updating suite variables in", fname
        update_suite_variables(fname, this_version)

    # Update KGO version, VN and CASIM, JULES, SOCRATES revision in rose-suite-conf
    print "Updating rose-suite.conf"
    update_rose_suite_conf("rose-suite.conf", this_version)

    # Update um_version_mod.F90
    print "Updating um_version_mod"
    update_um_version_mod("../src/control/top_level/um_version_mod.F90", this_version)

    # Update um_script_functions
    print "Updating um_script_functions"
    update_um_script_functions("../bin/um_script_functions", this_version)

    print "Updating stash tool"
    update_stash_tool("../admin/stash", this_version)

    print "This script will now add upgrade macros to bring apps up to"
    print "version vn%s."%(this_version)
    print "Please press CTRL-c if you do not wish to continue."
    ticket = raw_input("Please enter your ticket number: ")
    

    
    print "Generating upgrade macros"
    atmos_after, empty = ascertain_most_recent_after_tag("../rose-meta/um-atmos/version%s.py"%(verstr))
    generate_upgrade_macro("../rose-meta/um-atmos/version%s.py"%(verstr), atmos_after, this_version, ticket, empty)

    ummake_after, empty = ascertain_most_recent_after_tag("../rose-meta/um-fcm-make/version%s.py"%(verstr))
    generate_upgrade_macro("../rose-meta/um-fcm-make/version%s.py"%(verstr), ummake_after, this_version, ticket, empty)

    createbc_after, empty = ascertain_most_recent_after_tag("../rose-meta/um-createbc/version%s.py"%(verstr))
    generate_upgrade_macro("../rose-meta/um-createbc/version%s.py"%(verstr), createbc_after, this_version, ticket, empty)

    print "Initialising New Upgrade Macros"    
    initialise_new_upgrade_macro_file("../rose-meta/um-atmos/version%s.py"%(newverstr), this_version, "UM runtime")
    initialise_new_upgrade_macro_file("../rose-meta/um-fcm-make/version%s.py"%(newverstr), this_version, "UM fcm-make")
    initialise_new_upgrade_macro_file("../rose-meta/um-createbc/version%s.py"%(newverstr), this_version, "CreateBC")

    print "Adding import"
    add_import("../rose-meta/um-atmos", newverstr)
    add_import("../rose-meta/um-fcm-make", newverstr)
    add_import("../rose-meta/um-createbc", newverstr)
    
    something = raw_input("Press CTRL-c to abort, or enter to run update_all.py: ")
    print "Upgrading apps"
    basedir = os.path.abspath(os.path.join(os.getcwd(), '..'))
    execute_upgrade_macro(basedir, this_version)
    print "Complete."

    # Print rose-suite.conf variables for user to manually confirm
    print "\nPlease verify the following settings from rose-suite.conf:"
    ascertain_rose_suite_conf_settings('rose-suite.conf')
    print "Please also check the minimum Rose and Cylc versions in rose-stem/bin/rose_cylc_version_test.py"

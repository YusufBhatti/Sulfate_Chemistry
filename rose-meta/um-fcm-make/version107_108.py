import rose.upgrade
import re
import sys
import os
import glob
import fileinput

class UpgradeError(Exception):

    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__



class vn107_t2524(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2524 by Richard Hill."""

    BEFORE_TAG = "vn10.7"
    AFTER_TAG = "vn10.7_t2524"

    def upgrade(self, config, meta_config=None):
        """Check for models using the, now unsupported, OASIS3 coupler"""

	# Check if OASIS3 coupler is active. If it is then raise an exception.
	# This is chosen in preference to converting to oasis3-mct
	# since that may lead to the job being run with invalid namcouple
	# control files and resource requests.
	# Setting to none would similarly eventually lead to failure
	# some way down the line, so we nip things in the bud at the
	# earliest possible stage.
	coupler = self.get_setting_value(config, ["env", "COUPLER"])
	if coupler == "oasis3":
            raise UpgradeError('This suite uses the OASIS3 coupler which is no longer supported.\n' +
                               'Please edit your suite to select the OASIS3-MCT coupler.\n' +
                               'Please ensure you employ an appropriate namcouple control file.')

        # Make a list of keys to be checked for the presence of OASIS3.
        keys = ['keys_atmos_app', 'keys_atmos', 'keys_atmos_extra']

        # Set up the strings we're interested in.
        OASIS3_flag = 'OASIS3=oasis3'
	MCT_flag = 'MCT=mct'

        for key in keys:
            #  Get the key values:
            keylist = self.get_setting_value(config, ["env", key])

            if keylist is not None:
                # If both OASIS3 and MCT are present then 
                # this is a legitimate OASIS3-MCT model
                # and we just need to remove the defunct 
                # OASIS3 setting. If MCT is not present
                # we need to raise an exception.
                if (re.search(OASIS3_flag,keylist)):
                    if (re.search(MCT_flag,keylist)):
                        # Remove the, now redundant, OASIS3 key, leaving just 
                        # the MCT setting. Splitting and reforming the list
                        # should weed out any white space.          
			keylist = keylist.replace(OASIS3_flag,'') 
                        bits = keylist.split()
                        space = " "
                        keylist = space.join(bits)
                        self.change_setting_value(config, ['env', key], keylist)
                    else:
                        # This is an OASIS3 based model which we can't
                        # upgrade automatically
                        raise UpgradeError('The OASIS3 coupler is no longer supported.\n' +
		                           'Your CPP key settings indicate that you have OASIS3 activated.\n' +
                                           'Please edit your suite to select the OASIS3-MCT coupler.\n' +
                                           'Please ensure you employ an appropriate namcouple control file.')
	
        return config, self.reports


class vn107_t1681(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1681 by Steve Wardle"""

    BEFORE_TAG = "vn10.7_t2524"
    AFTER_TAG = "vn10.7_t1681"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        self.add_setting(config, ["env", "compile_sstpert_lib"],
                         "preprocess-sstpert_lib build-sstpert_lib")
        return config, self.reports


class vn107_t773(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #773 by Jonathan Wilkinson."""

    BEFORE_TAG = "vn10.7_t1681"
    AFTER_TAG = "vn10.7_t773"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app for scm and atmos configurations."""

        config_type = self.get_setting_value(config, [ 'env', 'config_type'] )

        # Test if this is a rose stem app or not. Rose-stem apps will have 
        # jules_sources = $HOST_SOURCE_JULES
        jules_sources = self.get_setting_value(config, [ 'env', 'jules_sources'])

        if jules_sources == "$HOST_SOURCE_JULES":
            # This is a rose-stem app. Add $BASE_CASIM_REV and $HOST_SOURCE_CASIM,
            # which are specific to rose-stem.
            self.add_setting(config, ['env', 'casim_rev'], '$BASE_CASIM_REV')
            self.add_setting(config, ['env', 'casim_sources'] ,'$HOST_SOURCE_CASIM')

        else:
            # This is not a rose-stem app. Set the version to the um10.7 keyword
            # and leave casim_sources empty on upgrade, allowing the user to
            # select their own casim_sources from branches etc
            self.add_setting(config, ["env", "casim_rev"], "um10.7")
            self.add_setting(config, ["env", "casim_sources"], "")


        if config_type == 'atmos' or config_type == 'scm':

            # Next find and update all ./file/*.cfg files within the fcm_make app directory
            cfg_files = glob.glob('./file/*.cfg')

            for fnm in cfg_files:
                # Only add casim_sources to file if it isn't already there.
                add_casim_sources = True
                for line in fileinput.input([fnm], inplace=True):
                    if '$casim_sources' in line:
                        add_casim_sources = False
                    sys.stdout.write(line)

                if add_casim_sources:
                    for line in fileinput.input([fnm], inplace=True):
                        if '$um_sources' in line:
                            line = line.rstrip() + '\n' + 'extract.location{diff}[casim] = $casim_sources' + '\n'
                        sys.stdout.write(line)

                    # Write out a message to report changes in this file
                    print '''\
[U] File upgrade %s-%s:
    file = %s
         Added text 'extract.location{diff}[casim] = $casim_sources'\
''' % (self.BEFORE_TAG, self.AFTER_TAG, fnm)

        else:
            # Not an atmos or scm app:
            print "/file/*.cfg unchanged as not a scm or atmos fcm-make configuration"

        return config, self.reports


class vn107_t3093(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3093 by Roddy Sharp & Steve Wardle"""

    BEFORE_TAG = "vn10.7_t773"
    AFTER_TAG = "vn10.7_t3093"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        self.remove_setting(config, ["env", "compile_packing_lib"])
        return config, self.reports

class vn108_t3162(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.7_t3093"
    AFTER_TAG = "vn10.8"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        socrates_rev = self.get_setting_value(config, ["env", "socrates_rev"])
        if socrates_rev != "$BASE_SOCRATES_REV":
            self.change_setting_value(config, ["env", "socrates_rev"],
                                      "um10.8", forced=True)
        jules_rev = self.get_setting_value(config, ["env", "jules_rev"])
        if jules_rev != "$BASE_JULES_REV":
            self.change_setting_value(config, ["env", "jules_rev"],
                                      "um10.8", forced=True)
        casim_rev = self.get_setting_value(config, ["env", "casim_rev"])
        if casim_rev != "$BASE_CASIM_REV":
            self.change_setting_value(config, ["env", "casim_rev"],
                                      "um10.8", forced=True)
        um_rev = self.get_setting_value(config, ["env", "um_rev"])
        if um_rev != "$BASE_UM_REV":
            self.change_setting_value(config, ["env", "um_rev"],
                                      "vn10.8", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$HOST_SOURCE_UM_BASE':
            pass
        elif (config_root == '$CONFIG_ROOT_PATH' and
              config_revision == '$CONFIG_REVISION'):
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.8'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.8',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn10.8_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn10.8_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path, forced=True)

        return config, self.reports


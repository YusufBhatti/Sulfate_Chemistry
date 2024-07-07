import rose.upgrade
import re
import sys
import copy

class UpgradeError(Exception):

    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__



class vn106_t2572(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.6"
    AFTER_TAG = "vn10.6.1"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        socrates_rev = self.get_setting_value(config, ["env", "socrates_rev"])
        if socrates_rev != "$BASE_SOCRATES_REV":
            self.change_setting_value(config, ["env", "socrates_rev"],
                                      "um10.6.1", forced=True)
        jules_rev = self.get_setting_value(config, ["env", "jules_rev"])
        if jules_rev != "$BASE_JULES_REV":
            self.change_setting_value(config, ["env", "jules_rev"],
                                      "um10.6.1", forced=True)
        um_rev = self.get_setting_value(config, ["env", "um_rev"])
        if um_rev != "$BASE_UM_REV":
            self.change_setting_value(config, ["env", "um_rev"],
                                      "vn10.6.1", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$HOST_SOURCE_UM_BASE':
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.6.1'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.6.1',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)

        return config, self.reports


class vn106_t2607(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2607 by Sam Cusworth."""

    BEFORE_TAG = "vn10.6.1"
    AFTER_TAG = "vn10.6_t2607"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        # Input your macro commands here
        self.remove_setting(config, ["env", "compile_fieldop"])
        return config, self.reports


class vn106_t2480(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2480 by Paul Cresswell."""

    BEFORE_TAG = "vn10.6_t2607"
    AFTER_TAG = "vn10.6_t2480"

    def upgrade(self, config, meta_config=None):
        """Warn anyone using a monsoon config that their suite needs updating.
"""
        platform = self.get_setting_value(config,
                                          ['env', 'platform_config_dir'])
        if 'monsoon-xc40-' in platform:
            warn_msg = """
  !!!                              WARNING!
  !!! From UM 10.7 onwards the central fcm-make files for the MONSooN
  !!! are configured for remote extract, with no mirror step.
  !!!
  !!! This means fcm_make2 tasks are no longer required for UM compilations.
  !!!
  !!! Please update your suite accordingly. See the UM release notes for
  !!! further details.
"""
            self.add_report(info=warn_msg, is_warning=True)

        return config, self.reports

class vn106_t2643(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2643 by Paul Cresswell."""

    BEFORE_TAG = "vn10.6_t2480"
    AFTER_TAG = "vn10.6_t2643"

    def examine_keys(self, keylist, expected):
       """Compare the supplied keys with those expected for a given exec"""
       # Get portio version (common to all execs)
       portio = None
       if re.search('C95_2A', keylist):
           portio = '2A'
       elif re.search('C95_2B', keylist):
           portio = '2B'

       # Timer is also common enough to be worth doing once here:
       timer = None
       if re.search('C97_1A', keylist):
           timer = '1A'
       elif re.search('C97_3A', keylist):
           timer = '3A'
       elif re.search('C97_4A', keylist):
           timer = '4A'

       # Check all the expected keys are present:
       complete = True
       for key in expected:
           if not re.search(key, keylist):
               complete = False

       # Remove all the expected keys and see if anything is left over:
       keys = copy.deepcopy(keylist)
       for key in expected:
           keys = re.sub(r'(?i)({0}\w*?)=\1'.format(key), '', keys)
       # Turn what remains into a list:
       extra_keys = keys.split()

       return complete, portio, timer, extra_keys

    def remove_valid_extras(self, keys):
        """Remove any extra keys that aren't mandatory but may be expected"""
        expected_extras = ['OASIS3', 'MCT', 'C_DP_HLM', 'DRHOOK', 'RECON']
        for key in expected_extras:
            if '{0}={1}'.format(key,key.lower()) in keys:
                keys.remove('{0}={1}'.format(key,key.lower()))
        return keys

    def remove_keylist(self, config, extra, name):
        """Remove the keylist and add back any extras separately"""
        self.remove_setting(config, ["env", "keys_{0}_app".format(name)])
        if len(extra) > 0:
            # Remove extras that might be there for a valid reason:
            extra = self.remove_valid_extras(extra)
            if len(extra) > 0:   # May have shrunk!
                self.add_setting(config, ["env", "keys_{0}_extra".format(name)],
                                 ' '.join(extra))
        return config

    def upgrade(self, config, meta_config=None):
        """Replace complete CPP key overrides with additions only."""

        config_type = self.get_setting_value(config, ["env", "config_type"])

        warn_cpp = False
        if config_type == 'atmos':
            atmos_keys = ['C84', 'C95', 'C96', 'C97', 'UM_JULES']
            recon_keys = ['RECON', 'C95', 'UM_JULES']

            keylist = self.get_setting_value(config, ["env", "keys_atmos_app"])
            if keylist is not None:
                warn_cpp = True
                complete, portio_atmos, timer, extra = self.examine_keys(keylist, atmos_keys)
                # Other checks not valid for all execs:
                if (re.search('OASIS3=oasis3',keylist)
                and re.search('MCT=mct',keylist)):
                    coupler = 'oasis3_mct'
                elif re.search('OASIS3=oasis3', keylist):
                    coupler = 'oasis3'
                else:
                    coupler = 'none'

                if not re.search('C_DP_HLM=c_dp_hlm', keylist):
                    solver = 'single'
                else:
                    solver = 'double'

                self.change_setting_value(config, ["env", "timer_version"],
                                          timer, forced=True)
                self.change_setting_value(config, ["env", "COUPLER"],
                                          coupler, forced=True)
                self.change_setting_value(config, ["env", "solver_precision"],
                                          solver, forced=True)

                # Does the keylist contain everything necessary?
                if complete:
                    config = self.remove_keylist(config, extra, 'atmos')

            else:
                portio_atmos = None

            # --- RECON ---
            keylist = self.get_setting_value(config, ["env", "keys_recon_app"])
            if keylist is not None:
                warn_cpp = True
                complete, portio_recon, _, extra = self.examine_keys(keylist, recon_keys)
                if re.search('RECON_SERIAL=recon_serial', keylist):
                    recon_mpi = 'serial'
                else:
                    recon_mpi = 'parallel'

                self.change_setting_value(config, ["env", "recon_mpi"],
                                          recon_mpi, forced=True)

                if complete:
                    config = self.remove_keylist(config, extra, 'recon')

            else:
                portio_recon = None

            # atmos portio settings take precedence, unless recon-only build:
            compile_atmos = self.get_setting_value(config, 
                                                   ["env", "compile_atmos"])
            compile_recon = self.get_setting_value(config, 
                                                   ["env", "compile_recon"])
            if (compile_atmos == '' and
                 compile_recon == 'preprocess-recon build-recon'):
                if portio_recon is not None:
                    self.change_setting_value(config, ["env", "portio_version"],
                                              portio_recon, forced=True)
            else:
                if portio_atmos is not None:
                    self.change_setting_value(config, ["env", "portio_version"],
                                              portio_atmos, forced=True)

        elif config_type == 'scm':
            scm_keys = ['SCMA', 'C84', 'C95', 'C96', 'C97', 'UM_JULES']
            keylist = self.get_setting_value(config, ["env", "keys_scm_app"])
            if keylist is not None:
                warn_cpp = True
                complete, portio, timer, extra = self.examine_keys(keylist, scm_keys)
                self.change_setting_value(config, ["env", "portio_version"],
                                          portio, forced=True)
                self.change_setting_value(config, ["env", "timer_version"],
                                          timer, forced=True)
                if complete:
                    config = self.remove_keylist(config, extra, 'scm')

        elif config_type == 'createbc':
            createbc_keys = ['C95', 'C97', 'UTILIO', 'UTILS']
            keylist = self.get_setting_value(config, ["env", "keys_createbc_app"])
            if keylist is not None:
                warn_cpp = True
                complete, portio, timer, extra = self.examine_keys(keylist, createbc_keys)
                self.change_setting_value(config, ["env", "portio_version"],
                                          portio, forced=True)
                self.change_setting_value(config, ["env", "timer_version"],
                                          timer, forced=True)
                if complete:
                    config = self.remove_keylist(config, extra, 'createbc')

        elif config_type == 'utils-mpp':
            crmstyle_keys = ['C95', 'C97', 'UTILIO']
            keylist = self.get_setting_value(config, ["env", "keys_crmstyle_coarse_grid_app"])
            if keylist is not None:
                warn_cpp = True
                complete, portio, timer, extra = self.examine_keys(keylist, crmstyle_keys)
                self.change_setting_value(config, ["env", "portio_version"],
                                          portio, forced=True)
                self.change_setting_value(config, ["env", "timer_version"],
                                          timer, forced=True)
                if complete:
                    config = self.remove_keylist(config, extra,
                                                 'crmstyle_coarse_grid')

        elif config_type == 'utils-serial':
            keydict = {'cumf' : ['C95', 'UTILIO', 'CONVIEEE', 'CUMF', 'FLDC'],
                       'pumf' : ['C95', 'UTILIO', 'PUMF'],
                       'fieldcalc' : ['C95', 'UTILIO', 'FLDCALC'],
                       'vomext' : ['C95', 'UTILIO', 'VOMEXT', 'UTILS'],
                       'fieldcos' : ['C95', 'FLDC', 'UTILIO', 'UTILS'],
                       'convieee' : ['C95', 'UTILIO', 'CONVIEEE', 'FLDC', 'IEEE'],
                       'pptoanc' : ['C95', 'UTILIO', 'UTILS', 'PPTOANC'],
                       'fieldop' : ['C95', 'UTILIO', 'UTILS', 'FLDOP'],
                       'fieldmod' : ['C95', 'UTILIO', 'UTILS', 'FLDMOD'],
                       'convpp' : ['C95', 'UTILIO', 'CONVPP', 'CONVPP64']
                      }

            for util, util_keys in keydict.items():
                keylist = self.get_setting_value(config,
                          ["env", "keys_{0}_app".format(util)])
                if keylist is not None:
                    warn_cpp = True
                    complete, portio, _, extra = self.examine_keys(keylist,
                                                 util_keys)
                    self.change_setting_value(config, ["env", "portio_version"],
                                               portio, forced=True)
                    if complete:
                        config = self.remove_keylist(config, extra, util)

        elif config_type == 'utils-serial-32B':
            keydict = {'cumf' : ['C95', 'UTILIO', 'CONVIEEE', 'CUMF', 'FLDC'],
                       'pumf' : ['C95', 'UTILIO', 'PUMF'],
                       'convpp' : ['C95', 'UTILIO', 'CONVPP', 'FLDC', 'UTILS']
                      }

            for util, util_keys in keydict.items():
                keylist = self.get_setting_value(config,
                          ["env", "keys_{0}_app".format(util)])
                if keylist is not None:
                    warn_cpp = True
                    complete, portio, _, extra = self.examine_keys(keylist,
                                                 util_keys)
                    self.change_setting_value(config, ["env", "portio_version"],
                                               portio, forced=True)
                    if complete:
                        config = self.remove_keylist(config, extra, util)

        if warn_cpp:
            warn_msg = """
  !!!                              WARNING!
  !!! A "keys_*_app" variable was found. Your CPP keys for preprocessing
  !!! may have been updated. Please check the preprocessing settings of your
  !!! fcm_make app are still set as intended.
"""
            self.add_report(info=warn_msg, is_warning=True)

        return config, self.reports


class vn107_t2744(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.6_t2643"
    AFTER_TAG = "vn10.7"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        socrates_rev = self.get_setting_value(config, ["env", "socrates_rev"])
        if socrates_rev != "$BASE_SOCRATES_REV":
            self.change_setting_value(config, ["env", "socrates_rev"],
                                      "um10.7", forced=True)
        jules_rev = self.get_setting_value(config, ["env", "jules_rev"])
        if jules_rev != "$BASE_JULES_REV":
            self.change_setting_value(config, ["env", "jules_rev"],
                                      "um10.7", forced=True)
        um_rev = self.get_setting_value(config, ["env", "um_rev"])
        if um_rev != "$BASE_UM_REV":
            self.change_setting_value(config, ["env", "um_rev"],
                                      "vn10.7", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$HOST_SOURCE_UM_BASE':
            pass
        elif (config_root == '$CONFIG_ROOT_PATH' and
              config_revision == '$CONFIG_REVISION'):
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.7'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.7',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn10.7_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn10.7_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path, forced=True)

        return config, self.reports


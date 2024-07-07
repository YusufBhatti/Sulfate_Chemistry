import rose.upgrade
import re
import fileinput
import sys

class UpgradeError(Exception):

      """Exception created when an upgrade fails."""
      
      def __init__(self, msg):
          self.msg = msg
      
      def __repr__(self):
          sys.tracebacklimit = 0
          return self.msg
          
      __str__ = __repr__


class vn90_t6074(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6074 by Steve Wardle."""

    BEFORE_TAG = "vn9.0"
    AFTER_TAG = "vn9.0_t6074"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""

        c98_keys = ["C98_1A=c98_1a","C98_0A=c98_0a"]
        
        env_vars = ["keys_atmos_app","keys_recon_app","keys_scm_app"]

        for env_var in env_vars:
            app_keys = self.get_setting_value(config,["env",env_var])
            if app_keys:
                app_keys = app_keys.split()
                for c98_key in c98_keys:
                    if c98_key in app_keys:
                        app_keys.remove(c98_key)
                self.change_setting_value(config, ["env",env_var],
                                          " ".join(app_keys))

        return config, self.reports


class vn90_t6195(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6195 by Paul Cresswell.
       Break up the config path into a number of separate options."""

    BEFORE_TAG = "vn9.0_t6074"
    AFTER_TAG = "vn9.0_t6195"

    def upgrade(self, config, meta_config=None):
        # Get the current config value:
        inc_config = self.get_setting_value(config, ["env", "include_config"])

        # Break it up into a number of separate options:
        searchObj = re.search(r'^(.*)/fcm-make/(.*)/um-(.*)-(.*)\.cfg(.*)$',
                              inc_config)

        root_path = searchObj.group(1)
        platform_config = searchObj.group(2)
        config_type = searchObj.group(3)
        optimisation = searchObj.group(4)
        config_rev = searchObj.group(5)

        self.add_setting(config, ["env", "config_root_path"], root_path)
        self.add_setting(config, ["env", "platform_config_dir"], 
                         platform_config)
        self.add_setting(config, ["env", "config_type"], config_type)
        self.add_setting(config, ["env", "optimisation_level"], optimisation)
        self.add_setting(config, ["env", "config_revision"], config_rev)

        # Add to the config file
        filename = './file/fcm-make.cfg'
        new_config = "include = %s/fcm-make/%s/um-%s-%s.cfg%s\n" % (
                     '$config_root_path', '$platform_config_dir',
                     '$config_type', '$optimisation_level', '$config_revision')
        for line in fileinput.input([filename], inplace=True):
            if re.match(r'\s*include\s*=', line):
                line = new_config
            sys.stdout.write(line)

        # Remove the old variable
        self.remove_setting(config, ["env", "include_config"])

        return config, self.reports


class vn90_t6274(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6274 by Paul Cresswell.
       Make both opt. level and OpenMP easily selectable from the GUI."""

    BEFORE_TAG = "vn9.0_t6195"
    AFTER_TAG = "vn9.0_t6274"

    def upgrade(self, config, meta_config=None):
        # Get the current config settings:
        config_type = self.get_setting_value(config, ["env", "config_type"])
        platform    = self.get_setting_value(config, 
                                             ["env", "platform_config_dir"])
        fcflags_omp = self.get_setting_value(config, ["env", "fcflags_omp"])
        ldflags_omp = self.get_setting_value(config, ["env", "ldflags_omp"])

        # Assume OpenMP on to begin
        openmp = 'openmp_on'

        # Work through config types to best-guess OpenMP settings:
        if re.match(r'utils-serial', config_type):  # serial and serial-32B
            openmp = 'openmp_off'
            if config_type == 'utils-serial' and platform == 'ncas-xc30-cce':
                openmp = 'openmp_on'
            if platform == 'meto-pwr7-xlf':
                openmp = 'openmp_on'

        # The SCM will not compile with OpenMP yet.
        if config_type == 'scm':
            openmp = 'openmp_off'
 
        # Now check to see if the user has overridden anything.
        warn_user = False
        blank_fcflags = False
        blank_ldflags = False
        if fcflags_omp is not None:
            warn_user = True
            # Use this flag to decide if the user is overriding OpenMP settings.
            if re.search(r'\S', fcflags_omp):
                openmp = 'openmp_on'
            else:
                openmp = 'openmp_off'
                blank_fcflags = True
 
        if ldflags_omp is not None:
            warn_user = True
            if not re.search(r'\S', ldflags_omp):
                blank_ldflags = True

        # If both flags are blank (and so openmp=off) assume we can safely
        # delete both variables. Otherwise leave them as ignored so they're
        # still available to the user if the upgrade settings are wrong.
        if blank_fcflags and blank_ldflags:
            self.remove_setting(config, ["env", "fcflags_omp"])
            self.remove_setting(config, ["env", "ldflags_omp"])
        else:
            self.ignore_setting(config, ["env", "fcflags_omp"])
            self.ignore_setting(config, ["env", "ldflags_omp"])
 
        if warn_user:
            warn_msg = """
     !!!!! OpenMP compilation and/or link settings have been disabled. !!!!!
     !!!!! Please check the new openmp setting is correct in your app. !!!!!
                       """
            self.add_report(info=warn_msg, is_warning=True)
 
        self.add_setting(config, ["env", "openmp"], openmp)
 
        return config, self.reports


class vn90_t6288(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6288 by Paul Cresswell.
       Make steplists easier to control from the GUI."""

    BEFORE_TAG = "vn9.0_t6274"
    AFTER_TAG = "vn9.0_t6288"

    def upgrade(self, config, meta_config=None):

        config_type = self.get_setting_value(config, ["env", "config_type"])
        platform = self.get_setting_value(config,
                                          ["env", "platform_config_dir"])

        # Read and store any steplist settings:
        steplist = self.get_setting_value(config, ["env", "steplist"])
        mirror_steplist = self.get_setting_value(config, 
                                                 ["env", "mirror_steplist"])

        steps = ''
        if steplist is None:
            steps = 'extract mirror'
        else:
            steps = steplist

        if mirror_steplist is not None:
            steps = steps + ' ' + mirror_steplist


        # Set defaults for everything (dict. of (value, variable) pairs):
        atmos_steps = {
                       'preprocess_atmos': 'preprocess-atmos',
                       'preprocess_recon': 'preprocess-recon',
                       'build_atmos':      'build-atmos',
                       'build_recon':      'build-recon',
                      }
        scm_steps = {
                     'preprocess_scm': 'preprocess-scm',
                     'build_scm':      'build-scm',
                    }
        utils_serial_steps = {
                              'preprocess_makebc':    'preprocess-makebc',
                              'preprocess_pumf':      'preprocess-pumf',
                              'preprocess_cumf':      'preprocess-cumf',
                              'preprocess_frames':    'preprocess-frames',
                              'preprocess_fieldcalc': 'preprocess-fieldcalc',
                              'preprocess_vomext':    'preprocess-vomext',
                              'preprocess_fieldcos':  'preprocess-fieldcos',
                              'preprocess_convieee':  'preprocess-convieee',
                              'build_makebc':    'build-makebc',
                              'build_pumf':      'build-pumf',
                              'build_cumf':      'build-cumf',
                              'build_frames':    'build-frames',
                              'build_fieldcalc': 'build-fieldcalc',
                              'build_vomext':    'build-vomext',
                              'build_fieldcos':  'build-fieldcos',
                              'build_convieee':  'build-convieee',
                             }
        utils_serial32_steps = {
                                'preprocess_convpp': 'preprocess-convpp',
                                'preprocess_cumf':   'preprocess-cumf',
                                'preprocess_pumf':   'preprocess-pumf',
                                'build_convpp': 'build-convpp',
                                'build_cumf':   'build-cumf',
                                'build_pumf':   'build-pumf',
                               }
        utils_serialmpp_steps = {
                                 'preprocess_crmstyle_coarse_grid': 'preprocess-crmstyle_coarse_grid',
                                 'build_crmstyle_coarse_grid': 'build-crmstyle_coarse_grid',
                                }

        # Set up different configs for dictionaries:
        if config_type == 'atmos':
            d = atmos_steps
        elif config_type == 'scm':
            d = scm_steps
        elif config_type == 'utils-serial':
            d = utils_serial_steps
        elif config_type == 'utils-serial-32B':
            d = utils_serial32_steps
        elif config_type == 'utils-mpp':
            d = utils_serialmpp_steps

        """Note that single-stage builds may produce "mirror=", which may need
           changing if the app is later adapted for a different platform.

           One edge case worth attempting to correct: if a 2-stage build app
           (i.e. with an fcm_make2 task) includes $steplist but NOT 
           $mirror_steplist, we will end up with an app which extracts and/or
           mirrors (or neither) but does not build anything! To correct this we
           add the missing steps where we know we can, and warn where we're not
           sure what's required.
        """

        if steplist is not None and mirror_steplist is None:

            # Two-stage build: add the missing steps.
            if (platform == 'meto-pwr7-xlf' 
             or platform == 'nci-x86-ifort'
             or platform == 'ncas-xc30-cce'):
    
                 for key,value in d.iteritems():
                     steps = steps + ' ' + value
    
            # One-stage build: do nothing.
            elif (platform == 'meto-x86-ifort'
               or platform == 'kma-xe6-cce'):
               pass

            # Don't know - we may have disabled the remote build.
            # Leave a message to warn the user.
            else:
                warn_msg = """
     !!!!! If you use an fcm_make2 step your remote make settings may !!!!!
     !!!!! have been disabled. Please restore any missing preprocess  !!!!!
     !!!!! and build steps manually in your app using the Rose GUI.   !!!!!
                           """
                self.add_report(info=warn_msg, is_warning=True)


        # Loop over the different configs and dictionaries:
        if steplist is not None or mirror_steplist is not None:
            # For each config type check if any steps are manually specified.
            keep = ''
            for key,value in d.iteritems():

                # Remember which, if any, steps we find:
                if re.search(r"%s" % value, steps):
                    keep = keep + value + ' '

            # Now blank all except these:
            for key,value in d.iteritems():
                if not re.search(r"%s" % value, keep):
                    atmos_steps[key] = ''

        # Now add the new variables:
        dict_list = [atmos_steps, scm_steps, utils_serial_steps, utils_serial32_steps, utils_serialmpp_steps]
        for d in dict_list:
            for key,value in d.iteritems():
                self.add_setting(config, ["env", key], value)

        # Add extract and mirror to the app:
        if re.search(r'extract', steps):
            self.add_setting(config, ["env", "extract"], "extract")
        else:
            self.add_setting(config, ["env", "extract"], "")

        if re.search(r'mirror', steps):
            self.add_setting(config, ["env", "mirror"], "mirror")
        else:
            self.add_setting(config, ["env", "mirror"], "")

        # Delete the old settings:
        self.remove_setting(config, ["env", "steplist"])
        self.remove_setting(config, ["env", "mirror_steplist"])

        return config, self.reports

class vn91_t6270(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn9.0_t6288"
    AFTER_TAG = "vn9.1"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        self.change_setting_value(config, ["env", "jules_rev"],
                                  "um9.1", forced=True)
        self.change_setting_value(config, ["env", "um_rev"],
                                  "vn9.1", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$SOURCE_UM_BASE':
            pass
        else:
            config_root = 'fcm:um_tr'
            config_revision = '@vn9.1'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn9.1',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn9.1_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn9.1_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path)

        return config, self.reports


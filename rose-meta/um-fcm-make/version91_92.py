import rose.upgrade
import re
import sys

class UpgradeError(Exception):

      """Exception created when an upgrade fails."""

      def __init__(self, msg):
          self.msg = msg

      def __repr__(self):
          sys.tracebacklimit = 0
          return self.msg

      __str__ = __repr__



class vn91_t5942(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5942 by Glenn Greed."""

    BEFORE_TAG = "vn9.1"
    AFTER_TAG = "vn9.1_t5942"

    def upgrade(self, config, meta_config=None):
        """Add the build and preprocess settings for additional execs."""
        # Input your macro commands here
        steps = {
                              'preprocess_pptoanc':    'preprocess-pptoanc',
                              'preprocess_fieldop':    'preprocess-fieldop',
                              'preprocess_fieldmod':   'preprocess-fieldmod',
                              'build_pptoanc':         'build-pptoanc',
                              'build_fieldop':         'build-fieldop',
                              'build_fieldmod':        'build-fieldmod',
                             }
        for key,value in steps.iteritems():
            self.add_setting(config, ["env", key], value)
        return config, self.reports


class vn91_t6298(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6298 by Paul Cresswell."""

    BEFORE_TAG = "vn9.1_t5942"
    AFTER_TAG = "vn9.1_t6298"

    def upgrade(self, config, meta_config=None):
        # These are now compulsory=true; values match central configs:
        self.add_setting(config, ["env", "um_rev"],    "head")
        self.add_setting(config, ["env", "jules_rev"], "head")

        # Also now compulsory=true; should make it easier to apply overrides:
        self.add_setting(config, ["env", "fcflags_overrides"], "")
        self.add_setting(config, ["env", "ldflags_overrides_prefix"], "")
        self.add_setting(config, ["env", "ldflags_overrides_suffix"], "")

        # Merge preprocess- and build- buttons into one.
        execs = ['atmos', 'recon', 'scm', 'makebc', 'pumf', 'cumf', 'frames',
                 'fieldcalc', 'vomext', 'fieldcos', 'convieee', 'pptoanc',
                 'fieldop', 'fieldmod', 'convpp', 'crmstyle_coarse_grid']

        partial_compiles = []
        full_compiles = []
        # Loop through config types:
        for name in execs:
            preprocess = self.get_setting_value(config,
                         ["env", "preprocess_%s" % name])
            build = self.get_setting_value(config,
                    ["env", "build_%s" % name])
            # But do we *actually* do these things?
            do_preprocess = self.get_setting_value(config,
                            ["env", "preprocess_%s" % name], no_ignore=True)
            do_build = self.get_setting_value(config,
                       ["env", "build_%s" % name], no_ignore=True)

            # Save anything that can't be captured by the new buttons:
            if do_preprocess and not do_build:
                partial_compiles.append(preprocess)
            elif do_build and not do_preprocess:
                partial_compiles.append(build)
            elif do_preprocess and do_build:
                # Save anything else that would need to appear in the steplist
                # alongside a partial compile.
                full_compiles.append(preprocess)
                full_compiles.append(build)

            # Decide on the new value and add the new variable:
            if preprocess and build:
                value = "preprocess-%s build-%s" % (name, name)
            else:
                value = ""
            self.add_setting(config, ["env", "compile_%s" % name], value)

            self.remove_setting(config, ["env", "preprocess_%s" % name])
            self.remove_setting(config, ["env", "build_%s" % name])

        # That covers most cases. Apps where something is built OR preprocessed
        # (but not both) now need their steplists overridden manually.

        # Decide how to split this up between the local and remote steplists
        # using the same logic as vn90_t6288.
        # This is only needed if there are *compile* steps to process manually
        # (extract and mirror settings are unaffected).

        if partial_compiles:

            # The compile steps that need resolving manually are:
            remote_steps = ( ' '.join(partial_compiles) + ' ' 
                           + ' '.join(full_compiles) )

            platform = self.get_setting_value(config,
                       ["env", "platform_config_dir"])
            steplist = ''
            mirror_steplist = ''
            if (platform == 'meto-pwr7-xlf'
             or platform == 'nci-x86-ifort'
             or platform == 'ncas-xc30-cce'):
                # Two-stage build required.
                # Local steplist is handled via the extract/mirror buttons.
                mirror_steplist = remote_steps
 
            elif (platform == 'meto-x86-ifort'
               or platform == 'kma-xe6-cce'):
                # One-stage build required.
                # We have to merge local and remote steplists.

                # See if we need to extract and mirror:
                extract = self.get_setting_value(config, ["env", "extract"],
                                                 no_ignore=True)
                mirror = self.get_setting_value(config, ["env", "mirror"],
                                                no_ignore=True)
                local_steplist = []
                if extract:
                    local_steplist.append(extract)
                if mirror:
                    local_steplist.append(mirror)

                local_steps = ' '.join(local_steplist)
                steplist = local_steps + ' ' + remote_steps
 
            else:
               # Don't know - we may have disabled a remote compile sub-step.
               # Leave a message to warn the user.
                warn_msg = """
     !!!!! Your compilation (make step) settings may have changed. !!!!!
     !!!!! Please use the steplist or mirror_steplist variables to !!!!!
     !!!!! restore any missing preprocess and build steps.         !!!!!
                           """
                self.add_report(info=warn_msg, is_warning=True)
 
            if steplist:
                self.add_setting(config, ["env", "steplist"], steplist)
            if mirror_steplist:
                self.add_setting(config, ["env", "mirror_steplist"],
                                 mirror_steplist)

        return config, self.reports


class vn91_t5830(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5830 by Paul Cresswell."""

    BEFORE_TAG = "vn9.1_t6298"
    AFTER_TAG = "vn9.1_t5830"

    def upgrade(self, config, meta_config=None):
        """Add a Dr Hook compile switch."""
        self.add_setting(config, ["env", "drhook"], "drhook_off")

        return config, self.reports

class vn92_t6571(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn9.1_t5830"
    AFTER_TAG = "vn9.2"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        self.change_setting_value(config, ["env", "jules_rev"],
                                  "um9.2", forced=True)
        self.change_setting_value(config, ["env", "um_rev"],
                                  "vn9.2", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$SOURCE_UM_BASE':
            pass
        else:
            config_root = 'fcm:um_tr'
            config_revision = '@vn9.2'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn9.2',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn9.2_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn9.2_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path)

        return config, self.reports


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



class vn109_t3496(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3496 by Stuart Whitehouse."""

    BEFORE_TAG = "vn10.9"
    AFTER_TAG = "vn10.9_t3496"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""

        self.remove_setting(config, ["env", "compile_fieldcos"])
        return config, self.reports


class vn109_t3822(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3822 by Paul Cresswell."""

    BEFORE_TAG = "vn10.9_t3496"
    AFTER_TAG = "vn10.9_t3822"

    def upgrade(self, config, meta_config=None):
        """Dummy upgrade macro to force re-evaluation of triggers."""

        return config, self.reports


class vn109_t3109(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3109 by Maff Glover."""

    BEFORE_TAG = "vn10.9_t3822"
    AFTER_TAG = "vn10.9_t3109"

    def upgrade(self, config, meta_config=None):

        """Upgrade optimisation level from HIGH to FAST."""

        warn = False

        platform = self.get_setting_value(config,
                                          ['env', 'platform_config_dir'])
        config_type = self.get_setting_value(config,
                                             ['env', 'config_type'])
        level = self.get_setting_value(config,
                                       ['env', 'optimisation_level'])

        warn_on_level = ('high' in level or re.match("\$", level))

        if 'meto-xc40-cce' in platform:

            if 'atmos' in config_type:

                if 'high' in level:

                    self.change_setting_value(config,
                                        ['env', 'optimisation_level'], 'fast')

                elif re.match("\$", level):
                    warn = True

            elif re.match("\$", config_type):

                if warn_on_level:
                    warn = True

        elif re.match("\$", platform):

            if 'atmos' in config_type or re.match("\$", config_type):

                if warn_on_level:
                    warn = True

        if (warn):
            warn_msg = """
!!!                          WARNING!
!!!  Some compilation options are controlled by suite variables
!!!  and the correct upgrade path cannot be determined.
!!!  If you wish to retain the current behaviour, any fcm-make tasks that
!!!  set all three of:
!!!    - platform_config_dir  =  'meto-xc40-cce'
!!!    - config_type          =  'atmos'
!!!    - optimisation         =  'high'
!!!  should be edited to change the optimisation setting to 'fast'.
!!!  Please refer to the fcm-make metadata and the UM 11.0 release notes
!!!  for more information.
!!!
"""
            self.add_report(info=warn_msg, is_warning=True)

        return config, self.reports



class vn110_t3821(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.9_t3109"
    AFTER_TAG = "vn11.0"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        socrates_rev = self.get_setting_value(config, ["env", "socrates_rev"])
        if socrates_rev != "$BASE_SOCRATES_REV":
            self.change_setting_value(config, ["env", "socrates_rev"],
                                      "um11.0", forced=True)
        jules_rev = self.get_setting_value(config, ["env", "jules_rev"])
        if jules_rev != "$BASE_JULES_REV":
            self.change_setting_value(config, ["env", "jules_rev"],
                                      "um11.0", forced=True)
        casim_rev = self.get_setting_value(config, ["env", "casim_rev"])
        if casim_rev != "$BASE_CASIM_REV":
            self.change_setting_value(config, ["env", "casim_rev"],
                                      "um11.0", forced=True)
        um_rev = self.get_setting_value(config, ["env", "um_rev"])
        if um_rev != "$BASE_UM_REV":
            self.change_setting_value(config, ["env", "um_rev"],
                                      "vn11.0", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$HOST_SOURCE_UM_BASE':
            pass
        elif (config_root == '$CONFIG_ROOT_PATH' and
              config_revision == '$CONFIG_REVISION'):
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn11.0'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn11.0',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn11.0_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn11.0_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path, forced=True)

        return config, self.reports


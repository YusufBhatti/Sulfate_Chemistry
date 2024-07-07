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


class vn105_t2063(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2063 by Paul Cresswell."""

    BEFORE_TAG = "vn10.5"
    AFTER_TAG = "vn10.5_t2063"

    def upgrade(self, config, meta_config=None):
        """Update meto-xc40 apps for single-step builds and give a warning."""
        platform = self.get_setting_value(config, ['env', 'platform_config_dir'])
        if 'meto-xc40-' in platform:
            # Update the prebuild path automatically, if possible.
            prebuild = self.get_setting_value(config, ['env', 'prebuild'])
            if (re.match(r'\s*/home/h01/frum/', prebuild) or
                re.match(r'\s*\$UMDIR/', prebuild) or
                re.match(r'\s*~frum/', prebuild)):
                prebuild = re.sub(r'^.*/cylc-run/', '/home/d05/frum/cylc-run/',
                                  prebuild)
                self.change_setting_value(config, ['env', 'prebuild'], prebuild)
                
            warn_msg = """
  !!!                              WARNING!
  !!! From UM 10.6 onwards the central fcm-make files for the Met Office
  !!! Cray XC40 are configured for remote extract, with no mirror step.
  !!!
  !!! This means fcm_make2 tasks are no longer required for UM compilations.
  !!!
  !!! Please update your suite accordingly. See the UM release notes for
  !!! further details.
"""
            self.add_report(info=warn_msg, is_warning=True)

        return config, self.reports


class vn106_t2411(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.5_t2063"
    AFTER_TAG = "vn10.6"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        socrates_rev = self.get_setting_value(config, ["env", "socrates_rev"])
        if socrates_rev != "$BASE_SOCRATES_REV":
            self.change_setting_value(config, ["env", "socrates_rev"],
                                      "um10.6", forced=True)
        jules_rev = self.get_setting_value(config, ["env", "jules_rev"])
        if jules_rev != "$BASE_JULES_REV":
            self.change_setting_value(config, ["env", "jules_rev"],
                                      "um10.6", forced=True)
        um_rev = self.get_setting_value(config, ["env", "um_rev"])
        if um_rev != "$BASE_UM_REV":
            self.change_setting_value(config, ["env", "um_rev"],
                                      "vn10.6", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$HOST_SOURCE_UM_BASE':
            pass
        elif (config_root == '$CONFIG_ROOT_PATH' and
              config_revision == '$CONFIG_REVISION'):
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.6'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.6',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn10.6_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn10.6_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path, forced=True)

        return config, self.reports


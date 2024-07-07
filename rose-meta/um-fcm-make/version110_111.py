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


class vn110_t3921(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3921 by Steve Wardle"""

    BEFORE_TAG = "vn11.0"
    AFTER_TAG = "vn11.0_t3921"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        self.add_setting(config, ["env", "compile_wafccb_lib"],
                         "preprocess-wafccb_lib build-wafccb_lib")
        return config, self.reports

class vn110_t3382(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3382 by Michele Guidolin."""

    BEFORE_TAG = "vn11.0_t3921"
    AFTER_TAG = "vn11.0_t3382"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        # Do nothing.
        return config, self.reports


class vn110_t3815(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3815 by Steve Wardle"""

    BEFORE_TAG = "vn11.0_t3382"
    AFTER_TAG = "vn11.0_t3815"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        self.remove_setting(config, ["env", "compile_fieldcalc"])
        return config, self.reports

class vn110_t4112(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn11.0_t3815"
    AFTER_TAG = "vn11.1"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        socrates_rev = self.get_setting_value(config, ["env", "socrates_rev"])
        if socrates_rev != "$BASE_SOCRATES_REV":
            self.change_setting_value(config, ["env", "socrates_rev"],
                                      "um11.1", forced=True)
        jules_rev = self.get_setting_value(config, ["env", "jules_rev"])
        if jules_rev != "$BASE_JULES_REV":
            self.change_setting_value(config, ["env", "jules_rev"],
                                      "um11.1", forced=True)
        casim_rev = self.get_setting_value(config, ["env", "casim_rev"])
        if casim_rev != "$BASE_CASIM_REV":
            self.change_setting_value(config, ["env", "casim_rev"],
                                      "um11.1", forced=True)
        um_rev = self.get_setting_value(config, ["env", "um_rev"])
        if um_rev != "$BASE_UM_REV":
            self.change_setting_value(config, ["env", "um_rev"],
                                      "vn11.1", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$HOST_SOURCE_UM_BASE':
            pass
        elif (config_root == '$CONFIG_ROOT_PATH' and
              config_revision == '$CONFIG_REVISION'):
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn11.1'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn11.1',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn11.1_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn11.1_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path, forced=True)

        return config, self.reports


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


class vn108_t3444(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3444 by Steve Wardle"""

    BEFORE_TAG = "vn10.8"
    AFTER_TAG = "vn10.8_t3444"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        self.remove_setting(config, ["env", "compile_pumf"])
        self.remove_setting(config, ["env", "compile_cumf"])
        self.remove_setting(config, ["env", "compile_convieee"])

        config_type = self.get_setting_value(config, ["env", "config_type"])
        if config_type == "utils-serial-32B":
            raise UpgradeError("""
            Unable to upgrade app; the app specifies "env:config_type" as
            "utils-serial-32B" but there are no longer any build configs
            for 32-bit serial utilities. Please either remove the app or
            modify its "env:config_type" to specify a 64-bit serial build
            ("utils-serial") instead.
            """)

        return config, self.reports

class vn108_t3355(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3355 by Richard Gilham"""

    BEFORE_TAG = "vn10.8_t3444"
    AFTER_TAG = "vn10.8_t3355"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        self.add_setting(config, ["env", "ls_precipitation_precision"],"double")
        return config, self.reports

class vn109_t3471(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.8_t3355"
    AFTER_TAG = "vn10.9"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        socrates_rev = self.get_setting_value(config, ["env", "socrates_rev"])
        if socrates_rev != "$BASE_SOCRATES_REV":
            self.change_setting_value(config, ["env", "socrates_rev"],
                                      "um10.9", forced=True)
        jules_rev = self.get_setting_value(config, ["env", "jules_rev"])
        if jules_rev != "$BASE_JULES_REV":
            self.change_setting_value(config, ["env", "jules_rev"],
                                      "um10.9", forced=True)
        casim_rev = self.get_setting_value(config, ["env", "casim_rev"])
        if casim_rev != "$BASE_CASIM_REV":
            self.change_setting_value(config, ["env", "casim_rev"],
                                      "um10.9", forced=True)
        um_rev = self.get_setting_value(config, ["env", "um_rev"])
        if um_rev != "$BASE_UM_REV":
            self.change_setting_value(config, ["env", "um_rev"],
                                      "vn10.9", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$HOST_SOURCE_UM_BASE':
            pass
        elif (config_root == '$CONFIG_ROOT_PATH' and
              config_revision == '$CONFIG_REVISION'):
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.9'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.9',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn10.9_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn10.9_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path, forced=True)

        return config, self.reports


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



class vn104_t1118(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1118 by Stuart Whitehouse."""

    BEFORE_TAG = "vn10.4"
    AFTER_TAG = "vn10.4_t1118"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        # Input your macro commands here
        self.remove_setting(config,["env","compile_makebc"])
        self.remove_setting(config,["env","compile_frames"])
        return config, self.reports


class vn104_t1636(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1636 by Stuart Whitehouse."""

    BEFORE_TAG = "vn10.4_t1118"
    AFTER_TAG = "vn10.4.1"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""

        um_rev = self.get_setting_value(config, ["env", "um_rev"])
        if um_rev != "$BASE_UM_REV":
            self.change_setting_value(config, ["env", "um_rev"],
                                      "vn10.4.1", forced=True)

        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$SOURCE_UM_BASE':
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.4.1'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.4.1',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn10.4.1_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn10.4.1_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path, forced=True)

        return config, self.reports


class vn105_t1963(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.4.1"
    AFTER_TAG = "vn10.5"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        socrates_rev = self.get_setting_value(config, ["env", "socrates_rev"])
        if socrates_rev != "$BASE_SOCRATES_REV":
            self.change_setting_value(config, ["env", "socrates_rev"],
                                      "um10.5", forced=True)
        jules_rev = self.get_setting_value(config, ["env", "jules_rev"])
        if jules_rev != "$BASE_JULES_REV":
            self.change_setting_value(config, ["env", "jules_rev"],
                                      "um10.5", forced=True)
        um_rev = self.get_setting_value(config, ["env", "um_rev"])
        if um_rev != "$BASE_UM_REV":
            self.change_setting_value(config, ["env", "um_rev"],
                                      "vn10.5", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$HOST_SOURCE_UM_BASE':
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.5'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.5',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn10.5_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn10.5_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path, forced=True)

        return config, self.reports


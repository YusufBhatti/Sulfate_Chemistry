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



class vn102_t759(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #759 by Steve Wardle"""

    BEFORE_TAG = "vn10.2"
    AFTER_TAG = "vn10.2_t759"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""

        # Remove the mule setting (added last release, but a more elegant
        # method for building the library then emerged)
        self.remove_setting(config, ['env', 'mule'])

        return config, self.reports


class vn103_t1119(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.2_t759"
    AFTER_TAG = "vn10.3"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        self.change_setting_value(config, ["env", "socrates_rev"],
                                  "um10.3", forced=True)
        self.change_setting_value(config, ["env", "jules_rev"],
                                  "um10.3", forced=True)
        self.change_setting_value(config, ["env", "um_rev"],
                                  "vn10.3", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$SOURCE_UM_BASE':
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.3'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.3',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn10.3_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn10.3_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path)

        return config, self.reports


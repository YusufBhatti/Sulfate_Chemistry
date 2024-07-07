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



class vn100_t2(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2 by Stuart Whitehouse."""

    BEFORE_TAG = "vn10.0"
    AFTER_TAG = "vn10.0_t2"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration to include CreateBC."""
        # In all cases other than CreateBC, this should be trigger ignored,
        # and no CreateBC builds exist as they're introduced in this ticket.

        self.add_setting(config, ["env", "compile_createbc"] ,
                     "preprocess-createbc build-createbc")

        return config, self.reports


class vn101_t327(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.0_t2"
    AFTER_TAG = "vn10.1"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        self.change_setting_value(config, ["env", "jules_rev"],
                                  "um10.1", forced=True)
        self.change_setting_value(config, ["env", "um_rev"],
                                  "vn10.1", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$SOURCE_UM_BASE':
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.1'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.1',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn10.1_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn10.1_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path)

        return config, self.reports


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


# !!! Note:
#     This macro is actually a duplicate of an earlier macro from the vn91-92
#     upgrade section.  A subtle bug caused things to conflict with an earlier
#     macro leaving apps in a partially upgraded state with some orphaned nodes
#
#     The ticket cited here fixes this bug in the original macro, but in cases
#     of user apps which already contain the introduced problems, this macro
#     aims to fix their apps.  In order to do this it sub-classes the original
#     conflicting macro (so it will re-run that macro again)
#
#     This should **not** be followed as an example for future macros unless
#     for a very, very good reason
#
from version91_92 import vn91_t6297
class vn92_t463(vn91_t6297):

    """Upgrade macro for ticket #463 by Steve Wardle."""

    BEFORE_TAG = "vn9.2"
    AFTER_TAG = "vn9.2_t463"

    # Intentionally blank (see the comment above)

class vn92_t463a(rose.upgrade.MacroUpgrade):

    """2nd upgrade macro for ticket #463 by Steve Wardle."""

    BEFORE_TAG = "vn9.2_t463"
    AFTER_TAG = "vn9.2_t463a"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Remove the idealise namelist entry (if it existed)
        self.remove_setting(config, ["namelist:nlcfiles", "idealise"])

        return config, self.reports


class vn92_t1389(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1389 (SRS) by Paul Cresswell."""

    BEFORE_TAG = "vn9.2_t463a"
    AFTER_TAG = "vn9.2_t1389"

    def upgrade(self, config, meta_config=None):
        """Add essential ioscntl inputs to all UM apps."""
        self.add_setting(config, ["namelist:ioscntl", "ios_spacing"], "0")
        self.add_setting(config, ["namelist:ioscntl", "ios_offset"], "0")
        self.add_setting(config, ["namelist:ioscntl", "ios_tasks_per_server"],
                         "1")

        # Fix invalid values of ios_spacing:
        ios_spacing = self.get_setting_value(config,
                                 ["namelist:ioscntl", "ios_spacing"])
        if int(ios_spacing) < 0:
           self.change_setting_value(config,
                                 ["namelist:ioscntl", "ios_spacing"], "0")

        return config, self.reports


class vn100_t40(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn9.2_t1389"
    AFTER_TAG = "vn10.0"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "10.0")
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMSTR"])
        if stashmaster_path:
           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",
                                     "vn10.0/ctldata",
                                      stashmaster_path)
           self.change_setting_value(config, ["env", "STASHMSTR"],
                                     stashmaster_path)
        return config, self.reports


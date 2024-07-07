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



class vn103_t1330(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1330 by Stuart Whitehouse."""

    BEFORE_TAG = "vn10.3"
    AFTER_TAG = "vn10.3_t1330"

    def upgrade(self, config, meta_config=None):
        """Upgrade a CreateBC app configuration."""
        self.add_setting(config, ["namelist:lbc_grid", 
            "write_header_only_once"], ".false.")
        return config, self.reports

class vn103_t1292(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1292 by Stuart Whitehouse."""

    BEFORE_TAG = "vn10.3_t1330"
    AFTER_TAG = "vn10.3_t1292"

    def upgrade(self, config, meta_config=None):
        """Upgrade a CreateBC app configuration."""
        self.add_setting(config, ["namelist:lbc_grid", 
                                  "frames_packing_option"], "-1")

        return config, self.reports

class vn104_t1450(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.3_t1292"
    AFTER_TAG = "vn10.4"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""

        return config, self.reports


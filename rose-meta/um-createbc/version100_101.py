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



class vn100_tXXXX(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #XXXX by <author>."""

    BEFORE_TAG = "vn10.0_t2"
    AFTER_TAG = "vn10.0_tXXXX"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        return config, self.reports


class vn101_t327(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.0_tXXXX"
    AFTER_TAG = "vn10.1"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""

        return config, self.reports


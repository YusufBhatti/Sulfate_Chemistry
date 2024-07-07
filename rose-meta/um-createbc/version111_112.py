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



class vn111_tXXXX(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #XXXX by <author>."""

    BEFORE_TAG = "vn11.1"
    AFTER_TAG = "vn11.1_tXXXX"

    def upgrade(self, config, meta_config=None):
        """Upgrade a CreateBC app configuration."""
        # Input your macro commands here
        return config, self.reports


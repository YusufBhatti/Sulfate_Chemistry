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



class vn110_t4112(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn11.0"
    AFTER_TAG = "vn11.1"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""

        return config, self.reports


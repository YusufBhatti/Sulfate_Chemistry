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


class vn107_t2998(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2998 by Adam Voysey (adamvoysey)."""

    BEFORE_TAG = "vn10.7"
    AFTER_TAG = "vn10.7_t2998"

    def upgrade(self, config, meta_config=None):
        """Upgrade a CreateBC app configuration."""
        self.add_setting(config, ["env", "DR_HOOK"], "false")
        self.add_setting(config, ["env", "DR_HOOK_CATCH_SIGNALS"], "0")
        self.add_setting(config, ["env", "DR_HOOK_IGNORE_SIGNALS"], "0")
        self.add_setting(config, ["env", "DR_HOOK_OPT"], "self,wallprof")
        self.add_setting(config, ["env", "DR_HOOK_PROFILE"],
                         "$CYLC_TASK_WORK_DIR/drhook.prof.%d")
        self.add_setting(config, ["env", "DR_HOOK_PROFILE_LIMIT"], "-10.0")
        self.add_setting(config, ["env", "DR_HOOK_PROFILE_PROC"], "-1")
        return config, self.reports

class vn108_t3162(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.7_t2998"
    AFTER_TAG = "vn10.8"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""

        return config, self.reports


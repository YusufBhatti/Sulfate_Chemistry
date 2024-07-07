import rose.upgrade
import re
import sys


class vn106_t2411(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.5"
    AFTER_TAG = "vn10.6"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""

        return config, self.reports


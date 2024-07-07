import rose.upgrade
import re
import sys


class vn105_t1963(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.4"
    AFTER_TAG = "vn10.5"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""

        return config, self.reports


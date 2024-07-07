import rose.upgrade
import fileinput
import os
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



class vn101_t464(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #464 by James Manners."""

    BEFORE_TAG = "vn10.1"
    AFTER_TAG = "vn10.1_t464"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        # Add SOCRATES source:
        self.add_setting(config, ['env', 'socrates_rev'] , 'um10.1')
        jules_sources = self.get_setting_value(config, [ 'env', 'jules_sources'])
        if jules_sources == '$SOURCE_JULES':
            self.add_setting(config, ['env', 'socrates_sources'] ,'$SOURCE_SOCRATES')
        else:
            self.add_setting(config, ['env', 'socrates_sources'], '')

        # Edit the fcm-make.cfg file to add an extract statement:
        filename = './file/fcm-make.cfg'
        # Make sure the fcm-make.cfg file is available for upgrading:
        if not os.path.exists(filename):
            raise UpgradeError('Cannot find file/fcm-make.cfg, please ensure the fcm_make app contains the necessary config file')

        # Insert the extra line:
        for line in fileinput.input([filename], inplace=True):
            if 'jules_sources' in line:
                line = line.rstrip() + '\nextract.location{diff}[socrates] = $socrates_sources\n'
            sys.stdout.write(line)

        # No means yet of reporting changes in file/, so we provide our own:
        print '''\
[U] File upgrade %s-%s: 
    file=fcm-make.cfg
        Added text 'extract.location{diff}[socrates] = $socrates_sources'\
''' % (self.BEFORE_TAG, self.AFTER_TAG)

        return config, self.reports


class vn101_t691(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #691 by Stuart Whitehouse."""

    BEFORE_TAG = "vn10.1_t464"
    AFTER_TAG = "vn10.1_t691"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        else:
           prebuild_path = re.sub(r'fcm_make_metolinux', r'fcm_make_meto_linux_ifort', prebuild_path)
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path)

        return config, self.reports


class vn101_t70(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #70 by Steve Wardle"""

    BEFORE_TAG = "vn10.1_t691"
    AFTER_TAG = "vn10.1_t70"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        # Add new compilation option for the packing library
        self.add_setting(config, ['env', 'compile_packing_lib'],
                         'preprocess-packing_lib build-packing_lib')
        # Add setting which triggers the inclusion of mule scripts
        self.add_setting(config, ['env', 'mule'], 'mule_off')

        return config, self.reports

class vn102_t718(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.1_t70"
    AFTER_TAG = "vn10.2"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        self.change_setting_value(config, ["env", "socrates_rev"],
                                  "um10.2", forced=True)
        self.change_setting_value(config, ["env", "jules_rev"],
                                  "um10.2", forced=True)
        self.change_setting_value(config, ["env", "um_rev"],
                                  "vn10.2", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$SOURCE_UM_BASE':
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.2'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.2',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn10.2_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn10.2_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path)

        return config, self.reports


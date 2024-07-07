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


class vn86_t5536(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5536 by Paul Cresswell.
       Add the ability to use prebuilds in Rose."""
    
    BEFORE_TAG = "vn8.6"
    AFTER_TAG = "vn8.6_t5536"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM fcm-make app configuration."""
        # Add a prebuild variable (left blank for the user to select one):
        self.add_setting(config, ["env", "prebuild"], "")

        # Edit the fcm-make.cfg file to add the necessary 'use' statement:

        filename = './file/fcm-make.cfg'
        # Make sure the fcm-make.cfg file is available for upgrading:
        if not os.path.exists(filename):
            raise UpgradeError('Cannot find file/fcm-make.cfg, please apply this upgrade macro directly to the fcm_make app')

        # Insert the extra line:
        for line in fileinput.input([filename], inplace=True):
            if fileinput.isfirstline():
                line = 'use = $prebuild\n\n' + line
            sys.stdout.write(line)

        # No means yet of reporting changes in file/, so we provide our own:
        print '''\
[U] File upgrade %s-%s: 
    file=fcm-make.cfg
        Added text 'use = $prebuild'\
''' % (self.BEFORE_TAG, self.AFTER_TAG)

        return config, self.reports

class vn90_t5995(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn8.6_t5536"
    AFTER_TAG = "vn9.0"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        um_version = self.get_setting_value(config, ["env", "um_rev"])
        if um_version:
            if um_version == "vn8.5" or um_version == "vn8.4":
                raise UpgradeError(
                    '!!!! Upgrade of rose app from UM version %s' % um_version  +
                    ' not supported. !!!!\nUpgrade macros are only available ' +
                    'from vn8.6 onwards.\nPlease upgrade equivalent job to '+
                    'vn8.6 using the UMUI. This may then be converted to Rose '+
                    'and upgraded to vn9.0.\n For simple apps it may be possible'+
                    ' to manually upgrade the app by hand.'
                    )
        self.change_setting_value(config, ["env", "jules_rev"],
                                  "um9.0", forced=True)
        self.change_setting_value(config, ["env", "um_rev"],
                                  "vn9.0", forced=True)

        include = self.get_setting_value(config, ["env", "include_config"])
        if include.endswith('$SOURCE_UM_REV'):
           # Do nothing - this is a rose-stem-style fcm_make app
           pass
        else:
            include = re.sub(r'.*/fcm-make', r'fcm:um_tr/fcm-make', include)
            include = re.sub(r'@.*', r'@vn9.0', include)
            if '@' not in include:
                include = include + '@vn9.0'
            self.add_report("env", "include_config", include, 
             info="Upgrading fcm_make config version to trunk@vn9.0",
             is_warning=True)

        self.change_setting_value(config, ["env", "include_config"],
                                      include)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn9.0_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn9.0_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path)

        return config, self.reports


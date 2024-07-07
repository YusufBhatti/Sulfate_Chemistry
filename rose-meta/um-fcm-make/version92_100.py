import rose.upgrade
import re
import sys
from collections import defaultdict

class UpgradeError(Exception):

      """Exception created when an upgrade fails."""

      def __init__(self, msg):
          self.msg = msg

      def __repr__(self):
          sys.tracebacklimit = 0
          return self.msg

      __str__ = __repr__



class vn92_t1(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1 by Paul Cresswell."""

    BEFORE_TAG = "vn9.2"
    AFTER_TAG = "vn9.2.1"

    def make_relative(self, project, sources):
        """Convert standard absolute paths into relative paths"""
        # Case-insensitive keywords:
        iproject = "(" + project + ")(?i)"
        sources = re.sub(r'fcm:'+iproject+'[-_](br)(?i)', 'branches', sources)
        sources = re.sub(r'fcm:'+iproject+'/branches',    'branches', sources)
        sources = re.sub(r'fcm:'+iproject+'[-_](tr)(?i)', 'trunk', sources)
        sources = re.sub(r'fcm:'+iproject+'/trunk',       'trunk', sources)
        # Met Office only:
        sources = re.sub(r'svn://fcm\d/'+project+'_svn/'+project+'/', '',
                         sources)

        return sources

    def upgrade(self, config, meta_config=None):
        """Redirect fcm_make app to build from the mirrored configs"""
        warn_msg = ''
        config_root_path = self.get_setting_value(config, 
                           ["env", "config_root_path"])
        config_revision = self.get_setting_value(config, 
                           ["env", "config_revision"])

        if re.match(r'\s*fcm:(um)(?i)[-_/]', config_root_path):
            config_root_path=re.sub(r'\s*fcm:(um)(?i)', 'fcm:um.xm',
                                    config_root_path)
            if not re.match(r'@(head)(?i)', config_revision):
                # Presumably not a branch, so:
                config_revision = '@vn9.2.1'

        # Met Office only:
        elif re.match(r'\s*svn://fcm\d/UM_svn/UM', config_root_path):
            config_root_path=re.sub(r'\s*svn://fcm\d/UM_svn/UM','fcm:um.xm',
                                    config_root_path)
            if not re.match(r'@(head)(?i)', config_revision):
                # Presumably not a branch, so:
                config_revision = '@vn9.2.1'

        # Last resort - set it to use the mirrored trunk.
        # If users need to edit this to use a branch, at least that's easy.
        else:
            if config_root_path != '$SOURCE_UM_BASE':
                config_root_path = 'fcm:um.xm_tr'
                config_revision = '@vn9.2.1'
                warn_msg = """
  !!!!! WARNING: config_root_path has been reset to use the new trunk.
                           """
            

        # Does nothing if the value hasn't changed:
        self.change_setting_value(config, ["env", "config_root_path"],
                                  config_root_path)
        self.change_setting_value(config, ["env", "config_revision"],
                                  config_revision)

        # The new variables *_project_location mean we can simplify base and
        # source variables to use relative paths, without worrying about
        # inter-dependencies.

        strings = [('UM', 'um_base'),       ('UM', 'um_sources'),
                   ('JULES', 'jules_base'), ('JULES', 'jules_sources')]
        extracts = defaultdict(list)
        for key, value in strings:
            extracts[key].append(value)

        for project in extracts:
            basename   = extracts[project][0]   # e.g. 'um_base'
            sourcename = extracts[project][1]   # e.g. 'um_sources'

            # Begin updating extract variables:
            # (*_base are compulsory=false variables)
            baseval = self.get_setting_value(config, ["env", basename])

            if baseval is not None:

                if (re.match(r'\s*fcm:('+project+'[-_/]tr)(?i)', baseval)
                 or re.match(r'\s*svn://fcm\d/'+project+'_svn/'+project+'/trunk',
                             baseval)):
                    # trunk is the default, remove:
                    self.remove_setting(config, ["env", basename])
    
                elif (re.match(r'\s*(fcm:'+project+'[-_/]br)(?i)', baseval)
                   or re.match(r'\s*svn://fcm\d/'+project+'_svn/'+project+'/branches',
                               baseval)):
                    # Keep branch, remove project:
                    baseval = self.make_relative(project, baseval)
                    self.change_setting_value(config, ["env",basename], baseval)
    
                elif re.match(r'\s*(fcm:|svn:|/)', baseval):
                    # Any remaining match is an unexpected keyword, URL or
                    # filepath. Don't second-guess; simply warn the user:
                    warn_msg = warn_msg + """
  !!!!! WARNING: %s may require attention.
                                          """ % basename

            # Allow everything to be specified by a relative path:
            # (*_sources are compulsory=true variables)
            sourcelist = self.get_setting_value(config, ["env", sourcename])
            sourcelist = self.make_relative(project, sourcelist)
            self.change_setting_value(config, ["env", sourcename], sourcelist)

            # Are there any non-standard URLs left that might need fixing?
            if re.search(r'(^| )(fcm:|svn:|/)', sourcelist):
                warn_msg = warn_msg + """
  !!!!! WARNING: %s may require attention.
                                      """ % sourcename

        # Print a general warning:
        warn_msg = """
  !!!!! All extract paths (um_base, um_sources, jules_base, jules_sources),
  !!!!! unless including a full path to a project or working copy, are now 
  !!!!! relative to the values of um_project_location & jules_project_location.
  !!!!! A relative path is one specified relative to the project location URL,
  !!!!! e.g. 
  !!!!! um_sources=branches/dev/... instead of um_sources=fcm:um.xm_br/dev/...
  !!!!! You are advised to use relative paths where possible.
  !!!!!
  !!!!! Additionally, any value of um_base or jules_base referring to that
  !!!!! project's trunk can be removed, as this is the default in all cases.
                   """ + warn_msg
        self.add_report(info=warn_msg, is_warning=True)

        return config, self.reports


class vn92_t18(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #18 by Paul Cresswell."""

    BEFORE_TAG = "vn9.2.1"
    AFTER_TAG = "vn9.2.2"

    def upgrade(self, config, meta_config=None):
        """Update the config rev to vn9.2.2, using external JULES"""
        config_revision = self.get_setting_value(config,
                           ["env", "config_revision"])

        # Replicate the structure of the vn9.2.1 macro above:
        if config_revision == '$SOURCE_UM_REV':
            # A rose-stem suite
            pass

        elif re.match(r'@(head)(?i)', config_revision):
            # Most likely a branch
            pass

        else:
            # Everything else
            self.change_setting_value(config, ["env", "config_revision"],
                                      "vn9.2.2")

        return config, self.reports
    

class vn100_t40(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn9.2.2"
    AFTER_TAG = "vn10.0"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        self.change_setting_value(config, ["env", "jules_rev"],
                                  "um10.0", forced=True)
        self.change_setting_value(config, ["env", "um_rev"],
                                  "vn10.0", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$SOURCE_UM_BASE':
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.0'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.0',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn10.0_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn10.0_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path)

        return config, self.reports


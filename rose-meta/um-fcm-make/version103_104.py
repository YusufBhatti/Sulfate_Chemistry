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



class vn103_t1199(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1199."""

    BEFORE_TAG = "vn10.3"
    AFTER_TAG = "vn10.3.1"

    def upgrade(self, config, meta_config=None):
        """Upgrade meto* fcm_make apps to vn10.3.1"""
        platform = self.get_setting_value(config, 
                                          ['env', 'platform_config_dir'])
        # Remove retired variables
        if platform == "meto-x86-ifort":
            self.remove_setting(config, ['env', 'diag_disable'])
            self.remove_setting(config, ['env', 'diag_enable'])
            self.remove_setting(config, ['env', 'diag_error'])

        # Upgrade MO platforms to next minor version
        if re.match('meto-', platform):

            self.change_setting_value(config, ["env", "um_rev"],
                                      "vn10.3.1", forced=True)

            config_revision = self.get_setting_value(config, 
                                      ['env', 'config_revision'])
            config_root = self.get_setting_value(config,
                                      ['env', 'config_root_path'])
            if config_root == '$SOURCE_UM_BASE':
                pass
            else:
                config_root = 'fcm:um.xm_tr'
                config_revision = '@vn10.3.1'
                self.add_report('env', 'config_root_path', config_root,
                    info='Upgrading fcm_make config version to trunk@vn10.3.1',
                    is_warning=True)
            self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)
            self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)

            # Update prebuilds:
            prebuild = self.get_setting_value(config,['env', 'prebuild'])
            if prebuild == r'$PREBUILD':
                pass
            elif re.search('/vn10.3_prebuilds', prebuild):
                prebuild = re.sub('/vn10.3_prebuilds',
                                  '/vn10.3.1_prebuilds', prebuild)
            elif re.search(r'/r\d+_', prebuild):
                prebuild = re.sub(r'/r\d+_', '/vn10.3.1_', prebuild)
            else:
                prebuild = ''
            self.change_setting_value(config, ['env', 'prebuild'], prebuild)

        return config, self.reports


class vn103_t1238(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1238 by Paul Cresswell."""

    BEFORE_TAG = "vn10.3.1"
    AFTER_TAG = "vn10.3_t1238"

    def upgrade(self, config, meta_config=None):
        """Replace CPP def strings with gui-friendly switches."""
        # Set defaults:
        coupler            = 'none'
        stash_version      = '1A'
        portio_version     = '2A'
        mpp_version        = '1C'
        timer_version      = '3A'
        solver_precision   = 'double'
        land_surface_model = 'jules'
        recon_mpi_model    = 'parallel'

        # Builds with different defaults:
        config_type = self.get_setting_value(config, ["env", "config_type"])
        if config_type == "scm":
            timer_version  = "4A"

        # Get all existing keys:
        keys = ['keys_atmos_app', 'keys_recon_app', 'keys_createbc_app',
                'keys_scm_app', 'keys_crmstyle_coarse_grid_app',
                'keys_makebc_app', 'keys_cumf_app', 'keys_pumf_app',
                'keys_frames_app', 'keys_fieldcalc_app', 'keys_vomext_app',
                'keys_fieldcos_app', 'keys_vomext_app', 'keys_fieldcos_app',
                'keys_convieee_app', 'keys_pptoanc_app', 'keys_fieldop',
                'keys_fieldmod', 'keys_convpp_app']

        for key in keys:
            #  Get the key values:
            keylist = self.get_setting_value(config, ["env", key])

            if keylist is not None:

                # Parse each string to set the new values:
                if re.search('C95_2A', keylist):
                    portio_version = '2A'
                elif re.search('C95_2B', keylist):
                    portio_version = '2B'
                if re.search('C97_1A', keylist):
                    timer_version = '1A'
                elif re.search('C97_3A', keylist):
                    timer_version = '3A'
                elif re.search('C97_4A', keylist):
                    timer_version = '4A'
    
                # Special/simple cases:
                if key == 'keys_atmos_app':
                    if (re.search('OASIS3=oasis3',keylist)
                    and re.search('MCT=mct',keylist)):
                        coupler = 'oasis3_mct'
                    elif re.search('OASIS3=oasis3', keylist):
                        coupler = 'oasis3'
    
                    if not re.search('C_DP_HLM=c_dp_hlm', keylist):
                        solver_precision = 'single'
    
                elif key == 'keys_recon_app':
                    if re.search('RECON_SERIAL', keylist):
                        recon_mpi_model = 'serial'

        # Update the config:
        self.add_setting(config, ["env", "coupler"], coupler)
        self.add_setting(config, ["env", "stash_version"], stash_version)
        self.add_setting(config, ["env", "portio_version"], portio_version)
        self.add_setting(config, ["env", "mpp_version"], mpp_version)
        self.add_setting(config, ["env", "timer_version"], timer_version)
        self.add_setting(config, ["env", "solver_precision"], solver_precision)
        self.add_setting(config, ["env", "land_surface_model"],
                         land_surface_model)
        self.add_setting(config, ["env", "recon_mpi_model"], recon_mpi_model)

        for key in keys:
            self.remove_setting(config, ["env", key])

        return config, self.reports


class vn103_t1314(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1314 by Paul Cresswell."""

    BEFORE_TAG = "vn10.3_t1238"
    AFTER_TAG = "vn10.3_t1314"

    def upgrade(self, config, meta_config=None):
        """Update values of OpenMP and Dr Hook variables."""
        openmp = self.get_setting_value(config, ["env", "openmp"])
        if openmp == 'openmp_off':
            self.change_setting_value(config, ["env", "openmp"], "false")
        elif openmp == 'openmp_on':
            self.change_setting_value(config, ["env", "openmp"], "true")
        else:
            # Cater for $OPENMP settings:
            warn_openmp = """
        !!! The value of the $openmp variable has changed as follows: !!!
        !!!   openmp_on   ->  true                                    !!!
        !!!   openmp_off  ->  false                                   !!!
        !!! Please update your suite with the new settings.           !!!
            """
            self.add_report(info=warn_openmp, is_warning=True)

        drhook = self.get_setting_value(config, ["env", "drhook"])
        self.remove_setting(config, ["env", "drhook"])
        if drhook == 'drhook_off':
            self.add_setting(config, ["env", "DR_HOOK"], "false")
        elif drhook == 'drhook_on':
            self.add_setting(config, ["env", "DR_HOOK"], "true")
        else:
            # Cater for $DRHOOK settings:
            self.add_setting(config, ["env", "DR_HOOK"], "$DRHOOK")
            warn_drhook = """
        !!! The value of the $drhook variable has changed as follows: !!!
        !!!   drhook_on   ->  true                                    !!!
        !!!   drhook_off  ->  false                                   !!!
        !!! Please update your suite with the new settings.           !!!
        !!! Additionally, please note that the variable name has      !!!
        !!! changed from $drhook to $DR_HOOK in the fcm-make app.     !!!
            """
            self.add_report(info=warn_drhook, is_warning=True)

        return config, self.reports


class vn103_t1556(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1556 by Paul Cresswell."""

    BEFORE_TAG = "vn10.3_t1314"
    AFTER_TAG = "vn10.3_t1556"

    def upgrade(self, config, meta_config=None):
        """Update variable names."""
        renames = ['coupler', 'atmos_exec', 'recon_exec', 'scm_exec',
                   'createbc_exec', 'crmstyle_coarse_grid_exec',
                   'convpp_exec', 'cumf_exec', 'pumf_exec', 'convieee_exec',
                   'fieldcalc_exec', 'fieldcos_exec', 'frames_exec',
                   'makebc_exec', 'vomext_exec', 'pptoanc_exec',
                   'fieldop_exec', 'fieldmod_exec' ]

        for var in renames:
            self.rename_setting(config, ["env", var], ["env", var.upper()])

        self.rename_setting(config, ["env", "recon_mpi_model"],
                                    ["env", "recon_mpi"])

        return config, self.reports

class vn104_t1450(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.3_t1556"
    AFTER_TAG = "vn10.4"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        self.change_setting_value(config, ["env", "socrates_rev"],
                                  "um10.4", forced=True)
        self.change_setting_value(config, ["env", "jules_rev"],
                                  "um10.4", forced=True)
        self.change_setting_value(config, ["env", "um_rev"],
                                  "vn10.4", forced=True)
        config_revision = self.get_setting_value(config, ['env', 'config_revision'])
        config_root = self.get_setting_value(config, ['env', 'config_root_path'])
        if config_root == '$SOURCE_UM_BASE':
            pass
        else:
            config_root = 'fcm:um.xm_tr'
            config_revision = '@vn10.4'
            self.add_report('env', 'config_root_path', config_root, 
                info='Upgrading fcm_make config version to trunk@vn10.4',
                is_warning=True)
        self.change_setting_value(config, ['env', 'config_revision'], config_revision, forced=True)
        self.change_setting_value(config, ['env', 'config_root_path'], config_root, forced=True)

        prebuild_path = self.get_setting_value(config, ['env', 'prebuild'])
        if prebuild_path == r'$PREBUILD':
            pass
        elif re.search(r'/vn\d+\.\d+_', prebuild_path):
            prebuild_path = re.sub(r'/vn\d+\.\d+_', r'/vn10.4_', prebuild_path)
        elif re.search(r'/r\d+_', prebuild_path):
            prebuild_path = re.sub(r'/r\d+_', r'/vn10.4_', prebuild_path)
        else:
            prebuild_path = ''
        self.change_setting_value(config, ['env', 'prebuild'], prebuild_path, forced=True)

        return config, self.reports


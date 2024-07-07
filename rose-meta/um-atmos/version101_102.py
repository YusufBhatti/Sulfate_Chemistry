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


class vn101_t384(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #384 by Roddy Sharp."""

    BEFORE_TAG = "vn10.1"
    AFTER_TAG = "vn10.1_t384"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Get setting values for select_input_fields from
        # recon which will be renammed select_output_fields
        select_input_fields_value = self.get_setting_value(config,
                                                  ["namelist:recon", "select_input_fields"])

        # remove the existing "select_input_fields" entry
        self.remove_setting(config,["namelist:recon", "select_input_fields"] )

        # add the replacement "select_output_fields" entry
        self.add_setting(config, ["namelist:recon", "select_output_fields"],
                         value=select_input_fields_value, forced=True)

        return config, self.reports

class vn101_t420(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #420 by Steve Wardle."""

    BEFORE_TAG = "vn10.1_t384"
    AFTER_TAG = "vn10.1_t420"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        msg = """
              !!!!! WARNING: You appear to have set the "last_field" or   !!!!!
              !!!!! "partial" value in the nlstcall_pp namelist/s:        !!!!!
              !!!!!  *  {0:s}
              !!!!! These values should not normally be set by the user   !!!!!
              !!!!! since they are part of the history mechanism, and as  !!!!!
              !!!!! of UM ticket #420 they are passed via a different     !!!!!
              !!!!! (history only) namelist.  They are being REMOVED from !!!!!
              !!!!! the app, please check your configuration carefully.   !!!!!
              """

        # Store the affected namelists in this list
        affected_namelists = []

        # Find all nlstcall_pp namelists in the app
        for section in config.value.keys():
            node = config.get([section])
            if (isinstance(node.value, dict) and
                re.match("namelist:nlstcall_pp\(.*?\)$", section)):

                # Check if the namelist contained the two variables which the
                # ticket aims to remove
                last_field = self.get_setting_value(config,
                                                    [section, "last_field"])
                partial    = self.get_setting_value(config,
                                                    [section, "partial"])
                if last_field is not None or partial is not None:
                    affected_namelists.append(section)

                # Remove the settings
                self.remove_setting(config, [section, "last_field"])
                self.remove_setting(config, [section, "partial"])

        # Issue the above warning for the user if any of their namelists were
        # affected, quoting any section names
        if len(affected_namelists) > 0:
            print msg.format(
                "\n              !!!!!  *  ".join(affected_namelists))

        return config, self.reports

# !!! Note:
#     This macro is actually a duplicate of an earlier macro from the vn91-92
#     upgrade section.  A subtle bug caused things to conflict with an earlier
#     macro leaving apps in a partially upgraded state with some orphaned nodes
#
#     The ticket cited here fixes this bug in the original macro, but in cases
#     of user apps which already contain the introduced problems, this macro
#     aims to fix their apps.  In order to do this it sub-classes the original
#     conflicting macro (so it will re-run that macro again)
#
#     This should **not** be followed as an example for future macros unless
#     for a very, very good reason
#
from version91_92 import vn91_t6297
class vn101_t463(vn91_t6297):

    """Upgrade macro for ticket #463 by Steve Wardle."""

    BEFORE_TAG = "vn10.1_t420"
    AFTER_TAG = "vn10.1_t463"

    # Intentionally blank (see the comment above)

class vn101_t463a(rose.upgrade.MacroUpgrade):

    """2nd upgrade macro for ticket #463 by Steve Wardle."""

    BEFORE_TAG = "vn10.1_t463"
    AFTER_TAG = "vn10.1_t463a"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Remove the idealise namelist entry (if it existed)
        self.remove_setting(config, ["namelist:nlcfiles", "idealise"])

        return config, self.reports


class vn101_t349(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #349 by <Warren Tennant>."""

    BEFORE_TAG = "vn10.1_t463a"
    AFTER_TAG = "vn10.1_t349"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Split tau into separate setting for SKEB and SPT
        tau=self.get_setting_value(config, ["namelist:run_stochastic","tau"])
        self.remove_setting(config, ["namelist:run_stochastic", "tau"])
        self.add_setting(config, ["namelist:run_stochastic", "tau_skeb"], tau)
        self.add_setting(config, ["namelist:run_stochastic", "tau_spt"], tau)
        # Re-name n1 and n2 to something more descriptive
        n1=self.get_setting_value(config, ["namelist:run_stochastic","n1"])
        n2=self.get_setting_value(config, ["namelist:run_stochastic","n2"])
        self.remove_setting(config, ["namelist:run_stochastic", "n1"])
        self.remove_setting(config, ["namelist:run_stochastic", "n2"])
        self.add_setting(config, ["namelist:run_stochastic", "stph_n1"], n1)
        self.add_setting(config, ["namelist:run_stochastic", "stph_n2"], n2)
        # Replace smagorinksy bi-harmonic logical with scheme select via integer
        l_skeb2_biharm=self.get_setting_value(config,["namelist:run_stochastic",
                                                      "l_skeb2_biharm"])
        self.remove_setting(config, ["namelist:run_stochastic", 
                                     "l_skeb2_biharm"])
        if l_skeb2_biharm == ".true.":
            self.add_setting(config, ["namelist:run_stochastic",
                                      "skeb2_sdisp"], "7")
        else:
            self.add_setting(config, ["namelist:run_stochastic",
                                      "skeb2_sdisp"], "6")
        return config, self.reports


class vn101_t68(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #68 by Claudio Sanchez."""

    BEFORE_TAG = "vn10.1_t349"
    AFTER_TAG = "vn10.1_t68"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_track"])
        self.add_setting(config, ["namelist:run_track", "l_hoskins"], ".true.")
        self.add_setting(config, ["namelist:run_track", "nbot_850"], "6")
        self.add_setting(config, ["namelist:run_track", "ntop_850"], "42")
        self.add_setting(config, ["namelist:run_track", "ntop_tc"], "63")
        self.add_setting(config, ["namelist:run_track", "nlevs_avg"], "3")
        self.add_setting(config, ["namelist:run_track", "sm"], "0.1")

        model_type = self.get_setting_value(config,
                                      ["namelist:model_domain", "model_type"])
        if model_type == "5":
            # For SCM, get the source from the SHARED namelist
            source_file = "file:SHARED"
        else:
            # For all other apps, from the main NAMELIST
            source_file = "file:NAMELIST"

        namelist_source = self.get_setting_value(config, [source_file, "source"])
        if namelist_source:
            namelist_source = re.sub(r'namelist:run_dust',
                                     r'namelist:run_dust namelist:run_track',
                                     namelist_source)
            self.change_setting_value(config, [source_file, "source"],
                                      namelist_source)

        return config, self.reports 

class vn101_t402(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #402 by malcolmbrooks."""

    BEFORE_TAG = "vn10.1_t68"
    AFTER_TAG = "vn10.1_t402"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_cloud", "l_pc2_check_init"], ".false.")
        return config, self.reports


class vn101_t179(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #179 by Claudio Sanchez."""

    BEFORE_TAG = "vn10.1_t402"
    AFTER_TAG = "vn10.1_t179"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_stochastic", "l_x_eq_sin_x"], ".true.")

        return config, self.reports
        
class vn101_t343(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #343 by paulearnshaw."""

    BEFORE_TAG = "vn10.1_t179"
    AFTER_TAG = "vn10.1_t343"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.remove_setting(config, ["namelist:idealise", "instability_diagnostics"])
        return config, self.reports
        

class vn101_t195(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #195 by Anne McCabe."""

    BEFORE_TAG = "vn10.1_t343"
    AFTER_TAG = "vn10.1_t195"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_stochastic", "i_rp_scheme"], "0")
        self.add_setting(config, ["namelist:run_stochastic", "rp2_decorr_ts"], "216000.0")
        self.add_setting(config, ["namelist:run_stochastic", "l_rp2_cycle_in"], ".false.")
        self.add_setting(config, ["namelist:run_stochastic", "l_rp2_cycle_out"], ".false.")
        self.add_setting(config, ["namelist:run_stochastic", "rp2_cycle_tm"], "43200")
        self.add_setting(config, ["namelist:run_stochastic", "x1r_rp"], "0.22")
        self.add_setting(config, ["namelist:run_stochastic", "x1r_rp_min"], "0.07")
        self.add_setting(config, ["namelist:run_stochastic", "x1r_rp_max"], "0.52")
        self.add_setting(config, ["namelist:run_stochastic", "ndrop_surf_rp"], "7.5e+07")
        self.add_setting(config, ["namelist:run_stochastic", "ndrop_surf_rp_min"], "2.0e+07")
        self.add_setting(config, ["namelist:run_stochastic", "ndrop_surf_rp_max"], "10.0e+07")
        self.add_setting(config, ["namelist:run_stochastic", "ec_auto_rp"], "0.55")
        self.add_setting(config, ["namelist:run_stochastic", "ec_auto_rp_min"], "0.01")
        self.add_setting(config, ["namelist:run_stochastic", "ec_auto_rp_max"], "0.6")
        self.add_setting(config, ["namelist:run_stochastic", "lam_meta_rp"], "1.0")
        self.add_setting(config, ["namelist:run_stochastic", "lam_meta_rp_min"], "0.2")
        self.add_setting(config, ["namelist:run_stochastic", "lam_meta_rp_max"], "3.0")
        return config, self.reports

class vn101_t338(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #338 by Rachel Stratton."""

    BEFORE_TAG = "vn10.1_t195"
    AFTER_TAG = "vn10.1_t338"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # adding run_dyntest to shared so that recon can read this
        shared_source = self.get_setting_value(config, ["file:SHARED","source"])
        if shared_source:
            shared_source = re.sub(r'namelist:run_electric',
                                   r'namelist:run_electric namelist:run_dyntest',
                                   shared_source)
            self.change_setting_value(config, ["file:SHARED","source"], shared_source)
        # Add l_idealised_data to run_dyntest
        l_idealised_data_local = self.get_setting_value(
                config, ["namelist:idealise", "l_idealised_data"])
        if l_idealised_data_local == ".true.":
            self.add_setting(config, ["namelist:run_dyntest", "l_idealised_data"], ".true.")
        else:
            self.add_setting(config, ["namelist:run_dyntest", "l_idealised_data"], ".false.")

        # Remove l_idealised_data from idealise
        self.remove_setting(config, ["namelist:idealise", "l_idealised_data"])
        return config, self.reports


class vn101_t163(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #163 by Joe Mancell."""

    BEFORE_TAG = "vn10.1_t338"
    AFTER_TAG = "vn10.1_t163"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Coupling settings
        l_oasis_timers    = self.get_setting_value(config, ["namelist:nlstcatm",
                                                            "l_oasis_timers"])
        oasis_couple_freq = self.get_setting_value(config, ["namelist:nlstcatm",
                                                            "oasis_couple_freq"])
        l_couple_master   = self.get_setting_value(config, ["namelist:nlstcatm",
                                                            "l_couple_master"])
        l_oasis_icecalve  = self.get_setting_value(config, ["namelist:nlstcatm",
                                                            "l_oasis_icecalve"])
        # All four namelist variables are compulsory = true, so can safely add them
        # to the new namelist as they will already by present in all metadata
        # compliant apps
        self.add_setting(config, ["namelist:coupling_control", 
                                  "l_oasis_timers"], l_oasis_timers)
        self.add_setting(config, ["namelist:coupling_control", 
                                  "oasis_couple_freq"], oasis_couple_freq)
        self.add_setting(config, ["namelist:coupling_control", 
                                  "l_couple_master"], l_couple_master)
        self.add_setting(config, ["namelist:coupling_control", 
                                  "l_oasis_icecalve"], l_oasis_icecalve)
        # AC assim settings
        a_assim_start_min = self.get_setting_value(config, ["namelist:nlstcatm",
                                                            "a_assim_start_min"])
        a_assim_end_min = self.get_setting_value(config, ["namelist:nlstcatm",
                                                          "a_assim_end_min"])
        # Add to ACP namelist
        self.add_setting(config, ["namelist:acp", 
                                  "a_assim_start_min"], a_assim_start_min)
        self.add_setting(config, ["namelist:acp", 
                                  "a_assim_end_min"], a_assim_end_min)
        # 360 day calendar
        lcal360 = self.get_setting_value(config, ["namelist:nlstcatm", "lcal360"])
        # Add to nlstcall
        self.add_setting(config, ["namelist:nlstcall", "lcal360"], lcal360)
        
        # Carbon settings
        i_co2_opt = self.get_setting_value(config, ["namelist:nlstcatm",
                                                            "i_co2_opt"])
        l_co2_emits = self.get_setting_value(config, ["namelist:nlstcatm",
                                                          "l_co2_emits"])
        # Add to carbon options namelist
        self.add_setting(config, ["namelist:carbon_options", 
                                  "i_co2_opt"], i_co2_opt)
        self.add_setting(config, ["namelist:carbon_options", 
                                  "l_co2_emits"], l_co2_emits)
        # Add new option to general phys namelist
        self.add_setting(config, ["namelist:gen_phys_inputs", 
                                  "l_leads_temp_prog"], ".false.")
        # Add new namelists to the appropiate place in the text files containing the 
        # namelists
        shared_source = self.get_setting_value(config, ["file:SHARED","source"])
        if shared_source:
           shared_source = re.sub(r'namelist:model_domain',
                                  r'namelist:carbon_options namelist:'+
                                  'coupling_control namelist:model_domain',
                                  shared_source)
           shared_source = re.sub(r'namelist:run_eng_corr',
                                  r'namelist:run_eng_corr namelist:nlstcall ' +
                                  r'namelist:gen_phys_inputs ',
                                  shared_source)
           shared_source = re.sub(r'namelist:nlstcatm ','', shared_source)
           self.change_setting_value(config, ["file:SHARED","source"], shared_source)
        namelist_source = self.get_setting_value(config, ["file:NAMELIST","source"])
        if namelist_source:
            namelist_source = re.sub(r'namelist:model_domain',
                                     r'namelist:carbon_options namelist:'+
                                     'coupling_control namelist:model_domain',
                                     namelist_source)
            namelist_source = re.sub(r'namelist:nlstcatm ','', namelist_source)
            self.change_setting_value(config, ["file:NAMELIST","source"], namelist_source)
        # Now that nlstcall is in SHARED there is no need for the SCM to have a separate
        # CNTLALL file
        self.remove_setting(config, ["file:CNTLALL"])
        # remove whole of nlstcatm namelist
        self.remove_setting(config, ["namelist:nlstcatm"])
                                        
        return config, self.reports


class vn101_t575(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #575 by johnmedwards."""

    BEFORE_TAG = "vn10.1_t163"
    AFTER_TAG = "vn10.1_t575"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Include the option to correct the evolution of canopy temperature,
        # but in a disabled state.
        self.add_setting(config, ["namelist:temp_fixes", "l_dtcanfix"], ".false.")
        return config, self.reports

class vn101_t257(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #257 by James Manners, Nathan Mayne,
       and Chris Smith"""

    BEFORE_TAG = "vn10.1_t575"
    AFTER_TAG = "vn10.1_t257"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Set rotation rate in planet_constants namelist only
        l_rotating = self.get_setting_value(config, ["namelist:idealise", "l_rotating"])
	self.remove_setting(config, ["namelist:idealise", "angular_velocity"])
	self.remove_setting(config, ["namelist:idealise", "l_rotating"])
        if l_rotating == ".false.":
           self.add_setting(config, ["namelist:planet_constants","l_set_planet_rotation"], ".true.")
           self.change_setting_value(config, ["namelist:planet_constants", "omega"], "0.0")
        else:
           self.add_setting(config, ["namelist:planet_constants","l_set_planet_rotation"], ".false.")

        # RUN_Dyn namelist
        self.add_setting(config, ["namelist:run_dyn","eg_vert_damp_coeff"], "0.05")
        self.add_setting(config, ["namelist:run_dyn","l_viscosity"], ".false.")
        self.add_setting(config, ["namelist:run_dyn","horiz_viscosity"], "0.0")
        self.add_setting(config, ["namelist:run_dyn","vert_viscosity"], "0.0")

        # RUN_dyntest namelist
        self.add_setting(config, ["namelist:run_dyntest","l_hydrostatic_eg"], ".false.")

        # RUN_diffusion
        self.add_setting(config, ["namelist:run_diffusion","l_diff_wind_ew_only"], ".false.")

        # New items for Idealise namelist (from planet_suite_mod)
        self.add_setting(config, ["namelist:idealise", "tforce_number"], "0")
        self.add_setting(config, ["namelist:idealise", "trelax_number"], "0")
        self.add_setting(config, ["namelist:idealise", "fric_number"], "0")
        self.add_setting(config, ["namelist:idealise", "nsteps_consv_print"], "0")
        self.add_setting(config, ["namelist:idealise", "initial_tolerance"], "1.0e-11")
        self.add_setting(config, ["namelist:idealise", "profile_filename"], "")
        return config, self.reports

class vn101_t427(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #427 by Jess Standen."""

    BEFORE_TAG = "vn10.1_t257"
    AFTER_TAG = "vn10.1_t427"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        spiral_old_setting = self.get_setting_value(
                                  config, ["namelist:recon", "lspiral_s"])
        coast_value="1"
        if spiral_old_setting == ".true.":
            coast_value="2"
        self.add_setting(config, ["namelist:recon", "coast_adj_method"],
                             coast_value)
        self.remove_setting(config, ["namelist:recon", "lspiral_s"])
        return config, self.reports

class vn101_t605(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #605 by <author>."""

    BEFORE_TAG = "vn10.1_t427"
    AFTER_TAG = "vn10.1_t605"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:jules_pftparm", "ief_io"], "25.0,8.00,16.00,24.00,20.00")
        self.add_setting(config, ["namelist:jules_pftparm", "tef_io"], "1.2,2.4,0.8,1.2,0.8")
        self.add_setting(config, ["namelist:jules_pftparm", "mef_io"], "0.9,1.8,0.6,0.9,0.57")
        self.add_setting(config, ["namelist:jules_pftparm", "aef_io"], "0.43,0.87,0.29,0.43,0.20")
        self.add_setting(config, ["namelist:jules_pftparm", "ci_st_io"], "33.46,33.46,34.26,29.98,34.26")
        self.add_setting(config, ["namelist:jules_pftparm", "gpp_st_io"], "1.29E-07,2.58E-08,2.07E-07,3.42E-07,1.68E-007")

        return config, self.reports
        
class vn101_t294(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #294 by Stuart Webster."""

    BEFORE_TAG = "vn10.1_t605"
    AFTER_TAG = "vn10.1_t294"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_gwd", "nsigma"], "2.5")
        return config, self.reports


class vn101_t462(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #462 by Nic Gedney."""

    BEFORE_TAG = "vn10.1_t294"
    AFTER_TAG = "vn10.1_t462"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here

        self.add_setting(config, ["namelist:jules_hydrology", "l_wetland_unfrozen"], ".false.")
 
        return config, self.reports

class vn101_t634(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #634 by Ian Boutle."""

    BEFORE_TAG = "vn10.1_t462"
    AFTER_TAG = "vn10.1_t634"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_radiation","l_fsd_eff_res"],".false.")
        return config, self.reports

class vn101_t243(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #243 by johnmedwards."""

    BEFORE_TAG = "vn10.1_t634"
    AFTER_TAG = "vn10.1_t243"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Input your macro commands here
        npft = int(self.get_setting_value(config, ["namelist:jules_surface_types", "npft"]))

        self.add_setting(config, ["namelist:jules_radiation", "l_embedded_snow"], ".false.")
        self.add_setting(config, ["namelist:jules_radiation", "l_mask_snow_orog"], ".false.")
        self.add_setting(config, ["namelist:jules_radiation", "wght_alb"], "0.0,0.5,0.0,0.5")

        self.add_setting(config, ["namelist:jules_snow", "l_et_metamorph"], ".false.")
        self.add_setting(config, ["namelist:jules_snow", "l_snow_nocan_hc"], ".false.")
        self.add_setting(config, ["namelist:jules_snow", "l_snow_infilt"], ".false.")
        self.add_setting(config, ["namelist:jules_snow", "i_snow_cond_parm"], "0")
        self.add_setting(config, ["namelist:jules_snow", "a_snow_et"], "0.0")
        self.add_setting(config, ["namelist:jules_snow", "b_snow_et"], "0.0")
        self.add_setting(config, ["namelist:jules_snow", "c_snow_et"], "0.0")
        self.add_setting(config, ["namelist:jules_snow", "rho_snow_et_crit"], "0.0")
        self.add_setting(config, ["namelist:jules_snow", "unload_rate_cnst"], ','.join(['0.0']*npft))
        self.add_setting(config, ["namelist:jules_snow", "unload_rate_u"], ','.join(['0.0']*npft))
        self.add_setting(config, ["namelist:jules_snow", "can_clump"], ','.join(['0.0']*npft))
        self.add_setting(config, ["namelist:jules_snow", "lai_alb_lim_sn"], ','.join(['0.5']*npft))
        self.add_setting(config, ["namelist:jules_snow", "n_lai_exposed"], ','.join(['0.0']*npft))

        self.add_setting(config, ["namelist:jules_pftparm", "lai_alb_lim_io"], ','.join(['0.5']*npft))

        self.add_setting(config, ["namelist:jules_vegetation", "l_vegcan_soilfx"], ".false.")

        return config, self.reports


class vn101_t75(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #75 by Jamie Rae."""

    BEFORE_TAG = "vn10.1_t243"
    AFTER_TAG = "vn10.1_t75"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        self.add_setting(config, ["namelist:jules_sea_seaice", "ahmax"], "0.3")
        self.add_setting(config, ["namelist:jules_sea_seaice", "albicei_cice"], "0.36")
        self.add_setting(config, ["namelist:jules_sea_seaice", "albicev_cice"], "0.78")
        self.add_setting(config, ["namelist:jules_sea_seaice", "albpondi_cice"], "0.07")
        self.add_setting(config, ["namelist:jules_sea_seaice", "albpondv_cice"], "0.27")
        self.add_setting(config, ["namelist:jules_sea_seaice", "albsnowi_cice"], "0.70")
        self.add_setting(config, ["namelist:jules_sea_seaice", "albsnowv_cice"], "0.98")
        self.add_setting(config, ["namelist:jules_sea_seaice", "dalb_mlt_cice"], "-0.075")
        self.add_setting(config, ["namelist:jules_sea_seaice", "dalb_mlts_i_cice"], "-0.15")
        self.add_setting(config, ["namelist:jules_sea_seaice", "dalb_mlts_v_cice"], "-0.1")
        self.add_setting(config, ["namelist:jules_sea_seaice", "dt_bare_cice"], "1.0")
        self.add_setting(config, ["namelist:jules_sea_seaice", "dt_snow_cice"], "0.0")
        self.add_setting(config, ["namelist:jules_sea_seaice", "l_saldep_freeze"], ".false.")
        self.add_setting(config, ["namelist:jules_sea_seaice", "l_sice_meltponds_cice"], ".false.")
        self.add_setting(config, ["namelist:jules_sea_seaice", "l_sice_swpen"], ".false.")
        self.add_setting(config, ["namelist:jules_sea_seaice", "nice_use"], "1")
        self.add_setting(config, ["namelist:jules_sea_seaice", "pen_rad_frac_cice"], "0.2")
        self.add_setting(config, ["namelist:jules_sea_seaice", "snowpatch"], "0.02")
        self.add_setting(config, ["namelist:jules_sea_seaice", "sw_beta_cice"], "0.6")
        alphab = self.get_setting_value(config, ["namelist:run_radiation", "alphab"])
        self.remove_setting(config, ["namelist:run_radiation", "alphab"])
        self.add_setting(config, ["namelist:jules_sea_seaice", "alphab"], alphab)
        alphac = self.get_setting_value(config, ["namelist:run_radiation", "alphac"])
        self.remove_setting(config, ["namelist:run_radiation", "alphac"])
        self.add_setting(config, ["namelist:jules_sea_seaice", "alphac"], alphac)
        alpham = self.get_setting_value(config, ["namelist:run_radiation", "alpham"])
        self.remove_setting(config, ["namelist:run_radiation", "alpham"])
        self.add_setting(config, ["namelist:jules_sea_seaice", "alpham"], alpham)
        dalb_bare_wet = self.get_setting_value(config, ["namelist:run_radiation", "dalb_bare_wet"])
        self.remove_setting(config, ["namelist:run_radiation", "dalb_bare_wet"])
        self.add_setting(config, ["namelist:jules_sea_seaice", "dalb_bare_wet"], dalb_bare_wet)
        dt_bare = self.get_setting_value(config, ["namelist:run_radiation", "dt_bare"])
        self.remove_setting(config, ["namelist:run_radiation", "dt_bare"])
        self.add_setting(config, ["namelist:jules_sea_seaice", "dt_bare"], dt_bare)
        dtice = self.get_setting_value(config, ["namelist:run_radiation", "dtice"])
        self.remove_setting(config, ["namelist:run_radiation", "dtice"])
        self.add_setting(config, ["namelist:jules_sea_seaice", "dtice"], dtice)
        pen_rad_frac = self.get_setting_value(config, ["namelist:run_radiation", "pen_rad_frac"])
        self.remove_setting(config, ["namelist:run_radiation", "pen_rad_frac"])
        self.add_setting(config, ["namelist:jules_sea_seaice", "pen_rad_frac"], pen_rad_frac)
        ssalphac = self.get_setting_value(config, ["namelist:run_radiation", "ssalphac"])
        self.remove_setting(config, ["namelist:run_radiation", "ssalphac"])
        self.add_setting(config, ["namelist:jules_sea_seaice", "ssalphac"], ssalphac)
        ssalpham = self.get_setting_value(config, ["namelist:run_radiation", "ssalpham"])
        self.remove_setting(config, ["namelist:run_radiation", "ssalpham"])
        self.add_setting(config, ["namelist:jules_sea_seaice", "ssalpham"], ssalpham)
        ssdtice = self.get_setting_value(config, ["namelist:run_radiation", "ssdtice"])
        self.remove_setting(config, ["namelist:run_radiation", "ssdtice"])
        self.add_setting(config, ["namelist:jules_sea_seaice", "ssdtice"], ssdtice)
        sw_beta = self.get_setting_value(config, ["namelist:run_radiation", "sw_beta"])
        self.remove_setting(config, ["namelist:run_radiation", "sw_beta"])
        self.add_setting(config, ["namelist:jules_sea_seaice", "sw_beta"], sw_beta)

        return config, self.reports

class vn101_t637(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #637 by Ian Boutle."""

    BEFORE_TAG = "vn10.1_t75"
    AFTER_TAG = "vn10.1_t637"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        l_init_tile_t = self.get_setting_value(config, ["namelist:recon","l_init_tile_t_zerofrac"])
        if l_init_tile_t == ".true.":
            self.add_setting(config, ["namelist:recon","l_use_zero_frac_tile_temp"],".false.")
        else:
            self.add_setting(config, ["namelist:recon","l_use_zero_frac_tile_temp"],".true.")
        self.remove_setting(config, ["namelist:recon","l_init_tile_t_zerofrac"])
        return config, self.reports

class vn101_t570(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #570 by Jane Mulcahy."""

    BEFORE_TAG = "vn10.1_t637"
    AFTER_TAG = "vn10.1_t570"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_ukca", "mode_incld_so2_rfrac"], "0.25")
        return config, self.reports

class vn101_t578(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #578 by Maggie Hendry."""

    BEFORE_TAG = "vn10.1_t570"
    AFTER_TAG = "vn10.1_t578"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Compulsory l_snow_tile_gbm_ancil_fix set to false
        self.add_setting(config, ["namelist:recon", "l_snow_tile_gbm_ancil_fix"], ".false.")
        return config, self.reports


class vn101_t576(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #576 by stephaniewoodward."""

    BEFORE_TAG = "vn10.1_t578"
    AFTER_TAG = "vn10.1_t576"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_ukca", "l_ukca_dust" ],
                         str(".false."))
        self.add_setting(config, ["namelist:run_ukca", "l_ukca_primdu" ],
                         str(".false."))
        return config, self.reports

class vn101_t425(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #425 by Roddy Sharp."""

    BEFORE_TAG = "vn10.1_t576"
    AFTER_TAG = "vn10.1_t425"

    def upgrade(self, config, meta_config=None):
        """Upgrade macro for ticket #425 by Roddy Sharp.
           Remove wet_levels from nlsizes namelist.
           Also check that wet_levels = model_levels and warn if not."""
        # find values of model_levels and wet_levels
        model_levels = self.get_setting_value(config, 
                                           ["namelist:nlsizes", "model_levels"])
        wet_levels   = self.get_setting_value(config, 
                                             ["namelist:nlsizes", "wet_levels"])
        if wet_levels is not None:
            if model_levels != wet_levels:
                warn_msg = ('Wet_Levels does not equal Model_Levels' +
                ' but code is going to assume it does ')
                self.add_report("namelist:nlsizes", "wet_levels", "None",
                                 warn_msg, is_warning=True)
            self.remove_setting(config, ["namelist:nlsizes", "wet_levels"])
        
        return config, self.reports

class vn101_t617(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #617 by Adam Voysey."""

    BEFORE_TAG = "vn10.1_t425"
    AFTER_TAG = "vn10.1_t617"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add new compulsory settings
        self.add_setting(config, ["namelist:io_control", "io_filesystem_profile"], "0")
        self.add_setting(config, ["namelist:lustre_control", "default_stripe_count"], "0")
        self.add_setting(config, ["namelist:lustre_control", "default_stripe_offset"], "-1")
        self.add_setting(config, ["namelist:lustre_control", "default_stripe_pattern"], "0")
        self.add_setting(config, ["namelist:lustre_control", "default_stripe_size"], "0")
        self.add_setting(config, ["namelist:lustre_control", "number_custom_files"], "0")
        self.add_setting(config, ["namelist:lustre_control_custom_files", "striped_file"], "'dummy_file'")
        self.add_setting(config, ["namelist:lustre_control_custom_files", "stripe_offset"], "-1")
        self.add_setting(config, ["namelist:lustre_control_custom_files", "stripe_pattern"], "0")
        self.add_setting(config, ["namelist:lustre_control_custom_files", "stripe_size"], "0")
        self.add_setting(config, ["namelist:lustre_control_custom_files", "stripe_count"], "0")
        # If using IOSCNTL, add new namelists
        IOSCNTL_file_list = self.get_setting_value(config, ["file:IOSCNTL", "source"])
        if IOSCNTL_file_list:
            IOSCNTL_file_list += (" (namelist:lustre_control)")
            IOSCNTL_file_list += (" (namelist:lustre_control_custom_files)")
            self.change_setting_value(config, ["file:IOSCNTL", "source"], IOSCNTL_file_list)

        return config, self.reports

class vn101_t262(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #262 by Adrian Lock."""

    BEFORE_TAG = "vn10.1_t617"
    AFTER_TAG = "vn10.1_t262"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_bl","dec_thres_cu"],"0.05")
        self.add_setting(config, ["namelist:run_bl","l_reset_dec_thres"],".false.")
        self.add_setting(config, ["namelist:run_cloud","forced_cu_fac"],"0.5")

        return config, self.reports


class vn101_t1389(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1389 by Paul Cresswell."""

    BEFORE_TAG = "vn10.1_t262"
    AFTER_TAG = "vn10.1_t1389"

    def upgrade(self, config, meta_config=None):
        """Add essential ioscntl inputs to all UM apps."""
        self.add_setting(config, ["namelist:ioscntl", "ios_spacing"], "0")
        self.add_setting(config, ["namelist:ioscntl", "ios_offset"], "0")
        self.add_setting(config, ["namelist:ioscntl", "ios_tasks_per_server"],
                         "1")

        # Fix invalid values of ios_spacing:
        ios_spacing = self.get_setting_value(config,
                                 ["namelist:ioscntl", "ios_spacing"])
        if int(ios_spacing) < 0:
           self.change_setting_value(config,
                                 ["namelist:ioscntl", "ios_spacing"], "0")

        return config, self.reports


class vn102_t718(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.1_t1389"
    AFTER_TAG = "vn10.2"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "10.2")
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMSTR"])
        if stashmaster_path:
           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",
                                     "vn10.2/ctldata",
                                      stashmaster_path)
           self.change_setting_value(config, ["env", "STASHMSTR"],
                                     stashmaster_path)
        return config, self.reports


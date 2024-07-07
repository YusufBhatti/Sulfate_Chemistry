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


class vn105_t1883(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1883 by Mohamed Zerroukat."""

    BEFORE_TAG = "vn10.5"
    AFTER_TAG = "vn10.5_t1883"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, [ "namelist:run_dyn", "l_conservation_moist_zlf" ],
                         str(".false."))
        self.add_setting(config, [ "namelist:run_dyn", "zlf_conservation_moist_option" ],
                         "1")
        self.add_setting(config, [ "namelist:run_dyn", "l_mono_moist_pmf" ],
                         str(".false."))
        self.add_setting(config, [ "namelist:run_dyn", "l_mono_moist_pmf_stringent" ],
                         str(".false."))
        return config, self.reports


class vn105_t1411(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1411 by NickSavage."""

    BEFORE_TAG = "vn10.5_t1883"
    AFTER_TAG = "vn10.5_t1411"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # add new logical l_ukca_ddep_lev1 to the UKCA namelist
        # with a default of false (no change to results)
        self.add_setting(config, ["namelist:run_ukca", "l_ukca_ddep_lev1"],
                         ".false.")
        return config, self.reports


class vn105_t1642(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1642 by Joe Mancell."""

    BEFORE_TAG = "vn10.5_t1411"
    AFTER_TAG = "vn10.5_t1642"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:recon_technical", 
                                  "output_field_list"], "-32768")
        select_output_fields = int(self.get_setting_value(
                                                      config, 
                                                      ["namelist:recon_technical",
                                                      "select_output_fields"]))
        l_basic_interp = ".false."
        select_output_fields_new_val = "0"
        if (select_output_fields == 1):
            # interp_all_basic
            select_output_fields_new_val = "3"
            l_basic_interp = ".true."
        elif (select_output_fields == 2):
            # interp_all_rebalance
            select_output_fields_new_val = "3"
            l_basic_interp = ".false."
        self.add_setting(config, ["namelist:recon_technical", 
                                  "l_basic_interp"], l_basic_interp)
        self.change_setting_value(config, ["namelist:recon_technical",
                                           "select_output_fields"], 
                                           select_output_fields_new_val)
        return config, self.reports
        
class vn105_t1708(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1708 by Glenn Greed."""

    BEFORE_TAG = "vn10.5_t1642"
    AFTER_TAG = "vn10.5_t1708"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        endgame = self.get_setting_value(config,
                       ['namelist:run_dyn', 'l_endgame'])

        if endgame == ".false." or endgame == None :
            raise UpgradeError("New Dynamics apps are not supported from UM10.6. "  
                               "The only exception is if this is an l_endgame=.false. Recon task for use with VAR. "
                               "For these, set l_endgame=ND_RCF prior to running this upgrade macro." )

        elif endgame == "ND_RCF" : 

           self.add_setting(config,['namelist:horizont', 'output_grid_stagger'],
                         value = '3')
        else :
            
           self.add_setting(config,['namelist:horizont', 'output_grid_stagger'],
                         value = '6')
            
        self.remove_setting(config, ["namelist:run_dyn", "l_endgame"])
        self.remove_setting(config, ["namelist:idealise", "mod_layers"])
        self.remove_setting(config, ["namelist:idealise", "mag"])
        self.remove_setting(config, ["namelist:run_dyn", "alpha_1"])
        self.remove_setting(config, ["namelist:run_dyn", "alpha_2"])
        self.remove_setting(config, ["namelist:run_dyn", "alpha_3"])
        self.remove_setting(config, ["namelist:run_dyn", "alpha_4"])
        self.remove_setting(config, ["namelist:run_dyn", "alpha_1_2"])
        self.remove_setting(config, ["namelist:run_dyn", "alpha_2_2"])
        self.remove_setting(config, ["namelist:run_dyn", "alpha_3_2"])
        self.remove_setting(config, ["namelist:run_dyn", "alpha_4_2"])
        self.remove_setting(config, ["namelist:run_dyn", "numcycles"])
        self.remove_setting(config, ["namelist:run_dyn", "l_new_tdisc"])
        self.remove_setting(config, ["namelist:run_dyn", "l_thmono_fixed"])
        self.remove_setting(config, ["namelist:run_dyn", "i_nd_solver_vn"])
        self.remove_setting(config, ["namelist:run_dyn", "l_transparent"])
        self.remove_setting(config, ["namelist:run_dyn", "gcr_tol2"])
        self.remove_setting(config, ["namelist:run_dyn", "gcr_use_residual_tol"])
        self.remove_setting(config, ["namelist:run_dyn", "l_gcr_cycle_opt"])
        self.remove_setting(config, ["namelist:run_dyn", "gcr_adi_add_full_soln"])
        self.remove_setting(config, ["namelist:run_dyn", "gcr_adi_pseudo_timestep"])
        self.remove_setting(config, ["namelist:run_dyntest", "l_adjust_wet"])
        self.remove_setting(config, ["namelist:run_dyntest", "l_free_slip"])
        self.remove_setting(config, ["namelist:run_dyntest", "uv_limit"])
        self.remove_setting(config, ["namelist:run_dyntest", "trap_option"])
        self.remove_setting(config, ["namelist:run_dyntest", "l_trap_theta"])
        self.remove_setting(config, ["namelist:run_dyntest", "l_trap_w"])
        self.remove_setting(config, ["namelist:run_dyntest", "cw_max"])
        self.remove_setting(config, ["namelist:run_dyntest", "cw_test_lev"])
        self.remove_setting(config, ["namelist:run_dyntest", "max_thinc"])
        self.remove_setting(config, ["namelist:run_sl", "l_2d_sl_geometry"])
        self.remove_setting(config, ["namelist:run_sl", "l_sl_halo_reprod"])
        self.remove_setting(config, ["namelist:run_sl", "instability_diagnostics"])

        return config, self.reports        

class vn105_t2077(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2077 by Adrian Lock."""

    BEFORE_TAG = "vn10.5_t1708"
    AFTER_TAG = "vn10.5_t2077"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:temp_fixes",
                                  "l_fix_dyndiag"], ".false.")
        return config, self.reports

class vn105_t971(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #971 by markrichardson."""

    BEFORE_TAG = "vn10.5_t2077"
    AFTER_TAG = "vn10.5_t971"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_ukca","i_mode_ncpc"],"4")
        return config, self.reports

class vn105_t2031(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2031 by Claudio Sanchez."""

    BEFORE_TAG = "vn10.5_t971"
    AFTER_TAG = "vn10.5_t2031"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_free_tracers", "l_pv_tracer"], ".false.")
        self.add_setting(config, ["namelist:run_free_tracers", "l_pv_dyn"], ".false.")
        self.add_setting(config, ["namelist:run_free_tracers", "l_pv_split_rad"], ".false.")
        self.add_setting(config, ["namelist:run_free_tracers", "l_pv_adv_only"], ".false.")
        return config, self.reports

class vn105_t397(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #397 by Rachel Stratton."""

    BEFORE_TAG = "vn10.5_t2031"
    AFTER_TAG = "vn10.5_t397"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.remove_setting(config, ["namelist:temp_fixes", "l_glue_conv5a"])
        return config, self.reports

class vn105_t1466(rose.upgrade.MacroUpgrade):
    
    """Upgrade macro for ticket #1466 by Alan J Hewitt."""
    
    BEFORE_TAG = "vn10.5_t397"
    AFTER_TAG  = "vn10.5_t1466"
    
    def upgrade(self, config, meta_config=None):
        """Add netcdf_varname to namelist:items. Then hash sorting."""
        # Input your macro commands here
        values = []
        for obi in config.get_value():
            # is it in items namelist?
            if re.search(r'namelist:items', obi):
                stash_req = self.get_setting_value(config,[obi, "stash_req"])
                source = self.get_setting_value(config,[obi, "source"])
                if stash_req is not None:
                    # obtain number of comma separated elements in stash_req
                    lst = stash_req.split(",")
                    number = len(lst)
                    if len(lst) > 1:
                        if source == 10:
                            self.add_setting(config,
                                             [obi, "netcdf_varname"],
                                             ",".join(["'unset'"]*number))
                        else:
                            # State (trigger ignored) is required since
                            #  netcdf_varname is triggered when source = 10
                            self.add_setting(config,
                                             [obi, "netcdf_varname"],
                                             ",".join(["'unset'"]*number),
                                             state=config.STATE_SYST_IGNORED)
                    else:
                        if source == 10:
                            self.add_setting(config,
                                             [obi, "netcdf_varname"],
                                             "'unset'")
                        else:
                            # State (trigger ignored) is required since
                            #  netcdf_varname is triggered when source = 10
                            self.add_setting(config,
                                             [obi, "netcdf_varname"],
                                             "'unset'",
                                             state=config.STATE_SYST_IGNORED)
        
        import os
        sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "vn10.5", "lib", "python", "macros"))
        
        import stash_indices
        
        config, reports = stash_indices.TidyStashTransform().transform(config)
        self.reports += reports
        
        return config, self.reports
        

class vn105_t1350(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1350 by <mohitdalvi>."""

    BEFORE_TAG = "vn10.5_t1466"
    AFTER_TAG = "vn10.5_t1350"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Introduce logical l_ukca_so2ems_expvolc to the 
        # Run_ukca namelist with default 'False' value
 
        self.add_setting(config, ["namelist:run_ukca","l_ukca_so2ems_expvolc"],
                          value=".false.", forced=True)
        return config, self.reports

class vn105_t1467(rose.upgrade.MacroUpgrade):
    
    """Add namelist:run_glomap_aeroclim for ticket #1467 by Alan J Hewitt."""
    
    BEFORE_TAG = "vn10.5_t1350"
    AFTER_TAG = "vn10.5_t1467"
    
    def upgrade(self, config, meta_config=None):
        """Add namelist:run_glomap_aeroclim"""
        # Input your macro commands here
        
        # Add new section for namelist:run_glomap_aeroclim
        self.add_setting(config, ["namelist:run_glomap_aeroclim"], value=None,
                         forced=True, state=None, comments=None, info=None)
        
        # l_glomap_mode_clim should default to false
        self.add_setting(config, ["namelist:run_glomap_aeroclim",
                                  "l_glomap_mode_clim"],
                         value=".false.", forced=True)
        
        # i_glomap_clim_setup defaults to 2
        self.add_setting(config, ["namelist:run_glomap_aeroclim",
                                  "i_glomap_clim_setup"],
                         value="2", forced=True)
        
        # Add namelist model_domain to the list of namelists that are output in
        # the SHARED file
        current_value = self.get_setting_value(config,["file:SHARED", "source"])
        if current_value != None:
            new_value = current_value.replace(' namelist:run_ukca ',
                                              ' namelist:run_glomap_aeroclim namelist:run_ukca ')
            
            self.change_setting_value(config,["file:SHARED", "source"],
                                      value=new_value, forced=True)
        
        return config, self.reports

class vn105_t2315(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2315 by J. M. Edwards."""

    BEFORE_TAG = "vn10.5_t1467"
    AFTER_TAG = "vn10.5_t2315"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:jules_snow", "i_grain_growth_opt"], "0")
        self.add_setting(config, ["namelist:jules_snow", "i_relayer_opt"], "0")
        return config, self.reports

class vn105_t1990(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1990 by James Manners."""

    BEFORE_TAG = "vn10.5_t2315"
    AFTER_TAG = "vn10.5_t1990"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        l_mr = self.get_setting_value(config, ["namelist:gen_phys_inputs","l_mr_physics1"])
        if l_mr==".true.":
            self.add_setting(config,["namelist:run_radiation","l_hydrostatic_mass"], ".false.")
            self.add_setting(config,["namelist:run_radiation","l_moist_heat_capacity"], ".true.")
        else:
            self.add_setting(config,["namelist:run_radiation","l_hydrostatic_mass"], ".true.")
            self.add_setting(config,["namelist:run_radiation","l_moist_heat_capacity"], ".false.")
        l_extra_top = self.get_setting_value(config, ["namelist:r2lwclnl","l_extra_top_lw"])
        self.add_setting(config,["namelist:run_radiation","l_extra_top"], l_extra_top)
        self.remove_setting(config,["namelist:r2lwclnl","l_extra_top_lw"])
        self.remove_setting(config,["namelist:r2lwclnl","l_extra_top_lw2"])
        self.remove_setting(config,["namelist:r2swclnl","l_extra_top_sw"])
        self.remove_setting(config,["namelist:r2swclnl","l_extra_top_sw2"])
        return config, self.reports

class vn105_t1704(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1704 by James Manners."""

    BEFORE_TAG = "vn10.5_t1990"
    AFTER_TAG = "vn10.5_t1704"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config,["namelist:radfcdia","co2_mmr_scl"],"0.0")
        self.add_setting(config,["namelist:radfcdia","co2_mmr_add"],"4.348e-04")
        self.add_setting(config,["namelist:radfcdia","n2o_mmr_scl"],"0.0")
        self.add_setting(config,["namelist:radfcdia","n2o_mmr_add"],"4.205e-07")
        self.add_setting(config,["namelist:radfcdia","ch4_mmr_scl"],"0.0")
        self.add_setting(config,["namelist:radfcdia","ch4_mmr_add"],"4.461e-07")
        self.add_setting(config,["namelist:radfcdia","o2_mmr_scl"],"0.0")
        self.add_setting(config,["namelist:radfcdia","o2_mmr_add"],"0.2314")
        self.add_setting(config,["namelist:radfcdia","cfc11_mmr_scl"],"0.0")
        self.add_setting(config,["namelist:radfcdia","cfc11_mmr_add"],"0.0")
        self.add_setting(config,["namelist:radfcdia","cfc12_mmr_scl"],"0.0")
        self.add_setting(config,["namelist:radfcdia","cfc12_mmr_add"],"0.0")
        self.add_setting(config,["namelist:radfcdia","cfc113_mmr_scl"],"0.0")
        self.add_setting(config,["namelist:radfcdia","cfc113_mmr_add"],"0.0")
        self.add_setting(config,["namelist:radfcdia","hcfc22_mmr_scl"],"0.0")
        self.add_setting(config,["namelist:radfcdia","hcfc22_mmr_add"],"0.0")
        self.add_setting(config,["namelist:radfcdia","hfc125_mmr_scl"],"0.0")
        self.add_setting(config,["namelist:radfcdia","hfc125_mmr_add"],"0.0")
        self.add_setting(config,["namelist:radfcdia","hfc134a_mmr_scl"],"0.0")
        self.add_setting(config,["namelist:radfcdia","hfc134a_mmr_add"],"0.0")
        self.add_setting(config,["namelist:radfcdia","c2c_o2"],".false.")
        return config, self.reports

class vn105_t2371(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2371 by Glenn Greed."""

    BEFORE_TAG = "vn10.5_t1704"
    AFTER_TAG = "vn10.5_t2371"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        valid_high_order_scheme = ['0','1','2','5','7','8','9']
        high_order_scheme_string = self.get_setting_value(config, ["namelist:run_sl","high_order_scheme"])
        high_order_scheme = high_order_scheme_string.split(',')
        for choice in high_order_scheme:
            if choice not in valid_high_order_scheme:
                raise UpgradeError (
                    'Unable to upgrade app as this uses an invalid or deprecated high_order_scheme\n' +
                    'Please select a valid high_order_scheme for ENDGame configurations')

        valid_monotone_scheme = ['0','1']
        monotone_scheme_string = self.get_setting_value(config, ["namelist:run_sl","monotone_scheme"])
        monotone_scheme = monotone_scheme_string.split(',')
        for choice in monotone_scheme:
            if choice not in valid_monotone_scheme:
                raise UpgradeError (
                    'Unable to upgrade app as this uses an invalid or deprecated monotone_scheme\n' +
                    'Please select a valid monotone_scheme for ENDGame configurations')
        return config, self.reports

class vn105_t667(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #667 by Maggie Hendry."""

    BEFORE_TAG = "vn10.5_t2371"
    AFTER_TAG = "vn10.5_t667"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Rename tile_map to tile_map_ids
        tile_map = self.get_setting_value(config, ["namelist:jules_surface_types", "tile_map"], no_ignore=True)
        if tile_map is not None:
            self.add_setting(config, ["namelist:jules_surface_types", "tile_map_ids"], tile_map)
        self.remove_setting(config, ["namelist:jules_surface_types", "tile_map"])
        return config, self.reports


class vn105_t2105(rose.upgrade.MacroUpgrade):
    
    """Upgrade macro for ticket #2105 by AdrianLock."""

    BEFORE_TAG = "vn10.5_t667"
    AFTER_TAG = "vn10.5_t2105"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config,["namelist:run_bl","z_nl_bl_levels"],"6000.0")
        self.remove_setting(config,["namelist:run_bl","nl_bl_levels"])
        return config, self.reports

class vn105_t2306(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2306 by Joe Mancell."""

    BEFORE_TAG = "vn10.5_t2105"
    AFTER_TAG = "vn10.5_t2306"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        lamipii = self.get_setting_value(config, ["namelist:ancilcta", 
                                                  "lamipii"])
        # Set both new namelist variables to the value of lamipii
        self.add_setting(config, ["namelist:ancilcta", 
                                  "l_amipii_ice_processing"], lamipii)
        self.add_setting(config, ["namelist:ancilcta", 
                                  "use_lookup_dates_anc_time_interp"], lamipii)

        self.remove_setting(config,["namelist:ancilcta","lamipii"])

        return config, self.reports


class vn105_t2174(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2174 by Steven Rumbold."""

    BEFORE_TAG = "vn10.5_t2306"
    AFTER_TAG = "vn10.5_t2174"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_ukca", "l_ukca_scale_seadms_ems"], ".false.")
        self.add_setting(config,["namelist:run_ukca","seadms_ems_scaling"],"1.0")
        return config, self.reports

class vn105_t2195(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2195 by Matt Woodhouse."""

    BEFORE_TAG = "vn10.5_t2174"
    AFTER_TAG = "vn10.5_t2195"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_ukca", "l_ukca_prim_moc"], ".false.")
        return config, self.reports


class vn105_t2359(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2359 by AdrianLock."""

    BEFORE_TAG = "vn10.5_t2195"
    AFTER_TAG = "vn10.5_t2359"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        lev_pert_theta = self.get_setting_value(config, ["namelist:run_stochastic","lev_pert_theta"])
        self.remove_setting(config,["namelist:run_stochastic","lev_pert_theta"])
        if lev_pert_theta == "20":
            self.add_setting(config,
                ["namelist:run_stochastic", "z_pert_theta"], "1500.0")
        else:
            self.add_setting(config,
                ["namelist:run_stochastic", "z_pert_theta"], "400.0")
        return config, self.reports


class vn105_t1130(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1130 by Stephen Haddad."""
    
    BEFORE_TAG = "vn10.5_t2359"
    AFTER_TAG = "vn10.5_t1130"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.remove_setting(config,['namelist:nlcfiles','rpseed'])
        self.remove_setting(config,['namelist:run_stochastic','l_stphseed_file'])
        cycle_in_active = \
            self.get_setting_value(config,
                                   ['namelist:run_stochastic',
                                   'l_rp2_cycle_in'],
                                   no_ignore=True)
        cycle_out_active = \
            self.get_setting_value(config,
                                   ['namelist:run_stochastic',
                                   'l_rp2_cycle_out'],
                                   no_ignore=True)
        rp2_seed_value = "'unset'"
        if cycle_in_active is not None and cycle_out_active is not None:
           if cycle_in_active == '.true.' and cycle_out_active == '.true.':
               rp2_seed_value="'$DATAM/seedfiles/$RUNID.seed'"
        self.add_setting(config,
                         ['namelist:nlcfiles','rp2_seed'],
                         value=rp2_seed_value)
        return config, self.reports


class vn105_t2324(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2324 by John Hemmings."""

    BEFORE_TAG = "vn10.5_t1130" 
    AFTER_TAG = "vn10.5_t2324"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        # This is an empty macro to force apps to be upgraded to the new
        # metadata version
        return config, self.reports


class vn105_t1958(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1958 by <Mohit Dalvi>."""

    BEFORE_TAG = "vn10.5_t2324"
    AFTER_TAG = "vn10.5_t1958"

    def upgrade(self, config, meta_config=None):
        """
           Add a logical to the temporary logicals namelist,
           with default value = FALSE
        """

        # Input your macro commands here
        self.add_setting(config, 
                         ["namelist:temp_fixes", "l_fix_ukca_deriv_init"],
                         ".false.")
        return config, self.reports

class vn106_t2411(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.5_t1958"
    AFTER_TAG = "vn10.6"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "10.6", forced=True)
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMSTR"])
        if stashmaster_path:
           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",
                                     "vn10.6/ctldata",
                                      stashmaster_path)
           self.change_setting_value(config, ["env", "STASHMSTR"],
                                     stashmaster_path, forced=True)
        return config, self.reports


import rose.upgrade
import re
import sys
import os

class UpgradeError(Exception):

    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


class vn106_t1232(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1232 by <Mohit Dalvi>."""

    BEFORE_TAG = "vn10.6"
    AFTER_TAG = "vn10.6_t1232"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
       
        """ 
         Add new items in run_ukca namelist to control the
          resetting of near-surface values of age-air tracer
         - i_ageair_reset_method: default=1 (by level current method)
         - max_ageair_reset_level: default=10 (currently hardwired in code)
         - max_ageair_reset_height: default=2000.0 (inactive if 
                 i_ageair_reset_method=1)                  
        """

        self.add_setting(config, ["namelist:run_ukca", 
                                 "i_ageair_reset_method"], "1")
        self.add_setting(config, ["namelist:run_ukca", 
                                 "max_ageair_reset_level"], "10")
        self.add_setting(config, ["namelist:run_ukca", 
                                 "max_ageair_reset_height"], "2000.0")

        return config, self.reports
        
class vn106_t2321(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2321 by Alex Archibald."""

    BEFORE_TAG = "vn10.6_t1232"
    AFTER_TAG = "vn10.6_t2321"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # test to see if background aerosol is turned on 
        is_background = self.get_setting_value(config, ["namelist:run_ukca", "l_ukca_use_background_aerosol"])
        if is_background:
            self.add_setting(config, ["namelist:run_ukca", "i_ukca_sad_months"], "0")
            self.add_setting(config, ["namelist:run_ukca", "i_ukca_sad_start_year"], "0")
        else:
        # the only other option is to use the SPARC climatology Sulfate_SAD_SPARC_1950-2100.asc
            self.add_setting(config, ["namelist:run_ukca", "i_ukca_sad_months"], "1812", forced=True)
            self.add_setting(config, ["namelist:run_ukca", "i_ukca_sad_start_year"], "1950", forced=True)
        return config, self.reports


class vn106_t2503(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2503 by Alejandro Bodas-Salcedo."""

    BEFORE_TAG = "vn10.6_t2321"
    AFTER_TAG = "vn10.6_t2503"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add a new real cosp_sr_cloud with default value of 3.0.
        self.add_setting(config, [ "namelist:run_cosp", "cosp_sr_cloud" ],
                         value = "3.0")
        return config, self.reports


class vn106_t2535(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2535 by Adrian Lock."""

    BEFORE_TAG = "vn10.6_t2503"
    AFTER_TAG = "vn10.6_t2535"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:jules_snow", "graupel_options"],
                         "0")
        return config, self.reports


class vn106_t1964(rose.upgrade.MacroUpgrade):
    
    """Upgrade macro for ticket #1964 by Andy Wiltshire/Alan Hewitt"""
    
    BEFORE_TAG = "vn10.6_t2535"
    AFTER_TAG = "vn10.6_t1964"
    
    # Parameters for section and item numbers
    DIAG_SECT = 3
    NPP       = 262
    PR        = 263
    NPPP      = 291
    PRP       = 292
    
    def upgrade(self, config, meta_config=None):
        for obj in config.get_value():
            # is it an streq namelist?
            if re.search(r'namelist:streq', obj):
                # get the section and item numbers as integers
                item = self.get_setting_value(config,[obj,'item'])
                section = self.get_setting_value(config,[obj,'isec'])
                # check if they are the items we need to remove
                if (int(section) == self.DIAG_SECT):
                    if (int(item) == self.NPP):
                        self.change_setting_value(config,[obj,'item'],'662')
                    if (int(item) == self.PR):
                        self.change_setting_value(config,[obj,'item'],'663')
                    if (int(item) == self.NPPP):
                        self.change_setting_value(config,[obj,'item'],'691')
                    if (int(item) == self.PRP):
                        self.change_setting_value(config,[obj,'item'],'692')
        
        # Hash sorting:
        # Load the path of the hashing macro.
        sys.path.insert(1,
                        (os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "vn10.6", "lib", "python", "macros")))
        import stash_indices
        # Remove stash_indices from sys.path to avoid namespace conflicts.
        sys.path.pop(1)
        # Finally we rehash stash indices.
        config, reports = stash_indices.TidyStashTransform().transform(config)
        self.reports += reports

        return config, self.reports


class vn106_t2572(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2572 by Stuart Whitehouse."""

    BEFORE_TAG = "vn10.6_t1964"
    AFTER_TAG = "vn10.6.1"

    def upgrade(self, config, meta_config=None):
        """Increment version number."""
        return config, self.reports

class vn106_t2271(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2271 by Andy Malcolm."""

    BEFORE_TAG = "vn10.6.1"
    AFTER_TAG = "vn10.6_t2271"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here

        self.rename_setting(config, ["namelist:run_ukca", "i_mode_ncpc"], ["namelist:tuning_segments", "ukca_mode_seg_size" ])

        return config, self.reports

class vn106_t2117(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2117 by Michele Guidolin."""

    BEFORE_TAG = "vn10.6_t2271"
    AFTER_TAG = "vn10.6_t2117"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, [ "namelist:nlst_mpp", "nb_swap_bounds_method" ], "2")
        return config, self.reports


class vn106_t2579(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2579 by Steve Wardle"""

    BEFORE_TAG = "vn10.6_t2117"
    AFTER_TAG = "vn10.6_t2579"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        stashm = self.get_setting_value(config, ["env","STASHMSTR"])
        if stashm is not None:
            self.add_setting(
                config, ["env","STASHMASTER"], stashm, forced=True)
        self.remove_setting(config, ["env", "STASHMSTR"])
        return config, self.reports

class vn106_t2592(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2592 by Mohamed Zerroukat."""

    BEFORE_TAG = "vn10.6_t2579"
    AFTER_TAG = "vn10.6_t2592"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.remove_setting(config, ["namelist:run_dyn", "l_mono_moist_pmf"])
        self.remove_setting(config, ["namelist:run_dyn", "l_mono_moist_pmf_stringent"])
        return config, self.reports


class vn106_t2594(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2594 by Chris Smith."""

    BEFORE_TAG = "vn10.6_t2592"
    AFTER_TAG = "vn10.6_t2594"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Delete variables from IDEALISE namelist that belonged to
        # New Dynamics only.
        self.remove_setting(config, ["namelist:idealise", "cool_rate"])
        self.remove_setting(config, ["namelist:idealise", "dmptim"])
        self.remove_setting(config, ["namelist:idealise", "grow_steps"])
        self.remove_setting(config, ["namelist:idealise", "hdmp"])
        self.remove_setting(config, ["namelist:idealise", "height_u_in"])
        self.remove_setting(config, ["namelist:idealise", "idl_interp_option"])
        self.remove_setting(config, ["namelist:idealise", "l_bomex"])
        self.remove_setting(config, ["namelist:idealise", "l_code_test"])
        self.remove_setting(config, ["namelist:idealise", "l_cyclone"])
        self.remove_setting(config, ["namelist:idealise", "l_damp"])
        self.remove_setting(config, ["namelist:idealise", "l_fix_orog_hgt_lbc"])
        self.remove_setting(config, ["namelist:idealise", "l_perturb_correlate_time"])
        self.remove_setting(config, ["namelist:idealise", "l_perturb_correlate_tq"])
        self.remove_setting(config, ["namelist:idealise", "l_perturb_correlate_vert"])
        self.remove_setting(config, ["namelist:idealise", "l_pforce"])
        self.remove_setting(config, ["namelist:idealise", "l_polar_wind_zero"])
        self.remove_setting(config, ["namelist:idealise", "l_pressure_balance"])
        self.remove_setting(config, ["namelist:idealise", "l_rotate_winds"])
        self.remove_setting(config, ["namelist:idealise", "l_sh_williamson"])
        self.remove_setting(config, ["namelist:idealise", "l_wind_balance"])
        self.remove_setting(config, ["namelist:idealise", "num_profile_data"])
        self.remove_setting(config, ["namelist:idealise", "num_qforce_levels"])
        self.remove_setting(config, ["namelist:idealise", "num_tforce_levels"])
        self.remove_setting(config, ["namelist:idealise", "num_uvforce_levels"])
        self.remove_setting(config, ["namelist:idealise", "num_uvprofile_data"])
        self.remove_setting(config, ["namelist:idealise", "orog_hgt_lbc"])
        self.remove_setting(config, ["namelist:idealise", "pforce_option"])
        self.remove_setting(config, ["namelist:idealise", "plat_size_x"])
        self.remove_setting(config, ["namelist:idealise", "plat_size_y"])
        self.remove_setting(config, ["namelist:idealise", "qforce_option"])
        self.remove_setting(config, ["namelist:idealise", "qprofile_data"])
        self.remove_setting(config, ["namelist:idealise", "r_plane"])
        self.remove_setting(config, ["namelist:idealise", "suhe_relax"])
        self.remove_setting(config, ["namelist:idealise", "t_horizfn_data"])
        self.remove_setting(config, ["namelist:idealise", "t_horizfn_number"])
        self.remove_setting(config, ["namelist:idealise", "tforce_option"])
        self.remove_setting(config, ["namelist:idealise", "tprofile_data"])
        self.remove_setting(config, ["namelist:idealise", "tropics_deg"])
        self.remove_setting(config, ["namelist:idealise", "u_ramp_end"])
        self.remove_setting(config, ["namelist:idealise", "u_ramp_start"])
        self.remove_setting(config, ["namelist:idealise", "ujet_lat"])
        self.remove_setting(config, ["namelist:idealise", "ujet_width"])
        self.remove_setting(config, ["namelist:idealise", "uprofile_data"])
        self.remove_setting(config, ["namelist:idealise", "uv_horizfn_number"])
        self.remove_setting(config, ["namelist:idealise", "uvforce_option"])
        self.remove_setting(config, ["namelist:idealise", "vprofile_data"])
        self.remove_setting(config, ["namelist:idealise", "witch_power"])
        self.remove_setting(config, ["namelist:idealise", "z_qforce_data"])
        self.remove_setting(config, ["namelist:idealise", "z_tforce_data"])
        self.remove_setting(config, ["namelist:idealise", "z_uvforce_data"])
        self.remove_setting(config, ["namelist:idealise", "z_uvprofile_data"])
        self.remove_setting(config, ["namelist:idealise", "zdmp"])
        self.remove_setting(config, ["namelist:idealise", "zprofile_data"])
        self.remove_setting(config, ["namelist:idealise", "zprofile_orog"])
        return config, self.reports


class vn106_t2449(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2449 by AdrianLock."""

    BEFORE_TAG = "vn10.6_t2594"
    AFTER_TAG = "vn10.6_t2449"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, 
                ["namelist:run_stochastic", "zmin_pert_theta"], "0.0")
        z_pert_theta = self.get_setting_value(config, ["namelist:run_stochastic","z_pert_theta"])
        self.remove_setting(config,["namelist:run_stochastic","z_pert_theta"])
        self.add_setting(config, 
                ["namelist:run_stochastic", "zmax_pert_theta"], z_pert_theta)
        return config, self.reports


class vn106_t2692(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2692 by Mohamed Zerroukat."""

    BEFORE_TAG = "vn10.6_t2449"
    AFTER_TAG = "vn10.6_t2692"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config,["namelist:run_dyn", "zlf_maximum_height"],"30000.0")
        return config, self.reports


class vn106_t2711(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2711 by Glenn Greed."""

    BEFORE_TAG = "vn10.6_t2692"
    AFTER_TAG = "vn10.6_t2711"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.remove_setting(config, ["namelist:run_diffusion","l_adjust_theta"])
        self.remove_setting(config, ["namelist:run_diffusion","adjust_theta_end"])
        self.remove_setting(config, ["namelist:run_diffusion","adjust_theta_start"])
        return config, self.reports
        
class vn106_t2380(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2380 by Mohit Dalvi."""

    BEFORE_TAG = "vn10.6_t2711"
    AFTER_TAG = "vn10.6_t2380"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Remove temporary logical l_fix_ukca_deriv_init
        # as 'True' condition is now accepted in the code

        self.remove_setting(config,
                         ["namelist:temp_fixes","l_fix_ukca_deriv_init"])

        return config, self.reports       

class vn106_t2591(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2591 by Ian Boutle."""

    BEFORE_TAG = "vn10.6_t2380"
    AFTER_TAG = "vn10.6_t2591"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config,["namelist:run_diffusion","hyd_mix_opt"],"0")
        return config, self.reports

class vn106_t2556(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2556 by Jonathan Wilkinson."""

    BEFORE_TAG = "vn10.6_t2591"
    AFTER_TAG = "vn10.6_t2556"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, 
                         ["namelist:temp_fixes", "l_fix_riming"],
                         ".false.")
        return config, self.reports


class vn106_t2139(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2139 by Glenn Greed."""

    BEFORE_TAG = "vn10.6_t2556"
    AFTER_TAG = "vn10.6_t2139"

    def _rename_section_fast(self, config, section, new_section):
        """Rename a setting in the configuration."""
        config.value[new_section] = config.value[section]
        config.unset([section])
        info = self.INFO_RENAMED_SECT.format(section, new_section)
        self.add_report(new_section, None, None, info)
        self.add_report(section, None, None, self.INFO_REMOVED)

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        domain_list = []
        time_list = []
        use_list = [] 
        streq_list = []
        other_list = []

        namelsts = self.get_setting_value(config, ["file:STASHC","source"],)
        if (namelsts):
            namelsts=namelsts.replace("namelist:streq", "namelist:umstash_streq")
            namelsts=namelsts.replace("namelist:domain", "namelist:umstash_domain")
            namelsts=namelsts.replace("namelist:time", "namelist:umstash_time")
            namelsts=namelsts.replace("namelist:use", "namelist:umstash_use")
            self.change_setting_value(config, ["file:STASHC","source"], 
                                                             namelsts)

        namelsts = self.get_setting_value(config, ["file:RECONA","source"],)
        if (namelsts):
            namelsts=namelsts.replace("namelist:vertical", "namelist:recon_vertical")
            self.change_setting_value(config, ["file:RECONA","source"], 
                                                             namelsts)

        namelsts = self.get_setting_value(config, ["file:SHARED","source"],)
        if (namelsts):
            namelsts=namelsts.replace("namelist:urban_switches", "namelist:jules_urban_switches")
            self.change_setting_value(config, ["file:SHARED","source"], 
                                                             namelsts)

        namelsts = self.get_setting_value(config, ["file:ATMOSCNTL","source"],)
        if (namelsts):
            namelsts=namelsts.replace("namelist:urban2t_param", "namelist:jules_urban2t_param")
            self.change_setting_value(config, ["file:ATMOSCNTL","source"], 
                                                             namelsts)

        # for SCM
        namelsts = self.get_setting_value(config, ["file:CNTLATM","source"],)
        if (namelsts):
            namelsts=namelsts.replace("namelist:urban2t_param", "namelist:jules_urban2t_param")
            self.change_setting_value(config, ["file:CNTLATM","source"], 
                                                             namelsts)

        # renaming of objects that may have closely linked names 
        # for safety needs to be done outside of the original loop through
        # the objects otherwise the list gets messed up.

        for obj in config.get_value():
           if re.search(r'namelist:use\(',obj) or obj == "namelist:use":
               use_list.append(obj)
           elif re.search(r'namelist:time\(',obj) or obj == "namelist:time":
               time_list.append(obj)
           elif re.search(r'namelist:streq\(',obj) or obj == "namelist:streq":
               streq_list.append(obj)
           elif re.search(r'namelist:domain\(',obj) or obj == "namelist:domain":
               domain_list.append(obj)
           elif re.search(r'namelist:vertical',obj):
               other_list.append(obj)
           elif re.search(r'namelist:urban2t_param',obj):
               other_list.append(obj)
           elif re.search(r'namelist:urban_switches',obj):
               other_list.append(obj)
 
        for obj in domain_list:
            new_obj = obj.replace("namelist:domain", "namelist:umstash_domain")
            self._rename_section_fast(config, obj,new_obj)

        for obj in use_list:
            new_obj = obj.replace("namelist:use", "namelist:umstash_use")
            self._rename_section_fast(config, obj,new_obj)

        for obj in streq_list:
            new_obj = obj.replace("namelist:streq", "namelist:umstash_streq")
            self._rename_section_fast(config, obj,new_obj)

        for obj in time_list:
            new_obj = obj.replace("namelist:time", "namelist:umstash_time")
            self._rename_section_fast(config, obj,new_obj)

        for obj in other_list:
            if re.search(r'namelist:vertical',obj):
                new_obj = obj.replace("namelist:vertical", "namelist:recon_vertical")
                self._rename_section_fast(config, obj,new_obj)
            elif re.search(r'namelist:urban2t_param',obj):
                new_obj = obj.replace("namelist:urban2t_param", "namelist:jules_urban2t_param")
                self._rename_section_fast(config, obj,new_obj)
            elif re.search(r'namelist:urban_switches',obj):
                new_obj = obj.replace("namelist:urban_switches", "namelist:jules_urban_switches")
                self._rename_section_fast(config, obj,new_obj)

        return config, self.reports

class vn106_t2584(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2584 by Douglas Clark."""

    BEFORE_TAG = "vn10.6_t2139"
    AFTER_TAG = "vn10.6_t2584"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        
        """ 
         Create new namelist jules_soil_biogeochem for soil biogeochemical options
         and populate it with new items and items moved from other namelists.
        """

        # Ensure that the setting exists for the new namelist
        self.add_setting(config, ["namelist:jules_soil_biogeochem"])

        # Check if l_triffid is selected.
        l_trif = self.get_setting_value(config, ["namelist:jules_vegetation", "l_triffid"])
        # Add the new compulsory namelist item.
        if l_trif == ".true.":
            self.add_setting(config, ["namelist:jules_soil_biogeochem", "soil_bgc_model"], "2")
        else:
            self.add_setting(config, ["namelist:jules_soil_biogeochem", "soil_bgc_model"], "1")

        # Move settings from jules_surface.
        self.rename_setting(config, ["namelist:jules_surface", "kaps"], ["namelist:jules_soil_biogeochem", "kaps"] )
        self.rename_setting(config, ["namelist:jules_surface", "kaps_roth"], ["namelist:jules_soil_biogeochem", "kaps_roth"] )
        self.rename_setting(config, ["namelist:jules_surface", "q10_soil"], ["namelist:jules_soil_biogeochem", "q10_soil"] )
        self.rename_setting(config, ["namelist:jules_surface", "bio_hum_cn"], ["namelist:jules_soil_biogeochem", "bio_hum_cn"] )
        self.rename_setting(config, ["namelist:jules_surface", "sorp"], ["namelist:jules_soil_biogeochem", "sorp"] )
        self.rename_setting(config, ["namelist:jules_surface", "n_inorg_turnover"], ["namelist:jules_soil_biogeochem", "n_inorg_turnover"] )
        self.rename_setting(config, ["namelist:jules_surface", "diff_n_pft"], ["namelist:jules_soil_biogeochem", "diff_n_pft"] )
        self.rename_setting(config, ["namelist:jules_surface", "l_layeredc"], ["namelist:jules_soil_biogeochem", "l_layeredc"] )
        self.rename_setting(config, ["namelist:jules_surface", "tau_resp"], ["namelist:jules_soil_biogeochem", "tau_resp"] )
        self.rename_setting(config, ["namelist:jules_surface", "tau_lit"], ["namelist:jules_soil_biogeochem", "tau_lit"] )

        # Move settings from jules_vegetation.
        self.rename_setting(config, ["namelist:jules_vegetation", "l_q10"], ["namelist:jules_soil_biogeochem", "l_q10"] )
        self.rename_setting(config, ["namelist:jules_vegetation", "l_soil_resp_lev2"], ["namelist:jules_soil_biogeochem", "l_soil_resp_lev2"] )
                
        # update the source line of the UM SHARED file to add jules_soil_biogeochem
        shared_source = self.get_setting_value(config, ["file:SHARED","source"])
        if shared_source:
            shared_source = re.sub(r'namelist:jules_vegetation',
                                   r'namelist:jules_vegetation namelist:jules_soil_biogeochem',
                                   shared_source)
            self.change_setting_value(config, ["file:SHARED","source"], shared_source)
                
        return config, self.reports

class vn106_t2526(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2526 by Rachel Stratton."""

    BEFORE_TAG = "vn10.6_t2584"
    AFTER_TAG = "vn10.6_t2526"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Removing 0A convection - replacing with 6A diagnosis with 
	# l_param_conv = .false.
        i_convection_local = self.get_setting_value(
                config, ["namelist:run_convection", "i_convection_vn"])
        if i_convection_local == "0":
	    # Changing to convection version 6A
            self.change_setting_value(config, ["namelist:run_convection", "i_convection_vn"],"6")
	    # Ensure calling of convection scheme is off 
	    # !! infront of l_param_conv should be removed by triggering
            self.change_setting_value(config, ["namelist:run_convection", "l_param_conv"],".false.")
        return config, self.reports


class vn106_t1233(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1233 by Sam Cusworth."""

    BEFORE_TAG = "vn10.6_t2526"
    AFTER_TAG = "vn10.6_t1233"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Loop through entries in config
        for obj in config.get_value():
            # We want to add items to domain namelists but we also need to
            # rehash the stash indices. The full rose triggering is done after
            # the macro is applied. If we naively add these new variables to
            # a namelist with a hashed name and rehash the name before the full
            # triggering then the upgraded apps will fail to validate.
            # To get around this we carefully add the new items with the state
            # set appropriately and do not rely on the triggering to set it.
            if re.search(r'namelist:umstash_domain\(', obj):
                timeseries = self.get_setting_value(config, [obj, 'ts'])
                if timeseries == '.true.':
                    l_spml = self.get_setting_value(config, [obj, 'l_spml_ts'])
                    if l_spml == '.true.':
                        # this is the extraordinary case- someone has added
                        # l_spml_ts=.true already.
                        logical_state = None
                        location_state = None
                    else:
                        # a regular (old format) time series
                        logical_state = None
                        location_state = config.STATE_SYST_IGNORED
                else:
                    # this is the most common case- no time series
                    logical_state = config.STATE_SYST_IGNORED
                    location_state = config.STATE_SYST_IGNORED

                self.add_setting(config, [obj, "l_spml_ts"],
                                 ".false.", state=logical_state)
                self.add_setting(config, [obj, "spml_ns"],
                                 "0", state=location_state)
                self.add_setting(config, [obj, "spml_ew"],
                                 "0", state=location_state)
                self.add_setting(config, [obj, "spml_bot"],
                                 "0", state=location_state)
                self.add_setting(config, [obj, "spml_top"],
                                 "0", state=location_state)

        # We need to rehash the stash indices using TidyStashTransform().
        # First we add the macro to the front of the path.
        sys.path.insert(1,
                        (os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "vn10.7", "lib", "python", "macros")))
        import stash_indices
        # We must reload the module to clear the vn10.6 stash_indices macro
        # from #1964 above, which is still held in cache.
        reload(stash_indices)
        # Then we quickly remove stash_indices from sys.path to avoid namespace
        # conflicts.
        sys.path.pop(1)
        # Finally we rehash stash indices.
        config, reports = stash_indices.TidyStashTransform().transform(config)
        self.reports += reports

        return config, self.reports

class vn106_t2509(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2509 by James Manners."""

    BEFORE_TAG = "vn10.6_t1233"
    AFTER_TAG = "vn10.6_t2509"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config,["namelist:planet_constants","planet_obs_lon"],"0.0")
        self.add_setting(config,["namelist:planet_constants","planet_obs_lat"],"0.0")
        self.add_setting(config,["namelist:r2lwclnl","l_so2_lw"],".false.")
        self.add_setting(config,["namelist:r2lwclnl","l_so2_lw2"],".false.")
        self.add_setting(config,["namelist:r2swclnl","l_so2_sw"],".false.")
        self.add_setting(config,["namelist:r2swclnl","l_so2_sw2"],".false.")
        self.remove_setting(config,["namelist:run_radiation","is_ncol"])
        self.add_setting(config,["namelist:run_radiation","so2mmr"],"0.0")
        return config, self.reports


class vn106_t2473(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2473 by Chris Smith."""

    BEFORE_TAG = "vn10.6_t2509"
    AFTER_TAG = "vn10.6_t2473"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Idealised initialisation and forcing options (bicyclic configuration)
        self.add_setting(config, ["namelist:idealise","mv_inc_data"],"0.0")
        self.add_setting(config, ["namelist:idealise","mv_inc_height"],"0.0")
        self.add_setting(config, ["namelist:idealise","mv_inc_time"],"0.0")
        self.add_setting(config, ["namelist:idealise","mv_init_data"],"0.0")
        self.add_setting(config, ["namelist:idealise","mv_init_field_type"],"20")
        self.add_setting(config, ["namelist:idealise","mv_init_height"],"0.0")
        self.add_setting(config, ["namelist:idealise","mv_relax_data"],"0.0")
        self.add_setting(config, ["namelist:idealise","mv_relax_height"],"0.0")
        self.add_setting(config, ["namelist:idealise","mv_relax_time"],"0.0")
        self.add_setting(config, ["namelist:idealise","mv_relax_timescale"],"0.0")
        self.add_setting(config, ["namelist:idealise","num_mv_inc_heights"],"0")
        self.add_setting(config, ["namelist:idealise","num_mv_inc_times"],"0")
        self.add_setting(config, ["namelist:idealise","num_mv_init_heights"],"0")
        self.add_setting(config, ["namelist:idealise","num_mv_relax_heights"],"0")
        self.add_setting(config, ["namelist:idealise","num_mv_relax_times"],"0")
        self.add_setting(config, ["namelist:idealise","num_theta_inc_heights"],"0")
        self.add_setting(config, ["namelist:idealise","num_theta_inc_times"],"0")
        self.add_setting(config, ["namelist:idealise","num_theta_init_heights"],"0")
        self.add_setting(config, ["namelist:idealise","num_theta_relax_heights"],"0")
        self.add_setting(config, ["namelist:idealise","num_theta_relax_times"],"0")
        self.add_setting(config, ["namelist:idealise","num_uv_inc_heights"],"0")
        self.add_setting(config, ["namelist:idealise","num_uv_inc_times"],"0")
        self.add_setting(config, ["namelist:idealise","num_uv_init_heights"],"0")
        self.add_setting(config, ["namelist:idealise","num_uv_relax_heights"],"0")
        self.add_setting(config, ["namelist:idealise","num_uv_relax_times"],"0")
        self.add_setting(config, ["namelist:idealise","pressure_balance"],"0")
        self.add_setting(config, ["namelist:idealise","theta_inc_data"],"0.0")
        self.add_setting(config, ["namelist:idealise","theta_inc_field_type"],"11")
        self.add_setting(config, ["namelist:idealise","theta_inc_height"],"0.0")
        self.add_setting(config, ["namelist:idealise","theta_inc_time"],"0.0")
        self.add_setting(config, ["namelist:idealise","theta_init_data"],"0.0")
        self.add_setting(config, ["namelist:idealise","theta_init_field_type"],"10")
        self.add_setting(config, ["namelist:idealise","theta_init_height"],"0.0")
        self.add_setting(config, ["namelist:idealise","theta_relax_data"],"0.0")
        self.add_setting(config, ["namelist:idealise","theta_relax_height"],"0.0")
        self.add_setting(config, ["namelist:idealise","theta_relax_time"],"0.0")
        self.add_setting(config, ["namelist:idealise","theta_relax_timescale"],"0.0")
        self.add_setting(config, ["namelist:idealise","u_inc_data"],"0.0")
        self.add_setting(config, ["namelist:idealise","u_init_data"],"0.0")
        self.add_setting(config, ["namelist:idealise","u_relax_data"],"0.0")
        self.add_setting(config, ["namelist:idealise","uv_inc_height"],"0.0")
        self.add_setting(config, ["namelist:idealise","uv_inc_time"],"0.0")
        self.add_setting(config, ["namelist:idealise","uv_init_height"],"0.0")
        self.add_setting(config, ["namelist:idealise","uv_relax_height"],"0.0")
        self.add_setting(config, ["namelist:idealise","uv_relax_time"],"0.0")
        self.add_setting(config, ["namelist:idealise","uv_relax_timescale"],"0.0")
        self.add_setting(config, ["namelist:idealise","v_inc_data"],"0.0")
        self.add_setting(config, ["namelist:idealise","v_init_data"],"0.0")
        self.add_setting(config, ["namelist:idealise","v_relax_data"],"0.0")
        return config, self.reports

class vn107_t2744(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.6_t2473"
    AFTER_TAG = "vn10.7"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "10.7", forced=True)
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMASTER"])
        if stashmaster_path:
           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",
                                     "vn10.7/ctldata",
                                      stashmaster_path)
           self.change_setting_value(config, ["env", "STASHMASTER"],
                                     stashmaster_path, forced=True)
        return config, self.reports


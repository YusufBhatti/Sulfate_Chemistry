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


class vn104_t1371(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1371 by Glenn Greed"""

    BEFORE_TAG = "vn10.4"
    AFTER_TAG = "vn10.4_t1371"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        stashc_source = self.get_setting_value(config, ["file:STASHC","source"])
        if stashc_source:
            self.change_setting_value(config, ["file:STASHC","source"],
                            stashc_source+" "+'(namelist:exclude_package(:))')
        # Input your macro commands here
        return config, self.reports


class vn104_t1437(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1437 by Ian Boutle."""

    BEFORE_TAG = "vn10.4_t1371"
    AFTER_TAG = "vn10.4_t1437"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        i_dump_year=self.get_setting_value(config,["namelist:headers","i_dump_year"])
        if i_dump_year != None:
              self.add_setting(config,["namelist:headers","i_override_date_time"],"1")
              self.add_setting(config,["namelist:headers","new_date_time"],str(i_dump_year)+",5*0")
        else:
              self.add_setting(config,["namelist:headers","i_override_date_time"],"0")
              self.add_setting(config,["namelist:headers","new_date_time"],"1984,3,10,0,0,0")
        self.remove_setting(config,["namelist:headers","i_dump_year"])
        l_pc2=self.get_setting_value(config,["namelist:run_cloud","l_pc2"])
        if l_pc2==".true.":
              self.add_setting(config,["namelist:run_cloud","i_cld_vn"],"2")
        else:
              self.add_setting(config,["namelist:run_cloud","i_cld_vn"],"1")
        self.remove_setting(config,["namelist:run_cloud","l_pc2"])
        i_bl_vn = int(self.get_setting_value(config,["namelist:run_bl","i_bl_vn"]))
        # beta_evap needs to be trigger ignored if running with a boundary layer scheme,
        #  otherwise it needs to be on.
        if i_bl_vn == 0:
            self.add_setting(config,["namelist:run_bl","beta_evap"],"0.0")
        else:
            self.add_setting(config,["namelist:run_bl","beta_evap"],"0.0", state=config.STATE_SYST_IGNORED)
        fric_number=self.get_setting_value(config,["namelist:idealise","fric_number"])
        self.add_setting(config,["namelist:run_bl","fric_number"],str(fric_number))
        self.remove_setting(config,["namelist:idealise","fric_number"])
        return config, self.reports

class vn104_t336(rose.upgrade.MacroUpgrade):
    """Upgrade macro for ticket #336 by Yongming Tang."""
    BEFORE_TAG = "vn10.4_t1437"
    AFTER_TAG = "vn10.4_t336"
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:jules_vegetation", "l_nitrogen"], ".false.")
        return config, self.reports  

class vn104_t1165(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1165 by Jonathan Wilkinson."""

    BEFORE_TAG = "vn10.4_t336"
    AFTER_TAG = "vn10.4_t1165"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        
        # This is an intentionally blank upgrade macro.
        # It prevents failures in Rose Stem metadata checking
        # by forcing changes to triggering to rose metadata in
        # this ticket to look at the head metadata and not the vn10.4
        # metadata. If this is not applied, Rose stem will complain
        # that the triggering is wrong and the Python stash testmask
        # will fall over.

        return config, self.reports

class vn104_t1598(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1598 by <Mohit Dalvi>."""

    BEFORE_TAG = "vn10.4_t1165"
    AFTER_TAG = "vn10.4_t1598"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Remove item 3-25 from UKCA coupling stash package
        #   if requested with 'UPUKCA' usage profile or 
        #    in the 'UKCA coupling macro package'
        itm_del = '25'
        sec_del = '3'
        to_del = []  # store items (as python objects) to delete
        for obj in config.get_value():
       
           # Check if STASH request
           if re.search(r'namelist:streq', obj):
             
              # get the section, item, Usage profile and package name 
              item = self.get_setting_value(config, [ obj, 'item'])
              section = self.get_setting_value(config, [ obj, 'isec'])
              uprof = self.get_setting_value(config, [ obj, 'use_name'])
              pkg = self.get_setting_value(config, [ obj, 'package'])

              # check if the stash code, usage profile or package names
              #  match
              if section == sec_del and item == itm_del:
                 if re.search(r'UPUKCA', uprof) or \
                    re.search(r'ukca coupling macro', pkg.lower()):
                    # append to the list to be deleted
                    to_del.append(obj)
       
        # loop over the selected streq items and actually delete them
        for delme in to_del:
           self.remove_setting(config, [delme],info='UKCA coupling item'+
               ' no longer required')
         
        return config, self.reports


class vn104_t1310(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1310 by <Mohit Dalvi>."""

    BEFORE_TAG = "vn10.4_t1598"
    AFTER_TAG = "vn10.4_t1310"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Add a scaling factor for Lightning emissions to
        #  the run_ukca namelist. Default value is 1.0
        self.add_setting(config, ["namelist:run_ukca", "lightnox_scale_fac"],
                        value="1.0", forced=True)

        return config, self.reports
        
class vn104_t1645(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1645 by Ian Boutle."""

    BEFORE_TAG = "vn10.4_t1310"
    AFTER_TAG = "vn10.4_t1645"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_bl","l_conv_tke"], ".false.")
        return config, self.reports
        
class vn104_t1167(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1167 by Michael Whitall."""

    BEFORE_TAG = "vn10.4_t1645"
    AFTER_TAG = "vn10.4_t1167"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        self.add_setting(config, ["namelist:temp_fixes",
                                  "l_fix_conv_precip_evap"], ".false.")

        return config, self.reports


class vn104_t1640(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1640 by Joe Mancell"""

    BEFORE_TAG = "vn10.4_t1167"
    AFTER_TAG = "vn10.4_t1640"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        recon_sci_items = ("q_min", "w_zero_start", "w_zero_end",
              "use_smc_stress", "polar_check", "l_canopy_snow_throughfall",
              "l_force_relayer", "coast_adj_method", "l_regularize_landice",
              "snow_landice_min", "snow_icefree_max", "l_adj_t_soil",
              "l_use_zero_frac_tile_temp", "l_snow_tile_gbm_ancil_fix")
        recon_tech_items = ("dump_pack", "input_dump_type", "reset_data_time",
              "l_trans", "var_recon", "select_output_fields", "ainitial",        
              "transp", "l_validity_lookup_u", "l_rcf_init_flexi")

        for item in recon_sci_items:
            item_val = self.get_setting_value(config, ["namelist:recon", item])
            self.add_setting(config, ["namelist:recon_science", item], item_val)
        for item in recon_tech_items:
            item_val = self.get_setting_value(config, ["namelist:recon", item])
            self.add_setting(config, ["namelist:recon_technical", item], item_val)
        self.remove_setting(config, ["namelist:recon",])
        recona_source = self.get_setting_value(config, ["file:RECONA","source"])
        if recona_source:
            recona_source = re.sub(r'namelist:recon', r'namelist:recon_technical '
                                   'namelist:recon_science', recona_source)
        self.change_setting_value(config, ["file:RECONA", "source"], recona_source)

        return config, self.reports

class vn104_t1375(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1375 by Robin Smith."""

    BEFORE_TAG = "vn10.4_t1640"
    AFTER_TAG = "vn10.4_t1375"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:jules_surface",
                                  "l_elev_land_ice"],
                         ".false.")
        self.add_setting(config, ["namelist:jules_soil",
                                  "dzsoil_elev"],
                         "0.")
        self.add_setting(config, ["namelist:jules_surface_types",
                                  "elev_ice"],
                         "0")
        
        return config, self.reports


class vn104_t1818(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1818 by J. M. Edwards."""

    BEFORE_TAG = "vn10.4_t1375"
    AFTER_TAG = "vn10.4_t1818"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:jules_radiation", "l_niso_direct"], ".false.")
        return config, self.reports


class vn104_t1836(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1836 by J.M. Edwards."""
    
    BEFORE_TAG = "vn10.4_t1818"
    AFTER_TAG = "vn10.4_t1836"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM app configuration."""
        self.add_setting(config, ["namelist:jules_sea_seaice", "l_iceformdrag_lupkes"], ".false.")
        self.add_setting(config, ["namelist:jules_sea_seaice", "l_stability_lupkes"], ".false.")
        self.add_setting(config, ["namelist:jules_sea_seaice", "h_freeboard_min"], "0.286")
        self.add_setting(config, ["namelist:jules_sea_seaice", "h_freeboard_max"], "0.534")
        self.add_setting(config, ["namelist:jules_sea_seaice", "beta_floe"], "1.0")
        self.add_setting(config, ["namelist:jules_sea_seaice", "d_floe_min"], "8.0")
        self.add_setting(config, ["namelist:jules_sea_seaice", "d_floe_max"], "300.0")
        self.add_setting(config, ["namelist:jules_sea_seaice", "ss_floe"], "0.5")
        self.add_setting(config, ["namelist:jules_sea_seaice", "ce_floe"], "0.222")
        return config, self.reports


class vn104_t1705(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1705 by AdrianLock."""

    BEFORE_TAG = "vn10.4_t1836"
    AFTER_TAG = "vn10.4_t1705"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config,
            ["namelist:run_bl", "l_new_kcloudtop"], ".false.")
        return config, self.reports


class vn104_t1244(rose.upgrade.MacroUpgrade):
    """Upgrade macro for ticket #1244 by David Walters."""

    BEFORE_TAG = "vn10.4_t1705"
    AFTER_TAG = "vn10.4_t1244"

    def upgrade(self, config, meta_config=None):
          """ Set land_field to -99 for apps reading a """
          """ land-sea mask from an ancillary file.    """

          # First, we check the value of select_output_fields.
          # If this is non-zero, then the reconfiguration is 
          # running an "all input fields" regridding, which ignores
          # the items namelist and hence won't read a land mask.
          select_output_fields = self.get_setting_value(
                config, ["namelist:recon_technical", "select_output_fields"])
          
          if select_output_fields == "0":
                
                # Loop over app namelist objects
                for obj in config.get_value():
                      # Test if this is an items namelist
                      if re.search(r'namelist:items', obj):
                            # Convert string of the stashcode list into 
                            # a list of stashcode strings
                            stashcodes = self.get_setting_value(config, [ obj, 'stash_req'])
                            source = self.get_setting_value(config, [obj, 'source'])
                            stashcodes = stashcodes.split(',')

                            # Test to see if the stashcodes contain the 
                            # land sea mask (30) and source is ancil (2)
                            # If this case is set we can set land_field to
                            # -99 as this will now be read from the
                            # ancillary by the reconfiguration code
                            stashcode_lsm="30"
                            source_ancil="2"
                            if ( stashcode_lsm in stashcodes and source == source_ancil ) :
                                  self.change_setting_value(
                                        config, ["namelist:nlsizes", "land_field"], "-99")
          return config, self.reports

class vn104_t1841(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1841 by nicksavage."""
    
    BEFORE_TAG = "vn10.4_t1244"
    AFTER_TAG = "vn10.4_t1841"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # this is an empty upgrade macro so that all apps will be 
        # forcibly updated to the new metadata version
        return config, self.reports


class vn104_t1882(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1882 by Anna Harper."""

    BEFORE_TAG = "vn10.4_t1841"
    AFTER_TAG = "vn10.4_t1882"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:jules_vegetation", "l_scale_resp_pm"], ".false.")
        self.add_setting(config, ["namelist:jules_vegetation", "l_stem_resp_fix"], ".false.")
        
        return config, self.reports

class vn104_t1638(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1638 by annemccabe."""

    BEFORE_TAG = "vn10.4_t1882"
    AFTER_TAG = "vn10.4_t1638"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.remove_setting(config,["namelist:run_stochastic","kay_gwave_min"])
        self.remove_setting(config,["namelist:run_stochastic","kay_gwave_max"])
        self.remove_setting(config,["namelist:run_stochastic","gwd_frc_min"])
        self.remove_setting(config,["namelist:run_stochastic","gwd_frc_max"])
        self.add_setting(config,
            ["namelist:temp_fixes","l_fix_rp_shock_amp"], ".false.")
        self.add_setting(config,
            ["namelist:run_stochastic","a_ent_shr_rp"], "5.0")
        self.add_setting(config,
            ["namelist:run_stochastic","a_ent_shr_rp_max"], "8.7")
        return config, self.reports


class vn104_t1729(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1729 by AdrianLock."""

    BEFORE_TAG = "vn10.4_t1638"
    AFTER_TAG = "vn10.4_t1729"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, 
            ["namelist:temp_fixes", "l_fix_ustar_dust"], ".false.")
        return config, self.reports

class vn104_t1713(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1713 by Michele Guidolin."""

    BEFORE_TAG = "vn10.4_t1729"
    AFTER_TAG = "vn10.4_t1713"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, [ "namelist:acp", "l_lhn_1a" ],
                         ".false.")
        return config, self.reports


class vn104_t1702(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1702 by Jane Mulcahy."""

    BEFORE_TAG = "vn10.4_t1713"
    AFTER_TAG = "vn10.4_t1702"

    def upgrade(self, config, meta_config=None):
        """Add new parameter knl to the JULES namelist."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_radiation","l_use_liu_spec"], ".false.")
        return config, self.reports


class vn104_t1782(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1782 by Karina Williams."""

    BEFORE_TAG = "vn10.4_t1702"
    AFTER_TAG = "vn10.4_t1782"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        npft = int(self.get_setting_value(config, ["namelist:jules_surface_types", "npft"]))
        self.add_setting(config, ["namelist:jules_pftparm", "can_struct_a_io"], ','.join(['1.0'] * npft))
       
        return config, self.reports


class vn104_t1376(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1376 by robinsmith."""

    BEFORE_TAG = "vn10.4_t1782"
    AFTER_TAG = "vn10.4_t1376"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:jules_snow",
                                  "aicemax"],
                         "0.78, 0.36")
        self.add_setting(config, ["namelist:jules_snow",
                                  "rho_firn_albedo"],
                         "550.")
        return config, self.reports


class vn104_t80(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #80 by andrewbushell."""

    BEFORE_TAG = "vn10.4_t1376"
    AFTER_TAG = "vn10.4_t80"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # RUN_gwd namelist - settings Standard 1A Default
        self.add_setting(config, ["namelist:run_gwd","i_ussp_vn"], "1")
        self.add_setting(config, ["namelist:run_gwd","wavelstar"], "4300.00")
        is_add_cgw = ".false."
        self.add_setting(config, ["namelist:run_gwd","l_add_cgw"],
                         value = is_add_cgw)

        # only needed IF l_add_cgw is true
        if (".true." in is_add_cgw):
            self.add_setting(config, ["namelist:run_gwd","cgw_scale_factor"], "1.0000")
        #
        # trigger ignore IF l_add_cgw switched off
        else:
            self.add_setting(config, ["namelist:run_gwd","cgw_scale_factor"], "0.0000")

        return config, self.reports

class vn104_t1869(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1869 by Anna Harper."""

    BEFORE_TAG = "vn10.4_t80"
    AFTER_TAG = "vn10.4_t1869"

    def upgrade(self, config, meta_config=None):
        """Add new parameter knl to the JULES namelist."""
        # Input your macro commands here
        npft = int(self.get_setting_value(config, ["namelist:jules_surface_types", "npft"]))
        crm = self.get_setting_value(config, ["namelist:jules_vegetation", "can_rad_mod"])
        kn_string = self.get_setting_value(config, ["namelist:jules_pftparm", "kn_io"])
        if kn_string is None:
            # if kn_io is not present, set a default:
            kn_string = ','.join(['0.20'] * npft)
        kn = kn_string.split(',')

        # Set kn and knl
        if crm == '6':
           # CRM6 uses knl, use previous value of Kn
           self.add_setting(config, ["namelist:jules_pftparm", "knl_io"], kn_string)
           self.change_setting_value(config, ["namelist:jules_pftparm", "kn_io"], ','.join(['0.78'] * npft), forced=True)
           if any(float(a) > 0.5 for a in kn ):
              kn_warning='''
                    Warning! A smaller value of knl is required with can_rad_mod=6 
                    than with previous canopy radiation schemes. Recommended value 
                    is 0.20.
                    '''
              self.add_report(info=kn_warning, is_warning=True)
        else:
           self.add_setting(config, ["namelist:jules_pftparm", "knl_io"], ','.join(['0.20'] * npft))
           # Does nothing if kn_io is already present:
           self.add_setting(config, ["namelist:jules_pftparm", "kn_io"], ','.join(['0.78'] * npft))

        # Add l_leaf_n_resp_fix
        self.add_setting(config, ["namelist:jules_vegetation", "l_leaf_n_resp_fix"], ".false.")
        return config, self.reports

class vn104_t1770(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1770 by Ben Johnson."""

    BEFORE_TAG = "vn10.4_t1869"
    AFTER_TAG = "vn10.4_t1770"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        vars_to_update = {
                 'ukcaaclw' : 'nml_ac_lw',
                 'ukcaacsw' : 'nml_ac_sw',
                 'ukcaanlw' : 'nml_an_lw',
                 'ukcaansw' : 'nml_an_sw',
                 'ukcacrlw' : 'nml_cr_lw',
                 'ukcacrsw' : 'nml_cr_sw' 
                         }
                           
        for item in vars_to_update.keys():
            value = self.get_setting_value(config, ["namelist:run_ukca", item])
            # Search for a path which looks like it's in $UMDIR/vn$VN
            search = 'ctldata/UKCA/radaer/' + vars_to_update[item]
            if search in value:
                new_value = re.sub(r'radaer/nml', r'radaer/ga3-7/nml' , value)
                self.change_setting_value(config, [ "namelist:run_ukca", item], new_value)
       
        return config, self.reports

class vn104_t1607(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1607 by James Manners."""

    BEFORE_TAG = "vn10.4_t1770"
    AFTER_TAG = "vn10.4_t1607"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:planet_constants","planet_albedo"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_emissivity"], "1.0")
        self.add_setting(config, ["namelist:r2lwclnl","l_h2o_lw"], ".true.")
        self.add_setting(config, ["namelist:r2lwclnl","l_h2o_lw2"], ".true.")
        self.add_setting(config, ["namelist:r2swclnl","l_h2o_sw"], ".true.")
        self.add_setting(config, ["namelist:r2swclnl","l_h2o_sw2"], ".true.")
        self.add_setting(config, ["namelist:r2lwclnl","l_o3_lw"], ".true.")
        self.add_setting(config, ["namelist:r2lwclnl","l_o3_lw2"], ".true.")
        self.add_setting(config, ["namelist:r2swclnl","l_o3_sw"], ".true.")
        self.add_setting(config, ["namelist:r2swclnl","l_o3_sw2"], ".true.")
        return config, self.reports


class vn104_t1346(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1346 by Andy Malcolm"""

    BEFORE_TAG = "vn10.4_t1607"
    AFTER_TAG = "vn10.4_t1346"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        self.rename_setting(config, ["namelist:run_bl", "bl_segment_size"], ["namelist:tuning_segments", "bl_segment_size" ])
        self.rename_setting(config, ["namelist:run_gwd", "gw_seg_size"], ["namelist:tuning_segments", "gw_seg_size" ])
        self.rename_setting(config, ["namelist:run_precip", "precip_segment_size"], ["namelist:tuning_segments", "precip_segment_size" ])
        self.rename_setting(config, ["namelist:run_convection", "a_convect_seg_size"], ["namelist:tuning_segments", "a_convect_seg_size" ])
        self.rename_setting(config, ["namelist:run_convection", "a_convect_segments"], ["namelist:tuning_segments", "a_convect_segments" ])
        self.rename_setting(config, ["namelist:run_radiation", "a_sw_seg_size"], ["namelist:tuning_segments", "a_sw_seg_size" ])
        self.rename_setting(config, ["namelist:run_radiation", "a_sw_segments"], ["namelist:tuning_segments", "a_sw_segments" ])
        self.rename_setting(config, ["namelist:run_radiation", "a_lw_seg_size"], ["namelist:tuning_segments", "a_lw_seg_size" ])
        self.rename_setting(config, ["namelist:run_radiation", "a_lw_segments"], ["namelist:tuning_segments", "a_lw_segments" ])

        # Insert it after iau_nl, which is read just before it in the code
        model=self.get_setting_value(config,["namelist:model_domain","model_type"])
        if model=="5":
            src_cntlatm = (self.get_setting_value(config, ["file:CNTLATM","source"])).split()
            index = src_cntlatm.index("namelist:iau_nl")
            src_cntlatm.insert(index+1,"namelist:tuning_segments")
            self.change_setting_value(config, ["file:CNTLATM","source"]," ".join(src_cntlatm))
        else:
            src_atmoscntl = (self.get_setting_value(config, ["file:ATMOSCNTL","source"])).split()
            index = src_atmoscntl.index("namelist:iau_nl")
            src_atmoscntl.insert(index+1,"namelist:tuning_segments")
            self.change_setting_value(config, ["file:ATMOSCNTL","source"]," ".join(src_atmoscntl))

        return config, self.reports


class vn104_t1676(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1676 by Marc Stringer."""

    BEFORE_TAG = "vn10.4_t1346"
    AFTER_TAG = "vn10.4_t1676"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here.
        # l_oasis_ocean is true if there is oasis3/-mct coupling 
        coupler=self.get_setting_value(config, ["env","COUPLER"])
        if re.search("oasis", coupler):
              # Coupler is set
              self.add_setting(config, ["namelist:coupling_control",
                                        "l_oasis_ocean"],".true.")
        else:
              # Coupler is not set
              self.add_setting(config, ["namelist:coupling_control",
                                        "l_oasis_ocean"],".false.")
        return config, self.reports

class vn105_t1963(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.4_t1676"
    AFTER_TAG = "vn10.5"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "10.5", forced=True)
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMSTR"])
        if stashmaster_path:
           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",
                                     "vn10.5/ctldata",
                                      stashmaster_path)
           self.change_setting_value(config, ["env", "STASHMSTR"],
                                     stashmaster_path, forced=True)
        return config, self.reports


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



class vn110_t3799 (rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3799 by Rachel Stratton."""

    BEFORE_TAG = "vn11.0"
    AFTER_TAG = "vn11.0_t3799"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
	# No commands just altering triggering
        return config, self.reports


class vn110_t3620(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3620 by Ian Boutle."""

    BEFORE_TAG = "vn11.0_t3799"
    AFTER_TAG = "vn11.0_t3620"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config,["namelist:run_convection","l_jules_flux"],
                         ".false.")
        return config, self.reports


class vn110_t3025(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3025 by Alan J Hewitt."""

    BEFORE_TAG = "vn11.0_t3620"
    AFTER_TAG = "vn11.0_t3025"

    def upgrade(self, config, meta_config=None):
        """
        Add new control logicals to namelist:run_glomap_aeroclim
        """
        
        # Input your macro commands here
        
        self.add_setting(config, ["namelist:run_glomap_aeroclim",
                                  "l_glomap_clim_arg_act"],
                         value=".false.")
        
        self.add_setting(config, ["namelist:run_glomap_aeroclim",
                                  "l_glomap_clim_aie1"],
                         value=".false.")
        
        self.add_setting(config, ["namelist:run_glomap_aeroclim",
                                  "l_glomap_clim_aie2"],
                         value=".false.")
        
        return config, self.reports


class vn110_t3803(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3803 by Michele Guidolin."""

    BEFORE_TAG = "vn11.0_t3025"
    AFTER_TAG = "vn11.0_t3803"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        l_swap_bounds_ddt = self.get_setting_value(config, ["namelist:nlst_mpp","l_swap_bounds_ddt"])

        if l_swap_bounds_ddt == ".false.":
            self.add_setting(config,[ "namelist:nlst_mpp", "swap_bounds_method_3d_large" ],"11") 
            self.add_setting(config,[ "namelist:nlst_mpp", "swap_bounds_method_3d_small" ],"10") 
        else:
            self.add_setting(config,[ "namelist:nlst_mpp", "swap_bounds_method_3d_large" ],"21") 
            self.add_setting(config,[ "namelist:nlst_mpp", "swap_bounds_method_3d_small" ],"20") 

        self.remove_setting(config,["namelist:nlst_mpp", "l_swap_bounds_ddt"])
        self.remove_setting(config,["namelist:nlst_mpp", "l_ddt_cornerless"])

        return config, self.reports


class vn110_t3806(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3806 by Luke Abraham."""

    BEFORE_TAG = "vn11.0_t3803"
    AFTER_TAG = "vn11.0_t3806"

    def upgrade(self, config, meta_config=None):
        """Rename l_ukca_linox_logp to l_ukca_linox_scaling
           as logical now also controls whether to use the 
           same NOx production efficiency for C2C and C2G 
           flashes, or not."""

        self.rename_setting(config, ['namelist:run_ukca', 'l_ukca_linox_logp'],
                            ['namelist:run_ukca', 'l_ukca_linox_scaling'])

        return config, self.reports

class vn110_t3056(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3056 by James Manners."""

    BEFORE_TAG = "vn11.0_t3806"
    AFTER_TAG = "vn11.0_t3056"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:planet_constants", "stellar_radius"], "6.957e+08")
        self.add_setting(config, ["namelist:r2swclnl", "l_spherical_solar"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl", "l_spherical_solar_2"], ".false.")
        return config, self.reports

class vn110_t3920(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3920 by Rachel Stratton."""

    BEFORE_TAG = "vn11.0_t3056"
    AFTER_TAG = "vn11.0_t3920"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
	# Need to know whether geostrophic forcing on
        l_geo_for = self.get_setting_value(config,
                                ["namelist:idealised", "l_geo_for"])
	
        if (l_geo_for==".true."):
	    # Need current u_geo and v_geo values to setup new variables
            uu = self.get_setting_value(config,
                                ["namelist:idealised", "u_geo"])
            vv = self.get_setting_value(config,
                                ["namelist:idealised", "v_geo"])
            self.remove_setting(config,["namelist:idealised", "u_geo"])
            self.remove_setting(config,["namelist:idealised", "v_geo"])
            # one time so forcing as original unchanging with time 
            self.add_setting(config, ["namelist:idealised", "num_uv_geo_times"],
                         "1")
            self.add_setting(config, ["namelist:idealised", "uv_geo_time"],
                         "0.0")
            # Need two heights 0.0 and model top to get same impact
	    # So going to use a top height of 100.km which is probably
	    # higher than all idealised tops. Most are ~40km
	    self.add_setting(config, ["namelist:idealised", "num_uv_geo_heights"],
                         "2")
            self.add_setting(config, ["namelist:idealised", "uv_geo_height"],
	                         "0.0,1.0e5")
            self.add_setting(config, ["namelist:idealised", "u_geo_data"],
	                         uu+","+uu)
            self.add_setting(config, ["namelist:idealised", "v_geo_data"],
	                         vv+","+vv)
	else:
            # No problems add new variables with null values
            self.remove_setting(config,["namelist:idealised", "u_geo"])
            self.remove_setting(config,["namelist:idealised", "v_geo"])
 	      	
            self.add_setting(config, ["namelist:idealised", "num_uv_geo_times"],
                         "0")
	    self.add_setting(config, ["namelist:idealised", "num_uv_geo_heights"],
                         "0")
            self.add_setting(config, ["namelist:idealised", "uv_geo_time"],
                         "0.0")
            self.add_setting(config, ["namelist:idealised", "uv_geo_height"],
	                         "0.0")
            self.add_setting(config, ["namelist:idealised", "u_geo_data"],
	                         "0.0")
            self.add_setting(config, ["namelist:idealised", "v_geo_data"],
	                         "0.0")
        return config, self.reports


class vn110_t1250(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1250 by Alan J Hewitt."""

    BEFORE_TAG = "vn11.0_t3920"
    AFTER_TAG = "vn11.0_t1250"

    def upgrade(self, config, meta_config=None):
        """
        Temporary logical to protect configurations from dry deposition changes.
        
        With new logical equals true, dry deposition velocities are set for
        HCl, HOCl, HBr, HOBr, H2SO4, MeOH and Sec_Org.
        
        With new logical equals true, dry deposition velocities for 9 tiles 
        will also be consistant with 13/17/27 tiles.
        
        New logical to be reviewed / retired in ticket #3997.
        """
        # Input your macro commands here
        
        self.add_setting(config,["namelist:temp_fixes","l_fix_improve_drydep"],
                         value=".false.")
        
        return config, self.reports


class vn110_t3018(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3018 by Chris Dearden."""

    BEFORE_TAG = "vn11.0_t1250"
    AFTER_TAG = "vn11.0_t3018"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, [ "namelist:run_ukca", "l_ukca_debug_asad" ],
                         str(".false."))
        return config, self.reports


class vn110_t4038(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4038 by Douglas Clark."""

    BEFORE_TAG = "vn11.0_t3018"
    AFTER_TAG = "vn11.0_t4038"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:temp_fixes", "l_fix_wind_snow"], ".false.")
        return config, self.reports


class vn110_t2293(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2293 by Gabriel Rooney."""

    BEFORE_TAG = "vn11.0_t4038"
    AFTER_TAG = "vn11.0_t2293"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_convection", "cnv_cold_pools"], "0")
        return config, self.reports


class vn110_t1373(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1373 by Gerd Folberth."""

    BEFORE_TAG = "vn11.0_t2293"
    AFTER_TAG = "vn11.0_t1373"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Input your macro commands here
        # Set default value of l_inferno from false to false
        l_inferno = self.get_setting_value(config, ["namelist:jules_vegetation", "l_inferno"])
        if l_inferno:
            self.change_setting_value(config, ["namelist:jules_vegetation", "l_inferno"], ".false.")
        else:
            self.add_setting(config,["namelist:jules_vegetation", "l_inferno"], ".false.")

        # Set default value of ignition_method to 1
        ignition_method = self.get_setting_value(config, ["namelist:jules_vegetation", "ignition_method"])
        if ignition_method:
            self.change_setting_value(config, ["namelist:jules_vegetation", "ignition_method"], "1")
        else:
            self.add_setting(config,["namelist:jules_vegetation", "ignition_method"], "1")

        # Set default value of l_trif_fire from false to false
        l_trif_fire = self.get_setting_value(config, ["namelist:jules_vegetation", "l_trif_fire"])
        if l_trif_fire:
            self.change_setting_value(config, ["namelist:jules_vegetation", "l_trif_fire"], ".false.")
        else:
            self.add_setting(config,["namelist:jules_vegetation", "l_trif_fire"], ".false.")

        """Upgrade the JULES pft parameters"""
        RMDI = str(-2**30)
        npft = int(self.get_setting_value(config, ["namelist:jules_surface_types", "npft"]))
        ncpft_in = self.get_setting_value(config, ["namelist:jules_surface_types", "ncpft"])

        if ncpft_in:
            ncpft = int(ncpft_in)
        else:
            ncpft = 0

        if npft == 5:
            # Add sensible values for 5 pfts
            fef_co2_io = ["1631", "1576", "1576", "1654", "1576"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "fef_co2_io"], "  ".join(fef_co2_io[:npft]))
            fef_co_io = ["100", "106", "106", "64", "106"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "fef_co_io"], "  ".join(fef_co_io[:npft]))
            fef_ch4_io = ["6.8", "4.8", "4.8", "2.4", "4.8"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "fef_ch4_io"], "  ".join(fef_ch4_io[:npft]))
            fef_nox_io = ["2.55", "3.24", "3.24", "2.49", "3.24"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "fef_nox_io"], "  ".join(fef_nox_io[:npft]))
            fef_so2_io = ["0.40", "0.40", "0.40", "0.48", "0.40"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "fef_so2_io"], "  ".join(fef_so2_io[:npft]))
            fef_oc_io = ["4.3", "9.1", "9.1", "3.2", "9.1"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "fef_oc_io"], "  ".join(fef_oc_io[:npft]))
            fef_bc_io = ["0.56", "0.56", "0.56", "0.47", "0.56"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "fef_bc_io"], "  ".join(fef_bc_io[:npft]))
            ccleaf_min_io = ["0.8", "0.8", "0.8", "0.8", "0.8", "0.56"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "ccleaf_min_io"], "  ".join(ccleaf_min_io[:npft]))
            ccleaf_max_io = ["1.0", "1.0", "1.0", "1.0", "1.0", "0.56"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "ccleaf_max_io"], "  ".join(ccleaf_max_io[:npft]))
            ccwood_min_io = ["0.0", "0.0", "0.0", "0.0", "0.0", "0.56"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "ccwood_min_io"], "  ".join(ccwood_min_io[:npft]))
            ccwood_max_io = ["0.4", "0.4", "0.4", "0.4", "0.4"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "ccwood_max_io"], "  ".join(ccwood_max_io[:npft]))
            avg_ba_io = ["0.6E6","0.6E6","1.4E6","1.4E6","1.2E6"] + ([RMDI] * (npft - 5))
            self.add_setting(config, ["namelist:jules_pftparm", "avg_ba_io"], "  ".join(avg_ba_io[:npft]))
        if npft == 9:
            # Add sensible values for 9 pfts
            if ncpft == 0:
                fef_co2_io = ["1643","1637","1643","1637","1489","1637","1686","1637","1489"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_co2_io"], "  ".join(fef_co2_io[:npft]))
                fef_co_io = ["93","89","93","89","127","89","63","89","127"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_co_io"], "  ".join(fef_co_io[:npft]))
                fef_ch4_io = ["5.07","3.92","5.07","3.92","5.96","3.92","1.94","3.92","5.96"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_ch4_io"], "  ".join(fef_ch4_io[:npft]))
                fef_nox_io = ["2.55","2.51","2.55","2.51","0.90","2.51","3.9","2.51","0.90"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_nox_io"], "  ".join(fef_nox_io[:npft]))
                fef_so2_io = ["0.40","0.40","0.40","0.40","0.40","0.40","0.48","0.40","0.40"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_so2_io"], "  ".join(fef_so2_io[:npft]))
                fef_oc_io = ["4.71","8.2","4.71","8.2","8.2","8.2","2.62","8.2","8.2"] + ([RMDI] * (npft - 5))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_oc_io"], "  ".join(fef_oc_io[:npft]))
                fef_bc_io = ["0.52","0.56","0.52","0.56","0.56","0.56","0.37","0.56","0.56"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_bc_io"], "  ".join(fef_bc_io[:npft]))
                ccleaf_min_io = ["0.8","0.8","0.8","0.8","0.8","0.8","0.8","0.8","0.8"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "ccleaf_min_io"], "  ".join(ccleaf_min_io[:npft]))
                ccleaf_max_io = ["1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "ccleaf_max_io"], "  ".join(ccleaf_max_io[:npft]))
                ccwood_min_io = ["0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "ccwood_min_io"], "  ".join(ccwood_min_io[:npft]))
                ccwood_max_io = ["0.4","0.4","0.4","0.4","0.4","0.4","0.4","0.4","0.4"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "ccwood_max_io"], "  ".join(ccwood_max_io[:npft]))
                avg_ba_io = ["0.6E6","0.6E6","0.6E6","0.6E6","0.6E6","1.4E6","1.4E6","1.2E6","1.2E6"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "avg_ba_io"], "  ".join(avg_ba_io[:npft]))
            elif ncpft == 4:
                fef_co2_io = ["1631", "1576", "1576", "1654", "1576","1576","1576","1654","1576"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_co2_io"], "  ".join(fef_co2_io[:npft]))
                fef_co_io = ["100", "106", "106", "64", "106","106","106","64","106"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_co_io"], "  ".join(fef_co_io[:npft]))
                fef_ch4_io = ["6.8", "4.8", "4.8", "2.4", "4.8","4.8","4.8","2.4","4.8"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_ch4_io"], "  ".join(fef_ch4_io[:npft]))
                fef_nox_io = ["2.55", "3.24", "3.24", "2.49", "3.24","3.24","3.24","2.49","3.24"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_nox_io"], "  ".join(fef_nox_io[:npft]))
                fef_so2_io = ["0.40", "0.40", "0.40", "0.48", "0.40","0.40","0.40","0.48","0.40"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_so2_io"], "  ".join(fef_so2_io[:npft]))
                fef_oc_io = ["4.3", "9.1", "9.1", "3.2", "9.1","9.1","9.1","3.2","9.1"] + ([RMDI] * (npft - 5))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_oc_io"], "  ".join(fef_oc_io[:npft]))
                fef_bc_io = ["0.56", "0.56", "0.56", "0.47", "0.56","0.56","0.56","0.47","0.56"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "fef_bc_io"], "  ".join(fef_bc_io[:npft]))
                ccleaf_min_io = ["0.8","0.8","0.8","0.8","0.8","0.8","0.8","0.8","0.8"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "ccleaf_min_io"], "  ".join(ccleaf_min_io[:npft]))
                ccleaf_max_io = ["1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "ccleaf_max_io"], "  ".join(ccleaf_max_io[:npft]))
                ccwood_min_io = ["0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0","0.0"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "ccwood_min_io"], "  ".join(ccwood_min_io[:npft]))
                ccwood_max_io = ["0.4","0.4","0.4","0.4","0.4","0.4","0.4","0.4","0.4"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "ccwood_max_io"], "  ".join(ccwood_max_io[:npft]))
                avg_ba_io = ["0.6E6","0.6E6","1.4E6","1.4E6","1.2E6","1.4E6","1.4E6","1.4E6","1.4E6"] + ([RMDI] * (npft - 9))
                self.add_setting(config, ["namelist:jules_pftparm", "avg_ba_io"], "  ".join(avg_ba_io[:npft]))

        else:
            # Set everything as RMDI
            fef_co2_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "fef_co2_io"], "  ".join(fef_co2_io[:npft]))
            fef_co_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "fef_co_io"], "  ".join(fef_co_io[:npft]))
            fef_ch4_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "fef_ch4_io"], "  ".join(fef_ch4_io[:npft]))
            fef_nox_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "fef_nox_io"], "  ".join(fef_nox_io[:npft]))
            fef_so2_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "fef_so2_io"], "  ".join(fef_so2_io[:npft]))
            fef_oc_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "fef_oc_io"], "  ".join(fef_oc_io[:npft]))
            fef_bc_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "fef_bc_io"], "  ".join(fef_bc_io[:npft]))
            ccleaf_min_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "ccleaf_min_io"], "  ".join(ccleaf_min_io[:npft]))
            ccleaf_max_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "ccleaf_max_io"], "  ".join(ccleaf_max_io[:npft]))
            ccwood_min_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "ccwood_min_io"], "  ".join(ccwood_min_io[:npft]))
            ccwood_max_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "ccwood_max_io"], "  ".join(ccwood_max_io[:npft]))
            avg_ba_io = [RMDI] * npft
            self.add_setting(config, ["namelist:jules_pftparm", "avg_ba_io"], "  ".join(avg_ba_io[:npft]))

        return config, self.reports

class vn110_t3963(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3963 by Chris Smith."""

    BEFORE_TAG = "vn11.0_t1373"
    AFTER_TAG = "vn11.0_t3963"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # === Read variables affecting dry air conservation ===
        model_type = self.get_setting_value(config,
                       ['namelist:model_domain', 'model_type'])
        l_fix_mass = self.get_setting_value(config,
                       ['namelist:run_dyn', 'l_fix_mass'])
        l_priestley_correct_thetav = self.get_setting_value(config,
                       ['namelist:run_sl', 'l_priestley_correct_thetav'])

        # === Replace l_fix_mass with conserve_dry_mass ===
        if (l_fix_mass == ".true."):
            if (model_type == "2"):
                conserve_dry_mass = "1"  # Constant factor for LAM
            else:                        # Linear factor for other domains
                if (l_priestley_correct_thetav == ".true."):
                    conserve_dry_mass = "3"  # Preserve gravitational energy
                else:
                    conserve_dry_mass = "2"  # Preserve internal and
                                             # gravitational energy
        else:
            conserve_dry_mass = "0"

        self.add_setting(config, ["namelist:run_dyn", 
                                  "conserve_dry_mass"], conserve_dry_mass)
        self.remove_setting(config, ["namelist:run_dyn", "l_fix_mass"])

        return config, self.reports


class vn110_t3545(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3545 by Rachel Stratton."""

    BEFORE_TAG = "vn11.0_t3963"
    AFTER_TAG = "vn11.0_t3545"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
	# Moving l_cartesian from recon_idealised to model_domain
	# First get value - note non-idealised values may be !!l_cartesian=.false.
        # or have no value.
        self.rename_setting(config,["namelist:recon_idealised", "l_cartesian"],
	                           ["namelist:model_domain", "l_cartesian"])
	
        return config, self.reports


class vn110_t3956(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3956 by Joe Mancell."""

    BEFORE_TAG = "vn11.0_t3545"
    AFTER_TAG = "vn11.0_t3956"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.remove_setting(config, ["namelist:temp_fixes", 
                                     "l_ignore_error_ancil_struct"])

        return config, self.reports


class vn110_t4029(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4029 by Alistair Sellar."""

    BEFORE_TAG = "vn11.0_t3956"
    AFTER_TAG = "vn11.0_t4029"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        return config, self.reports


class vn110_t3973(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3973 by Richard Barnes."""

    BEFORE_TAG = "vn11.0_t4029"
    AFTER_TAG = "vn11.0_t3973"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.remove_setting(config,["namelist:run_gwd", "l_ussp_opaque"])
        return config, self.reports

class vn110_t3089(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3089 by Rachel Stratton."""

    BEFORE_TAG = "vn11.0_t3973"
    AFTER_TAG = "vn11.0_t3089"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:recon_idealised", "l_saturate_bubble"], ".false.")
	
        return config, self.reports

class vn110_t4040(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4040 by David Walters and Stuart Whitehouse."""

    BEFORE_TAG = "vn11.0_t3089"
    AFTER_TAG = "vn11.0_t4040"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:nlstcall", "l_nrun_as_crun"], 
                         '.false.')
        return config, self.reports


class vn110_t3772(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3772 by Jonathan Wilkinson."""

    BEFORE_TAG = "vn11.0_t4040"
    AFTER_TAG = "vn11.0_t3772"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # No change, just need to increment metadata tag
        return config, self.reports


class vn110_t4018 (rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4018 by Alan J Hewitt."""

    BEFORE_TAG = "vn11.0_t3772"
    AFTER_TAG = "vn11.0_t4018"

    def upgrade(self, config, meta_config=None):
        """
        No commands just altering triggering.
        """
        # Input your macro commands here
	# No commands just altering triggering
        return config, self.reports


class vn110_t4023(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4023 by Michele Guidolin."""

    BEFORE_TAG = "vn11.0_t4018"
    AFTER_TAG = "vn11.0_t4023"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        
        nb_swap_bounds_method = self.get_setting_value(config, ["namelist:nlst_mpp","nb_swap_bounds_method"])
        
        # If the flag does not exist set the method to be "ddt_aao" (20)
        
        if nb_swap_bounds_method is None:
            self.add_setting(config,[ "namelist:nlst_mpp", "nb_swap_bounds_method" ],"20") 
            nb_swap_bounds_method="20"

        # Get the model domain
        model_domain = self.get_setting_value(config, ["namelist:model_domain", "model_type"])

        #  If it is a single column model (5) or site specific
        #  forecast model (6), use blocking (0) swap bounds method.
        if model_domain == "5" or  model_domain == "6":
            self.change_setting_value(config, ["namelist:nlst_mpp","nb_swap_bounds_method"], "0")
        else:
            # If the method is not a blocking (0), then set it to be "ddt_aao" (20)
            if nb_swap_bounds_method != "0":
                self.change_setting_value(config, ["namelist:nlst_mpp","nb_swap_bounds_method"], "20")

    

        return config, self.reports


class vn110_t920(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #920 by Adam Voysey (adamvoysey)."""

    BEFORE_TAG = "vn11.0_t4023"
    AFTER_TAG = "vn11.0_t920"

    def upgrade(self, config, meta_config=None):
        """Add io_alltoall_readflds"""

        self.add_setting(config, ["namelist:io_control", "io_alltoall_readflds"], ".false.")

        return config, self.reports


class vn110_t3628(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3628 by Claudio Sanchez."""

    BEFORE_TAG = "vn11.0_t920"
    AFTER_TAG = "vn11.0_t3628"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_free_tracers", "l_calc_pv_full"], ".false.")
        return config, self.reports


class vn110_t3757(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3757 by annemccabe."""

    BEFORE_TAG = "vn11.0_t3628"
    AFTER_TAG = "vn11.0_t3757"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        l_rp2 = self.get_setting_value(config,
                     ["namelist:run_stochastic", "l_rp2"])
        if (l_rp2 == '.true.'):
            self.change_setting_value(config,
                 ["namelist:nlcfiles","rp2_seed"], 
                 "'$DATAW/seedfiles/$RUNID.seed'")
        return config, self.reports

class vn110_t3783(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3783 by paulselwood."""

    BEFORE_TAG = "vn11.0_t3757"
    AFTER_TAG = "vn11.0_t3783"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:tuning_segments", "ussp_seg_size"], "32")
        return config, self.reports


class vn110_t3502(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3502 by Richard Gilham."""

    BEFORE_TAG = "vn11.0_t3783"
    AFTER_TAG = "vn11.0_t3502"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM app configuration."""

        #Simply add the three new qsat switches
        self.add_setting(config, ["namelist:gen_phys_inputs","l_new_qsat"], ".false.")
        return config, self.reports


class vn110_t4011(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4011 by <Richard Barnes>."""

    BEFORE_TAG = "vn11.0_t3502"
    AFTER_TAG = "vn11.0_t4011"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.remove_setting(config,["namelist:run_precip", "l_autoc_3b"])
        self.remove_setting(config,["namelist:run_precip", "l_cry_agg_dep"])
        self.remove_setting(config,["namelist:run_precip", "l_hallett_mossop"])
        self.remove_setting(config,["namelist:run_precip", "l_it_melting"])
        self.remove_setting(config,["namelist:run_precip", "aic"])
        self.remove_setting(config,["namelist:run_precip", "bic"])
        self.remove_setting(config,["namelist:run_precip", "lsp_eic"])
        self.remove_setting(config,["namelist:run_precip", "lsp_fic"])
        return config, self.reports

class vn110_t3902(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3902 by Harry Shepherd."""

    BEFORE_TAG = "vn11.0_t4011"
    AFTER_TAG = "vn11.0_t3902"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config,
                         ['namelist:nlstcall', 'lstashdumptimer'],
                         value = ".false.")
        return config, self.reports

class vn110_t3680(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3680 by annemccabe."""

    BEFORE_TAG = "vn11.0_t3902"
    AFTER_TAG = "vn11.0_t3680"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Move run_stochastic in the shared namelist to be read after the jules namelists
        namelsts = self.get_setting_value(config, ["file:SHARED","source"],)
        if (namelsts):
            namelsts=namelsts.replace("namelist:run_stochastic namelist:gen_phys_inputs","namelist:gen_phys_inputs")
            namelsts=namelsts.replace("namelist:run_electric","namelist:run_stochastic namelist:run_electric")
            self.change_setting_value(config, ["file:SHARED","source"], namelsts)
        # Add new entries to the run_stochastic namelist
        dz0v_dh_rp_default = self.get_setting_value(config,["namelist:jules_pftparm","dz0v_dh_io"])
        self.add_setting(config, ["namelist:run_stochastic", "dz0v_dh_rp"], dz0v_dh_rp_default)
        self.add_setting(config, ["namelist:run_stochastic", "dz0v_dh_rp_min"], dz0v_dh_rp_default)
        self.add_setting(config, ["namelist:run_stochastic", "dz0v_dh_rp_max"], dz0v_dh_rp_default)
        npft = int(self.get_setting_value(config,["namelist:jules_surface_types","npft"]))
        self.add_setting(config, ["namelist:run_stochastic", "lai_mult_rp_min"], "{0:}*0.5".format(npft))
        self.add_setting(config, ["namelist:run_stochastic", "lai_mult_rp_max"], "{0:}*1.5".format(npft))
        z0hm_pft_rp_default = self.get_setting_value(config,["namelist:jules_pftparm","z0hm_pft_io"])
        self.add_setting(config, ["namelist:run_stochastic", "z0hm_pft_rp"], z0hm_pft_rp_default)
        self.add_setting(config, ["namelist:run_stochastic", "z0hm_pft_rp_min"], z0hm_pft_rp_default)
        self.add_setting(config, ["namelist:run_stochastic", "z0hm_pft_rp_max"], z0hm_pft_rp_default)
        return config, self.reports


class vn110_t3826(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3826 by Breo Gomez"""

    BEFORE_TAG = "vn11.0_t3680"
    AFTER_TAG = "vn11.0_t3826"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        return config, self.reports

class vn110_t3027(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3027 by Adrian Lock."""

    BEFORE_TAG = "vn11.0_t3826"
    AFTER_TAG = "vn11.0_t3027"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        idl_flux = self.get_setting_value(config, ["namelist:idealised", "idlsurffluxseaoption"])
        if idl_flux == "0":
           self.add_setting(config,["namelist:run_bl","flux_bc_opt"],"0")
        else:
           self.add_setting(config,["namelist:run_bl","flux_bc_opt"],"1")
        return config, self.reports

class vn110_t4012(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4012 by Hiroshi Kusabiraki."""

    BEFORE_TAG = "vn11.0_t3027"
    AFTER_TAG = "vn11.0_t4012"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        npft = int(self.get_setting_value(config, ["namelist:jules_surface_types", "npft"]))
        self.add_setting(config, ["namelist:jules_vegetation", "l_vegdrag_pft"], ','.join(['.false.']*npft))
        self.add_setting(config, ["namelist:jules_vegetation", "l_rsl_scalar"], ".false.")
        self.add_setting(config, ["namelist:jules_vegetation", "cd_leaf"], "0.25")
        self.add_setting(config, ["namelist:jules_vegetation", "c1_usuh"], "0.32")
        self.add_setting(config, ["namelist:jules_vegetation", "c2_usuh"], "0.264")
        self.add_setting(config, ["namelist:jules_vegetation", "c3_usuh"], "8.0")
        self.add_setting(config, ["namelist:jules_vegetation", "stanton_leaf"], "0.3")
        return config, self.reports

class vn110_t2113(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2113 by Maggie Hendry."""

    BEFORE_TAG = "vn11.0_t4012"
    AFTER_TAG = "vn11.0_t2113"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Move l_urban2t switch from jules_urban_switches to jules_surface. If
        # it isn't present set to false.
        l_urban2t = self.get_setting_value(config,
                                           ["namelist:jules_urban_switches",
                                            "l_urban2t"])
        if not l_urban2t:
            l_urban2t = ".false."
        self.add_setting(config, ["namelist:jules_surface", "l_urban2t"],
                         l_urban2t)
        self.remove_setting(config, ["namelist:jules_urban_switches",
                                     "l_urban2t"])
        # All urban surface types changed to compulsory so add with out of
        # range value when trigger ignored to alert the user to set it to
        # something sensible when enabled and "not used" value otherwise.
        if l_urban2t == ".true.":
            # urban_canyon and urban_roof has always been required by the UM so
            # these should already be present.
            self.add_setting(config, ["namelist:jules_surface_types", "urban"],
                             "0")
        else:
            self.add_setting(config, ["namelist:jules_surface_types",
                                      "urban_canyon"], "0")
            self.add_setting(config, ["namelist:jules_surface_types",
                                      "urban_roof"], "0")
            # Although 'urban' is usually present, it is not incorrect not to
            # have it, therefore set it to "not used" value.
            self.add_setting(config, ["namelist:jules_surface_types", "urban"],
                             "-1")

        # Allow urban namelists to be optional as they are now triggered. Need
        # to check first that they are there.
        if any([re.search(r'^file:ATMOSCNTL', s) for s in
                config.value.keys()]):
            source = self.get_setting_value(config, ["file:ATMOSCNTL",
                                                     "source"])
            source = source.replace("namelist:jules_urban2t_param",
                                    "(namelist:jules_urban2t_param)")
            self.change_setting_value(config, ["file:ATMOSCNTL","source"],
                                      source)
        else:
            source = self.get_setting_value(config, ["file:CNTLATM",
                                                     "source"])
            source = source.replace("namelist:jules_urban2t_param",
                                    "(namelist:jules_urban2t_param)")
            self.change_setting_value(config, ["file:CNTLATM","source"],
                                      source)

        source = self.get_setting_value(config, ["file:SHARED","source"])
        source = source.replace("namelist:jules_urban_switches",
                                "(namelist:jules_urban_switches)")
        self.change_setting_value(config, ["file:SHARED","source"], source)

        #Add new namelist item
        self.add_setting(config, ["namelist:jules_urban_switches",
                                  "l_urban_empirical"], ".false.")

        # anthrop_heat_scale doesn't appear to always be present, although it
        # should be
        self.add_setting(config, ["namelist:jules_urban2t_param",
                                  "anthrop_heat_scale"],
                         "1.0")

        return config, self.reports


class vn110_t3810(rose.upgrade.MacroUpgrade): 
    """Upgrade macro for ticket #3810 by Maff Glover."""
 
    BEFORE_TAG = "vn11.0_t2113"
    AFTER_TAG = "vn11.0_t3810"
 
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:tuning_segments", "l_autotune_segments"],    ".false.")
        self.add_setting(config, ["namelist:tuning_segments", "l_autotune_verbose"],     ".false.")
        self.add_setting(config, ["namelist:tuning_segments", "autotune_trial_radius"],       "16")
        self.add_setting(config, ["namelist:tuning_segments", "autotune_init_time_coeff"],  "1.05")
        self.add_setting(config, ["namelist:tuning_segments", "autotune_cooling_fraction"],"0.001")
        self.add_setting(config, ["namelist:tuning_segments", "autotune_cooling_steps"],      "50")
        self.add_setting(config, ["namelist:tuning_segments", "autotune_num_leaders"],         "3")
        self.add_setting(config, ["namelist:tuning_segments", "autotune_min_seg_size"],        "4")
        self.add_setting(config, ["namelist:tuning_segments", "autotune_max_seg_size"],      "128")
        self.add_setting(config, ["namelist:tuning_segments", "autotune_granularity"],         "4")
        return config, self.reports



class vn110_t4112(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn11.0_t3810"
    AFTER_TAG = "vn11.1"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "11.1", forced=True)
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMASTER"])
        if stashmaster_path:
           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",
                                     "vn11.1/ctldata",
                                      stashmaster_path)
           self.change_setting_value(config, ["env", "STASHMASTER"],
                                     stashmaster_path, forced=True)
        return config, self.reports


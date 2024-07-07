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



class vn109_t3556(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3556 by Steve Wardle"""

    BEFORE_TAG = "vn10.9"
    AFTER_TAG = "vn10.9_t3556"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Note - this macro does nothing; it exists because the given
        # ticket updates the STASHmaster in HEAD but the apps on the trunk
        # are still vn10.9, meaning they don't pick up the changes.  So this
        # macro exists purely to bump the app versions
        return config, self.reports


class vn109_t3484(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3484 by Sam Cusworth."""

    BEFORE_TAG = "vn10.9_t3556"
    AFTER_TAG = "vn10.9_t3484"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # No change, just need to increment metadata tag
        return config, self.reports


class vn109_t3460(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3460 by <Rachel Stratton>."""

    BEFORE_TAG = "vn10.9_t3484"
    AFTER_TAG = "vn10.9_t3460"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:idealised", "num_w_force_times"],
                         "0")
        self.add_setting(config, ["namelist:idealised", "num_w_force_heights"],
                         "0")
        self.add_setting(config, ["namelist:idealised", "w_force_height"],
                         "0.0")
        self.add_setting(config, ["namelist:idealised", "w_force_time"],
                         "0.0")
        self.add_setting(config, ["namelist:idealised", "w_force_data"],
                         "0.0")
        return config, self.reports

class vn109_t3685(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3685 by Adrian Lock."""
 
    BEFORE_TAG = "vn10.9_t3460"
    AFTER_TAG = "vn10.9_t3685"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_murk", "murk_source_scale"],"1.0")
        return config, self.reports


class vn109_t3011(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3011 by J. M. Edwards."""

    BEFORE_TAG = "vn10.9_t3685"
    AFTER_TAG = "vn10.9_t3011"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:temp_fixes", "l_fix_albsnow_ts"], ".false.")
        self.add_setting(config, ["namelist:temp_fixes", "l_fix_rcf_mlsnow_icefreemax"], ".false.")
        self.add_setting(config, ["namelist:jules_snow", "i_basal_melting_opt"], "0")
        return config, self.reports


class vn109_t3672(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3672 by Rachel Stratton."""

    BEFORE_TAG = "vn10.9_t3011"
    AFTER_TAG = "vn10.9_t3672"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:idealised", "num_surface_flux_times"],
                         "0")
        self.add_setting(config, ["namelist:idealised", "surface_flux_time"],
                         "0.0")
        self.add_setting(config, ["namelist:idealised", "sh_flux"],
                         "0.0")
        self.add_setting(config, ["namelist:idealised", "lh_flux"],
                         "0.0")
        return config, self.reports

class vn109_t3542(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3542 by paulearnshaw."""

    BEFORE_TAG = "vn10.9_t3672"
    AFTER_TAG = "vn10.9_t3542"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Set new variables
        self.add_setting(config,
                         ["namelist:jules_radiation", "fixed_sea_albedo"],
                         "0.31", state=config.STATE_SYST_IGNORED)
        self.add_setting(config,
                         ["namelist:jules_sea_seaice", "hcap_sea"],
                         "0.0")
        self.add_setting(config,
                         ["namelist:jules_sea_seaice", "l_use_dtstar_sea"],
                         ".false.")
        self.add_setting(config,
                         ["namelist:jules_sea_seaice", "z0hsea"],
                         "4.0e-5")

        # Move beta_evap from UM boundary layer to JULES sea/seaice
        # - If beta_evap not set then default to 1.0 as this is the value for a null effect
        beta_evap = self.get_setting_value(config, ["namelist:run_bl", "beta_evap"], no_ignore=True)
        if not beta_evap:
            beta_evap = "1.0"
        self.add_setting(config, ["namelist:jules_sea_seaice", "beta_evap"], beta_evap)
        self.remove_setting(config, ["namelist:run_bl", "beta_evap"])

        return config, self.reports


class vn100_t247(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #247 by Cyril Morcrette."""
 
    BEFORE_TAG = "vn10.9_t3542"
    AFTER_TAG = "vn10.9_t247"
 
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM app configuration."""
        self.add_setting(config, ["namelist:run_cloud", "i_pc2_checks_cld_frac_method"], "2")
        self.add_setting(config, ["namelist:run_cloud", "l_simplify_pc2_init_logic"], ".false.")
        return config, self.reports

class vn109_t3758(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3758 by Ian Boutle."""

    BEFORE_TAG = "vn10.9_t247"
    AFTER_TAG = "vn10.9_t3758"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config,["namelist:run_cloud","ice_fraction_method"],"1")
        l_eacf = self.get_setting_value(config,["namelist:run_cloud","l_eacf"])
        if l_eacf == ".true.":
            self.add_setting(config,["namelist:run_cloud","i_eacf"],"1")
        else:
            self.add_setting(config,["namelist:run_cloud","i_eacf"],"0")
        self.remove_setting(config,["namelist:run_cloud","l_eacf"])
        return config, self.reports


class vn109_t3681(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3681 by Martin Willett."""

    BEFORE_TAG = "vn10.9_t3758"
    AFTER_TAG = "vn10.9_t3681"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config,
                         ["namelist:temp_fixes", "l_fix_conv_diags_var"],
                         ".false.")
        return config, self.reports


class vn109_t3576(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3576 by Richard Barnes."""

    BEFORE_TAG = "vn10.9_t3681"
    AFTER_TAG = "vn10.9_t3576"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Note - this macro does nothing; it exists because the given
        # ticket makes metadata changes that need to be picked up
        # by the apps
        return config, self.reports


class vn109_t3419(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3419 by Maggie Hendry."""

    BEFORE_TAG = "vn10.9_t3576"
    AFTER_TAG = "vn10.9_t3419"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:recon_science",
                                  "l_regularize_landice_all_lice_pts"],
                         ".false.")
        self.add_setting(config, ["namelist:recon_science",
                                  "l_regularize_landice_all_soil_pts"],
                         ".false.")
        return config, self.reports


class vn109_t3781(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3781 by Edward Comyn-Platt."""

    BEFORE_TAG = "vn10.9_t3419"
    AFTER_TAG = "vn10.9_t3781"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Add settings
        self.add_setting(config, ["namelist:jules_soil_biogeochem", "t0_ch4"], "273.15")
        self.add_setting(config, ["namelist:jules_soil_biogeochem", "const_ch4_npp"], "9.99e-3")
        self.add_setting(config, ["namelist:jules_soil_biogeochem", "const_ch4_resps"], "4.36e-3")
        self.add_setting(config, ["namelist:jules_soil_biogeochem", "q10_ch4_cs"], "3.7")
        self.add_setting(config, ["namelist:jules_soil_biogeochem", "q10_ch4_npp"], "1.5")
        self.add_setting(config, ["namelist:jules_soil_biogeochem", "q10_ch4_resps"], "1.5")

        # If l_triffid = .true. then const_ch4_cs = 5.41e-10
        l_triffid = self.get_setting_value(config,
                                ["namelist:jules_vegetation", "l_triffid"])
        if (l_triffid==".true."):
            self.add_setting(config,
                 ["namelist:jules_soil_biogeochem", "const_ch4_cs"], "5.41e-10")
        else:
            self.add_setting(config,
                 ["namelist:jules_soil_biogeochem", "const_ch4_cs"], "5.41e-12")

        return config, self.reports


class vn109_t3786(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3786 by sarahchadburn."""

    BEFORE_TAG = "vn10.9_t3781"
    AFTER_TAG = "vn10.9_t3786"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:jules_soil", "l_holdwater"], ".false.")
        return config, self.reports


class vn109_t3622(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3622 by Fiona O'Connor / Luke Abraham."""

    BEFORE_TAG = "vn10.9_t3786"
    AFTER_TAG = "vn10.9_t3622"

    def upgrade(self, config, meta_config=None):
        """
        Introduce logical to interpolate linearly in LOG(p)
        for the vertical redistribution of Lightning NOx.
        """

        self.add_setting(config, ["namelist:run_ukca",
                                  "l_ukca_linox_logp"], ".false.")

        return config, self.reports


class vn109_t2545(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2545 by paulearnshaw."""

    BEFORE_TAG = "vn10.9_t3622"
    AFTER_TAG = "vn10.9_t2545"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config,
                         ["namelist:temp_fixes", "l_fix_lsp_incs_to_spt"],
                         ".false.")
        return config, self.reports

class vn109_t2070(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2070 by paulearnshaw."""

    BEFORE_TAG = "vn10.9_t2545"
    AFTER_TAG = "vn10.9_t2070"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here

        self.add_setting(config, ["namelist:temp_fixes",
                                  "l_fix_ec_gen_hgt"], ".false.")

        return config, self.reports


class vn109_t3667(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3667 by Mohamed Zerroukat."""

    BEFORE_TAG = "vn10.9_t2070"
    AFTER_TAG = "vn10.9_t3667"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        l_conservation_moist_zlf = self.get_setting_value(config, ["namelist:run_dyn","l_conservation_moist_zlf"])
        if l_conservation_moist_zlf == ".false.":
           self.change_setting_value(config, ["namelist:run_dyn","zlf_conservation_moist_option"],"0")
        self.remove_setting(config,["namelist:run_dyn", "l_conservation_moist_zlf"])
        self.add_setting(config,["namelist:run_dyn", "zlf_conservation_tracers_option"],"0")
        self.add_setting(config,["namelist:run_dyn", "zlf_conservation_theta_option"],"0")
        return config, self.reports

class vn110_t3821(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.9_t3667"
    AFTER_TAG = "vn11.0"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "11.0", forced=True)
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMASTER"])
        if stashmaster_path:
           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",
                                     "vn11.0/ctldata",
                                      stashmaster_path)
           self.change_setting_value(config, ["env", "STASHMASTER"],
                                     stashmaster_path, forced=True)
        return config, self.reports


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


class vn107_t2843(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2843 by andymalcolm."""

    BEFORE_TAG = "vn10.7"
    AFTER_TAG = "vn10.7_t2843"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        
        ios_offset_num = self.get_setting_value(
               config, ["namelist:ioscntl", "ios_offset"])
        num=int(ios_offset_num)
        if (num < 0):
            self.change_setting_value(
                     config, ["namelist:ioscntl", "ios_offset"], "0")
        return config, self.reports


class vn107_t2745(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2745 by Chris Smith."""

    BEFORE_TAG = "vn10.7_t2843"
    AFTER_TAG = "vn10.7_t2745"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Delete New Dynamics settings
        self.remove_setting(config, ["namelist:idealise", "brunt_vaisala"])
        self.remove_setting(config, ["namelist:idealise", "hf"])
        self.remove_setting(config, ["namelist:idealise", "l_baroclinic"])
        self.remove_setting(config, ["namelist:idealise", "l_force"])
        l_trivial_trigs = self.get_setting_value(config,
                               ["namelist:idealise", "l_trivial_trigs"])
        if l_trivial_trigs == ".true.":
           self.change_setting_value(config,
                ["namelist:idealise", "l_cartesian"], ".true.")
        self.remove_setting(config, ["namelist:idealise", "l_trivial_trigs"])
        self.remove_setting(config, ["namelist:idealise", "q1"])
        self.remove_setting(config, ["namelist:idealise", "u_in"])
        self.remove_setting(config, ["namelist:idealise", "v_in"])
        self.remove_setting(config, ["namelist:idealise", "uvprofile_number"])
        self.remove_setting(config, ["namelist:idealise", "num_tforce_times"])
        self.remove_setting(config, ["namelist:idealise", "num_qforce_times"])
        self.remove_setting(config, ["namelist:idealise", "num_uvforce_times"])
        self.remove_setting(config, ["namelist:idealise", "num_pforce_times"])
        self.remove_setting(config, 
             ["namelist:idealise", "tforce_time_interval"])
        self.remove_setting(config, 
             ["namelist:idealise", "qforce_time_interval"])
        self.remove_setting(config, 
             ["namelist:idealise", "uvforce_time_interval"])
        self.remove_setting(config, 
             ["namelist:idealise", "pforce_time_interval"])
        self.remove_setting(config, ["namelist:idealise", "tforce_data"])
        self.remove_setting(config, ["namelist:idealise", "qforce_data"])
        self.remove_setting(config, ["namelist:idealise", "uforce_data"])
        self.remove_setting(config, ["namelist:idealise", "vforce_data"])
        self.remove_setting(config, ["namelist:idealise", "p_surface_data"])
        return config, self.reports


class vn107_t2861(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2861 by andymalcolm."""

    BEFORE_TAG = "vn10.7_t2745"
    AFTER_TAG = "vn10.7_t2861"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config,
                         ["namelist:nlstcall", "ltimers_user"],
                         ".false.")
        return config, self.reports

class vn107_t2730(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2730 by Ian Boutle."""

    BEFORE_TAG = "vn10.7_t2861"
    AFTER_TAG = "vn10.7_t2730"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        l_planet_grey_surface = self.get_setting_value(config, ["namelist:planet_constants", "l_planet_grey_surface"])
        if l_planet_grey_surface == ".true.":
            t_surf = self.get_setting_value(config, ["namelist:planet_constants", "planet_t_surface"])
            self.add_setting(config, ["namelist:items(12345)", "user_prog_rconst"],t_surf)
            self.add_setting(config, ["namelist:items(12345)", "source"],"6")
            self.add_setting(config, ["namelist:items(12345)", "stash_req"],"24")
            self.add_setting(config, ["namelist:items(12345)", "ancilfilename"],"''")
            self.add_setting(config, ["namelist:items(12345)", "domain"],"1")
            self.add_setting(config, ["namelist:items(12345)", "update_anc"],".false.")
            self.add_setting(config, ["namelist:items(12345)", "interval"],"0",state=config.STATE_SYST_IGNORED)
            self.add_setting(config, ["namelist:items(12345)", "netcdf_varname"],"''",state=config.STATE_SYST_IGNORED)
            self.add_setting(config, ["namelist:items(12345)", "period"],"1",state=config.STATE_SYST_IGNORED)
            self.add_setting(config, ["namelist:items(12345)", "user_prog_ancil_stash_req"],"0",state=config.STATE_SYST_IGNORED)            
        problem_number = self.get_setting_value(config, ["namelist:run_dyntest", "problem_number"])
        if problem_number == "5":
            t_surf = self.get_setting_value(config, ["namelist:idealise", "t_surface"])
            self.add_setting(config, ["namelist:items(12345)", "user_prog_rconst"],t_surf)
            self.add_setting(config, ["namelist:items(12345)", "source"],"6")
            self.add_setting(config, ["namelist:items(12345)", "stash_req"],"24")
            self.add_setting(config, ["namelist:items(12345)", "ancilfilename"],"''")
            self.add_setting(config, ["namelist:items(12345)", "domain"],"1")
            self.add_setting(config, ["namelist:items(12345)", "update_anc"],".false.")
            self.add_setting(config, ["namelist:items(12345)", "interval"],"0",state=config.STATE_SYST_IGNORED)
            self.add_setting(config, ["namelist:items(12345)", "netcdf_varname"],"''",state=config.STATE_SYST_IGNORED)
            self.add_setting(config, ["namelist:items(12345)", "period"],"1",state=config.STATE_SYST_IGNORED)
            self.add_setting(config, ["namelist:items(12345)", "user_prog_ancil_stash_req"],"0",state=config.STATE_SYST_IGNORED)            
        self.remove_setting(config, ["namelist:idealise", "theta_surface"])
        self.remove_setting(config, ["namelist:idealise", "newtonian_timescale"])
        self.remove_setting(config, ["namelist:planet_constants", "planet_t_surface"])

        # We need to rehash the stash indices using TidyStashTransform().
        # First we add the macro to the front of the path.
        sys.path.insert(1,
                        (os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "vn10.7", "lib", "python", "macros")))
        import stash_indices
        # We must reload the module
        reload(stash_indices)
        # Then we quickly remove stash_indices from sys.path to avoid namespace
        # conflicts.
        sys.path.pop(1)
        # Finally we rehash stash indices.
        config, reports = stash_indices.TidyStashTransform().transform(config)
        self.reports += reports

        return config, self.reports

class vn107_t1448(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1448 by Mohit Dalvi."""

    BEFORE_TAG = "vn10.7_t2730"
    AFTER_TAG = "vn10.7_t1448"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Remove items 38-40x from UKCA coupling stash package
        #   if requested with 'UPUKCA' usage profile or 
        #   'UKCA coupling macro' as package name
        itm_del = ['401', '402', '403', '404', '405', '406', '407']
        sec_del = '38'
        to_del = []  # store items (as python objects) to delete
        for obj in config.get_value():
       
           # Check if STASH request
           if re.search(r'namelist:umstash_streq', obj):
             
              # get the section, item, Usage profile and package name 
              item = self.get_setting_value(config, [ obj, 'item'])
              section = self.get_setting_value(config, [ obj, 'isec'])
              uprof = self.get_setting_value(config, [ obj, 'use_name'])
              pkg = self.get_setting_value(config, [ obj, 'package'])

              # check if the stash code, usage profile or package names
              #  match
              if section == sec_del and item in itm_del:
                 if re.search(r'UPUKCA', uprof) or \
                    re.search(r'ukca coupling macro', pkg.lower()):
                    # append to the list to be deleted
                    to_del.append(obj)
       
        # loop over the selected streq items and actually delete them
        for delme in to_del:
           self.remove_setting(config, [delme],info='UKCA coupling item'+
               ' no longer required')

        return config, self.reports


class vn107_t2988(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2988 by Glenn Greed."""

    BEFORE_TAG = "vn10.7_t1448"
    AFTER_TAG = "vn10.7_t2988"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        for section in config.value.keys():
            node = config.get([section])
            if (isinstance(node.value, dict) and
                    re.match("namelist:nlstcall_pp\(.*?\)$", section)):
                packing = self.get_setting_value(config, [section, "packing"])
                if (packing == "6"):
                    raise UpgradeError("WARNING! Your app requests Simple GRIB packing (6)\n"
                                       "This UM option is not supported and has been retired \n" 
                                       "Please revisit your app and consider an alternative packing option (packing)")

        ppxm_setting = self.get_setting_value(config, ["namelist:nlstcgen", "ppxm"])
        if ( ppxm_setting == "6" ):
            raise UpgradeError("WARNING! Your app requests Simple GRIB packing (6)\n"
                               "This UM option is not supported and has been retired \n" 
                               "Please revisit your app and consider an alternative packing option for means (ppxm)")

        return config, self.reports


class vn107_1908(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1908 by Ben Shipway."""

    BEFORE_TAG = "vn10.7_t2988"
    AFTER_TAG = "vn10.7_t1908"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        value='.false.'
        p1=self.get_setting_value(config, ["namelist:gen_phys_inputs", "l_mr_physics1"])
        p2=self.get_setting_value(config, ["namelist:gen_phys_inputs", "l_mr_physics2"])
        if (p1[:2]=='.t' or p1[:2]=='.T') and (p2[:2]=='.t' or p2[:2]=='.T'): value='.true.'
        if not p1==p2:
            warning='''
                    Warning! The option to mix use of mixing ratios and specific quantities
                    in physics1 and physics2 is no longer supported.  You job will default
                    to using specific quantities.
                    '''
            self.add_report(info=warning, is_warning=True)

        self.add_setting(config, ["namelist:gen_phys_inputs", "l_mr_physics"], value)
        self.remove_setting(config,["namelist:gen_phys_inputs", "l_mr_physics1"])
        self.remove_setting(config,["namelist:gen_phys_inputs", "l_mr_physics2"])
        return config, self.reports


class vn107_t2817(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2817 by Huw Lewis."""

    BEFORE_TAG = "vn10.7_t1908"
    AFTER_TAG = "vn10.7_t2817"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Just add the new compulsory namelist item with its default value
        self.add_setting(config, ["namelist:jules_hydrology", "l_spdmvar"], ".false.")
        self.add_setting(config, ["namelist:jules_hydrology", "slope_pdm_max"], "6.0")
        self.add_setting(config, ["namelist:jules_hydrology", "s_pdm"], "0.0")
        return config, self.reports

class vn107_t2772(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2772 by Rachel Stratton."""

    BEFORE_TAG = "vn10.7_t2817"
    AFTER_TAG = "vn10.7_t2772"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        # Only altered triggering of variables
        return config, self.reports


class vn107_t2819(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2819 by Luke Abraham."""

    BEFORE_TAG = "vn10.7_t2772"
    AFTER_TAG = "vn10.7_t2819"

    def upgrade(self, config, meta_config=None):
        """
        Add-in options for using a quasi-Newton method
        in the UKCA chemical solver.
        This can only be used by schemes which use the 
        ASAD chemical package with the Newton-Raphson 
        solver, i.e. i_ukca_chem=50,51,52,54.
        Default values:
         l_ukca_quasinewton       = .false.
         i_ukca_quasinewton_start = 2 (recommended value)
         i_ukca_quasinewton_end   = 3 (recommended value)
        Given that l_ukca_quasinewton = .false. the 
        values for the other variables do not matter
        here, so recommended ones are provided.
        The range for i_ukca_quasinewton_start and 
        i_ukca_quasinewton_end can only be between 2 
        and 50, although only 2 & 3 have been tested in
        detail (see ticket #2819).
        """

        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_quasinewton"],".false.")
        self.add_setting(config,["namelist:run_ukca",
                                 "i_ukca_quasinewton_start"],"2")
        self.add_setting(config,["namelist:run_ukca",
                                 "i_ukca_quasinewton_end"],"3")

        return config, self.reports


class vn107_t3080(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3080 by jamierae."""

    BEFORE_TAG = "vn10.7_t2819"
    AFTER_TAG = "vn10.7_t3080"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:temp_fixes", "l_fix_alb_ice_thick"], ".false.")
        return config, self.reports


class vn107_t3082(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3082 by Robin Smith."""

    BEFORE_TAG = "vn10.7_t3080"
    AFTER_TAG = "vn10.7_t3082"

    def upgrade(self, config, meta_config=None):

        """Make fsmc_p0 a PFT parameter available to the UM namelist"""

        npft = int(self.get_setting_value(config, ["namelist:jules_surface_types", "npft"]))
        self.add_setting(config, ["namelist:jules_pftparm", "fsmc_p0_io"],','.join(['0.0']*npft))

        return config, self.reports

class vn107_t3046(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3046 by J. M. Edwards"""
   
    BEFORE_TAG = "vn10.7_t3082"
    AFTER_TAG = "vn10.7_t3046"
   
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM app configuration."""
        self.add_setting(config, ["namelist:jules_sea_seaice", "i_high_wind_drag"], "0")
        self.add_setting(config, ["namelist:jules_sea_seaice", "cd_limit_sea"], "0.0024")
        return config, self.reports


class vn107_t2524(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2524 by Richard Hill."""

    BEFORE_TAG = "vn10.7_t3046"
    AFTER_TAG = "vn10.7_t2524"

    def upgrade(self, config, meta_config=None):
        """Remove all settings associated with pseudo-parallel OASIS3."""

        # Remove redundant l_couple_master namelist entry now that
	# OASI3 is no longer supported.
        self.remove_setting(config, ["namelist:coupling_control", "l_couple_master"])

        # Remove OASIS3 launch options environment variable, if set.
	prismlaunch = self.get_setting_value(config, ["env", "LAUNCH_MPI_PRISM"])
	if prismlaunch is not None:
            self.remove_setting(config, ["env", "LAUNCH_MPI_PRISM"])

	# Check if OASIS3 coupler is active. If it is then raise an exception.
	# This is chosen in preference to converting to oasis3-mct
	# since that may lead to the job being run with invalid namcouple
	# control files and resource requests.
	# Setting to none would similarly eventually lead to failure
	# some way down the line, so we nip things in the bud at the
	# earliest possible stage.
	coupler = self.get_setting_value(config, ["env", "COUPLER"])
	if coupler == "oasis3":
            raise UpgradeError('This suite uses the OASIS3 coupler which is no longer supported.\n' +
                               'Please edit your suite to select the OASIS3-MCT coupler.\n' +
                               'Please ensure you employ an appropriate namcouple control file.')

	return config, self.reports


class vn107_t3005(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3005 by Adrian Lock."""

    BEFORE_TAG = "vn10.7_t2524"
    AFTER_TAG = "vn10.7_t3005"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, 
                         ["namelist:temp_fixes", "l_fix_zh"],
                         ".false.")
        return config, self.reports

class vn107_t3076(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3076 by Mike Whitall."""

    BEFORE_TAG = "vn10.7_t3005"
    AFTER_TAG = "vn10.7_t3076"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, 
                         ["namelist:temp_fixes", "l_fix_ccb_cct"],
                         ".false.")
        return config, self.reports

class vn107_t205(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #205 by Samantha Smith."""

    BEFORE_TAG = "vn10.7_t3076"
    AFTER_TAG = "vn10.7_t205"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_precip", "l_orograin"], ".false.")
        self.add_setting(config, ["namelist:run_precip", "l_orograin_block"], ".false.")
        self.add_setting(config, ["namelist:run_precip", "l_orogrime"], ".false.")
        self.add_setting(config, ["namelist:run_precip", "fcrit"], "1.0")
        self.add_setting(config, ["namelist:run_precip", "nsigmasf"], "2.82843")
        self.add_setting(config, ["namelist:run_precip", "nscalesf"], "1.0")
        return config, self.reports


class vn107_t3104(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3104 by David Walters and Thomas Allen."""

    BEFORE_TAG = "vn10.7_t205"
    AFTER_TAG = "vn10.7_t3104"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Set maximum iteration count to previously hardwired value of 999
        # which was set within bicgstab - this now uses the related
        # ND variable.
        self.change_setting_value(
            config, ["namelist:run_dyn", "gcr_max_iterations"], "999")
        return config, self.reports


class vn107_t2870(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2870 by Rachel Stratton."""

    BEFORE_TAG = "vn10.7_t3104"
    AFTER_TAG = "vn10.7_t2870"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:idealise","l_ideal_2d"], ".false.")
        self.add_setting(config, ["namelist:idealise",
	                          "local_heat_option"], "0")
        self.add_setting(config, ["namelist:idealise",
	                          "local_heat_xoffset"], "0.0")
        self.add_setting(config, ["namelist:idealise",
	                          "local_heat_yoffset"], "0.0")
        self.add_setting(config, ["namelist:idealise",
	                          "local_heat_amp"], "0.0")
        self.add_setting(config, ["namelist:idealise",
	                          "local_heat_sigma"], "0.0")
        self.add_setting(config, ["namelist:idealise",
	                          "local_heat_base"], "0.0")
        self.add_setting(config, ["namelist:idealise",
	                          "local_heat_top"], "0.0")
        self.add_setting(config, ["namelist:idealise",
	                          "local_heat_period"], "0.0")
        return config, self.reports


class vn107_t2878(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2878 by Claudio Sanchez"""

    BEFORE_TAG = "vn10.7_t2870"
    AFTER_TAG = "vn10.7_t2878"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_free_tracers", "l_diab_tracer"], ".false.")
        self.add_setting(config, ["namelist:run_free_tracers", "l_diab_tr_rad"], ".false.")
        self.add_setting(config, ["namelist:run_free_tracers", "l_diab_tr_bl"], ".false.")
        return config, self.reports


class vn107_t773(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #773 by Jonathan Wilkinson."""

    BEFORE_TAG = "vn10.7_t2878"
    AFTER_TAG = "vn10.7_t773"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_precip", "casim_moments_choice"], "1")
        self.add_setting(config, ["namelist:run_precip", "l_casim"], ".false.")
        return config, self.reports

class vn107_t2995(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2995 by sarahchadburn."""

    BEFORE_TAG = "vn10.7_t773"
    AFTER_TAG = "vn10.7_t2995"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        
        l_npp = self.get_setting_value(config, ["namelist:jules_hydrology", "l_wetland_ch4_npp"])
        if l_npp == ".true.":
            self.add_setting(config, ["namelist:jules_soil_biogeochem", "ch4_substrate"], "2")
        else:
            self.add_setting(config, ["namelist:jules_soil_biogeochem", "ch4_substrate"], "1")

        self.add_setting(config, ["namelist:jules_soil_biogeochem", "l_ch4_interactive"], ".false.")
        self.add_setting(config, ["namelist:jules_soil_biogeochem", "l_ch4_tlayered"], ".false.")
        self.remove_setting(config, ["namelist:jules_hydrology", "l_wetland_ch4_npp"])

        return config, self.reports


class vn107_t2087(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2087 by Maggie Hendry."""

    BEFORE_TAG = "vn10.7_t2995"
    AFTER_TAG = "vn10.7_t2087"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add elev_rock to the namelist with a value of 0 or -1"""
        l_elev_land_ice = self.get_setting_value(config, ["namelist:jules_surface", "l_elev_land_ice"])
        if l_elev_land_ice == ".true.":
            # Add elev_rock with -1 so existing configurations pass fail-if
            self.add_setting(config, ["namelist:jules_surface_types", "elev_rock"], "-1" )
        else:
            # Add an out of range value otherwise to alert user to change the
            # value to something sensible
            self.add_setting(config, ["namelist:jules_surface_types", "elev_rock"], "0" )
        return config, self.reports

class vn107_t3048(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3048 by James Manners."""

    BEFORE_TAG = "vn10.7_t2087"
    AFTER_TAG = "vn10.7_t3048"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add namelist items for previously hardwired values
        self.add_setting(config,["namelist:run_radiation","aparam"],"0.07")
        self.add_setting(config,["namelist:run_radiation","bparam"],"-0.14")
        return config, self.reports


class vn107_t288(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #288 by Alan J Hewitt"""

    BEFORE_TAG = "vn10.7_t3048"
    AFTER_TAG = "vn10.7_t288"

    def upgrade(self, config, meta_config=None):
        """Add l_glomap_clim_radaer & l_glomap_clim_radaer_sustrat &
        gclmaclw & gclmacsw & gclmanlw & gclmansw & gclmcrlw & gclmcrsw &
        gclmprec"""
        # Input your macro commands here
        
        # l_glomap_clim_radaer should default to false
        self.add_setting(config, ["namelist:run_glomap_aeroclim",
                                  "l_glomap_clim_radaer"],
                         value=".false.", forced=True)
        
        # l_glomap_clim_radaer_sustrat should default to false
        self.add_setting(config, ["namelist:run_glomap_aeroclim",
                                  "l_glomap_clim_radaer_sustrat"],
                         value=".false.", forced=True)
        
        self.add_setting(config,['namelist:run_glomap_aeroclim','gclmaclw'],
                         value = "'unset'")
        
        self.add_setting(config,['namelist:run_glomap_aeroclim','gclmacsw'],
                         value = "'unset'")
        
        self.add_setting(config,['namelist:run_glomap_aeroclim','gclmanlw'],
                         value = "'unset'")
        
        self.add_setting(config,['namelist:run_glomap_aeroclim','gclmansw'],
                         value = "'unset'")
        
        self.add_setting(config,['namelist:run_glomap_aeroclim','gclmcrlw'],
                         value = "'unset'")
        
        self.add_setting(config,['namelist:run_glomap_aeroclim','gclmcrsw'],
                         value = "'unset'")
        
        self.add_setting(config,['namelist:run_glomap_aeroclim','gclmprec'],
                         value = "'unset'")
        
        return config, self.reports


class vn107_t2972(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2972 by harryshepherd."""

    BEFORE_TAG = "vn10.7_t288"
    AFTER_TAG = "vn10.7_t2972"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # The following run scripts have been removed in this ticket
        invalid_cmd = ["GC3-coupled", "um-coupled", "um-coupled-hybrid"]
        
        default_run_cmd = self.get_setting_value(
            config, ["command", "default"])
        if default_run_cmd in invalid_cmd:
            report_msg = "Unable to upgrade suite. This run command has been" \
                " removed from the UM at version 10.8. Please use the" \
                " coupled drivers in the MOCI repository to run a coupled" \
                " configuration"
            raise UpgradeError(report_msg)
        return config, self.reports


class vn107_t2860(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2860 by Anne McCabe."""

    BEFORE_TAG = "vn10.7_t2972"
    AFTER_TAG = "vn10.7_t2860"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Move run_stochastic from ATMOSCNTL to SHARED
        namelsts = self.get_setting_value(config, ["file:SHARED","source"],)
        if (namelsts):
            namelsts=namelsts.replace("namelist:gen_phys_inputs", "namelist:run_stochastic namelist:gen_phys_inputs")
            self.change_setting_value(config, ["file:SHARED","source"], 
                                                             namelsts)
        namelsts = self.get_setting_value(config, ["file:ATMOSCNTL","source"],)
        if (namelsts):
            namelsts=namelsts.replace("namelist:run_stochastic","")
            self.change_setting_value(config, ["file:ATMOSCNTL","source"], 
                                                             namelsts)
        return config, self.reports

class vn108_t3162(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.7_t2860"
    AFTER_TAG = "vn10.8"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "10.8", forced=True)
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMASTER"])
        if stashmaster_path:
           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",
                                     "vn10.8/ctldata",
                                      stashmaster_path)
           self.change_setting_value(config, ["env", "STASHMASTER"],
                                     stashmaster_path, forced=True)
        return config, self.reports


import os
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



class vn108_t3170(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3170 by Joe Mancell/Saeed Sadri."""

    BEFORE_TAG = "vn10.8"
    AFTER_TAG = "vn10.8_t3170"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Blank upgrade macro to force apps to use new metadata
        # associated with ticket #3170 which is only present in HEAD.
        # This is needed as the head of trunk
        # rose-stem apps currently point at vn10.8 rather than an
        # intermeditary version tag so will not pick up HEAD metadata.
        # This is an edge case scenario affecting only the head of trunk
        # rose stem suite.
        pass
        return config, self.reports

class vn108_t3354(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3354 by Ian Boutle."""

    BEFORE_TAG = "vn10.8_t3170"
    AFTER_TAG = "vn10.8_t3354"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.remove_setting(config, ["namelist:run_bl", "variable_ric"])
        self.remove_setting(config, ["namelist:run_bl", "entr_enhance_by_cu"])
        self.remove_setting(config, ["namelist:run_bl", "dec_thres_cloud"])
        self.remove_setting(config, ["namelist:run_bl", "my_prod_adj_fact"])
        alpha_cd=self.get_setting_value(config, ["namelist:run_bl", "alpha_cd"])
        alpha_cd = alpha_cd.split(',')
        alpha_cd = alpha_cd[:2]
        alpha_cd = ','.join(alpha_cd)
        self.add_setting(config, ["namelist:run_bl", "alpha_cd_in"],alpha_cd)
        self.remove_setting(config, ["namelist:run_bl", "alpha_cd"])
        return config, self.reports


class vn108_t1444(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1444 by Chris Smith."""

    BEFORE_TAG = "vn10.8_t3354"
    AFTER_TAG = "vn10.8_t1444"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Remove unused settings
        self.remove_setting(config, ["namelist:run_dyn", "l_mix_ratio"])
        self.remove_setting(config, ["namelist:idealise", "eccentricity"])
        # Move intrand_seed from run_dyn to run_dyntest
        intrand_seed = self.get_setting_value(config, 
                            ["namelist:run_dyn", "intrand_seed"])
        self.remove_setting(config, ["namelist:run_dyn", "intrand_seed"])
        self.add_setting(config, ["namelist:run_dyntest", "intrand_seed"],
                         intrand_seed)
        # alpha_relax_type = 0 is set in some reconfiguration apps
        # but this is an invalid setting (though not actually used
        # in reconfiguration).
        alpha_relax_type = self.get_setting_value(config, 
                                ["namelist:run_dyn", "alpha_relax_type"])
        if (alpha_relax_type == "0"):
            self.change_setting_value(config, 
                 ["namelist:run_dyn", "alpha_relax_type"], "1")

        return config, self.reports


class vn108_t2029(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2029 by Kirsty Hanley and Mike Whitall."""

    BEFORE_TAG = "vn10.8_t1444"
    AFTER_TAG = "vn10.8_t2029"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_diffusion", "l_leonard_term"],
                         ".false.")
        self.add_setting(config, ["namelist:run_diffusion", "leonard_kl"],
                         "1.0")
        return config, self.reports


class vn108_t739(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #739 by Jeff Cole."""

    BEFORE_TAG = "vn10.8_t2029"
    AFTER_TAG = "vn10.8_t739"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        atmoscntl_source = self.get_setting_value(config, ["file:ATMOSCNTL","source"])
        if atmoscntl_source:
            atmoscntl_source = re.sub(r'namelist:nlstcall_pp\(:\)', '', atmoscntl_source)
            atmoscntl_source += ' (namelist:nlstcall_pp(:)) (namelist:nlstcall_nc(:)) namelist:nlstcall_nc_options'
            self.change_setting_value(config, ["file:ATMOSCNTL","source"], atmoscntl_source)

        # Add default values for nlstcall_nc namelists
        self.add_setting(config, [ "namelist:nlstcall_nc(nc0)", "file_id" ], value = "'nc0'")
        self.add_setting(config, [ "namelist:nlstcall_nc(nc0)", "filename_base" ], value = "'$DATAM/${RUNID}a_pa%N.nc'")
        self.add_setting(config, [ "namelist:nlstcall_nc(nc0)", "reinit_unit" ], value = "2")
        self.add_setting(config, [ "namelist:nlstcall_nc(nc0)", "reinit_start" ], value = "0")
        self.add_setting(config, [ "namelist:nlstcall_nc(nc0)", "reinit_end" ], value = "-1")
        self.add_setting(config, [ "namelist:nlstcall_nc(nc0)", "reinit_step" ], value = "1")
        self.add_setting(config, [ "namelist:nlstcall_nc(nc0)", "l_compress" ], value = ".false.")
        self.add_setting(config, [ "namelist:nlstcall_nc(nc0)", "packing" ], value = "0")
        self.add_setting(config, [ "namelist:nlstcall_nc(nc0)", "nccomp_level" ], value = "1")
        self.add_setting(config, [ "namelist:nlstcall_nc(nc0)", "l_shuffle" ], value = ".true.")
        # Add default values for nlstcall_nc_options namelist
        self.add_setting(config, [ "namelist:nlstcall_nc_options", "l_netcdf" ], value = ".false.")
        self.add_setting(config, [ "namelist:nlstcall_nc_options", "ncformat" ], value = "4")
        self.add_setting(config, [ "namelist:nlstcall_nc_options", "ncdim_prec" ], value = "2")
        self.add_setting(config, [ "namelist:nlstcall_nc_options", "ncvar_prec" ], value = "2")
        self.add_setting(config, [ "namelist:nlstcall_nc_options", "l_checksum" ], value = ".false.")

        return config, self.reports


class vn108_t142(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #142 by Mohit Dalvi."""

    BEFORE_TAG = "vn10.8_t739"
    AFTER_TAG = "vn10.8_t142"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Remove option for 'new' emissions processing system
        # in UKCA as there is only one available now.
        # Also removes specification of level to inject SO2 emissions
        # from an ancil file.

        self.remove_setting(config,
                         ["namelist:run_ukca","l_ukca_new_emiss"])
        self.remove_setting(config,
                         ["namelist:run_ukca","i_so2_hi_level"])

        return config, self.reports

class vn108_t3209(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3209 by Luke Abraham."""

    BEFORE_TAG = "vn10.8_t142"
    AFTER_TAG = "vn10.8_t3209"

    def upgrade(self, config, meta_config=None):
        """
        Add column-call functionality to UKCA.
        This can only be used by schemes which use the
        ASAD chemical package with the Newton-Raphson
        solver, i.e. i_ukca_chem=50,51,52,54.
        Default value:
          l_ukca_asad_columns = .false.
        """
 
        self.add_setting(config,["namelist:run_ukca",
                                 "l_ukca_asad_columns"],".false.")

        return config, self.reports


class vn108_t3426(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3426 by Andy Wiltshire ."""

    BEFORE_TAG = "vn10.8_t3209"
    AFTER_TAG = "vn10.8_t3426"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:jules_vegetation",
                                  "l_trif_init_accum"],
                         ".true.")
        return config, self.reports


class vn108_t2405(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2405 by colinjohnson."""

    BEFORE_TAG = "vn10.8_t3426"
    AFTER_TAG = "vn10.8_t2405"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config,
                         ["namelist:temp_fixes", "l_fix_nacl_density"],
                         ".false.")
        return config, self.reports


class vn108_t3454(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3454 by Stuart Whitehouse."""

    BEFORE_TAG = "vn10.8_t2405"
    AFTER_TAG = "vn10.8_t3454"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        for obj in config.get_value():
            # is it an streq namelist?
            if re.search(r'namelist:umstash_streq', obj):
                item = self.get_setting_value(config,[obj,'item'])
                section = self.get_setting_value(config,[obj,'isec'])
                if section == '20':
                    if item =='64' or item == '064':
                        self.change_setting_value(config,[obj,'item'], '84')
                    elif item == '65' or item == '065':
                        self.change_setting_value(config,[obj,'item'], '85')
                    elif item == '66' or item == '066':
                        self.change_setting_value(config,[obj,'item'], '86')
                    elif item == '67' or item == '067':
                        self.change_setting_value(config,[obj,'item'], '87')
                        
        # Hash sorting:
        # Load the path of the hashing macro.
        sys.path.insert(1,
                        (os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "vn10.8", "lib", "python", "macros")))
        import stash_indices
        # Reload the module
        reload(stash_indices)
        # Remove stash_indices from sys.path to avoid namespace conflicts.
        sys.path.pop(1)
        # Finally we rehash the stash indices.
        config, reports = stash_indices.TidyStashTransform().transform(config)
        self.reports += reports

                        
        return config, self.reports

class vn108_t3400(rose.upgrade.MacroUpgrade):
	
    """Upgrade macro for ticket #3400 by <Mohit Dalvi>."""

    BEFORE_TAG = "vn10.8_t3454"
    AFTER_TAG = "vn10.8_t3400"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Introduce a logical to limit the formation
        # of Nitric Acid Trihydrate below a specified height

        self.add_setting(config,["namelist:run_ukca",
                        "l_ukca_limit_nat"],".false.")

        return config, self.reports


class vn108_t2587(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2587 by fraserdennison."""

    BEFORE_TAG = "vn10.8_t3400"
    AFTER_TAG = "vn10.8_t2587"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, [ "namelist:run_ukca", "i_ukca_solcyc"], "0")
        self.add_setting(config, [ "namelist:run_ukca", "jvsolar_file"], str("FJX_solcyc_May17.dat"))
        return config, self.reports


class vn108_t3451(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3451 by Stuart Whitehouse."""

    BEFORE_TAG = "vn10.8_t2587"
    AFTER_TAG = "vn10.8_t3451"

    WARNING = """WARNING! The list of namelists present in this file does not comply with
the metadata. The file contents have been updated to comply with the
metadata, but this may cause errors if the required namelists are not in
the app."""


    REQUIRED_NAMELISTS = {
        "file:ATMOSCNTL" : "namelist:configid namelist:nlstcgen namelist:nlst_mpp namelist:run_track namelist:run_calc_pmsl namelist:lbc_options namelist:run_nudging namelist:run_sl namelist:run_diffusion namelist:run_cosp namelist:radfcdia namelist:r2swclnl namelist:r2lwclnl namelist:clmchfcg namelist:acp namelist:acdiag namelist:jules_nvegparm namelist:jules_pftparm namelist:jules_triffid namelist:jules_elevate namelist:jules_urban2t_param namelist:iau_nl namelist:tuning_segments (namelist:nlstcall_pp(:)) (namelist:nlstcall_nc(:)) namelist:nlstcall_nc_options",
        "file:IDEALISE"  : "(namelist:idealise)",
        "file:IOSCNTL"   : "namelist:ioscntl namelist:io_control namelist:prnt_control (namelist:lustre_control) (namelist:lustre_control_custom_files)",
        "file:RECONA"    : "namelist:recon_technical namelist:recon_science namelist:recon_vertical namelist:horizont namelist:headers (namelist:trans(:))",
        "file:SHARED"    : "(namelist:nlcfiles) namelist:nlstcall namelist:ancilcta namelist:temp_fixes namelist:carbon_options namelist:coupling_control namelist:model_domain namelist:planet_constants namelist:run_dust namelist:run_glomap_aeroclim namelist:run_ukca namelist:run_gwd namelist:run_murk namelist:run_convection namelist:run_bl namelist:run_rivers namelist:run_precip namelist:run_radiation namelist:run_cloud namelist:run_aerosol namelist:lam_config namelist:run_ozone namelist:run_free_tracers namelist:run_eng_corr namelist:run_stochastic namelist:gen_phys_inputs namelist:run_dyn namelist:run_dyntest namelist:jules_surface_types namelist:jules_surface namelist:jules_radiation namelist:jules_hydrology namelist:jules_sea_seaice namelist:jules_soil namelist:jules_vegetation namelist:jules_soil_biogeochem namelist:jules_snow namelist:jules_urban_switches namelist:run_electric namelist:ancilcta (namelist:items(:))",
        "file:SIZES"     : "namelist:nlsizes",
        "file:STASHC"    : "(namelist:umstash_streq(:)) (namelist:umstash_domain(:)) (namelist:umstash_time(:)) (namelist:umstash_use(:)) (namelist:exclude_package(:))",
                         }

    def check_source(self, config, key, correct):
        source = self.get_setting_value(config, [key,"source"])
        if source != correct:
            self.change_setting_value(config, [key,"source"], correct)
            self.add_report(key, "source", source, self.WARNING, is_warning=True)
        return config


    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        for fname in self.REQUIRED_NAMELISTS.keys():
            config = self.check_source(config, fname, self.REQUIRED_NAMELISTS[fname])
        
        return config, self.reports


class vn108_t3476(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3476 by Stuart Whitehouse."""

    BEFORE_TAG = "vn10.8_t3451"
    AFTER_TAG = "vn10.8_t3476"

    def upgrade(self, config, meta_config=None):
        """
        Remove the uancnum namelist from any apps which still have it -
        it was deleted from the UM in #2016.
        """

        self.remove_setting(config, ["namelist:uancnum"])

        return config, self.reports


class vn108_t1824(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1824 by Paul Cresswell."""

    BEFORE_TAG = "vn10.8_t3476"
    AFTER_TAG = "vn10.8_t1824"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Split the idealise namelist in two.

        # Store if this is an idealised app or not, and what kind:
        is_idealised = False
        problem_number = self.get_setting_value(config,
                         ["namelist:idealise", "problem_number"])
        if problem_number != '0':
            is_idealised = True

        # The following three tests cases are mutually exclusive.
        # The switches for them are combined, preserving the order of priority.
        initial_profile = '0'
        if (self.get_setting_value(config,
           ['namelist:idealise', 'l_baro_inst']) == '.true.'):
            initial_profile = '1'
        elif (self.get_setting_value(config,
             ['namelist:idealise', 'l_solid_body']) == '.true.'):
            initial_profile = '2'
        elif (self.get_setting_value(config,
             ['namelist:idealise', 'l_deep_baro_inst']) == '.true.'):
            initial_profile = '3'
        self.add_setting(config, ['namelist:recon_idealised',
                                 'initial_profile'], initial_profile)

        # Rename l_initialise_data.
        self.rename_setting(config, ['namelist:idealise', 'l_initialise_data'],
                            ['namelist:recon_idealised', 'l_init_idealised'])

        # Now run through the rest of the items in namelist:idealise.
        # l_perturb and l_idl_bubble_saturate are to be deleted and so are
        # omitted from both lists.
        # eccentricity was retired by #1444 above.
        recon_items = ('l_cartesian', 'delta_xi1', 'delta_xi2', 'base_xi1',
          'base_xi2', 'grid_number', 'grid_flat', 'first_theta_height',
          'thin_theta_height', 'height_domain', 'big_layers', 'transit_layers',
          'big_factor', 'vert_grid_ratio', 'first_constant_r_rho_level_new',
          'surface_type', 'h_o', 'lambda_fraction', 'phi_fraction',
          'half_width_x', 'half_width_y',
          'tprofile_number', 'qprofile_number', 'pressure_balance',
          'dtheta_dz1', 'height_dz1', 'l_constant_dz',
          'aa_jet_m', 'aa_jet_n', 'aa_jet_u0', 'aa_jet_a',
          'grid_np_lon', 'grid_np_lat',
          'l_rotate_grid', 'l_baro_perturbed',
          't0_p', 't0_e', 'b_const', 'k_const',
          'l_perturb_t', 'perturb_magnitude_t', 'l_perturb_q',
          'perturb_magnitude_q', 'perturb_type', 'perturb_height',
          'idl_bubble_option', 'idl_bubble_max', 'idl_bubble_xoffset',
          'idl_bubble_yoffset', 'idl_bubble_width', 'idl_bubble_depth',
          'idl_bubble_height',
          'num_theta_init_heights', 'theta_init_field_type',
          'num_mv_init_heights', 'mv_init_field_type', 'num_uv_init_heights',
          'theta_init_height', 'mv_init_height', 'uv_init_height',
          'theta_init_data', 'mv_init_data', 'u_init_data', 'v_init_data',
          'initial_tolerance', 'profile_filename')

        model_items = ('l_shallow', 'l_const_grav', 'nxi1l',
          'nxi1v', 'nxi2l', 'nxi2v', 'delta_xi1_h', 'delta_xi1_l',
          'delta_xi2_h', 'delta_xi2_l', 'l_vert_coriolis', 'l_fixed_lbcs',
          'l_force_lbc', 'tstep_plot_frequency', 'tstep_plot_start',
          'roughlen_z0h', 't_surface', 'p_surface',
          'l_spec_z0', 'roughlen_z0m', 'f_plane', 'ff_plane',
          'idlsurffluxseaoption', 'idlsurffluxseaparams', 'l_geo_for', 'u_geo',
          'v_geo', 'num_theta_relax_heights', 'num_theta_relax_times',
          'num_mv_relax_heights', 'num_mv_relax_times', 'num_uv_relax_heights',
          'num_uv_relax_times', 'theta_relax_timescale', 'mv_relax_timescale',
          'uv_relax_timescale', 'theta_relax_height', 'theta_relax_time',
          'theta_relax_data', 'mv_relax_height', 'mv_relax_time',
          'mv_relax_data', 'uv_relax_height', 'uv_relax_time', 'u_relax_data',
          'v_relax_data', 'num_theta_inc_times', 'num_theta_inc_heights',
          'theta_inc_field_type', 'num_mv_inc_times', 'num_mv_inc_heights',
          'num_uv_inc_times', 'num_uv_inc_heights', 'theta_inc_time',
          'theta_inc_height', 'mv_inc_time', 'mv_inc_height', 'uv_inc_time',
          'uv_inc_height', 'theta_inc_data', 'mv_inc_data', 'u_inc_data',
          'v_inc_data', 'l_ideal_2d', 'local_heat_option',
          'local_heat_xoffset', 'local_heat_yoffset', 'local_heat_amp',
          'local_heat_sigma', 'local_heat_base', 'local_heat_top',
          'local_heat_period', 'tforce_number', 'trelax_number',
          'nsteps_consv_print', 'suhe_newtonian_timescale_ka',
          'suhe_newtonian_timescale_ks', 'suhe_pole_equ_deltat',
          'suhe_static_stab', 'suhe_sigma_cutoff', 'suhe_fric',
          'base_frictional_timescale', 'l_heldsuarez', 'l_heldsuarez1_drag')

        # Realistically there are only two options for qprofile_number when
        # problem_number=3 (idealised_problem): qp_dump and everything else.
        # We replace these with a logical to simplify the option structure.
        # qprofile_number is still used by problem_number=5 (idealised_planet).
        l_reset_mixing = '.false.'

        is_cartesian = False

        # Transfer items to recon_idealised
        for item in recon_items:
            value = self.get_setting_value(config, ["namelist:idealise", item])

            # Store if this is a Cartesian app:
            if item == 'l_cartesian':
                if value == '.true.':
                    is_cartesian = True

            # Deal with some retired values. These were previously allowed
            # by the metadata but had no code to support them.
            if item == 'grid_flat':
                if value in ['1', '2']:  # Retired linear options
                    value = '0'   # Only remaining 'linear' option.

            if item == 'grid_number':
                if value in ['22', '23']:  # Retired quadratic options
                    value = '21'  # Nearest remaining quad. heights structure

            if item == 'surface_type':
                if value in ['1', '2', '3', '4', '5'] : # Retired presets
                    value = '10'  # Revert to input dump

            if item == 'tprofile_number':
                if value in ['6', '7', '8']:  # Retired presets/by namelist
                    value = '10'  # Revert to input dump

            # If problem_number=2 this was previously trigger-ignored and
            # picking up a default, so we reset the value to that.
            if item == 'initial_tolerance':
                if problem_number == '2':
                    value = '1.0e-11'

            # qprofile_number has many options, some of which are only meant
            # for problem_number=3, some only for problem_number=5, but most of
            # which have the same outcome for any given app anyway. We split it
            # into two simpler variables, one for each problem type.
            # Both values are set regardless of the value of problem_number;
            # this ensures qprofile_number is valid even if trigger-ignored,
            # and both options have the 'nearest' setting to before.
            if item == 'qprofile_number':
                if value != '10':
                    l_reset_mixing = '.true.'
                if value not in ['10', '11']:
                    value = '0'

            if (item == 'num_mv_init_heights'
                or item == 'num_theta_init_heights'
                or item == 'num_uv_init_heights'):
                if value == '0':  # Previously allowed but invalid
                    value = '1'

            self.add_setting(config, ["namelist:recon_idealised", item], value)

        self.add_setting(config, ["namelist:recon_idealised",
                                  "l_reset_mixing"], l_reset_mixing)

        # Transfer the remaining items to the renamed idealised namelist.
        for item in model_items:
            value = self.get_setting_value(config, ["namelist:idealise", item])

            # An idealised, non-Cartesian app would previously have had these
            # values trigger-ignored. Now they will be active (an app
            # configured for atmos-only work cannot know if the input dump is
            # Cartesian), so set them to 'off':
            if is_idealised and not is_cartesian:
                if item == 'f_plane' or item == 'ff_plane':
                    value = '-90.0'
                if item == 'l_vert_coriolis':
                    value = '.false.'


            self.add_setting(config, ["namelist:idealised", item], value)

        # Delete the old idealise namelist.
        self.remove_setting(config, ["namelist:idealise"])

        # Insert the new namelists in the right files.
        recona = self.get_setting_value(config, ["file:RECONA","source"])
        if recona:
            recona = re.sub('\(namelist:trans\(:\)\)',
                            '(namelist:trans(:)) (namelist:recon_idealised)',
                            recona)
            self.change_setting_value(config, ["file:RECONA", "source"],
                                      recona)

        idealise = self.get_setting_value(config, ["file:IDEALISE","source"])
        if idealise:
            idealise = re.sub('\(namelist:idealise\)',
                              '(namelist:idealised)', idealise)
            self.change_setting_value(config, ["file:IDEALISE", "source"],
                                      idealise)

        return config, self.reports


class vn108_t1971(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1971 by Malcolm Brooks."""

    BEFORE_TAG = "vn10.8_t1824"
    AFTER_TAG = "vn10.8_t1971"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Adds pwsdiag_sfc_em
        self.add_setting(config, ["namelist:run_dust", "pwsdiag_sfc_em"],
                         "0.0")
        return config, self.reports


class vn108_t1711(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1711 by Adam Clayton."""

    BEFORE_TAG = "vn10.8_t1971"
    AFTER_TAG = "vn10.8_t1711"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:iau_nl",
                                  "l_iau_removess"], ".true.")
        return config, self.reports


class vn108_t2710(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #2710 by Adam Clayton."""

    BEFORE_TAG = "vn10.8_t1711"
    AFTER_TAG = "vn10.8_t2710"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:temp_fixes",
                                  "l_fix_iau_rim_density"], ".false.")
        return config, self.reports


class vn108_t1814(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket Easy Aerosol #1814 by Ben Johnson."""

    BEFORE_TAG = "vn10.8_t2710"
    AFTER_TAG = "vn10.8_t1814"

    def upgrade(self, config, meta_config=None):
        """Upgrade UM app configurations."""

        # Add namelist easyaerosol to the list of namelists that are output in the 
        # SHARED file                  
        current_value = self.get_setting_value(config,["file:SHARED", "source"])
        if current_value != None:
              new_value = current_value.replace(' namelist:run_dyntest ',
                                                ' namelist:run_dyntest namelist:easyaerosol ')

              self.change_setting_value(config,["file:SHARED", "source"],
                                        value=new_value, forced=True)

        self.add_setting(config, ["namelist:easyaerosol", "easyaerosol_dir"], "'unset'")
        self.add_setting(config, ["namelist:easyaerosol", "easyaerosol_files"], ",".join(["'unset'"]*7))
        self.add_setting(config, ["namelist:easyaerosol", "l_easyaerosol_autoconv"], ".false.")
        self.add_setting(config, ["namelist:easyaerosol", "l_easyaerosol_cdnc"], ".false.")
        self.add_setting(config, ["namelist:easyaerosol", "l_easyaerosol_lw"], ".false.")
        self.add_setting(config, ["namelist:easyaerosol", "l_easyaerosol_sw"], ".false.")
        self.add_setting(config, ["namelist:easyaerosol", "l_easyaerosol_zonal"], ".false.")
        self.add_setting(config, ["namelist:radfcdia", "c2c_easy_d"], ".false.")
        return config, self.reports

class vn109_t3471(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.8_t1814"
    AFTER_TAG = "vn10.9"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "10.9", forced=True)
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMASTER"])
        if stashmaster_path:
           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",
                                     "vn10.9/ctldata",
                                      stashmaster_path)
           self.change_setting_value(config, ["env", "STASHMASTER"],
                                     stashmaster_path, forced=True)
        return config, self.reports


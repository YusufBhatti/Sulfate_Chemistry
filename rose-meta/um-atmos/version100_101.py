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



class vn100_t45(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #45 by Steve Wardle."""

    BEFORE_TAG = "vn10.0"
    AFTER_TAG = "vn10.0_t45"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Note: this upgrade macro intentionally does nothing.
        # It is required because as part of #45 the transform
        # macro "stash_requests.StashProfilesRemoveUnused" will
        # be run on all apps in the trunk and this requires
        # an upgrade version to apply against (otherwise all
        # apps would have their metadata revert to "HEAD")

        return config, self.reports

class vn100_t47(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #47 by Warren Tennant."""

    BEFORE_TAG = "vn10.0_t45"
    AFTER_TAG = "vn10.0_t47"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_stochastic", "l_stphseed_file"], ".false.")
        self.add_setting(config, ["namelist:run_stochastic","rp2_callfreq"], "10800")
        self.add_setting(config, ["namelist:temp_fixes", "l_stph_rhcrit_unbias" ], ".false.")
        return config, self.reports

class vn100_t67(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #67 by Claudio Sanchez"""

    BEFORE_TAG = "vn10.0_t47"
    AFTER_TAG = "vn10.0_t67"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_stochastic", "sd_orog_thres"], "500.")
        self.add_setting(config, ["namelist:run_stochastic", "psif_orog_thres"], "0.5")
        return config, self.reports

class vn100_t137(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #137 by Glenn Greed.
       Add run_dyn namelist to shared file."""

    BEFORE_TAG = "vn10.0_t67"
    AFTER_TAG = "vn10.0_t137"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        shared_source = self.get_setting_value(config, ["file:SHARED","source"])
        if shared_source:
            shared_source = shared_source.replace('namelist:run_eng_corr',
                    'namelist:run_eng_corr namelist:run_dyn')
            self.change_setting_value(config, ["file:SHARED","source"], shared_source)
        l_endgame = self.get_setting_value(config, ["namelist:nlstcatm", "l_endgame"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_endgame"])
        self.add_setting(config, ["namelist:run_dyn","l_endgame"], l_endgame)
        return config, self.reports


class vn100_t134(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #134 by <Ricky Wong>."""

    BEFORE_TAG = "vn10.0_t137"
    AFTER_TAG = "vn10.0_t134"

    def upgrade(self, config, meta_config=None):
        """Make model_domain namelist and move 
           model_domain, l_regular from namelist nlstcatm
           to namelist model_domain. In addition rename model_domain
           to model_type in the move."""

        # Scm jobs are defined as model type number 5, we'll use this to test
        # if its a scm app
        mt_scm = 5

        # Add new section for namelist:model_domain
        self.add_setting(config, ["namelist:model_domain"], value=None,
                         forced=False, state=None,comments=None,info=None)

        # Get setting values for l_regular and model_domain from
        # nlstcatm which will be moved to the new namelist, model_domain
        model_domain_value = self.get_setting_value(config,
                                                  ["namelist:nlstcatm", "model_domain"])

        if model_domain_value != None:
            self.add_setting(config, ["namelist:model_domain", "model_type"],
                             value=model_domain_value, forced=True)

        l_regular_value = self.get_setting_value(config,
                                                  ["namelist:nlstcatm", "l_regular"])
        if l_regular_value != None:
            self.add_setting(config, ["namelist:model_domain", "l_regular"],
                             value=l_regular_value, forced=True)

        # Add namelist model_domain to the list of namelists that are output in the 
        # SHARED file (For Recon Jobs)
        current_value = self.get_setting_value(config,["file:SHARED", "source"])
        if current_value != None:
              new_value = current_value.replace(' namelist:temp_fixes ',
                                                ' namelist:temp_fixes namelist:model_domain ')

              self.change_setting_value(config,["file:SHARED", "source"],
                                        value=new_value, forced=True)

        # Add namelist model_domain to the list of namelists that are output in the 
        # NAMELIST file (for Non-Scm Apps)
        if model_domain_value != mt_scm:
            current_value = self.get_setting_value(config,["file:NAMELIST", "source"])
            if current_value != None:
                  new_value = current_value.replace(' namelist:temp_fixes ',
                                                    ' namelist:temp_fixes namelist:model_domain ')

                  self.change_setting_value(config,["file:NAMELIST", "source"],
                                            value=new_value, forced=True)

        # Remove model_domain and l_regular from the namelist nlstcatm
        self.remove_setting(config,["namelist:nlstcatm", "model_domain"] )
        self.remove_setting(config,["namelist:nlstcatm", "l_regular"] )
                   
        return config, self.reports

class vn100_t53(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #53 by Ian Boutle."""

    BEFORE_TAG = "vn10.0_t134"
    AFTER_TAG = "vn10.0_t53"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.remove_setting(config, ["namelist:run_bl","weightlouistolong"])
        self.remove_setting(config, ["namelist:run_bl","trweights1"])
        return config, self.reports


class vn100_t12(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #12 by richardbarnes."""

    BEFORE_TAG = "vn10.0_t53"
    AFTER_TAG = "vn10.0_t12"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        alphac = self.get_setting_value(config, ["namelist:run_stochastic", "alphac"])
        n1 = self.get_setting_value(config, ["namelist:run_stochastic", "n1"])
        n2 = self.get_setting_value(config, ["namelist:run_stochastic", "n2"])
        tau = self.get_setting_value(config, ["namelist:run_stochastic", "tau"])
        nsmooth = self.get_setting_value(config, ["namelist:run_stochastic", "nsmooth"])

        self.remove_setting(config, ["namelist:run_stochastic", "alphac"])
        self.remove_setting(config, ["namelist:run_stochastic", "n1"])
        self.remove_setting(config, ["namelist:run_stochastic", "n2"])
        self.remove_setting(config, ["namelist:run_stochastic", "tau"])
        self.remove_setting(config, ["namelist:run_stochastic", "nsmooth"])
        """Reset to current value, if any"""
        self.add_setting(config, ["namelist:run_stochastic", "alphac"], alphac)
        self.add_setting(config, ["namelist:run_stochastic", "n1"], n1)
        self.add_setting(config, ["namelist:run_stochastic", "n2"], n2)
        self.add_setting(config, ["namelist:run_stochastic", "tau"], tau)
        if nsmooth == None:
             nsmooth = "3"
        self.add_setting(config, ["namelist:run_stochastic", "nsmooth"], nsmooth)

        return config, self.reports

class vn100_t193(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #193 by Rachel Stratton."""

    BEFORE_TAG = "vn10.0_t12"
    AFTER_TAG = "vn10.0_t193"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        t_snow_melt_local = self.get_setting_value(
                config, ["namelist:run_convection", "t_melt_snow"])
        # If commented out this variable has a non valid value of 0 set when the variable
        # was first added to the meta data. Need to reset to a valid value otherwise
        # leave as it is. 
        # There should not be any non-valid commented out values in any apps.
        # Any app with a non-zero value commented out will keep this value even if it is 
        # non-valid as this will prompt a user to think of an appropriate value to set.
        # 
        # Note at UM10.0 none of the rose stem app have non-zero values
        # of t_melt_snow so all get reset to 274.15 
        if t_snow_melt_local == "0":
            self.change_setting_value(config, ["namelist:run_convection", "t_melt_snow"], "274.15")
        # Adding new switches
        self.add_setting(config, ["namelist:run_convection", "l_new_dd"], ".false.")
        self.add_setting(config, ["namelist:run_convection", "l_use_dd"], ".false.")
        self.add_setting(config, ["namelist:run_convection", "l_mom_dd"], ".false.")
        # Adding temporary fix
        self.add_setting(config, ["namelist:temp_fixes", "l_glue_conv5a"], ".false.")
        
        return config, self.reports


class vn100_t246(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #246 by Colin Johnson."""

    BEFORE_TAG = "vn10.0_t193"
    AFTER_TAG = "vn10.0_t246"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_ukca", "mode_activation_dryr"], "37.5")
        # Input your macro commands here
        return config, self.reports


class vn100_t181(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #181 by Maggie Hendry."""

    BEFORE_TAG = "vn10.0_t246"
    AFTER_TAG = "vn10.0_t181"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.remove_setting(config, ["namelist:urban_switches","l_urban_empirical"])
        return config, self.reports

class vn100_t234(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #234 by peterhill."""

    BEFORE_TAG = "vn10.0_t181"
    AFTER_TAG = "vn10.0_t234"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_radiation", "two_d_fsd_factor"], "1.414")
        return config, self.reports

class vn100_t111(rose.upgrade.MacroUpgrade):

   """Upgrade macro for ticket #111 by Mohit Dalvi."""

   BEFORE_TAG = "vn10.0_t234"
   AFTER_TAG = "vn10.0_t111"

   def upgrade(self, config, meta_config=None):
      """Add Logical for AgeAir to UKCA namelist"""

      # Check the UKCA scheme being used and set the 
      # new logical as True for Stratospheric schemes
      # and False for others
      ukca_chem = self.get_setting_value (config, 
             ["namelist:run_ukca", "i_ukca_chem"] )
    
      if int(ukca_chem) == 51 or int(ukca_chem) == 52:
         self.add_setting(config, 
            ["namelist:run_ukca", "l_ukca_ageair"],
             ".true.")
      else:
         self.add_setting(config, 
            ["namelist:run_ukca", "l_ukca_ageair"],
             ".false.")

      return config, self.reports

class vn100_t244(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #244 by Ian Boutle."""

    BEFORE_TAG = "vn10.0_t111"
    AFTER_TAG = "vn10.0_t244"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["env","ENS_MEMBER"], "0")
        self.add_setting(config, ["namelist:run_stochastic","i_pert_theta"], "0")
        self.add_setting(config, ["namelist:run_stochastic","mag_pert_theta"],"0.5")
        self.add_setting(config, ["namelist:run_stochastic","lev_pert_theta"],"10")
        self.add_setting(config, ["namelist:run_stochastic","npts_pert_theta"],"8")

        return config, self.reports

class vn100_t51(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #51 by Ian Boutle."""

    BEFORE_TAG = "vn10.0_t244"
    AFTER_TAG = "vn10.0_t51"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:recon","l_adj_t_soil"], ".false.")
        self.add_setting(config, ["namelist:recon","l_init_tile_t_zerofrac"], ".false.")
        self.add_setting(config, ["namelist:run_dyn","l_sl_bc_correction"], ".false.")
        return config, self.reports

class vn100_t197(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #197 by Irina Linova-Pavlova."""

    BEFORE_TAG = "vn10.0_t51"
    AFTER_TAG = "vn10.0_t197"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Request to calculate sections 34 or 50 on pressure levels
        self.add_setting(config, ["namelist:run_ukca", "l_ukca_chem_plev"],
                        ".false.")
        self.add_setting(config, ["namelist:run_ukca", "l_ukca_asad_plev"],
                        ".false.")

        return config, self.reports
        
class vn100_t159(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #159 by Chris Smith."""

    BEFORE_TAG = "vn10.0_t197"
    AFTER_TAG = "vn10.0_t159"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_sl","l_priestley_correct_thetav"],".false.")
        return config, self.reports
        

class vn100_t110(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #110 by richardbarnes."""

    BEFORE_TAG = "vn10.0_t159"
    AFTER_TAG = "vn10.0_t110"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Create "filler" arrays for each of the arrays
        full_levels = [["-32768.0" for _ in range(36) ] for _ in range(12)]
        full_rates  = [["-32768.0" for _ in range(36) ] for _ in range(12)]
        full_years  = [["-32768"   for _ in range(36) ] for _ in range(12)]
        full_nyears = ["0" for _ in range(12)]

        variables_with_defaults = [
            [ "clim_fcg_levls",  full_levels],
            [ "clim_fcg_nyears", full_nyears],
            [ "clim_fcg_rates",  full_rates],
            [ "clim_fcg_years",  full_years]
                                  ]

        index_to_species_mapping = [[0, "co2"],
                                    [1, "ch4"],
                                    [2, "n2o"],
                                    [3, "cfc11"],
                                    [4, "cfc12"],
                                    [5, "so4"],
                                    [7, "cfc113"],
                                    [8, "hcfc22"],
                                    [9, "hfc125"],
                                    [10,"hfc134a"],
                                    [11,"cfc114"]]

        namelist = "namelist:clmchfcg"

        l_clmchfcg = self.get_setting_value(config, [namelist,"l_clmchfcg"])

        namelist_contents = self.get_setting_value(config, [namelist])

        items_to_delete = []
        items_to_add = []

        # Pre-compile the search pattern for the array slices.  This pattern
        # will expect a 2d index/slice
        pattern_indices_2d = re.compile(
                            r"\(((?P<low_1>\d+)|)((?P<c_1>:)|)((?P<up_1>(\d+))|),"
                            r"((?P<low_2>\d+)|)((?P<c_2>:)|)((?P<up_2>(\d+))|)\)")
        pattern_indices_1d = re.compile(
                            r"\(((?P<low_1>\d+)|)((?P<c_1>:)|)((?P<up_1>(\d+))|)\)")
       
        # Process each instance of each variable in the config object
        for varname, full_array in variables_with_defaults:
            for obj in namelist_contents:
                if re.search(varname, obj):
                    # When we find one, add it to the list of things to delete...
                    items_to_delete.append([namelist, obj])

                    # If the object has brackets at the end (i.e. slice/index)
                    # which imply it is one of the 2d arrays
                    if re.search(r"\(.*,.*\)\s*$", obj):
                        result = re.search(pattern_indices_2d, obj)

                        # Process the first index
                        if result.group("low_1") is None:
                            low_1 = 0
                        else:
                            low_1 = int(result.group("low_1")) - 1
                        if result.group("c_1") is None:
                            # No slice specified, upper index == lower index
                            up_1 = low_1
                        else:
                            if result.group("up_1") is None:
                                # Slice with no upper bound
                                up_1 = 36 - 1
                            else:
                                # Slice with upper bound
                                up_1 = int(result.group("up_1")) - 1

                        # Process the second index
                        if result.group("low_2") is None:
                            low_2 = 0
                        else:
                            low_2 = int(result.group("low_2")) - 1
                        if result.group("c_2") is None:
                            # No slice specified, upper index == lower index
                            up_2 = low_2
                        else:
                            if result.group("up_2") is None:
                                # Slice with no upper bound
                                up_2 = 12 - 1
                            else:
                                # Slice with upper bound
                                up_2 = int(result.group("up_2")) - 1

                        # Now read in the values and put them into the array
                        values = self.get_setting_value(config, [namelist, obj])
                        values = values.split(",")

                        # If the array dimensions suggest that only the first
                        # element of the dimensions is specified, adjust it
                        expected_length = ((up_1 - low_1 + 1) *
                                           (up_2 - low_2 + 1))
                        if len(values) != expected_length:
                            if result.group("c_1") is None:
                                if len(values) % 36 == 0 and up_1 == 0:
                                    up_1 = 36 - 1

                        expected_length = ((up_1 - low_1 + 1) *
                                           (up_2 - low_2 + 1))
                        if len(values) != expected_length:
                            if result.group("c_2") is None:
                                if len(values) % 12 == 0 and up_2 == 0:
                                    up_2 = 12 - 1

                        # Now read the values into the right positions in the
                        # full array
                        val_idx = 0
                        for idx2 in range(low_2, up_2 + 1):
                            for idx1 in range(low_1, up_1 + 1):
                                full_array[idx2][idx1] = values[val_idx]
                                val_idx = val_idx + 1

                    # Alternatively, this is the 1d array case
                    elif re.search(r"\(.*\)\s*$", obj):
                        result = re.search(pattern_indices_1d, obj)

                        # Process the first index
                        low_1 = int(result.group("low_1")) - 1
                        if result.group("c_1") is None:
                            # No slice specified, upper index == lower index
                            up_1 = low_1
                        else:
                            if result.group("up_1") is None:
                                # Slice with no upper bound
                                up_1 = 12 - 1
                            else:
                                # Slice with upper bound
                                up_1 = int(result.group("up_1")) - 1

                        # Now read in the values and put them into the array
                        values = self.get_setting_value(config, [namelist, obj])
                        values = values.split(",")

                        val_idx = 0
                        for idx1 in range(low_1, up_1 + 1):
                            full_array[idx1] = values[val_idx]
                            val_idx = val_idx + 1

            # Now that the "full" array contains the data to be written out, we can
            # extract the appropriate columns

            # All the compilers we have access to (including CCE on the XC40) need 2D arrays
            # to have exact bounds, but 1D arrays don't, so we discard that information
            # in the interests of simplicity.

            # Generate the new variable name by converting the array location number to
            # a text suffix
            for out_index, species_name in index_to_species_mapping:
                new_var = "%s_%s"%(varname, species_name)
                if type(full_array[0]) is not list:
                    self.add_setting(config, [namelist, new_var], full_array[out_index])
                else:
                    self.add_setting(config, [namelist, new_var], ",".join(full_array[out_index]))
            
        # Delete the old items
        for item in items_to_delete:
            self.remove_setting(config, item)
        
        return config, self.reports


class vn100_t105(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #105 by Chris Smith."""

    BEFORE_TAG = "vn10.0_t110"
    AFTER_TAG = "vn10.0_t105"

    def upgrade(self, config, meta_config=None):
        """Interpolation options for tracers are same as for moisture."""
        # Monotone scheme
        ms_old=self.get_setting_value(
                         config, ["namelist:run_sl", "monotone_scheme"])
        # If full array specified, copy 2nd element to new 5th element,
        # giving tracers same interpolation options as moisture.
        if ',' in ms_old:
            ms_old=ms_old.split(',') 
            ms_old.append(ms_old[1])
            ms_new=','.join(ms_old)
        else:
            ms_new=ms_old

        self.change_setting_value(config, ["namelist:run_sl",
                                  "monotone_scheme"], str(ms_new))
         # High order scheme
        hs_old=self.get_setting_value(
                         config, ["namelist:run_sl", "high_order_scheme"])
        # If full array specified, copy 2nd element to new 5th element,
        # giving tracers same interpolation options as moisture.
        if ',' in hs_old:
            hs_old=hs_old.split(',') 
            hs_old.append(hs_old[1])
            hs_new=','.join(hs_old)
        else:
            hs_new=hs_old

        self.change_setting_value(config, ["namelist:run_sl",
                                  "high_order_scheme"], str(hs_new))
        return config, self.reports


class vn100_t113(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #113 by Steve Wardle"""

    BEFORE_TAG = "vn10.0_t105"
    AFTER_TAG = "vn10.0_t113"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Work out if this is an SCM app (these need different treatment below)
        model_domain = self.get_setting_value(config, ["namelist:model_domain", "model_type"])
        is_scm = (model_domain == "5")

        # The existing time convention dictates what sort of separator to use
        # as well as which of the "legacy" flags to use in the new filename
        # pattern
        time_convention = self.get_setting_value(
            config, ["namelist:nlstcall", "time_convention"])
        if time_convention == "'Absolute_Dstamp'":
            separator = "."
            legacy_flag = "%C"
        elif time_convention == "'Relative'":
            separator = "_"
            legacy_flag = "%N"
        elif time_convention == "'Sub-hourly'":
            separator = "_"
            legacy_flag = "%n"
        elif time_convention == "'Timestep'":
            separator = "_"
            legacy_flag = "%T"
        else:
            raise UpgradeError("Unrecognised time_convention: {0:s}"
                               .format(time_convention))
        # Time convention isn't required after this upgrade
        self.remove_setting(config, ["namelist:nlstcall", "time_convention"])

        # Re-format the nlstcall_pp namelists (output streams)
        for section in config.value.keys():
            node = config.get([section])
            if (isinstance(node.value, dict) and
                re.match("namelist:nlstcall_pp\(.*?\)$", section)):
                # Get the existing reinit letter value for the stream
                letter = self.get_setting_value(
                    config, [section, "reinit_letter"]).replace("'","")
                # And the current file id to set the filetype appropriately
                file_id = self.get_setting_value(
                    config, [section, "file_id"]).replace("'","")
                if file_id == "ppvar" or file_id == "ppbmc":
                    filetype = "c"
                else:
                    filetype = "p"

                # Generate an appropriate filename pattern, using the legacy
                # indicators to reproduce the same filename as before
                filename = ("'$DATAM/${RUNID}a" +
                            "{0:s}{1:s}{2:s}{3:s}'"
                            .format(separator,filetype,letter,legacy_flag))
                self.add_setting(config, [section, "filename_base"], filename)
                self.remove_setting(config, [section, "reinit_letter"])

        # Create the new dump filename pattern, again using legacy indicators
        self.add_setting(config, ["namelist:nlstcgen", "dump_filename_base"],
                         "'$DATAM/${RUNID}a"+
                         "{0:s}d%z{1:s}'".format(separator, legacy_flag))

        # As above but for the mean filename patterns, note that these all use
        # the "Absolute" style legacy convention, there are also 4 (one for
        # each of the meaning periods)
        mean_basename = "'$DATAM/${RUNID}a"+"{0:s}p%C'".format(separator)
        for i_mean in range(4):
            self.add_setting(config, ["namelist:nlstcgen",
                                      "mean_{0:d}_filename_base"
                                      .format(i_mean+1)], mean_basename)

        # The partial sum files don't use the generation in the same way, so
        # just create basenames for them
        self.add_setting(config, ["namelist:nlstcgen", "psum_filename_base"],
                         "'$DATAM/${RUNID}a_s'")

        # Add the new input for the STASH requests log
        if is_scm:
            self.add_setting(config, ["namelist:nlcfiles", "streqlog"], "'unset'")
        else:
            self.add_setting(config, ["namelist:nlcfiles", "streqlog"],
                            "'$DATAW/$RUNID.stash'")

        # Remove the nlcfiles entries which have been deleted (probably only
        # used rarely but explicitly removing them will protect apps)
        self.remove_setting(config, ["namelist:nlcfiles", "icfile"])
        self.remove_setting(config, ["namelist:nlcfiles", "aopsum3"])
        self.remove_setting(config, ["namelist:nlcfiles", "aopsum4"])
        self.remove_setting(config, ["namelist:nlcfiles", "aopstmp3"])
        self.remove_setting(config, ["namelist:nlcfiles", "aopstmp4"])        

        # Move l_postp to io_control namelist and ppxm to nlstcgen namelist
        l_postp = self.get_setting_value(config, ["namelist:nlstcall_qs", "l_postp"])
        self.add_setting(config, ["namelist:io_control", "l_postp"], l_postp)
        self.remove_setting(config, ["namelist:nlstcall_qs", "l_postp"])

        ppxm = self.get_setting_value(config, ["namelist:nlstcall_qs", "ppxm"])
        self.add_setting(config, ["namelist:nlstcgen", "ppxm"], ppxm)
        self.remove_setting(config, ["namelist:nlstcall_qs", "ppxm"])

        # Delete the nlstcall_qs namelist and remove it from the file
        self.remove_setting(config, ["namelist:nlstcall_qs"])
        namelist_file_list = self.get_setting_value(config, ["file:NAMELIST", "source"])
        if namelist_file_list:
            self.change_setting_value(config, ["file:NAMELIST", "source"],
                                    namelist_file_list.replace("namelist:nlstcall_qs ",""))

        # Get the existing value of RUNID, if it was set at the app level,
        # otherwise use the reference to it
        runid = self.get_setting_value(config, ["env", "RUNID"])
        if runid is None:
            runid = "${RUNID}"
        if not is_scm:
            self.add_setting(config, ["env", "RUNID"], "$RUNID")
        

        # Get the existing value of DATAM, if it was set at the app level,
        # otherwsie use the reference to it
        datam = self.get_setting_value(config, ["env", "DATAM"])
        if datam is None:
            datam = "$DATAM"

        # Likewise for DATAW
        dataw = self.get_setting_value(config, ["env", "DATAW"])
        if dataw is None:
            dataw = "$DATAW"

        # Add an entry to the app to ensure the DATAM directory gets created and
        # ensure it exists at the app level to improve the way this appears to
        # the user if they aren't passing it from the suite
        if not is_scm:
            self.add_setting(config, ["env", "DATAM"], datam)
            self.add_setting(config, ["file:$DATAM"])
            self.add_setting(config, ["file:$DATAM", "mode"], "mkdir")
            self.add_setting(config, ["file:$DATAM", "source"], "")

        # And the same for DATAW
        if not is_scm:
            self.add_setting(config, ["env", "DATAW"], dataw)
            self.add_setting(config, ["file:$DATAW"])
            self.add_setting(config, ["file:$DATAW", "mode"], "mkdir")
            self.add_setting(config, ["file:$DATAW", "source"], "")

        # Add the path to the history file
        if not is_scm:
            self.add_setting(config, ["env", "HISTORY"], datam + "/" + runid + ".xhist")

        # Get the current value of the PE output directory (if it is set)
        current_stdout_dir = self.get_setting_value(config,
                                                    ["env", "ATMOS_STDOUT_DIR"])
        if current_stdout_dir is None:
            stdout_dir = "pe_output"
        else:
            stdout_dir = current_stdout_dir

        # Add the full basename path of the PE output files, unless it would
        # result in the value being the default (see description of HISTORY)
        if runid != "atmos" and not is_scm:
            self.add_setting(config, ["env", "ATMOS_STDOUT_FILE"],
                             stdout_dir+"/"+runid+".fort6.pe")
        self.remove_setting(config, ["env", "ATMOS_STDOUT_DIR"])

        # Get the current value of the PE output directory (if it is set)
        current_stdout_dir = self.get_setting_value(config,
                                                    ["env", "RECON_STDOUT_DIR"])
        if current_stdout_dir is None:
            stdout_dir = "pe_output"
        else:
            stdout_dir = current_stdout_dir

        # Add the full basename path of the PE output files, unless it would
        # result in the value being the default (see description of HISTORY)
        if runid != "atmos" and not is_scm:
            self.add_setting(config, ["env", "RECON_STDOUT_FILE"],
                             stdout_dir+"/"+runid+".fort6.pe")
        self.remove_setting(config, ["env", "RECON_STDOUT_DIR"])

        # Try to make the user aware that some of the settings aren't actually
        # used anymore, and that there is a chance they may need to take a
        # little further action
        msg = """
        !!!!! NOTE: RUNID, DATAM and DATAW are no longer required   !!!!! 
        !!!!! or used directly by the UM.                           !!!!!
        !!!!!                                                       !!!!!
        !!!!! They are now merely environment variables defined by  !!!!!
        !!!!! the app/task/suite to provide common elements of some !!!!!
        !!!!! filenames.                                            !!!!!
        !!!!!                                                       !!!!!
        !!!!! It is now up to the suite to make sure that           !!!!!
        !!!!! DATAM/DATAW exist if they are required - this macro   !!!!!
        !!!!! has added the relevant sections to the apps to        !!!!!
        !!!!! attempt to always create DATAM and DATAW.             !!!!!
        !!!!!                                                       !!!!!
        !!!!! If your suite/app did not set DATAM, DATAW or RUNID   !!!!!
        !!!!! but relied on picking up the defaults in the um-atmos !!!!!
        !!!!! script, you will need to ensure these are set.        !!!!!
              """
        self.add_report(info=msg, is_warning=True)

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
class vn100_t463(vn91_t6297):

    """Upgrade macro for ticket #463 by Steve Wardle."""

    BEFORE_TAG = "vn10.0_t113"
    AFTER_TAG = "vn10.0_t463"

    # Intentionally blank (see the comment above)

class vn100_t463a(rose.upgrade.MacroUpgrade):

    """2nd upgrade macro for ticket #463 by Steve Wardle."""

    BEFORE_TAG = "vn10.0_t463"
    AFTER_TAG = "vn10.0_t463a"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Remove the idealise namelist entry (if it existed)
        self.remove_setting(config, ["namelist:nlcfiles", "idealise"])

        return config, self.reports


class vn100_t1389(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1389 by Paul Cresswell."""

    BEFORE_TAG = "vn10.0_t463a"
    AFTER_TAG = "vn10.0_t1389"

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


class vn101_t327(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn10.0_t1389"
    AFTER_TAG = "vn10.1"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "10.1")
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMSTR"])
        if stashmaster_path:
           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",
                                     "vn10.1/ctldata",
                                      stashmaster_path)
           self.change_setting_value(config, ["env", "STASHMSTR"],
                                     stashmaster_path)
        return config, self.reports


import rose.upgrade
import rose.namelist_dump
import re
import sys
import os
import tempfile

class UpgradeError(Exception):

      """Exception created when an upgrade fails."""

      def __init__(self, msg):
          self.msg = msg

      def __repr__(self):
          sys.tracebacklimit = 0
          return self.msg

      __str__ = __repr__


class vn91_t6314b(rose.upgrade.MacroUpgrade):
      
    """Upgrade macro for ticket #6314 for Joe Mancell"""

    BEFORE_TAG = "vn9.1"
    AFTER_TAG = "vn9.1_t6314"
    def upgrade(self, config, meta_config=None):
          """Upgrade a UM runtime app configuration."""
          # Add compulsory = true env vars
          self.add_setting(config,             
                           ['env', 'ATMOS_KEEP_MPP_STDOUT'], "false")
          self.add_setting(config,             
                           ['env', 'RECON_KEEP_MPP_STDOUT'], "false")
          self.add_setting(config,             
                           ['env', 'PRINT_STATUS'], "PrStatus_Normal")
          self.add_setting(config,             
                           ['env', 'RCF_PRINTSTATUS'], "PrStatus_Normal")
          self.add_setting(config,
                           ['env', 'RCF_TIMER'], "false")

          return config, self.reports

class vn91_t6367(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6367 by <Glenn Greed>."""

    BEFORE_TAG = "vn9.1_t6314"
    AFTER_TAG = "vn9.1_t6367"

    def upgrade(self, config, meta_config=None):
        """Make streq, domain time and use optional namelists."""
        # Check if the app file contains STASHC related namelists
        var_list = ("namelist:streq","namelist:domain",
                    "namelist:time", "namelist:use"  )

        namelists=[]
        for profile in var_list:
            stashc_streq = self.get_setting_value(config, 
                                        ["file:STASHC","source"])
            if stashc_streq:
                old_string=profile+"(:)"
                new_string="("+old_string+")"
                search_string=r""+re.escape(old_string)

                if (re.search(search_string, stashc_streq)):
                # Allow STASHC related namelists to be optional if not already
                    search_string=r""+re.escape(new_string)
                    if ( not (re.search(search_string, stashc_streq))):
                        stashc_streq = stashc_streq.replace(old_string,new_string)
                        self.change_setting_value(config, ["file:STASHC","source"],
                                                           stashc_streq)
                else:
                    self.change_setting_value(config, ["file:STASHC","source"],
                                                    stashc_streq+" "+new_string)
            else:
                # create list of namelists to remove as redundant if no STASHC
                for obj in config.get_value():
                    search_string=r""+re.escape(profile)
                    if re.search(search_string,obj):
                        namelists.append(obj)

        for name in namelists:
            self.remove_setting(config, [name])
                   
        return config, self.reports

class vn91_t6326(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6326 by Claudio Sanchez."""

    BEFORE_TAG = "vn9.1_t6367"
    AFTER_TAG = "vn9.1_t6326"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
	
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_stochastic", "l_spt_mse"], ".false.")
        return config, self.reports

class vn91_t5757(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5757 by James Manners."""

    BEFORE_TAG = "vn9.1_t6326"
    AFTER_TAG = "vn9.1_t5757"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_radiation","l_t_bdy_surf"], ".false.")
        return config, self.reports

class vn91_t6151(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6151 by SJ Swarbrick Group update ancil requests which
    relate to the same ancillary file."""

    BEFORE_TAG = "vn9.1_t5757"
    AFTER_TAG = "vn9.1_t6151"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        reqs_group = []
        reqs_done = []
        rem_list = []
        req_list = ''
        
        for obj in config.get_value():
        
          if re.search(r'namelist:items', obj):
          
          # Loop through items name lists, identify the ones which have updating on
          
            update_anc = self.get_setting_value(config, [obj, "update_anc"])
     
            if update_anc == ".true.":

             ancfile1   = self.get_setting_value(config, [obj, "ancilfilename"])
             domain     = self.get_setting_value(config, [obj, "domain"])
             interval   = self.get_setting_value(config, [obj, "interval"])
             period     = self.get_setting_value(config, [obj, "period"])
             source     = self.get_setting_value(config, [obj, "source"])
             stash_req1 = self.get_setting_value(config, [obj, "stash_req"])

             done = False
             for val in reqs_done:
               if stash_req1 == val:
                 done = True

             if done == True:
               continue
             else:

               reqs_group = []
               reqs_group.append(stash_req1)
               reqs_done.append(stash_req1)
              
               for obi in config.get_value():
               
               # Loop through items name lists again to construct request groupings
                 if re.search(r'namelist:items', obi):
                 
                   ancfile2   = self.get_setting_value(config, [obi, "ancilfilename"])
                   stash_req2 = self.get_setting_value(config, [obi, "stash_req"])

                   if ancfile2 == ancfile1 and stash_req2 != stash_req1:
                  
                     reqs_group.append(stash_req2)
                     reqs_done.append(stash_req2)
                     # Construct list of ancil requests to remove
                     rem_list.append(obi)
                  
                     req_list = ",".join(reqs_group)
                     
               if len(req_list) > 0:
               
               # Add grouped requests
               
                 self.remove_setting(config, [obj, "stash_req"])
                 self.add_setting   (config, [obj, "stash_req"], req_list)
                 req_list = ''
               
        # Remove ancil request no longer required
        rems = ",".join(rem_list)
        for nml in rem_list:
          self.remove_setting(config, [nml])
                                    
        return config, self.reports


class vn91_t5944(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5944 by Erica Neininger."""

    BEFORE_TAG = "vn9.1_t6151"
    AFTER_TAG = "vn9.1_t5944"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        autopp = self.get_setting_value(config, ["namelist:nlstcall_qs","autopp"])
        self.remove_setting(config, ["namelist:nlstcall_qs","autopp"])
        if autopp == "'Y'":
              self.add_setting(config, ["namelist:nlstcall_qs","l_postp"], '.true.')
        else:
              self.add_setting(config, ["namelist:nlstcall_qs","l_postp"], '.false.')

        self.remove_setting(config, ["namelist:nlstcall_qs","arch_sys"])
        self.remove_setting(config, ["namelist:nlstcall_qs","gddel"])
        self.remove_setting(config, ["namelist:nlstcall_qs","gpdel"])
        self.remove_setting(config, ["namelist:nlstcall_qs","gcmdel"])
        self.remove_setting(config, ["namelist:nlstcgen","archdump_freqim"])
        self.remove_setting(config, ["namelist:nlstcgen","archdump_offsetim"])
        self.remove_setting(config, ["namelist:nlstcgen","archppselim"])
        self.remove_setting(config, ["namelist:nlstcgen","jobrel_stepim"])
        self.remove_setting(config, ["namelist:nlstcgen","jobrel_offsetim"])
        return config, self.reports


class vn91_t6132(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6132 by Harry Shepherd (hshep)."""

    BEFORE_TAG = "vn9.1_t5944"
    AFTER_TAG = "vn9.1_t6132"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config,
                         ['namelist:nlstcgen', 'greg_monthly_dump'],
                         '.false.')
        self.add_setting(config,
                         ['namelist:nlstcgen', 'end_of_run_dump'],
                         '.false.')
        return config, self.reports


class vn91_t5704(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5704 by Glenn Greed."""

    BEFORE_TAG = "vn9.1_t6132"
    AFTER_TAG = "vn9.1_t5704"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.remove_setting(config, ["namelist:nlstcall","model_assim_mode"])
        self.remove_setting(config, ["namelist:nlstcall","run_assim_mode"])
        return config, self.reports

class vn91_t6256(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6256 by Joe Mancell."""

    BEFORE_TAG = "vn9.1_t5704"
    AFTER_TAG = "vn9.1_t6256"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        l_interp_input_only = self.get_setting_value(config, 
                                                     ["namelist:recon",
                                                      "l_interp_input_only"])
        var_recon = self.get_setting_value(config, ["namelist:recon",
                                                    "var_recon"])
        # If cannot find variable set to sensible default
        if l_interp_input_only is None:
            self.add_setting(config, ["namelist:recon", "select_input_fields"],
                             "0")
        elif re.search('true', l_interp_input_only, flags=re.IGNORECASE):
            self.add_setting(config, ["namelist:recon", "select_input_fields"],
                             "1")
        else:
            self.add_setting(config, ["namelist:recon", "select_input_fields"],
                             "0")
        self.remove_setting(config, ["namelist:recon", "l_interp_input_only"])

        # Add new flag to separate off some of the var_recon functionality
        if var_recon is None:
            self.add_setting(config, ["namelist:recon", "l_validity_lookup_u"],
                             ".false.")
        elif re.search('true', var_recon, flags=re.IGNORECASE):
            self.add_setting(config, ["namelist:recon", "l_validity_lookup_u"],
                             ".true.")
        else:
            self.add_setting(config, ["namelist:recon", "l_validity_lookup_u"],
                             ".false.")

        return config, self.reports


class vn91_t6463(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6463 by hshep."""

    BEFORE_TAG = "vn9.1_t6256"
    AFTER_TAG = "vn9.1_t6463"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here

        steps_per_period = self.get_setting_value( config, \
                            ['namelist:nlstcgen', 'steps_per_periodim'])
        secs_per_period = self.get_setting_value( config, \
                           ['namelist:nlstcgen', 'secs_per_periodim'])

        # return T3HRMN, T3HRLY, T3HRLY2 time profiles to how they were
        # at version 8.6. Make the T3HR profile more generic. It must
        # have a value of ifre=3hours, istr=3hours-1timestep, but to
        # accomplish this must have units of timesteps.
        for obj in config.get_value():
            if re.search(r'namelist:time', obj):
                i_name = self.get_setting_value(config,
                                                [obj,'tim_name'])
                #get rid of the quotes around the profile name
                i_name = i_name[1:-1]
                # recall below that a value of 2 when refering to time units
                # signifies hours, and a value of 1 signifies timesteps.
                if i_name == 'T3HRMN':
                    self.change_setting_value(config,
                                              [obj,'unt3'],
                                              '2')
                    self.change_setting_value(config,
                                              [obj,'ifre'],
                                              '3')
                    self.change_setting_value(config,
                                              [obj,'istr'],
                                              '3')
                elif i_name == 'T3HRLY':
                    self.change_setting_value(config,
                                              [obj,'unt3'],
                                              '2')
                    self.change_setting_value(config,
                                              [obj,'ifre'],
                                              '3')
                    self.change_setting_value(config,
                                              [obj,'istr'],
                                              '3')
                elif i_name == 'T3HRLY2':
                    self.change_setting_value(config,
                                              [obj,'unt3'],
                                              '2')
                    self.change_setting_value(config,
                                              [obj,'ifre'],
                                              '3')
                    self.change_setting_value(config,
                                              [obj,'istr'],
                                              '2')
                elif i_name == 'T3HR':
                    steps_per_day = int(steps_per_period) * \
                        (86400/int(secs_per_period))
                    steps_per_hour = steps_per_day / 24
                    self.change_setting_value(config,
                                              [obj,'ifre'],
                                              str(steps_per_hour*3))
                    self.change_setting_value(config,
                                              [obj,'istr'],
                                              str((steps_per_hour*3)-1))
                    self.change_setting_value(config,
                                              [obj,'unt3'],
                                              '1')
                else:
                    pass
        return config, self.reports


class vn91_t6395(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6395 by James Manners."""

    BEFORE_TAG = "vn9.1_t6463"
    AFTER_TAG = "vn9.1_t6395"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Get current settings
        l_eqt = self.get_setting_value(config, ["namelist:run_radiation","l_eqt"])
        sc = self.get_setting_value(config, ["namelist:run_radiation","sc"])
        # Remove retired items
        self.remove_setting(config, ["namelist:run_radiation","l_eqt"])
        self.remove_setting(config, ["namelist:run_radiation","sc"])
        self.remove_setting(config, ["namelist:planet_constants","l_planet_locked"])
        self.remove_setting(config, ["namelist:planet_constants","planet_lock_lon"])
        self.remove_setting(config, ["namelist:planet_constants","planet_solday"])
        self.remove_setting(config, ["namelist:planet_constants","planet_sidday"])
        # Add new items
        if l_eqt:
            self.add_setting(config, ["namelist:planet_constants","i_eqt"], "1")
        else:
            self.add_setting(config, ["namelist:planet_constants","i_eqt"], "0")
        self.add_setting(config, ["namelist:planet_constants","sc"], sc)
        self.add_setting(config, ["namelist:planet_constants","l_planet_intrinsic_flux"], ".false.")
        self.add_setting(config, ["namelist:planet_constants","planet_t_intrinsic"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","l_planet_aerosol"], ".false.")
        self.add_setting(config, ["namelist:planet_constants","planet_aerosol_mmr"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","l_fix_solang"], ".false.")
        self.add_setting(config, ["namelist:planet_constants","solar_zenith_angle"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","solar_azimuth_angle"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_ha"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_dha"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","sclht"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","lapse"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","omega"], "0.0")
        self.add_setting(config, ["namelist:r2lwclnl","l_na_lw"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_na_lw2"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_k_lw"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_k_lw2"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_na_sw"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_na_sw2"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_k_sw"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_k_sw2"], ".false.")
        return config, self.reports


class vn91_t6498(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6498 by Ian Boutle."""

    BEFORE_TAG = "vn9.1_t6395"
    AFTER_TAG = "vn9.1_t6498"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:temp_fixes", "l_iau_pc2check"], ".false.")
        return config, self.reports

class vn91_t6091(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6091 by Malcolm Brooks."""

    BEFORE_TAG = "vn9.1_t6498"
    AFTER_TAG = "vn9.1_t6091"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:temp_fixes", "l_fail_p_layers_inconsis"], ".false.")
        return config, self.reports

class vn91_t6514(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6514 by Malcolm Brooks."""

    BEFORE_TAG = "vn9.1_t6091"
    AFTER_TAG = "vn9.1_t6514"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:temp_fixes", "l_pc2_homog_turb_q_neg"], ".false.")
        return config, self.reports

class vn91_t6194(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6194 by Glenn Greed."""

    BEFORE_TAG = "vn9.1_t6514"
    AFTER_TAG = "vn9.1_t6194"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration, remove some temp logicals."""
        self.remove_setting(config,["namelist:temp_fixes", 
                                    "l_exclude_fix_theta_source"])
        self.remove_setting(config,["namelist:temp_fixes", 
                                    "l_p2t_weight_fix"])
        self.remove_setting(config,["namelist:temp_fixes", 
                                    "l_riverinland_fix"])
        self.remove_setting(config,["namelist:temp_fixes", 
                                    "l_error_ancil_struct"])
        self.add_setting( config,["namelist:temp_fixes", 
                                    "l_ignore_error_ancil_struct"],
                                    value = '.false.')  
        self.remove_setting(config,["namelist:temp_fixes", 
                                    "l_use_old_mass_fix"])
        return config, self.reports

class vn91_t6339(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6339 by Maggie Hendry."""

    BEFORE_TAG = "vn9.1_t6194"
    AFTER_TAG = "vn9.1_t6339"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Compulsory l_rcf_init_flexi set to false
        self.add_setting(config, ["namelist:recon", "l_rcf_init_flexi"], ".false.")
        return config, self.reports

class vn91_t6297(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6297 by Steve Wardle."""

    BEFORE_TAG = "vn9.1_t6339"
    AFTER_TAG = "vn9.1_t6297"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Entries from nlcfiles which have been removed from the code, and
        # aren't output stream filenames or possible old macro names
        delete_entries = ["mctl"    , "ictl"    , "hkfile"  , "ppxref"  ,
                          "config"  , "stashctl", "namelist", "output"  ,
                          "output2" , "xhist"   , "thist"   , "icecalve",
                          "cache1"  , "cache2"  , "aswap"   , "oswap"   ,
                          "arestart", "aopsum1" , "aopsum2" , "ssu"     ,
                          "ozone"   , "smcsnowd", "dsoiltmp", "soiltype",
                          "gentype" , "sstin"   , "sicein"  , "perturb" ,
                          "curntin" , "mask"    , "oinitial", "ostart"  ,
                          "orestart", "aopstmp1", "aopstmp2", "ocnanl"  ,
                          "atracer" , "otracer" , "wfin"    , "hfluxin" ,
                          "pmein"   , "icefin"  , "airtmp"  , "gensea"  ,
                          "fluxcorr", "swspectd", "bas_ind" , "slabhcon",
                          "dustsoil", "biomass" , "rivstor" , "rivchan" ,
                          "river2a" , "lwspectd", "dmsconc" , "orog"    ,
                          "olabcin" , "ocndepth", "murkfile", "sulpemis",
                          "usrancil", "usrmulti", "ousrancl", "ousrmult",
                          "so2natem", "chemoxid", "aerofcg" , "co2emits",
                          "tppsozon", "landfrac", "wlabcou1", "wlabcou2",
                          "wlabcou3", "wlabcou4", "ocffemis", "horzgrid",
                          "surfemis", "aircrems", "stratems", "extraems",
                          "radonems", "fracinit", "veginit" , "disturb" ,
                          "cached"  , "sootemis", "alabcou1", "alabcou2",
                          "alabcou3", "alabcou4", "alabcou5", "alabcou6",
                          "alabcou7", "alabcou8", "cariolo3", "arclbiog",
                          "arclbiom", "arclblck", "arclsslt", "arclsulp",
                          "arcldust", "arclocff", "arcldlta", "topmean" ,
                          "topstdev", "ukcastrd", "ukcasto3"]

        for entry in delete_entries:
              self.remove_setting(config, ["namelist:nlcfiles", entry])

        # Mappings for either the existing streams or possible old macros
        # Format is entry_name, old unit, old stream letter
        unit_mapping = {"60"  : "pp0", 
                        "61"  : "pp1", 
                        "62"  : "pp2", 
                        "63"  : "pp3", 
                        "64"  : "pp4", 
                        "65"  : "pp5", 
                        "66"  : "pp6", 
                        "67"  : "pp7", 
                        "68"  : "pp8", 
                        "69"  : "pp9", 
                        "150" : "ppvar", 
                        "151" : "pp10", 
                        "164" : "ppmbc", 
                        "27"  : "aomean", 
                        "81"  : "surgeou1", 
                        "82"  : "surgeout", 
                        "100" : "foamout1", 
                        "101" : "foamout2", 
                        "85"  : "wfout", 
                        "102" : "cxbkgerr", 
                        "83"  : "ppscreen", 
                        "84"  : "ppsmc", 
                        "86"  : "uarsout1", 
                        "87"  : "uarsout2", 
                        "103" : "rfmout", 
                        "89"  : "mosout", 
                        "93"  : "curntout", 
                        "94"  : "flxcrout", 
                        "88"  : "icefout", 
                        "92"  : "siceout", 
                        "91"  : "sstout"} 

        # Iterate through the use namelists in the app, recording which files
        # from the above list are in use, and update the entry to use the old
        # name as the new file id
        used_ids = []
        for section in config.value.keys():
              if not re.match("namelist:use\(.*\)$", section):
                    continue

              iunt = self.get_setting_value(config, [section, "iunt"])
              locn = self.get_setting_value(config, [section, "locn"])
              self.remove_setting(config, [section, "iunt"])

              # Locn determines the meaning of the unit number variable
              if locn == "3":
                    # For actual output stream files, provide the new id string
                    # with the correct mapping to the file unit given above

                    # This block was added post vn10.1, when a problem was
                    # discovered with this macro conflicting against another
                    # macro - this block allows the macro to be run again on
                    # an app which contains problems, by assuming that if iunt
                    # is not found, the app was already upgraded
                    if iunt is None:
                        # Get the id which was set by the previous upgrade
                        file_id = self.get_setting_value(config,
                                                         [section, "file_id"])
                        if file_id is not None:
                            file_id = file_id.strip("'")
                            # Get the filename currently in the app
                            current_fname = self.get_setting_value(config,
                                ["namelist:nlstcall_pp({0:s})".format(file_id),
                                 "filename"])
                            if current_fname == "":
                                # Overwrite it with the name from nlcfiles
                                # (only if it's blank, as the user might have
                                # manually fixed it)
                                new_fname = self.get_setting_value(config,
                                    ["namelist:nlcfiles", file_id])
                                if new_fname is None:
                                    new_fname = ""
                                self.change_setting_value(config,
                                ["namelist:nlstcall_pp({0:s})".format(file_id),
                                 "filename"], new_fname)
                        continue
                    
                    file_id = "'"+unit_mapping[iunt]+"'"
                    self.add_setting(config, [section, "file_id"], file_id)
                    self.add_setting(config, [section, "macrotag"], "0")
                    used_ids.append(file_id.strip("'"))
              else:
                    # Other uses are providing a stash macrotag to be referred

                    # to later.  These must stay as integers but the name has
                    # been changed to improve clarity
                    self.add_setting(config, [section, "file_id"], "''")
                    self.add_setting(config, [section, "macrotag"], iunt)



        # Copying the filenames and ids from nlcfiles to the existing nlstcall_pp
        # output streams 
        map_files = [("pp0",   "60",  "a"),
                     ("pp1",   "61",  "b"),
                     ("pp2",   "62",  "c"),
                     ("pp3",   "63",  "d"),
                     ("pp4",   "64",  "e"),
                     ("pp5",   "65",  "f"),
                     ("pp6",   "66",  "g"),
                     ("pp7",   "67",  "h"),
                     ("pp8",   "68",  "i"),
                     ("pp9",   "69",  "j"),
                     ("ppvar", "150", "a"),
                     ("pp10",  "151", "k"),
                     ("ppmbc", "164", "b")]

        # New names for the settings in the namelist (note the absence of i_ppa, 
        # since the archiving flag isn't used since #5944 it can be retired)
        map_settings = [("i_ppie", "reinit_end"),
                        ("i_ppis", "reinit_start"),
                        ("i_ppiu", "reinit_unit"),
                        ("i_ppos", "reserved_headers"),
                        ("i_ppx",  "packing"),
                        ("i_ppif", "reinit_step"),
                        ("i_ppi",  "l_reinit")]

        for nlcfiles_entry, namelist_index, reinit_letter in map_files:
              old_namelist = "namelist:nlstcall_pp({0:s})".format(namelist_index)
              new_namelist = "namelist:nlstcall_pp({0:s})".format(nlcfiles_entry)
              if nlcfiles_entry in used_ids:
                    # If the stream is in use by the app, get the filename from nlcfiles and
                    # copy it into a new nlstcall_pp namelist, named according to the id
                    self.add_setting(config, [new_namelist])
                    
                    filename = (
                          self.get_setting_value(config, ["namelist:nlcfiles", nlcfiles_entry]))
                    self.add_setting(config, [new_namelist, "filename"], filename)
                    self.add_setting(config, [new_namelist, "file_id"], "'"+nlcfiles_entry+"'")
                    self.add_setting(config, [new_namelist, "reinit_letter"], "'"+reinit_letter+"'")

                    # Rename the settings as per the code change
                    for old_setting, new_setting in map_settings:
                          value = self.get_setting_value(config, [old_namelist, old_setting])
                          self.add_setting(config, [new_namelist, new_setting], value)

              # Now remove the old namelist
              self.remove_setting(config, [old_namelist])
                    
              # Either way, cleanup the entry from nlcfiles
              self.remove_setting(config, ["namelist:nlcfiles", nlcfiles_entry])

        # To correspond with a line removed from the code, note that the ppvar stream
        # always had its reserved header setting overwritten to 8192.  So we'll change
        # the input value to that to preserve behaviour
        ppvar_nl = "namelist:nlstcall_pp(ppvar)"
        self.change_setting_value(config, [ppvar_nl, "reserved_headers"], "8192")

        # Adding output streams which were missing - these are mostly carried across
        # from what used to be macros in the UMUI
        new_streams = [("aomean",   ""),
                       ("surgeou1", ""),
                       ("surgeout", ""),
                       ("foamout1", "l"),
                       ("foamout2", "m"),
                       ("wfout",    ""),
                       ("cxbkgerr", "n"),
                       ("ppscreen", ""),
                       ("ppsmc",    ""),
                       ("uarsout1", ""),
                       ("uarsout2", ""),
                       ("rfmout",   ""),
                       ("upmean",   ""), 
                       ("mosout",   ""),
                       ("curntout", ""),
                       ("flxcrout", ""),
                       ("icefout",  ""),
                       ("siceout",  ""),
                       ("sstout",   "")]

        for file_id, reinit_letter in new_streams:
              if file_id in used_ids:
                    new_namelist = "namelist:nlstcall_pp({0:s})".format(file_id)
                    self.add_setting(config, [new_namelist, "file_id"], "'"+file_id+"'")
                    self.add_setting(config, [new_namelist, "filename"], "'$DATAW/$RUNID."+file_id+"'")
                    self.add_setting(config, [new_namelist, "l_reinit"], ".false.")
                    self.add_setting(config, [new_namelist, "reinit_end"], "0")
                    self.add_setting(config, [new_namelist, "reinit_step"], "0")
                    self.add_setting(config, [new_namelist, "reinit_start"], "0")
                    self.add_setting(config, [new_namelist, "reinit_unit"], "1")
                    self.add_setting(config, [new_namelist, "packing"], "1")
                    self.add_setting(config, [new_namelist, "reinit_letter"], "'"+reinit_letter+"'")

                    # As with ppvar above a couple of these had hard-coded values for
                    # the override which have been removed from the code, so make sure
                    # they get the correct values in the inputs
                    if file_id == "cxbkgerr":
                          self.add_setting(config, [new_namelist, "reserved_headers"], "12288")
                    elif file_id == "rfmout":
                          self.add_setting(config, [new_namelist, "reserved_headers"], "8192")
                    else:
                          self.add_setting(config, [new_namelist, "reserved_headers"], "0")

              self.remove_setting(config, ["namelist:nlcfiles", file_id])

        return config, self.reports


class vn91_t4251(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4251 by Warren Tennant."""

    BEFORE_TAG = "vn9.1_t6297"
    AFTER_TAG = "vn9.1_t4251"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:iau_nl", "l_iau_usetsoilperts"], ".false.")
        return config, self.reports


class vn91_t6480(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6480 by Cyril Morcrette (frme)"""

    BEFORE_TAG = "vn9.1_t4251"
    AFTER_TAG = "vn9.1_t6480"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here

        l_filter_cloud = self.get_setting_value(config, ["namelist:run_cloud", "l_filter_cloud"])
        if l_filter_cloud == ".true.":
            self.add_setting(config, ["namelist:run_cloud", "l_od_cld_filter"], ".true.")
        else:
            self.add_setting(config, ["namelist:run_cloud", "l_od_cld_filter"], ".false.")

        self.remove_setting(config, ["namelist:run_cloud", "l_filter_cloud"])
        self.add_setting(config, ["namelist:run_cloud", "l_ceil_cld_filter"], ".false.")
        self.add_setting(config, ["namelist:run_cloud", "l_sharpen_cbh_diags"], ".false.")

        return config, self.reports


class vn91_t6483(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6483 by Adrian Lock."""

    BEFORE_TAG = "vn9.1_t6480"
    AFTER_TAG = "vn9.1_t6483"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        l_subfilter_blend = self.get_setting_value(config, ["namelist:run_diffusion", "l_subfilter_blend"])
        if l_subfilter_blend == ".true.":
          self.add_setting(config, ["namelist:run_bl", "blending_option"], "1")
        else:
          self.add_setting(config, ["namelist:run_bl", "blending_option"], "0")
        self.remove_setting(config, ["namelist:run_diffusion", "l_subfilter_blend"])
        self.add_setting(config, ["namelist:run_bl", "keep_ri_fa"], "1")
        return config, self.reports


class vn91_t6357(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6357 by David Walters."""

    BEFORE_TAG = "vn9.1_t6483"
    AFTER_TAG = "vn9.1_t6357"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, [ "namelist:temp_fixes", "l_methox_fix" ],
                         str(".false."))
        return config, self.reports


class vn91_t5830(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5830 by Paul Cresswell."""

    BEFORE_TAG = "vn9.1_t6357"
    AFTER_TAG = "vn9.1_t5830"

    def upgrade(self, config, meta_config=None):
        """Add Dr Hook variables."""
        self.add_setting(config, ["env", "DR_HOOK"], "0")
        self.add_setting(config, ["env", "DR_HOOK_OPT"], "wallprof,self")
        self.add_setting(config, ["env", "DR_HOOK_PROFILE"],
                         "$CYLC_TASK_WORK_DIR/drhook.prof.%d")
        self.add_setting(config, ["env", "DR_HOOK_PROFILE_PROC"], "-1")
        self.add_setting(config, ["env", "DR_HOOK_PROFILE_LIMIT"], "-10.0")
        self.add_setting(config, ["env", "DR_HOOK_IGNORE_SIGNALS"], "0")
        self.add_setting(config, ["env", "DR_HOOK_CATCH_SIGNALS"], "0")

        return config, self.reports


class vn91_t6470(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6470 by J. M. Edwards."""

    BEFORE_TAG = "vn9.1_t5830"
    AFTER_TAG = "vn9.1_t6470"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        iseaz0t=self.get_setting_value(config,
            ["namelist:jules_sea_seaice","iseaz0t"])
        self.remove_setting(config,
            ["namelist:jules_sea_seaice","iseaz0t"])
#       Most, if not all, jobs should have iseaz0t=1 by now, but allow
#       for leagcy cases with 0.
        if iseaz0t == '1':
            self.add_setting(config,
                ["namelist:jules_sea_seaice","iseasurfalg"], "1")
        else:
            self.add_setting(config,
                ["namelist:jules_sea_seaice","iseasurfalg"], "0")

#       Add defaults based on COARE3.0.
        self.add_setting(config,
                ["namelist:jules_sea_seaice","l_10m_neut"], ".true.")
        self.add_setting(config,
                ["namelist:jules_sea_seaice","u10_min_coare"], "10")
        self.add_setting(config,
                ["namelist:jules_sea_seaice","u10_max_coare"], "18")
        self.add_setting(config,
                ["namelist:jules_sea_seaice","a_chrn_coare"], "0.000875")
        self.add_setting(config,
                ["namelist:jules_sea_seaice","b_chrn_coare"], "0.0022500")

        return config, self.reports

class vn91_t5503(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5503 by Gerd Folberth."""

    BEFORE_TAG = "vn9.1_t6470"
    AFTER_TAG = "vn9.1_t5503"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Set default value of l_bvoc_emis from false to false
        self.add_setting(config, ["namelist:jules_vegetation", "l_bvoc_emis"], ".false.")

        RMDI = str(-2**30)
        npft = int(self.get_setting_value(config, ["namelist:jules_surface_types", "npft"]))

        # Add sensible values for the first (up to) 5 pfts and missing data for the rest
        # This will trigger a sensible runtime error in standalone JULES to allow the user
        # to fix the problem
        ci_st_io = ["33.46", "33.46", "34.26", "29.98", "34.26"] + ([RMDI] * (npft - 5))
        self.add_setting(config, ["namelist:jules_pftparm", "ci_st_io"], "  ".join(ci_st_io[:npft]))

        gpp_st_io = ["1.29E-07", "2.58E-08", "2.07E-07", "3.42E-07", "1.68E-007"] + ([RMDI] * (npft - 5))
        self.add_setting(config, ["namelist:jules_pftparm", "gpp_st_io"], "  ".join(gpp_st_io[:npft]))

        self.add_setting(config, ["namelist:run_ukca", "l_ukca_ibvoc"], ".false.")
        
        return config, self.reports


class vn91_t6009(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6009 by Ben Johnson."""

    BEFORE_TAG = "vn9.1_t5503"
    AFTER_TAG = "vn9.1_t6009"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_ukca", "l_ukca_scale_biom_aer_ems"], ".false.")
        self.add_setting(config, ["namelist:run_ukca", "biom_aer_ems_scaling"], "1.0")
        self.add_setting(config, ["namelist:run_ukca", "l_ukca_scale_soa_yield"], ".false.")
        self.add_setting(config, ["namelist:run_ukca", "soa_yield_scaling"], "1.0")
        return config, self.reports

class vn91_t6277(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6277 by Jonathan Wilkinson/Kalli Furtado."""

    BEFORE_TAG = "vn9.1_t6009"
    AFTER_TAG = "vn9.1_t6277"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_precip", "l_subgrid_qcl_mp"], ".false.")
        self.add_setting(config, ["namelist:run_precip", "l_subgrid_cfl_mp_by_erosion"], ".false.")
        self.add_setting(config, ["namelist:run_precip", "l_mixed_phase_t_limit"], ".false.")
        self.add_setting(config, ["namelist:run_precip", "nbins_mp"], "100")
        self.add_setting(config, ["namelist:run_precip", "mp_t_limit"], "273.14")
        self.add_setting(config, ["namelist:run_precip", "mp_dz_scal"], "1.0")
        self.add_setting(config, ["namelist:run_precip", "mp_tau_d_lim"], "1200.0")
        self.add_setting(config, ["namelist:run_precip", "mp_czero"], "10.0")
        return config, self.reports

class vn91_t6507(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn9.1_t6277"
    AFTER_TAG = "vn9.1_t6507"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a JULES runtime app configuration."""
        
        # ADD the trait parameters to pft_params namelist

	self.add_setting(config, ["namelist:jules_pftparm", "kn_io"],   "5*0.78")
	self.add_setting(config, ["namelist:jules_pftparm", "q10_leaf_io"], "5*2.00")
	self.add_setting(config, ["namelist:jules_pftparm", "lma_io"], "0.0824,0.2263,0.0498,0.1370,0.0695")
	self.add_setting(config, ["namelist:jules_pftparm", "nmass_io"], "0.0210,0.0115,0.0219,0.0131,0.0219")
	self.add_setting(config, ["namelist:jules_pftparm", "vint_io"],"5.73,6.32,6.42,0.00,14.71")
	self.add_setting(config, ["namelist:jules_pftparm", "vsl_io"], "29.81,18.15,40.96,10.24,23.15")

	self.add_setting(config, ["namelist:jules_vegetation", "l_ht_compete"], ".false.")
	self.add_setting(config, ["namelist:jules_vegetation", "l_trait_phys"], ".false.")
	self.add_setting(config, ["namelist:jules_vegetation", "l_landuse"], ".false.")

	self.remove_setting(config, ["namelist:jules_surface", "q10_leaf"])

        return config, self.reports

class vn91_t6403(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6403 by Paul Cresswell, Carol Halliwell
       and Chris Smith."""

    BEFORE_TAG = "vn9.1_t6507"
    AFTER_TAG = "vn9.1_t6403"

    def upgrade(self, config, meta_config=None):
        """Import the idealise namelist into a UM app."""

        # Only get the namelist if it is active
        # (Lots of obsolete files hanging around in jobs and apps)
        idealise = self.get_setting_value(config, 
                   ['namelist:nlcfiles', 'idealise'], no_ignore=True)
        self.remove_setting(config, ['namelist:nlcfiles', 'idealise'])

        # Set up header for new file:
        self.add_setting(config, ["file:IDEALISE"])
        self.add_setting(config, ["file:IDEALISE", "source"],
                                  "(namelist:idealise)")

        if idealise and idealise != "'unset'":
            # Strip quotes:
            idealise = idealise.strip("'").strip('"')
            # Try and get the file:
            if os.path.isfile(idealise):

                # Create a tempfile to store roseified namelist:
                file = tempfile.NamedTemporaryFile()

                # Dump app data to file:
                rose.namelist_dump.namelist_dump([idealise], file.name, 'lower')

                if os.path.isfile(file.name):
                    for line in file:
                        if line.strip() == '[namelist:idealise]':
                            break
                    line = line.strip().strip('[').strip(']')
                    self.add_setting(config, [line])

                    for line in file:
                        line = line.strip().split('=', 1)
                        self.add_setting(config, 
                             ['namelist:idealise',line[0]],line[1])

                else:
                    raise UpgradeError (
                    'Cannot read temporary file for importing idealise namelist.')

            else:
                raise UpgradeError (
                'Cannot find the idealised namelist file.\n' + 
                'Please ensure the file is available on the local filesystem\n'+
                'and then edit namelist:nlcfiles=idealise with the full filepath.\n' +
                'If the namelist settings are not required you may set\n' +
                '"idealise=\'unset\'" instead.')

        # Fix weirdness from namelists stored as: item(1) = ..., item(2) = ...,
        # where the add_setting function will create a duplicate entry.
        # We expect this behaviour to be corrected in a later Rose release.
        # We can't cater for every scenario here (e.g. item(1,2) = ..., 
        # item(3,4) = ...) but this particular combination is known to occur.
        surf = [None, None, None, None]
        surf[0] = self.get_setting_value(config, 
                ["namelist:idealise", "idlsurffluxseaparams(1)"])
        surf[1] = self.get_setting_value(config, 
                ["namelist:idealise", "idlsurffluxseaparams(2)"])
        surf[2] = self.get_setting_value(config,
                ["namelist:idealise", "idlsurffluxseaparams(3)"])
        surf[3] = self.get_setting_value(config,
                ["namelist:idealise", "idlsurffluxseaparams(4)"])

        # Ensure the old version is gone...
        self.remove_setting(config, 
            ["namelist:idealise", "idlsurffluxseaparams(1)"])
        self.remove_setting(config,
            ["namelist:idealise", "idlsurffluxseaparams(2)"])
        self.remove_setting(config,
            ["namelist:idealise", "idlsurffluxseaparams(3)"])
        self.remove_setting(config,
            ["namelist:idealise", "idlsurffluxseaparams(4)"])

        # ...and if it was found like that, replace with the single-list form.
        if any(surf):
            for s in surf:
                if s is None:
                    s = 0.0
            self.add_setting(config, 
                  ["namelist:idealise", "idlsurffluxseaparams"],
                  surf[0]+","+surf[1]+","+surf[2]+","+surf[3])


        # Ensure every item is present with a default:
        self.add_setting(config, ["namelist:idealise", "aa_jet_a"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "aa_jet_m"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "aa_jet_n"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "aa_jet_u0"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "angular_velocity"],
                         "7.292116E-5")
        self.add_setting(config, ["namelist:idealise", "b_const"],
                         "2")
        self.add_setting(config, ["namelist:idealise", "base_frictional_timescale"],
                         "1.1574e-5")
        self.add_setting(config, ["namelist:idealise", "base_xi1"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "base_xi2"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "big_factor"],
                         "1.0")
        self.add_setting(config, ["namelist:idealise", "big_layers"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "brunt_vaisala"],
                         "0.01")
        self.add_setting(config, ["namelist:idealise", "cool_rate"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "delta_xi1"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "delta_xi1_h"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "delta_xi1_l"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "delta_xi2"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "delta_xi2_h"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "delta_xi2_l"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "diffuse_dse_ap1"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "diffuse_dse_ap2"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "diffuse_q_ap1"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "diffuse_q_ap2"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "diffuse_u_ap1"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "diffuse_u_ap2"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "diffuse_v_ap1"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "diffuse_v_ap2"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "dmptim"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "do_evap_ap1"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_evap_ap2"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_limit_surf_theta"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_precip_ap1"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_precip_ap2"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_radiation_lw_ap1"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_radiation_lw_ap2"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_radiation_sw_ap1"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_radiation_sw_ap2"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_surface_drag_ap1"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_surface_drag_ap2"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_surface_flux_ap1"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "do_surface_flux_ap2"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "dtheta_dz1"],
                         "0.01,0.01,0.01")
        self.add_setting(config, ["namelist:idealise", "eccentricity"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "f_plane"],
                         "-90.0")
        self.add_setting(config, ["namelist:idealise", "ff_plane"],
                         "-90.0")
        self.add_setting(config, ["namelist:idealise", "first_constant_r_rho_level_new"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "first_theta_height"],
                         "10.0")
        self.add_setting(config, ["namelist:idealise", "grid_flat"],
                         "3")
        self.add_setting(config, ["namelist:idealise", "grid_np_lat"],
                         "90.0")
        self.add_setting(config, ["namelist:idealise", "grid_np_lon"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "grid_number"],
                         "10")
        self.add_setting(config, ["namelist:idealise", "grow_steps"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "h_o"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "half_width_x"],
                         "2500000.0")
        self.add_setting(config, ["namelist:idealise", "half_width_y"],
                         "2500000.0")
        self.add_setting(config, ["namelist:idealise", "hdmp"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "height_domain"],
                         "10000.0")
        self.add_setting(config, ["namelist:idealise", "height_dz1"],
                         "0.0,0.0")
        self.add_setting(config, ["namelist:idealise", "height_u_in"],
                         "0.0,0.0,0.0")
        self.add_setting(config, ["namelist:idealise", "hf"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "idl_bubble_depth"],
                         "1000.0,1000.0,1000.0")
        self.add_setting(config, ["namelist:idealise", "idl_bubble_height"],
                         "1000.0,1000.0,1000.0")
        self.add_setting(config, ["namelist:idealise", "idl_bubble_max"],
                         "1.0,1.0,1.0")
        self.add_setting(config, ["namelist:idealise", "idl_bubble_option"],
                         "0,0,0")
        self.add_setting(config, ["namelist:idealise", "idl_bubble_width"],
                         "1000.0,1000.0,1000.0")
        self.add_setting(config, ["namelist:idealise", "idl_bubble_xoffset"],
                         "0.5,0.5,0.5")
        self.add_setting(config, ["namelist:idealise", "idl_bubble_yoffset"],
                         "0.5,0.5,0.5")
        self.add_setting(config, ["namelist:idealise", "idl_interp_option"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "idlsurffluxseaoption"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "idlsurffluxseaparams"],
                         "0.0,0.0,0.0,0.0")
        self.add_setting(config, ["namelist:idealise", "instability_diagnostics"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "k_const"],
                         "3")
        self.add_setting(config, ["namelist:idealise", "l_baro_inst"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_baro_perturbed"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_baroclinic"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_bomex"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_cartesian"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_code_test"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_const_grav"],
                         ".true.")
        self.add_setting(config, ["namelist:idealise", "l_constant_dz"],
                         ".true.")
        self.add_setting(config, ["namelist:idealise", "l_cyclone"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_damp"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_deep_baro_inst"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_fix_orog_hgt_lbc"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_fixed_lbcs"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_force_lbc"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_frierson"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_geo_for"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_heldsuarez"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_heldsuarez1_drag"],
                         ".true.")
        self.add_setting(config, ["namelist:idealise", "l_idealised_data"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_idl_bubble_saturate"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_initialise_data"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_perturb_correlate_time"],
                         ".true.")
        self.add_setting(config, ["namelist:idealise", "l_perturb_correlate_tq"],
                         ".true.")
        self.add_setting(config, ["namelist:idealise", "l_perturb_correlate_vert"],
                         ".true.")
        self.add_setting(config, ["namelist:idealise", "l_perturb_q"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_perturb_t"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_pforce"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_polar_wind_zero"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_pressure_balance"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_rotate_grid"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_rotate_winds"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_rotating"],
                         ".true.")
        self.add_setting(config, ["namelist:idealise", "l_sh_williamson"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_shallow"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_solid_body"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_spec_z0"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_trivial_trigs"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_vert_coriolis"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "l_wind_balance"],
                         ".false.")
        self.add_setting(config, ["namelist:idealise", "lambda_fraction"],
                         "0.5")
        self.add_setting(config, ["namelist:idealise", "mag"],
                         "1.0")
        self.add_setting(config, ["namelist:idealise", "mod_layers"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "newtonian_timescale"],
                         "3600.0")
        self.add_setting(config, ["namelist:idealise", "num_pforce_times"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "num_profile_data"],
                         "100")
        self.add_setting(config, ["namelist:idealise", "num_qforce_levels"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "num_qforce_times"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "num_tforce_levels"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "num_tforce_times"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "num_uvforce_levels"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "num_uvforce_times"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "num_uvprofile_data"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "nxi1l"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "nxi1v"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "nxi2l"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "nxi2v"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "orog_hgt_lbc"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "p_surface"],
                         "100000.0")
        self.add_setting(config, ["namelist:idealise", "p_surface_data"],
                         str( "100000.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "perturb_height"],
                         "0.0,0.0")
        self.add_setting(config, ["namelist:idealise", "perturb_magnitude_q"],
                         "0.5e-3")
        self.add_setting(config, ["namelist:idealise", "perturb_magnitude_t"],
                         "0.5")
        self.add_setting(config, ["namelist:idealise", "perturb_type"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "pforce_option"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "pforce_time_interval"],
                         "600.0")
        self.add_setting(config, ["namelist:idealise", "phi_fraction"],
                         "0.5")
        self.add_setting(config, ["namelist:idealise", "plat_size_x"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "plat_size_y"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "q1"],
                         "70.0")
        # Full array has 10,000 elements, this is enough for one time profile:
        self.add_setting(config, ["namelist:idealise", "qforce_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "qforce_option"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "qforce_time_interval"],
                         "600.0")
        self.add_setting(config, ["namelist:idealise", "qprofile_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "qprofile_number"],
                         "10")
        self.add_setting(config, ["namelist:idealise", "r_plane"],
                         "-90.0")
        self.add_setting(config, ["namelist:idealise", "roughlen_z0h"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "roughlen_z0m"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "suhe_fric"],
                         "2")
        self.add_setting(config, ["namelist:idealise", "suhe_newtonian_timescale_ka"],
                         "2.893e-7")
        self.add_setting(config, ["namelist:idealise", "suhe_newtonian_timescale_ks"],
                         "2.893e-6")
        self.add_setting(config, ["namelist:idealise", "suhe_pole_equ_deltat"],
                         "60.0")
        self.add_setting(config, ["namelist:idealise", "suhe_relax"],
                         "2")
        self.add_setting(config, ["namelist:idealise", "suhe_sigma_cutoff"],
                         "0.7")
        self.add_setting(config, ["namelist:idealise", "suhe_static_stab"],
                         "10.0")
        self.add_setting(config, ["namelist:idealise", "surface_type"],
                         "10")
        self.add_setting(config, ["namelist:idealise", "t0_e"],
                         "310.0")
        self.add_setting(config, ["namelist:idealise", "t0_p"],
                         "240.0")
        self.add_setting(config, ["namelist:idealise", "t_horizfn_data"],
                         "0.0,0.0,0.0,0.0")
        self.add_setting(config, ["namelist:idealise", "t_horizfn_number"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "t_surface"],
                         "280.0")
        # Full array has 10,000 elements, this is enough for one time profile:
        self.add_setting(config, ["namelist:idealise", "tforce_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "tforce_option"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "tforce_time_interval"],
                         "600.0")
        self.add_setting(config, ["namelist:idealise", "theta_surface"],
                         "280.0")
        self.add_setting(config, ["namelist:idealise", "thin_theta_height"],
                         "1.0")
        self.add_setting(config, ["namelist:idealise", "tprofile_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "tprofile_number"],
                         "10")
        self.add_setting(config, ["namelist:idealise", "transit_layers"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "tropics_deg"],
                         "30.0")
        self.add_setting(config, ["namelist:idealise", "tstep_plot_frequency"],
                         "1")
        self.add_setting(config, ["namelist:idealise", "tstep_plot_start"],
                         "-1")
        self.add_setting(config, ["namelist:idealise", "u_geo"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "u_in"],
                         "0.0,0.0,0.0,0.0")
        self.add_setting(config, ["namelist:idealise", "u_ramp_end"],
                         "-90.0")
        self.add_setting(config, ["namelist:idealise", "u_ramp_start"],
                         "0.0")
        # Full array has 10,000 elements, this is enough for one time profile:
        self.add_setting(config, ["namelist:idealise", "uforce_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "ujet_lat"],
                         "-90.0")
        self.add_setting(config, ["namelist:idealise", "ujet_width"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "uprofile_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "uv_horizfn_number"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "uvforce_option"],
                         "0")
        self.add_setting(config, ["namelist:idealise", "uvforce_time_interval"],
                         "600.0")
        self.add_setting(config, ["namelist:idealise", "uvprofile_number"],
                         "10")
        self.add_setting(config, ["namelist:idealise", "v_geo"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "v_in"],
                         "0.0,0.0,0.0,0.0")
        self.add_setting(config, ["namelist:idealise", "vert_grid_ratio"],
                         "1.0")
        # Full array has 10,000 elements, this is enough for one time profile:
        self.add_setting(config, ["namelist:idealise", "vforce_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "vprofile_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "witch_power"],
                         "1.5")
        self.add_setting(config, ["namelist:idealise", "z_qforce_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "z_tforce_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "z_uvforce_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "z_uvprofile_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "zdmp"],
                         "0.0")
        self.add_setting(config, ["namelist:idealise", "zprofile_data"],
                         str( "0.0," * 100).strip(','))
        self.add_setting(config, ["namelist:idealise", "zprofile_orog"],
                         "0.0")

        # Remove any items deleted from the namelist:
        self.remove_setting(config, ["namelist:idealise","l_balanced"])
        self.remove_setting(config, ["namelist:idealise","l_exact_profile"])
        self.remove_setting(config, ["namelist:idealise","l_force"])
        self.remove_setting(config, ["namelist:idealise","l_heldsuarez2_drag"])
        self.remove_setting(config, ["namelist:idealise","l_isothermal"])
        self.remove_setting(config, ["namelist:idealise","ring_height"])
        self.remove_setting(config, ["namelist:idealise","theta_pert"])
        self.remove_setting(config, ["namelist:idealise","trefer_number"])

        return config, self.reports


class vn91_t1389(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1389 (SRS) by Paul Cresswell."""

    BEFORE_TAG = "vn9.1_t6403"
    AFTER_TAG = "vn9.1_t1389"

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


class vn92_t6571(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn9.1_t1389"
    AFTER_TAG = "vn9.2"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "9.2")
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMSTR"])
        if stashmaster_path:
           stashmaster_path = re.sub("vn\d+\.\d+/ctldata",
                                     "vn9.2/ctldata",
                                      stashmaster_path)
           self.change_setting_value(config, ["env", "STASHMSTR"],
                                     stashmaster_path)
        return config, self.reports


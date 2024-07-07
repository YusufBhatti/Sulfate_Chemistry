import rose.upgrade
import os
import re
from rose.popen import RosePopener
import sys
import shlex
from collections import defaultdict


class UpgradeError(Exception): 

      """Exception created when an upgrade fails.""" 
           
      def __init__(self, msg): 
          self.msg = msg 
       
      def __repr__(self): 
          sys.tracebacklimit = 0 
          return self.msg 
           
      __str__ = __repr__ 


class vn90_t6006(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6006 by Glenn Greed."""

    BEFORE_TAG = "vn9.0"
    AFTER_TAG = "vn9.0_t6006"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        tim_check={'tim_name':"",
                   'ityp':'0',
                   'isam':'0',
                   'ioff':'0',
                   'intv':'0', 
                   'unt1':"'H '",
                   'unt2':"'H '",
                   'unt3':"'H '",
                   'iopt':'0',
                   'istr':'0',
                   'iend':'0',
                   'ifre':'0',
                   'itimes':'0',
                   'iser':'0',
                   'isdt':'0',
                   'iedt':'0'
                   }
        
        dom_check={'dom_name':"",
                   'iopl':'0',
                   'ilevs':'0',
                   'levb':'0',
                   'levt':'0',
                   'plt':'0',
                   'iopa':'0',
                   'imsk':'0',
                   'imn':'0',
                   'iwt':'0',
                   'ts':'.false.',
                   'levlst':'0',
                   'rlevlst':'0.0',
                   'inth':'0',
                   'isth':'0',
                   'iest':'0',
                   'iwst':'0',
                   'pslist':'0',
                   'tsnum':'0',
                   'tblim':'0',
                   'ttlim':'0',
                   'tblimr':'0.0',
                   'ttlimr':'0.0',
                   'tnlim':'0',
                   'tslim':'0',
                   'telim':'0',
                   'twlim':'0',
                   }

        for obj in config.get_value():
            if re.search(r'namelist:time\(',obj) or obj == "namelist:time":
                for item in tim_check.keys():
                    if not self.get_setting_value(
                                    config, [obj, item]):
                        self.add_setting(config, [obj, item], tim_check[item] )

            if re.search(r'namelist:domain\(',obj):
                timeseries=self.get_setting_value(
                                    config, [obj, "ts"])
                # will return none when run out of domain namelists
                if not timeseries:
                    break
                elif timeseries == "'Y'" :
                    self.change_setting_value(config,
                                              [obj, "ts"],
                                                  ".true.")
                elif timeseries == "'N'" :
                    self.change_setting_value(config,
                                              [obj, "ts"],
                                                  ".false.")
                for item in dom_check.keys():
                    if not self.get_setting_value(
                                    config, [obj, item]):
                        self.add_setting(config, [obj, item], dom_check[item] )
        return config, self.reports

class vn90_t6037(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6037 by Ian Boutle."""

    BEFORE_TAG = "vn9.0_t6006"
    AFTER_TAG = "vn9.0_t6037"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
#rhcrit parametrization settings
        l_rhcpt=self.get_setting_value(config, ["namelist:run_cloud","l_rhcpt"])
        self.remove_setting(config, ["namelist:run_cloud","l_rhcpt"])
        if l_rhcpt == ".true.":
            self.add_setting(config, ["namelist:run_cloud","i_rhcpt"], "1")
        else:
            self.add_setting(config, ["namelist:run_cloud","i_rhcpt"], "0")
#area cloud fraction settings
        l_cld_area=self.get_setting_value(config, ["namelist:run_cloud","l_cld_area"])
        self.remove_setting(config, ["namelist:run_cloud","l_cld_area"])
        self.remove_setting(config, ["namelist:run_cloud","l_acf_cusack"])
        if l_cld_area == ".false.":
            self.add_setting(config, ["namelist:run_cloud","i_cld_area"], "0")
        else:
            self.add_setting(config, ["namelist:run_cloud","i_cld_area"], "1")
#remove fixbug_pc2_checks
        self.remove_setting(config, ["namelist:run_cloud","i_fixbug_pc2_checks"])
#remove fixbug_pc2_qcl_incr
        self.remove_setting(config, ["namelist:run_cloud","l_fixbug_pc2_qcl_incr"])
#add forced_cu
        self.add_setting(config, ["namelist:run_cloud","forced_cu"],"0")
#add l_add_cca_to_mcica
        self.add_setting(config, ["namelist:run_cloud","l_add_cca_to_mcica"],".false.")
#give non-PC2 standard jobs some sensible PC2 settings so that it can
#be easily turned on in future
        l_pc2=self.get_setting_value(config, ["namelist:run_cloud","l_pc2"])
        if l_pc2 == ".false.":
            self.change_setting_value(config, ["namelist:run_cloud","l_ensure_min_in_cloud_qcf"],".true.")
            self.change_setting_value(config, ["namelist:run_cloud","l_fixbug_pc2_mixph"],".true.")
            self.change_setting_value(config, ["namelist:run_cloud","cff_spread_rate"],"2.0e-2")
            self.change_setting_value(config, ["namelist:run_cloud","l_micro_eros"],".true.")
            self.change_setting_value(config, ["namelist:run_cloud","i_pc2_erosion_method"],"3")
            self.change_setting_value(config, ["namelist:run_cloud","dbsdtbs_turb_0"],"7.5e-6")
#rename pc2_falliceshear_method
        falliceshear_method=self.get_setting_value(config, ["namelist:run_cloud","pc2_falliceshear_method"])
        self.remove_setting(config, ["namelist:run_cloud","pc2_falliceshear_method"])
        self.add_setting(config, ["namelist:run_cloud","falliceshear_method"], falliceshear_method)
        return config, self.reports

class vn90_t6042(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6042 by Harry Shepherd"""

    BEFORE_TAG = "vn9.0_t6037"
    AFTER_TAG = "vn9.0_t6042"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        l_ukca = self.get_setting_value(config,
                                        ['namelist:run_ukca',
                                         'l_ukca'])

        # if i_ukca_interval does not exist, this will return None
        i_ukca_interval = self.get_setting_value(config,
                                                 ['namelist:run_ukca',
                                                  'i_ukca_interval'])

        i_ukca_chem = self.get_setting_value(config,
                                             ['namelist:run_ukca',
                                              'i_ukca_chem'])

        # if the variable i_ukca_interval doesnt exist, we need to add
        # a fix for #5434, to add the variable when i_ukca_chem=1 as this
        # was missed. Set it to have a default value of 1
        if i_ukca_interval is None:
            if 'true' in l_ukca and int(i_ukca_chem) == 1:
                self.add_setting(config,
                                 ['namelist:run_ukca', 'i_ukca_interval'],
                                 value = '1')
        else:
            if 'false' in l_ukca:
                # change the i_ukca_interval variable from user ignored
                # to trigger ignored
                self.ignore_setting(config,
                                    ['namelist:run_ukca','i_ukca_interval'],
                                    state=config.STATE_SYST_IGNORED)
            elif 'true' in l_ukca:
                # enable the setting if ukca is enabled
                self.enable_setting(config,
                                    ['namelist:run_ukca','i_ukca_interval'])
            else:
                pass
        
        return config, self.reports

class vn90_t6045(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6045 by Jonathan Wilkinson."""

    BEFORE_TAG = "vn9.0_t6042"
    AFTER_TAG = "vn9.0_t6045"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        niter_bs = self.get_setting_value(
               config, ["namelist:run_precip", "niter_bs"])

        l_warm_new = self.get_setting_value(
               config, ["namelist:run_precip", "l_warm_new"])

        l_mcr_arcl = self.get_setting_value(
               config, ["namelist:run_precip", "l_mcr_arcl"])

        # Set a default value for c_r_correl where it is not
        # already being used (l_warm_new is false)
        if l_warm_new == ".false.":
            self.change_setting_value(
                   config, ["namelist:run_precip", "c_r_correl"], "0.9")

        # Set a default value for arcl_inhom_sc where it is not
        # already being used (l_mcr_arcl is false)
        if l_mcr_arcl == ".false.":
            self.change_setting_value(
                   config, ["namelist:run_precip", "arcl_inhom_sc"], "1.4")

        # Replace any value of niter_bs that is less than zero.
        # This should catch out where it is wrongly displayed as
        # 0 in rose apps, which is actually outside the range.
        if niter_bs.isdigit() and float(niter_bs) < 1:
            self.change_setting_value(
               config, ["namelist:run_precip", "niter_bs"], "1")


        return config, self.reports

class vn90_t5544(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5544 by Rachel Stratton."""

    BEFORE_TAG = "vn9.0_t6045"
    AFTER_TAG = "vn9.0_t5544"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        # Get values of convection scheme and cape_opt
        cape_opt_local = self.get_setting_value(
                config, ["namelist:run_convection", "cape_opt"])
        # Alter according to convection scheme
        # 0 - no changes required  
        # 4 - no changes required  - still uses 
        # 0 & 4 will add cldbase_opt_dp &  cldbase_opt_md which will get
        # commented out as not active
        # 5, 6 - leave cape_opt  as will get commented out add cldbase_opt_dp &
        #        cldbase_opt_md set to the same value as cape_opt had.
        self.add_setting(config, ["namelist:run_convection",
                                      "cldbase_opt_dp"],
                             str(cape_opt_local))
        self.add_setting(config, ["namelist:run_convection",
                                      "cldbase_opt_md"],
                             str(cape_opt_local))
        return config, self.reports

class vn90_t4481(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4481 by Jonathan Wilkinson."""

    BEFORE_TAG = "vn9.0_t5544"
    AFTER_TAG = "vn9.0_t4481"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:temp_fixes", "l_fix_drop_settle"], ".false.")
        return config, self.reports


class vn90_t6074(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6074 by Steve Wardle."""

    BEFORE_TAG = "vn9.0_t4481"
    AFTER_TAG = "vn9.0_t6074"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Remove C98
        self.remove_setting(config, ["namelist:umsections","section_c98"])

        return config, self.reports


class vn90_t6092(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6092 by Glenn Greed."""
    
    BEFORE_TAG = "vn9.0_t6074"
    AFTER_TAG = "vn9.0_t6092"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        tim_unit=['unt1',
                  'unt2',
                  'unt3']

        tim_replace={"'T '":'1',
                     "'T'":'1',
                     "'H '":'2',
                     "'H'":'2',
                     "'DA'":'3',
                     "'DU'":'4'}

        for obj in config.get_value():
            if re.search(r'namelist:time\(',obj):
                for tim in tim_unit:
                    current=self.get_setting_value(config, [obj, tim])
                    for item in tim_replace.keys():
                        if current == item:
                            self.change_setting_value(config,
                                                  [obj, tim],
                                                  tim_replace[item])
        return config, self.reports


class vn90_t6027(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6027 by Richard Barnes."""

    BEFORE_TAG = "vn9.0_t6092"
    AFTER_TAG = "vn9.0_t6027"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here

        # REMOVE CAT A / DUPLICATED ITEMS                   # 
        ##################################################### 
        self.remove_setting(config, ["namelist:nlstcatm", "l_auto_debias"])

        # L_ INTO DIFFERENT NAMELISTS                     # 
        ##################################################### 
        # Get l_settings from nlstcatm 
        l_use_arclbiom = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_arclbiom"])
        l_use_arclblck = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_arclblck"])
        l_use_arcldlta = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_arcldlta"])
        l_use_arcldust = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_arcldust"])
        l_use_arclocff = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_arclocff"])
        l_use_arclsslt = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_arclsslt"])
        l_use_arclsulp = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_arclsulp"])
        l_use_bmass_autoconv = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_bmass_autoconv"])
        l_use_nitrate_autoconv = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_nitrate_autoconv"])
        l_use_ocff_autoconv = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_ocff_autoconv"])
        l_use_seasalt_autoconv = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_seasalt_autoconv"])
        l_use_sulphate_autoconv = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_sulphate_autoconv"])

        # remove from existing namelist                     #
        #####################################################
        self.remove_setting(config, ["namelist:nlstcatm", "l_use_aod"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_use_arclbiom"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_use_arclblck"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_use_arcldlta"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_use_arcldust"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_use_arclocff"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_use_arclsslt"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_use_arclsulp"])
        self.remove_setting(config, ["namelist:nlstcatm",
                                     "l_use_bmass_autoconv"])
        self.remove_setting(config, ["namelist:nlstcatm",
                                     "l_use_nitrate_autoconv"])
        self.remove_setting(config, ["namelist:nlstcatm",
                                     "l_use_ocff_autoconv"])
        self.remove_setting(config, ["namelist:nlstcatm",
                                     "l_use_seasalt_autoconv"])
        self.remove_setting(config, ["namelist:nlstcatm",
                                     "l_use_soot_autoconv"])
        self.remove_setting(config, ["namelist:nlstcatm",
                                     "l_use_sulphate_autoconv"])

        # Adding variables to their new namelists 
        self.add_setting(config, ["namelist:run_radiation","l_use_arclbiom"],
        l_use_arclbiom)
        self.add_setting(config, ["namelist:run_radiation","l_use_arclblck"],
        l_use_arclblck)
        self.add_setting(config, ["namelist:run_radiation","l_use_arcldlta"],
        l_use_arcldlta)
        self.add_setting(config, ["namelist:run_radiation","l_use_arcldust"],
        l_use_arcldust)
        self.add_setting(config, ["namelist:run_radiation","l_use_arclocff"],
        l_use_arclocff)
        self.add_setting(config, ["namelist:run_radiation","l_use_arclsslt"],
        l_use_arclsslt)
        self.add_setting(config, ["namelist:run_radiation","l_use_arclsulp"],
        l_use_arclsulp)

        self.add_setting(config, ["namelist:run_precip","l_use_bmass_autoconv"],
        l_use_bmass_autoconv)
        self.add_setting(config, ["namelist:run_precip",
        "l_use_nitrate_autoconv"],l_use_nitrate_autoconv)
        self.add_setting(config, ["namelist:run_precip","l_use_ocff_autoconv"],
        l_use_ocff_autoconv)
        self.add_setting(config, ["namelist:run_precip",
        "l_use_seasalt_autoconv"],l_use_seasalt_autoconv)
        self.add_setting(config, ["namelist:run_precip",
        "l_use_sulphate_autoconv"],l_use_sulphate_autoconv)


        return config, self.reports


class vn90_t6076(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6076 by Paul Cresswell.
       Retire the namelist inputs for providing RUNID."""

    BEFORE_TAG = "vn9.0_t6027"
    AFTER_TAG = "vn9.0_t6076"

    def upgrade(self, config, meta_config=None):

        """If $RUNID is already defined, do not overwrite it.
           Otherwise, use the value from the namelist."""
        expt = re.sub(r"'", "", (self.get_setting_value(
            config, ["namelist:nlstcall", "expt_id"])))
        job = re.sub(r"'", "", (self.get_setting_value(
            config, ["namelist:nlstcall", "job_id"])))
        runid = expt + job

        # The SCM doesn't use RUNID at all:
        if config.get(["namelist:scm_cntl"]) is None:
          self.add_setting(config, ["env", "RUNID"], runid)

        self.remove_setting(config, ["namelist:nlstcall", "expt_id"])
        self.remove_setting(config, ["namelist:nlstcall", "job_id"])

        return config, self.reports

class vn90_t4177(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4177 by Glenn Greed."""

    BEFORE_TAG = "vn9.0_t6076"
    AFTER_TAG = "vn9.0_t4177"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        lbc_warn="""Upgrade macro has added  l_old_lbc_file=.false. by default.
                    It is not possible for this upgrade macro to determine if your
                    app is reading in LBCs with 10 (new) or 13 (old) compulsory
                    fields. If you are using older LBCs please set to true manually.
                 """
        self.add_report("namelist:lbc_options", "l_old_lbc_file", None,
                            info=lbc_warn,is_warning=True)
        self.add_setting(config, ["namelist:lbc_options","l_old_lbc_file"],".false.")
        return config, self.reports


class vn90_t5882(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5882 by Nick Savage."""

    BEFORE_TAG = "vn9.0_t4177"
    AFTER_TAG = "vn9.0_t5882"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        
        self.add_setting(config,["namelist:run_ukca","dir_strat_aer"],"''")
        self.add_setting(config,["namelist:run_ukca","dts0"],"0")
        self.add_setting(config,["namelist:run_ukca","jvspec_dir"],"''")
        self.add_setting(config,["namelist:run_ukca","l_mode_bln_on"],".false.")
        self.add_setting(config,["namelist:run_ukca","nit"],"0")

        return config, self.reports


class vn90_t6123(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6123 by Steve Wardle."""

    BEFORE_TAG = "vn9.0_t5882"
    AFTER_TAG = "vn9.0_t6123"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Only apply to SCM apps
        if config.get(["namelist:scm_cntl"]):
            src_cntlatm = (
                self.get_setting_value(config, ["file:CNTLATM","source"])).split()
            src_shared = (
                self.get_setting_value(config, ["file:SHARED","source"])).split()

            # Remove from CNTLATM (if it appears there)
            if "namelist:run_eng_corr" in src_cntlatm:
                src_cntlatm.remove("namelist:run_eng_corr")
                self.change_setting_value(
                    config, ["file:CNTLATM","source"]," ".join(src_cntlatm))
                
            # If already in SHARED, remove it to make sure it gets added
            # in the right place
            if "namelist:run_eng_corr" in src_shared:
                src_shared.remove("namelist:run_eng_corr")

            # Insert it after run_cloud, which is read just before it 
            # in the code
            index = src_shared.index("namelist:run_cloud")
            src_shared.insert(index+1,"namelist:run_eng_corr")
            self.change_setting_value(
                config, ["file:SHARED","source"]," ".join(src_shared))

        return config, self.reports

class vn90_t5849(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5849 by Sean Swarbrick."""

    BEFORE_TAG = "vn9.0_t6123"
    AFTER_TAG = "vn9.0_t5849"

# ROSE rationalisation of ancillary updating
# This macro merges each upanca name list into its related items name list,
#  and deletes the upanca name lists. 

# --------------------------------
    def refno_to_item(self,refno):
      """Convert ancil ref no to stash item"""
    
      # Array of stash items indexed by ancil ref number
    
      items = [30,  33,  34,  35,  36,  37,  60,  96,  23,  20,
               40,  41, 190,  43,  44,   0,  46,  47,   0,   0,
                0,   0,   0,   0,   0,  26,  31,  24,  32,  28,
               29,  93, 274, 275,  48,   9,   0,   0,  58,  59,
                0,   0,   0,  57,  90,  17,  18, 301, 302, 303,
              304, 305, 306, 307, 308, 309, 310, 311, 312, 313,
              314, 315, 316, 317, 318, 319, 320, 127, 128, 129,
                0, 121, 122, 123, 124, 125, 126, 251, 207,   0,
                0, 160, 216, 217, 218, 213, 219, 220, 223, 321,
              322, 323, 324, 325, 326, 327, 328, 329, 330, 331,
              332, 333, 334, 335, 336, 337, 338, 339, 340, 341,
              505, 418, 419, 420, 421, 422, 423, 424, 425, 426,
              130, 131, 132, 153, 151, 152,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                0,   0,   0,   0,   5,   6, 351, 352, 353, 354,
              355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 
              365, 366, 367, 368, 369, 370, 371, 480, 481, 482,
              483, 484, 485, 486, 487, 134, 135,   7,   0,   0,
                0,   0,   0, 243, 244, 245]

      # Get item number for given ancil ref number
    
      item = items[refno-1]

      return item
    
# --------------------------------------------    
    def merge_upanca_to_items(self, config, items_req, namelstnm):
      """Transfer upanca name list info to items name lists"""      

      array1 = []
    
      for obj in config.get_value():
      
        # Loop through the upanca name lists
        
        if re.search(r'namelist:upanca', obj):
        
          anc_ref_no = self.get_setting_value(config, [obj, "anc_ref_no"])
          interval   = self.get_setting_value(config, [obj, "interval"])
          period     = self.get_setting_value(config, [obj, "period"])

          # Convert anc_ref_no to related stash item

          item = self.refno_to_item(int(anc_ref_no))

          for obi in config.get_value():
          
            # Loop through items name lists
          
            if re.search(r'namelist:item', obi):

              stash_req = self.get_setting_value(config, [obi, "stash_req"])
              # check whether there are commas in stash_req
              if re.search(r',', stash_req):
                continue
              else:              
              # Find which items NL relates to the current upanca NL
              # Insert upanca info into items NL
              
                if stash_req == str(item):

                  self.add_setting(config, [obi, "interval"], interval)
                  self.add_setting(config, [obi, "period"], period)
                  self.add_setting(config, [obi, "update_anc"], ".true.")

                  # Get ancil file path from dictionary and add it in
                  for anc, stash in items_req.iteritems():
                  
                    stash_no = str(stash).replace("[","")
                    stash_no = stash_no.replace("]","")
                    stash_no = stash_no.replace("'","")
                    anc_line = "'"+anc+"'"

                    if re.search(r',', stash_no):
                      array1 = stash_no.split(',')
                      array2 = []
                      for i in array1:
                        array2.append(i.replace(" ",""))
                      for i in array2:                      
                        if stash_req == i:
                          self.add_setting(config,[obi, "ancilfilename"],anc_line)
                          break
                    else:
                      if stash_req == stash_no:
                        self.add_setting(config, [obi, "ancilfilename"],anc_line)
                        break              
                      
      return config      

# ----------------------------------
    def delete_upanca(self, config):
      """Delete upanca name lists from rose-app.conf"""

      num_nml = 0

      for obj in config.get_value():
      
        # Count upanca name lists
        
        if re.search(r'namelist:upanca', obj):

          num_nml = num_nml + 1
 
      nml = "namelist:upanca(0)"
          
      for i in range(num_nml):
      
        # delete upanca namelists
      
        nml = nml.replace(str(i),str(i+1))

        self.remove_setting(config, [nml])

      return config

# ----------------------------------------------

    def UM_version(self, config) :
        ''' a simple routine to identify the UM version 
            from the rose-app.conf '''

        meta = self.get_setting_value(config, ["meta"])
        pversion = re.compile('um-atmos/(?P<VERSION>vn[0-9].[0-9])')
        version = pversion.search(meta)

        if version:
           return version.group('VERSION')

        return None

# ----------------------------------------------

    def stash2ancilenv (self, vn,umdir) :
        '''Given a UM version, construct a dictionary of 
           fieldcodes to ancil env names'''

        # ANCILMASTER files need to be opened. os.path.join?
        anc_file=open(umdir+"/"+vn+"/ctldata/ANCILmaster/ANCILfiles_A", "r")
        anc_field=open(umdir+"/"+vn+"/ctldata/ANCILmaster/ANCILfields_A", "r")

        # initiate dictionaries
        ancil_filenumber={}
        ancil_fieldnumber={}
        stash2ancilenvname={}

        # read in the ancilfiles and construct env var to filenumber link
        for line in anc_file.readlines():

            # simply split up record and keep required 
            # ENV VAR and file number info.
            line = line.replace("\n", "")
            first_line=re.search('^1', line)

            # only interested in lines thst start with a 1 
            # and ignore end record.
            if first_line: 
                file_rec= line.split('|')

                if int(file_rec[1]) != -1:
                    envname=file_rec[3].replace(" ", "")
                    ancil_filenumber[int(file_rec[1])]=envname


        # read in ancil fields and construct env var to stashcodes
        for line in anc_field.readlines():

            # this time need to read in pair of record lines to 
            # obtain stashcode file number linkage.
            line = line.replace("\n", "")
            first_line=re.search('^1', line)
            second_line=re.search('^2', line)

            # take information from both lines in form of lists
            if first_line: 
                field_rec1= line.split('|')
                continue

            if second_line: 
                field_rec2= line.split('|')
                if int(field_rec2[1]) != 0:
                    ancil_fieldnumber[(1000*int(field_rec1[3])+
                                     int(field_rec1[4]))]=int(field_rec2[1])

        # finally create dictionary of stashcodes with ENV VAR for ancil
        for k,v in ancil_fieldnumber.iteritems():
            stash2ancilenvname[str(k)]=ancil_filenumber[ancil_fieldnumber[k]]

        return stash2ancilenvname


# ----------------------------------------------

    def ancilenvs (self, config, ancil_dir_file, initfilenv, umdir) :
        '''source the ancil versions files along with initifileenv 
           to set up the required dictionary that
           provides the actual filepaths for each UM ENV VAR'''

        # extract the env variables from the rose config app.
        for envvar in self.get_setting_value(config, [ 'env' ] ):
            value = self.get_setting_value(config, [ 'env', envvar])

            # ensure UMDIR is set to local form 
            value = value.replace("$UMDIR",umdir)
            value = value.replace("/projects/um1",umdir)

            os.environ[envvar] = value
         

        if initfilenv=="" and ancil_dir_file!="" :
            # set up bash command to source the file provided.
            command = ". %s ; env"%(ancil_dir_file) 
 
        elif initfilenv!="" and ancil_dir_file=="" :
            # set up bash command to source the file provided.
            command = ". %s ; env"%(initfilenv)  

        elif (ancil_dir_file!="" and initfilenv!="") :
            # set up bash command to source the two files provided.
            command = ". %s;. %s ; env"%(ancil_dir_file,initfilenv) 

        else :
            # not enough information has been provided to progress
            # if no items are requested in app, eg SCM then that is
            # not a problem, but if app does set these then fail
 
            for obj in config.get_value():

                if re.search(r'namelist:items', obj):
                   
                  raise UpgradeError (
                   'At the very least either ANCIL_VERSIONS or '+ 
                   'INITFILENV must be provided to proceed ' )                  

                else :
                    
                    return os.environ

        rc, stdout, stderr = self.run_command(command, shell=True)
 
        # create the dictionary of actual files paths linked to ENV VARs
        for line in stdout:
            (key, _, value) = line.partition("=")

            # ignore UM ENV VARs that are unset.
            if key=='':
                continue
            os.environ[key] = value

        return os.environ

# ----------------------------------------------

    def run_command(self, command, shell=False):
        '''Given a command as a string, run it and return the exit code, 
        standard out and standard error. The optional shell argument 
        allows a shell to be spawned to allow multiple commands to be run.'''

        import subprocess

        if shell:
            # Create the Popen object and connect out and err to pipes using
            # the shell=True option.
            p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE, shell=True)
        else:
            # Turn command into a list
            command_list = command.split()

            # Create the Popen object and connect out and err to pipes
            p = subprocess.Popen(command_list, stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)

        # Do the communiate and wait to get the results of the command
        stdout, stderr = p.communicate()

        # Reformat stdout
        stdout = ''.join(stdout)
        stdout = stdout.split("\n")

        # Get exit code
        rc = p.wait()

        return rc, stdout, stderr

# ----------------------------------------------

    def generate_stash_to_anc_dict(self, config, envvars, 
                                   stash2ancilenvname,hpc_file):
        '''Read in config object and extract the stashcode for any item
           that is to be sourced from an ancillary file. This then creates a 
           dictionary of stashcodes with ancillary files; the explicit paths.

           Inputs are the config object to be interogated, 
           the dictionary linking env vars to filenames and
           the dictionary linking stashcodes to ancil env var names '''
        items_req = {}
        namelstnm = []

        # find the upanca namelist objects.
        for obj in config.get_value():

            if re.search(r'namelist:upanca', obj):
               
                phash = re.compile('upanca\((?P<HASH>[a-g0-9]+)\)')
                namelstname = phash.search(obj) 

                anc_ref_no = self.get_setting_value(config, [obj, "anc_ref_no"])
               
                # If the above 'get' returns "None", we've reached the 
                # last item in the file, so end the loop
                if not anc_ref_no:
                    break

                stashcode = str(self.refno_to_item(int(anc_ref_no)))

                # create array of namelist names to reuse later.
                namelstnm.append(namelstname.group('HASH'))

                try:
                    ancfile=envvars[stash2ancilenvname[stashcode]]
                except KeyError as e:
                    missing_key = str(e)
                    raise UpgradeError("Unable to find %s in environment.\n"%(
                        missing_key) + "This app is not using std "+
                        "ancil_versions format. Manual editing of the app" +
                        " is required.\n" + 
                        "If not using a centrally controlled ancil versions" +
                        " file please ensure that the\nancil versions file" +
                        " and ancil filenames file are available on the file " +
                        "system\non which you are performing the app upgrade.")

                # if original app is a remote app need to ensure
                # the use of remote ancils path, installation specific.
                if hpc_file:
                   ancfile = ancfile.replace("/home/h01/frum","/projects/um1")

                if ancfile not in items_req:
                    items_req[ancfile]=[stashcode]
                else :
                    items_req[ancfile].append(stashcode)

        return items_req, namelstnm

# ----------------------------------------------

    def upgrade(self, config, meta_config=None):
        """For stash items with ancillary updating, insert upanca name list 
        information into the related items name list. Also delete the upanca
        name lists."""

        initfilenv='./file/INITFILEENV'
        roseapp='./rose-app.conf'

        # check that rose app is available.
        if os.path.isfile(roseapp):
            # file exists so continue
            pass
        else :
            # no file so exit with message.
            raise UpgradeError (
              'Cannot continue unable to find the rose app config file '
               +roseapp )                   


        # check whether initfilenv is available.
        if os.path.isfile(initfilenv):
            # file exists so continue
            pass
        else :
            # no file so set to empty string
            initfilenv=""

        # check whether $UMDIR is set in our local env or not.
        umdir=os.environ["UMDIR"]

        if not umdir:
            self.add_report("env", "UMDIR", None,
              '$UMDIR is not set - assuming /home/h01/frum', is_warning=True)
            umdir='/home/h01/frum'

        # This script will use vn8.6/9.0 ANCILmaster unless the rose app 
        # provided is clearly using an older UM version

        vn = self.UM_version(config)

        if not vn:
           self.add_report(None, None, None,
                  info="UM VN is not set - assuming vn9.0", is_warning=True)
           vn='vn9.0'

        # set up dictionary of   [ancil ENV names : stashcodes]
        # by using the ancilmaster file and fields tables.
        stash2ancilenvname = self.stash2ancilenv(vn,umdir)

        # identify the ancil versions files used by given UM job
        ancil_dir_file = self.get_setting_value(config, 
                                               ["env", "ANCIL_VERSIONS"])

        # if using MP single threads env then HPC (ugly and MO specific)
        MP_thread_file = self.get_setting_value(config, 
                                               ["env", "MP_SINGLE_THREAD"])

        # work needs to be done on the local machine so redirect
        # ancil_versions if appropriate and possible.
        ancil_dir_file_local=""
        hpc_app=False
        
        if ancil_dir_file:
            if (ancil_dir_file.find("/projects/um1") != -1 or 
                MP_thread_file != None ):
                hpc_app=True

            ancil_dir_file_local = ancil_dir_file.replace("/projects/um1",
                                                          umdir)

            # some apps are not using the central ancil_versions
            if ancil_dir_file.find("/data/um1/") != -1:
                error=('This app is not using std ancil_versions format. '+
                       ancil_dir_file+
                     '\nYou will have to update the app by hand!\n' + 
                        "ancil_versions format. Manual editing of the app "+
                     "is required.\n" + 
                     "If not using a centrally controlled ancil versions" +
                     " file please ensure that the\nancil versions file " +
                     " and ancil filenames file are available on the file " +
                     " system on which you are performing the app upgrade.")
                raise UpgradeError(error)

        # source the ancil versions files so to aid the creation of 
        # a dictionary of [ ancil ENV names: ancillary full paths]
        envvars = self.ancilenvs (config, ancil_dir_file_local,
                                  initfilenv,umdir)

        # set up a dictionary of requested [ ancillary full paths: stashcodes]
        items_req, namelstnm = self.generate_stash_to_anc_dict(config, 
                               envvars, stash2ancilenvname,hpc_app)


        #   Transfer upanca information to items namelists
        config = self.merge_upanca_to_items(config, items_req, namelstnm)


        #   Remove upanca from file:NAMELIST   
        namelsts = self.get_setting_value(config, ["file:NAMELIST","source"],)

        items_present = False
        for obi in config.get_value():
            if re.search(r'namelist:item', obi):
                items_present = True
                break

        
        if (namelsts):
          namelsts=namelsts.replace("namelist:upanca(:) ", "")
          if items_present:
              namelsts=namelsts.replace("namelist:ancilcta ",
                                        "namelist:ancilcta namelist:items(:) ")
          self.change_setting_value(config, ["file:NAMELIST","source"],namelsts)

#   Delete upanca name lists from rose-app.conf 
          config = self.delete_upanca(config)


        return config, self.reports

class vn90_t6117(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6117 by Stuart Whitehouse."""

    BEFORE_TAG = "vn9.0_t5849"
    AFTER_TAG = "vn9.0_t6117"

    def remove_namelist(self, config, fname, namelist):
        """Remove a namelist from a file."""
        
        source_line = self.get_setting_value(config, ["file:%s"%(fname),
                           "source"])
        if source_line == None:
            return
        
        source_line = re.sub(r'namelist:%s\s*'%(namelist), r'', source_line)
        
        self.change_setting_value(config, ["file:%s"%(fname), "source"], 
                                  source_line)

        return        

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        self.remove_setting(config, [ "namelist:intfcnsta" ])
        self.remove_namelist(config, 'NAMELIST', 'intfcnsta')
        self.remove_namelist(config, 'CNTLATM', 'intfcnsta')

        return config, self.reports

class vn90_t6023(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6023 by Andy Malcolm."""

    BEFORE_TAG = "vn9.0_t6117"
    AFTER_TAG = "vn9.0_t6023"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:ioscntl", "ios_print_start_time"],
                         ".false.")
        self.add_setting(config, ["namelist:io_control", "print_runtime_info"],
                         ".false.")
        return config, self.reports


class vn90_t5401(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5401 by <Warren Tennant>."""

    BEFORE_TAG = "vn9.0_t6023"
    AFTER_TAG = "vn9.0_t5401"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:temp_fixes", "l_roughnesslength_fix"], ".false.")
        return config, self.reports


class vn90_t6173(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6173 by Stuart Whitehouse."""

    BEFORE_TAG = "vn9.0_t5401"
    AFTER_TAG = "vn9.0_t6173"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        for obj in config.get_value():
            if re.search(r'namelist:time\(',obj) or obj == "namelist:time":
                istd = self.get_setting_value(config, [ obj, "istd"])
                if istd:
                    self.remove_setting(config, [obj, "istd"])
                    
                ietd = self.get_setting_value(config, [ obj, "ietd"])
                if ietd:
                    self.remove_setting(config, [obj, "ietd"])

        return config, self.reports


class vn90_t5864(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5864 by Joe Mancell."""

    BEFORE_TAG = "vn9.0_t6173"
    AFTER_TAG = "vn9.0_t5864"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.remove_setting(config, ["namelist:nlsizes", "pp_len_inthd"])
        self.remove_setting(config, ["namelist:nlsizes", "pp_len_realhd"])
        self.remove_setting(config, ["namelist:nlsizes", "tr_levels"])
        return config, self.reports


class vn90_t6152(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6152 by Ian Boutle."""

    BEFORE_TAG = "vn9.0_t5864"
    AFTER_TAG = "vn9.0_t6152"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.remove_setting(config, ["namelist:run_cloud","cloud_fraction_method"]) 
        return config, self.reports

class vn90_t6135(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6135 by Harry Shepherd (hshep)."""

    BEFORE_TAG = "vn9.0_t6152"
    AFTER_TAG = "vn9.0_t6135"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Remove the trailing spaces from the variable time_convention in
        # nlstcall to ensure compatibility with updated meta_data
        time_convention = self.get_setting_value(config,
                                                 ['namelist:nlstcall',
                                                  'time_convention'])
        # for some reason rstrip() isn't working in this scenario
        time_convention = time_convention.replace(' ','')
        self.change_setting_value(config,
                                  ['namelist:nlstcall','time_convention'],
                                  time_convention)
        return config, self.reports

class vn90_t6097(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6097 by Stuart Whitehouse."""

    BEFORE_TAG = "vn9.0_t6135"
    AFTER_TAG = "vn9.0_t6097"

    def remove_namelist(self, config, fname, namelist):
        """Remove a namelist from a file."""
        
        source_line = self.get_setting_value(config, ["file:%s"%(fname),
                           "source"])
        if source_line == None:
            return
        
        source_line = re.sub(r'namelist:%s\s*'%(namelist), r'', source_line)
        
        self.change_setting_value(config, ["file:%s"%(fname), "source"], 
                                  source_line)

        return        

    def do_var_stash_macro(self, config):
        # Copy the VAR stash macro settings
        
        self.add_setting(config, ['namelist:nlstcall_pp(150)'])

        # Start time
        varstart = self.get_setting_value(config, [ "namelist:nlstcall_var",
                                         "varstart"])
        if varstart:
            self.add_setting(config, ['namelist:nlstcall_pp(150)', 'i_ppis'],
                          str(varstart))

        # Period
        varintvl = self.get_setting_value(config, [ "namelist:nlstcall_var",
                                         "varintvl"])
        if varintvl:
            self.add_setting(config, ['namelist:nlstcall_pp(150)', 'i_ppif'],
                          str(varintvl))

        # units
        varunits = self.get_setting_value(config, [ "namelist:nlstcall_var",
                                         "varunits"])
        if varunits:
            self.add_setting(config, ['namelist:nlstcall_pp(150)', 'i_ppiu'],
                          str(varunits))

        # Packing
        packvar = self.get_setting_value(config, [ "namelist:nlstcall_var",
                                         "packvar"])
        if packvar == 'Y':
            self.add_setting(config, ['namelist:nlstcall_pp(150)', 'i_ppx'],
                          str(1))
        else:
            self.add_setting(config, ['namelist:nlstcall_pp(150)', 'i_ppx'],
                          str(0))

        # End time
        varend = self.get_setting_value(config, [ "namelist:nlstcall_var",
                                         "varend"])
        if varend:
            self.add_setting(config, ['namelist:nlstcall_pp(150)', 'i_ppie'],
                          str(varend))

        self.add_setting(config, ['namelist:nlstcall_pp(150)', 'i_ppg'], 
                          str("'N'"))

        varintfc = self.get_setting_value(config, [ "namelist:nlstcall_var",
                                         "varintfc"])
  
        if varintfc == '4':
            self.add_setting(config, ['namelist:nlstcall_pp(150)', 'i_ppi'], 
                             str("'Y'"))
        else:
            self.add_setting(config, ['namelist:nlstcall_pp(150)', 'i_ppi'], 
                             str("'N'"))
        
        self.add_setting(config, ['namelist:nlstcall_pp(150)', 'i_ppos'], 
                          str(0))

        # No archiving on this stream
        self.add_setting(config, ['namelist:nlstcall_pp(150)', 'i_ppa'], 
                         str("'N'"))

        return config

    def do_makebc_stash_macro(self, config):
        # Copy the MakeBC stash macro settings
        # Start time
        self.add_setting(config, ['namelist:nlstcall_pp(164)'])

        mbcstart = self.get_setting_value(config, [ "namelist:nlstcall_makebc",
                                         "mbc_strt"])
        if mbcstart:
            self.add_setting(config, ['namelist:nlstcall_pp(164)', 'i_ppis'],
                          str(mbcstart))
        
        # Period
        mbcintvl = self.get_setting_value(config, [ "namelist:nlstcall_makebc",
                                         "mbc_freq"])
        if mbcintvl:
            self.add_setting(config, ['namelist:nlstcall_pp(164)', 'i_ppif'],
                          str(mbcintvl))

        # units
        mbcunits = self.get_setting_value(config, [ "namelist:nlstcall_makebc",
                                         "mbc_unt"])
        if mbcunits:
            self.add_setting(config, ['namelist:nlstcall_pp(164)', 'i_ppiu'],
                          str(mbcunits))

        # Packing (not sure why this is in the VAR namelist!)
        packmbc = self.get_setting_value(config, [ "namelist:nlstcall_var",
                                         "packvarbc"])
        if packmbc == 'Y':
            self.add_setting(config, ['namelist:nlstcall_pp(164)', 'i_ppx'],
                          str(1))
        else:
            self.add_setting(config, ['namelist:nlstcall_pp(164)', 'i_ppx'],
                          str(0))

        # End time
        mbcend = self.get_setting_value(config, [ "namelist:nlstcall_makebc",
                                         "mbc_end"])
        if mbcend:
            self.add_setting(config, ['namelist:nlstcall_pp(164)', 'i_ppie'],
                          str(mbcend))

        self.add_setting(config, ['namelist:nlstcall_pp(164)', 'i_ppg'], 
                          str("'N'"))
        self.add_setting(config, ['namelist:nlstcall_pp(164)', 'i_ppos'], 
                          str(0))

        imkbc = self.get_setting_value(config, [ "namelist:nlstcall_makebc",
                                         "imkbc"])
        if imkbc == '1':
            self.add_setting(config, ['namelist:nlstcall_pp(164)', 'i_ppi'], 
                             str("'Y'"))
        else:
            self.add_setting(config, ['namelist:nlstcall_pp(164)', 'i_ppi'], 
                             str("'N'"))

        # No archiving on this stream
        self.add_setting(config, ['namelist:nlstcall_pp(164)', 'i_ppa'], 
                         str("'N'"))
        return config


    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Remove GRIB packing for meaning files
        self.remove_setting(config, ['namelist:nlstcall_qs', 'ppxg'])

        config = self.do_var_stash_macro(config)
        config = self.do_makebc_stash_macro(config)

        self.remove_setting(config, ['namelist:nlstcall_makebc'])
        self.remove_namelist(config, 'NAMELIST', 'nlstcall_makebc')
        self.remove_setting(config, ['namelist:nlstcall_var'])
        self.remove_namelist(config, 'NAMELIST', 'nlstcall_var')

        conversion_dict = {
           "''" : 0, "'H'" : 1,  "'DA'" : 2, "'T'" : 3, "'RM'" : 4,
                          }

        # Now loop over all nlstcall_pp namelists and:
        # * Remove GRIB packing
        # * Convert i_ppi to be a logical
        # * Convert i_ppu to be an integer
        for obj in config.get_value():
            if re.search(r'namelist:nlstcall_pp\(',obj):
                self.remove_setting(config, [ obj, 'i_ppg'])
                i_ppi = self.get_setting_value(config, [ obj, 'i_ppi'])
                if i_ppi == "'Y'":
                    self.change_setting_value(config, [obj,'i_ppi'], '.true.')
                else:
                    self.change_setting_value(config, [obj,'i_ppi'], '.false.')
                i_ppiu = self.get_setting_value(config, [ obj, 'i_ppiu'])
                self.change_setting_value(config, [ obj, 'i_ppiu'],
                                          str(conversion_dict[i_ppiu]))


        return config, self.reports        

class vn90_t6088(rose.upgrade.MacroUpgrade):

  """Upgrade macro for ticket #6088 by Sean J Swarbrick."""

  BEFORE_TAG = "vn9.0_t6097"
  AFTER_TAG = "vn9.0_t6088"

  def upgrade(self, config, meta_config=None):
    """Upgrade a UM runtime app configuration."""

    acp_dict = {'ac_obs_types': 0, 'ac_order': 0, 'l_lhn': ".false.", 'lac_mes': ".false.",
                'no_obs_files': 2, 'obs_format': 2, 'alpha_lhn': 0.333, 'diag_rdobs': 1,
                'epsilon_lhn': 0.333, 'fi_scale_lhn': 17000.0, 'lhn_range': 4, 'macdiag': 0,
                'nudge_lam': 0.0, 'relax_cf_lhn': 1.0, 'remove_neg_lh': ".false.", 
                'tgetoba': 0, 'tgetobb': 0, 'timea': 0, 'timeb': 0, 'use_conv_in_mops': ".true.",
                'npass_rf_lhn': 2, 'l_lhn_limit': ".false.", 'lhn_limit': 1.0, 
                'l_lhn_scale': ".true.",  'l_lhn_search': ".true.", 'l_lhn_fact': ".true.",
                'l_lhn_filt': ".true.", 'lhn_diag': ".true."}
    
    acdiag_dict = {'ldiagac': ".false.", 'lldac': ".false.", 'lrms': ".false.", 
                   'ltemp': ".false.", 'lverif': ".false."}
    
    if config.get(["namelist:acp(1)"]) is None:
    
      namelsts = self.get_setting_value(config, ["file:NAMELIST","source"],)    
      
      if not namelsts:
        pass
      else:      
        
        namelsts=namelsts.replace("namelist:clmchfcg ", "namelist:clmchfcg namelist:acp namelist:acdiag ")
        self.change_setting_value(config, ["file:NAMELIST","source"], namelsts)
       
        for item, value in acp_dict.iteritems():
         
          self.add_setting(config, ["namelist:acp", item], str(value)) 
          
        for item, value in acdiag_dict.iteritems():
        
          self.add_setting(config, ["namelist:acdiag", item], str(value))       
      
    else:    
      
      for item, value in acp_dict.iteritems():
      
        nl_value = self.get_setting_value(config, ["namelist:acp(1)", item])
        if not nl_value:
          pass
        else:
          acp_dict[item] = nl_value
        
      for item, value in acp_dict.iteritems():
      
        nl_value = self.get_setting_value(config, ["namelist:acp(2)", item])
        if not nl_value:
          pass
        else:
          acp_dict[item] = nl_value
   
      for item, value in acdiag_dict.iteritems():
      
        nl_value = self.get_setting_value(config, ["namelist:acdiag(1)", item])
        if not nl_value:
          pass
        else:
          acdiag_dict[item] = nl_value
          
      for item, value in acdiag_dict.iteritems():
      
        nl_value = self.get_setting_value(config, ["namelist:acdiag(2)", item])
        if not nl_value:
          pass
        else:
          acdiag_dict[item] = nl_value
     
     
      namelsts = self.get_setting_value(config, ["file:NAMELIST","source"],)
    
      if not namelsts:
        pass
      else:
      
        # convert acp(1), acp(2) to single NL acp
        # convert acdiag(1), acdiag(2) to single NL acdiag
        
        namelsts=namelsts.replace("namelist:acp(1)", "namelist:acp")
        namelsts=namelsts.replace("namelist:acp(2) ", "")
        
        namelsts=namelsts.replace("namelist:acdiag(1)", "namelist:acdiag")
        namelsts=namelsts.replace("namelist:acdiag(2) ", "")
        self.change_setting_value(config, ["file:NAMELIST","source"], namelsts)
     
        self.remove_setting(config, ["namelist:acp(1)"])
        self.remove_setting(config, ["namelist:acp(2)"])
        
        self.add_setting(config, ["namelist:acp"])
        
        for item, value in acp_dict.iteritems():
        
          self.add_setting(config, ["namelist:acp", item], str(value))
 
        self.remove_setting(config, ["namelist:acdiag(1)"])
        self.remove_setting(config, ["namelist:acdiag(2)"])         
        
        self.add_setting(config, ["namelist:acdiag"]) 
        
        for item, value in acdiag_dict.iteritems():
        
          self.add_setting(config, ["namelist:acdiag", item], str(value))
 
    return config, self.reports
    
class vn90_t5846(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5846 by Paul Selwood."""

    BEFORE_TAG = "vn9.0_t6088"
    AFTER_TAG = "vn9.0_t5846"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.remove_setting(config,["namelist:ustsnum"],
                            info="User STASHmasters are now defunct")
        # remove this namelist from the NAMELIST file
        namelist_source = self.get_setting_value(config, ["file:NAMELIST","source"])
        if namelist_source:
            namelist_source = re.sub(r'namelist:ustsnum ','',namelist_source)
            self.change_setting_value(config, ["file:NAMELIST","source"],namelist_source)
        # remove this namelist from the SHARED file
        shared_source = self.get_setting_value(config, ["file:SHARED","source"])
        if shared_source:
            shared_source = re.sub(r'namelist:ustsnum ','',shared_source)
            self.change_setting_value(config, ["file:SHARED","source"],shared_source)
        return config, self.reports

class vn90_t6085 (rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6085 by <Mohit Dalvi>."""

    BEFORE_TAG = "vn9.0_t5846"
    AFTER_TAG = "vn9.0_t6085"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

      # Introduce control variables (default 'OFF' state) to enable use of 
      # new Priestley conservation for moist/tracers 

        # Logicals for Priestley to tracers
        self.add_setting( config,
               ['namelist:run_sl', 'l_priestley_correct_tracers'],
                value = '.false.' )
                  
        self.add_setting( config,
               ['namelist:run_sl', 'tr_priestley_opt'],
                value = '0' )
                  
        self.add_setting( config,
               ['namelist:run_sl', 'l_tr_src_in_conserve'],
                value = '.false.' )
                
        # Logicals for Priestley to moisture
        
        self.add_setting( config,
               ['namelist:run_sl', 'l_priestley_correct_moist'],
                value = '.false.' )
                  
        self.add_setting( config,
               ['namelist:run_sl', 'moist_priestley_opt'],
                value = '0' )
                  
        self.add_setting( config,
               ['namelist:run_sl', 'l_moist_src_in_conserve'],
                value = '.false.' )
                  
        # Logicals for Priestley to UKCA
        self.add_setting( config,
               ['namelist:run_ukca', 'i_ukca_conserve_method'],
                value = '0' )
                  
        self.add_setting( config,
               ['namelist:run_ukca', 'i_ukca_hiorder_scheme'],
                value = '0' )
                  
        self.add_setting( config,
               ['namelist:run_ukca', 'l_ukca_src_in_conservation'],
                value = '.false.' )

        return config, self.reports

class vn90_t5694(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5694 by James Manners."""

    BEFORE_TAG = "vn9.0_t6085"
    AFTER_TAG = "vn9.0_t5694"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        #####################################################
        # planet_constants - new namelist                   #
        #####################################################
        self.add_setting(config, ["namelist:planet_constants"])
        self.add_setting(config, ["namelist:planet_constants","i_planet"], "3")
        self.add_setting(config, ["namelist:planet_constants","l_planet_orbit"], ".false.")
        self.add_setting(config, ["namelist:planet_constants","l_planet_locked"], ".false.")
        self.add_setting(config, ["namelist:planet_constants","planet_lock_lon"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","l_planet_grey_surface"], ".false.")
        self.add_setting(config, ["namelist:planet_constants","planet_t_surface"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","l_planet_g"], ".false.")
        self.add_setting(config, ["namelist:planet_constants","planet_radius"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_sidday"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_solday"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_epoch"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_e"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_de"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_lph"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_dlph"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_oblq"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_doblq"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_a"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_da"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_m"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","planet_dm"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","g"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","r"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","cp"], "0.0")
        self.add_setting(config, ["namelist:planet_constants","pref"], "0.0")
        namelist_source = self.get_setting_value(config, ["file:NAMELIST","source"])
        if namelist_source:
            namelist_source = re.sub(r'namelist:nlstcatm',
                                     r'namelist:nlstcatm namelist:planet_constants',
                                     namelist_source)
            self.change_setting_value(config, ["file:NAMELIST","source"], namelist_source)
        shared_source = self.get_setting_value(config, ["file:SHARED","source"])
        if shared_source:
            shared_source = re.sub(r'namelist:temp_fixes',
                                   r'namelist:temp_fixes namelist:planet_constants',
                                    shared_source)
            self.change_setting_value(config, ["file:SHARED","source"], shared_source)
        #####################################################
        # run_radiation - new logical                       #
        #####################################################
        self.add_setting(config, ["namelist:run_radiation","l_bs1999_abundances"], ".false.")
        #####################################################
        # r2lwclnl - new logicals                           #
        #####################################################
        self.add_setting(config, ["namelist:r2lwclnl","l_co_lw"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_co_lw2"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_nh3_lw"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_nh3_lw2"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_tio_lw"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_tio_lw2"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_vo_lw"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_vo_lw2"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_h2_lw"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_h2_lw2"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_he_lw"], ".false.")
        self.add_setting(config, ["namelist:r2lwclnl","l_he_lw2"], ".false.")
        #####################################################
        # r2swclnl - new logicals                           #
        #####################################################
        self.add_setting(config, ["namelist:r2swclnl","l_co_sw"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_co_sw2"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_nh3_sw"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_nh3_sw2"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_tio_sw"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_tio_sw2"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_vo_sw"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_vo_sw2"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_h2_sw"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_h2_sw2"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_he_sw"], ".false.")
        self.add_setting(config, ["namelist:r2swclnl","l_he_sw2"], ".false.")
        #############################################################
        # run_radiation - logicals previously missing from metadata #
        #############################################################
        self.add_setting(config, ['namelist:run_radiation','i_cloud_representation'], '2')
        self.add_setting(config, ['namelist:run_radiation','i_cloud_representation_2'], '2')

        return config, self.reports


class vn90_t6032(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6032 by James Manners."""

    BEFORE_TAG = "vn9.0_t5694"
    AFTER_TAG = "vn9.0_t6032"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Remove namelists r2lwncal and r2swncal
        namelsts = self.get_setting_value(config, ["file:NAMELIST","source"],)
        if (namelsts): 
            namelsts=namelsts.replace("namelist:r2lwncal", "")
            namelsts=namelsts.replace("namelist:r2swncal", "")
            self.change_setting_value(config, ["file:NAMELIST","source"], 
                                                             namelsts)
        namelsts = self.get_setting_value(config, ["file:CNTLATM","source"],)
        if (namelsts): 
            namelsts=namelsts.replace("namelist:r2lwncal", "")
            namelsts=namelsts.replace("namelist:r2swncal", "")
            self.change_setting_value(config, ["file:CNTLATM","source"], 
                                                             namelsts)
        self.remove_setting(config, ["namelist:r2lwncal"])
        self.remove_setting(config, ["namelist:r2swncal"])

        # Remove paths from spectral file names
        # Should also remove any trailing/leading quotes and whitespace
        Items2Split = [ ("namelist:r2lwclnl","spectral_file_lw" ) \
                      , ("namelist:r2lwclnl","spectral_file_lw2") \
                      , ("namelist:r2swclnl","spectral_file_sw" ) \
                      , ("namelist:r2swclnl","spectral_file_sw2") ]

        for namelist, item in Items2Split:
            # Get the namelist value
            sp_file = self.get_setting_value(config, [namelist,item])

            # Strip off any quotes/whitespace
            sp_file = sp_file.replace("'","")
            sp_file = sp_file.replace('"',"")

            # Get the filename (i.e. string after last '/')
            sp_file = sp_file.strip().split('/')[-1]

            # Fortran expects string to be quoted so finally
            # place the extracted string in quotes
            sp_file = "'%s'"%(sp_file)
            self.change_setting_value(config, [namelist,item], sp_file)


        # Add namelist items previously missing from metadata
        self.add_setting(config, ['namelist:r2swclnl','i_cnv_water_sw2'], '5')
        self.add_setting(config, ['namelist:run_radiation','hfc125mmr'], '0.0')
        self.add_setting(config, ['namelist:run_radiation','i_fsd_2'], '0')
        self.add_setting(config, ['namelist:run_radiation','i_cloud_representation'], '2')
        self.add_setting(config, ['namelist:run_radiation','i_cloud_representation_2'], '2')
        self.add_setting(config, ['namelist:run_radiation','i_inhom'], '0')
        self.add_setting(config, ['namelist:run_radiation','i_inhom_2'], '0')
        self.add_setting(config, ['namelist:run_radiation','i_overlap'], '0')
        self.add_setting(config, ['namelist:run_radiation','i_overlap_2'], '0')
        self.add_setting(config, ['namelist:run_radiation','l_quad_t_coast'], '.false.')
        self.add_setting(config, ['namelist:run_radiation','l_rad_snow_emis'], '.false.')
        self.add_setting(config, ['namelist:run_radiation','l_t_land_nosnow'], '.false.')
        self.add_setting(config, ['namelist:run_radiation','l_t_rad_solid'], '.false.')

        return config, self.reports

class vn90_t6093(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6093 by Steve Mullerworth."""
    BEFORE_TAG = "vn9.0_t6032"
    AFTER_TAG = "vn9.0_t6093"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Input your macro commands here
        self.remove_setting(config, ["namelist:nlstwritedata"])

        # Remove reference to namelist from UM jobs
        values = self.get_setting_value(config, ["file:NAMELIST", "source"])
        if values:
            items = shlex.split(self.get_setting_value(config, ["file:NAMELIST", "source"]))
            items = [item for item in items if item != "namelist:nlstwritedata"]
            new_setting = RosePopener.list_to_shell_str(items)
            self.change_setting_value(config, ["file:NAMELIST", "source"], new_setting)

        # Remove reference to namelist from SCM jobs
        values = self.get_setting_value(config, ["file:CNTLALL", "source"])
        if values:
            items = shlex.split(self.get_setting_value(config, ["file:CNTLALL", "source"]))
            items = [item for item in items if item != "namelist:nlstwritedata"]
            new_setting = RosePopener.list_to_shell_str(items)
            self.change_setting_value(config, ["file:CNTLALL", "source"], new_setting)

        return config, self.reports




class vn90_t5441(rose.upgrade.MacroUpgrade):
    """Upgrade macro for ticket #5441 by Matt Pryor."""
    
    BEFORE_TAG = "vn9.0_t6093"
    AFTER_TAG = "vn9.0_t5441"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        
        config = self.upgrade_hydrology(config)
        config = self.upgrade_radiation(config)
        config = self.upgrade_sea_seaice(config)
        config = self.upgrade_snow(config)
        config = self.upgrade_soil(config)
        config = self.upgrade_surface(config)
        config = self.upgrade_surface_types(config)
        config = self.upgrade_vegetation(config)
        
        self.remove_setting(config, ["namelist:jules_switches"])
        
        # Insert the new namelists into the correct place in NAMELIST and SHARED in the correct order
        # Basically, we replace the JULES_SWITCHES namelist with the new JULES namelists
        # We do this here instead of in individual macros to ensure we have the correct order
        for file in ["NAMELIST", "SHARED"]:
            source = self.get_setting_value(config, ["file:%s" % file, "source"])
            if source:
                # Remove jules_vegetation from its current position first
                source = re.sub(r'namelist:jules_vegetation', r'', source)
                source = re.sub(r'namelist:jules_switches',
                                r'namelist:jules_surface_types namelist:jules_surface namelist:jules_radiation namelist:jules_hydrology namelist:jules_sea_seaice namelist:jules_soil namelist:jules_vegetation namelist:jules_snow',
                                source)
                self.change_setting_value(config, ["file:%s" % file, "source"], source)
        
        return config, self.reports
    
    
    def upgrade_hydrology(self, config):
        # Ensure that the setting exists for the new namelist
        self.add_setting(config, ["namelist:jules_hydrology"])
        
        # Move values from switches to new hydrology namelist
        for member in ["l_hydrology", "l_top", "l_pdm", "l_baseflow_corr", "dz_pdm", "b_pdm"]:
            value = self.get_setting_value(config, ["namelist:jules_switches",member])
            self.remove_setting(config, ["namelist:jules_switches",member])
            if value:
                self.add_setting(config, ["namelist:jules_hydrology",member], value)

        ## ti_max, ti_wetl and zw_max are PARAMETERS being made into namelist variables, so
        ## they don't need to be upgraded, just added to the metadata

        return config
    
    
    def upgrade_radiation(self, config):
        # Ensure that the setting exists for the new namelist
        self.add_setting(config, ["namelist:jules_radiation"])
        
        # Move values from switches to new radiation namelist
        for member in ["l_spec_albedo", "l_spec_alb_bs", "l_snow_albedo", "l_albedo_obs",
                       "l_dolr_land_black", "l_spec_sea_alb", "l_sea_alb_var_chl", "i_sea_alb_method"]:
            value = self.get_setting_value(config, ["namelist:jules_switches",member])
            self.remove_setting(config, ["namelist:jules_switches",member])
            if value:
                self.add_setting(config, ["namelist:jules_radiation",member], value)

        return config
    
    
    def upgrade_sea_seaice(self, config):
        # Ensure that the setting exists for the new namelist
        self.add_setting(config, ["namelist:jules_sea_seaice"])
        
        # Move values from switches to new sea and sea-ice namelist
        for member in ["l_tstar_sice_new", "l_ssice_albedo", "l_sice_scattering",
                       "l_sice_meltponds", "l_sice_multilayers", "l_cice_alb",
                       "l_sice_heatflux", "l_ctile", "iseaz0t", "buddy_sea"]:
            value = self.get_setting_value(config, ["namelist:jules_switches",member])
            self.remove_setting(config, ["namelist:jules_switches",member])
            if value:
                self.add_setting(config, ["namelist:jules_sea_seaice",member], value)
            
        # l_sice_hadgem1a is no longer required
        self.remove_setting(config, ["namelist:jules_switches","l_sice_hadgem1a"])

        return config
    
    
    def upgrade_snow(self, config):
        # Ensure that the setting exists for the new namelist
        self.add_setting(config, ["namelist:jules_snow"])
        
        # Move nsmax from atm_sizes to new jules_snow namelist
        nsmax = self.get_setting_value(config, ["namelist:atm_sizes","nsmax"])
        self.remove_setting(config, ["namelist:atm_sizes","nsmax"])
        if nsmax:
            self.add_setting(config, ["namelist:jules_snow","nsmax"], nsmax)
            
        # Move switches from jules_switches to new jules_snow namelist
        for member in ["l_snowdep_surf", "l_rho_snow_corr", "frac_snow_subl_melt"]:
            value = self.get_setting_value(config, ["namelist:jules_switches",member])
            self.remove_setting(config, ["namelist:jules_switches",member])
            if value:
                self.add_setting(config, ["namelist:jules_snow",member], value)
            
        # Copy all items from jules_rad_param to jules_snow and remove jules_rad_param
        for member in ["r0", "rmax", "snow_ggr", "amax", "maskd", "dtland", "kland_numerator"]:
            value = self.get_setting_value(config, ["namelist:jules_rad_param",member])
            if value:
                self.add_setting(config, ["namelist:jules_snow",member], value)
                        
        self.remove_setting(config, ["namelist:jules_rad_param"])
        
        # Copy all items from jules_snow_param to jules_snow and remove jules_snow_param
        for member in ["rho_snow_const", "rho_snow_fresh", "snow_hcon", "snow_hcap",
                       "snowliqcap", "snowinterceptfact", "snowloadlai", "snowunloadfact", "cansnowpft"]:
            value = self.get_setting_value(config, ["namelist:jules_snow_param",member])
            if value:
                self.add_setting(config, ["namelist:jules_snow",member], value)

        # dzsnow_io is renamed to dzsnow in the new namelist
        dzsnow = self.get_setting_value(config, ["namelist:jules_snow_param","dzsnow_io"])
        if dzsnow:
            self.add_setting(config, ["namelist:jules_snow","dzsnow"], dzsnow)
        
        self.remove_setting(config, ["namelist:jules_snow_param"])
        
        # Update the source of the UM NAMELIST, SHARED and CNTLATM files to remove jules_rad_param and jules_snow_param
        for file in ["NAMELIST", "SHARED", "CNTLATM"]:
            source = self.get_setting_value(config, ["file:%s" % file,"source"])
            if source:
                source = re.sub(r'namelist:jules_rad_param', r'', source)
                source = re.sub(r'namelist:jules_snow_param', r'', source)
                self.change_setting_value(config, ["file:%s" % file,"source"], source)

        return config
    
    
    def upgrade_soil(self, config):
        # Ensure that the setting exists for the new namelist
        self.add_setting(config, ["namelist:jules_soil"])
        
        # Move switches from jules_switches to new jules_soil namelist
        for member in ["l_vg_soil", "l_dpsids_dsdz", "l_soil_sat_down", "soilhc_method"]:
            value = self.get_setting_value(config, ["namelist:jules_switches",member])
            self.remove_setting(config, ["namelist:jules_switches",member])
            if value:
                self.add_setting(config, ["namelist:jules_soil",member], value)
                
        # Move values from JULES_SOIL_PARAM to JULES_SOIL and remove it
        for member in ["zsmc", "zst", "confrac", "dzsoil_io"]:
            value = self.get_setting_value(config, ["namelist:jules_soil_param",member])
            if value:
                self.add_setting(config, ["namelist:jules_soil",member], value)
        self.remove_setting(config, ["namelist:jules_soil_param"])
        
        # Move cs_min from JULES_CSMIN to JULES_SOIL and remove it
        value = self.get_setting_value(config, ["namelist:jules_csmin","cs_min"])
        if value:
            self.add_setting(config, ["namelist:jules_soil","cs_min"], value)
        self.remove_setting(config, ["namelist:jules_csmin"])
        
        # Update the source of the UM NAMELIST, SHARED and CNTLATM files to remove jules_soil_param and jules_csmin
        for file in ["NAMELIST", "SHARED", "CNTLATM"]:
            source = self.get_setting_value(config, ["file:%s" % file,"source"])
            if source:
                source = re.sub(r'namelist:jules_csmin', r'', source)
                source = re.sub(r'namelist:jules_soil_param', r'', source)
                self.change_setting_value(config, ["file:%s" % file,"source"], source)
        
        return config
    
    
    def upgrade_surface(self, config):
        # Ensure that the setting exists for the new namelist
        self.add_setting(config, ["namelist:jules_surface"])
        
        # Move switches from jules_switches to new jules_surface namelist
        for member in ["l_flake_model", "l_epot_corr", "l_point_data", "l_aggregate",
                       "l_land_ice_imp", "l_anthrop_heat_src", "i_modiscopt", "all_tiles",
                       "cor_mo_iter", "iscrntdiag", "i_aggregate_opt", "formdrag", 
                       "fd_stab_dep", "isrfexcnvgust", "orog_drag_param"]:
            value = self.get_setting_value(config, ["namelist:jules_switches",member])
            self.remove_setting(config, ["namelist:jules_switches",member])
            if value:
                self.add_setting(config, ["namelist:jules_surface",member], value)
                
        # Move surface related values from JULES_SURF_PARAM to JULES_SURFACE
        for member in ["hleaf", "hwood", "beta1", "beta2", "fwe_c3", "fwe_c4", "q10_leaf",
                       "q10_soil", "kaps", "kaps_roth"]:
            value = self.get_setting_value(config, ["namelist:jules_surf_param",member])
            if value:
                self.add_setting(config, ["namelist:jules_surface",member], value)
        
        # Move sea and sea-ice related values from JULES_SURF_PARAM to JULES_SEA_SEAICE
        for member in ["z0miz", "z0sice", "z0h_z0m_miz", "z0h_z0m_sice", "emis_sea",
                       "emis_sice", "kappai", "kappai_snow", "kappa_seasurf",
                       "seasalinityfactor", "charnock"]:
            value = self.get_setting_value(config, ["namelist:jules_surf_param",member])
            if value:
                self.add_setting(config, ["namelist:jules_sea_seaice",member], value)
        
        # Remove JULES_SURF_PARAM
        self.remove_setting(config, ["namelist:jules_surf_param"])
        
        # Update the source of the UM NAMELIST, SHARED and CNTLATM files to remove jules_surf_param
        for file in ["NAMELIST", "SHARED", "CNTLATM"]:
            source = self.get_setting_value(config, ["file:%s" % file,"source"])
            if source:
                source = re.sub(r'namelist:jules_surf_param', r'', source)
                self.change_setting_value(config, ["file:%s" % file,"source"], source)
        
        return config
    
    
    def upgrade_surface_types(self, config):
        # Ensure that the setting exists for the new namelist
        self.add_setting(config, ["namelist:jules_surface_types"])
        
        # Move values from nstypes to new jules_surface_types namelist
        for member in ["npft", "nnvg", "urban", "lake", "soil", "ice"]:
            value = self.get_setting_value(config, ["namelist:jules_nstypes",member])
            self.remove_setting(config, ["namelist:jules_nstypes",member])
            if value:
                self.add_setting(config, ["namelist:jules_surface_types",member], value)
                
        # Remove the now redundant nstypes namelist
        self.remove_setting(config, ["namelist:jules_nstypes"])
                
        # Update the source of the UM NAMELIST and SIZES files to remove jules_nstypes
        for file in ["NAMELIST", "SIZES"]:
            source = self.get_setting_value(config, ["file:%s" % file,"source"])
            if source:
                source = re.sub(r'namelist:jules_nstypes', r'', source)
                self.change_setting_value(config, ["file:%s" % file,"source"], source)

        return config
    
    
    def upgrade_vegetation(self, config):
        # Ensure that the setting exists for the new namelist
        self.add_setting(config, ["namelist:jules_vegetation"])
        
        # Replace i_veg_vn with a suitable setting for l_triffid
        i_veg_vn = self.get_setting_value(config, ["namelist:jules_vegetation","i_veg_vn"])
        self.remove_setting(config, ["namelist:jules_vegetation","i_veg_vn"])
        l_triffid = ".true." if (i_veg_vn and int(i_veg_vn) == 2) else ".false."
        self.add_setting(config, ["namelist:jules_vegetation","l_triffid"], l_triffid)
        
        # Remove l_bvoc_emis and l_o3_damage from switches as they have no effect in the UM
        self.remove_setting(config, ["namelist:jules_switches","l_bvoc_emis"])
        self.remove_setting(config, ["namelist:jules_switches","l_o3_damage"])
        
        # Move other items from switches to jules_vegetation
        for member in ["can_model", "can_rad_mod", "ilayers"]:
            value = self.get_setting_value(config, ["namelist:jules_switches",member])
            self.remove_setting(config, ["namelist:jules_switches",member])
            if value:
                self.add_setting(config, ["namelist:jules_vegetation",member], value)
        
        # Copy any items from jules_seed and remove it
        for member in ["frac_min", "frac_seed"]:
            value = self.get_setting_value(config, ["namelist:jules_seed",member])
            if value:
                self.add_setting(config, ["namelist:jules_vegetation",member], value)
        self.remove_setting(config, ["namelist:jules_seed"])
        
        # Copy pow from jules_sigm and remove it
        value = self.get_setting_value(config, ["namelist:jules_sigm","pow"])
        if value:
            self.add_setting(config, ["namelist:jules_vegetation","pow"], value)
        self.remove_setting(config, ["namelist:jules_sigm"])
        
        # Update the source of the UM NAMELIST, SHARED and CNTLATM files to remove jules_seed and jules_sigm
        for file in ["NAMELIST", "SHARED", "CNTLATM"]:
            source = self.get_setting_value(config, ["file:%s" % file,"source"])
            if source:
                source = re.sub(r'namelist:jules_seed', r'', source)
                source = re.sub(r'namelist:jules_sigm', r'', source)
                self.change_setting_value(config, ["file:%s" % file,"source"], source)
        
        return config



class vn90_t6209(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6209 by Glenn Greed."""

    BEFORE_TAG = "vn9.0_t5441"
    AFTER_TAG = "vn9.0_t6209"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.remove_setting(config,["namelist:umsections"],
                            info="Cxx no longer runtime inputs")
        # remove this namelist from the NAMELIST file
        namelist_source = self.get_setting_value(config, ["file:NAMELIST","source"])
        if namelist_source:
            namelist_source = re.sub(r'namelist:umsections ','',namelist_source)
            self.change_setting_value(config, ["file:NAMELIST","source"],namelist_source)
        # remove this namelist from the SIZES file
        shared_source = self.get_setting_value(config, ["file:SIZES","source"])
        if shared_source:
            shared_source = re.sub(r' namelist:umsections',' ',shared_source)
            self.change_setting_value(config, ["file:SIZES","source"],shared_source)
        return config, self.reports    


class vn90_t5551(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5551 by Ruth Lewis."""

    BEFORE_TAG = "vn9.0_t6209"
    AFTER_TAG = "vn9.0_t5551"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:temp_fixes","l_emis_ssi_full"], ".false.") 
        return config, self.reports

class vn90_t5945(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5945 by Ian Boutle."""

    BEFORE_TAG = "vn9.0_t5551"
    AFTER_TAG = "vn9.0_t5945"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:jules_hydrology","l_var_rainfrac"],".false.")
        return config, self.reports

class vn90_t6121(rose.upgrade.MacroUpgrade):

    """Upgrade macro for tickets #6120 and #6121 by Adrian Lock."""

    BEFORE_TAG = "vn9.0_t5945"
    AFTER_TAG = "vn9.0_t6121"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here
        self.add_setting(config, ["namelist:run_cloud","forced_cu"],"0")
        self.add_setting(config, ["namelist:run_bl","kprof_cu"],"0")
        self.add_setting(config, ["namelist:run_bl","bl_res_inv"],"0")
        self.add_setting(config, ["namelist:run_bl","a_ent_shr_nml"],"5.0")
        self.add_setting(config, ["namelist:run_convection","cldbase_opt_sh"],"0")
        return config, self.reports

class vn90_t6022(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6022 by hadzc Claudio Sanchez."""

    BEFORE_TAG = "vn9.0_t6121"
    AFTER_TAG = "vn9.0_t6022"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_stochastic", "l_spt"], ".false.")
        self.add_setting(config, ["namelist:run_stochastic", "l_spt_cfl"], ".false.")
        self.add_setting(config, ["namelist:run_stochastic", "l_spt_rain"], ".false.")
        self.add_setting(config, ["namelist:run_stochastic", "l_spt_rad"], ".false.")
        self.add_setting(config, ["namelist:run_stochastic", "l_spt_gwd"], ".false.")
        self.add_setting(config, ["namelist:run_stochastic", "l_spt_conv"], ".false.")
        self.add_setting(config, ["namelist:run_stochastic", "l_spt_conv_mom"], ".false.")
        self.add_setting(config, ["namelist:run_stochastic", "l_spt_qcons"], ".false.")
        self.add_setting(config, ["namelist:run_stochastic", "spt_bot_tap_lev"], "9")
        self.add_setting(config, ["namelist:run_stochastic", "spt_botlev"], "15")
        self.add_setting(config, ["namelist:run_stochastic", "spt_toplev"], "41")
        self.add_setting(config, ["namelist:run_stochastic", "spt_top_tap_lev"], "45")
        self.add_setting(config, ["namelist:run_stochastic", "nsmooth_spt"], "3")
        self.add_setting(config, ["namelist:run_stochastic", "rain_std"], "1.73")
        self.add_setting(config, ["namelist:run_stochastic", "rad_std"], "1.73")
        self.add_setting(config, ["namelist:run_stochastic", "gwd_std"], "1.5")
        self.add_setting(config, ["namelist:run_stochastic", "conv_std"], "1.73")
        return config, self.reports

class vn90_t4471(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #XXXX by <author>."""

    BEFORE_TAG = "vn9.0_t6022"
    AFTER_TAG = "vn9.0_t4471"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Introduce two items if UKCA aerosol chemistry 
        # is selected but CLASSIC is Off

        l_ukca_achem = self.get_setting_value( config,
                     ["namelist:run_ukca", "l_ukca_chem_aero"])

        if l_ukca_achem  == ".true.":
                        # Check if Hi-level emissions are active 
           l_ClassSO2hi = self.get_setting_value( config, 
                     ["namelist:run_aerosol", "l_so2_hilem"])
           if l_ClassSO2hi == ".false.":
                 self.add_setting(config, 
                          ['namelist:run_ukca', 'i_so2_hi_level'],
                           value = '8')
    
         # Add DMS flux option   
           l_ClassDMS = self.get_setting_value( config, 
                     ["namelist:run_aerosol", "l_sulpc_dms"])
           if l_ClassDMS == ".false.":
                 self.add_setting(config, 
                          ['namelist:run_ukca', 'i_ukca_dms_flux'],
                          value = '1')
        else:
     
           self.add_setting(config,            
                          ['namelist:run_ukca', 'i_so2_hi_level'],
                           value = '0')
      
           self.add_setting(config,            
                          ['namelist:run_ukca', 'i_ukca_dms_flux'],
                           value = '0')

        return config, self.reports
    
class vn90_t6184(rose.upgrade.MacroUpgrade):
    """Upgrade macro for ticket #6184 by Steve Wardle."""
    
    BEFORE_TAG = "vn9.0_t4471"
    AFTER_TAG = "vn9.0_t6184"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Move nice to JULES seaice namelist
        nice = self.get_setting_value(config, ["namelist:atm_sizes","nice"])
        if nice:
              self.add_setting(config, ["namelist:jules_sea_seaice","nice"], nice)
        else:
              self.add_setting(config, ["namelist:jules_sea_seaice","nice"], 1)

        # Move nancil_lookupsa to the ancilcta namelist
        nancil_lookupsa = self.get_setting_value(config, ["namelist:atm_sizes","nancil_lookupsa"])
        if nancil_lookupsa:
              self.add_setting(config, ["namelist:ancilcta","nancil_lookupsa"], nancil_lookupsa)
        else:
              self.add_setting(config, ["namelist:ancilcta","nancil_lookupsa"], 0)

        # Move nrim_timesa to the lbc_options namelist
        nrim_timesa = self.get_setting_value(config, ["namelist:atm_sizes","nrim_timesa"])
        if nrim_timesa:
              self.add_setting(config, ["namelist:lbc_options","nrim_timesa"],nrim_timesa)
        else:
              self.add_setting(config, ["namelist:lbc_options","nrim_timesa"], 1)

        # Atm_sizes is no longer needed
        self.remove_setting(config, ["namelist:atm_sizes"])

        for nlfile in ["NAMELIST","SIZES"]:
              entry = "file:{0:s}".format(nlfile)
              source = self.get_setting_value(config, [entry,"source"])
              if source:
                    source = re.sub(r'namelist:atm_sizes\s*',r'', source)
                    self.change_setting_value(config, [entry,"source"], source)

        return config, self.reports
        
    
class vn90_t6211(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6211 by Joe Mancell"""
    BEFORE_TAG = "vn9.0_t6184"
    AFTER_TAG = "vn9.0_t6211"

    def upgrade(self, config, meta_config=None):
          """Upgrade a UM runtime app configuration."""
          # Retire inputs
          self.remove_setting(config, ["namelist:nlstcgen","llboutim"])
          self.remove_setting(config, ["namelist:nlstcgen","jobrel_stepim"])
          self.remove_setting(config, ["namelist:nlstcgen","jobrel_offsetim"])
          self.remove_setting(config, ["namelist:nlstcall","model_status"])
          dumpfreqim = self.get_setting_value(config, ["namelist:nlstcgen",
                                                       "dumpfreqim"])
          # Get irregular dump time array and convert to integer list
          dumptimesim = map(int, 
                            self.get_setting_value(config, 
                                                   ["namelist:nlstcgen",
                                                    "dumptimesim"]).split(","))
          dumptimesim_mask = [value > 0 for value in dumptimesim]
          if int(dumpfreqim) > 0 and not any(dumptimesim_mask):
                # Regular dumping
                i_dump_output = "2"
          elif int(dumpfreqim) == 0 and any(dumptimesim_mask):
                # Irregular dumping
                i_dump_output = "3"
          elif int(dumpfreqim) == 0 and not any(dumptimesim_mask):
                # No dumping
                i_dump_output = "1"
          else:
                i_dump_output = "1"
                dump_warn = ("Namelist values indicate UM run with both " +
                             "irregular and regular dumping. Dumping has been" +
                             " turned off on upgrade. Please configure " +
                             "dumping manually.")
                self.add_report("namelist:nlstcgen", "i_dump_freq", None,
                            info=dump_warn,is_warning=True)
          self.add_setting(config,  ["namelist:nlstcgen",
                                     "i_dump_output"], i_dump_output)
          # Get meaning period length array and convert to integer list
          meanfreqim = map(int, 
                           self.get_setting_value(config, 
                                                  ["namelist:nlstcgen",
                                                   "meanfreqim"]).split(","))
          meanfreqim_mask = [value > 0 for value in meanfreqim]
          if any(meanfreqim_mask):
                l_meaning_sequence = ".true."
          else:
                l_meaning_sequence = ".false."
          self.add_setting(config,  ["namelist:nlstcgen",
                                     "l_meaning_sequence"], l_meaning_sequence)

          # Add sensible default for new compulsory=true items
          surf_hgt_io = r'0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
          l_moruses_albedo       = ".false."
          l_moruses_emissivity   = ".false."
          l_moruses_macdonald    = ".false."
          l_moruses_rough        = ".false."
          l_moruses_storage      = ".false."
          l_moruses_storage_thin = ".false."
          l_urban2t              = ".false."
          l_urban_empirical      = ".false."
          a_assim_start_min      = "0"
          a_assim_end_min        = "0"
          io_timing              = "0"
          itemc                  = "0"                 
          sctnc                  = "0"    
          lev1                   = "0"    
          lev2                   = "0"    
          row1                   = "0"    
          row2                   = "0"    
          col1                   = "0"    
          col2                   = "0"    
          self.add_setting(config,  ["namelist:jules_elevate",
                                     "surf_hgt_io"],
                                      surf_hgt_io)
          self.add_setting(config,  ["namelist:urban_switches",
                                     "l_moruses_albedo"],
                                      l_moruses_albedo)
          self.add_setting(config,  ["namelist:urban_switches",
                                     "l_moruses_emissivity"],
                                      l_moruses_emissivity)
          self.add_setting(config,  ["namelist:urban_switches",
                                     "l_moruses_macdonald"],
                                      l_moruses_macdonald )
          self.add_setting(config,  ["namelist:urban_switches",
                                     "l_moruses_rough"],
                                      l_moruses_rough)
          self.add_setting(config,  ["namelist:urban_switches",
                                     "l_moruses_storage"],
                                      l_moruses_storage)
          self.add_setting(config,  ["namelist:urban_switches",
                                     "l_moruses_storage_thin"],
                                      l_moruses_storage_thin)
          self.add_setting(config,  ["namelist:urban_switches",
                                     "l_urban2t"],
                                      l_urban2t)
          self.add_setting(config,  ["namelist:urban_switches",
                                     "l_urban_empirical"],
                                      l_urban_empirical)
          self.add_setting(config,  ["namelist:nlstcatm",
                                     "a_assim_start_min"],
                                      a_assim_start_min)
          self.add_setting(config,  ["namelist:nlstcatm",
                                     "a_assim_end_min"],
                                      a_assim_end_min)
          self.add_setting(config,  ["namelist:io_control",
                                     "io_timing"],
                                      io_timing)
          # Check if the app file already contains a trans namelist
          if any([re.search(r'^namelist:trans', s) for s in config.value.keys()]):
                recona_source = self.get_setting_value(config, 
                                                       ["file:RECONA","source"])
                if recona_source:
                      # Convert single trans namelist to use duplicate notation
                      if not re.search('namelist:trans\(:\)', recona_source):
                            recona_source = recona_source.replace("namelist:trans",
                                                                  "namelist:trans(:)")
                      # Allow trans namelist to be optional
                      recona_source = recona_source.replace("namelist:trans(:)",
                                                            "(namelist:trans(:))")
                      self.change_setting_value(config, ["file:RECONA","source"],
                                                recona_source)
          else:
                # Add defaults
                self.add_setting(config,  ["namelist:trans(1)", "itemc"], itemc)
                self.add_setting(config,  ["namelist:trans(1)", "sctnc"], sctnc)
                self.add_setting(config,  ["namelist:trans(1)", "lev1"], lev1)
                self.add_setting(config,  ["namelist:trans(1)", "lev2"], lev2)
                self.add_setting(config,  ["namelist:trans(1)", "row1"], row1)
                self.add_setting(config,  ["namelist:trans(1)", "row2"], row2)
                self.add_setting(config,  ["namelist:trans(1)", "col1"], col1)
                self.add_setting(config,  ["namelist:trans(1)", "col2"], col2)
                recona_source = self.get_setting_value(config, 
                                                       ["file:RECONA","source"])
                if recona_source:
                      # Add namelist to the RECON file
                      recona_source = recona_source + " (namelist:trans(:))" 
                      self.change_setting_value(config, ["file:RECONA","source"],
                                                recona_source)
          # Deal with single trans namelist and convert to trans(1)
          if any([re.search(r'^namelist:trans$', s) for s in config.value.keys()]):
                itemc = self.get_setting_value(config, ["namelist:trans", "itemc"])
                sctnc = self.get_setting_value(config, ["namelist:trans", "sctnc"])
                lev1 = self.get_setting_value(config, ["namelist:trans", "lev1"])
                lev2 = self.get_setting_value(config, ["namelist:trans", "lev2"])
                row1 = self.get_setting_value(config, ["namelist:trans", "row1"])
                row2 = self.get_setting_value(config, ["namelist:trans", "row2"])
                col1 = self.get_setting_value(config, ["namelist:trans", "col1"])
                col2 = self.get_setting_value(config, ["namelist:trans", "col2"])
                # Add to new trans(1) section
                self.add_setting(config,  ["namelist:trans(1)", "itemc"], itemc)
                self.add_setting(config,  ["namelist:trans(1)", "sctnc"], sctnc)
                self.add_setting(config,  ["namelist:trans(1)", "lev1"], lev1)
                self.add_setting(config,  ["namelist:trans(1)", "lev2"], lev2)
                self.add_setting(config,  ["namelist:trans(1)", "row1"], row1)
                self.add_setting(config,  ["namelist:trans(1)", "row2"], row2)
                self.add_setting(config,  ["namelist:trans(1)", "col1"], col1)
                self.add_setting(config,  ["namelist:trans(1)", "col2"], col2)
                # Remove old trans section
                self.remove_setting(config, ["namelist:trans"])
          return config, self.reports

class vn90_t4300(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4300 by hadcj."""

    BEFORE_TAG = "vn9.0_t6211"
    AFTER_TAG = "vn9.0_t4300"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Input your macro commands here

        # Introduce two items if UKCA offline oxidants
        # chemistry is selected

        i_ukca_chemx = self.get_setting_value( config, 
                     ["namelist:run_ukca", "i_ukca_chem"]) 

        if i_ukca_chemx  == "54":
                     # check if offline oxidants chemistry selected
                 self.add_setting(config, 
                         ['namelist:run_ukca', 'ukca_offline_dir'], 
                           value = '')

                 self.add_setting(config, 
                         ['namelist:run_ukca', 'ukca_offline_files'], 
                           value = '')
        else:
                 self.add_setting(config,             
                         ['namelist:run_ukca', 'ukca_offline_dir'], 
                           state=config.STATE_SYST_IGNORED) 
       
                 self.add_setting(config,             
                         ['namelist:run_ukca', 'ukca_offline_files'], 
                           state=config.STATE_SYST_IGNORED) 
       
        return config, self.reports


class vn90_t6262(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6262 by Joe Mancell"""

    BEFORE_TAG = "vn9.0_t4300"
    AFTER_TAG = "vn9.0_t6262"
    def upgrade(self, config, meta_config=None):
          """Upgrade a UM runtime app configuration."""
          # Input your macro commands here
          # Currently all app have dumping period specified in timesteps
          self.add_setting(config,             
                           ['namelist:nlstcgen', 'dump_frequency_units'], "3")

          return config, self.reports

class vn90_t6246(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6246 by Joe Mancell"""
    BEFORE_TAG = "vn9.0_t6262"
    AFTER_TAG = "vn9.0_t6246"

    def upgrade(self, config, meta_config=None):
          """Upgrade a UM runtime app configuration."""
          items_defaults = {"ancilfilename" : "''",
                            "interval" : "1",
                            "period" : "1",
                            "update_anc" : '.false.',
                            "source" : "1",
                            "user_prog_rconst" : "0.0",
                            "user_prog_ancil_stash_req" : "",
                            "domain" : "1",
                            "stash_req" : ""}
          for key, value in items_defaults.iteritems():
                for obi in config.get_value():
                      # Loop through items name lists and add missing items
                      if re.search(r'namelist:item', obi):
                            self.add_setting(config, [obi, str(key)], str(value))
                            
          return config, self.reports

class vn90_t4634(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #4634 by Maggie Hendry."""

    BEFORE_TAG = "vn9.0_t6246"
    AFTER_TAG = "vn9.0_t4634"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Defaults removed for tile number so have to be added to namelist
        self.add_setting(config, ["namelist:jules_surface_types", "brd_leaf"], "1")
        self.add_setting(config, ["namelist:jules_surface_types", "ndl_leaf"], "2")
        self.add_setting(config, ["namelist:jules_surface_types", "c3_grass"], "3")
        self.add_setting(config, ["namelist:jules_surface_types", "c4_grass"], "4")
        self.add_setting(config, ["namelist:jules_surface_types", "shrub"],    "5")
        self.add_setting(config, ["namelist:jules_surface_types", "urban"],    "6")
        self.add_setting(config, ["namelist:jules_surface_types", "lake"],     "7")
        self.add_setting(config, ["namelist:jules_surface_types", "soil"],     "8")
        self.add_setting(config, ["namelist:jules_surface_types", "ice"],      "9")

        return config, self.reports


class vn90_t5903(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5903 by Paul Cresswell.
       Move INITFILEENV variables into namelists."""

    BEFORE_TAG = "vn9.0_t4634"
    AFTER_TAG = "vn9.0_t5903"

    def upgrade(self, config, meta_config=None):

        scm_app = False
        if config.get(["namelist:scm_cntl"]) is not None:
            scm_app = True

        # Insert the new namelist reads into the right files:
        new_namelist = 'namelist:nlcfiles'

        # NAMELIST shouldn't exist in SCM apps anyway, but just in case:
        if not scm_app:
            if config.get(["file:NAMELIST"]):
                src_namelist = (self.get_setting_value(config,
                                ["file:NAMELIST","source"])).split()
                src_namelist.insert(0, new_namelist)
                self.change_setting_value(config, ["file:NAMELIST", "source"],
                                          " ".join(src_namelist))

        # nlcfiles is compulsory in non-SCM apps but should be trigger-ignored
        # out of SCM apps; to appear consistent we add it to SCM apps as an
        # optional namelist.
        if config.get(["file:SHARED"]):
            src_shared = (self.get_setting_value(config,
                            ["file:SHARED","source"])).split()
            if scm_app:
                new_namelist = '(' + new_namelist + ')'
            src_shared.insert(0, new_namelist)
            self.change_setting_value(config,
                           ["file:SHARED", "source"], " ".join(src_shared))

        # Make a dict of variables to keep along with their target namelists
        # and values, in the form { env : [namelist, value] }

        var_list = [
                     ('AINITIAL', 'recon'),
                     ('ALABCIN1', 'nlcfiles'),
                     ('ALABCIN2', 'nlcfiles'),
                     ('ASTART',   'nlcfiles'),
                     ('ATMANL',   'nlcfiles'),
                     ('CURNTOUT', 'nlcfiles'),
                     ('CXBKGERR', 'nlcfiles'),
                     ('FLXCROUT', 'nlcfiles'),
                     ('FOAMOUT1', 'nlcfiles'),
                     ('FOAMOUT2', 'nlcfiles'),
                     ('IAU_inc',  'nlcfiles'),
                     ('ICEFOUT',  'nlcfiles'),
                     ('IDEALISE', 'nlcfiles'),
                     ('MOSOUT',   'nlcfiles'),
                     ('OBS01',    'nlcfiles'),
                     ('OBS02',    'nlcfiles'),
                     ('OBS03',    'nlcfiles'),
                     ('OBS04',    'nlcfiles'),
                     ('OBS05',    'nlcfiles'),
                     ('OCNANL',   'nlcfiles'),
                     ('PP0',      'nlcfiles'),
                     ('PP1',      'nlcfiles'),
                     ('PP2',      'nlcfiles'),
                     ('PP3',      'nlcfiles'),
                     ('PP4',      'nlcfiles'),
                     ('PP5',      'nlcfiles'),
                     ('PP6',      'nlcfiles'),
                     ('PP7',      'nlcfiles'),
                     ('PP8',      'nlcfiles'),
                     ('PP9',      'nlcfiles'),
                     ('PP10',     'nlcfiles'),
                     ('PPMBC',    'nlcfiles'),
                     ('PPSCREEN', 'nlcfiles'),
                     ('PPSMC',    'nlcfiles'),
                     ('PPVAR',    'nlcfiles'),
                     ('RFMOUT',   'nlcfiles'),
                     ('RPSEED',   'nlcfiles'),
                     ('SICEOUT',  'nlcfiles'),
                     ('SSTOUT',   'nlcfiles'),
                     ('SURGEOU1', 'nlcfiles'),
                     ('SURGEOUT', 'nlcfiles'),
                     ('TRANSP',   'recon'),
                     ('UARSOUT1', 'nlcfiles'),
                     ('UARSOUT2', 'nlcfiles'),
                     ('UKCAACLW', 'nlcfiles'),
                     ('UKCAACSW', 'nlcfiles'),
                     ('UKCAANLW', 'nlcfiles'),
                     ('UKCAANSW', 'nlcfiles'),
                     ('UKCACRLW', 'nlcfiles'),
                     ('UKCACRSW', 'nlcfiles'),
                     ('UKCAPREC', 'nlcfiles'),
                     ('VAR_GRID', 'nlsizes'),
                     ('VERT_LEV', 'nlsizes'),
                     ('WFOUT',    'nlcfiles'),
                   ]

        vars = defaultdict(list)
        for key,value in var_list:
            vars[key].append(value)

        # Need to check the env section of the app before the file itself
        for envvar in vars:
            value = self.get_setting_value(config, ["env", envvar])
            if value:
                vars[envvar].append(value)
                self.remove_setting(config, ["env", envvar])

        # Now scan the file to get the value of this item
        initfileenv = 'file/INITFILEENV'

        # Check whether INITFILEENV is present.
        file_exists = False
        if os.path.isfile(initfileenv):

            file_exists = True
            file = open(initfileenv, 'r')
            for line in file:
                for envvar in vars:
                    searchObj = re.search(r'^\s*export\s*%s=(.*)\s*#?' % envvar,
                                          line)
                    if searchObj:
                        # Check the user doesn't set var=$var in the file and
                        # pass a value from the app; in this case we want to
                        # keep the value from the app. (Values passed directly
                        # from the suite to the file, such as rose-stem does, 
                        # are fine.) Otherwise, take the value from the file.
                        if (len(vars[envvar]) > 1
                        and searchObj.group(1) == "$%s" % envvar):
                            pass
                        else:
                            vars[envvar].append((searchObj.group(1)))

        # Tidy the resulting values, creating an entry if none was found,
        # then add each one as a namelist item.
        # Always work with the last entry in each value-list because it may 
        # have been appended from both the app and the file,
        # in which case we want the file's value.
        for envvar in vars:
            if len(vars[envvar]) > 1:
                # Remove variable and remove padding:
                vars[envvar][-1] = re.sub(r'\$PREFIXW(?=\W)', '', 
                                          vars[envvar][-1]).strip()
            else:
                # Value did not exist in either app or file;
                # create a blank entry.
                vars[envvar].append('')

            # Give blank values an appropriate default.
            # We make a special exception for ASTART, which is required by
            # every recon and/or atmos app; if it's not set to anything
            # assume it's being passed down from the suite level.
            if vars[envvar][-1] == '':
                if envvar == 'ASTART':
                    vars[envvar][-1] = '$ASTART'
                else:
                    vars[envvar][-1] = 'unset'

            # Add quotes:
            value = "'" + vars[envvar][-1] + "'"

            self.add_setting(config, ["namelist:%s" % vars[envvar][0],
                                      envvar.lower()], value)

        # User may now delete INITFILEENV. We cannot do this in the macro
        # because it will happen even if the user rejects the upgrade.
        if file_exists:
           delete_msg = """
        !!!!!  file/INITFILEENV is no longer used by the UM.  !!!!! 
        !!!!! Please delete INITFILEENV from all UM 9.1 apps. !!!!!
                        """
           self.add_report(info=delete_msg, is_warning=True)
 
        return config, self.reports

class vn90_t6314a(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #6314 for Joe Mancell"""

    BEFORE_TAG = "vn9.0_t5903"
    AFTER_TAG = "vn9.0_t6314"
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


class vn90_t1389(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1389 (SRS) by Paul Cresswell."""

    BEFORE_TAG = "vn9.0_t6314"
    AFTER_TAG = "vn9.0_t1389"

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


class vn91_t6270(rose.upgrade.MacroUpgrade):
      
    BEFORE_TAG = "vn9.0_t1389"
    AFTER_TAG = "vn9.1"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        if self.get_setting_value(config, ["env", "VN"]):
            self.change_setting_value(config, ["env", "VN"],
                                      "9.1")
            # Update STASHMSTR
            stashmaster_path = self.get_setting_value(config, ["env", "STASHMSTR"])
            if stashmaster_path:
                stashmaster_path = re.sub("vn\d+\.\d+\/ctldata", 
                                          "vn9.1/ctldata", 
                                          stashmaster_path)
                self.change_setting_value(config, ['env', 'STASHMSTR'], 
                                        stashmaster_path)
        self.remove_setting(config, ["env", "ANCIL_VERSIONS"])
        self.remove_setting(config, ["env", "UM_ANCIL_FILENAMES"])

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



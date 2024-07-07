import rose.upgrade
import re
import os
import sys

class UpgradeError(Exception):

      """Exception created when an upgrade fails."""
      
      def __init__(self, msg):
          self.msg = msg
      
      def __repr__(self):
          sys.tracebacklimit = 0
          return self.msg
          
      __str__ = __repr__


# Please replace XXXX with your ticket number
class vn86_t5702(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5702 by Joe Mancell."""
    
    BEFORE_TAG = "vn8.6"
    AFTER_TAG = "vn8.6_t5702"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add your upgrade macro commands here
        # Get values of old logicals
        grib = self.get_setting_value(
                config, ["namelist:recon", "grib"])
        grib2ff = self.get_setting_value(
            config, ["namelist:recon", "grib2ff"])
        # Remove old logicals
        self.remove_setting(config, ["namelist:recon", "grib"])
        self.remove_setting(config, ["namelist:recon", "grib2ff"])
        # Calculate new integer based on old logicals
        if grib2ff == ".false." and grib == ".false.":
            input_dump_type_local = 1
        elif grib2ff == ".false." and grib == ".true.":
            input_dump_type_local = 2
        elif grib2ff == ".true." and grib == ".false.":
            input_dump_type_local = 3
        else:
            self.add_report("namelist:recon", "input_dump_type", None,
                            info="UPGRADE NOT POSSIBLE for input_dump_type, grib and grib2ff settings clash", 
                            is_warning=True)
            input_dump_type_local = None
        if input_dump_type_local is not None:
            self.add_setting(config, ["namelist:recon", 
                                      "input_dump_type"], 
                             str(input_dump_type_local))
        return config, self.reports
        

class vn86_t5331(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5331 by Joe Mancell."""
    
    BEFORE_TAG = "vn8.6_t5702"
    AFTER_TAG = "vn8.6_t5331"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add your upgrade macro commands here
        #####################################################
        # REMOVE CAT A / DUPLICATED ITEMS                   #
        #####################################################
        self.remove_setting(config, ["namelist:stshcomp", "ocaaa"])
        self.remove_setting(config, ["namelist:stshcomp", "run_target_end"])
        self.remove_setting(config, ["namelist:stshcomp", "lons"])
        self.remove_setting(config, ["namelist:stshcomp", "lats"])
        self.remove_setting(config, ["namelist:stshcomp", "swbnd"])
        self.remove_setting(config, ["namelist:stshcomp", "lwbnd"])
        self.remove_setting(config, ["namelist:stshcomp", "orogr"])
        self.remove_setting(config, ["namelist:stshcomp", "swmcr"])
        self.remove_setting(config, ["namelist:stshcomp", "meso"])
        self.remove_setting(config, ["namelist:stshcomp", "ocalb"])
        self.remove_setting(config, ["namelist:stshcomp", "lfloor"])
        self.remove_setting(config, ["namelist:stshcomp", "aobgrp"])
        self.remove_setting(config, ["namelist:stshcomp", "aobinc"])
        #####################################################
        # LAM_CONFIG  - NEW NAMELIST                        #
        #####################################################
        # Get lam config settings from stshcomp
        ewspacea = self.get_setting_value(config, ["namelist:stshcomp","ewspacea"])
        nsspacea = self.get_setting_value(config, ["namelist:stshcomp","nsspacea"])
        frstlata = self.get_setting_value(config, ["namelist:stshcomp","frstlata"])
        frstlona = self.get_setting_value(config, ["namelist:stshcomp","frstlona"])
        polelata = self.get_setting_value(config, ["namelist:stshcomp","polelata"])
        polelona = self.get_setting_value(config, ["namelist:stshcomp","polelona"])
        # Remove lam config variables from stshcomp
        self.remove_setting(config, ["namelist:stshcomp","ewspacea"])
        self.remove_setting(config, ["namelist:stshcomp","nsspacea"])
        self.remove_setting(config, ["namelist:stshcomp","frstlata"])
        self.remove_setting(config, ["namelist:stshcomp","frstlona"])
        self.remove_setting(config, ["namelist:stshcomp","polelata"])
        self.remove_setting(config, ["namelist:stshcomp","polelona"])
        # Add lam config variables to new namelist
        self.add_setting(config, ["namelist:lam_config"])
        if ewspacea:
            self.add_setting(config, ["namelist:lam_config","delta_lon"], ewspacea)
        if nsspacea:
            self.add_setting(config, ["namelist:lam_config","delta_lat"], nsspacea)
        if frstlata:
            self.add_setting(config, ["namelist:lam_config","frstlata"], frstlata)
        if frstlona:
            self.add_setting(config, ["namelist:lam_config","frstlona"], frstlona)
        if polelata:
            self.add_setting(config, ["namelist:lam_config","polelata"], polelata)
        if polelona:
            self.add_setting(config, ["namelist:lam_config","polelona"], polelona)
        #####################################################
        # OZONE  - NEW NAMELIST                             #
        #####################################################
        zon_av_ozone = self.get_setting_value(config, ["namelist:stshcomp","zonavozone"])
        self.remove_setting(config, ["namelist:stshcomp","zonavozone"])
        self.add_setting(config, ["namelist:run_ozone"])
        if zon_av_ozone:
            self.add_setting(config, ["namelist:run_ozone","zon_av_ozone"], zon_av_ozone)
        #####################################################
        # FREE TRACERS - NEW NAMELIST                       #
        #####################################################
        i_free_tracer = self.get_setting_value(config, ["namelist:stshcomp","tca"])
        i_free_tracer_lbc = self.get_setting_value(config, ["namelist:stshcomp","tca_lbc"])
        l_bl_tracer_mix =self.get_setting_value(config, ["namelist:nlstcatm","l_bl_tracer_mix"])
        self.remove_setting(config, ["namelist:stshcomp","tca"])
        self.remove_setting(config, ["namelist:stshcomp","tca_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_bl_tracer_mix"])
        if "1" in i_free_tracer:
            l_free_tracer = ".true."
        else:
            l_free_tracer = ".false."
        if "1" in i_free_tracer_lbc:
            l_free_tracer_lbc = ".true."
        else:
            l_free_tracer_lbc = ".false."
        self.add_setting(config, ["namelist:run_free_tracers"])
        if i_free_tracer:
            self.add_setting(config, ["namelist:run_free_tracers","i_free_tracer"], i_free_tracer)
        if i_free_tracer_lbc:
            self.add_setting(config, ["namelist:run_free_tracers","i_free_tracer_lbc"], i_free_tracer_lbc)
        if l_free_tracer:
            self.add_setting(config, ["namelist:run_free_tracers","l_free_tracer"], l_free_tracer)
        if l_free_tracer_lbc:
            self.add_setting(config, ["namelist:run_free_tracers","l_free_tracer_lbc"], l_free_tracer_lbc)
        if l_bl_tracer_mix:
            self.add_setting(config, ["namelist:run_free_tracers","l_bl_tracer_mix"], l_bl_tracer_mix)
        #####################################################
        # UPDATE SOURCE FOR NAMELIST TEXT FILES             #
        #####################################################
        # Add the new namelists to just after run_cloud in the SHARED file
        shared_source = self.get_setting_value(config, ["file:SHARED","source"])
        if shared_source:
            shared_source = re.sub(r'namelist:run_aerosol',
                                     r'namelist:run_aerosol namelist:lam_config namelist:run_ozone namelist:run_free_tracers',
                                     shared_source)
            self.change_setting_value(config, ["file:SHARED","source"], shared_source)
        # Add the new namelists to just after run_cloud in the NAMELIST file
        namelist_source = self.get_setting_value(config, ["file:NAMELIST","source"])
        if namelist_source:
            namelist_source = re.sub(r'namelist:run_aerosol',
                                     r'namelist:run_aerosol namelist:lam_config namelist:run_ozone namelist:run_free_tracers',
                                     namelist_source)
            self.change_setting_value(config, ["file:NAMELIST","source"], namelist_source)

        return config, self.reports


class vn86_t5367(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5367 by Harry Shepherd."""
    
    BEFORE_TAG = "vn8.6_t5331"
    AFTER_TAG = "vn8.6_t5367"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        #aliases for clarity in code
        gval = self.get_setting_value
        cntlgen_vars = {}
        original_vars = {'expt_id':'',
                         'job_id':'',
                         'run_resubmit_inc':'',
                         'model_status':'',
                         'model_basis_time':'',
                         'ancil_reftime':'',
                         'run_target_end':'',
                         'lclimrealyr':'',
                         'ltimer':'',
                         'time_convention':'',
                         'model_analysis_mins':'',
                         'model_assim_mode':'',
                         'run_assim_mode':'',
                         'control_resubmit':'',
                         'pp_pack_code':'',
                         'pp_len2_look':'',
                         'ft_steps':'',
                         'ft_select':'',
                         'ft_archsel':'',
                         'ft_firststep':'',
                         'ft_laststep':'',
                         'type_letter_2':'',
                         'num_albcs':''}
        
        #
        # We require two variables from cntlgen
        #
        cntlgen_vars['secs_per_period'] = self.get_setting_value(config,
                                               ['namelist:nlstcgen',
                                                'secs_per_periodim'])
        cntlgen_vars['steps_per_period'] = self.get_setting_value(config,
                                                ['namelist:nlstcgen',
                                                 'steps_per_periodim'])
        #
        # Read everything from cntlall
        #
        for original_var in original_vars:
            original_vars[original_var] = self.get_setting_value(config,
                                               ['namelist:nlstcall',
                                                original_var])
            # if the value is an array, we need to split this
            if ',' in original_vars[original_var]:
                original_vars[original_var] = \
                    original_vars[original_var].split(',')

        #
        # Caclulate the variables for the new namelists
        #
        output_vars = self.create_new_namelists(original_vars,cntlgen_vars)
        

        #
        # Remove defunct namelist variables
        #
        variables_to_remove = ['ft_archsel',
                               'ft_firststep',
                               'ft_laststep',
                               'ft_select',
                               'ft_steps',
                               'pp_len2_look',
                               'pp_pack_code',
                               'type_letter_2']
        for v_t_r in variables_to_remove:
            self.remove_setting(config, ["namelist:nlstcall", \
                                             v_t_r])

        #
        # Existing namelist variables that require changing
        #
        #
        # Swap run_assim_mode and model_assim_mod that have the values
        # 'None      ' (6 spaces) to None'
        r_a_m = self.get_setting_value(config,
                                       ['namelist:nlstcall',
                                        'run_assim_mode'])
        if 'None'.lower() in r_a_m.lower():
            self.change_setting_value(config,
                                      ['namelist:nlstcall',
                                       'run_assim_mode'],
                                      '\'None\'')
        
        m_a_m = self.get_setting_value(config,
                                       ['namelist:nlstcall',
                                        'model_assim_mode'])
        if 'None'.lower() in m_a_m.lower():
            self.change_setting_value(config,
                                      ['namelist:nlstcall',
                                       'model_assim_mode'],
                                      '\'None\'')

        #
        # Write our new file:namelist config
        #
        header = self.get_setting_value(config,
                                        ['file:NAMELIST', 'source'])
        try:
            header = header + ' ' + \
                'namelist:nlstcall_qs ' + \
                'namelist:nlstcall_makebc ' + \
                'namelist:nlstcall_var ' + \
                'namelist:nlstcall_pp(:)'
            self.change_setting_value(config,
                                      ['file:NAMELIST', 'source'],
                                      header)
        except TypeError:
            #dealing with the single column model
            pass
        #
        # Write our new namelists
        #
        # QS NAMELIST
        #
        qs_vars = ['arch_sys',
                   'autopp',
                   'gcmdel',
                   'gddel',
                   'gpdel',
                   'ppxg',
                   'ppxm']
        for qs_var in qs_vars:
            self.add_setting(config,
                             ['namelist:nlstcall_qs',
                              qs_var],
                             value = str(output_vars[qs_var]))
        #
        # MAKEBC NAMELIST
        #
        mkbc_vars = ['imkbc','mbc_end','mbc_freq','mbc_strt','mbc_unt']
        for mkbc_var in mkbc_vars:
            self.add_setting(config,
                             ['namelist:nlstcall_makebc',
                              mkbc_var],
                             value = str(output_vars[mkbc_var]))
        #
        # VAR NAMELIST
        #
        var_vars = ['varend',
                    'varintfc',
                    'varintm',
                    'varintvl',
                    'varstart',
                    'varunits',
                    'packvar',
                    'packvarbc']
        for var_var in var_vars:
            self.add_setting(config,
                             ['namelist:nlstcall_var',
                              var_var],
                             value = str(output_vars[var_var]))
        #
        # PP NAMELIST
        #
        pp_vars = ['ppa',
                   'ppg',
                   'ppi',
                   'ppie',
                   'ppif',
                   'ppis',
                   'ppiu',
                   'ppos',
                   'ppx']
        for i in xrange(0,11):
            if i == 10:
                j = '151'
            else:
                j = str(i+60)
            for pp_var in pp_vars:
                self.add_setting(config,
                                 ['namelist:nlstcall_pp('+j+')',
                                  'i_'+pp_var],
                                 value = str(output_vars[pp_var][i]))
            
        return config, self.reports

    def fromtime(self,outtime,period_h,spp):        
        '''
        Inverse of the UMUI function totime
        '''
        if spp != 0:
            outtime = int(outtime)
            period_h = int(period_h)
            return str(outtime*period_h / spp)
        else:
            #This is probably the single column model
            return '-32768'

    def create_new_namelists(self,in_vars,cg_vars):
        '''
        Inverse of the logic in src/control/top_level/pop_nlstcall_mod.F90
        to take the existing UM inputs from nlstcall and divide them out into
        the new namelists
        nlstcall -> unchanged variables
        nlstcall_qs -> questions that pertain to all the new namelists
        nlstcall_lbc(1) - nlstcall_lbc(8) -> lateral boundary conditions
        nlstcall_makebc -> makebc
        nlstcall_var -> var macro
        nlstcall_pp(60) - nlstcall(69) and nlstcall(151) -> pp files
        '''

        cg_vars['steps_per_period'] = int(cg_vars['steps_per_period'])
        out_vars = {'autopp':'\'N\'',
                'gcmdel':'\'N\'',
                'gddel':'\'N\'',
                'ppos':11*['0'],
                'ppx':11*['0'],
                'ppg':11*['\'N\''],
                'ppxm':'0',
                'ppxg':'\'N\'',
                'packvar': '\'N\'',
                'imkbc':'0',
                'packvarbc':'\'N\'',
                'gpdel':'\'N\'',
                'ppa':11*['\'N\''],
                'arch_sys':'0',
                'ppi':11*['\'N\''],
                'ppiu':11*['\'\''],
                'ppif':11*['0'],
                'ppis':11*['0'],
                'ppie':11*['0'],
                'ilmari':8*['\'N\''],
                'ocbila':8*['0'],
                'ilmarh':8*['0'],
                'ilmas':8*['0'],
                'ilmara':8*['\'N\''],
                'imkbc':'0',
                'mbc_unt':'\'\'',
                'mbc_freq':'0',
                'mbc_strt':'0',
                'mbc_end':'0',
                'varunits':'\'T\'',
                'varintvl':'0',
                'varstart':'0',
                'varend':'0',
                'varintfc':'3',
                'varintm':'0',
                }
    
    
        #
        # ft select
        #
        if in_vars['ft_select'][2].lower() == '\'y\'' or \
                in_vars['ft_select'][22].lower() == '\'y\'':
            out_vars['autopp'] = '\'Y\''
            out_vars['gddel'] = '\'Y\''
        if in_vars['ft_select'][7].lower() == 'y':
            out_vars['autopp'] = '\'Y\''
            out_vars['gcmdel'] = '\'Y\''

        #
        # pp_len_2
        #
        for i_pp in xrange(1,12):
            if i_pp == 11:
                i_u = 151
            else:
                i_u = i_pp + 59
            if (i_u >= 60 and i_u <= 69):
                out_vars['ppos'][i_pp-1] = in_vars['pp_len2_look'][i_u-20]
            elif i_u == 151:
                out_vars['ppos'][i_pp-1] = in_vars['pp_len2_look'][i_u-20]

        #
        # pp_pack_code
        #
        for i_pp in xrange(1,12):
            if i_pp == 11:
                i_u = 151
            else:
                i_u = i_pp + 59
            if int(in_vars['pp_pack_code'][i_u-20]) >= 100:
                out_vars['ppx'][i_pp-1] = \
                    str(int(in_vars['pp_pack_code'][i_u-20]) - 100)
                out_vars['ppg'][i_pp-1] = '\'Y\''
            else:
                out_vars['ppx'][i_pp-1] = in_vars['pp_pack_code'][i_u-20]
        #
        if int(in_vars['pp_pack_code'][7]) >= 100:
            out_vars['ppxg'] = '\'Y\''
        else:
            out_vars['ppxm'] = in_vars['pp_pack_code'][7]
        #
        if in_vars['pp_pack_code'][130] == '1':
            packvar = '\'Y\''
        if in_vars['pp_pack_code'][144] == '1':
            imkbc = '1'
            packvarbc = '\'Y\''

        #
        # PP files
        #
        for i_pp in xrange(1,12):
            if i_pp == 11:
                i_u = 151
            else:
                i_u = i_pp + 59
            if int(in_vars['ft_steps'][i_u-20]) != 0:
                out_vars['autopp'] = '\'Y\''
                out_vars['gpdel'] = '\'Y\''
                out_vars['ppi'][i_pp-1] = '\'Y\''
            if in_vars['ft_archsel'][i_u-20] == '\'Y\'':
                out_vars['ppa'][i_pp-1] = '\'Y\''
                out_vars['arch_sys'] = '2'
            if int(in_vars['ft_steps'][i_u-20]) < 0:
                out_vars['ppiu'][i_pp-1] = 'RM'
                out_vars['ppif'][i_pp-1] = str(-1*int(in_vars['ft_steps'] \
                                                          [i_u-20]))
                out_vars['ppis'][i_pp-1] = in_vars['ft_firststep'][i_u-20]
                out_vars['ppie'][i_pp-1] = in_vars['ft_laststep'][i_u-20]
            else:
                hours_per_period = int(float(cg_vars['secs_per_period']) / \
                                           3600.)
                out_vars['ppiu'][i_pp-1] = '\'H\''
                out_vars['ppif'][i_pp-1] = self.fromtime( \
                    in_vars['ft_steps'][i_u-20], hours_per_period, \
                        cg_vars['steps_per_period'])
                out_vars['ppis'][i_pp-1] = self.fromtime( \
                    in_vars['ft_firststep'][i_u-20],hours_per_period, \
                        cg_vars['steps_per_period'])
                out_vars['ppie'][i_pp-1] = self.fromtime( \
                    in_vars['ft_laststep'][i_u-20],hours_per_period, \
                        cg_vars['steps_per_period'])
        

        #
        # MakeBC
        #
        if in_vars['type_letter_2'][144] == '\'a\'' :
            out_vars['imkbc'] = 1
            hours_per_period = int(float(cg_vars['secs_per_period']) / 3600.)
            out_vars['mbc_unt'] = '\'T\''
            out_vars['mbc_freq'] = in_vars['ft_steps'][144]
            out_vars['mbc_strt'] = str(int(in_vars['ft_firststep'][144])+1)
            out_vars['mbc_end'] = in_vars['ft_laststep'][144]

        #
        # 4D var
        #
        if in_vars['ft_select'][130] == '\'N\'' and \
                in_vars['type_letter_2'][130] == '\'a\'':
            out_vars['varintfc'] = '4'
            out_vars['varintm'] = '1'
            out_vars['varunits'] = '\'T\''
            out_vars['varintvl'] = in_vars['ft_steps'][130]
            out_vars['varstart'] = str(int(in_vars['ft_firststep'][130]) + 1)
            out_vars['varend'] = str(int(in_vars['ft_laststep'][130]) + 1)
        
        out_vars = dict(out_vars.items()+in_vars.items())
        return out_vars

class vn86_t5850(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5850 by Glenn Greed.  
       Retire the need for the nsubmodl namelist"""
    
    BEFORE_TAG = "vn8.6_t5367"
    AFTER_TAG = "vn8.6_t5850"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add your upgrade macro commands here
        namelsts = self.get_setting_value(config, ["file:SIZES","source"],)

        if (namelsts):
            namelsts=namelsts.replace("namelist:nsubmodl", "")
            self.change_setting_value(config, ["file:SIZES","source"], 
                                                             namelsts)
        namelsts = self.get_setting_value(config, ["file:NAMELIST","source"],)
        if (namelsts): 
            namelsts=namelsts.replace("namelist:nsubmodl", "")
            self.change_setting_value(config, ["file:NAMELIST","source"], 
                                                             namelsts)
        self.remove_setting(config, ["namelist:nsubmodl"])
        return config, self.reports

class vn86_t5887(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5887 by Thomas Allen."""
    
    BEFORE_TAG = "vn8.6_t5850"
    AFTER_TAG = "vn8.6_t5887"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.add_setting(config, ["namelist:run_dyn", "damp_height"], "80000.0")
        return config, self.reports

class vn86_t5817(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5817 by Glenn Greed."""
    
    BEFORE_TAG = "vn8.6_t5887"
    AFTER_TAG = "vn8.6_t5817"

    def UM_version(self, config) :
        ''' a simple routine to identify the UM version 
            from the rose-app.conf '''

        meta = self.get_setting_value(config, ["meta"])
        pversion = re.compile('um-atmos/(?P<VERSION>vn[0-9].[0-9])')
        version = pversion.search(meta)

        if version:
           return version.group('VERSION')

        return None


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

                    error = ('At the very least either ANCIL_VERSIONS or '+ 
                             'INITFILENV must be provided to upgrade.')
                    raise UpgradeError(error) 

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

        # find the items namelist objects.
        for obj in config.get_value():

            if re.search(r'namelist:items', obj):
               
                phash = re.compile('items\((?P<HASH>[a-g0-9]+)\)')
                namelstname = phash.search(obj) 

                source = self.get_setting_value(config, [obj, "source"])
               
                # If the above 'get' returns "None", we've reached the 
                # last item in the file, so end the loop
                if not source:
                    break

                # If source is not equal to two, we don't care about it
                # as an ancillary is not being sourced.
                if source != "2":
                    continue

                item = self.get_setting_value(config, [obj, "item"])
                section = self.get_setting_value(config, [obj, "section"])
                stashcode = str(1000*int(section) + int(item))

                # create array of source=2 namelist names to reuse later.
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


    def reformat_items_namelist(self, config):
        '''Where ITEMS are requested but not sourced from an ancillary 
        file then reformat their requests to use stashcodes rather  
        than section and item.'''

        for obj in config.get_value():
            if re.search(r'namelist:items', obj):

                source = self.get_setting_value(config, [obj, "source"])
                domain = self.get_setting_value(config, [obj, "domain"])
                section = self.get_setting_value(config, [obj, "section"])
                item = self.get_setting_value(config, [obj, "item"])
                usect = self.get_setting_value(config, 
                                                [obj, "user_prog_ancil_sctnc"])
                uitem = self.get_setting_value(config, 
                                                [obj, "user_prog_ancil_itemc"])
                filen=self.get_setting_value(config, 
                                                [obj, "user_prog_ancil_file"])

            # If the above 'get' returns "None", we've reached the last item
            # in the file, so end the loop
                if not source:
                    break

            # If the source is an ancil, we'll deal with that later.
                if source == "2":
                    continue

               # If user progs also need to reformat those requests
                elif source == "7":
                    # namelists may not have section set
                    if not usect:
                        usect="0"
                    if not section:
                        section="0"
                    stashcode = str(1000*int(section) + int(item))
                    ustshcode = str(1000*int(usect) + int(uitem))
                    self.remove_setting(config, [obj, "item"] )
                    self.remove_setting(config, [obj, "section"] )
                    self.add_setting(config, [obj, "stash_req"], stashcode)
                    self.remove_setting(config, [obj, 
                                        "user_prog_ancil_sctnc"] )
                    self.remove_setting(config, [obj, 
                                        "user_prog_ancil_itemc"] )
                    self.add_setting(config, [obj, 
                                     "user_prog_ancil_stash_req"], ustshcode)
                    self.remove_setting(config, [obj, 
                                        "user_prog_ancil_file"] )
                    self.add_setting(config, [obj, 
                                     "ancilfilename"], filen)

                else:
                    # update the namelist items to new form.
                    stashcode = str(1000*int(section) + int(item))
                    self.remove_setting(config, [obj, "section"] )
                    self.remove_setting(config, [obj, "item"] )
                    self.add_setting(config, [obj, "stash_req"], stashcode)

        return config


    def add_ancil_items(self, items_req, config,namelstnm):
        '''Now for each item sourced from ancillary we need to add
        the new updated format including the ancillary paths. '''

        # remove all source 2 namelists using old form.
        for name in namelstnm:
            namelist = "namelist:items(%s)"%name
            self.remove_setting(config, [namelist] )

        # simply loop over the ancils required to add new namelists
        # reusing the names that have been removed. Thus ensuring
        # uniquness of names.
 
        next_item=0
        
        for anc, stash in items_req.iteritems():
            namelist = "namelist:items(%s)"%(namelstnm[next_item])
            self.add_setting(config, [namelist, "source"], "2")
            self.add_setting(config, [namelist, "domain"], "1")
            stashreq = ",".join(stash)
            self.add_setting(config, [namelist, "stash_req"], stashreq)
            self.add_setting(config, 
                            [namelist, "ancilfilename"], "'%s'"%(anc))
            self.remove_setting(config, [namelist, "section"] )
            self.remove_setting(config, [namelist, "item"] )

            next_item = next_item + 1

        return config



    def upgrade(self, config, meta_config=None):
        """Upgrade, where possible, UM runtime app configurations to use new 
           RCF items namelist formats.
           This is only possible where std ancil version files are used that
           appear on both the remote and local machines"""

        # This upgrade macro is required to at the very least to interogate
        # the app file to be updated and for most configurations
        # it also will need to look at the INITFILEENV file.
        # It will then generate a link (dictionary) between stash codes
        # and ancillary files with the help of the ANCILmaster.

        initfilenv='./file/INITFILEENV'
        roseapp='./rose-app.conf'

        # check the important rose app is available.
        if os.path.isfile(roseapp):
            # file exists so continue
            pass
        else :
            # no file so exit with message.
            raise UpgradeException ('Cannot continue unable to find the rose app config file '
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

        # This script will use vn8.6 ANCILmaster unless the rose app provided
        # is clearly using an older UM version

        vn = self.UM_version(config)

        if not vn:
            self.add_report("env", "VN", None, 
               'UM VN is not set - assuming vn8.6', is_warning=True)
            vn='vn8.6'


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

        # Interogate the provided config to look for any item 
        # requests that require sourcing from ancillary files
        # set up a dictionary of requested [ ancillary full paths: stashcodes]

        items_req,namelstnm = self.generate_stash_to_anc_dict(config, 
                              envvars, stash2ancilenvname,hpc_app)


        # Now reformat the app files, using the new ITEMS format.
        # firstly all items not soruced from an ancillary
        config = self.reformat_items_namelist(config)
        
        # secondly add the new namelists for ancillary sourcing.
        config = self.add_ancil_items(items_req, config, namelstnm)     
        
        return config, self.reports


class vn86_t5722(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5722 by Harry Shepherd"""
    
    BEFORE_TAG = "vn8.6_t5817"
    AFTER_TAG = "vn8.6_t5722"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # Get list of domain profiles (dom_name) used by diagnostic requests
        # with isec=2 and item in the list of items

        domain_profiles = []
        req_isec = 2
        req_item = [ 251,
                    252,
                    253,
                    254,
                    255,
                    256,
                    284,
                    285,
                    286,
                    287,
                    288,
                    289,
                    295,
                    296,
                    297,
                    298,
                    299,
                    300,
                    301,
                    302,
                    303,
                    304,
                    305,
                    421,
                    422,
                    423,
                    424,
                    425,
                    426,
                    427]

        for obj in config.get_value():
            if re.search(r'namelist:streq',obj):
                streq_isec = self.get_setting_value(config,
                                                    [obj,'isec'])
                # will return none when run out of streq namelists
                if not streq_isec:
                    break
                elif streq_isec == str(req_isec):
                    streq_item = self.get_setting_value(config,
                                                        [obj,'item'])
                    if streq_item in map(str,req_item):
                        domain_profiles.append(self.get_setting_value(config,
                                                        [obj,'dom_name']))
                else:
                    pass

        #remove duplicates from domain_profiles as this is required
        #for checking if we can upgrade
        #note a set by definitiion can not contain duplicates and we dont
        #care about order
        domain_profiles = list(set(domain_profiles))


        # check to see if the domain profiles are used by any other stash
        # request
        for obj in config.get_value():
            if re.search(r'namelist:streq',obj):
                i_dom = self.get_setting_value(config,
                                               [obj,'dom_name'])
                if i_dom in domain_profiles:
                    i_sec = self.get_setting_value(config,
                                                   [obj,'isec'])
                    i_item = self.get_setting_value(config,
                                                    [obj,'item'])
                    if (i_sec == str(req_isec) and \
                            i_item in map(str,req_item)):
                        # we can continue
                        pass
                    else:
                        #generate an error string that can be passed to
                        #std out
                        raise UpgradeError(
                            'Domain profile %s used by section %s' \
                            ' item %s, and so can not be automatically' \
                            ' upgraded' % (i_dom, i_sec, i_item))

        #perform the upgrade
        if len(domain_profiles) == 0:
            pass
        else:
            for obj in config.get_value():
                if re.search(r'namelist:domain',obj):
                    name = self.get_setting_value(config,
                                                  [obj,'dom_name'])
                    if (name in domain_profiles):
                        self.change_setting_value(config,
                                                  [obj, 'plt'],
                                                  '4')

        return config, self.reports 

class vn86_t5841(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5841 by Malcolm Brooks."""
    
    BEFORE_TAG = "vn8.6_t5722"
    AFTER_TAG = "vn8.6_t5841"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add your upgrade macro commands here
        self.add_setting(config, ["namelist:temp_fixes", 
                                  "l_fix_arcl_eg_levs"], ".false.")
        return config, self.reports


class vn86_t5448(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5448 by <SJ Swarbrick>."""
    
    BEFORE_TAG = "vn8.6_t5841"
    AFTER_TAG = "vn8.6_t5448"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add your upgrade macro commands here

        self.add_setting(config, ["namelist:run_diffusion", "hdiffopt"], "0")
        l_diffusion = self.get_setting_value(
                       config, ["namelist:run_diffusion", "l_diffusion"])
        l_subfilter_horiz = self.get_setting_value(
                       config, ["namelist:run_diffusion", "l_subfilter_horiz"])
        l_subfilter_vert = self.get_setting_value(
                       config, ["namelist:run_diffusion", "l_subfilter_vert"])
        l_subfilter_blend = self.get_setting_value(
                       config, ["namelist:run_diffusion", "l_subfilter_blend"])
        if l_diffusion == '.true.':
            self.change_setting_value(
                       config, ["namelist:run_diffusion", "hdiffopt"],"1")
        if l_subfilter_blend == '.true.' or l_subfilter_horiz == '.true.':
            self.change_setting_value(
                       config, ["namelist:run_diffusion", "hdiffopt"],"3")
        if l_subfilter_vert == '.true.':
            self.change_setting_value(
                       config, ["namelist:run_diffusion", "hdiffopt"],"3")

                     
        l_cdiffusion = self.get_setting_value(
                       config, ["namelist:run_diffusion", "l_cdiffusion"])
        if l_cdiffusion == '.true.':
            self.change_setting_value(
                       config, ["namelist:run_diffusion", "hdiffopt"],"2")
        
        
        l_vertical_diffusion = self.get_setting_value(
                    config, ["namelist:run_diffusion", "l_vertical_diffusion"])
        l_ramp = self.get_setting_value(
                    config, ["namelist:run_diffusion", "l_ramp"])
    
        
        self.add_setting(config, ["namelist:run_diffusion", "vdiffopt"], "0")
        if l_vertical_diffusion == '.true.':
          if l_ramp == '.false.':
            self.change_setting_value(
                          config, ["namelist:run_diffusion", "vdiffopt"],"1")
          elif l_ramp == '.true.':
            self.change_setting_value(
                          config, ["namelist:run_diffusion", "vdiffopt"],"2")
        
        
        self.add_setting(config, ["namelist:run_diffusion", "pofil_opt"], "0")
        model_domain = self.get_setting_value(
                                   config, ["namelist:nlstcatm", "model_domain"])
        l_endgame = self.get_setting_value(
                                   config, ["namelist:nlstcatm","l_endgame"])
        if model_domain == '1' and l_endgame == '.false.':
            self.change_setting_value(
                    config, ["namelist:run_diffusion", "pofil_opt"], "1")

        
        self.remove_setting(config, ["namelist:run_diffusion", "l_diffusion"])
        self.remove_setting(config, ["namelist:run_diffusion", "l_cdiffusion"])
        self.remove_setting(
                    config, ["namelist:run_diffusion", "l_vertical_diffusion"])
        self.remove_setting(config, ["namelist:run_diffusion", "l_ramp"])
        self.remove_setting(config, ["namelist:run_diffusion", "l_divdamp"])
        self.remove_setting(
                    config, ["namelist:run_diffusion", "div_damp_coefficient"])
        self.remove_setting(config, ["namelist:run_diffusion", "l_diff_ctl"])
        self.remove_setting(config, ["namelist:run_diffusion", "l_pfexner"]) 
        self.remove_setting(config, ["namelist:run_diffusion", "blockx_in"])
        self.remove_setting(config, ["namelist:run_diffusion", "blocky_in"])
        self.remove_setting(config, ["namelist:run_diffusion", 
                                                      "diffusion_coefficient_q"])
        self.remove_setting(config, ["namelist:run_diffusion", 
                                                 "diffusion_coefficient_thermo"])
        self.remove_setting(config, ["namelist:run_diffusion", 
                                                   "diffusion_coefficient_wind"])
        self.remove_setting(config, ["namelist:run_diffusion","diffusion_order_q"])
        self.remove_setting(config, ["namelist:run_diffusion",
                                                         "diffusion_order_thermo"])
        self.remove_setting(config, ["namelist:run_diffusion",
                                                           "diffusion_order_wind"])
        
        return config, self.reports


class vn86_t5193(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5193 by Richard Barnes."""
    
    BEFORE_TAG = "vn8.6_t5448"
    AFTER_TAG = "vn8.6_t5193"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add your upgrade macro commands here
        ##################################################### 
        # Remove UKCA tracer input from namelist            #
        ##################################################### 
        """Remove tc_ukca list from stshcomp namelist."""
        self.remove_setting(config, ["namelist:stshcomp", "tc_ukca"])
        ##################################################### 
        # UKCA LBC configuration - move namelist            #
        ##################################################### 
        tc_lbc_ukca = self.get_setting_value(config, ["namelist:stshcomp","tc_lbc_ukca"])
        self.remove_setting(config, ["namelist:stshcomp","tc_lbc_ukca"]) 
        self.add_setting(config, ["namelist:run_ukca","tc_lbc_ukca"], tc_lbc_ukca)

        # Remove stshcomp completely - no longer in UM
        self.remove_setting(config, ["namelist:stshcomp"]) 

        # update the source of the UM NAMELIST file to remove stshcomp
        namelist_source = self.get_setting_value(config, ["file:NAMELIST","source"])
        if namelist_source:
            namelist_source = re.sub(r'namelist:stshcomp','',namelist_source)
            self.change_setting_value(config, ["file:NAMELIST","source"], namelist_source)
        # update the source of the UM SIZES file to remove stshcomp
        sizes_source = self.get_setting_value(config, ["file:SIZES","source"])
        if sizes_source:
            sizes_source = re.sub(r'namelist:stshcomp','',sizes_source)
            self.change_setting_value(config, ["file:SIZES","source"], sizes_source)

        return config, self.reports


class vn86_t5893(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5893 by Joe Mancell."""
    
    BEFORE_TAG = "vn8.6_t5193"
    AFTER_TAG = "vn8.6_t5893"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add your upgrade macro commands here
        n_rims_to_do = self.get_setting_value(config, ["namelist:run_dyn",
                                                       "n_rims_to_do"])
        self.remove_setting(config, ["namelist:run_dyn", "n_rims_to_do"])
        if n_rims_to_do: 
            self.add_setting(config, ["namelist:lam_config", "n_rims_to_do"], 
                             n_rims_to_do)
                                                         
        return config, self.reports


class vn86_t5796(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5796 by Ian Boutle."""
    
    BEFORE_TAG = "vn8.6_t5893"
    AFTER_TAG = "vn8.6_t5796"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add your upgrade macro commands here
        bllevs=self.get_setting_value(config, ["namelist:nlsizes", "bl_levels"])
# changes to alpha_cd, removing 3 options and adding a list
        if bllevs == "13":
            self.add_setting(config, ["namelist:run_bl", "alpha_cd"], "2*2.0,1.5,10*1.0")
        else:
            self.add_setting(config, ["namelist:run_bl", "alpha_cd"], "2.0,%s*1.5" %(str(int(bllevs)-1)))
        self.remove_setting(config, ["namelist:run_bl", "alpha_cd_batches"])
        self.remove_setting(config, ["namelist:run_bl", "alpha_cd_items"])
        self.remove_setting(config, ["namelist:run_bl", "alpha_cd_vals"])
# adding appropriate defaults for TKE scheme
        self.change_setting_value(config, ["namelist:run_bl", "bdy_tke"], "3")
        self.change_setting_value(config, ["namelist:run_bl", "tke_levels"], "-1")
        self.change_setting_value(config, ["namelist:run_bl", "l_my_initialize"], ".true.")
        self.change_setting_value(config, ["namelist:run_bl", "my_ini_dbdz_min"], "1.0e-5")
        self.add_setting(config, ["namelist:run_bl", "l_adv_turb_field"], ".true.")
        self.change_setting_value(config, ["namelist:run_bl", "l_my_condense"], ".true.")
        self.change_setting_value(config, ["namelist:run_bl", "l_shcu_buoy"], ".true.")
        self.change_setting_value(config, ["namelist:run_bl", "shcu_levels"], "-1")
        self.change_setting_value(config, ["namelist:run_bl", "wb_ng_max"], "0.05")
        self.change_setting_value(config, ["namelist:run_bl", "my_lowest_pd_surf"], "2")
        self.change_setting_value(config, ["namelist:run_bl", "my_prod_adj_fact"], "%s*0.225" %(bllevs))
        self.change_setting_value(config, ["namelist:run_bl", "my_z_limit_elb"], "1.0e10")
        self.change_setting_value(config, ["namelist:run_bl", "tke_cm_mx"], "0.1")
        self.change_setting_value(config, ["namelist:run_bl", "tke_cm_fa"], "0.1")
        self.change_setting_value(config, ["namelist:run_bl", "tke_dlen"], "1")
# moving jules parameters from run_bl to jules_surf_param
        charnock=self.get_setting_value(config, ["namelist:run_bl", "charnock"])
        self.add_setting(config, ["namelist:jules_surf_param", "charnock"], charnock)
        self.remove_setting(config, ["namelist:run_bl", "charnock"])
        seasalinityfactor=self.get_setting_value(config, ["namelist:run_bl", "seasalinityfactor"])
        self.add_setting(config, ["namelist:jules_surf_param", "seasalinityfactor"], seasalinityfactor)
        self.remove_setting(config, ["namelist:run_bl", "seasalinityfactor"])
# deleting duplicated items in run_bl & jules_switches
        self.remove_setting(config, ["namelist:run_bl", "buddy_sea"])
        self.remove_setting(config, ["namelist:run_bl", "cor_mo_iter"])
        self.remove_setting(config, ["namelist:run_bl", "iseaz0t"])
# moving jules options from run_bl to jules_switches
        formdrag=self.get_setting_value(config, ["namelist:run_bl", "formdrag"])
        self.add_setting(config, ["namelist:jules_switches", "formdrag"], formdrag)
        self.remove_setting(config, ["namelist:run_bl", "formdrag"])
        orog_drag_param=self.get_setting_value(config, ["namelist:run_bl", "orog_drag_param"])
        self.add_setting(config, ["namelist:jules_switches", "orog_drag_param"], orog_drag_param)
        self.remove_setting(config, ["namelist:run_bl", "orog_drag_param"])
        fd_stab_dep=self.get_setting_value(config, ["namelist:run_bl", "fd_stab_dep"])
        self.add_setting(config, ["namelist:jules_switches", "fd_stab_dep"], fd_stab_dep)
        self.remove_setting(config, ["namelist:run_bl", "fd_stab_dep"])
        isrfexcnvgust=self.get_setting_value(config, ["namelist:run_bl", "isrfexcnvgust"])
        self.add_setting(config, ["namelist:jules_switches", "isrfexcnvgust"], isrfexcnvgust)
        self.remove_setting(config, ["namelist:run_bl", "isrfexcnvgust"])
        return config, self.reports

class vn86_t5857(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5857 by Joe Mancell."""
    
    BEFORE_TAG = "vn8.6_t5796"
    AFTER_TAG = "vn8.6_t5857"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        problem_number = self.get_setting_value(config, 
                                                ["namelist:nlstcatm", 
                                                 "problem_number"])
        npmsl_height = self.get_setting_value(config, 
                                              ["namelist:nlstcatm", 
                                               "npmsl_height"])
        l_pmsl_sor = self.get_setting_value(config, 
                                            ["namelist:nlstcatm", 
                                             "l_pmsl_sor"])
        l_use_methox = self.get_setting_value(config, 
                                              ["namelist:nlstcatm", 
                                               "l_use_methox"])
        l_mr_physics1 = self.get_setting_value(config, 
                                               ["namelist:nlstcatm", 
                                                "l_mr_physics1"])
        l_mr_physics2 = self.get_setting_value(config, 
                                               ["namelist:nlstcatm", 
                                                "l_mr_physics2"])
        l_int_uvw_lbc = self.get_setting_value(config, 
                                               ["namelist:nlstcatm", 
                                                "l_int_uvw_lbc"])
        l_lateral_boundary = self.get_setting_value(config, 
                                                    ["namelist:nlstcatm", 
                                                     "l_lateral_boundary"])
        rimweightsa = self.get_setting_value(config, 
                                             ["namelist:bouncnst", 
                                              "rimweightsa"])
        l_hydrology = self.get_setting_value(config,
                                             ["namelist:nlstcatm",
                                              "l_hydrology"])
        self.remove_setting(config, ["namelist:bouncnst"])
        self.add_setting(config, ["namelist:run_calc_pmsl"])
        self.add_setting(config, ["namelist:gen_phys_inputs"])
        self.add_setting(config, ["namelist:lbc_options"])
        if l_int_uvw_lbc:
            self.add_setting(config, ["namelist:lbc_options", 
                                      "l_int_uvw_lbc"], l_int_uvw_lbc)
            self.remove_setting(config, ["namelist:nlstcatm", "l_int_uvw_lbc"])
        if l_lateral_boundary:
            self.add_setting(config, ["namelist:lbc_options", 
                                      "l_lateral_boundary"], l_lateral_boundary)
            self.remove_setting(config, ["namelist:nlstcatm", 
                                         "l_lateral_boundary"])
        if rimweightsa:
            self.add_setting(config, ["namelist:lbc_options", 
                                      "rimweightsa"], rimweightsa)
        if problem_number:
            self.add_setting(config, ["namelist:run_dyntest", 
                                      "problem_number"], problem_number)
            self.remove_setting(config, ["namelist:nlstcatm", 
                                         "problem_number"])
        if npmsl_height:
            self.add_setting(config, ["namelist:run_calc_pmsl", 
                                      "npmsl_height"], npmsl_height)
            self.remove_setting(config, ["namelist:nlstcatm", 
                                         "npmsl_height"])
        if l_pmsl_sor:
            self.add_setting(config, ["namelist:run_calc_pmsl", 
                                      "l_pmsl_sor"], l_pmsl_sor)
            self.remove_setting(config, ["namelist:nlstcatm", 
                                         "l_pmsl_sor"])
        if l_use_methox:
            self.add_setting(config, ["namelist:gen_phys_inputs", 
                                      "l_use_methox"], l_use_methox)
            self.remove_setting(config, ["namelist:nlstcatm", "l_use_methox"]) 
        if l_mr_physics1:
            self.add_setting(config, ["namelist:gen_phys_inputs", 
                                      "l_mr_physics1"], l_mr_physics1)
            self.remove_setting(config, ["namelist:nlstcatm", "l_mr_physics1"]) 
        if l_mr_physics2:
            self.add_setting(config, ["namelist:gen_phys_inputs", 
                                      "l_mr_physics2"], l_mr_physics2)
            self.remove_setting(config, ["namelist:nlstcatm", "l_mr_physics2"]) 
        if l_hydrology:
            self.add_setting(config, ["namelist:jules_switches",
                                      "l_hydrology"], l_hydrology)
            self.remove_setting(config, ["namelist:nlstcatm", "l_hydrology"])
        namelist_source = self.get_setting_value(config, ["file:NAMELIST","source"])
        cntlatm_source = self.get_setting_value(config, ["file:CNTLATM","source"])
        if namelist_source:
            namelist_source = namelist_source.replace('namelist:run_stochastic',
                                                      'namelist:run_stochastic namelist:run_calc_pmsl namelist:gen_phys_inputs namelist:lbc_options')
            namelist_source = namelist_source.replace('namelist:bouncnst',
                                                      '')
            self.change_setting_value(config, ["file:NAMELIST","source"], namelist_source)
        if cntlatm_source:
            cntlatm_source = cntlatm_source.replace('namelist:run_eng_corr',
                                                      'namelist:run_eng_corr  namelist:gen_phys_inputs')
            self.change_setting_value(config, ["file:CNTLATM","source"], cntlatm_source)
        return config, self.reports

class vn86_t5778(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5778 by Erica Neininger."""
    
    BEFORE_TAG = "vn8.6_t5857"
    AFTER_TAG = "vn8.6_t5778"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Remove automatic resubmission controls
        self.remove_setting(config, ["namelist:nlstcall", "control_resubmit"])
        self.remove_setting(config, ["namelist:nlstcall", "run_resubmit_inc"])

        return config, self.reports


class vn86_t5868(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5868 by <Glenn Greed>.
       Simply retire l_oasis from UM inputs"""
    
    BEFORE_TAG = "vn8.6_t5778"
    AFTER_TAG = "vn8.6_t5868"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        self.remove_setting(config, ["namelist:nlstcatm", "l_oasis"])
        return config, self.reports


class vn86_t5832(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5832 by Erica Neininger."""

    BEFORE_TAG = "vn8.6_t5868"
    AFTER_TAG = "vn8.6_t5832"

    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Retire INITHIS #
        self.remove_setting(config,["namelist:nlcfiles"])
        self.remove_setting(config,["namelist:nlchistg"])
        self.remove_setting(config,["namelist:nlchisto"])
        self.remove_setting(config,["namelist:nlihistg"])
        self.remove_setting(config,["namelist:nlihisto"])
        self.remove_setting(config,["file:INITHIS"])
        return config, self.reports

class vn86_t3972(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #3972 by hshep Harry Shepherd."""
    
    BEFORE_TAG = "vn8.6_t5832"
    AFTER_TAG = "vn8.6_t3972"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""

        # check we have a time profile TDMPMN, and if not create one
        is_tdmpmn = False
        for obj in config.get_value():
            if re.search(r'namelist:time',obj):
                time_name = self.get_setting_value(config,
                                                   [obj,'tim_name'])
                # will return none when run out of streq namelists
                if not time_name:
                    break
                else:
                    if 'TDMPMN' in time_name:
                        is_tdmpmn = True
                        break
                    else:
                        pass

        if not is_tdmpmn:
            # generate an error string that can be passed to std out
            time_warn = 'Time profile TDMPMN does not exist, and so the job' \
                ' may not be automatically upgraded properly'
            self.add_report("namelist:time", "tim_name", "TDMPMN",
                   time_warn, is_warning=True)
        else:
            # change any river routing diagnostics with usage UPMEAN to use
            # the standard time profile
            for obj in config.get_value():
                if re.search(r'namelist:streq',obj):
                    i_sec = self.get_setting_value(config,
                                                   [obj,'isec'])
                    i_use = self.get_setting_value(config,
                                                   [obj,'use_name'])
                    i_time = self.get_setting_value(config,
                                                    [obj,'tim_name'])
                    i_item = self.get_setting_value(config,
                                                    [obj,'item'])
                    sys.stdout.write(i_sec + ' and ' + i_use + '\n')
                    if ('26' in i_sec) and ('UPMEAN' in i_use):
                        self.change_setting_value(config,
                                                  [obj, 'tim_name'],
                                                  '\'TDMPMN\'')
                    elif ('26' in i_sec) and ('T3HDAYM' in i_time):
                        self.change_setting_value(config,
                                                  [obj, 'tim_name'],
                                                  '\'TDAYM\'')
                    elif ('8' in i_sec) and ('245' in i_item) and \
                            ('T3HDMRV' in i_time):
                        self.change_setting_value(config,
                                                  [obj, 'tim_name'],
                                                  '\'TDMPMN\'')

        # Fix the coupling time profile to allow the coupled tests to run
        for obj in config.get_value():
            if re.search(r'namelist:time',obj):
                i_name = self.get_setting_value(config,
                                                [obj,'tim_name'])
                if 'T3HR' in i_name:
                    self.change_setting_value(config,
                                              [obj,'unt3'],
                                              '\'T \'')
                    self.change_setting_value(config,
                                              [obj,'ifre'],
                                              '9')
                    self.change_setting_value(config,
                                              [obj,'istr'],
                                              '8')

        return config, self.reports

class vn86_t5769(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5769 by Joe Mancell."""
    
    BEFORE_TAG = "vn8.6_t3972"
    AFTER_TAG = "vn8.6_t5769"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add your upgrade macro commands here
        # Grab vegetation scheme version and move it to new namelist
        # delete the old namelist
        i_veg_vn = self.get_setting_value(config, ["namelist:run_blveg","i_veg_vn"])
        self.remove_setting(config, ["namelist:run_blveg","i_veg_vn"])
        self.remove_setting(config, ["namelist:run_blveg"])
        if i_veg_vn:
            self.add_setting(config, ["namelist:jules_vegetation","i_veg_vn"],
                             i_veg_vn)
        # move items from nlstcatm to jules_vegetation
        l_nrun_mid_trif = self.get_setting_value(config, ["namelist:nlstcatm","l_nrun_mid_trif"])
        l_q10 = self.get_setting_value(config, ["namelist:nlstcatm","l_q10"])
        l_trif_eq = self.get_setting_value(config, ["namelist:nlstcatm","l_trif_eq"])
        l_phenol = self.get_setting_value(config, ["namelist:nlstcatm","l_phenol"])
        phenol_period = self.get_setting_value(config, ["namelist:nlstcatm","phenol_period"])
        triffid_period = self.get_setting_value(config, ["namelist:nlstcatm","triffid_period"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_nrun_mid_trif"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_trif_eq"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_phenol"])
        self.remove_setting(config, ["namelist:nlstcatm", "phenol_period"])
        self.remove_setting(config, ["namelist:nlstcatm", "triffid_period"])
        if l_nrun_mid_trif:
            self.add_setting(config, ["namelist:jules_vegetation","l_nrun_mid_trif"], l_nrun_mid_trif)
        if l_q10:
            self.add_setting(config, ["namelist:jules_vegetation","l_q10"], l_q10)
        if l_trif_eq:
            self.add_setting(config, ["namelist:jules_vegetation","l_trif_eq"], l_trif_eq)
        if l_phenol:
            self.add_setting(config, ["namelist:jules_vegetation","l_phenol"], l_phenol)
        if phenol_period:
            self.add_setting(config, ["namelist:jules_vegetation","phenol_period"], phenol_period)
        if triffid_period:
            self.add_setting(config, ["namelist:jules_vegetation","triffid_period"], triffid_period)
        # remove items from namelist
        self.remove_setting(config, ["namelist:nlstcatm", "l_disturb"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_veg_fracs"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_triffid"])
        self.remove_setting(config, ["namelist:jules_switches", "l_triffid"])
        self.remove_setting(config, ["namelist:jules_switches", "l_phenol"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_q10"])
        # update the source of the UM NAMELIST file to remove run_blveg
        # and to add jules_vegetation
        namelist_source = self.get_setting_value(config, ["file:NAMELIST","source"])
        if namelist_source:
            namelist_source = re.sub(r'namelist:run_blveg','',namelist_source)
            namelist_source = re.sub(r'namelist:urban2t_param',
                                     r'namelist:urban2t_param namelist:jules_vegetation',
                                     namelist_source)
            self.change_setting_value(config, ["file:NAMELIST","source"], namelist_source)
        # update the source of the SCM CNTLATM file to remove run_blveg
        # and to add jules_vegetation
        cntlatm_source = self.get_setting_value(config, ["file:CNTLATM","source"])
        if cntlatm_source:
            cntlatm_source = re.sub(r'namelist:run_blveg','',cntlatm_source)
            self.change_setting_value(config, ["file:CNTLATM","source"], cntlatm_source)
        # update the source line of the UM SHARED file to add jules_vegetation
        shared_source = self.get_setting_value(config, ["file:SHARED","source"])
        if shared_source:
            shared_source = re.sub(r'namelist:jules_snow_param',
                                   r'namelist:jules_snow_param namelist:jules_vegetation',
                                   shared_source)
            self.change_setting_value(config, ["file:SHARED","source"], shared_source)
        return config, self.reports


class vn86_t5843(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5843 by Richard Barnes."""
    
    BEFORE_TAG = "vn8.6_t5769"
    AFTER_TAG = "vn8.6_t5843"
    
    def upgrade(self, config, meta_config=None):
        """Upgrade a UM runtime app configuration."""
        # Add your upgrade macro commands here
        ##################################################### 
        # REMOVE CAT A / DUPLICATED ITEMS                   # 
        ##################################################### 
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_bmass_agd_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_bmass_cld_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_bmass_new_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_dms_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_dust_div1_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_dust_div2_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_dust_div3_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_dust_div4_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_dust_div5_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_dust_div6_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", "l_nh3_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_nitr_acc_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_nitr_diss_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_ocff_agd_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_ocff_cld_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_ocff_new_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_so2_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_so4_accu_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_so4_aitken_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_so4_diss_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_soot_agd_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_soot_cld_lbc_out"])
        self.remove_setting(config, ["namelist:nlstcatm", 
        "l_soot_new_lbc_out"])


        ##################################################### 
        # LBCs INTO DIFFERENT NAMELISTS                     # 
        ##################################################### 
        # Get _lbc and other settings from nlstcatm 
        l_bmass_agd_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_bmass_agd_lbc"])
        l_dms_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_dms_lbc"])
        l_dust_div1_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_dust_div1_lbc"])
        l_dust_div2_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_dust_div2_lbc"])
        l_dust_div3_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_dust_div3_lbc"])
        l_dust_div4_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_dust_div4_lbc"])
        l_dust_div5_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_dust_div5_lbc"])
        l_dust_div6_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_dust_div6_lbc"])
        l_nh3_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_nh3_lbc"])
        l_nitr_acc_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_nitr_acc_lbc"])
        l_ocff_agd_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_ocff_agd_lbc"])
        l_so2_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_so2_lbc"])
        l_soot_agd_lbc = self.get_setting_value(config, 
        ["namelist:nlstcatm","l_soot_agd_lbc"])
        l_consistent_cdnc = self.get_setting_value(config,
        ["namelist:nlstcatm","l_consistent_cdnc"])
        l_use_sulpc_direct = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_sulpc_direct"])
        l_use_sulpc_indirect_lw = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_sulpc_indirect_lw"])
        l_use_sulpc_indirect_sw = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_sulpc_indirect_sw"])
        l_use_seasalt_direct = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_seasalt_direct"])
        l_use_seasalt_indirect = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_seasalt_indirect"])
        l_use_biogenic = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_biogenic"])
        l_use_nitrate_direct = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_nitrate_direct"])
        l_use_nitrate_indirect = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_nitrate_indirect"])
        l_use_soot_direct = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_soot_direct"])
        l_use_soot_indirect = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_soot_indirect"])
        l_use_bmass_direct = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_bmass_direct"])
        l_use_bmass_indirect = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_bmass_indirect"])
        l_use_ocff_direct = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_ocff_direct"])
        l_use_ocff_indirect = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_ocff_indirect"])
        l_use_dust = self.get_setting_value(config,
        ["namelist:nlstcatm","l_use_dust"])


        # Remove _lbc and other variables from nlstcatm 
        self.remove_setting(config, ["namelist:nlstcatm","l_bmass_agd_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_bmass_cld_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_bmass_new_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_dms_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_dust_div1_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_dust_div2_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_dust_div3_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_dust_div4_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_dust_div5_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_dust_div6_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_nh3_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_nitr_acc_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_nitr_diss_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_ocff_agd_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_ocff_cld_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_ocff_new_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_so2_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_so4_accu_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_so4_aitken_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_so4_diss_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_soot_agd_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_soot_cld_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_soot_new_lbc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_consistent_cdnc"])
        self.remove_setting(config, ["namelist:nlstcatm","l_use_sulpc_direct"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_sulpc_indirect_lw"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_sulpc_indirect_sw"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_seasalt_direct"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_seasalt_indirect"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_biogenic"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_nitrate_direct"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_nitrate_indirect"])
        self.remove_setting(config, ["namelist:nlstcatm","l_use_soot_direct"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_soot_indirect"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_bmass_direct"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_bmass_indirect"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_ocff_direct"])
        self.remove_setting(config, ["namelist:nlstcatm",
        "l_use_ocff_indirect"])
        self.remove_setting(config, ["namelist:nlstcatm","l_use_dust"])

        # Adding _lbc and other variables to their new namelists 
        self.add_setting(config, ["namelist:run_dust","l_dust_div1_lbc"], 
        l_dust_div1_lbc)
        self.add_setting(config, ["namelist:run_dust","l_dust_div2_lbc"], 
        l_dust_div2_lbc)
        self.add_setting(config, ["namelist:run_dust","l_dust_div3_lbc"], 
        l_dust_div3_lbc)
        self.add_setting(config, ["namelist:run_dust","l_dust_div4_lbc"], 
        l_dust_div4_lbc)
        self.add_setting(config, ["namelist:run_dust","l_dust_div5_lbc"], 
        l_dust_div5_lbc)
        self.add_setting(config, ["namelist:run_dust","l_dust_div6_lbc"], 
        l_dust_div6_lbc)
        if l_bmass_agd_lbc: 
            self.add_setting(config, ["namelist:run_aerosol","l_bmass_lbc"], 
            l_bmass_agd_lbc)
        if l_dms_lbc:
            self.add_setting(config, ["namelist:run_aerosol","l_dms_lbc"], 
            l_dms_lbc)
        if l_nh3_lbc:
            self.add_setting(config, ["namelist:run_aerosol","l_nh3_lbc"], 
            l_nh3_lbc)
        if l_nitr_acc_lbc:
            self.add_setting(config, ["namelist:run_aerosol","l_nitr_lbc"], 
            l_nitr_acc_lbc)
        if l_ocff_agd_lbc:
            self.add_setting(config, ["namelist:run_aerosol","l_ocff_lbc"], 
            l_ocff_agd_lbc)
        if l_so2_lbc:
            self.add_setting(config, ["namelist:run_aerosol","l_so2_lbc"], 
            l_so2_lbc)
        if l_soot_agd_lbc:
            self.add_setting(config, ["namelist:run_aerosol","l_soot_lbc"], 
            l_soot_agd_lbc)
        if l_consistent_cdnc:
            self.add_setting(config, ["namelist:run_radiation",
            "l_consistent_cdnc"],l_consistent_cdnc)
        if l_use_sulpc_direct:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_sulpc_direct"],l_use_sulpc_direct)
        if l_use_sulpc_indirect_lw:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_sulpc_indirect_lw"],l_use_sulpc_indirect_lw)
        if l_use_sulpc_indirect_sw:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_sulpc_indirect_sw"],l_use_sulpc_indirect_sw)
        if l_use_seasalt_direct:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_seasalt_direct"],l_use_seasalt_direct)
        if l_use_seasalt_indirect:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_seasalt_indirect"],l_use_seasalt_indirect)
        if l_use_biogenic:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_biogenic"],l_use_biogenic)
        if l_use_nitrate_direct:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_nitrate_direct"],l_use_nitrate_direct)
        if l_use_nitrate_indirect:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_nitrate_indirect"],l_use_nitrate_indirect)
        if l_use_soot_direct:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_soot_direct"],l_use_soot_direct)
        if l_use_soot_indirect:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_soot_indirect"],l_use_soot_indirect)
        if l_use_bmass_direct:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_bmass_direct"],l_use_bmass_direct)
        if l_use_bmass_indirect:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_bmass_indirect"],l_use_bmass_indirect)
        if l_use_ocff_direct:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_ocff_direct"],l_use_ocff_direct)
        if l_use_ocff_indirect:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_ocff_indirect"],l_use_ocff_indirect)
        if l_use_dust:
            self.add_setting(config, ["namelist:run_radiation",
            "l_use_dust"],l_use_dust)

        return config, self.reports


class vn86_t5981(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5981 by Joe Mancell."""
    
    BEFORE_TAG = "vn8.6_t5843"
    AFTER_TAG = "vn8.6_t5981"
    
    def upgrade(self, config, meta_config=None):
        coupler = self.get_setting_value(config, ["env", "COUPLER"])
        if coupler is None:
            self.add_setting(config, ["env", "COUPLER"], "")

        return config, self.reports

class vn86_t5434(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #5981 by Joe Mancell."""

    BEFORE_TAG = "vn8.6_t5981"
    AFTER_TAG = "vn8.6_t5434"

    def upgrade(self, config, meta_config=None):

        cntlgen_vars, ukca_vars, = {}, {},

        is_ukca = self.get_setting_value(config,
                                         ['namelist:run_ukca',
                                          'l_ukca'])
        if '.false.' in is_ukca:
            #write a default value for the interval and ignore it
            self.add_setting(config,
                             ['namelist:run_ukca', 'i_ukca_interval'],
                             value = '1',
                             state = rose.config.ConfigNode.STATE_USER_IGNORED)
        else:
            #the value kcdt was hard wired into the UM
                kcdt = 3600
                #read in required cntlgen variables
                cntlgen_vars['secs_per_period'] = self.get_setting_value(config,
                                                   ['namelist:nlstcgen',
                                                    'secs_per_periodim'])
                cntlgen_vars['steps_per_period'] =self.get_setting_value(config,
                                                    ['namelist:nlstcgen',
                                                     'steps_per_periodim'])
                ukca_vars['i_ukca_chem'] = self.get_setting_value(config,
                                            ['namelist:run_ukca',
                                             'i_ukca_chem'])
                if ukca_vars['i_ukca_chem'] in ('11','13'):
                    interval = '1'
                    self.add_setting(config,
                                     ['namelist:run_ukca', 'i_ukca_interval'],
                                     value = interval)
                elif ukca_vars['i_ukca_chem'] in ('50','51','52','53'):
                    secs_per_step = int(cntlgen_vars['secs_per_period']) / \
                        int(cntlgen_vars['steps_per_period'])
                    interval = str(kcdt/secs_per_step)
                    #write the interval variable
                    self.add_setting(config,
                                     ['namelist:run_ukca', 'i_ukca_interval'],
                                     value = interval)

        return config, self.reports


class vn86_t1389(rose.upgrade.MacroUpgrade):

    """Upgrade macro for ticket #1389 (SRS) by Paul Cresswell."""

    BEFORE_TAG = "vn8.6_t5434"
    AFTER_TAG = "vn8.6_t1389"

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


class vn90_t5995(rose.upgrade.MacroUpgrade):

    BEFORE_TAG = "vn8.6_t1389"
    AFTER_TAG = "vn9.0"

    def upgrade(self, config, meta_config=None):
        """Upgrade configuration to next major version."""
        um_version = self.get_setting_value(config, ["env", "VN"])
        if um_version:
            if um_version != "8.6":
                raise UpgradeError(
                    '!!!! Upgrade of rose app from UM version %s' % um_version  +
                    ' not supported. !!!!\nUpgrade macros are only available ' +
                    'from vn8.6 onwards.\nPlease upgrade equivalent job to '+
                    'vn8.6 using the UMUI. This may then be converted to Rose '+
                    'and upgraded to vn9.0.'
                    )
            self.change_setting_value(config, ["env", "VN"],
                                      "9.0")
        # Update STASHMSTR
        stashmaster_path = self.get_setting_value(config, ["env", "STASHMSTR"])
        if stashmaster_path:
            stashmaster_path = re.sub("vn\d+\.\d+\/ctldata", 
                                      "vn9.0/ctldata", 
                                      stashmaster_path)
            self.change_setting_value(config, ['env', 'STASHMSTR'], 
                                      stashmaster_path)
        return config, self.reports


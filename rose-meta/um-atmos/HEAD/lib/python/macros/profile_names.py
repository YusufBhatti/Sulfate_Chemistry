#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
"""This module contains code to check STASH domain and time profiles have
names compatible with the CF metadata naming convention.

For netCDF output the STASH domain and time profile names are used as
netCDF dimension and variable names and therefore need to be compatible
with the CF metadata naming convention. This means names should begin 
with a letter and be composed of letters, digits, and underscores.

"""


import re
import rose.macro
import rose.config

class ProfileNameValidate(rose.macro.MacroBase):
    """Tests for STASH profile names incompatible with netCDF output streams"""

    def validate(self, config, meta_config=None):
        """Return a list of errors, if any."""
        self.reports = []
        
        # Patterns which identify the STASH requests, usage, time, domain and 
        # netCDF output stream namelists
        nc_stream_nml_pattern = "namelist:nlstcall_nc\(.*\)$"
        streq_nml_pattern = "namelist:umstash_streq\(.*?\)$"
        use_nml_pattern = "namelist:umstash_use\(.*?\)$"
        tim_nml_pattern = "namelist:umstash_time\(.*?\)$"
        dom_nml_pattern = "namelist:umstash_domain\(.*?\)$"

        # If profile name matches this pattern then it is a valid CF string
        # Need to include surrounding ''
        cf_name_pattern = "^'[a-zA-Z0-9_]+'$"
        # If the profile name first character matches this pattern then it is a
        # valid CF string
        cf_first_char_pattern = "[a-zA-Z]"

        # Initialise dictionaries to store a sub-dictionary for time and domain
        # profiles, the sub-dictionary will have time and domain names as keys,
        # 'use_names' will store the usage profile names and 
        # 'sections' will store the section names
        use_names = {}
        use_names["time"] = {}
        use_names["domain"] = {}
        sections = {}
        sections["time"] = {}
        sections["domain"] = {}

        # Initialise a dictionary to hold file_ids with usage profile names
        # as the keys
        file_ids = {}

        # Initialise a list to hold file_ids used in netCDF output streams
        nc_stream_ids = []

        # Extract a list of the apps sections and ensure it's in order
        sect_keys = config.value.keys()
        sorter = rose.config.sort_settings
        sect_keys.sort(sorter)

        # Populate the dictionaries/list
        for section in sect_keys:
            node = config.get([section])
            if (isinstance(node.value, dict)):
            
                if (re.match(streq_nml_pattern, section)):
                    tim_name = config.get_value([section, "tim_name"])
                    dom_name = config.get_value([section, "dom_name"])
                    use_name = config.get_value([section, "use_name"])

                    # tim_name and dom_name will be undefined if STASH request
                    # is set to ignored
                    if tim_name:
                        if tim_name not in use_names["time"].keys():
                            use_names["time"][tim_name] = set()
                        use_names["time"][tim_name].add(use_name)

                    if dom_name:
                        if dom_name not in use_names["domain"].keys():
                            use_names["domain"][dom_name] = set()
                        use_names["domain"][dom_name].add(use_name)

                if (re.match(use_nml_pattern, section)):
                    file_id = config.get_value([section, "file_id"])
                    use_name = config.get_value([section, "use_name"])

                    # Note: for the usage profiles it is valid for there to be
                    # no file_id defined (if the usage requests is sending to
                    # the dump or climate meaning output)
                    if file_id:
                        file_ids[use_name] = file_id

                if (re.match(tim_nml_pattern, section)):
                    tim_name = config.get_value([section, "tim_name"])
                    sections["time"][tim_name] = section

                if (re.match(dom_nml_pattern, section)):
                    dom_name = config.get_value([section, "dom_name"])
                    sections["domain"][dom_name] = section

                if (re.match(nc_stream_nml_pattern, section)):
                    file_id = config.get_value([section, "file_id"])
                    if file_id:
                        nc_stream_ids.append(file_id)

        for profile in ["time", "domain"]:

            # Loop over time and domain profile names
            for name in use_names[profile]:

                # Check profile name is defined and get section name
                if name in sections[profile].keys():
                    section = sections[profile][name]

                    # Loop over usage profiles this particular time or domain
                    # profile is output to
                    for use_name in use_names[profile][name]:

                        # Get file_id associated with usage profile, if any
                        if use_name in file_ids.keys():
                            file_id = file_ids[use_name]

                            # Only check validity of names used in netCDF output
                            if file_id in nc_stream_ids:

                                # Check profile name is composed of letters,  
                                # digits, and underscores only
                                if (not re.match(cf_name_pattern,name)):
                                    self.add_report(section, None, None,
                                        "Invalid {0:s} profile name {1:s}."
                                        " Can only be composed of letters,"
                                        " digits, and underscores"
                                        .format(profile, name))

                                # Check profile name begins with a letter
                                # First character seems to always be a 
                                # quote (') so check second character
                                if (not re.match(cf_first_char_pattern,name[1])):
                                    self.add_report(section, None, None,
                                        "Invalid {0:s} profile name {1:s}."
                                        " First character must be a letter"
                                        .format(profile, name))

                                # Only need to check the profile name once
                                break

        return self.reports

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************
"""This module contains code to ensure that the list of STASH codes in CreateBC
is in a sensible order.

CreateBC requires that certain STASH codes are followed by certain other ones
for the code to work correctly. Currently these requirements are:

 * The first STASH code must be orography (33)
 * u-wind (2) must be followed by v-wind (3)
 * Dust Mass-mixing-ratio bin 1 (431) must be followed by bin 2 (432)
   * If bin 3 (433) is present it must follow bin 2, and be followed by bins
     4, 5 and 6 (434, 435, 436) in that order.
 * Advected u-wind (256) must be followed by advected v-wind (257)     

The UM imposes a more rigorous restriction on STASH codes, requiring LBC fields
to be in a certain fixed order. The transformer macro alters the list of STASH 
codes to be in the correct order for the UM to be able to read the LBC file.

"""

import rose.macro

# List of potential STASH codes in LBC file, in order but excluding tracers
MASTER_STASH = [ 
            33,  # Orogaphy
            2,   # u-wind 
            3,   # v-wind
            150, # w-wind 
            253, # rho
            4,   # theta
            10,  # q
            254, # qcl
            12,  # qcf
            255, # exner
            256, # adv u-wind
            257, # adv v-wind
            258, # adv w-wind
            271, # qcf2
            272, # qrain
            273, # qgraup
            266, # bulk_cf
            267, # liquid_cf
            268, # frozen_cf
            90,  # total_aero
            431, # dust mmr bin 1
            432, # dust mmr bin 2
            433, # dust mmr bin 3
            434, # dust mmr bin 4
            435, # dust mmr bin 5
            436, # dust mmr bin 6
            101, # so2
            102, # dms
            103, # mmr so4_aitken
            104, # mmr so4_accum
            105, # mmr so4_diss
            107, # mmr nh3
            108, # mmr bc_fr
            109, # mmr bc_ag
            110, # mmr bc_cl
            111, # mmr smoke fr
            112, # mmr smoke ag
            113, # mmr smoke cl
            114, # mmr ocff_fr
            115, # mmr ocff_ag
            116, # mmr ocff_cl
            117, # mmr nitr acc
            118, # mmr nitr diss
            ]

NUM_FREE_TRACERS = 150
NUM_UKCA_TRACERS = 150
STASH_SEC_FREE_TRACER = 33
STASH_SEC_UKCA_TRACER = 34
STASH_SEC = "namelist:lbc_grid"
STASH_VAL = "stash_codes"


def _populate_stash():
    """Generate the ordered list of potential STASH codes by taking the 
    MASTER_STASH list and adding any potential tracers."""
    
    stash_codes = MASTER_STASH
    for i in range(1, NUM_FREE_TRACERS):
        stash_codes.append(i+1000*STASH_SEC_FREE_TRACER)
    for i in range(1, NUM_UKCA_TRACERS):
        stash_codes.append(i+1000*STASH_SEC_UKCA_TRACER)    
    return stash_codes


class StashOrderUMTransform(rose.macro.MacroBase):
    """Force stash codes to be in the order required by the UM."""
      
    def transform(self, config, meta_config=None):
        """Reorder stash_codes into the order the UM expects."""
        
        # Populate valid STASH codes
        self.stash_codes = _populate_stash()
        requested_stash = self.stash_codes
      
        # Get user-specified STASH codes
        node = config.get([STASH_SEC, STASH_VAL])
        if node is None or node.is_ignored():
            return config, self.reports
            
        # Convert STASH to integers            
        stash_codes = [ int(x) for x in node.value.split(",")]

        missing_codes = [x for x in stash_codes if x not in self.stash_codes]
        for i in missing_codes:
            self.add_report(STASH_SEC, STASH_VAL, i, "STASH code %s not "%(i)+
                  "valid content of an LBC file" )

        # Iterate until there are no further changes
        for i in range(0, len(self.stash_codes)):
            for i in self.stash_codes:
                if i not in stash_codes:
                    requested_stash.remove(i)

        # Convert STASH back to string
        requested_stash = [str(x) for x in requested_stash]

        # Set new value in config object and add user message
        config.set([STASH_SEC, STASH_VAL], ','.join(requested_stash))
        message = 'Reordering STASH to that expected by the UM'
        self.add_report(STASH_SEC, STASH_VAL, requested_stash, message)

        return config, self.reports


#!/usr/bin/env python2.7
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************

import os
import sys
import mule
import um_wafccb
import numpy as np

# Simple argument list - input file then output file
input_file = sys.argv[1]
output_file = sys.argv[2]

# Load the input file
ff = mule.FieldsFile.from_file(input_file)

# Extract the correct fields
cpnrt = []
bulk = []
conv = []
pres = []
for field in ff.fields:
    if field.lbrel in (2, 3):
        if (field.lbuser4 == 5205
                and field.lbproc == 128
                and field.lbft == 0):
            cpnrt.append(field)
        if (field.lbuser4 == 266
                and field.lbft == 0):
            bulk.append(field)
        if (field.lbuser4 == 5212
                and field.lbft == 0):
            conv.append(field)
        if (field.lbuser4 == 408
                and field.lbft == 0):
            pres.append(field)

# Get the dimensions (assume they are the same for all fields)
rows = pres[0].lbrow
cols = pres[0].lbnpt
levs = len(pres)

# Now need to construct the arrays to pass to the CB extension...

# CPNRT is easy as it's a single field... the code has the row/column
# dimensions the other way around to the way they arrive via Mule's
# reading so the array needs to be transposed.  Also note the ordering
# is set to "F" since this will be passed to Fortran ultimately and it
# allows us to keep the index ordering in both interfaces the same
cpnrt_array = np.empty((cols, rows), order="F")
cpnrt_array[:, :] = cpnrt[0].get_data().transpose()


# Define a quick sorting function
def level_sort(field):
    return field.lblev
# and use it to make sure the multi-level fields are in order
bulk = sorted(bulk, key=level_sort)
conv = sorted(conv, key=level_sort)
pres = sorted(pres, key=level_sort)

# Create some arrays to hold the data, as above note the ordering is set to
# "F" to avoid any confusion in the indices when crossing to Fortran
bulk_array = np.empty((cols, rows, levs), order="F")
conv_array = np.empty((cols, rows, levs), order="F")
pres_array = np.empty((cols, rows, levs), order="F")

# Populate the arrays, again noting the transpose for the same reason as the
# above
for ifield, fields in enumerate(zip(bulk, conv, pres)):
    bulk_array[:, :, ifield] = fields[0].get_data().transpose()
    conv_array[:, :, ifield] = fields[1].get_data().transpose()
    pres_array[:, :, ifield] = fields[2].get_data().transpose()

# Call the extension twice - once to return pressures and once to return
# ICAO heights
icao_out = False
cb_bottom_pres, cb_top_pres, cb_extent = um_wafccb.um_wafccb.wafccb(
    cpnrt_array, bulk_array, conv_array, pres_array, cpnrt[0].bmdi, icao_out)

icao_out = True
cb_bottom_icao, cb_top_icao, _ = um_wafccb.um_wafccb.wafccb(
    cpnrt_array, bulk_array, conv_array, pres_array, cpnrt[0].bmdi, icao_out)

# Use Mule to output a FieldsFile with the result fields in it

# Copy the input file headers
ff_out = ff.copy()

# Fieldcalc used to set more than this, but I'm not going to bother; I'll
# just set the STASH to make it possible to remember which is which!
data_out = [cb_top_pres, cb_bottom_pres, cb_extent,
            cb_top_icao, cb_bottom_icao]
stash_out = [20051, 20050, 20049, 20055, 20054]

for stash, data in zip(stash_out, data_out):
    # Copy the precip field's headers
    field_out = cpnrt[0].copy()
    # Remember to switch off packing as it will have the wrong accuracy
    field_out.lbuser4 = stash
    field_out.lbpack = 0
    # Simple array provider to shove the data into the field
    field_out.set_data_provider(mule.ArrayDataProvider(data))
    ff_out.fields.append(field_out)

# Write out the file
ff_out.to_file(output_file)

#!/usr/bin/env bash
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************

# This task is designed to install the files required by the next phase 
# of the group of CRUNs in the naming test.  It should be run before 
# the "archiving" script.  NOTE: this is not a generic script, as it
# expects the names of the files to conform to certain conventions.
#
# Required environment variables:
#
#   ATMOS_NEW  - should be set to the name of the new working directory
#                (i.e. the one which the suite will create for the next
#                 CRUN phase)
#   ATMOS_PREV - should be set to the name of the old working directory
#                (i.e. the CRUN phase which has just been run)
#
set -eux
# Ensure the new directory exists
mkdir -p ../$ATMOS_NEW

# Move the *last* dump file into the new directory - this is the dump that will
# be used to initialise the next phase of the set of CRUNs
mv $(ls ../$ATMOS_PREV/atmos_dump* | tail -n1) ../$ATMOS_NEW/

# Find all output stream files (be careful to omit any .arch files - these have
# a similar name and get caught in the expansion so must be filtered back out)
for file in $(ls ../$ATMOS_PREV/atmos_pp[0-9]_stream* | grep -v '.arch$') ; do
    # If the stream file isn't due to be "archived" move it into the working
    # directory for the next phase of the set of CRUNs
    if [ ! -e $file.arch ] ; then
        mv $file ../$ATMOS_NEW/
    fi
done

# Next find all of the partial sum, mean files and the history file
for file in $(ls ../$ATMOS_PREV/atmos_psum* \
                 ../$ATMOS_PREV/atmos_[0-9][0-9]dump_mean* \
                 ../$ATMOS_PREV/atmos.xhist | grep -v '.arch$') ; do
    # If they aren't due to be "archived" these should be *copied* (not moved!)
    # into the working directory for the next phase of the set of CRUNs. 
    # (The reason for copying rather than moving is that these files get
    #  overwritten as the model updates them - therefore to be able to restart
    #  the next phase of the CRUN a clean copy of these files must exist at the
    #  ending state of the previous phase)
    if [ ! -e $file.arch ] ; then
        cp $file ../$ATMOS_NEW/
    fi
done

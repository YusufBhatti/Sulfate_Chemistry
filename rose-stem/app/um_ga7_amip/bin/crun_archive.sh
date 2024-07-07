#!/usr/bin/env bash
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************

# Thsi task is designed to "archive" the files from a single phase 
# of the group of CRUNs in the naming test.  It should be run after 
# the "install" script (which installs the next phase's files).  
# NOTE: this is not a generic script, as it expects the names of 
# the files to conform to certain conventions.
#
# Required environment variables:
#
#   ATMOS_DIR   - should be set to the name of the CRUN working directory
#                 (i.e. the CRUN phase which has just been run)
#   ARCHIVE_DIR - should be set to the name of a directory where the 
#                 "archived" files should be copied
#   CRUN_NO     - should be set to the CRUN phase number (this is used to
#                 allow the script to treat the final phase differently)
#
set -eux
# Ensure the "archive" directory exists
mkdir -p $ARCHIVE_DIR

# Find any files in the CRUN task's working directory which have ".arch"
# files associated with them
for file in $(ls ../$ATMOS_DIR/*) ; do
    if [ -e $file.arch ] ; then
        # Move these files to the "archive" directory and remove the flag file
        mv ../$ATMOS_DIR/$file $ARCHIVE_DIR
        rm $file.arch
    fi
done

# Since dumps don't have the ".arch" files just move any found (note that this
# is why the script must run *after* the "install" script - if it ran before it
# would move the final dump to the archive directory before it was installed)
for file in $(ls ../$ATMOS_DIR/atmos_dump*) ; do
    mv $file $ARCHIVE_DIR
done

# For the final phase of the set of CRUNs, also pick up any leftover output
# stream files (some streams will continue indefinitely and therefore never
# produce a ".arch" file)
if [[ $CRUN_NO == 3 ]] ; then
   mv $(ls ../$ATMOS_DIR/atmos_pp[0-9]_stream*) $ARCHIVE_DIR 
fi


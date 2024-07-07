#!/usr/bin/env bash
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************

DESTINATION=$UM_INSTALL_DIR/vn$VN/$PLATFORM/utilities
COPY_CMD=${COPY_CMD:-cp}
echo Installing utilities to $DESTINATION

mkdir -p $DESTINATION
ERR_STATE=$?

if [[ $ERR_STATE != 0 ]] ; then
echo Unable to create directory $DESTINATION
exit 1
fi

for SOURCE in $SOURCES; do
$COPY_CMD $CYLC_SUITE_SHARE_DIR/$SOURCE/build-*/bin/* $DESTINATION
done

# We exit 0 here to ignore warnings about duplicated um_script_functions
exit 0

#!/usr/bin/env bash
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************

DESTINATION=$UM_INSTALL_DIR/vn$VN/ctldata
SOURCE=$CYLC_SUITE_SHARE_DIR/$TASK_DIR/extract
COPY_CMD=${COPY_CMD:-cp}
echo Installing ctldata to $DESTINATION

mkdir -p $DESTINATION
ERR_STATE=$?

if [[ $ERR_STATE != 0 ]] ; then
echo Unable to create directory $DESTINATION
exit 1
fi

$COPY_CMD -r $SOURCE/ctldata/* $DESTINATION
$COPY_CMD -r $SOURCE/stash/rose-meta/um-atmos/$INSTALL_META_VERSION/etc/stash/STASHmaster $DESTINATION
$COPY_CMD -r $SOURCE/stash/rose-meta/um-atmos/$INSTALL_META_VERSION/etc/STASH2CF $DESTINATION

# symlink spectral files as extract system does not create symlinks
cd $DESTINATION/spectral/ga3_0/
if [ ! -e "mcica_data" ] ; then
ln -s mcica_data_ga3_0 mcica_data
fi
cd $DESTINATION/spectral/ga3_1/
if [ ! -e "mcica_data" ] ; then
ln -s mcica_data_ga3_1 mcica_data
fi

# Check if ancil dir exists. If ancil dir does not exist and this is a central installation then
# report an error. Otherwise link in the ancil dir from the central installation. 
cd $UM_INSTALL_DIR
if [ ! -d "ancil" ] ; then
    if [ $UM_INSTALL_DIR = $UMDIR ] ; then
	# Send the error message to both stdout and stderr
	echo $UMDIR/ancil directory does not exist
	echo $UMDIR/ancil directory does not exist >&2
       exit 1
    else
	ln -s $UMDIR/ancil ancil
    fi
fi


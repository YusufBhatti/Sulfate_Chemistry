#!/usr/bin/env bash
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************

set -eux

DESTINATION=$UM_INSTALL_DIR/vn$VN/$PLATFORM
COPY_CMD=${COPY_CMD:-cp}
LIB_SUFFIX=${LIB_SUFFIX:-}

# Make the base destination directory
mkdir -p $DESTINATION
if [[ $? != 0 ]] ; then
    echo Unable to create directory $DESTINATION
    exit 1
fi

# Copy the sstpert library (if present)
SSTPERT_DESTINATION=$DESTINATION/sstpert${LIB_SUFFIX}
SSTPERT_BUILD=$CYLC_SUITE_SHARE_DIR/$TASK_DIR/build-sstpert_lib
if [[ -d $SSTPERT_BUILD ]] ; then
  echo "Installing UM sstpert library to $SSTPERT_DESTINATION/lib"
  mkdir -p $SSTPERT_DESTINATION/lib
  $COPY_CMD $SSTPERT_BUILD/lib/libum_sstpert.so $SSTPERT_DESTINATION/lib

  echo "Installing UM sstpert library includes to $SSTPERT_DESTINATION/include"
  mkdir -p $SSTPERT_DESTINATION/include
  EXTRACT_DIR=$CYLC_SUITE_SHARE_DIR/$TASK_DIR/extract/um/src/include/other
  $COPY_CMD $EXTRACT_DIR/sstpert.h \
            $SSTPERT_DESTINATION/include
fi

# Copy the wafccb library (if present)
WAFCCB_DESTINATION=$DESTINATION/wafccb${LIB_SUFFIX}
WAFCCB_BUILD=$CYLC_SUITE_SHARE_DIR/$TASK_DIR/build-wafccb_lib
if [[ -d $WAFCCB_BUILD ]] ; then
  echo "Installing UM wafccb library to $WAFCCB_DESTINATION/lib"
  mkdir -p $WAFCCB_DESTINATION/lib
  $COPY_CMD $WAFCCB_BUILD/lib/libum_wafccb.so $WAFCCB_DESTINATION/lib

  echo "Installing UM wafccb library includes to $WAFCCB_DESTINATION/include"
  mkdir -p $WAFCCB_DESTINATION/include
  EXTRACT_DIR=$CYLC_SUITE_SHARE_DIR/$TASK_DIR/extract/um/src/include/other
  $COPY_CMD $EXTRACT_DIR/wafccb.h \
            $WAFCCB_DESTINATION/include
fi
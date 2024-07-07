#!/usr/bin/env bash
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************

set -eux

SUFFIX=${SUFFIX:-}

# Switch to the extract directory
EXTRACT_DIR=$CYLC_SUITE_SHARE_DIR/$TASK_DIR/extract/mule
echo "[INFO] Switching to $EXTRACT_DIR"
cd $EXTRACT_DIR

if [ -z "${SHUMLIB_LIB:-}" ] ; then
    echo "[FAIL] Shumlib path not set"
    exit 1
fi

# Setup install script arguments
LIB_DEST=${LIB_DEST:-$UM_INSTALL_DIR/mule_lib}${SUFFIX}
BIN_DEST=${BIN_DEST:-$UM_INSTALL_DIR/mule_bin}${SUFFIX}

# Add the flags for the optional libraries based on whether
# they have been installed
OPTIONAL_FLAGS=
SSTPERT_LIB=$UM_INSTALL_DIR/vn$VN/$PLATFORM/sstpert${SUFFIX}
if [ -d $SSTPERT_LIB ] ; then
    OPTIONAL_FLAGS="$OPTIONAL_FLAGS --sstpert_lib $SSTPERT_LIB"
fi
WAFCCB_LIB=$UM_INSTALL_DIR/vn$VN/$PLATFORM/wafccb${SUFFIX}
if [ -d $WAFCCB_LIB ] ; then
    OPTIONAL_FLAGS="$OPTIONAL_FLAGS --wafccb_lib $WAFCCB_LIB"
fi

# Clear any previous installations
rm -rf $LIB_DEST $BIN_DEST

# Run the install script
echo "[INFO] Installing Mule project..."
./admin/install_mule.sh $OPTIONAL_FLAGS $LIB_DEST $BIN_DEST $SHUMLIB_LIB


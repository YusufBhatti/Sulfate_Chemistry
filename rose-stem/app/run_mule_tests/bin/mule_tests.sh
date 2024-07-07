#!/usr/bin/env bash
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************

set -eux

SUFFIX=${SUFFIX:-}
SSTPERT_INPUT=${SSTPERT_INPUT:-}
WAFCCB_INPUT=${WAFCCB_INPUT:-}

# The directories where mule should have been installed
LIB_DEST=${LIB_DEST:-$UM_INSTALL_DIR/mule_lib}${SUFFIX}
BIN_DEST=${BIN_DEST:-$UM_INSTALL_DIR/mule_bin}${SUFFIX}

# Python executable to use (this is copied from what is done by the
# install script, in case an alternative python install is needed
PYTHONEXEC=${PYTHONEXEC:-python2.7}

# Allow python to pickup the modules installed by this suite
export PYTHONPATH=$LIB_DEST
export PATH=$PATH:$BIN_DEST

# And finally run all 3 sets of testing (redirect stderr to stdout for these,
# since unittest insists on writing everything to stderr for some reason)
for module in um_packing mule um_utils ; do
    echo "[INFO] Test '$module' module..."
    $PYTHONEXEC -s -c "import $module ; print ($module)"
    $PYTHONEXEC -s -m unittest discover -v $module.tests 2>&1
done

# The SSTpert library test is a little different; it has no built-in tests
# so we instead use what is essentially a traditional KGO test
if [ -n "$SSTPERT_INPUT" ]; then
    ulimit -s unlimited
    mule-sstpert $SSTPERT_INPUT 2.5 201505141200 1 sstpert_out.ff
else
    echo "[WARN] mule-sstpert is not being tested."
    echo "[INFO] To turn on testing of mule-sstpert set SSTPERT_INPUT"
fi

# Similarly, the WAFC CB library test has no built in tests, and also no
# entry-point script (as SSTpert does above) so we'll instead run a simple
# script which is stored alongside this one
if [ -n "$WAFCCB_INPUT" ]; then
    $PYTHONEXEC $(dirname $0)/test_wafccb.py $WAFCCB_INPUT wafccb_out.ff
else    
    echo "[WARN] Mule's WAFC CB library is not being tested."
    echo "[INFO] To turn on testing of WAFC CB set WAFCCB_INPUT"
fi
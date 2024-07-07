#!/usr/bin/env bash
# *********************************COPYRIGHT************************************
# (C) Crown copyright Met Office. All rights reserved. 
# For further details please refer to the file COPYRIGHT.txt 
# which you should have received as part of this distribution. 
# *********************************COPYRIGHT************************************

# This script is designed to test an all-in-one upgrade of an app from a
# prior UM release to the head of the trunk. It runs the upgrade macros in
# one jump, whereas the normal rose-stem apps are upgraded incrementally 
# as each goes on the trunk. There are a few weird circumstances in which
# the incremental step-by-step approach and the all-in-one jump can give 
# different answers; this script is designed to detect those.

APPS=${OVERRIDE_APP_DIRS:-$APP_DIRS}
TRUNK=${OVERRIDE_TRUNK:-fcm:um.xm_tr}
METAPATH=$CYLC_TASK_WORK_DIR/../script_source/um/rose-meta
WCPATH=$CYLC_TASK_WORK_DIR/../script_source/um/rose-stem/app
retcode=0


# Loop over apps
for app in $APPS; do
  # Separate appname and version number
  arrapp=(${app/@/ })
  appname=${arrapp[0]}
  version=${arrapp[1]}
  mkdir -p "${appname}_${version}"
  cd "${appname}_${version}"

  # Retrieve the base app at a given UM version
  url="${TRUNK}/rose-stem/app/${appname}/rose-app.conf@${version}"
  fcm export --force -q $url
  if [[ $? != 0 ]]; then
    echo "[FAIL] Failed to export $url"
    echo "[FAIL] Failed to export $url" 1>&2
    retcode=1
    cd $CYLC_TASK_WORK_DIR
    continue
  fi

  # Work out the aftertag for this app
  aftertag=`rose app-upgrade -a -M ${METAPATH} | tail -n 1 | awk '{print $2}' | sed -e 's/_tXXXX//'`

  # Perform the upgrade
  echo "[INFO] Testing upgrade of $appname from $version to $aftertag"
  rose app-upgrade -q -a -M ${METAPATH} -y ${aftertag}

  # Does the upgraded app match the metadata? Only test UM/recon apps.
  rose macro --no-warn version -M ${METAPATH} -q --validate 1>&2
  if [[ $? == 0 ]]; then
    echo "[PASS] Upgraded app $appname from $version validates"
  else
    retcode=1
    echo "[FAIL] Upgraded app $appname from $version does not validate"
    echo "[FAIL] Upgraded app $appname from $version does not validate" 1>&2
  fi

  # Diff the upgrade app with the working copy version
  diff rose-app.conf $WCPATH/$appname/rose-app.conf 1>&2
  if [[ $? == 0 ]]; then
    echo "[PASS] Upgraded $appname from $version matches current working copy"
  else
    retcode=1
    echo "[FAIL] Upgraded $appname from $version does not match current working copy"
    echo "[FAIL] Upgraded $appname from $version does not match current working copy" 1>&2
  fi
  cd $CYLC_TASK_WORK_DIR
done

echo "[DONE] Script completed."
exit $retcode

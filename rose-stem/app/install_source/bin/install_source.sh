#!/usr/bin/env bash
# *****************************COPYRIGHT******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT******************************

# For each project, create a new repository in $HOME/source/$project/offline
# containing the source code obtained by Rose.

DESTINATION=$HOME/source

mkdir -p $DESTINATION
ERR_STATE=$?
if [[ $ERR_STATE != 0 ]] ; then
  echo Unable to create directory $DESTINATION
  exit 1
fi
echo Installing source in $DESTINATION

for project in um jules socrates casim um_aux mule; do

  mkdir -p $DESTINATION/$project

  cd $DESTINATION/$project

  rm -rf offline
  svnadmin create offline
  svn import -q $CYLC_TASK_WORK_DIR/$project file://$DESTINATION/$project/offline/trunk -m "Create $project offline trunk"

done

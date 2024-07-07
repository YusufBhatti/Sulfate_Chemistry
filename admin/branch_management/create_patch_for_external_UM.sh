#!/usr/bin/env bash
#
# Script to create a patch from the central Met Office UM repository for
# application to an external installation of the UM.
#
# Original Author - Roddy Sharp, Met Office
#
# To use this script.
# 1. Either copy the script giving it a name which clearly indicates the
#    purpose of the patch, or edit the line defining OUTDIR.
# 2. Change OUTLOC to a suitable location where you have write permissions.
# 3. Change the URL to the correct one for the changes you wish to replicate.
# 4. Change the two revision numbers. Remember to add at least 1 to the revision
#    that created your branch if you're sending a patch to create a branch.
#    The revision numbers can be Keywords if desired.
# 5. Check the fcm-import-patch for any reference to UKCA or isccp. If there is one, there
#    may be an error please contact the PUM manager, or the UM Systems development Team, or
#    the external collaborations team.

# OUTDIR is the directory you want the patch output to.
# If not supplied it defaults to fcm-makepatch-out
OUTDIR=`basename $0 | sed -e ' s/create_patch_// ' -e ' s/\.sh// ' `_patch
echo "\$OUTDIR is $OUTDIR"

# OUTLOC is the location you want the output directory in.
OUTLOC=$HOME

# SVN_URL is the URL of the trunk or branch that includes the changes
# you want the patch to be able to replicate.
SVN_URL="svn://fcm2/UM_svn/UM/trunk"

# START_REV is the revision you want to capture changes from. If you're creating
# a patch of a branch, it should at least be the revision number one greater
# than the revision that created the actual branch in the repository
START_REV=458

# END_REV is the revision you want to capture changes up to. 
# For example, the revision number of the last relevent change to your branch.
END_REV=vn6.4

#===============================================================================
#          ** Editing below this line should not be necessary **
#===============================================================================
echo "$(date): Creating patch ..."
echo "         ${OUTLOC}/$OUTDIR"
cd $OUTLOC
fcm mkpatch --organisation Met_Office \
            --exclude src/atmosphere/UKCA \
            --exclude src/atmosphere/short_wave_radiation/isccp.F90 \
            --exclude src/atmosphere/short_wave_radiation/isccp_cloud_types.F90 \
            --exclude src/atmosphere/short_wave_radiation/isccp_cloudtypes_fld.F90 \
            --revision ${START_REV}:${END_REV} \
            $SVN_URL \
            $OUTDIR

RET_CODE=$?

if [[ "$RET_CODE" == 0 ]] ; then
  echo "$(date): Patch created in $OUTDIR"
else
  echo "#--------------### WARNING ###-------------------#"
  echo "#    fcm mkpatch did not complete sucessfully    #"
  echo "#-------------#################------------------#"
  exit 1
fi

# 6. To apply the patch created the command is
#    <$OUTLOC>/<$OUTDIR>/fcm-import-patch <URL>
#    URL should be a valid repository URL or a working copy.
#    When URL is the repository it should include the path to what the patch is
#    being applied to. e.g. /UM/trunk, or /UM/branches/dev/<user>/<branchname>
#    In the case of applying the patch to a branch, the branch needs to have
#    been created first.
#    The fcm gui can be used to do this. e.g.
#
#    $ fcm gui
#   
#    [[= In the gui =
#    select branch tab
#    Type the URL in the relevant text box (include <project>/branches)
#    Type the branch name in the relevant text box
#    Type the correct revision (probably head)
#    Check the other details are correct
#    select 'run'
#    = End of gui use =]]
#
#    You could use the create_branch util if available or the command line
#    if you're comfortable with that.


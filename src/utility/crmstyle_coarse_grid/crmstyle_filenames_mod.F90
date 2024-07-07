! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Environmental filenames

MODULE crmstyle_filenames_mod

USE filenamelength_mod, ONLY:                                           &
  filenamelength

IMPLICIT NONE
SAVE

! Description:
!   Module containing information on filenames obtained from environmental
!  variables. Also contains program name and environment variables for PP/ff
!  files.
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CRMstyle_coarse_grid
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3  programming standards.
!
! Declarations:

INTEGER, PARAMETER :: max_ff = 25        ! maximum number of fieldsfiles

CHARACTER(LEN=*), PARAMETER :: ProgName = "crmstyle_coarse_grid"

! Environment variable names for pp/ff input files

CHARACTER(LEN=15), PARAMETER ::                                               &
  ff_env_name(max_ff)=(/"PP_INPUT_FILE01","PP_INPUT_FILE02","PP_INPUT_FILE03",&
                     "PP_INPUT_FILE04","PP_INPUT_FILE05","PP_INPUT_FILE06",   &
                     "PP_INPUT_FILE07","PP_INPUT_FILE08","PP_INPUT_FILE09",   &
                     "PP_INPUT_FILE10","PP_INPUT_FILE11","PP_INPUT_FILE12",   &
                     "PP_INPUT_FILE13","PP_INPUT_FILE14","PP_INPUT_FILE15",   &
                     "PP_INPUT_FILE16","PP_INPUT_FILE17","PP_INPUT_FILE18",   &
                     "PP_INPUT_FILE19","PP_INPUT_FILE20","PP_INPUT_FILE21",   &
                     "PP_INPUT_FILE22","PP_INPUT_FILE23","PP_INPUT_FILE24",   &
                     "PP_INPUT_FILE25"/)

! filenames

CHARACTER(LEN=filenamelength) :: input_file  = ""
CHARACTER(LEN=filenamelength) :: orogfile    = ""
CHARACTER(LEN=filenamelength) :: landseafile = ""
CHARACTER(LEN=filenamelength) :: lev_nl_file = ""
CHARACTER(LEN=filenamelength) :: pp_file(max_ff) = ""

! Output filenames
CHARACTER(LEN=filenamelength) :: all_file = ""
CHARACTER(LEN=filenamelength) :: acc_file = ""
CHARACTER(LEN=filenamelength) :: acu_file = ""
CHARACTER(LEN=filenamelength) :: bcu_file = ""
CHARACTER(LEN=filenamelength) :: wg1_file = ""
CHARACTER(LEN=filenamelength) :: ppd_file = ""
CHARACTER(LEN=filenamelength) :: nbd_file = ""
CHARACTER(LEN=filenamelength) :: nid_file = ""
CHARACTER(LEN=filenamelength) :: adu_file = ""
CHARACTER(LEN=filenamelength) :: acw_file = ""
CHARACTER(LEN=filenamelength) :: bcw_file = ""
CHARACTER(LEN=filenamelength) :: bcu_mask_file = ""

END MODULE crmstyle_filenames_mod

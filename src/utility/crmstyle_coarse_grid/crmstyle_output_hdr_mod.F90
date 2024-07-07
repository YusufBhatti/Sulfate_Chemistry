! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Holds Fields file Headers for output files

MODULE crmstyle_output_hdr_mod

USE IO_Mod, ONLY:         &
  UM_Header_type

IMPLICIT NONE

SAVE

! Description:
!   Module holding the output headers for the fieldsfiles
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CRMstyle_coarse_grid
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3  programming standards.
!
! Declarations:

TYPE(UM_Header_type) ::  &
  all_hdr                & ! UM Headers:  all_file
 ,acc_hdr                & ! UM Headers:  acc_file
 ,acu_hdr                & ! UM Headers:  acu_file
 ,bcu_hdr                & ! UM Headers:  bcu_file
 ,wg1_hdr                & ! UM Headers:  wg1_file
 ,ppd_hdr                & ! UM Headers:  ppd_file
 ,nbd_hdr                & ! UM Headers:  nbd_file
 ,nid_hdr                & ! UM Headers:  nid_file
 ,adu_hdr                & ! UM Headers:  adu_file
 ,acw_hdr                & ! UM Headers:  acw_file
 ,bcw_hdr                & ! UM Headers:  bcw_file
 ,ucu_hdr                & ! UM Headers:  ucu_file
 ,ppw_hdr                & ! UM Headers:  ppw_file
 ,nbw_hdr                & ! UM Headers:  nbw_file
 ,single_hdr             & ! UM Headers:  single_file
 ,bcu_mask_hdr             ! UM Headers:  bcu_mask_file

END MODULE crmstyle_output_hdr_mod

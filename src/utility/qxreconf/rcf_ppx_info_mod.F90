! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Defines the STASHmaster record

MODULE Rcf_Ppx_Info_Mod

! Description:
!   Defines the STASHmaster record format, the USTSNUM namelist
!   and other related variables
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


USE Submodel_Mod, ONLY:     &
    N_Internal_Model_Max

USE ppxlook_mod, ONLY: &
    ppxref_sections,          &
    ppxref_items

IMPLICIT NONE

! A parameter needed for the following definition.
! Number of packing codes in stashmaster record.
INTEGER, PARAMETER        :: Ppxref_pack_profs = 10

! Length of Option and Version Codes.
INTEGER, PARAMETER        :: STM_OptCodeLen = 30 !Must be multiple of 5
INTEGER, PARAMETER        :: STM_VerCodeLen = 20

! For broadcast: quantity of integer/character data in each STMrecord.
! Must be changed if STMrecord format changes.
INTEGER, PARAMETER        :: STM_IntDataLen  = 28 + Ppxref_pack_profs &
                                             + STM_OptCodeLen / 5
INTEGER, PARAMETER        :: STM_CharDataLen = 37

! Define STMrecord type - to hold STASHmaster records
! Note that as these are in a sequence and are integer followed by
! character, we can do some fairly efficient comms with them.
! Thus order *is* vital!
TYPE STM_record_type
  SEQUENCE
  INTEGER                 :: model
  INTEGER                 :: section
  INTEGER                 :: item
  INTEGER                 :: space_code
  INTEGER                 :: ptr_code
  INTEGER                 :: timavail_code
  INTEGER                 :: grid_type
  INTEGER                 :: lv_code
  INTEGER                 :: lb_code
  INTEGER                 :: lt_code
  INTEGER                 :: pt_code
  INTEGER                 :: pf_code
  INTEGER                 :: pl_code
  INTEGER                 :: lev_flag
  INTEGER                 :: opt_code( STM_OptCodeLen / 5 )
  INTEGER                 :: version_mask
  INTEGER                 :: halo_type
  INTEGER                 :: data_type
  INTEGER                 :: dump_packing
  INTEGER                 :: packing_acc(Ppxref_pack_profs)
  INTEGER                 :: rotate_code
  INTEGER                 :: field_code
  INTEGER                 :: user_code
  INTEGER                 :: lbvc_code
  INTEGER                 :: base_level
  INTEGER                 :: top_level
  INTEGER                 :: ref_lbvc_code
  INTEGER                 :: cf_levelcode
  INTEGER                 :: cf_fieldcode
  INTEGER                 :: RowIndex

  CHARACTER (LEN=36)      :: NAME
  CHARACTER (LEN=1)       :: OriginFlag
END TYPE STM_record_type

!-------------------------------------------------------------------
! Guts for the storage of STASHmaster information
!-------------------------------------------------------------------
INTEGER, TARGET, SAVE    :: ppxRecs = 0     ! No. of stash records

TYPE (STM_record_type), POINTER, SAVE   :: STM_record(:)

! Referencing for the above
INTEGER, TARGET, SAVE    :: ppxptr( N_Internal_Model_Max, 0:Ppxref_Sections , &
                                    Ppxref_Items ) = 0

END MODULE Rcf_Ppx_Info_mod

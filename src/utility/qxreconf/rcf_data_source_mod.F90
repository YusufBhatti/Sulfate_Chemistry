! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Definitions for data-type defining source of a field

MODULE Rcf_data_source_Mod

! Description:
!   Data module defining the data_source_type data-type and
!   related magic numbers.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE filenamelength_mod, ONLY: &
    filenamelength

IMPLICIT NONE

TYPE data_source_type

  INTEGER                          :: Source
  INTEGER                          :: Domain
  INTEGER                          :: Ancil_SctnC
  INTEGER                          :: Ancil_ItemC
  REAL                             :: RConst
  CHARACTER (LEN=filenamelength)   :: Ancil_File  ! user supplied
  CHARACTER (LEN=filenamelength)   :: ancilfile   ! central ancils

END TYPE data_source_type

! Almost magic number currently only used internally.
! defined for fields where there is logic in rcf_locate_alt_field
INTEGER, PARAMETER     :: Other_field           = -90

! Magic number for internal use only
INTEGER, PARAMETER     :: Already_Processed = -99

END MODULE Rcf_data_source_Mod

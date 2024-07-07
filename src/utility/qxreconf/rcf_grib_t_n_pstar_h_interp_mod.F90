! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Module used to hold data to allow recall later in execution tree.

MODULE Rcf_GRIB_T_n_Pstar_H_Interp_Mod

! Description:
!   Stores the horizontally interpolated versions of T and Pstar used
!   to generate heights when reconfiguring from ECMWF pressure levels.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE Rcf_Field_Type_Mod, ONLY:  &
  field_type

IMPLICIT NONE

TYPE (field_type), TARGET, SAVE :: grib_tv
TYPE (field_type), TARGET, SAVE :: GRIB_Pstar

REAL, TARGET, ALLOCATABLE, SAVE :: GRIB_Levels(:)

! Store ecmwf model level definitions
REAL, TARGET, ALLOCATABLE, SAVE :: ak(:)
REAL, TARGET, ALLOCATABLE, SAVE :: bk(:)
REAL, TARGET, ALLOCATABLE, SAVE :: akh(:)
REAL, TARGET, ALLOCATABLE, SAVE :: bkh(:)

END MODULE Rcf_GRIB_T_n_Pstar_H_Interp_Mod

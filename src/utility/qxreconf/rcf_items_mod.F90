! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Arrays for the items namelists

MODULE Rcf_Items_Mod

! Description:
!    Arrays for the items namelists
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE filenamelength_mod, ONLY: filenamelength
USE items_nml_mod,      ONLY: netcdf_varname_len


IMPLICIT NONE

INTEGER, ALLOCATABLE :: Source_Array(:) ! Specify source of output field
INTEGER, ALLOCATABLE :: Sctn_Array(:) ! STASH section number
INTEGER, ALLOCATABLE :: Item_Array(:) ! STASH item number
INTEGER, ALLOCATABLE :: area_array(:) ! Domain
INTEGER, ALLOCATABLE :: upas_array(:) ! User prognostic STASH section number
INTEGER, ALLOCATABLE :: upaa_array(:) ! User prognostic STASH item number
REAL,    ALLOCATABLE :: uprc_array(:) ! Constant value to reset field to
! Ancil filename
CHARACTER (LEN=filenamelength),     ALLOCATABLE :: upaf_array(:)
! NetCDF varname
CHARACTER (LEN=netcdf_varname_len), ALLOCATABLE :: upnv_array(:)

INTEGER, SAVE :: num_items            ! number of items in items namelists

END MODULE Rcf_Items_Mod

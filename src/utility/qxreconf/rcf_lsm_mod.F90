! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Land-Sea mask data

MODULE Rcf_Lsm_Mod

! Description:
!   Stores land sea masks and related sizes/data & coastal adjustment
!   gather indexes etc.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE items_nml_mod, ONLY: &
    input_dump

USE missing_data_mod, ONLY: &
    rmdi

IMPLICIT NONE

! Integer to indicate source of LSM. Default to input dump
! and change if specified by the items namelist.
INTEGER :: lsm_source = input_dump

! Real to store a constant value for the LSM.
REAL    :: lsm_fixed_value = rmdi

! Logical to indicate if LSM is present in output dump
LOGICAL :: l_lsm_out_present = .FALSE.

! Storage for land-sea masks
LOGICAL, ALLOCATABLE, TARGET, SAVE   :: glob_lsm_in(:)
LOGICAL, ALLOCATABLE, TARGET, SAVE   :: local_lsm_in(:)
LOGICAL, ALLOCATABLE, TARGET, SAVE   :: glob_lsm_out(:)
LOGICAL, ALLOCATABLE, TARGET, SAVE   :: local_lsm_out(:)
INTEGER, TARGET, SAVE                :: local_land_in
INTEGER, TARGET, SAVE                :: local_land_out
INTEGER, TARGET, SAVE                :: glob_land_in
INTEGER, TARGET, SAVE                :: glob_land_out

! Pointers to active lsm fields
LOGICAL, POINTER, SAVE       :: glob_atmos_landmask(:)
LOGICAL, POINTER, SAVE       :: local_atmos_landmask(:)
INTEGER, POINTER, SAVE       :: local_land_field
INTEGER, POINTER, SAVE       :: glob_land_field

! Storage for gather indexes etc for coastal adjustment
INTEGER, SAVE                :: N_Coastal_Points
INTEGER, SAVE                :: N_land_points_unres
INTEGER, SAVE                :: N_sea_points_unres
INTEGER, ALLOCATABLE, SAVE   :: Coast_Index_In(:)
INTEGER, ALLOCATABLE, SAVE   :: Coast_Index_Out(:)
INTEGER, ALLOCATABLE, SAVE   :: Index_Targ_Land(:)
INTEGER, ALLOCATABLE, SAVE   :: Index_Targ_Sea(:)
INTEGER, ALLOCATABLE, SAVE   :: Land_Unres_Index(:)
INTEGER, ALLOCATABLE, SAVE   :: Sea_Unres_Index(:)

! For spiral circle method
INTEGER, ALLOCATABLE, SAVE   :: Land_Unres_Constrain_Index(:)
INTEGER, ALLOCATABLE, SAVE   :: Land_Unres_NotConstrain_Index(:)
INTEGER, ALLOCATABLE, SAVE   :: Sea_Unres_NotConstrain_Index(:)

! Logical for coast adjustment
LOGICAL, SAVE                :: Cyclic

END MODULE Rcf_Lsm_Mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE stashmaster_constants_mod

! Import variables from UM modules and rename to aid readability
USE cppxref_mod, ONLY:                               &
    ! Field data type indicators (DataT)
    data_type_real => ppx_type_real,                 &
    data_type_integer => ppx_type_int,               &
    data_type_logical => ppx_type_log,               &
    ! Address for vertical co-ord type (LBVC)
    vertical_coord_hybrid_height => ppx_lbvc_hybrid, &
    ! Address for horiz grid type code (Grid)
    horiz_grid_type => ppx_grid_type,                &
    ! Address for level type (vert grid) code (LevelT)
    vert_level_type => ppx_lv_code,                  &
    ! Indicators for horiz grid type
    u_points => ppx_atm_cuall,                       &
    v_points => ppx_atm_cvall,                       &
    p_points => ppx_atm_tall,                        &
    ! Address of Meto8  field code (CFFF)
    meto8_code => ppx_meto8_fieldcode,               &
    ! Address of first and last level codes (LevelF  & LevelL)
    first_level_code => ppx_lb_code,                 &
    last_level_code => ppx_lt_code,                  &
    ! Address of halo code (Halo)
    halo_type_code => ppx_halo_type,                 &
    ! Address of PP field code         
    pp_field_code => ppx_field_code,                 &
    ! Address of packing accuracy
    packing_acc => ppx_packing_acc

USE stparam_mod, ONLY:                               &
    ! Indicators for level type
    theta_levels => st_levels_model_theta,           &
    rho_levels => st_levels_model_rho,               &
    single_level => st_levels_single

IMPLICIT NONE

! Description:
!   Contains the numeric values corresponding to quantities defined in the
!   STASHmaster file, documented in UMDP F3 and C4.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

END MODULE stashmaster_constants_mod

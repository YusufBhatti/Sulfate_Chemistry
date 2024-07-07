! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Module to define a grid data type.

MODULE Rcf_Grid_Type_Mod

! Description:
!   Defines the grid_type data type
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

TYPE grid_type

  ! General attributes
  LOGICAL      :: global       ! global == .true., LAM ==.false.
  LOGICAL      :: rotated      ! .true. if rotated LAM
  INTEGER      :: grid_stagger ! As read from fixed header - same values.


  ! Grid dimensions - horizontally
  ! 1st the Global variables (in MPP sense)
  INTEGER      :: glob_p_row_length   ! length of 'pressure' rows
  INTEGER      :: glob_u_row_length   ! length of U rows
  INTEGER      :: glob_v_row_length   ! length of V rows
  INTEGER      :: glob_r_row_length   ! length of River rows
  INTEGER      :: glob_p_rows         ! No. of 'pressure' rows
  INTEGER      :: glob_u_rows         ! No. of U rows
  INTEGER      :: glob_v_rows         ! No. of V rows
  INTEGER      :: glob_r_rows         ! No. of River rows
  INTEGER      :: glob_p_field        ! Total size of P field
  INTEGER      :: glob_u_field        ! Total size of U field
  INTEGER      :: glob_v_field        ! Total size of V field
  INTEGER      :: glob_r_field        ! Total size of River field
  INTEGER      :: glob_land_field     ! No. of land points in field

  ! Now the Local variables (in MPP sense)
  INTEGER      :: loc_p_row_length   ! length of 'pressure' rows
  INTEGER      :: loc_u_row_length   ! length of U rows
  INTEGER      :: loc_v_row_length   ! length of V rows
  INTEGER      :: loc_r_row_length   ! length of River rows
  INTEGER      :: loc_p_rows         ! No. of 'pressure' rows
  INTEGER      :: loc_u_rows         ! No. of U rows
  INTEGER      :: loc_v_rows         ! No. of V rows
  INTEGER      :: loc_r_rows         ! No. of River rows
  INTEGER      :: loc_p_field        ! Total size of P field
  INTEGER      :: loc_u_field        ! Total size of U field
  INTEGER      :: loc_v_field        ! Total size of V field
  INTEGER      :: loc_r_field        ! Total size of River field
  INTEGER      :: loc_land_field     ! No of land points in field

  ! Grid dimensions - vertically
  ! These are *all* Global...
  INTEGER      :: model_levels ! No. of pressure levels
  INTEGER      :: cloud_levels ! No. of cloud levels
  INTEGER      :: st_levels    ! No. of soil temp. levels
  INTEGER      :: sm_levels    ! No. of soil moisture levels
  INTEGER      :: bl_levels    ! No. of boundary-layer levels
  INTEGER      :: ozone_levels ! No. of ozone levels
  INTEGER      :: tr_levels    ! No. of tracer levels
  INTEGER      :: conv_levels  ! No. of convective levels

  ! Stuff for calculation of heights
  INTEGER       :: height_gen_method   ! method used to generate heights
  INTEGER       :: first_constant_r_rho_level
  REAL          :: z_top_of_model
  REAL, POINTER :: eta_theta_levels( : )
  REAL, POINTER :: eta_rho_levels( : )
  REAL, POINTER :: rhcrit( : )
  REAL, POINTER :: soil_depths( : )

  ! VarRes horizontal grid spacing
  REAL, POINTER ::  Lambda_p( : )
  REAL, POINTER ::  Lambda_u( : )
  REAL, POINTER ::  Phi_p( : )
  REAL, POINTER ::  Phi_v( : )

  ! Rotated pole information
  REAL :: lambda_npole
  REAL :: phi_npole

END TYPE grid_type

TYPE (grid_type), SAVE :: Input_grid
TYPE (grid_type), SAVE :: Output_grid
END MODULE Rcf_Grid_Type_Mod

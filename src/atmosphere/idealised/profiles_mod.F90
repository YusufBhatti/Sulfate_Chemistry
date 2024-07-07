! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE profiles_mod

! Purpose:
!   Contains profiles and associated variables for idealised configurations.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

USE missing_data_mod, ONLY: rmdi, imdi

IMPLICIT NONE

SAVE

! Parameters
! ======================
! Types of vertical coordinate
INTEGER, PARAMETER :: eta_coord     = 1
INTEGER, PARAMETER :: z_coord       = 2
INTEGER, PARAMETER :: p_coord       = 3

! Types of fields
INTEGER, PARAMETER :: theta_dry     = 10
INTEGER, PARAMETER :: temperature   = 11
INTEGER, PARAMETER :: dtheta_dz     = 12
INTEGER, PARAMETER :: Brunt_Vaisala = 13
INTEGER, PARAMETER :: mixing_ratio  = 20
INTEGER, PARAMETER :: rel_humidity  = 21

! Maximum size of arrays that can be read through namelist
INTEGER, PARAMETER :: num_data_max  = 100

! Profiles
! ========

! Derived type for carrying time-varying profile data
TYPE varying_profile
  INTEGER           :: n_times
  INTEGER           :: n_heights
  INTEGER           :: coord_type
  INTEGER           :: field_type
  REAL              :: timescale
  REAL, ALLOCATABLE :: tsec(:)
  REAL, ALLOCATABLE :: height(:)
  REAL, ALLOCATABLE :: vprof(:,:)
END TYPE varying_profile

! Initial Data Profiles
! =====================
TYPE(varying_profile) :: theta_init ! dry potential temperature
TYPE(varying_profile) :: mv_init    ! vapour mixing ratio
TYPE(varying_profile) :: u_init     ! xi1-component of wind
TYPE(varying_profile) :: v_init     ! xi2-component of wind

! Forcing Profiles
! ================

! Geostrophic forcing profiles
TYPE(varying_profile) :: u_geostrophic
TYPE(varying_profile) :: v_geostrophic

! Newton relaxation profiles
TYPE(varying_profile) :: theta_relax ! dry potential temperature
TYPE(varying_profile) :: mv_relax    ! vapour mixing ratio
TYPE(varying_profile) :: u_relax     ! xi1-component of wind
TYPE(varying_profile) :: v_relax     ! xi2-component of wind


! Increment profiles
TYPE(varying_profile) :: theta_inc
TYPE(varying_profile) :: mv_inc
TYPE(varying_profile) :: u_inc
TYPE(varying_profile) :: v_inc

! Subsidence
TYPE(varying_profile) :: w_subs

! Data Read from Namelist
! =======================

! Relaxation
INTEGER :: num_theta_relax_heights = imdi
INTEGER :: num_theta_relax_times   = imdi
INTEGER :: theta_relax_field_type  = imdi
INTEGER :: num_mv_relax_heights    = imdi
INTEGER :: num_mv_relax_times      = imdi
INTEGER :: num_uv_relax_heights    = imdi
INTEGER :: num_uv_relax_times      = imdi

! Timescale (in seconds) on which fields are relaxed
REAL    :: theta_relax_timescale   = rmdi
REAL    :: mv_relax_timescale      = rmdi
REAL    :: uv_relax_timescale      = rmdi

 
REAL    :: theta_relax_height(num_data_max) = rmdi
REAL    :: theta_relax_time(num_data_max)   = rmdi
REAL    :: theta_relax_data(num_data_max)   = rmdi
REAL    :: mv_relax_height(num_data_max)    = rmdi
REAL    :: mv_relax_time(num_data_max)      = rmdi
REAL    :: mv_relax_data(num_data_max)      = rmdi
REAL    :: uv_relax_height(num_data_max)    = rmdi
REAL    :: uv_relax_time(num_data_max)      = rmdi
REAL    :: u_relax_data(num_data_max)       = rmdi
REAL    :: v_relax_data(num_data_max)       = rmdi

! Increment
INTEGER :: num_theta_inc_times   = imdi
INTEGER :: num_theta_inc_heights = imdi
INTEGER :: theta_inc_field_type  = imdi
INTEGER :: num_mv_inc_times      = imdi
INTEGER :: num_mv_inc_heights    = imdi
INTEGER :: num_uv_inc_times      = imdi
INTEGER :: num_uv_inc_heights    = imdi
INTEGER :: num_w_force_times     = imdi
INTEGER :: num_w_force_heights   = imdi

REAL    :: theta_inc_height(num_data_max) = rmdi
REAL    :: theta_inc_time(num_data_max)   = rmdi
REAL    :: theta_inc_data(num_data_max)   = rmdi
REAL    :: mv_inc_height(num_data_max)    = rmdi
REAL    :: mv_inc_time(num_data_max)      = rmdi
REAL    :: mv_inc_data(num_data_max)      = rmdi
REAL    :: uv_inc_height(num_data_max)    = rmdi
REAL    :: uv_inc_time(num_data_max)      = rmdi
REAL    :: u_inc_data(num_data_max)       = rmdi
REAL    :: v_inc_data(num_data_max)       = rmdi
REAL    :: w_force_height(num_data_max)   = rmdi
REAL    :: w_force_time(num_data_max)     = rmdi
REAL    :: w_force_data(num_data_max)     = rmdi

! For time and vertical height varying geostrophic forcing
INTEGER :: num_uv_geo_times = imdi
INTEGER :: num_uv_geo_heights = imdi

REAL    :: uv_geo_height(num_data_max) = rmdi
REAL    :: uv_geo_time(num_data_max)   = rmdi
REAL    :: u_geo_data(num_data_max)    = rmdi
REAL    :: v_geo_data(num_data_max)    = rmdi


! Hydrostatic Reference Pressure
! ==============================

! Flag to indicate that at least one profile is specified against pressure
LOGICAL                   :: l_p_profile = .FALSE.
REAL, ALLOCATABLE, TARGET :: p_prof_theta(:)
REAL, ALLOCATABLE, TARGET :: p_prof_rho(:)

END MODULE profiles_mod

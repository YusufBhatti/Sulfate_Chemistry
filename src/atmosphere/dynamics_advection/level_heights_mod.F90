! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Constants necessary for model level heights in advection and other schemes.

MODULE level_heights_Mod

! Description:
! This module is used to hold levels values set in atmos_init.
!
! Method:
! The height levels arrays are calculated in routine control/top_level/
! set_levels called from atmos_init,
! and used in dynamics_advection and other routines.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

! Arguments

! Heights Arrays
! - eta values of theta levels
REAL, ALLOCATABLE, TARGET  ::  eta_theta_levels(:)
! - eta values of rho levels
REAL, ALLOCATABLE, TARGET  ::  eta_rho_levels  (:)
! - reference height of levels
REAL, ALLOCATABLE, TARGET  ::  z_ref_theta(:)
! - reference height of levels
REAL, ALLOCATABLE, TARGET  ::  z_ref_rho(:)
! - height of theta levels
REAL, ALLOCATABLE, TARGET  ::  r_theta_levels (:,:,:)
! - height of rho levels
REAL, ALLOCATABLE, TARGET  ::  r_rho_levels   (:,:,:)

REAL :: z_top_theta ! top of model as read in from dump header

! Next two variables are targets for use in  multivariate swap_bounds
REAL, ALLOCATABLE, TARGET ::  r_at_u (:,:,:)    ! height at u points on
                                                ! rho levels
REAL, ALLOCATABLE, TARGET ::  r_at_v (:,:,:)    ! height at v points on
                                                ! rho levels

REAL, ALLOCATABLE, TARGET :: r_at_u_w(:,:,:)  ! height at u points on rho levels
REAL, ALLOCATABLE, TARGET :: r_at_v_w(:,:,:)  ! height at v points on rho levels

! Heights of physics layer centres and boundaries above the centre of the planet
REAL, SAVE, ALLOCATABLE :: r_layer_centres(:,:,:)
REAL, SAVE, ALLOCATABLE :: r_layer_boundaries(:,:,:)

! Depth of physics layers
REAL, SAVE, ALLOCATABLE :: d_layer(:,:,:)

END MODULE level_heights_Mod

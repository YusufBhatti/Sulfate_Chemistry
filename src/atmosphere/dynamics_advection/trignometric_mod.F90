! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Constants necessary for Coriolis terms in advection scheme.

MODULE trignometric_Mod

! Description:
! This module is used to hold trignometric values set in atmos_init.
!
! Method:
! The trignometric arrays are calculated in routine control/top_level/
! set_trig called from atmos_init,
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

! Trignometric Arrays
! Most names are self-defining

! Variables that are targets, are so as they are targets in
! multivariate swap_bounds.

REAL, ALLOCATABLE, TARGET ::  cos_theta_latitude (:,:)
REAL, ALLOCATABLE, TARGET ::  sec_theta_latitude (:,:)
REAL, ALLOCATABLE, TARGET ::  FV_cos_theta_latitude (:,:)
REAL, ALLOCATABLE, TARGET ::  FV_sec_theta_latitude (:,:)
REAL, ALLOCATABLE  ::  sin_theta_latitude (:,:)
REAL, ALLOCATABLE  ::  tan_theta_latitude (:,:)
REAL, ALLOCATABLE  ::  sin_v_latitude (:,:)
REAL, ALLOCATABLE  ::  tan_v_latitude (:,:)
REAL, ALLOCATABLE, TARGET ::  cos_v_latitude (:,:)
REAL, ALLOCATABLE, TARGET ::  sec_v_latitude (:,:)
REAL, ALLOCATABLE  ::  cos_theta_longitude (:,:)
REAL, ALLOCATABLE  ::  sin_theta_longitude (:,:)
REAL, ALLOCATABLE  ::  cos_u_longitude (:,:)
REAL, ALLOCATABLE  ::  sin_u_longitude (:,:)
REAL, ALLOCATABLE  ::  true_latitude (:,:)
REAL, ALLOCATABLE  ::  true_longitude (:,:)

END MODULE trignometric_Mod

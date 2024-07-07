! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Holds Idealised diagnostic arrays

MODULE idealised_diag_mod

IMPLICIT NONE
SAVE

! Description:
!   Module containing idealised diagnostics
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: idealised
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 standards
!
! Declarations:

! 3d fields
REAL, ALLOCATABLE  ::        &
  dt_inc_ideal_um(:,:,:)     & ! dT on model levels (K/timestep)
 ,dq_inc_ideal_um(:,:,:)     & ! dq on model levels (kg/kg/timestep)
 ,du_inc_ideal_um(:,:,:)     & ! dU on model levels (m/s/timestep)
 ,dv_inc_ideal_um(:,:,:)     & ! dV on model levels (m/s/timestep)
 ,dtheta_inc_ideal_um(:,:,:)   ! dtheta on model levels (K/timestep)

! 2d fields
REAL, ALLOCATABLE  ::        &
  dcolqdt_ideal_um(:,:)      & ! rate of change of column water vapour (kg/m2/s)
 ,de_cvt_ideal_um(:,:)       & ! Change of col energy due to dT forcing J/m2
 ,de_u2_ideal_um(:,:)        & ! Change of col energy due to du forcing J/m2
 ,de_v2_ideal_um(:,:)          ! Change of col energy due to dv forcing J/m2

! copies of reference profiles
REAL, ALLOCATABLE  ::        &
  diag_theta_ref(:)          & ! theta reference profile
 ,diag_exner_ref(:)          & ! exner reference profile
 ,diag_rho_ref(:)            & ! density reference profile
 ,diag_u_ref(:)              & ! u wind reference profile
 ,diag_v_ref(:)              & ! u wind reference profile
 ,diag_q_ref(:)                ! water vapour reference profile

LOGICAL :: l_stored_ref      ! .true. if reference profiles stored

END MODULE idealised_diag_mod

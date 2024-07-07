! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Diagnostics necessary for AC assimilation scheme.


MODULE ac_diagnostics_Mod

! Description:
! This module is used to pass diagnostics from physics routines
! to the AC assimilation scheme.
!
! Method:
! By allocating and filling the arrays in routines control/top_level/
! diagnostics_adv/_conv/_lscld/_lsrain, and using them in ac_ctl.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

! Arguments

! Diagnostic Arrays

REAL, ALLOCATABLE  ::  lsrr(:)  ! large scale rain rate
REAL, ALLOCATABLE  ::  lssr(:)  ! large scale snow rate
REAL, ALLOCATABLE  ::  cvrr(:)  ! convective rain rate
REAL, ALLOCATABLE  ::  cvsr(:)  ! convective snow rate
REAL, ALLOCATABLE  ::  convcc(:,:)  ! convective cloud cover
REAL, ALLOCATABLE  ::  cf_lsc(:,:)  ! bulk cloud fraction after LSC
REAL, ALLOCATABLE  ::  tinc_cvn(:,:)! Temp incr across conv
REAL, ALLOCATABLE  ::  tinc_ppn(:,:)! Temp incr across LS precip
REAL, ALLOCATABLE  ::  qcl_lsc(:,:) ! cloud liq water after LSC
REAL, ALLOCATABLE  ::  qcl_adv(:,:) ! cloud liq water after advection

END MODULE ac_diagnostics_Mod

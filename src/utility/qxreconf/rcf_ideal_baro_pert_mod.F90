! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE rcf_ideal_baro_pert_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_IDEAL_BARO_PERT_MOD'

CONTAINS

FUNCTION baro_pert(x,y)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE conversions_mod, ONLY: pi

IMPLICIT NONE
!
! Description:
!
! Applies initial perturbation to zonal wind field.
!
! Method:
!
! Baroclinic instability test case (QJRMS 132, 2943--2975).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Idealised
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Arguments:
REAL, INTENT(IN) :: x,y

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BARO_PERT'

REAL, PARAMETER :: a = 6.371229e6
REAL, PARAMETER :: up = 1.0
REAL            :: r, xc,yc,rg
REAL            :: baro_pert

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

r = a / 10.0
xc = pi / 9.0
yc = 2 * pi / 9.0
rg = a * ACOS(SIN(yc)*SIN(y) + COS(yc) * COS(y) * COS(x - xc))

baro_pert = up * EXP(-(rg / r)**2)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION baro_pert
END MODULE rcf_ideal_baro_pert_mod

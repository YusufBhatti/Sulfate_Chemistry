! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE rcf_ideal_deep_baro_pert_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_IDEAL_DEEP_BARO_PERT_MOD'

CONTAINS

FUNCTION deep_baro_pert(x,y,option, dir)

USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE conversions_mod,     ONLY: pi
USE planet_constants_mod, ONLY: a=>planet_radius
USE umPrintMgr

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

! Arguments
INTEGER, INTENT(IN) :: option, dir
REAL,    INTENT(IN) :: x,y

REAL            :: r, xc,yc,r0,q
REAL            :: denom, drdx, drdy, psi0
REAL            :: deep_baro_pert

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEEP_BARO_PERT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

r0 = a / 10.0
xc = pi / 9.0
yc = 2 * pi / 9.0

r = a * ACOS(SIN(yc)*SIN(y) + COS(yc) * COS(y) * COS(x - xc))

SELECT CASE ( option )

CASE ( 1 ) ! JW06 pertubation

  IF ( dir == 1) THEN
    ! u
    deep_baro_pert = EXP(-(r / r0)**2)
  ELSE
    ! v
    deep_baro_pert = 0.0
  END IF

CASE ( 2 )
  r0 = a/6.0
  psi0 = -r0/2.0
  ! Stream function
  q = 0.0
  IF ( r <= r0 ) THEN
    q  = psi0*EXP( -2.0*r**2/(r - r0)**2 )
    q = -4.0*r/(r-r0)**2*(1.0 - r/(r - r0))*q
  END IF
  denom = 1.0/SQRT(1.0-COS(r/a)**2);

  IF (dir == 1) THEN
    ! u
    drdy = -a*(SIN(yc)*COS(y)-COS(yc)*SIN(y)*COS(x-xc))*denom
    deep_baro_pert = -1.0/a*drdy*q
  ELSE
    ! v
    drdx = a*(COS(yc)*COS(y)*SIN(x-xc))*denom;
    deep_baro_pert = 1.0/(a*COS(y))*drdx*q
  END IF

CASE ( 3 )
  r0 = a/6.0
  !     psi0 = -R0/2.0

  psi0 = -8.0*r0/(3.0*SQRT(3.0)*pi)
  ! Stream function
  q = 0.0
  IF ( r <= r0 ) q  = 0.5*pi*r/r0
  ! Cos^4 Stream function
  deep_baro_pert = -4.0*psi0*q/r*COS(q)**3*SIN(q)
  ! Cos^6 Stream function
  !     deep_baro_pert = -6.0*psi0*q/r*COS(q)**5*SIN(q)

  denom = 1.0/SQRT(1.0-COS(r/a)**2);

  IF (dir == 1) THEN
    ! u
    drdy = -a*(SIN(yc)*COS(y)-COS(yc)*SIN(y)*COS(x-xc))*denom
    deep_baro_pert =  -1.0/a*deep_baro_pert*drdy
  ELSE IF (dir == 2) THEN
    ! v
    drdx = a*(COS(yc)*COS(y)*SIN(x-xc))*denom;
    deep_baro_pert =  1.0/(a*COS(y))*deep_baro_pert*drdx
  ELSE
    ! p
    q = 0.5*pi
    IF ( r <= r0 ) q  = 0.5*pi*r/r0
    deep_baro_pert = -psi0*COS(q)**4
  END IF

CASE DEFAULT
  CALL umPrint('Invalid perturbation used in deep_baro_pert', &
      src='deep_baro_pert',stdErrorToo=.TRUE.)
END SELECT


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION deep_baro_pert
END MODULE rcf_ideal_deep_baro_pert_mod

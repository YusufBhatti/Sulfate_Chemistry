! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
MODULE rcf_ideal_baro_eta_conv_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_IDEAL_BARO_ETA_CONV_MOD'

CONTAINS

REAL FUNCTION baro_eta_conv(xi2,z)

USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE planet_constants_mod, ONLY: r, g
USE rcf_ideal_baro_geo_p_mod, ONLY: baro_geo_p
USE rcf_ideal_baro_T_p_mod, ONLY: baro_T_p

USE ereport_mod, ONLY: ereport
IMPLICIT NONE
!
! Description:
!
! Converts height coordinate (z) to normalised pressure coordinate (eta)
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
!
! Function arguments
REAL, INTENT(IN) :: xi2 ! Latitude (radians)
REAL, INTENT(IN) :: z   ! Distance above Earth's surface

! Local constants
INTEGER, PARAMETER :: maxits = 1000

! Local variables
LOGICAL :: conv
INTEGER :: itrn
INTEGER :: err_code
REAL    :: error
REAL    :: eta_n0, eta_n1
REAL    :: f, f_der,dn

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BARO_ETA_CONV'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

conv = .FALSE.
error = 1.0e-14
eta_n0 = 1.0e-7

DO itrn=0,maxits

  f= -g * z + baro_geo_p(xi2,eta_n0)
  f_der = -r / eta_n0 * baro_T_p(xi2,eta_n0)
  eta_n1 = eta_n0 - f / f_der
  dn = ABS(eta_n1 - eta_n0)
  IF (dn < error) THEN
    conv = .TRUE.
    EXIT
  END IF
  eta_n0 = eta_n1
END DO

IF (conv) THEN
  baro_eta_conv = eta_n1
ELSE
  err_code = 1
  CALL ereport("baro_eta_conv", err_code,                               &
               "did not converge" )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION baro_eta_conv
END MODULE rcf_ideal_baro_eta_conv_mod

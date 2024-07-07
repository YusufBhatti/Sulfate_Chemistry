! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_alpha_ramp_mod

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

INTEGER, SAVE :: alpha_relax_type = imdi
INTEGER, SAVE :: alpha_relax_int  = 1

INTEGER, PARAMETER :: alpha_relax_type_constant = 1
INTEGER, PARAMETER :: alpha_relax_type_sqrt     = 2
INTEGER, PARAMETER :: alpha_relax_type_coslaw   = 3
INTEGER, PARAMETER :: alpha_relax_type_step     = 4


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_ALPHA_RAMP_MOD'

CONTAINS

SUBROUTINE eg_alpha_ramp()

USE eg_alpha_mod
USE timestep_mod,    ONLY: timestep_number
USE um_parcore,      ONLY: mype
USE nlstcall_nrun_as_crun_mod, ONLY: l_nrun_as_crun
USE yomhook,         ONLY: lhook, dr_hook
USE parkind1,        ONLY: jprb, jpim
USE conversions_mod, ONLY: pi
USE umPrintMgr
USE ereport_mod,     ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
! Description:
!
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


LOGICAL, SAVE :: first_call = .TRUE.
REAL,    SAVE :: alpha_star,beta_star

CHARACTER(LEN=errormessagelength)            :: Cmessage
CHARACTER(LEN=15)             :: Routine = 'eg_ramp_alpha'
INTEGER                       :: icode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_ALPHA_RAMP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first_call) THEN

  alpha_star = alpha_u
  beta_star  = 1.0-alpha_u
  first_call = .FALSE.

  IF (l_nrun_as_crun) THEN
    alpha_relax_type = alpha_relax_type_constant
    cmessage = 'l_nrun_as_crun is TRUE:' //  &
               'Overriding value of alpha_relax_type (now set to 1=constant)'
    icode = -1
    CALL ereport(RoutineName, icode, cmessage)
  END IF

END IF

SELECT CASE (alpha_relax_type)

CASE (alpha_relax_type_sqrt)

  alpha_u     = beta_star*1.0/DBLE(timestep_number)**.5 + alpha_star

  alpha_v     = alpha_u
  alpha_w     = alpha_u
  alpha_rho   = alpha_u
  alpha_p     = alpha_u

  alpha_changed = .TRUE.

CASE (alpha_relax_type_constant)

  !   nothing to be done here, we simply use the alphas from the namelist

CASE (alpha_relax_type_coslaw)

  IF (timestep_number <  alpha_relax_int) THEN

    alpha_u   = beta_star*COS(pi*DBLE(timestep_number)/(2.0*DBLE(alpha_relax_int)))**2 + alpha_star

  ELSE

    alpha_u   = alpha_star

  END IF

  alpha_v     = alpha_u
  alpha_w     = alpha_u
  alpha_rho   = alpha_u
  alpha_p     = alpha_u
  alpha_changed = .TRUE.

CASE (alpha_relax_type_step)

  IF (timestep_number <= alpha_relax_int) THEN

    alpha_u   = 1.0

  ELSE

    alpha_u   = alpha_star

  END IF

  alpha_v     = alpha_u
  alpha_w     = alpha_u
  alpha_rho   = alpha_u
  alpha_p     = alpha_u

  alpha_changed = .TRUE.

CASE DEFAULT

  icode = 1
  cmessage = 'unknown alpha relax type'
  CALL Ereport(Routine,icode,cmessage)

END SELECT

IF ( alpha_changed ) THEN
  tau_u   = alpha_u
  tau_v   = alpha_v
  tau_w   = alpha_w
  tau_rho = alpha_rho
  tau_p   = alpha_p
END IF


IF ( mype == 0 .AND. PrintStatus > PrStatus_Normal) THEN
  WRITE(umMessage,FMT='(A,E22.15)') 'ALPHA:', alpha_u
  CALL umPrint(umMessage,src='eg_alpha_ramp_mod')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE
END MODULE

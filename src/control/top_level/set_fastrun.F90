! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine Interface:
SUBROUTINE set_fastrun

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE Control_Max_Sizes
USE nlstcall_mod, ONLY: l_fastrun
USE get_env_var_mod, ONLY: get_env_var

IMPLICIT NONE
!
!    Method:  Reset Rose GUI determined namelist item FASTRUN to 
!             allow IAU assimilation to be skipped

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

! System component covered: Control

! Declarations:

! Local scalars:
CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='SET_FASTRUN')
CHARACTER(LEN=80) :: fastrun    ! value of EnvVar FASTRUN

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Initialise string to false.
fastrun   = "false"

! FASTRUN
CALL get_env_var('FASTRUN',fastrun, allow_missing=.TRUE., allow_empty=.TRUE.)
IF (fastrun == 'true') THEN
  l_fastrun=.TRUE.
  IF (printstatus >= prstatus_oper) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A)')"FASTRUN=true"
      CALL umPrint(umMessage,src='set_fastrun')
    END IF
  END IF  ! PrintStatus test
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_fastrun

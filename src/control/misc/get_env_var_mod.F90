! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc.
!
! Description:
! Take the name of an environment variable as input and try to read it using
! the intrinsic subroutine get_environment_variable. All errors are handled
! internally.
!
! Method:
! The routine will abort if the variable:
! * is unset (i.e. missing)
! * is empty (i.e. zero length string)
! unless the appropriate optional arguments are provided, and ALWAYS if the
! variable:
! * returns a value too long to store
! * cannot be read for hardware reasons.
! The optional argument "length" returns:
! * -1 if the variable is undefined
! * an integer in the range [0:] containing the length of the value otherwise.
!
! NOTE: This routine cannot be called before gc_init_thread (etc.) because
!       ereport requires that MPI comms be set up first.

MODULE get_env_var_mod

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER, PRIVATE  :: ModuleName = 'GET_ENV_VAR_MOD'

CONTAINS

SUBROUTINE get_env_var(env_var, contents, allow_missing, allow_empty, length)

USE ereport_mod, ONLY: ereport
USE umprintmgr, ONLY: newline
USE errormessagelength_mod, ONLY: errormessagelength
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Arguments
CHARACTER(LEN=*), INTENT(IN)               :: env_var    ! Env var to be read
CHARACTER(LEN=*), INTENT(OUT)              :: contents   ! Contents of env var
LOGICAL, INTENT(IN), OPTIONAL              :: allow_missing
                                                         ! If .TRUE.,
                                                         ! allow success if env
                                                         ! var is undefined
LOGICAL, INTENT(IN), OPTIONAL              :: allow_empty! If .TRUE.,
                                                         ! allow success if env
                                                         ! var is empty
INTEGER, INTENT(OUT), OPTIONAL             :: length     ! Length of contents

! Local variables
INTEGER   :: local_length        ! Length return from get_environment_variable
LOGICAL   :: local_allow_empty   ! Local value for allow_empty argument
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message

CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'GET_ENV_VAR'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
INTEGER                       :: StatusReturn   ! Status return value from
                                                ! get_environment_variable
INTEGER                       :: icode          ! Error code passed to ereport

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(allow_empty)) THEN
  local_allow_empty = allow_empty
ELSE
  local_allow_empty = .FALSE.
END IF

CALL GET_ENVIRONMENT_VARIABLE(env_var, contents, local_length, StatusReturn)

IF (StatusReturn /= 0) THEN
  icode = 10
  SELECT CASE (StatusReturn)

  CASE (-1)
    cmessage = 'Environment variable ' // TRIM(env_var) // newline //        &
               'is too long to be stored in the variable provided.'
  CASE (1)
    cmessage = 'Environment variable ' // TRIM(env_var) // ' is not set'     &
               // newline //                                                 &
               'and the call to this routine did not permit a missing input.'
    ! Check whether we should allow the code to continue in this case:
    IF (PRESENT(allow_missing)) THEN
      IF (allow_missing) THEN
        icode = -10
        cmessage = 'Environment variable ' // TRIM(env_var) // ' is not set.'
        ! Overwrite the intrinsic's returned length to something more useful:
        local_length = -1
      END IF
    END IF

  CASE (2)
    cmessage = 'Environment variable ' // TRIM(env_var) //                   &
               ' could not be read;' // newline //                           &
               'the processor does not support environment variables.'
  END SELECT

  CALL ereport(RoutineName, icode, cmessage)

END IF

! Warn/abort if environment variable is empty:
IF (local_length == 0) THEN
  IF (local_allow_empty) THEN
    icode = -20
    cmessage = 'Contents of environment variable ' // TRIM(env_var)          &
               // newline // 'have zero length.'
  ELSE
    icode = 20
    cmessage = 'Contents of environment variable ' // TRIM(env_var)          &
               // newline // 'have zero length and the call to this routine' &
               // newline // 'did not permit an empty variable.'
  END IF
  CALL ereport(RoutineName, icode, cmessage)
END IF

IF (PRESENT(length)) THEN
  length = local_length
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE get_env_var
END MODULE get_env_var_mod

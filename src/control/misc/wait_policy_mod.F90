! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc

MODULE wait_policy_mod

!$ USE OMP_LIB
USE, INTRINSIC :: ISO_C_BINDING

USE umPrintMgr,             ONLY: umprint
USE get_env_var_mod,        ONLY: get_env_var
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,               ONLY: jpim, jprb
USE yomhook,                ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, PARAMETER                  :: default_policy    = 0
INTEGER, PARAMETER                  :: active_policy     = 1
INTEGER, PARAMETER                  :: passive_policy    = 2

INTEGER, SAVE                       :: original_policy   = active_policy

!DrHook-related parameters
INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_in     = 0
INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_out    = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='WAIT_POLICY_MOD'

CONTAINS


!------------------------------------------------------------------------------!
! Change the policy of the waiting OMP threads to the following options:
!   * default - as set by OMP_WAIT_POLICY
!   * active
!   * passive
! 
! At the moment this routine works only with a cray compiler 
! ATTENTION: this routine must not be called inside a parallel region.

SUBROUTINE init_wait_policy()

IMPLICIT NONE

INTEGER             :: length              ! length of env var contents
INTEGER             :: i
INTEGER             :: error_code
CHARACTER(LEN=32)   :: c_policy
CHARACTER(LEN=1)    :: ch
CHARACTER(LEN=errormessagelength) :: cmessage  ! used for ereport

LOGICAL             :: check_env

!DrHook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_WAIT_POLICY'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check the value of the OMP_WAIT_POLICY variable if the waiting
! policy can be modified at run time. At the moment this can only be
! done when a Cray compiler is used.
check_env = .FALSE.

!$ check_env = .TRUE.

#if ! defined (CRAY_FORTRAN)
check_env = .FALSE.
#endif

! Default value
original_policy = active_policy

IF (check_env) THEN
  CALL get_env_var('OMP_WAIT_POLICY',c_policy,allow_missing=.TRUE., &
    length=length)
  IF (length >  0) THEN

    ! Lower case string ASCIII 
    DO i = 1, length
      ch = c_policy(i:i) 
      IF (ch >= 'A' .AND. ch <= 'Z') THEN
        c_policy(i:i) = ACHAR(IACHAR(c_policy(i:i))+32)
      END IF
    END DO

    SELECT CASE(c_policy)
    CASE('active')
      original_policy = active_policy
      CALL umPrint('OMP wait policy is active',src='wait_policy_mod')
    CASE('passive')
      original_policy = passive_policy
      CALL umPrint('OMP wait policy is passive',src='wait_policy_mod')
    CASE DEFAULT
      error_code = 10
      cmessage='OMP_WAIT_POLICY value not identified: '//c_policy
      CALL ereport(RoutineName, error_code, cmessage)
    END SELECT

  ELSE
    CALL umPrint('OMP wait policy unknown, probably active', &
      src='wait_policy_mod')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_wait_policy


!------------------------------------------------------------------------------!
! Change the policy of the waiting OMP threads to the following options:
!   * default - as set by OMP_WAIT_POLICY
!   * active
!   * passive
! 
! At the moment this routine works only with a cray compiler 
! ATTENTION: this routine must not be called inside a parallel region.

SUBROUTINE set_wait_policy(p_type) BIND(c,NAME="set_wait_policy")

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN) :: p_type

! Local variables
INTEGER(KIND=C_INT64_T) :: policy_type

INTEGER :: error_code

!DrHook
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_WAIT_POLICY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


IF (p_type == default_policy) THEN
  policy_type = original_policy
ELSE
  policy_type = p_type
END IF

SELECT CASE(policy_type)
  
CASE(passive_policy)

#if defined (CRAY_FORTRAN)
!$ CALL cray_omp_set_wait_policy ( 'PASSIVE' ) 
#endif

CASE(active_policy)

#if defined (CRAY_FORTRAN)
!$ CALL cray_omp_set_wait_policy ( 'ACTIVE' ) 
#endif

CASE DEFAULT
  error_code = 10
  CALL ereport(RoutineName, error_code, 'Unknown wait policy type')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_wait_policy

END MODULE wait_policy_mod

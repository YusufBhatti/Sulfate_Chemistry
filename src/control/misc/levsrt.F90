! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE levsrt_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LEVSRT_MOD'

CONTAINS
!
! Subroutine Interface:
SUBROUTINE levsrt(TYPE,nlevs,il,rl)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
! Description:
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project

! Subroutine arguments:

!   Scalar arguments with intent(in):

CHARACTER(LEN=1) :: TYPE
INTEGER ::  nlevs

!   Array arguments with intent(inout):

REAL ::     rl(nlevs)
INTEGER ::  il(nlevs)

! Local variables:

LOGICAL ::  lswap
INTEGER ::  i
INTEGER ::  j
INTEGER ::  ilt
REAL ::     rlt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LEVSRT'

!- End of Header ----------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO i=1,nlevs
  lswap=.FALSE.
  DO j=1,nlevs-1
    IF (TYPE == 'I') THEN
      IF (il(j) >  il(j+1)) THEN
        lswap=.TRUE.
        ilt=il(j)
        il(j)=il(j+1)
        il(j+1)=ilt
      END IF
    ELSE
      IF (rl(j) <  rl(j+1)) THEN
        lswap=.TRUE.
        rlt=rl(j)
        rl(j)=rl(j+1)
        rl(j+1)=rlt
      END IF
    END IF
  END DO
  IF (.NOT. lswap) THEN
    GO TO 9999
  END IF
END DO

9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE levsrt
END MODULE levsrt_mod

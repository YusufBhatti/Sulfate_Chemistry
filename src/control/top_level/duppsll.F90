! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Determine whether a pseudo lev list is a duplicate of another one
! Subroutine Interface:

SUBROUTINE duppsll(ldupll,ndupll)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE version_mod, ONLY:                                            &
    nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
    nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
    nlevp, npslevp, npslistp
USE stextend_mod, ONLY: lenplst
USE cstash_mod,   ONLY: pslist_d, npslists

IMPLICIT NONE
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
! Global variables:


! Subroutine arguments:

!   Scalar arguments with intent(out):

LOGICAL :: ldupll
LOGICAL :: llocal
INTEGER :: ndupll

! Local scalars:

INTEGER :: i
INTEGER :: j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DUPPSLL'

!- End of Header -------------------------------------------------


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
ldupll=.FALSE.
ndupll=0
DO i=1,npslists-1
  IF (lenplst(i) == lenplst(npslists)) THEN
    llocal=.TRUE.
    DO j=1,lenplst(npslists)
      IF (pslist_d(j,i) /= pslist_d(j,npslists)) THEN
        llocal=.FALSE.
        EXIT
      END IF
    END DO

    IF (llocal) THEN
      ldupll=.TRUE.
      ndupll=i
      GO TO 9999
    END IF
  END IF
END DO

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE duppsll

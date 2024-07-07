! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Compress out unused pseudo levels lists
SUBROUTINE pslcom(nrecs)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE version_mod, ONLY:                                            &
    nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
    nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
    nlevp, npslevp, npslistp
USE stextend_mod, ONLY: list_s, lenplst
USE stparam_mod, ONLY: st_pseudo_out
USE cstash_mod,   ONLY: pslist_d, npslists

IMPLICIT NONE

! Subroutine arguments:

!   Scalar arguments with intent(in):

INTEGER :: nrecs

! Local Scalars:

INTEGER :: icount

! Local arrays:

INTEGER :: ipos (npslistp)  ! POSITION IN OLD LIST OF THE NEW
INTEGER :: ipos1(npslistp)  ! POSITION IN THE NEW LIST OF THE OLD

INTEGER :: i,j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PSLCOM'

!- End of Header ---------------------------------------------------


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
icount=0

DO i=1,npslists         ! LOOP DOWN THE LISTS
  IF (lenplst(i) /= 0) THEN   ! USED LIST
    icount=icount+1
    ipos(icount)=i
    ipos1(i)=icount
  ELSE
    ipos1(i)=0
  END IF
END DO

npslists=icount

DO    i=1,npslists
  lenplst(i)=lenplst(ipos(i))
  DO  j=1,npslevp
    pslist_d(j,i)=pslist_d(j,ipos(i))
  END DO
END DO

DO   i=npslists+1,npslistp
  lenplst(i)=0
  DO j=1,npslevp
    pslist_d(j,i)=0
  END DO
END DO

DO i=1,nrecs
  IF (list_s(st_pseudo_out,i) /= 0) THEN
    list_s(st_pseudo_out,i)=ipos1(list_s(st_pseudo_out,i))
  END IF
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pslcom

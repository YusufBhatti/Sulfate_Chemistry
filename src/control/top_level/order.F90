! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Bubble sort of STASH list by model/section/item/input code
!
! Subroutine Interface:

SUBROUTINE order(nrecs)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE version_mod, ONLY:                                            &
    nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
    nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
    nlevp, npslevp, npslistp
USE stextend_mod, ONLY: list_s
USE stparam_mod, ONLY: st_model_code, st_sect_no_code, st_item_code,&
                       st_input_code

IMPLICIT NONE

!  Description:
!    Bubble sort of STASH list by model/section/item/input code
!
!  Method:
!  Loops through records in LIST_S (STASH list) and compares values of
!  model number, then section, then item, then input code, then LIST_S
!  elements in order. If the values are not in descending order for any
!  one of these categories, the records are interchanged. A logical
!  LSWAP is used to determine when the sort is complete.
!    Called by DUPLIC, STPROC
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  Global variables
!
!  Subroutine arguments
!    Scalar argument with intent(in):
INTEGER :: nrecs       !  No. of records in STASH list

!  Local scalars:

INTEGER :: list_t !  Used for swapping values in consecutive records
LOGICAL :: lswap  !  If TRUE, order not complete
INTEGER :: i
INTEGER :: j
INTEGER :: irec
INTEGER :: irec1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ORDER'

!- End of Header ------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
lswap=.TRUE.

100 CONTINUE            ! return point for order loop

lswap=.FALSE.

!  Loop over records in LIST_S

DO irec=nrecs-1,1,-1
  irec1=irec+1

  !  Sort by internal model

  IF (list_s(st_model_code,irec) >                               &
     list_s(st_model_code,irec1)) THEN
    lswap=.TRUE.
    DO i=1,nelemp+1
      list_t=list_s(i,irec)
      list_s(i,irec)=list_s(i,irec1)
      list_s(i,irec1)=list_t
    END DO
  END IF

  !  Sort by section for same internal model

  IF ( (list_s(st_model_code,irec) ==                            &
       list_s(st_model_code,irec1)) .AND.                        &
      (list_s(st_sect_no_code,irec) >                           &
       list_s(st_sect_no_code,irec1)) ) THEN
    lswap=.TRUE.
    DO i=1,nelemp+1
      list_t=list_s(i,irec)
      list_s(i,irec)=list_s(i,irec1)
      list_s(i,irec1)=list_t
    END DO
  END IF

  !  Sort by item for same model, section

  IF ( (list_s(st_model_code,irec) ==                            &
       list_s(st_model_code,irec1)) .AND.                        &
      (list_s(st_sect_no_code,irec) ==                          &
       list_s(st_sect_no_code,irec1)) .AND.                      &
      (list_s(st_item_code,irec) >                              &
       list_s(st_item_code,irec1)) ) THEN
    lswap=.TRUE.
    DO i=1,nelemp+1
      list_t=list_s(i,irec)
      list_s(i,irec)=list_s(i,irec1)
      list_s(i,irec1)=list_t
    END DO
  END IF

  !  Sort by input code when model, section, item are correct

  IF ( (list_s(st_input_code,irec) <                             &
       list_s(st_input_code,irec1)) .AND.                        &
      (list_s(st_model_code,irec) ==                            &
       list_s(st_model_code,irec1)) .AND.                        &
      (list_s(st_item_code,irec) ==                             &
       list_s(st_item_code,irec1)) .AND.                         &
      (list_s(st_sect_no_code,irec) ==                          &
       list_s(st_sect_no_code,irec1)) ) THEN
    lswap=.TRUE.
    DO i=1,nelemp+1
      list_t=list_s(i,irec)
      list_s(i,irec)=list_s(i,irec1)
      list_s(i,irec1)=list_t
    END DO
  END IF

  !  Sort by all elements when model, section, item, input code are correct

  IF ( (list_s(st_input_code,irec) ==                            &
       list_s(st_input_code,irec1)) .AND.                        &
      (list_s(st_model_code,irec) ==                            &
       list_s(st_model_code,irec1)) .AND.                        &
      (list_s(st_item_code,irec) ==                             &
       list_s(st_item_code,irec1)) .AND.                         &
      (list_s(st_sect_no_code,irec) ==                          &
       list_s(st_sect_no_code,irec1)) ) THEN
    DO i=1,nelemp+1
      IF (list_s(i,irec) > list_s(i,irec1)) THEN
        lswap=.TRUE.
        DO j=1,nelemp+1
          list_t=list_s(j,irec)
          list_s(j,irec)=list_s(j,irec1)
          list_s(j,irec1)=list_t
        END DO
        EXIT
      ELSE IF (list_s(i,irec) < list_s(i,irec1)) THEN
        EXIT  ! Stop the sort here.
      END IF
    END DO
  END IF

END DO     !  Loop over records

IF (lswap) GO TO 100

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE order

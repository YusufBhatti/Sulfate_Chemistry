! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Update the dump array of a SCM diagnostic entry in SCMop

SUBROUTINE add2dump (x, nelements, d, SCMop, startperiod, endperiod)

! SCMop_type is defined in here...
USE scmoptype_defn
USE s_scmop_mod, ONLY:                                                     &
    t_inst, t_avg, t_max, t_min, t_acc, t_div, t_mult, t_acc_div           &
  , t_acc_mult, t_const

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE


! Description:
!   Take input array x and add to the dump array in the manner
!   specified by the time profile. This performs this timestep's
!   role in constructing the diagnostic.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran90

TYPE(SCMop_type) :: SCMop ! InOut The derived-type structure
                          !       containing all the diagnostic
                          !       information

INTEGER :: nelements      ! The no. of elements in array x

REAL :: x(nelements)      ! In The array from which the diagnostic
                          !    is constructed

INTEGER :: d              ! In The diagnostic index in SCMop

LOGICAL ::               &
  startperiod            &! In These flag the start and end
, endperiod               !    of the dumping period, when special
                          !    actions may need to be taken


REAL, POINTER :: dump(:)  ! The dump array in which the
                          ! diagnostic is constructed

INTEGER :: timprof        ! The time profile of the diagnostic

INTEGER :: nadd2dump      ! The number by which to divide
                          ! accumulated diag's in order to get the mean
INTEGER :: wnelements     ! The no. of elements in array wdump

REAL, POINTER :: wdump(:) ! A dump array of another
                          ! diagnostic which may be used in the
                          ! construction of this diagnostic

INTEGER :: wtimprof       ! The timeprofile of the other diag.
INTEGER :: wnadd2dump     ! As nadd2dump, but for the other diagnostic

INTEGER :: i              ! Counter

! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ADD2DUMP'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)


! We don't want lots of '%'s around, make some abbreviations...
dump      => SCMop%diag(d)%dump
nelements =  SCMop%nelements(d)
timprof   =  SCMop%timprof(d)
nadd2dump =  SCMop%nadd2dump(d)

IF (SCMop%wd(d) >  0) THEN
  wdump      => SCMop%diag(SCMop%wd(d))%dump
  wnelements =  SCMop%nelements(SCMop%wd(d))
  wtimprof   =  SCMop%timprof(SCMop%wd(d))
  wnadd2dump =  SCMop%nadd2dump(SCMop%wd(d))
ELSE
  NULLIFY(wdump)
  wnelements = 0
  wtimprof   = 0
  wnadd2dump = 0
END IF

!-T_CONST--A constant value---------------------------------------------
IF (timprof == t_const) THEN
  DO i=1, nelements
    dump(i) = x(i)
  END DO

  !-T_INST--Instaneous value----------------------------------------------
ELSE IF (timprof == t_inst) THEN
  IF (startperiod) THEN

    ! Set the dump to a default value. Shouldn't really be
    ! necessary, but handy for checking everything's working OK.
    DO i=1, nelements
      dump(i) = -999.0
    END DO
  END IF

  IF (endperiod) THEN

    ! The dump will ocurr at this timestep, so we want the value now.
    DO i=1, nelements
      dump(i) = x(i)
    END DO
  END IF

  !-T_ACC--Accumulated value----------------------------------------------
ELSE IF (timprof == t_acc) THEN
  IF (startperiod) THEN

    ! Take this timestep's values
    DO i=1, nelements
      dump(i) = x(i)
    END DO
  ELSE

    ! Add to the current values
    DO i=1, nelements
      dump(i) = dump(i) + x(i)
    END DO
  END IF

  !-T_AVG--Value averaged over a period----------------------------------
ELSE IF (timprof == t_avg) THEN
  IF (startperiod) THEN

    ! Take this timestep's values
    DO i=1, nelements
      dump(i) = x(i)
    END DO
  ELSE

    ! Add to the current values
    DO i=1, nelements
      dump(i) = dump(i) + x(i)
    END DO
  END IF

  ! And form average if it's the end of the period
  IF (endperiod) THEN
    DO i=1, nelements
      dump(i) = dump(i)/nadd2dump
    END DO
  END IF

  !-T_MIN--The minimum value over a time period---------------------------
ELSE IF (timprof == t_min) THEN
  IF (startperiod) THEN

    ! Take this timestep's values
    DO i=1, nelements
      dump(i) = x(i)
    END DO
  ELSE

    ! Overwrite the current values if smaller
    DO i=1, nelements
      IF (x(i) <  dump(i)) THEN
        dump(i) = x(i)
      END IF
    END DO
  END IF

  !-T_MAX--The maximum value over a time period---------------------------
ELSE IF (timprof == t_max) THEN
  IF (startperiod) THEN

    ! Take this timestep's values
    DO i=1, nelements
      dump(i) = x(i)
    END DO
  ELSE

    ! Overwrite the current value if bigger
    DO i=1, nelements
      IF (x(i) >  dump(i)) THEN
        dump(i) = x(i)
      END IF
    END DO
  END IF

  !-T_DIV-T_MULT--The accumulated value of one diagnostic multiplied or---
  !---------------divided by the value of another-------------------------
ELSE IF (timprof == t_div    .OR. timprof == t_mult .OR.               &
         timprof == t_acc_div .OR. timprof == t_acc_mult) THEN

  IF (startperiod) THEN
    DO i=1, nelements
      dump(i) = x(i)
    END DO
  ELSE

    ! Add to the current value
    DO i=1, nelements
      dump(i) = dump(i) + x(i)
    END DO
  END IF

  IF (endperiod) THEN
    IF (timprof == t_mult .OR. timprof == t_acc_mult) THEN

      ! Multiply the two diagnostics
      DO i=1, nelements
        dump(i) = dump(i)*wdump(MOD(i-1,wnelements)+1)
      END DO
    ELSE IF (timprof == t_div .OR. timprof == t_acc_div) THEN

      ! Divide the two diagnostics
      DO i=1, nelements
        IF (wdump(MOD(i-1,wnelements)+1) /= 0) THEN
          dump(i) = dump(i)/                                          &
                    wdump(MOD(i-1,wnelements)+1)
        ELSE IF (dump(i) /= 0) THEN
          WRITE(umMessage,'(A,1X,A,1X,I3,1X,A,1X,I6)')               &
            'Add2Dump WARNING:divide by zero avoided:'
          CALL umPrint(umMessage,src='add2dump')
        END IF
      END DO
    END IF

    ! Divide by dumping period?
    IF (timprof == t_div .OR. timprof == t_mult) THEN
      DO i=1, nelements
        dump(i) = dump(i)/nadd2dump
      END DO
    END IF
  END IF

  !-----------------------------------------------------------------------
ELSE
  WRITE(umMessage,'(A,I2)')                                           &
    'Add2Dump ERROR: Unknown temporal type: ', timprof
  CALL umPrint(umMessage,src='add2dump')
END IF
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE add2dump


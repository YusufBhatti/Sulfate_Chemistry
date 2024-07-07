! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE scmoutput_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SCMOUTPUT_MOD'

CONTAINS


!  Create a SCM output diagnostic

SUBROUTINE scmoutput                                                          &
  ( x, sname, lname, units, timprof, domprof, streams, sname2                 &
  , calling_routine )

USE global_scmop
USE scm_utils, ONLY: scm_trap_nan
USE s_scmop_mod, ONLY: stepnmbr, t_const

USE umPrintMgr,             ONLY: newline
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim


IMPLICIT NONE


! Description:
!   Create output diagnostic based on given inputs. This routine
!   cannot be called more than once per timestep with the same
!   "sname" unless the call is inside a sub-stepped part of the
!   model and the sub-stepping has been delimited by calls to
!   scm_substep_start and scm_substepping_end. The order of calls
!   to SCMoutput should not change between timesteps.

! Method:

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran77 with some bits of Fortran90

INTEGER :: domprof  ! In Domain profile for the diagnostic
REAL :: x  &        ! In Variable from which diagnostic will be constructed
    ( SCMop%d_rowa1(domprof):SCMop%d_rowa2(domprof),           &
      SCMop%d_rowb1(domprof):SCMop%d_rowb2(domprof),           &
      SCMop%d_lev1 (domprof):SCMop%d_lev2 (domprof) )

CHARACTER(LEN=*) :: &
  sname             &! In Short name for the diagnostic,
                     !    this should be unique
, lname             &! In Long name for the diagnostic
, units             &! In Units of the diagnostic
, sname2             ! In Short name of another, previously
                     !    defined diagnostic which will be used
                     !    in the construction of this one
                     !    according to the time profile

CHARACTER(LEN=*) :: &
  calling_routine    ! In Routine that has called scmoutput

INTEGER ::          &
  timprof           &! In The time profile for the diagnostic
, streams            ! In An encoded integer specifying
                     !    which output streams the diagnostic
                     !    is to go to

INTEGER :: d,n,i,j,k     ! General use

CHARACTER(LEN=30) :: sdum0

LOGICAL ::          &
  startperiod       &! Will be used to flag if we are at
, endperiod          ! the start or end of a dumping period

INTEGER ::                     &
  ndistinct_diags              &
, distinct_diags(maxnstreams)  &
, nThroughPeriod               &
, ntrad                        &
, ntrad1

LOGICAL ::       &
  call_add2dump  &
, order_changing

! Will hold the contents of SCMop%diag_mem while being
! re-allocated
INTEGER, ALLOCATABLE :: itemp(:,:)

! Flag for found NaN in the diagnostic field.
LOGICAL :: l_found_nan

! Inputs for ereport
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage


CHARACTER(LEN=*), PARAMETER :: RoutineName='SCMOUTPUT'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------


! Perform no action if SCMop is not turned on, except in the
! case of a constant diagnostic
IF (.NOT. SCMop%on .AND. timprof /= t_const) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Enforce a rule that sname must be at least one character long
IF (LEN(sname) == 0) THEN
  IF (SCMop%first_pass) THEN
    ! Raise a warning that the diagnostic will be ignored
    WRITE(cmessage,'(A)')                                                      &
      '=========================================================='//  newline//&
      '| ERROR scmoutput:'//                                          newline//&
      '| The variable sname, is a null string, this is not allowed.'//newline//&
      '| Diagnostic ignored:'//                                       newline//&
      '| '//TRIM(ADJUSTL(lname))//                                    newline//&
      '=========================================================='
    icode = -789
    CALL ereport(routinename, icode, cmessage)
  END IF
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Increment recorded no. of calls to this routine this timestep
IF (SCMop%on) SCMop%nSCMoutput = SCMop%nSCMoutput+1

! If the requested list of streams span a range of dumping
! periods then the given inputs will correspond to more than one
! diagnostic entry in SCMop. Thus we need to get
! distinct_diags(1:ndistinct_diags) - the indices of those
! entries. The routine getdistinctdiags can do this but, if this
! is not the first timestep, and assuming the order of the calls
! to this routine doesn't change between timesteps, we can just
! use our memory of a previous timestep instead.

order_changing = .FALSE.
IF (.NOT. SCMop%first_pass) THEN
  ! Use memory of previous timesteps to know which diagnostic
  ! entries this call pertains to.
  ndistinct_diags = SCMop%diag_mem(SCMop%nSCMoutput,1)
  DO n=1, ndistinct_diags
    distinct_diags(n) = SCMop%diag_mem(SCMop%nSCMoutput,2)+n-1

    ! Check we're right...
    IF (TRIM(sname(1:MIN(lsname,LEN(sname)))) /=                               &
        TRIM(SCMop%sname(distinct_diags(n)))) THEN
      WRITE(cmessage,'(A,I4,A)')                                               &
      '=========================================================='//  newline//&
      '| Warning scmoutput:'//                                        newline//&
      '| The order of the calls to scmoutput seems to be changing.'// newline//&
      '| On step ',stepnmbr(SCMop),' expected '                              //&
      TRIM(SCMop%sname(distinct_diags(n))) //                         newline//&
      '| but encountered ' // sname(1:MIN(lsname,LEN(sname))) //      newline//&
      '=========================================================='
      icode = -678
      CALL ereport(routinename, icode, cmessage)

      order_changing = .TRUE.

    ELSE IF (SCMop%substep_number /=                                           &
             SCMop%substep(distinct_diags(n))) THEN
      WRITE(cmessage,'(A,I5,A,I3,A,I3,A)')                                     &
      '=========================================================='//  newline//&
      '| ERROR scmoutput:'//                                          newline//&
      '| The order of the sub-steps seems to be changing.'//          newline//&
      '| On step ',stepnmbr(SCMop),' diagnostic '                            //&
      sname(1:MIN(lsname,LEN(sname))) //                              newline//&
      '| expected sub-step ', SCMop%substep(distinct_diags(n)),                &
      ' but got sub-step ',SCMop%substep_number,                      newline//&
      '| This does not make sense and is a sign of a potentially'//   newline//&
      '| serious problem.'//                                          newline//&
      '=========================================================='
      icode = 678
      CALL ereport(routinename, icode, cmessage)
      order_changing = .TRUE.
    END IF
  END DO
END IF

IF (SCMop%first_pass .OR. order_changing) THEN

  ! Either this is the first timestep or the order of the calls
  ! to SCMoutput is changing

  ! DEPENDS ON: getdistinctdiags
  CALL getdistinctdiags                                       &
    ( sname, lname, units, timprof, domprof, streams, sname2  &! In
    , ndistinct_diags, distinct_diags, SCMop )      ! Out,Out,InOut

  ! Don't want to do this next bit for constant diagnostics
  ! (which should be declared when the system is off)
  IF (SCMop%on) THEN

    ! Store the values of ndistinct_diags and distinct_diags(1)
    ! for future reference (to avoid unnecessary calls to
    ! getdistinctdiags). But is there enough space in the
    ! SCMop%diag_mem array?
    IF (SIZE(SCMop%diag_mem,1) == SCMop%nSCMoutput-1) THEN
      ! No. Make the array bigger...
      ! Allocate a temporary array with the same size and
      ! shape as diag_mem
      ALLOCATE(itemp(SCMop%nSCMoutput-1,2))

      ! Copy the contents of diag_mem into it
      itemp = SCMop%diag_mem

      ! Re-allocate diag_mem with a larger size
      DEALLOCATE(SCMop%diag_mem)
      ALLOCATE(SCMop%diag_mem(SCMop%nSCMoutput+49,2))

      ! Copy the original contents back in
      SCMop%diag_mem(1:SCMop%nSCMoutput-1,1:2) = itemp
    END IF

    SCMop%diag_mem(SCMop%nSCMoutput,1) = ndistinct_diags
    SCMop%diag_mem(SCMop%nSCMoutput,2) = distinct_diags(1)
  END IF
END IF

! From here on in, none of the input parameters are
! referred to at all, they have been distilled to
! distinct_diags(1:ndistinct_diags) and information
! in SCMop

ntrad  = SCMop%ntrad       ! No. of timesteps between calls to rad'n
ntrad1 = SCMop%ntrad1      ! Timestep containing 1st call to rad'n

DO n=1, ndistinct_diags
  d = distinct_diags(n)

  ! If this is not the first time we've seen this diagnostic,
  ! check the last time was the previous timestep
  IF (SCMop%lastencounter(d) >= 0 .AND.                                 &
      SCMop%lastencounter(d) /= stepnmbr(SCMop)-1) THEN

    WRITE(cmessage,'(A,5I5,A)')                                                &
      '=========================================================='//  newline//&
      '| ERROR scmoutput:'//                                          newline//&
      '| The last encounter with this diagnostic was not the last'//  newline//&
      '| timestep:'//                                                 newline//&
      '| ' // TRIM(ADJUSTL(sname)), SCMop%lastencounter(d),                    &
      stepnmbr(SCMop), (SCMop%daycount-1), SCMop%full_daysteps,                &
      SCMop%stepcount,                                                newline//&
      '=========================================================='
    icode = 987
    CALL ereport(routinename, icode, cmessage)
  END IF

  ! Record the fact that this diagnostic was seen by
  ! this routine on this timestep.
  SCMop%lastencounter(d) = stepnmbr(SCMop)

  ! Calculate how many timesteps we are through the current
  ! dumping period. 1=first time step, dump_step=last timestep
  ! (and so a dump will occur this timestep).
  nThroughPeriod = MOD(stepnmbr(SCMop)-1,SCMop%dump_step(d))+1

  ! Decide whether we are at the start of a dumping period, at
  ! the end of of a dumping period, and whether we need to call
  ! add2dump (using nThroughPeriod this is trivial for most
  ! diagnostics, but has added complications in the case of
  ! diagnostics only calculated on radiation timesteps)
  startperiod   = .FALSE.
  endperiod     = .FALSE.
  call_add2dump = .TRUE.

  IF (.NOT. SCMop%only_radsteps(d)) THEN
    ! Diagnostic d is a normal diagnostic, valid at every timestep
    IF (nThroughPeriod == 1) startperiod = .TRUE.
    IF (nThroughPeriod == SCMop%dump_step(d)) endperiod = .TRUE.

  ELSE
    ! Diagnostic d is based on a variable which only has
    ! valid values on radiation timesteps.

    IF (MOD(SCMop%stepcount-ntrad1,ntrad) /= 0) THEN

      ! This is not a radiation timestep, assume input
      ! array x contains nonsense information - do not
      ! call add2dump.
      call_add2dump = .FALSE.
    ELSE

      ! The criteria for startperiod and endperiod are now
      ! altered slightly, since startperiod must be true if
      ! this is the first radiation time step during this
      ! dumping period, and endperiod must be true if this is
      ! the last radiation time step during this dumping
      ! period.
      IF (nThroughPeriod-1 <  ntrad) startperiod = .TRUE.
      IF (SCMop%dump_step(d)-nThroughPeriod <  ntrad) endperiod = .TRUE.

      ! This is a radiation timestep and so we can call
      ! add2dump, but the no. of timesteps by which to divide
      ! in order to calculate the average (or whatever) is
      ! not dump_step, but the no. of times this part of the
      ! code has been reached during this dumping period,
      ! stored in SCMop%nadd2dump(d) (which for normal
      ! diagnostics is set to dump_step in newdiag)

      IF (startperiod) THEN
        SCMop%nadd2dump(d) = 1
      ELSE
        SCMop%nadd2dump(d) = SCMop%nadd2dump(d)+1
      END IF

    END IF            ! (mod(SCMop%stepcount-ntrad1,ntrad) /= 0)
  END IF             ! (.not.SCMop%only_radsteps(d))

  IF (call_add2dump) THEN

    ! Add this diag to the diagnostic dump
    ! (pass in first element of the 3-D array x, to use sequential memory
    ! access to map onto 1-D array dummy argument in add2dump).
    ! DEPENDS ON: add2dump
    CALL add2dump( x( SCMop%d_lev1(domprof),                           &
                      SCMop%d_rowb1(domprof),                          &
                      SCMop%d_rowa1(domprof) ),                        &
                   SCMop%nelements(d), d, SCMop, startperiod, endperiod )

    ! Check for NaNs in the data
    l_found_nan = .FALSE.
    DO k = SCMop%d_lev1(domprof), SCMop%d_lev2(domprof)
      DO j = SCMop%d_rowb1(domprof), SCMop%d_rowb2(domprof)
        DO i = SCMop%d_rowa1(domprof), SCMop%d_rowa2(domprof)
          IF (.NOT. l_found_nan ) THEN

            WRITE(sdum0,"(ES30.10)") x(i,j,k)

            IF ( INDEX(sdum0, '*'  ) /= 0 .OR.                        &
                 INDEX(sdum0, 'NaN') /= 0 .OR.                        &
                 INDEX(sdum0, 'nan') /= 0 .OR.                        &
                 INDEX(sdum0, 'NAN') /= 0      ) THEN
              l_found_nan = .TRUE.
              CALL scm_trap_nan (TRIM(ADJUSTL(sname)), calling_routine)
            END IF

          END IF
        END DO
      END DO
    END DO

  END IF

END DO

!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!----------------------------------------------------------------------
RETURN
END SUBROUTINE scmoutput


END MODULE scmoutput_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Create a new SCM diagnostic entry in SCMop.

SUBROUTINE newdiag                                                            &
  ( sname, lname, units, timprof, domprof, istrm, lnot_written, sname2, d     &
  , SCMop )

! SCMop_type is defined in here...
USE scmoptype_defn
USE s_scmop_mod, ONLY: stream, DoNotWrite                                   &
  , t_inst, t_avg, t_max, t_min, t_acc, t_div, t_mult, t_acc_div            &
  , t_acc_mult, t_const, only_radsteps

USE umPrintMgr, ONLY: umPrint, umMessage, newline
USE yomhook,    ONLY: lhook, dr_hook
USE parkind1,   ONLY: jprb, jpim

IMPLICIT NONE


! Description:
!   Create a new diagnostic entry in SCMop. Returns index of
!   newly created diagnostic, which is unchanged from its input
!   value if an error ocurrs and the entry could not be created.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran 90

TYPE(SCMop_type) :: SCMop       ! InOut The derived-type structure
                                !       containing all the diagnostic
                                !       information
CHARACTER (LEN=lsname) :: sname ! In Short name for the diagnostic
CHARACTER (LEN=*) :: lname      ! In Long name for the diagnostic
CHARACTER (LEN=*) :: units      ! In Units of the diagnostic
CHARACTER (LEN=*) :: sname2     ! In Short name of a previously
                                !    defined diagnostic used in the
                                !    construction of this one
INTEGER ::  &
  timprof   &! In Time profile for the diagnostic
, domprof   &! In Domain profile for the diagnostic
, istrm      ! In Stream to which diagnostic is to be sent

LOGICAL ::  &
  lnot_written ! In If true, diagnostic will not be written out

INTEGER ::  &
  d            ! InOut Index of the newly created
               !       diagnostic entry in SCMop. Unchanged from
               !       input value if entry could not be created.

INTEGER :: d_input,i,d1,d2,d3

LOGICAL ::                     &
  sname2_found                 &
, right_dumping_period

! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NEWDIAG'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

d_input = d ! record input value of d

IF (SCMop%nentries == SCMop%maxnentries) THEN

  ! No room at the inn, make the inn bigger.
  ! DEPENDS ON: expand_scmop
  CALL expand_scmop(SCMop)
END IF

! Increment the recorded number of entries
SCMop%nentries = SCMop%nentries+1

! d is the index of this new entry
d = SCMop%nentries

! Record the details of this entry in the respective arrays...
SCMop%sname(d)   = sname
SCMop%lname(d)   = lname
SCMop%units(d)   = units
SCMop%domprof(d) = domprof

IF (timprof <  only_radsteps) THEN

  ! This is a normal diagnostic - calculated on every timestep
  SCMop%timprof(d) = timprof
  SCMop%only_radsteps(d) = .FALSE.
ELSE

  ! This diagnostic is based on an array only valid on
  ! radiation timesteps
  SCMop%timprof(d) = timprof-only_radsteps
  SCMop%only_radsteps(d) = .TRUE.
END IF

SCMop%streams(d) = Stream(istrm)
IF (lnot_written) THEN
  SCMop%streams(d) = DoNotWrite(SCMop%streams(d))
END IF

! The dumping period of the diagnostic is that of the stream
! it's being sent to
SCMop%dump_step(d) = SCMop%strm(istrm)%dump_step
SCMop%nadd2dump(d) = SCMop%strm(istrm)%dump_step

! The dimensions of the diagnostic array can be obtained from
! the domain profile
d1 = SCMop%d_rowa2(domprof)-SCMop%d_rowa1(domprof)+1
d2 = SCMop%d_rowb2(domprof)-SCMop%d_rowb1(domprof)+1
d3 = SCMop%d_lev2 (domprof)-SCMop%d_lev1 (domprof)+1
SCMop%ncols(d) = d1
SCMop%nrows(d) = d2
SCMop%nlevs(d) = d3
SCMop%nelements(d) = d1*d2*d3

! Allocate the space for the dump array
ALLOCATE(SCMop%diag(d)%dump(d1*d2*d3))
! Initialise it for initialisation's sake
SCMop%diag(d)%dump = -999.0

! lastencounter will be set in SCMoutput
SCMop%lastencounter(d) = -1

! Set the substep number that this entry is being created for
SCMop%substep(d) = SCMop%substep_number

! wd will be the index of a diagnostic entry upon which this
! diagnostic depends (set to zero if sname2 is a null string)
SCMop%wd(d) = 0

IF (LEN(sname2) >  0) THEN
  ! Find the index of the weighting diagnostic
  sname2_found = .FALSE. ! (flags if at least found the right name)
  right_dumping_period = .FALSE. ! (" " " " right dumping period)

  DO i=1, d-1
    IF (SCMop%sname(i) == sname2) THEN

      ! We have found a diagnostic of the correct name
      sname2_found = .TRUE.

      ! But it must have the same dumping period or be constant
      IF (SCMop%dump_step(i) == SCMop%strm(istrm)%dump_step    &
        .OR.                                                   &
          SCMop%timprof(i) == t_const) THEN

        right_dumping_period = .TRUE.

        ! But it must also be defined on the same substep as
        ! the current entry, or have been defined outside of
        ! a sub-stepped part of the model.

        IF (SCMop%substep(i) == SCMop%substep_number .OR.     &
            SCMop%substep(i) == 0) THEN

          ! We have found an entry with the correct name & substep
          SCMop%wd(d) = i
          EXIT

        END IF
      END IF
    END IF
  END DO

  IF (SCMop%wd(d) == 0) THEN
    WRITE(umMessage,'(A)')                                            newline//&
      '=========================================================='//  newline//&
      '| ERROR newdiag:'//                                            newline//&
      '| The following error stems from a call to SCMoutput...'//     newline//&
      '| You have requested that diagnostic "'//TRIM(sname)//'"'//    newline//&
      '| be dependent on diagnostic "'//TRIM(sname2)//'"'//           newline//&
      '| (which is non-constant),'
    CALL umPrint(umMessage,src='newdiag')

    IF (.NOT. sname2_found) THEN
      WRITE(umMessage,'(A)')                                                   &
      '| but the latter diagnostic has not yet been defined.'//       newline//&
      '| Diagnostic "',TRIM(sname),'"'//                              newline//&
      '| will therefore not be calculated or output.'//               newline//&
      '|'//                                                           newline//&
      '| Note: This message will be repeated for every stream'//      newline//&
      '|       you have requested for this diagnostic'
      CALL umPrint(umMessage,src='newdiag')

    ELSE IF (.NOT. right_dumping_period) THEN
      WRITE(umMessage,'(A,I5,A,I5)')                                           &
      '| but you have requested that the former be sent to a'//       newline//&
      '| stream with a dumping period of ',                                    &
      SCMop%strm(istrm)%dump_step,', while the latter is'//           newline//&
      '| not. This is not permitted: if diagnostic "A" is to'//       newline//&
      '| depend on diagnostic "B", and diagnostic "B" is not'//       newline//&
      '| constant, then diagnostic "B" must be calculated for'//      newline//&
      '| every dumping period that is requested for diagnostic "A".'//newline//&
      '| This can be guaranteed by sending "B" to all the streams'//  newline//&
      '| which you have requested for "A". Diagnostic "'                     //&
      TRIM(ADJUSTL(sname)),'"'//                                      newline//&
      '| will therefore not be calculated or output with a dumping'// newline//&
      '| period of ',SCMop%strm(istrm)%dump_step
    ELSE IF (SCMop%substep_number /= 0) THEN
      WRITE(umMessage,'(A,I2)')                                                &
      '| but while that diagnostic appears to be defined within'//    newline//&
      '| a sub-stepped part of the code, no entry could be found'//   newline//&
      '| for it on the current substep of ',SCMop%substep_number
      CALL umPrint(umMessage,src='newdiag')
    ELSE
      WRITE(umMessage,'(A)')                                                   &
      '| but the latter is defined within a sub-stepped part of'//    newline//&
      '| the code while the former is not. This does not make'//      newline//&
      '| sense.'
      CALL umPrint(umMessage,src='newdiag')
    END IF

    WRITE(umMessage,'(A)')                                                     &
      '=========================================================='//  newline//&
      newline
    CALL umPrint(umMessage,src='newdiag')

    d = d_input
    SCMop%nentries = SCMop%nentries-1

    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN

  END IF
ELSE
  IF (timprof == t_div .OR. timprof == t_mult .OR.                             &
      timprof == t_acc_div .OR.                                                &
      timprof == t_acc_mult) THEN
    WRITE(umMessage,'(A,I3,A)')                                                &
      '=========================================================='//  newline//&
      '| ERROR newdiag:'//                                            newline//&
      '| You have requested this diagnostic to be dependent on'//     newline//&
      '| another, but have not specified which.'//                    newline//&
      '| ', TRIM(ADJUSTL(sname)), timprof,                            newline//&
      '=========================================================='
    CALL umPrint(umMessage,src='newdiag')
    d = d_input
    SCMop%nentries = SCMop%nentries-1

    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN

  END IF
END IF

! Assign this diagnostic an integer unique to its sname
! (do this by searching for a previously defined diagnostic
! with the same sname, if you find it give it that number,
! if you don't give it the highest number you came across
! plus one)
SCMop%sname_id(d) = 0

DO i=1, SCMop%nentries-1
  IF (SCMop%sname(i) == sname) THEN

    ! i has the same sname as d, give it the same sname_id
    SCMop%sname_id(d) = SCMop%sname_id(i)
    EXIT
  ELSE

    ! record the largest sname_id, but store as -ve to
    ! indicate that we have not found a match for sname
    SCMop%sname_id(d) = MIN(SCMop%sname_id(d),-SCMop%sname_id(i))
  END IF
END DO

IF (SCMop%sname_id(d) <= 0) THEN

  ! A diagnostic with the same sname was not found: assign a
  ! value one larger than the current largest value of sname_id
  SCMop%sname_id(d) = -SCMop%sname_id(d)+1

  ! Record the number of diagnostics that will be output (with
  ! no double counting from the same diagnostic being
  ! calculated with different dumping periods)
  IF (.NOT. lnot_written) SCMop%n_output = SCMop%n_output+1
END IF

! Record the number of diagnostics that will be sent to
! this stream
IF (.NOT. lnot_written) THEN
  SCMop%strm(istrm)%n_output = SCMop%strm(istrm)%n_output+1
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE newdiag


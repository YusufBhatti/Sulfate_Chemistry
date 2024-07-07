! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Obtain indices of diagnostic entries in SCMop given SCMoutput inputs

SUBROUTINE getdistinctdiags                                                   &
  ( sname_i, lname, units, timprof, domprof, streams, sname2_i                &
  , ndistinct_diags, distinct_diags, SCMop )

! SCMop_type is defined in here...
USE scmoptype_defn
USE s_scmop_mod,                                                            &
    ONLY: stream, streamlist, StreamIsOn, NotWritten, DoNotWrite            &
  , t_inst, t_avg, t_max, t_min, t_acc, t_div, t_mult, t_acc_div            &
  , t_acc_mult, t_const, only_radsteps

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE


! Description:
!   On the first call to this routine with a given sname_i, create
!   the corresponding diagnostic entries in SCMop and return their
!   indices in distinct_diags(1:ndistinct_diags). On subsequent
!   calls, simply look up the previously created entries and return
!   in distinct_diags(1:ndistinct_diags).

! Method:
!   For all streams to which the diagnostic is to be sent, look
!   for a corresponding entry in SCMop with the correct dump_step.
!   If one does not exist create it by calling NEWDIAG. Record its
!   index if it's different to the index found in any previous
!   cycle of the loop over streams. Thus build up a list of the
!   entries in SCMop to which this diagnostic is associated, and
!   return.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran90

  ! See SCMoutput for a description of the input variables
CHARACTER(LEN=*) :: &! In
  sname_i           &! In
, lname             &! In
, units             &! In
, sname2_i           ! In

INTEGER ::  &
  timprof   &! In
, domprof   &! In
, streams    ! In

INTEGER ::                    &
  ndistinct_diags             &! Out
, distinct_diags(maxnstreams)  ! Out

TYPE(SCMop_type) :: SCMop ! InOut The derived-type structure
                          !       containing all the diagnostic
                          !       information

! d will equal this if the entry could not be created
INTEGER, PARAMETER :: unset = 0

INTEGER :: i,j,d,initial_nentries

CHARACTER(LEN=lsname) :: sname,sname2

LOGICAL :: distinct

ndistinct_diags = 0

! Elsewhere we are going to rely on the "sname" for any
! diagnostic being exactly lsname characters long, so here we
! make them the right length.
sname  = sname_i
sname2 = sname2_i

IF (LEN(sname_i) >  lsname .AND. SCMop%first_pass) THEN

  ! Warn the user that we have done this
  WRITE(umMessage,*)'GetDistinctDiags WARNING: diagnostic name truncated: ' &
        ,sname_i(1:LEN(sname_i)),' -> ', sname
  CALL umPrint(umMessage,src='getdistinctdiags')
END IF

! If this is the first timestep, check this "sname" has not been
! used in a previous call to SCMoutput
IF (SCMop%first_pass) THEN
  DO i=1, SCMop%nentries
    IF (SCMop%sname(i) == sname .AND.                                       &
        SCMop%substep(i) == SCMop%substep_number) THEN
      WRITE(umMessage,*)'GetDistinctDiags ERROR: same sname used in >1 '// &
             'call to SCMoutput on same sub-step:',sname
      CALL umPrint(umMessage,src='getdistinctdiags')
      GO TO 9999
    END IF
  END DO
END IF

! Get list of streams to which the diagnostic is to be
! sent. Nominally this is determined by the variable, streams
! (the so-called "hard-wired" list of streams). But namelist
! information can alter this by choosing to ignore the
! hard-wired list and/or to add/remove diagnostics to/from
! the list. The two lists of extra diagnostics to be added to,
! and rejected from, stream X are in SCMop%strm(X)%accept_list and
! SCMop%strm(X)%reject_list respectively.

streamlist = 0 ! This will be the final list of streams encoded into
               ! one integer. For now it is set to "no streams".

! Loop over all streams.
DO j=1, maxnstreams

  ! If this stream is closed do nothing.
  IF (SCMop%strm(j)%switch == 0) THEN
    CYCLE ! Go to next value of j.
  END IF

  ! Shall we pay attention to the list of streams as
  ! specified in the call to SCMoutput?
  IF (SCMop%strm(j)%heed_hardwired /= 0) THEN

    ! Yes. Add stream j if requested in the call.
    IF (StreamIsOn(streams,j) .AND. .NOT. StreamIsOn(streamlist,j)) THEN
      streamlist = streamlist+Stream(j)
    END IF

    ! If it has been requested in the call that this
    ! diagnostic should not be written out, make sure this
    ! is so.
    IF (NotWritten(streams) .AND. .NOT. NotWritten(streamlist)) THEN
      streamlist = DoNotWrite(streamlist)
    END IF
  END IF

  ! Shall we pay attention to the list of diagnostics
  ! requested by namelist to be sent to this stream?
  IF (SCMop%strm(j)%heed_acceptlist /= 0) THEN
    ! Yes.

    ! If the diagnostic is not already being sent to this stream
    IF (.NOT. StreamIsOn(streamlist,j)) THEN

      ! then check if it's in the list of extra diagnostics
      ! to suck into this stream.

      ! Loop over the the names of diagnostics to be accepted
      DO i=1, SIZE(SCMop%strm(j)%accept_list)
        ! Do we have a name match?
        IF (TRIM(SCMop%strm(j)%accept_list(i)) == TRIM(sname)) THEN

          ! The name was in the list for stream j. Send
          ! this diagnostic to stream j.
          streamlist = streamlist+Stream(j)
          EXIT   ! this inner loop
        END IF
      END DO

    END IF

  END IF ! SCMop%strm(j)%heed_acceptlist /= 0

  ! Shall we pay attention to list of diagnostics requested
  ! by namelist to be prevented from going to this stream?
  IF (SCMop%strm(j)%heed_rejectlist /= 0) THEN
    ! Yes.

    ! If the diagnostic is being sent to this stream
    IF (StreamIsOn(streamlist,j)) THEN

      ! then check if it's in the list of diagnostics
      ! which are not to be sent to this stream
      DO i=1, SIZE(SCMop%strm(j)%reject_list)
        IF (TRIM(SCMop%strm(j)%reject_list(i)) == TRIM(sname)) THEN

          ! This diagnostic is not to be sent to this stream.
          streamlist = streamlist-Stream(j)
          EXIT
        END IF
      END DO

    END IF ! StreamIsOn(streamlist,j)

  END IF    ! SCMop%strm(j)%heed_rejectlist /= 0

END DO       ! j=1,maxnstreams

! From now on we will use streamlist instead of streams as the
! integer representing the list of streams to which the
! diagnostic is to be sent.

initial_nentries = SCMop%nentries

! Loop over all output streams
DO j=1, maxnstreams
  ! j represents a stream, is the diagnostic to go
  ! to this stream?
  IF (StreamIsOn(streamlist,j) .AND. SCMop%strm(j)%switch /= 0) THEN

    ! Yes, look to see if the entry already exists in SCMop
    d = unset ! (label d as unset for now)

    DO i=1, SCMop%nentries
      IF (SCMop%sname(i) == sname .AND.                                      &
          SCMop%substep(i) == SCMop%substep_number .AND.                    &
          SCMop%dump_step(i) == SCMop%strm(j)%dump_step) THEN
        ! This diagnostic exists.
        d = i
        EXIT
      END IF
    END DO

    IF (d == unset) THEN
      IF (SCMop%first_pass) THEN

        ! We didn't find an existing diagnostic fitting the
        ! inputs, create a new diagnostic entry in SCMop
        ! DEPENDS ON: newdiag
        CALL newdiag                                                        &
          ( sname, lname, units, timprof, domprof, j                        &
          , NotWritten(streamlist), sname2(1:lsname*MIN(LEN(sname2_i),1))   &
          , d, SCMop )
      ELSE

        ! Should not be having new diagnostics beyond
        ! the first timestep
        WRITE(umMessage,*) &
            'GetDistinctDiags ERROR: new diag after first timestep:',    &
            sname(1:LEN(sname)),                                         &
            SCMop%stepcount,SCMop%daycount,d
        CALL umPrint(umMessage,src='getdistinctdiags')
        GO TO 9999
      END IF
    ELSE
      IF (SCMop%first_pass) THEN

        ! d is the index of a diagnostic with the same sname
        ! as given in the input parameter list and the same
        ! dump_step as stream j. Send diagnostic d to stream
        ! j as well as it's existing streams then.
        IF (.NOT. StreamIsOn(SCMop%streams(d),j)) THEN

          SCMop%streams(d) = SCMop%streams(d) + Stream(j)

          ! Increment n_output to count the number we're
          ! going to output to this stream
          IF (.NOT. NotWritten(streamlist)) THEN
            SCMop%strm(j)%n_output = SCMop%strm(j)%n_output+1
          END IF

        ELSE

          ! Diagnostic d is already going to this stream,
          ! this should not ocurr.
          WRITE(umMessage,*)'GetDistinctDiags ERROR: Same diag. sent '//    &
                 'to same stream twice',d,sname(1:LEN(sname)),              &
                 j,SCMop%streams(d)
          CALL umPrint(umMessage,src='getdistinctdiags')
        END IF
      ELSE
        ! Check that diagnostic d is set up to go to
        ! this stream
        IF (.NOT. StreamIsOn(SCMop%streams(d),j)) THEN
          WRITE(umMessage,*)'GetDistinctDiags ERROR: the requested '//      &
                 'streams for this diagnostic have changed',                &
                 d,sname(1:LEN(sname)),j,SCMop%streams(d)
          CALL umPrint(umMessage,src='getdistinctdiags')
        END IF
      END IF
    END IF

    ! We should now have a value for d, but it may still be
    ! unset if an error ocurred in newdiag
    IF (d /= unset) THEN

      ! Is this diagnostic distinct? i.e. is the value of d at
      ! this point different than for any previous value of j
      ! in this loop?
      distinct = .TRUE.

      DO i=1, ndistinct_diags
        IF (distinct_diags(i) == d) THEN
          distinct = .FALSE.
          EXIT
        END IF
      END DO

      IF (distinct) THEN
        ndistinct_diags = ndistinct_diags+1
        distinct_diags(ndistinct_diags) = d

        ! Make a requirement that diagnostic entries must be
        ! created in order (see use of diag_mem in SCMoutput)
        IF (ndistinct_diags >  1) THEN
          IF (d /= distinct_diags(ndistinct_diags-1)+1) THEN
            WRITE(umMessage,*) &
                'GetDistinctDiags ERROR: non-consecutive diags',         &
                distinct_diags(1:ndistinct_diags)
            CALL umPrint(umMessage,src='getdistinctdiags')
          END IF
        END IF
      END IF
    END IF

  END IF                  ! (StreamIsOn(streams,j))
END DO                     ! j=1,maxnstreams

9999 CONTINUE

RETURN

END SUBROUTINE getdistinctdiags

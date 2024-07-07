! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Write out info about the output diagnostics for one or all streams

SUBROUTINE write_scumlist (SCMop, istrm)

! SCMop_type is defined in here...
USE scmoptype_defn
USE s_scmop_mod,  &
    ONLY: StreamIsOn, NotWritten, AnyStreamOn

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE


! Description:
!   For a given stream, writes information about each domain profile and
!   diagnostic that is active in that stream to the stream's output file
!   as part of its header. PV-wave can then read this info to make sense
!   of the rest of the data file. If the requested stream number is zero,
!   a new file will be opened and all the same information will be written
!   to it, but now pertaining to all streams. Such a file will list all
!   domain profiles and diagnostics, and will have additional comments.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran 90

TYPE(SCMop_type) :: &
  SCMop              ! In The derived-type structure containing
                     !    all the diagnostic information

INTEGER ::                    &
  istrm                       &! In The stream in question. Zero=all streams
, ndiags                      &! The number of diagnostic entries found in
                               ! SCMop which we will write about.
, diags(SCMop%nentries)       &! Contains in 1:niags the indices of the
                               ! entries in SCMop
, domprof_sparse(maxndomprof)  ! A sparse array indicating the domain
                               ! profiles possessed by the ndiags entries.

INTEGER :: i,j    ! Counters

! Character function that incorporates the substep number into
! the short name of a diagnostic. Somewhat longer than a normal
! short name.
CHARACTER (LEN=lsname+10) :: add_substep_to_sname,short_name

! The longest possible sname after the substep has been appended
INTEGER :: sname_length

INTEGER :: UNIT ! The unit to which we'll write
CHARACTER (LEN=100) :: FMT
CHARACTER (LEN=lsname) :: sname
INTEGER :: substep

! Make a list the diagnostics we're going to write and
! which domain profiles they use. ndiags will hold the total
! number of diagnostic entries in SCMop fitting the bill, and
! domprof_sparse will be a sparse array indicating which domain
! profiles are used by these diagnostics.
ndiags = 0
domprof_sparse = 0

! If a specific stream has been specified...
IF (istrm /= 0) THEN

  ! ... select only those diagnostics going to that stream
  DO i=1, SCMop%nentries
    IF (StreamIsOn(SCMop%streams(i),istrm) .AND. .NOT.              &
          NotWritten(SCMop%streams(i))) THEN
      ndiags = ndiags+1
      diags(ndiags) = i
      domprof_sparse(SCMop%domprof(i)) = 1
    END IF
  END DO

  ! Make a check just for the sake of it
  IF (ndiags /= SCMop%strm(istrm)%n_output) THEN
    WRITE(umMessage,*)'dump_streams_init ERROR: an inconsistency '//         &
           'has ocurred !',ndiags,SCMop%strm(istrm)%n_output,     &
           istrm,SCMop%nentries
    CALL umPrint(umMessage,src='write_scumlist')

    ! Switch the diagnostic system off so dodgy data is
    ! not mistakenly used in good faith.
    WRITE(umMessage,*)' -> switching diagnostic system off!'
    CALL umPrint(umMessage,src='write_scumlist')
    SCMop%on = .FALSE.
  END IF
ELSE
  ! If strm=0 then use all diagnostics which are going to
  ! at least one stream. But since the same diagnostic can
  ! go to several streams, and thus have several entries in
  ! SCMop, then make a condition that if the i'th entry has
  ! the same short name as the (i-1)'th, then skip the i'th
  ! to avoid multiple identical lines. We -do- want multiple
  ! lines for different sub-steps though, so allow identical
  ! consecutive snames if their sub-step numbers are different.

  sname = ''
  substep = -999

  DO i=1, SCMop%nentries
    IF (AnyStreamOn(SCMop%streams(i)) .AND.                        &
         (SCMop%sname(i) /= sname .OR.                            &
          SCMop%substep(i) /= substep)) THEN

      ndiags = ndiags+1
      diags(ndiags) = i
      domprof_sparse(SCMop%domprof(i)) = 1
      sname = SCMop%sname(i)
      substep = SCMop%substep(i)

    END IF
  END DO
END IF

! Determine the unit number to write to
IF (istrm /= 0) THEN

  ! Write to the unit of the specified stream. The file
  ! should already be open.
  UNIT = SCMop%strm(istrm)%op_unit
ELSE

  ! Write to a new file we'll attach to unit 10
  UNIT = 10
  OPEN (UNIT=unit,FILE='scumlist.dat')
END IF

!-------------------------------------------------------------
!     Write about the domain profiles
!-------------------------------------------------------------

  ! Write the number of domain profiles used by
  ! diagnostics in this stream
IF (istrm == 0) WRITE(UNIT,*)'No. of domain profiles:'
WRITE(UNIT,'(I3)') SUM(domprof_sparse)

! Write the format with which we will write some of
! the upcoming lines
FMT = '(I3,1X,A15,I3,1X,I3,1X,I3,1X,I3,1X,I3,1X,I3)'
IF (istrm == 0) WRITE(UNIT,'(A)')'Line format:'
WRITE(UNIT,'(A)') TRIM(FMT)

! Write info about each domain profile
IF (istrm == 0) WRITE(UNIT,*)'List of domain profiles:'

DO i=1, maxndomprof
  IF (domprof_sparse(i) == 1) THEN
    WRITE(UNIT,FMT)                                          &
      i,SCMop%d_name(i),                                     &
      SCMop%d_rowa1(i),SCMop%d_rowa2(i),                     &
      SCMop%d_rowb1(i),SCMop%d_rowb2(i),                     &
      SCMop%d_lev1(i),SCMop%d_lev2(i)
  END IF
END DO

!-------------------------------------------------------------
!     Write about the diagnostics themselves
!-------------------------------------------------------------

  ! Write the number
IF (istrm == 0) WRITE(UNIT,*)'No. of diagnostics:'
WRITE(UNIT,'(I3)') ndiags

! If there are multiple sub-steps then the short-name will
! have "_#N" tagged onto it, where N is the sub-step number.
! In this case we need to increase the amount of space made
! available for the short name.

IF (SCMop%num_substeps <= 1) THEN
  sname_length = lsname
ELSE
  sname_length = lsname                                      &
               + 2+(INT(LOG10(REAL(SCMop%num_substeps)))+1)
END IF

! Compose the format of the line that will describe
! each diagnostic
! NOTE: for the long-name field, ensure we write the original field-width
! llname_old=50, to preserve backwards compatibility with old scripts for
! reading SCM text-format output files (the long-name field-width llname
! can therefore be increased for netcdf output without breaking old
! text-format SCM-read scripts).
WRITE(FMT,'(A,I2, A,I2, A,I2, A)')'(I3' //                   &
     ',1X,A',sname_length,                                   &
     ',1X,A',llname_old,                                     &
     ',1X,A',lunits,                                         &
     ',I2)'
! Write the format to the file
IF (istrm == 0) WRITE(UNIT,'(A)')'Line format:'
WRITE(UNIT,'(A)') TRIM(FMT)

! Write a 1-line description of each diagnostic
! consisting of a unique integer ID, its short name,
! its long name, its units and its domain profile.
IF (istrm == 0) WRITE(UNIT,'(A)')'List of diagnostics '//    &
                '(i,sname,lname,units,domprof):'
DO i=1, ndiags
  j = diags(i)
  ! DEPENDS ON: add_substep_to_sname
  short_name = add_substep_to_sname(SCMop,j)
  WRITE(UNIT,FMT)                                            &
       SCMop%sname_id(j),short_name,                         &
       SCMop%lname(j)(1:llname_old), SCMop%units(j),         &
       SCMop%domprof(j)
END DO

IF (istrm == 0) THEN
  ! We need to close the new file we opened above
  CLOSE(UNIT)
END IF

RETURN
END SUBROUTINE write_scumlist


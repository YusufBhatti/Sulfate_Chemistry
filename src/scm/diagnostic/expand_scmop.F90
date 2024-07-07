! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Increase the size of the SCM diagnostic arrays in SCMop

SUBROUTINE expand_scmop (SCMop)
USE UM_types
! SCMop_type is defined in here...
USE scmoptype_defn
USE s_scmop_mod

USE umPrintMgr, ONLY: umPrint, umMessage, newline
USE yomhook,    ONLY: lhook, dr_hook
USE parkind1,   ONLY: jprb, jpim

IMPLICIT NONE


! Description:
!   Increases the maximum number of diagnostics entries
!   allowed in SCMop by reallocating the relevant arrays.

! Method:
!   On the first timestep, as calls to SCMoutput are being made,
!   memory has to be allocated to arrays to hold the resulting
!   information. Since it is not known at the start of the run how
!   many calls to SCMoutput there will be and what their input
!   parameters are, no memory is allocated at the outset and the
!   variable SCMop%maxnentries (the maximum no. of diagnostic
!   "entries" that the arrays in SCMop can handle before they run
!   out of space) is zero. In this case most of the
!   statements in this routine are ignored (since they start with
!   "if (maxnentries >  0)"), and the arrays are simply allocated
!   with a size equal to "chunk". On subsequent calls the contents
!   of the allocatable arrays are copied into temporary arrays,
!   de-allocated, re-allocated with their original size plus
!   "chunk", and then the data is copied back in. i.e. the arrays
!   are re-allocated with large sizes without losing their
!   information.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran 90

TYPE(SCMop_type) :: SCMop ! INOUT The derived-type structure
                          ! containing all the diagnostic
                          ! information

! Temporary arrays to hold data while SCMop arrays are being reallocated
INTEGER, ALLOCATABLE :: Ixxx(:)
LOGICAL, ALLOCATABLE :: Lxxx(:)
CHARACTER (LEN=llname), ALLOCATABLE :: Cxxx(:)
TYPE(allocatable_array), ALLOCATABLE :: Dxxx(:)

! The number by which to increment SCMop%maxnentries
INTEGER, PARAMETER :: chunk=100

! Holds SCMop%maxnentries
INTEGER :: maxnentries

! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXPAND_SCMOP'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (SCMop%nentries /= SCMop%maxnentries) THEN
  WRITE(umMessage,'(A,I3,A)')                                                 &
    '=========================================================='//   newline//&
    '| Warning expand_scmop:'//                                      newline//&
    '| SCMop is being expanded before it is full, this could'//      newline//&
    '| be dangerous.'//                                              newline//&
    '=========================================================='
  CALL umPrint(umMessage,src='expand_scmop')
END IF

maxnentries = SCMop%maxnentries

! Allocate the space required for the temporary arrays
IF (maxnentries >  0) ALLOCATE(Cxxx(maxnentries))
IF (maxnentries >  0) ALLOCATE(Ixxx(maxnentries))
IF (maxnentries >  0) ALLOCATE(Lxxx(maxnentries))
IF (maxnentries >  0) ALLOCATE(Dxxx(maxnentries))

! Increase the size of all the arrays in SCMop associated
! to specific diagnostics...

IF (maxnentries >  0) Cxxx = SCMop%sname
IF (maxnentries >  0) DEALLOCATE(SCMop%sname)
ALLOCATE(SCMop%sname(maxnentries+chunk))
IF (maxnentries >  0) SCMop%sname(1:maxnentries) = Cxxx

IF (maxnentries >  0) Cxxx = SCMop%lname
IF (maxnentries >  0) DEALLOCATE(SCMop%lname)
ALLOCATE(SCMop%lname(maxnentries+chunk))
IF (maxnentries >  0) SCMop%lname(1:maxnentries) = Cxxx

IF (maxnentries >  0) Cxxx = SCMop%units
IF (maxnentries >  0) DEALLOCATE(SCMop%units)
ALLOCATE(SCMop%units(maxnentries+chunk))
IF (maxnentries >  0) SCMop%units(1:maxnentries) = Cxxx

IF (maxnentries >  0) Ixxx = SCMop%domprof
IF (maxnentries >  0) DEALLOCATE(SCMop%domprof)
ALLOCATE(SCMop%domprof(maxnentries+chunk))
IF (maxnentries >  0) SCMop%domprof(1:maxnentries) = Ixxx

IF (maxnentries >  0) Ixxx = SCMop%timprof
IF (maxnentries >  0) DEALLOCATE(SCMop%timprof)
ALLOCATE(SCMop%timprof(maxnentries+chunk))
IF (maxnentries >  0) SCMop%timprof(1:maxnentries) = Ixxx

IF (maxnentries >  0) Ixxx = SCMop%streams
IF (maxnentries >  0) DEALLOCATE(SCMop%streams)
ALLOCATE(SCMop%streams(maxnentries+chunk))
IF (maxnentries >  0) SCMop%streams(1:maxnentries) = Ixxx

IF (maxnentries >  0) Ixxx = SCMop%dump_step
IF (maxnentries >  0) DEALLOCATE(SCMop%dump_step)
ALLOCATE(SCMop%dump_step(maxnentries+chunk))
IF (maxnentries >  0) SCMop%dump_step(1:maxnentries) = Ixxx

IF (maxnentries >  0) Ixxx = SCMop%nadd2dump
IF (maxnentries >  0) DEALLOCATE(SCMop%nadd2dump)
ALLOCATE(SCMop%nadd2dump(maxnentries+chunk))
IF (maxnentries >  0) SCMop%nadd2dump(1:maxnentries) = Ixxx

IF (maxnentries >  0) Lxxx = SCMop%only_radsteps
IF (maxnentries >  0) DEALLOCATE(SCMop%only_radsteps)
ALLOCATE(SCMop%only_radsteps(maxnentries+chunk))
IF (maxnentries >  0) SCMop%only_radsteps(1:maxnentries) = Lxxx

IF (maxnentries >  0) Ixxx = SCMop%ncols
IF (maxnentries >  0) DEALLOCATE(SCMop%ncols)
ALLOCATE(SCMop%ncols(maxnentries+chunk))
IF (maxnentries >  0) SCMop%ncols(1:maxnentries) = Ixxx

IF (maxnentries >  0) Ixxx = SCMop%nrows
IF (maxnentries >  0) DEALLOCATE(SCMop%nrows)
ALLOCATE(SCMop%nrows(maxnentries+chunk))
IF (maxnentries >  0) SCMop%nrows(1:maxnentries) = Ixxx

IF (maxnentries >  0) Ixxx = SCMop%nlevs
IF (maxnentries >  0) DEALLOCATE(SCMop%nlevs)
ALLOCATE(SCMop%nlevs(maxnentries+chunk))
IF (maxnentries >  0) SCMop%nlevs(1:maxnentries) = Ixxx

IF (maxnentries >  0) Ixxx = SCMop%nelements
IF (maxnentries >  0) DEALLOCATE(SCMop%nelements)
ALLOCATE(SCMop%nelements(maxnentries+chunk))
IF (maxnentries >  0) SCMop%nelements(1:maxnentries) = Ixxx

IF (maxnentries >  0) Ixxx = SCMop%sname_id
IF (maxnentries >  0) DEALLOCATE(SCMop%sname_id)
ALLOCATE(SCMop%sname_id(maxnentries+chunk))
IF (maxnentries >  0) SCMop%sname_id(1:maxnentries) = Ixxx

IF (maxnentries >  0) Ixxx = SCMop%wd
IF (maxnentries >  0) DEALLOCATE(SCMop%wd)
ALLOCATE(SCMop%wd(maxnentries+chunk))
IF (maxnentries >  0) SCMop%wd(1:maxnentries) = Ixxx

IF (maxnentries >  0) Ixxx = SCMop%lastencounter
IF (maxnentries >  0) DEALLOCATE(SCMop%lastencounter)
ALLOCATE(SCMop%lastencounter(maxnentries+chunk))
IF (maxnentries >  0) SCMop%lastencounter(1:maxnentries) = Ixxx

IF (maxnentries > 0) Ixxx = SCMop%substep
IF (maxnentries > 0) DEALLOCATE(SCMop%substep)
ALLOCATE(SCMop%substep(maxnentries+chunk))
IF (maxnentries > 0) SCMop%substep(1:maxnentries) = Ixxx

IF (maxnentries >  0) Dxxx = SCMop%diag
IF (maxnentries >  0) DEALLOCATE(SCMop%diag)
ALLOCATE(SCMop%diag(maxnentries+chunk))
IF (maxnentries >  0) SCMop%diag(1:maxnentries) = Dxxx

IF (maxnentries >  0) Ixxx = SCMop%netcdf_id
IF (maxnentries >  0) DEALLOCATE(SCMop%netcdf_id)
ALLOCATE(SCMop%netcdf_id(maxnentries+chunk))
IF (maxnentries >  0) SCMop%netcdf_id(1:maxnentries) = Ixxx

IF (maxnentries >  0) DEALLOCATE(Cxxx)
IF (maxnentries >  0) DEALLOCATE(Ixxx)
IF (maxnentries >  0) DEALLOCATE(Lxxx)
IF (maxnentries >  0) DEALLOCATE(Dxxx)

SCMop%maxnentries = SCMop%maxnentries+chunk

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE expand_scmop


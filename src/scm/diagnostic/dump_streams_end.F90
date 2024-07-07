! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Close SCM output files

SUBROUTINE dump_streams_end(SCMop)

USE netcdf
USE s_scmop_mod, ONLY: SCMop_type, maxnstreams, incdf

USE umPrintMgr, ONLY: umPrint, umMessage, newline
USE yomhook,    ONLY: lhook, dr_hook
USE parkind1,   ONLY: jprb, jpim

IMPLICIT NONE


! Description:
!   Close all open SCM output files

! Method:
!

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Language: Fortran90

TYPE(SCMop_type) :: SCMop ! InOut The derived-type structure
                          !       containing all the diagnostic
                          !       information

INTEGER :: n,STATUS

! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DUMP_STREAMS_END'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Loop over each output stream
DO n=1, maxnstreams

  ! Is this stream switched on and are any diagnostics being sent
  ! to it?
  IF (SCMop%strm(n)%switch /= 0 .AND.                              &
      SCMop%strm(n)%n_output > 0) THEN

    ! This stream was opened. Close it now.
    IF (SCMop%strm(n)%FORMAT >= 0 .AND.                            &
        SCMop%strm(n)%FORMAT <= 3) THEN
      CLOSE(SCMop%strm(n)%op_unit)

    ELSE IF (SCMop%strm(n)%FORMAT == 4) THEN
      ! Close the NetCDF file
      STATUS = Nf90_Close(INT(SCMop%strm(n)%op_unit,incdf))

    ELSE
      WRITE(umMessage,'(A)')                                                   &
        '=========================================================='//newline//&
        '| Dump_Streams_End ERROR:'                                 //newline//&
        '| Unknown format for stream ',n,                             newline//&
        '=========================================================='
      CALL umPrint(umMessage,src='dump_streams_end')
    END IF

  END IF
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE dump_streams_end


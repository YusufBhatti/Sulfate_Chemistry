! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE IOERROR----------------------------------------
!
!    Purpose: Prints out a message after using buffer in/out when
!             either a return code < 0.0 is encountered
!             by UNIT function or value returned by LENGTH
!             differs from length of I/O request.
!
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Dump I/O

SUBROUTINE ioerror(string,error,len_io,len_io_req)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE


INTEGER  :: len_io      ! Number of 64-bit words transferred as
                        !  registered  by LENGTH function
INTEGER  :: len_io_req  ! Number of 64-bit words requested for
                        !  transfer via BUFFER IN/OUT
CHARACTER(LEN=*) :: string

REAL     :: error       ! Error code returned by UNIT function

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOERROR'

! -------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
WRITE(umMessage,'('' **FATAL ERROR WHEN READING/WRITING MODEL DUMP**'')')
CALL umPrint(umMessage,src='ioerror')
WRITE(umMessage,'('' '',A)') TRIM(string)
CALL umPrint(umMessage,src='ioerror')
WRITE(umMessage,'('' Error code = '',F6.2)') error
CALL umPrint(umMessage,src='ioerror')
WRITE(umMessage,'('' Length requested            = '',I9)') len_io_req
CALL umPrint(umMessage,src='ioerror')
WRITE(umMessage,'('' Length actually transferred = '',I9)') len_io
CALL umPrint(umMessage,src='ioerror')
WRITE(umMessage,'('' Fatal error codes are as follows:'')')
CALL umPrint(umMessage,src='ioerror')
CALL umPrint(                                                    &
    ' -1.0 Mismatch between actual and requested data length')
WRITE(umMessage,'(''  0.0 End-of-file was read'')')
CALL umPrint(umMessage,src='ioerror')
WRITE(umMessage,'(''  1.0 Error occurred during read'')')
CALL umPrint(umMessage,src='ioerror')
WRITE(umMessage,'(''  2.0 Other disk malfunction'')')
CALL umPrint(umMessage,src='ioerror')
WRITE(umMessage,'(''  3.0 File does not exist'')')
CALL umPrint(umMessage,src='ioerror')
WRITE(umMessage,'('' ***********************************************'')')
CALL umPrint(umMessage,src='ioerror')

WRITE(0,'(//)')
WRITE(0,'('' **FATAL ERROR WHEN READING/WRITING MODEL DUMP**'')')
WRITE(0,'('' '',A)') TRIM(string)
WRITE(0,'('' Error code = '',F6.2)') error
WRITE(0,'('' Length requested            = '',I9)') len_io_req
WRITE(0,'('' Length actually transferred = '',I9)') len_io
WRITE(0,'('' Fatal Error codes are as follows:'')')
WRITE(0,'(A)')                                                    &
' -1.0 Mismatch between actual and requested data length'
WRITE(0,'(''  0.0 End-of-file was read'')')
WRITE(0,'(''  1.0 Error occurred during read'')')
WRITE(0,'(''  2.0 Other disk malfunction'')')
WRITE(0,'(''  3.0 File does not exist'')')
WRITE(0,'('' ***********************************************'')')

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ioerror

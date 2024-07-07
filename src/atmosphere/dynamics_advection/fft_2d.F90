! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine fft_2d
SUBROUTINE fft_2d(global_rows,global_row_length,spectra2          &
              ,spectra_im2,init,trigm,trign,work,ifail)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE

!
! Description:
!   Dummy subroutine.  To be replaced by a subroutine which calculates
!   2D FFTs.
!
! Method:
!   Passes the original field back without performing any calculations
!   on it.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
!  Dates should have leading zeroes in dd/mm/yy format
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Global variables (#include statements etc):

INTEGER ::                                                        &
 global_rows                                                      &
,global_row_length                                                &
,ifail

REAL ::                                                           &
 spectra2(global_row_length*global_rows)                          &
,spectra_im2(global_row_length*global_rows)                       &
,trign(2*global_row_length)                                       &
                               !  HOLDS TRIGONOMETRIC TERMS
,trigm(2*global_rows)                                             &
                               !  USED IN FFT'S
,work(2*global_row_length*global_rows)

CHARACTER(LEN=1) :: init

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FFT_2D'


!---------------------------------------------------------------------
! Section 1.
!---------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
WRITE(umMessage,*)'**************WARNING*****************'
CALL umPrint(umMessage,src='fft_2d')
WRITE(umMessage,*)'******FFT_2D is a dummy routine*******'
CALL umPrint(umMessage,src='fft_2d')
WRITE(umMessage,*)'*****original field is passed back****'
CALL umPrint(umMessage,src='fft_2d')
WRITE(umMessage,*)'**************************************'
CALL umPrint(umMessage,src='fft_2d')

!!    END OF ROUTINE FFT_2D
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fft_2d

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE PR_REHDA---------------------------------------
!
!    Purpose: Prints out real constants record and checks
!             validity of information.
!
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Documentation: Unified Model Documentation Paper No F3
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Dump I/O

SUBROUTINE pr_rehda(realhd,len_realhd)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE


INTEGER ::                                                        &
 len_realhd !IN Length of real header

REAL ::                                                           &
 realhd(len_realhd) !IN Real header

! Local variables:---------------------------------------------
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PR_REHDA'
!--------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
WRITE(umMessage,'('' '')')
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'('' REAL CONSTANTS'')')
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'('' --------------'')')
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'('' E-W grid spacing in degrees -'',e12.4)')             &
realhd(1)
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'('' N-S grid spacing in degress -'',e12.4)')             &
realhd(2)
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'('' Latitude of first row in degrees -'',e12.4)')        &
realhd(3)
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'(A,e12.4)')                                              &
' Longitude of first point in a row in degrees -',realhd(4)
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'(A,e12.4)')                                              &
' Real latitude of pseudo North Pole in degrees - ',realhd(5)
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'(A,e12.4)')                                              &
' Real longitude of pseudo North Pole in degrees - ',realhd(6)
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'(A,e12.4)')                                              &
' Grid orientation in degrees - ',realhd(7)
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'(26x,A)')                                                &
' Year        Day       Hour      Minute     Second'
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'(A,5e12.4)')                                             &
' Atmosphere time = ',(realhd(i),i=8,12)
CALL umPrint(umMessage,src='pr_rehda')
WRITE(umMessage,'('' Mass, energy, energy drift = '',3e12.4)')            &
realhd(19),realhd(20),realhd(21)
CALL umPrint(umMessage,src='pr_rehda')

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pr_rehda

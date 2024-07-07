! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    MODULE POSERROR_MOD ---------------------------------------
!
!    Purpose:
!             Prints out a message when position of a data block as
!             pointed to by fixed length header differs from actual
!             position in model dump.
!
!
!    Programming standard:
!             Unified Model Documentation Paper No 3
!
!             Code Owner: Please refer to the UM file CodeOwners.txt
!             This file belongs in section: Dump I/O

MODULE poserror_mod

IMPLICIT NONE

PRIVATE
PUBLIC poserror

CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName = 'POSERROR_MOD'

CONTAINS

SUBROUTINE poserror(string,start_block,head_pos,head_address)

USE parkind1,   ONLY: jprb, jpim
USE umPrintMgr, ONLY: umPrint, umMessage
USE yomhook,    ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: string
              !IN Description of block
INTEGER, INTENT(IN) :: start_block
              !IN Actual position of data block
INTEGER, INTENT(IN) :: head_pos
              !IN Position in FIXHD of pointer
INTEGER, INTENT(IN) :: head_address
              !IN Position in file pointed to by FIXHD(HEAD_POS)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='POSERROR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

WRITE(umMessage,'('' ******FATAL ERROR WHEN READING MODEL DUMP******'')')
CALL umPrint(umMessage,src='poserror')

WRITE(umMessage,'('' Conflict between start position of '',A)') string
CALL umPrint(umMessage,src='poserror')

WRITE(umMessage,'(A,i3,A,i9)')                                             &
    ' Block and Pointer in Fixed Length Header: fixhd(',                   &
    head_pos,') =',head_address
CALL umPrint(umMessage,src='poserror')

WRITE(umMessage,'('' Current position in file ='',i9,'' words in'')')      &
    start_block
CALL umPrint(umMessage,src='poserror')

WRITE(umMessage,'('' ***********************************************'')')
CALL umPrint(umMessage,src='poserror')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE poserror

END MODULE poserror_mod


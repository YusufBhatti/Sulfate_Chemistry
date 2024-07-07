#if ! defined(DRHOOK)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE yomhook

!
! Description:
!   Dummy module to replace the DrHook library.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dummy libraries
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.

USE parkind1, ONLY: jpim, jprb
IMPLICIT NONE

LOGICAL, PARAMETER :: lhook = .FALSE.

INTERFACE dr_hook
MODULE PROCEDURE dr_hook
MODULE PROCEDURE dr_hook_comm
END INTERFACE dr_hook

CONTAINS

SUBROUTINE dr_hook(NAME,code,handle)
IMPLICIT NONE

!Arguments
CHARACTER(LEN=*),   INTENT(IN)    :: NAME
INTEGER(KIND=jpim), INTENT(IN)    :: code
REAL(KIND=jprb),    INTENT(INOUT) :: handle

!Nothing to do

RETURN
END SUBROUTINE dr_hook

SUBROUTINE dr_hook_comm(NAME,code,handle,lcomm,comm)
IMPLICIT NONE

!Arguments
CHARACTER(LEN=*),   INTENT(IN)    :: NAME
INTEGER(KIND=jpim), INTENT(IN)    :: code
LOGICAL(KIND=jpim), INTENT(IN)    :: lcomm
INTEGER(KIND=jpim), INTENT(IN)    :: comm
REAL(KIND=jprb),    INTENT(INOUT) :: handle

!Nothing to do

RETURN
END SUBROUTINE dr_hook_comm

END MODULE yomhook
#endif

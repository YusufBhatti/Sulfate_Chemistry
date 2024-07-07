! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************

!
! Description:
!     Define a set of UM inclusive timers for coupling operations.
!     Adds a synchronisation before the timers in order
!     to obtain a unique measurement on all tasks and avoid load
!     imbalance.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!
!=====================================================================

MODULE oasis_timers

USE coupling_control_mod, ONLY: l_oasis_timers
 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

! Main control switch
LOGICAL :: l_synchronise=.TRUE.

! List of timer names for the UM

CHARACTER(LEN=*) :: InitName
PARAMETER (InitName='oasis3_inita2o')
CHARACTER(LEN=*) :: GetName
PARAMETER (GetName='oasis3_geto2a')
CHARACTER(LEN=*) :: PutName
PARAMETER (PutName='oasis3_puta2o')
CHARACTER(LEN=*) :: GridName
PARAMETER (GridName='oasis3_grid')
CHARACTER(LEN=*) :: ModelName
PARAMETER (ModelName='Atmosphere Time')
CHARACTER(LEN=*) :: CouplingName
PARAMETER (CouplingName='Coupling Time')

CONTAINS

SUBROUTINE initstarts
IMPLICIT NONE

IF (l_synchronise) CALL sync
CALL timer(InitName,3)
END SUBROUTINE initstarts

SUBROUTINE initends
IMPLICIT NONE

IF (l_synchronise) CALL sync
CALL timer(InitName,4)
END SUBROUTINE initends

SUBROUTINE getstarts
IMPLICIT NONE
IF (l_synchronise) CALL sync
CALL timer(GetName,3)
END SUBROUTINE getstarts

SUBROUTINE getends
IMPLICIT NONE

IF (l_synchronise) CALL sync
CALL timer(GetName,4)
END SUBROUTINE getends

SUBROUTINE putstarts
IMPLICIT NONE

IF (l_synchronise) CALL sync
CALL timer(PutName,3)
END SUBROUTINE putstarts

SUBROUTINE putends
IMPLICIT NONE

IF (l_synchronise) CALL sync
CALL timer(PutName,4)
END SUBROUTINE putends

SUBROUTINE gridstarts
IMPLICIT NONE

IF (l_synchronise) CALL sync
CALL timer(GridName,3)
END SUBROUTINE gridstarts

SUBROUTINE gridends
IMPLICIT NONE

IF (l_synchronise) CALL sync
CALL timer(GridName,4)
END SUBROUTINE gridends

SUBROUTINE sync
IMPLICIT NONE

INTEGER :: icode, icomm

CALL gc_get_communicator(icomm, icode)
CALL gcg_gsync(icomm, icode)
END SUBROUTINE sync

END MODULE oasis_timers

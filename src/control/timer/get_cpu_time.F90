! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE get_cpu_time_mod

IMPLICIT NONE

CONTAINS
!
! Gets the cpu time from the system
!
! Function Interface:
REAL FUNCTION get_cpu_time()

IMPLICIT NONE
!
! Description:
!   There used to be no Fortran standard for calculating CPU time.
!   This routine was therefore an in interface to propriatory library
!   functions, but now uses the Fortran 95 intrinsic.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Timer
!
! Code Description:
!   Language: Fortran 95
!
REAL :: cputime

CALL CPU_TIME(cputime)
get_cpu_time = cputime

RETURN
END FUNCTION get_cpu_time
END MODULE get_cpu_time_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module contains the interfaces to call c code within fortran.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C Code
!

MODULE io_timing_interfaces_mod

IMPLICIT NONE

INTERFACE
SUBROUTINE io_total_timings() BIND(c,NAME="io_total_timings")
END SUBROUTINE
END INTERFACE

END MODULE io_timing_interfaces_mod

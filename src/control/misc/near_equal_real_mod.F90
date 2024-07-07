! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Purpose:
!    Check whether two REAL values are nearly equal (to ~32bit).
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Misc


MODULE near_equal_real_mod

IMPLICIT NONE

CONTAINS

LOGICAL FUNCTION near_equal_real(p1, p2)

IMPLICIT NONE

REAL, INTENT(IN) :: p1, p2
near_equal_real = ((ABS(p1-p2))  >   (1.0e-6 * ABS(p1+p2)))
RETURN

END FUNCTION near_equal_real

END MODULE near_equal_real_mod

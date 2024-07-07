! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module: DRHOOK_CONTROL ---------------------------------------------
!
!  Purpose: Allows the value of lhook to be managed cleanly,
!           without a direct assignment in other code.  The dummy
!           DrHook code has lhook set as a parameter for performance
!           reasons, which prevents assignment.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Dummy libraries
!
!  Code description:
!    Language: Fortran 90.
!    This code is written to UM programming standards version 9.0.

MODULE drhook_control_mod

USE yomhook, ONLY: lhook

IMPLICIT NONE
PRIVATE

PUBLIC :: drhook_control_enable
PUBLIC :: drhook_control_disable

!----------------------------------------------------------------------
! Contained functions/subroutines
!----------------------------------------------------------------------
CONTAINS

!----------------------------------------------------------------------
! 
!----------------------------------------------------------------------

SUBROUTINE drhook_control_enable()
IMPLICIT NONE

#if defined(DRHOOK)
lhook = .TRUE.
#endif

END SUBROUTINE drhook_control_enable

!----------------------------------------------------------------------
! 
!----------------------------------------------------------------------

SUBROUTINE drhook_control_disable()
IMPLICIT NONE

#if defined(DRHOOK)
lhook = .FALSE.
#endif

END SUBROUTINE drhook_control_disable

END MODULE drhook_control_mod



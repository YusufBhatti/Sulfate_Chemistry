! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine: LEVEL_RLEVEL ------------------------------------------
!
!  Purpose: To return a real value even though the routine is called
!  with integer arguments.
!
!  Tested under compiler:   cft77
!  Tested under OS version: UNICOS 5.1
!
!  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!  Logical components covered: ...
!
!  Project task: ...
!
!  External documentation:
!
!  -------------------------------------------------------------------
!  Interface and arguments: ------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
SUBROUTINE level_rlevel(int_level,real_level,real_level_out)
IMPLICIT NONE
INTEGER ::                                                        &
     int_level              !    first dimension of the lookup
REAL ::                                                           &
     real_level,                                                  &
                            !    secnd dimension of the lookup
     real_level_out         !    secnd dimension of the lookup
!
real_level_out=real_level


RETURN
END SUBROUTINE level_rlevel

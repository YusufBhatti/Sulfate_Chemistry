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

MODULE fort2c_exceptions_interfaces

! DEPENDS ON: exceptions.o

USE, INTRINSIC :: ISO_C_BINDING, ONLY:                                         &
  C_LOC,                                                                       &
  C_CHAR,                                                                      &
  C_PTR,                                                                       &
  C_INT64_T,                                                                   &
  C_INT,                                                                       &
  C_FUNPTR,                                                                    &
  C_FUNLOC

IMPLICIT NONE

PRIVATE

PUBLIC ::                                                                      &
  signal_trap,                                                                 &
  signal_add_handler,                                                          &
  signal_command,                                                              &
  signal_rank,                                                                 &
  signal_controlled_exit,                                                      &
  signal_set_verbose,                                                          &
  signal_traceback,                                                            &
  signal_core,                                                                 &
  signal_off,                                                                  &
  signal_unregister_callbacks

!------------------------------------------------------------------------------!
! Signal handling method enumerator                                            !
!------------------------------------------------------------------------------!

ENUM, BIND(c)
ENUMERATOR ::                                                                  &
  signal_traceback,                                                            &
  signal_core,                                                                 &
  signal_off
END ENUM

!------------------------------------------------------------------------------!
! Callback handler abstract interface                                          !
!------------------------------------------------------------------------------!

ABSTRACT INTERFACE
SUBROUTINE handler_sub()

IMPLICIT NONE

END SUBROUTINE handler_sub
END INTERFACE

!------------------------------------------------------------------------------!
! C binding interfaces                                                         !
!------------------------------------------------------------------------------!

INTERFACE
SUBROUTINE signal_trap(option) BIND(c,NAME="signal_trap")

IMPORT :: C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN), VALUE :: option

END SUBROUTINE signal_trap
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE signal_rank(pe) BIND(c,NAME="signal_rank")

IMPORT :: C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN), VALUE :: pe

END SUBROUTINE signal_rank
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE signal_controlled_exit(e) BIND(c,NAME="signal_controlled_exit")

IMPORT :: C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN), VALUE :: e

END SUBROUTINE signal_controlled_exit
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE c_signal_add_handler(handler) BIND(c,NAME="signal_add_handler")

IMPORT :: C_FUNPTR

IMPLICIT NONE

TYPE(C_FUNPTR), INTENT(IN), VALUE :: handler

END SUBROUTINE c_signal_add_handler
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE signal_command(command, length) BIND(c,NAME="signal_command")

IMPORT :: C_INT64_T, C_CHAR

IMPLICIT NONE

CHARACTER(KIND=C_CHAR), INTENT(IN) :: command
INTEGER(KIND=C_INT64_T), INTENT(IN), VALUE :: length

END SUBROUTINE signal_command
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE signal_set_verbose() BIND(c,NAME="signal_set_verbose")

IMPLICIT NONE

END SUBROUTINE signal_set_verbose
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE signal_unregister_callbacks()                                       &
           BIND(c,NAME="signal_unregister_callbacks")

IMPLICIT NONE

END SUBROUTINE signal_unregister_callbacks
END INTERFACE

! -----------------------------------------------------------------------------!
CONTAINS
! -----------------------------------------------------------------------------!

SUBROUTINE signal_add_handler(handler)

  IMPLICIT NONE

  PROCEDURE(handler_sub) :: handler

  TYPE(C_FUNPTR) :: c_handler_ptr

  c_handler_ptr = C_FUNLOC(handler)

  CALL c_signal_add_handler(c_handler_ptr)

END SUBROUTINE signal_add_handler

! -----------------------------------------------------------------------------!

END MODULE fort2c_exceptions_interfaces

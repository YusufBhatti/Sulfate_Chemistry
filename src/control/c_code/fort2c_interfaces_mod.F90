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

MODULE fort2c_interfaces

! DEPENDS ON: portutils.o

USE, INTRINSIC :: ISO_C_BINDING, ONLY:                                         &
  C_BOOL,                                                                      &
  C_CHAR

IMPLICIT NONE

PRIVATE

PUBLIC ::                                                                      &
  get_file,                                                                    &
  um_sleep,                                                                    &
  logical_to_bool,                                                             &
  bool_to_logical

! -----------------------------------------------------------------------------!
! INTERFACE Blocks
! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE get_file (unit, filename, file_len, err) BIND(c,NAME="get_file")

IMPORT :: C_CHAR

IMPLICIT NONE

INTEGER ::  unit, file_len, err
CHARACTER(KIND=C_CHAR) :: filename(*)

END SUBROUTINE
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE um_sleep(sleep_time) BIND(c,NAME="um_sleep")
IMPLICIT NONE
INTEGER :: sleep_time
END SUBROUTINE
END INTERFACE

! -----------------------------------------------------------------------------!
CONTAINS
! -----------------------------------------------------------------------------!

ELEMENTAL FUNCTION bool_to_logical(cbool) RESULT(flog)
IMPLICIT NONE 
LOGICAL(KIND=C_BOOL), INTENT(IN) :: cbool
LOGICAL :: flog
flog = LOGICAL(cbool)
END FUNCTION bool_to_logical

! -----------------------------------------------------------------------------!

ELEMENTAL FUNCTION logical_to_bool(flog) RESULT(cbool)
IMPLICIT NONE 
LOGICAL, INTENT(IN) :: flog
LOGICAL(KIND=C_BOOL) :: cbool
cbool = LOGICAL(flog, KIND=C_BOOL)
END FUNCTION logical_to_bool

! -----------------------------------------------------------------------------!

END MODULE fort2c_interfaces


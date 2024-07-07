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

MODULE fort2c_getpos_interfaces

! DEPENDS ON: portio2a.o
! DEPENDS ON: portio2b.o

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT64_T

IMPLICIT NONE

PRIVATE

PUBLIC ::                                                                      &
  getpos,                                                                      &
  getpos8,                                                                     &
  getpos32,                                                                    &
  getpos64

! -----------------------------------------------------------------------------!
! INTERFACES                                                                   !
! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE c_getpos(unit, word_address) BIND(c,NAME="getpos")

IMPORT ::C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN) :: unit
INTEGER(KIND=C_INT64_T), INTENT(OUT) :: word_address

END SUBROUTINE
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE c_getpos8(unit, word_address) BIND(c,NAME="getpos8")

IMPORT ::C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN) :: unit
INTEGER(KIND=C_INT64_T), INTENT(OUT) :: word_address

END SUBROUTINE
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE c_getpos32(unit, word_address) BIND(c,NAME="getpos32")

IMPORT ::C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN) :: unit
INTEGER(KIND=C_INT64_T), INTENT(OUT) :: word_address

END SUBROUTINE
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE c_getpos64(unit, word_address) BIND(c,NAME="getpos64")

IMPORT ::C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN) :: unit
INTEGER(KIND=C_INT64_T), INTENT(OUT) :: word_address

END SUBROUTINE
END INTERFACE

! -----------------------------------------------------------------------------!

CONTAINS

! -----------------------------------------------------------------------------!
! Glue Routines                                                                !
! -----------------------------------------------------------------------------!

SUBROUTINE getpos(unit, word_address)

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit
INTEGER, INTENT(OUT) :: word_address

INTEGER(KIND=C_INT64_T) :: word_address_loc

CALL c_getpos(INT(unit,KIND=C_INT64_T), word_address_loc)

word_address = word_address_loc

END SUBROUTINE

! -----------------------------------------------------------------------------!

SUBROUTINE getpos8(unit, word_address)

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit
INTEGER, INTENT(OUT) :: word_address

INTEGER(KIND=C_INT64_T) :: word_address_loc

CALL c_getpos8(INT(unit,KIND=C_INT64_T), word_address_loc)

word_address = word_address_loc

END SUBROUTINE

! -----------------------------------------------------------------------------!

SUBROUTINE getpos32(unit, word_address)

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit
INTEGER, INTENT(OUT) :: word_address

INTEGER(KIND=C_INT64_T) :: word_address_loc

CALL c_getpos32(INT(unit,KIND=C_INT64_T), word_address_loc)

word_address = word_address_loc

END SUBROUTINE

! -----------------------------------------------------------------------------!

SUBROUTINE getpos64(unit, word_address)

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit
INTEGER, INTENT(OUT) :: word_address

INTEGER(KIND=C_INT64_T) :: word_address_loc

CALL c_getpos64(INT(unit,KIND=C_INT64_T), word_address_loc)

word_address = word_address_loc

END SUBROUTINE

! -----------------------------------------------------------------------------!

END MODULE

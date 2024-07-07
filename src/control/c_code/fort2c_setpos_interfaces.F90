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

MODULE fort2c_setpos_interfaces

! DEPENDS ON: portio2a.o
! DEPENDS ON: portio2b.o

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT64_T

IMPLICIT NONE

PRIVATE

PUBLIC ::                                                                      &
  setpos_single,                                                               &
  setpos8_single,                                                              &
  setpos32_single,                                                             &
  setpos64_single

! -----------------------------------------------------------------------------!
! INTERFACES                                                                   !
! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE c_setpos_single(unit, word_address, err)                            &
                           BIND(c,NAME="setpos_single")

IMPORT ::C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN)  :: unit, word_address
INTEGER(KIND=C_INT64_T), INTENT(OUT) :: err

END SUBROUTINE
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE c_setpos8_single(unit, word_address, err)                           &
                            BIND(c,NAME="setpos8_single")

IMPORT ::C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN)  :: unit, word_address
INTEGER(KIND=C_INT64_T), INTENT(OUT) :: err

END SUBROUTINE
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE c_setpos32_single(unit, word_address, err)                          &
                             BIND(c,NAME="setpos32_single")

IMPORT ::C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN)  :: unit, word_address
INTEGER(KIND=C_INT64_T), INTENT(OUT) :: err

END SUBROUTINE
END INTERFACE

! -----------------------------------------------------------------------------!

INTERFACE
SUBROUTINE c_setpos64_single(unit, word_address, err)                          &
                             BIND(c,NAME="setpos64_single")

IMPORT ::C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN)  :: unit, word_address
INTEGER(KIND=C_INT64_T), INTENT(OUT) :: err

END SUBROUTINE
END INTERFACE

! -----------------------------------------------------------------------------!

CONTAINS

! -----------------------------------------------------------------------------!
! Glue Routines                                                                !
! -----------------------------------------------------------------------------!

SUBROUTINE setpos_single(unit, word_address, err)

IMPLICIT NONE

INTEGER, INTENT(IN)  :: unit, word_address
INTEGER, INTENT(OUT) :: err

INTEGER(KIND=C_INT64_T) :: loc_err

CALL c_setpos_single(INT(unit,KIND=C_INT64_T),                                 &
                     INT(word_address,KIND=C_INT64_T),                         &
                     loc_err)

err = loc_err

END SUBROUTINE

! -----------------------------------------------------------------------------!

SUBROUTINE setpos8_single(unit, word_address, err)

IMPLICIT NONE

INTEGER, INTENT(IN)  :: unit, word_address
INTEGER, INTENT(OUT) :: err

INTEGER(KIND=C_INT64_T) :: loc_err

CALL c_setpos8_single(INT(unit,KIND=C_INT64_T),                                &
                      INT(word_address,KIND=C_INT64_T),                        &
                      loc_err)

err = loc_err

END SUBROUTINE

! -----------------------------------------------------------------------------!

SUBROUTINE setpos32_single(unit, word_address, err)

IMPLICIT NONE

INTEGER, INTENT(IN)  :: unit, word_address
INTEGER, INTENT(OUT) :: err

INTEGER(KIND=C_INT64_T) :: loc_err

CALL c_setpos32_single(INT(unit,KIND=C_INT64_T),                               &
                       INT(word_address,KIND=C_INT64_T),                       &
                       loc_err)

err = loc_err

END SUBROUTINE

! -----------------------------------------------------------------------------!

SUBROUTINE setpos64_single(unit, word_address, err)

IMPLICIT NONE

INTEGER, INTENT(IN)  :: unit, word_address
INTEGER, INTENT(OUT) :: err

INTEGER(KIND=C_INT64_T) :: loc_err

CALL c_setpos64_single(INT(unit,KIND=C_INT64_T),                               &
                       INT(word_address,KIND=C_INT64_T),                       &
                       loc_err)

err = loc_err

END SUBROUTINE

! -----------------------------------------------------------------------------!

END MODULE

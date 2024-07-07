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

MODULE fort2c_pio_byteswap_interfaces

! DEPENDS ON: pio_byteswap.o

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC, C_PTR, C_INT64_T

IMPLICIT NONE

PRIVATE

PUBLIC ::                                                                      &
  pio_byteswap,                                                                &
  endian_types,                                                                &
  bigEndian,                                                                   &
  littleEndian,                                                                &
  numEndians,                                                                  &
  get_machine_endianism

! -----------------------------------------------------------------------------!
! Endianism enumerator                                                         !
! -----------------------------------------------------------------------------!

ENUM, BIND(c)
ENUMERATOR ::                                                                  &
  bigEndian,                                                                   &
  littleEndian,                                                                &
  numEndians
END ENUM

INTEGER, PARAMETER :: endian_types = KIND(numEndians)

! -----------------------------------------------------------------------------!
! Interfaces to pio_byteswap                                                   !
! -----------------------------------------------------------------------------!

INTERFACE c_pio_byteswap

FUNCTION c_pio_byteswap(array, array_len, word_len)                            &
           BIND(c,NAME="pio_byteswap")

IMPORT :: C_INT64_T, C_PTR

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), VALUE :: array_len, word_len
INTEGER(KIND=C_INT64_T) :: c_pio_byteswap
TYPE(C_PTR), VALUE :: array

END FUNCTION c_pio_byteswap

END INTERFACE c_pio_byteswap

! -----------------------------------------------------------------------------!

INTERFACE pio_byteswap
MODULE PROCEDURE                                                               &
  pio_byteswap_64,                                                             &
  pio_byteswap_32
END INTERFACE pio_byteswap

! -----------------------------------------------------------------------------!
! Interfaces to get_machine_endianism                                          !
! -----------------------------------------------------------------------------!

INTERFACE get_machine_endianism

FUNCTION c_get_machine_endianism()                                             &
           BIND(c,NAME="get_machine_endianism")

IMPORT :: endian_types

IMPLICIT NONE

INTEGER(KIND=endian_types) :: c_get_machine_endianism

END FUNCTION c_get_machine_endianism

END INTERFACE get_machine_endianism

! -----------------------------------------------------------------------------!
CONTAINS
! -----------------------------------------------------------------------------!

FUNCTION pio_byteswap_64(array, array_len, word_len)

USE um_types, ONLY: integer64

IMPLICIT NONE

INTEGER(KIND=integer64), INTENT(IN) :: array_len, word_len
INTEGER(KIND=integer64) :: pio_byteswap_64
TYPE(C_PTR), INTENT(IN) :: array

pio_byteswap_64 = c_pio_byteswap(array,                                        &
                                 INT(array_len,kind=C_INT64_T),                &
                                 INT(word_len,kind=C_INT64_T))

END FUNCTION pio_byteswap_64

! -----------------------------------------------------------------------------!

FUNCTION pio_byteswap_32(array, array_len, word_len)

USE um_types, ONLY: integer32

IMPLICIT NONE

INTEGER(KIND=integer32), INTENT(IN) :: array_len, word_len
INTEGER(KIND=integer32) :: pio_byteswap_32
TYPE(C_PTR), INTENT(IN) :: array

pio_byteswap_32 = c_pio_byteswap(array,                                        &
                                 INT(array_len,kind=C_INT64_T),                &
                                 INT(word_len,kind=C_INT64_T))

END FUNCTION pio_byteswap_32

! -----------------------------------------------------------------------------!

END MODULE


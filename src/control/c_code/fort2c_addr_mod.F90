! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module contains the interfaces for interoperability of c and
! fortran addressing subroutines
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C Code
!

MODULE fort2c_addr_mod

! DEPENDS ON: c_address_routines.o

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INTPTR_T, C_SIZE_T

IMPLICIT NONE

PRIVATE

PUBLIC :: um_addr_diff, um_addr_size, um_addr_value

INTERFACE
SUBROUTINE um_addr_diff(object1,object2,diff) BIND(c,NAME="um_addr_diff")

IMPORT :: C_PTR, C_INTPTR_T

IMPLICIT NONE

TYPE(C_PTR), VALUE :: object1
TYPE(C_PTR), VALUE :: object2
INTEGER(KIND=C_INTPTR_T) :: diff

END SUBROUTINE
END INTERFACE

INTERFACE
FUNCTION um_addr_size() BIND(c,NAME="um_addr_size")

IMPORT :: C_SIZE_T

IMPLICIT NONE

INTEGER(KIND=C_SIZE_T) :: um_addr_size

END FUNCTION
END INTERFACE

! Return the integer value of the address of the given pointer.
INTERFACE
FUNCTION um_addr_value(object1) BIND(c,NAME="um_addr_value")

IMPORT :: C_PTR, C_INTPTR_T

IMPLICIT NONE

TYPE(C_PTR), VALUE       :: object1
INTEGER(KIND=C_INTPTR_T) :: um_addr_value  

END FUNCTION
END INTERFACE

END MODULE

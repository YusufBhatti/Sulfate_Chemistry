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

MODULE fort2c_portio_interfaces

#if defined(C95_2A)
! DEPENDS ON: portio2a.o
#elif defined(C95_2B)
! DEPENDS ON: portio2b.o
#endif

USE, INTRINSIC :: ISO_C_BINDING, ONLY:                                         &
  C_LOC,                                                                       &
  C_CHAR,                                                                      &
  C_PTR,                                                                       &
  C_INT64_T,                                                                   &
  C_INT8_T,                                                                    &
  C_BOOL

IMPLICIT NONE

PRIVATE

PUBLIC ::                                                                      &
  open_single,                                                                 &
  close_single,                                                                &
  sync_single,                                                                 &
  getextent,                                                                   &
  set_printstatus,                                                             &
  is_unit_open_local,                                                          &
  portioinit,                                                                  &
  portioshutdown,                                                              &
  portiodetachallhelpers,                                                      &
  portioaddhelper

!------------------------------------------------------------------------------!
! INTERFACE Blocks
!------------------------------------------------------------------------------!

INTERFACE open_single
SUBROUTINE c_open_single(unit, file_name, char_len, useintent,                 &
                         environ_var_flag, err) BIND(c,NAME="open_single")

IMPORT :: C_PTR, C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T) :: unit, char_len, useintent, environ_var_flag, err
TYPE(C_PTR), VALUE      :: file_name

END SUBROUTINE
MODULE PROCEDURE open_single_fort
END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE close_single
SUBROUTINE c_close_single(unit, file_name, char_len, environ_var_flag,         &
                          delete, err) BIND(c,NAME="close_single")

IMPORT :: C_PTR, C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T) :: unit, char_len, delete, environ_var_flag, err
TYPE(C_PTR), VALUE :: file_name

END SUBROUTINE
MODULE PROCEDURE close_single_fort
END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE
SUBROUTINE sync_single(unit, icode) BIND(c,NAME="sync_single")
IMPLICIT NONE
INTEGER :: unit, icode
END SUBROUTINE
END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE
SUBROUTINE getextent(unit, extent, icode) BIND(c,NAME="getextent")
IMPLICIT NONE
INTEGER :: unit, extent, icode
END SUBROUTINE
END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE
SUBROUTINE set_printstatus(printstatus_f) BIND(c,NAME="set_printstatus")
IMPLICIT NONE
INTEGER :: printstatus_f
END SUBROUTINE
END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE
SUBROUTINE c_is_unit_open_local(unit,ret_code) BIND(c,NAME="is_unit_open_local")

IMPORT :: C_BOOL, C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T) :: unit
LOGICAL(KIND=C_BOOL) :: ret_code

END SUBROUTINE
END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE
SUBROUTINE portioinit(wb_size, rb_size, rb_count, rb_update, rb_prefetch,      &
                      io_timing) BIND(c,NAME="portioinit")

IMPORT :: C_BOOL, C_INT64_T

IMPLICIT NONE

INTEGER(KIND=C_INT64_T) :: wb_size, rb_size, rb_count, rb_prefetch
LOGICAL(KIND=C_BOOL) :: rb_update, io_timing

END SUBROUTINE
END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE
SUBROUTINE portioshutdown() BIND(c,NAME="portioshutdown")

IMPLICIT NONE

END SUBROUTINE
END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE
SUBROUTINE portiodetachallhelpers() BIND(c,NAME="portiodetachallhelpers")

IMPLICIT NONE

END SUBROUTINE
END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE
SUBROUTINE portioaddhelper(role) BIND(c,NAME="portioaddhelper")

IMPLICIT NONE

INTEGER, INTENT(IN) :: role

END SUBROUTINE
END INTERFACE

!------------------------------------------------------------------------------!
CONTAINS
!------------------------------------------------------------------------------!

SUBROUTINE open_single_fort(unit, file_name, char_len, useintent,              &
                            environ_var_flag, err)

USE f_shum_string_conv_mod, ONLY: f_shum_f2c_string

IMPLICIT NONE

INTEGER :: unit, char_len, useintent, environ_var_flag, err
CHARACTER(LEN=*) :: file_name

INTEGER(KIND=C_INT64_T) :: err_loc
CHARACTER(KIND=C_CHAR,LEN=1), TARGET :: c_file_name(LEN_TRIM(file_name)+1)


c_file_name = f_shum_f2c_string(file_name)

CALL c_open_single(INT(unit,KIND=C_INT64_T), C_LOC(c_file_name),               &
                   INT(char_len,KIND=C_INT64_T),                               &
                   INT(useintent,KIND=C_INT64_T),                              &
                   INT(environ_var_flag,KIND=C_INT64_T), err_loc)


err = err_loc

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE close_single_fort(unit, file_name, char_len, environ_var_flag,      &
                             delete, err)



USE f_shum_string_conv_mod, ONLY: f_shum_f2c_string

IMPLICIT NONE

INTEGER :: unit, char_len, delete, environ_var_flag, err
CHARACTER(LEN=*) :: file_name

INTEGER(KIND=C_INT64_T) :: err_loc
CHARACTER(KIND=C_CHAR,LEN=1), TARGET :: c_file_name(LEN_TRIM(file_name)+1)


c_file_name = f_shum_f2c_string(file_name)

CALL c_close_single(INT(unit,KIND=C_INT64_T), C_LOC(c_file_name),              &
                    INT(char_len,KIND=C_INT64_T),                              &
                    INT(environ_var_flag,KIND=C_INT64_T),                      &
                    INT(delete,KIND=C_INT64_T), err_loc)


err = err_loc

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE is_unit_open_local(unit,open_status)

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit
LOGICAL, INTENT(OUT) :: open_status

LOGICAL(c_bool) :: c_open_status

CALL c_is_unit_open_local(INT(unit,KIND=C_INT64_T),c_open_status)

open_status = c_open_status

END SUBROUTINE

!------------------------------------------------------------------------------!

END MODULE fort2c_portio_interfaces

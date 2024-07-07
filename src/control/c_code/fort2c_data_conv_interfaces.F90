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

MODULE fort2c_data_conv_interfaces

! DEPENDS ON: pio_data_conv.o

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC, C_PTR, C_INT64_T
USE um_types, ONLY: integer64, integer32, real64, real32

IMPLICIT NONE

PRIVATE

PUBLIC ::                                                                      &
  ibm2ieee,                                                                    &
  ieee2ibm,                                                                    &
  ieee2ieg,                                                                    &
  ieg2ieee,                                                                    &
  data_types,                                                                  &
  typenull,                                                                    &
  typeless,                                                                    &
  integer_type,                                                                &
  real_type,                                                                   &
  complex_type,                                                                &
  logical_type,                                                                &
  character_type

! -----------------------------------------------------------------------------!

ENUM, BIND(c)
ENUMERATOR ::                                                                  &
  typenull,                                                                    &
  typeless,                                                                    &
  integer_type,                                                                &
  real_type,                                                                   &
  complex_type,                                                                &
  logical_type,                                                                &
  character_type
END ENUM

INTEGER, PARAMETER :: data_types = KIND(typenull)

! -----------------------------------------------------------------------------!
! Interfaces

! C Interfaces

INTERFACE
FUNCTION c_ieg2ieee (data_type, num, ieg_num_in, offset_in, cri_num_out,       &
                     stride, size_num_out, size_num_in) BIND(c,NAME="ieg2ieee")

IMPORT :: data_types, c_int64_t, c_ptr

IMPLICIT NONE

INTEGER(KIND=c_int64_t) ::                                                     &
  c_ieg2ieee

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=c_int64_t), INTENT(IN) ::                                         &
  num,                                                                         &
  offset_in,                                                                   &
  stride,                                                                      &
  size_num_out,                                                                &
  size_num_in

TYPE(c_ptr), INTENT(IN), VALUE ::                                              &
  ieg_num_in,                                                                  &
  cri_num_out

END FUNCTION
END INTERFACE

INTERFACE
FUNCTION c_ieee2ieg (data_type, num, ieg_num_out, offset_out, cri_num_in,      &
                     stride, size_num_in, size_num_out) BIND(c,NAME="ieee2ieg")

IMPORT :: data_types, c_int64_t, c_ptr

IMPLICIT NONE

INTEGER(KIND=c_int64_t) ::                                                     &
  c_ieee2ieg

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=c_int64_t), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

TYPE(c_ptr), INTENT(IN), VALUE ::                                              &
  cri_num_in,                                                                  &
  ieg_num_out

END FUNCTION
END INTERFACE

INTERFACE
FUNCTION c_ieee2ibm (data_type, num, ibm_num_out, offset_out, cri_num_in,      &
                     stride, size_num_in, size_num_out) BIND(c,NAME="ieee2ibm")

IMPORT :: data_types, C_INT64_T, c_ptr

IMPLICIT NONE

INTEGER(KIND=C_INT64_T) ::                                                     &
  c_ieee2ibm

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=C_INT64_T), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

TYPE(c_ptr), INTENT(IN), VALUE ::                                              &
  cri_num_in,                                                                  &
  ibm_num_out

END FUNCTION
END INTERFACE

INTERFACE
FUNCTION c_ibm2ieee (data_type, num, ibm_num_in, offset_in, cri_num_out,       &
                     stride, size_num_out, size_num_in) BIND(c,NAME="ibm2ieee")

IMPORT :: data_types, C_INT64_T, c_ptr

IMPLICIT NONE

INTEGER(KIND=C_INT64_T) ::                                                     &
  c_ibm2ieee

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=C_INT64_T), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_in,                                                                   &
  size_num_in,                                                                 &
  size_num_out

TYPE(c_ptr), INTENT(IN), VALUE ::                                              &
  ibm_num_in,                                                                  &
  cri_num_out

END FUNCTION
END INTERFACE

! Generic (Fortran) Interfaces

INTERFACE ieee2ieg
MODULE PROCEDURE                                                               &
  ieee2ieg_real,                                                               &
  ieee2ieg_real_scalar,                                                        &
  ieee2ieg_real_to_int,                                                        &
  ieee2ieg_int_to_real,                                                        &
  ieee2ieg_int_to_real_scalar,                                                 &
  ieee2ieg_64

END INTERFACE

INTERFACE ieee2ibm
MODULE PROCEDURE                                                               &
  ieee2ibm_scalar,                                                             &
  ieee2ibm_int_to_real_scalar,                                                 &
  ieee2ibm_int_to_real32_scalar,                                               &
  ieee2ibm_real_to_int,                                                        &
  ieee2ibm_64,                                                                 &
  ieee2ibm_real,                                                               &
  ieee2ibm_real_scalar,                                                        &
  ieee2ibm_32_scalar_32_to_64
END INTERFACE

INTERFACE ibm2ieee
MODULE PROCEDURE                                                               &
  ibm2ieee_scalar,                                                             &
  ibm2ieee_real_to_int_scalar,                                                 &
  ibm2ieee_real32_to_int_scalar,                                               &
  ibm2ieee_int_to_real,                                                        &
  ibm2ieee_real32_real64,                                                      &
  ibm2ieee_real32_to_real64_2d,                                                &
  ibm2ieee_64,                                                                 &
  ibm2ieee_32_scalar_64_to_32
END INTERFACE

INTERFACE ieg2ieee
MODULE PROCEDURE                                                               &
  ieg2ieee_64
END INTERFACE

CONTAINS

! -----------------------------------------------------------------------------!
! ieg2ieee

FUNCTION ieg2ieee_64 (data_type, num, ieg_num_in, offset_in, cri_num_out,      &
                      stride, size_num_out, size_num_in)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieg2ieee_64

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_in,                                                                   &
  size_num_in,                                                                 &
  size_num_out

INTEGER(kind=integer64), INTENT(INOUT), TARGET ::                              &
  cri_num_out(num)

INTEGER(kind=integer64), INTENT(IN), TARGET ::                                 &
  ieg_num_in(CEILING(num*REAL(size_num_in)/size_num_out))

ieg2ieee_64 = c_ieg2ieee(data_type, num, C_LOC(ieg_num_in), offset_in,         &
                         C_LOC(cri_num_out), stride, size_num_out, size_num_in)

END FUNCTION

!------------------------------------------------------------------------------!
! ieee2ieg

FUNCTION ieee2ieg_64 (data_type, num, ieg_num_out, offset_out,  cri_num_in,    &
                      stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ieg_64

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

INTEGER(kind=integer64), INTENT(INOUT), TARGET ::                              &
  ieg_num_out(num)

INTEGER(kind=integer64), INTENT(IN), TARGET ::                                 &
  cri_num_in(CEILING(num*REAL(size_num_in)/size_num_out))

ieee2ieg_64 = c_ieee2ieg(data_type, num, C_LOC(ieg_num_out), offset_out,       &
                         C_LOC(cri_num_in), stride, size_num_in, size_num_out)

END FUNCTION

FUNCTION ieee2ieg_real (data_type, num, ieg_num_out, offset_out, cri_num_in,   &
                        stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ieg_real

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

REAL(kind=real64), INTENT(IN), TARGET ::                                       &
  cri_num_in(CEILING(num*REAL(size_num_in)/size_num_out))

REAL(kind=real64), INTENT(INOUT), TARGET ::                                    &
  ieg_num_out(num)

ieee2ieg_real = c_ieee2ieg(data_type, num, C_LOC(ieg_num_out), offset_out,     &
                           C_LOC(cri_num_in), stride, size_num_in, size_num_out)

END FUNCTION

FUNCTION ieee2ieg_int_to_real (data_type, num, ieg_num_out, offset_out,        &
                               cri_num_in, stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ieg_int_to_real

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

REAL(kind=real64), INTENT(INOUT), TARGET ::                                    &
  ieg_num_out(num)

INTEGER(kind=integer64), INTENT(IN), TARGET ::                                 &
  cri_num_in(CEILING(num*REAL(size_num_in)/size_num_out))


ieee2ieg_int_to_real = c_ieee2ieg(data_type, num, C_LOC(ieg_num_out),          &
                                  offset_out, C_LOC(cri_num_in), stride,       &
                                  size_num_in, size_num_out)

END FUNCTION

FUNCTION ieee2ieg_real_to_int (data_type, num, ieg_num_out, offset_out,        &
                               cri_num_in, stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ieg_real_to_int

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

REAL(kind=real64), INTENT(IN), TARGET ::                                       &
  cri_num_in(CEILING(num*REAL(size_num_in)/size_num_out))

INTEGER(kind=integer64), INTENT(INOUT), TARGET ::                              &
  ieg_num_out(num)


ieee2ieg_real_to_int = c_ieee2ieg(data_type, num, C_LOC(ieg_num_out),          &
                                  offset_out, C_LOC(cri_num_in), stride,       &
                                  size_num_in, size_num_out)

END FUNCTION

FUNCTION ieee2ieg_real_scalar (data_type, num, ieg_num_out, offset_out,        &
                               cri_num_in,  stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ieg_real_scalar

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

REAL(kind=real64), INTENT(IN), TARGET ::                                       &
  cri_num_in

REAL(kind=real64), INTENT(INOUT), TARGET ::                                    &
  ieg_num_out

ieee2ieg_real_scalar = c_ieee2ieg(data_type, num, C_LOC(ieg_num_out),          &
                                  offset_out, C_LOC(cri_num_in), stride,       &
                                  size_num_in, size_num_out)

END FUNCTION

FUNCTION ieee2ieg_int_to_real_scalar(data_type, num, ieg_num_out, offset_out,  &
                                 cri_num_in,  stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ieg_int_to_real_scalar

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

INTEGER(kind=integer64), INTENT(IN), TARGET ::                                 &
  cri_num_in

REAL(kind=real64), INTENT(INOUT), TARGET ::                                    &
  ieg_num_out

ieee2ieg_int_to_real_scalar = c_ieee2ieg(data_type, num, C_LOC(ieg_num_out),   &
                                         offset_out, C_LOC(cri_num_in), stride,&
                                         size_num_in, size_num_out)

END FUNCTION

! -----------------------------------------------------------------------------!
! ieee2ibm

FUNCTION ieee2ibm_64 (data_type, num, ibm_num_out, offset_out, cri_num_in,     &
                      stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ibm_64

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

INTEGER(kind=integer64), INTENT(INOUT), TARGET ::                              &
  ibm_num_out(num)

INTEGER(kind=integer64), INTENT(IN), TARGET ::                                 &
  cri_num_in(CEILING(num*REAL(size_num_in)/size_num_out))

ieee2ibm_64 = c_ieee2ibm(data_type, num, C_LOC(ibm_num_out), offset_out,       &
                         C_LOC(cri_num_in), stride, size_num_in, size_num_out)

END FUNCTION


FUNCTION ieee2ibm_real_scalar (data_type, num, ibm_num_out, offset_out,        &
                               cri_num_in, stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ibm_real_scalar

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

REAL(kind=real64), INTENT(INOUT), TARGET ::                                    &
  ibm_num_out

REAL(kind=real64), INTENT(IN), TARGET ::                                       &
  cri_num_in

ieee2ibm_real_scalar = c_ieee2ibm(data_type, num, C_LOC(ibm_num_out),          &
                                  offset_out, C_LOC(cri_num_in), stride,       &
                                  size_num_in, size_num_out)

END FUNCTION

FUNCTION ieee2ibm_real (data_type, num, ibm_num_out, offset_out, cri_num_in,   &
                        stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ibm_real

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num, stride, offset_out, size_num_in, size_num_out

REAL(kind=real64), INTENT(INOUT), TARGET ::                                    &
  ibm_num_out(num)

REAL(kind=real64), INTENT(IN), TARGET ::                                       &
  cri_num_in(CEILING(num*REAL(size_num_in)/size_num_out))

ieee2ibm_real = c_ieee2ibm(data_type, num, C_LOC(ibm_num_out), offset_out,     &
                           C_LOC(cri_num_in), stride, size_num_in, size_num_out)

END FUNCTION

FUNCTION ieee2ibm_real_to_int (data_type, num, ibm_num_out, offset_out,        &
                               cri_num_in, stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ibm_real_to_int

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

INTEGER(kind=integer64), INTENT(INOUT), TARGET ::                              &
  ibm_num_out(num)

REAL(kind=real64), INTENT(IN), TARGET ::                                       &
  cri_num_in(CEILING(num*REAL(size_num_in)/size_num_out))

ieee2ibm_real_to_int = c_ieee2ibm(data_type, num, C_LOC(ibm_num_out),          &
                                  offset_out, C_LOC(cri_num_in), stride,       &
                                  size_num_in, size_num_out)

END FUNCTION

FUNCTION ieee2ibm_scalar (data_type, num, ibm_num_out, offset_out, cri_num_in, &
                          stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ibm_scalar

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

INTEGER(kind=integer64), INTENT(IN), TARGET ::                                 &
  cri_num_in

INTEGER(kind=integer64), INTENT(INOUT), TARGET ::                              &
  ibm_num_out


ieee2ibm_scalar = c_ieee2ibm(data_type, num, C_LOC(ibm_num_out), offset_out,   &
                             C_LOC(cri_num_in), stride, size_num_in,           &
                             size_num_out)

END FUNCTION

FUNCTION ieee2ibm_32_scalar_32_to_64 (data_type, num, ibm_num_out, offset_out, &
                            cri_num_in, stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer32) ::                                                     &
  ieee2ibm_32_scalar_32_to_64

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer32), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

INTEGER(kind=integer32), INTENT(IN), TARGET ::                                 &
  cri_num_in

INTEGER(kind=integer64), INTENT(INOUT), TARGET ::                              &
  ibm_num_out

ieee2ibm_32_scalar_32_to_64 = c_ieee2ibm(data_type, INT(num,kind=C_INT64_T),   &
                                         C_LOC(ibm_num_out),                   &
                                         INT(offset_out,kind=C_INT64_T),       &
                                         C_LOC(cri_num_in),                    &
                                         INT(stride,kind=C_INT64_T),           &
                                         INT(size_num_in,kind=C_INT64_T),      &
                                         INT(size_num_out,kind=C_INT64_T))

END FUNCTION

FUNCTION ieee2ibm_int_to_real_scalar (data_type, num, ibm_num_out, offset_out, &
                                  cri_num_in, stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(kind=integer64) ::                                                     &
  ieee2ibm_int_to_real_scalar

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(kind=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

INTEGER(kind=integer64), INTENT(IN), TARGET ::                                 &
  cri_num_in

REAL(kind=real64), INTENT(INOUT), TARGET ::                                    &
  ibm_num_out

ieee2ibm_int_to_real_scalar = c_ieee2ibm(data_type, num, C_LOC(ibm_num_out),   &
                                         offset_out, C_LOC(cri_num_in), stride,&
                                         size_num_in, size_num_out)

END FUNCTION

FUNCTION ieee2ibm_int_to_real32_scalar (data_type, num, ibm_num_out,           &
     offset_out,cri_num_in, stride, size_num_in, size_num_out)

IMPLICIT NONE

INTEGER(KIND=integer64) ::                                                     &
  ieee2ibm_int_to_real32_scalar

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_out,                                                                  &
  size_num_in,                                                                 &
  size_num_out

INTEGER(KIND=integer64), INTENT(IN), TARGET ::                                 &
  cri_num_in

REAL(KIND=real32), INTENT(INOUT), TARGET ::                                    &
  ibm_num_out

ieee2ibm_int_to_real32_scalar = c_ieee2ibm(data_type, num, C_LOC(ibm_num_out), &
                                         offset_out, C_LOC(cri_num_in), stride,&
                                         size_num_in, size_num_out)

END FUNCTION

! -----------------------------------------------------------------------------!
! ibm2ieee

FUNCTION ibm2ieee_64 (data_type, num, ibm_num_in, offset_in, cri_num_out,      &
                      stride, size_num_out, size_num_in)

IMPLICIT NONE

INTEGER(KIND=integer64) ::                                                     &
  ibm2ieee_64

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_in,                                                                   &
  size_num_in,                                                                 &
  size_num_out

INTEGER(KIND=integer64), INTENT(IN), TARGET ::                                 &
  ibm_num_in(CEILING(num*REAL(size_num_in)/size_num_out))

INTEGER(KIND=integer64), INTENT(INOUT), TARGET ::                              &
  cri_num_out(num)

ibm2ieee_64 = c_ibm2ieee(data_type, num, C_LOC(ibm_num_in), offset_in,         &
                         C_LOC(cri_num_out), stride, size_num_out, size_num_in)

END FUNCTION

FUNCTION ibm2ieee_scalar (data_type, num, ibm_num_in, offset_in, cri_num_out,  &
                          stride, size_num_out, size_num_in)

IMPLICIT NONE

INTEGER(KIND=integer64) ::                                                     &
  ibm2ieee_scalar

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_in,                                                                   &
  size_num_in,                                                                 &
  size_num_out

INTEGER(KIND=integer64), INTENT(IN), TARGET ::                                 &
  ibm_num_in

INTEGER(KIND=integer64), INTENT(INOUT), TARGET ::                              &
  cri_num_out

ibm2ieee_scalar = c_ibm2ieee(data_type, num, C_LOC(ibm_num_in), offset_in,     &
                             C_LOC(cri_num_out), stride, size_num_out,         &
                             size_num_in)

END FUNCTION

FUNCTION ibm2ieee_32_scalar_64_to_32 (data_type, num, ibm_num_in, offset_in,   &
                                 cri_num_out, stride, size_num_out, size_num_in)

IMPLICIT NONE

INTEGER(KIND=integer32) ::                                                     &
  ibm2ieee_32_scalar_64_to_32

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=integer32), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_in,                                                                   &
  size_num_in,                                                                 &
  size_num_out

INTEGER(kind=integer64), INTENT(IN), TARGET ::                                 &
  ibm_num_in

INTEGER(kind=integer32), INTENT(INOUT), TARGET ::                              &
  cri_num_out

ibm2ieee_32_scalar_64_to_32 = c_ibm2ieee(data_type, INT(num,kind=integer64),   &
                                         C_LOC(ibm_num_in),                    &
                                         INT(offset_in,kind=integer64),        &
                                         C_LOC(cri_num_out),                   &
                                         INT(stride,kind=integer64),           &
                                         INT(size_num_out,kind=integer64),     &
                                         INT(size_num_in,kind=integer64))

END FUNCTION

FUNCTION ibm2ieee_real32_real64 (data_type, num, ibm_num_in, offset_in,        &
                                 cri_num_out, stride, size_num_out, size_num_in)

IMPLICIT NONE

INTEGER(KIND=integer64) ::                                                     &
  ibm2ieee_real32_real64

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_in,                                                                   &
  size_num_in,                                                                 &
  size_num_out

REAL(kind=real32), INTENT(IN), TARGET ::                                       &
  ibm_num_in(CEILING(num*size_num_in*2.0/size_num_out))

REAL(kind=real64), INTENT(INOUT), TARGET ::                                    &
  cri_num_out(num)

ibm2ieee_real32_real64 = c_ibm2ieee(data_type, num, C_LOC(ibm_num_in),         &
                                    offset_in, C_LOC(cri_num_out), stride,     &
                                    size_num_out, size_num_in)

END FUNCTION

FUNCTION ibm2ieee_real_to_int_scalar (data_type, num, ibm_num_in, offset_in,   &
                                 cri_num_out, stride, size_num_out, size_num_in)

IMPLICIT NONE

INTEGER(KIND=integer64) ::                                                     &
  ibm2ieee_real_to_int_scalar

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_in,                                                                   &
  size_num_in,                                                                 &
  size_num_out

INTEGER(KIND=integer64), INTENT(INOUT), TARGET ::                              &
  cri_num_out

REAL(KIND=real64), INTENT(IN), TARGET ::                                       &
  ibm_num_in

ibm2ieee_real_to_int_scalar = c_ibm2ieee(data_type, num, C_LOC(ibm_num_in),    &
                                         offset_in, C_LOC(cri_num_out), stride,&
                                         size_num_out, size_num_in)

END FUNCTION

FUNCTION ibm2ieee_real32_to_int_scalar (data_type, num, ibm_num_in, offset_in, &
                                 cri_num_out, stride, size_num_out, size_num_in)

IMPLICIT NONE

INTEGER(KIND=integer64) ::                                                     &
  ibm2ieee_real32_to_int_scalar

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_in,                                                                   &
  size_num_in,                                                                 &
  size_num_out

INTEGER(KIND=integer64), INTENT(INOUT), TARGET ::                              &
  cri_num_out

REAL(KIND=real32), INTENT(IN), TARGET ::                                       &
  ibm_num_in

ibm2ieee_real32_to_int_scalar = c_ibm2ieee(data_type, num, C_LOC(ibm_num_in),  &
                                         offset_in, C_LOC(cri_num_out), stride,&
                                         size_num_out, size_num_in)

END FUNCTION

FUNCTION ibm2ieee_real32_to_real64_2d (data_type, num, ibm_num_in, offset_in,  &
                                 cri_num_out, stride, size_num_out, size_num_in)

IMPLICIT NONE

INTEGER(KIND=integer64) ::                                                     &
  ibm2ieee_real32_to_real64_2d

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_in,                                                                   &
  size_num_in,                                                                 &
  size_num_out

REAL(kind=real64), INTENT(INOUT), TARGET ::                                    &
  cri_num_out(:,:)

REAL(kind=real32), INTENT(IN), TARGET ::                                       &
  ibm_num_in(CEILING(num*size_num_in*2.0/size_num_out))

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Thererefore use address of (1,1)th element of the array

!NB: this means the user must ensure the array is contiguos in memory
!    and stride one

ibm2ieee_real32_to_real64_2d = c_ibm2ieee(data_type, num, C_LOC(ibm_num_in),   &
                                          offset_in, C_LOC(cri_num_out(1,1)),  &
                                          stride, size_num_out, size_num_in)

END FUNCTION

FUNCTION ibm2ieee_int_to_real (data_type, num, ibm_num_in, offset_in,          &
                               cri_num_out, stride, size_num_out, size_num_in)

IMPLICIT NONE

INTEGER(KIND=integer64) ::                                                     &
  ibm2ieee_int_to_real

INTEGER(KIND=data_types), INTENT(IN) ::                                        &
  data_type

INTEGER(KIND=integer64), INTENT(IN) ::                                         &
  num,                                                                         &
  stride,                                                                      &
  offset_in,                                                                   &
  size_num_in,                                                                 &
  size_num_out

REAL(KIND=real64), INTENT(INOUT), TARGET ::                                    &
  cri_num_out(num)

INTEGER(KIND=integer64), INTENT(IN), TARGET ::                                 &
  ibm_num_in(CEILING(num*REAL(size_num_in)/size_num_out))


ibm2ieee_int_to_real = c_ibm2ieee(data_type, num, C_LOC(ibm_num_in), offset_in,&
                                  C_LOC(cri_num_out), stride, size_num_out,    &
                                  size_num_in)

END FUNCTION

END MODULE


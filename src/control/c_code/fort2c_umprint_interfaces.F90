! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module contains the interfaces for interoperability of umPrint()
! between C and Fortran.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C Code
!

MODULE fort2c_umprint_interfaces

USE umprintmgr, ONLY:                                                          &
  umPrint,                                                                     &
  maxLineLen

USE f_shum_string_conv_mod, ONLY:                                              &
  f_shum_c2f_string

USE, INTRINSIC :: ISO_C_BINDING, ONLY:                                         &
  C_INT

IMPLICIT NONE

PRIVATE

! The following PARAMETERS are used by the C wrapper for umPrint() in order to
! encode information about the arguments it was called with. (They function as
! bit masks).

! They are not intended to be used directly.

INTEGER(KIND=c_int), PARAMETER ::                                              &
  one_bit = INT(b'1')

INTEGER(KIND=c_int), PARAMETER ::                                              &
  level_set      = ISHFT(one_bit, 0),                                          &
  pe_set         = ISHFT(one_bit, 1),                                          &
  model_set      = ISHFT(one_bit, 2),                                          &
  UsrPrefix_set  = ISHFT(one_bit, 3),                                          &
  HangIndent_set = ISHFT(one_bit, 4)

! ---------------------------------------------------------------------------- !
CONTAINS
! ---------------------------------------------------------------------------- !

! ---------------------------------------------------------------------------- !
! Functions to expose Fortran PARAMETERS to C code                             !
! ---------------------------------------------------------------------------- !

FUNCTION c_maxLineLen() BIND(c,NAME="c_maxLineLen")

  IMPLICIT NONE

  INTEGER(KIND=c_int) :: c_maxLineLen

  c_maxLineLen = maxLineLen

END FUNCTION

! ---------------------------------------------------------------------------- !

FUNCTION c_level_set() BIND(C,NAME="level_set")

  IMPLICIT NONE

  INTEGER(KIND=c_int) :: c_level_set

  c_level_set = level_set

END FUNCTION

! ---------------------------------------------------------------------------- !

FUNCTION c_pe_set() BIND(C,NAME="pe_set")

  IMPLICIT NONE

  INTEGER(KIND=c_int) :: c_pe_set

  c_pe_set = pe_set

END FUNCTION

! ---------------------------------------------------------------------------- !

FUNCTION c_model_set() BIND(C,NAME="model_set")

  IMPLICIT NONE

  INTEGER(KIND=c_int) :: c_model_set

  c_model_set = model_set

END FUNCTION

! ---------------------------------------------------------------------------- !

FUNCTION c_UsrPrefix_set() BIND(C,NAME="UsrPrefix_set")

  IMPLICIT NONE

  INTEGER(KIND=c_int) :: c_UsrPrefix_set

  c_UsrPrefix_set = UsrPrefix_set

END FUNCTION

! ---------------------------------------------------------------------------- !

FUNCTION c_HangIndent_set() BIND(C,NAME="HangIndent_set")

  IMPLICIT NONE

  INTEGER(KIND=c_int) :: c_HangIndent_set

  c_HangIndent_set = HangIndent_set

END FUNCTION

! ---------------------------------------------------------------------------- !
! umPrint Interface routines                                                   !
! ---------------------------------------------------------------------------- !

SUBROUTINE f_umPrint(message, src, level, pe, model, UsrPrefix, HangIndent,    &
                     stdErrorToo, setsum)                                      &
           BIND(c,NAME="f_umPrint")

USE, INTRINSIC :: ISO_C_BINDING, ONLY:                                         &
  C_CHAR, C_BOOL

IMPLICIT NONE

CHARACTER(KIND=C_CHAR,LEN=1), INTENT(IN) ::                                    &
  message(*),                                                                  &
  src(*),                                                                      &
  model(*),                                                                    &
  UsrPrefix(*),                                                                &
  HangIndent(*)

INTEGER(KIND=C_INT), INTENT(IN) ::                                             &
  level,                                                                       &
  pe,                                                                          &
  setsum

LOGICAL(KIND=C_BOOL), INTENT(IN) ::                                            &
  stdErrorToo

IF (IAND(setsum,level_set)==0) THEN

  CALL f_umPrint_opt_lev(line=TRIM(f_shum_c2f_string(message)),                &
                         src=TRIM(f_shum_c2f_string(src)),                     &
                         pe=pe,                                                &
                         model=model,                                          &
                         UsrPrefix=UsrPrefix,                                  &
                         HangIndent=HangIndent,                                &
                         stdErrorToo=LOGICAL(stdErrorToo),                     &
                         setsum=setsum)

ELSE

  CALL f_umPrint_opt_lev(line=TRIM(f_shum_c2f_string(message)),                &
                         src=TRIM(f_shum_c2f_string(src)),                     &
                         level=INT(level),                                     &
                         pe=pe,                                                &
                         model=model,                                          &
                         UsrPrefix=UsrPrefix,                                  &
                         HangIndent=HangIndent,                                &
                         stdErrorToo=LOGICAL(stdErrorToo),                     &
                         setsum=setsum)

END IF

END SUBROUTINE f_umPrint

! ---------------------------------------------------------------------------- !

SUBROUTINE f_umPrint_opt_lev(line, src, level, pe, model, UsrPrefix,           &
                             HangIndent, stdErrorToo, setsum)

USE, INTRINSIC :: ISO_C_BINDING, ONLY:                                         &
  C_CHAR

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) ::                                                &
  line,                                                                        &
  src

INTEGER, INTENT(IN), OPTIONAL ::                                               &
  level

INTEGER(KIND=C_INT), INTENT(IN), OPTIONAL ::                                   &
  pe

CHARACTER(KIND=C_CHAR,LEN=1), INTENT(IN), OPTIONAL ::                          &
  model(*),                                                                    &
  UsrPrefix(*),                                                                &
  HangIndent(*)

LOGICAL, INTENT(IN) ::                                                         &
  stdErrorToo

INTEGER(KIND=C_INT), INTENT(IN) ::                                             &
  setsum

IF (IAND(setsum,pe_set)==0) THEN

  CALL f_umPrint_opt_pe(line=line,                                             &
                        src=src,                                               &
                        level=level,                                           &
                        model=model,                                           &
                        UsrPrefix=UsrPrefix,                                   &
                        HangIndent=HangIndent,                                 &
                        stdErrorToo=stdErrorToo,                               &
                        setsum=setsum)

ELSE

  CALL f_umPrint_opt_pe(line=line,                                             &
                        src=src,                                               &
                        level=level,                                           &
                        pe=INT(pe),                                            &
                        model=model,                                           &
                        UsrPrefix=UsrPrefix,                                   &
                        HangIndent=HangIndent,                                 &
                        stdErrorToo=stdErrorToo,                               &
                        setsum=setsum)

END IF

END SUBROUTINE f_umPrint_opt_lev

! ---------------------------------------------------------------------------- !

SUBROUTINE f_umPrint_opt_pe(line, src, level, pe, model, UsrPrefix,            &
                            HangIndent, stdErrorToo, setsum)

USE, INTRINSIC :: ISO_C_BINDING, ONLY:                                         &
  C_CHAR

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) ::                                                &
  line,                                                                        &
  src

INTEGER, OPTIONAL, INTENT(IN) ::                                               &
  level,                                                                       &
  pe

CHARACTER(KIND=C_CHAR,LEN=1), INTENT(IN), OPTIONAL ::                          &
  model(*),                                                                    &
  UsrPrefix(*),                                                                &
  HangIndent(*)

LOGICAL, INTENT(IN) ::                                                         &
  stdErrorToo

INTEGER(KIND=C_INT), INTENT(IN) ::                                             &
  setsum

IF (IAND(setsum,model_set)==0) THEN

  CALL f_umPrint_opt_model(line=line,                                          &
                           src=src,                                            &
                           level=level,                                        &
                           pe=pe,                                              &
                           UsrPrefix=UsrPrefix,                                &
                           HangIndent=HangIndent,                              &
                           stdErrorToo=stdErrorToo,                            &
                           setsum=setsum)

ELSE

  CALL f_umPrint_opt_model(line=line,                                          &
                           src=src,                                            &
                           level=level,                                        &
                           pe=pe,                                              &
                           model=TRIM(f_shum_c2f_string(model)),               &
                           UsrPrefix=UsrPrefix,                                &
                           HangIndent=HangIndent,                              &
                           stdErrorToo=stdErrorToo,                            &
                           setsum=setsum)

END IF

END SUBROUTINE f_umPrint_opt_pe

! ---------------------------------------------------------------------------- !

SUBROUTINE f_umPrint_opt_model(line, src, level, pe, model, UsrPrefix,         &
                               HangIndent, stdErrorToo, setsum)

USE, INTRINSIC :: ISO_C_BINDING, ONLY:                                         &
  C_CHAR

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) ::                                                &
  line,                                                                        &
  src

INTEGER, OPTIONAL, INTENT(IN) ::                                               &
  level,                                                                       &
  pe

CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  ::                                     &
  model

CHARACTER(KIND=C_CHAR,LEN=1), INTENT(IN), OPTIONAL ::                          &
  UsrPrefix(*),                                                                &
  HangIndent(*)

LOGICAL, INTENT(IN) ::                                                         &
  stdErrorToo

INTEGER(KIND=C_INT), INTENT(IN) ::                                             &
  setsum

IF (IAND(setsum,UsrPrefix_set)==0) THEN

  CALL f_umPrint_opt_UsrPrefix(line=line,                                      &
                               src=src,                                        &
                               level=level,                                    &
                               pe=pe,                                          &
                               model=model,                                    &
                               HangIndent=HangIndent,                          &
                               stdErrorToo=stdErrorToo,                        &
                               setsum=setsum)

ELSE

  CALL f_umPrint_opt_UsrPrefix(line=line,                                      &
                               src=src,                                        &
                               level=level,                                    &
                               pe=pe,                                          &
                               model=model,                                    &
                               UsrPrefix=TRIM(f_shum_c2f_string(UsrPrefix)),   &
                               HangIndent=HangIndent,                          &
                               stdErrorToo=stdErrorToo,                        &
                               setsum=setsum)

END IF

END SUBROUTINE f_umPrint_opt_model

! ---------------------------------------------------------------------------- !

SUBROUTINE f_umPrint_opt_UsrPrefix(line, src, level, pe, model, UsrPrefix,     &
                                   HangIndent, stdErrorToo, setsum)

USE, INTRINSIC :: ISO_C_BINDING, ONLY:                                         &
  C_CHAR

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) ::                                                &
  line,                                                                        &
  src

INTEGER, OPTIONAL, INTENT(IN) ::                                               &
  level,                                                                       &
  pe

CHARACTER(LEN=*), OPTIONAL, INTENT(IN) ::                                      &
  model,                                                                       &
  UsrPrefix

CHARACTER(KIND=C_CHAR,LEN=1), INTENT(IN), OPTIONAL ::                          &
  HangIndent(*)

LOGICAL, INTENT(IN) ::                                                         &
  stdErrorToo

INTEGER(KIND=C_INT), INTENT(IN) ::                                             &
  setsum

IF (IAND(setsum,HangIndent_set)==0) THEN

  CALL umPrint(line=line,                                                      &
               src=src,                                                        &
               level=level,                                                    &
               pe=pe,                                                          &
               model=model,                                                    &
               UsrPrefix=UsrPrefix,                                            &
               stdErrorToo=stdErrorToo)

ELSE

  CALL umPrint(line=line,                                                      &
               src=src,                                                        &
               level=level,                                                    &
               pe=pe,                                                          &
               model=model,                                                    &
               UsrPrefix=UsrPrefix,                                            &
               HangIndent=TRIM(f_shum_c2f_string(HangIndent)),                 &
               stdErrorToo=stdErrorToo)

END IF

END SUBROUTINE f_umPrint_opt_UsrPrefix

! ---------------------------------------------------------------------------- !

END MODULE fort2c_umprint_interfaces

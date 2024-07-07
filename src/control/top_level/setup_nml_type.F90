! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
! Description:
!   Module containing routine to setup namelist types
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 95

MODULE setup_namelist

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SETUP_NAMELIST'

CONTAINS

SUBROUTINE setup_nml_type(no_of_types, mpl_nml_type,                   &
                          n_int_in, n_real_in, n_log_in, n_chars_in)

USE mpl, ONLY: mpl_integer, mpl_real, mpl_address_kind, mpl_logical, &
                mpl_character

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN)           :: no_of_types
INTEGER, OPTIONAL, INTENT(IN) :: n_int_in
INTEGER, OPTIONAL, INTENT(IN) :: n_real_in
INTEGER, OPTIONAL, INTENT(IN) :: n_log_in
INTEGER, OPTIONAL, INTENT(IN) :: n_chars_in
INTEGER, INTENT(OUT) :: mpl_nml_type

INTEGER :: oldtypes(0:no_of_types-1)
INTEGER :: blockcounts(0:no_of_types-1)
INTEGER (KIND=mpl_address_kind) :: offsets(0:no_of_types-1)
INTEGER (KIND=mpl_address_kind) :: extent,extent2, extent3
INTEGER (KIND=mpl_address_kind) :: lb
INTEGER :: icode
INTEGER :: counter
INTEGER :: n_int
INTEGER :: n_real
INTEGER :: n_log
INTEGER :: n_chars

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETUP_NML_TYPE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

counter=0

IF (PRESENT(n_int_in)) THEN
  n_int=n_int_in
  offsets    (counter) = 0
  oldtypes   (counter) = mpl_integer
  blockcounts(counter) = n_int
  CALL mpl_type_get_extent(mpl_integer, lb, extent, icode)
  counter=counter+1
ELSE
  extent=0
  n_int=0
END IF

IF (PRESENT(n_real_in)) THEN
  n_real=n_real_in
  offsets    (counter) = n_int * extent
  oldtypes   (counter) = mpl_real
  blockcounts(counter) = n_real
  CALL mpl_type_get_extent(mpl_real, lb, extent2, icode)
  counter=counter+1
ELSE
  extent2=0
  n_real=0
END IF

IF (PRESENT(n_log_in)) THEN
  n_log=n_log_in
  offsets    (counter) = n_int * extent + n_real * extent2
  oldtypes   (counter) = mpl_logical
  blockcounts(counter) = n_log
  CALL mpl_type_get_extent(mpl_logical, lb, extent3, icode)
  counter=counter+1
ELSE
  extent3=0
  n_log=0
END IF

IF (PRESENT(n_chars_in)) THEN
  n_chars=n_chars_in
  offsets(counter) = n_int * extent + n_real * extent2           &
                     + n_log * extent3
  oldtypes(counter) = mpl_character
  blockcounts(counter) = n_chars
ELSE
  n_chars=0
END IF

CALL mpl_type_create_struct(no_of_types, blockcounts, offsets,   &
                            oldtypes, mpl_nml_type, icode)
CALL mpl_type_commit(mpl_nml_type, icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE setup_nml_type

END MODULE setup_namelist

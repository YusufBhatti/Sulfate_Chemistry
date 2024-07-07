! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:
!
! Subroutine interface:
SUBROUTINE read_pp_header (ftin1,pp_int,pp_real,ibm_to_cray,      &
                           l_bit_32)

USE um_types
USE fort2c_data_conv_interfaces, ONLY: ibm2ieee, real_type, integer_type

IMPLICIT NONE
!
! Description:
! This reads the pp headers from a Fortan unformatted file.
!
! Method:
! Data is read either into a 32 bit word, or a 64 bit depending on
! the l_bit_32 or ibm_to_cray logicals. If data is in IBM format it
! is converted to IEEE before passing back. If data is IEEE 32 bit,
! it is converted to IEEE 64 bit.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
!
! Declarations:
!
! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER (KIND=integer64) :: ftin1

#if defined (NEC_FORTRAN_SX)
! The NEC cross-compiler requires "kind 8" logicals due to the way it
! enforcing type promotion.
LOGICAL (KIND=8)  :: ibm_to_cray
LOGICAL (KIND=8)  :: l_bit_32
#else
LOGICAL :: ibm_to_cray
LOGICAL :: l_bit_32
#endif

!   Array  arguments with intent(in):
INTEGER (KIND=integer64) :: pp_int(45)
REAL    (KIND=real64)    :: pp_real(19)

! local scalars
INTEGER (KIND=integer64) :: i
INTEGER (KIND=integer64) ::ier

! local arrays
INTEGER (KIND=integer64) :: pp_buffer(32)
INTEGER (KIND=integer32) :: pp_int_32(45)
REAL    (KIND=real32)    :: pp_real_32(19)

! End of header

IF (ibm_to_cray) THEN

  ! Read in the PP header
  READ(ftin1) pp_buffer

  ! Convert Integer part of header (Words 1-45)
  ier = ibm2ieee(integer_type,45,pp_buffer(1:23),0,pp_int(1:45),1,64,32)

  ! Convert Real part of header (Words 46-64)
  ier = ibm2ieee(real_type,19,pp_buffer(23:32),32,pp_real(1:19),1,64,32)

ELSE

  ! Read in the PP header
  IF (L_bit_32) THEN
    READ(ftin1) pp_int_32,pp_real_32
    DO i = 1, 45
      pp_int(i) = pp_int_32(i)
    END DO

    DO i = 1, 19
      pp_real(i) = pp_real_32(i)
    END DO
  ELSE
     ! read in straight 64 bit data
    READ(ftin1) pp_int, pp_real
  END IF  ! 32 bit

END IF


RETURN
END SUBROUTINE read_pp_header

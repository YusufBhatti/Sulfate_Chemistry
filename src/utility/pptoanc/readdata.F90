! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:
SUBROUTINE readdata(rows,columns,ftin1,ibm_to_cray,               &
                    len_extra,l_bit_32,l_skip,                    &
                    DATA, extra_data)

USE um_types
USE fort2c_data_conv_interfaces, ONLY: ibm2ieee, real_type

IMPLICIT NONE
!
! Description:
! Reads PP data from a Fortran unformatted file.
!
! Method:
! Data is read either into a 32 bit word, or a 64 bit depending on
! the l_bit_32 or ibm_to_cray logicals. If data is in IBM format it
! is converted to IEEE before passing back. If data is IEEE 32 bit,
! it is converted to IEEE 64 bit.
!
! If l_skip is true, the data is read, but nothing is passed back.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
!
! Declarations:
!
! Subroutine arguments
!   Scalar arguments with intent(in):


INTEGER (KIND=integer64):: rows
INTEGER (KIND=integer64):: columns
INTEGER (KIND=integer64):: ftin1
INTEGER (KIND=integer64):: len_extra

#if defined (NEC_FORTRAN_SX)
! The NEC cross-compiler requires kind=8 logicals due to the way it
! enforces type promotion.
LOGICAL (KIND=8) :: ibm_to_cray
LOGICAL (KIND=8) :: l_bit_32
LOGICAL (KIND=8) :: l_skip
#else
LOGICAL :: ibm_to_cray
LOGICAL :: l_bit_32
LOGICAL :: l_skip
#endif

REAL (KIND=real64) :: DATA(rows*columns)
REAL (KIND=real64) :: extra_data(len_extra+1)

! Local scalars:
INTEGER (KIND=integer64) ::  i
INTEGER (KIND=integer64) :: ierr

! local arrays:
REAL (KIND=real64) :: field2(rows*columns)
REAL (KIND=real64) :: extra_data2(len_extra+1)

REAL (KIND=real32) :: extra_data1(len_extra)
REAL (KIND=real32) :: field1(rows*columns)

! local parameters
INTEGER(KIND=integer64), PARAMETER :: i64b32  = 32
INTEGER(KIND=integer64), PARAMETER :: i64b64  = 64
INTEGER(KIND=integer64), PARAMETER :: i64unit =  1
INTEGER(KIND=integer64), PARAMETER :: i64zero =  0

! End of header

IF (ibm_to_cray .OR. l_bit_32) THEN     ! 32 bit read
  IF (len_extra >  0) THEN
    READ (ftin1) field1,(extra_data1(i),i=1,len_extra)
  ELSE
    READ (ftin1) field1
  END IF

  ! If not skipping, copy data to output in relevant
  ! format
  IF (.NOT. l_skip) THEN
    IF (ibm_to_cray) THEN

      ierr = ibm2ieee(real_type, rows*columns, field1(1:rows*columns),         &
                      i64zero, DATA(1:rows*columns), i64unit, i64b64, i64b32)

      ierr = ibm2ieee(real_type, rows*columns, extra_data1(1:rows*columns),    &
                      i64zero, extra_data(1:rows*columns), i64unit, i64b64,    &
                      i64b32)

    ELSE      ! straight 32 bit data

      DO i = 1, rows*columns
        DATA(i) = field1(i)
      END DO

      DO i = 1, len_extra
        extra_data(i) = extra_data1(i)
      END DO

    END IF ! ibm_to_cray
  END IF   ! l_skip


ELSE  ! 64 bit read
  IF (l_skip) THEN   ! read into junk array
    IF (len_extra >  0) THEN
      READ (ftin1) field2,(extra_data2(i),i=1,len_extra)
    ELSE
      READ (ftin1) field2
    END IF
  ELSE              ! read into output array
    IF (len_extra >  0) THEN
      READ (ftin1) DATA,(extra_data(i),i=1,len_extra)
    ELSE
      READ (ftin1) DATA
    END IF
  END IF
END IF

RETURN
END SUBROUTINE readdata
!
! Subroutine interface:

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:
SUBROUTINE calc_len_cfi(ftin2,cols_nowrap,len1_rowdepc,           &
     nlevels,len_cfi,fldsizelev,ibm_to_cray,                      &
     add_wrap_pts,l_bit_32,icode)

USE fort2c_data_conv_interfaces, ONLY: ibm2ieee, real_type

USE um_types, ONLY: real32, real64
USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage
IMPLICIT NONE
!
! Description:
!           this subroutine calculates the dimensions of the
!           compression arrays. It is a subset of the subroutine
!           calc_cfi_and_fld
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Subroutine arguments
!   Scalar arguments with intent(in):

INTEGER :: ftin2           ! (in) unit number for levels dataset
INTEGER :: cols_nowrap     ! (in) number of points east-west
                        !      (without wrap points)
INTEGER :: len1_rowdepc    ! (in) number of points north-south
INTEGER :: nlevels         ! (in) number of points in vertical


LOGICAL :: ibm_to_cray   ! T => input pp data is in IBM number
                      !      format and needs to be converted to
                      !      run on the Cray.
LOGICAL :: add_wrap_pts  ! T => add wrap points to the output file
LOGICAL :: l_bit_32

CHARACTER(LEN=80) :: levels

INTEGER :: icode      ! error code

! Array arguments with intent(in):

INTEGER :: len_cfi(3) ! (out) total number of sea segments
INTEGER :: fldsizelev(nlevels) !(out) no. of points on each level
                            ! of compressed field


! Local Scalars

INTEGER :: columns    ! no. of columns in levels dataset
INTEGER :: rows       ! no. of rows in levels dataset
INTEGER :: i,j,k      ! local loop indices
INTEGER :: ierr       ! return code from ibm2ieee
INTEGER :: COUNT      ! local counter for points in a sea segment
INTEGER :: seg_count  ! local counter for number of sea segments
INTEGER :: last_count ! local counter for calcultaing points in a
                   ! sea segment
! Local dynamic arrays:

INTEGER :: pp_int(45)     ! integer part of levels lookup header
REAL :: pp_real(19)       ! real part of levels lookup header

REAL(KIND=real32) :: levels_in(cols_nowrap*len1_rowdepc)
                        ! temp array for levels dataset to be
                        ! converted to cray number format
REAL(KIND=real64) :: levels_array(cols_nowrap,len1_rowdepc)
                        ! array of ocean levels

!- End of header

!  1. Take required dimensions from levels dataset

!  1.2 Obtain columns and rows by reading header

! DEPENDS ON: read_pp_header
CALL read_pp_header(ftin2,pp_int,pp_real,ibm_to_cray,l_bit_32)

rows = pp_int(18)
columns = pp_int(19)

WRITE(umMessage,*)'rows = ',rows
CALL umPrint(umMessage,src='calc_len_cfi')
WRITE(umMessage,*)'columns = ',columns
CALL umPrint(umMessage,src='calc_len_cfi')

!  1.3 Check the dimensions and read the levels_array.

IF (len1_rowdepc  /=  rows) THEN
  WRITE(umMessage,*)'wrong number of rows in SIZES namelist'
  CALL umPrint(umMessage,src='calc_len_cfi')
  WRITE(umMessage,*)'len1_rowdepc should equal rows in levels dataset'
  CALL umPrint(umMessage,src='calc_len_cfi')
  WRITE(umMessage,*)'resubmit'
  CALL umPrint(umMessage,src='calc_len_cfi')
  icode = 222
  GO TO 9999      ! Jump out
END IF

IF (cols_nowrap  /=  columns) THEN
  WRITE(umMessage,*)'wrong number of columns in SIZES namelist'
  CALL umPrint(umMessage,src='calc_len_cfi')
  WRITE(umMessage,*)'len1_coldepc should equal columns in levels dataset'
  CALL umPrint(umMessage,src='calc_len_cfi')
  WRITE(umMessage,*)'resubmit'
  CALL umPrint(umMessage,src='calc_len_cfi')
  icode = 223
  GO TO 9999      ! Jump out
END IF

!  Do number conversion if required

IF (ibm_to_cray) THEN
  READ(ftin2) levels_in

  ierr = ibm2ieee(real_type, rows*columns, levels_in(1:rows*columns), 0,       &
                  levels_array(1:columns,1:rows), 1, 64, 32)

  ierr = 0
ELSE
  READ(ftin2)levels_array
END IF

CLOSE(ftin2)

!  2. Calculate len_cfi and fldsizelev

!  2.1 Loop over the points in the field to calculate the number of
!  segments

COUNT=0
seg_count=0
last_count=0

DO k=1,nlevels
  DO j=1,rows

    IF (k  <=  levels_array(1,j)) THEN
      COUNT = COUNT + 1
      seg_count = seg_count + 1
    END IF

    DO i=2,cols_nowrap

      IF (k  <=  levels_array(i,j)) THEN
        COUNT = COUNT + 1
      END IF

      IF ((k  >   levels_array(i-1,j)) .AND.                          &
                       (k <= levels_array(i,j))) THEN
        seg_count=seg_count+1
      END IF

    END DO
  END DO

  fldsizelev(k) = COUNT - last_count
  WRITE(umMessage,*)'k = ',k
  CALL umPrint(umMessage,src='calc_len_cfi')
  WRITE(umMessage,*)'fldsizelev(k) = ',fldsizelev(k)
  CALL umPrint(umMessage,src='calc_len_cfi')

  last_count = COUNT

END DO

len_cfi(1) = seg_count
len_cfi(2) = seg_count
len_cfi(3) = rows * nlevels

WRITE(umMessage,*)'len_cfi(1) = ',len_cfi(1)
CALL umPrint(umMessage,src='calc_len_cfi')
WRITE(umMessage,*)'len_cfi(2) = ',len_cfi(2)
CALL umPrint(umMessage,src='calc_len_cfi')
WRITE(umMessage,*)'len_cfi(3) = ',len_cfi(3)
CALL umPrint(umMessage,src='calc_len_cfi')


9999 CONTINUE
RETURN
END SUBROUTINE calc_len_cfi

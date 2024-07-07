! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
! Subroutine interface:
SUBROUTINE calc_cfi_and_fld(ftin2,nlevels,len1_coldepc,           &
      cols_nowrap,len1_rowdepc,len1_flddepc,len2_flddepc,         &
      fields_const,fldsizelev,len_cfi,cfi1,cfi2,cfi3,compress,    &
      flddepc,ibm_to_cray,add_wrap_pts,imdi,l_bit_32,icode)

USE filenamelength_mod, ONLY:                                          &
    filenamelength

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage

USE fort2c_interfaces, ONLY: get_file

IMPLICIT NONE
!
! Description:
!     this subroutine calculates the compression arrays:
!     cfi1(len_cfi(1)), cfi2(len_cfi(2)) and cfi3(len_cfi(3))
!     using an array of numbers of ocean levels at each point:
!     levels_array(len1_coldepc,len1_rowdepc)
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
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER :: ftin2      ! (in) unit numbr for levels dataset
INTEGER :: nlevels    ! (in) number of points in vertical

INTEGER :: len1_coldepc ! (in) 1st dimension of col_dep_consts
INTEGER :: cols_nowrap  ! (in) no. of points east-west
INTEGER :: len1_rowdepc ! (in) 1st dimension of row_dep_consts
INTEGER :: len1_flddepc ! (in) 1st dimension of fields_const
INTEGER :: len2_flddepc ! (in) 2nd dimension of fields_const
INTEGER :: imdi         ! (in) integer missing data indicator
INTEGER :: icode        ! error code

LOGICAL :: compress     ! T => the dump is to be compressed
LOGICAL :: flddepc      ! T => fields_const are wanted in the dump
LOGICAL :: ibm_to_cray  ! T => input pp data is in IBM number
                     !      format and needs to be converted to
                     !      run on the Cray.

LOGICAL :: add_wrap_pts ! T => add wrap points to the output field
LOGICAL :: l_bit_32


!   Array arguments with intent(in):


INTEGER :: fldsizelev(nlevels) ! number of points on each compressed
                            ! level

INTEGER :: len_cfi(3) ! (in) total number of sea segments

INTEGER :: cfi1(len_cfi(1))  ! (out) index array for compressed array
INTEGER :: cfi2(len_cfi(2))  ! (out) index array for expanded array
INTEGER :: cfi3(len1_rowdepc,nlevels)  ! (out)
               ! contains number of first sea
               ! segment in each row at each levelc
               ! if there is a sea segment in the row
               ! contains number of next sea segment
               ! otherwise

REAL :: fields_const(len1_flddepc,len2_flddepc) ! (out) array for
                     ! fields of constants

! Local scalars :

INTEGER :: columns  ! no. of points east-west
INTEGER :: rows     ! no. of points north-south
INTEGER :: i,j,k    ! local loop indices
INTEGER :: COUNT    ! local counter for points in a sea segment
INTEGER :: seg_count! local counter for number of sea segments

LOGICAL :: l_skip ! If T, the data is read, but nothing is passed back.

CHARACTER(LEN=filenamelength) :: levels

! Local dynamic arrays :

INTEGER :: pp_int(45)

REAL :: pp_real(19)
REAL :: temp_levels_array(len1_coldepc,len1_rowdepc)
                            ! local array of ocean levels
REAL :: levels_array(cols_nowrap,len1_rowdepc)
                            ! local array of ocean levels
REAL :: dummy(1)


!- End of header

!  1. Read the fields_const from levels dataset

!  1.1 Read the data from levels dataset
CALL get_file(ftin2,levels,filenamelength,icode)
levels=TRIM(levels)
OPEN( UNIT=ftin2, FILE=levels, FORM="unformatted" )

! DEPENDS ON: read_pp_header
CALL read_pp_header(ftin2,pp_int,pp_real,ibm_to_cray,l_bit_32)

rows = pp_int(18)
columns = pp_int(19)

WRITE(umMessage,*)'rows = ',rows
CALL umPrint(umMessage,src='calc_cfi_and_fld')
WRITE(umMessage,*)'columns = ',columns
CALL umPrint(umMessage,src='calc_cfi_and_fld')

!  1.4 Read in levels_array and check the dataset is on the
!  same grid as the input pp fields.  If add_wrap_pts and flddepc
!  then add wrap points to the levels dataset.  The compression indices
!  are the same for an output dump with or without wrap points.

l_skip = .FALSE.
! DEPENDS ON: readdata
CALL readdata( rows, columns, ftin2, ibm_to_cray, 0,              &
               l_bit_32, l_skip, levels_array, dummy)
CLOSE(ftin2)

IF (add_wrap_pts) THEN

  IF (len1_rowdepc  /=  rows) THEN
    WRITE(umMessage,*)'wrong number of rows in SIZES namelist'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    WRITE(umMessage,*)'len1_rowdepc should equal rows in levels dataset'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    WRITE(umMessage,*)'resubmit'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    icode = 222
    GO TO 9999      ! Jump out
  END IF

  IF (len1_coldepc  /=  columns+2) THEN
    WRITE(umMessage,*)'wrong number of columns in SIZES namelist'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    WRITE(umMessage,*)'len1_coldepc should equal columns+2 in '//&
        'levels dataset'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    WRITE(umMessage,*)'resubmit'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    icode = 223
    GO TO 9999      ! Jump out
  END IF

  IF (len1_coldepc*len1_rowdepc  /=  len1_flddepc) THEN
    WRITE(umMessage,*)'len1_flddepc should equal '//&
        'len1_coldepc*len1_rowdepc'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    WRITE(umMessage,*)'resubmit'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    icode = 224
    GO TO 9999      ! Jump out
  END IF

  DO j = 1,rows
    DO i = 1,columns
      temp_levels_array(i,j) = levels_array(i,j)
    END DO
  END DO

  DO j = 1,rows
    temp_levels_array(columns+1,j)=temp_levels_array(1,j)
    temp_levels_array(columns+2,j)=temp_levels_array(2,j)
  END DO

  DO j = 1,rows
    DO i = 1,len1_coldepc
      fields_const(i+(j-1)*len1_coldepc,1) = temp_levels_array(i,j)
    END DO
  END DO

ELSE

  IF (len1_rowdepc  /=  rows) THEN
    WRITE(umMessage,*)'wrong number of rows in SIZES namelist'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    WRITE(umMessage,*)'len1_rowdepc should equal rows in levels dataset'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    WRITE(umMessage,*)'resubmit'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    icode = 222
    GO TO 9999      ! Jump out
  END IF

  IF (len1_coldepc  /=  columns) THEN
    WRITE(umMessage,*)'wrong number of columns in SIZES namelist'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    WRITE(umMessage,*)'len1_coldepc should equal columns in '//&
        'levels dataset'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    WRITE(umMessage,*)'resubmit'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    icode = 223
    GO TO 9999      ! Jump out
  END IF

  IF (len1_coldepc*len1_rowdepc  /=  len1_flddepc) THEN
    WRITE(umMessage,*)'len1_flddepc should equal '//&
        'len1_coldepc*len1_rowdepc'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    WRITE(umMessage,*)'resubmit'
    CALL umPrint(umMessage,src='calc_cfi_and_fld')
    icode = 224
    GO TO 9999      ! Jump out
  END IF

  DO j = 1,rows
    DO i = 1,len1_coldepc
      fields_const(i+(j-1)*len1_coldepc,1) = levels_array(i,j)
    END DO
  END DO

END IF

!  2.1 Initialise cfi3 array and create the compression indices.

IF (compress) THEN

  DO k=1,nlevels
    DO j=1,rows
      cfi3(j,k)=imdi
    END DO
  END DO

  COUNT=0
  seg_count=0

  DO k=1,nlevels
    DO j=1,rows
      !
      !     if the first element in a row is sea, a new segment is starting,
      !     so count and seg_count are both incremented, and cfi1 and
      !     cfi2 have new entries.  Columns is used here instead of
      !     len1_coldepc as the index to this array expects that in oa_pack.
      !
      IF (k <= levels_array(1,j)) THEN
        COUNT=count+1
        seg_count=seg_count+1
        cfi1(seg_count)=COUNT
        cfi2(seg_count)=1+(j-1)*columns+(k-1)*columns*rows
        cfi3(j,k)=seg_count
      END IF

      DO i=2,columns
        !
        !     if present point is sea, add one to count
        !
        IF (k <= levels_array(i,j)) THEN
          COUNT=count+1
        END IF
        !
        !     if present point is sea and previous point is land,
        !     a new segment is starting, so seg_count is incremented
        !     and cfi1 and cfi2 have new entries. Columns is used here instead
        !     of len1_coldepc as the index to this array expects that
        !     in oa_pack.
        !
        IF ((k >  levels_array(i-1,j)) .AND.                      &
                         (k <= levels_array(i,j))) THEN
          seg_count=seg_count+1
          cfi1(seg_count)=COUNT
          cfi2(seg_count)=i+(j-1)*columns+(k-1)*columns*rows
          !
          !     if cfi3(j,k) has not been reset,
          !     then the present segment must be the first in the row
          !
          IF (cfi3(j,k) == imdi) THEN
            cfi3(j,k)=seg_count
          END IF
        END IF
      END DO

      !
      !     if there is no sea segment in the row,
      !     then set cfi3 to seg_count+1
      !
      IF (cfi3(j,k) == imdi) THEN
        cfi3(j,k)=seg_count+1
      END IF

    END DO
  END DO

END IF
9999  CONTINUE
RETURN
END SUBROUTINE calc_cfi_and_fld

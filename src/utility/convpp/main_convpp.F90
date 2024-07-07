! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    PROGRAM MAIN_CONVPP --------------------------------------------
!
!    Purpose: Converts a UM file into PP format.
!
!    Programming standards:
!
!    Logical components covered:
!
!    System Tasks: F3,F4,F6
!
!    Documentation: UM Doc Paper F5
!
!    -----------------------------------------------------------------
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Small execs
PROGRAM main_convpp

USE filenamelength_mod, ONLY:                                    &
    filenamelength
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE io, ONLY: setpos, file_close, file_open, buffin
USE io_constants, ONLY: ioNoDelete, ioOpenReadOnly
USE ereport_mod, ONLY: ereport
USE UM_Config, ONLY: &
    appInit,          &
    appTerminate,     &
    exe_convpp
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE get_env_var_mod, ONLY: get_env_var
USE errormessagelength_mod, ONLY: errormessagelength
USE hostname_mod,           ONLY: get_hostname
USE ini_ppheader_mod, ONLY: ini_ppheader, ibm
IMPLICIT NONE


CHARACTER(LEN=filenamelength) :: file_in       ! Input file 
CHARACTER(LEN=filenamelength) :: file_out      ! Output file
CHARACTER(LEN=filenamelength) :: namelist_file ! Namelist

INTEGER :: namelist_unit

INTEGER ::                                                        &
 fixhd(256)                                                       &
                   !Space for fixed length header
,inthd(100)        !Space for integer header

INTEGER ::                                                        &
 len_fixhd                                                        &
                !Length of fixed length header on input file
,len_inthd                                                        &
                !Length of integer header on input file
,joc_no_seapts                                                    &
                !Number of points in compressed array
,len_ocfld                                                        &
                !Length of uncompressed ocean field
,len_realhd                                                       &
                !Length of real header on input file
,len1_levdepc                                                     &
                !1st dim of lev dependent consts on input file
,len2_levdepc                                                     &
                !2nd dim of lev dependent consts on input file
,len1_rowdepc                                                     &
                !1st dim of row dependent consts on input file
,len2_rowdepc                                                     &
                !2nd dim of row dependent consts on input file
,len1_coldepc                                                     &
                !1st dim of col dependent consts on input file
,len2_coldepc                                                     &
                !2nd dim of col dependent consts on input file
,len1_flddepc                                                     &
                !1st dim of field dependent consts on input file
,len2_flddepc                                                     &
                !2nd dim of field dependent consts on input file
,len_extcnst                                                      &
                !Length of extra consts on input file
,len_dumphist                                                     &
                !Length of history header on input file
,len_cfi1                                                         &
                !Length of index1 on input file
,len_cfi2                                                         &
                !Length of index2 on input file
,len_cfi3                                                         &
                !Length of index3 on input file
,len1_lookup                                                      &
                !1st dim of LOOKUP on input file
,len2_lookup                                                      &
                !2nd dim of LOOKUP on input file
,len_data                                                         &
                !Length of data on input file
,row_length                                                       &
                !No of points E-W on input file
,p_rows                                                           &
                !No of p-rows on input file
,p_field                                                          &
                !No of p-points per level on input file
,max_field_size !Maximum field size on file

INTEGER ::                                                        &
 len_io                                                           &
          !Length of I/O returned by BUFFER IN
,i                                                                &
          !Loop index
,nftin                                                            &
          !Unit number of input UM dump 1
,nftout                                                           &
          !Unit number of output file
,err                                                              &
          !Return code from OPEN

,icode    !Return code from setpos
REAL :: a    !BUFFER IN UNIT function
INTEGER :: iocode ! IO status 
CHARACTER(LEN=errormessagelength) :: iomessage ! IO message
CHARACTER(LEN=errormessagelength) :: cmessage ! Error message
CHARACTER(LEN=11), PARAMETER :: routinename='main_convpp'

CALL appInit(exe_convpp)

! Write out hostname
WRITE(umMessage,'(A,A)') 'Host is ',TRIM(get_hostname())
CALL umPrint(umMessage, src='main_convpp')

CALL get_env_var("NAMELIST", namelist_file)

CALL assign_file_unit(namelist_file, namelist_unit, handler="fortran")
OPEN(namelist_unit, FILE=namelist_file, FORM='FORMATTED', ACTION='READ')
CALL ini_ppheader(namelist_unit)
CALL release_file_unit(namelist_unit, handler="fortran")

!  1. Assign unit numbers

WRITE(umMessage,'(20x,''FILE STATUS'')')
CALL umPrint(umMessage,src='main_convpp')
WRITE(umMessage,'(20x,''==========='')')
CALL umPrint(umMessage,src='main_convpp')

CALL get_env_var("FILE1", file_in)

CALL get_env_var("FILE2", file_out)

CALL assign_file_unit(file_in,  nftin, handler="portio")
CALL assign_file_unit(file_out, nftout, handler="fortran")

CALL file_open(nftin, file_in, filenamelength, &
    read_write=ioOpenReadOnly,error=err)

! PP files are defined as 32-bit big-endian files. For fieldsfiles this would
! be handled by the C I/O routines, but as we're writing pp data here
! we manually set this here to ensure consistency with the file format standard. 
OPEN(nftout,FILE=file_out,FORM='UNFORMATTED',CONVERT='big_endian',           &
     IOSTAT=iocode, IOMSG=iomessage)
IF (iocode /= 0) THEN
  cmessage = 'Error opening file ' // TRIM(file_out) //': ' // TRIM(iomessage)
  CALL ereport(routinename, iocode, cmessage)
END IF         


!  2. Buffer in fixed length header record

CALL buffin(nftin,fixhd,256,len_io,a)

! Check for I/O errors
IF (a /= -1.0 .OR. len_io /= 256) THEN
  ! DEPENDS ON: ioerror
  CALL ioerror('buffer in of fixed length header of input dump',  &
  a,len_io,256)
  iocode = 1000
  CALL ereport('MAIN_CONVPP', iocode,                             &
   'Buffer in of fixed length header of input dump wrong size')

END IF

! Set missing data indicator to zero
DO  i=1,256
  IF (fixhd(i) <  0)fixhd(i)=0
END DO

! Input file dimensions (ensure sizes are >= 0 for allocating)
len_fixhd    = 256
len_inthd    = MAX(fixhd(101),0)
len_realhd   = MAX(fixhd(106),0)
len1_levdepc = MAX(fixhd(111),0)
len2_levdepc = MAX(fixhd(112),0)
len1_rowdepc = MAX(fixhd(116),0)
len2_rowdepc = MAX(fixhd(117),0)
len1_coldepc = MAX(fixhd(121),0)
len2_coldepc = MAX(fixhd(122),0)
len1_flddepc = MAX(fixhd(126),0)
len2_flddepc = MAX(fixhd(127),0)
len_extcnst  = MAX(fixhd(131),0)
len_dumphist = MAX(fixhd(136),0)
len_cfi1     = MAX(fixhd(141),0)
len_cfi2     = MAX(fixhd(143),0)
len_cfi3     = MAX(fixhd(145),0)
len1_lookup  = MAX(fixhd(151),0)
len2_lookup  = MAX(fixhd(152),0)
len_data     = MAX(fixhd(161),0)


!  3. Buffer in integer constants from dump

CALL buffin(nftin,inthd,len_inthd,len_io,a)

! Check for I/O errors
IF (a /= -1.0 .OR. len_io /= len_inthd) THEN
  ! DEPENDS ON: ioerror
  CALL ioerror('buffer in of integer constants in input dump',    &
  a,len_io,len_inthd)
  iocode = 1001
  CALL ereport('MAIN_CONVPP', iocode,                             &
   'Buffer in of integer constants in input dump wrong size')

END IF

! Set missing data indicator to zero
DO  i=1,len_inthd
  IF (inthd(i) <  0)inthd(i)=0
END DO

row_length=inthd(6)
p_rows=inthd(7)
p_field=row_length*p_rows

!  Extract maximum field size from LOOKUP header
! DEPENDS ON: find_max_field_size
CALL find_max_field_size(                                        &
     nftin,len1_lookup,len2_lookup,fixhd,max_field_size)

! Calculate sizes of compressed and uncompressed ocean fields
joc_no_seapts=inthd(11)
IF (fixhd(2) == 2) THEN
  len_ocfld=inthd(6)*inthd(7)*inthd(8)
ELSE
  len_ocfld=0
END IF
! Rewind file
CALL setpos(nftin,0,icode)

IF (fixhd(2) == 1) THEN
  ! Atmospheric dump

  IF (ibm) THEN
    ! DEPENDS ON: atmos_convpp64_ibm
    CALL atmos_convpp64_ibm &
    (len_fixhd,len_inthd,len_realhd,              &
    len1_levdepc,len2_levdepc,len1_rowdepc,       &
    len2_rowdepc,len1_coldepc,len2_coldepc,       &
    len1_flddepc,len2_flddepc,len_extcnst,        &
    len_dumphist,len_cfi1,len_cfi2,len_cfi3,      &
    len1_lookup,len2_lookup,len_data,p_field,     &
    nftin,nftout,max_field_size)


  ELSE
    ! DEPENDS ON: atmos_convpp64
    CALL atmos_convpp64 &
    (len_fixhd,len_inthd,len_realhd,              &
    len1_levdepc,len2_levdepc,len1_rowdepc,       &
    len2_rowdepc,len1_coldepc,len2_coldepc,       &
    len1_flddepc,len2_flddepc,len_extcnst,        &
    len_dumphist,len_cfi1,len_cfi2,len_cfi3,      &
    len1_lookup,len2_lookup,len_data,p_field,     &
    nftin,nftout,max_field_size)
  END IF
ELSE IF (fixhd(2) == 2) THEN
  iocode = 1000
  CALL ereport('MAIN_CONVPP', iocode, 'Ocean data no longer supported')
END IF

CALL file_close(nftin, file_in, filenamelength, delete=ioNoDelete)
CLOSE(nftout)
CALL release_file_unit(nftin, handler="portio")
CALL release_file_unit(nftout, handler="fortran")

CALL appTerminate()

END PROGRAM main_convpp

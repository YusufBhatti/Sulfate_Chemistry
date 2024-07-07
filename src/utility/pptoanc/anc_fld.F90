! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
! Subroutine interface:
SUBROUTINE anc_fld(ftin2,ftout,nolevsmax,number_of_codes,         &
  max_n_pp_files,len_cfi,fldsizelev,                              &
  fld_types,n_times,nlevels,n_pp_files,stash_code,field_code,     &
  nlevs_code,unit_no,n_freq_waves,n_dir_waves,len_intc,len_realc, &
  len_extra,len1_levdepc,len2_levdepc,len1_rowdepc,len2_rowdepc,  &
  len1_coldepc,len2_coldepc,len1_flddepc,len2_flddepc,            &
  len_extcnst,rmdi_input,                                         &
  add_wrap_pts,periodic,single_time,ibm_to_cray,compress,wave,    &
  levdepc,rowdepc,coldepc,flddepc,extcnst,pack32,pphead,          &
  grid_of_tracer,field_order,lwfio,L_bit_32,                      &
  len2_lookup_max,cols_nowrap,icode)

USE filenamelength_mod, ONLY:                                    &
      filenamelength
USE errormessagelength_mod, ONLY: errormessagelength
USE io
USE um_types
USE check_iostat_mod
USE writhead_mod
USE UM_ParParams
USE Atmos_Max_Sizes
USE lookup_addresses

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE missing_data_mod, ONLY: rmdi,imdi

USE fort2c_interfaces, ONLY: get_file

IMPLICIT NONE


!
! Description:
!          Main subroutine. Creates the ancillary/dump header,
!          lookup tables and writes the header using WRITHEAD.
!          Calls dataw which writes the data out
!
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
! Code Description:
!   Language: FORTRAN 90
!
! Subroutine arguments
!   Scalar arguments with intent(in):

INTEGER :: ftin2                  ! input unit for mask file used
                               ! for fields consts and
                               ! compression indices.
INTEGER :: ftout                  ! unit number for output ancillary
                               ! file

INTEGER :: nolevsmax              ! max number of levels; dimensions
                               !  fldsizelev array
INTEGER :: number_of_codes        ! max number of stash/field codes
INTEGER :: max_n_pp_files         ! max number of input pp files



INTEGER :: fld_types      ! number of field types in I/O files
INTEGER :: n_times          ! number of time periods in I/O files
INTEGER :: nlevels          ! number of levels (default = 1)
INTEGER :: n_pp_files       ! number of input pp files


INTEGER :: n_freq_waves  ! number of wave frequencies
INTEGER :: n_dir_waves   ! number of wave directions

INTEGER :: len_intc      !  Actual  length of integer constants array
INTEGER :: len_realc     !  Actual  length of real constants array

INTEGER :: len_extra     ! length of extra data (minimum value = 0)

INTEGER :: len1_levdepc   ! Actual 1st dimension of lev_dep_consts
INTEGER :: len2_levdepc   ! Actual 2nd dimension of lev_dep_consts

INTEGER :: len1_rowdepc   ! Actual 1st dimension of row_dep_consts
INTEGER :: len2_rowdepc   ! Actual 2nd dimension of row_dep_consts

INTEGER :: len1_coldepc   ! Actual 1st dimension of col_dep_consts
INTEGER :: len2_coldepc   ! Actual 2nd dimension of col_dep_consts

INTEGER :: len1_flddepc   ! Actual 1st dimension of fields_const
INTEGER :: len2_flddepc   ! Actual 2nd dimension of fields_const

INTEGER :: len_extcnst    ! Actual 1st dimension of fields_const

INTEGER :: len2_lookup_max ! maximum 2nd dimension of the lookup
                        ! table
INTEGER :: cols_nowrap     ! no. of columns in field without wrap
INTEGER :: icode           ! error code variable

REAL :: rmdi_input     ! real :: missing data indicator
                       ! in input pp field

LOGICAL :: add_wrap_pts   ! T => adds wrapping columns
                       !      e.g. for global grid
LOGICAL :: periodic       ! T => periodic in time
                       !      e.g. climate field
LOGICAL :: single_time    ! T => all fields input valid at one time
LOGICAL :: ibm_to_cray    ! T => input pp data is in IBM number
                       !      format and needs to be converted to
                       !      run on the Cray.
                       !      (Only use if running on Cray)
LOGICAL :: compress       ! T => fields are packed into ancillary
                       !      field compressed field indices are
                       !      calculated
LOGICAL :: wave           ! T => a wave dump is to be created
LOGICAL :: levdepc        ! T => if level dependent constants array
                       !      required
LOGICAL :: rowdepc        ! T => if row dependant constants are
                       !      required
LOGICAL :: coldepc        ! T => if column dependant constants are
                       !      required
LOGICAL :: flddepc        ! T => if fields of constants are
                       !      required
LOGICAL :: extcnst        ! T => if fields of constants are
                       !      required
LOGICAL :: pack32         ! T => use 32 bit Cray numbers
LOGICAL :: pphead         ! T => print out pp headers read in


LOGICAL :: field_order    ! T => input pp fields ordered by time.
                       !      i.e. different months in input
                       !         files, same fields in all files
                       ! F => inout pp fields ordered by fields.
                       !      i.e. different fields in input
                       !         files, all months in all files
LOGICAL :: lwfio          ! T => set the LBEGIN and LBNREC fields
                       !      in the LOOKUP Headers for VN 16
                       !      Type Dumpfiles.
                       ! F => Old dumpfiles

LOGICAL :: l_bit_32       ! T => pp file input is 32 bit
                       ! F => pp file input is 64 bit


!   Array  arguments with intent(in):

INTEGER :: len_cfi(3)          ! lengths of compressed field indices
INTEGER :: fldsizelev(nolevsmax)! size of packed field on each level
INTEGER :: stash_code(number_of_codes) ! array of stash codes
INTEGER :: field_code(number_of_codes) ! array of field codes
INTEGER :: nlevs_code(number_of_codes) ! array of levels depending
                                    ! on field code
INTEGER :: unit_no(number_of_codes)    ! array of unit numbers for
                                    ! input
LOGICAL :: grid_of_tracer(number_of_codes) ! T => fields are on a
                                        ! tracer grid


! Local parameters:

INTEGER :: len_look_user     ! No. of changes to the lookup table
                          ! made by the user
INTEGER :: len1_lookup       ! 1st dimension of the lookup
INTEGER :: len1_lookup_all   ! Dimension of the whole lookup array
INTEGER :: lfh               ! length of the fixed length header

INTEGER :: max_len_intc      ! Max dimension of integer constants
INTEGER :: max_len_realc     ! Max dimension of real constants
INTEGER :: max_len1_levdepc  ! Max 1st dimension of lev_dep_consts
INTEGER :: max_len2_levdepc  ! Max 2nd dimension of lev_dep_consts
INTEGER :: max_len1_rowdepc  ! Max 1st dimension of row_dep_consts
INTEGER :: max_len2_rowdepc  ! Max 2nd dimension of row_dep_consts
INTEGER :: max_len1_coldepc  ! Max 1st dimension of col_dep_consts
INTEGER :: max_len2_coldepc  ! Max 2nd dimension of col_dep_consts
INTEGER :: max_len_extcnst   ! Max dimension of extra_const

INTEGER :: len_dumphist      ! Actual dimension of dumphist

PARAMETER (len_look_user = 50)
PARAMETER (len1_lookup = 45)
PARAMETER (len1_lookup_all = 64)
PARAMETER (lfh=256)

PARAMETER (max_len_intc=40)
PARAMETER (max_len_realc=40)
PARAMETER (max_len1_levdepc=100)
PARAMETER (max_len2_levdepc=5)
PARAMETER (max_len1_rowdepc=rows_max)
PARAMETER (max_len2_rowdepc=5)
PARAMETER (max_len1_coldepc=row_length_max)
PARAMETER (max_len2_coldepc=5)
PARAMETER (max_len_extcnst=1000)


PARAMETER (len_dumphist=1)

! Local Scalars


INTEGER :: ftin1            ! unit number for input pp fields

INTEGER :: len_data         ! length of the data record
INTEGER :: start_block      ! position of start for WRITHEAD

INTEGER :: n_sea_points     ! number of sea points for wave dump

INTEGER :: rows             ! no. of rows in input pp field
INTEGER :: columns          ! no. of columns in input pp field

INTEGER :: i,j              ! loop counters
INTEGER :: levn             ! level number
INTEGER :: m,n              ! loop counters
INTEGER :: np               ! Number of points in pp field

INTEGER :: len2_lookup      ! 2nd dimension of lookup table
                         ! total number of fields in output file
INTEGER :: len2_step        ! calculation step for len2_lookup

INTEGER :: nlevs_this_code  ! # of levels for this field code and
                         ! limit for do -loop over levels(waves)

INTEGER :: fieldn           ! present field number
INTEGER :: fieldsize        ! size of field when it is stored
                         ! in output data set
INTEGER :: runtot           ! running total of start address
                         ! in data array of present field


INTEGER :: no_cmp           ! total no. of compressed points in
                         ! compressed array
INTEGER :: ipos             ! position counter

INTEGER ::                                                        &
 disk_address                                                     &
                                 ! Current rounded disk address
,number_of_data_words_on_disk                                     &
                                 ! Number of data words on disk
,number_of_data_words_in_memory  ! Number of Data Words in memory

LOGICAL :: tracer_grid   ! T => field is on a tracer grid
                      ! F =>field is on a velocity grid

LOGICAL :: t_compress    ! compress argument for DATA subroutine
                      ! used for wave dump LSmask set f whatever
                      ! compress is


CHARACTER(LEN=errormessagelength) :: cmessage    ! error message from WRITHEAD
CHARACTER(LEN=filenamelength) :: ancfile, levels

! Local arrays dimensioned by parameters:

      ! arrays to overwrite integer lookup tables

INTEGER :: ifld_int(len_look_user)  ! int lookup fields to change
INTEGER :: item_int(len_look_user)  ! item number to change
INTEGER :: ival_int(len_look_user)  ! integer value to use

! arrays to overwrite real lookup tables

INTEGER :: ifld_real(len_look_user)  ! lookup fields to change
INTEGER :: item_real(len_look_user)  ! item number to change
REAL :: rval_real(len_look_user)  ! real :: value to use

INTEGER :: fixhd(lfh)          ! fixed length header

INTEGER :: int_const(max_len_intc)  ! integer constants
REAL :: real_const(max_len_realc)   ! real constants

REAL :: lev_dep_consts(1+max_len1_levdepc*max_len2_levdepc)
REAL :: row_dep_consts(1+max_len1_rowdepc*max_len2_rowdepc)
REAL :: col_dep_consts(1+max_len1_coldepc*max_len2_coldepc)
REAL :: extra_const(max_len_extcnst)
REAL :: dumphist(len_dumphist)

! Local dynamic arrays:

INTEGER :: pp_int(45)
INTEGER :: lookup(45,len2_lookup_max)   ! Integer part of lookup
                                     ! table array

REAL :: pp_real(19)
REAL :: rlookup(46:64,len2_lookup_max)  ! Integer part of lookup
                                     ! table array

INTEGER :: lookup_all(len1_lookup_all,len2_lookup_max)
                              ! Whole lookup table array

INTEGER :: cfi1(len_cfi(1))   ! compressed field index
INTEGER :: cfi2(len_cfi(2))   ! arrays
INTEGER :: cfi3(len_cfi(3))

INTEGER :: n_pp_flds(max_n_pp_files)  ! Number of pp fields array

LOGICAL :: lsmask(len1_coldepc*len1_rowdepc) ! land sea mask
                                          ! for wave dump

LOGICAL :: l_skip ! If T, the data is read, but nothing is passed back.

REAL :: fields_const(len1_flddepc,len2_flddepc)
                            ! array for fields of constants
REAL :: dummy(1)
CHARACTER(LEN=errormessagelength) :: iomessage
INTEGER               :: ErrorStatus
                           ! Return code : 0 Normal Exit : >0 Error

! Function & Subroutine calls:
INTEGER :: find_namelist

!- End of header

NAMELIST /header_data/ fixhd,int_const,real_const,                &
                       lev_dep_consts,row_dep_consts,             &
                       col_dep_consts,extra_const,                &
                       ifld_int, item_int, ival_int,              &
                       ifld_real, item_real, rval_real

!  0. Preliminaries

!  0.1 Check sizes namelist dimensions do not exceed the maximum
!  dimensions and intialise arrays.


IF (len_intc >  max_len_intc) THEN
  WRITE(umMessage,*) ' len_intc in namelist is too big.'
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Max value allowed is ',max_len_intc
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Increase MAX_LEN_INTC in program'
  CALL umPrint(umMessage,src='anc_fld')
  icode = 1
  GO TO 9999   !  Return
END IF

IF (len_realc >  max_len_realc) THEN
  WRITE(umMessage,*) ' len_realc in namelist is too big.'
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Max value allowed is ',max_len_realc
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Increase MAX_LEN_REALC in program'
  CALL umPrint(umMessage,src='anc_fld')
  icode = 2
  GO TO 9999   !  Return
END IF

IF (len1_levdepc >  max_len1_levdepc) THEN
  WRITE(umMessage,*) ' len1_levdpec in namelist is too big.'
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Max value allowed is ',max_len1_levdepc
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Increase MAX_LEN1_LEVDEPC in program'
  CALL umPrint(umMessage,src='anc_fld')
  icode = 3
  GO TO 9999   !  Return
END IF

IF (len2_levdepc >  max_len2_levdepc) THEN
  WRITE(umMessage,*) ' len2_levdpec in namelist is too big.'
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Max value allowed is ',max_len2_levdepc
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Increase MAX_LEN2_LEVDEPC in program'
  CALL umPrint(umMessage,src='anc_fld')
  icode = 4
  GO TO 9999   !  Return
END IF

IF (len1_rowdepc >  max_len1_rowdepc) THEN
  WRITE(umMessage,*) ' len1_rowdpec in namelist is too big.'
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Max value allowed is ',max_len1_rowdepc
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Increase MAX_LEN1_ROWDEPC in program'
  CALL umPrint(umMessage,src='anc_fld')
  icode = 5
END IF

IF (len2_rowdepc >  max_len2_rowdepc) THEN
  WRITE(umMessage,*) ' len2_rowdpec in namelist is too big.'
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Max value allowed is ',max_len2_rowdepc
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Increase MAX_LEN2_ROWDEPC in program'
  CALL umPrint(umMessage,src='anc_fld')
  icode = 6
  GO TO 9999   !  Return
END IF

IF (len1_coldepc >  max_len1_coldepc) THEN
  WRITE(umMessage,*) ' len1_coldpec in namelist is too big.'
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Max value allowed is ',max_len1_coldepc
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Increase MAX_LEN1_COLDEPC in program'
  CALL umPrint(umMessage,src='anc_fld')
  icode = 7
  GO TO 9999   !  Return
END IF

IF (len2_coldepc >  max_len2_coldepc) THEN
  WRITE(umMessage,*) ' len2_coldepc in namelist is too big.'
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Max value allowed is ',max_len2_coldepc
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Increase MAX_LEN2_COLDEPC in program'
  CALL umPrint(umMessage,src='anc_fld')
  icode = 8
  GO TO 9999   !  Return
END IF

IF (len_extcnst >  max_len_extcnst) THEN
  WRITE(umMessage,*) ' len_extcnst in namelist is too big.'
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Max value allowed is ',max_len_extcnst
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) ' Increase MAX_LEN_EXTCNST in program'
  CALL umPrint(umMessage,src='anc_fld')
  icode = 11
  GO TO 9999   !  Return
END IF

! Initialise namelist arrays.
DO n=1,max_len2_levdepc
  DO m=1,max_len1_levdepc
    lev_dep_consts(m+(n-1)*max_len1_levdepc)=0.0
  END DO
END DO

DO n=1,max_len2_rowdepc
  DO m=1,max_len1_rowdepc
    row_dep_consts(m+(n-1)*max_len1_rowdepc)=0.0
  END DO
END DO

DO n=1,max_len2_coldepc
  DO m=1,max_len1_coldepc
    col_dep_consts(m+(n-1)*max_len1_coldepc)=0.0
  END DO
END DO

DO n=1,max_len_extcnst
  extra_const(n)=0.0
END DO

DO n=1,len_dumphist
  dumphist(n)=0.0
END DO

!
!  0.2 Read wave land/sea mask
!
IF (wave .AND. compress) THEN

  WRITE(umMessage,*) 'reading in landsea mask for waves from pp dataset'
  CALL umPrint(umMessage,src='anc_fld')

  CALL get_file(ftin2,levels,filenamelength,icode)
  levels=TRIM(levels)
  OPEN(UNIT=ftin2, FILE=levels, FORM="unformatted",               &
       CONVERT="big_endian", IOSTAT=ErrorStatus, IOMSG=iomessage)
  IF (ErrorStatus /= 0) THEN
    cmessage = 'Error opening file ' // TRIM(levels) //': ' //    &
               TRIM(iomessage)
    CALL ereport('ANC_FLD', ErrorStatus, cmessage)
  END IF          

  ! DEPENDS ON: read_pp_header
  CALL read_pp_header(ftin2,pp_int,pp_real,ibm_to_cray,           &
                      l_bit_32)

  l_skip = .FALSE.
  ! DEPENDS ON: readdata
  CALL readdata( len1_rowdepc,len1_coldepc,ftin2,ibm_to_cray,     &
                 0, l_bit_32, l_skip, lsmask, dummy)

  CLOSE(ftin2)

  ! reset so sea points are true
  !
  n_sea_points=0

  DO i=1,len1_coldepc*len1_rowdepc

    lsmask(i)= .NOT. lsmask(i)

    IF (lsmask(i)) THEN
      n_sea_points=n_sea_points+1
    END IF

  END DO

  WRITE(umMessage,*)'n_sea_points set to ', n_sea_points
  CALL umPrint(umMessage,src='anc_fld')
  fldsizelev(1)=n_sea_points

END IF     ! wave .and. compress

!  0.3 Calculate number of fields (len2_lookup), which depends on
!      nlevs_code.

len2_lookup = 0
len2_step = 0

DO i =1,fld_types
  len2_step = n_times * nlevs_code(i)
  len2_lookup = len2_lookup + len2_step
END DO

WRITE(umMessage,*)' '
CALL umPrint(umMessage,src='anc_fld')
WRITE(umMessage,*)'len2_lookup = ',len2_lookup
CALL umPrint(umMessage,src='anc_fld')

!  0.4 Initialise arrays

DO n=1,max_n_pp_files
  n_pp_flds(n)=0
END DO

DO n = 1, len_look_user
  ifld_int(n) = imdi
  item_int(n) = imdi
  ival_int(n) = imdi
  ifld_real(n) = imdi
  item_real(n) = imdi
  rval_real(n) = rmdi
END DO

DO i = 1,len_dumphist
  dumphist(i)= 0.0
END DO

!  1. Read headers of all input data

!  1.0 Open the UM/ancillary file

CALL get_file(ftout,ancfile,filenamelength,icode)
CALL file_open (ftout,ancfile,filenamelength,1,1,icode)
IF (icode >  0) THEN
  WRITE(umMessage,*) ' Problem with opening Ancillary File on Unit, '
  CALL umPrint(umMessage,src='anc_fld')
  icode = 2
  GO TO 9999  ! Return
END IF

! note these values are used in pp_table so need to be set here
! before the do-loops
!     * set default lev_dep_consts for WAVE frequency using the
!       factor 1.1
!     * as set in real_const(15) by namelist
!     * need to set frmin as lev_dep_consts(1) in namelist input
!     * pick up the CO value from namelist input HEADER value

IF (wave) THEN

  REWIND(5)
  ! DEPENDS ON: find_namelist
  i=find_namelist(5,"header_data")

  IF (i == 0) THEN
    READ (UNIT=5, NML=header_data, IOSTAT=ErrorStatus, IOMSG=iomessage)
    CALL check_iostat(errorstatus, "namelist HEADER_DATA", iomessage)
  ELSE
    WRITE(umMessage,*)'Cannot find namelist HEADER_DATA'
    CALL umPrint(umMessage,src='anc_fld')
  END IF

  WRITE(umMessage,*)'real_const(15) is',real_const(15)
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*)'lev_dep_consts(1,1)=',lev_dep_consts(1)
  CALL umPrint(umMessage,src='anc_fld')

  DO m=2,len1_levdepc
    lev_dep_consts(m) = real_const(15)*lev_dep_consts(m-1)
  END DO

END IF

!  1.1 Read through all data sets calculating the ancillary
!     file headers.  Loop over n_times then fld_types.

runtot=1  ! points to start point in data array for next field
fieldn=0  ! field number counter

DO n=1,n_times

  DO m=1,fld_types

    fieldn = fieldn + 1

    !  1.2 do steps which are independent of the level first
    !  Read the pp header for each field

    IF (n_pp_files == 1) THEN
      ftin1= unit_no(1)
    ELSE
      IF (field_order) THEN
        ftin1= unit_no(n)
      ELSE
        ftin1= unit_no(m)
      END IF
    END IF

    ! DEPENDS ON: read_pp_header
    CALL read_pp_header (ftin1,pp_int,pp_real,ibm_to_cray,l_bit_32)

    n_pp_flds(ftin1-29) = n_pp_flds(ftin1-29)+1
    WRITE(umMessage,*) 'Field No ',fieldn,' read in. From PP File ', &
        ftin1-29,' Field No ',n_pp_flds(ftin1-29)
    CALL umPrint(umMessage,src='anc_fld')


    IF (pphead) THEN
      WRITE(umMessage,*) 'pp_int for field ',fieldn
      CALL umPrint(umMessage,src='anc_fld')
      DO j=1,145
        WRITE(umMessage,*) pp_int(j)
        CALL umPrint(umMessage,src='anc_fld')
      END DO
      WRITE(umMessage,*) 'pp_real for field ',fieldn
      CALL umPrint(umMessage,src='anc_fld')
      DO j=1,19
        WRITE(umMessage,*) pp_real(j)
        CALL umPrint(umMessage,src='anc_fld')
      END DO
    END IF


    !  1.3 Extract the data dimensions and determine whether field
    !       is on tracer or velocity grid and how many levels this
    !       field has.

    rows          = pp_int(lbrow)
    columns       = pp_int(lbnpt)
    np            = pp_int(lblrec)
    len_extra     = MAX(pp_int(lbext), 0)  ! extra data in pp-field
    pp_int(lbext )= 0                      ! get rid of extra data

    DO i = 1, number_of_codes
      IF ( pp_int(item_code)  ==  stash_code(i) ) THEN

        !         Get grid and number of levels for this stash code
        tracer_grid = grid_of_tracer(i)
        nlevs_this_code = nlevs_code(i)
        WRITE(umMessage,*) ' PP Code ',pp_int(lbfc),' tracer_grid ',         &
        tracer_grid,' nlevs_this_code ',nlevs_this_code
        CALL umPrint(umMessage,src='anc_fld')

        !         Check field code set ; if not, set from FIELD_CODE
        IF (pp_int(lbfc) /= field_code(i)) THEN
          WRITE(umMessage,*) 'Field No ',fieldn,' Field code',               &
          ' incorrect or not set. Reset from ',pp_int(lbfc),                 &
          ' to ',field_code(i)
          CALL umPrint(umMessage,src='anc_fld')
          pp_int(lbfc) =  field_code(i)
        END IF

        GO TO 8
      END IF
    END DO

    WRITE(umMessage,*) ' WARNING from subroutine ANC_FLD'
    CALL umPrint(umMessage,src='anc_fld')
    WRITE(umMessage,*) ' Stash code ', pp_int(item_code),' in PP Header ',   &
      ' was not found in STASH_CODE of CODES namelist.'
    CALL umPrint(umMessage,src='anc_fld')

    8     CONTINUE

    !  1.4 Start loop over levels

    DO levn= 1, nlevs_this_code

      IF (levn  /=  1) THEN

        fieldn = fieldn + 1

        ! DEPENDS ON: read_pp_header
        CALL read_pp_header(ftin1,pp_int,pp_real,ibm_to_cray,l_bit_32)

        n_pp_flds(ftin1-29) = n_pp_flds(ftin1-29)+1

        WRITE(umMessage,*) 'Field No ',fieldn,' read in. From PP File '      &
       ,ftin1-29, ' Field No ',n_pp_flds(ftin1-29)
        CALL umPrint(umMessage,src='anc_fld')

        IF (pphead) THEN
          WRITE(umMessage,*) 'pp_int for field ',fieldn
          CALL umPrint(umMessage,src='anc_fld')
          DO j=1,45
            WRITE(umMessage,*) pp_int(j)
            CALL umPrint(umMessage,src='anc_fld')
          END DO
          WRITE(umMessage,*) 'pp_real for field ',fieldn
          CALL umPrint(umMessage,src='anc_fld')
          DO j=1,19
            WRITE(umMessage,*) pp_real(j)
            CALL umPrint(umMessage,src='anc_fld')
          END DO
        END IF

      END IF   ! levn  /=  1

      !  1.5 Set t_compress depending on wave and nlevs_this_code.
      !      Don't compress fields of only one level.

      IF (compress .AND. nlevs_this_code  /=  1) THEN

        t_compress = .TRUE.

      ELSE IF (wave .AND. pp_int(lbfc) == 38) THEN

        WRITE(umMessage,*)'re-setting compress false for lsmask'
        CALL umPrint(umMessage,src='anc_fld')
        WRITE(umMessage,*)'in ppheader section of anc_fld'
        CALL umPrint(umMessage,src='anc_fld')
        t_compress=.FALSE.

      ELSE IF (wave .AND.  pp_int(lbfc) /= 38                            &
                         .AND. nlevs_this_code  ==  1) THEN

        t_compress = .TRUE.

      ELSE

        t_compress = .FALSE.

      END IF

      !  1.6 Calculate fieldsize

      IF (add_wrap_pts) THEN
        IF (t_compress) THEN
          IF (.NOT. wave) THEN
            fieldsize= fldsizelev(levn)
          ELSE
            fieldsize=n_sea_points
            IF (pp_int(lbfc) == 38) THEN  ! LS mask data field
              WRITE(umMessage,*)'fieldsize set uncomp for lsmask'
              CALL umPrint(umMessage,src='anc_fld')
              fieldsize=rows*columns
            END IF
          END IF
        ELSE
          fieldsize=rows*(columns+2)
        END IF
      ELSE
        IF (t_compress) THEN
          IF (.NOT. wave) THEN
            fieldsize= fldsizelev(levn)
          ELSE
            fieldsize=n_sea_points
            IF (pp_int(lbfc) == 38) THEN  ! LS mask data field
              WRITE(umMessage,*)'fieldsize set uncomp for lsmask'
              CALL umPrint(umMessage,src='anc_fld')
              fieldsize=rows*columns
            END IF
          END IF
        ELSE
          fieldsize=rows*columns
        END IF
      END IF

      !  1.7 Set the fixed header, integer and real constants

      IF (fieldn == 1) THEN

        icode = 0

        !  Calculate no_cmp

        no_cmp = 0
        DO i = 1, nlevels   ! do not use levn in this loop
          no_cmp = no_cmp + fldsizelev(i)
        END DO

        !
        ! note - anc_head only uses compress if .not. wave
        !
        ! DEPENDS ON: anc_head
        CALL anc_head(pp_int,pp_real,rows,columns,fieldsize,len2_lookup,  &
         fld_types,n_times,nlevels,n_freq_waves,n_dir_waves,no_cmp,       &
         len1_levdepc,len2_levdepc,len1_rowdepc,len2_rowdepc,len1_coldepc,&
         len2_coldepc,len1_flddepc,len2_flddepc,len_extcnst,len_cfi,      &
         tracer_grid,add_wrap_pts,periodic,single_time,ibm_to_cray,       &
         t_compress,levdepc,rowdepc,coldepc,flddepc,extcnst,wave,         &
         lfh,fixhd,len_intc,int_const,len_realc,real_const,icode)

        IF (icode >  0) GO TO 9999  !  Error detected ; Return

        ! initialise header for length of data in ancillary file
        ! (reset to zero from mdi value set in anc_head)
        fixhd(161) = 0

      END IF    ! fieldn  ==  1

      ! accumulate indicator of length of data in ancillary file
      fixhd(161)=fixhd(161)+fieldsize

      !  1.8 Set the lookup table for this field

      !CMH note for waves - to use frequency information in lev-dep-consts
      !CMH need to set before call pp_table as well as in proper place.
      !CMH so do before the loops

      ! DEPENDS ON: pp_table
      CALL pp_table(pp_int,pp_real,len2_lookup,lookup,rlookup,          &
       fieldsize,fieldn,levn,m,runtot,number_of_codes,field_code,       &
       stash_code,add_wrap_pts,t_compress,pack32,wave,len1_levdepc,     &
       len2_levdepc,lev_dep_consts,len_realc,real_const,icode)

      IF (icode >  0) GO TO 9999  !  Error detected ; Return

      !  1.9 Read past the data part of this pp field
      l_skip = .TRUE.
      ! DEPENDS ON: readdata
      CALL readdata(rows,columns,ftin1,ibm_to_cray,len_extra,l_bit_32,  &
                    l_skip, dummy, dummy)

    END DO ! end of loop over levels

  END DO ! end of loop over fld_types

END DO ! end of loop over times

WRITE(umMessage,*) '==================================='
CALL umPrint(umMessage,src='anc_fld')
WRITE(umMessage,*) fieldn,' PP fields have been read in'
CALL umPrint(umMessage,src='anc_fld')
WRITE(umMessage,*) '==================================='
CALL umPrint(umMessage,src='anc_fld')

!  2.  If flddepc=t or compress = t:  Read in fields of constant and
!      compressed field indices from levels dataset

!  2.1 For ocean dumps, create compressed field indices and
!      fields_const

icode = 0

IF ((compress .OR. flddepc) .AND. .NOT. wave) THEN

  ! DEPENDS ON: calc_cfi_and_fld
  CALL calc_cfi_and_fld(ftin2,nlevels,len1_coldepc,               &
    cols_nowrap,len1_rowdepc,len1_flddepc,len2_flddepc,           &
    fields_const,fldsizelev,len_cfi,cfi1,cfi2,cfi3,compress,      &
    flddepc,ibm_to_cray,add_wrap_pts,imdi,l_bit_32,icode)

  IF (icode  /=  0) THEN
    GO TO 9999
  END IF

END IF

!
!  3. Over-ride elements in header arrays
!
!      Arrays that can be over-ridden are fixed length, integer
!      constants, real constants and level dependent constants.

!     * set default lev_dep_consts for WAVE frequency using the
!       factor 1.1
!     * as set in real_const(15) by namelist
!     * need to set frmin as lev_dep_consts(1) in namelist input
!     * pick up the CO value from namelist input HEADER value

!  3.1 Read in the header_data namelist which includes the
!     lev_dep_consts,row_dep_consts,col_dep_consts,extra_consts

REWIND(5)
! DEPENDS ON: find_namelist
i=find_namelist(5,"header_data")

IF (i == 0) THEN
  READ (UNIT=5, NML=header_data, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist HEADER_DATA", iomessage)
ELSE
  WRITE(umMessage,*)'Cannot find namelist HEADER_DATA'
  CALL umPrint(umMessage,src='anc_fld')
END IF
WRITE(umMessage,*) ' '
CALL umPrint(umMessage,src='anc_fld')

!  For wave dumps calculate the lev_dep_consts again

IF (wave) THEN

  DO m=2,len1_levdepc
    !CC      FR(M) = CO*FR(M-1)
    lev_dep_consts(m) = real_const(15)*lev_dep_consts(m-1)
  END DO

END IF

!  3.2 Amend the lookup tables

DO i = 1, len_look_user

  IF ( ifld_int(i)  /=  imdi ) THEN

    IF ( ifld_int(i)  ==  0 ) THEN
      DO j = 1, len2_lookup
        lookup( item_int(i) , j ) = ival_int(i)
      END DO
    ELSE
      lookup( item_int(i) , ifld_int(i) ) = ival_int(i)
    END IF

  END IF    ! ifld_int(i)  /=  imdi

  IF ( ifld_real(i)  /=  imdi ) THEN

    IF ( ifld_real(i)  ==  0 ) THEN
      DO j = 1, len2_lookup
        rlookup( item_real(i) , j ) = rval_real(i)
      END DO
    ELSE
      rlookup( item_real(i) , ifld_real(i) ) = rval_real(i)
    END IF

  END IF    ! ifld_real(i)  /=  imdi

END DO  ! i = 1, len_look_user

!  3.3 Print out headers to screen

IF (pphead) THEN

  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) 'fixhd'
  CALL umPrint(umMessage,src='anc_fld')
  DO j=1,161
    WRITE(umMessage,*) fixhd(j)
    CALL umPrint(umMessage,src='anc_fld')
  END DO
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) 'int_const ; length ',len_intc
  CALL umPrint(umMessage,src='anc_fld')
  DO j=1,len_intc
    WRITE(umMessage,*) int_const(j)
    CALL umPrint(umMessage,src='anc_fld')
  END DO
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='anc_fld')
  WRITE(umMessage,*) 'real_const ; length ',len_realc
  CALL umPrint(umMessage,src='anc_fld')
  DO j=1,len_realc
    WRITE(umMessage,*) real_const(j)
    CALL umPrint(umMessage,src='anc_fld')
  END DO

  IF (levdepc) THEN
    WRITE(umMessage,*) ' '
    CALL umPrint(umMessage,src='anc_fld')
    WRITE(umMessage,*) 'level dependent constants '
    CALL umPrint(umMessage,src='anc_fld')
    DO j=1,len2_levdepc
      ipos=(j-1)*len1_levdepc
      WRITE(umMessage,*) 'variable ',j
      CALL umPrint(umMessage,src='anc_fld')
      DO i=1,len1_levdepc
        WRITE(umMessage,*) lev_dep_consts(ipos+i)
        CALL umPrint(umMessage,src='anc_fld')
      END DO
    END DO
  END IF
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='anc_fld')

  IF (rowdepc) THEN
    WRITE(umMessage,*) ' '
    CALL umPrint(umMessage,src='anc_fld')
    WRITE(umMessage,*) 'row dependent constants '
    CALL umPrint(umMessage,src='anc_fld')
    DO j=1,len2_rowdepc
      ipos=(j-1)*len1_rowdepc
      WRITE(umMessage,*) 'variable ',j
      CALL umPrint(umMessage,src='anc_fld')
      DO i=1,len1_rowdepc
        WRITE(umMessage,*) row_dep_consts(ipos+i)
        CALL umPrint(umMessage,src='anc_fld')
      END DO
    END DO
  END IF
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='anc_fld')

  IF (coldepc) THEN
    WRITE(umMessage,*) ' '
    CALL umPrint(umMessage,src='anc_fld')
    WRITE(umMessage,*) 'column dependent constants '
    CALL umPrint(umMessage,src='anc_fld')
    DO j=1,len2_coldepc
      ipos=(j-1)*len1_coldepc
      WRITE(umMessage,*) 'variable ',j
      CALL umPrint(umMessage,src='anc_fld')
      DO i=1,len1_coldepc
        WRITE(umMessage,*) col_dep_consts(ipos+i)
        CALL umPrint(umMessage,src='anc_fld')
      END DO
    END DO
  END IF
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='anc_fld')

  IF (flddepc) THEN
    WRITE(umMessage,*) ' '
    CALL umPrint(umMessage,src='anc_fld')
    WRITE(umMessage,*) 'fields constants'
    CALL umPrint(umMessage,src='anc_fld')
    WRITE(umMessage,*)' len1_flddepc = ',len1_flddepc
    CALL umPrint(umMessage,src='anc_fld')
    WRITE(umMessage,*)' len2_flddepc = ',len2_flddepc
    CALL umPrint(umMessage,src='anc_fld')
  END IF
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='anc_fld')

  IF (extcnst) THEN
    WRITE(umMessage,*) ' '
    CALL umPrint(umMessage,src='anc_fld')
    WRITE(umMessage,*) 'extra constants '
    CALL umPrint(umMessage,src='anc_fld')
    DO i=1,len_extcnst
      WRITE(umMessage,*) extra_const(i)
      CALL umPrint(umMessage,src='anc_fld')
    END DO
  END IF
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='anc_fld')

END IF
!
!  4. Write out header data to ancillary field file
!

!  4.1 Write out fixed, integer and real constants headers.
!  Set values for use in WRITHEAD and convert lookup and rlookup
!  into one array lookup_all using subroutine conv_real

len_data=fixhd(161)

! DEPENDS ON: conv_real
CALL conv_real(rlookup,lookup_all,len2_lookup)

DO i=1,len2_lookup
  lookup_all(1:45,i) = lookup(1:45,i)
END DO

! If logical lwfio (set in namelist LOGICALS) is true then set the
! LBEGIN and LBNREC fields in the LOOKUP Headers for VN 16 Type
! Dumpfiles.
IF (lwfio) THEN

  ! DEPENDS ON: set_dumpfile_address
  CALL set_dumpfile_address(fixhd,lfh,                           &
   lookup_all,len1_lookup_all,len2_lookup,                       &
   number_of_data_words_in_memory,                               &
   number_of_data_words_on_disk,                                 &
   disk_address)

END IF

!  4.2 Use WRITHEAD to write the headers and constants

CALL writhead(ftout,fixhd,lfh,int_const,len_intc,                &
real_const,len_realc,lev_dep_consts,len1_levdepc,len2_levdepc,   &
row_dep_consts,len1_rowdepc,len2_rowdepc,col_dep_consts,         &
len1_coldepc,len2_coldepc,fields_const,len1_flddepc,             &
len2_flddepc,extra_const,len_extcnst,dumphist,len_dumphist,      &
cfi1,len_cfi(1),cfi2,len_cfi(2),cfi3,len_cfi(3),lookup_all,      &
len1_lookup_all,len2_lookup,len_data,                            &
umFortranIntegerSize()*8,                                        &
start_block,icode,cmessage)

!  5.0  Write out (rest of) data to ancillary file

!  5.1  Return to start of pp input pp fields files

WRITE(umMessage,*) ' '
CALL umPrint(umMessage,src='anc_fld')
DO n=1,n_pp_files
  REWIND 29+n
  WRITE(umMessage,*) ' Rewinding PP file on Unit No ',29+n
  CALL umPrint(umMessage,src='anc_fld')
END DO
WRITE(umMessage,*) ' '
CALL umPrint(umMessage,src='anc_fld')

!  5.2 Start loop over fields

fieldn=0   ! field number counter

DO n=1,n_times

  DO m=1,fld_types

    !  5.3 Do steps which are independant of level first

    IF (n_pp_files == 1) THEN
      ftin1= unit_no(1)
    ELSE
      IF (field_order) THEN
        ftin1= unit_no(n)
      ELSE
        ftin1= unit_no(m)
      END IF
    END IF

    fieldn = fieldn + 1

    !  5.4 Read pp header and determine length of field to
    !        be output to ancillary file

    ! DEPENDS ON: read_pp_header
    CALL read_pp_header (ftin1,pp_int,pp_real,ibm_to_cray,l_bit_32)

    ! (* extract the data dimensions and tracer/velocity grid type *)

    rows        = pp_int(lbrow)
    columns     = pp_int(lbnpt)
    len_extra   = MAX(pp_int(lbext), 0)  ! extra data in pp-field

    DO i = 1, number_of_codes
      IF ( pp_int(lbfc)  ==  field_code(i) ) THEN
        tracer_grid = grid_of_tracer(i)
        nlevs_this_code=nlevs_code(i)
        GO TO 30
      END IF
    END DO

    30    CONTINUE

    !  5.5 Start loop over levels

    DO levn = 1, nlevs_this_code

      IF (levn  /=  1) THEN
        fieldn = fieldn + 1
        ! DEPENDS ON: read_pp_header
        CALL read_pp_header(ftin1,pp_int,pp_real,ibm_to_cray,l_bit_32)
      END IF

      !  5.6 Set t_compress which depends on wave and nlevs_this_code

      IF (compress .AND. nlevs_this_code  /=  1) THEN

        t_compress = .TRUE.

      ELSE IF (wave .AND. pp_int(lbfc) == 38) THEN

        WRITE(umMessage,*)'re-setting compress false for lsmask'
        CALL umPrint(umMessage,src='anc_fld')
        WRITE(umMessage,*)'in ppheader section of anc_fld'
        CALL umPrint(umMessage,src='anc_fld')
        t_compress = .FALSE.

      ELSE

        t_compress = .FALSE.

      END IF

      !  5.7 Calculate fieldsize

      IF (add_wrap_pts) THEN
        IF (t_compress) THEN
          IF (.NOT. wave) THEN
            fieldsize= fldsizelev(levn)
          ELSE
            fieldsize=n_sea_points
            IF (pp_int(lbfc) == 38) THEN  ! LS mask data field
              WRITE(umMessage,*)'fieldsize set uncomp for lsmask'
              CALL umPrint(umMessage,src='anc_fld')
              fieldsize=rows*columns
            END IF
          END IF
        ELSE
          fieldsize=rows*(columns+2)
        END IF
      ELSE
        IF (t_compress) THEN
          IF (.NOT. wave) THEN
            fieldsize= fldsizelev(levn)
          ELSE
            fieldsize=n_sea_points
            IF (pp_int(lbfc) == 38) THEN  ! LS mask data field
              WRITE(umMessage,*)'fieldsize set uncomp for lsmask'
              CALL umPrint(umMessage,src='anc_fld')
              fieldsize=rows*columns
            END IF
          END IF
        ELSE
          fieldsize=rows*columns
        END IF

      END IF

      !  5.8 Call subroutine data to write the fields to the dump/ancillary
      !      file.

      ! DEPENDS ON: dataw
      CALL dataw(rows,columns,fieldsize,nlevels,levn,len_extra,         &
       fieldn,len1_lookup_all,lookup_all,fixhd,                         &
       len_cfi, cfi1, cfi2, cfi3,fldsizelev,ftin1,ftout,                &
       tracer_grid,add_wrap_pts,ibm_to_cray,t_compress,rmdi_input,      &
       wave,lsmask, l_bit_32,                                           &
       icode)

    END DO
  END DO
END DO

WRITE(umMessage,*) '========================================'
CALL umPrint(umMessage,src='anc_fld')
WRITE(umMessage,*) fieldn,' fields written to Ancillary File'
CALL umPrint(umMessage,src='anc_fld')
WRITE(umMessage,*) '========================================'
CALL umPrint(umMessage,src='anc_fld')

9999  CONTINUE
RETURN
END SUBROUTINE anc_fld

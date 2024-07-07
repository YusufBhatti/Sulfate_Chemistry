! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

PROGRAM pptoanc
USE io
USE filenamelength_mod, ONLY:                                    &
    filenamelength
USE errormessagelength_mod, ONLY: errormessagelength
USE check_iostat_mod
USE ereport_mod, ONLY: ereport
USE UM_Config, ONLY: &
    appInit,          &
    appTerminate,     &
    exe_pptoanc
USE lookup_addresses
USE submodel_mod, ONLY:                                                        &
    n_internal_model_max, n_internal_model
USE ppxlook_mod, ONLY: read_atmos_stashmaster

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE missing_data_mod, ONLY: rmdi

USE fort2c_interfaces, ONLY: get_file
USE get_env_var_mod,   ONLY: get_env_var
USE hostname_mod,      ONLY: get_hostname

IMPLICIT NONE

!
! Routine: pptoanc -------------------------------------------------
!
! Description:
!  To create ancillary fields from pp fields.
!  pp fields are output in same order as they are input
!
! Method:
!
!     unit                    Description
!    ftin1=30 onwards  INPUT  pp files (unit #s provided by user)
!    ftin2=11          INPUT  levels dataset (only used if compress=t
!                             or flddepc=t)
!    ftout=10          OUTPUT ancillary file
!
! Use  namelists to set:
!
!
!      sizes          field_types,n_times,n_levels,n_pp_files,
!                     n_freq_waves,n_dir_waves,stash_code,field_code,
!                     nlevs_code,unit_no,len_intc,len_realc,
!                     len1_levdepc,len2_levdepc,len1_rowdepc,
!                     len2_rowdepc,len1_coldepc,len2_coldepc,
!                     len1_flddepc,len2_flddepc,len_extcnst,rmdi_input
!
!      logicals       add_wrap_pts,periodic,single_time,ibm_to_cray,
!                     compress,wave,levdepc,rowdepc,coldepc,
!                     flddepc,extcnst,pack32,pphead,grid_of_tracer,
!                     field_order
!
!      first_vt       fvhh,fvdd,fvmm,fvyy
!                     (first validity time)
!
!      last_vt        lvhh,lvdd,lvmm,lvyy
!                     (last validity time)
!
!      interval       year360,year_unspec,ivhh,ivdd,ivmm,ivyy
!                     (time interval between validity times)
!
!      header_data    fixhd, int_const,real_const,lev_dep_consts,
!                     row_dep_consts,col_dep_consts,extra_const,
!                     ifld_int, item_int, ival_int,
!                     ifld_real, item_real, rval_real
!
!  The last 6 variables are arrays allowing up to n= len_look_user
!  changes to the lookup tables. All arrays are initiated as missing
!  data. Changes are made in the order from n=1 to len_look_user. If
!  ifld_int(n) or ifld_real(n) is not a missing data value, changes
!  are made to the integer or real lookup tables.
!
!                for integer parts of lookup tables
!                     ifld_int(n)   field number to change
!                                   0 => all lookup tables
!                     item_int(n)   item number to change
!                     ival_int(n)   integer value to use
!
!                for real parts of lookup tables
!                     ifld_real(n)   field number to change
!                                   0 => all lookup tables
!                     item_real(n)   item number to change
!                     rval_real(n)   real value to use
!
!    The following elements are particularly worth checking
!                     fixhd(4)    grid type code (global is default)
!                     fixhd(8)    calendar type (Gregorian is default)
!                     fixhd(12)   UM version number
!
!                     fixhd and lookup for dates of validity
!
!  The stash_code and field_code in the namelist sizes must be in the
!  same order as they are in the input pp fields
!
!  Do not put field types requiring different compression indices
!  in one ancillary file.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:

! Routine arguments
!   Scalar arguments

INTEGER :: n_stash_codes    ,                                        &
                           ! counter for number of stash codes
        n_unit_no        ,                                        &
                           ! counter for number of unit numbers
        len2_lookup_max  ,                                        &
                           ! 2nd dimension for lookup array
                           ! in ancfld(maximum)
        cols_nowrap      ,                                        &
                           ! no. of columns east-west without
                           ! wrap_points
        n,i              ,                                        &
                           ! loop counter
        icode
                           ! error exit condition code

!  Define variables from SIZES namelist

INTEGER :: field_types   ,                                           &
                         ! number of field types in I/O files
        n_times       ,                                           &
                         ! number of time periods in I/O files
        nlevels       ,                                           &
                         ! number of levels (default = 1)
        n_pp_files    ,                                           &
                         ! number of input pp files
        n_freq_waves  ,                                           &
                         ! number of wave frequencies
        n_dir_waves   ,                                           &
                         ! number of wave directions
        len_intc      ,                                           &
                         ! dimension for integer constants
        len_realc     ,                                           &
                         ! dimension for real constants
        len_extra     ,                                           &
                         ! dimension for extra data
        len1_levdepc  ,                                           &
                         ! dimension for levdepc array
        len2_levdepc  ,                                           &
                         ! 2nd dimension for levdepc array
        len1_rowdepc  ,                                           &
                         ! dimension for rowdepc array
        len2_rowdepc  ,                                           &
                         ! 2nd dimension for rowdepc array
        len1_coldepc  ,                                           &
                         ! dimension for coldepc array
        len2_coldepc  ,                                           &
                         ! 2nd dimension for coldepc array
        len1_flddepc  ,                                           &
                         ! dimension for flddepc array
        len2_flddepc  ,                                           &
                         ! 2nd dimension for flddepc array
        len_extcnst      ! dimension for extcnst array

REAL :: rmdi_input          ! real missing data indicator
                         ! in input pp field

!  Define variables from LOGICALS namelist

LOGICAL :: add_wrap_pts ,                                            &
                       ! T => adds wrapping columns
                       !      e.g. for global grid
        periodic     ,                                            &
                       ! T => periodic in time
                       !      e.g. climate field
        single_time  ,                                            &
                       ! T => all fields input valid at one time
        ibm_to_cray  ,                                            &
                       ! T => input pp data is in IBM number
                       !      format and needs to be converted to
                       !      run on the Cray.
                       !      (Only use if running on Cray)
        compress     ,                                            &
                       ! T => fields are packed into ancillary
                       !      field compressed field indices are
                       !      calculated
        wave         ,                                            &
                       ! T => a wave dump is to be created
        levdepc      ,                                            &
                       ! T => if level dependent constants array
                       !      required
        rowdepc      ,                                            &
                       ! T => if row dependant constants are
                       !      required
        coldepc      ,                                            &
                       ! T => if column dependant constants are
                       !      required
        flddepc      ,                                            &
                       ! T => if fields of constants are
                       !      required
        extcnst      ,                                            &
                       ! T => if fields of constants are
                       !      required

        pack32       ,                                            &
                       ! T => use 32 bit Cray numbers
        pphead       ,                                            &
                       ! T => print out pp headers read in

        field_order ,                                             &
                       ! T => input pp fields ordered by time.
                       !      i.e. different months in input
                       !        files, same fields in all files
                       ! F => inout pp fields ordered by fields.
                       !      i.e. different fields in input
                       !        files, all months in all files

        lwfio          ! T => set the LBEGIN and LBNREC fields
                       !      in the LOOKUP Headers for VN 16
                       !      Type Dumpfiles.
                       ! F => Old dumpfiles

CHARACTER(LEN=filenamelength) :: namelst
CHARACTER(LEN=errormessagelength) :: cmessage

! Parameters:
INTEGER :: ftin2                  ! input unit for mask file used
PARAMETER (ftin2=11)           ! for level dependent consts and
                               ! compression indices.
                               ! Only used when compress is T
INTEGER :: ftout
PARAMETER (ftout=10)           ! unit number for output ancillary
                               ! file
INTEGER :: nolevsmax
PARAMETER (nolevsmax=200)      ! max number of levels; dimensions
                               !  fldsizelev array
INTEGER :: number_of_codes
PARAMETER (number_of_codes=200)! max number of stash/field codes

INTEGER :: max_n_pp_files
PARAMETER (max_n_pp_files=50)  ! max number of input pp files

INTEGER, PARAMETER :: Max_Filename_Len = 256

! Array arguments:

INTEGER :: len_cfi(3)  ,                                             &
                           ! lengths of compressed field indices
        fldsizelev(nolevsmax)  ! size of packed field
                               ! on each level

LOGICAL :: grid_of_tracer(number_of_codes) ! T => fields are on a
                                        ! tracer grid

! Define variables from SIZES namelist

INTEGER :: stash_code(number_of_codes),                              &
                                    ! array of stash codes
        field_code(number_of_codes),                              &
                                    ! array of field codes
        nlevs_code(number_of_codes),                              &
                                    ! array of levels depending
                                    ! on field code
        unit_no(number_of_codes)    ! array of unit numbers for
                                    ! input

CHARACTER(LEN=Max_Filename_Len) :: pp_file_name
CHARACTER(LEN=1)                :: c_bit_32
CHARACTER(LEN=errormessagelength) :: iomessage
INTEGER                         :: err
INTEGER                         :: errorstatus
INTEGER                         :: bit_32
INTEGER                         :: length     ! length of returned string
LOGICAL                         :: L_bit_32   ! is this 32 bit?
! Function & Subroutine calls:
INTEGER :: find_namelist

!- End of header

NAMELIST /sizes/ field_types,n_times,nlevels,n_pp_files,          &
 n_freq_waves,n_dir_waves,stash_code,field_code,nlevs_code,       &
 unit_no,len_intc,len_realc,len1_levdepc,len2_levdepc,            &
 len1_rowdepc,len2_rowdepc,len1_coldepc,len2_coldepc,             &
 len1_flddepc,len2_flddepc,len_extcnst,rmdi_input

NAMELIST /logicals/ add_wrap_pts,periodic,single_time,            &
  ibm_to_cray,compress,wave,levdepc,rowdepc,coldepc,flddepc,      &
  extcnst,pack32,pphead,grid_of_tracer,field_order,lwfio


CALL appInit(exe_pptoanc)

! Write out hostname
WRITE(umMessage,'(A,A)') 'Host is ',TRIM(get_hostname())
CALL umPrint(umMessage, src='pptoanc')


!  1 Set values

!  1.0 Set default values for SIZES NAMELIST


field_types  = 2
n_times      = 12
nlevels      = 1
n_pp_files   = 1
n_freq_waves = 1
n_dir_waves  = 1
len_intc     = 40
len_realc    = 40
len1_levdepc = 1
len2_levdepc = 1
len1_rowdepc = 1
len2_rowdepc = 1
len1_coldepc = 1
len2_coldepc = 1
len1_flddepc = 1
len2_flddepc = 1
len_extcnst  = 1

rmdi_input   = rmdi

!  1.1 Initialise arrays in SIZES NAMELIST

DO n=1,number_of_codes
  field_code(n)=-99
  stash_code(n)=-99
  nlevs_code(n)=1
  unit_no(n)=-99
END DO

!  1.2 Open UNIT05 containing namelists and read in SIZES NAMELIST

CALL get_file(5,namelst,filenamelength,icode)
OPEN (UNIT=5,FILE=namelst, ACTION='READ')

REWIND(5)
! DEPENDS ON: find_namelist
i=find_namelist(5,"sizes")

IF (i == 0) THEN
  READ (UNIT=5, NML=sizes, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist SIZES", iomessage)
ELSE
  WRITE(umMessage,*)'Cannot find namelist SIZES'
  CALL umPrint(umMessage,src='pptoanc')
END IF

CALL print_nlist_sizes()

!  1.3 Check that n_pp_files is not greater than max_n_pp_files

IF (n_pp_files <= 0 .OR. n_pp_files >  max_n_pp_files) THEN
  CALL umPrint('',src='pptoanc')
  WRITE(umMessage,*) ' N_PP_FILES must in range 1-',max_n_pp_files
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) ' N_PP_FILES must in range 1-',number_of_codes
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) ' Resubmit job with new value for N_PP_FILES'
  CALL umPrint(umMessage,src='pptoanc')
  GO TO 9999   !  Return
END IF

!  1.4 Check that n_times and number of field_types is not greater
!      than number_of_codes

IF (n_times  >   number_of_codes   .OR.                           &
   field_types  >   number_of_codes ) THEN
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) ' ** WARNING ** WARNING ** '
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) ' N_TIMES = ',n_times,' or FIELD_TYPES = ',         &
      field_types,' greater than NUMBER_OF_CODES = ',number_of_codes
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) ' Dimension of UNIT_NO may be too small if used.'
  CALL umPrint(umMessage,src='pptoanc')
END IF

!  1.5 Count the number of stash codes and check they are not
!      greater than number of field_types

n_stash_codes = 0
DO n=1,number_of_codes
  IF (stash_code(n) >= 0) THEN
    n_stash_codes = n_stash_codes + 1
  END IF
END DO

IF (n_stash_codes /= field_types) THEN
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) ' Wrong number of stash codes provided.'
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) n_stash_codes,' stash codes in namelist.'
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) field_types  ,' stash codes expected.'
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) ' Rerun with correct no of stash codes'
  CALL umPrint(umMessage,src='pptoanc')
  GO TO 9999   !  Return
ELSE
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) n_stash_codes,' stash codes in SIZES namelist.'
  CALL umPrint(umMessage,src='pptoanc')
END IF

IF (nlevels  >   nolevsmax) THEN
  WRITE(umMessage,*) 'parameter nolevsmax is smaller than nlevels'
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) 'increase nolevsmax in program create'
  CALL umPrint(umMessage,src='pptoanc')
  GO TO 9999   !  Jump out
END IF

IF (rmdi_input  ==  rmdi) THEN
  WRITE(umMessage,*) 'rmdi_input should equal rmdi in input pp field'
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) 'WARNING !!! '
  CALL umPrint(umMessage,src='pptoanc')
  WRITE(umMessage,*) &
      'if not, RESUBMIT with the correct rmdi_input in SIZES namelist.'
  CALL umPrint(umMessage,src='pptoanc')
END IF

!
!  1.6 Set default values for LOGICALS NAMELIST
!
add_wrap_pts    = .FALSE.
periodic        = .FALSE.
single_time     = .FALSE.
ibm_to_cray     = .FALSE.
compress        = .FALSE.
wave            = .FALSE.
levdepc         = .FALSE.
rowdepc         = .FALSE.
coldepc         = .FALSE.
flddepc         = .FALSE.
extcnst         = .FALSE.
pack32          = .FALSE.
pphead          = .FALSE.
field_order     = .TRUE.
lwfio           = .TRUE.

!  1.7 Initialise array in LOGICAL NAMELIST

DO n = 1, number_of_codes
  grid_of_tracer(n)=.TRUE.
END DO

!  1.8 Read in the LOGICALS NAMELIST

REWIND(5)
! DEPENDS ON: find_namelist
i=find_namelist(5,"logicals")

IF (i == 0) THEN
  READ (UNIT=5, NML=logicals, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist LOGICALS", iomessage)
ELSE
  WRITE(umMessage,*)'Cannot find namelist LOGICALS'
  CALL umPrint(umMessage,src='pptoanc')
END IF

CALL print_nlist_logicals()

!  1.9 Count number of unit numbers needed which depends on field_order,
!      n_times and field_types

n_unit_no = 0
DO n=1,number_of_codes
  IF (unit_no(n) >  0) THEN
    n_unit_no = n_unit_no + 1
  END IF
END DO


IF (n_unit_no >  0) THEN
  DO n=1,n_unit_no
    IF (unit_no(n) <  30 .OR. unit_no(n) >  29+n_pp_files) THEN
      WRITE(umMessage,*) ' '
      CALL umPrint(umMessage,src='pptoanc')
      WRITE(umMessage,*) ' Unit no out of range in UNIT_NO :',unit_no(n)
      CALL umPrint(umMessage,src='pptoanc')
      WRITE(umMessage,*) ' Range is 30-',29+n_pp_files
      CALL umPrint(umMessage,src='pptoanc')
      WRITE(umMessage,*) ' Rerun with correct unit numbers'
      CALL umPrint(umMessage,src='pptoanc')
      GO TO 9999   !  Return
    END IF
  END DO
ELSE              ! n_unit_no
  DO n=1,max_n_pp_files
    unit_no(n)=29+n
  END DO
END IF

!  2 Set dimensions

!  2.0 If data are to be compressed calculate the lengths of compression
!     indices and number of points in field on each level using the
!     levels dataset

IF (add_wrap_pts) THEN
  cols_nowrap = len1_coldepc-2
ELSE
  cols_nowrap = len1_coldepc
END IF

IF (compress .AND. .NOT. wave) THEN  

  ! DEPENDS ON: calc_len_cfi
  CALL calc_len_cfi(ftin2,cols_nowrap,len1_rowdepc,              &
         nlevels,len_cfi,fldsizelev,ibm_to_cray,add_wrap_pts,    &
         l_bit_32,icode)

  IF (icode  /=  0) THEN
    GO TO 9999        ! jump out
  END IF

ELSE                 ! .not. compress

  len_cfi(1) = 1
  len_cfi(2) = 1
  len_cfi(3) = 1

END IF               ! compress

!
!  2.1 Calculate len2_lookup_max which depends on wave dimensions
!
icode = 0

IF (wave) THEN
  len2_lookup_max = field_types*n_times*nlevels                   &
   + (n_freq_waves*n_dir_waves -1)*n_times

ELSE
  len2_lookup_max = field_types*n_times*nlevels
END IF

WRITE(umMessage,*)'len2_lookup_max set to ',len2_lookup_max
CALL umPrint(umMessage,src='pptoanc')
WRITE(umMessage,*)' '
CALL umPrint(umMessage,src='pptoanc')
!
!  3 Read STASHmaster files

!  3.1 Initialise N_INTERNAL_MODEL

n_internal_model=n_internal_model_max

CALL read_atmos_stashmaster()

!  3.6 Find and open all pp files - check if 32 bit
!      data is required
DO i = 1, n_pp_files
  CALL Get_File(29+i, pp_file_name, Max_Filename_Len, err)
  pp_file_name = TRIM(pp_file_name)
  OPEN(UNIT=29+i, FILE=pp_file_name, FORM="unformatted", ACTION='READ', &
       CONVERT="big_endian" , IOSTAT=ErrorStatus, IOMSG=iomessage)
  IF (ErrorStatus /= 0) THEN
    cmessage ='Error opening file ' // TRIM(pp_file_name) //': '//&
              TRIM(iomessage)
    CALL ereport('PPTOANC', ErrorStatus, cmessage)
  END IF         
END DO

! Get the "32 bit" flag from environment
CALL get_env_var('BIT32', c_bit_32, allow_missing=.TRUE., allow_empty=.TRUE., &
                 length=length)
IF (length > 0) THEN
  c_bit_32 = TRIM(c_bit_32)
  READ(c_bit_32,'(i1)') bit_32
  IF (bit_32 == 1) THEN
    L_bit_32 = .TRUE.
  ELSE
    L_bit_32 = .FALSE.
  END IF
ELSE
  l_bit_32 = .FALSE.
END IF


!  4 Call main subroutine

! DEPENDS ON: anc_fld
CALL anc_fld(ftin2,ftout,nolevsmax,number_of_codes,               &
  max_n_pp_files,len_cfi,fldsizelev,                              &
  field_types,n_times,nlevels,n_pp_files,stash_code,field_code,   &
  nlevs_code,unit_no,n_freq_waves,n_dir_waves,len_intc,len_realc, &
  len_extra,len1_levdepc,len2_levdepc,len1_rowdepc,len2_rowdepc,  &
  len1_coldepc,len2_coldepc,len1_flddepc,len2_flddepc,            &
  len_extcnst,rmdi_input,                                         &
  add_wrap_pts,periodic,single_time,ibm_to_cray,compress,wave,    &
  levdepc,rowdepc,coldepc,flddepc,extcnst,pack32,pphead,          &
  grid_of_tracer,field_order,lwfio,L_bit_32,                      &
!    #  len2_lookup_max,cols_nowrap,icode)
        len2_lookup_max,cols_nowrap,icode)

IF (icode  >   0) THEN
  GO TO 9999   !  Jump out
END IF

!  5 Tidy up at end of program

!  5.1 Close ancillary file and pp files

CALL file_close (ftout,'ANCFILE',7,1,0,icode)
IF (icode  >   0) THEN
  WRITE(umMessage,*) ' Problem with FILE_CLOSE for unit no ',ftout
  CALL umPrint(umMessage,src='pptoanc')
  GO TO 9999   !  Jump out
END IF

DO i = i, n_pp_files
  CLOSE(i)
END DO


!     ===========================================================
!  5.2  NORMAL and ABNORMAL COMPLETION.
!     ===========================================================
9999  CONTINUE

IF (icode == 0) THEN
  CALL umPrint(' ',src='pptoanc')
  CALL umPrint('Program completed normally',src='pptoanc')
  WRITE(umMessage,'(A,I6)') 'Return code = ',icode
  CALL umPrint(umMessage,src='pptoanc')
  CALL umPrint('',src='pptoanc')
ELSE
  CALL umPrint( 'PPTOANC Program',src='pptoanc')
  CALL umPrint( 'Return code has been set in program',src='pptoanc')
  WRITE(umMessage,'(A,I6)') 'Return code = ',icode
  CALL umPrint(umMessage,src='pptoanc')
  CALL umPrint( ' ',src='pptoanc')
  CALL umPrint( 'Program aborted',src='pptoanc')
  CALL ereport('PPTOANC', icode, 'Return code error')
END IF


CALL appTerminate()

CONTAINS

SUBROUTINE print_nlist_sizes()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist sizes', &
    src='pptoanc')

WRITE(lineBuffer,*)' field_types = ',field_types
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' n_times = ',n_times
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' nlevels = ',nlevels
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' n_pp_files = ',n_pp_files
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' n_freq_waves = ',n_freq_waves
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' n_dir_waves = ',n_dir_waves
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' stash_code = ',stash_code
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' field_code = ',field_code
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' nlevs_code = ',nlevs_code
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' unit_no = ',unit_no
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' len_intc = ',len_intc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' len_realc = ',len_realc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' len1_levdepc = ',len1_levdepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' len2_levdepc = ',len2_levdepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' len1_rowdepc = ',len1_rowdepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' len2_rowdepc = ',len2_rowdepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' len1_coldepc = ',len1_coldepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' len2_coldepc = ',len2_coldepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' len1_flddepc = ',len1_flddepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' len2_flddepc = ',len2_flddepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' len_extcnst = ',len_extcnst
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' rmdi_input = ',rmdi_input
CALL umPrint(lineBuffer,src='pptoanc')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='pptoanc')

END SUBROUTINE print_nlist_sizes

SUBROUTINE print_nlist_logicals()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist logicals', &
    src='pptoanc')

WRITE(lineBuffer,*)' add_wrap_pts = ',add_wrap_pts
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' periodic = ',periodic
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' single_time = ',single_time
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' ibm_to_cray = ',ibm_to_cray
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' compress = ',compress
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' wave = ',wave
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' levdepc = ',levdepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' rowdepc = ',rowdepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' coldepc = ',coldepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' flddepc = ',flddepc
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' extcnst = ',extcnst
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' pack32 = ',pack32
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' pphead = ',pphead
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' grid_of_tracer = ',grid_of_tracer
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' field_order = ',field_order
CALL umPrint(lineBuffer,src='pptoanc')
WRITE(lineBuffer,*)' lwfio = ',lwfio
CALL umPrint(lineBuffer,src='pptoanc')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='pptoanc')

END SUBROUTINE print_nlist_logicals



END PROGRAM pptoanc

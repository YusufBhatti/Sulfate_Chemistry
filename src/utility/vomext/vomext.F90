! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! *********************************************************************
!  Program VOM_Extract : Extracts Vertical Profiles from UM data.
!
PROGRAM VOM_Extract
USE io
USE ereport_mod, ONLY: ereport
USE filenamelength_mod, ONLY:                                               &
    filenamelength
USE file_manager, ONLY: &
    assign_file_unit, release_file_unit, get_file_unit_by_id
USE UM_Config, ONLY:                                                        &
    appInit,                                                                 &
    appTerminate,                                                            &
    exe_vomext
USE umPrintMgr
USE errormessagelength_mod, ONLY: errormessagelength
USE hostname_mod,      ONLY: get_hostname
USE readflds_mod, ONLY: readflds
 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE
!
! Description :
!    Program to read in model data from UM Fieldsfiles and extrcat
!    required data to run the 3DVOM model.
!
! Method :
!    1. Namelist is read in - this will contain a list of lat/long of
!       points and times for which data is to be extracted for.
!       Maximums allowed : 50 profiles, 10 times.
!    2. The lookup headers for all the input fieldsfiles are read in.
!       Maximum of 20 fieldsfiles allowed.
!    3. For each lat/long the 'nearest' row/col on the model grid is
!       determined. Method used is very basic. The algorithm could be
!       improved to determine the nearest point more accurately. It
!       has only been tested with Global Model Level fieldsfiles. For
!       rotated grids, new code would be needed to convert the namelist
!       lat/long to the rotated lat/long.
!    4. Determine pointers to the data in the fieldsfiles. currently,
!       6 pointers :
!         - Orography               Stash Code 33
!         - U/V                     Stash Code 2/3
!         - Theta                   Stash Code 4
!         - Density                 Stash Code 253
!         - Pressure on Rho levels  Stash Code 407
!    5. Data is read in, extracted and written out.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
! Language : FORTRAN 90
! This code is written to UMDP3 programming standards.
!
! Declarations :

CHARACTER(LEN=errormessagelength) :: cmessage          ! Error Message
CHARACTER(LEN=*), PARAMETER :: RoutineName = "VOM_Extract"

INTEGER, PARAMETER :: max_n_times    = 10
INTEGER, PARAMETER :: max_n_profiles = 50
INTEGER, PARAMETER :: max_n_levs     = 71  ! ND has 70 levels for all fields
                                           ! EG 71 for theta level prognotics

!     ---------------
!     Namelist UMPROF
!     ---------------

REAL :: Latitude  (max_n_profiles) !  Latitude of profile
REAL :: Longitude (max_n_profiles) !  Longitude of profile
INTEGER :: FC_Time (max_n_times)   !  Forecast time data required
NAMELIST /umprof/ Latitude, Longitude, FC_Time

CALL appInit(exe_vomext)

! Write out hostname
WRITE(umMessage,'(A,A)') 'Host is ',TRIM(get_hostname())
CALL umPrint(umMessage, src='vomext')

CALL Timer ( RoutineName, 1 )
WRITE(umMessage,*) ' ###################################### '
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' Running VOM Extract Utility to extract '
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' vertical profiles from UM data         '
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' ###################################### '
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' '
CALL umPrint(umMessage,src='vomext')

CALL UM_Profiles()

CALL Timer ( RoutineName, 2 )

WRITE(umMessage,*) ' '
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' #######################################'
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' VOM Extract program completed normally.'
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' #######################################'
CALL umPrint(umMessage,src='vomext')

CALL appTerminate()

CONTAINS


!  Subroutine UM_Profiles

! Subroutine Interface :
SUBROUTINE UM_Profiles

USE ereport_mod, ONLY: ereport
USE missing_data_mod, ONLY: rmdi
USE get_env_var_mod, ONLY: get_env_var
IMPLICIT NONE
!
! Description :
!     Routine for extract UM data for 3DVOM model
!
! Method :
!     See main program.
!
!
! Language : FORTRAN 90
! This code is written to UMDP3 programming standards.
!
! Declarations :

!   Local variables

INTEGER :: icode                  !  Error code
INTEGER :: length                 !  Length of returned string
CHARACTER (LEN=errormessagelength) :: cmessage    !  Error Message
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'UM_Profiles'

INTEGER :: len_env    = 6  !  Length of env. variable
INTEGER :: env_var    = 0  !  Indicator that filename is in env var
INTEGER :: read_only  = 0  !  Input Fieldsfiles   - read only
INTEGER :: read_write = 1  !  Profile output file - read & write
CHARACTER(LEN=6) :: env            !  Env Variable for input fieldsfiles

INTEGER :: NL_unit_no         ! Unit No for Namelist File
INTEGER :: prof_unit_no       ! Unit No for Profile data
INTEGER :: FF_unit_no         ! Unit No for FFs (31-)

INTEGER :: Len_FixHd = 256
!     Integer :: FixHd(Len_FixHd)  ! Compiler doesn't like this line ?
INTEGER :: FixHd (256)

INTEGER :: Len_IHead = 10
INTEGER :: Len_RHead = 10
INTEGER :: IntHead (10)
REAL    :: RealHead(10)

INTEGER :: ReturnCode
INTEGER :: i,j,k,n

REAL :: orog  (max_n_profiles)
REAL :: u     (max_n_levs,max_n_profiles)
REAL :: v     (max_n_levs,max_n_profiles)
REAL :: p     (max_n_levs,max_n_profiles)
REAL :: rho   (max_n_levs,max_n_profiles)
REAL :: theta (max_n_levs,max_n_profiles)
CHARACTER (LEN=errormessagelength) :: iomessage
CHARACTER (LEN=filenamelength) :: FileName
!  FileName for Namelist and Profile output file

INTEGER :: Len1_Lookup
INTEGER :: Len2_Lookup
INTEGER :: Len_Lookup
INTEGER, ALLOCATABLE :: Lookup (:,:)
INTEGER, ALLOCATABLE :: Lookup_FF_No (:)

INTEGER :: n_lookups (0:20)
INTEGER :: LookUp_Start_Address(20)
INTEGER :: i_FF, n_FF
INTEGER :: ihead

REAL, ALLOCATABLE :: Data_In (:,:,:)
REAL, ALLOCATABLE :: srce_lat (:)
REAL, ALLOCATABLE :: srce_long(:)

REAL :: a_io
INTEGER :: len_io
REAL    :: prof_lat
REAL    :: prof_long
INTEGER :: iprof
INTEGER :: n_profiles
INTEGER :: n_times
INTEGER :: itime
INTEGER :: irow (max_n_profiles)
INTEGER :: icol (max_n_profiles)

INTEGER :: ipt_orog
INTEGER :: ipt_u
INTEGER :: ipt_v
INTEGER :: ipt_theta
INTEGER :: ipt_density
INTEGER :: ipt_rho_p

INTEGER :: ix
INTEGER :: iy

!     Validity Time
INTEGER :: vt_yy
INTEGER :: vt_mm
INTEGER :: vt_dd
INTEGER :: vt_hh

INTEGER :: row_length
INTEGER :: rows
INTEGER :: nlevs

! ENDGame input FF
LOGICAL :: l_EG_FF
INTEGER :: kstart  ! first level for reading theta data.


!   -------------------------
!   Get filename for Namelist
!   -------------------------

CALL get_env_var ('NAMELIST', FileName)

WRITE(umMessage,*) ' NameList File : ',FileName(1:LEN_TRIM(FileName))
CALL umPrint(umMessage,src='vomext')

!   ----------------------
!   Open the Namelist File
!   ----------------------

CALL assign_file_unit(FileName, NL_unit_no, handler="fortran")

OPEN ( UNIT=NL_unit_no, FILE=FileName,ACTION='READ',                       &
       IOSTAT= icode, IOMSG=iomessage )

WRITE(umMessage,*) ' Namelist file opened : icode ',icode
CALL umPrint(umMessage,src='vomext')

IF (icode /= 0) THEN
  CMessage = 'Error in opening namelist file: ' // TRIM(iomessage)

  CALL Ereport ( RoutineName, icode, CMessage )
END IF

!   -----------------
!   Read the Namelist
!   -----------------

Latitude (:)  = rmdi
Longitude (:) = rmdi
fc_time (:)   = -1

READ (NL_unit_no, NML=umprof)

CALL print_nlist_umprof()

!   -------------------
!   Close Namelist File
!   -------------------

WRITE(umMessage,*) ' closing namelist file'
CALL umPrint(umMessage,src='vomext')
CLOSE ( NL_unit_no, IOSTAT=icode, IOMSG=iomessage )
IF (icode /= 0) THEN
  CMessage = 'Error in closing namelist file: '// TRIM(iomessage)
  CALL Ereport ( RoutineName, Icode, Cmessage )
END IF
CALL release_file_unit(NL_unit_no, handler="fortran")

!   ------------------------------------
!   Work out how many profiles requested
!   ------------------------------------

n_profiles = 0
DO i = 1, max_n_profiles
  IF ( Latitude(i) /= rmdi ) THEN
    n_profiles = n_profiles + 1
  END IF
END DO
WRITE(umMessage,*) ' n_profiles ',n_profiles
CALL umPrint(umMessage,src='vomext')

!   ------------------------------------------
!   Work out how many forecast times requested
!   ------------------------------------------

n_times = 0
DO i = 1,max_n_times
  IF ( FC_time(i) /= -1 ) THEN
    n_times = n_times + 1
  END IF
END DO
WRITE(umMessage,*) ' n_times ',n_times
CALL umPrint(umMessage,src='vomext')

!   ------------------------------------
!   Get filename for profile output file
!   ------------------------------------

CALL get_env_var ('VOMFIL1', FileName)
WRITE(umMessage,*) ' Output Profile File : ',                              &
    FileName(1:LEN_TRIM(FileName))
CALL umPrint(umMessage,src='vomext')

!   ----------------------------
!   Open the profile output file
!   ----------------------------

!   Cannot use file_open, must use fortran open

CALL assign_file_unit(FileName, prof_unit_no, handler="fortran")

OPEN ( UNIT=prof_unit_no, FILE=FileName,                                   &
    STATUS='unknown', IOSTAT= icode, IOMSG=iomessage )

WRITE(umMessage,*) ' output file opened : icode ',icode
CALL umPrint(umMessage,src='vomext')

IF (icode /= 0) THEN
  CMessage = 'Error in opening profile output file: '// TRIM(iomessage)

  CALL Ereport ( RoutineName, icode, CMessage )
END IF

!   ----------------------------------
!   Work out how many FFs to be opened
!   ----------------------------------

n_FF = 0
DO i_FF = 1, 20

  env = 'FILE  '
  WRITE (env(5:6),'(I2.2)') i_FF

  FileName = 'NotSet'
  CALL get_env_var (env, FileName, allow_missing=.TRUE., length=length)
  IF (length < 1) THEN

    ! No such environment variable
    CYCLE

  ELSE

    ! Env Var found
    WRITE(umMessage,*) ' Env Var ',env,' found.'
    CALL umPrint(umMessage,src='vomext')
    WRITE(umMessage,*) ' Filename ',Filename
    CALL umPrint(umMessage,src='vomext')
    n_FF = n_FF + 1

  END IF

END DO
WRITE(umMessage,*) ' n_FF = ',n_FF
CALL umPrint(umMessage,src='vomext')

IF (n_FF == 0) THEN
  ICode = 13
  CMessage = 'No input FieldsFiles available?'

  CALL Ereport ( RoutineName, icode, CMessage )
END IF

!   --------------------
!   Open the fieldsfiles
!   --------------------

n_lookups (:) = 0

DO i_FF = 1, n_FF

  env = 'FILE  '
  WRITE (env(5:6),'(I2.2)') i_ff
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='vomext')

  !     Test env var exists

  CALL get_env_var (env, FileName)

  CALL assign_file_unit(FileName, ff_unit_no, handler="portio", id=env)

  CALL File_Open (ff_unit_no, filename, filenamelength,&
                  read_write=read_only,error=icode)
  IF (icode /= 0) THEN
    CMessage = 'Error in opening fieldfiles'

    CALL Ereport ( RoutineName, icode, CMessage )
  END IF

  !     ------------------------
  !     Read in the fixed header
  !     ------------------------

  ! DEPENDS ON: Read_FLH
  CALL Read_FLH ( ff_unit_no,fixhd,len_fixhd,icode,cmessage )
  IF (icode > 0) THEN
    CMessage = 'Error in Read_FLH'

    CALL Ereport ( RoutineName, icode, CMessage )
  END IF

  !     WRITE (6,*) ' Fixed Header for FFs '
  !     DO i=1,256
  !     WRITE (6,*) i,fixhd(i)
  !     END DO

  IF (fixhd(5) == 3 .AND. fixhd(9) == 3) THEN
    WRITE(umMessage,*) 'This is a ND FieldsFile'
    CALL umPrint(umMessage,src='vomext')
    l_EG_FF=.FALSE.
  ELSE IF (fixhd(5) == 3 .AND. fixhd(9) == 6) THEN
    WRITE(umMessage,*) 'This is an EG FieldsFile'
    CALL umPrint(umMessage,src='vomext')
    l_EG_FF=.TRUE.
  ELSE
    CMessage= 'Invalid input file type'
    Returncode = 10
    CALL Ereport( RoutineName, ReturnCode, CMessage )
  END IF

  !     -------------------------------
  !     Move to start of Integer Header
  !     -------------------------------

  CALL SetPos ( FF_unit_no, FixHd(100)-1, ICode )

  IF (ICode /= 0) THEN
    CMessage = 'Error in SetPos for Integer Header'

    CALL Ereport ( RoutineName, icode, CMessage )
  END IF

  !     ----------------------
  !     Read in Integer Header
  !     ----------------------

  CALL BuffIn ( FF_unit_no,                                                &
      IntHead, Len_IHead, Len_io, a_io )

  IF (a_io /= -1.0 .OR. len_io /= len_ihead) THEN
    ! DEPENDS ON: IOERROR
    CALL ioerror('buffer in of Integer header',a_io,len_io,                &
        len_lookup)
    cmessage = 'I/O ERROR with buffin of Integer Header'
    ICode    = 11

    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  !     ---------------------------------------
  !     Get grid dimensions from Integer Header
  !     ---------------------------------------

  row_length = IntHead (6)
  rows       = IntHead (7)
  nlevs      = IntHead (8)
  WRITE(umMessage,*) ' Model p grid : ',row_length,' x ',rows
  CALL umPrint(umMessage,src='vomext')
  WRITE(umMessage,*) ' No of model plevels in FFs ',nlevs
  CALL umPrint(umMessage,src='vomext')

  !     ----------------------------
  !     Move to start of Real Header
  !     ----------------------------

  CALL SetPos ( FF_unit_no, FixHd(105)-1, ICode )

  IF (ICode /= 0) THEN
    CMessage = 'Error in SetPos for Real Header'

    CALL Ereport ( RoutineName, icode, CMessage )
  END IF

  !     -------------------
  !     Read in Real Header
  !     -------------------

  CALL BuffIn ( FF_unit_no,                                                &
      RealHead, Len_RHead, Len_io, a_io )

  IF (a_io /= -1.0 .OR. len_io /= len_rhead) THEN
    ! DEPENDS ON: IOERROR
    CALL ioerror('buffer in of Real header',a_io,len_io,                   &
        len_lookup)
    cmessage = 'I/O ERROR with buffin of Real Header'
    ICode    = 11

    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  !     ---------------------------------------
  !     Get lookup dimensions from fixed header
  !     ---------------------------------------

  LookUp_Start_Address(i_FF) = FixHd(150)
  Len1_Lookup = FixHd(151)
  Len2_Lookup = FixHd(152)
  Len_Lookup  = Len1_Lookup * Len2_Lookup
  WRITE(umMessage,*) ' len1_lookup ',len1_lookup,                          &
      ' len2_lookup ',len2_lookup,                                         &
      ' len_lookup  ',len_lookup
  CALL umPrint(umMessage,src='vomext')

  !     -------------------------------
  !     Allocate space for LookUp table
  !     -------------------------------

  ALLOCATE ( LookUp (Len1_Lookup, Len2_Lookup) )

  !     -----------------------------
  !     Move to start of LookUp table
  !     -----------------------------

  CALL SetPos ( FF_unit_no, FixHd(150)-1, ICode )

  IF (ICode /= 0) THEN
    CMessage = 'Error in SetPos for LookUp table'

    CALL Ereport ( RoutineName, icode, CMessage )
  END IF

  !     --------------------
  !     Read in LookUP table
  !     --------------------

  CALL BuffIn ( FF_unit_no,                                                &
      LookUp, Len_Lookup, Len_io, a_io )

  IF (a_io /= -1.0 .OR. len_io /= len_lookup) THEN
    ! DEPENDS ON: IOERROR
    CALL ioerror('buffer in of lookup header',a_io,len_io,                 &
        len_lookup)
    cmessage = 'I/O ERROR with buffin of Lookup Header'
    ICode    = 11

    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  !     -----------------------------
  !     Find number of LookUp headers
  !     -----------------------------

  DO i = 1, FixHd(152)
    IF ( Lookup (42,i) == -99 ) THEN
      EXIT
    ELSE
      n_lookups(i_FF) = n_lookups(i_FF) +1
    END IF
  END DO

  !     Update total no of LookUp Entries
  n_lookups(0) = n_lookups(0) + n_lookups(i_FF)

  WRITE(umMessage,*) ' n_lookups ',n_lookups(i_FF)
  CALL umPrint(umMessage,src='vomext')

  !     Deallocate LookUp Table for File i_FF
  DEALLOCATE (LookUp)

END DO   !  Loop over i_FF

WRITE(umMessage,*) ' total n_lookups ',n_lookups(0)
CALL umPrint(umMessage,src='vomext')

!   ----------------------------------------
!   Total number of LookUp entries now known
!   Allocate space for all LookUp entries
!   ----------------------------------------

ALLOCATE ( LookUp (64, n_lookups(0) ) )
ALLOCATE ( LookUp_FF_No (n_lookups(0)) )

ihead = 0

DO i_FF = 1, n_FF

  env = 'FILE  '
  WRITE (env(5:6),'(I2.2)') i_ff
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='vomext')

  ff_unit_no = get_file_unit_by_id(env, handler="portio")

  !     -----------------------------
  !     Move to start of LookUp table
  !     -----------------------------

  CALL SetPos ( FF_unit_no, LookUp_Start_Address(i_FF)-1, ICode )

  IF (ICode /= 0) THEN
    CMessage = 'Error in SetPos for LookUp table'

    CALL Ereport ( RoutineName, icode, CMessage )
  END IF

  !     --------------------
  !     Read in LookUP table
  !     --------------------

  Len_LookUp = 64 * n_lookups(i_FF)

  CALL BuffIn ( FF_unit_no,                                                &
      LookUp (1:,ihead+1:), Len_Lookup, Len_io, a_io )

  IF (a_io /= -1.0 .OR. len_io /= len_lookup) THEN
    ! DEPENDS ON: IOERROR
    CALL ioerror('buffer in of lookup header',a_io,len_io,                 &
        len_lookup)
    cmessage = 'I/O ERROR with buffin of Lookup Header'
    ICode    = 11

    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  DO i = 1, n_lookups (i_FF)
    LookUp_FF_No (ihead+i) = ff_unit_no
  END DO

  ihead = ihead + n_lookups (i_FF)

END DO

!   -----------------------------
!   Allocate space for input data
!   -----------------------------
! Here we allocate Data_In array to be large enough to read in the
! the largest grids in the horizontal and vertical planes.
! The v grid has a different number of rows to the theta and u grids.
! The theta grid for EG has a zeroth level not present in ND.

ALLOCATE ( Srce_Lat (rows) )
ALLOCATE ( Srce_Long(row_length) )

! EG u and theta rows [1:rows], v [1,rows+1]
! EG theta levels [0:nlevs]
! ND u and theta rows [1:rows], v [1,rows-1]
! ND theta levels [1:nlevs]

  ALLOCATE ( Data_In  (row_length, rows, nlevs) )

!   -----------------------------------------------
!   Compute icol and jrow from lat/long in namelist
!   -----------------------------------------------

! ASIDE: The following attempts to calculate which row and col pt is 'near'
! to requested latitude and longitude.
! It is flawed in that an assumption is made that row x and col y
! are indeed available. But close to the North Pole this may not be the case.
! If one requests ND North Pole data, row=rows, we do not have suitable v data.
! While for EG it will be impossible to request data close to the North Pole
! as the p grid does not encompass either Pole.
! To fix these problems implies a major rewrite of this utitility including
! reading in BZY, BDY, BZX and BDX for each field from the lookup headers
! and interpolation between the p,u and v grids to ensure we have values for
! the whole domain.

WRITE(umMessage,*) ' Source Grid details'
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' No of columns   ',IntHead(6)
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' No of rows      ',IntHead(7)
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' deltax          ',RealHead(1)
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' deltay          ',RealHead(2)
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' first latitude  ',RealHead(3)
CALL umPrint(umMessage,src='vomext')
WRITE(umMessage,*) ' first longitude ',RealHead(4)
CALL umPrint(umMessage,src='vomext')

DO iprof = 1, n_profiles

  WRITE(umMessage,*) ' latitude  in namelist ',latitude(iprof)
  CALL umPrint(umMessage,src='vomext')
  WRITE(umMessage,*) ' longitude in namelist ',longitude(iprof)
  CALL umPrint(umMessage,src='vomext')

  prof_lat  = latitude (iprof)
  prof_long = longitude(iprof)

  !       WRITE (6,*) ' prof_lat  ',prof_lat
  !       WRITE (6,*) ' prof_long ',prof_long

  DO i = 1, Inthead(6)
    srce_long(i) = RealHead(4) + (i-1)*RealHead(1)
  END DO

  DO j = 1, IntHead(7)
    srce_lat(j) = RealHead(3) + (j-1)*RealHead(2)
  END DO

  ix = 0
  DO i = 1, IntHead(6)
    !         if (i < 10) WRITE (6,*) i,srce_long(i),prof_long
    IF ( srce_long(i) <= prof_long ) THEN
      ix = i
    ELSE
      EXIT
    END IF
  END DO
  !       WRITE (6,*) ' ix ',ix,' long ',srce_long(ix)

  iy = 0
  DO j = 1, IntHead(7)
    IF ( srce_lat(j) <= prof_lat ) THEN
      iy = j
    ELSE
      EXIT
    END IF
  END DO
  !       WRITE (6,*) ' iy ',iy,' lat ',srce_lat(iy)

  icol (iprof) = ix
  irow (iprof)  = iy
  WRITE(umMessage,*) ' icol set to ',icol (iprof)
  CALL umPrint(umMessage,src='vomext')
  WRITE(umMessage,*) ' irow set to ',irow (iprof)
  CALL umPrint(umMessage,src='vomext')

END DO

!   ---------------
!   Loop over times
!   ---------------

DO itime = 1, n_times

  !     ----------------------------
  !     Find pointers to data fields
  !     ----------------------------

  IF (itime == 1) THEN
    ipt_orog = 0
    DO i = 1, n_lookups (0)
      IF ( Lookup(42,i) == 33 .AND.                                        &
          LookUp(14,i) == FC_Time(itime) ) THEN
        ipt_orog = i
        EXIT
      END IF
    END DO
  END IF

  ipt_u = 0
  DO i = 1, n_lookups (0)
    IF ( Lookup(42,i) == 2 .AND.                                           &
        LookUp(14,i) == FC_Time(itime) ) THEN
      ipt_u = i
      EXIT
    END IF
  END DO

  ipt_v = 0
  DO i = 1, n_lookups (0)
    IF ( Lookup(42,i) == 3 .AND.                                           &
        LookUp(14,i) == FC_Time(itime) ) THEN
      ipt_v = i
      EXIT
    END IF
  END DO

  ipt_theta = 0
  DO i = 1, n_lookups (0)
    IF ( Lookup(42,i) == 4 .AND.                                           &
        LookUp(14,i) == FC_Time(itime) ) THEN
      ipt_theta = i
      EXIT
    END IF
  END DO

  ipt_density = 0
  DO i = 1, n_lookups (0)
    IF ( Lookup(42,i) == 253 .AND.                                         &
        LookUp(14,i) == FC_Time(itime) ) THEN
      ipt_density = i
      EXIT
    END IF
  END DO

  ipt_rho_p = 0
  DO i = 1, n_lookups (0)
    IF ( Lookup(42,i) == 407 .AND.                                         &
        LookUp(14,i) == FC_Time(itime) ) THEN
      ipt_rho_p = i
      EXIT
    END IF
  END DO

  WRITE(umMessage,*) ' ipt_orog    ',ipt_orog
  CALL umPrint(umMessage,src='vomext')
  WRITE(umMessage,*) ' ipt_u       ',ipt_u
  CALL umPrint(umMessage,src='vomext')
  WRITE(umMessage,*) ' ipt_v       ',ipt_v
  CALL umPrint(umMessage,src='vomext')
  WRITE(umMessage,*) ' ipt_theta   ',ipt_theta
  CALL umPrint(umMessage,src='vomext')
  WRITE(umMessage,*) ' ipt_density ',ipt_density
  CALL umPrint(umMessage,src='vomext')
  WRITE(umMessage,*) ' ipt_rho_p   ',ipt_rho_p
  CALL umPrint(umMessage,src='vomext')

  !     ----------------------------------------------
  !     Extract validity time from LookUp header for U
  !     ----------------------------------------------

  IF (ipt_u > 0) THEN
    vt_yy = LookUp (1, ipt_u)
    vt_mm = LookUp (2, ipt_u)
    vt_dd = LookUp (3, ipt_u)
    vt_hh = LookUp (4, ipt_u)
  END IF
  WRITE(umMessage,*) ' VT : ',vt_yy,vt_mm,vt_dd,vt_hh
  CALL umPrint(umMessage,src='vomext')

  !     ---------------------
  !     Get FieldFiles Number
  !     ---------------------

  FF_Unit_No = LookUp_FF_No (ipt_u)
  WRITE(umMessage,*) ' FF_Unit_No : ',FF_Unit_No
  CALL umPrint(umMessage,src='vomext')

  !     -----------------
  !     Read in orography
  !     -----------------

  IF (itime == 1) THEN   !   Read in orography for first time only

    WRITE(umMessage,*) ' Reading in orography'
    CALL umPrint(umMessage,src='vomext')

    CALL ReadFlds(FF_unit_no,                                              &
        1,                                                                 &
        ipt_orog,                                                          &
        lookup,                                                            &
        data_in(:,:,:),                                                 &
        fixhd,                                                             &
        1,                                                                 &
        icode,                                                             &
        CMessage)

    IF (icode /= 0) THEN
      WRITE(umMessage,'(A)') 'Error in ReadFlds ?'
      CALL umPrint(umMessage,src='vomext')

      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

    !     -----------------
    !     Extract Orography
    !     -----------------

    DO iprof = 1, n_profiles
      orog (iprof) = data_in ( icol(iprof), irow(iprof), 1 )
      WRITE(umMessage,*) ' orog ',orog(iprof)
      CALL umPrint(umMessage,src='vomext')
    END DO

  END IF

  !     --------------
  !     Read in U data
  !     --------------

  WRITE(umMessage,*) ' Reading in U'
  CALL umPrint(umMessage,src='vomext')

    CALL ReadFlds(FF_unit_no,                                              &
        nlevs,                                                             &
        ipt_u,                                                             &
        lookup,                                                            &
        data_in(:,:,:),                                                 &
        fixhd,                                                             &
        1,                                                                 &
        icode,                                                             &
        CMessage)

    IF (icode /= 0) THEN
      WRITE(umMessage,'(A)') 'Error in ReadFlds ?'
      CALL umPrint(umMessage,src='vomext')

      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

  !     ----------------
  !     Extract U data
  !     ----------------

  DO iprof = 1, n_profiles
    DO k=1,nlevs
      u(k,iprof) = data_in ( icol(iprof), irow(iprof), k )
      !       WRITE (6,*) ' u data ',k,u(k)
    END DO
  END DO

  !     --------------
  !     Read in V data
  !     --------------

  WRITE(umMessage,*) ' Reading in V'
  CALL umPrint(umMessage,src='vomext')

  IF (l_eg_ff) THEN
    DEALLOCATE (Data_In)
    ALLOCATE ( Data_In  (row_length, rows+1, nlevs) )
  ELSE
    DEALLOCATE (Data_In)
    ALLOCATE ( Data_In  (row_length, rows-1, nlevs) )
  END IF

    CALL ReadFlds(FF_unit_no,                                              &
        nlevs,                                                             &
        ipt_v,                                                             &
        lookup,                                                            &
        data_in(:,:,:),                                                 &
        fixhd,                                                             &
        1,                                                                 &
        icode,                                                             &
        CMessage)

    IF (icode /= 0) THEN
      WRITE(umMessage,'(A)') 'Error in ReadFlds ?'
      CALL umPrint(umMessage,src='vomext')

      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

  !     ----------------
  !     Extract V data
  !     ----------------

  DO iprof = 1, n_profiles
    DO k=1,nlevs
      v(k,iprof) = data_in ( icol(iprof), irow(iprof), k )
      !       WRITE (6,*) ' v data ',k,v(k)
    END DO
  END DO

  !     --------------
  !     Read in P data
  !     --------------

  WRITE(umMessage,*) ' Reading in P'
  CALL umPrint(umMessage,src='vomext')

  DEALLOCATE (Data_In)
  ALLOCATE ( Data_In  (row_length, rows, nlevs) )

    CALL ReadFlds(FF_unit_no,                                              &
        nlevs,                                                             &
        ipt_rho_p,                                                         &
        lookup,                                                            &
        data_in(:,:,:),                                                 &
        fixhd,                                                             &
        1,                                                                 &
        icode,                                                             &
        CMessage)

    IF (icode /= 0) THEN
      WRITE(umMessage,'(A)') 'Error in ReadFlds ?'
      CALL umPrint(umMessage,src='vomext')

      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

  !     --------------
  !     Extract P data
  !     --------------

  DO iprof = 1, n_profiles
    DO k=1,nlevs
      p(k,iprof) = data_in ( icol(iprof), irow(iprof), k )
      !       WRITE (6,*) ' p data ',k,p(k)
    END DO
  END DO

  !     -------------
  !     Read in Theta
  !     -------------

  WRITE(umMessage,*) ' Reading in theta'
  CALL umPrint(umMessage,src='vomext')

  kstart=1
  IF (l_eg_ff) THEN
    kstart=0
    DEALLOCATE (Data_In)
    ALLOCATE ( Data_In  (row_length, rows, 0:nlevs) )
  END IF

    CALL ReadFlds(FF_unit_no,                                              &
        nlevs-kstart+1,                                                    &
        ipt_theta,                                                         &
        lookup,                                                            &
        data_in(:,:,kstart:),                                            &
        fixhd,                                                             &
        1,                                                                 &
        icode,                                                             &
        CMessage)

    IF (icode /= 0) THEN
      WRITE(umMessage,'(A)') 'Error in ReadFlds ?'
      CALL umPrint(umMessage,src='vomext')

      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

  !     ------------------
  !     Extract Theta data
  !     ------------------
  ! output 3dvom expects only 1 to nlevs.

  DO iprof = 1, n_profiles
    DO k=1,nlevs
      theta(k,iprof) = data_in ( icol(iprof), irow(iprof), k )
      !       WRITE (6,*) ' theta data ',k,theta(k)
    END DO
  END DO

  !     ---------------
  !     Read in density
  !     ---------------

  WRITE(umMessage,*) ' Reading in density'
  CALL umPrint(umMessage,src='vomext')

  IF (l_eg_ff) THEN
    DEALLOCATE (Data_In)
    ALLOCATE ( Data_In  (row_length, rows, nlevs) )
  END IF

    CALL ReadFlds(FF_unit_no,                                              &
        nlevs,                                                             &
        ipt_density,                                                       &
        lookup,                                                            &
        data_in(1:,1:,1:),                                                 &
        fixhd,                                                             &
        1,                                                                 &
        icode,                                                             &
        CMessage)

    IF (icode /= 0) THEN
      WRITE(umMessage,'(A)') 'Error in ReadFlds ?'
      CALL umPrint(umMessage,src='vomext')

      CALL Ereport ( RoutineName, icode, CMessage )
    END IF


  !     --------------------
  !     Extract Density data
  !     --------------------

  DO iprof = 1, n_profiles
    DO k=1,nlevs
      rho(k,iprof) = data_in ( icol(iprof), irow(iprof), k )
      !       WRITE (6,*) ' rho data ',k,rho(k)
    END DO
  END DO

  !     --------------------------
  !     Write out the profile data
  !     --------------------------

  DO i = 1, n_profiles  !  Points

    WRITE (prof_unit_no,1)                                                 &
        'YEAR MONTH DAY HOUR         FORECAST TIME (H)'
    WRITE (prof_unit_no,2) vt_yy,vt_mm,vt_dd,vt_hh,fc_time(itime)
    WRITE (prof_unit_no,3)                                                 &
        'LATITUDE LONGITUDE  ROW  COLUMN   HEIGHT  LEVELS'
    WRITE (prof_unit_no,4)                                                 &
        latitude(i),longitude(i),irow(i),icol(i),orog(i),nlevs
    WRITE (prof_unit_no,5)                                                 &
        ' Level    U         V        Theta      P        Rho'

    DO k = 1, nlevs     !  Levels

      !         WRITE (6,6)&
      !         z(k),u(k),v(k),theta(k),p(k),rho(k)&
      !         k,u(k),v(k),theta(k),p(k),rho(k)&
      !         k,u(k,i),v(k,i),theta(k,i),p(k,i),rho(k,i)

      WRITE (prof_unit_no,6)                                               &
          k,u(k,i),v(k,i),theta(k,i),p(k,i),rho(k,i)

    END DO

  END DO

END DO    !  End of loop over forecast times

!   -------------------------
!   Deallocate space for data
!   -------------------------

DEALLOCATE (Data_In)
DEALLOCATE (Srce_Lat)
DEALLOCATE (Srce_Long)

DEALLOCATE (LookUp)
DEALLOCATE (LookUp_FF_No)

!   ---------------------
!   Close the fieldsfiles
!   ---------------------

WRITE(umMessage,*) ' closing fieldsfiles'
CALL umPrint(umMessage,src='vomext')
DO i_FF = 1, n_FF

  env = 'FILE  '
  WRITE (env(5:6),'(I2.2)') i_ff
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='vomext')

  ff_unit_no = get_file_unit_by_id(env, handler="portio")

  CALL File_Close (ff_unit_no,env,len_env,env_var,0,icode)
  IF (icode /= 0) THEN
    CMessage = 'Error in closing fieldsfiles'
    CALL Ereport ( RoutineName, Icode, Cmessage )
  END IF
  CALL release_file_unit(ff_unit_no, handler="portio")
END DO

!   ---------------------
!   Close the output file
!   ---------------------

!   Cannot use file_close, use fortran close

CLOSE ( prof_unit_no, IOSTAT=icode, IOMSG=iomessage )
IF (icode /= 0) THEN
  CMessage = 'Error in closing output file: '// TRIM(iomessage)
  CALL Ereport ( RoutineName, Icode, Cmessage )
END IF
CALL release_file_unit(prof_unit_no, handler="fortran")

!   ----------------------------------------------------
!   These FORMAT statements must not be changed without
!   making corresponding changes in the 3DVOM model
!   ----------------------------------------------------

1   FORMAT (a45)
2   FORMAT (i4,1x,i2,4x,i2,2x,i2,10x,i2)
3   FORMAT (a48)
4   FORMAT (f7.3,2x,f7.3,2x,i5,2x,i5,2x,f8.3,2x,i3)
5   FORMAT (a52)
6   FORMAT (i3,1x,2(f10.4,1x),1x,f8.3,1x,f10.2,2x,f17.2)

RETURN

END SUBROUTINE UM_Profiles

SUBROUTINE print_nlist_umprof()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist umprof',                              &
    src='vomext')

WRITE(lineBuffer,*)' Latitude = ',Latitude
CALL umPrint(lineBuffer,src='vomext')
WRITE(lineBuffer,*)' Longitude = ',Longitude
CALL umPrint(lineBuffer,src='vomext')
WRITE(lineBuffer,*)' FC_Time = ',FC_Time
CALL umPrint(lineBuffer,src='vomext')

CALL umPrint('- - - - - - end of namelist - - - - - -',                  &
    src='vomext')

END SUBROUTINE print_nlist_umprof

END PROGRAM VOM_Extract

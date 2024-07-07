! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Initialise SCM diagnostic output files

SUBROUTINE dump_streams_init                                                  &
  ( SCMop, row_length, rows, bl_levels                                        &
  , cloud_levels, ozone_levels, st_levels, sm_levels, ntiles, year_init       &
  , month_init, day_init, hour_init, min_init, sec_init, timestep, ndayin     &
  , nminin, nsecin, sec_day, tot_nsteps, a_sw_radstep_prog                    &
  , a_sw_radstep_diag, z_top_of_model, first_constant_r_rho_level             &
  , eta_theta, eta_rho, orog, r_theta_levels, r_rho_levels, netcdf_chunksize )

USE netcdf
USE UM_types

USE scm_utils, ONLY:                                                        &
    scm_timestep, scm_timestep_count, time_info, scm_nc_time_id             &
  , scm_nc_date_id, scm_nc_hour_id, scm_nc_local_id

USE nlsizes_namelist_mod, ONLY: model_levels
! SCMop_type is defined in here...
USE scmoptype_defn
USE s_scmop_mod,                                                            &
    ONLY: StreamIsOn, NotWritten
USE um_version_mod, ONLY: um_version_char
USE umPrintMgr,  ONLY: umPrint, umMessage, newline
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE


! Description:
!   Perform all steps necessary for subsequent writing of output stream data :
!   write the headers of the data files and any auxilliary files that may
!   accompany them

! Method:

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran77 with some bits of Fortran90

! NOTE: on the NEC this routine may be compiled in 32 bits because requests
!       to create 32 bit integers (necessary to communicate with the NetCDF
!       library routines - see incdf) are ignored by the compiler when
!       compiling at 64 bits. Thus all input and output variables are
!       explicitly declared as 64 bit using the KIND types i64, r64 and l64
!       delcared in scmoptype_defn module

TYPE(SCMop_type) :: SCMop ! InOut The derived-type structure containing
                          !       all the diagnostic information

! In Horizontal and vertical model domain...
INTEGER(i64) :: &
  row_length    &
, rows          &
, bl_levels     &
, cloud_levels  &
, ozone_levels  &
, st_levels     &
, sm_levels     &
, ntiles

INTEGER(i64) :: &
  year_init     &! In Initial year
, month_init    &! In Initial month
, day_init      &! In Initial day
, hour_init     &! In Initial hour
, min_init      &! In Initial minute
, sec_init       ! In Initial second

REAL(r64) ::    &
  timestep       ! In Timestep

INTEGER(i64) ::     &
  ndayin            &! In No. of days requested in run
, nminin            &! In No. of minutes requested in run
, nsecin            &! In No. of seconds requested in run
, sec_day           &! In No. of seconds in a day (Why this
                     !    isn't in an include file I don't know)
, tot_nsteps        &! In Total no. of steps to be done
, a_sw_radstep_prog &! In No. of timesteps between prognostic
                     !    calls to radiation (3C/3Z)
, a_sw_radstep_diag  ! In No. of timesteps between diagnostic
                     !    calls to radiation (3C/3Z)

REAL(r64) ::        &
   z_top_of_model    ! In The height of the "top of the
                     !    atmosphere"

INTEGER(i64) ::     &! In Lowest rho level that has constant height
   first_constant_r_rho_level

REAL(r64) ::                  &
   eta_theta(model_levels+1)  &! In The etas of the theta and rho levels
 , eta_rho(model_levels)      &!
 , orog(row_length,rows)       ! In Orography height

! IN The physical heights of the theta and rho levels...
REAL(r64) ::                                      &
   r_theta_levels(row_length,rows,0:model_levels) &
 , r_rho_levels(row_length,rows,model_levels)

! IN The ChunkSize input to NF90_Create(), the routine that
! creates a NetCDF file. Controls a space versus time trade-off:
! memory allocated in the netcdf library versus number of system
! calls.
INTEGER(i64) :: netcdf_chunksize

CHARACTER (LEN=100) :: FMT        ! Holds format specifiers
CHARACTER (LEN=100) :: basetime   ! Time units

INTEGER ::               &
  i, n, m                &! Counters
, nhourin                 ! No. of hours requested in run

INTEGER(i64) ::          &! 64-bit counter for passing to other
  n64                     ! routines compiled in 64 bit

INTEGER :: UNIT

INTEGER(i64) ::          &
  diags(SCMop%nentries)  &
, ndiags

REAL :: ndump

! Copy of netcdf_chunksize at the correct precision for passing
! into NetCCDF routines. May be altered by NF90_Create.
INTEGER(incdf) :: netcdf_chunksize_io

! NetCDF identifiers for each dimension of each domain
! profile. Integers that are passed into NetCDF library routines
! should have the same precision as those used internally within the
! library (which may have been compiled at 32 bit), hence the use of
! incdf.
INTEGER (incdf) ::                 &
  netcdf_dimension_time            &
, netcdf_dimension(3,maxndomprof)

! Other NetCDF-format related variables
INTEGER (incdf) :: STATUS,NcID
CHARACTER (LEN=5) :: c_unit

! Character function that incorporates the substep number into
! the short name of a diagnostic. Somewhat longer than a normal
! short name.
CHARACTER (LEN=lsname+10) :: add_substep_to_sname

! Temporary variables to hold output from add_substep_to_sname
CHARACTER (LEN=lsname+10) :: short_name
CHARACTER (LEN=lsname+10), ALLOCATABLE :: short_names(:)
CHARACTER (LEN=100) :: start_date_str

! For calls to ereport
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage


! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DUMP_STREAMS_INIT'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Write out the UM version number
WRITE(umMessage,'(A)')                                                        &
  '==================='//                                            newline//&
  'UM Version='//um_version_char//                                   newline//&
  '==================='
CALL umPrint(umMessage,src='dump_streams_init')

! Give a summary of the output streams
IF (SCMop%n_output /= 0) THEN

  WRITE(umMessage,'(A)')                                                      &
    '|-----------------------------------------------------------------'    //&
    newline // '| Summary of open streams:' // newline                      //&
    '|-----------------------------------------------------------------'    //&
    newline                                                                 //&
    '|No.| Unit | Time between dumps | No. dgnstcs | Format | Filename'     //&
    newline                                                                 //&
    '|   |      | (steps) |(minutes) |             |        |'              //&
    newline                                                                 //&
    '|---+------+---------+----------+-------------+--------+----------'
  CALL umPrint(umMessage,src='dump_streams_init')


  DO i=1, maxnstreams
    IF (SCMop%strm(i)%switch /= 0) THEN

      ! Write the stream output unit into a character variable,
      ! unless this is a NetCDF stream, in which case the unit
      ! is not used.
      IF (SCMop%strm(i)%FORMAT /= 4) THEN
        WRITE(c_unit,'(I5)') SCMop%strm(i)%op_unit
      ELSE
        c_unit='n/a'
        c_unit=ADJUSTR(c_unit)
      END IF

      WRITE(umMessage,'(A)')                                                  &
            '|   |      |         |          |             |        |'

      WRITE(umMessage(3:4),   '(I2)')    i
      WRITE(umMessage(7:11),  '(A5)')    c_unit
      WRITE(umMessage(15:21), '(I7)')    SCMop%strm(i)%dump_step
      WRITE(umMessage(25:32), '(F8.02)') SCMop%strm(i)%dump_step*timestep/60.0
      WRITE(umMessage(35:46), '(I12)')   SCMop%strm(i)%n_output
      WRITE(umMessage(49:55), '(I7)')    SCMop%strm(i)%FORMAT
      WRITE(umMessage(58:),   '(A)')     TRIM(SCMop%strm(i)%filename)

      CALL umPrint(ADJUSTL(umMessage),src='dump_streams_init')

    END IF
  END DO

  WRITE(umMessage,'(A)')                                                      &
    ' |-----------------------------------------------------------------'   //&
    newline  
  CALL umPrint(umMessage,src='dump_streams_init')


ELSE
  WRITE(umMessage,'(A)') 'Dump_Streams_Init: No diagnostics to be written.'
  CALL umPrint(umMessage,src='dump_streams_init')

  ! Nothing to do if there's no diagnostics to be written
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN

END IF

! Write the scumlist file
n64 = 0
! DEPENDS ON: write_scumlist
CALL write_scumlist (SCMop, n64)

! write the no. of levels, etc into a file
! DEPENDS ON: write_domain_data
CALL write_domain_data                                           &
  ( row_length, rows, model_levels, z_top_of_model               &
  , first_constant_r_rho_level, orog, eta_theta,eta_rho )

! Calculate nhourin from nminin and nsecin
nhourin = INT(nminin/60.0+nsecin/3600.0)

! Write the headers of the output data files
DO n=1, maxnstreams
  IF (SCMop%strm(n)%switch /= 0 .AND.                             &
      SCMop%strm(n)%n_output >  0) THEN
    ! This stream is open and there are diagnostics
    ! to be written out from it
    UNIT = SCMop%strm(n)%op_unit

    ! Check this unit is free (not applicable for NetCDF files)
    DO m=1, n-1
      IF (SCMop%strm(m)%switch  /= 0 .AND. SCMop%strm(m)%FORMAT  /= 4 .AND.   &
          SCMop%strm(m)%n_output > 0 .AND. SCMop%strm(m)%op_unit == UNIT) THEN
        WRITE(cmessage,'(A,3(I3),A)')                                         &
          '=========================================================='      //&
          newline                                                           //&
          '| Problem in dump_streams_init:' // newline                      //&
          '| Two output files with the same unit numbers: ',n,m,UNIT,         &
          newline                                                           //&
          '=========================================================='      //&
          newline                                                           //&
          'Output redirected to unit 103!'

        ! Raise a warning
        icode = -600
        CALL ereport( routinename, icode, cmessage )

        UNIT = 103
        SCMop%strm(n)%op_unit = UNIT
      END IF
    END DO

    ! The header depends on the format

    !-------------------------------------------------------------------
    ! Format=0: designed for reading by PV-wave routine
    !           scmread2.pro
    !-------------------------------------------------------------------
    IF (SCMop%strm(n)%FORMAT == 0) THEN

      ! Open the stream's output file
      OPEN(UNIT=unit,FILE=SCMop%strm(n)%filename)

      ! Write the UM version and the stream format
      WRITE(UNIT,'(A)')'vn'//um_version_char//',old_format'

      ! Write the total no. of days, hours, dumps per day and
      ! steps, plus the timestep, the dumping period of this
      ! stream and the no. of diagnostics per dump
      ndump = sec_day/(SCMop%strm(n)%dump_step*timestep)
      WRITE(UNIT,                                              &
           '(I4,1X,I4,1X,F7.3,1X,I7,1X,F9.2,1X,I9,1X,I4)')     &
            ndayin, nhourin, ndump, tot_nsteps, timestep,      &
            SCMop%strm(n)%dump_step,                           &
            SCMop%strm(n)%n_output

      ! Make a list the diagnostics we're going to write
      ndiags = 0
      DO i=1, SCMop%nentries
        IF (StreamIsOn(SCMop%streams(i),n) .AND. .NOT.           &
            NotWritten(SCMop%streams(i))) THEN

          ndiags = ndiags+1
          diags(ndiags) = i
        END IF
      END DO

      IF (ndiags /= SCMop%strm(n)%n_output) THEN
        WRITE(cmessage,'(A,4(I4))')                            &
          'Dump_Streams_Init ERROR: an inconsistency '       //&
          'has ocurred !', ndiags, SCMop%strm(n)%n_output,     &
           n, SCMop%nentries
        ! Raise a fatal error
        icode = 345
        CALL ereport( routinename, icode, cmessage )
      END IF

      ! Write that list to the file
      WRITE(FMT,'(A,I5,A)') '(',ndiags,'(I5,1X))'
      WRITE(UNIT,FMT)(SCMop%sname_id(diags(i)),i=1,ndiags)

      !-------------------------------------------------------------------
      ! Format=1 : designed for easy visual perusal.
      !-------------------------------------------------------------------
    ELSE IF (SCMop%strm(n)%FORMAT == 1) THEN

      ! Open the stream's output file
      OPEN(UNIT=unit,FILE=SCMop%strm(n)%filename)

      WRITE(UNIT,'(A,A)')'UM version=',um_version_char

      ! Write the total no. of days, hours, dumps per day and
      ! steps, plus the timestep, the dumping period of this
      ! stream and the no. of diagnostics per dump
      ndump = sec_day/(SCMop%strm(n)%dump_step*timestep)
      WRITE(UNIT,'(A)')'General run information...'
      WRITE(UNIT,'(A,I4,1X,A,I4,1X,A,F7.3,1X,A,I7)')                  &
           'No. of days=',ndayin,'No. of hrs=',nhourin,               &
           'No. of dumps=',ndump,'Total number of steps=',            &
            tot_nsteps
      WRITE(UNIT,'(A,F9.2,1X,A,I9,1X,A,I4)')                          &
           'Timestep=',timestep,                                      &
           'No. of steps between dumps=',                             &
           SCMop%strm(n)%dump_step,                                   &
           'No. of diags outputting to this stream=',                 &
           SCMop%strm(n)%n_output

      ! list the diagnostics we're going to write
      ndiags = 0
      DO i=1, SCMop%nentries
        IF (StreamIsOn(SCMop%streams(i),n) .AND. .NOT.                  &
          NotWritten(SCMop%streams(i))) THEN

          ndiags = ndiags+1
          diags(ndiags) = i
        END IF
      END DO

      IF (ndiags /= SCMop%strm(n)%n_output) THEN
        WRITE(cmessage,'(A,4(I4))')                            &
          'Dump_Streams_Init ERROR: an inconsistency '       //&
          'has ocurred !', ndiags, SCMop%strm(n)%n_output,     &
           n, SCMop%nentries
        ! Raise a fatal error
        icode = 345
        CALL ereport( routinename, icode, cmessage )
      END IF

      WRITE(UNIT,'(A)')'Diagnostics in this file... '//               &
           '(subset of entries in scumlist file)'

      ! Ensure that the sub-step is appended to the short
      ! name of each diagnostic entry, if required.

      ALLOCATE(short_names(ndiags))

      DO i=1, ndiags
        ! DEPENDS ON: add_substep_to_sname
        short_names(i) = add_substep_to_sname(SCMop,diags(i))
      END DO ! i

      ! Write the information
      WRITE(UNIT,'(I3,1X,A,1X,A)')                                    &
           (SCMop%sname_id(diags(i)),short_names(i),                  &
            SCMop%lname(diags(i)),i=1,ndiags)

      DEALLOCATE(short_names)

      !-------------------------------------------------------------------
      ! Format=2 : format required for FSSI database
      !-------------------------------------------------------------------
    ELSE IF (SCMop%strm(n)%FORMAT == 2) THEN

      ! Initialisation done in DUMP_STREAMS for the time being.

      !-------------------------------------------------------------------
      ! Format=3 : new format designed for reading by PV-wave
      !            routine scmread2.pro
      !-------------------------------------------------------------------
    ELSE IF (SCMop%strm(n)%FORMAT == 3) THEN

      ! Open the stream's output file
      OPEN(UNIT=unit,FILE=SCMop%strm(n)%filename)

      ! Write the UM version, plus any extra strings necessary
      ! for PV-Wave to know what to expect in the file.
      WRITE(UNIT,'(A)')'vn'//um_version_char//',start_time_present'

      ! Write information about the size and shape of the domain
      WRITE(UNIT,'(I4,I4,8I4)') row_length,rows,                      &
           model_levels,model_levels,bl_levels,cloud_levels,          &
           ozone_levels,st_levels,sm_levels,ntiles

      ! Write the starting time of the run as specified in
      ! the INDATA namelist
      WRITE(UNIT,'( 6(I0,1x) )')                                      &
           year_init, month_init, day_init,                           &
           hour_init, min_init,   sec_init

      ! Write information pertaining to the length of
      ! the model run and time in general
      WRITE(UNIT,'(I9,1x,F9.2,1x,I9,1x,I4,1x,I4,1x,I5,1x,I5,1x,I5)')  &
           tot_nsteps,timestep,SCMop%strm(n)%dump_step,               &
           a_sw_radstep_prog,a_sw_radstep_diag,ndayin,nminin,nsecin

      ! Write information about the levels and orography...
      WRITE(FMT,'(A,I3,A)') '(',model_levels+1 ,'(1PE18.10E3," "))'
      WRITE(UNIT,FMT) eta_theta

      WRITE(FMT,'(A,I3,A)') '(',model_levels   ,'(1PE18.10E3," "))'
      WRITE(UNIT,FMT) eta_rho

      WRITE(FMT,'(A,I6,A)') '(',row_length*rows,'(1PE18.10E3," "))'
      WRITE(UNIT,FMT) orog

      WRITE(FMT,'(A,I6,A)') '(',row_length*rows*(model_levels+1)      &
                           ,'(1PE18.10E3," "))'
      WRITE(UNIT,FMT) r_theta_levels

      WRITE(FMT,'(A,I6,A)') '(',row_length*rows*(model_levels  )      &
                           ,'(1PE18.10E3," "))'

      WRITE(UNIT,FMT) r_rho_levels

      ! Write data into the file pertaining to the
      ! individual diagnostics that it will contain.

      n64 = n
      ! DEPENDS ON: write_scumlist
      CALL write_scumlist(SCMop,n64)

      !-------------------------------------------------------------------
      ! Format=4 : NetCDF
      !-------------------------------------------------------------------
    ELSE IF (SCMop%strm(n)%FORMAT == 4) THEN

      ! Create the NetCDF file. Chunksize is a parameter,
      ! the size of which may greatly affect the I/O speed
      ! of the NetCDF file writing.
      netcdf_chunksize_io = netcdf_chunksize

      STATUS = Nf90_Create(SCMop%strm(n)%filename,                    &
      Nf90_Clobber,NcID,Chunksize = netcdf_chunksize_io)

      ! Write the start-date to a string for outputting as a global attribute
      WRITE(start_date_str,'( 6(I0,1x) )')                            &
           year_init, month_init, day_init,                           &
           hour_init, min_init,   sec_init

      ! Set global attributes
      STATUS = nf90_put_att                  &
             ( ncid   = ncid                 &
             , varid  = nf90_global          &
             , NAME   = "Start date"         &
             , values = TRIM(start_date_str) )

      ! Declare the three spatial dimensions of each domain
      ! profile
      DO i=1, maxndomprof

        ! Check this domain has been defined.
        IF (SCMop%d_lev1(i) /= -1) THEN
          STATUS = Nf90_Def_Dim(                                      &
            ncid  = NcID,                                             &
            NAME  = TRIM(ADJUSTL(SCMop%d_name(i)))//'_i',             &
            LEN   = INT(SCMop%d_rowa2(i)-SCMop%d_rowa1(i)+1,incdf),   &
            dimid = netcdf_dimension(1,i))

          STATUS = Nf90_Def_Dim(                                      &
            ncid  = NcID,                                             &
            NAME  = TRIM(ADJUSTL(SCMop%d_name(i)))//'_j',             &
            LEN   = INT(SCMop%d_rowb2(i)-SCMop%d_rowb1(i)+1,incdf),   &
            dimid = netcdf_dimension(2,i))

          STATUS = Nf90_Def_Dim(                                      &
            ncid  = NcID,                                             &
            NAME  = TRIM(ADJUSTL(SCMop%d_name(i)))//'_k',             &
            LEN   = INT(SCMop%d_lev2(i)-SCMop%d_lev1(i)+1,incdf),     &
            dimid = netcdf_dimension(3,i))
        END IF
      END DO

      WRITE(basetime,"(a,i4.4,a,i2.2,a,i2.2,1x,i2.2,a,i2.2,a,i2.2)")  &
        "seconds since ", year_init, "-", month_init, "-", day_init   &
                        , hour_init, ":", min_init,   ":", sec_init

      ! Declare the time dimension
      STATUS = Nf90_Def_Dim               &
        ( ncid   = NcID                   &
        , NAME   = "time"                 &
        , LEN    = nf90_unlimited         &
        , dimid  = netcdf_dimension_time  )

      ! Define a seconds variable
      STATUS = nf90_def_var               &
        ( ncid   = ncid                   &
        , NAME   = "time"                 &
        , xtype  = nf90_float             &
        , dimids = netcdf_dimension_time  &
        , varid  = scm_nc_time_id         )

      ! Define seconds variable attributes
      STATUS = nf90_put_att               &
        ( ncid   = ncid                   &
        , NAME   = "units"                &
        , varid  = scm_nc_time_id         &
        , values = basetime               )

      STATUS = nf90_put_att               &
        ( ncid   = ncid                   &
        , NAME   = "long_name"            &
        , varid  = scm_nc_time_id         &
        , values = "Number of seconds since start of simulations" )

      ! Define a date variable
      STATUS = nf90_def_var               &
        ( ncid   = ncid                   &
        , NAME   = "date"                 &
        , xtype  = nf90_int               &
        , dimids = netcdf_dimension_time  &
        , varid  = scm_nc_date_id         )

      ! Define date variable attributes
      STATUS = nf90_put_att               &
        ( ncid   = ncid                   &
        , NAME   = "units"                &
        , varid  = scm_nc_date_id         &
        , values = "yyyymmdd"             )

      STATUS = nf90_put_att               &
        ( ncid   = ncid                   &
        , NAME   = "long_name"            &
        , varid  = scm_nc_date_id         &
        , values = "Date string"          )


      ! Define each diagnostic that will go to this stream
      DO i=1, SCMop%nentries

        ! Is this entry to be sent to this stream, and
        ! not been flagged as one not to output?
        IF (StreamIsOn(SCMop%streams(i),n) .AND. .NOT.                  &
            NotWritten(SCMop%streams(i))) THEN

          ! Get the short name with the substep number
          ! appended if necessary.
          n64 = i
          short_name = add_substep_to_sname(SCMop,n64)

          ! Define this variable
          STATUS = Nf90_Def_Var(                                      &
            ncid   = NcID,                                            &
            NAME   = short_name,                                      &
            xtype  = Nf90_Float,                                      &
            dimids = (/                                               &
              netcdf_dimension(1,SCMop%domprof(i)),                   &
              netcdf_dimension(2,SCMop%domprof(i)),                   &
              netcdf_dimension(3,SCMop%domprof(i)),                   &
              netcdf_dimension_time/),                                &
            varid  = SCMop%netcdf_id(i))

          IF (STATUS /= nf90_NoErr) THEN
            WRITE(cmessage,'(A,A)')                                   &
              'Error defining NetCDF variable: ', short_name
            ! Raise a fatal error
            icode = 456
            CALL ereport( routinename, icode, cmessage )
          END IF

          ! Define its attributes
          STATUS = Nf90_Put_Att(                                      &
            ncid   = NcID,                                            &
            varid  = SCMop%netcdf_id(i),                              &
            NAME   = "units",                                         &
            values = SCMop%units(i))

          STATUS = Nf90_Put_Att(                                      &
            ncid   = NcID,                                            &
            varid  = SCMop%netcdf_id(i),                              &
            NAME   = "long_name",                                     &
            values = SCMop%lname(i))
        END IF
      END DO

      ! End of file definition
      STATUS = Nf90_EndDef(ncid = NcID)

      ! Record the file ID so it can be used later for output.
      SCMop%strm(n)%op_unit = NcID

    ELSE
      WRITE(cmessage,'(A,I3,A)')                                              &
        ' =========================================================='//       &
                                                                     newline//&
        ' | ERROR in dump_streams_init:'//                           newline//&
        ' | Unknown format for output stream ',n,                    newline//&
        ' =========================================================='
      icode = 567
      CALL ereport( routinename, icode, cmessage )
    END IF

  END IF                  ! SCMop%strm(n)%switch /= 0
END DO                     ! do n=1,maxnstreams

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE dump_streams_init


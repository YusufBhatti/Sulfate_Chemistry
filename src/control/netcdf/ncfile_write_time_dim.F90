! *****************************COPYRIGHT*******************************
! 
! Copyright 2017-2018 University of Reading
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice, 
! this list of conditions and the following disclaimer in the documentation 
! and/or other materials provided with the distribution.
! 
! 3. Neither the name of the copyright holder nor the names of its contributors 
! may be used to endorse or promote products derived from this software without 
! specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: NetCDF output
!
!  Purpose: Create the netCDF dimensions, coordinate variables, and bounds
!           variables, plus associated attributes needed for time coordinates.
!           Update cell_methods attributes for STASH variable.


MODULE ncfile_write_time_dim_mod

USE netcdf,                 ONLY: nf90_max_name ! External netCDF library module
USE yomhook,                ONLY: lhook,dr_hook
USE parkind1,               ONLY: jprb,jpim
USE umprintmgr,             ONLY: umprint,ummessage,printstatus,prdiag
USE ereport_mod,            ONLY: ereport
USE file_manager,           ONLY: um_file_type
USE submodel_mod,           ONLY: atmos_im
USE model_time_mod,         ONLY: stepim
USE nlstgen_mod,            ONLY: &
    steps_per_periodim,secs_per_periodim,dumpfreqim
USE nlstcall_mod,           ONLY: model_basis_time,lcal360
USE stash_array_mod,        ONLY: stlist,sttabl,nsttims
USE cv_param_mod,           ONLY: a_conv_step
USE jules_vegetation_mod,   ONLY: phenol_period,triffid_period
USE river_inputs_mod,       ONLY: river_step
USE asad_mod,               ONLY: interval
USE conversions_mod,        ONLY: isec_per_min,isec_per_hour,isec_per_day
USE profilename_length_mod, ONLY: profilename_length
USE errormessagelength_mod, ONLY: errormessagelength
USE reinit_file_times_mod,  ONLY: reinit_file_times
USE rad_input_mod,          ONLY: &
    a_lw_radstep_diag,a_sw_radstep_diag, &
    a_lw_radstep_prog,a_sw_radstep_prog
USE stparam_mod,            ONLY: &
    st_diag_address,st_freq_code,st_start_time_code,st_end_time_code, &
    st_accum_code,st_time_mean_code,st_max_code,st_min_code, &
    st_end_of_list,st_replace_code,st_output_ts0_code, &
    st_proc_no_code,st_input_code,st_time_unit2_code,st_time_unit3_code
USE cstash_mod,             ONLY: &
    itim_b,timpro, &
    stsh_timesteps,stsh_seconds,stsh_minutes,stsh_hours,stsh_days,stsh_dumps
USE umnetcdf_mod,           ONLY: &
    nc_create_dim,nc_create_var,nc_put_att, &
    check_cf_name,nc_dim_type_real,nc_max_attlen

IMPLICIT NONE

! The time units which can appear in the 'units' and 'cell_methods' attributes
CHARACTER(LEN=*), PARAMETER :: timeunit(4) = ['seconds','minutes', &
                                              'hours  ','days   ']

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NCFILE_WRITE_TIME_DIM_MOD'

PRIVATE add_interval
PRIVATE iunit_from_ts

CONTAINS

SUBROUTINE ncfile_write_time_dim(um_file,time_avail,section,ix, &
                                 cell_methods,dimid)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER         , INTENT(IN)    :: time_avail  ! Time availability code
INTEGER         , INTENT(IN)    :: section     ! STASH section number
INTEGER         , INTENT(IN)    :: ix          ! STASH item index
CHARACTER(LEN=*), INTENT(INOUT) :: cell_methods! CF cell_methods attribute
INTEGER         , INTENT(OUT)   :: dimid       ! NetCDF dimension ID

!
!  Local variables
!
INTEGER :: idiag              ! STASH diagnostic address
INTEGER :: input_code         ! When -ve gives index of parent record in stlist
INTEGER :: proccode           ! Time processing code for data
INTEGER :: sampling_freq      ! Sampling period of data
INTEGER :: starttime          ! Start time of diagnostic output
INTEGER :: endtime            ! End time of diagnostic output
INTEGER :: output_freq        ! Frequency of diagnostic output
INTEGER :: times_list(nsttims)! Used to store output times if given as a list
INTEGER :: ntimes_list        ! Number of times stored in times_list
INTEGER :: step               ! Timestep of model
INTEGER :: river_step_ts      ! River routing period in timesteps
INTEGER :: triffid_period_ts  ! Triffid period in timesteps
INTEGER :: phenol_period_ts   ! Phenol period in timesteps
INTEGER :: dimid_array(1)     ! used to convert dimid value into 1-D array
INTEGER :: dimid_bnds         ! dimension ID for bounds variable
INTEGER :: dimid_time_bnds(2) ! dimension IDs for time bounds variable
INTEGER :: varid              ! variable ID for time coordinates
INTEGER :: varid_bnds         ! variable ID for bounds variable
INTEGER :: i                  ! Loop counter
INTEGER :: timeunit_len       ! Length of timeunit string
INTEGER :: iunit              ! \ indices for timeunit, 1=seconds, 2=minutes,
INTEGER :: iunit2             ! /                       3=hours, 4=days
INTEGER :: filetime1          ! First timestep in current output file
INTEGER :: filetime2          ! Last timestep in current output file
INTEGER :: n1                 ! First output time index in current output file
INTEGER :: n2                 ! Last output time index in current output file
INTEGER :: ntimes             ! Number of output times in current output file
INTEGER :: unitin             ! Time unit specified in STASH
INTEGER :: ErrorStatus        ! Error code

LOGICAL :: lexists            ! TRUE if time dimension already exists
LOGICAL :: lstep0             ! TRUE if first output time is at timestep 0
LOGICAL :: loutput_ts0        ! TRUE if output of diagnostics at timestep 0
                              ! is enabled

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NCFILE_WRITE_TIME_DIM'
CHARACTER(LEN=errormessagelength) :: message ! Error message string
CHARACTER(LEN=profilename_length) :: timname ! Time profile name
CHARACTER(LEN=nf90_max_name)      :: dimname ! Dimension name
CHARACTER(LEN=nc_max_attlen)      :: stdname ! Standard name for dimension
CHARACTER(LEN=nc_max_attlen)      :: units   ! Units name for dimension

! Dr-Hook
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

dimid = 0

! Get time name for this STASH diagnostic
idiag = stlist(st_diag_address,ix)
timname = timpro(itim_b(idiag))

! Make sure timname has no characters disallowed by the CF standard
CALL check_cf_name(timname)

! Get STASH time information from stlist
input_code = stlist(st_input_code,ix)
IF (input_code < 0) THEN
  proccode = stlist(st_proc_no_code,-input_code)
ELSE
  proccode = stlist(st_proc_no_code,ix)
END IF
starttime = stlist(st_start_time_code,ix)
endtime = stlist(st_end_time_code,ix)
output_freq = stlist(st_freq_code,ix)
loutput_ts0 = stlist(st_output_ts0_code,ix) /= 0

! Current model time-step
step = stepim(atmos_im)

IF (proccode == st_accum_code .OR. proccode == st_time_mean_code .OR. &
    proccode == st_max_code   .OR. proccode == st_min_code) THEN
  sampling_freq = stlist(st_freq_code,-input_code)
ELSE
  sampling_freq = -1
END IF

IF (printstatus >= prdiag) THEN
  CALL umPrint(TRIM(RoutineName)//' variables:',src=RoutineName)
  WRITE(umMessage,'(A,2(I3),2(I6))') &
            'time_avail,proccode,step,sampling_freq = ', &
             time_avail,proccode,step,sampling_freq
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,3(I6))') &
            'starttime,endtime,output_freq = ', &
             starttime,endtime,output_freq
  CALL umPrint(umMessage,src=RoutineName)
END IF

! river_step is in seconds - calculate in timesteps
river_step_ts = river_step*(REAL(steps_per_periodim(atmos_im)) &
                /REAL(secs_per_periodim(atmos_im)))

! triffid_period and phenol_period are in days - calculate in timesteps
triffid_period_ts = triffid_period*(REAL(steps_per_periodim(atmos_im)) &
                    /REAL(secs_per_periodim(atmos_im)))*isec_per_day
phenol_period_ts = phenol_period*(REAL(steps_per_periodim(atmos_im)) &
                   /REAL(secs_per_periodim(atmos_im)))*isec_per_day

! For certain diagnostics the start time can be altered,
! these will need a new NetCDF time dimension
IF (time_avail == 2 .AND. a_lw_radstep_diag /= 1) THEN
  IF (a_lw_radstep_diag == a_sw_radstep_diag .AND. &
      a_lw_radstep_diag == a_lw_radstep_prog .AND. &
      a_lw_radstep_diag == a_sw_radstep_prog) THEN
    dimname = TRIM(timname)//'_rad'
  ELSE IF (a_lw_radstep_diag == a_sw_radstep_diag) THEN
    dimname = TRIM(timname)//'_rad_diag'
  ELSE IF (a_lw_radstep_diag == a_lw_radstep_prog) THEN
    dimname = TRIM(timname)//'_lwrad'
  ELSE
    dimname = TRIM(timname)//'_lwrad_diag'
  END IF
ELSE IF (time_avail == 3 .AND. a_sw_radstep_diag /= 1) THEN
  IF (a_sw_radstep_diag == a_lw_radstep_diag .AND. &
      a_sw_radstep_diag == a_sw_radstep_prog .AND. &
      a_sw_radstep_diag == a_lw_radstep_prog) THEN
    dimname = TRIM(timname)//'_rad'
  ELSE IF (a_sw_radstep_diag == a_lw_radstep_diag) THEN
    dimname = TRIM(timname)//'_rad_diag'
  ELSE IF (a_sw_radstep_diag == a_sw_radstep_prog) THEN
    dimname = TRIM(timname)//'_swrad'
  ELSE
    dimname = TRIM(timname)//'_swrad_diag'
  END IF
ELSE IF (time_avail == 13 .AND. a_conv_step /= 1) THEN
  dimname = TRIM(timname)//'_conv'
ELSE IF (time_avail == 14 .AND. phenol_period_ts /= 1) THEN
  dimname = TRIM(timname)//'_phenol'
ELSE IF (time_avail == 15 .AND. triffid_period_ts /= 1) THEN
  dimname = TRIM(timname)//'_triffid'
ELSE IF (time_avail == 16 .AND. river_step_ts /= 1) THEN
  dimname = TRIM(timname)//'_river'
ELSE IF (time_avail == 17 .AND. interval /= 1) THEN
  dimname = TRIM(timname)//'_ukca'
ELSE IF (time_avail == 18 .AND. a_lw_radstep_prog /= 1) THEN
  IF (a_lw_radstep_prog == a_sw_radstep_prog .AND. &
      a_lw_radstep_prog == a_lw_radstep_diag .AND. &
      a_lw_radstep_prog == a_sw_radstep_diag) THEN
    dimname = TRIM(timname)//'_rad'
  ELSE IF (a_lw_radstep_prog == a_sw_radstep_prog) THEN
    dimname = TRIM(timname)//'_rad_prog'
  ELSE IF (a_lw_radstep_prog == a_lw_radstep_diag) THEN
    dimname = TRIM(timname)//'_lwrad'
  ELSE
    dimname = TRIM(timname)//'_lwrad_prog'
  END IF
ELSE IF (time_avail == 19 .AND. a_sw_radstep_prog /= 1) THEN
  IF (a_sw_radstep_prog == a_lw_radstep_prog .AND. &
      a_sw_radstep_prog == a_sw_radstep_diag .AND. &
      a_sw_radstep_prog == a_lw_radstep_diag) THEN
    dimname = TRIM(timname)//'_rad'
  ELSE IF (a_sw_radstep_prog == a_lw_radstep_prog) THEN
    dimname = TRIM(timname)//'_rad_prog'
  ELSE IF (a_sw_radstep_prog == a_sw_radstep_diag) THEN
    dimname = TRIM(timname)//'_swrad'
  ELSE
    dimname = TRIM(timname)//'_swrad_prog'
  END IF
ELSE
  dimname = timname
END IF

stdname = "time"

! If list of times specified extract them and set lstep0
IF (output_freq < 0) THEN
  DO i=1,nsttims
    times_list(i) = sttabl(i,-output_freq)
    IF (times_list(i) == st_end_of_list) THEN
      ntimes_list = i-1
      EXIT
    END IF
  END DO
  lstep0 = times_list(1) == 0
ELSE
  lstep0 = starttime == 0
END IF

! Only sections 0,15,16,30,33,34 have diagnostics output at timestep 0,
! so there could possibly be 2 coordinates, with and without timestep 0
IF (lstep0 .AND. loutput_ts0 .AND. &
   (section ==  0 .OR. section == 15 .OR. section == 16 .OR. &
    section == 30 .OR. section == 33 .OR. section == 34)) THEN
  dimname = TRIM(dimname)//'_0'
END IF

! Get time unit information for output times
unitin = stlist(st_time_unit3_code,ix)
SELECT CASE (unitin)
CASE (stsh_seconds)
  iunit = 1
CASE (stsh_minutes)
  iunit = 2
CASE (stsh_hours)
  iunit = 3
CASE (stsh_days)
  iunit = 4
CASE (stsh_dumps)
  CALL iunit_from_ts(dumpfreqim(atmos_im),atmos_im,iunit)
CASE (stsh_timesteps)
  IF (output_freq > 0) THEN
    CALL iunit_from_ts(output_freq,atmos_im,iunit)
  ELSE IF (output_freq < 0) THEN
    ! Loop over output times to work out longest possible time unit
    CALL iunit_from_ts(times_list(1),atmos_im,iunit)
    IF (iunit > 1) THEN
      DO i=2,ntimes_list
        CALL iunit_from_ts(times_list(i),atmos_im,iunit2)
        IF (iunit2 < iunit) THEN
          iunit = iunit2
        END IF
        IF (iunit == 1) EXIT ! Seconds is shortest possible unit
      END DO
    END IF
  ELSE
    errorstatus = 101
    CALL ereport(RoutineName,ErrorStatus,'Error time frequency is 0')
  END IF
CASE DEFAULT
  errorstatus = 102
  WRITE(message,'(A,I3)') &
    'Unknown unit code for output times, unitin = ',unitin
  CALL ereport(RoutineName,ErrorStatus,message)
END SELECT
timeunit_len = LEN_TRIM(timeunit(iunit))

! Construct units attribute
units = timeunit(iunit)(1:timeunit_len)
WRITE(units(timeunit_len+1:timeunit_len+26), &
      '(A7,I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') &
      ' since ', &
      model_basis_time(1),'-',model_basis_time(2),'-', &
      model_basis_time(3),' ',model_basis_time(4),':', &
      model_basis_time(5),':',model_basis_time(6)

! Work out number of time points for this NetCDF dimension

! Get range of timesteps valid for current output file
CALL reinit_file_times(um_file,filetime1,filetime2)

IF (output_freq < 0) THEN
  ! Count number of output times valid for current output file
  ntimes = 0
  DO i=1,ntimes_list
    IF (times_list(i) >= filetime1 .AND. &
        times_list(i) <= filetime2) THEN
      ntimes = ntimes+1
    END IF
  END DO
ELSE IF ((endtime > 0 .AND. endtime < filetime1) .OR. &
         starttime > filetime2) THEN
  ! Output times outside range of current output file
  ntimes = 0
ELSE
  ! Calculate first and last output time indices for current output file
  ! Then use these values to calculate the number of output times
  IF (starttime < filetime1) THEN
    n1 = CEILING(REAL(filetime1-starttime)/REAL(output_freq))+1
  ELSE
    n1 = 1
  END IF
  IF (endtime > 0 .AND. endtime < filetime2) THEN
    n2 = FLOOR(REAL(endtime-starttime)/REAL(output_freq))+1
  ELSE
    n2 = FLOOR(REAL(filetime2-starttime)/REAL(output_freq))+1
  END IF
  IF (printstatus >= prdiag) THEN
    WRITE(umMessage,'(A,I4,A,I4)')'n1 = ',n1,' n2 = ',n2
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  ntimes = n2-n1+1
END IF

! There is one less time value if diagnostic has output at timestep 0 in
! curent output file, but either this has been disabled or is in a STASH section
! which doesn't allow output at timestep 0, i.e. not one of STASH sections 
! 0,15,16,30,33,34
IF (ntimes > 0 .AND. step == 0 .AND. lstep0 .AND. &
    (.NOT. loutput_ts0 .OR. &
     (section /=  0 .AND. section /= 15 .AND. section /= 16 .AND. &
      section /= 30 .AND. section /= 33 .AND. section /= 34))) THEN
  ntimes = ntimes-1
END IF

IF (printstatus >= prdiag) THEN
  WRITE(umMessage,'(A,3(I6))') &
            'filetime1,filetime2,ntimes = ', &
             filetime1,filetime2,ntimes
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(6(A))') &
            'dimname = ',TRIM(dimname), &
            ' : stdname = ',TRIM(stdname), &
            ' : units = ',TRIM(units)
  CALL umPrint(umMessage,src=RoutineName)
END IF

IF (ntimes == 0) THEN
  ! Dimension range outside bounds of current output file so finish
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Get cell_methods attribute information for the time dimension
SELECT CASE (proccode)
CASE (st_replace_code)
  IF (cell_methods == ' ') THEN
    cell_methods = TRIM(dimname)//': point'
  ELSE
    cell_methods = TRIM(cell_methods)//' '//TRIM(dimname)//': point'
  END IF
CASE (st_accum_code)
  IF (cell_methods == ' ') THEN
    cell_methods = TRIM(dimname)//': sum'
  ELSE
    cell_methods = TRIM(cell_methods)//' '//TRIM(dimname)//': sum'
  END IF
CASE (st_time_mean_code)
  IF (cell_methods == ' ') THEN
    cell_methods = TRIM(dimname)//': mean'
  ELSE
    cell_methods = TRIM(cell_methods)//' '//TRIM(dimname)//': mean'
  END IF
CASE (st_max_code)
  IF (cell_methods == ' ') THEN
    cell_methods = TRIM(dimname)//': maximum'
  ELSE
    cell_methods = TRIM(cell_methods)//' '//TRIM(dimname)//': maximum'
  END IF
CASE (st_min_code)
  IF (cell_methods == ' ') THEN
    cell_methods = TRIM(dimname)//': minimum'
  ELSE
    cell_methods = TRIM(cell_methods)//' '//TRIM(dimname)//': minimum'
  END IF
END SELECT

! Get time unit information for sampling frequency
IF (proccode == st_accum_code .OR. proccode == st_time_mean_code .OR. &
    proccode == st_max_code   .OR. proccode == st_min_code) THEN

  unitin = stlist(st_time_unit2_code,ix)
  SELECT CASE (unitin)
  CASE (stsh_seconds)
    iunit = 1
  CASE (stsh_minutes)
    iunit = 2
  CASE (stsh_hours)
    iunit = 3
  CASE (stsh_days)
    iunit = 4
  CASE (stsh_timesteps)
    CALL iunit_from_ts(sampling_freq,atmos_im,iunit)
  CASE (stsh_dumps)
    CALL iunit_from_ts(dumpfreqim(atmos_im),atmos_im,iunit)
  CASE DEFAULT
    errorstatus = 103
    WRITE(message,'(A,I3)') &
      'Unknown unit code for sampling period, unitin = ',unitin
    CALL ereport(RoutineName,ErrorStatus,message)
  END SELECT
  CALL add_interval(cell_methods,sampling_freq,atmos_im,iunit)
END IF
IF (printstatus >= prdiag) THEN
  WRITE(umMessage,'(A,A)')'cell_methods = ',TRIM(cell_methods)
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Create or read NetCDF dimensions
CALL nc_create_dim(um_file,dimname,ntimes,dimid,lexists)

IF (lexists) THEN
  ! Time dimension already exists so finish
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

dimid_array(1) = dimid

! Create NetCDF time dimension variable
CALL nc_create_var(um_file,dimname,nc_dim_type_real,varid,dimid_array)

! Save stash index value for later use
um_file % nc_meta % stash_index_dim(dimid) = ix

! Create or read bounds dimension if needed
IF (proccode == st_accum_code .OR. proccode == st_time_mean_code .OR. &
    proccode == st_max_code   .OR. proccode == st_min_code) THEN
  CALL nc_create_dim(um_file,'bounds2',2,dimid_bnds)
  dimid_time_bnds(1) = dimid_bnds
  dimid_time_bnds(2) = dimid
END IF

! Create NetCDF dimension variable attributes
CALL nc_put_att(um_file,varid,'standard_name',stdname)
CALL nc_put_att(um_file,varid,'units',units)
CALL nc_put_att(um_file,varid,'axis','T')
IF (lcal360) THEN
  CALL nc_put_att(um_file,varid,'calendar','360_day')
ELSE
  CALL nc_put_att(um_file,varid,'calendar','gregorian')
END IF
IF (proccode == st_accum_code .OR. proccode == st_time_mean_code .OR. &
    proccode == st_max_code   .OR. proccode == st_min_code) THEN
  CALL nc_put_att(um_file,varid,'bounds',TRIM(dimname)//'_bounds')
  CALL nc_create_var(um_file,TRIM(dimname)//'_bounds',nc_dim_type_real, &
                     varid_bnds,dimid_time_bnds)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_write_time_dim

! Work out suitable unit for given timestep. Chose longest time unit which can
! represent timestep as an integer value. Return value is an index to timeunit

SUBROUTINE iunit_from_ts(ntimestep,model,iunit)

IMPLICIT NONE

INTEGER,          INTENT(IN)   :: ntimestep   ! number of timesteps
INTEGER,          INTENT(IN)   :: model       ! internal model
INTEGER,          INTENT(OUT)  :: iunit       ! integer code for unit in
                                              ! increasing order of unit length

!
!  Local variables
!
INTEGER :: ts_sec   ! timestep value in seconds

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'IUNIT_FROM_TS'

! Dr-Hook
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Convert timestep value to seconds
ts_sec = ntimestep*secs_per_periodim(model)/steps_per_periodim(model)

IF (MOD(ts_sec,isec_per_min) /= 0) THEN
  ! Timestep not a whole number of minutes
  iunit = 1 ! seconds
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (MOD(ts_sec,isec_per_hour) /= 0) THEN
  ! Timestep not a whole number of hours
  iunit = 2 ! minutes
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (MOD(ts_sec,isec_per_day) /= 0) THEN
  ! Timestep not a whole number of days
  iunit = 3 ! hours
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Timestep is a whole number of days
iunit = 4 ! days

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE iunit_from_ts

! Add time interval to cell_methods attribute

SUBROUTINE add_interval(cell_methods,ntimestep,model,iunit)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(INOUT) :: cell_methods! cell_methods attribute string
INTEGER,          INTENT(IN)    :: ntimestep   ! number of timesteps
INTEGER,          INTENT(IN)    :: model       ! internal model
INTEGER,          INTENT(IN)    :: iunit       ! integer code for unit in
                                               ! increasing order of unit length

!
!  Local variables
!
INTEGER      :: interval      ! timestep value in seconds/minutes/hours/days
                              ! depending on value of iunit
INTEGER      :: timeunit_len  ! actual length of timeunit

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ADD_INTERVAL'

! Dr-Hook
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Convert timestep value to seconds
interval = ntimestep*secs_per_periodim(model)/steps_per_periodim(model)

IF (iunit == 2) THEN
  interval = interval/isec_per_min  ! convert to minutes
ELSE IF (iunit == 3) THEN
  interval = interval/isec_per_hour ! convert to hours
ELSE IF (iunit == 4) THEN
  interval = interval/isec_per_day  ! convert to days
END IF

timeunit_len = LEN_TRIM(timeunit(iunit))
IF (interval == 1) THEN
  ! Don't use plural form of timeunit
  timeunit_len = timeunit_len-1
END IF

WRITE(cell_methods,'(SS,A,I0,A)')TRIM(cell_methods)//' (interval: ',interval, &
                                 ' '//timeunit(iunit)(1:timeunit_len)//')'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE add_interval

END MODULE ncfile_write_time_dim_mod


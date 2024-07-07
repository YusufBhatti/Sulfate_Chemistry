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
!  Purpose: Writes data for coordinate variables and bounds variables
!           created in ncfile_write_time_dim


MODULE ncfile_write_time_var_mod

USE netcdf,               ONLY: nf90_max_name ! External netCDF library module
USE yomhook,              ONLY: lhook,dr_hook
USE parkind1,             ONLY: jprb,jpim
USE umprintmgr,           ONLY: umprint,ummessage,printstatus,prdiag
USE ereport_mod,          ONLY: ereport
USE file_manager,         ONLY: um_file_type
USE submodel_mod,         ONLY: atmos_im
USE model_time_mod,       ONLY: stepim
USE nlstgen_mod,          ONLY: steps_per_periodim,secs_per_periodim
USE conversions_mod,      ONLY: isec_per_min,isec_per_hour,rsec_per_day
USE stash_array_mod,      ONLY: stlist,sttabl,nsttims
USE reinit_file_times_mod,ONLY: reinit_file_times
USE stparam_mod,          ONLY: &
    st_freq_code,st_start_time_code,st_end_time_code,st_proc_no_code, &
    st_input_code,st_sect_no_code,st_period_code,st_end_of_list, &
    st_accum_code,st_time_mean_code,st_max_code,st_min_code, &
    st_output_ts0_code,st_infinite_time
USE umnetcdf_mod,         ONLY: &
    nc_inq_ndims,nc_inq_dim,nc_get_varid,nc_get_att,nc_put_var,nc_max_attlen

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NCFILE_WRITE_TIME_VAR_MOD'

CONTAINS

SUBROUTINE ncfile_write_time_var(um_file)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file

!
!  Local variables
! 
INTEGER :: ndims           ! Number of dimensions in NetCDF file
INTEGER :: dimlen          ! Dimension size
INTEGER :: varid           ! Variable id of dimension
INTEGER :: ix              ! STASH item index
INTEGER :: i               ! NetCDF dimension index
INTEGER :: itime           ! Time index
INTEGER :: jtime           ! Another time index
INTEGER :: itime1          ! Valid first time index value
INTEGER :: starttime       ! Start time of diagnostic output
INTEGER :: endtime         ! End time of diagnostic output
INTEGER :: output_freq     ! Frequency of diagnostic output
INTEGER :: period          ! Number of time steps in time processing period
INTEGER :: input_code      ! When -ve gives index of parent record in stlist
INTEGER :: proccode        ! Time processing code for data
INTEGER :: section         ! STASH section number
INTEGER :: filetime1       ! First timestep in current output file
INTEGER :: filetime2       ! Last timestep in current output file
INTEGER :: timeval         ! Time value from output times list
INTEGER :: n1              ! First output time index in current output file
INTEGER :: n2              ! Last output time index in current output file
INTEGER :: start(1)        ! Start position of dimension output data
INTEGER :: count1(1)       ! Number of items of dimension output data
INTEGER :: start_bnds(2)   ! Start position of dimension bounds output data
INTEGER :: count_bnds(2)   ! Number of items of dimension bounds output data
INTEGER :: step            ! Timestep of model
INTEGER :: ErrorStatus     ! Error code

REAL              :: fac1               ! Coordinate time units per period
REAL              :: fac2               ! Reciprocal of time steps per period
REAL, ALLOCATABLE :: dimdata(:)         ! Variable to contain dimension values
REAL, ALLOCATABLE :: dimdata_bnds(:,:)  ! Variable to contain dimension 
                                        ! boundary values

LOGICAL :: lexists                      ! TRUE if netCDF attribute exists
LOGICAL :: lstep0                       ! TRUE if first timestep value is 0
LOGICAL :: loutput_ts0                  ! TRUE if output of diagnostics
                                        ! at timestep 0 is enabled

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NCFILE_WRITE_TIME_VAR'
CHARACTER(LEN=nf90_max_name) :: dimname ! NetCDF dimension name
CHARACTER(LEN=1)             :: axis    ! Dimension axis type
CHARACTER(LEN=nc_max_attlen) :: units   ! Units name for dimension


! Dr-Hook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Current model time-step
step = stepim(atmos_im)

! Get numnber of dimensions in NetCDF file
CALL nc_inq_ndims(um_file,ndims)

! Loop over each dimension, i is netCDF dimid
DO i=1,ndims

  ! Get dimension name and length for this dimid
  CALL nc_inq_dim(um_file,i,dimname,dimlen)

  ! Not interested in bounds dimension
  IF (dimname(1:6) == 'bounds' .OR. dimlen == 0) CYCLE

  ! Get id of variable with same name as dimension
  CALL nc_get_varid(um_file,dimname,varid)

  ! Get axis attribute
  CALL nc_get_att(um_file,varid,'axis',axis,lexists)
  IF (.NOT. lexists) CYCLE

  ! This routine is for processing time dimensions, so only select those
  IF (axis == 'T') THEN

    ! Get stash index value
    ix = um_file % nc_meta % stash_index_dim(i)

    IF (printstatus >= prdiag) THEN
      WRITE(umMessage,'(A,A,A,I6,A,I4,A)') &
            'Time dimension name ',TRIM(dimname), &
            ' of size ',dimlen,' on unit ',um_file % UNIT,' found'
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,A,A,I5)') &
                'STASH index for time dimension ',TRIM(dimname),' is ',ix
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    ! Get dimension units to calculate conversion factors from timesteps
    CALL nc_get_att(um_file,varid,'units',units,lexists)

    ! Calculate time unit conversion factors
    IF (units(1:6) == 'second') THEN
      fac1 = REAL(secs_per_periodim(atmos_im))
    ELSE IF (units(1:6) == 'minute') THEN
      fac1 = REAL(secs_per_periodim(atmos_im))/REAL(isec_per_min)
    ELSE IF (units(1:4) == 'hour') THEN
      fac1 = REAL(secs_per_periodim(atmos_im))/REAL(isec_per_hour)
    ELSE IF (units(1:3) == 'day') THEN
      fac1 = REAL(secs_per_periodim(atmos_im))/rsec_per_day
    ELSE
      ErrorStatus = -102
      CALL ereport(RoutineName,ErrorStatus, &
                   'Unrecognised time dimension unit '//TRIM(units))
      CYCLE
    END IF
    fac2 = 1.0/REAL(steps_per_periodim(atmos_im))

    ! Intialise start and count variables needed to write netCDF data
    start(1) = 1
    count1(1) = dimlen
    start_bnds(1) = 1
    start_bnds(2) = 1
    count_bnds(1) = 2
    count_bnds(2) = dimlen

    ALLOCATE(dimdata(dimlen))

    ! Get STASH and time information and meaning options
    section = stlist(st_sect_no_code,ix)
    starttime = stlist(st_start_time_code,ix)
    endtime = stlist(st_end_time_code,ix)
    output_freq = stlist(st_freq_code,ix)
    loutput_ts0 = stlist(st_output_ts0_code,ix) /= 0
    IF (output_freq < 0) THEN
      lstep0 = sttabl(1,-output_freq) == 0
    ELSE
      lstep0 = starttime == 0
    END IF

    ! First time index value depends on whether timestep 0 has been selected
    ! and is available for selection
    IF (lstep0 .AND. step == 0) THEN
      IF (section ==  0 .OR. section == 15 .OR. section == 16 .OR. &
          section == 30 .OR. section == 33 .OR. section == 34) THEN
        IF (loutput_ts0) THEN
          itime1 = 1
        ELSE
          itime1 = 2
        END IF
      ELSE
        itime1 = 2
      END IF
    ELSE
      itime1 = 1
    END IF

    ! Get range of timesteps valid for current output file
    CALL reinit_file_times(um_file,filetime1,filetime2)

    IF (output_freq < 0) THEN
      ! Extract time values valid for current output file
      jtime = 1
      DO itime=itime1,nsttims
        timeval = sttabl(itime,-output_freq)
        IF (timeval >= filetime1 .AND. timeval <= filetime2) THEN
          dimdata(jtime) = fac1*timeval*fac2
          jtime = jtime+1
        END IF
        IF (timeval == st_end_of_list) EXIT
      END DO
    ELSE
      ! Calculate first and last output time indices for current output file
      ! Then use these values to calculate the time values
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
      DO itime=n1+itime1-1,n2
        dimdata(itime-n1-itime1+2) = fac1*(starttime+(itime-1)*output_freq)*fac2
      END DO
    END IF

    ! Write time dimension data
    CALL nc_put_var(um_file,varid,dimdata(1:dimlen),start,count1)

    ! Get processing code
    input_code = stlist(st_input_code,ix)
    IF (input_code < 0) THEN
      proccode = stlist(st_proc_no_code,-input_code)
    ELSE
      proccode = stlist(st_proc_no_code,ix)
    END IF

    IF (proccode == st_accum_code .OR. proccode == st_time_mean_code .OR. &
        proccode == st_max_code   .OR. proccode == st_min_code) THEN

      ! Get dimension bounds variable id
      CALL nc_get_varid(um_file,TRIM(dimname)//'_bounds',varid)
 
      ALLOCATE(dimdata_bnds(2,dimlen))

      ! Calculate time dimension boundary values
      period = stlist(st_period_code,ix)
      DO itime=1,dimlen
        IF (period == st_infinite_time) THEN
          dimdata_bnds(1,itime) = 0.0
        ELSE
          dimdata_bnds(1,itime) = dimdata(itime) - fac1*period*fac2
        END IF
        dimdata_bnds(2,itime) = dimdata(itime)
      END DO

      ! Write dimension bounds information to NetCDF file
      CALL nc_put_var(um_file,varid,dimdata_bnds,start_bnds,count_bnds)

      DEALLOCATE(dimdata_bnds)
   
    END IF

    DEALLOCATE(dimdata)

  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_write_time_var

END MODULE ncfile_write_time_var_mod

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
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: NetCDF output
!
!  Purpose: Initialises variables needed to write data on a CRUN

MODULE init_nc_crun_mod

! External netCDF library module
USE netcdf,                ONLY: nf90_max_name
USE yomhook,               ONLY: lhook,dr_hook
USE parkind1,              ONLY: jprb,jpim
USE ereport_mod,           ONLY: ereport
USE file_manager,          ONLY: um_file_type
USE umprintmgr,            ONLY: umprint,ummessage
USE model_time_mod,        ONLY: stepim
USE nlstgen_mod,           ONLY: steps_per_periodim,secs_per_periodim
USE conversions_mod,       ONLY: isec_per_min,isec_per_hour,rsec_per_day
USE submodel_mod,          ONLY: atmos_im
USE stash_array_mod,       ONLY: stlist,stindex,nitems,nsects
USE stextend_mod,          ONLY: in_s
USE stparam_mod,           ONLY: st_output_code,st_output_type,st_netcdf
USE umnetcdf_mod,          ONLY: &
    nc_file_open,nc_mode_open,nc_time_index,nc_dim_used, &
    nc_get_varid,nc_inq_var_dimid,nc_inq_dimname,nc_inq_dimlen, &
    nc_get_att,nc_get_var,set_varid_val,nc_max_attlen,nc_max_var_dims

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_NC_CRUN_MOD'

PRIVATE ncfile_get_last_time_index

CONTAINS

SUBROUTINE init_nc_crun(um_file)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file

!
!  Local variables
!
INTEGER :: dimids(nc_max_var_dims) ! NetCDF dimension IDs array
INTEGER :: stash_count ! Count of the number of fields of a particular
                       ! STASH section/item for a given netCDF file
INTEGER :: im_index    ! internal model index
INTEGER :: idim1       ! Index over netCDF dimensions
INTEGER :: dimid       ! NetCDF dimension ID
INTEGER :: ndims       ! Number of dimensions for STASH item
INTEGER :: time_index  ! Time index of current timestep in netCDF file
INTEGER :: varid       ! NetCDF variable id for STASH item
INTEGER :: dimvarid    ! NetCDF dimension variable id
INTEGER :: i           ! loop index over STASH item codes
INTEGER :: j           ! loop index over STASH section codes
INTEGER :: ix          ! STASH item index

LOGICAL :: lexists     ! TRUE if netCDF attribute already exists

CHARACTER(LEN=*), PARAMETER  :: RoutineName = "INIT_NC_CRUN"
CHARACTER(LEN=nf90_max_name) :: varname ! STASH variable name
CHARACTER(LEN=nf90_max_name) :: dimname ! STASH dimension name
CHARACTER(LEN=1)             :: axis    ! Dimension axis type

! Dr-Hook
REAL(KIND=jprb)              :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

WRITE(umMessage,'(A,I3)')TRIM(RoutineName)//': Opening netCDF file '// &
          TRIM(um_file % filename)//' on unit ', um_file % UNIT
CALL umPrint(umMessage,src=RoutineName)

CALL nc_file_open(um_file,nc_mode_open)

! Extract information needed to write STASH from netCDF file
im_index = 1

DO j=0,nsects
  DO i=1,nitems
    IF (in_s(1,atmos_im,j,i) >= 1) THEN

      ! Reset value of stash_count
      stash_count = 0

      DO ix=stindex(1,i,j,im_index), &
            stindex(1,i,j,im_index)+stindex(2,i,j,im_index)-1

        IF (-stlist(st_output_code,ix) == um_file % UNIT .AND. &
             stlist(st_output_type,ix) == st_netcdf) THEN

          stash_count = stash_count + 1
          nc_time_index(ix) = 0
          nc_dim_used(:,ix) = .FALSE.

          ! Create netCDF variable name from the model ID, 
          ! STASH section, STASH item and STASH count
          IF (stash_count == 1) THEN
            WRITE(varname,'(SS,A,I2.2,A,I2.2,A,I3.3)') &
                  'STASH_m',atmos_im,'s',j,'i',i
          ELSE
            WRITE(varname,'(SS,A,I2.2,A,I2.2,A,I3.3,A,I0)') &
                  'STASH_m',atmos_im,'s',j,'i',i,"_",stash_count
          END IF

          ! Get varid of variable with this varname
          CALL nc_get_varid(um_file,varname,varid)
          IF (varid == -1) THEN
            ! Variable doesn't exist in file for this period so skip
            WRITE(umMessage,'(A,A,A,A,A)') &
                      'Variable ',TRIM(varname), &
                      ' in file ',TRIM(um_file % filename), &
                      ' does not exist in this period, skipping'
            CALL umPrint(umMessage,src=RoutineName)
            CYCLE
          END IF

          ! Save varid value for later use
          CALL set_varid_val(ix,varid)

          ! Get dimension ids for this variable
          CALL nc_inq_var_dimid(um_file,varid,dimids,ndims)

          DO idim1=1,ndims
            dimid = dimids(idim1)
            ! Get dimension name
            CALL nc_inq_dimname(um_file,dimid,dimname)
              
            ! Get dimension variable id
            CALL nc_get_varid(um_file,dimname,dimvarid)
              
            ! Get dimension axis attribute
            CALL nc_get_att(um_file,dimvarid,'axis',axis,lexists)

            ! set nc_dim_used to .TRUE. if dimension used by this variable
            SELECT CASE (axis)
            CASE ('X')
              nc_dim_used(1,ix) = .TRUE.
            CASE ('Y')
              nc_dim_used(2,ix) = .TRUE.
            CASE ('Z')
              nc_dim_used(3,ix) = .TRUE.
            CASE ('T')
              nc_dim_used(5,ix) = .TRUE.
              ! calculate time index for this timestep
              CALL ncfile_get_last_time_index(um_file,dimvarid,dimid,time_index)
              nc_time_index(ix) = time_index
            CASE DEFAULT
              ! Pseudo dimension
              nc_dim_used(4,ix) = .TRUE.
            END SELECT
          END DO

        END IF
      END DO
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_nc_crun

! Work out time index, for this variable, of the current time step
SUBROUTINE ncfile_get_last_time_index(um_file,varid,dimid,time_index)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)  :: varid      ! Netcdf time dimension variable id
INTEGER, INTENT(IN)  :: dimid      ! Time dimension netCDF id
INTEGER, INTENT(OUT) :: time_index ! Time index of current timestep
                                   ! in netCDF file

!
!  Local variables
!
INTEGER           :: step             ! Model step number 
INTEGER           :: steps_per_period ! Number of time steps in defining period 
INTEGER           :: secs_per_period  ! Number of seconds in defining period 
INTEGER           :: dimlen           ! Length of time dimension
INTEGER           :: i                ! Index over time values
INTEGER           :: start(1)         ! Start position of time dimension
INTEGER           :: count1(1)        ! Number of items in time dimension
INTEGER           :: dimdata_step     ! Time dimension values
                                      ! converted into number of timesteps
INTEGER           :: ErrorStatus      ! Error return code

LOGICAL           :: lexists          ! TRUE if netCDF attribute already exists

REAL              :: fac1             ! Real value of steps_per_period
REAL              :: fac2             ! Number of defining periods in time unit
REAL, ALLOCATABLE :: dimdata(:)       ! Variable to contain dimension values

CHARACTER(LEN=*), PARAMETER  :: RoutineName = "NCFILE_GET_LAST_TIME_INDEX"
CHARACTER(LEN=nc_max_attlen) :: units ! Units name for time dimension

! Dr-Hook
REAL(KIND=jprb)             :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set timestep information
step = stepim(atmos_im)
secs_per_period = secs_per_periodim(atmos_im)
steps_per_period = steps_per_periodim(atmos_im)

! Get time dimension units to calculate conversion factors to timesteps
CALL nc_get_att(um_file,varid,'units',units,lexists)

fac1 = REAL(steps_per_period)
IF (units(1:6) == 'second') THEN
  fac2 = 1.0/REAL(secs_per_period)
ELSE IF (units(1:6) == 'minute') THEN
  fac2 = REAL(isec_per_min)/REAL(secs_per_period)
ELSE IF (units(1:4) == 'hour') THEN
  fac2 = REAL(isec_per_hour)/REAL(secs_per_period)
ELSE IF (units(1:3) == 'day') THEN
  fac2 = rsec_per_day/REAL(secs_per_period)
ELSE
  ErrorStatus = 1
  CALL ereport(RoutineName,ErrorStatus, &
       'Error unrecognised time dimension unit '//TRIM(units))
END IF

! Get number of time dimension values
CALL nc_inq_dimlen(um_file,dimid,dimlen) 

ALLOCATE(dimdata(dimlen))

start(1)=1
count1(1)=dimlen

! Get time dimension values
CALL nc_get_var(um_file,varid,dimdata,start,count1)

! Calculate index of last time value
IF (step < NINT(fac1*dimdata(1)*fac2)) THEN
  time_index = 0
ELSE IF (step > NINT(fac1*dimdata(dimlen)*fac2)) THEN
  time_index = dimlen
ELSE
  DO i=1,dimlen
    dimdata_step = NINT(fac1*dimdata(i)*fac2)
    IF (dimdata_step == step) THEN
      time_index=i
      EXIT
    ELSE IF (dimdata_step > step) THEN
      time_index=i-1
      EXIT
    END IF
  END DO
END IF

DEALLOCATE(dimdata)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_get_last_time_index

END MODULE init_nc_crun_mod

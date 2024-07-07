! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE rcf_ideal_file_forcing_mod

  USE planet_constants_mod
  USE atm_fields_bounds_mod
  USE model_domain_mod
  USE horiz_grid_mod

  IMPLICIT NONE
  ! Description: 
  !  Module containing subroutines to read in a NETCDF
  !  format T-P profile
  !
  ! Method:
  !  
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

  ! P-T profile from file
  REAL, ALLOCATABLE :: file_temperature(:), file_pressure(:)

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
    ModuleName='RCF_IDEAL_FILE_FORCING_MOD'

CONTAINS

  SUBROUTINE read_file_temp(filename)
  
    USE netcdf
    USE ereport_mod, ONLY: ereport
    USE errormessagelength_mod, ONLY: errormessagelength
    USE um_types, ONLY: integer32

    IMPLICIT NONE
  
    ! Input
    CHARACTER(LEN=*) :: filename
    
    ! Local
    INTEGER(KIND=integer32) :: status, ncid
    CHARACTER(LEN=errormessagelength) :: cmessage
    INTEGER(KIND=integer32) :: varid_t, varid_p, dimid(1)
    INTEGER(KIND=integer32) :: file_num_levels
        
    status = nf90_open(TRIM(filename), nf90_nowrite, ncid)
    IF(status /= nf90_noerr) THEN
      cmessage = 'Error in loading file ' // filename
      CALL ereport('READ_FILE_TEMP',status,cmessage)
    END IF
    
    ! Get variable IDs
    status = nf90_inq_varid(ncid, 'temperature', varid_t)
    IF(status /= nf90_noerr) THEN
      cmessage = 'Error reading temperature'
      CALL ereport('READ_FILE_TEMP',status,cmessage)
    END IF
    status = nf90_inq_varid(ncid, 'pressure_si', varid_p)
    IF(status /= nf90_noerr) THEN
      ! Fail if pressure in SI units cannot be found.
      status   = 1
      cmessage = 'Error reading pressure_si'
      CALL ereport('READ_FILE_TEMP',status,cmessage)
    END IF
    
    ! Get ID of vertical dimension
    status = nf90_inquire_variable(ncid, varid_t,                       &
        dimids=dimid)
    IF(status /= nf90_noerr) THEN
      cmessage = 'Error getting ID of vertical dimension'
      CALL ereport('READ_FILE_TEMP',status,cmessage)
    END IF
    ! Get length of vertical dimension
    status = nf90_inquire_dimension(ncid, dimid(1),                     &
        len=file_num_levels)
    IF(status /= nf90_noerr) THEN
      cmessage = 'Error getting length of vertical dimension'
      CALL ereport('READ_FILE_TEMP',status,cmessage)
    END IF
    
    ! Allocate temperature and pressure arrays
    IF (ALLOCATED(file_temperature)) THEN
      DEALLOCATE(file_temperature)
    END IF
    IF (ALLOCATED(file_pressure)) THEN
      DEALLOCATE(file_pressure)
    END IF
    ALLOCATE(file_temperature(file_num_levels),                         &
        file_pressure(file_num_levels))
    
    ! Read temperature and pressure
    status = nf90_get_var(ncid, varid_t, file_temperature)
    IF(status /= nf90_noerr) THEN
      cmessage = 'Error reading temperature'
      CALL ereport('READ_FILE_TEMP',status,cmessage)
    END IF
    status = nf90_get_var(ncid, varid_p, file_pressure)
    IF(status /= nf90_noerr) THEN
      cmessage = 'Error reading pressure'
      CALL ereport('READ_FILE_TEMP',status,cmessage)
    END IF
    
    ! Close file
    status = nf90_close(ncid)
    IF(status /= nf90_noerr) THEN
      cmessage = 'Error closing file ' // filename
      CALL ereport('READ_FILE_TEMP',status,cmessage)
    END IF
      
  END SUBROUTINE read_file_temp

! ----------------------------------------------------------------------
! FUNCTIONS: Interpolate in provided P-T profile
! ----------------------------------------------------------------------

  REAL FUNCTION file_temp(exner_theta_levels)
    ! Function which returns the night side equilibrium temperature
    ! for HD209458b, smoothed or original version
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE
    
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='FILE_TEMP'
    
    ! Input
    REAL, INTENT (IN) :: exner_theta_levels

    ! Local variables
    REAL, ALLOCATABLE :: log_file_pressure(:)
    REAL :: pressure

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! First construct the pressure variable
    pressure=(exner_theta_levels**recip_kappa)*p_zero
   
    ! Avoid temporary memory copies 
    ALLOCATE(log_file_pressure(SIZE(file_pressure)))
    log_file_pressure(:) = LOG10(file_pressure(:))

    ! Interpolate in log file_pressure to find temperature
    file_temp = interp(log_file_pressure, file_temperature,          &
        LOG10(pressure))

    DEALLOCATE(log_file_pressure)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  
CONTAINS
  
  FUNCTION interp(x ,y, xi) RESULT(yi)
  ! Linear interpolation function
  
    IMPLICIT NONE

    REAL, INTENT(IN) ::                                                 &
        x(:)                                                            &
!         Array with x-values
      , y(:)                                                            &
!         Array with y-values
      , xi
!         Value at which the y-coordinate is wanted

    REAL ::                                                             &
        yi
!         Value at xi.

    INTEGER ::                                                          &
        x_len                                                           &
!         Length of x-array
      , i
!         Loop index


    x_len = SIZE(x)

    IF (xi < x(1)) THEN
      yi = y(1)
      RETURN
    ELSE IF (xi > x(x_len)) THEN
      yi = y(x_len)
      RETURN
    ELSE
      DO i=2,x_len
        IF (xi <= x(i)) THEN
          yi = (y(i) - y(i-1))/(x(i) - x(i-1))*(xi - x(i-1)) + y(i-1)
          RETURN
        END IF
      END DO
    END IF
    
  END FUNCTION interp
  
  END FUNCTION file_temp

END MODULE rcf_ideal_file_forcing_mod

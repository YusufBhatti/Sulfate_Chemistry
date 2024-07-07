! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the points to calculate SW radiation
!
! Purpose: This routine mainly deals with the spatial degradation and
!          the use of segments
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!-----------------------------------------------------------------------
!
MODULE prelim_swrad_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'PRELIM_SWRAD_MOD'
CONTAINS

SUBROUTINE prelim_swrad(error_code,                                     &
! Parallel variables
       at_extremity, n_proc,                                            &
! Model Dimensions
       row_length,rows,n_rad_layers,                                    &
! Model Switches
       l_rad_deg, l_subsample, l_geostationary, l_spherical_solar,      &
! Timestepping information
       timestep_number,a_sw_radstep,                                    &
! Satelitte Geometry
       min_view_lon, max_view_lon,                                      &
       min_view_lat, max_view_lat,                                      &
! Number of call
       j_sw,                                                            &
! Other variables
       true_latitude,true_longitude, seconds_since_midnight,            &
       tot_daylight_points,daylight_points,day_fraction, day_frac_sph,  &
       list_daylight_points, rad_mask, switch,                          &
       first_data_interp, diag_row_list, diag_col_list)


USE um_parparams, ONLY: pnorth, psouth

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

INTEGER :: n_proc
!     Total number of processors
INTEGER :: error_code

LOGICAL :: at_extremity(4)
!     Indicates if this processor is at north, south
!     east or west of the processor grid

INTEGER :: row_length
INTEGER :: rows
INTEGER :: n_rad_layers

LOGICAL :: L_rad_deg
!     Controls the spatial degradation of E-S code

LOGICAL :: L_subsample
!     Flag to apply spatial subsampling (for satellite footprints)
!     in radiation.
LOGICAL :: L_geostationary
!     Flag to signal that a geostationary satellite is assumed.
LOGICAL :: l_spherical_solar
!     Use spherical geometry for the direct beam

LOGICAL :: rad_mask(row_length, rows)
!     A mask which ensures a chequerboard pattern of radiation
!     calculations over the whole domain (not just one PE)
LOGICAL :: switch(row_length,rows)

INTEGER :: j_sw

INTEGER :: a_sw_radstep
INTEGER :: timestep_number


INTEGER :: daylight_points
INTEGER :: tot_daylight_points
!     Total number of daylight points in whole domain
INTEGER :: List_daylight_points(row_length*rows)


INTEGER :: first_data_interp
!     The first data point, in data co-ords, which needs to be
!     interpolated on a PE (for use in interpolation routine)

INTEGER :: diag_row_list(row_length*rows)
!     List of row indices of points where diagnostics are
!     calculated
INTEGER :: diag_col_list(row_length*rows)
!     List of column indices of points where diagnostics are
!     calculated

REAL :: min_view_lon
!     Minimum longitude of viewing domain
REAL :: max_view_lon
!     Maximum longitude of viewing domain
REAL :: min_view_lat
!     Minimum latitude of viewing domain
REAL :: max_view_lat
!     Maximum latitude of viewing domain

REAL :: seconds_since_midnight
REAL :: true_longitude(row_length, rows)
REAL :: true_latitude(row_length, rows)

REAL :: day_fraction(row_length,rows)
REAL :: day_frac_sph(row_length,rows, 0:n_rad_layers+1)

! Local Variables

LOGICAL :: l_viewed
LOGICAL,PARAMETER :: l_sw = .TRUE.  ! called from the SW code.

INTEGER :: i,j

INTEGER :: interp_index
!         Variable which alternates between 0 and 1 every other
!         radiation timestep

! External routines

LOGICAL :: in_footprint

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRELIM_SWRAD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (L_rad_deg) THEN
  !
  interp_index=MOD((timestep_number/a_sw_radstep),2)
  !
  IF ( (interp_index == 0) .eqv. rad_mask(1,1)) THEN
    first_data_interp = 1
  ELSE
    first_data_interp = 0
  END IF
  DO j=1, rows
    DO i=1, row_length
      switch(i,j) = ( (day_fraction(i,j) > 0.0) .AND.      &
                    (rad_mask(i,j) .eqv.                   &
                    (interp_index == 0) ) )
    END DO
  END DO
  !
ELSE   ! Spatial degradation is switched off
  IF (l_spherical_solar) THEN
    ! SW calculations performed where the top layer is lit
    DO j = 1, rows
      DO i = 1, row_length
        switch(i,j) = day_frac_sph(i,j,n_rad_layers) >  0.0
      END DO
    END DO
  ELSE
    DO j = 1, rows
      DO i = 1, row_length
        switch(i,j) = day_fraction(i,j) >  0.0
      END DO
    END DO
  END IF
END IF   ! If l_rad_deg
!

daylight_points = 0
DO j = 1, rows
  DO i = 1, row_length
    IF (L_subsample) THEN
      ! DEPENDS ON: in_footprint
      l_viewed = in_footprint(j_sw, l_sw,                  &
        L_geostationary,                                   &
        min_view_lon, max_view_lon,                        &
        min_view_lat, max_view_lat,                        &
        true_latitude(i,j),                                &
        true_longitude(i,j),                               &
        seconds_since_midnight)
    ELSE
      l_viewed = .TRUE.
    END IF
    IF (switch(i,j) .AND. l_viewed) THEN
      daylight_points = daylight_points + 1
      list_daylight_points(daylight_points) =              &
                                  (j-1)*row_length + i
      !
      ! The following arrays are 2-D analogues of the above
      ! used for diagnostic purposes from 5.3.
      !
      diag_row_list(daylight_points) = j
      diag_col_list(daylight_points) = i
    END IF
  END DO
END DO
!
! Add up the total number of daylight points in the whole domain using
! the GC routine ISUM (integer sum).
!
IF (L_rad_deg) THEN
  tot_daylight_points=daylight_points
  CALL gc_isum(1,n_proc,Error_code,tot_daylight_points)
  IF ( error_code /= 0 ) THEN

    CALL Ereport('NI_rad_ctl', Error_code,                 &
                 'Unable to determine total lit points.')
  END IF
ELSE
  tot_daylight_points=0
END IF   ! If l_rad_deg
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
!
END SUBROUTINE prelim_swrad
END MODULE prelim_swrad_mod

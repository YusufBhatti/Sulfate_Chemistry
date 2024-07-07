! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: This part of the code has been removed from RAD_CTL2
!              in order to make that deck more readable. This part of
!              the code was contained in RAD_CTL2 up to UM 6.1
!
! Purpose: This routine mainly deals with the spatial degradation and
!          the use of segments
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!-----------------------------------------------------------------------
!
MODULE prelim_lwrad_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'PRELIM_LWRAD_MOD'
CONTAINS

SUBROUTINE prelim_lwrad(                                           &
! Parallel variables
       at_extremity, n_proc,                                         &
! Model Dimensions
       row_length,rows,model_levels,                                 &
! Model Switches
       l_rad_deg, l_extra_top_lw, l_subsample, l_geostationary,      &
! Time stepping Information
       timestep_number,a_lw_radstep,                                 &
! ancillary fields and fields needed to be kept from timestep to
! timestep
       lw_incs,                                                      &
! Satelitte Geometry
       min_view_lon, max_view_lon,                                   &
       min_view_lat, max_view_lat,                                   &
! Number of call
       j_lw,                                                         &
! Other variables
       true_latitude,true_longitude, seconds_since_midnight,         &
       rad_mask,list_lw_points,first_data_interp,                    &
       first_row,last_row,diag_row_list,diag_col_list,               &
       lw_points, olr, lw_down, lwsea, top_absorption)

! Declare variables

USE um_parparams, ONLY: pnorth, psouth

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

INTEGER :: N_proc
!     Total number of processors

LOGICAL :: at_extremity(4)
!     Indicates if this processor is at north, south,
!     east or west of the processor grid.

INTEGER :: row_length
INTEGER :: rows
INTEGER :: model_levels

LOGICAL :: L_rad_deg
!     Controls the spatial degradation of E-S code

LOGICAL :: L_subsample
!     Flag to apply spatial subsampling (for satellite footprints)
!     in radiation.
LOGICAL :: L_geostationary
!     Flag to signal that a geostationary satellite is assumed.

LOGICAL :: L_extra_top_lw
!     Flag to use an extra top layer in radiative calculations

LOGICAL :: rad_mask(row_length, rows)
!     A mask which ensures a chequerboard pattern of radiation
!     calculations over the whole domain (not just one PE)

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


REAL :: LW_incs(row_length, rows, 0:model_levels)
REAL :: olr(row_length, rows)
REAL :: lw_down(row_length, rows)
REAL :: LWsea(row_length, rows)
REAL :: top_absorption(row_length, rows)

INTEGER :: j_lw

INTEGER :: timestep_number
INTEGER :: a_lw_radstep

INTEGER :: list_lw_points(row_length*rows)
INTEGER :: lw_points
!     Variables to enable scatter/gather of LW (like SW):
!     creates a LIST of points to do a calculation, and
!     also facilitates segmenting
INTEGER :: first_data_interp
!     The first data point, in data co-ords, which needs to be
!     interpolated on a PE (for use in interpolation routine)
INTEGER :: first_row,last_row

INTEGER :: diag_row_list(row_length*rows)
!     List of row indices of points where diagnostics are
!     calculated
INTEGER :: diag_col_list(row_length*rows)
!     List of column indices of points where diagnostics are
!     calculated

! Local Variables

LOGICAL :: l_viewed
LOGICAL,PARAMETER :: l_sw = .FALSE.  ! called from the SW code.

INTEGER :: i,j

INTEGER :: interp_index
!     Variable which alternates between 0 and 1 every other
!     radiation timestep

! External Routines

LOGICAL :: in_footprint

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRELIM_LWRAD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
lw_points = 0
IF (L_rad_deg) THEN
  interp_index=MOD((timestep_number/a_lw_radstep),2)
  !
  IF ( ( interp_index == 0 ) .eqv. rad_mask(1,1) ) THEN
    first_data_interp = 1
  ELSE
    first_data_interp = 0
  END IF
  !
  ! Calculations need only be performed at one point on polar
  ! rows, so special treatment is required in these cases.
  !
  first_row = 1
  last_row  = rows
  !
  DO j=first_row, last_row
    DO i=1, row_length
      IF ( rad_mask(i,j) .eqv. ( interp_index == 0 ) ) THEN
        IF (L_subsample) THEN
          ! DEPENDS ON: in_footprint
          l_viewed = in_footprint(j_lw, l_sw,               &
            L_geostationary,                                &
            min_view_lon,max_view_lon,                      &
            min_view_lat,max_view_lat,                      &
            true_latitude(i,j),                             &
            true_longitude(i,j),                            &
            seconds_since_midnight)
        ELSE
          l_viewed=.TRUE.
        END IF
        IF (l_viewed) THEN
          lw_points                 = lw_points + 1
          list_lw_points(lw_points) = i+(j-1)*row_length
        END IF
      END IF
    END DO
  END DO
  !

  !
  ! Initialize all output fields.
  !
  olr(:,:)       = 0.0
  lw_down(:,:)   = 0.0
  LWsea(:,:)     = 0.0
  LW_incs(:,:,:) = 0.0
  IF (l_extra_top_lw) THEN
    top_absorption(:,:)=0.0
  END IF
  !
  !
ELSE  ! SPATIAL DEGRADATION IS SWITCHED OFF
  !
  ! Calculations need only be performed at one point on polar
  ! rows, so special treatment is required in these cases.
  !
  first_row = 1
  last_row  = rows
  !
  DO j=first_row, last_row
    DO i=1, row_length
      IF (L_subsample) THEN
        ! DEPENDS ON: in_footprint
        l_viewed = in_footprint(j_lw, l_sw,                     &
            L_geostationary,                                    &
            min_view_lon,max_view_lon,                          &
            min_view_lat,max_view_lat,                          &
            true_latitude(i,j),                                 &
            true_longitude(i,j),                                &
            seconds_since_midnight)
      ELSE
        l_viewed=.TRUE.
      END IF
      IF (l_viewed) THEN
        lw_points                 = lw_points + 1
        list_lw_points(lw_points) = i+(j-1)*row_length
      END IF
    END DO
  END DO
  !
  !
END IF  ! If L_rad_deg block
!
! Infer the row and column indices of the points where
! calculations are required. This is very mildly inefficient,
! but is cleaner than repeating code throughout the preceding
! block.
!
DO j=1, lw_points
  diag_row_list(j) = (LIST_LW_points(j) - 1)/row_length + 1
  diag_col_list(j) = LIST_LW_points(j) - row_length                &
       * (diag_row_list(j) - 1)
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
!
END SUBROUTINE prelim_lwrad
END MODULE prelim_lwrad_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_calc_coords_mod

IMPLICIT NONE

!  Subroutine rcf_calc_coords

! Description:
!   Calculate the co-ordinates of each (ENDGame) grid cell at p and u points.
!   Suitable for use with both longitude-latitude and Cartesian co-ordinates.
!   All 4 output arguments are OPTIONAL; request only arrays that are needed.
!   No special treatment is given to bi/cyclic LAM grids.
!
!   CAUTION:
!   In the reconfiguration xi1_u_global is defined over 1:glob_u_row_length,
!   but the equivalent model array runs from 0:glob_u_row_length. Therefore
!   what the reconfiguration stores as xi1_u_global(1,1) would be
!   xi1_u_global(0,1) in the model, and the reconfiguration does not store the
!   value of what would be xi1_u_global(glob_u_row_length,1) in the model.
!   Likewise xi2_v_global is zero-indexed in the model but starts at 1 here.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CALC_COORDS_MOD'

CONTAINS

SUBROUTINE rcf_calc_coords (hdr, grid, xi1_p_global, xi2_p_global,    &
                            xi1_u_global, xi2_v_global)

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE rcf_grid_type_mod, ONLY: &
    grid_type

USE rcf_headaddress_mod, ONLY:                   &
    fh_gridstagger_endgame,     FH_CoordSystem,  &
    FH_CoordSystem_Cartesian

USE latlon_eq_rotation_mod, ONLY: &
    rotate_eq_to_latlon

USE rcf_headaddress_mod, ONLY:                   &
  rc_polelong,                  rc_longspacing,  &
  rc_latspacing,                rc_polelat,      &
  rc_firstlong,                 rc_firstlat

USE lam_config_inputs_mod, ONLY: &
    frstlona,    frstlata,       &
    delta_lat

USE ereport_mod, ONLY: &
    ereport

USE conversions_mod, ONLY: &
    pi, pi_over_180

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE(um_header_type), INTENT(IN) :: hdr
TYPE (grid_type), INTENT(IN)     :: grid
REAL, INTENT(OUT), OPTIONAL      :: xi1_p_global(grid % glob_p_row_length, &
                                                 grid % glob_p_rows)
                                                 ! lon at global p-points
REAL, INTENT(OUT), OPTIONAL      :: xi2_p_global(grid % glob_p_row_length, &
                                                 grid % glob_p_rows)
                                                 ! lat at global p-points
REAL, INTENT(OUT), OPTIONAL      :: xi1_u_global(grid % glob_u_row_length, &
                                                 grid % glob_u_rows)
                                                 ! lon at global u-points
REAL, INTENT(OUT), OPTIONAL      :: xi2_v_global(grid % glob_v_row_length, &
                                                 grid % glob_v_rows)
                                                 ! lat at global v-points

! Local variables
INTEGER                           :: i,j
INTEGER                           :: errorstatus
CHARACTER (LEN=errormessagelength):: cmessage       ! used for EReport
CHARACTER (LEN=*), PARAMETER      :: routinename='RCF_CALC_COORDS'

! Temporary arrays for optional outputs that are always required internally:
REAL :: xi1_p_global_tmp(grid % glob_p_row_length, grid % glob_p_rows)
                                       ! lon at global p-points
REAL :: xi2_p_global_tmp(grid % glob_p_row_length, grid % glob_p_rows)
                                       ! lat at global p-points

INTEGER(KIND=jpim), PARAMETER     :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER     :: zhook_out = 1
REAL(KIND=jprb)                   :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (grid % grid_stagger /= fh_gridstagger_endgame) THEN
  ! Print error if grid staggering is not catered for in this routine.
  cmessage = 'Grid staggering method is not supported.'
  errorstatus= 10
  CALL ereport( routinename, errorstatus, cmessage )
END IF

!--------------------------------------------------------------------
! Set up a latitude/longitude field for the grid.
!--------------------------------------------------------------------

DO j = 1, grid % glob_p_rows
  DO i = 1, grid % glob_p_row_length

    xi1_p_global_tmp(i,j) = hdr % realc(rc_firstlong) +      &
                            (i-0.5) * hdr % realc(rc_longspacing)

    xi2_p_global_tmp(i,j) = hdr % realc(rc_firstlat) +       &
                            (j-0.5) * hdr % realc(rc_latspacing)

  END DO
END DO

!--------------------------------------------------------------------
! For rotated grids, get true lats and longs.
!--------------------------------------------------------------------

IF (grid % rotated) THEN

  ! Use the same arrays to store true lats & longs

  CALL rotate_eq_to_latlon (xi2_p_global_tmp, xi1_p_global_tmp, &
                            xi2_p_global_tmp, xi1_p_global_tmp, &
                            hdr % realc(rc_polelat),            &
                            hdr % realc(rc_polelong),           &
                            grid % glob_p_field)

END IF

! Check the grid isn't Cartesian:
IF (hdr % fixhd( FH_CoordSystem ) /= FH_CoordSystem_Cartesian) THEN
  ! Convert to radians
  DO j = 1, grid % glob_p_rows
    DO i = 1, grid % glob_p_row_length
      xi1_p_global_tmp(i,j) = xi1_p_global_tmp(i,j) * pi_over_180
      xi2_p_global_tmp(i,j) = xi2_p_global_tmp(i,j) * pi_over_180
    END DO
  END DO
END IF

! Copy the temporary arrays back to the dummy arguments if required:
IF (PRESENT(xi1_p_global)) THEN
  DO j = 1, grid % glob_p_rows
    DO i = 1, grid % glob_p_row_length
      xi1_p_global(i,j) = xi1_p_global_tmp(i,j)
    END DO
  END DO
END IF

IF (PRESENT(xi2_p_global)) THEN
  DO j = 1, grid % glob_p_rows
    DO i = 1, grid % glob_p_row_length
      xi2_p_global(i,j) = xi2_p_global_tmp(i,j)
    END DO
  END DO
END IF


IF (PRESENT(xi1_u_global)) THEN

  DO i = 2, grid % glob_u_row_length
    DO j = 1, grid % glob_u_rows
      xi1_u_global(i,j) = 0.5*(xi1_p_global_tmp(i-1,j) + xi1_p_global_tmp(i,j))
    END DO
  END DO

  IF (grid % global) THEN
    i = 1
    DO j = 1, grid % glob_u_rows
      xi1_u_global(i,j) = 0.5 * (xi1_p_global_tmp(i,j) +                     &
                           xi1_p_global_tmp(grid%glob_p_row_length,j) -2.0*pi)
    END DO
  ELSE
    i = 1
    DO j = 1, grid % glob_u_rows
      xi1_u_global(i,j) = frstlona
    END DO
  END IF

END IF

IF (PRESENT(xi2_v_global)) THEN

  DO i = 1, grid % glob_v_row_length
    DO j = 2, grid % glob_v_rows-1

      xi2_v_global(i,j) = 0.5*(xi2_p_global_tmp(i,j-1) + xi2_p_global_tmp(i,j))

    END DO
  END DO

  IF (grid % global) THEN
    DO i = 1, grid % glob_v_row_length
      j=1
      xi2_v_global(i,j) = -pi/2.0

      j = grid % glob_v_rows
      xi2_v_global(i,j) = pi/2.0
    END DO
  ELSE
    DO i = 1, grid % glob_v_row_length
      j = 1
      xi2_v_global(i,j) = frstlata

      j = grid % glob_v_rows
      xi2_v_global(i,j) = xi2_v_global(i,j-1) + delta_lat
    END DO
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_calc_coords
END MODULE rcf_calc_coords_mod

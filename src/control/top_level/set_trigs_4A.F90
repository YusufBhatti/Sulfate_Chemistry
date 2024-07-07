! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Initialises trigonometric and Coriolis dynamical arrays (ENDGAME version)
!
! Subroutine Interface:
SUBROUTINE set_trigs_4A(                                          &
                     row_length, rows, n_rows,model_levels,       &
                     base_lambda, base_phi,                       &
                     delta_lambda, delta_phi,                     &
                     base_lambda_wk, base_phi_wk,                 &
                     delta_lambda_wk, delta_phi_wk,               &
                     f_plane_rad, ff_plane_rad,                   &
                     two_Omega,                                   &
                     lat_rot_NP, long_rot_NP,                     &
                     lat_rot_NP_deg, long_rot_NP_deg,             &
                     lat_n, lat_s, long_e, long_w,                &
                     long_e_model, long_w_model, rmdi)

USE trignometric_mod
USE rimtypes
USE dyn_coriolis_mod
USE dyn_var_res_Mod,    ONLY: glambda_p, phi_p
USE rot_coeff_mod,      ONLY: rot_coeff1, rot_coeff2
USE conversions_mod,    ONLY: pi_over_180

USE swap_bounds_2d_mv_mod, ONLY: swap_bounds_2d_mv
USE swapable_field_mod, ONLY: swapable_field_pointer_type

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE eg_alpha_mod
USE eg_alpha_ramp_mod
USE horiz_grid_mod

USE Control_Max_Sizes
USE UM_ParVars
USE UM_ParParams, ONLY: pnorth, psouth
USE Field_Types,  ONLY: fld_type_p, fld_type_u, fld_type_v
USE ereport_mod, ONLY: ereport

USE atm_fields_bounds_mod
USE coriolis_mod

USE latlon_eq_rotation_mod, ONLY: rotate_eq_to_latlon, eq_latlon_vector_coeffs

USE model_domain_mod, ONLY: l_regular, model_type, mt_global
USE idealise_run_mod, ONLY: f_plane, ff_plane, l_shallow, l_vert_coriolis

IMPLICIT NONE

!
! Description:
!   Initialises trigonometric and Coriolis dynamic arrays
!
! Method:
!   Sets values from grid information
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

! Input variables.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_TRIGS_4A'

INTEGER  ::  row_length
INTEGER  ::  rows
INTEGER  ::  n_rows
INTEGER  ::  model_levels

REAL  :: base_lambda
REAL  :: base_phi
REAL  :: delta_lambda
REAL  :: delta_phi
REAL  :: f_plane_rad
REAL  :: ff_plane_rad
REAL  :: lat_rot_NP_deg
REAL  :: long_rot_NP_deg
REAL  :: lat_rot_NP
REAL  :: long_rot_NP
REAL  :: base_lambda_wk ! Longitude of first theta point in degs
REAL  :: base_phi_wk    ! Latitude of first theta point in degs
REAL  :: delta_lambda_wk ! EW (x) grid spacing in degrees
REAL  :: delta_phi_wk    ! NS (y) grid spacing in degrees
REAL  :: rmdi            ! For unset (missing) data
REAL  :: two_Omega

REAL :: lat_n
REAL :: lat_s
REAL :: long_e
REAL :: long_w
REAL :: long_w_model
REAL :: long_e_model

! Local variables.
INTEGER  ::  i, j
INTEGER  ::  gi, gj
INTEGER  ::  j0, j1

REAL  ::  f1_temp        ! f1 term for f-plane or ff-plane
REAL  ::  f2_temp        ! f2 term for f-plane or ff-plane
REAL  ::  f3_temp        ! f3 term for f-plane or ff-plane
REAL  ::  temp1 ! Used in calculation of Coriolis terms (LAM only)
REAL  ::  temp2 ! Used in calculation of Coriolis terms (LAM only)

REAL  ::  latitude_rot(row_length, n_rows)
REAL  ::  longitude_rot(row_length, n_rows)
REAL  ::  true_latitude_b(row_length, n_rows)
                                       ! true latitude on B-grid
REAL  ::  true_longitude_b(row_length,n_rows)
                                       ! true longitude on B-grid


INTEGER :: i_field

TYPE(swapable_field_pointer_type) :: fields_to_swap(10)

! -----------------------------------------
! 1. set trig fields
! -----------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Allocate latitude arrays for trignometric_mod module
IF (.NOT. ALLOCATED(cos_theta_latitude)) THEN
  ALLOCATE (cos_theta_latitude (tdims_s%i_start:tdims_s%i_end,      &
                                tdims_s%j_start:tdims_s%j_end))
END IF
IF (.NOT. ALLOCATED(sec_theta_latitude)) THEN
  ALLOCATE (sec_theta_latitude (tdims_s%i_start:tdims_s%i_end,      &
                                tdims_s%j_start:tdims_s%j_end))
END IF
IF (.NOT. ALLOCATED(FV_cos_theta_latitude)) THEN
  ALLOCATE (FV_cos_theta_latitude (tdims_s%i_start:tdims_s%i_end,   &
                                   tdims_s%j_start:tdims_s%j_end))
END IF
IF (.NOT. ALLOCATED(FV_sec_theta_latitude)) THEN
  ALLOCATE (FV_sec_theta_latitude (tdims_s%i_start:tdims_s%i_end,   &
                                   tdims_s%j_start:tdims_s%j_end))
END IF
IF (.NOT. ALLOCATED(sin_theta_latitude)) THEN
  ALLOCATE (sin_theta_latitude (row_length, rows))
END IF
IF (.NOT. ALLOCATED(tan_theta_latitude)) THEN
  ALLOCATE (tan_theta_latitude (row_length, rows))
END IF
IF (.NOT. ALLOCATED(sin_v_latitude)) THEN
  ALLOCATE (sin_v_latitude (vdims%i_start:vdims%i_end,              &
                            vdims%j_start:vdims%j_end))
END IF
IF (.NOT. ALLOCATED(tan_v_latitude)) THEN
  ALLOCATE (tan_v_latitude (vdims%i_start:vdims%i_end,              &
                            vdims%j_start:vdims%j_end))
END IF
IF (.NOT. ALLOCATED(cos_v_latitude)) THEN
  ALLOCATE (cos_v_latitude (vdims_s%i_start:vdims_s%i_end,          &
                            vdims_s%j_start:vdims_s%j_end))
END IF
IF (.NOT. ALLOCATED(sec_v_latitude)) THEN
  ALLOCATE (sec_v_latitude (vdims_s%i_start:vdims_s%i_end,          &
                            vdims_s%j_start:vdims_s%j_end))
END IF

! Allocate longitude and 'true' arrays for trignometric_mod module
IF (.NOT. ALLOCATED(cos_theta_longitude)) THEN
  ALLOCATE (cos_theta_longitude (row_length, rows))
END IF
IF (.NOT. ALLOCATED(sin_theta_longitude)) THEN
  ALLOCATE (sin_theta_longitude (row_length, rows))
END IF
IF (.NOT. ALLOCATED(cos_u_longitude)) THEN
  ALLOCATE (cos_u_longitude (udims%i_start:udims%i_end,             &
                             udims%j_start:udims%j_end))
END IF
IF (.NOT. ALLOCATED(sin_u_longitude)) THEN
  ALLOCATE (sin_u_longitude (udims%i_start:udims%i_end,             &
                             udims%j_start:udims%j_end))
END IF


IF (.NOT. cartesian_grid) THEN
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      cos_theta_latitude(i,j) = COS(xi2_p(j))
      sin_theta_latitude(i,j) = SIN(xi2_p(j))
      FV_cos_theta_latitude(i,j) = cos_theta_latitude(i,j)
    END DO
  END DO
ELSE
  cos_theta_latitude    = 1.0
  sin_theta_latitude    = 0.0
  FV_cos_theta_latitude = 1.0
END IF



DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    sec_theta_latitude(i, j) = 1.0 / cos_theta_latitude(i, j)
    tan_theta_latitude(i, j) = sin_theta_latitude(i,j) /            &
                                 cos_theta_latitude(i,j)
  END DO
END DO


DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    FV_sec_theta_latitude(i, j) = 1.0/FV_cos_theta_latitude(i,j)
  END DO
END DO

IF (.NOT. cartesian_grid) THEN

  DO j = vdims%j_start,vdims%j_end
    DO i = vdims%i_start,vdims%i_end
      sin_v_latitude(i,j) = SIN(xi2_v(j))
      cos_v_latitude(i,j) = COS(xi2_v(j))
      sec_v_latitude(i,j) = 1.0/cos_v_latitude(i,j)
      tan_v_latitude(i,j) = sin_v_latitude(i,j) /                   &
                              cos_v_latitude(i,j)
    END DO
  END DO

  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      sin_theta_longitude(i,j) = SIN(xi1_p(i))
      cos_theta_longitude(i,j) = COS(xi1_p(i))
    END DO
  END DO

  DO j = udims%j_start,udims%j_end
    DO i = udims%i_start,udims%i_end
      sin_u_longitude(i,j)     = SIN(xi1_u(i))
      cos_u_longitude(i,j)     = COS(xi1_u(i))
    END DO
  END DO

ELSE

  DO j = vdims%j_start,vdims%j_end
    DO i = vdims%i_start,vdims%i_end
      sin_v_latitude(i,j) = 0.0
      cos_v_latitude(i,j) = 1.0
      sec_v_latitude(i,j) = 1.0
      tan_v_latitude(i,j) = 0.0
    END DO
  END DO

  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      sin_theta_longitude(i,j) = 0.0
      cos_theta_longitude(i,j) = 1.0
    END DO
  END DO

  DO j = udims%j_start,udims%j_end
    DO i = udims%i_start,udims%i_end
      sin_u_longitude(i,j)     = 0.0
      cos_u_longitude(i,j)     = 1.0
    END DO
  END DO

END IF

! -----------------------------------------
! 2. set Coriolis components
! -----------------------------------------

! ENDGAME only
!----------------------------------------------------------------------
!  Initialize fi_star coriolis terms: fi_star(k) = fi where fi is the
!  actual Coriolis term. Will be converted to fi_star after the
!  computing the metric terms.
!----------------------------------------------------------------------

IF (model_type == mt_global) THEN

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      IF ( L_shallow ) THEN
        f1_comp(i,j) = 0.0
        f2_comp(i,j) = 0.0 !cos_theta_latitude(i,j)
        f3_comp(i,j) = two_omega*SIN(xi2_p(j)) !sin_theta_latitude(i,j)
      ELSE
        f1_comp(i,j) = 0.0
        f2_comp(i,j) = two_omega*COS(xi2_p(j)) !cos_theta_latitude(i,j)
        f3_comp(i,j) = two_omega*SIN(xi2_p(j)) !sin_theta_latitude(i,j)
      END IF

    END DO
  END DO

ELSE

  ! limited area model
  temp1 = COS( lat_rot_NP)
  temp2 = SIN( lat_rot_NP)

  IF ( cartesian_grid .OR. f_plane  >   -89.0) THEN

    f_plane_rad = f_plane * pi_over_180
    ff_plane_rad = ff_plane * pi_over_180
    f1_temp = - two_omega * temp1 * SIN(ff_plane_rad)
    f2_temp = two_omega * (COS(f_plane_rad) * temp2                 &
               - SIN(f_plane_rad) * temp1 * COS(ff_plane_rad))
    f3_temp = two_omega * ( SIN(f_plane_rad) * temp2                &
               + COS(f_plane_rad) * temp1 * COS(ff_plane_rad))
    IF ( L_vert_Coriolis) THEN
      f1_temp = 0.0
      f2_temp = 0.0
      f3_star = two_omega * ( SIN(f_plane_rad) * temp2              &
              + COS(f_plane_rad) * temp1 * COS(ff_plane_rad))
    END IF  !  L_vert_Coriolis

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        f1_comp(i,j) = f1_temp
        f2_comp(i,j) = f2_temp
        f3_comp(i,j) = f3_temp
      END DO
    END DO

  ELSE      !  cartesian_grid = .false. or f_plane <= -89

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        f1_comp(i,j) = - two_omega * temp1 *                        &
                         sin_theta_longitude(i,j)
        f2_comp(i,j) = two_omega*(cos_theta_latitude(i,j)*          &
                           temp2 -sin_theta_latitude(i,j)*temp1     &
                                  *cos_theta_longitude(i,j) )
        f3_comp(i,j) = two_omega*(sin_theta_latitude(i,j)*          &
                           temp2 +cos_theta_latitude(i,j)*temp1     &
                                  *cos_theta_longitude(i,j) )
      END DO
    END DO

  END IF !  cartesian_grid or f_plane > -89

END IF ! If (model_type == mt_global)
! End of ENDGAME-only

! Allocate arrays for dyn_coriolis_mod module
IF (.NOT. ALLOCATED(f1_at_v)) THEN
  ALLOCATE (f1_at_v (vdims_s%i_start:vdims_s%i_end,                 &
                     vdims_s%j_start:vdims_s%j_end ))
END IF
IF (.NOT. ALLOCATED(f2_at_u)) THEN
  ALLOCATE (f2_at_u (udims_s%i_start:udims_s%i_end,                 &
                     udims_s%j_start:udims_s%j_end ))
END IF
IF (.NOT. ALLOCATED(f3_at_u)) THEN
  ALLOCATE (f3_at_u (udims_s%i_start:udims_s%i_end,                 &
                     udims_s%j_start:udims_s%j_end))
END IF
IF (.NOT. ALLOCATED(f3_at_v)) THEN
  ALLOCATE (f3_at_v (vdims_s%i_start:vdims_s%i_end,                 &
                     vdims_s%j_start:vdims_s%j_end))
END IF
IF (.NOT. ALLOCATED(true_latitude)) THEN
  ALLOCATE (true_latitude (row_length, rows))
END IF
IF (.NOT. ALLOCATED(true_longitude)) THEN
  ALLOCATE (true_longitude (row_length, rows))
END IF

IF (model_type == mt_global) THEN

  DO j = vdims%j_start,vdims%j_end
    DO i = vdims%i_start,vdims%i_end
      f1_at_v(i,j) = 0.0
      f3_at_v(i,j) = two_omega * sin_v_latitude(i,j)
    END DO
  END DO
  DO j = udims%j_start,udims%j_end
    DO i = udims%i_start,udims%i_end
      f2_at_u(i,j) = two_omega * cos_theta_latitude(i+1,j)
      f3_at_u(i,j) = two_omega * sin_theta_latitude(i+1,j)
    END DO
  END DO
  IF (l_regular) THEN
    DO j = 1, rows
      DO i = 1, row_length
        gi = datastart(1) + i - 1
        true_longitude(i,j) = (gi-1)*delta_lambda
      END DO
    END DO
  ELSE !  variable resolution
    DO j = 1, rows
      DO i = 1, row_length
        gi = datastart(1) + i - 1
        true_longitude(i,j) = glambda_p(gi) - base_lambda
      END DO
    END DO
  END IF !  l_regular

  true_latitude = ASIN(sin_theta_latitude)

  ! set polar values of true longitude to be all the same

  IF (at_extremity(PNorth)) THEN
    DO i = 1, row_length
      true_longitude(i,rows) = 0.0
    END DO
  END IF

  IF (at_extremity(PSouth)) THEN
    DO i = 1, row_length
      true_longitude(i,1) = 0.0
    END DO
  END IF

  ! get parameters for common latlonmax within mppac
  ! these need to be in degrees
  IF (l_regular) THEN
    lat_n = Base_phi_wk + (datastart(2)+rows-2) * Delta_phi_wk
    long_e = Base_lambda_wk +                                         &
                (datastart(1)+row_length-2) * Delta_lambda_wk
  ELSE !  variable resolution
    lat_n =  phi_p(1,rows)  / pi_over_180
    gi = datastart(1) + row_length - 1
    long_e = glambda_p(gi) / pi_over_180
  END IF !  l_regular

  lat_s = Base_phi_wk
  long_w       = Base_lambda_wk
  long_w_model = long_w
  long_e_model = long_e

  IF (long_w >  180.0) long_w = long_w-360.0
  IF (long_e >  180.0) long_e = long_e-360.0

  ! Allocate arrays for rot_coeff_mod module (not used for global)
  IF (.NOT. ALLOCATED(rot_coeff1)) THEN
    ALLOCATE (rot_coeff1 ( 1, 1 ))
  END IF
  IF (.NOT. ALLOCATED(rot_coeff2)) THEN
    ALLOCATE (rot_coeff2 ( 1, 1 ))
  END IF

ELSE     ! limited area model

  temp1 = COS(lat_rot_NP)
  temp2 = SIN(lat_rot_NP)

  IF ( cartesian_grid .OR. f_plane  >   -89.0) THEN
    f1_temp = - two_omega * temp1 * SIN(ff_plane_rad)
    f2_temp = two_omega * (COS(f_plane_rad) * temp2                 &
                 - SIN(f_plane_rad) * temp1 * COS(ff_plane_rad))
    f3_temp = two_omega * ( SIN(f_plane_rad) * temp2                &
               + COS(f_plane_rad) * temp1 * COS(ff_plane_rad))
    IF (L_vert_Coriolis) THEN
      f1_temp = 0.0
      f2_temp = 0.0
      f3_temp = two_omega * ( SIN(f_plane_rad) * temp2              &
                 + COS(f_plane_rad) * temp1 * COS(ff_plane_rad))
    END IF  !  L_vert_Coriolis

    DO j = vdims%j_start,vdims%j_end
      DO i = vdims%i_start,vdims%i_end
        f1_at_v(i,j) = f1_temp
        f3_at_v(i,j) = f3_temp
      END DO
    END DO
    DO j = udims%j_start,udims%j_end
      DO i = udims%i_start,udims%i_end
        f2_at_u(i,j) = f2_temp
        f3_at_u(i,j) = f3_temp
      END DO
    END DO

  ELSE

    DO j = vdims%j_start,vdims%j_end
      DO i = vdims%i_start,vdims%i_end
        f1_at_v(i,j) = - two_omega * temp1 *                        &
                         sin_theta_longitude(i,MIN(j+1,rows))
        f3_at_v(i,j) = two_omega * ( sin_v_latitude(i,j) * temp2    &
                                  +cos_v_latitude(i,j) * temp1      &
                                  *cos_theta_longitude(i,MIN(j+1,rows)))
      END DO
    END DO

    DO j = udims%j_start,udims%j_end
      DO i = udims%i_start,udims%i_end
        f2_at_u(i,j) = two_omega * (                                &
                             cos_theta_latitude(i+1,j) * temp2 -    &
                             sin_theta_latitude(i+1,j) * temp1 *    &
                                cos_u_longitude(i,j) )
        f3_at_u(i,j) = two_omega * (                                &
                             sin_theta_latitude(i+1,j) * temp2 +    &
                             cos_theta_latitude(i+1,j) * temp1 *    &
                                cos_u_longitude(i,j) )
      END DO
    END DO

  END IF

  ! calculate true longitude etc. in radians
  IF (cartesian_grid) THEN

    ! Allocate arrays for rot_coeff_mod module
    IF (.NOT. ALLOCATED(rot_coeff1)) THEN
      ALLOCATE (rot_coeff1 ( row_length, n_rows ))
    END IF
    IF (.NOT. ALLOCATED(rot_coeff2)) THEN
      ALLOCATE (rot_coeff2 ( row_length, n_rows ))
    END IF

    DO j = 1, rows
      DO i = 1, row_length
        true_latitude(i,j)    = 0.0
        true_longitude(i,j)   = 0.0
        true_latitude_b(i,j)  = 0.0
        true_longitude_b(i,j) = 0.0
        longitude_rot(i,j)    = 0.0
        latitude_rot(i,j)     = 0.0

        ! Rotation coefficients for wind.
        ! These should never be needed for Cartesian LAMs, as e.g. the
        ! w_eqtoll routines shouldn't be called - but here are some sensible
        ! defaults.
        rot_coeff1(i,j) = 1.0
        rot_coeff2(i,j) = 0.0
      END DO
    END DO

  ELSE

    DO j = 1, rows
      gj = datastart(2) + j - 1
      DO i = 1, row_length
        gi = datastart(1) + i - 1
        longitude_rot(i,j) = (base_lambda + (gi-1) * delta_lambda)  &
                             / pi_over_180
        latitude_rot(i,j) = (base_phi + (gj-1) * delta_phi)         &
                             / pi_over_180
      END DO
    END DO

    CALL rotate_eq_to_latlon(latitude_rot, longitude_rot,             &
                             true_latitude, true_longitude,           &
                             lat_rot_NP_deg, long_rot_NP_deg,         &
                             rows*row_length)

    DO j = 1, rows
      DO i = 1, row_length
        true_longitude(i,j) = true_longitude(i,j) * pi_over_180
        true_latitude(i,j) = true_latitude(i,j) * pi_over_180
      END DO
    END DO

    ! get parameters for common latlonmax within mppac
    ! these need to be in degrees
    lat_n = latitude_rot(1,rows)
    lat_s = latitude_rot(1,1)

    long_w       = longitude_rot(1,1)
    long_w_model = long_w
    long_e       = longitude_rot(row_length,1)
    long_e_model = long_e

    IF (long_w >  180.0) long_w = long_w-360.0
    IF (long_e >  180.0) long_e = long_e-360.0
    IF (long_w_model >= 360.0)long_w_model = long_w_model-360.0
    IF (long_e_model >= 360.0)long_e_model = long_e_model-360.0

    ! calculate lat/longitude for points on equatorial grid for B grid
    DO j = 1, n_rows
      gj = datastart(2) + j - 1
      DO i = 1, row_length
        gi = datastart(1) + i - 1
        longitude_rot(i,j) = (base_lambda + (gi-.5) *               &
                                   delta_lambda) / pi_over_180
        latitude_rot(i,j) = (base_phi + (gj-.5) * delta_phi)        &
                                                / pi_over_180
      END DO
    END DO

    CALL rotate_eq_to_latlon(latitude_rot, longitude_rot,             &
                             true_latitude_b, true_longitude_b,       &
                             lat_rot_NP_deg, long_rot_NP_deg,         &
                             n_rows*row_length)

    ! Calculate rotation coefficients for wind

    ! Allocate arrays for rot_coeff_mod module
    IF (.NOT. ALLOCATED(rot_coeff1)) THEN
      ALLOCATE (rot_coeff1 ( row_length, n_rows ))
    END IF
    IF (.NOT. ALLOCATED(rot_coeff2)) THEN
      ALLOCATE (rot_coeff2 ( row_length, n_rows ))
    END IF

    CALL eq_latlon_vector_coeffs(rot_coeff1, rot_coeff2,              &
                                 true_longitude_b, longitude_rot,     &
                                 lat_rot_NP_deg, long_rot_NP_deg,     &
                                 n_rows*row_length)

  END IF ! cartesian_grid

END IF ! model_type


! ------------------------------------------------
! 3. swap_bounds for fields which have halos
! ------------------------------------------------

i_field = 0

i_field = i_field + 1
fields_to_swap(i_field) % field_2d   => fv_cos_theta_latitude(:,:)
fields_to_swap(i_field) % field_type = fld_type_p
fields_to_swap(i_field) % levels     = 1
fields_to_swap(i_field) % rows       = rows
fields_to_swap(i_field) % vector     = .FALSE.

i_field = i_field + 1
fields_to_swap(i_field) % field_2d   => cos_theta_latitude(:,:)
fields_to_swap(i_field) % field_type = fld_type_p
fields_to_swap(i_field) % levels     = 1
fields_to_swap(i_field) % rows       = rows
fields_to_swap(i_field) % vector     = .FALSE.

i_field = i_field + 1
fields_to_swap(i_field) % field_2d   => sec_theta_latitude(:,:)
fields_to_swap(i_field) % field_type = fld_type_p
fields_to_swap(i_field) % levels     = 1
fields_to_swap(i_field) % rows       = rows
fields_to_swap(i_field) % vector     = .FALSE.

i_field = i_field + 1
fields_to_swap(i_field) % field_2d   => fv_sec_theta_latitude(:,:)
fields_to_swap(i_field) % field_type = fld_type_p
fields_to_swap(i_field) % levels     = 1
fields_to_swap(i_field) % rows       = rows
fields_to_swap(i_field) % vector     = .FALSE.

i_field = i_field + 1
fields_to_swap(i_field) % field_2d   => cos_v_latitude(:,:)
fields_to_swap(i_field) % field_type = fld_type_v
fields_to_swap(i_field) % levels     = 1
fields_to_swap(i_field) % rows       = n_rows
fields_to_swap(i_field) % vector     = .FALSE.

i_field = i_field + 1
fields_to_swap(i_field) % field_2d   => sec_v_latitude(:,:)
fields_to_swap(i_field) % field_type = fld_type_v
fields_to_swap(i_field) % levels     = 1
fields_to_swap(i_field) % rows       = n_rows
fields_to_swap(i_field) % vector     = .FALSE.

i_field = i_field + 1
fields_to_swap(i_field) % field_2d   => f3_at_u(:,:)
fields_to_swap(i_field) % field_type = fld_type_u
fields_to_swap(i_field) % levels     = 1
fields_to_swap(i_field) % rows       = rows
fields_to_swap(i_field) % vector     = .FALSE.

i_field = i_field + 1
fields_to_swap(i_field) % field_2d   => f3_at_v(:,:)
fields_to_swap(i_field) % field_type = fld_type_v
fields_to_swap(i_field) % levels     = 1
fields_to_swap(i_field) % rows       = n_rows
fields_to_swap(i_field) % vector     = .FALSE.

i_field = i_field + 1
fields_to_swap(i_field) % field_2d   => f1_at_v(:,:)
fields_to_swap(i_field) % field_type = fld_type_v
fields_to_swap(i_field) % levels     = 1
fields_to_swap(i_field) % rows       = n_rows
fields_to_swap(i_field) % vector     = .FALSE.

i_field = i_field + 1
fields_to_swap(i_field) % field_2d   => f2_at_u(:,:)
fields_to_swap(i_field) % field_type = fld_type_u
fields_to_swap(i_field) % levels     = 1
fields_to_swap(i_field) % rows       = rows
fields_to_swap(i_field) % vector     = .FALSE.

CALL swap_bounds_2d_mv(fields_to_swap, i_field, row_length,         &
                        offx, offy)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_trigs_4A

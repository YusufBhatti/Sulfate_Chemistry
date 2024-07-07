! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_Friction_Suarez_Held

SUBROUTINE IDL_Friction_Suarez_Held                               &
                        (row_length, rows, n_rows                 &
,                        model_levels, timestep                   &
,                        off_x, off_y                             &
, global_row_length,n_proc, n_procy, proc_row_group,at_extremity  &
,                        friction_level                           &
,                        base_frictional_timescale                &
,                        SuHe_sigma_cutoff, SuHe_fric             &
,                        p, p_star, u, v, R_u, R_v )

! Purpose:
!          Provides temperature relaxation for Suarez Held test problem.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

USE um_parparams, ONLY: nodomain, pnorth, peast, psouth, pwest

USE atm_fields_bounds_mod, ONLY:                                  &
    udims, vdims, wdims, tdims, pdims,                            &
    udims_s, vdims_s, tdims_s, pdims_s

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Field_Types
USE ereport_mod, ONLY: ereport
USE p_to_u_mod, ONLY: p_to_u
USE p_to_v_mod, ONLY: p_to_v
USE polar_row_mean_mod, ONLY: polar_row_mean
USE u_to_p_mod, ONLY: u_to_p
USE v_to_p_mod, ONLY: v_to_p
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

USE model_domain_mod, ONLY: model_type, mt_global

USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

! Parallel setup variables

INTEGER ::            &
   off_x              & ! Size of small halo in i
 , off_y              & ! Size of small halo in j
 , global_row_length  & ! number of points on a row
 , proc_row_group     & ! Group id for processors on the same row
 , n_proc             & ! Total number of processors
 , n_procy              ! Number of processors in latitude


LOGICAL ::  at_extremity(4)
                        ! Indicates if this processor is at north,
                        ! south, east or west of the processor grid

INTEGER ::                                                        &
  row_length                                                      &
                   ! number of points on a row
, rows                                                            &
                   ! number of rows in a u field
, n_rows                                                          &
                     ! number of rows in a v field
, model_levels
                   ! number of model levels


REAL ::                                                           &
  timestep

REAL ::                                                           &
  p(pdims_s%i_start:pdims_s%i_end,                                &
    pdims_s%j_start:pdims_s%j_end,                                &
         model_levels)                                            &
, p_star(row_length, rows)                                        &
, u( udims_s%i_start:udims_s%i_end,                               &
     udims_s%j_start:udims_s%j_end,                               &
        model_levels)                                             &
, v( vdims_s%i_start:vdims_s%i_end,                               &
     vdims_s%j_start:vdims_s%j_end,                               &
        model_levels)                                             &
, R_u(udims_s%i_start:udims_s%i_end,                              &
      udims_s%j_start:udims_s%j_end,                              &
        model_levels)                                             &
, R_v(vdims_s%i_start:vdims_s%i_end,                              &
      vdims_s%j_start:vdims_s%j_end,                              &
        model_levels)

! Suarez Held variables
REAL ::                                                           &
  friction_level(model_levels)                                    &
, base_frictional_timescale                                       &
, SuHe_sigma_cutoff

INTEGER ::                                                        &
  SuHe_fric       ! Switch to choose form of friction

! local variables
INTEGER ::                                                        &
  i, j, k                                                         &
, j0, j1                                                          &
, levels_h
INTEGER :: ErrorStatus

REAL ::                                                           &
  temp

REAL ::                                                          &
 work_u_to_p(row_length, rows, model_levels)                     &
, work_v_to_p(row_length, rows, model_levels)                     &
, work_p_to_v(row_length, n_rows, model_levels)                   &
, work_p_halo(pdims_s%i_start:pdims_s%i_end,                      &
             pdims_s%j_start:pdims_s%j_end,                      &
                model_levels)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'IDL_FRICTION_SUAREZ_HELD'
CHARACTER(LEN=errormessagelength)             :: cmessage


! No External routines

! ----------------------------------------------------------------------
! Section 0. Initialise number of working levels levels_h
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!  Expect boundary -layer levels to be no more than half model_levels
levels_h = model_levels / 2

! ----------------------------------------------------------------------
! Section 1.  Calculate friction terms
! ----------------------------------------------------------------------
IF (SuHe_fric  ==  1) THEN
  ! ----------------------------------------------------------------------
  ! Section 1.1 Original simple friction applied at u,v points
  ! ----------------------------------------------------------------------
  j0 = 1
  j1 = rows
  IF (at_extremity(PSouth)) j0 = 2
  IF (at_extremity(PNorth)) j1 = rows - 1

  ! add on to increment field
  DO k = 1, levels_h
    DO j = j0, j1
      DO i = 1, row_length
        R_u(i,j,k) = R_u(i,j,k) - timestep *                      &
                          friction_level(k) * u(i,j,k)
      END DO
    END DO
  END DO

  ! add on to increment field
  DO k = 1, levels_h
    DO j = 1, n_rows
      DO i = 1, row_length
        R_v(i,j,k) = R_v(i,j,k) -  timestep *                     &
                           friction_level(k) * v(i,j,k)
      END DO
    END DO
  END DO

ELSE IF (SuHe_fric  ==  2) THEN
  ! ----------------------------------------------------------------------
  ! Section 1.2 Simple friction evaluated at p points then
  !             interpolated to u,v points
  ! ----------------------------------------------------------------------

  CALL u_to_p(u,                                                    &
                    udims_s%i_start,udims_s%i_end,                  &
                    udims_s%j_start,udims_s%j_end,                  &
                    pdims%i_start,pdims%i_end,                      &
                    pdims%j_start,pdims%j_end,                      &
                    levels_h, at_extremity,work_u_to_p)

  CALL v_to_p(v,                                                    &
                    vdims_s%i_start,vdims_s%i_end,                  &
                    vdims_s%j_start,vdims_s%j_end,                  &
                    pdims%i_start,pdims%i_end,                      &
                    pdims%j_start,pdims%j_end,                      &
                    levels_h, at_extremity,work_v_to_p)

 
  ! Calculate k_v*u and k_v*v
  DO k = 1, levels_h
    DO j = 1, rows
      DO i = 1, row_length
        temp = (p(i,j,k)/p_star(i,j) - SuHe_sigma_cutoff) /       &
                          (1.0 - SuHe_sigma_cutoff)
        work_u_to_p(i,j,k) = base_frictional_timescale *          &
                            MAX(0.0, temp) * work_u_to_p(i,j,k)
        work_v_to_p(i,j,k) = base_frictional_timescale *          &
                            MAX(0.0, temp) * work_v_to_p(i,j,k)
      END DO
    END DO
  END DO

  !
  ! copy u-tendencies into haloed arrays
  DO k = 1, levels_h
    DO j = 1, rows
      DO i = 1, row_length
        work_p_halo(i,j,k) = work_u_to_p(i,j,k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: swap_bounds
  CALL Swap_Bounds(work_p_halo, row_length, rows,                 &
           levels_h, off_x, off_y, fld_type_p,swap_field_is_scalar)
  ! interpolate to u grid
  CALL p_to_u(work_p_halo,                                    &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               udims%i_start,udims%i_end,                     &
               udims%j_start,udims%j_end,                     &
               udims%k_start,levels_h,work_u_to_p)

  ! add on to increment field
  DO k = 1, levels_h
    DO j = 1, rows
      DO i = 1, row_length
        R_u(i,j,k) = R_u(i,j,k) - work_u_to_p(i,j,k) * timestep
      END DO
    END DO
  END DO

  ! copy v-tendencies into haloed arrays
  DO k = 1, levels_h
    DO j = 1, rows
      DO i = 1, row_length
        work_p_halo(i,j,k) = work_v_to_p(i,j,k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: swap_bounds
  CALL Swap_Bounds(work_p_halo, row_length, rows,                 &
           levels_h, off_x, off_y, fld_type_p,swap_field_is_scalar)
  ! interpolate to v grid
  CALL p_to_v(work_p_halo,                                    &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               vdims%i_start,vdims%i_end,                     &
               vdims%j_start,vdims%j_end,                     &
               vdims%k_start,levels_h,work_p_to_v)


  ! set polar rows to common mean value.
  CALL polar_row_mean(                                        &
                  work_p_to_v,                                &
                  vdims%i_start,vdims%i_end,                  &
                  vdims%j_start,vdims%j_end,                  &
                  vdims%k_start,levels_h,                     &
                  global_row_length,                          &
                  n_proc, n_procy, proc_row_group,            &
                  at_extremity)


  ! add on to increment field
  DO k = 1, levels_h
    DO j = 1, n_rows
      DO i = 1, row_length
        R_v(i,j,k) = R_v(i,j,k) - work_p_to_v(i,j,k) * timestep
      END DO
    END DO
  END DO

ELSE
  ErrorStatus = 1
  cmessage = "Value of SuHe_fric not supported - put correct "// &
             "value of 1 or 2 in idealised namelist"
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF !  SuHe_fric  ==  1

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE IDL_Friction_Suarez_Held

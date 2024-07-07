! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_sl_diab_trac_mod
IMPLICIT NONE
CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName = 'EG_SL_DIAB_TRAC_MOD'
CONTAINS

SUBROUTINE eg_sl_diab_trac(                                           &
             row_length, rows, n_rows, model_levels,                  &
             g_i_pe, high_order_scheme, monotone_scheme,              &
             l_high, l_mono,                                          &
             ! dtheta fields
             dtheta_0, dtheta_bl, dtheta_bl_mix, dtheta_bl_LH,        &
             dtheta_conv, dtheta_mic, dtheta_rad,                     &
             dtheta_SW, dtheta_LW, dtheta_slow, dtheta_cld,           & 
             ! errotstatus in case there are problems 
             ! with No of tracers
             errorstatus)

! Purpose:
!          Performs semi-Lagrangian advection of the theta-tracer variables
!          on theta points. 
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics diagnostics
!
! Code Description: Loads in dtheta arrays and specific model domain
!                   variables. Copy all the dtheta into the superarray
!                   super_dtheta_in, calls eg_interpolation_eta to
!                   interpolate these to the departure points of theta,
!                   whose values are hold in superarray super_dtheta_out.

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE atm_fields_bounds_mod,     ONLY: tdims , tdims_s, tdims_l

USE departure_pts_mod,         ONLY:                                  &
    depart_xi1_w, depart_xi2_w, depart_xi3_w

USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE Field_Types

USE nlsizes_namelist_mod,      ONLY: global_rows, global_row_length

USE um_parcore,                ONLY: mype, nproc

USE um_parvars,                ONLY: nproc_x,nproc_y,                 &
                                     at_extremity,gc_proc_row_group,  &
                                     gc_proc_col_group,               &
                                     halo_i, halo_j,                  &
                                     offx, offy, datastart

USE cv_run_mod,                ONLY: l_param_conv

USE cloud_inputs_mod,          ONLY: i_cld_vn
USE pc2_constants_mod,         ONLY: i_cld_pc2


USE level_heights_mod,         ONLY: eta_theta_levels

USE model_domain_mod,          ONLY: model_type, mt_lam

USE dynamics_input_mod,        ONLY: l_sl_bc_correction

USE free_tracers_inputs_mod,   ONLY: l_diab_tr_bl, l_diab_tr_rad,     &
                                     num_dtheta

! Print out and error statements
USE errormessagelength_mod,    ONLY: errormessagelength
USE ereport_mod,               ONLY: ereport
USE umPrintMgr

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.
INTEGER, INTENT(IN) ::                                                &
! number of points on a row
  row_length                                                          &
! number of rows.
, rows                                                                &
! number of v-rows.
, n_rows                                                              &
! number of model levels.
, model_levels                                                        &
! processor on my processor-row holding a given value in i direction
, g_i_pe(1-halo_i:global_row_length+halo_i)                           &
! a code saying which high order scheme to use for interpolation
, high_order_scheme                                                   &
! a code saying which monotone scheme to use for interpolation
, monotone_scheme

LOGICAL, INTENT(IN) ::                                                &
! True, if high order interpolation required
  l_high                                                              &
 ! True, if interpolation required to be monotone
, l_mono

! Prognostic fields for the diabatic-tracer scheme
REAL, INTENT(INOUT) ::                                                &
                    dtheta_0    (tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end)           &

,                   dtheta_bl    (tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end)           &

,                   dtheta_bl_mix(tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end)          &

,                   dtheta_bl_LH(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end)           &

,                   dtheta_conv (tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end)           &


,                   dtheta_mic  (tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end)           &

,                   dtheta_rad  (tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end)           &

,                   dtheta_SW   (tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end)           &

,                   dtheta_LW   (tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end)           &

,                   dtheta_slow  (tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end)           &

,                   dtheta_cld  (tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end)


INTEGER, INTENT(OUT) :: errorstatus   ! Non-zero on exit if error detected.

! Counting variable
INTEGER :: array_size_count
! Loop indices
INTEGER :: i, j, k
! Linear interpolation is used at departure points in this layer and below.
! (Optional argument for subroutine eg_interpolation_eta.)
INTEGER :: k_int_linear

! Work array containing all the prognostics to be advected
! (same dimensions as "work" in eg_sl_theta)
REAL :: super_dtheta_in                                         &
                          (tdims_l%i_start:tdims_l%i_end,       &
                           tdims_l%j_start:tdims_l%j_end,       &
                           tdims_l%k_start:tdims_l%k_end,       &
                           num_dtheta )

! Work array containing all the prognostics after advection
! (same dimensions as the output r_theta_d in eg_sl_theta)
REAL :: super_dtheta_out                                        &
                          (tdims_s%i_start:tdims_s%i_end,       &
                           tdims_s%j_start:tdims_s%j_end,       &
                           tdims_s%k_start:tdims_s%k_end,       &
                           num_dtheta )

!Print out error statements
INTEGER :: icode    ! error code for ereport
! local temporary arrays
CHARACTER(LEN=errormessagelength)       :: cmessage      ! out error message
CHARACTER(LEN=*), PARAMETER  :: routinename='EG_SL_DIAB_TRAC'

! Dr Hook stuff
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Copy fields into a super-array to pass to the advection routine
! ----------------------------------------------------------------------
array_size_count = 0

array_size_count = array_size_count + 1
DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      super_dtheta_in (i,j,k,array_size_count) = dtheta_0(i,j,k)
    END DO
  END DO
END DO

IF (l_diab_tr_bl) THEN 
  ! BL mixing
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        super_dtheta_in (i,j,k,array_size_count) = dtheta_bl_mix(i,j,k)
      END DO
    END DO
  END DO
  
  ! Latent heat contribution
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        super_dtheta_in (i,j,k,array_size_count) = dtheta_bl_LH(i,j,k)
      END DO
    END DO
  END DO
ELSE 
  ! All BL combined
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        super_dtheta_in (i,j,k,array_size_count) = dtheta_bl(i,j,k)
      END DO
    END DO
  END DO
END IF  !end if over l_diab_tr_bl


IF (l_param_conv) THEN
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        super_dtheta_in (i,j,k,array_size_count) = dtheta_conv(i,j,k)
      END DO
    END DO
  END DO
END IF 

array_size_count = array_size_count + 1
DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      super_dtheta_in (i,j,k,array_size_count) = dtheta_mic(i,j,k)
    END DO
  END DO
END DO

IF (l_diab_tr_rad) THEN
  ! SW rad
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        super_dtheta_in (i,j,k,array_size_count) = dtheta_SW(i,j,k)
      END DO
    END DO
  END DO
  ! LW rad
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        super_dtheta_in (i,j,k,array_size_count) = dtheta_LW(i,j,k)
      END DO
    END DO
  END DO
ELSE 
  ! rad combining SW + LW
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        super_dtheta_in (i,j,k,array_size_count) = dtheta_rad(i,j,k)
      END DO
    END DO
  END DO
END IF !end if over l_diab_tr_rad


!Slow physics
array_size_count = array_size_count + 1
DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      super_dtheta_in (i,j,k,array_size_count) = dtheta_slow(i,j,k)
    END DO
  END DO
END DO

!Cloud rebalancing
IF (i_cld_vn == i_cld_pc2) THEN
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        super_dtheta_in (i,j,k,array_size_count) = dtheta_cld(i,j,k)
      END DO
    END DO
  END DO
END IF 

! Stop the model if num_dtheta is not equal to array_size_count
IF (num_dtheta /= array_size_count) THEN
  WRITE(umMessage,'(A)')'Number of dtheta tracers \= Number of tracers advected'
  CALL umPrint(umMessage,src='eg_sl_diab_trac')
  icode   = 1
  WRITE (cmessage,'(A)')'num_dtheta \= array_size_count.' &
          // 'Check value of num_dtheta in diab_tracer_sav'
  CALL ereport(routinename, icode, cmessage)
END IF

! ----------------------------------------------------------------------
! If any fields to handle, advect them!
! ----------------------------------------------------------------------
IF (num_dtheta>=1) THEN

  ! Swap haloes...
  ! DEPENDS ON: swap_bounds
  CALL Swap_Bounds(super_dtheta_in ,                                     &
                   tdims_l%i_len - 2*tdims_l%halo_i,                     &
                   tdims_l%j_len - 2*tdims_l%halo_j,                     &
                   tdims_l%k_len*array_size_count,                       &
                   tdims_l%halo_i, tdims_l%halo_j, fld_type_p,.FALSE.)

  ! Set layers over which linear interpolation is used
  IF (l_sl_bc_correction) THEN
    k_int_linear=2
  ELSE
    k_int_linear=1
  END IF

  ! Interpolate fields to the departure points...
  CALL eg_interpolation_eta_pmf(                                         &
                eta_theta_levels,fld_type_w,                             &
                num_dtheta,                                              &
                row_length, rows, model_levels + 1,                      &
                rows,                                                    &
                row_length, rows, model_levels + 1,                      &
                high_order_scheme, monotone_scheme,                      &
                l_high, l_mono,                                          &
                depart_xi3_w,depart_xi1_w,depart_xi2_w,                  &
                mype, nproc, nproc_x, nproc_y,                           &
                halo_i, halo_j,                                          &
                global_row_length, datastart, at_extremity,              &
                g_i_pe, gc_proc_row_group, gc_proc_col_group,            &
                offx, offy, errorstatus,                                 &
                super_dtheta_in, super_dtheta_out,                       &
                k_int_linear_in=k_int_linear)

END IF  ! num_dtheta>=1


! ----------------------------------------------------------------------
! Copy the super-array fields back into their respective variables
! ----------------------------------------------------------------------
array_size_count = 0

array_size_count = array_size_count + 1
DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_0(i,j,k) = super_dtheta_out(i,j,k,array_size_count)
    END DO
  END DO
END DO

IF (l_diab_tr_bl) THEN 
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_bl_mix(i,j,k) = super_dtheta_out(i,j,k,array_size_count)
      END DO
    END DO
  END DO

  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_bl_LH(i,j,k) = super_dtheta_out(i,j,k,array_size_count)
      END DO
    END DO
  END DO

ELSE 
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_bl(i,j,k) = super_dtheta_out(i,j,k,array_size_count)
      END DO
    END DO
  END DO
END IF !end if over l_diab_tr_bl

IF (l_param_conv) THEN
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_conv(i,j,k) = super_dtheta_out(i,j,k,array_size_count)
      END DO
    END DO
  END DO
END IF 

array_size_count = array_size_count + 1
DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_mic(i,j,k) = super_dtheta_out(i,j,k,array_size_count)
    END DO
  END DO
END DO

IF (l_diab_tr_rad) THEN 
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_SW(i,j,k) = super_dtheta_out(i,j,k,array_size_count)
      END DO
    END DO
  END DO
  
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_LW(i,j,k) = super_dtheta_out(i,j,k,array_size_count)
      END DO
    END DO
  END DO  
ELSE 
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_rad(i,j,k) = super_dtheta_out(i,j,k,array_size_count)
      END DO
    END DO
  END DO
END IF !end if over l_diab_tr_rad

array_size_count = array_size_count + 1
DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_slow(i,j,k) = super_dtheta_out(i,j,k,array_size_count)
    END DO
  END DO
END DO

IF (i_cld_vn == i_cld_pc2) THEN
  array_size_count = array_size_count + 1
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_cld(i,j,k) = super_dtheta_out(i,j,k,array_size_count)
      END DO
    END DO
  END DO
END IF 

! END of routine.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_sl_diab_trac

END MODULE eg_sl_diab_trac_mod

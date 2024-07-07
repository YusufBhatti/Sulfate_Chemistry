! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sl_moisture_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SL_MOISTURE_MOD'

CONTAINS
SUBROUTINE eg_sl_moisture(                                            &
              moisture_array_size,                                    &
              row_length, rows, n_rows, model_levels, halo_i,         &
              halo_j, offx, offy,datastart,g_i_pe,high_order_scheme,  &
              monotone_scheme,  l_high, l_mono, L_mcr_rain,           &
              l_mcr_cf2,l_mcr_graup,                                  &
              r_m_v, r_m_cl, r_m_cf, r_m_r, r_m_gr, r_m_cf2,          &
              cf_bulk, cf_liquid, cf_frozen,exner_theta_levels,       &
              exner_star,r_m_v_d, r_m_cl_d, r_m_cf_d,r_m_r_d,         &
              r_m_gr_d, r_m_cf2_d,cf_star, cfl_star, cff_star,        &
              conv_prog_1, conv_prog_2, conv_prog_3, conv_prog_precip,&
              rho_n, rho_np1, error_code )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE um_parcore,        ONLY: mype, nproc
USE nlsizes_namelist_mod,  ONLY: global_row_length
USE um_parvars,        ONLY: nproc_x,nproc_y, at_extremity,           &
                             gc_proc_row_group, gc_proc_col_group

USE timestep_mod,      ONLY: timestep
USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels,        &
                              xi3_at_theta=>r_theta_levels,           &
                              xi3_at_rho=>r_rho_levels,               &
                              xi3_at_u=>r_at_u, xi3_at_v=>r_at_v

USE cloud_inputs_mod, ONLY: l_fixbug_pc2_mixph, i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2
USE cv_run_mod, ONLY: l_conv_prog_group_1, l_conv_prog_group_2,       &
                      l_conv_prog_group_3, l_conv_prog_precip
USE atm_fields_bounds_mod
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE horiz_grid_mod
USE metric_terms_mod
USE departure_pts_mod
USE ereport_mod, ONLY: ereport
USE Field_Types
USE dynamics_input_mod, ONLY: l_sl_bc_correction,                &
                              zlf_conservation_moist_option

USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE enforce_mono_mz_mod
USE eg_mass_conserv_mod

USE atm_step_local, ONLY: L_print_L2norms, cycleno
USE sl_moist_norm_mod
USE lam_config_inputs_mod, ONLY: n_rims_to_do
USE turb_diff_mod, ONLY: norm_lev_start, norm_lev_end
USE umPrintMgr, ONLY: umMessage, umPrint
USE eg_zlf_conservation_mod, ONLY: eg_zlf_conservation
USE eg_zlf_mod, ONLY: zlf_cfl_top_level_moist

IMPLICIT NONE
!
! Description:
!   Find departure point timelevel n dependent quantity R_theta_d.
!
!
! Method: ENDGame formulation version 1.01,
!         section 7.3.

!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SL_MOISTURE'



! Model dimensions

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels

! MPP options

INTEGER, INTENT(IN) ::                                                &
  halo_i,                                                             &
                     ! Size of halo in i.
  halo_j,                                                             &
                     ! Size of halo in j.
  offx,                                                               &
                     ! Size of small halo in i
  offy,                                                               &
                     ! Size of small halo in j.
  datastart(3),                                                       &
                     ! First gridpoints held by this processor.
  g_i_pe(1-halo_i:global_row_length+halo_i)
                     ! processor on my processor-row
                     ! holding a given value in i direction


! Loop index bounds for arrays defined on p, u, v points respectively


INTEGER, INTENT(IN) :: moisture_array_size

! Integer parameters for advection

INTEGER, INTENT(IN) ::                                                &
  high_order_scheme,                                                  &
                     ! a code saying which high order scheme to
                     ! use. 1 = tensor tri-cubic lagrange order
                     ! (j,i,k) no other options available at
                     ! present.
  monotone_scheme
                     ! a code saying which monotone scheme to use.
                     ! 1 = tri-linear order (j,i,k)
                     ! no other options available at present.

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                                &
  l_high,                                                             &
                   ! True, if high order interpolation required.
  l_mono
                   ! True, if interpolation required to be monotone.

LOGICAL, INTENT(IN) :: l_mcr_rain, l_mcr_cf2, l_mcr_graup


! Timelevel n arrival point quantities

REAL, INTENT(IN) ::                                                   &
  r_m_v(row_length,rows,0:model_levels),                              &
  r_m_cl(row_length,rows,0:model_levels),                             &
  r_m_cf(row_length,rows,0:model_levels),                             &
  r_m_r(row_length,rows,0:model_levels),                              &
  r_m_gr(row_length,rows,0:model_levels),                             &
  r_m_cf2(row_length,rows,0:model_levels)

REAL, INTENT(OUT) ::                                                  &
  exner_star(1-offx:row_length+offx,1-offy:rows+offy, model_levels)

REAL, INTENT (INOUT) ::                                               &
        cf_star (row_length, rows, 0:model_levels),                   &
        cfl_star(row_length, rows, 0:model_levels),                   &
        cff_star(row_length, rows, 0:model_levels)

REAL, INTENT (INOUT) :: cf_bulk(  tdims_l%i_start:tdims_l%i_end,      &
                                  tdims_l%j_start:tdims_l%j_end,      &
                                  tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: cf_liquid(tdims_l%i_start:tdims_l%i_end,      &
                                  tdims_l%j_start:tdims_l%j_end,      &
                                  tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: cf_frozen(tdims_l%i_start:tdims_l%i_end,      &
                                  tdims_l%j_start:tdims_l%j_end,      &
                                  tdims_l%k_start:tdims_l%k_end)

REAL, INTENT (INOUT) :: conv_prog_1(      tdims_s%i_start:tdims_s%i_end,  &
                                          tdims_s%j_start:tdims_s%j_end,  &
                                          tdims_s%k_start:tdims_s%k_end )
REAL, INTENT (INOUT) :: conv_prog_2(      tdims_s%i_start:tdims_s%i_end,  &
                                          tdims_s%j_start:tdims_s%j_end,  &
                                          tdims_s%k_start:tdims_s%k_end )
REAL, INTENT (INOUT) :: conv_prog_3(      tdims_s%i_start:tdims_s%i_end,  &
                                          tdims_s%j_start:tdims_s%j_end,  &
                                          tdims_s%k_start:tdims_s%k_end )
REAL, INTENT (INOUT) :: conv_prog_precip( tdims_s%i_start:tdims_s%i_end,  &
                                          tdims_s%j_start:tdims_s%j_end,  &
                                          tdims_s%k_start:tdims_s%k_end )

REAL :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end)


INTEGER :: error_code   ! Non-zero on exit if error detected.

! Timelevel n departure point quantities

REAL, INTENT(OUT) ::                                                  &
  r_m_v_d(row_length,rows,0:model_levels),                            &
  r_m_cl_d(row_length,rows,0:model_levels),                           &
  r_m_cf_d(row_length,rows,0:model_levels),                           &
  r_m_r_d(row_length,rows,0:model_levels),                            &
  r_m_gr_d(row_length,rows,0:model_levels),                           &
  r_m_cf2_d(row_length,rows,0:model_levels)

REAL, INTENT(IN) ::   rho_n(pdims_s%i_start:pdims_s%i_end,            &
                            pdims_s%j_start:pdims_s%j_end,            &
                            pdims_s%k_start:pdims_s%k_end)
REAL, INTENT(IN) :: rho_np1(pdims_s%i_start:pdims_s%i_end,            &
                            pdims_s%j_start:pdims_s%j_end,            &
                            pdims_s%k_start:pdims_s%k_end)

! Local variables

INTEGER :: i,j,k, number_of_inputs, number_of_inputs_zlf,  &
                  number_of_inputs_nozlf
INTEGER :: k_int_linear ! Linear interpolation is used at departure
                        ! points in this layer and below.
                        ! (Optional argument for subroutine
                        !  eg_interpolation_eta.)

LOGICAL :: l_high_in   ! local setting for use in CALLs
LOGICAL :: l_mono_in   ! local setting for use in CALLs

! tmp & dummy arrays


REAL :: super_array(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,    &
                0:model_levels, moisture_array_size )

REAL :: super_array_out(1-offx:row_length+offx,1-offy:rows+offy,        &
                0:model_levels, moisture_array_size)

INTEGER :: counter


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!
! add all the fields that should not be conserved or treated with
! ZLF first in the "super_array" , then in the second part you
! add the tracers that are treated with ZLF (or conserved)
! 
! Also specify the integer "number_of_inputs_nozlf" that split the 
! super_array into no-zlf and zlf
!  1: number_of_inputs_nozlf                    ---> No    ZLF
!  number_of_inputs_nozlf+1:moisture_array_size --> apply ZLF
!

counter = 0
number_of_inputs_zlf = counter 
IF ( i_cld_vn == i_cld_pc2 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)  &
!$OMP& SHARED(model_levels,pdims,super_array,counter,exner_theta_levels)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        super_array(i,j,k,counter+1) = exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF

!
! add the tracers to be conserved with ZLF and set the splitting
! identifier "number_of_inputs_nozlf = counter"
!

number_of_inputs_nozlf = counter
 
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array,counter,r_m_v,r_m_cl,r_m_cf)
DO k=0, model_levels
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      super_array(i,j,k,counter+1)  =  r_m_v(i,j,k)
      super_array(i,j,k,counter+2)  = r_m_cl(i,j,k)
      super_array(i,j,k,counter+3)  = r_m_cf(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
counter = counter + 3

IF ( i_cld_vn == i_cld_pc2 ) THEN
  IF (l_fixbug_pc2_mixph) THEN
    ! new method, advect the mixed phase cloud fraction
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array,counter,cf_liquid,  &
!$OMP& cfl_star,cf_frozen,cff_star,cf_bulk,cf_star)
    DO k=1, model_levels
      DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
          super_array(i,j,k,counter+1) = cf_liquid(i,j,k) +               &
                                         cfl_star (i,j,k) +               &
                                         cf_frozen(i,j,k) +               &
                                         cff_star (i,j,k) -               &
                                         cf_bulk(i,j,k) -                 &
                                         cf_star (i,j,k)
          super_array(i,j,k,counter+2) = cf_liquid(i,j,k) + cfl_star(i,j,k)
          super_array(i,j,k,counter+3) = cf_frozen(i,j,k) + cff_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ! original method, advect the bulk cloud fraction
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)  &
!$OMP& SHARED(model_levels,pdims,super_array,counter,cf_liquid,  &
!$OMP& cfl_star,cf_frozen,cff_star,cf_bulk,cf_star)
    DO k=1, model_levels
      DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
          super_array(i,j,k,counter+1) = cf_bulk(i,j,k)   + cf_star(i,j,k)
          super_array(i,j,k,counter+2) = cf_liquid(i,j,k) + cfl_star(i,j,k)
          super_array(i,j,k,counter+3) = cf_frozen(i,j,k) + cff_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      super_array(i,j,0,counter+1) = super_array(i,j,1,counter+1)
      super_array(i,j,0,counter+2) = super_array(i,j,1,counter+2)
      super_array(i,j,0,counter+3) = super_array(i,j,1,counter+3)
    END DO
  END DO
  counter = counter + 3  
END IF

IF ( L_mcr_rain ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array,counter,r_m_r)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        super_array(i,j,k,counter+1)   = r_m_r(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF

IF ( l_mcr_cf2 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array,counter,r_m_cf2)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        super_array(i,j,k,counter+1) = r_m_cf2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF

IF ( l_mcr_graup ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array,counter,r_m_gr)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        super_array(i,j,k,counter+1) = r_m_gr(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF

IF ( l_conv_prog_group_1 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array,counter,conv_prog_1)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        super_array(i,j,k,counter+1) = conv_prog_1(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF
IF ( l_conv_prog_group_2 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array,counter,conv_prog_2)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        super_array(i,j,k,counter+1) = conv_prog_2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF
IF ( l_conv_prog_group_3 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array,counter,conv_prog_3)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        super_array(i,j,k,counter+1) = conv_prog_3(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF
IF ( l_conv_prog_precip ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array,counter,conv_prog_precip)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        super_array(i,j,k,counter+1) = conv_prog_precip(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF


! DEPENDS ON: swap_bounds
CALL Swap_Bounds(                                                       &
             super_array, row_length, rows,                             &
             moisture_array_size*(model_levels+1),                      &
             halo_i, halo_j, fld_type_p, swap_field_is_scalar)

number_of_inputs = moisture_array_size


! Set layers over which linear interpolation is used
IF (l_sl_bc_correction) THEN
  k_int_linear=2
ELSE
  k_int_linear=1
END IF

IF ( zlf_conservation_moist_option > 0 ) THEN 

  number_of_inputs_zlf   = number_of_inputs - number_of_inputs_nozlf  

  CALL eg_zlf_conservation(super_array, super_array_out,            &
           number_of_inputs, number_of_inputs_zlf,                  &
           row_length, rows, n_rows, model_levels, halo_i, halo_j,  &
           offx, offy, datastart, g_i_pe, high_order_scheme,        &
           monotone_scheme,  l_high, l_mono,                        &
           zlf_conservation_moist_option, zlf_cfl_top_level_moist,  &
           rho_n, rho_np1, error_code                               )


ELSE ! use the old non conservative advection

  CALL eg_interpolation_eta_pmf(                                   &
                  eta_theta_levels,fld_type_w,                     &
                  number_of_inputs,                                &
                  row_length, rows, model_levels+1,                &
                  rows,                                            &
                  row_length, rows, model_levels+1,                &
                  high_order_scheme, monotone_scheme,              &
                  l_high, l_mono, depart_xi3_w, depart_xi1_w,      &
                  depart_xi2_w, mype, nproc, nproc_x, nproc_y,     &
                  halo_i, halo_j,                                  &
                  global_row_length, datastart, at_extremity,      &
                  g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                  offx, offy, error_code,                          &
                  super_array, super_array_out,                    &
                  k_int_linear_in=k_int_linear)


END IF ! endif for use ZLF

! Unpack super array

counter = 0
IF ( i_cld_vn == i_cld_pc2 ) THEN 
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)  &
!$OMP& SHARED(model_levels,pdims,super_array_out,counter,exner_star)
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        exner_star(i,j,k) = super_array_out(i,j,k,counter+1)
      END DO
    END DO
  END DO
  counter = counter + 1
END IF

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array_out,counter,       &
!$OMP& r_m_v_d,r_m_cl_d,r_m_cf_d)
DO k=0, model_levels
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      r_m_v_d(i,j,k)  = super_array_out(i,j,k,counter+1)
      r_m_cl_d(i,j,k) = super_array_out(i,j,k,counter+2)
      r_m_cf_d(i,j,k) = super_array_out(i,j,k,counter+3)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
counter = counter + 3

IF ( i_cld_vn == i_cld_pc2 ) THEN
  IF (l_fixbug_pc2_mixph) THEN
    ! new method, advect the mixed phase cloud fraction
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        cf_star(i,j,0) = super_array_out(i,j,1,counter+2)             &
                        +super_array_out(i,j,1,counter+3)             &
                        -super_array_out(i,j,1,counter+1)
        cfl_star(i,j,0)= super_array_out(i,j,1,counter+2)
        cff_star(i,j,0)= super_array_out(i,j,1,counter+3)
      END DO
    END DO
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)     &
!$OMP& SHARED(model_levels,pdims,super_array_out,counter,cf_star,   &
!$OMP& cfl_star,cff_star)
    DO k=1, model_levels
      DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
          cf_star(i,j,k) = super_array_out(i,j,k,counter+2)          &
                          +super_array_out(i,j,k,counter+3)          &
                          -super_array_out(i,j,k,counter+1)
          cfl_star(i,j,k)   = super_array_out(i,j,k,counter+2)
          cff_star(i,j,k)   = super_array_out(i,j,k,counter+3)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ! original method, advect the bulk cloud fraction
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        cf_star(i,j,0) = super_array_out(i,j,1,counter+1)
        cfl_star(i,j,0)= super_array_out(i,j,1,counter+2)
        cff_star(i,j,0)= super_array_out(i,j,1,counter+3)
      END DO
    END DO
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)     &
!$OMP& SHARED(model_levels,pdims,super_array_out,counter,cf_star,   &
!$OMP& cfl_star,cff_star)
    DO k=1, model_levels
      DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
          cf_star(i,j,k)    = super_array_out(i,j,k,counter+1)
          cfl_star(i,j,k)   = super_array_out(i,j,k,counter+2)
          cff_star(i,j,k)   = super_array_out(i,j,k,counter+3)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF !l_fixbug_pc2_mixph
  counter = counter + 3
END IF

IF ( L_mcr_rain ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array_out,counter,r_m_r_d)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        r_m_r_d(i,j,k) = super_array_out(i,j,k,counter+1)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF

IF ( l_mcr_cf2 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array_out,counter,r_m_cf2_d)   
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        r_m_cf2_d(i,j,k) = super_array_out(i,j,k,counter+1)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF

IF ( l_mcr_graup ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array_out,counter,r_m_gr_d)   
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        r_m_gr_d(i,j,k) = super_array_out(i,j,k,counter+1)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF

IF ( l_conv_prog_group_1 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array_out,counter,conv_prog_1)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        conv_prog_1(i,j,k) = super_array_out(i,j,k,counter+1)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF
IF ( l_conv_prog_group_2 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array_out,counter,conv_prog_2)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        conv_prog_2(i,j,k) = super_array_out(i,j,k,counter+1)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF
IF ( l_conv_prog_group_3 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array_out,counter,conv_prog_3)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        conv_prog_3(i,j,k) = super_array_out(i,j,k,counter+1)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF
IF ( l_conv_prog_precip ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,pdims,super_array_out,counter,conv_prog_precip)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        conv_prog_precip(i,j,k) = super_array_out(i,j,k,counter+1)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  counter = counter + 1
END IF


IF ( L_print_L2norms ) THEN
  WRITE(umMessage,'(A, I2, A)') ' ** cycleno =  ', cycleno              &
                           , ' **   L2 norms after eg_sl_moisture  **'
  CALL umPrint(umMessage,src='EG_SL_MOISTURE')
  CALL sl_moist_norm(                                                   &
                     norm_lev_start, norm_lev_end, n_rims_to_do,        &
                     exner_star, r_m_v_d, r_m_cl_d, r_m_cf_d,           &
                     r_m_r_d, r_m_gr_d, r_m_cf2_d,                      &
                     cf_star, cfl_star, cff_star,                       &
                    .TRUE., .FALSE., .FALSE. )
END IF !  L_print_L2norms

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_sl_moisture
END MODULE eg_sl_moisture_mod

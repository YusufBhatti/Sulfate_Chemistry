! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_moisture_pseudo_lbflux_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='EG_MOISTURE_PSEUDO_LBFLUX_MOD'

CONTAINS
SUBROUTINE eg_moisture_pseudo_lbflux(                                 &
              moisture_array_size,                                    &
              row_length, rows, n_rows, model_levels, halo_i,         &
              halo_j, offx, offy,datastart,g_i_pe,high_order_scheme,  &
              monotone_scheme,  l_high, l_mono, L_mcr_rain,           &
              l_mcr_cf2,l_mcr_graup,                                  &
              r_m_v, r_m_cl, r_m_cf, r_m_r, r_m_gr, r_m_cf2,          &
              cf_bulk, cf_liquid, cf_frozen,exner_theta_levels,       &
              cf_star, cfl_star, cff_star,                            &
              rho_n, rho_np1, l_conserv,                              &
              pseudo_lbflux_moisture, number_qs,                      &
              error_code )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE um_parcore,      ONLY: mype, nproc
USE um_parvars,      ONLY: nproc_x,nproc_y,                           &
                              at_extremity,gc_proc_row_group,         &
                              gc_proc_col_group

USE nlsizes_namelist_mod,  ONLY: global_row_length                       

USE timestep_mod,      ONLY: timestep
USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels,       &
                              xi3_at_theta=>r_theta_levels,           &
                              xi3_at_rho=>r_rho_levels,               &
                              xi3_at_u=>r_at_u, xi3_at_v=>r_at_v

USE cloud_inputs_mod, ONLY: l_fixbug_pc2_mixph, i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2
USE atm_fields_bounds_mod, ONLY : pdims, tdims_s, tdims_l
USE eg_interpolation_eta_pseudo_lbflux_mod, ONLY: &
                                     eg_interpolation_eta_pseudo_lbflux
USE horiz_grid_mod
USE metric_terms_mod
USE departure_pts_mod
USE mpp_conf_mod,     ONLY: swap_field_is_scalar
USE ereport_mod, ONLY: ereport
USE Field_Types

IMPLICIT NONE
!
! Description:
!   Compute pseudo lateral boundary flux (PLF) of moisture.
!
! Method: ENDGame formulation version 4.xx
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dynamics advection
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.


! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_MOISTURE_PSEUDO_LBFLUX'



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

LOGICAL, INTENT(IN) ::                                                &
  l_conserv
                   ! True, if LAM-conservation is required.

LOGICAL, INTENT(IN) :: l_mcr_rain, l_mcr_cf2, l_mcr_graup


! Timelevel n arrival point quantities

REAL, INTENT(IN) ::                                                   &
  r_m_v(row_length,rows,0:model_levels),                              &
  r_m_cl(row_length,rows,0:model_levels),                             &
  r_m_cf(row_length,rows,0:model_levels),                             &
  r_m_r(row_length,rows,0:model_levels),                              &
  r_m_gr(row_length,rows,0:model_levels),                             &
  r_m_cf2(row_length,rows,0:model_levels)

REAL, INTENT(IN) ::                                                   &
    rho_n(1-offx:row_length+offx,1-offy:rows+offy,1:model_levels),    &
  rho_np1(1-offx:row_length+offx,1-offy:rows+offy,1:model_levels)

REAL, INTENT (IN) ::                                               &
        cf_star (row_length, rows, 0:model_levels),                &
        cfl_star(row_length, rows, 0:model_levels),                &
        cff_star(row_length, rows, 0:model_levels)

REAL, INTENT (INOUT) :: cf_bulk  (tdims_l%i_start:tdims_l%i_end,      &
                        tdims_l%j_start:tdims_l%j_end,                &
                        tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: cf_liquid(tdims_l%i_start:tdims_l%i_end,      &
                        tdims_l%j_start:tdims_l%j_end,                &
                        tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: cf_frozen(tdims_l%i_start:tdims_l%i_end,      &
                        tdims_l%j_start:tdims_l%j_end,                &
                        tdims_l%k_start:tdims_l%k_end)

REAL :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end)


INTEGER :: error_code   ! Non-zero on exit if error detected.

! Pseudo lateral boundary flux
INTEGER, INTENT(IN) :: number_qs
REAL, INTENT(OUT) :: pseudo_lbflux_moisture(number_qs)

! Local variables

INTEGER :: i,j,k, number_of_inputs

! tmp & dummy arrays


REAL :: super_array(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,      &
                0:model_levels, moisture_array_size )

REAL :: super_array_out(1-offx:row_length+offx,1-offy:rows+offy,          &
                0:model_levels, moisture_array_size)


INTEGER :: nqs


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(NONE)                 &
!$OMP& PRIVATE(k, j, i)                         &
!$OMP& SHARED(model_levels, pdims,              &
!$OMP&        super_array, r_m_v, r_m_cl, r_m_cf)
DO k=0, model_levels
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      super_array(i,j,k,1)  = r_m_v(i,j,k)
      super_array(i,j,k,2)  = r_m_cl(i,j,k)
      super_array(i,j,k,3)  = r_m_cf(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

nqs = 3
IF ( i_cld_vn == i_cld_pc2 ) THEN
  IF (l_fixbug_pc2_mixph) THEN
    ! new method, advect the mixed phase cloud fraction
!$OMP PARALLEL DO DEFAULT(NONE)                                 &
!$OMP& PRIVATE(k, j, i)                                         &
!$OMP& SHARED(model_levels, pdims,                              &
!$OMP&        nqs, super_array, cf_liquid, cfl_star, cf_frozen, &
!$OMP&        cff_star, cf_bulk, cf_star, exner_theta_levels)
    DO k=1, model_levels
      DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
          super_array(i,j,k,nqs+1) = cf_liquid(i,j,k) +                &
                                       cfl_star (i,j,k) +                &
                                       cf_frozen(i,j,k) +                &
                                       cff_star (i,j,k) -                &
                                       cf_bulk(i,j,k) -                  &
                                       cf_star (i,j,k)
          super_array(i,j,k,nqs+2) = cf_liquid(i,j,k) + cfl_star(i,j,k)
          super_array(i,j,k,nqs+3) = cf_frozen(i,j,k) + cff_star(i,j,k)
          super_array(i,j,k,nqs+4) = exner_theta_levels(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ! original method, advect the bulk cloud fraction
!$OMP PARALLEL DO DEFAULT(NONE)                                &
!$OMP& PRIVATE(k, j, i)                                        &
!$OMP& SHARED(model_levels, pdims,                             &
!$OMP&        nqs, super_array, cf_bulk, cf_star, cf_liquid,   &
!$OMP&        cfl_star, cf_frozen, cff_star, exner_theta_levels)
    DO k=1, model_levels
      DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end
          super_array(i,j,k,nqs+1) = cf_bulk(i,j,k)   + cf_star(i,j,k)
          super_array(i,j,k,nqs+2) = cf_liquid(i,j,k) + cfl_star(i,j,k)
          super_array(i,j,k,nqs+3) = cf_frozen(i,j,k) + cff_star(i,j,k)
          super_array(i,j,k,nqs+4) = exner_theta_levels(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

!$OMP PARALLEL DO DEFAULT(NONE)              &
!$OMP& PRIVATE(j, i)                         &
!$OMP& SHARED(pdims, nqs,                    &
!$OMP&        super_array, exner_theta_levels)
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      super_array(i,j,0,nqs+1) = super_array(i,j,1,nqs+1)
      super_array(i,j,0,nqs+2) = super_array(i,j,1,nqs+2)
      super_array(i,j,0,nqs+3) = super_array(i,j,1,nqs+3)
      super_array(i,j,0,nqs+4) = exner_theta_levels(i,j,0)
    END DO
  END DO
!$OMP END PARALLEL DO
  nqs = nqs + 4
END IF

IF ( L_mcr_rain ) THEN
!$OMP PARALLEL DO DEFAULT(NONE)      &
!$OMP& PRIVATE(k, j, i)              &
!$OMP& SHARED(model_levels, pdims,   &
!$OMP&        nqs, super_array, r_m_r)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        super_array(i,j,k,nqs+1)   = r_m_r(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  nqs = nqs + 1
END IF

IF ( l_mcr_cf2 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE)        &
!$OMP& PRIVATE(k, j, i)                &
!$OMP& SHARED(model_levels, pdims,     &
!$OMP&        nqs, super_array, r_m_cf2)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        super_array(i,j,k,nqs+1) = r_m_cf2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  nqs = nqs + 1
END IF

IF ( l_mcr_graup ) THEN
!$OMP PARALLEL DO DEFAULT(NONE)       &
!$OMP& PRIVATE(k, j, i)               &
!$OMP& SHARED(model_levels, pdims,    &
!$OMP&        nqs, super_array, r_m_gr)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        super_array(i,j,k,nqs+1) = r_m_gr(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  nqs = nqs + 1
END IF


! DEPENDS ON: swap_bounds
CALL Swap_Bounds(                                               &
     super_array, row_length, rows,                             &
     moisture_array_size*(model_levels+1),                      &
     halo_i, halo_j, fld_type_p, swap_field_is_scalar )

number_of_inputs = moisture_array_size

pseudo_lbflux_moisture(:) = 0.0
CALL eg_interpolation_eta_pseudo_lbflux(                            &
                   eta_theta_levels,fld_type_w,                     &
                   number_of_inputs,                                &
                   row_length, rows, model_levels+1, rows,          &
                   row_length, rows, model_levels+1,                &
                   high_order_scheme, monotone_scheme,              &
                   l_high, l_mono, depart_xi3_w, depart_xi1_w,      &
                   depart_xi2_w, mype, nproc, nproc_x, nproc_y,     &
                   halo_i, halo_j,                                  &
                   global_row_length, datastart, at_extremity,      &
                   g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                   offx, offy, offx ,offy, error_code,              &
                   super_array, super_array_out,                    &
                   rho_n, rho_np1, l_conserv,                       &
                   pseudo_lbflux_moisture)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_moisture_pseudo_lbflux
END MODULE eg_moisture_pseudo_lbflux_mod

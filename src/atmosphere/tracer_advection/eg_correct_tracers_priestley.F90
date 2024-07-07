! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE eg_correct_tracers_priestley_mod
IMPLICIT NONE

! Description:
!
!            This routine enforces the mass conservation of tracers,
!            using a modified version of the Priestley(1993) scheme.
!            It should work with either mixing ratios and dry density
!            or with specific humdities and wet density.
!            If L_conserve_tracers=.FALSE., then the routine make
!            sure that the tracers are strictly positive only and
!            doesn't impose any mass conservation constraint.
!
! Method: ENDGame formulation version 4.00 (Feb 2014)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Tracer Advection
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='EG_CORRECT_TRACERS_PRIESTLEY_MOD'

CONTAINS
SUBROUTINE eg_correct_tracers_priestley(                          &
                         super_array_size,                        &
                         super_tracer_phys1, super_tracer_phys2,  &
                         rho_n, rho_np1,                          &
                         co2, L_CO2_interactive,                  &
                         murk, L_murk_advect,                     &
                         soot_new, soot_agd, soot_cld, L_soot,    &
                         bmass_new, bmass_agd, bmass_cld,         &
                         L_biomass,                               &
                         ocff_new, ocff_agd, ocff_cld, l_ocff,    &
                         dust_div1,dust_div2,dust_div3,           &
                         dust_div4,dust_div5,dust_div6,           &
                         l_dust,                                  &
                         so2, so4_aitken, so4_accu,               &
                         so4_diss, nh3, dms,                      &
                         L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms,   &
                         nitr_acc, nitr_diss, L_nitrate,          &
                         l_use_cariolle, ozone_tracer,            &
                         tracers, tr_vars,                        &
                         tr_ukca, tracer_ukca,                    &
                         L_conserve_tracers,                      &
                         g_i_pe,depart_lambda,depart_phi,depart_r )


USE um_parvars,     ONLY: nproc_x, nproc_y,                       &
                          at_extremity,gc_proc_row_group,         &
                          gc_proc_col_group, datastart
USE um_parcore,     ONLY: mype, nproc

USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows


USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels
USE eg_group_tracers_mod
USE eg_ungroup_tracers_mod
USE dust_parameters_mod, ONLY: L_twobin_dust
USE atm_fields_bounds_mod
USE Field_Types
USE priestley_algorithm_mod
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE sl_input_mod, ONLY: tr_priestley_opt, L_tr_src_in_conserve
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE eg_helmholtz_mod,      ONLY: ec_vol
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE
INTEGER, INTENT(IN) :: g_i_pe(1-(tdims%i_start-tdims_l%i_start):  &
              global_row_length+(tdims%i_start-tdims_l%i_start) )

REAL, INTENT(IN)    :: depart_lambda( wdims%i_start:wdims%i_end,  &
                                      wdims%j_start:wdims%j_end,  &
                                      wdims%k_start:wdims%k_end), &
!
                             depart_phi( wdims%i_start:wdims%i_end,     &
                                         wdims%j_start:wdims%j_end,     &
                                         wdims%k_start:wdims%k_end),    &
!
                             depart_r( wdims%i_start:wdims%i_end,       &
                                       wdims%j_start:wdims%j_end,       &
                                       wdims%k_start:wdims%k_end)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CORRECT_TRACERS_PRIESTLEY'

INTEGER, INTENT(IN)  :: super_array_size
INTEGER, INTENT(IN)  :: tr_ukca, tr_vars
LOGICAL, INTENT(IN)  :: L_CO2_interactive
LOGICAL, INTENT(IN)  :: L_murk_advect, L_Soot
LOGICAL, INTENT(IN)  :: L_biomass, L_ocff, l_dust, L_sulpc_so2
LOGICAL, INTENT(IN)  :: L_sulpc_nh3, l_sulpc_dms, l_use_cariolle
LOGICAL, INTENT(IN)  :: L_nitrate, L_conserve_tracers

REAL, INTENT(INOUT) ::  co2      (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  murk     (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  soot_new (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  soot_agd (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  soot_cld (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  so2      (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  so4_aitken                              &
                                 (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  so4_accu (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  so4_diss (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  nh3      (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  dms      (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  dust_div1(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  dust_div2(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  dust_div3(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  dust_div4(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  dust_div5(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  dust_div6(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  bmass_new(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  bmass_agd(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  bmass_cld(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  ocff_new (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  ocff_agd (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  ocff_cld (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  nitr_acc (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  nitr_diss(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  ozone_tracer                            &
                                 (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) :: tracer_ukca                              &
                                 (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end,&
                                  tr_ukca)
REAL, INTENT(INOUT) :: tracers                                  &
                                 (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end,&
                                  tr_vars)

! Note super_tracer_phys1 should have up-to-date halos
REAL, INTENT(IN) :: super_tracer_phys1                          &
               (tdims_l%i_start:tdims_l%i_end,                  &
                tdims_l%j_start:tdims_l%j_end,                  &
                tdims_l%k_start:tdims_l%k_end,                  &
                super_array_size)

REAL, INTENT(IN) :: super_tracer_phys2                          &
               (tdims%i_start:tdims%i_end,                      &
                tdims%j_start:tdims%j_end,                      &
                tdims%k_start:tdims%k_end,                      &
                super_array_size)

REAL, INTENT(IN)   :: rho_n                                     &
               (pdims_s%i_start:pdims_s%i_end,                  &
                pdims_s%j_start:pdims_s%j_end,                  &
                pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(IN)   :: rho_np1                                   &
               (pdims_s%i_start:pdims_s%i_end,                  &
                pdims_s%j_start:pdims_s%j_end,                  &
                pdims_s%k_start:pdims_s%k_end)

! Local Variables

LOGICAL :: l_high_in   ! local setting for use in CALLs
LOGICAL :: l_mono_in   ! local setting for use in CALLs

REAL ::    super_array (tdims_s%i_start:tdims_s%i_end,          &
                        tdims_s%j_start:tdims_s%j_end,          &
                        tdims_s%k_start:tdims_s%k_end,          &
                        super_array_size)

REAL :: alfa_za(pdims%k_start:pdims%k_end)
REAL :: beta_za(pdims%k_start:pdims%k_end)

REAL :: psi_n(tdims%i_start:tdims%i_end,                      &
              tdims%j_start:tdims%j_end,                      &
              tdims%k_start:tdims%k_end)
REAL :: psi_np1(tdims%i_start:tdims%i_end,                    &
              tdims%j_start:tdims%j_end,                      &
              tdims%k_start:tdims%k_end)

INTEGER :: tracers_switches(super_array_size)
INTEGER :: ist,iend, number_qs, total_number_tracers
INTEGER :: i,j,k,kk, mpierr
INTEGER :: i_halo_big, j_halo_big, i_halo_small, j_halo_small

REAL,    ALLOCATABLE :: qs_n(:,:,:,:)
REAL,    ALLOCATABLE :: qs_n2(:,:,:,:)
REAL,    ALLOCATABLE :: qs_np1(:,:,:,:)
REAL,    ALLOCATABLE :: qs_np1_low(:,:,:,:)
REAL,    ALLOCATABLE :: qs_s(:,:,:,:)
REAL,    ALLOCATABLE :: qs_s2(:,:,:,:)
REAL,    ALLOCATABLE :: qsmin(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!================================================================
! pack/group all the tracers sources into a superarray
!================================================================

CALL eg_group_tracers(                                          &
                    super_array_size,                           &
                    total_number_tracers,                       &
                    super_array,                                &
                    L_CO2_interactive, co2,                     &
                    L_Murk_advect, murk,                        &
                    L_Soot, soot_new, soot_agd, soot_cld,       &
                    l_sulpc_so2, so2, SO4_aitken, so4_accu,     &
                                 so4_diss,                      &
                    L_sulpc_nh3, nh3,                           &
                    L_sulpc_dms, dms,                           &
                    L_dust, dust_div1, dust_div2, dust_div3,    &
                            dust_div4, dust_div5, dust_div6,    &
                    L_biomass, bmass_new, bmass_agd, bmass_cld, &
                    L_ocff, ocff_new, ocff_agd, ocff_cld,       &
                    L_nitrate, nitr_acc, nitr_diss,             &
                    l_use_cariolle, ozone_tracer,               &
                    tracers, tr_vars,                           &
                    tr_ukca, tracer_ukca,                       &
                    L_twobin_dust, tracers_switches             )


IF ( L_conserve_tracers ) THEN

  i_halo_big   = tdims%i_start-tdims_l%i_start
  i_halo_small = tdims%i_start-tdims_s%i_start
  j_halo_big   = tdims%j_start-tdims_l%j_start
  j_halo_small = tdims%j_start-tdims_s%j_start

  number_qs =  total_number_tracers

  ALLOCATE(                                                     &
         qs_n(tdims_s%i_start:tdims_s%i_end,                    &
              tdims_s%j_start:tdims_s%j_end,                    &
              tdims_s%k_start:tdims_s%k_end,number_qs))
  ALLOCATE(                                                     &
         qs_n2(tdims_l%i_start:tdims_l%i_end,                   &
               tdims_l%j_start:tdims_l%j_end,                   &
               tdims_l%k_start:tdims_l%k_end,number_qs))
  ALLOCATE(                                                     &
       qs_np1(tdims_s%i_start:tdims_s%i_end,                    &
              tdims_s%j_start:tdims_s%j_end,                    &
              tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(                                                     &
       qs_np1_low(tdims_s%i_start:tdims_s%i_end,                &
                  tdims_s%j_start:tdims_s%j_end,                &
                  tdims_s%k_start:tdims_s%k_end,number_qs))
  ALLOCATE(                                                     &
         qs_s(tdims_s%i_start:tdims_s%i_end,                    &
              tdims_s%j_start:tdims_s%j_end,                    &
              tdims_s%k_start:tdims_s%k_end,number_qs))
  ALLOCATE(                                                     &
         qs_s2(tdims_s%i_start:tdims_s%i_end,                   &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(                                                    &
      qsmin(number_qs))


  qsmin = 0.0
  ist   = 1
  iend  = total_number_tracers

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP PRIVATE(kk, k, i, j)                                             &
!$OMP SHARED(tdims_s, tdims_l, rho_n, ec_vol, alfa_za, beta_za, psi_n, &
!$OMP        rho_np1, qs_np1, super_array, qs_n2, super_tracer_phys1,  &
!$OMP        qs_n, qs_s, number_qs, ist, iend, pdims, psi_np1, tdims,  &
!$OMP        eta_theta_levels, super_tracer_phys2, eta_rho_levels)

!$OMP DO SCHEDULE(STATIC)
  DO kk = 1, number_qs
    DO k = tdims%k_start, tdims%k_end    
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end    
          qs_np1(i,j,k,kk) = super_array(i,j,k,kk)
        END DO
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  !======================================================================
  ! Put fields into one superarray conistent with the field indexing
  !======================================================================

!$OMP DO SCHEDULE(STATIC)
  DO kk = ist, iend
    DO k = tdims_l%k_start, tdims_l%k_end    
      DO j = tdims_l%j_start, tdims_l%j_end
        DO i = tdims_l%i_start, tdims_l%i_end    
            qs_n2(i,j,k,kk) = super_tracer_phys1(i,j,k,kk)
        END DO
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO kk = ist, iend
    DO k = tdims%k_start, tdims%k_end    
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end    
            qs_n(i,j,k,kk)  = super_tracer_phys1(i,j,k,kk)
        END DO
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO kk = ist, iend
    DO k = tdims%k_start, tdims%k_end    
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end    
            qs_s(i,j,k,kk)  = super_tracer_phys2(i,j,k,kk)
        END DO
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  ! Enforce the boundary condition field(surface)=field(level 1)

!$OMP DO SCHEDULE(STATIC)
  DO kk = 1, number_qs   
    DO j = tdims_l%j_start, tdims_l%j_end
      DO i = tdims_l%i_start, tdims_l%i_end   
        qs_n2(i,j,tdims_l%k_start,kk) =  qs_n2(i,j,tdims_l%k_start+1,kk)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO kk = 1, number_qs   
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end   
        qs_n(i,j,tdims_s%k_start,kk)  =  qs_n(i,j,tdims_s%k_start+1,kk)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  !=========================================================================
  ! Compute the vertical linear averaging weights used to compute tracer-mass
  !=========================================================================

  ! Set values for bottom-most level, assuming tr(:,:,0) = tr(:,:,1)
!$OMP SINGLE
  alfa_za(pdims%k_start) = 1.0
  beta_za(pdims%k_start) = 0.0
!$OMP END SINGLE NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = pdims%k_start + 1, pdims%k_end
    alfa_za(k) = (   eta_rho_levels(k) - eta_theta_levels(k-1) ) /    &
                ( eta_theta_levels(k) - eta_theta_levels(k-1) )
    beta_za(k) = 1.0 - alfa_za(k)
  END DO
!$OMP END DO NOWAIT

  ! Explicitly define psi_n and psi_n1 arrays. 

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end              
      psi_n  (i,j,tdims%k_start) = 0.0
      psi_np1(i,j,tdims%k_start) = 0.0
    END DO
  END DO
!$OMP END DO 
  !=========================================================================
  ! Compute the 3D array psi(i,j,k) = rho*volume*average
  !   so mass_tracer = SUM[psi*tracer] = SUM[av(tracer)*rho*vol]
  !=========================================================================

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start + 1, tdims%k_end - 1

    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        psi_n(i,j,k)  =  rho_n(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )    &
                      +  rho_n(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)

        psi_np1(i,j,k) = rho_np1(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )  &
                       + rho_np1(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)
      END DO
    END DO
  END DO ! k
!$OMP END DO NOWAIT

  k = tdims%k_end

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      psi_n(i,j,k)   =   rho_n(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
      psi_np1(i,j,k) = rho_np1(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !=====================================================================
  ! Call 3-tri-linear interpolation to compute qs_np1_low
  !=====================================================================

  l_high_in=.FALSE.
  l_mono_in=.TRUE.
  CALL eg_interpolation_eta_pmf(                               &
                eta_theta_levels,fld_type_w, number_qs,        &
                tdims%i_len,                                   &
                tdims%j_len,                                   &
                tdims%k_len,                                   &
                tdims%j_len,                                   &
                tdims%i_len,                                   &
                tdims%j_len,                                   &
                tdims%k_len,0, 1,                              &
                l_high_in,l_mono_in,                           &
                depart_r, depart_lambda, depart_phi,           &
                mype, nproc, nproc_x, nproc_y,                 &
                i_halo_big, j_halo_big,                        &
                global_row_length, datastart, at_extremity,    &
                g_i_pe, gc_proc_row_group, gc_proc_col_group,  &
                i_halo_small, j_halo_small,                    &
                mpierr,                                        &
                qs_n2, qs_np1_low                              )

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP PRIVATE(kk, k, j, i)                                             &
!$OMP SHARED(number_qs, tdims, qs_np1, qs_s2, qs_s, qs_np1_low,        &
!$OMP        L_tr_src_in_conserve) 

  IF ( .NOT. L_tr_src_in_conserve ) THEN  
    DO kk = 1, number_qs
      DO k = tdims%k_start, tdims%k_end
!$OMP DO SCHEDULE(STATIC)  
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            qs_np1(i,j,k,kk) = qs_np1(i,j,k,kk) - qs_s(i,j,k,kk)
            qs_s2(i,j,k,kk) = 0.0
          END DO
        END DO
!$OMP END DO NOWAIT
      END DO
    END DO
  ELSE

    DO kk = 1, number_qs
      DO k = tdims%k_start, tdims%k_end
!$OMP DO SCHEDULE(STATIC)   
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            qs_np1_low(i,j,k,kk) = qs_np1_low(i,j,k,kk) + qs_s(i,j,k,kk)
            qs_s2(i,j,k,kk) =  qs_s(i,j,k,kk)
          END DO
        END DO
!$OMP END DO NOWAIT
      END DO ! k
    END DO
  END IF

    DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qs_np1_low(i,j,tdims%k_start,kk) = qs_np1_low(i,j,tdims%k_start+1,kk)
          qs_np1(i,j,tdims%k_start,kk) =     qs_np1(i,j,tdims%k_start+1,kk)
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
!$OMP END PARALLEL

  !=====================================================================
  ! Apply conservation to qs_np1
  !=====================================================================

  IF ( tr_priestley_opt == 1 ) THEN

    CALL priestley_algorithm( qs_n, qs_np1, qs_np1_low, qs_s2,          &
                              psi_n, psi_np1,                           &
                              tdims_s%i_start, tdims_s%i_end,           &
                              tdims_s%j_start, tdims_s%j_end,           &
                              tdims_s%k_start, tdims_s%k_end,           &
                              tdims%i_start  , tdims%i_end,             &
                              tdims%j_start  , tdims%j_end,             &
                              tdims%k_start  , tdims%k_end,             &
                              number_qs                                 )
  ELSE

    CALL priestley_algorithm2(qs_n, qs_np1, qs_np1_low, qs_s2,         &
                              psi_n, psi_np1,                          &
                              tdims_s%i_start, tdims_s%i_end,          &
                              tdims_s%j_start, tdims_s%j_end,          &
                              tdims_s%k_start, tdims_s%k_end,          &
                              tdims%i_start  , tdims%i_end,            &
                              tdims%j_start  , tdims%j_end,            &
                              tdims%k_start  , tdims%k_end,            &
                              number_qs                                )

  END IF

  IF ( .NOT. L_tr_src_in_conserve ) THEN  
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP PRIVATE(kk, k, j, i)                                             &
!$OMP SHARED(number_qs, tdims, qs_np1, qs_s)
    DO kk = 1, number_qs
      DO k = tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            qs_np1(i,j,k,kk) = qs_np1(i,j,k,kk) + qs_s(i,j,k,kk)
          END DO
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  !=========================================================================
  ! Impose minimum values on the corrected qs while maintaining conservation
  !=========================================================================

  CALL impose_minima_with_conservation(qs_np1, qsmin, psi_np1,          &
                                       tdims_s%i_start, tdims_s%i_end,  &
                                       tdims_s%j_start, tdims_s%j_end,  &
                                       tdims_s%k_start, tdims_s%k_end,  &
                                       tdims%i_start  , tdims%i_end,    &
                                       tdims%j_start  , tdims%j_end,    &
                                       tdims%k_start  , tdims%k_end,    &
                                       number_qs                        )

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP PRIVATE(kk, k, j, i)                                             &
!$OMP SHARED(number_qs, tdims, super_array, qs_np1)
  DO kk = 1, number_qs
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          super_array(i,j,k,kk) = qs_np1(i,j,k,kk)
        END DO
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO


  DEALLOCATE( qsmin      )
  DEALLOCATE( qs_s       )
  DEALLOCATE( qs_s2      )
  DEALLOCATE( qs_np1     )
  DEALLOCATE( qs_np1_low )
  DEALLOCATE( qs_n       )
  DEALLOCATE( qs_n2      )

ELSE

  !=========================================================================
  ! if L_conserve_tracers==.false. then simply make sure that the tracers
  ! are strictly positive
  !=========================================================================

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) COLLAPSE(2)           &
!$OMP PRIVATE(kk, k, j, i)                                             &
!$OMP SHARED(tdims_s, super_array, super_array_size)
  DO kk = 1, super_array_size
    DO k = tdims_s%k_start, tdims_s%k_end
      DO j = tdims_s%j_start, tdims_s%j_end
        DO i = tdims_s%i_start, tdims_s%i_end
          super_array(i,j,k,kk) = MAX(super_array(i,j,k,kk), 0.0 )
        END DO
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! L_conserve_tracers

! unpack/ungroup the superarray. Copy the superarray fields into their
! appropriate arrays (i.e.,reverse the operation in "eg_group_tracers")

CALL eg_ungroup_tracers(                                        &
                    super_array_size,                           &
                    super_array,                                &
                    L_CO2_interactive, co2,                     &
                    L_Murk_advect, murk,                        &
                    L_Soot, soot_new, soot_agd, soot_cld,       &
                    l_sulpc_so2, so2, SO4_aitken, so4_accu,     &
                                 so4_diss,                      &
                    L_sulpc_nh3, nh3,                           &
                    L_sulpc_dms, dms,                           &
                    L_dust, dust_div1, dust_div2, dust_div3,    &
                            dust_div4, dust_div5, dust_div6,    &
                    L_biomass, bmass_new, bmass_agd, bmass_cld, &
                    L_ocff, ocff_new, ocff_agd, ocff_cld,       &
                    L_nitrate, nitr_acc, nitr_diss,             &
                    l_use_cariolle, ozone_tracer,               &
                    tracers, tr_vars,                           &
                    tr_ukca, tracer_ukca,                       &
                    L_twobin_dust                               )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_correct_tracers_priestley

END MODULE eg_correct_tracers_priestley_mod

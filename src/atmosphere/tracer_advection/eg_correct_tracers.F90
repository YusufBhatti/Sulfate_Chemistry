! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE eg_conserv_tracers_mod
IMPLICIT NONE

! Description:
!
!            This routine enforces the mass conservation of tracers.
!            It should work with either mixing ratios and dry density
!            or with specific humdities and wet density.
!            If L_conserve_tracers=.FALSE., then the routine make
!            sure that the tracers are strictly positive only and
!            doesn't impose any mass conservation constraint.
!
! Method: ENDGame formulation version 3.02
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Tracer Advection
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_CONSERV_TRACERS_MOD'

CONTAINS
SUBROUTINE eg_correct_tracers_fix(                              &
                         mype, super_array_size,                &
                         super_tracer_phys1, super_tracer_phys2,&
                         rho_n, rho_np1,                        &
                         co2, L_co2_interactive,                &
                         murk, L_murk_advect,                   &
                         soot_new, soot_agd, soot_cld, L_soot,  &
                         bmass_new, bmass_agd, bmass_cld,       &
                         L_biomass,                             &
                         ocff_new, ocff_agd, ocff_cld, l_ocff,  &
                         dust_div1,dust_div2,dust_div3,         &
                         dust_div4,dust_div5,dust_div6,         &
                         L_dust,                                &
                         so2, so4_aitken, so4_accu,             &
                         so4_diss, nh3, dms,                    &
                         L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms, &
                         nitr_acc, nitr_diss, L_nitrate,        &
                         L_use_cariolle, ozone_tracer,          &
                         tracers, tr_vars,                      &
                         tr_ukca, tracer_ukca,                  &
                         L_conserve_tracers,                    &
                         L_conserv_smooth_lap                   )

USE eg_group_tracers_mod
USE eg_ungroup_tracers_mod
USE eg_check_conserv_mod
USE eg_mass_conserv_mod
USE dust_parameters_mod, ONLY: L_twobin_dust
USE atm_fields_bounds_mod
USE umPrintMgr
USE Field_Types
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CORRECT_TRACERS_FIX'

INTEGER, INTENT(IN)  :: mype
INTEGER, INTENT(IN)  :: super_array_size
INTEGER, INTENT(IN)  :: tr_ukca, tr_vars
LOGICAL, INTENT(IN)  :: L_CO2_interactive
LOGICAL, INTENT(IN)  :: L_murk_advect, L_Soot
LOGICAL, INTENT(IN)  :: L_biomass, L_ocff, l_dust, L_sulpc_so2
LOGICAL, INTENT(IN)  :: L_sulpc_nh3, l_sulpc_dms, l_use_cariolle
LOGICAL, INTENT(IN)  :: L_nitrate, L_conserve_tracers
LOGICAL, INTENT(IN)  :: L_conserv_smooth_lap

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

REAL ::    super_array (tdims_s%i_start:tdims_s%i_end,          &
                        tdims_s%j_start:tdims_s%j_end,          &
                        tdims_s%k_start:tdims_s%k_end,          &
                        super_array_size)

INTEGER :: ist,iend, number_qs, total_number_tracers
INTEGER :: i, mpierr

REAL,    ALLOCATABLE :: qs_n(:,:,:,:)
REAL,    ALLOCATABLE :: qs_np1(:,:,:,:)
REAL,    ALLOCATABLE :: qs_s(:,:,:,:)
REAL,    ALLOCATABLE :: qsmin(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! pack/group all the tracers sources into a superarray

CALL eg_group_tracers_fix(                                      &
                    super_array_size,                           &
                    total_number_tracers,                       &
                    super_array,                                &
                    L_CO2_interactive, co2,                     &
                    L_Murk_advect, murk,                        &
                    L_Soot, soot_new, soot_agd, soot_cld,       &
                    L_sulpc_so2, so2, so4_aitken, so4_accu,     &
                                 so4_diss,                      &
                    L_sulpc_nh3, nh3,                           &
                    L_sulpc_dms, dms,                           &
                    L_dust, dust_div1, dust_div2, dust_div3,    &
                            dust_div4, dust_div5, dust_div6,    &
                    L_biomass, bmass_new, bmass_agd, bmass_cld, &
                    L_ocff, ocff_new, ocff_agd, ocff_cld,       &
                    L_nitrate, nitr_acc, nitr_diss,             &
                    L_use_cariolle, ozone_tracer,               &
                    tracers, tr_vars,                           &
                    tr_ukca, tracer_ukca,                       &
                    L_twobin_dust                               )


IF (L_conserve_tracers ) THEN

  number_qs =  total_number_tracers

  ALLOCATE(                                                     &
         qs_n(tdims_s%i_start:tdims_s%i_end,                    &
              tdims_s%j_start:tdims_s%j_end,                    &
              tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(                                                     &
       qs_np1(tdims_s%i_start:tdims_s%i_end,                    &
              tdims_s%j_start:tdims_s%j_end,                    &
              tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(                                                     &
         qs_s(tdims_s%i_start:tdims_s%i_end,                    &
              tdims_s%j_start:tdims_s%j_end,                    &
              tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(                                                    &
      qsmin(number_qs))


  qsmin     = 0.0
  qs_np1    = super_array(:,:,:,1:number_qs)

  ist  = 1
  iend = total_number_tracers

  ! put the sources into one superarray conistent with the field
  ! superarray

  qs_n(:,:,:,ist:iend) = super_tracer_phys1(                    &
                       tdims_s%i_start:tdims_s%i_end,           &
                       tdims_s%j_start:tdims_s%j_end,           &
                       tdims_s%k_start:tdims_s%k_end, ist:iend )

  qs_s( tdims%i_start:tdims%i_end,                              &
        tdims%j_start:tdims%j_end,                              &
        tdims%k_start:tdims%k_end, ist:iend ) =                 &
  super_tracer_phys2(                                           &
           tdims%i_start:tdims%i_end,                           &
           tdims%j_start:tdims%j_end,                           &
           tdims%k_start:tdims%k_end, ist:iend      )

  ! Maybe this swap_bounds is not necessary, but it's here for safety.
  ! If all the tracers are definitely swap_bounded before then this call
  ! can be removed

  CALL Swap_Bounds(qs_np1,                                          &
       tdims%i_len,tdims%j_len,                                     &
       number_qs*tdims_s%k_len,                                     &
       tdims%i_start-tdims_s%i_start,tdims%j_start-tdims_s%j_start, &
       fld_type_p, swap_field_is_scalar  )

  ! print mass conservation error before correction if needed

  IF (printstatus > prstatus_normal) THEN
    CALL umPrint(                                                    &
         'Error in mass conservation for tracers before correction', &
         src='eg_correct_tracers',pe=0)
    CALL eg_check_mass_conservation_fix(rho_n, rho_np1, qs_n,   &
                                    qs_np1, qs_s,               &
                        number_qs, mype,'eg_correct_tracers_fix')
  END IF

  ! enfore mass conservation

  CALL eg_mass_conservation_fix(rho_n, rho_np1, qs_n, qs_np1,   &
                            qs_s, qsmin, number_qs,             &
                            L_conserv_smooth_lap                )

  ! print mass conservation error after correction if needed

  IF (printstatus > prstatus_normal) THEN
    CALL umPrint( &
        'Error in mass conservation for tracers after correction', &
        src='eg_correct_tracers',pe=0)
    CALL eg_check_mass_conservation_fix( rho_n, rho_np1,  &
                                     qs_n, qs_np1, qs_s,  &
                                      number_qs, mype,    &
                                 'eg_correct_tracers_fix' )
  END IF

  super_array(:,:,:,1:number_qs) = qs_np1

  ! may be there is no need for this swap-bound if the tracers
  ! are swap-bounded later on -- but it's safer to leave it

  CALL Swap_Bounds(super_array,                                    &
      tdims%i_len,tdims%j_len,                                     &
      number_qs*(tdims%k_len),                                     &
      tdims%i_start-tdims_s%i_start,tdims%j_start-tdims_s%j_start, &
      fld_type_p, swap_field_is_scalar   )

  DEALLOCATE( qsmin     )
  DEALLOCATE( qs_s      )
  DEALLOCATE( qs_np1    )
  DEALLOCATE( qs_n      )

ELSE

  ! if L_conserve_tracers==.false. then simply make sure that the tracers
  ! are strictly positive

  super_array = MAX(super_array, 0.0 )

END IF ! L_conserve_tracers

! unpack/ungroup the superarray. Copy the superarray fields into their
! appropriate arrays (i.e.,reverse the operation in "eg_group_tracers")

CALL eg_ungroup_tracers(                                        &
                    total_number_tracers,                       &
                    super_array,                                &
                    L_co2_interactive, co2,                     &
                    L_Murk_advect, murk,                        &
                    L_Soot, soot_new, soot_agd, soot_cld,       &
                    L_sulpc_so2, so2, so4_aitken, so4_accu,     &
                                 so4_diss,                      &
                    L_sulpc_nh3, nh3,                           &
                    L_sulpc_dms, dms,                           &
                    L_dust, dust_div1, dust_div2, dust_div3,    &
                            dust_div4, dust_div5, dust_div6,    &
                    L_biomass, bmass_new, bmass_agd, bmass_cld, &
                    L_ocff, ocff_new, ocff_agd, ocff_cld,       &
                    L_nitrate, nitr_acc, nitr_diss,             &
                    L_use_cariolle, ozone_tracer,               &
                    tracers, tr_vars,                           &
                    tr_ukca, tracer_ukca,                       &
                    L_twobin_dust                               )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_correct_tracers_fix

!======================================================================
!
! Version to be used if l_fix_conserv = .FALSE.
!
SUBROUTINE eg_correct_tracers(                                  &
                         mype, super_array_size,                &
                         super_tracer_phys1, super_tracer_phys2,&
                         rho_n, rho_np1,                        &
                         co2, L_CO2_interactive,                &
                         murk, L_murk_advect,                   &
                         soot_new, soot_agd, soot_cld, L_soot,  &
                         bmass_new, bmass_agd, bmass_cld,       &
                         L_biomass,                             &
                         ocff_new, ocff_agd, ocff_cld, l_ocff,  &
                         dust_div1,dust_div2,dust_div3,         &
                         dust_div4,dust_div5,dust_div6,         &
                         L_dust,                                &
                         so2, so4_aitken, so4_accu,             &
                         so4_diss, nh3, dms,                    &
                         L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms, &
                         nitr_acc, nitr_diss, L_nitrate,        &
                         L_use_cariolle, ozone_tracer,          &
                         tracers, tr_vars,                      &
                         tr_ukca, tracer_ukca,                  &
                         L_conserve_tracers,                    &
                         L_conserv_smooth_lap                   )

USE eg_group_tracers_mod
USE eg_ungroup_tracers_mod
USE eg_check_conserv_mod
USE eg_mass_conserv_mod
USE dust_parameters_mod, ONLY: L_twobin_dust
USE atm_fields_bounds_mod
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE umPrintMgr
USE Field_Types
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CORRECT_TRACERS'

INTEGER, INTENT(IN)  :: mype
INTEGER, INTENT(IN)  :: super_array_size
INTEGER, INTENT(IN)  :: tr_ukca, tr_vars
LOGICAL, INTENT(IN)  :: L_CO2_interactive
LOGICAL, INTENT(IN)  :: L_murk_advect, L_Soot
LOGICAL, INTENT(IN)  :: L_biomass, L_ocff, l_dust, L_sulpc_so2
LOGICAL, INTENT(IN)  :: L_sulpc_nh3, l_sulpc_dms, l_use_cariolle
LOGICAL, INTENT(IN)  :: L_nitrate, L_conserve_tracers
LOGICAL, INTENT(IN)  :: L_conserv_smooth_lap

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

REAL ::    super_array (tdims_s%i_start:tdims_s%i_end,          &
                        tdims_s%j_start:tdims_s%j_end,          &
                        tdims_s%k_start:tdims_s%k_end,          &
                        super_array_size)

INTEGER :: tracers_switches(super_array_size)
INTEGER :: number_qs, total_number_tracers
INTEGER :: i, j, k, l, mpierr

REAL,    ALLOCATABLE :: qs_n(:,:,:,:)
REAL,    ALLOCATABLE :: qs_np1(:,:,:,:)
REAL,    ALLOCATABLE :: qs_s(:,:,:,:)
REAL,    ALLOCATABLE :: qsmin(:)
INTEGER, ALLOCATABLE :: qswitches(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! pack/group all the tracers sources into a superarray

CALL eg_group_tracers(                                          &
                    super_array_size,                           &
                    total_number_tracers,                       &
                    super_array,                                &
                    L_co2_interactive, co2,                     &
                    L_Murk_advect, murk,                        &
                    L_Soot, soot_new, soot_agd, soot_cld,       &
                    L_sulpc_so2, so2, so4_aitken, so4_accu,     &
                                 so4_diss,                      &
                    L_sulpc_nh3, nh3,                           &
                    L_sulpc_dms, dms,                           &
                    L_dust, dust_div1, dust_div2, dust_div3,    &
                            dust_div4, dust_div5, dust_div6,    &
                    L_biomass, bmass_new, bmass_agd, bmass_cld, &
                    L_ocff, ocff_new, ocff_agd, ocff_cld,       &
                    L_nitrate, nitr_acc, nitr_diss,             &
                    L_use_cariolle, ozone_tracer,               &
                    tracers, tr_vars,                           &
                    tr_ukca, tracer_ukca,                       &
                    L_twobin_dust, tracers_switches             )


IF (L_conserve_tracers ) THEN

  number_qs =  total_number_tracers

  ALLOCATE(                                                     &
         qs_n(tdims_s%i_start:tdims_s%i_end,                    &
              tdims_s%j_start:tdims_s%j_end,                    &
              tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(                                                     &
       qs_np1(tdims_s%i_start:tdims_s%i_end,                    &
              tdims_s%j_start:tdims_s%j_end,                    &
              tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(                                                     &
         qs_s(tdims_s%i_start:tdims_s%i_end,                    &
              tdims_s%j_start:tdims_s%j_end,                    &
              tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(                                                    &
      qsmin(number_qs))

  ALLOCATE(                                                     &
        qswitches(number_qs))

  qswitches = tracers_switches(1:number_qs)
  qsmin     = 0.0
!$OMP  PARALLEL DEFAULT(NONE)                                   &
!$OMP  SHARED( number_qs, tdims_s, tdims, qs_np1, super_array,  &
!$OMP          qs_n, super_tracer_phys1, qs_s,                  &
!$OMP          super_tracer_phys2 )                             &
!$OMP  PRIVATE( i, j, k, l )
  DO l = 1, number_qs
!$OMP  DO SCHEDULE(STATIC)
    DO k = tdims_s%k_start, tdims_s%k_end
      DO j = tdims_s%j_start, tdims_s%j_end
        DO i = tdims_s%i_start, tdims_s%i_end
          qs_np1(i,j,k,l) = super_array(i,j,k,l)
  
! put the sources into one superarray conistent with the field
  ! superarray
          qs_n(i,j,k,l) = super_tracer_phys1(i,j,k,l)

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
    ! super_tracer_phys2 has no halos
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qs_s(i,j,k,l) = super_tracer_phys2(i,j,k,l)
        END DO
      END DO
    END DO
!$OMP END DO 


  ! Enforce the boundary condition filed(surface)=field(level 1)
  !
  ! These are the assumption used in computing the total mass.
  ! In some array they are forced earlier on and in some arrays
  ! they are not --- and this is why they are forced here
  ! so this routine doesn't miss-diagnose the total mass of tracers

!$OMP  DO SCHEDULE(STATIC)
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        qs_n(i,j,tdims_s%k_start,l)   = qs_n(i,j,tdims_s%k_start+1,l)
        qs_np1(i,j,tdims_s%k_start,l) = qs_np1(i,j,tdims_s%k_start+1,l)
        qs_s(i,j,tdims_s%k_start,l)   = qs_s(i,j,tdims_s%k_start+1,l) 
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL

  ! maybe this swap_bounds is not necessary, but it's here for safety after
  ! the boundary conditions are explicitly imposed -- because some tracers
  ! don't have anything at level 0 higher-up and this swap-bounds to make
  ! sure all the points have sensible data so mass is not missdiagnosed.

  CALL Swap_Bounds(qs_np1,                                          &
       tdims%i_len,tdims%j_len,                                     &
       number_qs*tdims%k_len,                                       &
       tdims%i_start-tdims_s%i_start,tdims%j_start-tdims_s%j_start, &
       fld_type_p, swap_field_is_scalar    )

  ! print mass conservation error before correction if needed

  IF (printstatus > prstatus_normal) THEN
    CALL umPrint(                                                   &
        'Error in mass conservation for tracers before correction', &
        src='eg_correct_tracers',pe=0)
    CALL eg_check_mass_conservation(rho_n, rho_np1, qs_n,       &
                                    qs_np1, qs_s,               &
                 qswitches, number_qs, mype,'eg_correct_tracers')
  END IF

  ! enfore mass conservation

  CALL eg_mass_conservation(rho_n, rho_np1, qs_n, qs_np1,       &
                            qs_s, qsmin,qswitches, number_qs,   &
                            L_conserv_smooth_lap                )

  ! print mass conservation error after correction if needed

  IF (printstatus > prstatus_normal) THEN
    CALL umPrint( &
        'Error in mass conservation for tracers after correction', &
        src='eg_correct_tracers',pe=0)
    CALL eg_check_mass_conservation( rho_n, rho_np1,      &
        qs_n, qs_np1, qs_s,  &
        qswitches, number_qs, mype, &
        'eg_correct_tracers' )
  END IF

!$OMP  PARALLEL DEFAULT(NONE)                                   &
!$OMP& SHARED( qs_np1, super_array, tdims_s, number_qs )        &
!$OMP& PRIVATE( i, j, k, l )
  DO l = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims_s%k_start, tdims_s%k_end
      DO j = tdims_s%j_start, tdims_s%j_end
        DO i = tdims_s%i_start, tdims_s%i_end
          super_array(i,j,k,l) = qs_np1(i,j,k,l)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL

  ! may be there is no need for this swap-bounds if the tracers
  ! are swap-bounded later on -- but it's safer to leave it

  CALL Swap_Bounds(super_array,                                    &
      tdims%i_len,tdims%j_len,                                     &
      number_qs*tdims%k_len,                                       &
      tdims%i_start-tdims_s%i_start,tdims%j_start-tdims_s%j_start, &
      fld_type_p, swap_field_is_scalar  )

  DEALLOCATE( qswitches )
  DEALLOCATE( qsmin     )
  DEALLOCATE( qs_s      )
  DEALLOCATE( qs_np1    )
  DEALLOCATE( qs_n      )

ELSE

  ! if L_conserve_tracers==.false. then simply make sure that the tracers
  ! are strictly positive

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)               &
!$OMP& SHARED( super_array_size, tdims_s, super_array )         &
!$OMP& PRIVATE( i, j, k, l ) COLLAPSE(2)
  DO l = 1, super_array_size
    DO k = tdims_s%k_start, tdims_s%k_end
      DO j = tdims_s%j_start, tdims_s%j_end
        DO i = tdims_s%i_start, tdims_s%i_end
          super_array(i,j,k,l) = MAX(super_array(i,j,k,l), 0.0)
        END DO
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! L_conserve_tracers

! unpack/ungroup the superarray. Copy the superarray fields into their
! appropriate arrays (i.e.,reverse the operation in "eg_group_tracers")

CALL eg_ungroup_tracers(                                        &
                    total_number_tracers,                       &
                    super_array,                                &
                    L_co2_interactive, co2,                     &
                    L_Murk_advect, murk,                        &
                    L_Soot, soot_new, soot_agd, soot_cld,       &
                    L_Sulpc_so2, so2, so4_aitken, so4_accu,     &
                                 so4_diss,                      &
                    L_sulpc_nh3, nh3,                           &
                    L_sulpc_dms, dms,                           &
                    L_dust, dust_div1, dust_div2, dust_div3,    &
                            dust_div4, dust_div5, dust_div6,    &
                    L_biomass, bmass_new, bmass_agd, bmass_cld, &
                    L_ocff, ocff_new, ocff_agd, ocff_cld,       &
                    L_nitrate, nitr_acc, nitr_diss,             &
                    L_use_cariolle, ozone_tracer,               &
                    tracers, tr_vars,                           &
                    tr_ukca, tracer_ukca,                       &
                    L_twobin_dust                               )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_correct_tracers

END MODULE eg_conserv_tracers_mod

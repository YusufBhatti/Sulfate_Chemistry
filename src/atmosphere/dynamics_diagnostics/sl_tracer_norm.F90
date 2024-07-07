! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Dynamics Diagnostics
!
! Subroutine Interface:
MODULE sl_tracer_norm_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SL_TRACER_NORM_MOD'

CONTAINS

SUBROUTINE sl_tracer_norm(                                              &
                          level_start, level_end, rims,                 &
                          super_array_size, co2, murk,                  &
                          dust_div1, dust_div2, dust_div3,              &
                          dust_div4, dust_div5, dust_div6,              &
                          soot_new, soot_agd, soot_cld,                 &
                          bmass_new, bmass_agd, bmass_cld,              &
                          ocff_new, ocff_agd, ocff_cld,                 &
                          so2, so4_aitken, so4_accu,                    &
                          so4_diss, nh3, dms,                           &
                          nitr_acc, nitr_diss, ozone_tracer,            &
                          tracers, tr_vars,                             &
                          tracers_ukca, tr_ukca,                        &
                          L_do_rims, L_do_halos, L_print_pe )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod,  ONLY: tdims_s
USE array_norm_mod, ONLY: array_norm
USE pr_block4_mod, ONLY: dom_no, dom_eo, dom_s, dom_w
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels,        &
                                 global_row_length, global_rows
USE um_parcore, ONLY: mype
USE um_parvars, ONLY: offx, offy
USE carbon_options_mod, ONLY: l_co2_interactive
USE dust_parameters_mod, ONLY: l_dust, l_twobin_dust
USE murk_inputs_mod, ONLY: l_murk_advect
USE run_aerosol_mod, ONLY: l_sulpc_so2, l_soot, l_ocff, l_biomass,      &
                           l_sulpc_nh3, l_nitrate, l_sulpc_dms
USE rad_input_mod, ONLY: l_use_cariolle

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr, ONLY: umPrint, umMessage, umPrintFlush

!
! Description: Printing of norms for tracer arrays after sl_tracer1_4A
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SL_TRACER_NORM'

! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_do_rims
LOGICAL :: L_do_halos
LOGICAL :: L_print_pe   ! print diagnostics on all pe's if true

INTEGER, INTENT(IN) :: super_array_size
INTEGER, INTENT(IN) :: level_start, level_end, rims

INTEGER ::                                                        &
  tr_vars                                                         &
, tr_ukca              ! No of UKCA tracers

REAL, INTENT(IN) :: co2                                         &
                  (tdims_s%i_start:tdims_s%i_end,               &
                   tdims_s%j_start:tdims_s%j_end,               &
                   tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: murk                                       &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: soot_new                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: soot_agd                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: soot_cld                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so2                                     &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so4_aitken                              &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so4_accu                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so4_diss                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: nh3                                     &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dms                                     &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div1                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div2                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div3                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div4                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div5                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div6                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: bmass_new                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: bmass_agd                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: bmass_cld                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN)  :: ocff_new                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: ocff_agd                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: ocff_cld                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: nitr_acc                                &
                 (tdims_s%i_start:tdims_s%i_end,                &
                  tdims_s%j_start:tdims_s%j_end,                &
                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN)  :: nitr_diss                               &
                 (tdims_s%i_start:tdims_s%i_end,                &
                  tdims_s%j_start:tdims_s%j_end,                &
                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  ::                                         &
 tracers(tdims_s%i_start:tdims_s%i_end,                         &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end,tr_vars)           &
, tracers_ukca(tdims_s%i_start:tdims_s%i_end,                   &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end,tr_ukca)           &
! Add cariolle specific parameters for ozone tracer
            , ozone_tracer(tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end)

! Arguments with INTENT OUT. ie: Output variables.

! Local Variables.
INTEGER :: k, num, counter
INTEGER :: levels      !  model_levels + 1 
INTEGER :: start_level, end_level
INTEGER :: k0, s0      !  level zero pointers
INTEGER :: rim_n, rim_s, rim_e, rim_w

! Local array for storing two_norms
REAL ::  two_norm(super_array_size)

! ----------------------------------------------------------------------
! 1.0 Start of subroutine code
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( L_do_rims ) THEN
  rim_n = 1
  rim_s = 1
  rim_e = 1
  rim_w = 1
ELSE ! L_do_rims
  rim_n = rims
  rim_s = rims
  rim_e = rims
  rim_w = rims
!  set up block norm if dom_w and dom_s /= 1
  IF ( dom_w /= 1 .AND. dom_s /= 1 ) THEN
    rim_s = dom_s
    rim_w = dom_w
    rim_n = global_rows - dom_no
    rim_e = global_row_length - dom_eo
  END IF ! dom_w /= 1 .AND. dom_s /= 1
END IF ! L_do_rims

levels = model_levels + 1
start_level = level_start
end_level = level_end
s0 = 0
IF ( level_start == level_end ) THEN
  IF ( level_start < 1 ) THEN
    start_level = 1
    end_level   = 1
  ELSE IF ( level_start == model_levels) THEN
    start_level = 1
    s0 = 1
    levels = model_levels
  END IF ! level_start == model_levels 
ELSE     !    0 < level_start < model_levels
  s0 = 1
  levels = model_levels
END IF   !   level_start == level_end

IF ( level_start == level_end ) THEN

  DO k =  start_level, end_level
    k0 = k
    IF( level_start < 1 ) k0 = 0

    num = 0
    IF (l_CO2_interactive) THEN
      num = num + 1
      CALL array_norm(                                               &
                      co2(1-offx,1-offy,s0), row_length, rows,       &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! l_CO2_interactive
    IF (l_Soot) THEN
      num = num + 1
      CALL array_norm(                                               &
                      soot_new(1-offx,1-offy,s0), row_length, rows,  &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      soot_agd(1-offx,1-offy,s0), row_length, rows,  &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      soot_cld(1-offx,1-offy,s0), row_length, rows,  &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! l_Soot
    IF (l_Biomass) THEN
      num = num + 1
      CALL array_norm(                                               &
                      bmass_new(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      bmass_agd(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      bmass_cld(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! l_Biomass

    IF (l_Sulpc_so2) THEN
      num = num + 1
      CALL array_norm(                                               &
                      so4_aitken(1-offx,1-offy,s0), row_length, rows,&
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      so4_accu(1-offx,1-offy,s0), row_length, rows,  &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      so4_diss(1-offx,1-offy,s0), row_length, rows,  &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      so2(1-offx,1-offy,s0), row_length, rows,       &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      IF (L_sulpc_nh3) THEN
        num = num + 1
        CALL array_norm(                                               &
                        nh3(1-offx,1-offy,s0), row_length, rows,       &
                        levels, k, k,                                  &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
      END IF  ! l_Sulpc_nh3
      IF (L_sulpc_dms) THEN
        num = num + 1
        CALL array_norm(                                               &
                        dms(1-offx,1-offy,s0), row_length, rows,       &
                        levels, k, k,                                  &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
      END IF  ! l_Sulpc_dms
    END IF  ! l_Sulpc_so2

    IF (l_dust) THEN
      num = num + 1
      CALL array_norm(                                             &
                      dust_div1(1-offx,1-offy,s0), row_length, rows,&
                      levels, k, k,                                &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                             &
                      dust_div2(1-offx,1-offy,s0), row_length, rows,&
                      levels, k, k,                                &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      IF (.NOT. l_twobin_dust) THEN
        num = num + 1
        CALL array_norm(                                            &
                        dust_div3, row_length, rows,                &
                        levels, k, k,                               &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
        num = num + 1
        CALL array_norm(                                            &
                        dust_div4, row_length, rows,                &
                        levels, k, k,                               &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
        num = num + 1
        CALL array_norm(                                            &
                        dust_div5, row_length, rows,                &
                        levels, k, k,                               &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
        num = num + 1
        CALL array_norm(                                            &
                        dust_div6, row_length, rows,                &
                        levels, k, k,                               &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
      END IF ! .NOT. l_twobin_dust
    END IF ! L_DUST

    IF (L_ocff) THEN
      num = num + 1
      CALL array_norm(                                              &
                      ocff_new(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                              &
                      ocff_agd(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                              &
                      ocff_cld(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF !   L_ocff

    IF (l_use_cariolle) THEN
      num = num + 1
      CALL array_norm(                                             &
                      ozone_tracer, row_length, rows,              &
                      levels, k, k,                                &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF ! L_use_cariolle

    IF (L_nitrate) THEN
      num = num + 1
      CALL array_norm(                                              &
                      nitr_acc(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                              &
                      nitr_diss(1-offx,1-offy,s0), row_length, rows,&
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! L_nitrate

    IF ( tr_vars > 0 ) THEN
      DO counter = 1, tr_vars
        num = num + 1
        CALL array_norm(                                             &
                        tracers(1-offx,1-offy,s0,counter),             &
                        row_length, rows, levels, k, k,              &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
      END DO !  counter = 1, tr_vars
    END IF   ! tr_vars > 0

    IF ( tr_ukca > 0 ) THEN
      DO counter = 1, tr_ukca
        num = num + 1
        CALL array_norm(                                             &
                        tracers_ukca(1-offx,1-offy,s0,counter),        &
                        row_length, rows, levels, k, k,              &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
      END DO !  counter = 1, tr_ukca
    END IF   ! tr_ukca > 0

    IF ( l_murk_advect ) THEN
        num = num + 1
        CALL array_norm(                                            &
                        murk(1-offx,1-offy,s0), row_length, rows,   &
                        levels, k, k,                               &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! l_murk_advect

    IF ( L_print_pe .OR. mype == 0 ) THEN
      num = 0
      IF (l_CO2_interactive) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' co2 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END IF  ! l_CO2_interactive
      IF (l_Soot) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' soot_new Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' soot_agd Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' soot_cld Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END IF  ! l_Soot
      IF (l_Biomass) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' bmass_new Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' bmass_agd Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' bmass_cld Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END IF  ! l_Biomass

      IF (l_Sulpc_so2) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' so4_aitken Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' so4_accu Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' so4_diss Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' so2 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        IF (L_sulpc_nh3) THEN
          num = num + 1
          WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' nh3 Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER_NORM')
        END IF  ! l_Sulpc_nh3
        IF (L_sulpc_dms) THEN
          num = num + 1
          WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' dms Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER_NORM')
        END IF  ! l_Sulpc_dms
      END IF  ! l_Sulpc_so2

      IF (l_dust) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' dust_div1 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' dust_div2 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        IF (.NOT. l_twobin_dust) THEN
          num = num + 1
          WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' dust_div3 Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER_NORM')
          num = num + 1
          WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' dust_div4 Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER_NORM')
          num = num + 1
          WRITE(umMessage, '(A, I5, A, E23.16)')                       &
           ' Level ', k0, ' dust_div5 Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER_NORM')
          num = num + 1
          WRITE(umMessage, '(A, I5, A, E23.16)')                       &
           ' Level ', k0, ' dust_div6 Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER_NORM')
        END IF ! .NOT. l_twobin_dust
      END IF ! L_DUST

      IF (L_ocff) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' ocff_new Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' ocff_agd Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' ocff_cld Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END IF !   L_ocff

      IF (l_use_cariolle) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' ozone_tracer Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END IF ! L_use_cariolle

      IF (L_nitrate) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' nitr_acc Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' nitr_diss Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END IF  ! L_nitrate

      IF ( tr_vars > 0 ) THEN
        DO counter = 1, tr_vars
          num = num + 1
          WRITE(umMessage, '(A, I5, A, I2, A, E23.16)')                &
                ' Level ', k0, ' tracer' , counter,                      & 
                ' Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER_NORM')
        END DO !  counter = 1, tr_vars
      END IF   ! tr_vars > 0

      IF ( tr_ukca > 0 ) THEN
        DO counter = 1, tr_ukca
          num = num + 1
          WRITE(umMessage, '(A, I5, A, I2, A, E23.16)')                &
                ' Level ', k0, ' tracer_ukca' , counter,                 & 
                ' Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER_NORM')
        END DO !  counter = 1, tr_ukca
      END IF   ! tr_ukca > 0

      IF ( l_murk_advect ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' murk Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END IF  ! l_murk_advect

    END IF ! L_print_pe .OR. mype == 0

  END DO !  k =  start_level, end_level

ELSE   !  3d norms from  level_start to level_end

  num = 0
  IF (l_CO2_interactive) THEN
    num = num + 1
    CALL array_norm(                                               &
                    co2(1-offx,1-offy,s0), row_length, rows,       &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF  ! l_CO2_interactive

  IF (l_Soot) THEN
    num = num + 1
    CALL array_norm(                                                &
                    soot_new(1-offx,1-offy,s0), row_length, rows,   &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    soot_agd(1-offx,1-offy,s0), row_length, rows,   &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    soot_cld(1-offx,1-offy,s0), row_length, rows,   &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF  ! l_Soot
 
  IF (l_Biomass) THEN
    num = num + 1
    CALL array_norm(                                                &
                    bmass_new(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    bmass_agd(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    bmass_cld(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF  ! l_Biomass

  IF (l_Sulpc_so2) THEN
    num = num + 1
    CALL array_norm(                                               &
                    so4_aitken(1-offx,1-offy,s0), row_length, rows, &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    so4_accu(1-offx,1-offy,s0), row_length, rows,   &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    so4_diss(1-offx,1-offy,s0), row_length, rows,   &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    so2(1-offx,1-offy,s0), row_length, rows,        &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    IF (L_sulpc_nh3) THEN
      num = num + 1
      CALL array_norm(                                               &
                      nh3(1-offx,1-offy,s0), row_length, rows,        &
                      levels, start_level, level_end,                &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! l_Sulpc_nh3
    IF (L_sulpc_dms) THEN
      num = num + 1
      CALL array_norm(                                               &
                      dms(1-offx,1-offy,s0), row_length, rows,       &
                      levels, start_level, level_end,                &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! l_Sulpc_dms
  END IF  ! l_Sulpc_so2

  IF (l_dust) THEN
    num = num + 1
    CALL array_norm(                                                &
                    dust_div1(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    dust_div1(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    IF (.NOT. l_twobin_dust) THEN
      num = num + 1
      CALL array_norm(                                                &
                      dust_div3(1-offx,1-offy,s0), row_length, rows,  &
                      levels, start_level, level_end,                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                                &
                      dust_div3(1-offx,1-offy,s0), row_length, rows,  &
                      levels, start_level, level_end,                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                                &
                      dust_div5(1-offx,1-offy,s0), row_length, rows,  &
                      levels, start_level, level_end,                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                                &
                      dust_div6(1-offx,1-offy,s0), row_length, rows,  &
                      levels, start_level, level_end,                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF ! .NOT. l_twobin_dust
  END IF ! L_DUST

  IF (L_ocff) THEN
    num = num + 1
    CALL array_norm(                                               &
                    ocff_new(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    ocff_agd(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    ocff_cld(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF !   L_ocff

  IF (l_use_cariolle) THEN
    num = num + 1
    CALL array_norm(                                                  &
                    ozone_tracer(1-offx,1-offy,s0), row_length, rows, &
                    levels, start_level, level_end,                   &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,          &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF ! L_use_cariolle

  IF (L_nitrate) THEN
    num = num + 1
    CALL array_norm(                                                &
                    nitr_acc(1-offx,1-offy,s0), row_length, rows,   &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    nitr_diss(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF  ! L_nitrate

  IF ( tr_vars > 0 ) THEN
    DO counter = 1, tr_vars
      num = num + 1
      CALL array_norm(                                             &
                      tracers(1-offx,1-offy,s0,counter),             &
                      row_length, rows, model_levels,              &
                      level_start, level_end,                      &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END DO !  counter = 1, tr_vars
  END IF   ! tr_vars > 0

  IF ( tr_ukca > 0 ) THEN
    DO counter = 1, tr_ukca
      num = num + 1
      CALL array_norm(                                             &
                      tracers_ukca(1-offx,1-offy,s0,counter),        &
                      row_length, rows, levels,                    &
                      level_start, level_end,                      &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END DO !  counter = 1, tr_ukca
  END IF   ! tr_ukca > 0

  IF ( l_murk_advect ) THEN
      num = num + 1
      CALL array_norm(                                             &
                      murk(1-offx,1-offy,s0), row_length, rows,    &
                      levels, start_level, level_end,              &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF  ! l_murk_advect

  IF ( L_print_pe .OR. mype == 0 ) THEN
    num = 0
    IF (l_CO2_interactive) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' co2 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
    END IF  ! l_CO2_interactive
    IF (l_Soot) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' soot_new Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' soot_agd Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' soot_cld Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
    END IF  ! l_Soot

    IF (l_Biomass) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' bmass_new Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' bmass_agd Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' bmass_cld Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
    END IF  ! l_Biomass

    IF (l_Sulpc_so2) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' so4_aitken Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' so4_accu Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' so4_diss Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' so2 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      IF (L_sulpc_nh3) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,        &
                  ' nh3 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END IF  ! l_Sulpc_nh3
      IF (L_sulpc_dms) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,        &
                  ' dms Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END IF  ! l_Sulpc_dms
    END IF  ! l_Sulpc_so2

    IF (l_dust) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div1 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div2 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      IF (.NOT. l_twobin_dust) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')            &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div3 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div4 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div5 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div6 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END IF ! .NOT. l_twobin_dust
    END IF ! L_DUST

    IF (L_ocff) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' ocff_new Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' ocff_agd Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' ocff_cld Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
    END IF !   L_ocff

    IF (l_use_cariolle) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' ozone_tracer Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
    END IF ! L_use_cariolle

    IF (L_nitrate) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' nitr_acc Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' nitr_diss Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
    END IF  ! L_nitrate

    IF ( tr_vars > 0 ) THEN
      DO counter = 1, tr_vars
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, I2, A, E23.16)')     &
                  ' Levels', start_level, ' to', level_end,      &
                  ' tracer' , counter, ' Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END DO !  counter = 1, tr_vars
    END IF   ! tr_vars > 0

    IF ( tr_ukca > 0 ) THEN
      DO counter = 1, tr_ukca
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, I2, A, E23.16)')     &
                  ' Levels', start_level, ' to', level_end,      &
           ' tracer_ukca' , counter, ' Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER_NORM')
      END DO !  counter = 1, tr_ukca
    END IF   ! tr_ukca > 0

    IF ( l_murk_advect ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' murk Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER_NORM')
    END IF  ! l_murk_advect

  END IF ! L_print_pe .OR. mype == 0

END IF ! level_start == level_end

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE sl_tracer_norm

SUBROUTINE sl_tracer1_norm(                                            &
                           super_array_size,                           &
                           level_start, level_end, rims,               &
                           co2, aerosol, dust_div1, dust_div2,         &
                           dust_div3, dust_div4, dust_div5, dust_div6, &
                           soot_new, soot_agd, soot_cld,               &
                           bmass_new, bmass_agd, bmass_cld,            &
                           ocff_new, ocff_agd, ocff_cld,               &
                           so2, so4_aitken, so4_accu, so4_diss,        &
                           nh3, nitr_acc, nitr_diss,                   &
                           tracers_ukca, tr_ukca, ozone_tracer,        &
                           L_do_rims, L_do_halos, L_print_pe )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE carbon_options_mod, ONLY: l_co2_interactive
USE dust_parameters_mod, ONLY: l_dust, l_twobin_dust
USE murk_inputs_mod, ONLY: l_murk_advect
USE run_aerosol_mod, ONLY: l_sulpc_so2, l_soot, l_ocff, l_biomass,      &
                           l_sulpc_nh3, l_nitrate
USE rad_input_mod, ONLY: l_use_cariolle

USE atm_fields_bounds_mod,  ONLY: tdims_s
USE array_norm_mod
USE pr_block4_mod
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels,        &
                                 global_row_length, global_rows
USE um_parcore, ONLY: mype
USE um_parvars, ONLY: offx, offy

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr, ONLY: umPrint, umMessage, umPrintFlush
!
! Description: Printing of norms for tracer arrays after sl_tracer1_4A
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SL_TRACER1_NORM'

! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_do_rims
LOGICAL :: L_do_halos
LOGICAL :: L_print_pe   ! print diagnostics on all pe's if true

INTEGER, INTENT(IN) :: super_array_size
INTEGER, INTENT(IN) :: level_start, level_end, rims

INTEGER ::  tr_ukca        ! No of UKCA tracers

REAL, INTENT(IN) :: co2                                      &
                  (tdims_s%i_start:tdims_s%i_end,               &
                   tdims_s%j_start:tdims_s%j_end,               &
                   tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: aerosol                                    &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: soot_new                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: soot_agd                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: soot_cld                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so2                                     &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so4_aitken                              &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so4_accu                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so4_diss                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: nh3                                     &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div1                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div2                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div3                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div4                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div5                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div6                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: bmass_new                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: bmass_agd                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: bmass_cld                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN)  :: ocff_new                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: ocff_agd                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: ocff_cld                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: nitr_acc                                &
                 (tdims_s%i_start:tdims_s%i_end,                &
                  tdims_s%j_start:tdims_s%j_end,                &
                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN)  :: nitr_diss                               &
                 (tdims_s%i_start:tdims_s%i_end,                &
                  tdims_s%j_start:tdims_s%j_end,                &
                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  ::                                         &
  tracers_ukca(tdims_s%i_start:tdims_s%i_end,                   &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end,tr_ukca)           &

! Add cariolle specific parameters for ozone tracer
            , ozone_tracer(tdims_s%i_start:tdims_s%i_end,       &
                           tdims_s%j_start:tdims_s%j_end,        &
                           tdims_s%k_start:tdims_s%k_end)

! Arguments with INTENT OUT. ie: Output variables.

! Local Variables.
INTEGER :: k, num, counter
INTEGER :: levels      !  model_levels + 1 
INTEGER :: start_level, end_level
INTEGER :: k0, s0      !  level zero pointers
INTEGER :: rim_n, rim_s, rim_e, rim_w

! Local array for storing two_norms
REAL ::  two_norm(super_array_size)

! ----------------------------------------------------------------------
! 1.0 Start of subroutine code
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( L_do_rims ) THEN
  rim_n = 1
  rim_s = 1
  rim_e = 1
  rim_w = 1
ELSE ! L_do_rims
  rim_n = rims
  rim_s = rims
  rim_e = rims
  rim_w = rims
!  set up block norm if dom_w and dom_s /= 1
  IF ( dom_w /= 1 .AND. dom_s /= 1 ) THEN
    rim_s = dom_s
    rim_w = dom_w
    rim_n = global_rows - dom_no
    rim_e = global_row_length - dom_eo
  END IF ! dom_w /= 1 .AND. dom_s /= 1
END IF ! L_do_rims

levels = model_levels + 1
start_level = level_start
end_level = level_end
s0 = 0
IF ( level_start == level_end ) THEN
  IF ( level_start < 1 ) THEN
    start_level = 1
    end_level   = 1
  ELSE IF ( level_start == model_levels) THEN
    start_level = 1
    s0 = 1
    levels = model_levels
  END IF ! level_start == model_levels 
ELSE     !    0 < level_start < model_levels
  s0 = 1
  levels = model_levels
END IF   !   level_start == level_end

IF ( level_start == level_end ) THEN

  DO k =  start_level, end_level
    k0 = k
    IF( level_start < 1 ) k0 = 0

    num = 0
    IF (l_CO2_interactive) THEN
      num = num + 1
      CALL array_norm(                                               &
                      co2(1-offx,1-offy,s0), row_length, rows,       &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! l_CO2_interactive
    IF (l_Soot) THEN
      num = num + 1
      CALL array_norm(                                               &
                      soot_new(1-offx,1-offy,s0), row_length, rows,  &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      soot_agd(1-offx,1-offy,s0), row_length, rows,  &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      soot_cld(1-offx,1-offy,s0), row_length, rows,  &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! l_Soot
    IF (l_Biomass) THEN
      num = num + 1
      CALL array_norm(                                               &
                      bmass_new(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      bmass_agd(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      bmass_cld(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! l_Biomass

    IF (l_Sulpc_so2) THEN
      num = num + 1
      CALL array_norm(                                               &
                      so4_aitken(1-offx,1-offy,s0), row_length, rows,&
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      so4_accu(1-offx,1-offy,s0), row_length, rows,  &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      so4_diss(1-offx,1-offy,s0), row_length, rows,  &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                               &
                      so2(1-offx,1-offy,s0), row_length, rows,       &
                      levels, k, k,                                  &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      IF (L_sulpc_nh3) THEN
        num = num + 1
        CALL array_norm(                                              &
                        nh3(1-offx,1-offy,s0), row_length, rows,      &
                        levels, k, k,                                 &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
      END IF  ! l_Sulpc_nh3
    END IF  ! l_Sulpc_so2

    IF (l_dust) THEN
      num = num + 1
      CALL array_norm(                                              &
                      dust_div1(1-offx,1-offy,s0), row_length, rows,&
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                              &
                      dust_div2(1-offx,1-offy,s0), row_length, rows,&
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      IF (.NOT. l_twobin_dust) THEN
        num = num + 1
        CALL array_norm(                                            &
                        dust_div3, row_length, rows,                &
                        levels, k, k,                               &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
        num = num + 1
        CALL array_norm(                                            &
                        dust_div4, row_length, rows,                &
                        levels, k, k,                               &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
        num = num + 1
        CALL array_norm(                                            &
                        dust_div5, row_length, rows,                &
                        levels, k, k,                               &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
        num = num + 1
        CALL array_norm(                                            &
                        dust_div6, row_length, rows,                &
                        levels, k, k,                               &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
      END IF ! .NOT. l_twobin_dust
    END IF ! L_DUST

    IF (L_ocff) THEN
      num = num + 1
      CALL array_norm(                                              &
                      ocff_new(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                              &
                      ocff_agd(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                              &
                      ocff_cld(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF !   L_ocff

    IF (l_use_cariolle) THEN
      num = num + 1
      CALL array_norm(                                             &
                      ozone_tracer, row_length, rows,              &
                      levels, k, k,                                &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF ! L_use_cariolle

    IF (L_nitrate) THEN
      num = num + 1
      CALL array_norm(                                              &
                      nitr_acc(1-offx,1-offy,s0), row_length, rows, &
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                              &
                      nitr_diss(1-offx,1-offy,s0), row_length, rows,&
                      levels, k, k,                                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! L_nitrate

    IF ( tr_ukca > 0 ) THEN
      DO counter = 1, tr_ukca
        num = num + 1
        CALL array_norm(                                            &
                        tracers_ukca(1-offx,1-offy,0,counter),        &
                        row_length, rows, levels, k, k,             &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
      END DO !  counter = 1, tr_ukca
    END IF   ! tr_ukca > 0

    IF ( l_murk_advect ) THEN
        num = num + 1
        CALL array_norm(                                             &
                        aerosol(1-offx,1-offy,s0), row_length, rows, &
                        levels, k, k,                                &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! l_murk_advect

    IF ( L_print_pe .OR. mype == 0 ) THEN
      num = 0
      IF (l_CO2_interactive) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' co2 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      END IF  ! l_CO2_interactive
      IF (l_Soot) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' soot_new Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' soot_agd Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' soot_cld Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      END IF  ! l_Soot
      IF (l_Biomass) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' bmass_new Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' bmass_agd Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' bmass_cld Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      END IF  ! l_Biomass

      IF (l_Sulpc_so2) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' so4_aitken Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' so4_accu Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' so4_diss Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' so2 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        IF (L_sulpc_nh3) THEN
          num = num + 1
          WRITE(umMessage, '(A, I5, A, E23.16)')                     &
              ' Level ', k0, ' nh3 Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        END IF  ! l_Sulpc_nh3
      END IF  ! l_Sulpc_so2

      IF (l_dust) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' dust_div1 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' dust_div2 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        IF (.NOT. l_twobin_dust) THEN
          num = num + 1
          WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' dust_div3 Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER1_NORM')
          num = num + 1
          WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' dust_div4 Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER1_NORM')
          num = num + 1
          WRITE(umMessage, '(A, I5, A, E23.16)')                       &
           ' Level ', k0, ' dust_div5 Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER1_NORM')
          num = num + 1
          WRITE(umMessage, '(A, I5, A, E23.16)')                       &
           ' Level ', k0, ' dust_div6 Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        END IF ! .NOT. l_twobin_dust
      END IF ! L_DUST

      IF (L_ocff) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' ocff_new Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' ocff_agd Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' ocff_cld Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      END IF !   L_ocff

      IF (l_use_cariolle) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' ozone_tracer Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      END IF ! L_use_cariolle

      IF (L_nitrate) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' nitr_acc Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
            ' Level ', k0, ' nitr_diss Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      END IF  ! L_nitrate

      IF ( tr_ukca > 0 ) THEN
        DO counter = 1, tr_ukca
          num = num + 1
          WRITE(umMessage, '(A, I5, A, I2, A, E23.16)')               &
                ' Level ', k0, ' tracer_ukca' , counter,                & 
                ' Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        END DO !  counter = 1, tr_ukca
      END IF   ! tr_ukca > 0

      IF ( l_murk_advect ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' aerosol Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      END IF  ! l_murk_advect

    END IF ! L_print_pe .OR. mype == 0

  END DO !  k = 0, model_levels

ELSE   !  3d norms from level_start to level_end

  num = 0
  IF (l_CO2_interactive) THEN
    num = num + 1
    CALL array_norm(                                                &
                    co2(1-offx,1-offy,s0), row_length, rows,        &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF  ! l_CO2_interactive

  IF (l_Soot) THEN
    num = num + 1
    CALL array_norm(                                                &
                    soot_new(1-offx,1-offy,s0), row_length, rows,   &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    soot_agd(1-offx,1-offy,s0), row_length, rows,   &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    soot_cld(1-offx,1-offy,s0), row_length, rows,   &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF  ! l_Soot
 
  IF (l_Biomass) THEN
    num = num + 1
    CALL array_norm(                                                &
                    bmass_new(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    bmass_agd(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    bmass_cld(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF  ! l_Biomass

  IF (l_Sulpc_so2) THEN
    num = num + 1
    CALL array_norm(                                                &
                    so4_aitken(1-offx,1-offy,s0), row_length, rows, &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    so4_accu(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    so4_diss(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    so2(1-offx,1-offy,s0), row_length, rows,       &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    IF (L_sulpc_nh3) THEN
      num = num + 1
      CALL array_norm(                                               &
                      nh3(1-offx,1-offy,s0), row_length, rows,       &
                      levels, start_level, level_end,                &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF  ! l_Sulpc_nh3
  END IF  ! l_Sulpc_so2

  IF (l_dust) THEN
    num = num + 1
    CALL array_norm(                                                &
                    dust_div1(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                 &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    dust_div1(1-offx,1-offy,s0), row_length, rows,   &
                    levels, start_level, level_end,                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,         &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    IF (.NOT. l_twobin_dust) THEN
      num = num + 1
      CALL array_norm(                                                &
                      dust_div3(1-offx,1-offy,s0), row_length, rows,  &
                      levels, start_level, level_end,                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                                &
                      dust_div4(1-offx,1-offy,s0), row_length, rows,  &
                      levels, start_level, level_end,                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                                &
                      dust_div5(1-offx,1-offy,s0), row_length, rows,  &
                      levels, start_level, level_end,                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
      num = num + 1
      CALL array_norm(                                                &
                      dust_div6(1-offx,1-offy,s0), row_length, rows,  &
                      levels, start_level, level_end,                 &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF ! .NOT. l_twobin_dust
  END IF ! L_DUST

  IF (L_ocff) THEN
    num = num + 1
    CALL array_norm(                                               &
                    ocff_new(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    ocff_agd(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    ocff_cld(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF !   L_ocff

  IF (l_use_cariolle) THEN
    num = num + 1
    CALL array_norm(                                                  &
                    ozone_tracer(1-offx,1-offy,s0), row_length, rows, &
                    levels, start_level, level_end,                   &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,          &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF ! L_use_cariolle

  IF (L_nitrate) THEN
    num = num + 1 
    CALL array_norm(                                               &
                    nitr_acc(1-offx,1-offy,s0), row_length, rows,  &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    nitr_diss(1-offx,1-offy,s0), row_length, rows, &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF  ! L_nitrate

  IF ( tr_ukca > 0 ) THEN
    DO counter = 1, tr_ukca
      num = num + 1
      CALL array_norm(                                             &
                      tracers_ukca(1-offx,1-offy,s0,counter),        &
                      row_length, rows, model_levels,              &
                      level_start, level_end,                      &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END DO !  counter = 1, tr_ukca
  END IF   ! tr_ukca > 0

  IF ( l_murk_advect ) THEN
      num = num + 1
      CALL array_norm(                                             &
                      aerosol(1-offx,1-offy,s0), row_length, rows, &
                      levels, start_level, level_end,              &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF  ! l_murk_advect

  IF ( L_print_pe .OR. mype == 0 ) THEN
    num = 0
    IF (l_CO2_interactive) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' co2 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
    END IF  ! l_CO2_interactive
    IF (l_Soot) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' soot_new Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' soot_agd Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' soot_cld Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
    END IF  ! l_Soot

    IF (l_Biomass) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' bmass_new Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' bmass_agd Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' bmass_cld Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
    END IF  ! l_Biomass

    IF (l_Sulpc_so2) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' so4_aitken Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' so4_accu Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' so4_diss Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' so2 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      IF (L_sulpc_nh3) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,        &
                  ' nh3 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      END IF  ! l_Sulpc_nh3
    END IF  ! l_Sulpc_so2

    IF (l_dust) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div1 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div2 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      IF (.NOT. l_twobin_dust) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')            &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div3 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div4 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div5 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dust_div6 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      END IF ! .NOT. l_twobin_dust
    END IF ! L_DUST

    IF (L_ocff) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' ocff_new Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' ocff_agd Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' ocff_cld Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
    END IF !   L_ocff

    IF (l_use_cariolle) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' ozone_tracer Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
    END IF ! L_use_cariolle

    IF (L_nitrate) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' nitr_acc Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' nitr_diss Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
    END IF  ! L_nitrate

    IF ( tr_ukca > 0 ) THEN
      DO counter = 1, tr_ukca
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, I2, A, E23.16)')     &
                  ' Levels', start_level, ' to', level_end,      &
           ' tracer_ukca' , counter, ' Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_TRACER1_NORM')
      END DO !  counter = 1, tr_ukca
    END IF   ! tr_ukca > 0

    IF ( l_murk_advect ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' aerosol Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_TRACER1_NORM')
    END IF  ! l_murk_advect

  END IF ! L_print_pe .OR. mype == 0

END IF ! _level_start == level_end

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE sl_tracer1_norm

SUBROUTINE free_tracer_norm(                                            &
                            free_tracers, tr_vars,                       &
                            level_start, level_end, rims,                &
                            L_do_rims, L_do_halos, L_print_pe )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod,  ONLY: tdims_s
USE array_norm_mod
USE pr_block4_mod
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels,        &
                                 global_row_length, global_rows
USE um_parcore, ONLY: mype
USE um_parvars, ONLY: offx, offy

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr, ONLY: umPrint, umMessage, umPrintFlush

!
! Description: Printing of free tracer norms 
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FREE_TRACER_NORM'

! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_do_rims
LOGICAL :: L_do_halos
LOGICAL :: L_print_pe   ! print diagnostics on all pe's if true

INTEGER, INTENT(IN) ::  level_start, level_end, rims
INTEGER, INTENT(IN) ::  tr_vars        ! No of free tracers

REAL, INTENT(IN)  :: free_tracers                                       &
                     (tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end, tr_vars)

! Local Variables.
INTEGER :: k, num, counter
INTEGER :: levels      !  model_levels + 1 
INTEGER :: start_level, end_level
INTEGER :: k0, s0      !  level zero pointers
INTEGER :: rim_n, rim_s, rim_e, rim_w

! Local array for storing two_norms
REAL ::  two_norm(tr_vars)

! ----------------------------------------------------------------------
! 1.0 Start of subroutine code
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( tr_vars > 0 ) THEN

IF ( L_do_rims ) THEN
  rim_n = 1
  rim_s = 1
  rim_e = 1
  rim_w = 1
ELSE ! L_do_rims
  rim_n = rims
  rim_s = rims
  rim_e = rims
  rim_w = rims
!  set up block norm if dom_w and dom_s /= 1
  IF ( dom_w /= 1 .AND. dom_s /= 1 ) THEN
    rim_s = dom_s
    rim_w = dom_w
    rim_n = global_rows - dom_no
    rim_e = global_row_length - dom_eo
  END IF ! dom_w /= 1 .AND. dom_s /= 1
END IF ! L_do_rims

levels = model_levels + 1
start_level = level_start
end_level = level_end
s0 = 0
IF ( level_start == level_end ) THEN
  IF ( level_start < 1 ) THEN
    start_level = 1
    end_level   = 1
  ELSE IF ( level_start == model_levels) THEN
    start_level = 1
    s0 = 1
    levels = model_levels
  END IF ! level_start == model_levels 
ELSE     !    0 < level_start < model_levels
  s0 = 1
  levels = model_levels
END IF   !   level_start == level_end

IF ( level_start == level_end ) THEN

  DO k =  start_level, end_level
    k0 = k
    IF( level_start < 1 ) k0 = 0

      num = 0
      DO counter = 1, tr_vars
        num = num + 1
        CALL array_norm(                                             &
                        free_tracers(1-offx,1-offy,s0,counter),        &
                        row_length, rows, levels, k, k,              &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
      END DO !  counter = 1, tr_vars

      IF ( L_print_pe .OR. mype == 0 ) THEN

        num = 0
        DO counter = 1, tr_vars
          num = num + 1
          WRITE(umMessage, '(A, I5, A, I2, A, E23.16)')               &
                ' Level ', k0, ' tracer' , counter,                     & 
                ' Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='FREE_TRACER_NORM')
        END DO !  counter = 1, tr_vars

      END IF ! L_print_pe .OR. mype == 0

    END DO !  k = 0, model_levels

  ELSE   !  3d norms from start_level to end_level  

    num = 0
    DO counter = 1, tr_vars
      num = num + 1
      CALL array_norm(                                            &
                      free_tracers(1-offx,1-offy,s0,counter),       &
                      row_length, rows, model_levels,             &
                      level_start, level_end,                     &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END DO !  counter = 1, tr_vars

    IF ( L_print_pe .OR. mype == 0 ) THEN
      num = 0
      DO counter = 1, tr_vars
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, I2, A, E23.16)')     &
                  ' Levels', start_level, ' to', level_end,      &
                  ' tracer' , counter, ' Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='FREE_TRACER_NORM')
      END DO !  counter = 1, tr_vars

    END IF ! L_print_pe .OR. mype == 0

  END IF ! start_level == end_level

END IF   ! tr_vars > 0

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE free_tracer_norm

SUBROUTINE ukca_tracer_norm(                                            &
                            tracer_ukca, tr_ukca,                       &
                            level_start, level_end, rims,                &
                            L_do_rims, L_do_halos, L_print_pe )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod
USE array_norm_mod
USE pr_block4_mod
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels,        &
                                 global_row_length, global_rows
USE um_parcore, ONLY: mype
USE um_parvars, ONLY: offx, offy

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr, ONLY: umPrint, umMessage, umPrintFlush

!
! Description: Printing of ukca_tracer norms
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_TRACER_NORM'

! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_do_rims
LOGICAL :: L_do_halos
LOGICAL :: L_print_pe   ! print diagnostics on all pe's if true

INTEGER, INTENT(IN) ::  level_start, level_end, rims
INTEGER, INTENT(IN) ::  tr_ukca        ! No of ukca tracers

REAL, INTENT(IN)  :: tracer_ukca                                        &
                     (tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end, tr_ukca)

! Local Variables.
INTEGER :: k, num, counter
INTEGER :: levels      !  model_levels + 1 
INTEGER :: start_level, end_level
INTEGER :: k0, s0      !  level zero pointers
INTEGER :: rim_n, rim_s, rim_e, rim_w

! Local array for storing two_norms
REAL ::  two_norm(tr_ukca)

! ----------------------------------------------------------------------
! 1.0 Start of subroutine code
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( tr_ukca > 0 ) THEN

IF ( L_do_rims ) THEN
  rim_n = 1
  rim_s = 1
  rim_e = 1
  rim_w = 1
ELSE ! L_do_rims
  rim_n = rims
  rim_s = rims
  rim_e = rims
  rim_w = rims
!  set up block norm if dom_w and dom_s /= 1
  IF ( dom_w /= 1 .AND. dom_s /= 1 ) THEN
    rim_s = dom_s
    rim_w = dom_w
    rim_n = global_rows - dom_no
    rim_e = global_row_length - dom_eo
  END IF ! dom_w /= 1 .AND. dom_s /= 1
END IF ! L_do_rims

levels = model_levels + 1
start_level = level_start
end_level = level_end
s0 = 0
IF ( level_start == level_end ) THEN
  IF ( level_start < 1 ) THEN
    start_level = 1
    end_level   = 1
  ELSE IF ( level_start == model_levels) THEN
    start_level = 1
    s0 = 1
    levels = model_levels
  END IF ! level_start == model_levels 
ELSE     !    0 < level_start < model_levels
  s0 = 1
  levels = model_levels
END IF   !   level_start == level_end

IF ( level_start == level_end ) THEN

  DO k =  start_level, end_level
    k0 = k
    IF( level_start < 1 ) k0 = 0

      num = 0
      DO counter = 1, tr_ukca
        num = num + 1
        CALL array_norm(                                             &
                        tracer_ukca(1-offx,1-offy,s0,counter),         &
                        row_length, rows, levels, k, k,              &
                        offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                        rim_n, rim_s, rim_e, rim_w, two_norm(num))
      END DO !  counter = 1, tr_ukca

      IF ( L_print_pe .OR. mype == 0 ) THEN
        num = 0
        DO counter = 1, tr_ukca
          num = num + 1
          WRITE(umMessage, '(A, I5, A, I2, A, E23.16)')               &
                ' Level ', k0, ' tracer_ukca' , counter,                & 
                ' Two_Norm =' , two_norm(num)
          CALL umPrint(umMessage, src='UKCA_TRACER_NORM')
        END DO !  counter = 1, tr_ukca
      END IF ! L_print_pe .OR. mype == 0

    END DO !  k = 0, model_levels

  ELSE   !  3d norms from  level_start to level_end

    num = 0
    DO counter = 1, tr_ukca
      num = num + 1
      CALL array_norm(                                            &
                      tracer_ukca(1-offx,1-offy,s0,counter),        &
                      row_length, rows, levels,                   &
                      level_start, level_end,                     &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END DO !  counter = 1, tr_ukca

    IF ( L_print_pe .OR. mype == 0 ) THEN
      num = 0
      DO counter = 1, tr_ukca
        num = num + 1
        WRITE(umMessage, '(A, I5, A, I5, A, I2, A, E23.16)')     &
                  ' Levels', start_level, ' to', level_end,      &
                  ' ukca_tracer' , counter, ' Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='UKCA_TRACER_NORM')
      END DO !  counter = 1, tr_ukca

    END IF ! L_print_pe .OR. mype == 0

  END IF ! level_start == level_end

END IF   ! tr_ukca > 0

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_tracer_norm

SUBROUTINE dms_norm(                                                 &
                    dms, level_start, level_end, rims,               &
                    L_do_rims, L_do_halos, L_print_pe )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod
USE array_norm_mod
USE pr_block4_mod
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels,        &
                                 global_row_length, global_rows
USE um_parcore, ONLY: mype
USE um_parvars, ONLY: offx, offy

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr, ONLY: umPrint, umMessage, umPrintFlush
!
! Description: Printing of norms for dms only
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DMS_NORM'

! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_do_rims
LOGICAL :: L_do_halos
LOGICAL :: L_print_pe   ! print diagnostics on all pe's if true

INTEGER, INTENT(IN) :: level_start, level_end, rims

REAL, INTENT(IN)  :: dms                                        &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

! Local Variables.
INTEGER :: k
INTEGER :: levels      !  model_levels + 1 
INTEGER :: start_level, end_level
INTEGER :: k0, s0      !  level zero pointers
INTEGER :: rim_n, rim_s, rim_e, rim_w

REAL ::  two_norm

! ----------------------------------------------------------------------
! 1.0 Start of subroutine code
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( L_do_rims ) THEN
  rim_n = 1
  rim_s = 1
  rim_e = 1
  rim_w = 1
ELSE ! L_do_rims
  rim_n = rims
  rim_s = rims
  rim_e = rims
  rim_w = rims
!  set up block norm if dom_w and dom_s /= 1
  IF ( dom_w /= 1 .AND. dom_s /= 1 ) THEN
    rim_s = dom_s
    rim_w = dom_w
    rim_n = global_rows - dom_no
    rim_e = global_row_length - dom_eo
  END IF ! dom_w /= 1 .AND. dom_s /= 1
END IF ! L_do_rims

levels = model_levels + 1
start_level = level_start
end_level = level_end
s0 = 0
IF ( level_start == level_end ) THEN
  IF ( level_start < 1 ) THEN
    start_level = 1
    end_level   = 1
  ELSE IF ( level_start == model_levels) THEN
    start_level = 1
    s0 = 1
    levels = model_levels
  END IF ! level_start == model_levels 
ELSE     !    0 < level_start < model_levels
  s0 = 1
  levels = model_levels
END IF   !   level_start == level_end

IF ( level_start == level_end ) THEN

  DO k =  start_level, end_level
    k0 = k
    IF( level_start < 1 ) k0 = 0

      CALL array_norm(                                                 &
                      dms(1-offx,1-offy,s0), row_length, rows,         &
                      levels, k, k,                                    &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,         &
                      rim_n, rim_s, rim_e, rim_w, two_norm )
      IF ( L_print_pe .OR. mype == 0 ) THEN
        WRITE(umMessage, '(A, I5, A, E23.16)')                         &
                    ' Level ', k0, ' dms Two_Norm =' , two_norm
        CALL umPrint(umMessage, src='DMS_NORM')
      END IF ! L_print_pe .OR. mype == 0

  END DO !  k = 0, model_levels

ELSE   !  3d norms from start_level to end_level  

    CALL array_norm(                                               &
                    dms(1-offx,1-offy,s0), row_length, rows,       &
                    levels, start_level, level_end,                &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm )
    IF ( L_print_pe .OR. mype == 0 ) THEN
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                ' Levels', start_level, ' to', level_end,        &
                ' dms Two_Norm =' , two_norm
      CALL umPrint(umMessage, src='DMS_NORM')
    END IF ! L_print_pe .OR. mype == 0

END IF ! level_start == level_end

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE dms_norm
END MODULE sl_tracer_norm_mod

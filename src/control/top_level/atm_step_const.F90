! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Initialisation of constants for Atm_Step
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

MODULE atm_step_const_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ATM_STEP_CONST_MOD'

CONTAINS

SUBROUTINE atm_step_consts(                                            &
                           row_length, rows,                           &
                           model_levels, bl_levels, first_flat_level )

USE atm_step_local, ONLY: tr_size, l_tracer,                           &
                          first_constant_r_rho_level,                  &
                          first_constant_r_rho_level_m1,               &
                          check_bottom_levels,                         &
                          solver_tolerance, couple_app,                & 
                          i_start, i_stop, j_start, j_stop,            &
                          j_begin, j_end, lambda_start,                &
                          lbc_size, ij_u, ij_v, L_do_rims

USE carbon_options_mod, ONLY: l_co2_interactive
USE run_aerosol_mod, ONLY: l_sulpc_so2, l_soot, l_ocff, l_biomass,  &
                           l_nitrate
USE murk_inputs_mod, ONLY: l_murk

USE dust_parameters_mod, ONLY: l_dust
USE rad_input_mod, ONLY: l_use_cariolle, lrad_ccrad

USE dump_headers_mod, ONLY: a_realhd, rh_z_top_theta,                    &
                            a_inthd, ih_1_c_rho_level

USE atm_fields_bounds_mod, ONLY: tdims

USE mphys_bypass_mod, ONLY: l_crystals, mphys_mod_top,                   &
                            qcf2_idims_start, qcf2_idims_end,            &
                            qcf2_jdims_start, qcf2_jdims_end,            &
                            qcf2_kdims_end
USE mphys_inputs_mod, ONLY: l_psd, l_mcr_qcf2
USE cv_run_mod,       ONLY: l_ccrad

USE dynamics_testing_mod, ONLY: L_Physics, L_idealised_data,            &
                                problem_number
USE problem_mod, ONLY: standard

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE UM_ParVars, ONLY: halo_i, halo_j,                                   &
                      nproc_x, nproc_y, datastart, at_extremity,        &
                      g_lasize, g_pe_index_EW, g_pe_index_NS
USE UM_ParParams, ONLY: halo_type_no_halo,                              &
                        PNorth, PEast, PSouth, PWest

USE UM_ParParams, ONLY: halo_type_no_halo,                              &
                       PNorth, PEast, PSouth, PWest
USE Field_Types, ONLY: fld_type_p

USE dynamics_input_mod, ONLY: L_new_tdisc,                              &
                              numcycles, innits, tol_sc_fact,           &
                              gcr_use_residual_tol,                     &
                              GCR_tol_res, GCR_tol_abs

USE missing_data_mod, ONLY: rmdi
USE sl_param_mod, ONLY: g_row_length, g_rows, g_i_pe, g_j_pe, flat_level
USE rot_coeff_mod, ONLY: lam_max_cfl
USE sl_input_mod, ONLY: high_order_scheme, L_mono, L_high,              &
                        interp_vertical_search_tol
USE eg_sisl_consts_mod, ONLY: eg_sisl_consts

USE UM_parcore, ONLY: mype, nproc
USE highos_mod,  ONLY: quinticLagrange

USE nlsizes_namelist_mod, ONLY: global_row_length, global_rows,         &
                                theta_off_size,                         &
                                tr_levels, tr_ukca, tr_vars

USE model_domain_mod, ONLY: l_regular, model_type, mt_global, mt_lam,   &
                            mt_cyclic_lam, mt_bi_cyclic_lam

USE umPrintMgr, ONLY: umPrint, umMessage, printstatus, prstatus_oper

IMPLICIT NONE

! Input variables.
INTEGER, INTENT(IN)  ::  row_length
INTEGER, INTENT(IN)  ::  rows
INTEGER, INTENT(IN)  ::  model_levels
INTEGER, INTENT(IN)  ::  bl_levels
INTEGER, INTENT(IN)  ::  first_flat_level

! Loop indices.
INTEGER  ::  i, j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATM_STEP_CONSTS'

! ------------------------------------------------------------------
! Set non-initialised variables to one
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! 1. ENDGame
! ----------------------------------------------------------------------
!     ENDGame dynamical core
l_new_tdisc = .TRUE.

! ----------------------------------------------------------------------
! 2. Set solver_tolerance
! ----------------------------------------------------------------------
IF ( gcr_use_residual_tol ) THEN
  solver_tolerance = gcr_tol_res * tol_sc_fact**(numcycles*innits-1)
ELSE
  solver_tolerance = gcr_tol_abs * tol_sc_fact**(numcycles*innits-1)
END IF

! ----------------------------------------------------------------------
! 3. Set constants for atm_step
! ----------------------------------------------------------------------
lbc_size       = 1
ij_v           = 1
ij_u           = 1

! Tracers
tr_size = (tr_levels * theta_off_size) ! size of a single tracer

L_Tracer = ( tr_vars > 0 .OR. L_soot .OR. L_CO2_interactive .OR.  &
             L_sulpc_SO2 .OR. L_murk .OR. L_dust .OR. L_biomass .OR.  &
             l_use_cariolle .OR. L_ocff .OR. L_nitrate .OR.           &
             tr_ukca > 0 )

! Set microphysics switches and constants not in modules via
! mphys_bypass_mod

l_crystals      = l_mcr_qcf2 .OR. (.NOT. l_psd)
mphys_mod_top   = a_realhd(rh_z_top_theta)

IF (l_mcr_qcf2) THEN

  qcf2_idims_start = tdims%i_start
  qcf2_idims_end   = tdims%i_end
  qcf2_jdims_start = tdims%j_start
  qcf2_jdims_end   = tdims%j_end
  qcf2_kdims_end   = tdims%k_end

ELSE ! l_mcr_qcf2

  qcf2_idims_start = 1
  qcf2_idims_end   = 1
  qcf2_jdims_start = 1
  qcf2_jdims_end   = 1
  qcf2_kdims_end   = 1

END IF ! l_mcr_qcf2

! Set radiation switches in rad_switches_mod

lrad_ccrad        = l_ccrad            ! cntlatm

! ----------------------------------------------------------------------
! 4. Set grid information
! ----------------------------------------------------------------------

first_constant_r_rho_level = a_inthd(ih_1_c_rho_level)

flat_level = first_constant_r_rho_level

IF (first_constant_r_rho_level == 1 ) THEN
  first_constant_r_rho_level_m1 = 1
ELSE
  first_constant_r_rho_level_m1 = first_constant_r_rho_level - 1
END IF

interp_vertical_search_tol = MIN(interp_vertical_search_tol,          &
                                 model_levels/2)

check_bottom_levels = MAX( bl_levels, interp_vertical_search_tol )
check_bottom_levels = MIN(check_bottom_levels, model_levels)

lambda_start = datastart(1) - halo_i

! ----------------------------------------------------------------------
! 5. Set mpp arrays from UM_ParVars information
! ----------------------------------------------------------------------

ALLOCATE ( g_row_length(0:nproc-1) )
ALLOCATE ( g_rows(0:nproc-1) )
DO i= 0,nproc-1
  g_row_length(i) = g_lasize(1,fld_type_p,halo_type_no_halo,i)
  g_rows      (i) = g_lasize(2,fld_type_p,halo_type_no_halo,i)
END DO ! i processors

!  i_start, i_stop, j_start, j_stop define the solution domain
!  j_begin, j_end exclude the poles
i_start = 1
i_stop = row_length
j_start = 1
j_stop = rows
j_begin = 1
j_end = rows

! ----------------------------------------------------------------------
! 6. Set logical switch for diagnostic printing of l2norms
! ----------------------------------------------------------------------

L_do_rims = .FALSE.
!  mt_bi_cyclic_lam keeps defaults; other domains change as below
IF (model_type == mt_global) THEN
  L_do_rims = .TRUE.
END IF  ! model_type  ==  mt_global

IF (model_type == mt_lam) THEN
  IF (at_extremity(PWest)) i_start = 2
  IF (at_extremity(PEast)) i_stop = row_length - 1
END IF  ! model_type  ==  mt_lam

! ----------------------------------------------------------------------
! 7. Check large halos are sufficient for chosen interpolation 
! ----------------------------------------------------------------------

DO i = 1, 3
  IF ( high_order_scheme(i) == quinticLagrange ) THEN
    j = 3
    IF ( halo_j  <  5 .AND. mype == 0 ) THEN
      WRITE(umMessage,'(A, I4, A)') ' ** ERROR ** halo_j = ', halo_j    &
                    , ' is too small for quintic interpolation.'
      CALL umPrint(umMessage,src='ATM_STEP_CONSTS')
    END IF
  ELSE  ! cubic interpolation
    j = 2
    IF ( halo_j  <  4 .AND. mype == 0 ) THEN
      WRITE(umMessage,'(A, I4, A)') ' ** ERROR ** halo_j = ', halo_j    &
                    , ' is too small for cubic interpolation.'
      CALL umPrint(umMessage,src='ATM_STEP_CONSTS')
    END IF
  END IF
END DO

! ----------------------------------------------------------------------
! 8.  Calculate max cfl possible in LAM
!     Y-direction also applies in global model 
!     though variable is not used in code.
! ----------------------------------------------------------------------

lam_max_cfl(1) = halo_i - j
lam_max_cfl(2) = halo_j - j

! ----------------------------------------------------------------------
! 9. Error trapping for advection choices
! ----------------------------------------------------------------------

DO i=1,3
  IF (.NOT. L_mono(i) .AND. .NOT. L_high(i)) THEN
    IF (i == 1) THEN
      WRITE(umMessage,'(A)')                            &
      'WARNING Advection choices incompatible for theta'
      CALL umPrint(umMessage,src='ATM_STEP_CONSTS')
    END IF
    IF (i == 2) THEN
      WRITE(umMessage,'(A)')                            &
      'WARNING Advection choices incompatible for moisture'
      CALL umPrint(umMessage,src='ATM_STEP_CONSTS')
    END IF
    IF (i == 3) THEN
      WRITE(umMessage,'(A)')                            &
      'WARNING Advection choices incompatible for winds'
      CALL umPrint(umMessage,src='ATM_STEP_CONSTS')
    END IF
    WRITE(umMessage,'(A)') 'Both L_high and L_mono switches set to false '
    CALL umPrint(umMessage,src='ATM_STEP_CONSTS')

  END IF
END DO

! ----------------------------------------------------------------------
! 10. Set mpp arrays
! ----------------------------------------------------------------------

ALLOCATE ( g_i_pe(1-halo_i:global_row_length+halo_i) )
ALLOCATE ( g_j_pe(1-halo_j:global_rows+halo_j) )

DO i = 1-halo_i, 0
  g_i_pe(i) = 0
END DO ! i

DO i = 1,global_row_length
  g_i_pe(i) = g_pe_index_EW(i)
END DO ! i row points

DO i = global_row_length+1, global_row_length+halo_i
  g_i_pe(i) = nproc_x - 1
END DO ! i

DO j = 1-halo_j, 0
  g_j_pe(j) = 0
END DO ! j

DO j = 1, global_rows
  g_j_pe(j) = g_pe_index_NS(j)
END DO ! j col points

DO j = global_rows+1, global_rows+halo_j
  g_j_pe(j) = nproc_y - 1
END DO ! j

IF ( Problem_number /= standard ) THEN
  IF (.NOT. L_Physics .AND. mype  ==  0) THEN
    CALL umPrint('If source of dump is a full model run with orography',&
        src='ATM_STEP_CONSTS')
    CALL umPrint( 'then there is no guarantee that run will work.',     &
        src='ATM_STEP_CONSTS')
    CALL umPrint( 'since physics is OFF',src='ATM_STEP_CONSTS')
  END IF !.not. L_Physics.and. mype  ==  0
END IF ! Problem_number /= standard

! ----------------------------------------------------------------------
!  11. Set SISL consts 
! ----------------------------------------------------------------------

CALL eg_sisl_consts()

! ----------------------------------------------------------------------
!  12. Message if physics coupling is disabled
! ----------------------------------------------------------------------

IF (printstatus > prstatus_oper) THEN
  IF (.NOT. couple_app) CALL umPrint('coupling disabled!',               &
                                      src='ATM_STEP_CONSTS')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE atm_step_consts

END MODULE atm_step_const_mod


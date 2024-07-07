! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Initialisation section for Atm_Step
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE atm_step_init (           &
lambda_a, g_row_length, g_rows, g_i_pe, g_j_pe, flux_e, flux_h, ustar_in, &
z0m_scm, z0h_scm, errorstatus )

USE atm_step_local

USE carbon_options_mod, ONLY: l_co2_interactive
USE run_aerosol_mod, ONLY: l_sulpc_so2, l_soot, l_ocff, l_biomass,  &
                           l_nitrate
USE murk_inputs_mod, ONLY: l_murk

USE dust_parameters_mod, ONLY: l_dust
USE rad_input_mod
USE diag_print_mod

USE mphys_bypass_mod, ONLY: l_crystals, mphys_mod_top,                   &
                            qcf2_idims_start, qcf2_idims_end,            &
                            qcf2_jdims_start, qcf2_jdims_end,            &
                            qcf2_kdims_end
USE mphys_inputs_mod, ONLY: l_psd, l_mcr_qcf2
USE cv_run_mod,       ONLY: l_ccrad, l_pc2_diag_sh

USE dyn_var_res_mod

USE timestep_mod, ONLY: timestep_number
USE turb_diff_mod, ONLY: print_step, l_diag_l2helm, l_diag_l2norms,      &
                         first_norm_print
USE dump_headers_mod, ONLY: ih_1_c_rho_level, rh_z_top_theta, a_inthd,   &
                            a_realhd
USE atm_fields_bounds_mod, ONLY: tdims
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE UM_ParVars
USE UM_ParCore,   ONLY: nproc
USE UM_ParParams, ONLY: halo_type_no_halo, peast, pwest
USE Field_Types,  ONLY: fld_type_p
USE Control_Max_Sizes

USE sl_input_mod, ONLY: interp_vertical_search_tol

USE dynamics_testing_mod, ONLY: l_idealised_data
USE idealise_run_mod, ONLY: l_spec_z0
USE bl_option_mod, ONLY: flux_bc_opt, interactive_fluxes
USE cderived_mod, ONLY: delta_lambda, delta_phi, base_lambda
USE horiz_grid_mod, ONLY: cartesian_grid
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE missing_data_mod, ONLY: rmdi


USE nlsizes_namelist_mod, ONLY:                                            &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,        &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,         &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst,     &
    a_len_inthd, a_len_realhd, bl_levels, global_row_length, global_rows,  &
    len1_lookup, len_dumphist, len_fixhd,                                  &
    mpp_len1_lookup, row_length, rows, theta_off_size, tr_levels, tr_ukca, &
    tr_vars

USE model_domain_mod, ONLY: l_regular, model_type, mt_global, mt_lam           &
                          , mt_cyclic_lam, mt_bi_cyclic_lam

USE pws_diags_mod, ONLY: pws_diags_alloc

IMPLICIT NONE

! Subroutine arguments
INTEGER :: g_row_length(0:nproc-1)       ! Table of number of points on a row
INTEGER :: g_rows(0:nproc-1)             ! Table number of rows in theta field
INTEGER :: g_i_pe(1-halo_i:global_row_length+halo_i)
                    ! proc on my proc-row holding a given value in i direction
INTEGER :: g_j_pe(1-halo_j:global_rows+halo_j)
                    ! proc on my proc-col holding a given value in j direction

REAL :: flux_e(row_length, rows)
REAL :: flux_h(row_length, rows)
REAL :: ustar_in(row_length, rows)
REAL :: z0m_scm(row_length, rows)
REAL :: z0h_scm(row_length, rows)

REAL :: lambda_a (row_length) ! delta_lambda term for polar wind

! Local variables
INTEGER :: errorstatus
CHARACTER (LEN=errormessagelength) :: cmessage      ! used for ereport

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATM_STEP_INIT'

! --------------------------------------------------------------------------
! Set non-initialised variables to one
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
lbc_size       = 1
ij_v           = 1
ij_u           = 1

! DEPENDS ON: atm_step_timestep_init
CALL atm_step_timestep_init

IF ((model_type == mt_global) .AND. (cartesian_grid)) THEN
  ErrorStatus=123

  CALL Ereport("ATM_STEP", ErrorStatus,                           &
             "Unable to run global model with a Cartesian grid.")
END IF

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

! set up logical switch for diagnostic printing of l2norms
L_print_L2norms = .FALSE.
L_print_L2helm = .FALSE.
itemp = timestep_number - first_norm_print
L_do_rims = .FALSE.
IF ( model_type == mt_global ) L_do_rims = .TRUE.

IF ( print_step > 0 .AND. itemp >=0) THEN
  IF ( MOD( itemp , print_step ) == 0 ) THEN
    IF ( L_diag_L2norms ) THEN
      IF ( itemp / print_step  > 11 ) THEN
        WRITE(umMessage,*)'l2norms printing too often, limited to 12 occasions'
        CALL umPrint(umMessage,src='atm_step_init')
      ELSE
        L_print_L2norms = .TRUE.
      END IF  ! itemp / print_step > 11
      IF ( itemp / print_step  < 4 ) THEN
        ! Diagnostic print limited to 4 occasions
        L_print_diag = L_diag_L2norms
      END IF  ! itemp / print_step < 4
    END IF ! L_diag_L2norms
    IF (  L_diag_L2helm ) THEN
      IF ( itemp / print_step  > 11 ) THEN
        WRITE(umMessage,*)'l2norms printing too often, limited to 12 occasions'
        CALL umPrint(umMessage,src='atm_step_init')
      ELSE
        L_print_L2helm = .TRUE.
      END IF  ! itemp / print_step > 11
    END IF !  L_diag_L2helm
  END IF !  mod( itemp , print_step ) == 0)
END IF ! print_step > 0

! grid information
first_constant_r_rho_level = a_inthd(ih_1_c_rho_level)

IF (first_constant_r_rho_level  ==  1 ) THEN
  first_constant_r_rho_level_m1 = first_constant_r_rho_level
ELSE
  first_constant_r_rho_level_m1 = first_constant_r_rho_level - 1
END IF

gi = 1

lambda_a(:) = rmdi

lambda_start = datastart(1) - halo_i

! Set mpp arrays from UM_ParVars information
DO i= 0,nproc-1
  g_row_length(i) = g_lasize(1,fld_type_p,halo_type_no_halo,i)
  g_rows      (i) = g_lasize(2,fld_type_p,halo_type_no_halo,i)
END DO ! i processors

DO i= 1-halo_i, 0
  g_i_pe(i) = 0
END DO ! i

DO i= 1,global_row_length
  g_i_pe(i) = g_pe_index_EW(i)
END DO ! i row points

DO i= global_row_length+1, global_row_length+halo_i
  g_i_pe(i) = nproc_x-1
END DO ! i

DO j= 1-halo_j, 0
  g_j_pe(j) = 0
END DO ! j

DO j= 1,global_rows
  g_j_pe(j) = g_pe_index_NS(j)
END DO ! j row points

DO j= global_rows+1, global_rows+halo_j
  g_j_pe(j) = nproc_y-1
END DO ! j


! End Set mpp arrays

! set number of levels to check for trajectory inside orography at
! bottom of model, the same parameter is also used inside the
! interpolation routine.

check_bottom_levels = MAX(bl_levels, interp_vertical_search_tol )

! Define number of work arrays needed in sl_full_wind routine depending
! on model type choice.
IF (model_type == mt_global .OR.                              &
    model_type == mt_cyclic_lam .OR.                          &
    model_type == mt_bi_cyclic_lam) THEN
  n_Y_arrays    = 1
  n_Yw_arrays   = 1
  n_Yd_arrays   = 1
  n_Ydw_arrays  = 1
ELSE IF (model_type == mt_lam) THEN
  n_Y_arrays    = 3
  n_Yw_arrays   = 2
  n_Yd_arrays   = 3
  n_Ydw_arrays  = 2
END IF

!  i_start, i_stop, j_start, j_stop define the solution domain
!  j_begin, j_end exclude the poles
i_start = 1
i_stop = row_length
j_start = 1
j_stop = rows
j_begin = 1
j_end = rows

!  mt_bi_cyclic_lam keeps defaults; other domains change as below

IF (model_type == mt_lam) THEN
  IF (at_extremity(PWest)) i_start = 2
  IF (at_extremity(PEast)) i_stop = row_length - 1
END IF  ! model_type  ==  mt_lam

! Check settings if not using idealised forcing
IF (.NOT. L_idealised_data) THEN
  !  Only interactive_fluxes currently works for standard UM
  IF (flux_bc_opt /= interactive_fluxes) THEN 
    ErrorStatus=1
    WRITE(cmessage,'(A)')                                               &
       'The chosen flux_bc_opt option only works with l_idealised_data'
    CALL Ereport(RoutineName, ErrorStatus, cmessage)
  END IF
  ! Set variables for roughness length to false/0.0 if not using forcing
  L_spec_z0     = .FALSE.
  z0m_scm(:,:) = 0.0
  z0h_scm(:,:) = 0.0
END IF
!  Set dummy variables for SCM Diagnostics to false for full UM
L_SCMDiags(1:nSCMDpkgs)  = .FALSE.

! Allocate PWS diagnostic arrays if required
CALL pws_diags_alloc( &
          errorstatus)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE atm_step_init


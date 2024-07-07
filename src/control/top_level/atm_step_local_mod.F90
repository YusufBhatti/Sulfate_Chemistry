! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Declares variables and allocatable arrays which are local to atm_step
!  (and below)
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

MODULE atm_step_local

USE def_easyaerosol, ONLY: t_easyaerosol_rad, t_easyaerosol_cdnc

IMPLICIT NONE

LOGICAL :: l_do_rims

INTEGER :: cycleno ! Number of cycles (iterations) for iterative SISL

INTEGER :: tr_size ! size of a single tracer - FLUME variable

LOGICAL :: l_tracer ! T if *any* tracer variables present

LOGICAL :: l_print_l2norms ! diagnostic printing of l2norms

INTEGER :: first_constant_r_rho_level    ! 1st rho level with r constant
INTEGER :: first_constant_r_rho_level_m1 ! Max (1,first_constant_r_rho_level)

INTEGER :: i_start
INTEGER :: i_stop
INTEGER :: j_start
INTEGER :: j_stop
INTEGER :: j_begin
INTEGER :: j_end

INTEGER :: lambda_start ! pointer for start of lambda_p/lambda_u on this pe

! used in interpolation code - no. of levels to check for whether the departure
! point lies inside the orography.
INTEGER :: check_bottom_levels

INTEGER :: rhc_row_length
INTEGER :: rhc_rows

!  Dummy variables for SCM Diagnostics,
!  for passing down to atmos_physics1 and 2 routines
INTEGER, PARAMETER :: nscmdpkgs = 15
LOGICAL :: l_scmdiags(nscmdpkgs)
LOGICAL :: l_flux_bc    ! T if prescribed surface fluxes to be used

LOGICAL, PARAMETER :: set_halos = .TRUE.
LOGICAL, PARAMETER :: do_not_set_halos = .FALSE.
                   ! logicals used to set up the calls to tr_set_phys.

LOGICAL, PARAMETER :: couple_app=.TRUE.   ! this parameter can be used
                   ! to force a dynamics only run even with the full physics
                   ! suite executing. For testing purposes. If false it will
                   ! zero all source terms coming from atmos_physics1 &
                   ! atmos_physics2 (as far as dynamics is concerned - mostly).
                   ! It could be set to l_physics, however that would exclude
                   ! the eg_idl_forcing contribution which is enabled by
                   ! disabling physics in the namelist and setting the
                   ! appropriate IDEALISED namelist logical.

! time-stepping weights. Values may be different
! at different cycles namelist controlled
REAL :: alpha1, alpha2, alpha3, alpha4

REAL :: max_q_star, min_q_star    ! local diagnostics

REAL :: increment_factor ! For calculating value of LBC at next TL
REAL :: solver_tolerance

! Variables required for call to SET_LATERAL_BOUNDARIES
INTEGER :: lbc_size        ! size of a single level of LBC
LOGICAL :: l_do_halos      ! update the halos?
LOGICAL :: l_do_boundaries ! update the boundaries?

! -------------------------------------------------

INTEGER :: nd_o3  ! Total size of ozone array supplied

! Code to do with tropopause diagnostics from O3_to_3D
! Declare logicals for Stash Diagnostics used in call to O3_to_3D
LOGICAL :: l_o3_trop_level   ! stash code 2,280
LOGICAL :: l_o3_trop_height  ! stash code 2,281
LOGICAL :: l_t_trop_level    ! stash code 2,282
LOGICAL :: l_t_trop_height   ! stash code 2,283

INTEGER :: land_pts_trif ! For dimensioning variables in NI_bl_ctl
INTEGER :: npft_trif     ! Depending on whether TRIFFID is in use

REAL, POINTER :: cs(:,:)  => NULL() ! soil carb content & accum soil respiration
REAL, POINTER :: rsa(:,:) => NULL() ! (target varies according to l_TRIFFID)

INTEGER :: dim_cs1=1 ! soil C dimension: 1 for single, 4 for RothC
INTEGER :: dim_cs2=1 ! soil C dimension: 1 for single, LAND_FIELD for RothC

INTEGER :: co2_dim_len ! For dimension 3-D CO2 field to be passed
INTEGER :: co2_dim_row !     to NI_bl_ctl
INTEGER :: co2_dim_lev !     and NI_rad_ctl

! Array dimensions for sea-salt aerosol
INTEGER :: salt_dim1
INTEGER :: salt_dim2
INTEGER :: salt_dim3

! Array dimensions for Aero_Ctl
INTEGER :: aero_dim1
INTEGER :: aero_dim2
INTEGER :: aero_dim3

! Declare the tropopause variables output for O3_to_3D as allocatable
REAL, ALLOCATABLE:: o3_trop_level(:,:)
REAL, ALLOCATABLE:: t_trop_level(:,:)
REAL, ALLOCATABLE:: o3_trop_height(:,:)
REAL, ALLOCATABLE:: t_trop_height(:,:)

!3d ozone (expanded from zonal) for radiation
REAL, ALLOCATABLE:: ozone3D(:,:,:)

! Local Arrays to store microphysics fields
REAL, ALLOCATABLE, TARGET :: qcf2_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: qrain_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: qgraup_star(:,:,:)

! Local arrays for when using mixing ratios
REAL, ALLOCATABLE :: mix_v(:,:,:), mix_cl(:,:,:), mix_cf(:,:,:)
REAL, ALLOCATABLE :: mix_cf2(:,:,:), mix_rain(:,:,:), mix_graup(:,:,:)

!    Local arrays for when using mixing ratios
REAL, ALLOCATABLE, TARGET :: mix_v_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: mix_cl_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: mix_cf_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: mix_cf2_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: mix_rain_star(:,:,:)
REAL, ALLOCATABLE, TARGET :: mix_graup_star(:,:,:)

!  stashworki = stashwork for section i
REAL, ALLOCATABLE:: stashwork1(:),stashwork2(:),stashwork3(:), &
                    stashwork4(:),stashwork5(:),stashwork6(:), &
                    stashwork8(:),stashwork9(:),stashwork12(:),&
                    stashwork13(:),stashwork14(:),stashwork17(:), &
                    stashwork19(:),stashwork26(:),stashwork30(:), &
                    stashwork10(:),stashwork18(:),stashwork35(:), &
                    stashwork39(:),stashwork21(:),stashwork20(:), &
                    stashwork53(:),stashwork54(:)

! Local dynamic arrays for phys1 and phys2 increments for moisture:
REAL, ALLOCATABLE :: cf_phys1(:,:,:), cfl_phys1(:,:,:), cff_phys1(:,:,:)

! local dynamic arrays for PC2

REAL, ALLOCATABLE :: t_inc_pres(:,:,:)
      ! Temperature increment due to adiabatic temperature forcing (K)
REAL, ALLOCATABLE :: q_inc_pres(:,:,:)
      ! Vapour increment due to adiabatic temperature forcing (kg kg-1)
REAL, ALLOCATABLE :: qcl_inc_pres(:,:,:)
      ! Liquid increment due to adiabatic temperature forcing (kg kg-1)
REAL, ALLOCATABLE :: qcf_inc_pres(:,:,:)
      ! Ice increment due to adiabatic temperature forcing (kg kg-1)
REAL, ALLOCATABLE :: cf_inc_pres(:,:,:)
      ! Total cloud fraction increment due to adiabatic forcing (no units)
REAL, ALLOCATABLE :: cfl_inc_pres(:,:,:)
      ! Liquid cloud fraction increment due to adiabatic forcing (no units)
REAL, ALLOCATABLE :: cff_inc_pres(:,:,:)
      ! Ice cloud fraction increment due to adiabatic forcing (no units)
REAL, ALLOCATABLE :: t_dini(:,:,:)
      ! Temperature increment due to initiation (K)
REAL, ALLOCATABLE :: q_dini(:,:,:)
      ! Vapour increment due to initiation (kg kg-1)
REAL, ALLOCATABLE :: qcl_dini(:,:,:)
      ! Liquid increment due to initiation (kg kg-1)
REAL, ALLOCATABLE :: qcf_dini(:,:,:)
      ! Ice increment due to initiation (kg kg-1)
REAL, ALLOCATABLE :: cf_dini(:,:,:)
      ! Total cloud fraction increment due to initiation (no units)
REAL, ALLOCATABLE :: cfl_dini(:,:,:)
      ! Liquid cloud fraction increment due to initiation(no units)
REAL, ALLOCATABLE :: cff_dini(:,:,:)
      ! Ice cloud fraction increment due to initiation (no units)

REAL, ALLOCATABLE :: rhts(:,:,:), qtts(:,:,:), tlts(:,:,:), ptts(:,:,:)

! Extra variables needed for cycling
! Vars ending in _phys1 are copies holding the value the original
! variable had after exiting phys1.
! Vars ending in _np1 are tn+1 estimates holding the value the original
! variable had at the end of the last cycle (provided that CycleNo>1).
! obtained from the last
! cycle when CycleNo>1.

REAL, ALLOCATABLE :: r_u_phys1(:,:,:), r_v_phys1(:,:,:)
REAL, ALLOCATABLE :: thetastar_phys1(:,:,:), qstar_phys1(:,:,:)
REAL, ALLOCATABLE :: qclstar_phys1(:,:,:), qcfstar_phys1(:,:,:)
REAL, ALLOCATABLE :: qcf2_star_phys1(:,:,:), qrain_star_phys1(:,:,:)
REAL, ALLOCATABLE :: qgraup_star_phys1(:,:,:), ti_phys1(:,:,:)
REAL, ALLOCATABLE :: cca_phys1(:,:,:), area_cld_frac_phys1(:,:,:)
REAL, ALLOCATABLE :: e_trb_phys1(:,:,:), tsq_trb_phys1(:,:,:),     &
                     qsq_trb_phys1(:,:,:), cov_trb_phys1(:,:,:)
REAL, ALLOCATABLE :: bulk_cld_frac_phys1(:,:,:), bulk_cld_liq_phys1(:,:,:)
REAL, ALLOCATABLE :: bulk_cld_fr_phys1(:,:,:)

REAL, SAVE, ALLOCATABLE, TARGET :: wetrho_r_sq_np1(:,:,:)

REAL, SAVE, ALLOCATABLE, TARGET :: u_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: v_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: w_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: dryrho_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: rho_np1(:,:,:)

REAL, SAVE, ALLOCATABLE, TARGET :: etadot_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: thetav_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: exner_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: exner_surf_np1(:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: m_v_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: m_cl_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: m_cf_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: m_r_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: m_gr_np1(:,:,:)
REAL, SAVE, ALLOCATABLE, TARGET :: m_cf2_np1(:,:,:)


REAL, ALLOCATABLE:: z0msea_phys1(:,:), zh_phys1(:,:), ddmfx_phys1(:,:), &
                                     ti_gb_phys1(:,:)
REAL, ALLOCATABLE:: deep_flag_phys1(:,:), past_precip_phys1(:,:),&
                                     past_conv_ht_phys1(:,:)
REAL, ALLOCATABLE:: conv_prog_1_phys1(:,:,:)
REAL, ALLOCATABLE:: conv_prog_2_phys1(:,:,:)
REAL, ALLOCATABLE:: conv_prog_3_phys1(:,:,:)
REAL, ALLOCATABLE:: conv_prog_precip_phys1(:,:,:)
REAL, ALLOCATABLE:: t_land_ctile_phys1(:,:)
REAL, ALLOCATABLE:: t_sea_ctile_phys1(:,:)
REAL, ALLOCATABLE:: t_sice_ctile_phys1(:,:,:)
REAL, ALLOCATABLE:: t_surf_phys1(:,:), t_sf_tile_phys1(:,:)
REAL, ALLOCATABLE:: snow_tile_phys1(:,:), dolr_phys1(:,:)
REAL, ALLOCATABLE:: rho_lbc_real_tend(:,:)

INTEGER, ALLOCATABLE :: ccb_phys1(:,:), cct_phys1(:,:)

! Variables required for ice category selection
REAL, POINTER :: p_ti(:,:,:)        => NULL()
REAL, POINTER :: p_ice_fract(:,:,:) => NULL()
REAL, POINTER :: p_ice_thick(:,:,:) => NULL()
REAL, POINTER :: p_tstar_sice(:,:,:) => NULL()
REAL, POINTER :: p_snodep_sice(:,:,:) => NULL()
REAL, POINTER :: p_ice_thick_rad(:,:,:) => NULL()
REAL, POINTER :: p_ice_fract_rad(:,:,:) => NULL()

! increment diagnostics:
REAL, ALLOCATABLE :: u_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: v_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: t_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: q_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: qcl_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: qcf_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: qrain_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: qgraup_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: qcf2_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: cf_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: cfl_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: cff_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: w_incr_diagnostic(:,:,:)

! Allocatable arrays for use in AC_CTL call
REAL, ALLOCATABLE :: work_q(:,:,:), work_qcl(:,:,:), work_qcf(:,:,:)

! Workspace defined as allocatable arrays, since they each communicate
! fields between near adjacent calls and only need to use memory for
! a subset of the total routine.
REAL, ALLOCATABLE :: exner_prime(:,:,:) ! soln to helmholtz solver
REAL, ALLOCATABLE :: dtheta_dr_term(:,:,:)
REAL, ALLOCATABLE :: depart_lambda(:,:,:), depart_phi(:,:,:)
REAL, ALLOCATABLE :: depart_r_theta(:,:,:), depart_r_w(:,:,:)

REAL, ALLOCATABLE :: rhcpt(:,:,:)

! ----------------------------------------------------------------

LOGICAL ::  gather   ! Convert to sea_ice points within BL code
PARAMETER ( gather = .TRUE. ) ! (was l_compress_seaice)

INTEGER :: i, j, k, l   ! loop counters
INTEGER :: ij           ! 2D array index for offx/y variables
INTEGER :: lbc_pt       ! pointer to LBC array
INTEGER :: ij_u         ! 2D array index for U offx/y variables
INTEGER :: ij_v         ! 2D array index for V offx/y variables
INTEGER :: ji           ! 2D array index for halo_i/j variables
INTEGER :: icount       ! diagnostic count
INTEGER :: j_ti, j_ice_fract, j_ice_thick

INTEGER :: n_y_arrays    ! = 1 for global, 3 for LAM
INTEGER :: n_yw_arrays   ! = 1 for global, 2 for LAM
INTEGER :: n_yd_arrays   ! = 1 for global, 3 for LAM
INTEGER :: n_ydw_arrays  ! = 1 for global, 2 for LAM
!
INTEGER :: i_field       ! fields in multivariate swapbounds

INTEGER :: info          ! icode return from umFlush

REAL  ::   constant
REAL  ::   h_print     ! needed for printing idealised orography

LOGICAL :: l_update_lbcs, l_apply_lbcs, l_balance, gcr_zero_guess_it

INTEGER :: itemp
INTEGER :: gi

! Oxidant mass-mixing ratios and concentrations, for use in sulphur
! cycle.
REAL, ALLOCATABLE :: o3_mmr(:,:,:), hno3_mmr(:,:,:), h2o2_mmr(:,:,:)
REAL, ALLOCATABLE :: oh_conc(:,:,:), ho2_conc(:,:,:)

LOGICAL :: l_physics_store

! Local variables for using the aerosol climatology for NWP

! Internal array of mass-mixing ratios
REAL, ALLOCATABLE :: arcl(:,:,:,:)

! Number of requested species within the climatology
INTEGER :: n_arcl_species

! Corresponding number of requested components
INTEGER :: n_arcl_compnts

! Local structures for using EasyAerosol climatologies
TYPE (t_easyaerosol_rad)  :: easyaerosol_sw
TYPE (t_easyaerosol_rad)  :: easyaerosol_lw
TYPE (t_easyaerosol_cdnc) :: easyaerosol_cdnc

! Declare allocatable arrays for passing cloud fractions
! to LS_ACF_Brooks
REAL, ALLOCATABLE:: &
      cf_bulk_nohalo(:,:,:), cf_liquid_nohalo(:,:,:), cf_frozen_nohalo(:,:,:)


LOGICAL, SAVE :: first_atmstep_call = .TRUE.
LOGICAL, SAVE :: iau_in_initial     = .FALSE. ! Flag to test within atm_step if IAU ran in initial

END MODULE atm_step_local


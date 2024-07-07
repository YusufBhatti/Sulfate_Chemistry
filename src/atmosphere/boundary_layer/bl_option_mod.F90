! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

!  Data module for switches/options concerned with the BL scheme.
! Description:
!   Module containing runtime options/data used by the boundary
!   layer scheme

! Method:
!   Switches and associated data values used by the boundary layer
!   scheme are defined here and assigned default values. These may
!   be overridden by namelist input.

!   Any routine wishing to use these options may do so with the 'USE'
!   statement.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Boundary Layer

! Code Description:
!   Language: FORTRAN 95

MODULE bl_option_mod

USE missing_data_mod, ONLY: rmdi, imdi
USE visbty_constants_mod, ONLY: calc_prob_of_vis
USE mym_option_mod, ONLY: bdy_tke, tke_levels, l_local_above_tkelvs,    &
      l_my_initialize, l_my_ini_zero, my_ini_dbdz_min,                  &
      l_adv_turb_field, l_my_condense, l_shcu_buoy, shcu_levels,        &
      wb_ng_max, my_lowest_pd_surf, l_my_prod_adj,                      &
      my_z_limit_elb, l_print_max_tke, tke_cm_mx, tke_cm_fa, tke_dlen
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength
USE chk_opts_mod, ONLY: chk_var, def_src

IMPLICIT NONE

!=======================================================================
! BL namelist options - ordered as in meta-data sort-key, with namelist
! option and possible values
!=======================================================================

! 01 Switch to determine version of boundary layer scheme
INTEGER :: i_bl_vn = imdi
! 0 => no boundary layer scheme
INTEGER, PARAMETER :: i_bl_vn_0  = 0
! 1 => 1A prognostic TKE based turbulent closure
INTEGER, PARAMETER :: i_bl_vn_1a = 1
! 2 => 9B non-local scheme with revised K-prof depths
INTEGER, PARAMETER :: i_bl_vn_9b = 2
! 3 => 9C revised entr fluxes plus new scalar flux-grad
INTEGER, PARAMETER :: i_bl_vn_9c = 3

! 02 Stable boundary layer stability function option
INTEGER :: sbl_op= imdi
! Options for stable boundary layers
INTEGER, PARAMETER ::  long_tails           = 0
INTEGER, PARAMETER ::  sharpest             = 1
INTEGER, PARAMETER ::  sharp_sea_long_land  = 2
INTEGER, PARAMETER ::  mes_tails            = 3
INTEGER, PARAMETER ::  louis_tails          = 4
INTEGER, PARAMETER ::  depth_based          = 5
INTEGER, PARAMETER ::  sharp_sea_mes_land   = 6
INTEGER, PARAMETER ::  lem_stability        = 7
INTEGER, PARAMETER ::  sharp_sea_louis_land = 8
INTEGER, PARAMETER ::  equilibrium_sbl      = 9

! 03 Sets transition Ri for Sharpest function
REAL :: ritrans = rmdi  ! 1/g0 so 0.1 for original Sharpest
                        ! (closer to 0.2 makes it sharper)

! 04 Options for convective BL stability functions
INTEGER :: cbl_op = imdi
! UM_std   (=0) => Standard UM
INTEGER, PARAMETER ::  um_std     = 0
! neut_cbl(=1) => Keep neutral stability for Ri<0
INTEGER, PARAMETER ::  neut_cbl   = 1
! LEM_conv (=2) => Conventional LEM
INTEGER, PARAMETER ::  lem_conven = 2
! LEM_std (=3) => Standard LEM
INTEGER, PARAMETER ::  lem_std    = 3

! 05 Minimum value of the mixing length (m)
REAL :: lambda_min_nml = rmdi

! 06 Switch to mix qcf in the vertical
LOGICAL :: l_bl_mix_qcf = .FALSE.

! 07 Switch for free atmospheric mixing options
INTEGER :: local_fa= imdi
! to_sharp_across_1km (=1) => smoothly switch to sharpest across 1km AGL
INTEGER, PARAMETER :: to_sharp_across_1km = 1
! ntml_level_corrn (=2) => correct level for ntml_local (downwards)
INTEGER, PARAMETER :: ntml_level_corrn    = 2
! free_trop_layers (=3) => as "ntml_level_corrn" but also diagnose
! FA turbulent layer depths
INTEGER, PARAMETER :: free_trop_layers    = 3

! 08 Switch to keep local mixing in free atmosphere
INTEGER :: Keep_Ri_FA = imdi
! Only reduce local K to zero at discontinuous inversions
INTEGER, PARAMETER :: except_disc_inv = 2

! 09 SBL mixing dependent on sg orography
INTEGER :: sg_orog_mixing = imdi
!   1 - extending length of SHARPEST tail following McCabe and Brown (2007)
INTEGER, PARAMETER :: extended_tail = 1
!   2 - include subgrid drainage shear in Ri and diffusion coefficients,
INTEGER, PARAMETER :: sg_shear = 2
!   3 - as 2 + orographically enhanced lambda_m (note lambda_h enhancement
INTEGER, PARAMETER :: sg_shear_enh_lambda = 3
! not included in bdy_expl2 in error but now operational in UKV) and smooth
! decrease to lambda_min above

! 10 Switch to apply heating source from turbulence dissipation
INTEGER :: fric_heating = imdi
! OFF (=0) => not used
! ON  (=1) => used

! 11 Options for surface fluxes solver (explicit vs implicit)
INTEGER :: flux_bc_opt = imdi
! Options for surface fluxes solver (explicit vs implicit)
INTEGER, PARAMETER :: interactive_fluxes  = 0
! use land surface model to calculate fluxes interactively with PBL
INTEGER, PARAMETER :: specified_fluxes_only  = 1
! only specify sensible and latent surface fluxes, flux_h and flux_e
INTEGER, PARAMETER :: specified_fluxes_tstar = 2
! as (1) plus specify surface T, tstar (aka t_surf in atmos_phys2)
INTEGER, PARAMETER :: specified_fluxes_cd    = 3
! as (2) plus use specified ustar to give CD (still implicit solver)

!=======================================================================
! Options specific to 9B/9C (non-local) schemes
!=======================================================================

! 01 Height of number of levels for non_local scheme
REAL :: z_nl_bl_levels= rmdi

! 02 Switch to use zi/L in the diagnosis of shear-driven boundary layers:
INTEGER :: idyndiag = imdi
! OFF (=0) => not used
! DynDiag_ZL (=1) => over sea uses zi/L but gives problems when
! BL_LEVELS is high
INTEGER, PARAMETER :: DynDiag_ZL = 1
! DynDiag_ZL_corrn (=2) => as 1 but copes with high BL_LEVELS
INTEGER, PARAMETER :: DynDiag_ZL_corrn = 2
! DynDiag_ZL_CuOnly (=3) => as 2 but only applied to points diagnosed
! with Cumulus and strictly for sea points (fland<0.01, cf 0.5)
INTEGER, PARAMETER :: DynDiag_ZL_CuOnly = 3
! DynDiag_Ribased (=4) => as 3 but also overrides Cumulus diagnosis if
!     ZH(Ri) > ZLCL+zhloc_depth_fac*(ZHPAR-ZLCL)
! Note that here Ri accounts for gradient adjustment by the non-local scheme.
INTEGER, PARAMETER :: DynDiag_Ribased = 4

! 03 For idyndiag options DynDiag_Ribased: the fractional height into
! the cloud layer reached by the local Ri-based BL depth calculation
REAL :: zhloc_depth_fac = rmdi ! suggested 0.5

! 04
LOGICAL :: l_use_surf_in_ri = .FALSE.
                       ! dbdz in Ri calculation on theta level 1 (ri(k=2))
                       ! false => use dbdz between levels 1 and 2
                       ! true  => use Tstar to give centred difference

! 05 Switch for revised flux-gradient relationships (9C Only)
INTEGER :: flux_grad= imdi
! Flux gradients as in Lock et al. (2000)
INTEGER, PARAMETER :: Locketal2000   = 0
! Flux gradients as in Lock et al (2000) but using coefficients from
! Holtslag and Boville (1993)
INTEGER, PARAMETER :: HoltBov1993 = 1
! Flux gradients as in Lock and Whelan (2006)
INTEGER, PARAMETER :: LockWhelan2006 = 2

! 06 Method of doing entrainment at dsc top (9C only)
INTEGER :: entr_smooth_dec = imdi
! 0 - old method, include surface terms only in weakly coupled,
! no cumulus situations
! 1 - taper method - no hard limit but reduce surface terms according
! to svl difference

! 06a Decoupling thresholds for cloudy boundary layers and cumulus sub-cloud
! layers (larger makes decoupling less likely) (9C Only)
REAL :: dec_thres_cu    = rmdi ! suggested 0.05 (GA7)

! 07 Switch to reset the decoupling threshold to dec_thres_cloud for
! the Ktop iteration (in case it has been set to cu or clear values)
LOGICAL :: l_reset_dec_thres = .FALSE.

! 08 Switch to relax the requirements for diagnosis of Sc over Cu to be
! for all cumulus, not just shallow (9C only)
INTEGER :: relax_sc_over_cu = imdi

! 09 Options for non-local mixing across the LCL in Cu (9C only)
INTEGER :: kprof_cu  = imdi
! klcl_entr (=1) => set K-profile at the LCL using the standard CBL
! entrainment parametrization
INTEGER, PARAMETER :: klcl_entr = 1
! buoy_integ (=2) => use buoyancy flux integration, as used in decoupling
INTEGER, PARAMETER :: buoy_integ = 2

! 10 Switch for parametrizing fluxes across a resolved inversion (9C only)
INTEGER :: bl_res_inv= imdi
! OFF (=0) => not used
! ON  (=1) => used

! 11 Options for blending with 3D Smagorinsky
INTEGER :: blending_option = imdi
! grey_zone_cu (=1) => blending everywhere with grey zone Cu param
INTEGER, PARAMETER :: grey_zone_cu = 1
! blend_except_cu (=2) => blending everywhere except in cumulus layers
INTEGER, PARAMETER :: blend_except_cu = 2

! 12 empirical parameter scaling ustar entrainment term (9C only)
REAL :: a_ent_shr_nml = rmdi ! suggested 1.6

! 13 Improved method to calculate cloud-top radiative flux jump
LOGICAL :: l_new_kcloudtop = .FALSE.

! 14 Include convective effects in TKE diagnostic
LOGICAL :: l_conv_tke = .FALSE.
!=======================================================================
! Implicit solver options
!=======================================================================

! 01 Time weights for boundary layer levels
REAL :: alpha_cd_in (2) = rmdi
REAL, ALLOCATABLE :: alpha_cd(:)

! Parameters for uncond stable numerical solver
! 02 Puns : used in an unstable BL column
REAL :: puns = rmdi ! suggested 0.5
! 03 Pstb : used in an stable BL column
REAL :: pstb = rmdi ! suggested 2.0

!=======================================================================
! a few idealised switches available with the 0A scheme
!=======================================================================
! Choice of Rayleigh friction
INTEGER :: fric_number = imdi
! Options for Rayleigh friction
INTEGER, PARAMETER :: fric_none = 0
INTEGER, PARAMETER :: fric_HeldSuarez_1 = 1
INTEGER, PARAMETER :: fric_HeldSuarez_2 = 2

!=======================================================================
! Options not in namelist, Integers
!=======================================================================

! Generic switch options
INTEGER, PARAMETER :: off = 0  ! Switch disabled
INTEGER, PARAMETER :: on  = 1  ! Switch enabled

! Number of levels for non_local sheme
INTEGER :: nl_bl_levels= imdi

! Options for non-gradient stress following
INTEGER, PARAMETER :: BrownGrant97 = 1
! Brown and Grant (1997), version 2 including a limit on its size
INTEGER, PARAMETER :: BrownGrant97_limited = 2
! Switch for non-gradient stress
INTEGER, PARAMETER :: ng_stress = BrownGrant97_limited

! Options for Prandtl number (in local Ri scheme)
INTEGER, PARAMETER ::  LockMailhot2004 = 1
! Switch for Prandtl number options
INTEGER, PARAMETER :: Prandtl= LockMailhot2004

! Switch for shear-dominated b.l.
! Not a parameter as needs switching off with 1A scheme
INTEGER :: ishear_bl = on

! Switch on the non-local scheme
! Not a parameter as needs switching off for 3D Smag
INTEGER :: non_local_bl = on

! Switch to use implicit weights of 1 for tracers (if on)
INTEGER, PARAMETER :: trweights1 = on

! Options for blending height level
! Implicitly couple at model level one
INTEGER, PARAMETER :: blend_level1  = 0
! Implicitly couple to a specified (fixed) blending height
INTEGER, PARAMETER :: blend_levelk  = 1
! Switch to determine blending height for implicit coupling
INTEGER, PARAMETER :: blend_height_opt = blend_level1

! Switch to enhance entrainment in decoupled stratocu over cu (9C Only)
! OFF (=0) => not used
! Buoyrev_feedback (=1) => buoyancy reversal feedback enhanced by Cu
! for 0<D<0.1
INTEGER, PARAMETER :: Buoyrev_feedback = 1
INTEGER, PARAMETER :: entr_enhance_by_cu = Buoyrev_feedback

! 03 Switch to allow different critical Richardson numbers in the diagnosis
! of BL depth
! 0 => RiC=1 everywhere
! 1 => RiC=0.25 for SHARPEST, 1 otherwise
INTEGER, PARAMETER :: Variable_RiC = on

!=======================================================================
! Options not in namelist, Reals
!=======================================================================
! Weighting of the Louis tail towards long tails:
REAL, PARAMETER :: WeightLouisToLong = 0.0
! 0.0 = Louis tails
! 1.0 = Long tails

! Maximum value of the mixing length (m)
REAL, PARAMETER :: lambda_max_nml = 500.0

! fraction of BL depth for mixing length
REAL, PARAMETER :: lambda_fac = 0.15

! Maximum implied stress gradient across the boundary layer, used to limit
! the explicit stress applied in non-local scalings (m/s2)
REAL, PARAMETER :: max_stress_grad = 0.05

! Max height above LCL for K profile in cumulus
REAL, PARAMETER :: max_cu_depth = 250.0

! timescale for drainage flow development (s)
REAL, PARAMETER :: t_drain = 1800.0

! Horizontal scale for drainage flows (m)
! Currently hard-wired for UKV run over UK
REAL, PARAMETER :: h_scale = 1500.0

! Maximum allowed value of the Prandtl number with the LockMailhot2004 option
REAL, PARAMETER :: pr_max  = 5.0

! Parameters governing the speed of transition from 1D BL to Smagorinsky.
! beta = 0.15 matches Honnert et al well in the BL but a faster transition
! seems appropriate in the free atmosphere above.
REAL, PARAMETER :: beta_bl = 0.15
REAL, PARAMETER :: beta_fa = 1.0

! Constants in linear part of the blending weight function
! linear0 and linear1 are the values of delta/z_scale where linear factor
! equals 0 and 1 respectively
REAL, PARAMETER :: linear1=0.25
REAL, PARAMETER :: linear0=4.0
REAL, PARAMETER :: rlinfac=1.0/(linear0-linear1)

! Critical Ri when SHARPEST stability functions used
REAL, PARAMETER :: RiCrit_sharp = 0.25

! Parameter used in gradient adjustment
REAL, PARAMETER :: a_grad_adj = 3.26

! Parameter used in gradient adjustment
REAL, PARAMETER :: max_t_grad = 1.0e-3

! Powers to use in prescription of equilibrium profiles of stress and
! buoyancy flux in Equilib. SBL model
REAL, PARAMETER :: Muw_SBL= 1.5
REAL, PARAMETER :: Mwt_SBL= 1.0

! Pre-calculated powers for speeding up calculations
REAL, PARAMETER :: one_third = 1.0/3.0
REAL, PARAMETER :: two_thirds = 2.0 * one_third

! Maximum value of TKE to me used in TKE diagnostic
REAL, PARAMETER :: max_tke = 5.0

! Blending height (m) used to calculate the level for surface coupling
REAL, PARAMETER :: h_blend_fix = 0

! Decoupling thresholds for cloudy boundary layers and cumulus sub-cloud
! layers (larger makes decoupling less likely) (9C Only)
REAL, PARAMETER :: dec_thres_cloud = 0.1

!=======================================================================
! Options not in namelist, LOGICAL
!=======================================================================
! Switch for coupled gradient method in Equilibrium SBL model
LOGICAL, PARAMETER :: L_SBLco = .TRUE.

! LambdaM=2*LambdaH (operational setting)
LOGICAL, PARAMETER :: l_lambdam2     = .FALSE.

! Lambdas NOT reduced above NTML_LOCAL+1
LOGICAL, PARAMETER :: l_full_lambdas = .FALSE.

! logical for whether to skip calculations based on start of timestep
! quantities when using semi-lagrangian cycling with Endgame
! This is hard-wired to true in readlsta_4a so Endgame always uses it,
! otherwise it is false
! N.B. results should bit-compare whether this logical is true or false
! and so changing it will be a good test of whether new code has been
! added correctly
LOGICAL :: l_quick_ap2 = .FALSE.

! Logical for whether to compute the explicit momentum fluxes on the
! p-grid on all model-levels.  These might be wanted for a new convection
! scheme, or just for diagnostic purposes.
LOGICAL :: l_calc_tau_at_p = .FALSE.

!=======================================================================
!run_bl namelist
!=======================================================================

NAMELIST/run_bl/ i_bl_vn, sbl_op, cbl_op,                               &
    l_use_surf_in_ri, lambda_min_nml, ritrans,                          &
    local_fa, Keep_Ri_FA, l_bl_mix_qcf, l_conv_tke,                     &
    sg_orog_mixing, fric_heating, calc_prob_of_vis, z_nl_bl_levels,     &
    idyndiag, zhloc_depth_fac, flux_grad, entr_smooth_dec, dec_thres_cu,&
    relax_sc_over_cu, kprof_cu, bl_res_inv, blending_option,            &
    a_ent_shr_nml, alpha_cd_in, puns, pstb, l_reset_dec_thres,          &
    bdy_tke, tke_levels, l_local_above_tkelvs, l_my_initialize,         &
    l_my_ini_zero, my_ini_dbdz_min, l_adv_turb_field, l_my_condense,    &
    l_shcu_buoy, shcu_levels, wb_ng_max, my_lowest_pd_surf,             &
    l_my_prod_adj, my_z_limit_elb, l_print_max_tke, tke_cm_mx,          &
    tke_cm_fa, tke_dlen, l_new_kcloudtop, fric_number, flux_bc_opt

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='BL_OPTION_MOD'

CONTAINS

SUBROUTINE check_run_bl()

! Description:
!   Subroutine to apply logic checks and set control variables based on the
!   options selected in the run_bl namelist.

USE cv_run_mod,   ONLY: i_convection_vn, i_convection_vn_5a,            &
                        i_convection_vn_6a, i_cv_llcs
USE ereport_mod,  ONLY: ereport

USE mym_option_mod, ONLY: deardorff, mymodel25, mymodel3, no_pd_surf,   &
    businger, bh1991, my_length, ddf_length, non_local_like_length

IMPLICIT NONE

CHARACTER (LEN=errormessagelength) :: cmessage      ! used for ereport
INTEGER                            :: icode         ! used for ereport
CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'CHECK_RUN_BL'
REAL(KIND=jprb)                    :: zhook_handle



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

! Check sensible BL version is chosen
CALL chk_var( i_bl_vn, 'i_bl_vn'                                        &
            , [i_bl_vn_0, i_bl_vn_1a, i_bl_vn_9b, i_bl_vn_9c] )

IF (ANY(alpha_cd_in == rmdi)) THEN
  icode = 50
  WRITE(cmessage,'(A)')                                                 &
       'Not enough values of alpha_cd_in have been specified. ' //      &
       'Please specify 2 values'
  CALL ereport(RoutineName, icode, cmessage )
END IF

! Force any globally required defaults
IF (i_bl_vn /= i_bl_vn_9c) THEN
  kprof_cu = off
  cmessage = 'kprof_cu=off set since boundary layer version 9C not chosen'
  icode=-100
  CALL ereport(RoutineName, icode, cmessage)
END IF

! If BL scheme is on...
IF (i_bl_vn /= i_bl_vn_0) THEN

  ! Check options valid for all versions
  CALL chk_var( sbl_op, 'sbl_op',                                       &
              [ long_tails,  sharpest,    sharp_sea_long_land, mes_tails&
              , louis_tails, depth_based, sharp_sea_mes_land,  lem_stability&
              , sharp_sea_louis_land, equilibrium_sbl ] )

  CALL chk_var( ritrans,          'ritrans',          '[0.0:1.0]' )
  CALL chk_var( cbl_op,           'cbl_op',                             &
       [ um_std,neut_cbl,lem_conven,lem_std] )
  CALL chk_var( lambda_min_nml,   'lambda_min_nml',   '[0.0:9999.0]' )
  CALL chk_var( local_fa,         'local_fa',                           &
       [ off, to_sharp_across_1km, ntml_level_corrn, free_trop_layers ] )
  CALL chk_var( Keep_Ri_FA,       'Keep_Ri_FA',   [off,on,except_disc_inv] )
  CALL chk_var( sg_orog_mixing,   'sg_orog_mixing',                     &
       [ off,extended_tail,sg_shear,sg_shear_enh_lambda ] )
  CALL chk_var( fric_heating,     'fric_heating', [off,on] )
  CALL chk_var( calc_prob_of_vis, 'calc_prob_of_vis', '[0.0:1.0]' )
  CALL chk_var(flux_bc_opt,'flux_bc_opt',                                 &
       [interactive_fluxes,specified_fluxes_only,specified_fluxes_tstar,  &
        specified_fluxes_cd])

  ! Check options only valid for 9B/9C schemes
  IF (i_bl_vn == i_bl_vn_9b .OR. i_bl_vn == i_bl_vn_9c) THEN

    IF (i_convection_vn == i_cv_llcs) THEN
      non_local_bl = off
      cmessage = 'non_local_bl=off since LLCS chosen'
      icode=-100
      CALL ereport(RoutineName, icode, cmessage)
    END IF

    CALL chk_var( idyndiag, 'idyndiag',                                 &
         [ off, DynDiag_ZL, DynDiag_ZL_corrn, DynDiag_ZL_CuOnly         &
         , DynDiag_Ribased ] )

    IF ( idyndiag == DynDiag_Ribased )                                  &
         CALL chk_var( zhloc_depth_fac, 'zhloc_depth_fac', '[0.0:1.0]' )

    CALL chk_var( blending_option, 'blending_option',                   &
         [ off, grey_zone_cu, blend_except_cu] )

    ! Check options only valid of 9C scheme
    IF (i_bl_vn == i_bl_vn_9c) THEN

      CALL chk_var( flux_grad, 'flux_grad',                             &
           [ Locketal2000, HoltBov1993, LockWhelan2006 ] )
      CALL chk_var( entr_smooth_dec,  'entr_smooth_dec',  [off,on] )
      CALL chk_var( relax_sc_over_cu, 'relax_sc_over_cu', [off,on] )
      CALL chk_var( kprof_cu,         'kprof_cu', [off,klcl_entr,buoy_integ] )

      IF ( kprof_cu == buoy_integ )                                     &
           CALL chk_var( dec_thres_cu,     'dec_thres_cu',     '[0.0:1.0]' )

      IF ( i_convection_vn == i_convection_vn_5a .OR.                   &
           i_convection_vn == i_convection_vn_6a )                      &
           CALL chk_var( bl_res_inv, 'bl_res_inv', [off,on] )

      CALL chk_var( a_ent_shr_nml,    'a_ent_shr_nml',    '[0.0:10.0]' )


    END IF ! 9C only

  ELSE ! We're running the 1A scheme

    ! Force any 1A defaults
    ishear_bl = off
    cmessage  = 'ishear_bl=off set since boundary layer version 1A chosen'
    icode=-100
    CALL ereport(RoutineName, icode, cmessage)

    blending_option = off
    cmessage = 'blending_option=off set since boundary layer version 1A chosen'
    icode=-100
    CALL ereport(RoutineName, icode, cmessage)

    ! Check options only valid for 1A scheme
    CALL chk_var( bdy_tke, 'bdy_kde', [deardorff,mymodel25,mymodel3] )
    CALL chk_var( my_ini_dbdz_min, 'my_ini_dbdz_min', '[1.0e-10:1.0e10]' )

    IF (l_shcu_buoy) CALL chk_var( wb_ng_max, 'wb_ng_max', '[0.0:9999.0]' )

    CALL chk_var( my_lowest_pd_surf, 'my_lowest_pd_surf',               &
         [ no_pd_surf, businger, bh1991 ] )
    CALL chk_var( my_z_limit_elb, 'my_z_limit_elb', '[0.0:1.0e10]' )

    IF (bdy_tke == deardorff) THEN
      CALL chk_var( tke_cm_mx, 'tke_cm_mx', '[0.0:9999.0]' )
      CALL chk_var( tke_cm_fa, 'tke_cm_fa', '[0.0:9999.0]' )
      CALL chk_var( tke_dlen,  'tke_dlen',                              &
           [ my_length, ddf_length, non_local_like_length ] )
    END IF

  END IF ! version specific

  ! Check implicit solver options
  CALL chk_var( puns, 'puns', '[0.0:100.0]' )
  CALL chk_var( pstb, 'pstb', '[0.0:100.0]' )


ELSE ! i_bl_vn == 0

  CALL chk_var(fric_number,'fric_number',                               &
               [fric_none,fric_HeldSuarez_1,fric_HeldSuarez_2])

END IF ! i_bl_vn /= 0

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_bl

SUBROUTINE print_nlist_run_bl()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_BL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_bl',                             &
    src='bl_option_mod')

WRITE(lineBuffer,*) 'i_bl_vn = ',i_bl_vn
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'sbl_op = ',sbl_op
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'cbl_op = ',cbl_op
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_use_surf_in_ri = ',l_use_surf_in_ri
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'ritrans = ',ritrans
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'lambda_min_nml = ',lambda_min_nml
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_bl_mix_qcf = ',l_bl_mix_qcf
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'local_fa = ',local_fa
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'Keep_Ri_FA = ',Keep_Ri_FA
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'sg_orog_mixing = ',sg_orog_mixing
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'fric_heating = ',fric_heating
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'calc_prob_of_vis = ',calc_prob_of_vis
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,'(A,ES12.4)') 'z_nl_bl_levels = ',z_nl_bl_levels
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'idyndiag = ',idyndiag
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'zhloc_depth_fac = ',zhloc_depth_fac
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'flux_grad = ',flux_grad
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'entr_smooth_dec = ',entr_smooth_dec
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'dec_thres_cu = ',dec_thres_cu
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_reset_dec_thres = ',l_reset_dec_thres
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_new_kcloudtop = ',l_new_kcloudtop
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'relax_sc_over_cu = ',relax_sc_over_cu
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'kprof_cu = ',kprof_cu
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'bl_res_inv = ',bl_res_inv
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'a_ent_shr_nml = ',a_ent_shr_nml
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(linebuffer,*) 'blending_option = ',blending_option
CALL umprint(linebuffer,src='bl_option_mod')
WRITE(linebuffer,*) 'l_conv_tke = ',l_conv_tke
CALL umprint(linebuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'alpha_cd_in = ',alpha_cd_in
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'puns = ',puns
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'pstb = ',pstb
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,'(A,I4)') 'flux_bc_opt = ',flux_bc_opt
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'bdy_tke = ',bdy_tke
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'tke_levels = ',tke_levels
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_local_above_tkelvs = ',l_local_above_tkelvs
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_my_initialize = ',l_my_initialize
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_my_ini_zero = ',l_my_ini_zero
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'my_ini_dbdz_min = ',my_ini_dbdz_min
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_adv_turb_field = ',l_adv_turb_field
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_my_condense = ',l_my_condense
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_shcu_buoy = ',l_shcu_buoy
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'shcu_levels = ',shcu_levels
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'wb_ng_max = ',wb_ng_max
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'my_lowest_pd_surf = ',my_lowest_pd_surf
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_my_prod_adj = ',l_my_prod_adj
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'my_z_limit_elb = ',my_z_limit_elb
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'l_print_max_tke = ',l_print_max_tke
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'tke_cm_mx = ',tke_cm_mx
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'tke_cm_fa = ',tke_cm_fa
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*) 'tke_dlen = ',tke_dlen
CALL umPrint(lineBuffer,src='bl_option_mod')
WRITE(lineBuffer,*)'fric_number = ',fric_number
CALL umPrint(lineBuffer,src='bl_option_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -',                 &
    src='bl_option_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_bl

#if !defined(LFRIC)
SUBROUTINE read_nml_run_bl(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: errorstatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_BL'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 21
INTEGER, PARAMETER :: n_real = 14 + 2
INTEGER, PARAMETER :: n_log = 13

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_bl_vn
  INTEGER :: sbl_op
  INTEGER :: cbl_op
  INTEGER :: local_fa
  INTEGER :: Keep_Ri_FA
  INTEGER :: sg_orog_mixing
  INTEGER :: fric_heating
  INTEGER :: idyndiag
  INTEGER :: flux_grad
  INTEGER :: entr_smooth_dec
  INTEGER :: relax_sc_over_cu
  INTEGER :: kprof_cu
  INTEGER :: bl_res_inv
  INTEGER :: blending_option
  INTEGER :: bdy_tke
  INTEGER :: tke_levels
  INTEGER :: shcu_levels
  INTEGER :: my_lowest_pd_surf
  INTEGER :: tke_dlen
  INTEGER :: fric_number
  INTEGER :: flux_bc_opt
  REAL :: calc_prob_of_vis
  REAL :: zhloc_depth_fac
  REAL :: dec_thres_cu
  REAL :: a_ent_shr_nml
  REAL :: alpha_cd_in(2)
  REAL :: puns
  REAL :: pstb
  REAL :: my_ini_dbdz_min
  REAL :: wb_ng_max
  REAL :: my_z_limit_elb
  REAL :: tke_cm_mx
  REAL :: tke_cm_fa
  REAL :: ritrans
  REAL :: lambda_min_nml
  REAL :: z_nl_bl_levels
  LOGICAL :: l_bl_mix_qcf
  LOGICAL :: l_local_above_tkelvs
  LOGICAL :: l_my_initialize
  LOGICAL :: l_my_ini_zero
  LOGICAL :: l_adv_turb_field
  LOGICAL :: l_my_condense
  LOGICAL :: l_shcu_buoy
  LOGICAL :: l_my_prod_adj
  LOGICAL :: l_print_max_tke
  LOGICAL :: l_reset_dec_thres
  LOGICAL :: l_use_surf_in_ri
  LOGICAL :: l_new_kcloudtop
  LOGICAL :: l_conv_tke
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,          &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=run_bl, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_BL", iomessage)

  my_nml % i_bl_vn            = i_bl_vn
  my_nml % sbl_op             = sbl_op
  my_nml % cbl_op             = cbl_op
  my_nml % local_fa           = local_fa
  my_nml % Keep_Ri_FA         = Keep_Ri_FA
  my_nml % sg_orog_mixing     = sg_orog_mixing
  my_nml % fric_heating       = fric_heating
  my_nml % idyndiag           = idyndiag
  my_nml % flux_grad          = flux_grad
  my_nml % entr_smooth_dec    = entr_smooth_dec
  my_nml % relax_sc_over_cu   = relax_sc_over_cu
  my_nml % kprof_cu           = kprof_cu
  my_nml % bl_res_inv         = bl_res_inv
  my_nml % blending_option    = blending_option
  my_nml % bdy_tke            = bdy_tke
  my_nml % tke_levels         = tke_levels
  my_nml % shcu_levels        = shcu_levels
  my_nml % my_lowest_pd_surf  = my_lowest_pd_surf
  my_nml % tke_dlen           = tke_dlen
  my_nml % fric_number        = fric_number
  my_nml % flux_bc_opt        = flux_bc_opt
  ! end of integers
  my_nml % calc_prob_of_vis  = calc_prob_of_vis
  my_nml % zhloc_depth_fac   = zhloc_depth_fac
  my_nml % dec_thres_cu      = dec_thres_cu
  my_nml % a_ent_shr_nml     = a_ent_shr_nml
  my_nml % alpha_cd_in       = alpha_cd_in
  my_nml % puns              = puns
  my_nml % pstb              = pstb
  my_nml % my_ini_dbdz_min   = my_ini_dbdz_min
  my_nml % wb_ng_max         = wb_ng_max
  my_nml % my_z_limit_elb    = my_z_limit_elb
  my_nml % tke_cm_mx         = tke_cm_mx
  my_nml % tke_cm_fa         = tke_cm_fa
  my_nml % ritrans           = ritrans
  my_nml % lambda_min_nml    = lambda_min_nml
  my_nml % z_nl_bl_levels    = z_nl_bl_levels
  ! end of reals
  my_nml % l_bl_mix_qcf         = l_bl_mix_qcf
  my_nml % l_local_above_tkelvs = l_local_above_tkelvs
  my_nml % l_my_initialize      = l_my_initialize
  my_nml % l_my_ini_zero        = l_my_ini_zero
  my_nml % l_adv_turb_field     = l_adv_turb_field
  my_nml % l_my_condense        = l_my_condense
  my_nml % l_shcu_buoy          = l_shcu_buoy
  my_nml % l_my_prod_adj        = l_my_prod_adj
  my_nml % l_print_max_tke      = l_print_max_tke
  my_nml % l_reset_dec_thres    = l_reset_dec_thres
  my_nml % l_use_surf_in_ri     = l_use_surf_in_ri
  my_nml % l_new_kcloudtop      = l_new_kcloudtop
  my_nml % l_conv_tke           = l_conv_tke
  ! end of logicals

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  i_bl_vn            = my_nml % i_bl_vn
  sbl_op             = my_nml % sbl_op
  cbl_op             = my_nml % cbl_op
  local_fa           = my_nml % local_fa
  Keep_Ri_FA         = my_nml % Keep_Ri_FA
  sg_orog_mixing     = my_nml % sg_orog_mixing
  fric_heating       = my_nml % fric_heating
  idyndiag           = my_nml % idyndiag
  flux_grad          = my_nml % flux_grad
  entr_smooth_dec    = my_nml % entr_smooth_dec
  relax_sc_over_cu   = my_nml % relax_sc_over_cu
  kprof_cu           = my_nml % kprof_cu
  bl_res_inv         = my_nml % bl_res_inv
  blending_option    = my_nml % blending_option
  bdy_tke            = my_nml % bdy_tke
  tke_levels         = my_nml % tke_levels
  shcu_levels        = my_nml % shcu_levels
  my_lowest_pd_surf  = my_nml % my_lowest_pd_surf
  tke_dlen           = my_nml % tke_dlen
  fric_number        = my_nml % fric_number
  flux_bc_opt        = my_nml % flux_bc_opt
  ! end of integers
  calc_prob_of_vis  = my_nml % calc_prob_of_vis
  zhloc_depth_fac   = my_nml % zhloc_depth_fac
  dec_thres_cu      = my_nml % dec_thres_cu
  a_ent_shr_nml     = my_nml % a_ent_shr_nml
  alpha_cd_in       = my_nml % alpha_cd_in
  puns              = my_nml % puns
  pstb              = my_nml % pstb
  my_ini_dbdz_min   = my_nml % my_ini_dbdz_min
  wb_ng_max         = my_nml % wb_ng_max
  my_z_limit_elb    = my_nml % my_z_limit_elb
  tke_cm_mx         = my_nml % tke_cm_mx
  tke_cm_fa         = my_nml % tke_cm_fa
  ritrans           = my_nml % ritrans
  lambda_min_nml    = my_nml % lambda_min_nml
  z_nl_bl_levels    = my_nml % z_nl_bl_levels
  ! end of reals
  l_bl_mix_qcf         = my_nml % l_bl_mix_qcf
  l_local_above_tkelvs = my_nml % l_local_above_tkelvs
  l_my_initialize      = my_nml % l_my_initialize
  l_my_ini_zero        = my_nml % l_my_ini_zero
  l_adv_turb_field     = my_nml % l_adv_turb_field
  l_my_condense        = my_nml % l_my_condense
  l_shcu_buoy          = my_nml % l_shcu_buoy
  l_my_prod_adj        = my_nml % l_my_prod_adj
  l_print_max_tke      = my_nml % l_print_max_tke
  l_reset_dec_thres    = my_nml % l_reset_dec_thres
  l_use_surf_in_ri     = my_nml % l_use_surf_in_ri
  l_new_kcloudtop      = my_nml % l_new_kcloudtop
  l_conv_tke           = my_nml % l_conv_tke
  ! end of logicals

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_bl
#endif

END MODULE bl_option_mod

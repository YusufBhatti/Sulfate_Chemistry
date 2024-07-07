! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Main routine for IAU scheme

SUBROUTINE iau (                                  &
                 l_mr_physics,                    & ! in
                 frac_typ,                        & ! in
                 snowdepth_z, nsnow,              & ! in
                 vol_smc_wilt, vol_smc_crit,      & ! in
                 vol_smc_sat,                     & ! in
                 u,      v,          w,           & ! inout
                 theta,  rho,        murk,        & ! inout
                 q,      qCL,        qCF,         & ! inout
                 TStar,  TStar_tile, TSoil, smcl, & ! inout
                 tsnow, tstar_anom,               & ! inout
                 p_star, p,                       & ! inout
                 p_theta_levels,                  & ! inout
                 exner,                           & ! inout
                 exner_theta_levels,              & ! inout
                 snow_depth,                      & ! inout
                 area_cloud_fraction,             & ! inout
                 bulk_cloud_fraction,             & ! inout
                 cloud_fraction_liquid,           & ! inout
                 cloud_fraction_frozen,           & ! inout
                 dust_div1, dust_div2, dust_div3, & ! inout
                 dust_div4, dust_div5, dust_div6, & ! inout
                 ozone_tracer, SO2, tracer_ukca )   ! inout

! Description:
!
!   Main routine for the Incremental Analysis Update (IAU) scheme, as described
!   in UMDP 31.
!
!   Model fields can be updated in two ways:
!
!     1. Directly, via corresponding increment fields in a series of IAU
!        increment files.
!     2. Indirectly, via the increments to other variables.
!
!   Currently, direct incrementing is supported for the following fields:
!
!     u, v, w, theta, rho, exner, p, q, qCL, qCF, murk (aerosol), ozone,
!     sfctemp, TSoil, soilm, dust_div1-6
!
!   theta, rho and exner may also be updated indirectly via the following
!   relationships:
!
!     theta - hydrostatic equation
!     rho   - equation of state
!     exner - expression in terms of p
!
!   These three options are activated via namelist switches. When one of these
!   is activated, direct updates to the corresponding field are switched off.
!
!   As well as the above set of fields, the IAU increment files may also
!   contain increments to the total specific humidity qT. When qT increments
!   are available, there are two options: (1) use the routine Var_DiagCloud to
!   diagnose increments to q, qCL and (optionally) qCF, or (2) add qT' to q.
!   Direct and indirect increments to q, qCL or qCF are not allowed on
!   the same timestep.
!
!   When the PC2 cloud scheme is in use, direct updates to q, qCL or qCF
!   are followed up by adjustments to q, qCL and theta to ensure compatibilty
!   with the scheme.
!
!   There is also an option to pass the level-one temperature increments on to
!   the surface temperature fields (TStar and TStar_tile) at land points, and
!   the top-level soil temperature (TSoil(1)) field.
!
!   Sfctemp, TSoil and soilm are used as initial perturbations in ensemble
!   forecasts. Sfctemp is passed to TStar (at all points) and to TStar_tile
!   and TSoil(1) over land. Soilm and TSoil are passed to all soil levels.
!   These fields are automatically processed if present in the increment file.
!
!   At the end of the routine two options are available to limit humidity and
!   cloud. The first option is to remove supersaturation wrt water. The second
!   is to apply limits to tropospheric and non-troposphere humidities via the
!   routine QLimits, which also removes cloud outside the troposphere. The
!   frequency with which QLimits is called is controlled via a namelist
!   variable. For example, it can be called on every IAU call for which
!   increments are to be added, or only at the end of the IAU insertion period.
!
!   Finally, in order to minimise problems with the spin-down in precipitation
!   following the addition of an analysis increment, options are available to
!   modify the output of the cloud incrementing operator (Var_DiagCloud).
!   Also, the removal of supersaturation noted above can be done with respect
!   to a scaled value of qsat.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE swap_bounds_mv_mod, ONLY: swap_bounds_mv
USE submodel_mod, ONLY: atmos_im

USE atm_fields_bounds_mod, ONLY:   &
    pdims, pdims_s,                 &
    tdims, tdims_s, tdims_l,        &
    udims, udims_s, udims_l,        &
    vdims, vdims_s, vdims_l,        &
    wdims, wdims_s, wdims_l

USE planet_constants_mod, ONLY:    &
    kappa,                          &
    cp,                             &
    c_virtual,                      &
    p_zero,                         &
    r,                              &
    recip_kappa,                    &
    recip_epsilon,                  &
    g

USE Control_Max_Sizes, ONLY:       &
    max_121_rows,                   &
    max_look,                       &
    max_n_intf_a,                   &
    max_req_thpv_levs,              &
    max_sponge_width,               &
    max_updiff_levels

USE dyn_coriolis_mod, ONLY:        &
    f3_at_v

USE ereport_mod, ONLY:             &
    ereport

USE filenamelength_mod, ONLY :     &
    filenamelength

USE IAU_mod, ONLY:                 &
    IAU_NumFldCodes,                &
    IAU_FldCodes,                   &
    CallFreq_Never,                 &
    CallFreq_EveryCallIncsToAdd,    &
    CallFreq_EndInc1Window,         &
    CallFreq_LastIAUCall,           &
    IAU_incs,                       &
    IAU_LocFldLens,                 &
    q_min,                          &
    IAU_LastCallTS,                 &
    ukca_O3_index,                  &
    ukca_NO2_index,                 &
    Num_IAU_incs,                   &
    L_IAU_CalcExnerIncs,            &
    L_IAU_CalcThetaIncs,            &
    L_IAU_CalcRhoIncs,              &
    L_IAU_IncTStar,                 &
    L_IAU_IncTStar_tile,            &
    L_IAU_IncTSurf_Snow,            &
    L_IAU_RemoveSS,                 &
    L_IAU_ApplyQTCorrections,       &
    L_IAU_Add_qT_prime_to_q,        &
    L_IAU_IncrementIce,             &
    L_IAU_ScaleCloud,               &
    L_IAU_LimitUpperThetaIncs,      &
    L_IAU_SetOzoneMin,              &
    L_IAU_SMCChecks,                &
    L_IAU_TSoilLandIceMask,         &
    IAU_LimitUpperThetaIncs_pBound, &
    IAU_LimitUpperThetaIncs_maxInc, &
    IAU_QLimitsCallFreq,            &
    L_IAU_QLimitsDiags,             &
    L_IAU_RmNonTropQIncs,           &
    IAU_trop_min_RH,                &
    IAU_nonTrop_max_q,              &
    IAU_nonTrop_min_q,              &
    IAU_nonTrop_max_RH,             &
    IAU_trop_min_p,                 &
    IAU_trop_max_PV,                &
    IAU_nonTrop_max_p,              &
    L_IAU_IgnoreTopLevThetaIncs

USE io_constants, ONLY:             &
    ioOpenReadOnly,                 &
    ioNameProvided,                 &
    ioNoDelete

USE model_file, ONLY:               &
    model_file_open,                &
    model_file_close

USE file_manager, ONLY:             &
    assign_file_unit, release_file_unit

USE level_heights_mod, ONLY:       &
    r_theta_levels,                 &
    r_rho_levels,                   &
    r_at_u,                         &
    r_at_v

USE cloud_inputs_mod, ONLY:        &
    rhcrit,                         &
    i_cld_area,                     &
    i_cld_vn,                       &
    l_pc2_check_init
USE pc2_constants_mod, ONLY: acf_off, acf_cusack, acf_brooks, i_cld_pc2

USE ancilcta_namelist_mod, ONLY:    &
    l_sstanom

USE jules_surface_types_mod, ONLY: &
    ntype, lake, ice

USE jules_surface_mod, ONLY: &
    l_aggregate

USE jules_snow_mod, ONLY:               &
    nsmax, dzsnow
 
USE jules_soil_mod, ONLY:               &
    dzsoil
 
USE parkind1, ONLY:                &
    jprb,                           &
    jpim

USE umPrintMgr, ONLY:              &
    umPrint,                        &
    umMessage,                      &
    PrintStatus,                    &
    PrStatus_Normal
USE science_fixes_mod, ONLY: l_iau_pc2check
USE swapable_field_mod, ONLY:      &
    swapable_field_pointer_type

USE trignometric_mod, ONLY:        &
    sec_v_latitude,                 &
    tan_v_latitude,                 &
    sec_theta_latitude

USE cderived_mod, ONLY: delta_lambda, delta_phi

USE UM_ParVars, ONLY:              &
    halo_i,                         &
    halo_j,                         &
    offx,                           &
    offy

USE visbty_constants_mod, ONLY:    &
    aero0,                          &
    aeromax

USE water_constants_mod, ONLY:     &
    lc,                             &
    tm,                             &
    rho_water

USE yomhook, ONLY:                  &
    lhook,                          &
    dr_hook

USE lookup_addresses

USE qlimits_mod, ONLY: qlimits
USE readiaufield_mod, ONLY: readiaufield
USE var_diagcloud_mod, ONLY: var_diagcloud

USE atm_fields_mod, ONLY: m_v, m_cl, m_cf                          &
                          ,m_cf2, m_r, m_gr, thetav, theta3d => theta, &
                          qgraup,qrain,qcf2

USE mphys_inputs_mod,  ONLY:  l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup

USE update_rho_mod,  ONLY: update_dryrho_mix_tl_n
USE eg_q_to_mix_mod, ONLY: eg_q_to_mix

USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE field_types, ONLY:              &
    fld_type_p,                     &
    fld_type_u,                     &
    fld_type_v


USE nlsizes_namelist_mod, ONLY:                                          &
    bl_levels, global_row_length, global_rows, land_field, model_levels, &
    n_rows, ntiles, row_length, rows, sm_levels, st_levels,              &
    theta_field_size, v_field_size, tr_ukca

USE model_time_mod, ONLY: &
    secs_per_stepim, stepim

USE land_soil_dimensions_mod, ONLY: land_index

USE errormessagelength_mod, ONLY: errormessagelength

USE um_stashcode_mod, ONLY: stashcode_O3, stashcode_NO2, stashcode_SO2

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_wat_new     => qsat_wat,                       &
                    l_new_qsat_acassim !Currently defaults to FALSE

USE initial_pc2_check_mod, ONLY: initial_pc2_check
USE pc2_assim_mod, ONLY: pc2_assim
IMPLICIT NONE

! DEPENDS ON: qsat_wat
! DEPENDS ON: expand21
! DEPENDS ON: pack21
! DEPENDS ON: calc_p_star
! DEPENDS ON: calc_exner_at_theta
! DEPENDS ON: calc_p_from_exner
! DEPENDS ON: calc_pv_at_theta


! Subroutine arguments:

LOGICAL, INTENT(IN) :: l_mr_physics ! Use mixing ratio code?

REAL,    INTENT(IN) :: frac_typ ( land_field, ntype )

REAL,    INTENT(IN) :: snowdepth_z ( land_field, ntype )
!                      ! Depth of snow in m
REAL,    INTENT(IN) :: nsnow ( land_field, ntype )
!                      ! Number of snow layers

REAL,    INTENT(INOUT), TARGET ::           &
  u      ( udims_s%i_start : udims_s%i_end, &
           udims_s%j_start : udims_s%j_end, &
           udims_s%k_start : udims_s%k_end )

REAL,    INTENT(INOUT), TARGET ::           &
  v      ( vdims_s%i_start : vdims_s%i_end, &
           vdims_s%j_start : vdims_s%j_end, &
           vdims_s%k_start : vdims_s%k_end )

REAL,    INTENT(INOUT) ::                   &
  w      ( wdims_s%i_start : wdims_s%i_end, &
           wdims_s%j_start : wdims_s%j_end, &
           wdims_s%k_start : wdims_s%k_end )

REAL,    INTENT(INOUT), TARGET ::           &
  theta  ( tdims_s%i_start : tdims_s%i_end, &
           tdims_s%j_start : tdims_s%j_end, &
           tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT), TARGET ::           &
  rho    ( pdims_s%i_start : pdims_s%i_end, &
           pdims_s%j_start : pdims_s%j_end, &
           pdims_s%k_start : pdims_s%k_end )

REAL,    INTENT(INOUT) ::                   &
  murk   ( tdims_s%i_start : tdims_s%i_end, &
           tdims_s%j_start : tdims_s%j_end, &
           tdims_s%k_start : tdims_s%k_end )

! Note that level-zero humidity fields are not currently updated for ENDGame.
! We may add support for this later.
REAL,    INTENT(INOUT) ::                   &
  q      ( tdims_l%i_start : tdims_l%i_end, &
           tdims_l%j_start : tdims_l%j_end, &
           tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                   &
  qCL    ( tdims_l%i_start : tdims_l%i_end, &
           tdims_l%j_start : tdims_l%j_end, &
           tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                   &
  qCF    ( tdims_l%i_start : tdims_l%i_end, &
           tdims_l%j_start : tdims_l%j_end, &
           tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(IN)    :: vol_smc_wilt (land_field), &
                          vol_smc_crit (land_field), &
                          vol_smc_sat  (land_field)

REAL,    INTENT(INOUT) :: TStar      ( theta_field_size )
REAL,    INTENT(INOUT) :: tstar_anom ( theta_field_size )

REAL,    INTENT(INOUT) :: TStar_tile ( land_field, ntiles )

REAL,    INTENT(INOUT) :: TSoil      ( land_field, st_levels ), &
                          smcl       ( land_field, sm_levels )

REAL,    INTENT(INOUT) :: tsnow ( land_field, ntype, nsmax )
!                         ! Temperature of snow layers

REAL,    INTENT(INOUT) :: p_star     ( pdims_s%i_start : pdims_s%i_end,  &
                                       pdims_s%j_start : pdims_s%j_end )

REAL,    INTENT(INOUT) ::                                  &
  p                     ( pdims_s%i_start : pdims_s%i_end, &
                          pdims_s%j_start : pdims_s%j_end, &
                          pdims_s%k_start : pdims_s%k_end + 1 )

REAL,    INTENT(INOUT) ::                                  &
  p_theta_levels        ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT), TARGET ::                          &
  exner                 ( pdims_s%i_start : pdims_s%i_end, &
                          pdims_s%j_start : pdims_s%j_end, &
                          pdims_s%k_start : pdims_s%k_end + 1 )

REAL,    INTENT(INOUT) ::                                  &
  exner_theta_levels    ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT) :: snow_depth ( theta_field_size )

REAL,    INTENT(INOUT) :: area_cloud_fraction          &
                        ( tdims%i_start : tdims%i_end, &
                          tdims%j_start : tdims%j_end, &
                                      1 : tdims%k_end )

REAL,    INTENT(INOUT) ::                                  &
  bulk_cloud_fraction   ( tdims_l%i_start : tdims_l%i_end, &
                          tdims_l%j_start : tdims_l%j_end, &
                          tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                                  &
  cloud_fraction_liquid ( tdims_l%i_start : tdims_l%i_end, &
                          tdims_l%j_start : tdims_l%j_end, &
                          tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                                  &
  cloud_fraction_frozen ( tdims_l%i_start : tdims_l%i_end, &
                          tdims_l%j_start : tdims_l%j_end, &
                          tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                                  &
  dust_div1             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )
REAL,    INTENT(INOUT) ::                                  &
  dust_div2             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )
REAL,    INTENT(INOUT) ::                                  &
  dust_div3             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )
REAL,    INTENT(INOUT) ::                                  &
  dust_div4             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )
REAL,    INTENT(INOUT) ::                                  &
  dust_div5             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )
REAL,    INTENT(INOUT) ::                                  &
  dust_div6             ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT) ::                                  &
  ozone_tracer          ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT) ::                                  &
  SO2                   ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT) ::                                  &
  tracer_ukca           ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end, &
                          tr_ukca )

! Local constants:

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'IAU'

REAL, PARAMETER :: CC_tol = 1.0e-12 ! qCL/qCF tolerance for cloud clearing
REAL, PARAMETER :: oz_min = 1.0e-8  ! Minimum value allowed for ozone
                                    ! after addition of ozone increments

REAL :: p_zero_recip ! = 1.0/p_zero
REAL :: exner_exp    ! = 1.0/kappa - 1.0

! PC2 options.
! Seek advice from the PC2 team before altering these parameters from .TRUE.;
! you may need to have put in place large amounts of extra code first.
LOGICAL, PARAMETER :: L_pc2_cond = .TRUE. ! Do condensation and liquid cloud
                                          ! fraction changes?
LOGICAL, PARAMETER :: L_pc2_cfl  = .TRUE. ! Do liquid cloud fraction changes?
                                          ! This requires that the
                                          ! condensation as a result of
                                          ! assimilation is calculated either
                                          ! directly or via the estimate
                                          ! selected by setting l_pc2_cond to
                                          ! .TRUE..
                                          ! (Note that one must not run with
                                          ! condensation changes on but the
                                          ! liquid cloud fraction changes off.)
LOGICAL, PARAMETER :: L_pc2_cff  = .TRUE. ! Do ice cloud fraction changes?
                                          ! This requires that the ice
                                          ! increment from assimilation is
                                          ! calculated directly.

! Local variables:

INTEGER :: i, j, k, js
INTEGER :: TS_len_secs
INTEGER :: Sec
INTEGER :: IncNum
INTEGER :: FieldNum
INTEGER :: WeightNum
INTEGER :: tile_num
INTEGER :: field_size
INTEGER :: flds_StartAddr
INTEGER :: s_addr
INTEGER :: e_addr
INTEGER :: addr
INTEGER :: Code
INTEGER :: LocFldLen
INTEGER :: buflen      ! global field size read in
INTEGER :: ICode
INTEGER :: i_field     ! Counter for swap_bounds
INTEGER :: IAU_unit    ! File unit number

REAL :: pBound
REAL :: maxInc
REAL :: delta_r
REAL :: thetaV_in_int
REAL :: thetaV_new_int
REAL :: Weight
REAL :: Weight_1
REAL :: Weight_2
REAL :: Term_1
REAL :: Term_2
REAL :: vol_smc

REAL :: tile_frac      ! Fraction of current surface tile
REAL :: wgt_snow(nsmax)! Weighting factor for each snow level
REAL :: wgt_soil_snow  ! Weighting soil below snow
REAL :: z_snow         ! Depth of middle of snow layer

REAL :: t1_in     ( tdims%i_start : tdims%i_end, &
                    tdims%j_start : tdims%j_end )
REAL :: t1_inc    ( theta_field_size )   ! work array for theta(1) perts
REAL :: p_tmp     ( pdims%i_start : pdims%i_end, &
                    pdims%j_start : pdims%j_end )
REAL :: q_sat_wat ( tdims%i_start : tdims%i_end, &
                    tdims%j_start : tdims%j_end )

! Array for holding a single IAU increment field:
REAL, ALLOCATABLE :: IAUIncFld ( : )

! Work arrays required for diagnosing humidity/cloud incs from qT incs:
REAL, ALLOCATABLE :: p_theta_levels_old        (:,:,:)
REAL, ALLOCATABLE :: qT                        (:,:,:)
REAL, ALLOCATABLE :: qT_plus                   (:,:,:)
REAL, ALLOCATABLE :: VarTemp                   (:)
REAL, ALLOCATABLE :: AreaDivBulk               (:,:,:)
REAL, ALLOCATABLE :: cloud_fraction_liquid_0   (:)
REAL, ALLOCATABLE :: cloud_fraction_liquid_plus(:)
REAL, ALLOCATABLE :: qCL_0                     (:)
REAL, ALLOCATABLE :: qCL_plus                  (:)
REAL, ALLOCATABLE :: Var_qT                    (:)
REAL, ALLOCATABLE :: initialScaling            (:,:,:)
REAL, ALLOCATABLE :: scalingCoeff              (:,:,:)
REAL, ALLOCATABLE :: qcl_Inc                   (:,:,:)
REAL, ALLOCATABLE :: Cl_Inc                    (:,:,:)
REAL, ALLOCATABLE :: Cf_Inc                    (:,:,:)
REAL, ALLOCATABLE :: Cb_Inc                    (:,:,:)
REAL, ALLOCATABLE :: Var_BGT                   (:)
REAL, ALLOCATABLE :: cloud_fraction_frozen_0   (:)
REAL, ALLOCATABLE :: cloud_fraction_frozen_plus(:)
REAL, ALLOCATABLE :: qCF_0                     (:)
REAL, ALLOCATABLE :: qCF_max                   (:)
REAL, ALLOCATABLE :: qCF_plus                  (:)

! Work arrays required for PC2 calculations:
REAL, ALLOCATABLE :: q_work   (:,:,:)
REAL, ALLOCATABLE :: qcl_work (:,:,:)
REAL, ALLOCATABLE :: qcf_in   (:,:,:)
REAL, ALLOCATABLE :: t_work   (:,:,:)
REAL, ALLOCATABLE :: p_work   (:,:,:)
REAL, ALLOCATABLE :: delta_q  (:,:,:)
REAL, ALLOCATABLE :: delta_qcl(:,:,:)
REAL, ALLOCATABLE :: delta_qcf(:,:,:)
REAL, ALLOCATABLE :: delta_t  (:,:,:)
REAL, ALLOCATABLE :: delta_p  (:,:,:)

! Work arrays required for call to QLimits:
REAL, ALLOCATABLE :: q_in       (:,:,:)
REAL, ALLOCATABLE :: pv_at_theta(:,:,:)

! Work arrays for sfctemp and soilm perturbations
REAL,ALLOCATABLE :: t2_inc    (:)
REAL,ALLOCATABLE :: s1_inc    (:)
REAL,ALLOCATABLE :: ts1_inc   (:)

! Other large work arrays to be allocated on demand:
REAL, ALLOCATABLE :: exner_in  (:,:,:)
REAL, ALLOCATABLE :: thetaV_in (:,:,:)
REAL, ALLOCATABLE :: thetaV_new(:,:,:)

TYPE(swapable_field_pointer_type) :: fields_to_swap(4)

LOGICAL :: IncsThisTS
LOGICAL :: qTIncsThisTS
LOGICAL :: CallQLimitsThisTS
LOGICAL :: FileOpen
LOGICAL :: ReadIAUField_FirstCallThisInc

LOGICAL :: l_include_halos = .FALSE.

CHARACTER(LEN=4)   :: IncNumStr
CHARACTER(LEN=9)   :: TimeStr
CHARACTER(LEN=errormessagelength) :: CMessage

CHARACTER(filenamelength) :: FilePath

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header ---------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

p_zero_recip = 1.0/p_zero
exner_exp    = 1.0/kappa - 1.0

! Array for holding a single IAU increment field:
ALLOCATE ( IAUIncFld ( v_field_size ) )

! Timestep length in seconds:
TS_len_secs = secs_per_stepim(atmos_im)

IF (PrintStatus >= PrStatus_Normal) THEN
  Sec = stepim(atmos_im) * TS_len_secs
  WRITE (TimeStr,'(I3,":",I2.2,":",I2.2)')  &
    Sec/3600, MOD(Sec/60, 60), MOD(Sec, 60)
  WRITE(umMessage,'(A)')                                                &
        'IAU: Entering IAU. Time since basis time (hr:mm:ss): '//TimeStr
  CALL umPrint(umMessage,src='iau')
END IF

!-------------------------------------------------------------------------------
! [1]: Return if there are no increments to add on this timestep.
!-------------------------------------------------------------------------------

IncsThisTS = .FALSE.
DO IncNum = 1, Num_IAU_incs
  IF (IAU_incs(IncNum) % FieldsToAdd) THEN
    IF (stepim(atmos_im) >= IAU_incs(IncNum) % InsertionStartTS .AND. &
        stepim(atmos_im) <= IAU_incs(IncNum) % InsertionEndTS) THEN
      WeightNum = stepim(atmos_im) - IAU_incs(IncNum) % InsertionStartTS + 1
      IF (IAU_incs(IncNum) % Weights(WeightNum) /= 0.0) THEN
        IncsThisTS = .TRUE.
        EXIT
      END IF
    END IF
  END IF
END DO

IF (.NOT. IncsThisTS) THEN
  IF (PrintStatus >= PrStatus_Normal) THEN
    WRITE(umMessage,'(A)') 'IAU: No increments to add on this timestep'
    CALL umPrint(umMessage,src='iau')
    WRITE(umMessage,'(A)')                                              &
         'IAU: Exiting IAU. Time since basis time (hr:mm:ss): '//TimeStr
    CALL umPrint(umMessage,src='iau')
  END IF
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

!-------------------------------------------------------------------------------
! [2]: If exner increments are to be calculated from p increments, make sure p
!      is consistent with exner.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcExnerIncs) THEN

  p(1:row_length,1:rows,:) = exner(1:row_length,1:rows,:)**recip_kappa &
                           * p_zero

END IF

!-------------------------------------------------------------------------------
! [3]: Determine whether there are qT increments to be added on this timestep.
!-------------------------------------------------------------------------------

qTIncsThisTS = .FALSE.
DO IncNum = 1, Num_IAU_incs
  IF (IAU_incs(IncNum) % Contains_qT) THEN
    IF (stepim(atmos_im) >= IAU_incs(IncNum) % InsertionStartTS .AND. &
        stepim(atmos_im) <= IAU_incs(IncNum) % InsertionEndTS) THEN
      qTIncsThisTS = .TRUE.
    END IF
  END IF
END DO

!-------------------------------------------------------------------------------
! [4]: Determine whether QLimits will be called on this timestep.
!-------------------------------------------------------------------------------

SELECT CASE (IAU_QLimitsCallFreq)
CASE (CallFreq_Never)
  CallQLimitsThisTS = .FALSE.
CASE (CallFreq_EveryCallIncsToAdd)
  CallQLimitsThisTS = .TRUE.
CASE (CallFreq_EndInc1Window)
  CallQLimitsThisTS = ( IAU_incs(1) % FieldsToAdd .AND. &
                        stepim(atmos_im) == IAU_incs(1) % InsertionEndTS )
CASE (CallFreq_LastIAUCall)
  CallQLimitsThisTS = ( stepim(atmos_im) == IAU_LastCallTS )
END SELECT

!-------------------------------------------------------------------------------
! [5]: Save fields required for calculating indirect field increments.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! [5.1]: Exner and thetaV.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcThetaIncs .OR. L_IAU_CalcRhoIncs) THEN

  ALLOCATE (exner_in (pdims%i_start:pdims%i_end, &
                      pdims%j_start:pdims%j_end, &
                      pdims%k_start:pdims%k_end+1))
  ALLOCATE (thetaV_in(tdims%i_start:tdims%i_end, &
                      tdims%j_start:tdims%j_end, &
                      tdims%k_start:tdims%k_end)  )

  exner_in (:,:,:)             = exner(pdims%i_start:pdims%i_end, &
                                       pdims%j_start:pdims%j_end,: )

  thetaV_in(:,:,tdims%k_start:tdims%k_end)  =                     &
                                 theta(tdims%i_start:tdims%i_end, &
                                       tdims%j_start:tdims%j_end, &
                                       tdims%k_start:tdims%k_end) &
                               * ( 1.0 + C_Virtual                &
                               *   q  (tdims%i_start:tdims%i_end, &
                                       tdims%j_start:tdims%j_end, &
                                       tdims%k_start:tdims%k_end) )

END IF

!-------------------------------------------------------------------------------
! [5.2]: Level 1 temperature.
!-------------------------------------------------------------------------------

! Level 1 temperature:
IF (L_IAU_IncTStar) THEN

  t1_in(:,:) = exner_theta_levels(tdims%i_start:tdims%i_end,    &
                                  tdims%j_start:tdims%j_end, 1) &
             * theta             (tdims%i_start:tdims%i_end,    &
                                  tdims%j_start:tdims%j_end, 1)

END IF

!-------------------------------------------------------------------------------
! [5.3]: Fields required for processing qT increments.
!-------------------------------------------------------------------------------

IF (qTIncsThisTS) THEN

  ! NOTE : The zeroth theta level is being excluded in the initial
  !        implementation of ENDGAME IAU (UM 8.2) - to be reconsidered later.
  !        Moisture and temperature arrays therefore start from level 1.
  ALLOCATE (qT     (tdims%i_start:tdims%i_end, &
                    tdims%j_start:tdims%j_end, &
                                1:tdims%k_end))
  ALLOCATE (qT_plus(tdims%i_start:tdims%i_end, &
                    tdims%j_start:tdims%j_end, &
                                1:tdims%k_end))

  qT(:,:,:) = q  (tdims%i_start:tdims%i_end, &
                  tdims%j_start:tdims%j_end, &
                              1:tdims%k_end) &
            + qCL(tdims%i_start:tdims%i_end, &
                  tdims%j_start:tdims%j_end, &
                              1:tdims%k_end)

  IF (L_IAU_IncrementIce .AND. L_IAU_ApplyQTCorrections)   &
    qT(:,:,:) = qT(:,:,:) + qCF(tdims%i_start:tdims%i_end, &
                                tdims%j_start:tdims%j_end, &
                                            1:tdims%k_end)

  ! Initialise qT_plus:
  qT_plus(:,:,:) = qT(:,:,:)

  IF (.NOT. L_IAU_Add_qT_prime_to_q) THEN

    field_size = tdims%i_end * tdims%j_end * tdims%k_end

    ALLOCATE (p_theta_levels_old(tdims%i_start:tdims%i_end, &
                                 tdims%j_start:tdims%j_end, &
                                             1:tdims%k_end))

    ALLOCATE (VarTemp(field_size))

    p_theta_levels_old(:,:,:) = p_theta_levels(tdims%i_start:tdims%i_end, &
                                               tdims%j_start:tdims%j_end, &
                                                           1:tdims%k_end)

    VarTemp(:) = RESHAPE ( exner_theta_levels(tdims%i_start:tdims%i_end, &
                                              tdims%j_start:tdims%j_end, &
                                                          1:tdims%k_end) &
                         * theta             (tdims%i_start:tdims%i_end, &
                                              tdims%j_start:tdims%j_end, &
                                                          1:tdims%k_end) &
                         - qCL               (tdims%i_start:tdims%i_end, &
                                              tdims%j_start:tdims%j_end, &
                                                          1:tdims%k_end) &
                         * Lc/Cp                                         &
                         , (/field_size/) )

    IF (L_IAU_IncrementIce) THEN
      ALLOCATE (Var_BGT(field_size))
      Var_BGT(:) = RESHAPE ( exner_theta_levels(tdims%i_start:tdims%i_end, &
                                                tdims%j_start:tdims%j_end, &
                                                            1:tdims%k_end) &
                           * theta             (tdims%i_start:tdims%i_end, &
                                                tdims%j_start:tdims%j_end, &
                                                            1:tdims%k_end) &
                           , (/field_size/) )
    END IF

  END IF ! (.NOT.L_IAU_Add_qT_prime_to_q)

END IF ! (qTIncsThisTS)

!-------------------------------------------------------------------------------
! [5.4]: Fields required for PC2 calculations.
!-------------------------------------------------------------------------------

IF (.NOT. qTIncsThisTS .AND. &
    i_cld_vn==i_cld_pc2 .AND. (L_pc2_cond .OR. L_pc2_cfl .OR. L_pc2_cff)) THEN

  ALLOCATE (q_work  (tdims%i_start:tdims%i_end, &
                     tdims%j_start:tdims%j_end, &
                                 1:tdims%k_end))
  ALLOCATE (qcl_work(tdims%i_start:tdims%i_end, &
                     tdims%j_start:tdims%j_end, &
                                 1:tdims%k_end))
  ALLOCATE (qcf_in  (tdims%i_start:tdims%i_end, &
                     tdims%j_start:tdims%j_end, &
                                 1:tdims%k_end))
  ALLOCATE (p_work  (tdims%i_start:tdims%i_end, &
                     tdims%j_start:tdims%j_end, &
                                 1:tdims%k_end))
  ALLOCATE (t_work  (tdims%i_start:tdims%i_end, &
                     tdims%j_start:tdims%j_end, &
                                 1:tdims%k_end))

  ! Calculate increments to vapour, cloud, temperature and pressure:
  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        q_work  (i,j,k) = q  (i,j,k)
        qcl_work(i,j,k) = qcl(i,j,k)
        qcf_in  (i,j,k) = qcf(i,j,k)
        p_work  (i,j,k) = p_theta_levels    (i,j,k)
        t_work  (i,j,k) = exner_theta_levels(i,j,k) * theta(i,j,k)
      END DO
    END DO
  END DO

END IF

!-------------------------------------------------------------------------------
! [5.5]: Fields required for QLimits call.
!-------------------------------------------------------------------------------

IF (CallQLimitsThisTS) THEN

  ALLOCATE (q_in(tdims%i_start:tdims%i_end, &
                 tdims%j_start:tdims%j_end, &
                             1:tdims%k_end))

  q_in(:,:,:) = q(tdims%i_start:tdims%i_end, &
                  tdims%j_start:tdims%j_end, 1: )

END IF

!-------------------------------------------------------------------------------
! [6]: Loop over IAU increment files, and perform direct field updates.
!-------------------------------------------------------------------------------

DO IncNum = 1, Num_IAU_incs

  ! Get increment number as a string:
  WRITE (IncNumStr,'(I4)') IncNum
  IncNumStr = ADJUSTL(IncNumStr)
  ! Cycle if there are no increments to add on this timestep;
  IF (          .NOT. IAU_incs(IncNum) % FieldsToAdd      .OR. &
      stepim(atmos_im) < IAU_incs(IncNum) % InsertionStartTS .OR. &
      stepim(atmos_im) > IAU_incs(IncNum) % InsertionEndTS) CYCLE

  ! Get and report weight:
  WeightNum = stepim(atmos_im) - IAU_incs(IncNum) % InsertionStartTS + 1
  Weight    = IAU_incs(IncNum) % Weights(WeightNum)
  IF (PrintStatus >= PrStatus_Normal) THEN
    WRITE(umMessage,'(A,ES17.10,A)')                                    &
      ' IAU: Increments to be added with weight ', Weight,              &
      ' for inc no.'// TRIM(IncNumStr)
    CALL umPrint(umMessage,src='iau')
  END IF

  FileOpen = .FALSE.

  !-----------------------------------------------------------------------------
  ! [6.1]: If required, allocate storage space for increments.
  !-----------------------------------------------------------------------------

  IF (     IAU_incs(IncNum) % StoreFields .AND. &
      .NOT. IAU_incs(IncNum) % FieldsStored) THEN

    IF (IAU_incs(IncNum) % StorePacked) THEN
      ALLOCATE (IAU_incs(IncNum) % flds32bit(IAU_incs(IncNum) % flds_len))
    ELSE
      ALLOCATE (IAU_incs(IncNum) % flds     (IAU_incs(IncNum) % flds_len))
    END IF

  END IF

  !-----------------------------------------------------------------------------
  ! [6.2]: Loop over fields in the increment file.
  !-----------------------------------------------------------------------------

  flds_StartAddr = 1 ! Field start address in flds or flds32bit array

  ReadIAUField_FirstCallThisInc = .TRUE.

  DO FieldNum = 1, IAU_incs(IncNum) % Len2Lookup

    Code    = IAU_incs(IncNum) % Lookup(item_code, FieldNum)
    k       = IAU_incs(IncNum) % Lookup(lblev,     FieldNum)
    buflen  = IAU_incs(IncNum) % Lookup(lblrec,    FieldNum)

    ! Get local field length:
    LocFldLen = 0
    DO i = 1, IAU_NumFldCodes
      IF (IAU_FldCodes(i) == Code) LocFldLen = IAU_LocFldLens(i)
    END DO

    IF (LocFldLen == 0) CYCLE ! Unused field

    !---------------------------------------------------------------------------
    ! [6.2.1]: Obtain increment field, and save for next call if necessary.
    !---------------------------------------------------------------------------

    IF (IAU_incs(IncNum) % FieldsStored) THEN ! Field already in memory

      s_addr = flds_StartAddr
      e_addr = flds_StartAddr + LocFldLen - 1

      IF (IAU_incs(IncNum) % StorePacked) THEN
        ! Unpack field:
        CALL expand21 ( LocFldLen,                                   & ! in
                        IAU_incs(IncNum) % flds32bit(s_addr:e_addr), & ! in
                        IAUIncFld(1:LocFldLen) )                       ! out
      ELSE
        IAUIncFld(1:LocFldLen) = IAU_incs(IncNum) % flds(s_addr:e_addr)
      END IF

    ELSE ! Read field from file

      FilePath = IAU_incs(IncNum) % FilePath

      IF (.NOT. FileOpen) THEN
        CALL assign_file_unit(FilePath, IAU_unit, handler="portio")
        CALL model_file_open (IAU_unit, FilePath, read_write=ioOpenReadOnly)
        FileOpen = .TRUE.
      END IF

      CALL ReadIAUField (                                &
                          IAU_unit,                      & ! in
                          IAU_incs(IncNum) % FixHd,      & ! in
                          IAU_incs(IncNum) % Len1Lookup, & ! in
                          IAU_incs(IncNum) % Len2Lookup, & ! in
                          IAU_incs(IncNum) % Lookup,     & ! inout
                          IncNum,                        & ! in
                          FieldNum,                      & ! in
                          LocFldLen,                     & ! in
                          ReadIAUField_FirstCallThisInc, & ! in
                          IAUIncFld(1:LocFldLen) )         ! out

      ReadIAUField_FirstCallThisInc = .FALSE.

      IF (IAU_incs(IncNum) % StoreFields) THEN

        s_addr = flds_StartAddr
        e_addr = flds_StartAddr + LocFldLen - 1

        IF (IAU_incs(IncNum) % StorePacked) THEN
          CALL pack21 ( LocFldLen,              &                     ! in
                        IAUIncFld(1:LocFldLen), &                     ! in
                        IAU_incs(IncNum) % flds32bit(s_addr:e_addr) ) ! out
        ELSE
          IAU_incs(IncNum) % flds(s_addr:e_addr) = IAUIncFld(1:LocFldLen)
        END IF

      END IF

    END IF

    !---------------------------------------------------------------------------
    ! [6.2.2]: Add the increment to the relevant model field, applying limits
    !          where necessary.
    !---------------------------------------------------------------------------

    addr = 1

    ! u:
    IF (Code == 2 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          u(i,j,k) = u(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! v:
    IF (Code == 3 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          v(i,j,k) = v(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! w:
    IF ( Code == 150 .AND. &
         (k == 9999 .OR. (k >= 1 .AND. k <= model_levels)) ) THEN
      IF (k == 9999) k = 0
      DO j = wdims%j_start, wdims%j_end
        DO i = wdims%i_start, wdims%i_end
          w(i,j,k) = w(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! theta:
    IF (Code == 4 .AND. k >= 1 .AND. k <= model_levels) THEN
      ! ** Want to remove the following extra test at some point:
      IF (.NOT. (L_IAU_IgnoreTopLevThetaIncs .AND. k == model_levels)) THEN

        ! Apply limit to size of upper-level theta increments?
        IF (L_IAU_LimitUpperThetaIncs) THEN

          ! Shorten variable names:
          pBound = IAU_LimitUpperThetaIncs_pBound
          maxInc = IAU_LimitUpperThetaIncs_maxInc

          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end

              IF (p(i,j,k) < pBound .AND. &
                  ABS(Weight * IAUIncFld(addr)) > maxInc) THEN
                WRITE(umMessage,'(A,ES11.4,A,I5,A,ES11.4,A)')           &
                      'IAU: Theta inc (', Weight * IAUIncFld(addr),     &
                      ',) restricted at level ', k, ' pressure ',       &
                      p(i,j,k), ' for inc no.'// TRIM(IncNumStr)
                CALL umPrint(umMessage,src='iau')

                IF (Weight * IAUIncFld(addr) < 0) THEN
                  theta(i,j,k) = theta(i,j,k) - maxInc
                ELSE
                  theta(i,j,k) = theta(i,j,k) + maxInc
                END IF
              ELSE
                theta(i,j,k) = theta(i,j,k) + Weight * IAUIncFld(addr)
              END IF
              addr = addr + 1
            END DO
          END DO

        ELSE

          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              theta(i,j,k) = theta(i,j,k) + Weight * IAUIncFld(addr)
              addr = addr + 1
            END DO
          END DO

        END IF ! (L_IAU_LimitUpperThetaIncs)

      END IF
    END IF

    ! rho:
    IF (Code == 253 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          rho(i,j,k) = rho(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! exner:
    IF (Code == 255 .AND. k >= 1 .AND. k <= model_levels+1) THEN
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          exner(i,j,k) = exner(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! p:
    IF (Code == 407 .AND. k >= 1 .AND. k <= model_levels+1) THEN
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          p(i,j,k) = p(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! q:
    IF (Code == 10 .AND. k >= 1 .AND. k <  model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q(i,j,k) = q(i,j,k) + Weight * IAUIncFld(addr)
          ! Apply lower limit:
          IF (q(i,j,k) < q_min) q(i,j,k) = q_min
          addr = addr + 1
        END DO
      END DO
    END IF

    ! qCL:
    IF (Code == 254 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qCL(i,j,k) = qCL(i,j,k) + Weight * IAUIncFld(addr)
          ! Check for non-positive qCLs:
          IF (qCL(i,j,k) <= 0.0) THEN
            qCL(i,j,k) = 0.0
            cloud_fraction_liquid(i,j,k) = 0.0
            IF (qCF(i,j,k) <= CC_tol) THEN
              area_cloud_fraction(i,j,k) = 0.0
              bulk_cloud_fraction(i,j,k) = 0.0
            END IF
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF

    ! qCF:
    IF (Code == 12 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qCF(i,j,k) = qCF(i,j,k) + Weight * IAUIncFld(addr)
          ! Check for non-positive qCFs:
          IF (qCF(i,j,k) <= 0.0) THEN
            qCF(i,j,k) = 0.0
            cloud_fraction_frozen(i,j,k) = 0.0
            IF (qCL(i,j,k) <= CC_tol) THEN
              area_cloud_fraction(i,j,k) = 0.0
              bulk_cloud_fraction(i,j,k) = 0.0
            END IF
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF

    ! qT:
    IF ((Code == 16207 .OR. Code == 18001) .AND. &
        k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qT_plus(i,j,k) = qT_plus(i,j,k) + Weight * IAUIncFld(addr)
          addr = addr + 1
          IF (L_IAU_ApplyQTCorrections) THEN
            ! Apply lower limit:
            IF (qT_plus(i,j,k) < q_min) qT_plus(i,j,k) = q_min
          END IF
        END DO
      END DO
    END IF

    ! Aerosol:
    IF (Code == 90 .AND. k >= 1 .AND. k <= bl_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          murk(i,j,k) = murk(i,j,k) + Weight * IAUIncFld(addr)
          ! Apply upper (AEROMAX) and lower (AERO0) limits:
          murk(i,j,k) = MIN(MAX(aero0,murk(i,j,k)),aeromax)
          addr = addr + 1
        END DO
      END DO
    END IF

    ! dust_div1:
    IF (Code == 431 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div1(i,j,k) = dust_div1(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div1(i,j,k) <= 0.0) THEN
            dust_div1(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF
    ! dust_div2:
    IF (Code == 432 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div2(i,j,k) = dust_div2(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div2(i,j,k) <= 0.0) THEN
            dust_div2(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF
    ! dust_div3:
    IF (Code == 433 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div3(i,j,k) = dust_div3(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div3(i,j,k) <= 0.0) THEN
            dust_div3(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF
    ! dust_div4:
    IF (Code == 434 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div4(i,j,k) = dust_div4(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div4(i,j,k) <= 0.0) THEN
            dust_div4(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF
    ! dust_div5:
    IF (Code == 435 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div5(i,j,k) = dust_div5(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div5(i,j,k) <= 0.0) THEN
            dust_div5(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF
    ! dust_div6:
    IF (Code == 436 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dust_div6(i,j,k) = dust_div6(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive dust
          IF (dust_div6(i,j,k) <= 0.0) THEN
            dust_div6(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF

    ! Ozone tracer:
    IF (Code == 480 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ozone_tracer(i,j,k) = ozone_tracer(i,j,k) + Weight * IAUIncFld(addr)
          ! check for non-positive ozone tracer:
          IF (ozone_tracer(i,j,k) <= 0.0) THEN
            ozone_tracer(i,j,k) = 0.0
            IF (L_IAU_SetOzoneMin) ozone_tracer(i,j,k) = oz_min
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF

    ! SfcTemp: SST (sea) & TStar/TSoil (land)
    IF (Code == 24 .AND. k == 9999) THEN

      IF (PrintStatus >= PrStatus_Normal) THEN
        WRITE(umMessage,'(A)')                                          &
              'IAU: Adding surface temp perts to TStar + TSoil(1)'
        CALL umPrint(umMessage,src='iau')
      END IF
      ALLOCATE (t2_inc(theta_field_size))

      DO i = 1, theta_field_size
        t2_inc(i) = Weight * IAUIncFld(addr)
        addr = addr + 1
      END DO

      ! Add perts to TStar (Land and Sea)
      DO i = 1, theta_field_size
        TStar(i) = TStar(i) + t2_inc(i)
      END DO
      ! Add perts to tstar_anom if SST anomaly field is being used,
      ! otherwise anomaly update in REPLANCA will overwrite pert values
      IF (l_sstanom) THEN
        DO i = 1, theta_field_size
          tstar_anom(i) = tstar_anom(i) + t2_inc(i)
        END DO
      END IF

      ! Assuming an uncertainty in the snow depth, we perturb TSoil and
      ! TStar_tiles up to double the depth of the threshold used for
      ! analysis incr: i.e. 0.1 kg/m**2 instead of 0.05 kg/m**2

      ! Note on the multilayer snow scheme. The above threshold is
      ! roughly equivalent to a depth of 1 mm of fresh snow. This is
      ! sufficiently small that is seems unnecessary to make any
      ! extra adjustments for the multilayer scheme at this point.
      ! However, should the treatment of perturbations over snow
      ! be reconsidered, thought should be given to the updating
      ! of Tsnow in a manner similar to that used for assimilation
      ! of increments in section 7.8 below.

      ! Add sfctemp perts to TSoil (land)
      DO i = 1, land_field

        IF ( snow_depth(land_index(i)) <= 0.1 ) THEN
          TSoil(i, 1) = TSoil(i, 1) + t2_inc(land_index(i))
        END IF

      END DO

      ! Also add increments to TStar_tile:
      DO tile_num = 1, ntiles
        DO i = 1, land_field

          IF (nsmax == 0) THEN
            IF ( snow_depth(land_index(i)) <= 0.1 ) THEN
              TStar_tile(i,tile_num) = TStar_tile(i,tile_num) &
                                     + t2_inc(land_index(i))
            END IF
          ELSE
!           As for the 0-layer scheme, but to not allow the
!           accumulation of increments on empty tiles, which
!           will cause problems elsewhere.
!           The tile will match the surface type unless we are
!           aggregating, when there will be just one tile with
!           a fraction of 1.
            tile_frac = frac_typ(i,tile_num)
            IF (l_aggregate) tile_frac = 1.0
            IF ( ( snow_depth(land_index(i)) <= 0.1 ) .AND. &
                 ( tile_frac > 0.0 ) ) THEN
              TStar_tile(i,tile_num) = TStar_tile(i,tile_num) &
                                     + t2_inc(land_index(i))
            END IF
          END IF

        END DO
      END DO

      DEALLOCATE (t2_inc)

    END IF
    
    ! Deep Soil Temperature (st_levels)
    IF (Code == 20 .AND. k >= 1 .AND. k <= st_levels) THEN


      IF (PrintStatus >= PrStatus_Normal) THEN
        WRITE(umMessage,'(A)')                                          &
              'IAU: Adding deep soil temperature perts to TSoil'
        CALL umPrint(umMessage,src='iau')
      END IF


      ! Check whether data on land-points or full 2-D grid
      field_size = global_row_length * global_rows
      IF (buflen == field_size) THEN
        ALLOCATE (ts1_inc(theta_field_size))

        DO i = 1, theta_field_size
          ts1_inc(i) = Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO

        DO i = 1, land_field
          IF ( L_IAU_TSoilLandIceMask ) THEN
            ! Mask out TSoil perturbations on land-ice points, defined as
            ! points where the critical volumetric SMC is non-positive:
            IF ( vol_smc_crit(i) <= 0. ) ts1_inc(land_index(i)) = 0.0
          END IF
          IF ( snow_depth(land_index(i)) <= 0.05 ) THEN
            TSoil(i,k) = TSoil(i,k) + ts1_inc(land_index(i))
          ELSE IF (L_IAU_IncTSurf_Snow) THEN
              ! Limit the increment to prevent melting of snow
            IF (TSoil(i,k) < tm) TSoil(i,k) = MIN(TSoil(i,k)          &
                         + ts1_inc(land_index(i)), tm)
          END IF    
        END DO

        DEALLOCATE (ts1_inc)

      ELSE
        ! Land-points only
        DO i = 1, land_field
          IF ( L_IAU_TSoilLandIceMask ) THEN
            ! Mask out TSoil perturbations on land-ice points, defined as
            ! points where the critical volumetric SMC is non-positive:
            IF ( vol_smc_crit(i) <= 0.0 ) IAUIncFld(addr) = 0.0
          END IF
          IF ( snow_depth(land_index(i)) <= 0.05 ) THEN
            TSoil(i,k) = TSoil(i,k) + Weight * IAUIncFld(addr)
          ELSE IF (L_IAU_IncTSurf_Snow) THEN
              ! Limit the increment to prevent melting of snow
            IF (TSoil(i,k) < tm) TSoil(i,k) = MIN(TSoil(i,k)          &
                         + Weight * IAUIncFld(addr), tm)
          END IF    
          addr = addr + 1
        END DO

      END IF
    END IF

    ! Soil Moisture (sm_levels)
    IF (Code == 9 .AND. k >= 1 .AND. k <= sm_levels) THEN

      IF (PrintStatus >= PrStatus_Normal) THEN
        WRITE(umMessage,'(A)')                                          &
              'IAU: Adding soil moisture perts to smcl'
        CALL umPrint(umMessage,src='iau')
      END IF

      ! Check whether data on land-points or full 2-D grid
      field_size = global_row_length * global_rows
      IF (buflen == field_size) THEN
        ALLOCATE (s1_inc(theta_field_size))

        DO i = 1, theta_field_size
          s1_inc(i) = Weight * IAUIncFld(addr)
          addr = addr + 1
        END DO

        DO i = 1, land_field
          IF (L_IAU_SMCChecks) THEN
            ! Only apply smc perturbations at land points that
            ! (a) have a snow_depth <= 0.05 kg/m**2, and
            ! (b) are not land-ice points; i.e. critical volumetric smc > zero.
            IF (snow_depth(land_index(i)) <= 0.05 .AND. &
                vol_smc_crit(i) > 0.0) THEN
              smcl(i,k) = smcl(i,k) + s1_inc(land_index(i))
              ! Now apply limits so that
              ! wilting point * 0.1 <= volumetric SMC <= saturation.
              vol_smc = smcl(i,k) / (rho_water * dzsoil(k))
              IF (vol_smc > vol_smc_sat(i))      vol_smc = vol_smc_sat(i)
              IF (vol_smc < vol_smc_wilt(i)*0.1) vol_smc = vol_smc_wilt(i)*0.1
              smcl(i,k) = rho_water * dzsoil(k) * vol_smc
            END IF
          ELSE
            smcl(i,k) = smcl(i,k) + s1_inc(land_index(i))
          END IF
        END DO

        IF (ALLOCATED(s1_inc)) DEALLOCATE (s1_inc)

      ELSE
        ! Land-points only
        DO i = 1, land_field
          IF (L_IAU_SMCChecks) THEN
            ! Only apply smc perturbations at land points that
            ! (a) have a snow_depth <= 0.05 kg/m**2, and
            ! (b) are not land-ice points; i.e. critical volumetric smc > zero.
            IF (snow_depth(land_index(i)) <= 0.05 .AND. &
                vol_smc_crit(i) > 0.0) THEN
              smcl(i,k) = smcl(i,k) + Weight * IAUIncFld(addr)
              ! Now apply limits so that
              ! wilting point * 0.1 <= volumetric SMC <= saturation.
              vol_smc = smcl(i,k) / (rho_water * dzsoil(k))
              IF (vol_smc > vol_smc_sat(i))      vol_smc = vol_smc_sat(i)
              IF (vol_smc < vol_smc_wilt(i)*0.1) vol_smc = vol_smc_wilt(i)*0.1
              smcl(i,k) = rho_water * dzsoil(k) * vol_smc
            END IF
          ELSE
            smcl(i,k) = smcl(i,k) + Weight * IAUIncFld(addr)
          END IF
          addr = addr + 1
        END DO

      END IF

    END IF

    ! SO2:
    IF (Code == stashcode_SO2 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          SO2(i,j,k) = SO2(i,j,k) + Weight * IAUIncFld(addr) 
          ! Remove negatives:
          IF (SO2(i,j,k) < 0.0) THEN
            SO2(i,j,k) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF  

    ! UKCA tracer O3:
    IF (Code == stashcode_O3 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          tracer_ukca(i,j,k,ukca_O3_index) = tracer_ukca(i,j,k,ukca_O3_index) &
                                           + Weight * IAUIncFld(addr)
          ! Remove negatives:
          IF (tracer_ukca(i,j,k,ukca_O3_index) < 0.0) THEN
            tracer_ukca(i,j,k,ukca_O3_index) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF

    ! UKCA tracer NO2:
    IF (Code == stashcode_NO2 .AND. k >= 1 .AND. k <= model_levels) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          tracer_ukca(i,j,k,ukca_NO2_index) =                        &
                              tracer_ukca(i,j,k,ukca_NO2_index)      &
                              + Weight * IAUIncFld(addr)
          ! Remove negatives:
          IF (tracer_ukca(i,j,k,ukca_NO2_index) < 0.0) THEN
            tracer_ukca(i,j,k,ukca_NO2_index) = 0.0
          END IF
          addr = addr + 1
        END DO
      END DO
    END IF

    !---------------------------------------------------------------------------
    ! [6.2.3]: Update start address in flds or flds32bit array.
    !---------------------------------------------------------------------------

    flds_StartAddr = flds_StartAddr + LocFldLen

  END DO ! FieldNum

  !-----------------------------------------------------------------------------
  ! [6.3]: Set FieldsStored flag.
  !-----------------------------------------------------------------------------

  IF (IAU_incs(IncNum) % StoreFields) IAU_incs(IncNum) % FieldsStored = .TRUE.

  !-----------------------------------------------------------------------------
  ! [6.4]: Close increment file.
  !-----------------------------------------------------------------------------

  IF (FileOpen) THEN
    CALL model_file_close (IAU_unit, FilePath, delete=ioNoDelete, error=ICode)
    CALL release_file_unit (IAU_unit, handler="portio")
    IF (ICode > 0) THEN
      CMessage(1:80) = 'Error closing IAU increment no.'// TRIM(IncNumStr) //':'
      CMessage(81:)  = TRIM(FilePath)
      CALL EReport (RoutineName, ICode, CMessage)
    END IF
    FileOpen = .FALSE.
  END IF

END DO ! IncNum

DEALLOCATE(IAUIncFld)

!-------------------------------------------------------------------------------
! [7]: Perform indirect field updates.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! [7.1]: exner/p.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcExnerIncs) THEN
  ! Calculate exner from p:
  exner(1:row_length,1:rows,:) = ( p(1:row_length,1:rows,:) &
                                 * p_zero_recip             &
                                 )**kappa
ELSE                          ! Calculate p from exner
  ! Calculate p from exner:
  p(1:row_length,1:rows,:) =                                &
            p_zero * exner(1:row_length,1:rows,:)**recip_kappa
END IF

! Update pstar:
CALL Calc_P_star ( r_theta_levels, r_rho_levels,   & ! in
                   p, rho,                         & ! in
                   row_length, rows, model_levels, & ! in
                   offx, offy, halo_i, halo_j,     & ! in
                   p_star )                          ! out

! Update exner on theta levels:
CALL Calc_Exner_at_theta ( r_theta_levels, r_rho_levels,   & ! in
                           exner,                          & ! in
                           row_length, rows, model_levels, & ! in
                           offx, offy, halo_i, halo_j,     & ! in
                           exner_theta_levels,             & ! out
                           l_include_halos )                 ! in

! Get p on theta levels from exner on theta levels:
CALL Calc_P_from_Exner ( p_theta_levels,                 & ! out
                         row_length, rows,               & ! in
                         tdims_s%k_len,                  & ! in
                         offx, offy,                     & ! in
                         exner_theta_levels, l_include_halos)      ! in

!-------------------------------------------------------------------------------
! [7.2]: thetaV.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcThetaIncs .OR. L_IAU_CalcRhoIncs) THEN

  ALLOCATE (thetaV_new(tdims%i_start:tdims%i_end, &
                       tdims%j_start:tdims%j_end, &
                                   1:tdims%k_end))

  ! Add hydrostatic thetaV increments:
  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        IF (k /= model_levels) THEN
          delta_r = r_rho_levels(i,j,k+1) &
                  - r_rho_levels(i,j,k  )
        ELSE
          delta_r = 2.0                     &
                  * ( r_theta_levels(i,j,k) &
                    - r_rho_levels  (i,j,k) &
                    )
        END IF

        thetaV_new(i,j,k) = thetaV_in(i,j,k)       &
                          - g/cp * delta_r         &
                          * ( 1.0                  &
                            / ( exner    (i,j,k+1) &
                              - exner    (i,j,k  ) &
                            )                      &
                            - 1.0                  &
                            / ( exner_in(i,j,k+1)  &
                              - exner_in(i,j,k  )  &
                            ) )

      END DO
    END DO
  END DO

END IF

!-------------------------------------------------------------------------------
! [7.3]: rho.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcRhoIncs) THEN

  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        ! Interpolate old and new thetaVs to rho points:
        IF (k /= 1) THEN

          Weight_1 = r_rho_levels  (i,j,k  ) &
                   - r_theta_levels(i,j,k-1)
          Weight_2 = r_theta_levels(i,j,k  ) &
                   - r_rho_levels  (i,j,k  )

          thetaV_in_int  = ( Weight_1 * thetaV_in (i,j,k  ) &
                           + Weight_2 * thetaV_in (i,j,k-1) &
                           )                                &
                         / ( Weight_1 + Weight_2 )

          thetaV_new_int = ( Weight_1 * thetaV_new(i,j,k  ) &
                           + Weight_2 * thetaV_new(i,j,k-1) &
                           )                                &
                         / ( Weight_1 + Weight_2 )

        ELSE

          thetaV_in_int  = thetaV_in (i,j,k)
          thetaV_new_int = thetaV_new(i,j,k)

        END IF

        Term_1 = exner   (i,j,k)**exner_exp / thetaV_new_int
        Term_2 = exner_in(i,j,k)**exner_exp / thetaV_in_int

        rho(i,j,k) = rho(i,j,k)             &
                   + ( Term_1 - Term_2 )    &
                   * r_rho_levels(i,j,k)**2 &
                   * p_zero / (kappa * cp)

      END DO
    END DO
  END DO

END IF

!-------------------------------------------------------------------------------
! [7.4]: theta.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcThetaIncs) THEN

  theta(tdims%i_start:tdims%i_end,          &
        tdims%j_start:tdims%j_end,          &
                    1:tdims%k_end) =        &
             thetaV_new(:,:,1:tdims%k_end)  &
          / ( 1.0 + C_Virtual               &
             * q(tdims%i_start:tdims%i_end, &
                 tdims%j_start:tdims%j_end, &
                             1:tdims%k_end) &
             )

  theta(tdims%i_start:tdims%i_end, &
        tdims%j_start:tdims%j_end, &
                      tdims%k_end+1:) = thetaV_new(:,:,tdims%k_end+1:)

END IF

!-------------------------------------------------------------------------------
! [7.5]: Deallocate work arrays.
!-------------------------------------------------------------------------------

IF (L_IAU_CalcThetaIncs .OR. L_IAU_CalcRhoIncs) THEN
  DEALLOCATE (exner_in )
  DEALLOCATE (thetaV_in)
  DEALLOCATE (thetaV_new)
END IF

!-------------------------------------------------------------------------------
! [7.6]: Humidity and cloud updates via qT increments.
!-------------------------------------------------------------------------------

IF (qTIncsThisTS .AND. .NOT. L_IAU_Add_qT_prime_to_q) THEN

  ALLOCATE (AreaDivBulk(tdims%i_start:tdims%i_end, &
                        tdims%j_start:tdims%j_end, &
                                    1:tdims%k_end))

  ALLOCATE (cloud_fraction_liquid_0   (field_size))
  ALLOCATE (cloud_fraction_liquid_plus(field_size))
  ALLOCATE (qCL_0                     (field_size))
  ALLOCATE (qCL_plus                  (field_size))
  ALLOCATE (Var_qT                    (field_size))

  !----------------------------------------------------------------------------
  ! [7.6.1]: obtain initial and incremented values from Var scheme
  !-----------------------------------------------------------------------------

  ! Use Var cloud diagnostic to obtain  initial T, from qT & Tl
  Var_qT  = RESHAPE(qT(tdims%i_start:tdims%i_end, &
                       tdims%j_start:tdims%j_end, &
                                   1:tdims%k_end), (/field_size/))

  CALL Var_DiagCloud(                                       &
      field_size,                                           &
      RESHAPE(p_theta_levels_old, (/field_size/)),          &
      RESHAPE(SPREAD(spread(                                &
              rhcrit(1:model_levels),1,rows),1,row_length), &
              (/field_size/)),                              &
      L_IAU_IncrementIce,                                   &
      Var_qT,                                               &
      CMessage,                                             &
      ICode,                                                &
      tl=VarTemp)

  IF (ICode /= 0) CALL EReport (RoutineName, ICode, CMessage)

  ! (VarTemp now holds initial T values)

  IF (L_IAU_IncrementIce) THEN

    ALLOCATE (cloud_fraction_frozen_0   (field_size))
    ALLOCATE (cloud_fraction_frozen_plus(field_size))
    ALLOCATE (qCF_0                     (field_size))
    ALLOCATE (qCF_max                   (field_size))
    ALLOCATE (qCF_plus                  (field_size))

    ! Use Var cloud diagnostic to obtain  initial qCL,qCF,Cl & Cf
    qcl_0 = RESHAPE(qCL (tdims%i_start:tdims%i_end, &
                         tdims%j_start:tdims%j_end, &
                                     1:tdims%k_end), (/field_size/))
    qcf_0 = RESHAPE(qCF (tdims%i_start:tdims%i_end, &
                         tdims%j_start:tdims%j_end, &
                                     1:tdims%k_end), (/field_size/))

    CALL Var_DiagCloud(                                       &
        field_size,                                           &
        RESHAPE(p_theta_levels_old, (/field_size/)),          &
        RESHAPE(SPREAD(spread(                                &
                rhcrit(1:model_levels),1,rows),1,row_length), &
                (/field_size/)),                              &
        L_IAU_IncrementIce,                                   &
        Var_qT,                                               &
        CMessage,                                             &
        ICode,                                                &
        cl    = cloud_fraction_liquid_0,                      &
        qCL   = qCL_0,                                        &
        BGqCL = RESHAPE(qCL (tdims%i_start:tdims%i_end,       &
                             tdims%j_start:tdims%j_end,       &
                                         1:tdims%k_end),      &
                        (/field_size/)),                      &
        qCF   = qCF_0,                                        &
        cf    = cloud_fraction_frozen_0,                      &
        BGqCF = RESHAPE(qCF (tdims%i_start:tdims%i_end,       &
                             tdims%j_start:tdims%j_end,       &
                                         1:tdims%k_end),      &
                        (/field_size/)),                      &
        BG_t  = Var_BGT,                                      &
        t     = VarTemp)

  ELSE

    ! Use Var cloud diagnostic to obtain  initial qCL & Cl

    CALL Var_DiagCloud(                                       &
        field_size,                                           &
        RESHAPE(p_theta_levels_old, (/field_size/)),          &
        RESHAPE(SPREAD(spread(                                &
                rhcrit(1:model_levels),1,rows),1,row_length), &
                (/field_size/)),                              &
        L_IAU_IncrementIce,                                   &
        Var_qT,                                               &
        CMessage,                                             &
        ICode,                                                &
        cl    = cloud_fraction_liquid_0,                      &
        qCL   = qCL_0,                                        &
        t     = VarTemp)

  END IF

  IF (ICode /= 0) CALL EReport (RoutineName, ICode, CMessage)

  ! Assign VarTemp to hold incremented T values:
  IF (L_IAU_IncrementIce) THEN

    ! Assign Var_BGT to store the T associated with qcl_0
    ! Use qCL_Plus as a temporary variable
    ! This should be done in the qcl / Cl version too but
    !   will go in later as a bug-fix so that bit-comparison
    !   can be maintained when the code for qcf incrementing is
    !   implemented but switched off
    qCL_Plus = VarTemp
    VarTemp  = VarTemp + (                              &
         RESHAPE(                                       &
         exner_theta_levels(tdims%i_start:tdims%i_end,  &
                            tdims%j_start:tdims%j_end,  &
                                        1:tdims%k_end)  &
         * theta(tdims%i_start:tdims%i_end,             &
                            tdims%j_start:tdims%j_end,  &
                                        1:tdims%k_end), &
        (/field_size/)) - Var_BGT)
    Var_BGT  = qCL_Plus

  ELSE

    VarTemp = RESHAPE(                                 &
         exner_theta_levels(tdims%i_start:tdims%i_end, &
                            tdims%j_start:tdims%j_end, &
                                        1:tdims%k_end) &
         * theta(tdims%i_start:tdims%i_end,            &
                 tdims%j_start:tdims%j_end,            &
                             1:tdims%k_end),           &
        (/field_size/))

  END IF

  IF (L_IAU_IncrementIce) THEN

    ! Calculate maximum qcf values (for use later)
    qCL_Plus = 0.0
    qCF_Plus = 0.0
    Var_qT = RESHAPE(qT_plus(tdims%i_start:tdims%i_end,  &
                             tdims%j_start:tdims%j_end,  &
                                         1:tdims%k_end)  &
                     +  qCF (tdims%i_start:tdims%i_end,  &
                             tdims%j_start:tdims%j_end,  &
                                         1:tdims%k_end), &
                     (/field_size/))
    qCF_max = Var_qT

    CALL Var_DiagCloud(                                       &
        field_size,                                           &
        RESHAPE(p_theta_levels                                &
                (tdims%i_start:tdims%i_end,                   &
                 tdims%j_start:tdims%j_end,                   &
                             1:tdims%k_end),                  &
                (/field_size/)),                              &
        RESHAPE(SPREAD(spread(                                &
                rhcrit(1:model_levels),1,rows),1,row_length), &
                (/field_size/)),                              &
        L_IAU_IncrementIce,                                   &
        Var_qT,                                               &
        CMessage,                                             &
        ICode,                                                &
        qCL    = qCL_Plus,                                    &
        BGqCL  = RESHAPE(qCL(tdims%i_start:tdims%i_end,       &
                             tdims%j_start:tdims%j_end,       &
                                         1:tdims%k_end),      &
                         (/field_size/)),                     &
        qCF    = qCF_Plus,                                    &
        qCFmax = qCF_max,                                     &
        BGqCF  = RESHAPE(qCF(tdims%i_start:tdims%i_end,       &
                             tdims%j_start:tdims%j_end,       &
                                         1:tdims%k_end),      &
                         (/field_size/)),                     &
        BG_t   = Var_BGT,                                     &
        t      = VarTemp)

    IF (ICode /= 0) CALL EReport (RoutineName, ICode, CMessage)

    qCL_Plus  = qCL_0
    qCF_Plus  = qCF_0
    Var_qT = RESHAPE(qT     (tdims%i_start:tdims%i_end,  &
                             tdims%j_start:tdims%j_end,  &
                                         1:tdims%k_end)  &
                     +  qCF (tdims%i_start:tdims%i_end,  &
                             tdims%j_start:tdims%j_end,  &
                                         1:tdims%k_end), &
                     (/field_size/))

    WHERE (qCF_max <= 0.0)

      qcf_max = Var_qT

    END WHERE

  END IF

  Var_qT = RESHAPE(qT_plus(tdims%i_start:tdims%i_end, &
                        tdims%j_start:tdims%j_end,    &
                                    1:tdims%k_end),   &
                   (/field_size/))

  IF (L_IAU_IncrementIce) THEN

    ! Use Var cloud diagnostic to obtain  final qCL,qCF,Cl & Cf
    CALL Var_DiagCloud(                                       &
        field_size,                                           &
        RESHAPE(p_theta_levels                                &
                (tdims%i_start:tdims%i_end,                   &
                 tdims%j_start:tdims%j_end,                   &
                             1:tdims%k_end),                  &
                (/field_size/)),                              &
        RESHAPE(SPREAD(spread(                                &
                rhcrit(1:model_levels),1,rows),1,row_length), &
                (/field_size/)),                              &
        L_IAU_IncrementIce,                                   &
        Var_qT,                                               &
        CMessage,                                             &
        ICode,                                                &
        qCL    = qCL_Plus,                                    &
        BGqCL  = RESHAPE(qCL(tdims%i_start:tdims%i_end,       &
                             tdims%j_start:tdims%j_end,       &
                                         1:tdims%k_end),      &
                         (/field_size/)),                     &
        qCF    = qCF_Plus,                                    &
        BGqCF  = RESHAPE(qCF(tdims%i_start:tdims%i_end,       &
                             tdims%j_start:tdims%j_end,       &
                                         1:tdims%k_end),      &
                         (/field_size/)),                     &
        qCFmax = qCF_max,                                     &
        cl     = cloud_fraction_liquid_plus,                  &
        cf     = cloud_fraction_frozen_plus,                  &
        BG_t   = Var_BGT,                                     &
        t      = VarTemp)

  ELSE

    ! Use Var cloud diagnostic to obtain  final qCL & Cl
    CALL Var_DiagCloud(                                       &
        field_size,                                           &
        RESHAPE(p_theta_levels                                &
                (tdims%i_start:tdims%i_end,                   &
                 tdims%j_start:tdims%j_end,                   &
                             1:tdims%k_end),                  &
                (/field_size/)),                              &
        RESHAPE(SPREAD(spread(                                &
                rhcrit(1:model_levels),1,rows),1,row_length), &
                (/field_size/)),                              &
        L_IAU_IncrementIce,                                   &
        Var_qT,                                               &
        CMessage,                                             &
        ICode,                                                &
        cl     = cloud_fraction_liquid_plus,                  &
        qCL    = qCL_Plus,                                    &
        t      = VarTemp)

  END IF

  IF (ICode /= 0) CALL EReport (RoutineName, ICode, CMessage)

  !-----------------------------------------------------------------------------
  ! [7.6.2]: apply increments
  !-----------------------------------------------------------------------------

  ! qcl

  ALLOCATE (qcl_Inc(tdims%i_start:tdims%i_end, &
                    tdims%j_start:tdims%j_end, &
                                1:tdims%k_end))

  IF (L_IAU_ScaleCloud) THEN

    ! Scale qcl increments to be physical values
    ALLOCATE (initialScaling(tdims%i_start:tdims%i_end, &
                             tdims%j_start:tdims%j_end, &
                                         1:tdims%k_end))
    ALLOCATE (scalingCoeff  (tdims%i_start:tdims%i_end, &
                             tdims%j_start:tdims%j_end, &
                                         1:tdims%k_end))
    WHERE (RESHAPE(qCL_0,(/row_length,rows,model_levels/)) >     &
          qCL(tdims%i_start:tdims%i_end,                         &
              tdims%j_start:tdims%j_end,                         &
                          1:tdims%k_end) .AND.                   &
      RESHAPE(qCL_plus - qCL_0,(/row_length,rows,model_levels/)) &
      < 0.0)
      qcl_Inc = RESHAPE(qCL_plus - qCL_0, &
                        (/row_length, rows, model_levels/))
      WHERE (qCL(tdims%i_start:tdims%i_end, &
                 tdims%j_start:tdims%j_end, &
                             1:tdims%k_end) > 0.0)
        initialScaling = -qCL(tdims%i_start:tdims%i_end,  &
                              tdims%j_start:tdims%j_end,  &
                                          1:tdims%k_end)* &
            (1.0 - EXP(qcl_Inc /                          &
                       qCL(tdims%i_start:tdims%i_end,     &
                           tdims%j_start:tdims%j_end,     &
                                       1:tdims%k_end)))
        scalingCoeff = -qCL(tdims%i_start:tdims%i_end,             &
                            tdims%j_start:tdims%j_end,             &
                                        1:tdims%k_end)*            &
        (1.0-EXP(-RESHAPE(qCL_0,(/row_length,rows,model_levels/))/ &
        qCL(tdims%i_start:tdims%i_end,                             &
            tdims%j_start:tdims%j_end,                             &
                        1:tdims%k_end)))
        ! The following calculation should remain >0 so no error
        ! trap needed
        scalingCoeff =                                      &
            (scalingCoeff +                                 &
             qCL(tdims%i_start:tdims%i_end,                 &
                 tdims%j_start:tdims%j_end,                 &
                             1:tdims%k_end)) /              &
            (scalingCoeff +                                 &
             RESHAPE(qCL_0,(/row_length,rows,model_levels/)))
      ELSEWHERE
        initialScaling = 0.0
        scalingCoeff = 0.0
      END WHERE
      qcl_Inc = initialScaling - &
                (initialScaling - qcl_Inc) * scalingCoeff
    ELSEWHERE
      qcl_Inc = RESHAPE(qCL_plus - qCL_0, &
                        (/row_length, rows, model_levels/))
    END WHERE

  ELSE

    qcl_Inc = RESHAPE(qCL_plus - qCL_0, &
                      (/row_length, rows, model_levels/))

  END IF

  ! Update:
  qCL(tdims%i_start:tdims%i_end,   &
      tdims%j_start:tdims%j_end,   &
                  1:tdims%k_end) = &
  qCL(tdims%i_start:tdims%i_end,   &
      tdims%j_start:tdims%j_end,   &
                  1:tdims%k_end) + qcl_Inc(:,:,:)

  ! qcf
  IF (L_IAU_IncrementIce) THEN

    qCF(tdims%i_start:tdims%i_end,   &
        tdims%j_start:tdims%j_end,   &
                    1:tdims%k_end) = &
    qCF(tdims%i_start:tdims%i_end,   &
        tdims%j_start:tdims%j_end,   &
                    1:tdims%k_end) + &
        RESHAPE(qCF_plus - qCF_0,    &
                (/row_length, rows, model_levels/))

  END IF

  ! Store input values of area / bulk cloud fractions
  IF (i_cld_area == acf_cusack .OR. i_cld_area == acf_brooks) THEN

    WHERE (bulk_cloud_fraction(tdims%i_start:tdims%i_end, &
                               tdims%j_start:tdims%j_end, &
                                           1:tdims%k_end) &
           <= 0.0)

      AreaDivBulk = 1.0

    ELSEWHERE

      AreaDivBulk =                                      &
        area_cloud_fraction(tdims%i_start:tdims%i_end,   &
                            tdims%j_start:tdims%j_end,   &
                                        1:tdims%k_end) / &
        bulk_cloud_fraction(tdims%i_start:tdims%i_end,   &
                            tdims%j_start:tdims%j_end,   &
                                        1:tdims%k_end)

    END WHERE

  END IF

  ! Cl
  IF (L_IAU_ScaleCloud) THEN

    ! Scale Cl increments to be physical values
    ALLOCATE (Cl_Inc(tdims%i_start:tdims%i_end,  &
                     tdims%j_start:tdims%j_end,  &
                                 1:tdims%k_end))
    ALLOCATE (Cf_Inc(tdims%i_start:tdims%i_end,  &
                     tdims%j_start:tdims%j_end,  &
                                 1:tdims%k_end))
    WHERE (RESHAPE(cloud_fraction_liquid_0,              &
                  (/row_length,rows,model_levels/)) >    &
      cloud_fraction_liquid(tdims%i_start:tdims%i_end,   &
                            tdims%j_start:tdims%j_end,   &
                                        1:tdims%k_end)   &
      .AND. RESHAPE(cloud_fraction_liquid_plus -         &
                    cloud_fraction_liquid_0,             &
                    (/row_length, rows, model_levels/)) < 0.0)
      Cl_Inc = RESHAPE(cloud_fraction_liquid_plus - &
                       cloud_fraction_liquid_0,     &
                       (/row_length, rows, model_levels/))
      WHERE (cloud_fraction_liquid                 &
                       (tdims%i_start:tdims%i_end, &
                        tdims%j_start:tdims%j_end, &
                                    1:tdims%k_end) > 0.0)
        initialScaling = - cloud_fraction_liquid       &
                           (tdims%i_start:tdims%i_end, &
                            tdims%j_start:tdims%j_end, &
                                        1:tdims%k_end) &
                         * (1.0 - EXP(Cl_Inc /         &
                            cloud_fraction_liquid      &
                          (tdims%i_start:tdims%i_end,  &
                           tdims%j_start:tdims%j_end,  &
                                       1:tdims%k_end)))
        scalingCoeff =                                        &
            - cloud_fraction_liquid                           &
                       (tdims%i_start:tdims%i_end,            &
                        tdims%j_start:tdims%j_end,            &
                                    1:tdims%k_end) *          &
              (1.0 - EXP(-RESHAPE(cloud_fraction_liquid_0,    &
                          (/row_length,rows,model_levels/)) / &
              cloud_fraction_liquid                           &
                       (tdims%i_start:tdims%i_end,            &
                        tdims%j_start:tdims%j_end,            &
                                    1:tdims%k_end)))
        ! The following calculation should remain >0 so
        ! no error trap needed
        scalingCoeff =                                       &
            (scalingCoeff +                                  &
             cloud_fraction_liquid                           &
                       (tdims%i_start:tdims%i_end,           &
                        tdims%j_start:tdims%j_end,           &
                                    1:tdims%k_end)) /        &
            (scalingCoeff + RESHAPE(cloud_fraction_liquid_0, &
                            (/row_length,rows,model_levels/)))
      ELSEWHERE
        initialScaling = 0.0
        scalingCoeff = 0.0
      END WHERE
      Cl_Inc = initialScaling - (initialScaling - Cl_Inc) * &
               scalingCoeff
    ELSEWHERE(RESHAPE(cloud_fraction_liquid_0,            &
                      (/row_length,rows,model_levels/)) < &
              cloud_fraction_liquid                       &
              (tdims%i_start:tdims%i_end,                 &
               tdims%j_start:tdims%j_end,                 &
                           1:tdims%k_end) .AND.           &
              RESHAPE(cloud_fraction_liquid_plus -        &
                      cloud_fraction_liquid_0,            &
                      (/row_length, rows, model_levels/)) > 0.0)
      Cl_Inc = RESHAPE(cloud_fraction_liquid_plus - &
                       cloud_fraction_liquid_0,     &
                       (/row_length, rows, model_levels/))
      WHERE ((1.0 - cloud_fraction_liquid       &
                    (tdims%i_start:tdims%i_end, &
                     tdims%j_start:tdims%j_end, &
                                 1:tdims%k_end)) > 0.0)
        initialScaling = (1.0 - cloud_fraction_liquid   &
                          (tdims%i_start:tdims%i_end,   &
                           tdims%j_start:tdims%j_end,   &
                                       1:tdims%k_end))* &
                         (1.0 - EXP(-Cl_Inc / (1.0 -    &
                           cloud_fraction_liquid        &
                         (tdims%i_start:tdims%i_end,    &
                          tdims%j_start:tdims%j_end,    &
                                      1:tdims%k_end))))
        scalingCoeff = (1.0 - cloud_fraction_liquid                &
            (tdims%i_start:tdims%i_end,                            &
             tdims%j_start:tdims%j_end,                            &
                         1:tdims%k_end)) *                         &
            (1.0 - EXP(-(1.0 - RESHAPE(cloud_fraction_liquid_0,    &
                               (/row_length,rows,model_levels/)))/ &
            (1.0 - cloud_fraction_liquid                           &
                   (tdims%i_start:tdims%i_end,                     &
                    tdims%j_start:tdims%j_end,                     &
                                1:tdims%k_end))))
        scalingCoeff =                                  &
            ((1.0 - cloud_fraction_liquid               &
                    (tdims%i_start:tdims%i_end,         &
                     tdims%j_start:tdims%j_end,         &
                                 1:tdims%k_end)) -      &
             scalingCoeff) / ((1.0 -                    &
             RESHAPE(cloud_fraction_liquid_0,           &
                     (/row_length,rows,model_levels/))) &
             - scalingCoeff)
      ELSEWHERE
        initialScaling = 0.0
        scalingCoeff = 0.0
      END WHERE
      Cl_Inc = initialScaling + (Cl_Inc - initialScaling) * &
               scalingCoeff
    ELSEWHERE
      Cl_Inc = RESHAPE(cloud_fraction_liquid_plus - &
                       cloud_fraction_liquid_0,     &
                       (/row_length, rows, model_levels/))
    END WHERE
    ! Reuse cloud_fraction_liquid_0 to store initial values for
    ! use in Cb incrementing
    cloud_fraction_liquid_0 = RESHAPE(cloud_fraction_liquid &
                         (tdims%i_start:tdims%i_end,                  &
                          tdims%j_start:tdims%j_end,                  &
                                      1:tdims%k_end), (/field_size/))
    cloud_fraction_liquid(tdims%i_start:tdims%i_end,     &
                          tdims%j_start:tdims%j_end,     &
                                      1:tdims%k_end) =   &
        cloud_fraction_liquid(tdims%i_start:tdims%i_end, &
                              tdims%j_start:tdims%j_end, &
                                          1:tdims%k_end) &
        + Cl_Inc

  ELSE

    cloud_fraction_liquid(tdims%i_start:tdims%i_end,     &
                          tdims%j_start:tdims%j_end,     &
                                      1:tdims%k_end) =   &
        cloud_fraction_liquid(tdims%i_start:tdims%i_end, &
                              tdims%j_start:tdims%j_end, &
                                          1:tdims%k_end) &
        + RESHAPE(cloud_fraction_liquid_plus -           &
                  cloud_fraction_liquid_0,               &
                  (/row_length, rows, model_levels/))

  END IF

  ! Set Cl <= 1.0
  WHERE (                                              &
      cloud_fraction_liquid(tdims%i_start:tdims%i_end, &
                            tdims%j_start:tdims%j_end, &
                                        1:tdims%k_end) &
      > 1.0)

    cloud_fraction_liquid(tdims%i_start:tdims%i_end, &
                          tdims%j_start:tdims%j_end, &
                                      1:tdims%k_end) &
        = 1.0

  END WHERE

  ! Set qcl & Cl >= 0.0
  WHERE (                                              &
      cloud_fraction_liquid(tdims%i_start:tdims%i_end, &
                            tdims%j_start:tdims%j_end, &
                                        1:tdims%k_end) &
      <= 0.0                                           &
      .OR. qCL(tdims%i_start:tdims%i_end,              &
               tdims%j_start:tdims%j_end,              &
                           1:tdims%k_end) <= CC_tol)

    qCL(tdims%i_start:tdims%i_end, &
        tdims%j_start:tdims%j_end, &
                    1:tdims%k_end) = 0.0
    cloud_fraction_liquid(tdims%i_start:tdims%i_end, &
                          tdims%j_start:tdims%j_end, &
                                      1:tdims%k_end) &
        = 0.0

  END WHERE

  IF (L_IAU_IncrementIce) THEN

    ! Cf
    IF (L_IAU_ScaleCloud) THEN

      ! Scale Cf increments to be physical values
      WHERE (RESHAPE(cloud_fraction_frozen_0,              &
                    (/row_length,rows,model_levels/)) >    &
        cloud_fraction_frozen(tdims%i_start:tdims%i_end,   &
                              tdims%j_start:tdims%j_end,   &
                                          1:tdims%k_end)   &
        .AND. RESHAPE(cloud_fraction_frozen_plus -         &
                      cloud_fraction_frozen_0,             &
                      (/row_length,rows,model_levels/)) < 0.0)
        Cf_Inc = RESHAPE(cloud_fraction_frozen_plus - &
                         cloud_fraction_frozen_0,     &
                         (/row_length,rows,model_levels/))
        WHERE (cloud_fraction_frozen       &
               (tdims%i_start:tdims%i_end, &
                tdims%j_start:tdims%j_end, &
                            1:tdims%k_end) > 0.0)
          initialScaling = - cloud_fraction_frozen       &
                             (tdims%i_start:tdims%i_end, &
                              tdims%j_start:tdims%j_end, &
                                          1:tdims%k_end) &
                           * (1.0 - EXP(Cf_Inc /         &
                              cloud_fraction_frozen      &
                            (tdims%i_start:tdims%i_end,  &
                             tdims%j_start:tdims%j_end,  &
                                         1:tdims%k_end)))
          scalingCoeff =                                        &
              - cloud_fraction_frozen                           &
                (tdims%i_start:tdims%i_end,                     &
                 tdims%j_start:tdims%j_end,                     &
                             1:tdims%k_end) *                   &
                (1.0 - EXP(-RESHAPE(cloud_fraction_frozen_0,    &
                            (/row_length,rows,model_levels/)) / &
                cloud_fraction_frozen                           &
                (tdims%i_start:tdims%i_end,                     &
                 tdims%j_start:tdims%j_end,                     &
                             1:tdims%k_end)))
          ! The following calculation should remain >0 so
          ! no error trap needed
          scalingCoeff =                                       &
              (scalingCoeff +                                  &
               cloud_fraction_frozen                           &
               (tdims%i_start:tdims%i_end,                     &
                tdims%j_start:tdims%j_end,                     &
                            1:tdims%k_end)) /                  &
              (scalingCoeff + RESHAPE(cloud_fraction_frozen_0, &
                              (/row_length,rows,model_levels/)))
        ELSEWHERE
          initialScaling = 0.0
          scalingCoeff = 0.0
        END WHERE
        Cf_Inc = initialScaling - (initialScaling - Cf_Inc) * &
                 scalingCoeff
      ELSEWHERE(RESHAPE(cloud_fraction_frozen_0,            &
                        (/row_length,rows,model_levels/)) < &
                cloud_fraction_frozen                       &
                (tdims%i_start:tdims%i_end,                 &
                 tdims%j_start:tdims%j_end,                 &
                             1:tdims%k_end) .AND.           &
                RESHAPE(cloud_fraction_frozen_plus -        &
                        cloud_fraction_frozen_0,            &
                        (/row_length,rows,model_levels/)) > 0.0)
        Cf_Inc = RESHAPE(cloud_fraction_frozen_plus - &
                         cloud_fraction_frozen_0,     &
                         (/row_length,rows,model_levels/))
        WHERE ((1.0 - cloud_fraction_frozen       &
                      (tdims%i_start:tdims%i_end, &
                       tdims%j_start:tdims%j_end, &
                                   1:tdims%k_end)) > 0.0)
          initialScaling = (1.0 - cloud_fraction_frozen &
                      (tdims%i_start:tdims%i_end,       &
                       tdims%j_start:tdims%j_end,       &
                                   1:tdims%k_end))*     &
                           (1.0 - EXP(-Cf_Inc / (1.0 -  &
                             cloud_fraction_frozen      &
                      (tdims%i_start:tdims%i_end,       &
                       tdims%j_start:tdims%j_end,       &
                                   1:tdims%k_end))))
          scalingCoeff = (1.0 - cloud_fraction_frozen              &
              (tdims%i_start:tdims%i_end,                          &
               tdims%j_start:tdims%j_end,                          &
                           1:tdims%k_end)) *                       &
              (1.0 - EXP(-(1.0 -RESHAPE(cloud_fraction_frozen_0,   &
                           (/row_length,rows,model_levels/))) /    &
              (1.0 - cloud_fraction_frozen                         &
                     (tdims%i_start:tdims%i_end,                   &
                      tdims%j_start:tdims%j_end,                   &
                                  1:tdims%k_end))))
          scalingCoeff =                             &
              ((1.0 - cloud_fraction_frozen          &
                      (tdims%i_start:tdims%i_end,    &
                       tdims%j_start:tdims%j_end,    &
                                   1:tdims%k_end)) - &
               scalingCoeff) / ((1.0 -          &
               RESHAPE(cloud_fraction_frozen_0, &
               (/row_length,rows,model_levels/))) - scalingCoeff)
        ELSEWHERE
          initialScaling = 0.0
          scalingCoeff = 0.0
        END WHERE
        Cf_Inc = initialScaling + (Cf_Inc - initialScaling) * &
                 scalingCoeff
      ELSEWHERE
        Cf_Inc = RESHAPE(cloud_fraction_frozen_plus - &
                         cloud_fraction_frozen_0,     &
                         (/row_length,rows,model_levels/))
      END WHERE
      ! Reuse cloud_fraction_frozen_0 to store initial values for
      ! use in Cb incrementing
      cloud_fraction_frozen_0 = RESHAPE(cloud_fraction_frozen &
          (tdims%i_start:tdims%i_end,                         &
           tdims%j_start:tdims%j_end,                         &
                       1:tdims%k_end), (/field_size/))
      cloud_fraction_frozen(tdims%i_start:tdims%i_end,   &
                            tdims%j_start:tdims%j_end,   &
                                        1:tdims%k_end) = &
        cloud_fraction_frozen(tdims%i_start:tdims%i_end, &
        tdims%j_start:tdims%j_end,                       &
                    1:tdims%k_end)                       &
        + Cf_Inc

    ELSE

      cloud_fraction_frozen(tdims%i_start:tdims%i_end,     &
                            tdims%j_start:tdims%j_end,     &
                                        1:tdims%k_end) =   &
          cloud_fraction_frozen(tdims%i_start:tdims%i_end, &
                                tdims%j_start:tdims%j_end, &
                                            1:tdims%k_end) &
          + RESHAPE(cloud_fraction_frozen_plus -           &
           cloud_fraction_frozen_0,                        &
           (/row_length, rows, model_levels/))

    END IF

    ! Set Cf <= 1.0
    WHERE (                                            &
      cloud_fraction_frozen(tdims%i_start:tdims%i_end, &
                            tdims%j_start:tdims%j_end, &
                                        1:tdims%k_end) &
        > 1.0)

      cloud_fraction_frozen(tdims%i_start:tdims%i_end, &
                            tdims%j_start:tdims%j_end, &
                                        1:tdims%k_end) = 1.0

    END WHERE

    ! Set qcf & Cf >= 0.0
    WHERE (                                              &
        cloud_fraction_frozen(tdims%i_start:tdims%i_end, &
                              tdims%j_start:tdims%j_end, &
                                          1:tdims%k_end) &
        <= 0.0                                           &
        .OR. qCF(tdims%i_start:tdims%i_end,              &
                 tdims%j_start:tdims%j_end,              &
                             1:tdims%k_end) <= CC_tol)

      qCF(tdims%i_start:tdims%i_end, &
          tdims%j_start:tdims%j_end, &
                      1:tdims%k_end) = 0.0
      cloud_fraction_frozen(tdims%i_start:tdims%i_end, &
                            tdims%j_start:tdims%j_end, &
                                        1:tdims%k_end) &
          = 0.0

    END WHERE

  END IF

  ! Set other cloud variables to zero if qCL & qCF <= CC_tol
  WHERE (   qCL(tdims%i_start:tdims%i_end,           &
                tdims%j_start:tdims%j_end,           &
                            1:tdims%k_end) <= CC_tol &
      .AND. qCF(tdims%i_start:tdims%i_end,           &
                tdims%j_start:tdims%j_end,           &
                            1:tdims%k_end) <= CC_tol)

    qCF(tdims%i_start:tdims%i_end, &
        tdims%j_start:tdims%j_end, &
                    1:tdims%k_end) = 0.0
    cloud_fraction_frozen(tdims%i_start:tdims%i_end, &
                          tdims%j_start:tdims%j_end, &
                                      1:tdims%k_end) &
        = 0.0
    area_cloud_fraction(tdims%i_start:tdims%i_end, &
                        tdims%j_start:tdims%j_end, &
                                    1:tdims%k_end) &
        = 0.0
    bulk_cloud_fraction(tdims%i_start:tdims%i_end, &
                        tdims%j_start:tdims%j_end, &
                                    1:tdims%k_end) &
        = 0.0

  END WHERE

  IF (L_IAU_IncrementIce) THEN

    WHERE (   qCL(tdims%i_start:tdims%i_end,           &
                  tdims%j_start:tdims%j_end,           &
                              1:tdims%k_end) <= CC_tol &
        .AND. qCF(tdims%i_start:tdims%i_end,           &
                  tdims%j_start:tdims%j_end,           &
                              1:tdims%k_end) <= CC_tol)

      qCL(tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end, &
                                    1:tdims%k_end) = 0.0
      cloud_fraction_liquid(tdims%i_start:tdims%i_end, &
                            tdims%j_start:tdims%j_end, &
                                        1:tdims%k_end) = 0.0

    END WHERE

  END IF

  ! Apply stored area / bulk cloud fractions to output
  IF (L_IAU_ScaleCloud) THEN

    ! Scale Cb increments to be physical values
    ! Use Cb_inc and Cf_Inc as temporary storage for incremented
    ! and unincremented MIN(Cl+Cf, 1.0)
    ALLOCATE (Cb_Inc (tdims%i_start:tdims%i_end, &
                      tdims%j_start:tdims%j_end, &
                                  1:tdims%k_end))
    Cb_Inc = MIN(cloud_fraction_liquid               &
                       (tdims%i_start:tdims%i_end,   &
                        tdims%j_start:tdims%j_end,   &
                                    1:tdims%k_end) + &
                 cloud_fraction_frozen               &
                       (tdims%i_start:tdims%i_end,   &
                        tdims%j_start:tdims%j_end,   &
                                    1:tdims%k_end), 1.0)
    IF (L_IAU_IncrementIce) THEN
      Cf_Inc = RESHAPE(MIN(cloud_fraction_liquid_0 +     &
                           cloud_fraction_frozen_0,1.0), &
                       (/row_length,rows,model_levels/))
    ELSE
      Cf_Inc = MIN(RESHAPE(cloud_fraction_liquid_0,        &
                       (/row_length,rows,model_levels/)) + &
                       cloud_fraction_frozen               &
                       (tdims%i_start:tdims%i_end,         &
                        tdims%j_start:tdims%j_end,         &
                                    1:tdims%k_end), 1.0)
    END IF
    WHERE (Cf_Inc >                                       &
          bulk_cloud_fraction(tdims%i_start:tdims%i_end, &
                              tdims%j_start:tdims%j_end, &
                                          1:tdims%k_end) &
          .AND. Cb_Inc - Cf_Inc < 0.0)
      Cb_Inc = Cb_Inc - Cf_Inc
      WHERE (bulk_cloud_fraction         &
             (tdims%i_start:tdims%i_end, &
              tdims%j_start:tdims%j_end, &
                          1:tdims%k_end) > 0.0)
        initialScaling = - bulk_cloud_fraction           &
                           (tdims%i_start:tdims%i_end,   &
                            tdims%j_start:tdims%j_end,   &
                                        1:tdims%k_end) * &
                         (1.0 - EXP(Cb_Inc /             &
                                    bulk_cloud_fraction  &
                           (tdims%i_start:tdims%i_end,   &
                            tdims%j_start:tdims%j_end,   &
                                        1:tdims%k_end)))
        scalingCoeff  = - bulk_cloud_fraction            &
                             (tdims%i_start:tdims%i_end, &
                              tdims%j_start:tdims%j_end, &
                                          1:tdims%k_end) &
                           * (1.0 - EXP(-Cf_Inc /        &
                             bulk_cloud_fraction         &
                          (tdims%i_start:tdims%i_end,    &
                           tdims%j_start:tdims%j_end,    &
                                       1:tdims%k_end)))
        ! The following calculation should remain >0 so no error
        !trap needed
        scalingCoeff =                                    &
         (scalingCoeff +                                  &
          bulk_cloud_fraction(tdims%i_start:tdims%i_end,  &
                              tdims%j_start:tdims%j_end,  &
                                          1:tdims%k_end)) &
          / (scalingCoeff + Cf_Inc  )
      ELSEWHERE
        initialScaling = 0.0
        scalingCoeff = 0.0
      END WHERE
      Cb_Inc = initialScaling - (initialScaling - Cb_Inc) * &
                                scalingCoeff
    ELSEWHERE(Cf_Inc <                          &
              bulk_cloud_fraction               &
              (tdims%i_start:tdims%i_end,       &
               tdims%j_start:tdims%j_end,       &
                           1:tdims%k_end) .AND. &
              Cb_Inc - Cf_Inc > 0.0)
      Cb_Inc = Cb_Inc - Cf_Inc
      WHERE ((1.0 - bulk_cloud_fraction         &
                    (tdims%i_start:tdims%i_end, &
                     tdims%j_start:tdims%j_end, &
                                 1:tdims%k_end)) > 0.0)
        initialScaling = (1.0 - bulk_cloud_fraction             &
                          (tdims%i_start:tdims%i_end,           &
                           tdims%j_start:tdims%j_end,           &
                                       1:tdims%k_end))          &
                       * (1.0 - EXP(-Cb_Inc / (1.0 -            &
                                            bulk_cloud_fraction &
                         (tdims%i_start:tdims%i_end,            &
                          tdims%j_start:tdims%j_end,            &
                                      1:tdims%k_end))))
        scalingCoeff = (1.0 - bulk_cloud_fraction             &
                          (tdims%i_start:tdims%i_end,         &
                           tdims%j_start:tdims%j_end,         &
                                       1:tdims%k_end))        &
                        * (1.0 - EXP(-(1.0 - Cf_Inc) / (1.0 - &
                          bulk_cloud_fraction                 &
                         (tdims%i_start:tdims%i_end,          &
                          tdims%j_start:tdims%j_end,          &
                                      1:tdims%k_end))))
        scalingCoeff =                             &
            ((1.0 - bulk_cloud_fraction            &
                    (tdims%i_start:tdims%i_end,    &
                     tdims%j_start:tdims%j_end,    &
                                 1:tdims%k_end)) - &
             scalingCoeff) /                       &
            ((1.0 - Cf_Inc  ) - scalingCoeff)
      ELSEWHERE
        initialScaling = 0.0
        scalingCoeff = 0.0
      END WHERE
      Cb_Inc = initialScaling + (Cb_Inc - initialScaling) * &
                                scalingCoeff
    ELSEWHERE
      Cb_Inc = Cb_Inc - Cf_Inc
    END WHERE
    bulk_cloud_fraction(tdims%i_start:tdims%i_end,     &
                        tdims%j_start:tdims%j_end,     &
                                    1:tdims%k_end) =   &
        bulk_cloud_fraction(tdims%i_start:tdims%i_end, &
                            tdims%j_start:tdims%j_end, &
                                        1:tdims%k_end) &
        + Cb_Inc

  ELSE

    bulk_cloud_fraction(tdims%i_start:tdims%i_end,      &
                        tdims%j_start:tdims%j_end,      &
                                    1:tdims%k_end) =    &
       MIN(1.0,MAX(0.0,                                 &
        cloud_fraction_liquid                           &
                            (tdims%i_start:tdims%i_end, &
                             tdims%j_start:tdims%j_end, &
                                         1:tdims%k_end) &
        +                                               &
        cloud_fraction_frozen                           &
                            (tdims%i_start:tdims%i_end, &
                             tdims%j_start:tdims%j_end, &
                                         1:tdims%k_end)))

  END IF

  IF (i_cld_area == acf_off) THEN

    area_cloud_fraction  (tdims%i_start:tdims%i_end,   &
                          tdims%j_start:tdims%j_end,   &
                                      1:tdims%k_end) = &
      bulk_cloud_fraction(tdims%i_start:tdims%i_end,   &
                          tdims%j_start:tdims%j_end,   &
                                      1:tdims%k_end)

  ELSE

    area_cloud_fraction(tdims%i_start:tdims%i_end,     &
                        tdims%j_start:tdims%j_end,     &
                                    1:tdims%k_end) =   &
        MIN(1.0,                                       &
        AreaDivBulk *                                  &
        bulk_cloud_fraction(tdims%i_start:tdims%i_end, &
                            tdims%j_start:tdims%j_end, &
                                        1:tdims%k_end))

  END IF

  ! q
  IF (L_IAU_IncrementIce) THEN

    IF (L_IAU_ScaleCloud) THEN

      q(tdims%i_start:tdims%i_end,                    &
        tdims%j_start:tdims%j_end,                    &
                    1:tdims%k_end) =                  &
          q(tdims%i_start:tdims%i_end,                &
            tdims%j_start:tdims%j_end,                &
                        1:tdims%k_end)     +          &
          (qT_plus - qT)                            - &
          qcl_Inc -                                   &
          RESHAPE(qCF_plus - qCF_0,                   &
                  (/row_length, rows, model_levels/))

    ELSE

      q(tdims%i_start:tdims%i_end,                      &
        tdims%j_start:tdims%j_end,                      &
                    1:tdims%k_end) =                    &
          q(tdims%i_start:tdims%i_end,                  &
            tdims%j_start:tdims%j_end,                  &
                        1:tdims%k_end)     +            &
          (qT_plus - qT)                              - &
          RESHAPE(qCL_plus - qCL_0,                     &
                  (/row_length, rows, model_levels/)) - &
          RESHAPE(qCF_plus - qCF_0,                     &
                  (/row_length, rows, model_levels/))

    END IF

  ELSE

    IF (L_IAU_ScaleCloud) THEN

      q(tdims%i_start:tdims%i_end,                    &
        tdims%j_start:tdims%j_end,                    &
                    1:tdims%k_end) =                  &
          q(tdims%i_start:tdims%i_end,                &
            tdims%j_start:tdims%j_end,                &
                        1:tdims%k_end)     +          &
          (qT_plus - qT)                            - &
          qcl_Inc

    ELSE

      q(tdims%i_start:tdims%i_end,                    &
        tdims%j_start:tdims%j_end,                    &
                    1:tdims%k_end) =                  &
          q(tdims%i_start:tdims%i_end,                &
            tdims%j_start:tdims%j_end,                &
                        1:tdims%k_end)     +          &
          (qT_plus - qT)                            - &
          RESHAPE(qCL_plus - qCL_0,                   &
                  (/row_length, rows, model_levels/))

    END IF

  END IF

  ! Apply lower limit
  WHERE (q(tdims%i_start:tdims%i_end, &
           tdims%j_start:tdims%j_end, &
                       1:tdims%k_end) < q_min)

    q(tdims%i_start:tdims%i_end, &
      tdims%j_start:tdims%j_end, &
                  1:tdims%k_end) = q_min

  END WHERE

  !-----------------------------------------------------------------------------
  ! [7.6.3]: Deallocate work arrays.
  !-----------------------------------------------------------------------------

  DEALLOCATE (AreaDivBulk)
  DEALLOCATE (cloud_fraction_liquid_0)
  DEALLOCATE (cloud_fraction_liquid_plus)
  DEALLOCATE (p_theta_levels_old)
  DEALLOCATE (qCL_0)
  DEALLOCATE (qCL_plus)
  DEALLOCATE (Var_qT)
  DEALLOCATE (VarTemp)

  IF (L_IAU_IncrementIce) THEN
    DEALLOCATE (cloud_fraction_frozen_0)
    DEALLOCATE (cloud_fraction_frozen_plus)
    DEALLOCATE (qCF_0)
    DEALLOCATE (qCF_max)
    DEALLOCATE (qCF_plus)
    DEALLOCATE (Var_BGT)
  END IF

  IF (L_IAU_ScaleCloud) THEN
    DEALLOCATE (initialScaling)
    DEALLOCATE (scalingCoeff)
    DEALLOCATE (Cl_Inc)
    DEALLOCATE (Cf_Inc)
    DEALLOCATE (Cb_Inc)
  END IF

  DEALLOCATE (qcl_Inc)

END IF ! (qTIncsThisTS .AND. .NOT.L_IAU_Add_qT_prime_to_q)

IF (qTIncsThisTS) THEN

  ! Update q using qT_prime when Var_DiagCloud is bypassed:
  IF (L_IAU_Add_qT_prime_to_q) THEN
    q(tdims%i_start:tdims%i_end,       &
      tdims%j_start:tdims%j_end,       &
      1:tdims%k_end) =                 &
    qT_plus(tdims%i_start:tdims%i_end, &
            tdims%j_start:tdims%j_end, &
            1:tdims%k_end)  -          &
    qcl(tdims%i_start:tdims%i_end,     &
        tdims%j_start:tdims%j_end,     &
        1:tdims%k_end)

    IF (L_IAU_IncrementIce .AND. L_IAU_ApplyQTCorrections) THEN
      q(tdims%i_start:tdims%i_end, &
        tdims%j_start:tdims%j_end, &
        1:tdims%k_end) =           &
      q(tdims%i_start:tdims%i_end, &
        tdims%j_start:tdims%j_end, &
        1:tdims%k_end) -           &
      qcf(tdims%i_start:tdims%i_end, &
          tdims%j_start:tdims%j_end, &
          1:tdims%k_end)
    END IF

  END IF

  DEALLOCATE (qT)
  DEALLOCATE (qT_plus)

END IF

!-------------------------------------------------------------------------------
! [7.7]: PC2 updates to water vapour, cloud and temperature.
!-------------------------------------------------------------------------------

IF (.NOT. qTIncsThisTS .AND. &
   i_cld_vn==i_cld_pc2 .AND. (L_pc2_cond .OR. L_pc2_cfl .OR. L_pc2_cff)) THEN

  ALLOCATE (delta_q  (tdims%i_start:tdims%i_end, &
                      tdims%j_start:tdims%j_end, &
                                  1:tdims%k_end))
  ALLOCATE (delta_qcl(tdims%i_start:tdims%i_end, &
                      tdims%j_start:tdims%j_end, &
                                  1:tdims%k_end))
  ALLOCATE (delta_qcf(tdims%i_start:tdims%i_end, &
                      tdims%j_start:tdims%j_end, &
                                  1:tdims%k_end))
  ALLOCATE (delta_p  (tdims%i_start:tdims%i_end, &
                      tdims%j_start:tdims%j_end, &
                                  1:tdims%k_end))
  ALLOCATE (delta_t  (tdims%i_start:tdims%i_end, &
                      tdims%j_start:tdims%j_end, &
                                  1:tdims%k_end))

  ! Calculate increments to vapour, cloud, pressure and temperature.
  ! (The _work arrays currently contain values on entry to the routine.)
  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i=  tdims%i_start, tdims%i_end
        delta_q  (i,j,k) = q  (i,j,k)                - q_work  (i,j,k)
        delta_qcl(i,j,k) = qcl(i,j,k)                - qcl_work(i,j,k)
        delta_qcf(i,j,k) = qcf(i,j,k)                - qcf_in  (i,j,k)
        delta_p  (i,j,k) = p_theta_levels    (i,j,k) - p_work  (i,j,k)
        delta_t  (i,j,k) = exner_theta_levels(i,j,k) &
                         * theta(i,j,k)              - t_work  (i,j,k)
      END DO
    END DO
  END DO

  ! Note that the cloud_fraction arrays have halos, so their dimensions have
  ! to be specified on passing into pc2_assim, in order to exclude the halos.
  ! Also the bottom (zero) level is excluded. 

  CALL pc2_assim ( REAL(TS_len_secs),                            & ! in
                   l_pc2_cfl,                                    & ! in
                   l_pc2_cff,                                    & ! in
                   l_mr_physics,                                 & ! in
                   t_work,                                       & ! inout
      bulk_cloud_fraction  (tdims%i_start:tdims%i_end,           & 
                            tdims%j_start:tdims%j_end,           &  
                                        1:tdims%k_end),          & ! inout
      cloud_fraction_liquid(tdims%i_start:tdims%i_end,           & 
                            tdims%j_start:tdims%j_end,           &  
                                        1:tdims%k_end),          & ! inout
      cloud_fraction_frozen(tdims%i_start:tdims%i_end,           & 
                            tdims%j_start:tdims%j_end,           &  
                                        1:tdims%k_end),          & ! inout
                   q_work,                                       & ! inout
                   qcl_work,                                     & ! inout
                   qcf_in,                                       & ! in
                   p_work,                                       & ! in
                   delta_t,                                      & ! in
                   delta_q,                                      & ! in
                   delta_qcl,                                    & ! in
                   delta_qcf,                                    & ! in
                   delta_p                                       & ! in
                  )

  ! q_work, qcl_work and t_work have now been updated by the forcing and
  ! condensation terms together. cloud_fraction_liquid (and _frozen and _bulk)
  ! has been updated by the condensation. qcf_work and p_work are not updated.

  ! Now copy q_work, qcl_work and t_work variables back into the INOUT
  ! variables if we haven't a confident direct estimate of condensation from
  ! another source:
  IF (L_pc2_cond) THEN
    DO k =             1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i= tdims%i_start, tdims%i_end
          q    (i,j,k) = q_work  (i,j,k)
          qcl  (i,j,k) = qcl_work(i,j,k)
          ! Remember the INOUT variable is theta, not temperature:
          theta(i,j,k) = t_work  (i,j,k) / exner_theta_levels(i,j,k)
        END DO
      END DO
    END DO
  END IF

  DEALLOCATE (q_work)
  DEALLOCATE (qcl_work)
  DEALLOCATE (qcf_in)
  DEALLOCATE (p_work)
  DEALLOCATE (t_work)
  DEALLOCATE (delta_q)
  DEALLOCATE (delta_qcl)
  DEALLOCATE (delta_qcf)
  DEALLOCATE (delta_p)
  DEALLOCATE (delta_t)

ELSE IF (i_cld_vn==i_cld_pc2 .AND. (l_pc2_check_init .OR. l_iau_pc2check)) THEN
  ! check that the cloud fraction and condensate values produced by the
  ! IAU are consistent

  CALL initial_pc2_check(exner_theta_levels, p_theta_levels, theta,        &
                         bulk_cloud_fraction, cloud_fraction_liquid,       &
                         cloud_fraction_frozen, q, qcl, qcf, l_mr_physics)

END IF ! PC2 section

!-------------------------------------------------------------------------------
! [7.8]: TSoil(1), TStar and TStar_tile.
!-------------------------------------------------------------------------------

IF (L_IAU_IncTStar) THEN

  IF (PrintStatus >= PrStatus_Normal) THEN
    WRITE(umMessage,'(A)')                                              &
          'IAU: Adding t1 incs to TStar and TSoil(1)'
    CALL umPrint(umMessage,src='iau')
  END IF

  ! Get increment to level-one temperature:
  addr = 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      t1_inc(addr) = exner_theta_levels(i,j,1) &
                   * theta             (i,j,1) &
                   - t1_in(i,j)
      addr = addr + 1
    END DO
  END DO

  ! Separate the case of the multilayer snow scheme.   
  !   
  IF ( nsmax > 0 ) THEN   
   
    !   Multilayer snow scheme.   
   
    !   NOTES:   
    !   1.) The logical switches on the surface temperature present 
    !       in the original case of the zero-layer scheme,
    !       L_IAU_IncTStar_tile and L_IAU_IncTSurf_Snow, have   
    !       been deleted, since it is assumed that we will always want   
    !       to apply these increments in the case of multilayer snow
    !       without maintaining historical functionality.
    !   2.) The GBM TStar is not incremented, since it is assumed that   
    !       it will be recalculated from the tiled temperatures.
    !   3.) The tests that the skin and snow temperatures are below freezing   
    !       before incrementing may be redundant, but are retained for now.   

    !   Calculate the weighting factor for each snow layer.
    !   We expect the first snow layer to be so thin that 
    !   we shall always want to increment it. Typically,
    !   the phase of the diurnal cycle will be reversed
    !   at a depth of about 0.3 m. The weight applied to
    !   the increment is therefore decreased smoothly from
    !   half this depth (diffusive timescale 4 times smaller)
    !   to 0 at 0.3 m.
    !
    wgt_snow(1) = 1.0
    z_snow = 0.5 * dzsnow(1)
    DO js = 2, nsmax
      z_snow = z_snow + 0.5 * ( dzsnow(js-1) + dzsnow(js) )
      wgt_snow(js) = (0.3 - z_snow) / (0.3 - 0.15)
      wgt_snow(js) = MAX( 0.0, MIN( 1.0, wgt_snow(js) ) )
    END DO
   
    !   Increment the tiled surface temperatures and the soil temperature.
    DO tile_num = 1, ntiles
      IF ( tile_num /= lake ) THEN
        DO i = 1, land_field

          !         Set the tile fraction to the fraction of the
          !         surface type unless aggregating: tiles and types
          !         will match if not.
          tile_frac = frac_typ(i,tile_num)
          IF (l_aggregate) tile_frac = 1.0


          !         Increment the soil temperature if there is no snow and   
          !         the snow temperature if there is.   

          IF ( NINT(nsnow(i,tile_num)) == 0 ) THEN 

            IF ( tile_frac > 0.0 ) &
              TStar_tile(i,tile_num) = TStar_tile(i,tile_num) &
                                     + t1_inc(land_index(i))
            IF ( tile_num /= lake ) & 
              TSoil(i,1) = TSoil(i,1) &   
                + tile_frac * t1_inc(land_index(i))   

          ELSE   

            !       Apply the full increment to the skin temperature
            IF (tile_frac > 0.0) THEN
              IF (TStar_tile(i,tile_num) < tm) & 
                TStar_tile(i,tile_num) = MIN(TStar_tile(i,tile_num)  & 
                     + t1_inc(land_index(i)), tm) 

!             Apply a decreasing increment through the snow pack.
              DO js = 1, nsmax
                IF (js <= nsnow(i,tile_num)) THEN
                  IF (tsnow(i,tile_num,js) < tm)                       &
                    tsnow(i,tile_num,js) = MIN( tsnow(i,tile_num,js) + &
                         wgt_snow(js) * t1_inc(land_index(i)), tm)
                END IF
              END DO

!             Apply an increment to the soil temperature if the snow
!             is very thin. The threshold is partly historical and
!             a proper treatment would take full account of the physical
!             timescales in the DA interval.
              IF (TSoil(i,1) < tm) THEN
                wgt_soil_snow = MAX(0.0, 1.0-snowdepth_z(i,tile_num)/0.04)
                TSoil(i,1) = MIN(tm, TSoil(i,1) &   
                  + wgt_soil_snow * tile_frac * t1_inc(land_index(i)) )
              END IF
            END IF

          END IF

        END DO
      END IF
    END DO

  ELSE   
    !   
    !   Zero-layer snow scheme   
    !   
    !   Add increments to TSoil(1) and TStar where   
    !   snow_depth <= 0.05 kg/m**2 or explicitly requested:   
    DO i = 1, land_field   
   
      IF ( snow_depth(land_index(i)) <= 0.05 ) THEN   
        TSoil(i,1) = TSoil(i,1) + t1_inc(land_index(i))   
        TStar(land_index(i)) =  TStar(land_index(i)) &   
                           + t1_inc(land_index(i))   
      ELSE IF (L_IAU_IncTSurf_Snow) THEN   
        !       Limit the increment not to allow melting.   
        IF (TSoil(i,1) < tm) TSoil(i,1) = MIN(TSoil(i,1)       &   
                             + t1_inc(land_index(i)), tm)   
        IF (TStar(land_index(i)) < tm) TStar(land_index(i)) =  &   
                             MIN(TStar(land_index(i)) &   
                             + t1_inc(land_index(i)), tm)   
      END IF   
   
    END DO   
   
    ! If using a tile scheme, also add increments to TStar_tile:   
    IF (L_IAU_IncTStar_tile) THEN   
   
      DO tile_num = 1, ntiles   
        IF ( tile_num /= lake ) THEN 
          DO i = 1, land_field   
     
            IF ( snow_depth(land_index(i)) <= 0.05 ) THEN   
              TStar_tile(i,tile_num) = TStar_tile(i,tile_num) &   
                                     + t1_inc(land_index(i))   
            ELSE IF (L_IAU_IncTSurf_Snow) THEN   
              IF (TStar_tile(i,tile_num) < tm) TStar_tile(i,tile_num) &   
                       = MIN(TStar_tile(i,tile_num) &   
                       + t1_inc(land_index(i)), tm)   
            END IF   
   
          END DO   
        END IF   
      END DO   
   
    END IF   
    
  END IF

END IF ! (L_IAU_IncTStar)

!-------------------------------------------------------------------------------
! [8]: Remove supersaturation wrt water.
!-------------------------------------------------------------------------------

IF (L_IAU_RemoveSS) THEN

  IF (PrintStatus >= PrStatus_Normal) THEN
    WRITE(umMessage,'(A)') 'IAU: Removing supersaturation wrt water'
    CALL umPrint(umMessage,src='iau')
  END IF

  DO k = 1, tdims%k_end

    ! Get temperature and pressure:
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ! Use t1_in work array to hold temp on this level:
        t1_in(i,j) = exner_theta_levels(i,j,k) * theta(i,j,k)
        p_tmp(i,j) = p_theta_levels    (i,j,k)
      END DO
    END DO

    IF ( l_new_qsat_acassim ) THEN
      CALL qsat_wat_new(q_sat_wat, t1_in, p_tmp, row_length, rows)
    ELSE
      ! Calculate saturated specific humidity wrt water:
      CALL qsat_wat (q_sat_wat, t1_in, p_tmp, row_length*rows)
    END IF

    ! Limit q to q_sat:
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q(i,j,k) = MIN(q(i,j,k), q_sat_wat(i,j))
      END DO
    END DO

  END DO

END IF ! L_IAU_RemoveSS

!-------------------------------------------------------------------------------
! [9]: Apply bounds to tropospheric and non-tropospheric humidities.
!-------------------------------------------------------------------------------

IF (CallQLimitsThisTS) THEN

  IF (PrintStatus >= PrStatus_Normal) THEN
    WRITE(umMessage,'(A)') 'IAU: Calling QLimits'
    CALL umPrint(umMessage,src='iau')
  END IF

  ! Need to update halos ready for calculation of pv_at_theta:
  i_field = 1
  fields_to_swap(i_field) % field      => u(:,:,:)
  fields_to_swap(i_field) % field_type =  fld_type_u
  fields_to_swap(i_field) % levels     =  udims%k_len
  fields_to_swap(i_field) % rows       =  rows
  fields_to_swap(i_field) % vector     = .TRUE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field      => v(:,:,:)
  fields_to_swap(i_field) % field_type =  fld_type_v
  fields_to_swap(i_field) % levels     =  vdims%k_len
  fields_to_swap(i_field) % rows       =  n_rows
  fields_to_swap(i_field) % vector     = .TRUE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field      => theta(:,:,:)
  fields_to_swap(i_field) % field_type =  fld_type_p
  fields_to_swap(i_field) % levels     =  tdims%k_len
  fields_to_swap(i_field) % rows       =  rows
  fields_to_swap(i_field) % vector     = .FALSE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field      => rho(:,:,:)
  fields_to_swap(i_field) % field_type =  fld_type_p
  fields_to_swap(i_field) % levels     =  pdims%k_len
  fields_to_swap(i_field) % rows       =  rows
  fields_to_swap(i_field) % vector     = .FALSE.

  CALL swap_bounds_mv (fields_to_swap, i_field, row_length, offx, offy)

  ALLOCATE (pv_at_theta(tdims%i_start:tdims%i_end, &
                        tdims%j_start:tdims%j_end, &
                                    1:tdims%k_end))

  CALL Calc_PV_at_theta ( u, v, theta, rho,             & ! in
                          r_theta_levels, r_rho_levels, & ! in
                          r_at_u, r_at_v,               & ! in
                          sec_v_latitude,               & ! in
                          tan_v_latitude,               & ! in
                          sec_theta_latitude,           & ! in
                          f3_at_v,                      & ! in
                          delta_lambda, delta_phi,      & ! in
                          pv_at_theta )                   ! out

  CALL QLimits ( L_IAU_QLimitsDiags,    & ! in
                 L_IAU_RmNonTropQIncs,  & ! in
                 IAU_trop_min_RH,       & ! in
                 IAU_nonTrop_max_q,     & ! in
                 IAU_nonTrop_min_q,     & ! in
                 IAU_nonTrop_max_RH,    & ! in
                 IAU_trop_min_p,        & ! in
                 IAU_trop_max_PV,       & ! in
                 IAU_nonTrop_max_p,     & ! in
                 q_in,                  & ! in
                 pv_at_theta,           & ! in
                 theta,                 & ! in
                 p,                     & ! in
                 p_theta_levels,        & ! in
                 exner_theta_levels,    & ! in
                 q,                     & ! inout
                 qCL,                   & ! inout
                 qCF,                   & ! inout
                 area_cloud_fraction,   & ! inout
                 bulk_cloud_fraction,   & ! inout
                 cloud_fraction_liquid, & ! inout
                 cloud_fraction_frozen )  ! inout

  DEALLOCATE (pv_at_theta)
  DEALLOCATE (q_in)

END IF ! (CallQLimitsThisTS)

!-------------------------------------------------------------------------------
! [10]: IAU may have changed wet density. The dry density has to see this
!       for conversions at the end to be accurate.
!-------------------------------------------------------------------------------

CALL eg_q_to_mix                                                &
          (tdims_l,tdims_s,                                     &
           q, qcl, qcf,                                         &
           qcf2, qrain, qgraup,                                 &
           l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
           m_v, m_cl, m_cf                                      &
          ,m_cf2, m_r, m_gr, swap_in=.FALSE.)

CALL update_dryrho_mix_tl_n(rho)

DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      thetav(i,j,k) = theta3d(i,j,k)*(1.0+m_v(i,j,k)*recip_epsilon)
    END DO
  END DO
END DO

CALL swap_bounds(thetav,                                               &
                 tdims_s%i_len - 2*tdims_s%halo_i,                     &
                 tdims_s%j_len - 2*tdims_s%halo_j,                     &
                 tdims_s%k_len,                                        &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p,swap_field_is_scalar)

!-------------------------------------------------------------------------------
! [11]: Update halos as necessary.
!
!       In general we don't update the halo values of fields that have been
!       modified by this routine, on the assumption that they will be updated
!       by other parts of the code before they are used. This is risky, but
!       reduces the runtime.
!
!       However, there are a few fields for which halo updates have been shown
!       to be necessary to make results insensitive to processor configuration,
!       so this section performs the necessary swap_bounds calls.
!-------------------------------------------------------------------------------

i_field = 0

IF (.NOT.CallQLimitsThisTS) THEN
  i_field = i_field + 1
  fields_to_swap(i_field) % field      => u(:,:,:)
  fields_to_swap(i_field) % field_type =  fld_type_u
  fields_to_swap(i_field) % levels     =  udims%k_len
  fields_to_swap(i_field) % rows       =  rows
  fields_to_swap(i_field) % vector     = .TRUE.
END IF

IF (.NOT.CallQLimitsThisTS) THEN
  i_field = i_field + 1
  fields_to_swap(i_field) % field      => v(:,:,:)
  fields_to_swap(i_field) % field_type =  fld_type_v
  fields_to_swap(i_field) % levels     =  vdims%k_len
  fields_to_swap(i_field) % rows       =  n_rows
  fields_to_swap(i_field) % vector     = .TRUE.
END IF

i_field = i_field + 1
fields_to_swap(i_field) % field      => exner(:,:,:)
fields_to_swap(i_field) % field_type =  fld_type_p
fields_to_swap(i_field) % levels     =  pdims%k_len + 1
fields_to_swap(i_field) % rows       =  rows
fields_to_swap(i_field) % vector     = .FALSE.

CALL swap_bounds_mv (fields_to_swap(1:i_field), \
                     i_field, row_length, offx, offy)

!-------------------------------------------------------------------------------
! [12]: Loop through IAU structures and deallocate arrays that are no longer
!       needed.
!-------------------------------------------------------------------------------

DO IncNum = 1, Num_IAU_incs

  IF (stepim(atmos_im) == IAU_incs(IncNum) % InsertionEndTS .AND. &
                      IAU_incs(IncNum) % FieldsToAdd) THEN

    DEALLOCATE (IAU_incs(IncNum) % Lookup)
    DEALLOCATE (IAU_incs(IncNum) % Weights)

    IF (IAU_incs(IncNum) % FieldsStored) THEN
      IF (IAU_incs(IncNum) % StorePacked) THEN
        DEALLOCATE (IAU_incs(IncNum) % flds32bit)
      ELSE
        DEALLOCATE (IAU_incs(IncNum) % flds)
      END IF
      IAU_incs(IncNum) % FieldsStored = .FALSE.
    END IF

  END IF

END DO

!-------------------------------------------------------------------------------
! [13]: Exit the routine.
!-------------------------------------------------------------------------------

IF (PrintStatus >= PrStatus_Normal) THEN
  WRITE(umMessage,'(A)')                                                &
        'IAU: Exiting IAU. Time since basis time (hr:mm:ss): '//TimeStr
  CALL umPrint(umMessage,src='iau')
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE iau


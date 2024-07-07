! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection

MODULE diagnostics_conv_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DIAGNOSTICS_CONV_MOD'
CONTAINS

SUBROUTINE diagnostics_conv(                                      &
                       row_length, rows, at_extremity             &
,                      u, v, p, r_u, r_v                          &
,                      dubydt_conv_p, dvbydt_conv_p               &
,                      dubydt_conv_u, dvbydt_conv_v               &
,                      theta_star, q_star                         &
,                      exner_theta_levels                         &
,                      ls_rain, ls_snow                           &
,                      ccw, conv_rain, conv_snow                  &
,                      cca_2d, cca, ccb, cct                      &
,                      cu_over_orog,cape_undilute,cin_undilute    &
,                      lcbase, lctop, lcca, n_cca_levels          &
,                      l_dust                                     &
,                      conv_type                                  &
,                      timestep,                                  &
     stashwork                                                    &
     )

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                          &
   udims, vdims, udims_s, vdims_s, tdims, tdims_s, pdims_s

! Convective diagnostic output arrays
USE cv_diagnostic_array_mod, ONLY:                                      &
        precip_deep, precip_shall ,precip_mid ,precip_cong ,cape_out    &
       ,deep_ind ,shallow_ind ,congestus_ind ,congestus_ind2 ,mid_ind   &
       ,ntml_diag ,ntpar_diag ,freeze_diag ,kterm_diag                  &
       ,wstar_up_diag ,wstar_dn_diag ,mb1_diag ,mb2_diag ,wqt_cb        &
       ,wthetal_cb ,wqt_inv ,wthetal_inv ,sh_top ,sh_base ,cg_top       &
       ,cg_base ,cg_term ,cape_ts_diag ,ind_cape_reduced_diag           &
       ,conscav_dust ,conscav_so4ait ,conscav_so4acc ,conscav_so4dis    &
       ,conscav_agedsoot ,conscav_agedbmass ,conscav_agedocff           &
       ,conscav_nitracc ,conscav_nitrdiss ,conwash_so2 ,conwash_nh3     &
       ,uw_dp, vw_dp ,uw_shall ,vw_shall ,uw_mid ,vw_mid ,up_flux       &
       ,dwn_flux ,entrain_up ,detrain_up ,entrain_dwn ,detrain_dwn      &
       ,up_flux_half ,T_incr_diag_conv ,q_incr_diag_conv                &
       ,qcl_incr_diag_conv ,qcf_incr_diag_conv                          &
       ,cf_liquid_incr_diag_conv ,cf_frozen_incr_diag_conv              &
       ,bulk_cf_incr_diag_conv,t_incr_conv_only,q_incr_conv_only        &
       ,mf_deep ,mf_congest ,mf_shall ,mf_midlev                        &
       ,dt_deep ,dt_congest ,dt_shall ,dt_midlev                        &
       ,dq_deep ,dq_congest ,dq_shall ,dq_midlev                        &
       ,du_deep ,du_congest ,du_shall ,du_midlev                        &
       ,dv_deep ,dv_congest ,dv_shall ,dv_midlev                        &
       ,wqt_flux_sh ,wql_flux_sh ,wthetal_flux_sh ,wthetav_flux_sh      &
       ,conv_rain_3d ,conv_snow_3d                                      &
       ,qcl_incr_inhom_diag ,qcf_incr_inhom_diag                        &
       ,bulk_cf_incr_inhom_diag ,cf_liquid_incr_inhom_diag              &
       ,cf_frozen_incr_inhom_diag, deep_cfl_limited, mid_cfl_limited    &
       ,deep_tops, w2p, dt_dd, dq_dd, du_dd, dv_dd, area_ud, area_dd

! ----------------------------------------------------------------------+-------
USE science_fixes_mod,  ONLY: l_fix_conv_diags_var
USE cv_param_mod,  ONLY: max_mf_fall
USE ac_diagnostics_mod, ONLY:                                            &
   cvrr, cvsr, convcc, tinc_cvn
USE acp_namel_mod, ONLY: l_ac
USE pws_diags_mod, ONLY: pws_precip_sym_conv_rain, pws_precip_sym_conv_snow, &
      pws_conv_icao_base, pws_conv_icao_top
USE um_stashcode_mod, ONLY: stashcode_pws_sec, stashcode_pws_precip_sym, &
      stashcode_pws_conv_cld_dep, stashcode_pws_conv_icao_base,          &
      stashcode_pws_conv_icao_top
USE icao_ht_fc_mod, ONLY: ICAOHeight_FC

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE cv_run_mod,  ONLY:                                             &
   tower_factor, l_dcpl_cld4pc2, l_3d_cca

USE dust_parameters_mod, ONLY: l_twobin_dust, ndiv,                &
                               tot_dust_dep_flux
USE nlsizes_namelist_mod,   ONLY: model_levels

USE ereport_mod, ONLY: ereport
USE icao_ht_mod, ONLY: icao_ht
USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE missing_data_mod, ONLY: rmdi
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Calculates convection-related diagnostics (held in STASH section 5).
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Intercepted arrays and diagnostic arrays are input from the
!   convection routines called previously. Each diagnostic is simply
!   copied into the STASHwork array to be passed on to STASH for
!   output processing.

!   Diagnostics currently available (in numerical order):
!   Item  Description
!     10  q at end of timestep
!    140  QCL incr from conv_inhom, positive part
!    141  QCL incr from conv_inhom, negative part
!    142  QCF incr from conv_inhom, positive part
!    143  QCF incr from conv_inhom, negative part
!    146  CFL incr from conv_inhom, positive part
!    147  CFL incr from conv_inhom, negative part
!    148  CFF incr from conv_inhom, positive part
!    149  CFF incr from conv_inhom, negative part
!    150  QCL incr from hom_conv
!    151  QCL incr from hom_conv, positive part
!    152  QCL incr from hom_conv, negative part
!    156  CFL incr from hom_conv
!    157  CFL incr from hom_conv, positive part
!    158  CFL incr from hom_conv, negative part
!    161  Temperature increment (convection only)
!    162  Humidity increment (convection only)
!    163  Liquid Cloud Condensate increment (inhomog)
!    164  Frozen Cloud Condensate increment (inhomog)
!    172  Total  Cloud Volume increment (inhomog)
!    173  Liquid Cloud Volume increment (inhomog)
!    174  Frozen Cloud Volume increment (inhomog)
!    175  du/dt from downdraughts (m/s/s) on p-grid
!    176  dv/dt from downdraughts (m/s/s) on p-grid
!    181  Temperature increment (may include PC2 homo with PC2 on)
!    182  Humidity increment (may include PC2 homo with PC2 on)
!    183  Liquid Cloud Condensate increment
!    184  Frozen Cloud Condensate increment
!    185  u wind increment
!    186  v wind increment
!    187  Temperature increment in no shallow convection columns
!    188  Humidity increment in no shallow convection columns
!
!    192  Total  Cloud Volume increment
!    193  Liquid Cloud Volume increment
!    194  Frozen Cloud Volume increment
!    196  (w)^2 of parcel
!    197   w of parcel where w^2>0.0
!    198  dT/dt from downdraughts  (K/s)
!    199  dq/dt from downdraughts   (kg/kg/s)
!
!    201  Convective rainfall
!    202  Convective snowfall
!    205  Convective rainfall rates
!    206  Convective snowfall rates
!    207  Pressure at convective cloud base
!    208  Pressure at convective cloud top
!    209  Temperature at end of timestep
!    210  Convective cloud base (ICAO standard atmosphere heights)
!    211  Convective cloud top (ICAO standard atmosphere heights)
!    212  Convective cloud amount
!    213  Convective cloud liquid water
!    214  Total rain rate
!    215  Total snow rate
!    216  Total precipitation rate
!    217  CAPE
!    218  Lowest convective cloud base
!    219  Lowest convective cloud top
!    220  Lowest convective cloud amount
!    221
!    222  Pressure at lowest convective cloud base
!    223  Pressure at lowest convective cloud top
!    224  Lowest convective cloud base(ICAO standard atmosphere heights)
!    225  Lowest convective cloud top (ICAO standard atmosphere heights)
!    226  Total precipitation
!    227  3D Convective rainfall rate
!    228  3D Convective snowfall rate
!    229  fractional area of updraughts
!    230  fractional area of downdraughts
!    231  Actual CAPE timescale used for deep convection
!    232  indicator of reduced CAPE timescale
!    233  undilute parcel CAPE
!    234  undilute parcel CIN
!    235  u latest
!    236  v latest
!    237  NH3 SCAVENGED BY CONV PPN KG/M2/SEC
!    238  SO2 SCAVENGED BY CONV PPN KG/M2/SEC
!    239  SO4 AIT SCAVNGD BY CONV PPN KG/M2/S
!    240  SO4 ACC SCAVNGD BY CONV PPN KG/M2/S
!    241  SO4 DIS SCAVNGD BY CONV PPN KG/M2/S
!    242  Aged soot scavenged by convective precipitation kg/m2/s
!    243  Aged biomass scavenged by convective precipitation kg/m2/s
!    244  Aged fossil-fuel OC scavenged by convective precipitation kg/m2/s
!    246  Component of UPDRAUGHT MASS FLUX ON HALF (RHO) LEVELS (Pa/s)
!         needed for 4D-VAR
!    247  Accumulation-mode ammonium nitrate scavenged by convective
!         precip kg[N]/m2/s
!    248  Dissolved ammonium nitrate scavenged by convective precip kg[N]/m2/s
!    249  UPDRAUGHT MASS FLUX ON HALF (RHO) LEVELS (Pa/s)
!    250  UPDRAUGHT MASS FLUX (Pa/s)
!    251  DOWNDRAUGHT MASS FLUX
!    252  UPDRAUGHT ENTRAINMENT RATE
!    253  UPDRAUGHT DETRAINMENT RATE
!    254  DOWNDRAUGHT ENTRAINMENT RATE
!    255  DOWNDRAUGHT DETRAINMENT RATE
!    256  U INCREMENT (P GRID)
!    257  V INCREMENT (P GRID)
!    258  DEEP UW STRESS (P GRID)
!    259  DEEP VW STRESS (P GRID)
!    260  SHALLOW UW STRESS (P GRID)
!    261  SHALLOW UW STRESS (P GRID)
!    262  2D Convective cloud amount
!    263  Mid UW stress (P grid)
!    264  Mid VW stress (P grid)
!    267  CFL limited deep convection indicator
!    268  CFL limited mid convection indicator
!    269  Deep cumulus indicator
!    270  SHALLOW CUMULUS INDICATOR
!    271  CU OVER OROGRAPHY INDICATOR
!    272  Mid level onvection indicator
!    273  NTML
!    274  NTPAR
!    275  Freezing level number
!    276  KTERM deep
!    277  total precipitation from deep convection
!    278  total precipitation from shallow convection
!    279  total precipitation from mid-level convection
!    280  total precipitation from congestus convection
!    281  Dust wet dep flux conv precip div 1
!    282  Dust wet dep flux conv precip div 2
!    283  Dust wet dep flux conv precip div 3
!    284  Dust wet dep flux conv precip div 4
!    285  Dust wet dep flux conv precip div 5
!    286  Dust wet dep flux conv precip div 6
! NEW DIAGNOSTICS FOR 5A turbulence part of convection scheme only
!    290  w'qt' flux
!    291  w'ql' flux
!    292  w'thetal' flux
!    293  w'thetav' flux
!    300  wstar_dn
!    301  wstar_up
!    302  cloud base mass flux  from wstar_dn
!    303  alternative cloud base mass flux
!    304  cloud base w'qt'
!    305  cloud base w'thetal'
!    306  inversion w'qt'
!    307  inversion w'thetal'
!    308  base of shallow convection (m)
!    309  top of shallow convection (m)

!    310  Congestus indicator
!    311  Congestus indicator 2
!    312  Termination level for congestus
!    313  height of top of congestus
!    314  height of base of congestus

!    319  Frequency deep terminates on a given model level (4A as well as 5A)
!    320  mass flux deep convection
!    321  mass flux congestus convection
!    322  mass flux shallow convection
!    323  mass flux mid_lev convection
!    324  dT deep convection
!    325  dT congestus convection
!    326  dT shallow convection
!    327  dT mid_lev convection
!    328  dq deep convection
!    329  dq congestus convection
!    330  dq shallow convection
!    331  dq mid_lev convection
!    332  du deep convection
!    333  du congestus convection
!    334  du shallow convection
!    335  du mid_lev convection
!    336  dv deep convection
!    337  dv congestus convection
!    338  dv shallow convection
!    339  dv mid_lev convection

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards Version 8.3.
!------------------------------------------------------------------------------

LOGICAL, INTENT(IN) ::  &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid

! Arguments with Intent IN. ie: Input variables.
INTEGER, INTENT(IN) ::  &
  row_length            & ! number of points on a row
, rows                  & ! number of rows in a theta field
, n_cca_levels            ! number of convective cloud levels

LOGICAL, INTENT(IN) ::  &
  l_dust

REAL, INTENT(IN) :: &
  timestep            ! timestep


! Primary Arrays used in all models

REAL, INTENT(IN) ::                       &
    u(udims_s%i_start:udims_s%i_end,      & ! U (m/s)
      udims_s%j_start:udims_s%j_end,      &
      udims_s%k_start:udims_s%k_end)      &
  , v(vdims_s%i_start:vdims_s%i_end,      & ! V (m/s)
      vdims_s%j_start:vdims_s%j_end,      &
      vdims_s%k_start:vdims_s%k_end)      &
  , p(pdims_s%i_start:pdims_s%i_end,      & ! pressure (Pa)
      pdims_s%j_start:pdims_s%j_end,      &
      pdims_s%k_start:pdims_s%k_end)      &
   ,r_u(udims_s%i_start:udims_s%i_end,    & ! r_u (m/s)
        udims_s%j_start:udims_s%j_end,    &
        udims_s%k_start:udims_s%k_end)    &
  , r_v(vdims_s%i_start:vdims_s%i_end,    & ! r_v (m/s)
        vdims_s%j_start:vdims_s%j_end,    &
        vdims_s%k_start:vdims_s%k_end)

! Convective momentum transport tendencies on p-grid
REAL, INTENT(IN) :: dubydt_conv_p ( pdims_s%i_start:pdims_s%i_end,  &
                                    pdims_s%j_start:pdims_s%j_end,  &
                                    pdims_s%k_start:pdims_s%k_end )
REAL, INTENT(IN) :: dvbydt_conv_p ( pdims_s%i_start:pdims_s%i_end,  &
                                    pdims_s%j_start:pdims_s%j_end,  &
                                    pdims_s%k_start:pdims_s%k_end )
! Convective momentum transport tendencies on natives grids
REAL, INTENT(IN) :: dubydt_conv_u ( udims%i_start:udims%i_end,      &
                                    udims%j_start:udims%j_end,      &
                                    udims%k_start:udims%k_end )
REAL, INTENT(IN) :: dvbydt_conv_v ( vdims%i_start:vdims%i_end,      &
                                    vdims%j_start:vdims%j_end,      &
                                    vdims%k_start:vdims%k_end )


REAL, INTENT(IN) ::                                               &
  theta_star ( tdims%i_start: tdims%i_end,                        &
               tdims%j_start: tdims%j_end,                        &
               tdims%k_start: tdims%k_end )                       &
, q_star     ( tdims%i_start: tdims%i_end,                        &
               tdims%j_start: tdims%j_end,                        &
               tdims%k_start: tdims%k_end )                       &
, ls_rain(row_length, rows)                                       &
, ls_snow(row_length, rows)                                       &
, ccw(row_length, rows, model_levels)                             &
, conv_rain(row_length, rows)                                     &
, conv_snow(row_length, rows)                                     &
, cca(row_length, rows, n_cca_levels)                             &
, cca_2d(row_length, rows)                                        &
 ,exner_theta_levels(tdims_s%i_start:tdims_s%i_end,  & ! Exner on theta level
                     tdims_s%j_start:tdims_s%j_end,  & !
                     tdims_s%k_start:tdims_s%k_end)

! Type diagnostics
REAL, INTENT(IN)     ::          &
  cu_over_orog(row_length,rows)  & ! Indicator for cumulus over steep
                                   ! orography. Indicator set to 1.0 if true,
                                   !   0.0 if false. Exclusive.
, cape_undilute(row_length,rows) & ! CAPE from undilute parcel ascent
, cin_undilute(row_length,rows)    ! CIN from undilute parcel ascent

INTEGER, INTENT(IN) ::                                            &
  ccb(row_length, rows)                                           &
, cct(row_length, rows)                                           &
, lcbase(row_length, rows)                                        &
, lctop(row_length, rows)
REAL, INTENT(IN) ::                                               &
  lcca(row_length, rows)

INTEGER, INTENT(IN) :: conv_type(row_length, rows)

! Diagnostics info
REAL, INTENT(INOUT) ::                                            &
 stashwork(*)     ! STASH workspace

!-------------------------------------------------------------------------
! Local variables
INTEGER ::           &
  i, j, k, ji, itmp  & ! loop counters
 ,icode              & ! Return code  =0 Normal exit  >1 Error
 ,item               & ! STASH item
 ,im_index             ! internal model index

INTEGER, PARAMETER :: &
  sect = 5              ! STASH section for convection

INTEGER ::    &
  u_rows      &  ! number of u_rows
 ,v_rows         ! number of v_rows

LOGICAL :: l_tddf       ! Indicate if total dust dep flux is required

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*) :: RoutineName
PARAMETER ( RoutineName='DIAGNOSTICS_CONV')

! max_2d_cca is the maximum aount of 2d cca permitted
REAL, PARAMETER :: max_2d_cca = 0.5

REAL :: icao_height(row_length,rows)

! PC2 hom_conv diags (homog response to convection + erosion)
REAL   ::                                                         &
  qcl_hom_conv_incr_diag(row_length,rows,model_levels)            &
, cf_liquid_hom_conv_incr_diag(row_length,rows,model_levels)

REAL :: tmp       (row_length,rows)               ! Temporary array for 2d
REAL :: tmp3d     (row_length,rows,model_levels)  ! Temporary array for 3d
REAL :: tmp3d_u   (udims%i_start:udims%i_end,   & ! Temporary array for u wind
                   udims%j_start:udims%j_end,   &
                   udims%k_start:udims%k_end)
REAL :: tmp3d_v   (vdims%i_start:vdims%i_end,   & ! Temporary array for v wind
                   vdims%j_start:vdims%j_end,   &
                   vdims%k_start:vdims%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode = 0 ! Initialise error status
im_index = 1

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output.
! ----------------------------------------------------------------------
! Work out number of rows for U and V grids - will depend on whether
! ENDGame or not.

u_rows = udims%j_len
v_rows = vdims%j_len

! ----------------------------------------------------------------------
! 5 161 Temperature increments from convection only
! ----------------------------------------------------------------------
item = 161
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        t_incr_conv_only,                                         &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(t_incr_conv_only) 161"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  !
! ----------------------------------------------------------------------
! 5 162 q increments from convection only
! ----------------------------------------------------------------------
item = 162
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        q_incr_conv_only,                                         &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(q_incr_conv_only) 162"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  !

! ----------------------------------------------------------------------
! DIAG.05163 Copy qCL INC: conv inhom to stashwork
! ----------------------------------------------------------------------
item = 163
! Diag05163_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        qcl_incr_inhom_diag,                                      &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(Qcl INC: conv inhom)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05163_if1

! ----------------------------------------------------------------------
! DIAG.05140 Copy positive part of conv inhom qcl incr to stashwork
! ----------------------------------------------------------------------
item = 140
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,tmp3d,qcl_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        tmp3d(i,j,k) = MAX(0.0,qcl_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d,     &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(Qcl INC+: convection)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.05141 Copy negative part of conv inhom qcl incr to stashwork
! ----------------------------------------------------------------------
item = 141
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,tmp3d,qcl_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        tmp3d(i,j,k) = MIN(0.0,qcl_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d,     &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(Qcl INC-: convection)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.05164 Copy qCF INC: conv inhom to stashwork
! ----------------------------------------------------------------------
item = 164
! Diag05164_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        qcf_incr_inhom_diag,                                      &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(Qcf INC: conv inhom)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05164_if1

! ----------------------------------------------------------------------
! DIAG.05142 Copy positive part of conv inhom qcf incr to stashwork
! ----------------------------------------------------------------------
item = 142
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,tmp3d,qcf_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        tmp3d(i,j,k) = MAX(0.0,qcf_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d,             &
         row_length,rows,model_levels,0,0,0,0, at_extremity,              &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
         stash_levels,num_stash_levels+1,                                 &
         atmos_im,sect,item,                                              &
         icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(Qcf INC+: conv inhom)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.05143 Copy negative part of conv inhom qcf incr to stashwork
! ----------------------------------------------------------------------
item = 143
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,tmp3d,qcf_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        tmp3d(i,j,k) = MIN(0.0,qcf_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d,           &
         row_length,rows,model_levels,0,0,0,0, at_extremity,            &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,            &
         stash_levels,num_stash_levels+1,                               &
         atmos_im,sect,item,                                            &
         icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(Qcf INC-: conv inhom)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.05172 Copy BULK_CF INC: conv inhom to stashwork
! ----------------------------------------------------------------------
item = 172
! Diag05172_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        bulk_cf_incr_inhom_diag,                                  &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(BCF INC: conv inhom)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05172_if1

! ----------------------------------------------------------------------
! DIAG.05173 Copy CF_LIQUID INC: conv inhom to stashwork
! ----------------------------------------------------------------------
item = 173
! Diag05173_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        cf_liquid_incr_inhom_diag,                                &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(CFl INC: conv inhom)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05173_if1

! ----------------------------------------------------------------------
! DIAG.05146 Copy positive part of conv inhom cfl incr to stashwork
! ----------------------------------------------------------------------
item = 146
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,tmp3d,cf_liquid_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        tmp3d(i,j,k) = MAX(0.0,cf_liquid_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d,     &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(CFl INC+: conv inhom)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.05147 Copy negative part of conv inhom cfl incr to stashwork
! ----------------------------------------------------------------------
item = 147
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,tmp3d,cf_liquid_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        tmp3d(i,j,k) = MIN(0.0,cf_liquid_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d,     &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(CFl INC-: conv inhom)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.05174 Copy CF_FROZEN INC: conv inhom to stashwork
! ----------------------------------------------------------------------
item = 174
! Diag05174_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        cf_frozen_incr_inhom_diag,                                &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(CFf INC: conv inhom)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05174_if1

! ----------------------------------------------------------------------
! DIAG.05148 Copy positive part of conv inhom cff incr to stashwork
! ----------------------------------------------------------------------
item = 148
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,tmp3d,cf_frozen_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        tmp3d(i,j,k) = MAX(0.0,cf_frozen_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d,     &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(CFf INC+: conv inhom)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.05149 Copy negative part of conv inhom cff incr to stashwork
! ----------------------------------------------------------------------
item = 149
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,tmp3d,cf_frozen_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        tmp3d(i,j,k) = MIN(0.0,cf_frozen_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d,     &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(CFf INC-: conv inhom)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! increment diagnostics= modified - previous

item = 181  ! temperature increment
IF (icode <= 0 .AND. (sf(item,sect) .OR. l_ac)) THEN

  IF (.NOT. ALLOCATED(tinc_cvn)) THEN
    ALLOCATE ( tinc_cvn(row_length*rows,model_levels) )
    tinc_cvn(1,1) = rmdi
  END IF
END IF

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       t_incr_diag_conv,                                         &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 181)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i,ji)                                                       &
!$OMP SHARED(model_levels,rows,row_length,tinc_cvn,t_incr_diag_conv)
  DO k = 1,model_levels
    DO j = 1,rows
      DO i = 1,row_length
        ji = (j-1)*row_length+i
        tinc_cvn(ji,k) = t_incr_diag_conv(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF  !  sf(item,sect)

item = 182  ! humidity increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       q_incr_diag_conv,                                         &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 182)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.05183 Copy qCL INC: convection to stashwork
! ----------------------------------------------------------------------
item = 183
! Diag05183_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        qcl_incr_diag_conv,                                       &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(Qcl INC: convection)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05183_if1

! ----------------------------------------------------------------------
! DIAG.05184 Copy qCF INC: convection to stashwork
! ----------------------------------------------------------------------
item = 184
! Diag05184_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        qcf_incr_diag_conv,                                       &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(Qcf INC: convection)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05184_if1

! ----------------------------------------------------------------------
! DIAG.05185 Copy U Wind INC: convection to stashwork
! ----------------------------------------------------------------------
item = 185  ! u wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(udims,tmp3d_u,dubydt_conv_u,timestep)
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        tmp3d_u(i,j,k) = dubydt_conv_u(i,j,k) * timestep
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d_u,  &
       row_length,u_rows,model_levels,0,0,0,0, at_extremity,     &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 185)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

item = 186  ! v wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(vdims,tmp3d_v,dvbydt_conv_v,timestep)
  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        tmp3d_v(i,j,k) = dvbydt_conv_v(i,j,k) * timestep
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d_v,  &
       row_length,v_rows,model_levels,0,0,0,0, at_extremity,     &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 186)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)
! ----------------------------------------------------------------------
! convection source terms in non shallow regions
! ----------------------------------------------------------------------
item = 187  ! potential temperature/temperature increment
IF (icode <= 0 .AND. sf(item,sect)) THEN
  IF (l_fix_conv_diags_var) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          tmp3d(i,j,k) = t_incr_conv_only(i,j,k)/exner_theta_levels(i,j,k)    &
                                  *(1.0-shallow_ind(i,j))
        END DO
      END DO
    END DO
  ELSE
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,tmp3d,t_incr_diag_conv,shallow_ind)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          tmp3d(i,j,k) = t_incr_diag_conv(i,j,k)*(1.0-shallow_ind(i,j))
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),tmp3d,     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d(item 187)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

item = 188  ! humidity increment
IF (icode <= 0 .AND. sf(item,sect)) THEN
  IF (l_fix_conv_diags_var) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          tmp3d(i,j,k) = q_incr_conv_only(i,j,k)*(1.0-shallow_ind(i,j))
        END DO
      END DO
    END DO
  ELSE
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,tmp3d,q_incr_diag_conv,shallow_ind)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          tmp3d(i,j,k) = q_incr_diag_conv(i,j,k)*(1.0-shallow_ind(i,j))
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  CALL copydiag_3d (stashwork(si(item,sect,im_index)),tmp3d,     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode > 0) THEN
    cmessage=": error in copydiag_3d(item 188)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.05192 Copy BULK_CF INC: convection to stashwork
! ----------------------------------------------------------------------
item = 192
! Diag05192_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        bulk_cf_incr_diag_conv,                                   &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(BCF INC: convection)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05192_if1

! ----------------------------------------------------------------------
! DIAG.05193 Copy CF_LIQUID INC: convection to stashwork
! ----------------------------------------------------------------------
item = 193
! Diag05193_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        cf_liquid_incr_diag_conv,                                 &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(CFl INC: convection)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05193_if1

! ----------------------------------------------------------------------
! DIAG.05194 Copy CF_FROZEN INC: convection to stashwork
! ----------------------------------------------------------------------
item = 194
! Diag05194_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        cf_frozen_incr_diag_conv,                                 &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(CFf INC: convection)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05194_if1



!-----------------------------------------------------------------------
! W-EQN DIAGS
!-----------------------------------------------------------------------

! ----------------------------------------------------------------------
! DIAG.05196 Convective W^2 of parcel
! ----------------------------------------------------------------------
item = 196

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,sect,im_index)),            &
       w2p,                                                      &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diagnostics_conv  : error in copydiag_3d(w2p: convection)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.05197 Convective W of parcel (where w2p > 0.0)
! ----------------------------------------------------------------------
item = 197
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  tmp3d(:,:,:) = 0.0
  WHERE (w2p(:,:,:) > 0.0) tmp3d(:,:,:) = SQRT(w2p(:,:,:))

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,sect,im_index)), tmp3d,     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diagnostics_conv  : error in copydiag_3d(wp: convection)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

!-----------------------------------------------------------------------
! W-EQN DIAGS
!-----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Section: PC2 (hom_conv) homog response to conv plus erosion
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! DIAG.05150 Copy hom_conv qcl incr to stashwork
! ----------------------------------------------------------------------
item = 150
! Diag05150_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,qcl_hom_conv_incr_diag,             &
!$OMP        qcl_incr_diag_conv,qcl_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        qcl_hom_conv_incr_diag(i,j,k)=                            &
              qcl_incr_diag_conv(i,j,k) - qcl_incr_inhom_diag(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  !
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        qcl_hom_conv_incr_diag,                                   &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)
  !
  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(Qcl INC: hom_conv)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05150_if1

! ----------------------------------------------------------------------
! DIAG.05151 Copy positive part of hom_conv qcl incr to stashwork
! ----------------------------------------------------------------------
item = 151
! Diag05151_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,qcl_hom_conv_incr_diag,             &
!$OMP        qcl_incr_diag_conv,qcl_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        qcl_hom_conv_incr_diag(i,j,k)=MAX(0.0,                    &
              qcl_incr_diag_conv(i,j,k)                           &
              -qcl_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  !
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        qcl_hom_conv_incr_diag,                                   &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)
  !
  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(Qcl INC: hom_conv pos)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05151_if1

! ----------------------------------------------------------------------
! DIAG.05152 Copy negative part of hom_conv qcl incr to stashwork
! ----------------------------------------------------------------------
item = 152
! Diag05152_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,qcl_hom_conv_incr_diag,             &
!$OMP        qcl_incr_diag_conv,qcl_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        qcl_hom_conv_incr_diag(i,j,k)=MIN(0.0,                    &
              qcl_incr_diag_conv(i,j,k)                           &
              -qcl_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  !
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        qcl_hom_conv_incr_diag,                                   &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        icode,cmessage)
  !
  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(Qcl INC: hom_conv neg)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05152_if1

! ----------------------------------------------------------------------
! DIAG.05156 Copy hom_conv cfl incr to stashwork
! ----------------------------------------------------------------------
item = 156
! Diag05156_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,cf_liquid_hom_conv_incr_diag,       &
!$OMP        cf_liquid_incr_diag_conv,cf_liquid_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        cf_liquid_hom_conv_incr_diag(i,j,k)=                      &
              cf_liquid_incr_diag_conv(i,j,k)                     &
              -cf_liquid_incr_inhom_diag(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  !
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
       cf_liquid_hom_conv_incr_diag,                              &
       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
       stash_levels,num_stash_levels+1,                           &
       atmos_im,sect,item,                                        &
       icode,cmessage)
  !
  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(CFl INC: hom_conv)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05156_if1

! ----------------------------------------------------------------------
! DIAG.05157 Copy positive part of hom_conv cfl incr to stashwork
! ----------------------------------------------------------------------
item = 157
! Diag05157_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,cf_liquid_hom_conv_incr_diag,       &
!$OMP        cf_liquid_incr_diag_conv,cf_liquid_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        cf_liquid_hom_conv_incr_diag(i,j,k)=MAX(0.0,              &
              cf_liquid_incr_diag_conv(i,j,k)                     &
              -cf_liquid_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  !
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
       cf_liquid_hom_conv_incr_diag,                              &
       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
       stash_levels,num_stash_levels+1,                           &
       atmos_im,sect,item,                                        &
       icode,cmessage)
  !
  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(CFl INC: hom_conv pos)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05157_if1

! ----------------------------------------------------------------------
! DIAG.05158 Copy negative part of hom_conv cfl incr to stashwork
! ----------------------------------------------------------------------
item = 158
! Diag05158_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,cf_liquid_hom_conv_incr_diag,       &
!$OMP        cf_liquid_incr_diag_conv,cf_liquid_incr_inhom_diag)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        cf_liquid_hom_conv_incr_diag(i,j,k)=MIN(0.0,              &
              cf_liquid_incr_diag_conv(i,j,k)                     &
              -cf_liquid_incr_inhom_diag(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  !
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
       cf_liquid_hom_conv_incr_diag,                              &
       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
       stash_levels,num_stash_levels+1,                           &
       atmos_im,sect,item,                                        &
       icode,cmessage)
  !
  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(CFl INC: hom_conv neg)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag05157_if1

! ----------------------------------------------------------------------
! Downdraughts
! ----------------------------------------------------------------------
item = 198  !  dT/dt from downdraughts

IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
       dt_dd,                                                     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
       stash_levels,num_stash_levels+1,                           &
       atmos_im,sect,item,                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(for dT_dd)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

item = 199  !  dq/dt from downdraughts

IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
       dq_dd,                                                     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
       stash_levels,num_stash_levels+1,                           &
       atmos_im,sect,item,                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(for dq_dd)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

item = 175  !  du/dt from downdraughts on p_grid

IF (icode  <=  0  .AND.  sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
       du_dd,                                                     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
       stash_levels,num_stash_levels+1,                           &
       atmos_im,sect,item,                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(for du_dd)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

item = 176  !  dv/dt from downdraughts on p_grid

IF (icode  <=  0  .AND.  sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
       dv_dd,                                                     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
       stash_levels,num_stash_levels+1,                           &
       atmos_im,sect,item,                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(for dv_dd)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! Section  Convective Rain
! ----------------------------------------------------------------------
! Item 201 Convective rainfall,resolve to accumulate over timestep

item = 201 ! convective rainfall accumulation
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(201,5,im_index)),                   &
       conv_rain*timestep,row_length,rows,                       &
       0,0,0,0, at_extremity,atmos_im,5,201,                     &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag(for 5201)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! Section  Convective Snow
! ----------------------------------------------------------------------
! Item 202 Convective Snowfall, resolve to accumulate overtimestep.

item = 202 ! convective snowfall accumulation
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(202,5,im_index)),                   &
       conv_snow*timestep,row_length,rows,                       &
       0,0,0,0, at_extremity,atmos_im,5,202,                     &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag(for 5202)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! Section  Convective rainfall rates
! ----------------------------------------------------------------------
! C Item 205 Convective rainfall rates

! this is also used by stash 20014 precip_symbol
IF (icode <=0 .AND. (sf(stashcode_pws_precip_sym, stashcode_pws_sec))) THEN
  pws_precip_sym_conv_rain(:,:) = conv_rain(:,:)
END IF

item = 205 ! convective rainfall rate
IF (icode <= 0 .AND. (sf(item,sect) .OR. l_ac)) THEN

  IF (.NOT. ALLOCATED(cvrr)) THEN
    ALLOCATE ( cvrr(row_length*rows) )
    cvrr(1) = rmdi
  END IF
END IF

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(205,5,im_index)),conv_rain,         &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,205,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag(for 5205)"
    CALL ereport(routinename,icode,cmessage)
  END IF

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i,ji)                                                         &
!$OMP SHARED(rows,row_length,cvrr,conv_rain)
  DO j = 1,rows
    DO i = 1,row_length
      ji = (j-1)*row_length+i
      cvrr(ji) = conv_rain(i,j)
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

! ----------------------------------------------------------------------
! Section  206 Convective snowfall rates
! ----------------------------------------------------------------------
! C Item 206 Convective snowfall rates

! this is also used by stash 20014 precip_symbol
IF (icode <=0 .AND. (sf(stashcode_pws_precip_sym, stashcode_pws_sec))) THEN
  pws_precip_sym_conv_snow(:,:) = conv_snow(:,:)
END IF

item = 206 ! convective snowfall rate
IF (icode <= 0 .AND. (sf(item,sect) .OR. l_ac)) THEN

  IF (.NOT. ALLOCATED(cvsr)) THEN
    ALLOCATE ( cvsr(row_length*rows) )
  END IF
END IF

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(206,5,im_index)),conv_snow,         &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,206,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl  : error in copydiag(for 5206)"
    CALL ereport(routinename,icode,cmessage)
  END IF

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i,ji)                                                         &
!$OMP SHARED(rows,row_length,cvsr,conv_snow)
  DO j = 1,rows
    DO i = 1,row_length
      ji = (j-1)*row_length+i
      cvsr(ji) = conv_snow(i,j)
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

! ----------------------------------------------------------------------
! Section Convective Cloud Amount
! ----------------------------------------------------------------------
! Item 212

item = 212 ! convective cloud amount
IF (icode <= 0 .AND. sf(item,sect)) THEN

  IF (l_3d_cca) THEN
    IF (l_ac) THEN
      IF (.NOT. ALLOCATED(convcc)) THEN
        ALLOCATE ( convcc(row_length*rows,n_cca_levels) )
        convcc(1,1) = rmdi
      END IF
    END IF

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i,ji)                                                       &
!$OMP SHARED(n_cca_levels,rows,row_length,tmp3d,cca,l_ac,convcc)
    DO k = 1, n_cca_levels
      DO j = 1, rows
        DO i = 1, row_length
          tmp3d(i,j,k) = cca(i,j,k)
          IF (l_ac) THEN
            ji = (j-1)*row_length+i
            convcc(ji,k) = tmp3d(i,j,k)
          END IF
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d, &
         row_length,rows,n_cca_levels,0,0,0,0, at_extremity,    &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,    &
         stash_levels,num_stash_levels+1,                       &
         atmos_im,sect,item,                                    &
         icode,cmessage)
  ELSE
    IF (l_ac) THEN
      IF (.NOT. ALLOCATED(convcc)) THEN
        ALLOCATE ( convcc(row_length*rows,model_levels) )
      END IF
    END IF

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i,ji)                                                       &
!$OMP SHARED(model_levels,rows,row_length,ccb,cct,tmp3d,cca,convcc,l_ac)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          IF (k  >=  ccb(i,j) .AND.                          &
              k  <   cct(i,j)) THEN
            tmp3d(i,j,k) = cca(i,j,1)
          ELSE
            tmp3d(i,j,k) = 0.0
          END IF
          IF (l_ac) THEN
            ji = (j-1)*row_length+i
            convcc(ji,k) = tmp3d(i,j,k)
          END IF
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork(si(item,sect,im_index)), tmp3d, &
         row_length,rows,model_levels,0,0,0,0, at_extremity,    &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,    &
         stash_levels,num_stash_levels+1,                       &
         atmos_im,sect,item,                                    &
         icode,cmessage)
  END IF

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 212)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section Convective Cloud Liquid Water
! ----------------------------------------------------------------------
! Item 213

item = 213 ! convective cloud liquid water
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), ccw,      &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 213)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section Total precipitation
! ----------------------------------------------------------------------
! Item 226

item = 226 ! total precipitation accumulation
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,tmp,ls_rain,ls_snow,conv_rain,conv_snow,         &
!$OMP        timestep)
  DO j = 1, rows
    DO i = 1, row_length
      tmp(i,j) = ( ls_rain(i,j) + ls_snow(i,j)                   &
                 + conv_rain(i,j) + conv_snow(i,j) )             &
                       * timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(226,5,im_index)),tmp,               &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,226,                                           &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5226)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! Section Total Rain rate
! ----------------------------------------------------------------------
! Item 214 ls_rain +conv_rain

item = 214 ! total rain rate
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,tmp,ls_rain,conv_rain)
  DO j = 1, rows
    DO i = 1, row_length
      tmp(i,j) = (ls_rain(i,j) +conv_rain(i,j) )
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(214,5,im_index)),tmp,               &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,214,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5214)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section Total Snow rate
! ----------------------------------------------------------------------
! Item 215 ls_snow + conv_snow

item = 215 ! total snow rate
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,tmp,ls_snow,conv_snow)
  DO j = 1, rows
    DO i = 1, row_length
      tmp(i,j) = (ls_snow(i,j) + conv_snow(i,j) )
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(215,5,im_index)),tmp,               &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,215,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5215)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! Section Total Precipitation rate
! ----------------------------------------------------------------------
! Item 216 ls_rain + conv_rain + ls_snow + conv_snow

item = 216 ! total precipitation rate
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,tmp,ls_rain,conv_rain,ls_snow,conv_snow)
  DO j = 1, rows
    DO i = 1, row_length
      tmp(i,j) = (ls_rain(i,j) + conv_rain(i,j) +      &
                  ls_snow(i,j) + conv_snow(i,j))
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(216,5,im_index)),tmp,               &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,216,                                           &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5216)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! Calculate pressure at convective cloud base
! ----------------------------------------------------------------------

IF (icode <= 0 .AND. (sf(207,sect) .OR. sf(210,sect))) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,tmp,p,ccb)
  DO j = 1, rows
    DO i = 1, row_length
      IF (ccb(i,j)  /=  0 ) THEN
        tmp(i,j) = p(i,j,ccb(i,j))
      ELSE
        tmp(i,j) =  rmdi
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO
END IF


! ----------------------------------------------------------------------
! Section pressure at convective cloud base
! ----------------------------------------------------------------------
! Item 207 Pressure at convective cloud base

item = 207 ! Pressure at convective cloud base
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(207,5,im_index)),tmp,               &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,207,                                           &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5207)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Calculate height at convective cloud base
!  (ICAO standard atmosphere heights)
! ----------------------------------------------------------------------

IF (icode <= 0 .AND. sf(210,sect)) THEN
  CALL icao_ht(tmp,row_length,rows,icao_height)

! this is also used by stash 20012 convective cloud depth
! and by stash 20039 conv cld base icao ht
! but use ICAOHeight_FC routine derived from fieldcalc
  IF (sf(stashcode_pws_conv_cld_dep,stashcode_pws_sec) .OR. &
      sf(stashcode_pws_conv_icao_base,stashcode_pws_sec)) THEN
    CALL ICAOHeight_FC(tmp,pws_conv_icao_base,row_length,rows)
  END IF

END IF

! ----------------------------------------------------------------------
! Section ICAO height at convective cloud base
! ----------------------------------------------------------------------
! Item 210 ICAO height at convective cloud base

item = 210 ! ICAO height at convective cloud base
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(210,5,im_index)),icao_height,       &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,210,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5210)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Calculate pressure at convective cloud top
! ----------------------------------------------------------------------

IF (icode <= 0 .AND. (sf(208,sect) .OR. sf(211,sect))) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,tmp,p,cct)
  DO j = 1, rows
    DO i = 1, row_length
      IF (cct(i,j)  /=  0 ) THEN
        tmp(i,j) = p(i,j,cct(i,j))
      ELSE
        tmp(i,j) =  rmdi
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! ----------------------------------------------------------------------
! Section pressure at convective cloud top
! ----------------------------------------------------------------------
! Item 208 Pressure at convective cloud top

item = 208 ! Pressure at convective cloud top
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(208,5,im_index)),tmp,               &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,208,                                           &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5208)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! Calculate height at convective cloud top
!  (ICAO standard atmosphere heights)
! ----------------------------------------------------------------------

IF (icode <= 0 .AND. sf(211,sect)) THEN
  CALL icao_ht(tmp,row_length,rows,icao_height)

! this is also used by stash 20012 convective cloud depth
! and by stash 20040 conv cld top icao ht
! but use ICAOHeight_FC routine derived from fieldcalc
  IF (sf(stashcode_pws_conv_cld_dep,stashcode_pws_sec) .OR. &
      sf(stashcode_pws_conv_icao_top,stashcode_pws_sec)) THEN
    CALL ICAOHeight_FC(tmp,pws_conv_icao_top,row_length,rows)
  END IF

END IF

! ----------------------------------------------------------------------
! Section ICAO height at convective cloud top
! ----------------------------------------------------------------------
! Item 211 ICAO height at convective cloud top

item = 211 ! ICAO height at convective cloud top
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(211,5,im_index)),icao_height,       &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,211,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5211)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  CAPE
! ----------------------------------------------------------------------
! Item 217 CAPE

item=217
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(217,5,im_index)),cape_out,           &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,217,                                            &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5217)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section Lowest Convective Cloud Base
! ----------------------------------------------------------------------
! Item 218 Lowest convective cloud base

item = 218 ! Lowest convective cloud base
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(218,5,im_index)),lcbase,            &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,218,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5218)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section Lowest Convective Cloud Top
! ----------------------------------------------------------------------
! Item 219 Lowest convective cloud top

item = 219 ! Lowest convective cloud top
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(219,5,im_index)),lctop,             &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,219,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5219)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section Lowest Convective Cloud Amount
! ----------------------------------------------------------------------
! Item 220 Lowest convective cloud amount

item = 220 ! Lowest convective cloud amount
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(220,5,im_index)),lcca,              &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,220,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5220)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Calculate pressure at lowest convective cloud base
! ----------------------------------------------------------------------

IF (icode <= 0 .AND. (sf(222,sect) .OR. sf(224,sect))) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,tmp,lcbase,p)
  DO j = 1, rows
    DO i = 1, row_length
      IF (lcbase(i,j)  /=  0) THEN
        tmp(i,j) = p(i,j,lcbase(i,j))
      ELSE
        tmp(i,j) = rmdi
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! ----------------------------------------------------------------------
! Section pressure at lowest convective cloud base
! ----------------------------------------------------------------------
! Item 222 Pressure at lowest convective cloud base

item = 222 ! Pressure at lowest convective cloud base
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(222,5,im_index)),tmp,               &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,222,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5222)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Calculate height at convective cloud base
!  (ICAO standard atmosphere heights)
! ----------------------------------------------------------------------

IF (icode <= 0 .AND. sf(224,sect)) THEN
  CALL icao_ht(tmp,row_length,rows,icao_height)
END IF

! ----------------------------------------------------------------------
! Section ICAO height at lowest convective cloud base
! ----------------------------------------------------------------------
! Item 224 ICAO height at convective cloud base

item = 224 ! ICAO height at convective cloud base
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(224,5,im_index)),icao_height,       &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,224,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5224)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Calculate pressure at convective cloud top
! ----------------------------------------------------------------------

IF (icode <= 0 .AND. (sf(223,sect) .OR. sf(225,sect))) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,tmp,lctop,p)
  DO j = 1, rows
    DO i = 1, row_length
      IF (lctop(i,j)  /=  0) THEN
        tmp(i,j) = p(i,j,lctop(i,j))
      ELSE
        tmp(i,j) = rmdi
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! ----------------------------------------------------------------------
! Section pressure at lowest convective cloud top
! ----------------------------------------------------------------------
! Item 223 Pressure at lowest convective cloud top

item = 223 ! Pressure at lowest convective cloud top
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(223,5,im_index)),tmp,               &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,223,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5223)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Calculate ICAO height at convective cloud top
! ----------------------------------------------------------------------

IF (icode <= 0 .AND. sf(225,sect)) THEN
  CALL icao_ht(tmp,row_length,rows,icao_height)
END IF

! ----------------------------------------------------------------------
! Section ICAO height at lowest convective cloud top
! ----------------------------------------------------------------------
! Item 225 ICAO height at convective cloud top

item = 225 ! ICAO height at convective cloud top
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(225,5,im_index)),icao_height,       &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,5,225,                                           &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag(item 5225)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

!-----------------------------------------------------------------------
! Section 3D Convective rainfall rate
! Item 227
! ----------------------------------------------------------------------

item = 227                         ! 3D convective rainfall rate

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       conv_rain_3d,                                             &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 227)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section 3D Convective snowfall rate
! Item 228
! ----------------------------------------------------------------------

item = 228                        ! 3D convective snowfall rate

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       conv_snow_3d,                                             &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 228)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF
! ----------------------------------------------------------------------
! Fractional areas
! ----------------------------------------------------------------------

item = 229                        ! fractional updraught area

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       area_ud,                                                  &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 229)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

item = 230                        ! fractional downdraught area

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       area_dd,                                                  &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 230)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section 2D convective cloud amount
! ----------------------------------------------------------------------
item = 262 ! 2D convective cloud amount

!-----------------------------------------------------------------------
! Not sure if This switch is required as it should not
! affect evolution, i.e. a diagnostic only.
!-----------------------------------------------------------------------

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! This diagnostic appears to try and back out the CCA_2D from
  ! the original code before the anvil (and presumably any
  ! radiative cca tunings) was applied, it appears to assume
  ! that convection has occurred if cca at cct > cca at ccb. It
  ! is messy, confusing and unesseccary, may mean holding one
  ! more single level diag.

  IF (l_dcpl_cld4pc2) THEN
    ! This will change cca_2d results, and any associated
    ! diagnostics. Simpler, essentially it is the cca_2d, that
    ! was used to generate 3d-cca profiles. However this will
    ! differ from the original in the case of shallow convection
    ! in that the cca_2d will be the cca_2d used to generate
    ! the shallow cloud. (shallow cloud may not be uniform with
    ! height)

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,tmp,cca_2d)
    DO j=1, rows
      DO i=1, row_length
        tmp(i,j) = MIN(max_2d_cca, cca_2d(i,j))
      END DO
    END DO
!$OMP END PARALLEL DO

  ELSE
    ! Difficult to say what was intended by this diagnostic. It
    ! is backed out from the 3d cca field and backs out cca_2d
    ! from the cca at cloud base, which would  have had the
    ! tower_factor applied.

    ! It looks to take the cca at cloud base before the anvil
    ! was applied (or any tunings). It looks to be similar to
    ! what the diag lcca is doing but appears to be returning
    ! the upper convective layer, certainly the description

    ! 2D convective cloud amount is vague.

    ! Because it it only using cca, ccb, cct it is working of
    ! the highest cloud layer.
    IF (l_3d_cca) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,tmp,ccb,cct,cca,tower_factor)
      DO j=1, rows
        DO i=1, row_length
          ! Convert 3d cca to 2d cca
          IF ( ccb(i,j)  ==  0 .AND. cct(i,j)  ==  0 ) THEN
            tmp(i,j) = 0.0
          ELSE IF ( cca(i,j,cct(i,j)-1)  >   cca(i,j,ccb(i,j))) THEN
            IF (tower_factor  >   0.0) THEN
              ! Back out the original 2D convective cloud amount
              ! before the anvil parametrization was applied
              tmp(i,j) = cca(i,j,ccb(i,j))  / tower_factor
            ELSE
              ! Tower_factor is zero (likely to be a PC2 run)
              ! and hence assume that no convective cloud is
              ! required.
              tmp(i,j) = 0.0
            END IF  ! tower_factor  >   0.0
          ELSE
            tmp(i,j)=cca(i,j,ccb(i,j))
          END IF


          !   Ensure that the maximum 2d cca isn't breached
          IF ( tmp(i,j)  >   max_2d_cca ) THEN
            tmp(i,j) = max_2d_cca
          END IF
        END DO
      END DO
!$OMP END PARALLEL DO
    ELSE
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,tmp,cca)
      DO j=1, rows
        DO i=1, row_length
          ! Copy 2d cca
          tmp(i,j) = cca(i,j,1)
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF

  END IF ! l_dcpl_cld4pc2

  ! Item 262 2D convective cloud amount
  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),tmp,           &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,sect,item,                                       &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag (5262)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! Section 1.2  Tendency diagnostics
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section Temperature after convection.
! ----------------------------------------------------------------------
! Item 209 Temperature

item = 209 ! temperature
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,tmp3d,theta_star,exner_theta_levels)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        tmp3d(i,j,k) = theta_star(i,j,k)*exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(209,5,im_index)),tmp3d,          &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,209,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,209,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item 209)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section q after convection.
! ----------------------------------------------------------------------
! Item 010 q

item =  10 ! q
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(010,5,im_index)),                &
       q_star(:,:,1:model_levels),                               &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,010,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,010,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item 010)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  Item 231 - CAPE timescale used if deep convection
! ----------------------------------------------------------------------

item=231
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),cape_ts_diag,      &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag (5231)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
!  Item 232 - indicator of reduced CAPE timescale
! ----------------------------------------------------------------------

item=232
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),ind_cape_reduced_diag,  &
       row_length,rows,0,0,0,0, at_extremity,                          &
       atmos_im,5,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag (5232)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  233 undilute CAPE
! ----------------------------------------------------------------------
! Item 233 undilute CAPE

item=233
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),cape_undilute,     &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,                                           &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag (5233)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  234 undilute CIN
! ----------------------------------------------------------------------
! Item 234 undilute CIN

item=234
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),cin_undilute,      &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,                                           &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl: error in copydiag (5234)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section u latest.
! ----------------------------------------------------------------------
! Item 235 u

item = 235 ! u
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(udims,tmp3d_u,u,r_u)
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        tmp3d_u(i,j,k) = u(i,j,k) + r_u(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(235,5,im_index)), tmp3d_u,       &
       row_length,u_rows,model_levels,0,0,0,0, at_extremity,     &
       stlist(1,stindex(1,235,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,235,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item 235)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section v latest
! ----------------------------------------------------------------------
! Item 236 v

item = 236 ! u
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(vdims,tmp3d_v,v,r_v)
  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        tmp3d_v(i,j,k) = v(i,j,k) + r_v(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(236,5,im_index)), tmp3d_v,       &
       row_length,v_rows,model_levels,0,0,0,0, at_extremity,     &
       stlist(1,stindex(1,236,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,236,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item 236)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! deep convec. massflux components (for 4D-Var)
! ----------------------------------------------------------------------
! item 246 rho-level massflux in a) no shallow regions and
!                                b) predominantly entraining regions
!                                   (mass flux is growing with height).
! The condition on the mass flux increasing with height is relaxed 
! if l_fix_conv_diags_var=TRUE.

item=246 ! component of half level massflux
IF (icode == 0 .AND. sf(item,sect)) THEN
  IF (l_fix_conv_diags_var) THEN
    DO k = 1, model_levels-1
      DO j = 1, rows
        DO i = 1, row_length
          IF ((1.0-max_mf_fall)*up_flux_half(i,j,k) <               &
               up_flux_half(i,j,k+1)) THEN
            tmp3d(i,j,k) = up_flux_half(i,j,k)*(1.0-shallow_ind(i,j))
          ELSE
            tmp3d(i,j,k) =0.0
          END IF
        END DO
      END DO
    END DO
  ELSE !l_fix_conv_diags_var=F
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,up_flux_half,tmp3d,shallow_ind)
    DO k = 1, model_levels-1
      DO j = 1, rows
        DO i = 1, row_length
          IF (up_flux_half(i,j,k) < up_flux_half(i,j,k+1)) THEN
            tmp3d(i,j,k) = up_flux_half(i,j,k)*(1.0-shallow_ind(i,j))
          ELSE
            tmp3d(i,j,k) =0.0
          END IF
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF  !l_fix_conv_diags_var
  
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(model_levels,rows,row_length,tmp3d)
  DO j = 1, rows
    DO i = 1, row_length
      tmp3d(i,j,model_levels) = 0.0
    END DO
  END DO
!$OMP END PARALLEL DO

  CALL copydiag_3d(stashwork(si(item,5,im_index)),tmp3d,         &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="conv_ctl  : error IN copydiag_3d(mass flux component 246)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! item 249 updraught mass flux on half levels

item=249 ! updraught mass flux on half levels
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(249,5,im_index)),up_flux_half,   &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,249,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,249,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error IN copydiag_3d(up mass flux half levs)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! item 250 updraught mass flux

item=250 ! updraught mass flux
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(250,5,im_index)),up_flux,        &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,250,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,250,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(up mass flux)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! item 251 downdraught mass flux

item=251
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(251,5,im_index)),dwn_flux,       &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,251,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,251,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dwn mass flux)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! item 252 updraught entrainment rate

item=252
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(252,5,im_index)),entrain_up,     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,252,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,252,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(up entrainment)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! item 253 updraught detrainment rate

item=253
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(253,5,im_index)),detrain_up,     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,253,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,253,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(up detrainment)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! item 254 downdraught entrainment rate

item=254
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(254,5,im_index)),entrain_dwn,    &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,254,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,254,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dwn entrainment)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! item 255 downdraught detrainment rate

item=255
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(255,5,im_index)),detrain_dwn,    &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,255,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,255,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dwn detrainment)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! item 256 time rate of change of u on the p-grid

item=256
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(256,5,im_index)),dubydt_conv_p,  &
       row_length,rows,model_levels,0,0,                         &
       udims_s%halo_i, udims_s%halo_j, at_extremity,             &
       stlist(1,stindex(1,256,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,256,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dudt (p grid)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! item 257 time rate of change of v on the p-grid

item=257
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(257,5,im_index)),dvbydt_conv_p,  &
       row_length,rows,model_levels,0,0,                         &
       vdims_s%halo_i, vdims_s%halo_j, at_extremity,             &
       stlist(1,stindex(1,257,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,257,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dvdt (p grid)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! item 258 u-stress for deep convection (p-grid)

item=258
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(258,5,im_index)),uw_dp,          &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,258,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,258,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(uw deep)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! item 259 v-stress for deep convection (p-grid)

item=259
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(259,5,im_index)),vw_dp,          &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,259,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,259,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(vw deep)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! item 260 u stress for shallow convection (p-grid)

item=260
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(260,5,im_index)),uw_shall,       &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,260,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,260,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(uw shallow)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! item 261 v-stress for shallow convection (p-grid)

item=261
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(261,5,im_index)),vw_shall,       &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,261,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,261,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(vw shallow)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! item 263 u stress for mid-level convection (p-grid)

item=263
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),uw_mid,        &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(uw mid)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! item 264 v-stress for mid-level convection (p-grid)

item=264
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),vw_mid,        &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(vw mid)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

!----------------------------------------------------------------------
!  Item 267 CFL limited deep convection indicator
! ----------------------------------------------------------------------

item=267
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),deep_cfl_limited,  &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5267)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

!----------------------------------------------------------------------
! Item 268  CFL limited mid convection indicator
! ----------------------------------------------------------------------

item=268
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),mid_cfl_limited,   &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5268)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

!----------------------------------------------------------------------
!  Deep convection indicator
! ----------------------------------------------------------------------
! Item 269 Deep convection diagnostic

item=269
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),deep_ind,          &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5269)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

!----------------------------------------------------------------------
!  Shallow convection indicator
! ----------------------------------------------------------------------
! Item 270 Shallow convection diagnostic

item=270
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(270,5,im_index)),shallow_ind,           &
       row_length,rows,0,0,0,0, at_extremity,                        &
       atmos_im,5,270,                                               &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5270)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  Convection over orography indicator
! ----------------------------------------------------------------------
! Item 271 Convection over orography indicator

item=271
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(271,5,im_index)),cu_over_orog,       &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,271,                                            &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5271)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  Item 272 - indicator for Mid level convection
! ----------------------------------------------------------------------

item=272
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(272,5,im_index)),mid_ind,            &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,272,                                            &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5272)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 273 - NTML level number for top of mixing level
! ----------------------------------------------------------------------

item=273
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(273,5,im_index)),ntml_diag,          &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,273,                                            &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5273)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 274 -  NTPAR model level number for top of initial parcel ascent
! ----------------------------------------------------------------------

item=274
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(274,5,im_index)),ntpar_diag,         &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5274)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  Item 275 - Freeze_lev
! ----------------------------------------------------------------------

item=275
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(275,5,im_index)),freeze_diag,        &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,275,                                            &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5275)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 276 - kterm deep convection
! ----------------------------------------------------------------------

item=276
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(276,5,im_index)),kterm_diag,         &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,276,                                            &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5276)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 277 -  Total precipitation from deep convection
! ----------------------------------------------------------------------

item=277
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(277,5,im_index)),precip_deep,        &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,277,                                            &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5277)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 278 -  Total precipitation from shallow convection
! ----------------------------------------------------------------------

item=278
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(278,5,im_index)),precip_shall,       &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,278,                                            &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5278)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 279 -  Total precipitation from mid_level convection
! ----------------------------------------------------------------------

item=279
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(279,5,im_index)),precip_mid,         &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,279,                                            &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5279)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 280 -  Total precipitation from congestus convection
! ----------------------------------------------------------------------

item=280
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),precip_cong,       &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag(5280)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  Diagnostics only available from turbulence convection scheme
! ----------------------------------------------------------------------
! Item 290 -  w'qt' flux from turbulent convection
! ----------------------------------------------------------------------

item=290
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),wqt_flux_sh,   &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(wqt_flux)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 291 -  w'ql' flux from turbulent convection
! ----------------------------------------------------------------------

item=291
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),wql_flux_sh,   &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(wql_flux)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 292 -  w'thetal' flux from turbulent convection
! ----------------------------------------------------------------------

item=292
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),wthetal_flux_sh,&
       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                           &
       atmos_im,5,item,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(wthetal_flux)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 293 -  w'thetav' flux from turbulent convection
! ----------------------------------------------------------------------

item=293
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),wthetav_flux_sh,&
       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                           &
       atmos_im,5,item,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(wthetav_flux)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 300 -  subcloud layer convective velocity scale
! ----------------------------------------------------------------------

item=300
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),wstar_dn_diag,     &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5300)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 301 -  cumulus layer convective velocity scale
! ----------------------------------------------------------------------

item=301
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),wstar_up_diag,     &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5301)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 302 -   cloud base mass flux 1
! ----------------------------------------------------------------------

item=302
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),mb1_diag,          &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5302)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 303 -   cloud base mass flux 2
! ----------------------------------------------------------------------

item=303
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),mb2_diag,          &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5303)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 304 - wqt at cloud base
! ----------------------------------------------------------------------

item=304
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),wqt_cb,            &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5304)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 305 - wthetal at cloud base
! ----------------------------------------------------------------------

item=305
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),wthetal_cb,        &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5305)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 306 - wqt at inversion
! ----------------------------------------------------------------------

item=306
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),wqt_inv,           &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5306)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 307 - wthetal at inversion
! ----------------------------------------------------------------------

item=307
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),wthetal_inv,       &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5307)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 308 - height of top of shallow convection
!            Diagnostic = height * 1 if shallow convection otherwise
!            zero
! ----------------------------------------------------------------------

item=308
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),sh_top,            &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5308)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 309 - height of base of shallow convection
!            Diagnostic = height * 1 if shallow convection otherwise
!            zero
! ----------------------------------------------------------------------

item=309
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),sh_base,           &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5309)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  Item 310 - Congestus convection indicator
! ----------------------------------------------------------------------

item=310
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),congestus_ind,     &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5310)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  Item 311 - Congestus convection indicator  2
! ----------------------------------------------------------------------

item=311
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),congestus_ind2,    &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5311)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  Item 312 - Congestus termination level
! ----------------------------------------------------------------------

item=312
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),cg_term,           &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5312)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  Item 313 - Congestus top height
! ----------------------------------------------------------------------

item=313
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),cg_top,            &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5313)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
!  Item 314 - Congestus base height
! ----------------------------------------------------------------------

item=314
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,5,im_index)),cg_base,           &
       row_length,rows,0,0,0,0, at_extremity,                     &
       atmos_im,5,item,icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag (5314)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF
! ----------------------------------------------------------------------
! Item 319 - deep tops
!            Frequency deep convection terminates in model level k
!            Designed for use with some form of stash meaning to get
!            an ideal of the distribution of levels over which deep
!            convection terminates.
! ----------------------------------------------------------------------

item=319
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),deep_tops,     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(deep_tops)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 320 - mass flux deep
! ----------------------------------------------------------------------

item=320
IF (icode == 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),mf_deep,       &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(mf_deep)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 321 - mass flux congestus
! ----------------------------------------------------------------------

item=321
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),mf_congest,    &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(mf_congest)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 322 - mass flux shallow
! ----------------------------------------------------------------------

item=322
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),mf_shall,      &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(mf_shall)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! Item 323 - mass flux mid -level
! ----------------------------------------------------------------------

item=323
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),mf_midlev,     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(mf_midlev)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 324 - DT deep
! ----------------------------------------------------------------------

item=324
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dt_deep,       &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dt_deep)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 325 - dt congestus
! ----------------------------------------------------------------------

item=325
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dt_congest,    &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dt_congest)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 326 - dT shallow
! ----------------------------------------------------------------------

item=326
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dt_shall,      &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dt_shall)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! Item 327 - dT mid -level
! ----------------------------------------------------------------------

item=327
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dt_midlev,     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dt_midlev)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 328 - dq deep
! ----------------------------------------------------------------------

item=328
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dq_deep,       &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dq_deep)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 329 - dq congestus
! ----------------------------------------------------------------------

item=329
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dq_congest,    &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dq_congest)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 330 - dq shallow
! ----------------------------------------------------------------------

item=330
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dq_shall,      &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dq_shall)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 331 - dq mid -level
! ----------------------------------------------------------------------

item=331
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dq_midlev,     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dq_midlev)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 332 - du deep
! ----------------------------------------------------------------------

item=332
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),du_deep,       &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(du_deep)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 333 - du congestus
! ----------------------------------------------------------------------

item=333
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),du_congest,    &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(du_congest)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! Item 334 - du shallow
! ----------------------------------------------------------------------

item=334
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),du_shall,      &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(du_shall)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 335 - du mid -level
! ----------------------------------------------------------------------

item=335
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),du_midlev,     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(du_midlev)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 336 - dv deep
! ----------------------------------------------------------------------

item=336
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dv_deep,       &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dv_deep)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 337 - dv congestus
! ----------------------------------------------------------------------

item=337
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dv_congest,    &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dv_congest)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 338 - dv shallow
! ----------------------------------------------------------------------

item=338
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dv_shall,      &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dv_shall)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Item 339 - dv mid -level
! ----------------------------------------------------------------------

item=339
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(item,5,im_index)),dv_midlev,     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,5,item,                                          &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="conv_ctl  : error in copydiag_3d(dv_midlev)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section Mineral Dust scavenging diagnostics
! ----------------------------------------------------------------------

IF (l_dust) THEN

  ! Determine if total dust deposition flux is required
  l_tddf = sf(440,3)

  item = 281
  IF ( (icode <= 0) .AND. (l_tddf .OR. sf(item,sect)) ) THEN

    !  Convert scavenged dust per timestep to flux per sec
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_dust,timestep)
    DO j=1,rows
      DO i=1,row_length
        conscav_dust(i,j,1) = conscav_dust(i,j,1)/timestep
      END DO
    END DO
!$OMP END PARALLEL DO

    ! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(item,sect,im_index)),conscav_dust(1,1,1),&
        row_length,rows,0,0,0,0,at_extremity,                           &
        atmos_im,sect,item,                                             &
        icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": ERROR IN COPYDIAG(ITEM 281)"//cmessage
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF


  item = 282
  IF ( (icode <= 0) .AND. (l_tddf .OR. sf(item,sect)) ) THEN

    !  Convert scavenged dust per timestep to flux per sec

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_dust,timestep)
    DO j=1,rows
      DO i=1,row_length
        conscav_dust(i,j,2) = conscav_dust(i,j,2)/timestep
      END DO
    END DO
!$OMP END PARALLEL DO

    ! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(item,sect,im_index)),conscav_dust(1,1,2),&
        row_length,rows,0,0,0,0,at_extremity,                           &
        atmos_im,sect,item,                                             &
        icode,cmessage)

    IF (icode >  0) THEN
      cmessage=": ERROR IN COPYDIAG(ITEM 282)"//cmessage
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF

  IF (.NOT. l_twobin_dust) THEN      ! Six-bin dust scheme only

    item = 283
    IF ( (icode <= 0) .AND. (l_tddf .OR. sf(item,sect)) ) THEN

      !  Convert scavenged dust per timestep to flux per sec
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_dust,timestep)
      DO j=1,rows
        DO i=1,row_length
          conscav_dust(i,j,3) = conscav_dust(i,j,3)/timestep
        END DO
      END DO
!$OMP END PARALLEL DO

      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),conscav_dust(1,1,3),&
          row_length,rows,0,0,0,0,at_extremity,                           &
          atmos_im,sect,item,                                             &
          icode,cmessage)

      IF (icode >  0) THEN
        cmessage=": ERROR IN COPYDIAG(ITEM 283)"//cmessage
        CALL ereport(routinename,icode,cmessage)
      END IF

    END IF


    item = 284
    IF ( (icode <= 0) .AND. (l_tddf .OR. sf(item,sect)) ) THEN

      !  Convert scavenged dust per timestep to flux per sec

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_dust,timestep)
      DO j=1,rows
        DO i=1,row_length
          conscav_dust(i,j,4) = conscav_dust(i,j,4)/timestep
        END DO
      END DO
!$OMP END PARALLEL DO

      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),conscav_dust(1,1,4),&
          row_length,rows,0,0,0,0,at_extremity,                           &
          atmos_im,sect,item,                                             &
          icode,cmessage)

      IF (icode >  0) THEN
        cmessage=": ERROR IN COPYDIAG(ITEM 284)"//cmessage
        CALL ereport(routinename,icode,cmessage)
      END IF

    END IF


    item = 285
    IF ( (icode <= 0) .AND. (l_tddf .OR. sf(item,sect)) ) THEN

      !  Convert scavenged dust per timestep to flux per sec

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_dust,timestep)
      DO j=1,rows
        DO i=1,row_length
          conscav_dust(i,j,5) = conscav_dust(i,j,5)/timestep
        END DO
      END DO
!$OMP END PARALLEL DO

      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),conscav_dust(1,1,5),&
          row_length,rows,0,0,0,0,at_extremity,                           &
          atmos_im,sect,item,                                             &
          icode,cmessage)

      IF (icode >  0) THEN
        cmessage=": ERROR IN COPYDIAG(ITEM 285)"//cmessage
        CALL ereport(routinename,icode,cmessage)
      END IF

    END IF


    item = 286
    IF ( (icode <= 0) .AND. (l_tddf .OR. sf(item,sect)) ) THEN

      !  Convert scavenged dust per timestep to flux per sec
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_dust,timestep)
      DO j=1,rows
        DO i=1,row_length
          conscav_dust(i,j,6) = conscav_dust(i,j,6)/timestep
        END DO
      END DO
!$OMP END PARALLEL DO

      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),conscav_dust(1,1,6),&
          row_length,rows,0,0,0,0,at_extremity,                           &
          atmos_im,sect,item,                                             &
          icode,cmessage)

      IF (icode >  0) THEN
        cmessage=": ERROR IN COPYDIAG(ITEM 286)"//cmessage
        CALL ereport(routinename,icode,cmessage)
      END IF

    END IF

    ! Add dust conv dep to total dust dep flux
    IF ( l_tddf ) THEN

      ! If large scale precipitation is called, this array should
      ! already be allocated.
      IF (.NOT. ALLOCATED(tot_dust_dep_flux)) THEN
        ALLOCATE(tot_dust_dep_flux(row_length,rows))
        tot_dust_dep_flux(:,:) = 0.0
      END IF

      DO k=1, ndiv
        DO j=1,rows
          DO i=1,row_length
            tot_dust_dep_flux(i,j) = tot_dust_dep_flux(i,j) +              &
                                     conscav_dust(i,j,k)
          END DO
        END DO
      END DO !ndiv

    END IF

  END IF ! l_twobin_dust
END IF !L_DUST

! ----------------------------------------------------------------------
! Section Sulphur Cycle scavenging diagnostics
! ----------------------------------------------------------------------

item = 237                      !wet scav flux N in NH3
IF ( sf(item,sect) ) THEN

  !  Convert scavenged nh3 per timestep to flux per sec

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conwash_nh3,timestep)
  DO j=1,rows
    DO i=1,row_length
      conwash_nh3(i,j) = conwash_nh3(i,j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),conwash_nh3,    &
        row_length,rows,0,0,0,0,at_extremity,                     &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 237)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

item = 238                      !wet scav flux SO2
IF (icode <= 0 .AND. sf(item,sect)) THEN

  !  Convert scavenged so2 per timestep to flux per sec

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conwash_so2,timestep)
  DO j=1,rows
    DO i=1,row_length
      conwash_so2(i,j) = conwash_so2(i,j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),conwash_so2,    &
        row_length,rows,0,0,0,0,at_extremity,                     &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 238)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

item = 239                      !wet scav flux SO4_AIT
IF (icode <= 0 .AND. sf(item,sect)) THEN

  !  Convert scavenged so4_ait per timestep to flux per sec

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_so4ait,timestep)
  DO j=1,rows
    DO i=1,row_length
      conscav_so4ait(i,j) = conscav_so4ait(i,j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),conscav_so4ait, &
        row_length,rows,0,0,0,0,at_extremity,                     &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 239)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

item = 240                      !wet scav flux SO4_ACC
IF (icode <= 0 .AND. sf(item,sect)) THEN

  !  Convert scavenged so4_acc per timestep to flux per sec
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_so4acc,timestep)
  DO j=1,rows
    DO i=1,row_length
      conscav_so4acc(i,j) = conscav_so4acc(i,j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),conscav_so4acc, &
        row_length,rows,0,0,0,0,at_extremity,                     &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 240)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

item = 241                      !wet scav flux SO4_DIS
IF (icode <= 0 .AND. sf(item,sect)) THEN

  !  Convert scavenged so4_dis per timestep to flux per sec

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_so4dis,timestep)
  DO j=1,rows
    DO i=1,row_length
      conscav_so4dis(i,j) = conscav_so4dis(i,j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),conscav_so4dis, &
        row_length,rows,0,0,0,0,at_extremity,                     &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 241)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF
  !
END IF

! ----------------------------------------------------------------------
! Section Soot scheme scavenging diagnostics
! ----------------------------------------------------------------------

item = 242                      !wet scav flux soot
IF (icode <= 0 .AND. sf(item,sect)) THEN

  !  Convert scavenged soot per timestep to flux per sec

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_agedsoot,timestep)
  DO j=1,rows
    DO i=1,row_length
      conscav_agedsoot(i,j) = conscav_agedsoot(i,j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),                &
        conscav_agedsoot,                                         &
        row_length,rows,0,0,0,0,at_extremity,                     &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 242)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section Biomass aerosol scheme scavenging diagnostics
! ----------------------------------------------------------------------

item = 243                      !wet scav flux biomass
IF (icode <= 0 .AND. sf(item,sect)) THEN

  !  Convert scavenged biomass per timestep to flux per sec

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_agedbmass,timestep)
  DO j=1,rows
    DO i=1,row_length
      conscav_agedbmass(i,j) = conscav_agedbmass(i,j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),                &
        conscav_agedbmass,                                        &
        row_length,rows,0,0,0,0,at_extremity,                     &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 243)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

! ----------------------------------------------------------------------
! Section Fossil-fuel OC aerosol scheme scavenging diagnostics
! ----------------------------------------------------------------------

item = 244                     !wet scav flux ocff
IF (icode <= 0 .AND. sf(item,sect)) THEN

  !  Convert scavenged OCFF per timestep to flux per sec
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_agedocff,timestep)
  DO j=1,rows
    DO i=1,row_length
      conscav_agedocff(i,j) = conscav_agedocff(i,j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),                &
        conscav_agedocff,                                         &
        row_length,rows,0,0,0,0,at_extremity,                     &
        atmos_im,sect,item,                                       &
        icode,cmessage)

  IF (icode > 0) THEN
    cmessage=": error in copydiag(item 244)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

!
! ----------------------------------------------------------------------
! Section Ammonium nitrate scheme scavenging diagnostics
! ----------------------------------------------------------------------
!
item = 247                      !wet scav flux accumulation nitrate
IF (icode <= 0 .AND. sf(item,sect)) THEN
  !
  !  Convert scavenged accumulation nitrate per timestep to flux per sec

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_nitracc,timestep)
  DO j=1,rows
    DO i=1,row_length
      conscav_nitracc(i,j) = conscav_nitracc(i,j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),                &
        conscav_nitracc,                                          &
        row_length,rows,0,0,0,0,at_extremity,                     &
        atmos_im,sect,item,                                       &
        icode,cmessage)
  !
  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 247)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF
  !
END IF
!
item = 248                      !wet scav flux dissolved nitrate
IF (icode <= 0 .AND. sf(item,sect)) THEN
  !
  !  Convert scavenged dissolved nitrate per timestep to flux per sec

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(rows,row_length,conscav_nitrdiss,timestep)
  DO j=1,rows
    DO i=1,row_length
      conscav_nitrdiss(i,j) = conscav_nitrdiss(i,j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO
  !
  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),                &
        conscav_nitrdiss,                                         &
        row_length,rows,0,0,0,0,at_extremity,                     &
        atmos_im,sect,item,                                       &
        icode,cmessage)
  !
  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 248)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF
  !
END IF
!
!
!======================================================================
! New TCS (6A) warm diagnostics
!======================================================================
!----------------------------------------------------------------------
!  conv_type indicators
! ----------------------------------------------------------------------

! Item 400 conv_type all
item=400
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag ( stashwork(si(item,sect,im_index)), REAL(conv_type) &
     , row_length, rows, 0, 0, 0, 0, at_extremity                    &
     , atmos_im, sect, item, icode, cmessage )

  IF (icode >  0) THEN
    cmessage="diagnostics_conv: error in copydiag(5400)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF


! Item 401-404 conv_type 1-4
DO itmp=1,4
  item = 400+itmp
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    WHERE (conv_type == itmp)
      tmp=1.0
    ELSEWHERE
      tmp=0.0
    END WHERE

    ! DEPENDS ON: copydiag
    CALL copydiag ( stashwork(si(item,sect,im_index)), tmp          &
       , row_length, rows, 0, 0, 0, 0, at_extremity                 &
       , atmos_im, sect, item, icode, cmessage )

    IF (icode >  0) THEN
      cmessage="diagnostics_conv: error in copydiag(5401-404)"
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF
END DO


! Item 405-408 precip from conv_type 1-4
DO itmp=1,4
  item=404+itmp
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    WHERE (conv_type == itmp)
      tmp=precip_shall
    ELSEWHERE
      tmp=0.0
    END WHERE

    ! DEPENDS ON: copydiag
    CALL copydiag ( stashwork(si(item,sect,im_index)), tmp          &
       , row_length, rows, 0, 0, 0, 0, at_extremity                 &
       , atmos_im, sect, item, icode, cmessage )
    IF (icode >  0) THEN
      cmessage="diagnostics_conv: error in copydiag(5405-408)"
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF
END DO


! Item 409-412 dt from conv_type 1-4
DO itmp=1,4
  item=408+itmp
  IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,conv_type,itmp,tmp3d,dt_shall)
    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          IF (conv_type(i,j) == itmp) THEN
            tmp3d(i,j,k)=dt_shall(i,j,k)
          ELSE
            tmp3d(i,j,k)=0.0
          END IF
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d ( stashwork(si(item,sect,im_index)), tmp3d     &
       , row_length, rows, model_levels, 0, 0, 0, 0, at_extremity   &
       , stlist(1,stindex(1,item,sect,im_index)), len_stlist        &
       , stash_levels, num_stash_levels+1                           &
       , atmos_im, sect, item, icode, cmessage )

    IF (icode >  0) THEN
      cmessage="diagnostics_conv: error in copydiag(5409-412)"
      CALL ereport(routinename,icode,cmessage)
    END IF
  END IF
END DO


! Item 413-416 dq from conv_type 1-4
DO itmp=1,4
  item=412+itmp
  IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,conv_type,itmp,tmp3d,dq_shall)
    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          IF (conv_type(i,j) == itmp) THEN
            tmp3d(i,j,k)=dq_shall(i,j,k)
          ELSE
            tmp3d(i,j,k)=0.0
          END IF
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d ( stashwork(si(item,sect,im_index)), tmp3d     &
       , row_length, rows, model_levels, 0, 0, 0, 0, at_extremity   &
       , stlist(1,stindex(1,item,sect,im_index)), len_stlist        &
       , stash_levels, num_stash_levels+1                           &
       , atmos_im, sect, item, icode, cmessage )
    IF (icode >  0) THEN
      cmessage="diagnostics_conv: error in copydiag_3d(5413-416)"
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF
END DO


! Item 417-420 mass flux from conv_type 1-4
DO itmp=1,4
  item=416+itmp
  IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,conv_type,itmp,tmp3d,mf_shall)
    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          IF (conv_type(i,j) == itmp) THEN
            tmp3d(i,j,k)=mf_shall(i,j,k)
          ELSE
            tmp3d(i,j,k)=0.0
          END IF
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d ( stashwork(si(item,sect,im_index)), tmp3d     &
       , row_length, rows, model_levels, 0, 0, 0, 0, at_extremity   &
       , stlist(1,stindex(1,item,sect,im_index)), len_stlist        &
       , stash_levels, num_stash_levels+1                           &
       , atmos_im, sect, item, icode, cmessage )
    IF (icode >  0) THEN
      cmessage="diagnostics_conv: error in copydiag_3d(5417-420)"
      CALL ereport(routinename,icode,cmessage)
    END IF
  END IF
END DO

! Item 421-424 height from conv_type 1-4
DO itmp=1,4
  item=420+itmp
  IF (icode <= 0 .AND. sf(item,sect)) THEN

    WHERE (conv_type == itmp)
      tmp=ntpar_diag
    ELSEWHERE
      tmp=0.0
    END WHERE

    ! DEPENDS ON: copydiag
    CALL copydiag ( stashwork(si(item,sect,im_index)), tmp          &
       , row_length, rows, 0, 0, 0, 0, at_extremity                 &
       , atmos_im, sect,item, icode, cmessage )
    IF (icode >  0) THEN
      cmessage="diagnostics_conv: error in copydiag(5421-424)"
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF
END DO


! Item 425-428 wthetal from conv_type 1-4
DO itmp=1,4
  item=424+itmp
  IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,conv_type,itmp,tmp3d,               &
!$OMP        wthetal_flux_sh)
    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          IF (conv_type(i,j) == itmp) THEN
            tmp3d(i,j,k)=wthetal_flux_sh(i,j,k)
          ELSE
            tmp3d(i,j,k)=0.0
          END IF
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d ( stashwork(si(item,sect,im_index)), tmp3d     &
       , row_length, rows, model_levels, 0, 0, 0, 0, at_extremity   &
       , stlist(1,stindex(1,item,sect,im_index)), len_stlist        &
       , stash_levels, num_stash_levels+1                           &
       , atmos_im, sect, item, icode, cmessage )

    IF (icode >  0) THEN
      cmessage="diagnostics_conv: error in copydiag_3d(5425-428)"
      CALL ereport(routinename,icode,cmessage)
    END IF
  END IF
END DO


! Item 429-432 wqt from conv_type 1-4
DO itmp=1,4
  item=428+itmp
  IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,conv_type,itmp,tmp3d,wqt_flux_sh)
    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          IF (conv_type(i,j) == itmp) THEN
            tmp3d(i,j,k)=wqt_flux_sh(i,j,k)
          ELSE
            tmp3d(i,j,k)=0.0
          END IF
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d ( stashwork(si(item,sect,im_index)), tmp3d     &
       , row_length, rows, model_levels, 0, 0, 0, 0, at_extremity   &
       , stlist(1,stindex(1,item,sect,im_index)), len_stlist        &
       , stash_levels, num_stash_levels+1                           &
       , atmos_im, sect, item, icode, cmessage )

    IF (icode >  0) THEN
      cmessage="diagnostics_conv: error in copydiag_3d(5429-432)"
      CALL ereport(routinename,icode,cmessage)
    END IF
  END IF
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_conv
END MODULE diagnostics_conv_mod

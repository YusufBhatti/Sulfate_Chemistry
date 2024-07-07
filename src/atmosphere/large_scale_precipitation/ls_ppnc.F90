! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     SUBROUTINES LS_PPN and LS_PPNC------------------------------------
!    Purpose:
!            LS_PPN and LS_PPNC:
!             Calculate large-scale (dynamical) precipitation.
!             LS_PPNC is the gather/scatter routine which then
!             calls LSP_ICE.
!    Note: in all cases, level counters (incl subscripts) run from 1
!          (lowest model layer) to LEVELS (topmost  model layer)
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Documentation: UM Documentation Paper 26.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation
MODULE ls_ppnc_mod

  ! General atmosphere modules
USE planet_constants_mod, ONLY: g

USE murk_inputs_mod, ONLY: l_murk

! Microphysics modules
USE mphys_diags_mod,       ONLY: l_aggfr_diag,                                &
                                 l_psdep_diag, l_psaut_diag,                  &
                                 l_psacw_diag, l_psacr_diag,                  &
                                 l_psaci_diag, l_psmlt_diag,                  &
                                 l_psmltevp_diag, l_praut_diag,               &
                                 l_pracw_diag, l_prevp_diag,                  &
                                 l_pgaut_diag, l_pgacw_diag,                  &
                                 l_pgacs_diag, l_pgmlt_diag,                  &
                                 l_pifrw_diag, l_piprm_diag,                  &
                                 l_pidep_diag, l_piacw_diag,                  &
                                 l_piacr_diag, l_pimlt_diag,                  &
                                 l_pimltevp_diag, l_pifall_diag,              &
                                 l_psfall_diag, l_prfall_diag,                &
                                 l_pgfall_diag, l_plset_diag,                 &
                                 l_plevpset_diag, l_pifrr_diag,               &
                                 l_sfwater_diag, l_sfrain_diag,               &
                                 l_sfsnow_diag,                               &
                                 l_vtbranch_diag, l_vm_cry_diag,              &
                                 l_vm_agg_diag,                               &
                                 l_vm_used_diag,                              &
                                 psdep, psaut, frac_agg,                      &
                                 psacw, psacr, psaci, psmlt,                  &
                                 psmltevp, praut, pracw, prevp,               &
                                 pgaut, pgacw, pgacs, pgmlt, pifrw,           &
                                 piprm, pidep, piacw, piacr, pimlt,           &
                                 pimltevp, pifall, psfall, prfall,            &
                                 pgfall, plset, plevpset, pifrr,              &
                                 vm_cry, vm_agg, vtbranch_flag,               &
                                 vm_used, dbz_tot, dbz_g, dbz_i,              &
                                 dbz_i2, dbz_l, dbz_r,                        &
                                 sfwater, sfrain, sfsnow

USE mphys_inputs_mod,      ONLY: l_mcr_qcf2, l_mcr_qgraup,                    &
                                 l_mcr_qrain, l_fsd_generator

USE tuning_segments_mod,      ONLY: precip_segment_size

USE mphys_bypass_mod,      ONLY: l_ref_diag

USE cloud_inputs_mod,      ONLY:  i_cld_vn
USE pc2_constants_mod,     ONLY:  i_cld_pc2

! Grid bounds module

USE atm_fields_bounds_mod, ONLY: tdims, pdims

USE trignometric_mod,      ONLY: cos_theta_latitude

USE level_heights_mod,     ONLY: r_theta_levels

! Dr Hook modules
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

! Large scale precipitation modules
USE lsp_ice_mod,           ONLY: lsp_ice
USE lsp_scav_mod,          ONLY: lsp_scav

! FSD parameters module
USE fsd_parameters_mod,    ONLY: f_arr

! Use in KIND for large scale precip, used for compressed variables passed down
! from here
USE um_types,              ONLY: real_lsprec, real64

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1
CHARACTER(LEN=*),   PARAMETER, PRIVATE :: ModuleName='LS_PPNC_MOD'

!-----------------------------------------------------------------------
!  Enumerators with the indices for the super arrays
!
!  The variables in an enumerator have sequentially increasing integer
!  values. These variables represent the indices of the various data
!  in the super array used by a segment.
!
!  How to add a variable passed from ls_ppn to ls_ppnc.
! 
!  1) Add  the original field as argument into the ls_ppnc routine.
!
!  2) Add an extra element in the enumerator structure before the
!  last_index_??? element and after the element set as "= 1". The
!  enumerator that needs to be changed depends on the type of the
!  variable. For example, to add the logical variable "var", the
!  "var_i" element need to be added after "bland_i" into the logical
!  enumerator. The size of the logical super-array will increase
!  automatically since the value of last_index_??? is the total number
!  of elements in the enumerator.
!
!  3) Update the gather/scatter routines to take as argument the
!  original field and the new compressed array. Update the lsp
!  subroutines to use the new compressed sub-array. In the previous
!  example, the new var compressed array located in the logical
!  super-array will be passed to the sub-routines using the code
!  "spl(:,var_i)".
!
!  At the moment there are 2 super arrays: spr for REAL values and spl
!  for logical values. If a new type of compressed arrays is needed,
!  for example INTEGERS, then a new super array "spi" and a new
!  integer enumerator should be crated.
!  
!-----------------------------------------------------------------------

! Logical super-array enumerator
ENUM, BIND(C) 
ENUMERATOR ::                                                              &

  bland_i = 1,                                                             &
              ! gathered land/sea mask
  last_index_log
              ! Used to find the total number of variables in the
              ! supper array, i.e. this value minus 1.
END ENUM

! Real super-array enumerator
ENUM, BIND(C)
ENUMERATOR ::                                                              &
 cf_i = 1,                                                                 &
                      ! gathered Cloud fraction.
 q_i,                                                                      &
                       ! gathered Specific humidity (kg water/kg air).
 qcf_i,                                                                    &
                       ! gathered Cloud ice (kg per kg air).
 qcl_i,                                                                    &
                       ! gathered Cloud liquid water (kg per kg air).
 qcf2_i,                                                                   &
                        ! gathered cloud ice2 (kg per kg air).
 qrain_i,                                                                  &
                         ! gathered rain (kg per kg air).
 qgraup_i,                                                                 &
                          ! gathered graupel (kg per kg air).
 t_i,                                                                      &
                       ! gathered Temperature (K).
 uk_i,                                                                     &
                       ! gathered u wind on level k
 vk_i,                                                                     &
                       ! gathered v wind on level k
 ukp1_i,                                                                   &
                       ! gathered u wind on level k+1
 vkp1_i,                                                                   &
                       ! gathered v wind on level k+1
 r_theta_levels_i,                                                         &
                       ! gathered distance from centre of Earth
 r_theta_surf_i,                                                           &
                       ! ...and at surface
 g_cos_theta_latitude_i,                                                   &
                       ! gathered cos of latitude
 aero_i,                                                                   &
                       ! gathered Aerosol.
 lsrain_i,                                                                 &
                 !gathered Surface rainfall rate (kg per sq m per s).
 lssnow_i,                                                                 &
                 !gathered Surface snowfall rate (kg per sq m per s).
 lssnow2_i,                                                                &
                  !gathered layer snowfall rate (kg per sq m per s).
 lsgraup_i,                                                                &
                  !gathered layer graupel fall rate (kg/sq m/s)
 droplet_flux_i,                                                           &
                       ! gathered water droplet flux / kg m-2 s-1
 cttemp_i,                                                                 &
                            !gathered ice cloud top temperature.
 rainfrac_i,                                                               &
                            !gathered rain fraction.
 rainfrac_impr_i,                                                          &
                            !gathered improved rain fraction
 frac_ice_above_i,                                                         &
                            !gathered fraction of ice in layer above
 frac_agg_i,                                                               &
                       ! gathered aggregate fraction
 cfl_i,                                                                    &
                       ! gathered Cloud liquid fraction.
 cff_i,                                                                    &
                       ! gathered Cloud ice fraction.
 vfall_i,                                                                  &
                       ! gathered fall velocity (m per s).
 vfall2_i,                                                                 &
                        ! gathered fall velocity for qcf2 (m per s).
 vfall_rain_i,                                                             &
                        ! gathered fall velocity for qcf2 (m per s).
 vfall_graup_i,                                                            &
                        ! gathered fall velocity for qcf2 (m per s).
 rhc_i,                                                                    &
                        ! gathered RH_crit value at points.
 n_drop_tpr_i,                                                             &
 n_drop_out_i,                                                             &
                        ! gathered droplet numbers
 land_fract_i,                                                             &
                        ! gathered land fraction
 hmteff_i,                                                                 &
                        ! gathered effective mountain height
 zb_i,                                                                     &
                        ! gathered blocked layer depth for seederfeeder
 f_arr1_i,                                                                 &
 f_arr2_i,                                                                 &
 f_arr3_i,                                                                 &

  rhodz_i,                                                                 &
                    ! WORK Used for air mass p.u.a. in successive
!                            layers.
  deltaz_i,                                                                &
                       ! Thickness of layer (m)
   rhodz_dry_i,                                                            &
  ! Dry air density * layer thickness (kg m-2)
  rhodz_moist_i,                                                           &
                       ! Moist air density * layer thickness (kg m-2)

   p_i,                                                                    &    
                  ! WORK Used for pressure at successive levels.



! Microphysical process rate diagnostics (compressed arrays)
  psdep_i,                                                                 &
                  ! Deposition of vapour to snow aggregates
  psaut_i,                                                                 &
                  ! Autoconversion of aggregates from crystals
  psacw_i,                                                                 &
                  ! Accretion of liq. water by snow aggregates
  psacr_i,                                                                 &
                  ! Collection of rain by snow aggregates
  psaci_i,                                                                 &
                  ! Collection of ice crystals by aggregates
  psmlt_i,                                                                 &
                  ! Melting of snow aggregates
  psmltevp_i,                                                              &
                  ! Evaporation of melting aggregates
  praut_i,                                                                 &
                  ! Autoconversion of cloud drops to rain
  pracw_i,                                                                 &
                  ! Accretion of liq. water by rain
  prevp_i,                                                                 &
                  ! Evaporation of rain
  pgaut_i,                                                                 &
                  ! Autoconversion of graupel from aggregates
  pgacw_i,                                                                 &
                  ! Accretion of liq. water by graupel
  pgacs_i,                                                                 &
                  ! Collection of snow aggregates by graupel
  pgmlt_i,                                                                 &
                  ! Melting of graupel
  pifrw_i,                                                                 &
                  ! Homogeneous freezing nucleation
  pifrr_i,                                                                 &
                  ! Homogeneous freezing of rain
  piprm_i,                                                                 &
                  ! Heterogeneous (primary) nucleation
  pidep_i,                                                                 &
                  ! Deposition of vapour to ice crystals
  piacw_i,                                                                 &
                  ! Accretion of liq. water by ice crystals
  piacr_i,                                                                 &
                  ! Collection of rain by ice crystals
  pimlt_i,                                                                 &
                  ! Melting of ice crystals
  pimltevp_i,                                                              &
                  ! Evaporation of melting ice crystals
  pifall_i,                                                                &
                  ! Sedimentation of ice crystals
  psfall_i,                                                                &
                  ! Sedimentation of aggregates
  prfall_i,                                                                &
                  ! Sedimentation of rain
  pgfall_i,                                                                &
                  ! Sedimentation of graupel
  plset_i,                                                                 &
                  ! Droplet settling of liquid water
  plevpset_i,                                                              &
                  ! Evaporated settled droplets
  vm_agg_i,                                                                &
               ! Mass-weighted fallspeed with aggregate parameters
  vm_cry_i,                                                                &
               ! Mass-weighted fallspeed with crystal parameters
  vtbranch_flag_i,                                                         &
               ! Flag indicating which vt-D is used
  vm_used_i,                                                               &
               ! Mass-weighted fallspeed used

  sfwater_i,                                                               &
               ! Seeder feeder orographic water mixing ratio
  sfrain_i,                                                                &
               ! Seeder feeder orographic rain production
  sfsnow_i,                                                                &
               ! Seeder feeder orographic snow production

  dbz_tot_i,                                                               &
               ! Total reflectivity (dBZ)
  dbz_g_i,                                                                 &
              ! Graupel reflectivity (dBZ) 
  dbz_i_i,                                                                 &
              ! Ice Agg. reflectivity (dBZ) 
  dbz_i2_i,                                                                &
              ! Ice Cry. reflectivity (dBZ) 
  dbz_l_i,                                                                 &
              ! Cloud liquid reflectivity (dBZ) 
  dbz_r_i,                                                                 &
              ! Rain reflectivity (dBZ)

  last_index_real
              ! Used to find the total number of variables in the
              ! supper array, i.e. this value minus 1.
END ENUM

! A LSP segment data is composed of:
! 1) segment_size * lsp_num_real  REAL    values 
! 2) segment_size * lsp_num_log   LOGICAL values 


! The total number of variables in the LSP super array.
INTEGER, PARAMETER :: lsp_num_real = last_index_real - 1
INTEGER, PARAMETER :: lsp_num_log  = last_index_log  - 1

CONTAINS

!-----------------------------------------------------------------------
!  Main routine
!-----------------------------------------------------------------------

SUBROUTINE ls_ppnc( level, ix, n,                                             &
 lsrain,lssnow,lssnow2,lsgraup, droplet_flux,                                 &
 cf,cfl,cff,                                                                  &
 qcf,qcl,t, qcf2,qrain,qgraup,                                                &
 n_drop_tpr, n_drop_out,                                                      &
 aerosol, land_fract,                                                         &
 hmteff, zb,                                                                  &
!--------------------------------------------------------
! Layer thicknesses and variables passed from layer above
!--------------------------------------------------------
 q, p_theta_levels, layer_thickness,                                          &
 deltaz, rhodz_dry, rhodz_moist,                                              &
 rhc_row_length, rhc_rows, bland, rhcrit,                                     &
 vfall, vfall2, vfall_rain, vfall_graup,                                      &
 frac_ice_above, cttemp, rainfrac, rainfrac_impr, niters_mp,                  &
 uk, vk, ukp1, vkp1, onelevel,                                                &
!--------------------------------------------------------
! COSP flag
!--------------------------------------------------------
 l_cosp_lsp                                                                   &
 )


IMPLICIT NONE

INTEGER ::                                                                    &
 level,                                                                       &
              ! level number
 n,                                                                           &
              ! IN Number of points where pptn non-zero from above
!                    or where CF>CFMIN
   niters_mp,                                                                 &
                ! Number of total iterations of microphysics
   rhc_row_length,rhc_rows,                                                   &
   ix (  tdims%i_len  *                                                       &
         tdims%j_len  , 2 )
                                ! IN gather/scatter index
REAL ::                                                                       &
 cf (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                   &
                            ! IN Cloud fraction.
 cfl(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                   &
                            ! IN Cloud liquid fraction.
 cff(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                   &
                            ! IN Cloud ice fraction.
! CF, CFL and CFF are IN/OUT if the PC2 cloud scheme is in use.

    p_theta_levels( tdims%i_start:tdims%i_end,                                &
                    tdims%j_start:tdims%j_end ),                              &
    layer_thickness(tdims%i_start:tdims%i_end,                                &
                    tdims%j_start:tdims%j_end ),                              &
                                          ! IN thickness of layer (Pa)
    deltaz( tdims%i_start:tdims%i_end,                                        &
            tdims%j_start:tdims%j_end ),                                      &
                                          ! IN thickness of layer (m)
    rhodz_dry(tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end ),                                    &
                                          ! Dry air density
                                          ! * layer thickness (kg m-2)
    rhodz_moist(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end ),                                  &
                                          ! Moist air density
                                          ! * layer thickness (kg m-2)

    rhcrit(rhc_row_length,rhc_rows)
!                       IN Critical humidity for cloud formation.


LOGICAL, INTENT(IN) ::                                                        &
     bland(tdims%i_start:tdims%i_end,                                         &
           tdims%j_start:tdims%j_end )
                        !IN Land/sea mask

LOGICAL, INTENT(IN) ::  l_cosp_lsp
                        ! Switch for COSP LS diagnostics

REAL, INTENT(INOUT) ::                                                        &
     q(tdims%i_start:tdims%i_end,                                             &
       tdims%j_start:tdims%j_end),                                            &
                          ! INOUT Specific humidity (kg water/kg air).
      qcf(tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end),                                         &
                          ! INOUT Cloud ice (kg per kg air).
      qcl(tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end),                                         &
                          ! INOUT Cloud liquid water (kg per kg air).
      qcf2(tdims%i_start:tdims%i_end,                                         &
           tdims%j_start:tdims%j_end),                                        &
                          ! INOUT Cloud ice2 (kg per kg air).
      qrain(tdims%i_start:tdims%i_end,                                        &
            tdims%j_start:tdims%j_end),                                       &
                          ! INOUT Rain water (kg per kg air).
      qgraup(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end),                                      &
                          ! INOUT Graupel water (kg per kg air).
      t(tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end),                                           &
                          ! INOUT Temperature (K).
      aerosol(tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end),                                     &
                          ! INOUT Aerosol (K).
      lsrain(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end),                                      &
                          !INOUT Surface rainfall rate (kg m^-2 s^-1).
      lssnow(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end),                                      &
                          !INOUT Surface snowfall rate (kg m^-2 s^-1).
      lssnow2(tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end),                                     &
                          !INOUT layer snowfall rate (kg m^-2 s^-1).
      lsgraup(tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end),                                     &
                          !INOUT layer graupelfall rate (kg m^-2 s^-1)
      droplet_flux(tdims%i_start:tdims%i_end,                                 &
                   tdims%j_start:tdims%j_end),                                &
                          !INOUT water droplet flux / kg m^-2 s^-1
      cttemp(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end),                                      &
                          ! INOUT Ice cloud top temperature (K)
      rainfrac(tdims%i_start:tdims%i_end,                                     &
               tdims%j_start:tdims%j_end),                                    &
                          ! INOUT Rain fraction.
      frac_ice_above(tdims%i_start:tdims%i_end,                               &
                     tdims%j_start:tdims%j_end),                              &
                          ! INOUT Ice fraction from layer above
 vfall(tdims%i_start:tdims%i_end,                                             &
       tdims%j_start:tdims%j_end),                                            &
                                   ! INOUT fall velocity of ice (m per
 vfall2(tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end),                                           &
                                    ! INOUT fall vel. of rain (m/s)
 vfall_rain(tdims%i_start:tdims%i_end,                                        &
            tdims%j_start:tdims%j_end),                                       &
                                    ! INOUT fall vel. of rain (m/s)
 vfall_graup(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end)
                                    ! INOUT fall vel. of graupel (m/s)

REAL, INTENT(IN) :: rainfrac_impr( tdims%i_start:tdims%i_end,                 &
                                   tdims%j_start:tdims%j_end )
                                    ! Improved rain fraction

REAL,INTENT(IN) ::                                                            &
! U and V wind at levels k and k+1 on P grid, so defined using pdims.
   uk  (pdims%i_start : pdims%i_end,                                          &
        pdims%j_start : pdims%j_end),                                         &
   ukp1(pdims%i_start : pdims%i_end,                                          &
        pdims%j_start : pdims%j_end),                                         &
   vk  (pdims%i_start : pdims%i_end,                                          &
        pdims%j_start : pdims%j_end),                                         &
   vkp1(pdims%i_start : pdims%i_end,                                          &
        pdims%j_start : pdims%j_end)
! Level we are interested in for r_theta_level
INTEGER, INTENT(IN) :: onelevel

REAL ::                                                                       &
     land_fract(tdims%i_start:tdims%i_end,                                    &
              tdims%j_start:tdims%j_end)

REAL,INTENT(IN) ::                                                            &
     hmteff(tdims%i_start:tdims%i_end,                                        &
            tdims%j_start:tdims%j_end)
REAL,INTENT(IN) ::                                                            &
     zb(tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end)

REAL ::                                                                       &
                    !, Intent(IN)
      n_drop_tpr( tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )
                             ! Tapered droplet number
REAL ::                                                                       &
                    ! INTENT(INOUT)
      n_drop_out( tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )
                    ! droplet number from autoconversion

!    Workspace usage ---------------------------------------------------

! Super arrays containing all the data used by a segment.
REAL (KIND=real_lsprec), ALLOCATABLE :: spr(:,:)
LOGICAL,                 ALLOCATABLE :: spl(:,:)



INTEGER :: jj,i1,i2,ip

REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LS_PPNC'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(jj, i1, i2, ip, spr, spl)

!$OMP DO SCHEDULE(DYNAMIC)
DO jj = 1, n, precip_segment_size

  i1 = jj
  i2 = MIN(jj+precip_segment_size-1, n)
  ip = i2-i1+1

!-----------------------------------------------------------------------
!  Allocate super arrays
!-----------------------------------------------------------------------

  ALLOCATE( spl (ip, lsp_num_log)    )
  ALLOCATE( spr (ip, lsp_num_real)   )
  
!-----------------------------------------------------------------------
!  Gather variables using index
!-----------------------------------------------------------------------

  CALL ls_ppnc_gather( level, ix, ip, i1,                                     &
    lsrain,lssnow,lssnow2,lsgraup, droplet_flux,                              &
    cf,cfl,cff,                                                               &
    qcf,qcl,t, qcf2,qrain,qgraup,                                             &
    n_drop_tpr, n_drop_out,                                                   &
    aerosol, land_fract,                                                      &
    hmteff, zb,                                                               &
    q, p_theta_levels, layer_thickness,                                       &
    deltaz, rhodz_dry, rhodz_moist,                                           &
    rhc_row_length, rhc_rows, bland, rhcrit,                                  &
    vfall, vfall2, vfall_rain, vfall_graup,                                   &
    frac_ice_above, cttemp, rainfrac, rainfrac_impr,                          &
    uk, vk, ukp1, vkp1, onelevel,                                             &
    l_cosp_lsp,                                                               &
    spr(:,p_i),              spr(:,rhodz_i),       spr(:,deltaz_i),           &
    spr(:,rhodz_dry_i),      spr(:,rhodz_moist_i), spr(:,rhc_i),              &
    spr(:,land_fract_i),     spr(:,hmteff_i),      spr(:,zb_i),               &
    spr(:,qcf_i),            spr(:,qcl_i),         spr(:,q_i),                &
    spr(:,qcf2_i),           spr(:,qrain_i),                                  &
    spr(:,qgraup_i),         spr(:,n_drop_tpr_i),                             &
    spr(:,n_drop_out_i),     spr(:,lsrain_i),      spr(:,vfall_rain_i),       &
    spr(:,lssnow_i),         spr(:,vfall_i),       spr(:,lssnow2_i),          &
    spr(:,vfall2_i),         spr(:,lsgraup_i),     spr(:,vfall_graup_i),      &
    spr(:,droplet_flux_i),   spr(:,frac_ice_above_i),                         &
    spr(:,frac_agg_i),       spr(:,cttemp_i),      spr(:,rainfrac_i),         &
    spr(:,rainfrac_impr_i),  spr(:,t_i),           spr(:,cf_i),               &
    spr(:,cfl_i),            spr(:,cff_i),         spl(:,bland_i),            &
    spr(:,psdep_i),          spr(:,psaut_i),                                  &
    spr(:,psacw_i),          spr(:,psacr_i),       spr(:,psaci_i),            &
    spr(:,psmlt_i),          spr(:,psmltevp_i),    spr(:,praut_i),            &
    spr(:,pracw_i),          spr(:,prevp_i),       spr(:,pgaut_i),            &
    spr(:,pgacw_i),          spr(:,pgacs_i),       spr(:,pgmlt_i),            &
    spr(:,pifrw_i),          spr(:,pifrr_i),       spr(:,piprm_i),            &
    spr(:,pidep_i),          spr(:,piacw_i),       spr(:,piacr_i),            &
    spr(:,pimlt_i),          spr(:,pimltevp_i),    spr(:,pifall_i),           &
    spr(:,psfall_i),         spr(:,prfall_i),      spr(:,pgfall_i),           &
    spr(:,plset_i),          spr(:,plevpset_i),    spr(:,dbz_tot_i),          &
    spr(:,dbz_g_i),          spr(:,dbz_i_i),       spr(:,dbz_i2_i),           &
    spr(:,dbz_l_i),          spr(:,dbz_r_i),       spr(:,sfwater_i),          &
    spr(:,sfrain_i),         spr(:,sfsnow_i),      spr(:,uk_i),               &
    spr(:,vk_i),             spr(:,ukp1_i),                                   & 
    spr(:,vkp1_i),           spr(:,r_theta_levels_i),                         &
    spr(:,g_cos_theta_latitude_i),                 spr(:,r_theta_surf_i),     &
    spr(:,f_arr1_i),         spr(:,f_arr2_i),      spr(:,f_arr3_i),           &
    spr(:,vm_cry_i),         spr(:,vm_agg_i),      spr(:,vtbranch_flag_i),    &
    spr(:,vm_used_i),        spr(:,aero_i)                                    &
    )

!-----------------------------------------------------------------------
! Call subroutine lsp_ice to control transfers between different
! microphysical species
!-----------------------------------------------------------------------

  CALL lsp_ice(                                                               &
    spr(:,p_i), spr(:,rhodz_i), spr(:,deltaz_i), spr(:,rhodz_dry_i),          &
    spr(:,rhodz_moist_i), ip, spr(:,rhc_i),                                   &
    spr(:,land_fract_i), spr(:,hmteff_i), spr(:,zb_i), spr(:,qcf_i),          &
    spr(:,qcl_i), spr(:,q_i), spr(:,qcf2_i), spr(:,qrain_i),                  &
    spr(:,qgraup_i), spr(:,n_drop_tpr_i), spr(:,n_drop_out_i),                &
    spr(:,lsrain_i), spr(:,vfall_rain_i), spr(:,lssnow_i),                    &
    spr(:,vfall_i), spr(:,lssnow2_i), spr(:,vfall2_i),                        &
    spr(:,lsgraup_i), spr(:,vfall_graup_i), spr(:,droplet_flux_i),            &
    spr(:,frac_ice_above_i), spr(:,frac_agg_i), spr(:,cttemp_i),              &
    spr(:,rainfrac_i), spr(:,rainfrac_impr_i), spr(:,t_i),                    &
    spr(:,cf_i), spr(:,cfl_i), spr(:,cff_i), spl(:,bland_i),                  &
    spr(:,psdep_i), spr(:,psaut_i), spr(:,psacw_i), spr(:,psacr_i),           &
    spr(:,psaci_i), spr(:,psmlt_i), spr(:,psmltevp_i),                        &
    spr(:,praut_i), spr(:,pracw_i), spr(:,prevp_i), spr(:,pgaut_i),           &
    spr(:,pgacw_i), spr(:,pgacs_i), spr(:,pgmlt_i), spr(:,pifrw_i),           &
    spr(:,pifrr_i), spr(:,piprm_i), spr(:,pidep_i), spr(:,piacw_i),           &
    spr(:,piacr_i), spr(:,pimlt_i), spr(:,pimltevp_i),                        &
    spr(:,pifall_i), spr(:,psfall_i), spr(:,prfall_i),                        &
    spr(:,pgfall_i), spr(:,plset_i), spr(:,plevpset_i),                       &
    spr(:,dbz_tot_i), spr(:,dbz_g_i), spr(:,dbz_i_i),                         &
    spr(:,dbz_i2_i), spr(:,dbz_l_i), spr(:,dbz_r_i),                          & 
    spr(:,sfwater_i), spr(:,sfrain_i), spr(:,sfsnow_i), niters_mp,            &
    spr(:,uk_i), spr(:,vk_i), spr(:,ukp1_i), spr(:,vkp1_i),                   &
    spr(:,r_theta_levels_i), spr(:,g_cos_theta_latitude_i),                   &
    spr(:,r_theta_surf_i), spr(:,f_arr1_i), spr(:,f_arr2_i),                  &
    spr(:,f_arr3_i), spr(:,vm_cry_i), spr(:,vm_agg_i),                        &
    spr(:,vtbranch_flag_i), spr(:,vm_used_i)                                  &
    )

!-----------------------------------------------------------------------
! Lose Murk aerosol by scavenging: call LSP_SCAV
!-----------------------------------------------------------------------

  IF (l_murk) THEN
    CALL lsp_scav(ip,spr(:,lsrain_i), spr(:,lssnow_i), spr(:,droplet_flux_i), &
                  spr(:,aero_i))
  END IF

!-----------------------------------------------------------------------
! Scatter back arrays which will have been changed.
!-----------------------------------------------------------------------

  CALL ls_ppnc_scatter( level, ix, ip, i1,                                    &
    lsrain,lssnow,lssnow2,lsgraup, droplet_flux,                              &
    cf,cfl,cff,                                                               &
    qcf,qcl,t, qcf2,qrain,qgraup,                                             &
    n_drop_tpr, n_drop_out,                                                   &
    aerosol, land_fract,                                                      &
    hmteff, zb,                                                               &
    q, p_theta_levels, layer_thickness,                                       &
    deltaz, rhodz_dry, rhodz_moist,                                           &
    rhc_row_length, rhc_rows, bland, rhcrit,                                  &
    vfall, vfall2, vfall_rain, vfall_graup,                                   &
    frac_ice_above, cttemp, rainfrac, rainfrac_impr,                          &
    uk, vk, ukp1, vkp1, onelevel,                                             &
    l_cosp_lsp,                                                               &
    spr(:,p_i),              spr(:,rhodz_i),       spr(:,deltaz_i),           &
    spr(:,rhodz_dry_i),      spr(:,rhodz_moist_i), spr(:,rhc_i),              &
    spr(:,land_fract_i),     spr(:,hmteff_i),      spr(:,zb_i),               &
    spr(:,qcf_i),            spr(:,qcl_i),         spr(:,q_i),                &
    spr(:,qcf2_i),           spr(:,qrain_i),                                  &
    spr(:,qgraup_i),         spr(:,n_drop_tpr_i),                             &
    spr(:,n_drop_out_i),     spr(:,lsrain_i),      spr(:,vfall_rain_i),       &
    spr(:,lssnow_i),         spr(:,vfall_i),       spr(:,lssnow2_i),          &
    spr(:,vfall2_i),         spr(:,lsgraup_i),     spr(:,vfall_graup_i),      &
    spr(:,droplet_flux_i),   spr(:,frac_ice_above_i),                         &
    spr(:,frac_agg_i),       spr(:,cttemp_i),      spr(:,rainfrac_i),         &
    spr(:,rainfrac_impr_i),  spr(:,t_i),           spr(:,cf_i),               &
    spr(:,cfl_i),            spr(:,cff_i),         spl(:,bland_i),            &
    spr(:,psdep_i),          spr(:,psaut_i),                                  &
    spr(:,psacw_i),          spr(:,psacr_i),       spr(:,psaci_i),            &
    spr(:,psmlt_i),          spr(:,psmltevp_i),    spr(:,praut_i),            &
    spr(:,pracw_i),          spr(:,prevp_i),       spr(:,pgaut_i),            &
    spr(:,pgacw_i),          spr(:,pgacs_i),       spr(:,pgmlt_i),            &
    spr(:,pifrw_i),          spr(:,pifrr_i),       spr(:,piprm_i),            &
    spr(:,pidep_i),          spr(:,piacw_i),       spr(:,piacr_i),            &
    spr(:,pimlt_i),          spr(:,pimltevp_i),    spr(:,pifall_i),           &
    spr(:,psfall_i),         spr(:,prfall_i),      spr(:,pgfall_i),           &
    spr(:,plset_i),          spr(:,plevpset_i),    spr(:,dbz_tot_i),          &
    spr(:,dbz_g_i),          spr(:,dbz_i_i),       spr(:,dbz_i2_i),           &
    spr(:,dbz_l_i),          spr(:,dbz_r_i),       spr(:,sfwater_i),          &
    spr(:,sfrain_i),         spr(:,sfsnow_i),      spr(:,uk_i),               &
    spr(:,vk_i),             spr(:,ukp1_i),                                   & 
    spr(:,vkp1_i),           spr(:,r_theta_levels_i),                         &
    spr(:,g_cos_theta_latitude_i),                 spr(:,r_theta_surf_i),     &
    spr(:,f_arr1_i),         spr(:,f_arr2_i),      spr(:,f_arr3_i),           &
    spr(:,vm_cry_i),         spr(:,vm_agg_i),      spr(:,vtbranch_flag_i),    &
    spr(:,vm_used_i),        spr(:,aero_i)                                    &
    )
    

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------
  
  DEALLOCATE( spr )
  DEALLOCATE( spl )

END DO ! Loop over points
!$OMP END DO 

!$OMP END PARALLEL 


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ls_ppnc


!-----------------------------------------------------------------------
!  Gather variables, i.e initialise the data
!-----------------------------------------------------------------------

SUBROUTINE ls_ppnc_gather( level, ix, ip, i1,                                 &
  lsrain,lssnow,lssnow2,lsgraup, droplet_flux,                                &
  cf,cfl,cff,                                                                 &
  qcf,qcl,t, qcf2,qrain,qgraup,                                               &
  n_drop_tpr, n_drop_out,                                                     &
  aerosol, land_fract,                                                        &
  hmteff, zb,                                                                 &
!--------------------------------------------------------
! Layer thicknesses and variables passed from layer above
!--------------------------------------------------------
  q, p_theta_levels, layer_thickness,                                         &
  deltaz, rhodz_dry, rhodz_moist,                                             &
  rhc_row_length, rhc_rows, bland, rhcrit,                                    &
  vfall, vfall2, vfall_rain, vfall_graup,                                     &
  frac_ice_above, cttemp, rainfrac, rainfrac_impr,                            &
  uk, vk, ukp1, vkp1, onelevel,                                               &
!--------------------------------------------------------
! COSP flag
!--------------------------------------------------------
  l_cosp_lsp,                                                                 &
!--------------------------------------------------------
! Compressed arrays
!--------------------------------------------------------
  p_c, rhodz_c, deltaz_c, rhodz_dry_c,                                        &
  rhodz_moist_c, rhc_c, land_fract_c,                                         &
  hmteff_c, zb_c, qcf_c, qcl_c, q_c, qcf2_c, qrain_c,                         &
  qgraup_c, n_drop_tpr_c, n_drop_out_c,                                       &
  lsrain_c,vfall_rain_c, lssnow_c,vfall_c,                                    &
  lssnow2_c,vfall2_c, lsgraup_c,                                              &
  vfall_graup_c, droplet_flux_c,                                              &
  frac_ice_above_c, frac_agg_c,                                               &
  cttemp_c, rainfrac_c, rainfrac_impr_c,                                      &
  t_c,cf_c,cfl_c, cff_c, bland_c,                                             &
  psdep_c,  psaut_c,  psacw_c, psacr_c,                                       &
  psaci_c,  psmlt_c,  psmltevp_c,                                             &
  praut_c,  pracw_c,  prevp_c,                                                &
  pgaut_c,  pgacw_c,  pgacs_c,  pgmlt_c,                                      &
  pifrw_c,  pifrr_c,  piprm_c,  pidep_c,                                      &
  piacw_c,  piacr_c,  pimlt_c,  pimltevp_c,                                   &
  pifall_c, psfall_c, prfall_c, pgfall_c,                                     &
  plset_c,  plevpset_c, dbz_tot_c, dbz_g_c,                                   &
  dbz_i_c,  dbz_i2_c, dbz_l_c,  dbz_r_c,                                      &
  sfwater_c, sfrain_c, sfsnow_c,                                              &
  uk_c, vk_c, ukp1_c, vkp1_c,                                                 &
  r_theta_levels_c, g_cos_theta_latitude_c,                                   &
  r_theta_surf_c,                                                             &
  f_arr1_c, f_arr2_c, f_arr3_c,                                               &
  vm_cry_c, vm_agg_c, vtbranch_flag_c, vm_used_c, aero_c                      &
 )


IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                        &
 level,                                                                       &
              ! level number
 ip,                                                                          &
              ! Number of points
 i1,                                                                          &
              ! Starting point
 rhc_row_length,                                                              &
 rhc_rows,                                                                    &
 ix (  tdims%i_len * tdims%j_len , 2 )
              ! Gather/scatter index

REAL, INTENT(IN) ::                                                           &
  cf (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                  &
              ! Cloud fraction.
  cfl(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                  &
              ! Cloud liquid fraction.
  cff(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                  &
              ! Cloud ice fraction.
  p_theta_levels( tdims%i_start:tdims%i_end,                                  &
                  tdims%j_start:tdims%j_end ),                                &
              !    
  layer_thickness(tdims%i_start:tdims%i_end,                                  &
                  tdims%j_start:tdims%j_end ),                                &
              ! Thickness of layer (Pa)
  deltaz( tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end ),                                        &
              ! Thickness of layer (m)
  rhodz_dry(tdims%i_start:tdims%i_end,                                        &
            tdims%j_start:tdims%j_end ),                                      &
              ! Dry air density
              ! * layer thickness (kg m-2)
  rhodz_moist(tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end ),                                    &
              ! Moist air density
              ! * layer thickness (kg m-2)
  rhcrit(rhc_row_length,rhc_rows)
              ! Critical humidity for cloud formation.

LOGICAL, INTENT(IN) ::                                                        &
  bland(tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end )
              ! Land/sea mask

LOGICAL, INTENT(IN) ::  l_cosp_lsp
              ! Switch for COSP LS diagnostics

REAL, INTENT(IN) ::                                                           &
  q(tdims%i_start:tdims%i_end,                                                &
    tdims%j_start:tdims%j_end),                                               &
              ! Specific humidity (kg water/kg air).
  qcf(tdims%i_start:tdims%i_end,                                              &
      tdims%j_start:tdims%j_end),                                             &
              ! Cloud ice (kg per kg air).
  qcl(tdims%i_start:tdims%i_end,                                              &
      tdims%j_start:tdims%j_end),                                             &
              ! Cloud liquid water (kg per kg air).
  qcf2(tdims%i_start:tdims%i_end,                                             &
       tdims%j_start:tdims%j_end),                                            &
              ! Cloud ice2 (kg per kg air).
  qrain(tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end),                                           &
              ! Rain water (kg per kg air).
  qgraup(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end),                                          &
              ! Graupel water (kg per kg air).
  t(tdims%i_start:tdims%i_end,                                                &
    tdims%j_start:tdims%j_end),                                               &
              ! Temperature (K).
  aerosol(tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end),                                         &
              ! Aerosol (K).
  lsrain(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end),                                          &
              ! Surface rainfall rate (kg m^-2 s^-1).
  lssnow(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end),                                          &
              ! Surface snowfall rate (kg m^-2 s^-1).
  lssnow2(tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end),                                         &
              ! Layer snowfall rate (kg m^-2 s^-1).
  lsgraup(tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end),                                         &
              ! Layer graupelfall rate (kg m^-2 s^-1)
  droplet_flux(tdims%i_start:tdims%i_end,                                     &
               tdims%j_start:tdims%j_end),                                    &
              ! Water droplet flux / kg m^-2 s^-1
  cttemp(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end),                                          &
              ! Ice cloud top temperature (K)
  rainfrac(tdims%i_start:tdims%i_end,                                         &
           tdims%j_start:tdims%j_end),                                        &
              ! Rain fraction.
  frac_ice_above(tdims%i_start:tdims%i_end,                                   &
                 tdims%j_start:tdims%j_end),                                  &
              ! Ice fraction from layer above
  vfall(tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end),                                           &
              ! Fall velocity of ice (m per
  vfall2(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end),                                          &
              ! Fall vel. of rain (m/s)
  vfall_rain(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end),                                      &
              ! Fall vel. of rain (m/s)
  vfall_graup(tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end)
              ! Fall vel. of graupel (m/s)

REAL, INTENT(IN) ::                                                           &
  rainfrac_impr(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end)
              ! Improved rain fraction

REAL,INTENT(IN) ::                                                            &
  ! U and V wind at levels k and k+1 on P grid, so defined using pdims.
  uk  (pdims%i_start : pdims%i_end,                                           &
       pdims%j_start : pdims%j_end),                                          &
  ukp1(pdims%i_start : pdims%i_end,                                           &
       pdims%j_start : pdims%j_end),                                          &
  vk  (pdims%i_start : pdims%i_end,                                           &
       pdims%j_start : pdims%j_end),                                          &
  vkp1(pdims%i_start : pdims%i_end,                                           &
       pdims%j_start : pdims%j_end)

! Level we are interested in for r_theta_level
INTEGER, INTENT(IN) :: onelevel

REAL,INTENT(IN) ::                                                            &
  land_fract(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end)

REAL,INTENT(IN) ::                                                            &
  hmteff(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end)
REAL,INTENT(IN) ::                                                            &
  zb(tdims%i_start:tdims%i_end,                                               &
     tdims%j_start:tdims%j_end)

REAL,INTENT(IN) ::                                                            &
  n_drop_tpr(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end )
              ! Tapered droplet number
REAL, INTENT(IN) ::                                                           &
  n_drop_out( tdims%i_start : tdims%i_end,                                    &
                  tdims%j_start : tdims%j_end )
              ! Droplet number from autoconversion

! Compressed arrays ---------------

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  cf_c(ip),                                                                   &
  q_c(ip),                                                                    &
  qcf_c(ip),                                                                  &
  qcl_c(ip),                                                                  &
  qcf2_c(ip),                                                                 &
  qrain_c(ip),                                                                &
  qgraup_c(ip),                                                               &
  t_c(ip),                                                                    &
  uk_c(ip),                                                                   &
  vk_c(ip),                                                                   &
  ukp1_c(ip),                                                                 &
  vkp1_c(ip),                                                                 &
  r_theta_levels_c(ip),                                                       &
  r_theta_surf_c(ip),                                                         &
  g_cos_theta_latitude_c(ip),                                                 &
  aero_c(ip),                                                                 &
  lsrain_c(ip),                                                               &
  lssnow_c(ip),                                                               &
  lssnow2_c(ip),                                                              &
  lsgraup_c(ip),                                                              &
  droplet_flux_c(ip),                                                         &
  cttemp_c(ip),                                                               &
  rainfrac_c(ip),                                                             &
  rainfrac_impr_c(ip),                                                        &
  frac_ice_above_c(ip),                                                       &
  frac_agg_c(ip),                                                             &
  cfl_c(ip),                                                                  &
  cff_c(ip),                                                                  &
  vfall_c(ip),                                                                &
  vfall2_c(ip),                                                               &
  vfall_rain_c(ip),                                                           &
  vfall_graup_c(ip),                                                          &
  rhc_c(ip),                                                                  &
  n_drop_tpr_c(ip),                                                           &
  n_drop_out_c(ip),                                                           &
  land_fract_c(ip),                                                           &
  hmteff_c(ip),                                                               &
  zb_c(ip),                                                                   &
  f_arr1_c(ip),                                                               &
  f_arr2_c(ip),                                                               &
  f_arr3_c(ip)

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  rhodz_c(ip),                                                                &
  deltaz_c(ip),                                                               &
  rhodz_dry_c(ip),                                                            &
  rhodz_moist_c(ip),                                                          &
  p_c(ip)  

LOGICAL, INTENT(INOUT) :: bland_c(ip)          ! gathered land/sea mask

! Microphysical process rate diagnostics (compressed arrays)
REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  psdep_c(ip),                                                                &
  psaut_c(ip),                                                                &
  psacw_c(ip),                                                                &
  psacr_c(ip),                                                                &
  psaci_c(ip),                                                                &
  psmlt_c(ip),                                                                &
  psmltevp_c(ip),                                                             &
  praut_c(ip),                                                                &
  pracw_c(ip),                                                                &
  prevp_c(ip),                                                                &
  pgaut_c(ip),                                                                &
  pgacw_c(ip),                                                                &
  pgacs_c(ip),                                                                &
  pgmlt_c(ip),                                                                &
  pifrw_c(ip),                                                                &
  pifrr_c(ip),                                                                &
  piprm_c(ip),                                                                &
  pidep_c(ip),                                                                &
  piacw_c(ip),                                                                &
  piacr_c(ip),                                                                &
  pimlt_c(ip),                                                                &
  pimltevp_c(ip),                                                             &
  pifall_c(ip),                                                               &
  psfall_c(ip),                                                               &
  prfall_c(ip),                                                               &
  pgfall_c(ip),                                                               &
  plset_c(ip),                                                                &
  plevpset_c(ip),                                                             &
  vm_agg_c(ip),                                                               &
  vm_cry_c(ip),                                                               &
  vtbranch_flag_c(ip),                                                        &
  vm_used_c(ip)

REAL (KIND=real_lsprec), INTENT(INOUT)  ::                                    &
  sfwater_c(ip),                                                              &
  sfrain_c(ip),                                                               &
  sfsnow_c(ip)

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  dbz_tot_c(ip),                                                              &
  dbz_g_c(ip),                                                                &
  dbz_i_c(ip),                                                                &
  dbz_i2_c(ip),                                                               &
  dbz_l_c(ip),                                                                &
  dbz_r_c(ip)                                                                 

! Local variables

REAL (KIND=real_lsprec) :: p1upong      
              ! One upon g (sq seconds per m).

INTEGER :: multrhc
              ! Zero if (rhc_row_length*rhc_rows) le 1, else 1

INTEGER :: i,ii,ij,                                                           &
              ! Loop counters: I - horizontal field index; 
           irhi,irhj
              !  IRHI,IRHJ-indices for RHcrit.
           

REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LS_PPNC_GATHER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Multiple RHC rows?

IF ( (rhc_row_length * rhc_rows)  >   1) THEN
  multrhc = 1
ELSE
  multrhc = 0
END IF

!p1upong = 1.0_real_lsprec / g


!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO i=1, ip

  qcf2_c(i)   = 0.0_real_lsprec
  qrain_c(i)  = 0.0_real_lsprec
  qgraup_c(i) = 0.0_real_lsprec
  aero_c(i)   = 0.0_real_lsprec

  ii   =  ix(i+i1-1,1)
  ij   =  ix(i+i1-1,2)
  irhi = (multrhc * (ii - 1)) + 1
  irhj = (multrhc * (ij - 1)) + 1

  lsrain_c(i)          = REAL(lsrain(ii,ij),real_lsprec)
  lssnow_c(i)          = REAL(lssnow(ii,ij),real_lsprec)
  lssnow2_c(i)         = REAL(lssnow2(ii,ij),real_lsprec)
  lsgraup_c(i)         = REAL(lsgraup(ii,ij),real_lsprec)
  droplet_flux_c(i)    = REAL(droplet_flux(ii,ij),real_lsprec)
  cttemp_c(i)          = REAL(cttemp(ii,ij),real_lsprec)
  rainfrac_c(i)        = REAL(rainfrac(ii,ij),real_lsprec)
  rainfrac_impr_c(i)   = REAL(rainfrac_impr(ii,ij),real_lsprec)
  frac_ice_above_c(i)  = REAL(frac_ice_above(ii,ij),real_lsprec)
  p_c(i)               = REAL(p_theta_levels(ii,ij),real_lsprec)
! rhodz_c(i)           = REAL(-p1upong*layer_thickness(ii,ij),real_lsprec)
  deltaz_c(i)          = REAL(deltaz(ii,ij),real_lsprec)
  uk_c(i)              = REAL(uk(ii,ij),real_lsprec)
  vk_c(i)              = REAL(vk(ii,ij),real_lsprec)
  ukp1_c(i)            = REAL(ukp1(ii,ij),real_lsprec)
  vkp1_c(i)            = REAL(vkp1(ii,ij),real_lsprec)
  r_theta_levels_c(i)  = REAL(r_theta_levels(ii,ij,onelevel),real_lsprec)
  g_cos_theta_latitude_c(i)  = REAL(cos_theta_latitude(ii,ij),real_lsprec)
  rhodz_dry_c(i)   = REAL(rhodz_dry(ii,ij),real_lsprec)
  rhodz_moist_c(i) = REAL(rhodz_moist(ii,ij),real_lsprec)
  bland_c(i)       = bland(ii,ij) !This is a LOGICAL
  cf_c(i)          = REAL(cf(ii,ij),real_lsprec)
  cfl_c(i)         = REAL(cfl(ii,ij),real_lsprec)
  cff_c(i)         = REAL(cff(ii,ij),real_lsprec)
  qcf_c(i)         = REAL(qcf(ii,ij),real_lsprec)
  qcl_c(i)         = REAL(qcl(ii,ij),real_lsprec)
  q_c(i)           = REAL(q(ii,ij),real_lsprec)
  t_c(i)           = REAL(t(ii,ij),real_lsprec)
  n_drop_tpr_c(i)  = REAL(n_drop_tpr(ii, ij),real_lsprec)

  IF (l_mcr_qcf2)   qcf2_c(i)   = REAL(qcf2(ii,ij),real_lsprec)
  IF (l_mcr_qrain)  qrain_c(i)  = REAL(qrain(ii,ij),real_lsprec)
  IF (l_mcr_qgraup) qgraup_c(i) = REAL(qgraup(ii,ij),real_lsprec)
  IF (l_murk)       aero_c(i)   = REAL(aerosol(ii,ij),real_lsprec)

  !--------------------------------------------------------------------
  ! Pass through transfer diagnostics. Note that the 3D arrays are all
  ! initialized to zero in microphys_ctl, so the 1D versions below
  ! should just pick up zeros on the first microphysics iteration and
  ! if there are any subsequent iterations, then they should pick up
  ! the value of the diagnostic from previous iterations. Increments
  ! to each diagnostic are then added in the various transfer routines
  !--------------------------------------------------------------------

  IF (l_psdep_diag)    psdep_c(i) = REAL(psdep(ii,ij, level),real_lsprec)
  IF (l_psaut_diag)    psaut_c(i) = REAL(psaut(ii,ij, level),real_lsprec)
  IF (l_psacw_diag)    psacw_c(i) = REAL(psacw(ii,ij, level),real_lsprec)
  IF (l_psacr_diag)    psacr_c(i) = REAL(psacr(ii,ij, level),real_lsprec)
  IF (l_psaci_diag)    psaci_c(i) = REAL(psaci(ii,ij, level),real_lsprec)
  IF (l_psmlt_diag)    psmlt_c(i) = REAL(psmlt(ii,ij, level),real_lsprec)

  IF (l_psmltevp_diag) psmltevp_c(i) = REAL(psmltevp(ii,ij, level),         &
    real_lsprec)
  IF (l_praut_diag)    praut_c(i)    = REAL(praut(ii,ij, level),real_lsprec)
  IF (l_pracw_diag)    pracw_c(i)    = REAL(pracw(ii,ij, level),real_lsprec)
  IF (l_prevp_diag)    prevp_c(i)    = REAL(prevp(ii,ij, level),real_lsprec)
  IF (l_pgaut_diag)    pgaut_c(i)    = REAL(pgaut(ii,ij, level),real_lsprec)
  IF (l_pgacw_diag)    pgacw_c(i)    = REAL(pgacw(ii,ij, level),real_lsprec)
  IF (l_pgacs_diag)    pgacs_c(i)    = REAL(pgacs(ii,ij, level),real_lsprec)
  IF (l_pgmlt_diag)    pgmlt_c(i)    = REAL(pgmlt(ii,ij, level),real_lsprec)
  IF (l_pifrw_diag)    pifrw_c(i)    = REAL(pifrw(ii,ij, level),real_lsprec)
  IF (l_piprm_diag)    piprm_c(i)    = REAL(piprm(ii,ij, level),real_lsprec)
  IF (l_pidep_diag)    pidep_c(i)    = REAL(pidep(ii,ij, level),real_lsprec)
  IF (l_piacw_diag)    piacw_c(i)    = REAL(piacw(ii,ij, level),real_lsprec)
  IF (l_piacr_diag)    piacr_c(i)    = REAL(piacr(ii,ij, level),real_lsprec)
  IF (l_pimlt_diag)    pimlt_c(i)    = REAL(pimlt(ii,ij, level),real_lsprec)
  IF (l_pimltevp_diag) pimltevp_c(i) = REAL(pimltevp(ii,ij, level),         &
    real_lsprec)
  IF (l_pifall_diag)   pifall_c(i)   = REAL(pifall(ii,ij, level),real_lsprec)
  IF (l_psfall_diag)   psfall_c(i)   = REAL(psfall(ii,ij, level),real_lsprec)
  IF (l_prfall_diag)   prfall_c(i)   = REAL(prfall(ii,ij, level),real_lsprec)
  IF (l_pgfall_diag)   pgfall_c(i)   = REAL(pgfall(ii,ij, level),real_lsprec)
  IF (l_plset_diag)    plset_c(i)    = REAL(plset(ii,ij, level),real_lsprec)
  IF (l_plevpset_diag) plevpset_c(i) = REAL(plevpset(ii,ij, level),         &
    real_lsprec)
  IF (l_pifrr_diag)    pifrr_c(i)    = REAL(pifrr(ii, ij, level),real_lsprec)

  IF (l_sfwater_diag) sfwater_c(i) = REAL(sfwater(ii,ij, level),real_lsprec)
  IF (l_sfrain_diag) sfrain_c(i)   = REAL(sfrain(ii,ij, level),real_lsprec)
  IF (l_sfsnow_diag) sfsnow_c(i)   = REAL(sfsnow(ii,ij, level),real_lsprec)

  n_drop_out_c(i)  = REAL(n_drop_out(ii,ij),real_lsprec)
  vfall_c(i)       = REAL(vfall(ii,ij),real_lsprec)
  vfall2_c(i)      = REAL(vfall2(ii,ij),real_lsprec)
  vfall_rain_c(i)  = REAL(vfall_rain(ii,ij),real_lsprec)
  vfall_graup_c(i) = REAL(vfall_graup(ii,ij),real_lsprec)
  rhc_c(i)         = REAL(rhcrit(irhi,irhj),real_lsprec)

  land_fract_c(i) = REAL(land_fract(ii,ij),real_lsprec)
  hmteff_c(i)     = REAL(hmteff(ii,ij),real_lsprec)
  zb_c(i)         = REAL(zb(ii,ij),real_lsprec)

  IF ((l_fsd_generator) .AND. (ALLOCATED(f_arr))) THEN
    f_arr1_c(i) = REAL(f_arr(1,ii,ij,level),real_lsprec)
    f_arr2_c(i) = REAL(f_arr(2,ii,ij,level),real_lsprec)
    f_arr3_c(i) = REAL(f_arr(3,ii,ij,level),real_lsprec)
  END IF

END DO 

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ls_ppnc_gather


!-----------------------------------------------------------------------
!  Scatter variables, i.e post-process the data
!-----------------------------------------------------------------------

SUBROUTINE ls_ppnc_scatter( level, ix, ip, i1,                                &
  lsrain,lssnow,lssnow2,lsgraup, droplet_flux,                                &
  cf,cfl,cff,                                                                 &
  qcf,qcl,t, qcf2,qrain,qgraup,                                               &
  n_drop_tpr, n_drop_out,                                                     &
  aerosol, land_fract,                                                        &
  hmteff, zb,                                                                 &
!--------------------------------------------------------
! Layer thicknesses and variables passed from layer above
!--------------------------------------------------------
  q, p_theta_levels, layer_thickness,                                         &
  deltaz, rhodz_dry, rhodz_moist,                                             &
  rhc_row_length, rhc_rows, bland, rhcrit,                                    &
  vfall, vfall2, vfall_rain, vfall_graup,                                     &
  frac_ice_above, cttemp, rainfrac, rainfrac_impr,                            &
  uk, vk, ukp1, vkp1, onelevel,                                               &
!--------------------------------------------------------
! COSP flag
!--------------------------------------------------------
  l_cosp_lsp,                                                                 &
!--------------------------------------------------------
! Compressed arrays
!--------------------------------------------------------
  p_c, rhodz_c, deltaz_c, rhodz_dry_c,                                        &
  rhodz_moist_c, rhc_c, land_fract_c,                                         &
  hmteff_c, zb_c, qcf_c, qcl_c, q_c, qcf2_c, qrain_c,                         &
  qgraup_c, n_drop_tpr_c, n_drop_out_c,                                       &
  lsrain_c,vfall_rain_c, lssnow_c,vfall_c,                                    &
  lssnow2_c,vfall2_c, lsgraup_c,                                              &
  vfall_graup_c, droplet_flux_c,                                              &
  frac_ice_above_c, frac_agg_c,                                               &
  cttemp_c, rainfrac_c, rainfrac_impr_c,                                      &
  t_c,cf_c,cfl_c, cff_c, bland_c,                                             &
  psdep_c,  psaut_c,  psacw_c, psacr_c,                                       &
  psaci_c,  psmlt_c,  psmltevp_c,                                             &
  praut_c,  pracw_c,  prevp_c,                                                &
  pgaut_c,  pgacw_c,  pgacs_c,  pgmlt_c,                                      &
  pifrw_c,  pifrr_c,  piprm_c,  pidep_c,                                      &
  piacw_c,  piacr_c,  pimlt_c,  pimltevp_c,                                   &
  pifall_c, psfall_c, prfall_c, pgfall_c,                                     &
  plset_c,  plevpset_c, dbz_tot_c, dbz_g_c,                                   &
  dbz_i_c,  dbz_i2_c, dbz_l_c,  dbz_r_c,                                      &
  sfwater_c, sfrain_c, sfsnow_c,                                              &
  uk_c, vk_c, ukp1_c, vkp1_c,                                                 &
  r_theta_levels_c, g_cos_theta_latitude_c,                                   &
  r_theta_surf_c,                                                             &
  f_arr1_c, f_arr2_c, f_arr3_c,                                               &
  vm_cry_c, vm_agg_c, vtbranch_flag_c, vm_used_c, aero_c                      &
 )


IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                        &
 level,                                                                       &
              ! level number
 ip,                                                                          &
              ! Number of points
 i1,                                                                          &
              ! Starting point
 rhc_row_length,                                                              &
 rhc_rows,                                                                    &
 ix (  tdims%i_len * tdims%j_len , 2 )
              ! Gather/scatter index

REAL, INTENT(INOUT) ::                                                        &
  cf (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                  &
              ! Cloud fraction.
  cfl(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                  &
              ! Cloud liquid fraction.
  cff(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),                  &
              ! Cloud ice fraction.
  p_theta_levels( tdims%i_start:tdims%i_end,                                  &
                  tdims%j_start:tdims%j_end ),                                &
              !    
  layer_thickness(tdims%i_start:tdims%i_end,                                  &
                  tdims%j_start:tdims%j_end ),                                &
              ! Thickness of layer (Pa)
  deltaz( tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end ),                                        &
              ! Thickness of layer (m)
  rhodz_dry(tdims%i_start:tdims%i_end,                                        &
            tdims%j_start:tdims%j_end ),                                      &
              ! Dry air density
              ! * layer thickness (kg m-2)
  rhodz_moist(tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end ),                                    &
              ! Moist air density
              ! * layer thickness (kg m-2)
  rhcrit(rhc_row_length,rhc_rows)
              ! Critical humidity for cloud formation.

LOGICAL, INTENT(IN) ::                                                        &
  bland(tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end )
              ! Land/sea mask

LOGICAL, INTENT(IN) ::  l_cosp_lsp
              ! Switch for COSP LS diagnostics

REAL, INTENT(INOUT) ::                                                        &
  q(tdims%i_start:tdims%i_end,                                                &
    tdims%j_start:tdims%j_end),                                               &
              ! Specific humidity (kg water/kg air).
  qcf(tdims%i_start:tdims%i_end,                                              &
      tdims%j_start:tdims%j_end),                                             &
              ! Cloud ice (kg per kg air).
  qcl(tdims%i_start:tdims%i_end,                                              &
      tdims%j_start:tdims%j_end),                                             &
              ! Cloud liquid water (kg per kg air).
  qcf2(tdims%i_start:tdims%i_end,                                             &
       tdims%j_start:tdims%j_end),                                            &
              ! Cloud ice2 (kg per kg air).
  qrain(tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end),                                           &
              ! Rain water (kg per kg air).
  qgraup(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end),                                          &
              ! Graupel water (kg per kg air).
  t(tdims%i_start:tdims%i_end,                                                &
    tdims%j_start:tdims%j_end),                                               &
              ! Temperature (K).
  aerosol(tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end),                                         &
              ! Aerosol (K).
  lsrain(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end),                                          &
              ! Surface rainfall rate (kg m^-2 s^-1).
  lssnow(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end),                                          &
              ! Surface snowfall rate (kg m^-2 s^-1).
  lssnow2(tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end),                                         &
              ! Layer snowfall rate (kg m^-2 s^-1).
  lsgraup(tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end),                                         &
              ! Layer graupelfall rate (kg m^-2 s^-1)
  droplet_flux(tdims%i_start:tdims%i_end,                                     &
               tdims%j_start:tdims%j_end),                                    &
              ! Water droplet flux / kg m^-2 s^-1
  cttemp(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end),                                          &
              ! Ice cloud top temperature (K)
  rainfrac(tdims%i_start:tdims%i_end,                                         &
           tdims%j_start:tdims%j_end),                                        &
              ! Rain fraction.
  frac_ice_above(tdims%i_start:tdims%i_end,                                   &
                 tdims%j_start:tdims%j_end),                                  &
              ! Ice fraction from layer above
  vfall(tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end),                                           &
              ! Fall velocity of ice (m per
  vfall2(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end),                                          &
              ! Fall vel. of rain (m/s)
  vfall_rain(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end),                                      &
              ! Fall vel. of rain (m/s)
  vfall_graup(tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end)
              ! Fall vel. of graupel (m/s)

REAL, INTENT(IN) ::                                                           &
  rainfrac_impr(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end)
              ! Improved rain fraction

REAL,INTENT(IN) ::                                                            &
  ! U and V wind at levels k and k+1 on P grid, so defined using pdims.
  uk  (pdims%i_start : pdims%i_end,                                           &
       pdims%j_start : pdims%j_end),                                          &
  ukp1(pdims%i_start : pdims%i_end,                                           &
       pdims%j_start : pdims%j_end),                                          &
  vk  (pdims%i_start : pdims%i_end,                                           &
       pdims%j_start : pdims%j_end),                                          &
  vkp1(pdims%i_start : pdims%i_end,                                           &
       pdims%j_start : pdims%j_end)

! Level we are interested in for r_theta_level
INTEGER, INTENT(IN) :: onelevel

REAL,INTENT(INOUT) ::                                                         &
  land_fract(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end)

REAL,INTENT(IN) ::                                                            &
  hmteff(tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end)

REAL,INTENT(IN) ::                                                            &
  zb(tdims%i_start:tdims%i_end,                                               &
     tdims%j_start:tdims%j_end)

REAL,INTENT(INOUT) ::                                                         &
  n_drop_tpr(tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end )
              ! Tapered droplet number
REAL, INTENT(INOUT) ::                                                        &
  n_drop_out( tdims%i_start : tdims%i_end,                                    &
                  tdims%j_start : tdims%j_end )
              ! Droplet number from autoconversion

! Compressed arrays ---------------

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  cf_c(ip),                                                                   &
  q_c(ip),                                                                    &
  qcf_c(ip),                                                                  &
  qcl_c(ip),                                                                  &
  qcf2_c(ip),                                                                 &
  qrain_c(ip),                                                                &
  qgraup_c(ip),                                                               &
  t_c(ip),                                                                    &
  uk_c(ip),                                                                   &
  vk_c(ip),                                                                   &
  ukp1_c(ip),                                                                 &
  vkp1_c(ip),                                                                 &
  r_theta_levels_c(ip),                                                       &
  r_theta_surf_c(ip),                                                         &
  g_cos_theta_latitude_c(ip),                                                 &
  aero_c(ip),                                                                 &
  lsrain_c(ip),                                                               &
  lssnow_c(ip),                                                               &
  lssnow2_c(ip),                                                              &
  lsgraup_c(ip),                                                              &
  droplet_flux_c(ip),                                                         &
  cttemp_c(ip),                                                               &
  rainfrac_c(ip),                                                             &
  rainfrac_impr_c(ip),                                                        &
  frac_ice_above_c(ip),                                                       &
  frac_agg_c(ip),                                                             &
  cfl_c(ip),                                                                  &
  cff_c(ip),                                                                  &
  vfall_c(ip),                                                                &
  vfall2_c(ip),                                                               &
  vfall_rain_c(ip),                                                           &
  vfall_graup_c(ip),                                                          &
  rhc_c(ip),                                                                  &
  n_drop_tpr_c(ip),                                                           &
  n_drop_out_c(ip),                                                           &
  land_fract_c(ip),                                                           &
  hmteff_c(ip),                                                               &
  zb_c(ip),                                                                   &
  f_arr1_c(ip),                                                               &
  f_arr2_c(ip),                                                               &
  f_arr3_c(ip)

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  rhodz_c(ip),                                                                &
  deltaz_c(ip),                                                               &
  rhodz_dry_c(ip),                                                            &
  rhodz_moist_c(ip),                                                          &
  p_c(ip)  

LOGICAL, INTENT(IN) :: bland_c(ip)          ! gathered land/sea mask

! Microphysical process rate diagnostics (compressed arrays)
REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  psdep_c(ip),                                                                &
  psaut_c(ip),                                                                &
  psacw_c(ip),                                                                &
  psacr_c(ip),                                                                &
  psaci_c(ip),                                                                &
  psmlt_c(ip),                                                                &
  psmltevp_c(ip),                                                             &
  praut_c(ip),                                                                &
  pracw_c(ip),                                                                &
  prevp_c(ip),                                                                &
  pgaut_c(ip),                                                                &
  pgacw_c(ip),                                                                &
  pgacs_c(ip),                                                                &
  pgmlt_c(ip),                                                                &
  pifrw_c(ip),                                                                &
  pifrr_c(ip),                                                                &
  piprm_c(ip),                                                                &
  pidep_c(ip),                                                                &
  piacw_c(ip),                                                                &
  piacr_c(ip),                                                                &
  pimlt_c(ip),                                                                &
  pimltevp_c(ip),                                                             &
  pifall_c(ip),                                                               &
  psfall_c(ip),                                                               &
  prfall_c(ip),                                                               &
  pgfall_c(ip),                                                               &
  plset_c(ip),                                                                &
  plevpset_c(ip),                                                             &
  vm_agg_c(ip),                                                               &
  vm_cry_c(ip),                                                               &
  vtbranch_flag_c(ip),                                                        &
  vm_used_c(ip)

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  sfwater_c(ip),                                                              &
  sfrain_c(ip),                                                               &
  sfsnow_c(ip)

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  dbz_tot_c(ip),                                                              &
  dbz_g_c(ip),                                                                &
  dbz_i_c(ip),                                                                &
  dbz_i2_c(ip),                                                               &
  dbz_l_c(ip),                                                                &
  dbz_r_c(ip)                                                                  

! Local variables

INTEGER :: i,ii,ij
              ! Loop counters: I - horizontal field index; 
           

REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LS_PPNC_SCATTER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! assumption that all pairs ii,ij generated for each i are distinct
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO i=1, ip

  ii         = ix(i+i1-1,1)
  ij         = ix(i+i1-1,2)

  t(ii,ij)   = REAL(t_c(i))
  q(ii,ij)   = REAL(q_c(i))
  qcf(ii,ij) = REAL(qcf_c(i))
  qcl(ii,ij) = REAL(qcl_c(i))

  IF (l_mcr_qcf2) THEN
    qcf2(ii,ij) = REAL(qcf2_c(i))
  END IF

  IF (l_mcr_qrain) THEN
    qrain(ii,ij) = REAL(qrain_c(i))
  END IF

  IF (l_mcr_qgraup) THEN
    qgraup(ii,ij) = REAL(qgraup_c(i))
  END IF

  IF (l_murk) THEN
    aerosol(ii,ij) = REAL(aero_c(i))
  END IF

  n_drop_out(ii,ij)     = REAL(n_drop_out_c(i))
  lsrain(ii,ij)         = REAL(lsrain_c(i))
  lssnow(ii,ij)         = REAL(lssnow_c(i))
  lssnow2(ii,ij)        = REAL(lssnow2_c(i))
  lsgraup(ii,ij)        = REAL(lsgraup_c(i))
  droplet_flux(ii,ij)   = REAL(droplet_flux_c(i))
  cttemp(ii,ij)         = REAL(cttemp_c(i))
  rainfrac(ii,ij)       = REAL(rainfrac_c(i))
  frac_ice_above(ii,ij) = REAL(frac_ice_above_c(i))

  IF (i_cld_vn == i_cld_pc2) THEN
    cff(ii,ij)  = REAL(cff_c(i))
    cfl(ii,ij)  = REAL(cfl_c(i))
    cf(ii,ij)   = REAL(cf_c(i))
  END IF  ! i_cld_pc2

  vfall(ii,ij)       = REAL(vfall_c(i))
  vfall2(ii,ij)      = REAL(vfall2_c(i))
  vfall_rain(ii,ij)  = REAL(vfall_rain_c(i))
  vfall_graup(ii,ij) = REAL(vfall_graup_c(i))

  ! Only store process rates in array for diagnostic
  ! if a particular diagnostic is requested,
  ! otherwise overwriting will occur
  ! (space for the 3D array in microphys_ctl is only allocated
  ! if the diagnostic is active, to save memory)
  IF (l_aggfr_diag .OR. l_cosp_lsp)                                         &
    frac_agg(ii,ij, level) = REAL(frac_agg_c(i))
  IF (l_psdep_diag)    psdep(ii,ij, level)    = REAL(psdep_c(i))
  IF (l_psaut_diag)    psaut(ii,ij, level)    = REAL(psaut_c(i))
  IF (l_psacw_diag)    psacw(ii,ij, level)    = REAL(psacw_c(i))
  IF (l_psacr_diag)    psacr(ii,ij, level)    = REAL(psacr_c(i))
  IF (l_psaci_diag)    psaci(ii,ij, level)    = REAL(psaci_c(i))
  IF (l_psmlt_diag)    psmlt(ii,ij, level)    = REAL(psmlt_c(i))
  IF (l_psmltevp_diag) psmltevp(ii,ij, level) = REAL(psmltevp_c(i))
  IF (l_praut_diag)    praut(ii,ij, level)    = REAL(praut_c(i))
  IF (l_pracw_diag)    pracw(ii,ij, level)    = REAL(pracw_c(i))
  IF (l_prevp_diag)    prevp(ii,ij, level)    = REAL(prevp_c(i))
  IF (l_pgaut_diag)    pgaut(ii,ij, level)    = REAL(pgaut_c(i))
  IF (l_pgacw_diag)    pgacw(ii,ij, level)    = REAL(pgacw_c(i))
  IF (l_pgacs_diag)    pgacs(ii,ij, level)    = REAL(pgacs_c(i))
  IF (l_pgmlt_diag)    pgmlt(ii,ij, level)    = REAL(pgmlt_c(i))
  IF (l_pifrw_diag)    pifrw(ii,ij, level)    = REAL(pifrw_c(i))
  IF (l_piprm_diag)    piprm(ii,ij, level)    = REAL(piprm_c(i))
  IF (l_pidep_diag)    pidep(ii,ij, level)    = REAL(pidep_c(i))
  IF (l_piacw_diag)    piacw(ii,ij, level)    = REAL(piacw_c(i))
  IF (l_piacr_diag)    piacr(ii,ij, level)    = REAL(piacr_c(i))
  IF (l_pimlt_diag)    pimlt(ii,ij, level)    = REAL(pimlt_c(i))
  IF (l_pimltevp_diag) pimltevp(ii,ij, level) = REAL(pimltevp_c(i))
  IF (l_pifall_diag)   pifall(ii,ij, level)   = REAL(pifall_c(i))
  IF (l_psfall_diag)   psfall(ii,ij, level)   = REAL(psfall_c(i))
  IF (l_prfall_diag)   prfall(ii,ij, level)   = REAL(prfall_c(i))
  IF (l_pgfall_diag)   pgfall(ii,ij, level)   = REAL(pgfall_c(i))
  IF (l_plset_diag)    plset(ii,ij, level)    = REAL(plset_c(i))
  IF (l_plevpset_diag) plevpset(ii,ij, level) = REAL(plevpset_c(i))
  IF (l_pifrr_diag)    pifrr(ii, ij, level)   = REAL(pifrr_c(i))
  IF (l_vtbranch_diag) vtbranch_flag(ii, ij, level)                         &
    = REAL(vtbranch_flag_c(i))
  IF (l_vm_cry_diag)   vm_cry(ii, ij, level)  = REAL(vm_cry_c(i))
  IF (l_vm_agg_diag)   vm_agg(ii, ij, level)  = REAL(vm_agg_c(i))
  IF (l_vm_used_diag)  vm_used(ii, ij, level) = REAL(vm_used_c(i))

  ! Seeder feeder diagnostics
  IF (l_sfwater_diag)  sfwater(ii,ij, level)  = REAL(sfwater_c(i))
  IF (l_sfrain_diag)   sfrain(ii,ij, level)   = REAL(sfrain_c(i))
  IF (l_sfsnow_diag)   sfsnow(ii,ij, level)   = REAL(sfsnow_c(i))

  ! Radar reflectivity diagnostics
  IF (l_ref_diag) THEN
    dbz_tot(ii, ij, level) = REAL(dbz_tot_c(i))
    dbz_g(ii, ij, level)   = REAL(dbz_g_c(i))
    dbz_i(ii, ij, level)   = REAL(dbz_i_c(i))
    dbz_i2(ii, ij, level)  = REAL(dbz_i2_c(i))
    dbz_r(ii, ij, level)   = REAL(dbz_r_c(i))
    dbz_l(ii, ij, level)   = REAL(dbz_l_c(i))
  END IF

END DO 


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ls_ppnc_scatter


END MODULE ls_ppnc_mod

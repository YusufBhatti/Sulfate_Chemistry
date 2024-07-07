! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! ---------------------------------------------------------------------
! Purpose: Interface level to tr_mix for tracer variables that
!          require boundary layer mixing and/or dry deposition
!          with JULES surface scheme.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Boundary Layer

! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

! ---------------------------------------------------------------------
MODULE bl_trmix_dd_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'BL_TRMIX_DD_MOD'
CONTAINS

SUBROUTINE bl_trmix_dd(                                                 &
! IN arguments
        bl_levels, dtrdz_charney_grid, tr_vars,                         &
        alpha_cd, rhokh_mix,  p_star,                                   &
! IN Emissions fields
        dust_flux, co2_emits, co2flux, npp, resp_s,                     &
! IN variables needed for tr_mix
        kent, we_lim, t_frac, zrzi,                                     &
        kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl, zhsc, z_uv,   &
! IN for dry deposition of tracers
        rho_aresist, aresist,                                           &
        r_b_dust,tstar, land_points, land_index, ice_frac,              &
! IN variables for sresfact
        ntiles, tile_pts, tile_index, tile_frac,                        &
        canopy, catch, snow_tile, gc,                                   &
        aresist_tile, flandg,                                           &
! INOUT Fields to mix
        murk, free_tracers, ozone_tracer, drydep_str, resist_b,         &
        resist_b_tile,                                                  &
! INOUT Mineral Dust
        dust_div1,dust_div2,dust_div3,                                  &
        dust_div4,dust_div5,dust_div6,                                  &
! INOUT Sulphur cycle
        so2, dms, so4_aitken, so4_accu, so4_diss, nh3,                  &
! INOUT Soot cycle
        soot_new, soot_aged, soot_cld,                                  &
! INOUT Biomass aerosol
        bmass_new, bmass_agd, bmass_cld,                                &
! INOUT Fossil-fuel organic carbon aerosol
        ocff_new, ocff_agd, ocff_cld,                                   &
! INOUT Ammonium nitrate aerosol
        nitr_acc, nitr_diss,                                            &
! INOUT Carbon cycle
        co2, co2_flux_tot, land_co2,                                    &
! STASH related variables
        stashwork                                                       &
        )

USE atm_fields_bounds_mod, ONLY:                                        &
 pdims, tdims, tdims_s, pdims_s
USE atm_step_local, ONLY: dim_cs2
USE c_bm_bdy_mod, ONLY: resb_freshbmass, resb_agedbmass,                &
                        resb_bmassincloud, ress_bmass
USE c_nitrate_bdy_mod, ONLY: resb_nitr_acc, resb_nitr_diss,             &
                             ress_nitr_acc, ress_nitr_diss
USE c_ocff_bdy_mod, ONLY: resb_freshocff, resb_agedocff,                &
                          resb_ocffincloud, ress_ocff
USE c_st_bdy_mod, ONLY: resb_freshsoot, resb_agedsoot,                  &
                        resb_sootincloud, ress_soot
USE c_sulbdy_mod, ONLY: resb_so2, resb_nh3, resb_so4_ait,               &
                        resb_so4_acc, resb_so4_dis, ress_so2,           &
                        ress_nh3, ress_so4_ait, ress_so4_acc,           &
                        ress_so4_dis, cond_lim, r_snow, asnow

USE ccarbon, ONLY: m_co2, m_carbon
USE conversions_mod, ONLY: rsec_per_day
USE dust_parameters_mod, ONLY: ndiv, rhop, drep, l_twobin_dust,         &
                               l_dust, tot_dust_dep_flux
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE free_tracers_inputs_mod, ONLY: l_bl_tracer_mix
USE jules_surface_types_mod, ONLY: ntype
USE missing_data_mod, ONLY: rmdi
USE model_domain_mod, ONLY: model_type, mt_single_column
USE murk_inputs_mod, ONLY: l_murk_advect
USE um_parvars, ONLY: at_extremity
USE rad_input_mod, ONLY: l_use_cariolle
USE run_aerosol_mod, ONLY: l_sulpc_so2,                                 &
  l_sulpc_nh3, l_sulpc_dms, l_soot, l_biomass, l_ocff, l_nitrate
USE stash_array_mod, ONLY:                                              &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE submodel_mod,  ONLY: atmos_im
USE s_scmop_mod,   ONLY: default_streams,                               &
                         t_avg, d_bl, scmdiag_bl
USE scmoutput_mod, ONLY: scmoutput

USE tr_mix_mod, ONLY: tr_mix

USE carbon_options_mod, ONLY: l_co2_interactive, l_co2_emits

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE jules_vegetation_mod, ONLY: l_triffid

USE prognostics, ONLY: triffid_co2_gb

USE sresfact_mod, ONLY: sresfact
IMPLICIT NONE

! Model dimensions
INTEGER, INTENT(IN) ::                                                  &
  bl_levels,                                                            &
  tr_vars
                        ! number of free tracer variables
! Model parameters
REAL, INTENT(IN) ::                                                     &
  alpha_cd(bl_levels),                                                  &
  rhokh_mix (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             bl_levels),                                                &
  we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                                    !rho*entrainment rate implied by
                                    !placing of subsidence
  zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),          &
                                    !(z-z_base)/(z_i-z_base)
  t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                                    !a fraction of the timestep
  we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
                                    !rho*entrainment rate implied by
                                    !  placing of subsidence
  zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                                    !(z-z_base)/(z_i-z_base)
  t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
                                    !a fraction of the timestep
  z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,             &
         bl_levels),                                                    &
                                    !Z_uv(*,K) is height of half
                                    !    level k-1/2.
  zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                    !Top of decoupled layer
  zhnl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                    !Top of surface mixed layer
  dtrdz_charney_grid(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end,bl_levels)
                                    ! For tracer mixing

INTEGER, INTENT(IN) ::                                                  &
  kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                    !grid-level of SML inversion
  kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                    !grid-level of DSC inversion
  land_points,                                                          &
                                    !No.of land points being processed
  land_index(land_points),                                              &
                                    !set from land_sea_mask
! For sresfact
    ntiles,                                                             &
                                      !No.of land-surface tiles
    tile_pts(ntype),                                                    &
                                      !No.of tile points
    tile_index(land_points,ntype) !Index of tile points

! Required fields (input)
REAL, INTENT(IN) ::                                                     &
  p_star(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! Emissions fields (input)
REAL, INTENT(IN) ::                                                     &
 dust_flux(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,ndiv)
                                   ! dust emission flux

REAL, INTENT(IN) ::                                                     &
  co2_emits (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                      ! anthro CO2 emissions
                                      !  (Kg CO2 m-2 s-1)

  co2flux   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                      ! ocean  CO2 emissions
                                      !  (Kg CO2 m-2 s-1)

! For dry deposition of tracers (input)
REAL, INTENT(IN) ::                                                     &
  rho_aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                     !RHOSTAR*CD_STD*VSHR
                                     !for CLASSIC aerosol scheme
  aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                     !1/(CD_STD*VSHR)
  r_b_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,ndiv),   &
                                        !surf layer res for dust
  tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
                             ! temperature
  ice_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                     !sea_ice fractn (ice/qrclim.ice)
! For sresfact
    tile_frac(land_points,ntiles),                                      &
                                       !fractional coverage for each
!                                       surface tile
    canopy(land_points,ntiles),                                         &
                                       !Surface/canopy water (kg/m2)
    catch(land_points,ntiles),                                          &
                                       !Surface/canopy water capacity
!                                       of snow-free land tiles (kg/m2)
    snow_tile(land_points,ntiles),                                      &
                                       !Snow on tiles (kg/m2).
    gc(land_points,ntiles),                                             &
                                       !Stomatal conductance to evapn
!                                       for land tiles (m/s)
    aresist_tile(land_points,ntiles),                                   &
                                       ! 1/(CD_STD*VSHR) on land tiles
                                       !for CLASSIC aerosol scheme
    flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),&
                                     !land fraction of grid box
!                                       (0 if all sea, 1 if all land)
    npp(land_points),                                                   &
                                       ! net primary productivity
                                       !  (Kg C m-2 s-1)
    resp_s(dim_cs2)
                                       ! soil respiration
                                       !  (Kg C m-2 s-1)

! Tracer fields for mixing
REAL, INTENT(INOUT) ::                                                  &
  murk(tdims_s%i_start:tdims_s%i_end,                                   &
       tdims_s%j_start:tdims_s%j_end, tdims_s%k_start:tdims_s%k_end),   &
  free_tracers(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end,                           &
               tr_vars)

REAL, INTENT(INOUT) ::                                                  &
  resist_b_tile(land_points,ntiles),                                    &
                                     !(1/CH-1/CD_STD)/VSHR on land tiles
                                       !for CLASSIC aerosol scheme
  drydep_str(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                     !(1/CH-1/CD_STD)/VSHR
  co2_flux_tot(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                     ! total CO2 flux
                                     !  (Kg CO2 m-2 s-1)
  land_co2(land_points)             ! terrestrial CO2 flux
                                     !  (Kg CO2 m-2 s-1)

REAL, INTENT(INOUT) ::                                                  &
  dust_div1   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  dust_div2   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  dust_div3   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  dust_div4   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  dust_div5   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  dust_div6   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::                                                  &
  so2         (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  dms         (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  so4_aitken  (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  so4_accu    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  so4_diss    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  nh3         (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::                                                  &
  soot_new    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  soot_aged   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  soot_cld    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  bmass_new   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  bmass_agd   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  bmass_cld   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  ocff_new    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  ocff_agd    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  ocff_cld    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  nitr_acc    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  nitr_diss   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  co2         (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  ozone_tracer(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::  stashwork(*)

! Local variables

LOGICAL, PARAMETER :: soluble = .TRUE.
LOGICAL, PARAMETER :: insoluble = .FALSE.

INTEGER, PARAMETER :: sect = 3               ! BL section
INTEGER            :: item                   ! stash item code
INTEGER            :: im_index               ! internal model
INTEGER            :: n_tracer               ! looper
INTEGER            :: errorstatus            ! error reporting
INTEGER            :: idiv            ! loop counter, dust divs
INTEGER            :: lev1            ! =1 no. levs for vgrav calc
CHARACTER (LEN=2)  :: cdiv            ! character string for idiv

REAL        :: zero_field(pdims%i_start:pdims%i_end,                    &
                          pdims%j_start:pdims%j_end) ! field of 0s
REAL        :: tr_flux( pdims%i_start:pdims%i_end,                      &
                        pdims%j_start:pdims%j_end,bl_levels)
                                                     ! Fluxes returned
                                                     ! from mixing

CHARACTER (LEN=*), PARAMETER :: RoutineName='BL_TRMIX_DD'
CHARACTER (LEN=errormessagelength)           :: cmessage

INTEGER :: i, j, k, l                !loop counters

REAL ::                                                                 &
  work_2d(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                      !2d work array for stash
  str_resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                      !Rb for S Cycle dry dep
  str_resist_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                      !Rs for S Cycle dry dep
  res_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                      !Ra/(Ra+Rb+Rs) for dry dep
 res_factor_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
                                      !Ra/(Ra+Rb+Rs) mean over land
  resist_s(pdims%i_end*pdims%j_end),                                    &
                                      !stomatal resistance
  vstokes1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                               !dust settling velocity, lowest layer
  work1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels,ndiv),                                                &
                                            ! workspace
  work2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
  work3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!---------------------------------------------------------------------
! Initial setup
!---------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
zero_field(:,:) = 0.0
errorstatus = 0
im_index = 1
!  Initialise RES_FACTOR array to zero to prevent use of unset values
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    res_factor(i,j) = 0.0
  END DO
END DO

!---------------------------------------------------------------------
! Mixing for Aerosol
!---------------------------------------------------------------------
IF ( l_murk_advect ) THEN

  CALL tr_mix(                                                          &
  ! IN fields
           bl_levels,alpha_cd,rhokh_mix(:,:,2:), rhokh_mix(:,:,1),      &
           dtrdz_charney_grid, zero_field, zero_field,                  &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
  ! INOUT / OUT fields
           murk(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                1:bl_levels ),tr_flux,drydep_str                        &
           )

  !   Output BL flux diagnostic
  IF ( model_type == mt_single_column ) THEN
    CALL SCMoutput(tr_flux,                                             &
                   'murk_flux','BL flux of MURK','m/s',                 &
                   t_avg,d_bl,default_streams,'',RoutineName)
  ELSE
    item = 129
    IF ( sf(item,sect) ) THEN      ! diagnostic flux required
      ! DEPENDS ON: copydiag_3d
      CALL copydiag_3d( stashwork(si(item,sect,im_index)),              &
           tr_flux,                                                     &
           pdims%i_end,pdims%j_end,bl_levels,0,0,0,0, at_extremity,     &
           stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,sect,item,                                          &
           errorstatus ,cmessage)
    END IF
    IF (errorstatus /=0) THEN
      CALL ereport( RoutineName, errorstatus, cmessage )
    END IF
  END IF  ! model_type

END IF  ! ( l_murk_advect )

!------------------------------------------------------------------
! Mixing for mineral dust, including dry deposition
!------------------------------------------------------------------

IF (l_dust) THEN

  IF ( sf(440,3) ) THEN
    ! Allocate total dust depostion flux in the event that both
    ! sections 4 and 5 have not been called (otherwise, either
    ! should have already allocated this array)
    IF (.NOT. ALLOCATED(tot_dust_dep_flux)) THEN
      ALLOCATE(tot_dust_dep_flux(pdims%i_len,pdims%j_len))
      tot_dust_dep_flux(:,:) = 0.0
    END IF
  END IF

  lev1 = 1

  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        work1(i,j,k,1)=dust_div1(i,j,k)
        work1(i,j,k,2)=dust_div2(i,j,k)
        IF (.NOT. l_twobin_dust) THEN
          work1(i,j,k,3)=dust_div3(i,j,k)
          work1(i,j,k,4)=dust_div4(i,j,k)
          work1(i,j,k,5)=dust_div5(i,j,k)
          work1(i,j,k,6)=dust_div6(i,j,k)
        END IF
      END DO !i
    END DO !j
  END DO !BL_LEVELS

  DO idiv = 1, ndiv

    ! DEPENDS ON: vgrav
    CALL vgrav(                                                         &
  lev1,drep(idiv),                                                      &
  rhop,p_star,tstar,                                                    &
  vstokes1,work2,work3                                                  &
  )
    !     Do not allow tracer to settle directly from level 1
    !     but treat this in parallel with tracer mixing here.

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j)=aresist(i,j) * ( vstokes1(i,j) + 1.0/           &
             (aresist(i,j)+r_b_dust(i,j,idiv)+                          &
             aresist(i,j)*r_b_dust(i,j,idiv)*vstokes1(i,j)) )
      END DO !i
    END DO !j

    CALL tr_mix(                                                        &
    ! IN fields
           bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,           &
           dtrdz_charney_grid,                                          &
           dust_flux(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end, idiv), res_factor,      &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
    ! INOUT / OUT fields
           work1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,   &
                  1:bl_levels,idiv ),tr_flux,drydep_str                 &
           )

    !     Output BL flux diagnostic
    IF ( model_type == mt_single_column ) THEN
      WRITE(cdiv,'(I2.2)')idiv
      CALL SCMoutput(tr_flux,                                           &
           'dust_flux'//cdiv,'BL flux of dust category '//cdiv,'m/s',   &
           t_avg,d_bl,default_streams,'',RoutineName)
    ELSE
      item = 440 + idiv      !dust dry dep flux, from lowest layer
      IF ( sf(item,sect) ) THEN
        !         Change sign of dry dep flux (otherwise negative)
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            work_2d(i,j) = -drydep_str(i,j)
          END DO
        END DO
        ! DEPENDS ON: copydiag
        CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,        &
                      pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,     &
                      atmos_im,sect,item,                               &
                      errorstatus,cmessage)
        IF (errorstatus /=0) THEN
          CALL ereport( RoutineName, errorstatus, cmessage )
        END IF
      END IF  ! stashflag

      ! Add dust dry dep from lev1 to total dust dep field. (Note that
      ! drydep_str needs a change of sign, hence why there is a
      ! subtraction below and not an addition.)
      IF ( sf(440,3) ) THEN
        tot_dust_dep_flux(pdims%i_start:pdims%i_end,                    &
                          pdims%j_start:pdims%j_end) =                  &
          tot_dust_dep_flux(pdims%i_start:pdims%i_end,                  &
                            pdims%j_start:pdims%j_end) -                &
          drydep_str(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end)
      END IF
    END IF  ! model_type

  END DO !NDIV

  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        dust_div1(i,j,k)=work1(i,j,k,1)
        dust_div2(i,j,k)=work1(i,j,k,2)
        IF (.NOT. l_twobin_dust) THEN
          dust_div3(i,j,k)=work1(i,j,k,3)
          dust_div4(i,j,k)=work1(i,j,k,4)
          dust_div5(i,j,k)=work1(i,j,k,5)
          dust_div6(i,j,k)=work1(i,j,k,6)
        END IF
      END DO !i
    END DO !j
  END DO !BL_LEVELS

END IF !L_DUST
!------------------------------------------------------------------
! Mixing for Sulphur Cycle variables, including dry deposition
!------------------------------------------------------------------

IF (l_sulpc_so2) THEN

  !  Calculate resistance factors for species to be dry deposited.

  ! Note that for JULES with coastal tiling, at coastal points
  ! (i.e. those containing both land and sea):
  ! RESIST_B here contains  values approriate for SEA only (i.e. no
  ! meaning including land tiles has been done in SF_EXCH), but
  ! ARESIST and RHO_ARESIST contain grid box mean values calculated
  ! in SF_EXCH taking account of land/sea fractions.

  !   For RESIST_B check to eliminate possible negative values
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF (resist_b(i,j)  <   0.0) THEN
        resist_b(i,j) = 0.0
      END IF
    END DO
  END DO

  !   For STR_RESIST_S values depend on surface type (land, sea,
  !    snow,ice) as well as tracer identity.

  !   Initialise 'stomatal' resistance array to zero:

  DO k = 1, pdims%i_end*pdims%j_end
    resist_s(k)=0.0
  END DO

  ! CODE FOR SO2

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) <  1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_so2*resist_b(i,j)
        str_resist_s(i,j) = ress_so2*resist_s(k)

        !  Need to set STR_RESIST_S = R_SNOW for SO2 over sea ice
        !  (0 is acceptable over sea).
        IF (ice_frac(i,j)  >   0.0 ) THEN
          str_resist_s(i,j)=r_snow
        END IF

        !  Calculate RES_FACTOR for SO2 (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for SO2 to calculate correct RES_FACTOR_LAND values
  ! if land present:
  IF (land_points >  0) THEN

    CALL sresfact (land_points, land_index,                             &
                 ntiles, tile_index, tile_pts, soluble,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_so2, ress_so2,             &
                 resist_b_tile, res_factor_land)

    ! Recalculate RES_FACTOR for SO2 combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix SO2   (Note emiss added in AERO_CTL via TRSRCE call)

  CALL tr_mix(                                                          &
  ! IN fields
           bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,          &
           dtrdz_charney_grid, zero_field, res_factor,                  &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
  ! INOUT / OUT fields
           so2(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               1:bl_levels ), tr_flux, drydep_str                       &
           )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 270                      !dry dep flux SO2
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! ( sf(item,sect) )
  END IF  ! ( model_type /= mt_single_column )

  ! CODE FOR NH3 (if present)

  IF (l_sulpc_nh3) THEN

    !  Calculate RES_FACTOR for NH3 to allow dry deposition in same way
    !  as for SO2 (including code for snow and ice)

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        k=(j-1)*pdims%i_end + i

        IF (flandg(i,j) <  1.0) THEN
          !  For grid boxes containing sea:
          str_resist_b(i,j) = resb_nh3*resist_b(i,j)
          str_resist_s(i,j) = ress_nh3*resist_s(k)

          !  Need to set STR_RESIST_S = R_SNOW for SO2 over sea ice
          !  (0 is acceptable over sea).
          IF (ice_frac(i,j)  >   0.0 ) THEN
            str_resist_s(i,j)=r_snow
          END IF

          !  Calculate RES_FACTOR for NH3  (only correct for sea points)

          res_factor(i,j) = aresist(i,j) /                              &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

        END IF

      END DO
    END DO

    ! Call SRESFACT for NH3 to calculate correct RES_FACTOR_LAND values
    ! if land present:
    IF (land_points >  0) THEN

      CALL sresfact (land_points, land_index,                           &
                 ntiles, tile_index, tile_pts, soluble,                 &
                 canopy, catch, gc, tile_frac, snow_tile,               &
                 aresist, aresist_tile, resb_nh3, ress_nh3,             &
                 resist_b_tile, res_factor_land)

      ! Recalculate RES_FACTOR for NH3 combining sea and land values
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +         &
                         flandg(i,j)*res_factor_land(i,j)
        END DO
      END DO

    END IF

    ! Mix NH3   (Note emiss added in AERO_CTL via TRSRCE call)

    CALL tr_mix(                                                        &
    ! IN fields
           bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,          &
           dtrdz_charney_grid, zero_field, res_factor,                  &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
    ! INOUT/OUT fields
           nh3(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,    &
               1:bl_levels ), tr_flux, drydep_str                       &
             )

    !     Output BL flux diagnostic
    IF ( model_type /= mt_single_column ) THEN
      item = 300                      !dry dep flux NH3
      IF ( sf(item,sect) ) THEN
        !         Change sign of dry dep flux (otherwise negative)
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            work_2d (i,j) = -drydep_str(i,j)
          END DO
        END DO
        ! DEPENDS ON: copydiag
        CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,        &
                      pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,     &
                      atmos_im,sect,item,                               &
                      errorstatus,cmessage)
        IF (errorstatus /=0) THEN
          CALL ereport( RoutineName, errorstatus, cmessage )
        END IF
      END IF  ! ( sf(item,sect) )
    END IF  ! ( model_type /= mt_single_column )

  END IF                        ! End If nh3 condition

  ! CODE FOR SO4_AIT

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end+ i

      IF (flandg(i,j) <  1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_so4_ait*resist_b(i,j)
        str_resist_s(i,j) = ress_so4_ait*resist_s(k)

        !  Calculate RES_FACTOR for SO4_AIT  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for SO4_AIT to calculate correct RES_FACTOR_LAND values
  ! if land present:
  IF (land_points >  0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_so4_ait, ress_so4_ait,       &
               resist_b_tile, res_factor_land)

    ! Recalculate RES_FACTOR for SO4_AIT combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix so4_aitken

  CALL tr_mix(                                                          &
  ! IN fields
         bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,             &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,  &
  ! INOUT / OUT fields
         so4_aitken(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    1:bl_levels ), tr_flux, drydep_str                  &
         )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 271                      !dry dep flux SO4_AIT
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! ( sf(item,sect) )
  END IF  ! ( model_type /= mt_single_column )

  ! CODE FOR SO4_ACC

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) <  1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_so4_acc*resist_b(i,j)
        str_resist_s(i,j) = ress_so4_acc*resist_s(k)

        !  Calculate RES_FACTOR for SO4_ACC  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for SO4_ACC to calculate correct RES_FACTOR_LAND values
  ! if land present:
  IF (land_points >  0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_so4_acc, ress_so4_acc,       &
               resist_b_tile, res_factor_land)

    ! Recalculate RES_FACTOR for SO4_ACC combining sea and land values
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix so4_accu

  CALL tr_mix(                                                          &
  ! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,  &
  ! INOUT / OUT fields
         so4_accu(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, &
                  1:bl_levels ), tr_flux, drydep_str                    &
         )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 272                      !dry dep flux SO4_ACC
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! ( sf(item,sect) )
  END IF  ! ( model_type /= mt_single_column )

  ! CODE FOR SO4_DIS

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) <  1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_so4_dis*resist_b(i,j)
        str_resist_s(i,j) = ress_so4_dis*resist_s(k)

        !  Calculate RES_FACTOR for SO4_DIS  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for SO4_ACC to calculate correct RES_FACTOR_LAND values
  ! if land present:
  IF (land_points >  0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_so4_dis, ress_so4_dis,       &
               resist_b_tile, res_factor_land)

    ! Recalculate RES_FACTOR for SO4_DIS combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix so4_diss

  CALL tr_mix(                                                          &
  ! IN fields
           bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,          &
           dtrdz_charney_grid, zero_field, res_factor,                  &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
  ! INOUT / OUT fields
           so4_diss(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    1:bl_levels ), tr_flux, drydep_str                  &
           )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 273                      !dry dep flux SO4_DIS
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! ( sf(item,sect) )
  END IF  ! ( model_type /= mt_single_column )

  ! CODE FOR DMS (if present)

  IF (l_sulpc_dms) THEN

    ! (Note no dry dep for DMS, and emiss added in AERO_CTL via TRSRCE call)

    ! Mix DMS

    CALL tr_mix(                                                        &
    ! IN fields
           bl_levels,alpha_cd, rhokh_mix(:,:,2:), rhokh_mix(:,:,1),     &
           dtrdz_charney_grid, zero_field, zero_field,                  &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
    ! INOUT / OUT fields
           dms(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               1:bl_levels ), tr_flux, drydep_str                       &
             )

  END IF                          ! End If dms condition

END IF                            ! END IF l_sulpc_so2

!---------------------------------------------------------------------
! Mixing for Carbon cycle
!---------------------------------------------------------------------
IF ( l_co2_interactive ) THEN

  DO i = pdims%i_start, pdims%i_end
    DO j = pdims%j_start, pdims%j_end
      co2_flux_tot(i,j)=0.0

      !  (i) CO2 emissions from ancillary file.
      IF (l_co2_emits) THEN
        IF ( co2_emits(i,j)  /=  rmdi ) THEN
          co2_flux_tot(i,j)=co2_flux_tot(i,j) + co2_emits(i,j)
        END IF
      END IF

      !  (ii) CO2 flux from ocean. (+ve implies air to sea)
      IF ( co2flux(i,j)  /=  rmdi ) THEN
        co2_flux_tot(i,j)=co2_flux_tot(i,j) - co2flux(i,j)
      END IF

    END DO
  END DO

  !  (iii) CO2 flux from land processes. (+ve implies biosphere to atmos)
  !        Scale land_co2 by the land fraction before adding to the atmosphere
  !
  !        If TRIFFID is switched on, the CO2 flux from additional processes
  !        (land use and nitrogen limitation) needs to be added to the
  !        atmospheric CO2 tracer. This flux is stored as a prognostic so
  !        will either have been read in from a start dump, or updated by
  !        TRIFFID directly. It has units of kgC/m2/yr, so needs to be
  !        converted to the units of land_co2, kgCO2/m2/s (note that TRIFFID is
  !        constrained to a 360 day calendar).
  !           Although calculated on TRIFFID timesteps only, the flux will then
  !        be applied to the atmosphere, at a constant rate, on every
  !        atmosphere model timestep between TRIFFID calls.
  !

  DO l = 1, land_points

    j = (land_index(l)-1)/pdims%i_end + 1
    i =  land_index(l) - (j-1)*pdims%i_end
    IF (l_triffid) THEN
      land_co2(l) = (resp_s(l) - npp(l)                       +         &
                     triffid_co2_gb(l) / (rsec_per_day * 360.0 ))       &
                                   * m_co2 / m_carbon
    ELSE
      land_co2(l) = (resp_s(l) - npp(l)) * m_co2 / m_carbon
    END IF
    co2_flux_tot(i,j) = co2_flux_tot(i,j) + land_co2(l) * flandg(i,j)
  END DO

  CALL tr_mix(                                                          &
  ! IN fields
           bl_levels,alpha_cd, rhokh_mix(:,:,2:), rhokh_mix(:,:,1),     &
           dtrdz_charney_grid, co2_flux_tot, zero_field,                &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
  ! INOUT / OUT fields
           co2(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,    &
               1:bl_levels ), tr_flux, drydep_str                       &
           )

END IF   ! C-cycle

!------------------------------------------------------------------
! Mixing for Soot variables, including dry deposition
!------------------------------------------------------------------

! Note: Emissions of soot are dealt with in aero_ctl. Here we consider
! just boundary layer mixing and dry deposition of the three modes
! of soot.

IF (l_soot) THEN    ! If soot modelling is included

  !  Calculate resistance factors for species to be dry deposited.

  ! Note that for JULES with coastal tiling, at coastal points
  ! (i.e. those containing both land and sea):
  ! RESIST_B here contains  values approriate for SEA only (i.e. no
  ! meaning including land tiles has been done in SF_EXCH), but
  ! ARESIST and RHO_ARESIST contain grid box mean values calculated
  ! in SF_EXCH taking account of land/sea fractions.

  !   For RESIST_B check to eliminate possible negative values
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF (resist_b(i,j)  <   0.0) THEN
        resist_b(i,j) = 0.0
      END IF
    END DO
  END DO

  !   For STR_RESIST_S values depend on surface type (land, sea,
  !    snow,ice) as well as tracer identity.

  !   Initialise 'stomatal' resistance array to zero:

  DO k = 1, pdims%i_end*pdims%j_end
    resist_s(k)=0.0
  END DO

  ! CODE FOR FRESH SOOT

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) <  1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_freshsoot*resist_b(i,j)
        str_resist_s(i,j) = ress_soot*resist_s(k)

        !  Calculate RES_FACTOR for fresh soot  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for fresh soot to calculate correct RES_FACTOR_LAND
  ! values if land present:
  IF (land_points >  0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_freshsoot, ress_soot,        &
               resist_b_tile, res_factor_land)

    ! Recalculate RES_FACTOR for fresh soot combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix fresh soot

  CALL tr_mix(                                                          &
  ! IN fields
           bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,          &
           dtrdz_charney_grid, zero_field, res_factor,                  &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
  ! INOUT / OUT fields
           soot_new(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    1:bl_levels ), tr_flux, drydep_str                  &
           )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 301                      !dry dep flux fresh soot
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! SF(item,sect)
  END IF  ! ( model_type /= mt_single_column )

  ! CODE FOR AGED SOOT

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) <  1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_agedsoot*resist_b(i,j)
        str_resist_s(i,j) = ress_soot*resist_s(k)

        !  Calculate RES_FACTOR for aged soot  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for aged soot to calculate correct RES_FACTOR_LAND
  ! values if land present:
  IF (land_points >  0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_agedsoot, ress_soot,         &
               resist_b_tile,  res_factor_land)

    ! Recalculate RES_FACTOR for aged soot combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix aged soot

  CALL tr_mix(                                                          &
  ! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,  &
  ! INOUT / OUT fields
         soot_aged(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
                   1:bl_levels ), tr_flux, drydep_str                   &
         )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 302                      !dry dep flux aged soot
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! SF(item,sect)
  END IF  ! ( model_type /= mt_single_column )

  ! CODE FOR CLOUD SOOT

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) <  1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_sootincloud*resist_b(i,j)
        str_resist_s(i,j) = ress_soot*resist_s(k)

        !  Calculate RES_FACTOR for cloud soot  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for cloud soot to calculate correct RES_FACTOR_LAND
  ! values if land present:
  IF (land_points >  0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_sootincloud, ress_soot,      &
               resist_b_tile,  res_factor_land)

    ! Recalculate RES_FACTOR for cloud soot combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix cloud soot

  CALL tr_mix(                                                          &
  ! IN fields
           bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,          &
           dtrdz_charney_grid, zero_field, res_factor,                  &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
  ! INOUT / OUT fields
           soot_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    1:bl_levels ), tr_flux, drydep_str                  &
           )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 303                   !dry (occult) dep flux cloud soot
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! SF(item,sect)
  END IF  ! ( model_type /= mt_single_column )

END IF  ! L_SOOT Soot modelling included

! END of soot section

!------------------------------------------------------------------
! Mixing for Biomass aerosol, including dry deposition
!------------------------------------------------------------------

! Note: Emissions of biomass aerosol are dealt with in aero_ctl.
! Here we consider just boundary layer mixing and dry deposition.

IF (l_biomass) THEN    ! If biomass aerosol is included

  !  Calculate resistance factors for species to be dry deposited.

  ! Note that for JULES with coastal tiling, at coastal points
  ! (i.e. those containing both land and sea):
  ! RESIST_B here contains  values approriate for SEA only (i.e. no
  ! meaning including land tiles has been done in SF_EXCH), but
  ! ARESIST and RHO_ARESIST contain grid box mean values calculated
  ! in SF_EXCH taking account of land/sea fractions.

  !   For RESIST_B check to eliminate possible negative values
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF (resist_b(i,j)  <   0.0) THEN
        resist_b(i,j) = 0.0
      END IF
    END DO
  END DO

  !   For STR_RESIST_S values depend on surface type (land, sea,
  !    snow,ice) as well as tracer identity.

  !   Initialise 'stomatal' resistance array to zero:

  DO k = 1, pdims%i_end*pdims%j_end
    resist_s(k)=0.0
  END DO

  ! CODE FOR FRESH BIOMASS

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) <  1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_freshbmass*resist_b(i,j)
        str_resist_s(i,j) = ress_bmass*resist_s(k)

        ! Calculate RES_FACTOR for fresh biomass  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for fresh biomass to calculate correct RES_FACTOR_LAND
  ! values if land present:
  IF (land_points >  0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_freshbmass, ress_bmass,      &
               resist_b_tile, res_factor_land)

    ! Recalculate RES_FACTOR for fresh biomass combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix fresh biomass smoke

  CALL tr_mix(                                                          &
  ! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,  &
  ! INOUT / OUT fields
         bmass_new(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,&
                   1:bl_levels ), tr_flux, drydep_str                   &
         )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 396                      !dry dep flux fresh biomass
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! SF(item,sect)
  END IF  ! ( model_type /= mt_single_column )

  ! CODE FOR AGED BIOMASS

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) <  1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_agedbmass*resist_b(i,j)
        str_resist_s(i,j) = ress_bmass*resist_s(k)

        ! Calculate RES_FACTOR for aged biomass  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for aged biomass to calculate correct RES_FACTOR_LAND
  ! values if land present:
  IF (land_points >  0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_agedbmass, ress_bmass,       &
               resist_b_tile, res_factor_land)

    ! Recalculate RES_FACTOR for aged smoke combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix aged smoke

  CALL tr_mix(                                                          &
  ! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,  zhnl ,zhsc, z_uv, &
  ! INOUT / OUT fields
         bmass_agd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
                   1:bl_levels ), tr_flux, drydep_str                   &
         )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 397                      !dry dep flux aged biomass
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! SF(item,sect)
  END IF  ! ( model_type /= mt_single_column )

  ! CODE FOR CLOUD BIOMASS

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) <  1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_bmassincloud*resist_b(i,j)
        str_resist_s(i,j) = ress_bmass*resist_s(k)

        ! Calculate RES_FACTOR for cloud smoke  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for cloud smoke to calculate correct RES_FACTOR_LAND
  ! values if land present:
  IF (land_points >  0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_bmassincloud, ress_bmass,    &
               resist_b_tile,  res_factor_land)

    ! Recalculate RES_FACTOR for cloud smoke combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix cloud smoke

  CALL tr_mix(                                                          &
  ! IN fields
         bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,             &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,  &
  ! INOUT / OUT fields
         bmass_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
                   1:bl_levels ), tr_flux, drydep_str                   &
         )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 398                   !dry (occult) dep flux cloud smoke
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! SF(item,sect)
  END IF  ! ( model_type /= mt_single_column )

END IF  ! L_BIOMASS Biomass smoke modelling included

! END of biomass section

!------------------------------------------------------------------
! Mixing for Fossil-fuel organic carbon, including dry deposition
!------------------------------------------------------------------

! Note: Emissions of ocff are dealt with in aero_ctl. Here we consider
! just boundary layer mixing and dry deposition of the three modes
! of ocff.

IF (l_ocff) THEN    ! If ocff modelling is included

  !  Calculate resistance factors for species to be dry deposited.

  ! Note that for JULES with coastal tiling, at coastal points
  ! (i.e. those containing both land and sea):
  ! RESIST_B here contains  values approriate for SEA only (i.e. no
  ! meaning including land tiles has been done in SF_EXCH), but
  ! ARESIST and RHO_ARESIST contain grid box mean values calculated
  ! in SF_EXCH taking account of land/sea fractions.

  !   For RESIST_B check to eliminate possible negative values
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF (resist_b(i,j)  <  0.0) THEN
        resist_b(i,j) = 0.0
      END IF
    END DO
  END DO

  !   For STR_RESIST_S values depend on surface type (land, sea,
  !    snow,ice) as well as tracer identity.

  !   Initialise 'stomatal' resistance array to zero:

  DO k = 1, pdims%i_end*pdims%j_end
    resist_s(k)=0.0
  END DO

  ! CODE FOR FRESH OCFF

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) < 1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_freshocff*resist_b(i,j)
        str_resist_s(i,j) = ress_ocff*resist_s(k)

        !  Calculate RES_FACTOR for fresh ocff  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for fresh ocff to calculate correct RES_FACTOR_LAND
  ! values if land present:
  IF (land_points > 0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_freshocff, ress_ocff,        &
               resist_b_tile,  res_factor_land)

    ! Recalculate RES_FACTOR for fresh ocff combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix fresh ocff

  CALL tr_mix(                                                          &
  ! IN fields
           bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,          &
           dtrdz_charney_grid, zero_field, res_factor,                  &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
  ! INOUT / OUT fields
           ocff_new(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    1:bl_levels ), tr_flux, drydep_str                  &
           )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 407                        !dry dep flux fresh ocff
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! SF(item,sect)
  END IF  ! ( model_type /= mt_single_column )

  ! CODE FOR AGED OCFF

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) < 1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_agedocff*resist_b(i,j)
        str_resist_s(i,j) = ress_ocff*resist_s(k)

        !  Calculate RES_FACTOR for aged ocff  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for aged ocff to calculate correct RES_FACTOR_LAND
  ! values if land present:
  IF (land_points > 0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_agedocff, ress_ocff,         &
               resist_b_tile, res_factor_land)

    ! Recalculate RES_FACTOR for aged ocff combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix aged OCFF

  CALL tr_mix(                                                          &
  ! IN fields
         bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,             &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,  &
  ! INOUT / OUT fields
         ocff_agd(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, &
                  1:bl_levels ), tr_flux, drydep_str                    &
         )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 408                  ! dry dep flux aged ocff
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! SF(item,sect)
  END IF  ! ( model_type /= mt_single_column )

  ! CODE FOR CLOUD OCFF

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) < 1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_ocffincloud*resist_b(i,j)
        str_resist_s(i,j) = ress_ocff*resist_s(k)

        !  Calculate RES_FACTOR for cloud ocff  (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for cloud ocff to calculate correct RES_FACTOR_LAND
  ! values if land present:
  IF (land_points > 0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_ocffincloud, ress_ocff,      &
               resist_b_tile, res_factor_land)

    ! Recalculate RES_FACTOR for cloud ocff combining sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix cloud ocff

  CALL tr_mix(                                                          &
  ! IN fields
         bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,             &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,  &
  ! INOUT / OUT fields
         ocff_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  1:bl_levels ), tr_flux, drydep_str                    &
         )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 409                 !MSR dry (occult) dep flux cloud ocff
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! SF(item,sect)
  END IF  ! ( model_type /= mt_single_column )

END IF  ! L_ocff ocff modelling included

! END of fossil-fuel organic carbon section

!------------------------------------------------------------------
! Mixing for ammonium nitrate, including dry deposition
!------------------------------------------------------------------

! Note:  There are no direct emissions of ammonium nitrate.  Here
! we consider just boundary layer mixing and dry deposition of the
! two modes of ammonium nitrate.

IF (l_nitrate) THEN  ! If ammonium nitrate modelling is included

  !  Calculate resistance factors for species to be dry deposited.

  ! Note that for JULES with coastal tiling, at coastal points
  ! (i.e. those containing both land and sea):
  ! RESIST_B here contains  values approriate for SEA only (i.e. no
  ! meaning including land tiles has been done in SF_EXCH), but
  ! ARESIST and RHO_ARESIST contain grid box mean values calculated
  ! in SF_EXCH taking account of land/sea fractions.

  !   For RESIST_B check to eliminate possible negative values
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF (resist_b(i,j)  <  0.0) THEN
        resist_b(i,j) = 0.0
      END IF
    END DO
  END DO

  !   For STR_RESIST_S values depend on surface type (land, sea,
  !    snow,ice) as well as tracer identity.

  !   Initialise 'stomatal' resistance array to zero:

  DO k = 1, pdims%i_end*pdims%j_end
    resist_s(k)=0.0
  END DO

  ! CODE FOR ACCUMULATION-MODE AMMONIUM NITRATE

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) < 1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_nitr_acc * resist_b(i,j)
        str_resist_s(i,j) = ress_nitr_acc * resist_s(k)

        !  Calculate RES_FACTOR for accumulation nitrate
        ! (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for accumulation nitrate to calculate correct RES_FACTOR_LAND
  ! values if land present:
  IF (land_points > 0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_nitr_acc, ress_nitr_acc,     &
               resist_b_tile,  res_factor_land)

    ! Recalculate RES_FACTOR for accumulation nitrate combining
    ! sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix accumulation nitrate

  CALL tr_mix(                                                          &
  ! IN fields
           bl_levels,alpha_cd,rhokh_mix(:,:,2:), rho_aresist,           &
           dtrdz_charney_grid, zero_field, res_factor,                  &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
  ! INOUT / OUT fields
           nitr_acc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    1:bl_levels ), tr_flux, drydep_str                  &
           )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 274             !dry dep flux accumulation nitrate
    IF ( sf(item,sect) ) THEN
      !    Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! SF(item,sect)
  END IF  ! ( model_type /= mt_single_column )

  ! CODE FOR DISSOLVED AMMONIUM NITRATE

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k=(j-1)*pdims%i_end + i

      IF (flandg(i,j) < 1.0) THEN
        !  For grid boxes containing sea:
        str_resist_b(i,j) = resb_nitr_diss * resist_b(i,j)
        str_resist_s(i,j) = ress_nitr_diss * resist_s(k)

        !  Calculate RES_FACTOR for dissolved nitrate
        ! (only correct for sea points)

        res_factor(i,j) = aresist(i,j) /                                &
           (aresist(i,j)+str_resist_b(i,j)+str_resist_s(i,j))

      END IF

    END DO
  END DO

  ! Call SRESFACT for dissolved nitrate to calculate correct RES_FACTOR_LAND
  ! values if land present:
  IF (land_points > 0) THEN

    CALL sresfact (land_points, land_index,                             &
               ntiles, tile_index, tile_pts, insoluble,                 &
               canopy, catch, gc, tile_frac, snow_tile,                 &
               aresist, aresist_tile, resb_nitr_diss, ress_nitr_diss,   &
               resist_b_tile, res_factor_land)

    ! Recalculate RES_FACTOR for dissolved nitrate combining
    ! sea and land values

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        res_factor(i,j) = (1.0-flandg(i,j))*res_factor(i,j) +           &
                         flandg(i,j)*res_factor_land(i,j)
      END DO
    END DO

  END IF

  ! Mix dissolved nitrate

  CALL tr_mix(                                                          &
  ! IN fields
         bl_levels,alpha_cd, rhokh_mix(:,:,2:), rho_aresist,            &
         dtrdz_charney_grid, zero_field, res_factor,                    &
         kent, we_lim, t_frac, zrzi,                                    &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,  &
  ! INOUT / OUT fields
         nitr_diss(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
                   1:bl_levels ), tr_flux, drydep_str                   &
         )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 275             !dry dep flux dissolved nitrate
    IF ( sf(item,sect) ) THEN
      !       Change sign of dry dep flux (otherwise negative)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          work_2d(i,j) = -drydep_str(i,j)
        END DO
      END DO
      ! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(item,sect,im_index)),work_2d,          &
                    pdims%i_end,pdims%j_end,0,0,0,0,at_extremity,       &
                    atmos_im,sect,item,                                 &
                    errorstatus,cmessage)
      IF (errorstatus /=0) THEN
        CALL ereport( RoutineName, errorstatus, cmessage )
      END IF
    END IF  ! SF(item,sect)
  END IF  ! ( model_type /= mt_single_column )

END IF   ! L_nitrate

! END of ammonium nitrate section.

!---------------------------------------------------------------------
! Mixing for Cariolle Ozone tracer
!---------------------------------------------------------------------
IF ( l_use_cariolle ) THEN

  CALL tr_mix(                                                          &
  ! IN fields
          bl_levels,alpha_cd, rhokh_mix(:,:,2:), rhokh_mix(:,:,1),      &
          dtrdz_charney_grid, zero_field, zero_field,                   &
          kent, we_lim, t_frac, zrzi,                                   &
          kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv, &
  ! INOUT / OUT fields
          ozone_tracer(tdims%i_start:tdims%i_end,                       &
                       tdims%j_start:tdims%j_end,1:bl_levels ),         &
          tr_flux, drydep_str                                           &
          )

  !   Output BL flux diagnostic
  IF ( model_type /= mt_single_column ) THEN
    item = 480
    IF ( sf(item,sect) ) THEN      ! diagnostic flux required
      ! DEPENDS ON: copydiag_3d
      CALL copydiag_3d( stashwork(si(item,sect,im_index)),              &
           tr_flux,                                                     &
           pdims%i_end,pdims%j_end,bl_levels,0,0,0,0, at_extremity,     &
           stlist(1,stindex(1,item,sect,im_index)),len_stlist,          &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,sect,item,                                          &
           errorstatus ,cmessage)
    END IF
    IF (errorstatus /=0) THEN
      CALL ereport( RoutineName, errorstatus, cmessage )
    END IF
  END IF  ! ( model_type /= mt_single_column )

END IF  ! ( l_use_cariolle )

!---------------------------------------------------------------------
! Mixing for Free Tracers
! Will mimic 4.5 here, but note that free tracer levels exist only
! at the top of the model so this may be incorrect if
! model_levels /= tr_levels!
!---------------------------------------------------------------------

IF ( l_bl_tracer_mix ) THEN

  DO n_tracer = 1, tr_vars

    CALL tr_mix(                                                        &
    ! IN fields
           bl_levels,alpha_cd,rhokh_mix(:,:,2:), rhokh_mix(:,:,1),      &
           dtrdz_charney_grid, zero_field, zero_field,                  &
           kent, we_lim, t_frac, zrzi,                                  &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhnl ,zhsc, z_uv,&
    ! INOUT / OUT fields
           free_tracers(tdims%i_start:tdims%i_end,                      &
                        tdims%j_start:tdims%j_end,1:bl_levels,n_tracer),&
           tr_flux, drydep_str                                          &
           )

    !     Output BL flux diagnostic
    IF ( model_type == mt_single_column ) THEN
      WRITE(cdiv,'(I2.2)') n_tracer
      CALL SCMoutput(tr_flux,                                           &
           'tracer_flux'//cdiv,'BL flux of tracer '//cdiv,'m/s',        &
           t_avg,d_bl,default_streams,'',RoutineName)
    ELSE
      item = 99 + n_tracer
      IF ( sf(item,sect) ) THEN      ! diagnostic flux required
        ! DEPENDS ON: copydiag_3d
        CALL copydiag_3d( stashwork(si(item,sect,im_index)),            &
             tr_flux,                                                   &
             pdims%i_end,pdims%j_end,bl_levels,0,0,0,0, at_extremity,   &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             errorstatus, cmessage)
        IF (errorstatus /=0) THEN
          CALL ereport( RoutineName, errorstatus, cmessage )
        END IF
      END IF  ! ( sf(item,sect) )
    END IF  ! model_type

  END DO  ! n_tracer = 1, tr_vars
END IF  ! ( l_bl_tracer_mix )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE bl_trmix_dd
END MODULE bl_trmix_dd_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Top level module for the new UKCA emission system. It replaces
!    ukca_emission_ctl. It contains the subroutine ukca_new_emiss_ctl,
!    which is called from UKCA_MAIN1 at each timestep.
!
!  Method:
!  1) Call ukca_emiss_init to check the NetCDF emission files and
!     initialise the emissions structure (only in the first time step).
!  2) At each time step ukca_emiss_update is called to update
!     emissions (if necessary).
!  3) Look for online emissions (e.g. lightning NOx and CH4 from wetlands).
!  4) Call base_emiss_factors to update the values in the emissions
!     structure (depending on what the emission field is mapped to
!     and what is currently being emitted), and check whether to apply
!     a diurnal cycle to isoprene emissions
!  5) Call ukca_add_emiss to inject emissions and do tracer mixing
!  6) Call ukca_emiss_diags to produce emission diagnostics
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! --------------------------------------------------------------------------
!
MODULE ukca_emiss_ctl_mod

USE stash_array_mod, ONLY: sf

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EMISS_CTL_MOD'

CONTAINS

SUBROUTINE ukca_emiss_ctl (                                             &
  row_length,  rows, bl_levels,                                         &
  global_row_length, global_rows, n_tracers, iyear, imonth, iday, ihour,&
  timestep_number, timestep, l_first, land_points, land_index,          &
  conv_cloud_base, conv_cloud_top,                                      &
  delta_lambda, delta_phi, r_theta_levels, FV_cos_theta_latitude,       &
  cos_zenith_angle, int_zenith_angle, sin_theta_latitude,               &
  tan_theta_latitude, true_longitude, tropopause_height, f3u_by_omega,  &
  land_fraction, r_rho_levels, p_theta_levels, t_theta_levels,          &
  p_layer_boundaries,                                                   &
  totnodens, volume, mass, interf_z, ls_mask,  theta, q, qcl, qcf,      &
  exner_rho_levels, rho_r2, kent, kent_dsc, rhokh_mix, dtrdz_charney_grid,&
  we_lim, t_frac, zrzi, we_lim_dsc, t_frac_dsc, zrzi_dsc, ml_depth,     &
  zhsc, z_half, ch4_wetl_emiss, seaice_frac, area, dust_flux, u_scalar_10m,&
  tstar, dms_sea_conc, chloro_sea, tracers,                             &
  len_stashwork38, stashwork38, len_stashwork50, stashwork50)

USE ukca_emiss_mod,       ONLY: num_em_flds, emissions,                 &
                                ukca_emiss_init, ukca_emiss_update,     &
                                inox_light, ich4_wetl, biogenic_isop,   &
                                ic5h8_ibvoc, ic10h16_ibvoc,             &
                                ich3oh_ibvoc, ich3coch3_ibvoc,          &
                                ukca_emiss_update_mode, ukca_emiss_mode_map,&
                                idms_seaflux, iseasalt_first,           &
                                ipmoc_first, idust_first
USE ukca_emiss_mode_mod,  ONLY: aero_ems_species

USE ukca_emiss_factors,   ONLY: vertical_emiss_factors, base_emiss_factors
USE ukca_add_emiss_mod,   ONLY: ukca_add_emiss
USE ukca_d1_defs,         ONLY: ukca_diag_sect, em_chem_spec
USE ukca_emiss_diags_mod, ONLY: ukca_emiss_diags
USE ukca_emiss_diags_mode_mod, ONLY: ukca_emiss_diags_mode
USE ukca_mode_setup,      ONLY: nmodes, mm, component, cp_cl, cp_du
USE ukca_mode_verbose_mod, ONLY: verbose => glob_verbose
USE ukca_prim_du_mod,     ONLY: ukca_prim_du
USE um_stashcode_mod,     ONLY: stashcode_glomap_sec
USE planet_constants_mod, ONLY: two_omega
USE bvoc_vars,            ONLY: isoprene_gb, terpene_gb, methanol_gb, acetone_gb
USE jules_vegetation_mod, ONLY: l_bvoc_emis
USE conversions_mod,      ONLY: recip_pi_over_180
USE ukca_constants,       ONLY: L_ukca_diurnal_isopems,                 &
                                m_no2, zboltz, ra, avc
USE ukca_option_mod,      ONLY: l_ukca_qch4inter, l_ukca_prescribech4,  &
                                jpctr, l_ukca_mode, l_ukca_aerchem,     &
                                l_ukca_nr_aqchem, l_ukca_ibvoc,         &
                                l_ukca_primss, l_ukca_primdu,           &
                                l_ukca_primbcoc, l_ukca_prim_moc,       &
                                l_ukca_offline, l_ukca_offline_be,      &
                                i_ukca_dms_flux, l_ukca_scale_seadms_ems,&
                                seadms_ems_scaling
USE dms_flux_mod_4A,      ONLY: dms_flux_4A
USE asad_mod,             ONLY: advt
USE asad_chem_flux_diags
USE run_aerosol_mod,      ONLY: l_sulpc_dms
USE ereport_mod,          ONLY: ereport
USE dust_parameters_mod,  ONLY: ndiv
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,           ONLY: umPrint, umMessage, PrintStatus,        &
                                PrStatus_Diag, PrStatus_Normal
USE nlsizes_namelist_mod, ONLY: model_levels
USE parkind1,             ONLY: jpim, jprb      ! DrHook
USE yomhook,              ONLY: lhook, dr_hook  ! DrHook

USE ukca_diurnal_isop_ems_mod, ONLY: ukca_diurnal_isop_ems
USE ukca_light_ctl_mod, ONLY: ukca_light_ctl
USE ukca_prim_ss_mod, ONLY: ukca_prim_ss
USE ukca_prim_moc_mod, ONLY: ukca_prim_moc
IMPLICIT NONE


! Subroutine arguments

! Iput arguments with info on model dimensions
INTEGER,          INTENT(IN)    :: row_length
INTEGER,          INTENT(IN)    :: rows
INTEGER,          INTENT(IN)    :: bl_levels
INTEGER,          INTENT(IN)    :: global_row_length
INTEGER,          INTENT(IN)    :: global_rows

! Input arguments to get nr tracers & current model time
INTEGER, INTENT(IN) :: n_tracers       ! nr tracers
INTEGER, INTENT(IN) :: iyear           ! current yr, mon, day, hr
INTEGER, INTENT(IN) :: imonth
INTEGER, INTENT(IN) :: iday
INTEGER, INTENT(IN) :: ihour
INTEGER, INTENT(IN) :: timestep_number ! nr. timesteps since start

REAL,    INTENT(IN) :: timestep        ! timestep length (sec)
LOGICAL, INTENT(IN) :: l_first         ! T if first call to UKCA

! Input arguments used if iBVOC emissions from JULES
INTEGER, INTENT(IN) :: land_points
INTEGER, INTENT(IN) :: land_index(land_points)

! Input arguments to apply diurnal cycle to isoprene emissions,
! get SO2 emiss from explosive volcanic eruptions, calculate
! lightning emissions of NOx, get vertical profiles or calculate
! surf_area
INTEGER, INTENT(IN) :: conv_cloud_base    (1:row_length, 1:rows)
INTEGER, INTENT(IN) :: conv_cloud_top     (1:row_length, 1:rows)

REAL, INTENT(IN) :: delta_lambda
REAL, INTENT(IN) :: delta_phi
REAL, INTENT(IN) :: r_theta_levels        (1:row_length, 1:rows, 0:model_levels)
REAL, INTENT(IN) :: FV_cos_theta_latitude (1:row_length, 1:rows)
REAL, INTENT(IN) :: cos_zenith_angle      (1:row_length, 1:rows)
REAL, INTENT(IN) :: int_zenith_angle      (1:row_length, 1:rows)
REAL, INTENT(IN) :: sin_theta_latitude    (1:row_length, 1:rows)
REAL, INTENT(IN) :: tan_theta_latitude    (1:row_length, 1:rows)
REAL, INTENT(IN) :: true_longitude        (1:row_length, 1:rows)
REAL, INTENT(IN) :: tropopause_height     (1:row_length, 1:rows)
REAL, INTENT(IN) :: f3u_by_omega          (1:row_length, 1:rows)
REAL, INTENT(IN) :: land_fraction         (1:row_length, 1:rows)
REAL, INTENT(IN) :: r_rho_levels          (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: p_theta_levels        (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: p_layer_boundaries    (1:row_length, 1:rows, 0:model_levels)
REAL, INTENT(IN) :: totnodens             (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: volume                (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: mass                  (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: interf_z              (1:row_length, 1:rows, 0:model_levels)

LOGICAL, INTENT(IN) :: ls_mask            (1:row_length, 1:rows)

! Input arguments needed by UKCA_ADD_EMISS to call TSRCE
REAL, INTENT(IN) :: theta             (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: q                 (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: qcl               (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: qcf               (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: exner_rho_levels  (1:row_length, 1:rows, 1:model_levels+1)
REAL, INTENT(IN) :: rho_r2            (1:row_length, 1:rows, 1:model_levels)

! Input arguments needed by UKCA_ADD_EMISS to call TR_MIX
INTEGER, INTENT(IN) :: kent       (1:row_length, 1:rows)
INTEGER, INTENT(IN) :: kent_dsc   (1:row_length, 1:rows)

! Input arguments needed by UKCA_ADD_EMISS to call UKCA_PRIM_SS, UKCA_PRIM_MOC,
! UKCA_PRIM_DU and DMS_FLUX_4A.
REAL, INTENT(IN) :: seaice_frac    (1:row_length, 1:rows)
REAL, INTENT(IN) :: dust_flux(row_length,rows,ndiv) ! dust emiss (kg m-2 s-1)
REAL, INTENT(IN) :: u_scalar_10m   (:,:)
REAL, INTENT(IN) :: t_theta_levels (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: tstar          (1:row_length, 1:rows)
REAL, INTENT(IN) :: dms_sea_conc   (1:row_length, 1:rows)
REAL, INTENT(IN) :: chloro_sea     (1:row_length, 1:rows) ! surface chl-a kg/m3

! Area of grid cell, required by ukca_emiss_diags_mode
REAL, INTENT(IN) :: area(row_length, rows, model_levels)

REAL, INTENT(IN) :: rhokh_mix          (1:row_length, 1:rows, 1:bl_levels)
REAL, INTENT(IN) :: dtrdz_charney_grid (1:row_length, 1:rows, 1:bl_levels)

REAL, INTENT(IN) :: we_lim     (1:row_length, 1:rows, 1:3)
REAL, INTENT(IN) :: t_frac     (1:row_length, 1:rows, 1:3)
REAL, INTENT(IN) :: zrzi       (1:row_length, 1:rows, 1:3)

REAL, INTENT(IN) :: we_lim_dsc (1:row_length, 1:rows, 1:3)
REAL, INTENT(IN) :: t_frac_dsc (1:row_length, 1:rows, 1:3)
REAL, INTENT(IN) :: zrzi_dsc   (1:row_length, 1:rows, 1:3)

REAL, INTENT(IN) :: ml_depth   (1:row_length, 1:rows)
REAL, INTENT(IN) :: zhsc       (1:row_length, 1:rows)
REAL, INTENT(IN) :: z_half     (1:row_length, 1:rows,1:bl_levels)

! Length of diagnostics arrays
INTEGER, INTENT(IN) :: len_stashwork38, len_stashwork50

! Wetland emissions of methane - scaled by land_fraction while adding
! to emissions structure. Any negative values are removed.
REAL, INTENT(INOUT) :: ch4_wetl_emiss (1:row_length,1:rows)

! Tracer mass mixing ratios
REAL, INTENT(INOUT) :: tracers (1:row_length, 1:rows, 1:model_levels,   &
                                1:n_tracers)
! Diagnostics arrays
REAL, INTENT(INOUT) :: stashwork38(len_stashwork38)
REAL, INTENT(INOUT) :: stashwork50(len_stashwork50)

! Local variables
INTEGER                 :: i, j, k, l, n, ilev
INTEGER                 :: ii, jj, ll        ! loop counters for iBVOC mapping
INTEGER                 :: section           ! stash section
REAL                    :: base_scaling      ! factor to convert emiss units
REAL                    :: vert_fact_3d (1:row_length, 1:rows, 1:model_levels)
                                     ! 3-D vertical profiles of emissions
REAL                    :: aer_emmas (row_length, rows, model_levels, nmodes)
                                     ! aerosol mass emissions by mode
REAL                    :: aer_emnum (row_length, rows, model_levels, nmodes)
                                     ! aerosol number emissions by mode
REAL                    :: chla(row_length, rows)
                                     ! chlorophyll-a, mg.m-3
REAL                    :: mass_pmoc (row_length, rows, 1) 
                                     ! primary marine organic carbon emission
REAL                    :: dummy (row_length, rows, model_levels, nmodes)
                                          ! for interface to ukca_primdu
INTEGER                 :: emiss_levs     ! num levels in emiss field
INTEGER                 :: icp            ! component
LOGICAL                 :: lmode_emiss(nmodes)   ! modes which are emitted into
                                                 ! for a given emission field
REAL                    :: f_dms_sea(row_length, rows)  ! DMS flux over sea

REAL, SAVE, ALLOCATABLE :: surf_area      (:,:) ! gridbox area

LOGICAL                 :: gridbox_emiss        ! True for emissions reported
                                                ! per grid box (e.g. lightning)

REAL :: isoprene_2D (row_length, rows)       ! 2D isoprene ems if iBVOC emiss
REAL :: terpene_2D  (row_length, rows)       ! 2D terpene ems
REAL :: methanol_2D (row_length, rows)       ! 2D methanol ems
REAL :: acetone_2D  (row_length, rows)       ! 2D acetone ems

REAL :: tmp_in_em_field  (row_length, rows)  ! Climatological and daily varying
REAL :: tmp_out_em_field (row_length, rows)  ! isoprene emission field
LOGICAL, SAVE    :: testdcycl = .FALSE.      ! True only for debugging
                                             ! UKCA_DIURNAL_ISOP_EMS

LOGICAL  :: l_3dsource = .FALSE.             ! source emissions are 3D

INTEGER, SAVE, ALLOCATABLE :: ils_mask (:,:)  ! Land/sea mask (1/0) used
                                              ! for lightning NOx

! Arrays of NOx lightning emiss with different units:
! kg(N)/grid box/s and kg(NO2)/kg(air)/s
REAL :: lightningem_n_gridbox  (1:row_length, 1:rows, 1:model_levels)
REAL :: lightningem_no2_to_air (1:row_length, 1:rows, 1:model_levels)

! Number density of air (per cm3)
REAL :: aird(row_length,rows,model_levels)

!  Molar mass of dry air (kg/mol) =avc*zboltz/ra
REAL    :: mm_da

INTEGER                            :: errcode    ! Error code for ereport
CHARACTER (LEN=errormessagelength) :: cmessage   ! Error message
INTEGER            :: ierr       ! Req. for ASAD diags (do not initialise here)
INTEGER, SAVE      :: inox       ! Index of NOx tracer for ASAD 3D diags

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_CTL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialisation
errcode  = 0

isoprene_2D (:,:) = 0.0      ! 2D isoprene interactive emiss
terpene_2D  (:,:) = 0.0      ! 2D terpene  interactive emiss
methanol_2D (:,:) = 0.0      ! 2D methanol interactive emiss
acetone_2D  (:,:) = 0.0      ! 2D acetone  interactive emiss

IF ( l_first ) THEN

  IF (l_ukca_ibvoc .AND. (.NOT. l_bvoc_emis)) THEN
    errcode  = 1
    cmessage = 'UKCA cannot use iBVOC emissions because they are not' //&
               ' activated in JULES'
    CALL ereport (ModuleName//':'//RoutineName, errcode, cmessage)
  END IF

  ! Calculate the gridbox surface area (in m^2). This will be used to convert
  ! any emissions given as 'kg gridcell-1 s-1' to 'kg m-2 s-1'
  ALLOCATE (surf_area (row_length, rows))
  DO j = 1, rows
    DO i = 1, row_length
      surf_area(i,j) = r_theta_levels (i, j, 0)                         &
                     * r_theta_levels (i, j, 0)                         &
                     * delta_lambda * delta_phi                         &
                     * FV_cos_theta_latitude (i, j)
    END DO
  END DO

  ! Read emission NetCDF files and look for the the emission fields
  ! in them to allocate all variables in the emissions structure.
  CALL ukca_emiss_init (row_length,        rows,          model_levels, &
                        global_row_length, global_rows)

  ! Get index for NOx tracer (needed only if ASAD diagnostics
  ! for NOx lightning emiss are required)
  DO k = 1, jpctr
    SELECT CASE (advt(k))
    CASE ('NOx       ')
      inox = k
    CASE ('NO        ')
      inox = k
    END SELECT
  END DO

  IF (inox == -99 .AND. .NOT. (l_ukca_offline .OR. l_ukca_offline_be)) THEN
    errcode  = 1
    cmessage = 'Did not find NO or NOx tracer'
    CALL ereport (ModuleName//':'//RoutineName, errcode, cmessage)
  END IF

  ! Set up integer land/sea mask (needed for lightning emissions)
  ALLOCATE (ils_mask (row_length, rows))
  DO j = 1, rows
    DO i = 1, row_length
      IF (ls_mask (i, j)) THEN
        ils_mask (i, j) = 1
      ELSE
        ils_mask (i, j) = 0
      END IF
    END DO
  END DO

  ! Get scaling factors indicating how to to spread emissions
  ! over different vertical levels. Done it only once here and
  ! stored in the emissions structure, then valid for all time steps.
  DO l = 1, num_em_flds
    CALL vertical_emiss_factors (row_length, rows, model_levels,        &
            emissions(l)%lowest_lev, emissions(l)%highest_lev,          &
            interf_z (1:row_length, 1:rows, 0:model_levels),            &
            emissions(l)%vert_fact,  emissions(l)%var_name,             &
            vert_fact_3d (1:row_length, 1:rows, 1:model_levels))

    ALLOCATE (emissions(l)%vert_scaling_3d (1:row_length, 1:rows,       &
                                            1:model_levels))

    emissions(l)%vert_scaling_3d (:,:,:) = vert_fact_3d (:,:,:)
  END DO

END IF  ! IF l_first


! Check if it is time to update the emiss fields that
! we read from NetCDF files (depending on time step
! and update frequency). If needed then update the files.

CALL ukca_emiss_update (row_length, rows, model_levels,                 &
                        global_row_length, global_rows,                 &
                        timestep_number, timestep,                      &
                        l_first)

! ----------------------------------------------------------------------
! Deal with online emissions, which are not read from NetCDF files and
! are always updated at each time step. For the moment only NO2 from
! lightning, CH4 from wetlands and sea-air DMS flux are considered here.
!
IF (.NOT. (l_ukca_offline .OR. l_ukca_offline_be)) THEN
  ! Call the routine to diagnose NO2 lightning emissions.
  lightningem_n_gridbox  = 0.0
  lightningem_no2_to_air = 0.0

  CALL ukca_light_ctl (                                                 &
        rows, row_length, delta_lambda, delta_phi, model_levels,        &
        conv_cloud_base    (1:row_length, 1:rows),                      &
        conv_cloud_top     (1:row_length, 1:rows),                      &
        ils_mask           (1:row_length, 1:rows),                      &
        recip_pi_over_180 *                                        &
             ASIN (f3u_by_omega (1:row_length, 1:rows) ),               &
        surf_area          (1:row_length, 1:rows),                      &
        r_theta_levels     (1:row_length, 1:rows, 0:model_levels),      &
        r_rho_levels       (1:row_length, 1:rows, 1:model_levels),      &
        p_theta_levels     (1:row_length, 1:rows, 1:model_levels),      &
        p_layer_boundaries (1:row_length, 1:rows, 0:model_levels),      &
        lightningem_n_gridbox  (1:row_length, 1:rows, 1:model_levels),  &
        lightningem_no2_to_air (1:row_length, 1:rows, 1:model_levels))

  ! ASAD diagnostics for lightning NOx emissions, in kg(NO2)/kg(air)/s. Note:
  !  lightningem_no2_to_air ==> lightning NOx emissions, in kg(NO2)/kg(air)/s
  !  lightning_emissions    ==> type flag used in ASAD diagnostic routines
  IF (L_asad_use_chem_diags .AND. L_asad_use_light_ems) THEN
    CALL asad_3D_emissions_diagnostics (                                &
                row_length,                                             &
                rows,                                                   &
                model_levels,                                           &
                inox,                                                   &
                lightningem_no2_to_air,                                 &
                surf_area,                                              &
                totnodens,                                              &
                volume,                                                 &
                mass,                                                   &
                m_no2,                                                  &
                timestep,                                               &
                lightning_emissions,                                    &
                ierr)
  END IF

  ! Fill in values of lightning NOx in the emissions structure
  ! and indicate that the emission field is updated.
  ! For this we need to use lightningem_n_gridbox, because units are
  ! defined as 'kg(N) gridbox-1 s-1' in the subroutine ukca_emiss_init.
  emissions(inox_light)%values (:,:,:) = lightningem_n_gridbox (:,:,:) ! 3D
  emissions(inox_light)%l_update       = .TRUE.

END IF  ! .NOT. l_ukca_offline

! Set wetland methane emissions to zero over non-land surfaces
WHERE (ch4_wetl_emiss < 0.0) ch4_wetl_emiss = 0.0

! Fill in CH4 wetland emissions (if used) at surface level in the
! emissions structure and indicate that they are updated.
IF (L_ukca_qch4inter .AND. (.NOT. L_ukca_prescribech4)) THEN
  ! Scale wetland emission rates by the land fraction
  emissions(ich4_wetl)%values (:,:,1) = ch4_wetl_emiss (:,:) *          &
                                         land_fraction(:,:)
  emissions(ich4_wetl)%l_update       = .TRUE.
END IF

! Interactive BVOC emissions
IF (l_ukca_ibvoc .AND. l_bvoc_emis) THEN
  ! Regrid iBVOC emissions from landpoints (in JULES) to 2D-grid.
  DO ll = 1, land_points
    jj = (land_index (ll) - 1)     / row_length + 1
    ii =  land_index (ll) - (jj-1) * row_length

    isoprene_2D (ii,jj) = isoprene_gb (ll)
    terpene_2D  (ii,jj) = terpene_gb  (ll)
    methanol_2D (ii,jj) = methanol_gb (ll)
    acetone_2D  (ii,jj) = acetone_gb  (ll)
  END DO

  ! Apply coastal tile correction to BVOC emission fluxes
  isoprene_2D (:,:) = isoprene_2D (:,:) * land_fraction(:,:)
  terpene_2D  (:,:) = terpene_2D  (:,:) * land_fraction(:,:)
  methanol_2D (:,:) = methanol_2D (:,:) * land_fraction(:,:)
  acetone_2D  (:,:) = acetone_2D  (:,:) * land_fraction(:,:)

  IF (ic5h8_ibvoc > 0) THEN
    emissions(ic5h8_ibvoc)%values (:,:,1) = isoprene_2D (:,:)  ! surface emiss
    emissions(ic5h8_ibvoc)%l_update       = .TRUE.
  END IF

  IF (ic10h16_ibvoc > 0) THEN
    emissions(ic10h16_ibvoc)%values (:,:,1) = terpene_2D (:,:)
    emissions(ic10h16_ibvoc)%l_update       = .TRUE.
  END IF

  IF (ich3oh_ibvoc > 0) THEN
    emissions(ich3oh_ibvoc)%values (:,:,1) = methanol_2D  (:,:)
    emissions(ich3oh_ibvoc)%l_update       = .TRUE.

  END IF

  IF (ich3coch3_ibvoc > 0) THEN
    emissions(ich3coch3_ibvoc)%values (:,:,1) = acetone_2D (:,:)
    emissions(ich3coch3_ibvoc)%l_update       = .TRUE.
  END IF
END IF


! Call dms flux routine if CLASSIC OFF
IF (ANY(em_chem_spec == 'DMS     ') .AND. (.NOT. l_sulpc_dms)) THEN
  CALL dms_flux_4A(row_length, rows,                                    &
           u_scalar_10m, tstar, land_fraction, dms_sea_conc,            &
           i_ukca_dms_flux, f_dms_sea)

  IF (l_ukca_scale_seadms_ems) THEN
    ! Multiplication by sea DMS emission factor.
    f_dms_sea = seadms_ems_scaling*f_dms_sea
  END IF

  emissions(idms_seaflux)%values(:,:,1) = f_dms_sea(:,:) *              &
        (1.0 - land_fraction(:,:)) * (1.0 - seaice_frac(:,:))
  emissions(idms_seaflux)%l_update = .TRUE.
END IF     ! If 'DMS  ' & CLASSIC OFF

! ----------------------------------------------------------------------
! Update all fields in the emissions super array:
! * Do conversions so that emissions are given as 'kg(tracer) m-2 s-1'
! * Update isoprene emissions if they are diurnally varying
DO l = 1, num_em_flds
  ! Only apply conversion factors to each emission field if it
  ! is updated in the current time step
  ! Note that l_update is always false for online emissions, incl mode ems
  IF (emissions(l)%l_update) THEN
    ! Get conversion factors for the emission field. A base_scaling
    ! different from 1 will be applied if what is being emitted does
    ! not exactly coincide with the tracer the emission field is
    ! mapped to (e.g. if an NO emission field is expressed as kg of
    ! nitrogen, or if a VOC emission is expressed as kg of carbon).
    CALL base_emiss_factors (emissions(l)%tracer_name,                  &
                             emissions(l)%units,                        &
                             emissions(l)%var_name,                     &
                             emissions(l)%std_name,                     &
                             emissions(l)%lng_name,                     &
                             base_scaling, gridbox_emiss)

    ! Get emissions in kg of emitted tracer
    emissions(l)%values (:,:,:) = emissions(l)%values (:,:,:) * base_scaling

    ! Debug update times only for first emission field
    IF (l ==1 .AND. PrintStatus >= PrStatus_Diag) THEN
      WRITE (umMessage, '(A,1x,I6)')                                    &
         RoutineName//': Update first emission field for timestep ',&
         timestep_number
      CALL umPrint(umMessage,src=ModuleName//':'//RoutineName)
    END IF

    ! If emissions are in kg gridbox-1 s-1 then convert to kg m-2 s-1
    IF (gridbox_emiss) THEN
      IF (emissions(l)%three_dim) THEN      ! 3D emiss
        DO ilev = 1, model_levels
          emissions(l)%values (:,:,ilev) =                              &
          emissions(l)%values (:,:,ilev) / surf_area (:,:)
        END DO
      ELSE                                  ! surface or single-level emiss
        emissions(l)%values (:,:,1) =                                   &
        emissions(l)%values (:,:,1) / surf_area (:,:)
      END IF
    END IF

    ! In the case of biogenic emissions of isoprene, need to store the
    ! latest emission values updated in another variable to avoid re-
    ! applying the hourly scaling with the call to UKCA_DIURNAL_ISOP_EMS
    IF (emissions(l)%tracer_name == 'C5H8      '      .AND.             &
        emissions(l)%hourly_fact == 'diurnal_isopems' .AND.             &
        L_ukca_diurnal_isopems                        .AND.             &
        .NOT. L_ukca_ibvoc) THEN
      biogenic_isop = emissions(l)%values (:,:,1)  ! always 2D field
    END IF

  END IF  ! IF l_update ?


  ! If L_ukca_diurnal_isopems is set to true then
  ! use a diurnally varying isoprene emission field.
  ! This is only applied to biogenically emitted isoprene
  ! (currently with houly scaling = 'diurnal_isopems')
  ! but not to isoprene emiss from other sources.
  IF (emissions(l)%tracer_name == 'C5H8      '      .AND.               &
      emissions(l)%hourly_fact == 'diurnal_isopems' .AND.               &
      L_ukca_diurnal_isopems                        .AND.               &
      .NOT. L_ukca_ibvoc) THEN

    tmp_in_em_field (:,:) = biogenic_isop (:,:)

    CALL ukca_diurnal_isop_ems (row_length, rows,                       &
         tmp_in_em_field,  cos_zenith_angle,                            &
         int_zenith_angle,                                              &
         sin_theta_latitude, FV_cos_theta_latitude,                     &
         tan_theta_latitude, timestep, tmp_out_em_field,                &
         testdcycl)

    ! Update the emission field
    emissions(l)%values(:,:,1) = tmp_out_em_field (:,:)
  END IF

  ! Some basic checks for first time step if print extra diagnostics
  IF (l_first .AND. PrintStatus > PrStatus_Normal) THEN
    WRITE  (umMessage,'(A,I3,1X,A,1X,A)') RoutineName//':', l,    &
      TRIM (emissions(l)%var_name), TRIM (emissions(l)%tracer_name)
    CALL umPrint(umMessage,src=RoutineName)
    WRITE (umMessage,'(A,I3,1X,A)') RoutineName//': ', l,          &
      TRIM (emissions(l)%std_name)
    CALL umPrint(umMessage,src=RoutineName)
    WRITE (umMessage,'(A,I3,1X,A)') RoutineName//': ', l,          &
      TRIM (emissions(l)%lng_name)
    CALL umPrint(umMessage,src=RoutineName)
    WRITE (umMessage,'(A,I3,1X,A,1X,L7,1X,F10.5)') RoutineName//': ', l,&
      TRIM (emissions(l)%units), gridbox_emiss, base_scaling
    CALL umPrint(umMessage,src=RoutineName)
    WRITE (umMessage,'(A,I3,1X,A,1X,I4,1X,A,3X,I1)') RoutineName//': ', l,&
      'update_freq:', emissions(l)%update_freq,                         &
      'update_type:', emissions(l)%update_type
    CALL umPrint(umMessage,src=RoutineName)
    WRITE (umMessage,'(A)')  RoutineName//': ---------------------------'
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  ! Uncomment lines below for short simulations if needed:
  ! Check whether each emission value is updated for each timestep.
  ! Result should always be T for online emissions, and T for NetCDF
  ! emissions only if model time is multiple of the update frequency
  ! for the given emission field.
  ! IF (PrintStatus >= PrStatus_Normal) THEN
  !   WRITE (umMessage,'(A,I3,X,A,X,A,X,I6,A,X,L)') 'UKCA_NEW_EMISS_CTL - ', l,&
  !     TRIM (emissions(l)%var_name), 'updated for tstep', timestep_number, &
  !     "?", emissions(l)%l_update
  !   CALL umPrint(umMessage,src='ukca_new_emiss_ctl')
  ! END IF

END DO  ! loop through number of emission fields

! --------------------------------
! Aerosol emissions into modes:
! Online emissions first (sea salt, dust) which return number and mass for each
! mode, and these are mapped onto a seperate emissions structure for each
! mode/moment.
! Then offline emissions, which are provided as total mass and are projected
! onto number and mass for each mode by ukca_emiss_update_mode, before being
! mapped to emissions structures for each mode/moment in the same way as the
! online emissions. This has to be done after the scaling factors are applied
! to the offline emission structures above.

! Online emissions for sea salt
IF (l_ukca_primss) THEN
  aer_emmas(:,:,:,:) = 0.0
  aer_emnum(:,:,:,:) = 0.0
  aird(:,:,:) = p_theta_levels(:,:,:)/(t_theta_levels(:,:,:)*           &
                zboltz*1e6)           ! No. density (/cm^3)
  CALL ukca_prim_ss(row_length, rows, model_levels,                     &
                    verbose, land_fraction, seaice_frac, mass, u_scalar_10m,&
                    aird, aer_emmas, aer_emnum)

  ! Convert mass emission to kg(component)/m2/s
  mm_da = avc*zboltz/ra  ! Molar mass of dry air (kg/mol)
  aer_emmas = aer_emmas * mm(cp_cl) / mm_da
  lmode_emiss(:) = component(:,cp_cl)  ! Seasalt emiss covers all ss modes.
  emiss_levs = 1                       ! Surface emissions.

  ! Map num and mass arrays from ukca_prim_ss to corresponding emissions
  ! structures.
  CALL ukca_emiss_mode_map(row_length, rows, model_levels,              &
                           aer_emmas, aer_emnum, emiss_levs, cp_cl,     &
                           lmode_emiss, iseasalt_first)
END IF

! Online emissions for primary marine organic carbon.
! Requires sea-salt ems to be turned on
IF (l_ukca_prim_moc .AND. l_ukca_primbcoc .AND. l_ukca_primss) THEN
! Use the seasalt emissions and chlorophyll concentration to calculate
! marine organic primary emissions. aer_emass/num are the input seasalt
! emissions for each mode, while mass_pmoc is the total mass of marine
! organic emission. This is then distributed across modes by
! ukca_emission_update_mode, and mapped into the emissions super-array by
! ukca_emiss_mode_map. It is done this way to take advantage of
! ukca_emiss_update_mode to do the vol/number calculation and use the
! frac/size configuration information in ukca_emiss_mode_mod.

  ! Note that chloro_sea is not masked and contains land values.
  ! Land masking is accounted for implicitly here as a result of using
  ! sea-salt emissions as input to ukca_prim_moc
  chla(:,:) = chloro_sea(:,:) * 1e6 ! chlorophyll-a concentration, mg.m-3

  CALL ukca_prim_moc(row_length, rows, model_levels,                    &
                      aer_emmas, aer_emnum, chla, u_scalar_10m,         &
                      mass_pmoc)

  ! re-initialisation to zero after calculation of PMOC flux
  aer_emmas(:,:,:,:) = 0.0
  aer_emnum(:,:,:,:) = 0.0

  CALL ukca_emiss_update_mode(row_length, rows, model_levels,           &
                              "PMOC      ", l_3dsource, mass_pmoc,      &
                              aer_emmas, aer_emnum, lmode_emiss,        &
                              emiss_levs,icp)

  CALL ukca_emiss_mode_map(row_length, rows, model_levels,              &
                           aer_emmas, aer_emnum, emiss_levs, icp,       &
                           lmode_emiss, ipmoc_first)

ELSE IF (l_ukca_prim_moc .AND. (.NOT. l_ukca_primss)) THEN
  cmessage = 'Prim. marine OC requires sea-salt emissions to be turned on'
  errcode = 1
  CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)

END IF

! Online emissions for dust
IF (l_ukca_primdu) THEN
  aer_emmas(:,:,:,:) = 0.0
  aer_emnum(:,:,:,:) = 0.0
  aird(:,:,:) = p_theta_levels(:,:,:)/(t_theta_levels(:,:,:)*           &
                zboltz*1e6)           ! No. density (/cm^3)
  CALL ukca_prim_du(row_length, rows, model_levels,                     &
                    verbose, dust_flux,                                 &
                    aer_emnum, dummy, aer_emmas)
  lmode_emiss(:) = component(:,cp_du)  ! Dust emiss covers all du modes.
  emiss_levs = 1                       ! Surface emissions.

  ! Map num and mass arrays from ukca_prim_du to corresponding emissions
  ! structures.
  CALL ukca_emiss_mode_map(row_length, rows, model_levels,              &
                           aer_emmas, aer_emnum, emiss_levs, cp_du,     &
                           lmode_emiss, idust_first)
END IF ! l_ukca_primdu


! Map offline aerosol emissions into mode emissions structures (see block
! comment above).
DO l = 1, num_em_flds
  IF (ANY(aero_ems_species == emissions(l)%tracer_name) .AND. l_ukca_mode) THEN
    aer_emmas(:,:,:,:) = 0.0
    aer_emnum(:,:,:,:) = 0.0
    CALL ukca_emiss_update_mode(row_length, rows, model_levels, &
                                emissions(l)%tracer_name, &
                                emissions(l)%three_dim, &
                                emissions(l)%values(:,:,:), &
                                aer_emmas, aer_emnum, lmode_emiss, emiss_levs,&
                                icp)
    CALL ukca_emiss_mode_map(row_length, rows, model_levels,            &
                             aer_emmas, aer_emnum, emiss_levs, icp,     &
                             lmode_emiss, l)
  END IF
END DO
! End of aerosol emission into modes
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Inject emissions and do tracer mixing
CALL ukca_add_emiss (                                                   &
  row_length, rows, bl_levels,                                          &
  n_tracers, iyear, imonth, iday, ihour, timestep,                      &
  delta_lambda, delta_phi,                                              &
  sin_theta_latitude    (1:row_length, 1:rows),                         &
  true_longitude        (1:row_length, 1:rows),                         &
  tropopause_height     (1:row_length, 1:rows),                         &
  r_theta_levels        (1:row_length, 1:rows, 0:model_levels),         &
  theta                 (1:row_length, 1:rows, 1:model_levels),         &
  q                     (1:row_length, 1:rows, 1:model_levels),         &
  qcl                   (1:row_length, 1:rows, 1:model_levels),         &
  qcf                   (1:row_length, 1:rows, 1:model_levels),         &
  exner_rho_levels      (1:row_length, 1:rows, 1:model_levels+1),       &
  rho_r2                (1:row_length, 1:rows, 1:model_levels),         &
  kent                  (1:row_length, 1:rows),                         &
  kent_dsc              (1:row_length, 1:rows),                         &
  rhokh_mix             (1:row_length, 1:rows, 1:bl_levels),            &
  dtrdz_charney_grid    (1:row_length, 1:rows, 1:bl_levels),            &
  we_lim                (1:row_length, 1:rows, 1:3),                    &
  t_frac                (1:row_length, 1:rows, 1:3),                    &
  zrzi                  (1:row_length, 1:rows, 1:3),                    &
  we_lim_dsc            (1:row_length, 1:rows, 1:3),                    &
  t_frac_dsc            (1:row_length, 1:rows, 1:3),                    &
  zrzi_dsc              (1:row_length, 1:rows, 1:3),                    &
  ml_depth              (1:row_length, 1:rows),                         &
  zhsc                  (1:row_length, 1:rows),                         &
  z_half                (1:row_length, 1:rows, 1:bl_levels),            &
  surf_area             (1:row_length, 1:rows),                         &
  totnodens             (1:row_length, 1:rows, 1:model_levels),         &
  volume                (1:row_length, 1:rows, 1:model_levels),         &
  mass                  (1:row_length, 1:rows, 1:model_levels),         &
  tracers               (1:row_length, 1:rows, 1:model_levels,          &
                         1:n_tracers),                                  &
  len_stashwork50, stashwork50)

! ----------------------------------------------------------------------
! Call the emission diagnostics code if any of the diagnostics present
! in the routine GET_EMDIAG_STASH has been selected via stash.
! Chemical diagnostics first, then aerosol diagnostics.
section = UKCA_diag_sect

IF (sf(156,section) .OR. sf(157,section) .OR. sf(158,section) .OR.      &
    sf(159,section) .OR. sf(160,section) .OR. sf(161,section) .OR.      &
    sf(162,section) .OR. sf(163,section) .OR. sf(164,section) .OR.      &
    sf(165,section) .OR. sf(166,section) .OR. sf(167,section) .OR.      &
    sf(168,section) .OR. sf(169,section) .OR. sf(170,section) .OR.      &
    sf(171,section) .OR. sf(172,section) .OR. sf(211,section) .OR.      &
    sf(212,section) .OR. sf(213,section) .OR. sf(214,section) .OR.      &
    sf(215,section) .OR. sf(216,section) .OR. sf(217,section) ) THEN

  CALL ukca_emiss_diags (row_length, rows, model_levels,                &
    len_stashwork50, stashwork50)
END IF


! Aerosol diagnostics
section = stashcode_glomap_sec

IF (sf(201,section) .OR. sf(202,section) .OR. sf(203,section) .OR.      &
    sf(204,section) .OR. sf(205,section) .OR. sf(206,section) .OR.      &
    sf(207,section) .OR. sf(208,section) .OR. sf(209,section) .OR.      &
    sf(210,section) .OR. sf(211,section) .OR. sf(212,section) .OR.      &
    sf(213,section)  ) THEN

  CALL ukca_emiss_diags_mode (row_length, rows, model_levels, area,     &
    len_stashwork38, stashwork38)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_emiss_ctl

END MODULE ukca_emiss_ctl_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE aero_ctl_mod_4A

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='AERO_CTL_MOD_4A'

CONTAINS

SUBROUTINE aero_ctl_4A(                                                        &
  ! Parallel variables
  halo_i, halo_j, off_x, off_y, global_row_length, global_rows,                &
  proc_row_group, proc_col_group, at_extremity, n_proc, n_procx,               &
  n_procy, neighbour, g_rows, g_row_length, me,                                &
  ! model dimensions
  row_length, rows, n_rows, land_points,                                       &
  bl_levels, n_cca_levels,                                                     &
  theta_field_size,                                                            &
  salt_dim1, salt_dim2, salt_dim3,                                             &
  aero_dim1, aero_dim2, aero_dim3,                                             &
  ! Co-ordinate information
  delta_lambda, delta_phi,                                                     &
  lat_rot_np, long_rot_np,                                                     &
  ! Time stepping information
  val_year, val_month, val_day_number, val_day, val_hour, val_minute,          &
  val_second, timestep_number,                                                 &
  previous_time,                                                               &
  call_chem_freq,                                                              &
  ! Trig arrays
  sin_theta_longitude, cos_theta_longitude,                                    &
  fv_cos_theta_latitude,                                                       &
  ! Grid-dependent arrays
  f3_at_u, true_longitude, true_latitude,                                      &
  ! Data fields IN
  u, v, tstar, tstar_sea,                                                      &
  theta, q, qcl, qcf,                                                          &
  rho, land_mask, fland_ctile,                                                 &
  p_theta_levels, exner_rho_levels, exner_theta_levels,                        &
  ice_fract, snow_depth,                                                       &
  cloudf_halos,                                                                &
  oh_conc, h2o2_lmt, ho2_conc, o3, hno3,                                       &
  so2_surfem, so2_hilem, so2_natem,                                            &
  dms_em_ancil, dms_conc, nh3_em,                                              &
  soot_surem, soot_hilem, bmass_surem, bmass_hilem, ocff_surem,                &
  ocff_hilem, bmass_hilem_h1, bmass_hilem_h2, land_index,                      &
  ! Logicals IN
  l_sulpc_so2, l_sulpc_dms, l_sulpc_ozone,                                     &
  l_sulpc_so2_o3_nonbuffered, l_sulpc_nh3,                                     &
  l_sulpc_online_oxidants,                                                     &
  l_sulpc_2_way_coupling,                                                      &
  l_sulphate_ccn, l_seasalt_ccn, l_soot,                                       &
  l_biomass, l_biomass_ccn, l_ocff, l_ocff_ccn,                                &
  l_nitrate, l_nitrate_ccn,                                                    &
  l_so2_surfem, l_so2_hilem, l_so2_natem, l_dms_em,                            &
  l_dms_em_inter, i_dms_flux,                                                  &
  l_nh3_em, l_ctile,                                                           &
  l_soot_surem, l_soot_hilem, l_bmass_surem, l_bmass_hilem,                    &
  l_ocff_surem, l_ocff_hilem, l_use_biogenic,                                  &
  l_use_seasalt_direct, l_use_seasalt_indirect,                                &
  l_use_seasalt_autoconv, l_use_seasalt_pm, l_dust,                            &
  ! Data fields IN/OUT
  so2, dms,                                                                    &
  so4_ait, so4_acc, so4_dis,                                                   &
  h2o2_mxr, nh3,                                                               &
  soot_new, soot_agd, soot_cld,                                                &
  bmass_new, bmass_agd, bmass_cld,                                             &
  ocff_new, ocff_agd, ocff_cld,                                                &
  biogenic, nitr_acc, nitr_diss,                                               &
  ! Data fields IN
  dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6,            &
  ! Data fields OUT
  ! Diagnostic info
  stashwork17,                                                                 &
  error_code)

!----------------------------------------------------------------------
! Purpose: Interface for Aerosol Modelling, to include Sulphur Cycle,
!          soot, biomass burning aerosol, OCFF and nitrate modelling
!          Note that dust is simply passed in for use in PM diagnostics
!          Main mineral dust routines are called from boundary layer
!          routines

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols

! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP3 v8 programming standards

! System components covered:

! System task:

! Documentation: Not yet available

!-----------------------------------------------------------------

USE planet_constants_mod, ONLY: r

USE level_heights_mod, ONLY:                          &
    r_theta_levels, r_rho_levels

USE nlsizes_namelist_mod, ONLY: model_levels

USE timestep_mod, ONLY: timestep

USE run_aerosol_mod, ONLY: so2_high_level, soot_high_level,  &
    bmass_high_level_1, bmass_high_level_2, ocff_high_level, &
    l_temporal_emi, l_bmass_hilem_variable

USE calc_pm_diags_mod, ONLY: calc_pm_diags

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE stash_array_mod, ONLY: sf
USE pws_diags_mod, ONLY: pws_dustconc_tot, pws_dustconc_surf
USE um_stashcode_mod, ONLY: stashcode_pws_sec, stashcode_aerosol_sec, &
    stashcode_pws_dustconc_surf, stashcode_pws_dustconc_5000,         &
    stashcode_aero_total_dust

USE uv_p_pnts_mod,   ONLY: uv_p_pnts
USE dms_flux_mod_4A, ONLY: dms_flux_4A
USE set_seasalt_mod, ONLY: set_seasalt_4A

USE age_aerosol_mod,      ONLY: age_aerosol
USE aero_nuclscav_mod,    ONLY: aero_nuclscav
USE bmass_variable_hilem_mod, ONLY: calc_interf_z, get_injection_fraction, &
    inject_variable_emiss
USE diagnostics_aero_mod, ONLY: diagnostics_aero
USE nitrate_mod,          ONLY: nitrate
USE sootdiffscav_mod,     ONLY: sootdiffscav
USE sulphr_mod,           ONLY: sulphr
USE number_droplet_mod,   ONLY: number_droplet
USE trsrce_mod,           ONLY: trsrce
USE day_of_week_mod,      ONLY: day_of_week
USE conversions_mod,      ONLY: pi
USE umPrintMgr
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE solang_mod, ONLY: solang
USE solpos_mod, ONLY: solpos
IMPLICIT NONE


!
! Arguments with intent IN:

! Parallel setup variables
INTEGER :: halo_i                     ! Size of halo in i direction.
INTEGER :: halo_j                     ! Size of halo in j direction.
INTEGER :: off_x                      ! Size of small halo in i
INTEGER :: off_y                      ! Size of small halo in j.
INTEGER :: global_row_length          ! number of points on a row
INTEGER :: proc_row_group             ! Group id for processors on the same row
INTEGER :: proc_col_group             ! Group id for processors on the same col
INTEGER :: global_rows                ! NUMBER OF global rows
INTEGER :: n_proc                     ! Total number of processors
INTEGER :: n_procx                    ! Number of processors in longitude
INTEGER :: n_procy                    ! Number of processors in latitude
! Array with the Ids of the four neighbours in the horizontal plane
INTEGER :: neighbour(4)
INTEGER :: g_rows (0:n_proc-1)
INTEGER :: g_row_length (0:n_proc-1)
INTEGER :: me                         ! My processor number

! Indicates if this processor is at north, south east or west of the processor grid
LOGICAL  ::  at_extremity(4)


! Model dimensions
INTEGER :: row_length
INTEGER :: rows
INTEGER :: n_rows
INTEGER :: land_points
INTEGER :: bl_levels
INTEGER :: n_cca_levels        ! No. conv cloud levels (1 if 2D, nlevs if 3D)
INTEGER :: theta_field_size
! dimensions of seasalt array
INTEGER :: salt_dim1
INTEGER :: salt_dim2
INTEGER :: salt_dim3
INTEGER :: aero_dim1           ! row_length or 1
INTEGER :: aero_dim2           ! rows or 1
INTEGER :: aero_dim3           ! model_levels or 1

! Co-ordinate arrays
REAL :: delta_lambda
REAL :: delta_phi
REAL :: lat_rot_np
REAL :: long_rot_np

! Trig arrays
REAL :: cos_theta_longitude (row_length, rows)
REAL :: sin_theta_longitude (row_length, rows)
! Grid-dependent arrays
REAL :: f3_at_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y)
REAL :: true_longitude(row_length, rows) ! in radians
REAL :: true_latitude(row_length, rows)  ! in radians
REAL :: fv_cos_theta_latitude (1-off_x:row_length+off_x,1-off_y:rows+off_y)

! Time stepping information
INTEGER :: val_year
INTEGER,INTENT(IN) :: val_month
INTEGER :: val_day_number
INTEGER,INTENT(IN) :: val_day
INTEGER :: val_hour
INTEGER :: val_minute
INTEGER :: val_second
INTEGER :: timestep_number
INTEGER :: previous_time(7)
! frequency of calling chemistry per atmos phys timestep
INTEGER :: call_chem_freq


! Diagnostics info
REAL :: stashwork17(*)  ! STASH workspace for section 17 (Aero_Ctl_4A)

REAL :: u    (1-off_x:row_length+off_x, 1-off_y:rows+off_y,    model_levels)
REAL :: v    (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,  model_levels)
REAL :: theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    model_levels)
REAL :: q    (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,model_levels)
REAL :: qcl  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,model_levels)
REAL :: qcf  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,model_levels)
REAL :: rho  (1-off_x:row_length+off_x,   1-off_y:rows+off_y,  model_levels)
REAL :: p_theta_levels(1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)
REAL :: exner_rho_levels  (1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels+1)
REAL :: exner_theta_levels(1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)
REAL :: ice_fract (row_length, rows)
REAL :: snow_depth(row_length, rows)
REAL :: land_fract(row_length, rows)
REAL :: tstar(row_length, rows)
REAL :: tstar_sea(row_length, rows)
REAL :: fland_ctile(land_points)
REAL :: cloudf_halos(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
                     model_levels)
REAL :: oh_conc(row_length,rows,model_levels)
REAL :: ho2_conc(row_length,rows,model_levels)
REAL :: h2o2_lmt(row_length,rows,model_levels)
REAL :: o3(row_length,rows,model_levels)
REAL :: so2_surfem(row_length,rows)                    ! SO2 emiss at surface
REAL :: so2_hilem(row_length,rows)                     ! SO2 emiss at chimney level
REAL :: so2_natem(row_length,rows,model_levels)        ! Volcanic SO2 emiss
REAL :: dms_em_ancil(row_length,rows)                  ! Ancillary DMS emiss (surf)
REAL :: dms_conc(row_length,rows)                      ! Seawater DMS conc (n mol l-1)
REAL :: nh3_em(row_length,rows)                        ! NH3 emiss (surface)
REAL :: soot_surem(row_length,rows)                    ! SOOT emiss at surface
REAL :: soot_hilem(row_length,rows)                    ! SOOT emiss at chimney lev
REAL :: bmass_surem(row_length,rows)                   ! Biomass surface emissions
REAL :: bmass_hilem(row_length,rows)                   ! Biomass high level emissions
REAL :: ocff_surem(row_length,rows)                    ! OCFF emiss at surface
REAL :: ocff_hilem(row_length,rows)                    ! OCFF emiss at chimney lev
REAL :: bmass_hilem_h1(row_length,rows) ! min height for biomass hilev emiss (m)
REAL :: bmass_hilem_h2(row_length,rows) ! max height for biomass hilev emiss (m)

INTEGER :: land_index(land_points)

LOGICAL :: l_sulpc_so2                  ! T if Sulphur Cycle required
LOGICAL :: l_sulpc_dms                  ! T if DMS chemistry required
LOGICAL :: l_sulpc_ozone                ! T if O3 oxidn required

! l_sulpc_so2_o3_nonbuffered T if SO2+O3 reaction is NOT to be buffered by NH3.
LOGICAL :: l_sulpc_so2_o3_nonbuffered

! l_sulpc_nh3 T if NH3 buffering required (always T if L_SULPC_OZONE is T)
LOGICAL :: l_sulpc_nh3
LOGICAL :: l_sulpc_online_oxidants      ! T if oxidants from UKCA are used.

! l_sulpc_2_way_coupling T if coupling to UKCA is 2-way,
! i.e. if depleted oxidants are passed back to UKCA
LOGICAL :: l_sulpc_2_way_coupling
LOGICAL :: l_sulphate_ccn               !T if sulphate used for CCN
LOGICAL :: l_seasalt_ccn                !T if sea-salt used for CCN
LOGICAL :: l_biomass_ccn                !T if biomass used for CCN
LOGICAL :: l_ocff_ccn                   !T if OCFF used for CCN
LOGICAL :: l_nitrate_ccn                !T if nitrate used for CCN
LOGICAL :: l_ctile                      !T if coastal tiling is on
LOGICAL :: l_soot                       !T if SOOT modelling required
LOGICAL :: l_biomass                    !T if biomass modelling reqd
LOGICAL :: l_ocff                       !T if OCFF modelling required
LOGICAL :: l_nitrate                    !T if nitrate modelling required
LOGICAL :: l_so2_surfem                 !T if surface SO2 ems present
LOGICAL :: l_so2_hilem                  !T if high lev SO2 em present
LOGICAL :: l_so2_natem                  !T if volcanic SO2 ems present
LOGICAL :: l_dms_em                     !T if DMS emiss present (surf)
LOGICAL :: l_dms_em_inter               !T if interactive DMS emiss

! Switch to determine which scheme to use for interactive DMS emissions
INTEGER :: i_dms_flux

LOGICAL :: l_nh3_em                     ! T if NH3 emiss present (surf)
LOGICAL :: l_soot_surem                 ! T if surface SOOT ems present
LOGICAL :: l_soot_hilem                 ! T if high lev SOOT ems presnt
LOGICAL :: l_bmass_surem                ! T if sfc biomass ems present
LOGICAL :: l_bmass_hilem                ! T if hi lev bmass ems present
LOGICAL :: l_ocff_surem                 ! T if surface OCFF ems present
LOGICAL :: l_ocff_hilem                 ! T if high lev OCFF ems presnt
LOGICAL :: l_use_biogenic               ! T if using biogenics for CCN
LOGICAL :: land_mask(row_length,rows)   ! T IF LAND, F IF SEA
LOGICAL :: l_dust                       ! T if mineral dust used
LOGICAL :: l_use_seasalt_direct         ! T if SS dir. rad. effect.
LOGICAL :: l_use_seasalt_indirect       ! T if SS 1st indir. effect
LOGICAL :: l_use_seasalt_autoconv       ! T if SS 2nd indir. effect
LOGICAL :: l_use_seasalt_pm

! Arguments with intent IN/OUT:
! mmr S in SO2
REAL :: so2      (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
                  model_levels)
! mmr S in DMS
REAL :: dms      (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
                  model_levels)
! mmr S in AIT
REAL,TARGET :: so4_ait  (1-off_x:row_length+off_x, 1-off_y:rows+off_y,  &
                         model_levels)
! mmr S in ACC
REAL,TARGET :: so4_acc  (1-off_x:row_length+off_x, 1-off_y:rows+off_y,  &
                         model_levels)
! mmr S in DIS
REAL,TARGET :: so4_dis  (1-off_x:row_length+off_x, 1-off_y:rows+off_y,  &
                         model_levels)
! mmr N in NH3
REAL :: nh3      (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
                  model_levels)
! mmr H2O2
REAL :: h2o2_mxr (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
                  model_levels)
! mmr fresh soot
REAL :: soot_new (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
                  model_levels)
! mmr aged soot
REAL :: soot_agd (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
                  model_levels)
! mmr soot in cloud
REAL :: soot_cld (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
                  model_levels)
! mmr fresh smoke
REAL :: bmass_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
                  model_levels)
! mmr aged smoke
REAL,TARGET :: bmass_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,  &
                         model_levels)
! mmr cloud smoke
REAL,TARGET :: bmass_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,  &
                         model_levels)
! mmr fresh OCFF
REAL :: ocff_new (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
                  model_levels)
! mmr aged OCFF
REAL,TARGET :: ocff_agd (1-off_x:row_length+off_x, 1-off_y:rows+off_y,  &
                         model_levels)
! mmr OCFF in cloud
REAL,TARGET :: ocff_cld (1-off_x:row_length+off_x, 1-off_y:rows+off_y,  &
                         model_levels)

REAL,TARGET :: biogenic(row_length, rows, model_levels)    ! mmr biogenics
REAL :: hno3(row_length,rows,model_levels)                  ! mmr HNO3

! mmr accum mode nitrate
REAL,TARGET :: nitr_acc (1-off_x:row_length+off_x, 1-off_y:rows+off_y,  &
                         model_levels)
! mmr dissolved nitrate
REAL,TARGET :: nitr_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,  &
                         model_levels)

! Arguments with intent IN:
! Dust in 6 divisions
REAL, INTENT(IN) :: dust_div1(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)
REAL, INTENT(IN) :: dust_div2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)
REAL, INTENT(IN) :: dust_div3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)
REAL, INTENT(IN) :: dust_div4(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)
REAL, INTENT(IN) :: dust_div5(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)
REAL, INTENT(IN) :: dust_div6(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
  model_levels)

! arguments with intent OUT:
INTEGER ::   error_code

! Variables with intent OUT (diagnostics):
REAL :: msa(row_length,rows,model_levels)             ! mmr S in MSA
REAL :: f_dms_to_so2(row_length,rows,model_levels)    ! frac oxid DMS to SO2
REAL :: f_dms_to_so4(row_length,rows,model_levels)    ! frac oxid DMS to SO4
REAL :: f_dms_to_msa(row_length,rows,model_levels)    ! frac oxid DMS to MSA
REAL :: nh3_dep(row_length,rows,model_levels)         ! NH3 depleted
REAL :: deltas_dry(row_length,rows,model_levels)      ! SO2 dry ox per ts
REAL :: deltas_wet(row_length,rows,model_levels)      ! SO2 wet ox by H2O2
REAL :: deltas_wet_o3(row_length,rows,model_levels)   ! SO2 wet ox by O3
REAL :: deltas_tot(row_length,rows,model_levels)      ! total SO2 ox per ts
REAL :: deltas_dms(row_length,rows,model_levels)      !DMS dry ox per ts
! SO4_DIS released by evapn of cloud droplets to SO4_ACC per ts
REAL :: deltas_evap(row_length,rows,model_levels)
! SO4_ACC transfd by nucleation to SO4_DIS per ts
REAL :: deltas_nucl(row_length,rows,model_levels)
!SO4_AIT transfd to SO4_DIS by diffusion per ts
REAL :: deltas_diffuse(row_length,rows,model_levels)
REAL :: deltas_merge(row_length,rows,model_levels)
REAL :: deltas_coag(row_length,rows,model_levels)
REAL :: psi(row_length,rows,model_levels)
! Flux thro' reaction HNO3 + NH3 --> NH4NO3
REAL :: delta_n_chem(row_length,rows,model_levels)
! nitr_diss released by evaporation of cloud droplets to nitr_diss per ts
REAL :: delta_n_evap(row_length,rows,model_levels)
! nitr_acc transferred by nucleation to nitr_diss by nucleation per ts.
REAL :: delta_n_nuc(row_length,rows,model_levels)

! Local variables

INTEGER :: i,j,k,n
INTEGER :: i_start              ! Row start point for polar row tidying
INTEGER :: iday_week            ! Day of week: 1-7
INTEGER :: min_index, max_index ! variables to check array bounds
INTEGER :: t_local              ! local hour
INTEGER :: icode                ! code to pass to ereport
INTEGER , PARAMETER :: n_month = 12 ! Total number of months in one year
INTEGER , PARAMETER :: n_hour = 24  ! Total number of hours in one day

REAL    :: chemstep             ! chemistry timestep

LOGICAL :: l_planet_obs = .FALSE. !  Calculate angles towards distant observer

! -----------------------------------------------------------------------
! Hourly factors of emissions for Europe. Calculated by TNO for the
! MACC project.
!
! Source:
! Hugo Denier van der Gon, Carlijn Hendriks, Jeroen Kuenen, Arjo Segers,
! Antoon Visschedijk: Description of current temporal emission patterns
! and sensitivity of predicted AQ for temporal emission patterns.
! EU FP7 MACC deliverable report D_D-EMIS_1.3, TNO report, Dec 2011.
! https://gmes-atmosphere.eu/documents/deliverables/d-emis/
!         MACC_TNO_del_1_3_v2.pdf
!
! The following hourly factors have been taken from Table on pg 22 of
! the report for SNAP sector 7 (traffic)
REAL, PARAMETER :: hourly_scaling(0:23)=                                  &
  (/ 0.19,0.09,0.06,0.05,0.09,0.22,0.86,1.84,1.86,1.41,1.24,1.2,          &
     1.32,1.44,1.45,1.59,2.03,2.08,1.51,1.06,0.74,0.62,0.61,0.44 /)

! The following hourly factors have been calculated from monthly
! average daily solar cycles for 54.3 N, 0.0W, 2014
REAL, PARAMETER :: hourly_scaling_nh3 (0:n_hour-1, n_month) = RESHAPE ((/ &
  0.43,0.43,0.43,0.43,0.43,0.43,0.43,0.43,0.43,1.01,2.39,2.93,            & !Jan
  3.10,3.01,2.60,1.54,0.45,0.43,0.43,0.43,0.43,0.43,0.43,0.43,            & 
  0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.68,1.87,2.41,2.65,            & !Feb
  2.73,2.70,2.55,2.18,1.30,0.34,0.33,0.33,0.33,0.33,0.33,0.33,            & 
  0.26,0.26,0.26,0.26,0.26,0.26,0.26,0.69,1.59,2.01,2.21,2.30,            & !Mar
  2.33,2.32,2.24,2.09,1.77,1.02,0.27,0.26,0.26,0.26,0.26,0.26,            & 
  0.23,0.23,0.23,0.23,0.23,0.23,0.69,1.39,1.72,1.89,1.98,2.02,            & !Apr
  2.04,2.03,1.98,1.89,1.73,1.41,0.73,0.23,0.23,0.23,0.23,0.23,            & 
  0.21,0.21,0.21,0.21,0.21,0.51,1.15,1.47,1.65,1.76,1.82,1.85,            & !May
  1.86,1.85,1.81,1.75,1.64,1.45,1.10,0.44,0.21,0.21,0.21,0.21,            & 
  0.21,0.21,0.21,0.21,0.21,0.72,1.19,1.45,1.60,1.69,1.74,1.77,            & !Jun
  1.78,1.77,1.74,1.69,1.60,1.45,1.20,0.74,0.21,0.21,0.21,0.21,            & 
  0.22,0.22,0.22,0.22,0.22,0.52,1.11,1.43,1.61,1.71,1.78,1.81,            & !Jul
  1.83,1.82,1.79,1.73,1.63,1.47,1.18,0.66,0.22,0.22,0.22,0.22,            & 
  0.23,0.23,0.23,0.23,0.23,0.24,0.83,1.39,1.67,1.83,1.92,1.97,            & !Aug
  1.98,1.97,1.93,1.84,1.70,1.43,0.91,0.27,0.23,0.23,0.23,0.23,            & 
  0.27,0.27,0.27,0.27,0.27,0.27,0.32,1.17,1.74,2.01,2.15,2.22,            & !Sep
  2.24,2.21,2.13,1.97,1.65,0.98,0.30,0.27,0.27,0.27,0.27,0.27,            & 
  0.32,0.32,0.32,0.32,0.32,0.32,0.32,0.47,1.52,2.16,2.44,2.57,            & !Oct
  2.59,2.52,2.33,1.91,0.96,0.33,0.32,0.32,0.32,0.32,0.32,0.32,            &
  0.42,0.42,0.42,0.42,0.42,0.42,0.42,0.42,0.59,1.87,2.63,2.93,            & !Nov
  2.98,2.82,2.34,1.17,0.42,0.42,0.42,0.42,0.42,0.42,0.42,0.42,            & 
  0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.86,2.44,3.07,            & !Dec
  3.23,3.02,2.29,0.57,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50 /),         & 
  (/n_hour, n_month/))

! Daily factors of emissions derived for Great Britain, based on
! calculations by TNO for the MACC project (Denier van der Gon et
! al., 2011) for SNAP sector 7 (traffic)
! Days are: Sun, Mon, ..., Sat:
REAL, PARAMETER :: weekly_scaling(7)=                                          &
  (/0.886409, 1.040343, 1.017710, 1.033410, 1.037143, 1.105577, 0.879408/) 

! the actual scaling applied is a function of longitude, use a 2D array 
! to hold the factors to multiply the 2D emissions by. Can't use 
! vector here as true longitude in a LAM is not fixed for constant
! values of i.
REAL :: emiscale(row_length,rows)
REAL :: emiscale_nh3(row_length,rows)

! Diagnostic increments generated by each call of SULPHR
REAL :: msa_inc    (row_length,rows,model_levels)
REAL :: nh3_dep_inc(row_length,rows,model_levels)

REAL :: t     (row_length,rows,model_levels)   ! Temp (K) on theta levels
REAL :: t_surf(row_length,rows)                ! Surface temperature (K)

! For calculation of cos zenith angle
REAL :: cosza2d(row_length,rows)      !cos zenith angle
REAL :: sindec                        ! sin(solar declination)
REAL :: sol_azimuth(row_length,rows)
REAL :: cosz_beg(row_length,rows)
REAL :: cosz_end(row_length,rows)
REAL :: scs                           ! solar constant scaling factor
REAL :: seconds_since_midnight
REAL :: eq_time                       ! The Equation of time
REAL :: day_fraction(row_length,rows)
REAL :: sindec_obs, eqt_obs           ! dummy variables here

! For sea-salt aerosol and DMS emissions calculations
REAL :: u_1(aero_dim1,aero_dim2)                         ! level 1 wind u
REAL :: v_1(aero_dim1,aero_dim2)                         ! level 1 wind v
REAL :: windspeed_1(aero_dim1,aero_dim2)                 ! level 1 windspeed
REAL :: windspeed_10m(aero_dim1,aero_dim2)               ! 10 m windspeed
! Sea-salt film aerosol.
REAL,TARGET :: sea_salt_film(salt_dim1,salt_dim2,salt_dim3)
! Sea-salt jet aerosol.
REAL,TARGET :: sea_salt_jet(salt_dim1,salt_dim2,salt_dim3)
REAL :: height_theta(aero_dim1,aero_dim2,aero_dim3)      ! Theta level centre heights
REAL :: height_rho_1(aero_dim1,aero_dim2)                ! Rho level 1 centre heights
REAL :: dms_em_inter(aero_dim1,aero_dim2)                ! Interactive DMS emissns
REAL :: dms_em_combined(row_length,rows)                 ! Combined DMS emissns

! Increments from soot ageing and diffusional scavenging
! Increment from soot ageing
REAL :: delta_agesoot(row_length,rows,model_levels)
! Increment from soot diff. scav.
REAL :: delta_sootdiffscav(row_length,rows,model_levels)

! Increments from biomass smoke ageing and nucleation scavenging
! Increment from smoke ageing
REAL :: delta_agebmass(row_length,rows,model_levels)
! Increment from smoke nucl. scav.
REAL :: delta_bmassnuclscav(row_length,rows,model_levels)
! For emissions of biomass smoke
! Fraction to inject in each layer
REAL :: bmass_frac_em(row_length,rows,model_levels)
! Model level of biomass emissions
INTEGER ::  biomass_level

! Quantity of biomass smoke emitted per model level
REAL    ::  biomass_per_level(row_length,rows)

! Altitude of interfaces (metres) above ground level
REAL, SAVE, ALLOCATABLE :: interf_z (:,:,:)

! Increments from fossil-fuel organic carbon ageing and nucl. scavenging
REAL :: delta_ageocff(row_length, rows, model_levels)      ! Increment from OCFF ageing
REAL :: delta_ocffnuclscav(row_length, rows, model_levels) ! Increment from OCFF nucl. scav.


! For cloud fraction without halos
REAL :: cloudf(row_length,rows,model_levels)

LOGICAL,PARAMETER :: l_sulpc_newdms= .TRUE.       ! TRUE if new DMS scheme required
LOGICAL, SAVE     :: l_initialise_injection_levels = .TRUE. ! for BB elev emiss

! For droplet number calculation
REAL,TARGET :: dummy(row_length, rows, model_levels)
REAL :: n_droplet(row_length, rows, model_levels)
REAL :: rho_air(row_length, rows, model_levels)
REAL, POINTER :: ss_film(:,:,:)
REAL, POINTER :: ss_jet(:,:,:)
REAL, POINTER :: sulp_acc(:,:,:)
REAL, POINTER :: sulp_dis(:,:,:)
REAL, POINTER :: agd_bmass(:,:,:)
REAL, POINTER :: cld_bmass(:,:,:)
REAL, POINTER :: agd_ocff(:,:,:)
REAL, POINTER :: cld_ocff(:,:,:)
REAL, POINTER :: acc_nitrate(:,:,:)
REAL, POINTER :: dis_nitrate(:,:,:)
REAL, POINTER :: biogenic_ccn(:,:,:)

! For PM10 & PM2.5 calculation
REAL :: pm10       (row_length, rows, model_levels)
REAL :: pm2p5      (row_length, rows, model_levels)
REAL :: pm10_so4   (row_length, rows, model_levels)
REAL :: pm2p5_so4  (row_length, rows, model_levels)
REAL :: pm10_bc    (row_length, rows, model_levels)
REAL :: pm2p5_bc   (row_length, rows, model_levels)
REAL :: pm10_bb    (row_length, rows, model_levels)
REAL :: pm2p5_bb   (row_length, rows, model_levels)
REAL :: pm10_ocff  (row_length, rows, model_levels)
REAL :: pm2p5_ocff (row_length, rows, model_levels)
REAL :: pm10_soa   (row_length, rows, model_levels)
REAL :: pm2p5_soa  (row_length, rows, model_levels)
REAL :: pm10_ss    (row_length, rows, model_levels)
REAL :: pm2p5_ss   (row_length, rows, model_levels)
REAL :: conc_dust  (row_length, rows, model_levels)
REAL :: conc_dust_surf  (row_length, rows)
REAL :: pm10_dust  (row_length, rows, model_levels)
REAL :: pm2p5_dust (row_length, rows, model_levels)
REAL :: pm10_nitr  (row_length, rows, model_levels)
REAL :: pm2p5_nitr (row_length, rows, model_levels)

REAL, PARAMETER :: z0_sea=2.5e-04              ! Roughness length over sea (m)
REAL, PARAMETER :: p1 = LOG(10.0/z0_sea)

! Is the input "sulphate" aerosol in the form of
! ammonium sulphate (T) or just sulphur (F)? Used in the
! same way for the nitrate aerosol, in form of ammonium
! nitrate (T) or just nitrogen (F).
LOGICAL :: l_nh42so4 = .FALSE.

CHARACTER(LEN=errormessagelength) :: cmessage        ! Error message
CHARACTER (LEN=* ), PARAMETER ::  RoutineName='AERO_CTL_4A'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! Set up global land fraction field

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(k, j, i)

IF (l_ctile) THEN

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      land_fract(i,j) = 0.0
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO k = 1,land_points
    j = (land_index(k)-1)/row_length + 1
    i = land_index(k) - (j-1)*row_length
    land_fract(i,j) = fland_ctile(k)
  END DO
!$OMP END DO NOWAIT

ELSE

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      IF (land_mask(i,j)) THEN
        land_fract(i,j) = 1.0
      ELSE
        land_fract(i,j) = 0.0
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT

END IF


! Copy cloud fraction into an array without halos and set to zero if
! smaller than machine precision

!$OMP DO SCHEDULE(STATIC)
DO k = 1,model_levels
  DO j = 1,rows
    DO i = 1,row_length
      cloudf(i,j,k)=cloudf_halos(i,j,k)
      IF (cloudf(i,j,k) <  EPSILON(cloudf(i,j,k)))                             &
        cloudf(i,j,k)=0.0
      IF (cloudf(i,j,k) >  1.0) THEN
        cloudf(i,j,k) = 1.0
      END IF
    END DO
  END DO
END DO
!$OMP END DO NOWAIT


!$OMP END PARALLEL

! If sea-salt aerosol or interactive DMS emissions are to be
! used then calculate surface windspeed and height information

IF (l_seasalt_ccn .OR. l_use_seasalt_pm .OR. l_dms_em_inter) THEN

  ! calculate level 1 winds on p points using function
  CALL uv_p_pnts(u(1:row_length,1:rows,1),v(1:row_length,1:n_rows,1), &
                 cos_theta_longitude,sin_theta_longitude,             &
                 global_row_length,proc_row_group,u_1,v_1)

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(k, j, i)

  ! calculate level 1 windspeed (NB on rho levels)
  ! Height of theta levels above the surface (used in set_seasalt_4A)
!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      windspeed_1(i,j)=SQRT(u_1(i,j)**2+v_1(i,j)**2)
    END DO
  END DO
!$OMP END DO NOWAIT

  ! Height of theta levels above the surface (used in set_seasalt_4A)
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, aero_dim3
    DO j = 1, aero_dim2
      DO i = 1, aero_dim1
        height_theta(i,j,k) = r_theta_levels(i,j,k)                           &
          -r_theta_levels(i,j,0)
      END DO
    END DO
  END DO
!$OMP END DO

  ! Height of rho level 1 above the surface
!$OMP DO SCHEDULE(STATIC)
  DO j = 1, aero_dim2
    DO i = 1, aero_dim1
      height_rho_1(i,j) = r_rho_levels(i,j,1) - r_theta_levels(i,j,0)
    END DO
  END DO
!$OMP END DO NOWAIT

  ! Adjust level 1 windspeed to be 10m windspeed
!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      windspeed_10m(i,j)=windspeed_1(i,j)*p1/(LOG(height_theta(i,j,1)/z0_sea))
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

END IF ! L_seasalt_CCN .OR. L_DMS_em_inter

IF (l_seasalt_ccn .OR. l_use_seasalt_pm) THEN
  CALL set_seasalt_4A(windspeed_10m, height_theta, land_fract, ice_fract,      &
    row_length, rows, model_levels,                                            &
    salt_dim1, salt_dim2, salt_dim3,                                           &
    bl_levels, sea_salt_film, sea_salt_jet)
ELSE
  sea_salt_film(1,1,1) = 0.0
  sea_salt_jet(1,1,1) = 0.0
END IF

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(k, j, i)

! Calculate T from theta if soot, S cycle, biomass or OCFF switched on
! and droplet numbers for diffusional scavenging.

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      t(i,j,k)=theta(i,j,k)*exner_theta_levels(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

IF (l_sulpc_so2 .OR. l_soot .OR. l_biomass .OR. l_ocff) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        rho_air(i,j,k)=p_theta_levels(i,j,k)/                      &
          (r * t(i,j,k))

      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

END IF ! L_SULPC_SO2 .OR. L_SOOT .OR. L_BIOMASS .OR. L_OCFF

!$OMP END PARALLEL

IF (l_sulpc_so2 .OR. l_soot .OR. l_biomass .OR. l_ocff) THEN

  IF (l_seasalt_ccn) THEN
    ss_film=>sea_salt_film(1:row_length, 1:rows, :)
    ss_jet=>sea_salt_jet(1:row_length, 1:rows, :)
  ELSE
    ss_film=>dummy
    ss_jet=>dummy
  END IF

  IF (l_sulphate_ccn) THEN
    sulp_acc=>so4_acc(1:row_length, 1:rows, :)
    sulp_dis=>so4_dis(1:row_length, 1:rows, :)
  ELSE
    sulp_acc=>dummy
    sulp_dis=>dummy
  END IF

  IF (l_biomass_ccn) THEN
    agd_bmass=>bmass_agd(1:row_length, 1:rows, :)
    cld_bmass=>bmass_cld(1:row_length, 1:rows, :)
  ELSE
    agd_bmass=>dummy
    cld_bmass=>dummy
  END IF

  IF (l_ocff_ccn) THEN
    agd_ocff=>ocff_agd(1:row_length, 1:rows, :)
    cld_ocff=>ocff_cld(1:row_length, 1:rows, :)
  ELSE
    agd_ocff=>dummy
    cld_ocff=>dummy
  END IF

  IF (l_nitrate_ccn) THEN
    acc_nitrate=>nitr_acc(1:row_length, 1:rows, :)
    dis_nitrate=>nitr_diss(1:row_length, 1:rows, :)
  ELSE
    acc_nitrate=>dummy
    dis_nitrate=>dummy
  END IF

  IF (l_use_biogenic) THEN
    biogenic_ccn=>biogenic(1:row_length, 1:rows, :)
  ELSE
    biogenic_ccn=>dummy
  END IF

  CALL number_droplet(                                                   &
    1, row_length,                                                       &
    1, rows,                                                             &
    1, model_levels,                                                     &
    1, model_levels,                                                     &
    l_sulphate_ccn,l_nh42so4,                                            &
    sulp_acc,sulp_dis,                                                   &
    l_seasalt_ccn,ss_film,ss_jet,                                        &
    l_use_biogenic,biogenic_ccn,                                         &
    l_biomass_ccn,agd_bmass,cld_bmass,                                   &
    l_ocff_ccn,agd_ocff,cld_ocff,                                        &
    l_nitrate_ccn,acc_nitrate,dis_nitrate,                               &
    rho_air,snow_depth,land_fract,                                       &
    n_droplet )
END IF

IF (l_initialise_injection_levels) THEN
  ! If injection heights of biomass burning emissions are grid cell dependent
  ! then calculate altitude above ground for the interfaces (in metres).
  ! This will be used to find the model levels which correspond to the
  ! minimum and maximum injection heights.
  IF (l_bmass_hilem_variable) THEN
    ! Fill in array with altitudes of the interfaces at theta levels
    ALLOCATE (interf_z (row_length, rows, 0:model_levels))
    CALL calc_interf_z (row_length, rows, model_levels, interf_z)
  END IF ! l_bmass_hilem_variable
END IF   ! l_initialise_injection_levels

! calculate scaling for diurnal and weekly cycles, only if this option is on
IF (l_temporal_emi) THEN

  ! Find which day of the week it is
  iday_week=day_of_week(val_day,val_month,val_year)
  
  ! Check that won't get array bounds exception from the code below.
  ! Note t_local is used to index an array indexed from 0:23.
  ! If t_local is less than 0 we add 24, if t_local is more than
  ! 23 we subtract 24. Therefore valid values are from -24 to +47
  
  min_index =  NINT (MINVAL(true_longitude)*24.0/(2.0*pi)) + val_hour
  max_index =  NINT (MAXVAL(true_longitude)*24.0/(2.0*pi)) + val_hour
  
  IF (min_index < -24 .OR. max_index > 47) THEN
      cmessage='Invalid index for emissions t-profile'
      icode = 1
      CALL ereport (routinename, icode, cmessage)
  END IF 

  DO j = 1, rows
    DO i = 1, row_length
      ! Take nearest integer so that at 0Z all cells of -7.5 to +7.5
      ! degrees use the hourly scaling indexed at 0
      t_local = NINT (true_longitude(i,j)*24.0/(2.0*pi)) + val_hour
      !
      ! Ensure that t_local is bounded to be between 0 and 23
      IF (t_local >= 24) t_local = t_local - 24
      IF (t_local <   0) t_local = t_local + 24

      ! calculate scaling for diurnal cycle as function of local time
      ! scaling is product of weekly and hourly factors
      emiscale (i, j) = hourly_scaling (t_local)*weekly_scaling(iday_week)
      emiscale_nh3 (i, j) = hourly_scaling_nh3 (t_local, val_month)
    END DO
  END DO

ELSE

  ! If temporal scaling is off, factor is one
  emiscale(:, :) = 1.0
  emiscale_nh3(:, :) = 1.0
END IF

! Soot cycle code.
! Calculate emissions of soot from surface and chimney level, ageing
! of soot and diffusional scavenging of aged soot to soot-in-cloud.

IF (l_soot) THEN  ! If soot modelling is included

  IF (l_soot_surem) THEN  ! If surface soot emissions are included
    CALL trsrce(                                                               &
      rows, row_length, off_x, off_y, halo_i, halo_j,                          &
      theta, q , qcl , qcf , exner_rho_levels, rho,                            &
      soot_new(:,:,1), soot_surem(:, :)*emiscale(:, :), 1,                     &
      timestep, 1, 1, 0.0 )
  END IF  ! L_SOOT_SUREM

  IF (l_soot_hilem) THEN  ! If chimney level soot emissions
    ! are included
    CALL trsrce(                                                               &
      rows, row_length, off_x, off_y, halo_i, halo_j,                          &
      theta, q , qcl , qcf , exner_rho_levels, rho,                            &
      soot_new(:,:,soot_high_level), soot_hilem,                               &
      soot_high_level, timestep, 1, 1, 0.0 )
  END IF  ! L_SOOT_HILEM

  ! Calculate quantity of fresh soot converted to aged soot

  CALL age_aerosol(                                                            &
    row_length, rows, off_x, off_y,                                            &
    model_levels,                                                              &
    soot_new,                                                                  &
    delta_agesoot, 'soot'                                                      &
    )

  ! Calculate quantity of aged soot scavenged to cloud soot

  CALL sootdiffscav(                                                           &
    rows, row_length, off_x, off_y, halo_i, halo_j,                            &
    cloudf, qcl, qcf, p_theta_levels, t,                                       &
    n_droplet,                                                                 &
    soot_agd, soot_cld,                                                        &
    delta_sootdiffscav                                                         &
    )

  ! Update soot arrays with increments from ageing and
  ! diffusional scavenging

  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        soot_new(i,j,k)=soot_new(i,j,k) - delta_agesoot(i,j,k)
        soot_agd(i,j,k)=soot_agd(i,j,k) + delta_agesoot(i,j,k)                 &
          - delta_sootdiffscav(i,j,k)
        soot_cld(i,j,k)=soot_cld(i,j,k)                                        &
          + delta_sootdiffscav(i,j,k)
      END DO
    END DO
  END DO

END IF  ! L_SOOT

! End of soot cycle code

! Biomass aerosol code.
! Calculate emissions of smoke from surface and high level, ageing
! of smoke and nucleation scavenging of aged smoke to smoke-in-cloud.

IF (l_biomass) THEN  ! If biomass smoke modelling is included

  IF (l_bmass_surem) THEN  ! If surface smoke emissions included
    CALL trsrce(                                                               &
      rows, row_length, off_x, off_y, halo_i, halo_j,                          &
      theta, q , qcl , qcf , exner_rho_levels, rho,                            &
      bmass_new(:,:,1), bmass_surem, 1,                                        &
      timestep, 1, 1, 0.0  )
  END IF  ! L_BMASS_SUREM

  IF (l_bmass_hilem) THEN  ! If high level smoke emissions are included

    ! Different injection heights for each grid cell
    IF (l_bmass_hilem_variable) THEN
      ! Get 2-D varying injection layers
      CALL get_injection_fraction (row_length, rows, model_levels,             &
        bmass_hilem_h1 (1:row_length, 1:rows),                                 &
        bmass_hilem_h2 (1:row_length, 1:rows),                                 &
        interf_z       (1:row_length, 1:rows, 0:model_levels),                 &
        bmass_frac_em  (1:row_length, 1:rows, 1:model_levels))

      ! Update fresh biomass burning aerosol mmr (bmass_new) with high-level
      ! emissions (bmass_hilem) within the 2-D varying model layers 
      CALL inject_variable_emiss (                                             &
        row_length, rows, model_levels, off_x, off_y, halo_i, halo_j,          &
        timestep, bmass_hilem, bmass_frac_em,                                  &
        theta, rho, q, qcl, qcf, exner_rho_levels, bmass_new)

    ! Fixed values for first and last layers of injection
    ELSE
      ! Check that the range of emission levels is correctly specified
      ! and, if so, emit equal quantities of smoke on all model levels
      ! between bmass_high_level_1 and bmass_high_level_2.

      IF (bmass_high_level_1  >   bmass_high_level_2 .OR.                      &
        bmass_high_level_1  >   model_levels .OR.                              &
        bmass_high_level_2  >   model_levels) THEN
        WRITE(umMessage,'(A)') 'Aero_Ctl_4A: Invalid range of biomass emission '
        CALL umPrint(umMessage,src=routinename)
        WRITE(umMessage,'(A)') 'levels specified.'
        CALL umPrint(umMessage,src=routinename)
        WRITE(umMessage,'(A,I5)') 'Lowest level:  ', bmass_high_level_1
        CALL umPrint(umMessage,src=routinename)
        WRITE(umMessage,'(A,I5)') 'Highest level: ', bmass_high_level_2
        CALL umPrint(umMessage,src=routinename)

        icode    = 2
        cmessage = 'Inconsistency found for minimum and maximum levels '    // &
                   'of injection for biomass burning aerosol emissions'
        CALL ereport (routinename, icode, cmessage)
      ELSE

        DO j=1,rows
          DO i=1,row_length
            biomass_per_level(i,j) = bmass_hilem(i,j)/                         &
              ((bmass_high_level_2 - bmass_high_level_1) + 1)
          END DO
        END DO

        DO biomass_level = bmass_high_level_1, bmass_high_level_2

          CALL trsrce(                                                         &
            rows, row_length, off_x, off_y, halo_i, halo_j,                    &
            theta, q , qcl , qcf , exner_rho_levels, rho,                      &
            bmass_new(:,:,biomass_level), biomass_per_level,                   &
            biomass_level, timestep, 1, 1, 0.0 )

        END DO  ! biomass_level

      END IF  ! Test on emissions levels

    END IF  ! Injection within variable or fixed model levels

  END IF  ! L_BMASS_HILEM

  ! Calculate quantity of fresh smoke converted to aged smoke

  CALL age_aerosol(                                                            &
    row_length, rows, off_x, off_y,                                            &
    model_levels,                                                              &
    bmass_new,                                                                 &
    delta_agebmass,  'bmas'                                                    &
    )

  ! Calculate quantity of aged smoke scavenged to cloud smoke

  CALL aero_nuclscav(                                                          &
    rows, row_length, off_x, off_y, halo_i, halo_j,                            &
    cloudf, qcl, qcf,                                                          &
    bmass_agd, bmass_cld,                                                      &
    delta_bmassnuclscav                                                        &
    )

  ! Update smoke arrays with increments from ageing and
  ! nucleation scavenging

  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        bmass_new(i,j,k)=bmass_new(i,j,k) - delta_agebmass(i,j,k)
        !    Simulate the condensation of VOCs onto aged biomass by
        !    increasing the mass transferred upon ageing by (8.75/5.4),
        !    the ratio of the fraction of BC in fresh (8.75%) and aged
        !    (5.4%) biomass aerosol:
        bmass_agd(i,j,k)=bmass_agd(i,j,k)                                      &
          + ((8.75/5.4)*delta_agebmass(i,j,k))                                 &
          - delta_bmassnuclscav(i,j,k)
        bmass_cld(i,j,k)=bmass_cld(i,j,k)                                      &
          + delta_bmassnuclscav(i,j,k)
      END DO
    END DO
  END DO

END IF  ! L_BIOMASS

! End of biomass aerosol code


IF ( l_sulpc_so2 ) THEN

  ! Calculate cos zenith angle for SULPH2 by calling SOLPOS and SOLANG
  ! Calculate number of seconds since midnight to the beginning of
  ! the timetsep.

  seconds_since_midnight = REAL( previous_time(4) * 3600                       &
    + previous_time(5) * 60  + previous_time(6))

  CALL solpos (previous_time(7), previous_time(1),                             &
    seconds_since_midnight, timestep,l_planet_obs ,                            &
    eq_time, sindec, scs, sindec_obs, eqt_obs)

  CALL solang(                                                                 &
                                ! input constants
    sindec, seconds_since_midnight,                                            &
    timestep, eq_time,                                                         &
                                ! row and column dependent constants
    true_latitude,                                                             &
    true_longitude,                                                            &
                                ! size variables
    row_length*rows,                                                           &
                                ! output fields
    day_fraction, cosza2d, sol_azimuth, cosz_beg, cosz_end )


  ! Calculate interactive DMS emissions if required

  IF (l_dms_em_inter) THEN

    DO j = 1, rows
      DO i = 1, row_length
        IF (l_ctile) THEN
          t_surf(i, j) = tstar_sea(i, j)
        ELSE
          t_surf(i, j) = tstar(i, j)
        END IF
      END DO
    END DO

    CALL dms_flux_4A (                                                         &
      ! size variables
      row_length, rows,                                                        &
      ! input fields
      windspeed_10m,                                                           &
      t_surf,                                                                  &
      land_fract,                                                              &
      dms_conc,                                                                &
      ! switch
      i_dms_flux,               &
      ! output field
      dms_em_inter )

  END IF


  ! zero output diagnostics for full model timestep
  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        msa(i,j,k) = 0.0
        nh3_dep(i,j,k) = 0.0
      END DO
    END DO
  END DO

  ! Calculate length of chemistry timestep and use it to control input
  ! of emissions and S Chemistry

  chemstep=timestep/call_chem_freq

  DO n=1,call_chem_freq

    ! Call TRSRCE to insert emissions

    IF (l_so2_surfem) THEN         ! Insert surface SO2 emiss
      CALL trsrce(                                                             &
        rows, row_length, off_x, off_y, halo_i, halo_j,                        &
        theta, q , qcl , qcf , exner_rho_levels, rho,                          &
        so2(:,:,1), so2_surfem(:, :)*emiscale(:, :), 1,                        &
        chemstep, 1, 1, 0.0 )
    END IF

    IF (l_so2_hilem) THEN          ! Insert chimney SO2 emiss
      CALL trsrce(                                                             &
        rows, row_length, off_x, off_y, halo_i, halo_j,                        &
        theta, q , qcl , qcf , exner_rho_levels, rho,                          &
        so2(:,:,so2_high_level), so2_hilem, so2_high_level,                    &
        chemstep, 1, 1, 0.0 )
    END IF

    IF (l_so2_natem) THEN          ! Insert volcanic SO2 emiss
      DO k= 1,model_levels
        CALL trsrce(                                                           &
          rows, row_length, off_x, off_y, halo_i, halo_j,                      &
          theta, q , qcl , qcf , exner_rho_levels, rho,                        &
          so2(:,:,k), so2_natem(:,:,k), k,                                     &
          chemstep, 1, 1, 0.0 )
      END DO
    END IF

    IF (l_dms_em_inter) THEN
      DO j = 1, rows
        DO i = 1, row_length
          dms_em_combined(i, j)                                                &
            = (land_fract(i, j) * dms_em_ancil(i, j))                          &
            + ( ((1.0 - land_fract(i, j)) * dms_em_inter(i, j))                &
            * (1.0 - ice_fract(i, j)) )
        END DO
      END DO
    ELSE
      IF (l_dms_em) THEN     ! Just copy over the standard ancil
        DO j = 1, rows
          DO i = 1, row_length
            dms_em_combined(i, j) = dms_em_ancil(i, j)
          END DO
        END DO
      END IF
    END IF

    IF (l_dms_em) THEN             ! Insert DMS emiss (surface)
      CALL trsrce(                                                             &
        rows, row_length, off_x, off_y, halo_i, halo_j,                        &
        theta, q , qcl , qcf , exner_rho_levels, rho,                          &
        dms(:,:,1), dms_em_combined, 1,                                        &
        chemstep, 1, 1, 0.0 )
    END IF

    IF (l_nh3_em) THEN             ! Insert NH3 emiss (surface)
      CALL trsrce(                                                             &
        rows, row_length, off_x, off_y, halo_i, halo_j,                        &
        theta, q , qcl , qcf , exner_rho_levels, rho,                          &
        nh3(:,:,1), nh3_em(:, :)*emiscale_nh3(:, :), 1,                        &
        chemstep, 1, 1, 0.0 )
    END IF


    CALL sulphr(                                                               &
                                ! Arguments IN
      halo_i, halo_j, off_x, off_y,                                            &
      row_length, rows,                                                        &
      theta_field_size,                                                        &
      chemstep,                                                                &
      cloudf, cosza2d,                                                         &
      p_theta_levels, t, q, qcl, qcf,                                          &
      oh_conc, h2o2_lmt, ho2_conc, o3,                                         &
      n_droplet,                                                               &
      l_sulpc_dms, l_sulpc_newdms,                                             &
      l_sulpc_ozone,                                                           &
      l_sulpc_so2_o3_nonbuffered,                                              &
      l_sulpc_nh3,                                                             &
      l_sulpc_online_oxidants,                                                 &
      l_sulpc_2_way_coupling,                                                  &
                                ! Arguments IN/OUT
      so2, dms, so4_ait, so4_acc, so4_dis,                                     &
      nh3, h2o2_mxr,                                                           &
                                ! Arguments OUT (diagnostics)
      msa_inc,                                                                 &
      nh3_dep_inc,                                                             &
      f_dms_to_so2, f_dms_to_so4, f_dms_to_msa,                                &
      deltas_dry, deltas_wet, deltas_wet_o3,                                   &
      deltas_tot, deltas_dms,                                                  &
      deltas_evap, deltas_nucl, deltas_diffuse,                                &
      deltas_merge,                                                            &
      deltas_coag, psi )

    ! Add diagnostic increments from SULPHR to total for model timestep

    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          msa(i,j,k)=msa(i,j,k)+msa_inc(i,j,k)
          nh3_dep(i,j,k)=nh3_dep(i,j,k)+nh3_dep_inc(i,j,k)
        END DO
      END DO
    END DO

    IF (l_nitrate) THEN
      ! Call nitrate chemistry routine
      CALL nitrate(                                                            &
                                ! Arguments IN
        halo_i, halo_j, off_x, off_y,                                          &
        row_length, rows,                                                      &
        chemstep, cloudf, p_theta_levels, t, q, qcl, qcf,                      &
                                ! Arguments IN/OUT
        hno3, nh3, nitr_acc, nitr_diss,                                        &
                                ! Arguments OUT
        delta_n_chem, delta_n_evap, delta_n_nuc )
      !
    END IF           ! End L_NITRATE condition

  END DO             ! End CALL_CHEM_FREQ loop

END IF               ! End L_SULPC_SO2 test


! OCFF cycle code.
! Calculate emissions of ocff from surface and chimney level, ageing
! of ocff and nucleation scavenging of aged ocff to ocff-in-cloud.

IF (l_ocff) THEN  ! If ocff modelling is included

  IF (l_ocff_surem) THEN  ! If surface ocff emissions are included
    CALL trsrce(                                                               &
      rows, row_length, off_x, off_y, halo_i, halo_j,                          &
      theta, q , qcl , qcf , exner_rho_levels, rho,                            &
      ocff_new(:,:,1), ocff_surem(:, :)*emiscale(:, :), 1,                     &
      timestep, 1, 1, 0.0 )
  END IF  ! L_OCFF_SUREM

  IF (l_ocff_hilem) THEN  ! If chimney level ocff emissions
    ! are included
    CALL trsrce(                                                               &
      rows, row_length, off_x, off_y, halo_i, halo_j,                          &
      theta, q , qcl , qcf , exner_rho_levels, rho,                            &
      ocff_new(:,:,ocff_high_level), ocff_hilem,                               &
      ocff_high_level, timestep, 1, 1, 0.0 )
  END IF  ! L_OCFF_HILEM

  ! Calculate quantity of fresh ocff converted to aged ocff

  CALL age_aerosol(                                                            &
    row_length, rows, off_x, off_y,                                            &
    model_levels,                                                              &
    ocff_new,                                                                  &
    delta_ageocff, 'ocff'                                                      &
    )

  ! Calculate quantity of aged ocff scavenged to cloud ocff

  CALL aero_nuclscav(                                                          &
    rows, row_length, off_x, off_y, halo_i, halo_j,                            &
    cloudf, qcl, qcf,                                                          &
    ocff_agd, ocff_cld,                                                        &
    delta_ocffnuclscav                                                         &
    )

  ! Update ocff arrays with increments from ageing and
  ! nucleation scavenging

  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        ocff_new(i,j,k)=ocff_new(i,j,k) - delta_ageocff(i,j,k)
        ocff_agd(i,j,k)=ocff_agd(i,j,k) + delta_ageocff(i,j,k)                 &
          - delta_ocffnuclscav(i,j,k)
        ocff_cld(i,j,k)=ocff_cld(i,j,k)                                        &
          + delta_ocffnuclscav(i,j,k)
      END DO
    END DO
  END DO

END IF  ! L_OCFF

! End of ocff cycle code



! PM concentration code: If any of the diagnostics on PM10 and PM2.5
! or total mass concentrations (item numbers 220-235, Sect. 17) is requested
! then call the subroutine that considers all aerosol species and modes

IF (sf(220,17) .OR. sf(221,17) .OR. sf(222,17) .OR.                            &
    sf(223,17) .OR. sf(224,17) .OR. sf(225,17) .OR.                            &
    sf(226,17) .OR. sf(227,17) .OR. sf(228,17) .OR.                            &
    sf(229,17) .OR. sf(230,17) .OR. sf(231,17) .OR.                            &
    sf(232,17) .OR. sf(233,17) .OR. sf(234,17) .OR.                            &
    sf(235,17) .OR. sf(236,17) .OR. sf(237,17) .OR.                            &
    sf(257,17) ) THEN

  CALL calc_pm_diags (                                                         &
    off_x, off_y,                                                              &
    row_length, rows,                                                          &
    model_levels,                                                              &
    salt_dim1, salt_dim2, salt_dim3,                                           &
    l_sulpc_so2, l_soot, l_biomass,                                            &
    l_ocff, l_use_biogenic, l_dust, l_nitrate,                                 &
    l_use_seasalt_pm, p_theta_levels, t,                                       &
    so4_ait, so4_acc, soot_new, soot_agd, bmass_new, bmass_agd,                &
    ocff_new, ocff_agd, biogenic, sea_salt_film, sea_salt_jet,                 &
    dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6,          &
    nitr_acc, pm10, pm2p5,                                                     &
    pm10_so4, pm2p5_so4,                                                       &
    pm10_bc, pm2p5_bc,                                                         &
    pm10_bb, pm2p5_bb,                                                         &
    pm10_ocff, pm2p5_ocff,                                                     &
    pm10_soa, pm2p5_soa,                                                       &
    pm10_ss, pm2p5_ss,                                                         &
    conc_dust, conc_dust_surf, pm10_dust, pm2p5_dust,                          &
    pm10_nitr, pm2p5_nitr)

! Copy total dust concs for PWS diagnostics section
  IF ( sf(stashcode_aero_total_dust , stashcode_aerosol_sec) .AND.   &
      (sf(stashcode_pws_dustconc_surf,stashcode_pws_sec) .OR.    &
       sf(stashcode_pws_dustconc_5000,stashcode_pws_sec)) ) THEN
    pws_dustconc_tot(:,:,:) = conc_dust(:,:,:)
    pws_dustconc_surf(:,:) = conc_dust_surf(:,:)
  END IF
END IF


! Call diagnostic routine

CALL diagnostics_aero(                                                         &
  row_length, rows,                                                            &
  n_rows, global_row_length, global_rows,                                      &
  halo_i, halo_j, off_x, off_y, me,                                            &
  n_proc, n_procx, n_procy,                                                    &
  g_rows, g_row_length,                                                        &
  at_extremity,                                                                &
  l_sulpc_so2, l_dms_em,                                                       &
  l_sulpc_dms, l_sulpc_newdms,                                                 &
  l_sulpc_ozone, l_sulpc_nh3,                                                  &
  l_soot,                                                                      &
  msa, nh3_dep,                                                                &
  dms_em_combined,                                                             &
  deltas_dms,                                                                  &
  f_dms_to_so2,                                                                &
  f_dms_to_so4,                                                                &
  f_dms_to_msa,                                                                &
  deltas_dry,                                                                  &
  deltas_wet,                                                                  &
  deltas_wet_o3,                                                               &
  deltas_evap,                                                                 &
  deltas_nucl,                                                                 &
  deltas_diffuse,                                                              &
  deltas_coag,                                                                 &
  deltas_merge,                                                                &
  delta_n_chem,                                                                &
  delta_n_evap,                                                                &
  delta_n_nuc,                                                                 &
  psi,                                                                         &
  pm10, pm2p5,                                                                 &
  pm10_so4, pm2p5_so4,                                                         &
  pm10_bc, pm2p5_bc,                                                           &
  pm10_bb, pm2p5_bb,                                                           &
  pm10_ocff, pm2p5_ocff,                                                       &
  pm10_soa, pm2p5_soa,                                                         &
  pm10_ss, pm2p5_ss,                                                           &
  conc_dust, pm10_dust, pm2p5_dust,                                            &
  pm10_nitr, pm2p5_nitr,                                                       &
  stashwork17)

l_initialise_injection_levels = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE aero_ctl_4A


END MODULE aero_ctl_mod_4A


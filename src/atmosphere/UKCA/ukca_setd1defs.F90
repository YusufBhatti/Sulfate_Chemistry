! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Subroutine to define fields required from D1 by UKCA.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!  Set D1 section and item codes depending on the selected
!  chemistry, aerosol and other options.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
!  Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_setd1defs_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_SETD1DEFS_MOD'

CONTAINS

SUBROUTINE ukca_setd1defs(row_length, rows, n_rows, bl_levels,   &
                          tr_levels ,land_pts, sm_levels,        &
                          ntiles, tr_ukca, all_ntp)

USE ukca_d1_defs
USE ukca_tracer_stash,    ONLY: a_max_ukcavars
USE ukca_nmspec_mod, ONLY: nm_spec
USE ukca_option_mod, ONLY: L_ukca_rado3, L_ukca_radch4,          &
                           L_ukca_mode, L_ukca_chem, L_ukca_dust,&
                           L_ukca_qch4inter, i_ukca_photol,      &
                           L_ukca_arg_act, L_ukca_radn2o,        &
                           L_ukca_ageair, L_ukca_achem,          &
                           L_ukca_aerchem, L_ukca_tropisop,      &
                           L_ukca_trop, l_ukca_prim_moc,         &
                           L_ukca_stratcfc, L_ukca_trophet,      &
                           L_ukca_strattrop, L_ukca_raq,         &
                           L_ukca_strat, jpctr, l_ukca,          &
                           tr_ukca_a,                            &
                           l_ukca_offline, l_ukca_offline_be,    &
                           l_ukca_nr_aqchem,                     &
                           l_ukca_classic_hetchem,               &
                           l_ukca_raqaero, l_ukca_chem_aero
USE ukca_photo_scheme_mod,   ONLY: i_ukca_fastjx
USE ukca_ntp_mod,            ONLY: ntp_type, dim_ntp
USE asad_mod,                ONLY: advt, nadvt, speci
USE dust_parameters_mod,     ONLY: ndiv
USE carbon_options_mod,      ONLY: l_co2_interactive
USE jules_surface_types_mod, ONLY: ntype, npft
USE um_stashcode_mod,        ONLY: stashcode_glomap_sec
USE parkind1,                ONLY: jprb, jpim
USE yomhook,                 ONLY: lhook, dr_hook
USE ereport_mod,             ONLY: ereport
USE umPrintMgr
USE spec_sw_lw, ONLY: sw_spectrum
USE run_aerosol_mod, ONLY:                                       &
    l_sulpc_so2, l_so2_surfem, l_so2_hilem, l_so2_natem,         &
    l_nh3, l_dms_em, l_nh3_em, l_sulpc_dms, l_soot, l_ocff,      &
    l_use_seasalt_sulpc
USE set_rad_steps_mod, ONLY: l_rad_step_prog
USE rad_input_mod,     ONLY: l_use_biogenic,                     &
                             l_use_seasalt_direct,               & 
                             l_use_seasalt_indirect,             &
                             l_use_arclsulp
USE mphys_inputs_mod,  ONLY: l_use_seasalt_autoconv
USE cstash_mod,        ONLY: modl_b, isec_b, item_b, ndiag
USE jules_sea_seaice_mod,  ONLY: l_ctile
USE submodel_mod,          ONLY: submodel_for_sm, atmos_im
USE ozone_inputs_mod,      ONLY: zon_av_ozone
USE missing_data_mod,      ONLY: rmdi, imdi
USE nlsizes_namelist_mod,  ONLY: model_levels

USE errormessagelength_mod,  ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length     ! length of row
INTEGER, INTENT(IN) :: rows           ! number of rows
INTEGER, INTENT(IN) :: n_rows         ! number of rows on v grid
INTEGER, INTENT(IN) :: bl_levels      ! number of levels in BL
INTEGER, INTENT(IN) :: tr_levels      ! number of tracer levels
INTEGER, INTENT(IN) :: land_pts       ! no of land points
INTEGER, INTENT(IN) :: sm_levels      ! no of soil moisture levels
INTEGER, INTENT(IN) :: ntiles         ! no of land tile types
INTEGER, INTENT(IN) :: tr_ukca        ! no of activated tracers
! structure holding non-transported prognostics
TYPE(ntp_type), INTENT(IN) :: all_ntp(dim_ntp) 

INTEGER :: i,j,idiag                  ! counters
INTEGER :: errcode                    ! Variable passed to ereport

CHARACTER (LEN=errormessagelength) :: cmessage        ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SETD1DEFS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Set max fluxdiags to zero initially
nmax_strat_fluxdiags = 0
nmax_mode_diags      = 0
n_nonchem_tracers    = 0

! as lbc_spec and lbc_mmr are passed to emiss_ctl they must be allocated
ALLOCATE(lbc_mmr(n_boundary_vals))
ALLOCATE(lbc_spec(n_boundary_vals))

lbc_mmr(:) = rmdi
!  Species with potential lower boundary conditions (same for all flavours)
lbc_spec =                                                        &
     (/'N2O       ','CF2Cl2    ',                                 &
       'CFCl3     ','MeBr      ',                                 &
       'H2        ','CH4       ',                                 &
       'COS       '/)

IF (L_UKCA_trop) THEN

  ! Standard tropospheric chemistry for B-E solver
  ! ==============================================
  n_chem_emissions = 8
  n_3d_emissions = 1       ! aircraft NOX
  n_aero_tracers = 0
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  n_chem_tracers =  26
  nr_therm       = 102        ! thermal reactions
  nr_phot        = 27         ! photolytic ---"---
  nmax_strat_fluxdiags = n_chem_tracers
  em_chem_spec =                                                 &
  (/'NO        ','CH4       ','CO        ','HCHO      ',         &
    'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',         &
    'NO_aircrft'/)
ELSE IF (L_ukca_tropisop .AND. L_ukca_achem) THEN

  ! Std tropospheric chemistry + MIM with aerosol scheme (N-R)
  ! ==========================================================
  n_chem_emissions = 19    ! 2D emission fields
  n_3d_emissions = 4       ! SO2_nat, BC & OC biomass, aircraft NOX
  n_aero_tracers = 9       ! DMS, SO2... aerosol precursor species
  n_chem_tracers = 51
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  nr_therm       = 113
  nr_phot        = 37
  nmax_strat_fluxdiags = n_chem_tracers
  ! Table refers to emissions, more species are emitted using surrogates
  em_chem_spec =                                                  &
  (/'NO        ','CH4       ','CO        ','HCHO      ',          &
    'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',          &
    'C5H8      ','BC_fossil ','BC_biofuel','OC_fossil ',          &
    'OC_biofuel','Monoterp  ','NVOC      ','SO2_low   ',          &
    'SO2_high  ','NH3       ','DMS       ','SO2_nat   ',          &
    'BC_biomass','OC_biomass','NO_aircrft'/)
ELSE IF (L_ukca_tropisop .AND. .NOT. L_ukca_achem) THEN

  ! Std tropospheric chemistry with MIM isoprene scheme (N-R)
  ! =========================================================
  n_chem_emissions = 9
  n_3d_emissions = 1       ! aircraft NOX
  n_aero_tracers = 0
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  n_chem_tracers = 49
  nr_therm       = 132
  nr_phot        = 35
  nmax_strat_fluxdiags = n_chem_tracers
  em_chem_spec =                                                 &
  (/'NO        ','CH4       ','CO        ','HCHO      ',         &
    'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',         &
    'C5H8      ','NO_aircrft'/)
ELSE IF (L_ukca_aerchem) THEN

  ! Std trop chem with SO2, DMS, NH3, and monoterpene (BE)
  ! ======================================================
  n_chem_emissions = 18       ! Surface/ high-level emissions
  n_3d_emissions = 4          ! SO2_nat, aircraft NOX, OC & BC Biomass
  n_chem_tracers = 26         ! advected chemical tracers
  n_aero_tracers =  7         ! advected aerochem ---"---
  nr_therm       = 137        ! thermal reactions
  nr_phot        = 27         ! photolytic ---"---
  nmax_strat_fluxdiags = n_chem_tracers
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  em_chem_spec =                                                  &
  (/'NO        ','CH4       ','CO        ','HCHO      ',          &
    'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',          &
    'C5H8      ','BC_fossil ','BC_biofuel','OC_fossil ',          &
    'Monoterp  ','MeOH      ','SO2_low   ','SO2_high  ',          &
    'NH3       ','DMS       ','SO2_nat   ','BC_biomass',          &
    'OC_biomass','NO_aircrft'/)
ELSE IF (l_ukca_raq) THEN

  ! Regional air quality chemistry (RAQ), based on STOCHEM
  ! ========================================================
  n_chem_emissions  = 16
  n_3d_emissions    = 1       ! aircraft NOx
  n_chem_tracers    = 40      ! advected chemical tracers
  n_aero_tracers    = 0       ! advected aerochem tracers
  nr_therm          = 192     ! thermal reactions
  nr_phot           = 23      ! photolytic reacs
  nmax_strat_fluxdiags = n_chem_tracers
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  em_chem_spec =                                                  &
    (/'NO        ','CH4       ','CO        ','HCHO      ',        &
      'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',        &
      'C5H8      ','C4H10     ','C2H4      ','C3H6      ',        &
      'TOLUENE   ','oXYLENE   ','CH3OH     ','H2        ',        &
      'NO_aircrft' /)
ELSE IF (l_ukca_raqaero) THEN

! Regional air quality chemistry plus aerosols RAQ-AERO
! ========================================================
  n_chem_emissions  = 20
  n_3d_emissions    = 1       ! aircraft NOx
  n_chem_tracers    = 40      ! advected chemical tracers
  n_aero_tracers    = 8       ! advected aerochem tracers
  nr_therm          = 197     ! thermal reactions
  nr_phot           = 23      ! photolytic reacs 
  nmax_strat_fluxdiags = n_chem_tracers
  
  IF (l_ukca_mode) THEN
    ! BC & OC fossil, biofuel, & biomass emissions
    n_chem_emissions = n_chem_emissions + 6
  END IF
  
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  em_chem_spec(1:21) =                                            &
    (/'NO        ','CH4       ','CO        ','HCHO      ',        & ! 4
      'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',        & ! 8
      'C5H8      ','C4H10     ','C2H4      ','C3H6      ',        & ! 12
      'TOLUENE   ','oXYLENE   ','CH3OH     ','Monoterp  ',        & ! 16
      'SO2_low   ','SO2_high  ','DMS       ','NH3       ',        & ! 20
      'NO_aircrft'                                                & ! 21
      /)

  IF (l_ukca_mode) THEN
    ! The fossil/biofuel pairs will be 2D;
    ! the biomass emissions will also be 2D
    em_chem_spec(22:23) = (/ 'BC_fossil ', 'BC_biofuel' /)
    em_chem_spec(24:25) = (/ 'OC_fossil ', 'OC_biofuel' /)
    em_chem_spec(26:27) = (/ 'OC_biomass', 'BC_biomass' /)
  END IF


ELSE IF (l_ukca_offline_be) THEN

! Offline oxidants scheme with aerosol chemistry
! ==============================================
  n_chem_emissions = 8     ! 2D emission fields
  n_3d_emissions = 3       ! SO2_nat, BC & OC biomass
  n_aero_tracers = 7       ! DMS, SO2... aerosol precursor species
  n_chem_tracers = 0
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  nr_therm       = 11      ! ratb + ratt
  nr_phot        = 0
  nmax_strat_fluxdiags = n_chem_tracers

! Table refers to emissions, more species may be emitted using surrogates
  em_chem_spec =                                                  &
  (/'BC_fossil ','BC_biofuel','OC_fossil ','OC_biofuel',          &
    'Monoterp  ','SO2_low   ','SO2_high  ','DMS       ',          &
    'SO2_nat   ','BC_biomass','OC_biomass'/)

ELSE IF (L_UKCA_strat .OR. L_ukca_strattrop .OR. L_UKCA_stratcfc) &
  THEN

  ! Stratospheric chemistry
  ! =======================
  n_nonchem_tracers  = 1            ! Passive O3 included by default

  IF (L_ukca_strat) THEN
    IF (.NOT. L_ukca_achem) THEN ! NOT using aerosol chemistry
       ! emissions:
      n_chem_emissions = 4
      n_3d_emissions = 1       ! aircraft NOX
      n_aero_tracers =  0
      ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
      em_chem_spec =                                          &
           (/'NO        ','CH4       ',                       &
             'CO        ','HCHO      ',                       &
             'NO_aircrft'/)
      ! tracers and reactions:
      n_chem_tracers  = 38   ! CCMVal !!No H2OS, but does have H2O
      nr_therm       = 135
      nr_phot        = 34
    ELSE ! USING AEROSOL CHEMISTRY
       ! emissions:
      n_chem_emissions = 7       ! em_chem_spec below
      n_3d_emissions   = 2       ! volc SO2 & aircraft NOX
      ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
      em_chem_spec =                                          &
           (/'NO        ','CH4       ',                       &
             'CO        ','HCHO      ',                       &
             'SO2_low   ','SO2_high  ',                       &
             'DMS       ','SO2_nat   ',                       &
             'NO_aircrft'/)
      ! tracers and reactions:
      n_chem_tracers  = 46
      n_aero_tracers  = 8
      nr_therm        = 149
      nr_phot         = 38
    END IF
  ELSE IF (L_ukca_strattrop .AND. .NOT. L_ukca_achem) THEN
    n_chem_emissions = 10
    n_3d_emissions = 1       ! aircraft NOX
    n_aero_tracers =  0
    ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
    em_chem_spec =                                             &
        (/'NO        ','CH4       ','CO        ','HCHO      ', &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ', &
          'C5H8      ','NVOC      ','NO_aircrft'/)
    n_chem_tracers = 71         ! No chem tracers
    nr_therm       = 220        ! thermal reactions
    nr_phot        = 55         ! photolytic (ATA)

  ELSE IF (L_ukca_strattrop .AND. L_ukca_achem) THEN
    n_chem_emissions = 19      ! em_chem_spec below
    n_3d_emissions   = 4       ! BC, OC, volc SO2 & aircraft NOX
    ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
    em_chem_spec =                                             &
        (/'NO        ','CH4       ','CO        ','HCHO      ', &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ', &
          'C5H8      ','BC_fossil ','BC_biofuel','OC_fossil ', &
          'OC_biofuel','Monoterp  ','NVOC      ','SO2_low   ', &
          'SO2_high  ','NH3       ','DMS       ','SO2_nat   ', &
          'BC_biomass','OC_biomass','NO_aircrft'/)
    n_aero_tracers = 12
    n_chem_tracers = 72         ! No chem tracers !!!update 71 to 72 re MSIA inclusion, LER UC 15 Feb 2019
    IF (L_ukca_trophet) THEN
      nr_therm     = 241        ! thermal reactions
    ELSE
      nr_therm     = 239        ! thermal reactions
    END IF
    nr_phot        = 59         ! photolytic (ATA)

  ELSE IF (L_ukca_stratcfc) THEN
    n_chem_tracers = 43
  END IF
  nr_therm       = 102        !
  nr_phot        = 27         !
  nmax_strat_fluxdiags = n_chem_tracers
ELSE IF (l_ukca_offline) THEN

  ! Offline oxidants scheme with aerosol chemistry
  ! ==============================================
  n_chem_emissions = 8     ! 2D emission fields
  n_3d_emissions = 3       ! SO2_nat, BC & OC biomass
  n_aero_tracers = 7       ! DMS, SO2... aerosol precursor species
  n_chem_tracers = 0
  ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
  nr_therm       = 11      ! ratb + ratt, unused ?
  nr_phot        = 0
  nmax_strat_fluxdiags = n_chem_tracers

  ! Table refers to emissions, more species may be emitted using surrogates
  em_chem_spec =                                                  &
  (/'BC_fossil ','BC_biofuel','OC_fossil ','OC_biofuel',          &
    'Monoterp  ','SO2_low   ','SO2_high  ','DMS       ',          &
    'SO2_nat   ','BC_biomass','OC_biomass'/)
ELSE

  ! No chemistry
  ! ============
  n_chem_emissions = 0
  n_chem_tracers   = 0
END IF

IF (L_ukca_ageair) THEN
  ! Include Age of air tracer
  n_nonchem_tracers = n_nonchem_tracers + 1
END IF

IF (L_UKCA_dust) THEN
  n_dust_emissions = 6  ! Use 6 bin dust emissions from UM (Woodward,2001)
  n_dust_tracers   = 0  ! these currently included in n_mode_tracers
  ALLOCATE(em_dust_spec(n_dust_emissions))
  em_dust_spec(1:n_dust_emissions)=                              &
  (/'Dust_div_1','Dust_div_2','Dust_div_3','Dust_div_4',         &
    'Dust_div_5','Dust_div_6'/)
ELSE
  n_dust_emissions = 0
  n_dust_tracers   = 0
END IF

IF (l_ukca_mode) THEN
  IF (.NOT. (l_ukca_aerchem .OR. l_ukca_nr_aqchem .OR.            &
             l_ukca_offline_be .OR. l_ukca_raqaero)) THEN
    cmessage=' l_ukca_aerchem or l_ukca_tropisop or'//            &
             ' l_ukca_offline or l_ukca_offline_be or'//          &
             ' l_ukca_raqaero is needed for MODE'
    errcode=1
    CALL ereport('UKCA_SETD1DEFS',errcode,cmessage)
  END IF
  n_MODE_emissions = 0          ! See aerosol chemistry section
  ! n_mode_tracers is now calculated in UKCA_INIT
  nmax_mode_diags  = 359
ELSE
  n_MODE_emissions = 0
  IF (l_sulpc_so2 .AND. l_use_arclsulp) THEN
     cmessage=' cannot use CLASSIC (l_sulpc_so2) with'//          &
              ' aerosol climatologies (l_use_arclsulp)'
     errcode=1
     CALL ereport('UKCA_SETD1DEFS',errcode,cmessage)
  END IF
END IF

n_use_tracers   = n_chem_tracers + n_mode_tracers +               &
                  n_aero_tracers + n_dust_tracers +               &
                  n_nonchem_tracers
n_use_emissions = n_chem_emissions + n_mode_emissions +           &
                  n_3d_emissions   

n_in_progs     =  40     ! max no of progs reqd other than tracers/ems
n_in_diags0    =   4     ! max no of diags (sect 0) reqd
n_in_diags1    =   4     ! max no of diags (sect 1) reqd
n_in_diags2    =   1     ! max no of diags (sect 2) reqd
n_in_diags3    =  16     ! max no of diags (sect 3) reqd
n_in_diags4    =  10     ! max no of diags (sect 4) reqd
n_in_diags5    =   4     ! max no of diags (sect 5) reqd
n_in_diags8    =   1     ! max no of diags (sect 8) reqd
n_in_diags15   =   1     ! max no of diags (sect 15) reqd
n_in_diags30   =   1     ! max no of diags (sect 30) reqd
n_in_diags38   = nmax_mode_diags + nmax_strat_fluxdiags
!                              ! max no UKCA diags (sect 38) reqd

! If radiative feedback is specified, check that the tracer array addresses
! correspond to the named tracer
IF (L_ukca_rado3) THEN
  IF (nm_spec(i_ukca_grg_o3) /= 'O3        ' .AND.                &
      ANY(advt(:) /= 'O3        ')) THEN
    errcode = i_ukca_grg_o3
    cmessage = ' Tracer address for O3 radiation feedback '//     &
               ' does not correspond with O3'
    WRITE(umMessage,'(A80,I5,2A12)') cmessage,i_ukca_grg_o3,      &
               nm_spec(i_ukca_grg_o3), advt(i_ukca_grg_o3)
    CALL umPrint(umMessage,src='ukca_setd1defs')
    CALL ereport('UKCA_SETD1DEFS',errcode,cmessage)
  END IF
END IF
IF (L_ukca_radch4) THEN
  IF (nm_spec(i_ukca_grg_ch4) /= 'CH4       ' .AND.               &
      ANY(advt(:) /= 'CH4       ')) THEN
    errcode = i_ukca_grg_ch4
    cmessage = ' Tracer address for CH4 radiation feedback '//    &
               ' does not correspond with CH4'
    WRITE(umMessage,'(A80,I5,2A12)') cmessage,i_ukca_grg_ch4,     &
               nm_spec(i_ukca_grg_ch4), advt(i_ukca_grg_ch4)
    CALL umPrint(umMessage,src='ukca_setd1defs')
    CALL ereport('UKCA_SETD1DEFS',errcode,cmessage)
  END IF
END IF
IF (L_ukca_radn2o) THEN
  IF (nm_spec(i_ukca_grg_n2o) /= 'N2O       ' .AND.               &
      ANY(advt(:) /= 'N2O       ')) THEN
    errcode = i_ukca_grg_n2o
    cmessage = ' Tracer address for N2O radiation feedback '//    &
               ' does not correspond with CH4'
    WRITE(umMessage,'(A80,I5,2A12)') cmessage,i_ukca_grg_ch4,     &
               nm_spec(i_ukca_grg_ch4), advt(i_ukca_grg_ch4)
    CALL umPrint(umMessage,src='ukca_setd1defs')
    CALL ereport('UKCA_SETD1DEFS',errcode,cmessage)
  END IF
END IF

!     Specify the section and item codes of prognostics and diagnostics
!     required from D1.  Set array dimensions, but ignore halos as these
!     are set in the call to UKCA_SET_ARRAY_BOUNDS from the halosize
!     array.  Set %prognostic and %required t/f.
!     Other components of UkcaD1Codes are read in from D1 address array
!     in UKCA_MAIN1.

Nukca_D1items = n_use_tracers  +                                  &
                n_dust_emissions + dim_ntp       +                & 
                n_in_progs     + n_in_diags0     +                &
                n_in_diags1    + n_in_diags2     +                &
                n_in_diags3    +                                  &
                n_in_diags4    + n_in_diags5     +                &
                n_in_diags8    + n_in_diags15    +                &
                n_in_diags30   +                                  &
                n_in_diags38

ALLOCATE(UkcaD1Codes(Nukca_D1items))

UkcaD1Codes(:)%section=imdi
UkcaD1Codes(:)%item=imdi
UkcaD1Codes(:)%n_levels=imdi
UkcaD1Codes(:)%address=imdi
UkcaD1Codes(:)%length=imdi
UkcaD1Codes(:)%halo_type=imdi
UkcaD1Codes(:)%grid_type=imdi
UkcaD1Codes(:)%field_type=imdi
UkcaD1Codes(:)%len_dim1=imdi
UkcaD1Codes(:)%len_dim2=imdi
UkcaD1Codes(:)%len_dim3=imdi
UkcaD1Codes(:)%prognostic=.TRUE.
UkcaD1Codes(:)%required=.FALSE.

! If using Newton-Raphson solver, re-order tracers so that they
!  exist in the same order as in ASAD

! Only retain names for active tracers and set required logical

IF (printstatus >= Prstatus_oper) THEN
  WRITE(umMessage,'(A50)') ' UKCA: The following tracers were selected:'
  CALL umPrint(umMessage,src='ukca_setd1defs')
END IF

j = 0

! ----------------------------------------------------------------------
! Non transported prognostics - in UKCA section but not tracers
! ----------------------------------------------------------------------

! Loop over all_ntp array. Assume that always prognostic
! and always have same dimensions as for tracers
DO i=1, dim_ntp

  UkcaD1Codes(j+1)%section    = all_ntp(i)%section
  UkcaD1Codes(j+1)%item       = all_ntp(i)%item
  UkcaD1Codes(j+1)%prognostic = .TRUE.
  UkcaD1Codes(j+1)%required   = all_ntp(i)%l_required
  UkcaD1Codes(j+1)%len_dim1   = row_length
  UkcaD1Codes(j+1)%len_dim2   = rows
  UkcaD1Codes(j+1)%len_dim3   = tr_levels
  UkcaD1Codes(j+1)%name       = all_ntp(i)%varname

  ! Increment j before next entry
  j = j+1

END DO

! ----------------------------------------------------------------------
! Prognostics from section 0
! ----------------------------------------------------------------------
! 

IF (L_UKCA_dust) THEN
  DO i=1,n_dust_emissions
    UkcaD1Codes(j+i)%section    = 3
    UkcaD1Codes(j+i)%item       = 400+i     ! Dust emissions sec 3, items 401-6
    UkcaD1Codes(j+i)%len_dim1   = row_length        
    UkcaD1Codes(j+i)%len_dim2   = rows             
    UkcaD1Codes(j+i)%required   = .TRUE.           
    UkcaD1Codes(j+i)%prognostic = .FALSE.           
  END DO
END IF

j = j + n_dust_emissions

! Prognostic fields

UkcaD1Codes(j+1:j+n_in_progs)%section    = 0
UkcaD1Codes(j+1:j+n_in_progs)%prognostic = .TRUE.
UkcaD1Codes(j+1:j+n_in_progs)%required   = .TRUE.  ! default option
UkcaD1Codes(j+1)%item=4              ! Potential Temperature
UkcaD1Codes(j+1)%len_dim1=row_length
UkcaD1Codes(j+1)%len_dim2=rows
UkcaD1Codes(j+1)%len_dim3=model_levels
UkcaD1Codes(j+2)%item=9              ! Soil Moisture
UkcaD1Codes(j+2)%len_dim1=land_pts
UkcaD1Codes(j+2)%len_dim2=sm_levels
UkcaD1Codes(j+3)%item=10             ! Q
UkcaD1Codes(j+3)%len_dim1=row_length
UkcaD1Codes(j+3)%len_dim2=rows
UkcaD1Codes(j+3)%len_dim3=model_levels
UkcaD1Codes(j+4)%item=12             ! QCF
UkcaD1Codes(j+4)%len_dim1=row_length
UkcaD1Codes(j+4)%len_dim2=rows
UkcaD1Codes(j+4)%len_dim3=model_levels
UkcaD1Codes(j+5)%item=16             ! Conv cloud liquid water path
UkcaD1Codes(j+5)%len_dim1=row_length
UkcaD1Codes(j+5)%len_dim2=rows
UkcaD1Codes(j+6)%item=24             ! Surface temperature
UkcaD1Codes(j+6)%len_dim1=row_length
UkcaD1Codes(j+6)%len_dim2=rows
UkcaD1Codes(j+7)%item=25             ! Boundary layer height
UkcaD1Codes(j+7)%len_dim1=row_length
UkcaD1Codes(j+7)%len_dim2=rows
UkcaD1Codes(j+8)%item=26             ! Roughness length
UkcaD1Codes(j+8)%len_dim1=row_length
UkcaD1Codes(j+8)%len_dim2=rows
UkcaD1Codes(j+9)%item=30            ! Land sea mask
UkcaD1Codes(j+9)%len_dim1=row_length
UkcaD1Codes(j+9)%len_dim2=rows
UkcaD1Codes(j+10)%item=31            ! Sea ice fraction
UkcaD1Codes(j+10)%len_dim1=row_length
UkcaD1Codes(j+10)%len_dim2=rows
UkcaD1Codes(j+11)%item=60            ! Climatological ozone
IF (zon_av_ozone) THEN
  UkcaD1Codes(j+11)%len_dim1=1
ELSE
  UkcaD1Codes(j+11)%len_dim1=row_length
END IF
UkcaD1Codes(j+11)%len_dim2=rows
UkcaD1Codes(j+11)%len_dim3=model_levels
UkcaD1Codes(j+12)%item=96            ! Surface chlorophyll
UkcaD1Codes(j+12)%len_dim1=row_length
UkcaD1Codes(j+12)%len_dim2=rows
IF (.NOT. l_ukca_prim_moc) UkcaD1Codes(j+12)%required=.FALSE.
UkcaD1Codes(j+13)%item=103           ! SO4 Aitken Mode
UkcaD1Codes(j+13)%len_dim1=row_length
UkcaD1Codes(j+13)%len_dim2=rows
UkcaD1Codes(j+13)%len_dim3=tr_levels
IF ((.NOT. l_sulpc_so2) .OR. l_use_arclsulp)                          &
                                    UkcaD1Codes(j+13)%required=.FALSE.
UkcaD1Codes(j+14)%item=104           ! SO4 accumulation Mode
UkcaD1Codes(j+14)%len_dim1=row_length
UkcaD1Codes(j+14)%len_dim2=rows
UkcaD1Codes(j+14)%len_dim3=tr_levels
IF ((.NOT. l_sulpc_so2) .OR. l_use_arclsulp)                          &
                                    UkcaD1Codes(j+14)%required=.FALSE.
UkcaD1Codes(j+15)%item=150           ! W component of wind
UkcaD1Codes(j+15)%len_dim1=row_length
UkcaD1Codes(j+15)%len_dim2=rows
UkcaD1Codes(j+15)%len_dim3=model_levels+1
IF (.NOT. L_ukca_arg_act) UkcaD1Codes(j+15)%required=.FALSE.
UkcaD1Codes(j+16)%item=211           ! Conv cloud amount
UkcaD1Codes(j+16)%len_dim1=row_length
UkcaD1Codes(j+16)%len_dim2=rows
UkcaD1Codes(j+16)%len_dim3=model_levels
UkcaD1Codes(j+17)%item=216           ! Fraction of surface types
UkcaD1Codes(j+17)%len_dim1=land_pts
UkcaD1Codes(j+17)%len_dim2=ntype
UkcaD1Codes(j+18)%item=217           ! LAI of PFTs
UkcaD1Codes(j+18)%len_dim1=land_pts
UkcaD1Codes(j+18)%len_dim2=npft
UkcaD1Codes(j+19)%item=218           ! Canopy heights of PFTs
UkcaD1Codes(j+19)%len_dim1=land_pts
UkcaD1Codes(j+19)%len_dim2=npft
UkcaD1Codes(j+20)%item=229           ! Canopy water content on tiles
UkcaD1Codes(j+20)%len_dim1=land_pts
UkcaD1Codes(j+20)%len_dim2=ntiles
UkcaD1Codes(j+21)%item=233           ! Surface temperature on tiles
UkcaD1Codes(j+21)%len_dim1=land_pts
UkcaD1Codes(j+21)%len_dim2=ntiles
UkcaD1Codes(j+22)%item=234           ! Surface roughness lengths on t
UkcaD1Codes(j+22)%len_dim1=land_pts
UkcaD1Codes(j+22)%len_dim2=ntiles
UkcaD1Codes(j+23)%item=240           ! Snow depth on tiles
UkcaD1Codes(j+23)%len_dim1=land_pts
UkcaD1Codes(j+23)%len_dim2=ntiles
UkcaD1Codes(j+24)%item=253           ! Rho_r2
UkcaD1Codes(j+24)%len_dim1=row_length
UkcaD1Codes(j+24)%len_dim2=rows
UkcaD1Codes(j+24)%len_dim3=model_levels
UkcaD1Codes(j+25)%item=254           ! QCL
UkcaD1Codes(j+25)%len_dim1=row_length
UkcaD1Codes(j+25)%len_dim2=rows
UkcaD1Codes(j+25)%len_dim3=model_levels
UkcaD1Codes(j+26)%item=255           ! Exner pressure on rho levels
UkcaD1Codes(j+26)%len_dim1=row_length
UkcaD1Codes(j+26)%len_dim2=rows
UkcaD1Codes(j+26)%len_dim3=model_levels+1
UkcaD1Codes(j+27)%item=265           ! Area cloud fraction in each la
UkcaD1Codes(j+27)%len_dim1=row_length
UkcaD1Codes(j+27)%len_dim2=rows
UkcaD1Codes(j+27)%len_dim3=model_levels
UkcaD1Codes(j+28)%item=266           ! Bulk Cloud fraction
UkcaD1Codes(j+28)%len_dim1=row_length
UkcaD1Codes(j+28)%len_dim2=rows
UkcaD1Codes(j+28)%len_dim3=model_levels
UkcaD1Codes(j+29)%item=267           ! Cloud Liquid fraction
UkcaD1Codes(j+29)%len_dim1=row_length
UkcaD1Codes(j+29)%len_dim2=rows
UkcaD1Codes(j+29)%len_dim3=model_levels
UkcaD1Codes(j+30)%item=505           ! Land fraction
UkcaD1Codes(j+30)%len_dim1=land_pts
IF (.NOT. l_ctile) UkcaD1Codes(j+30)%required = .FALSE.
UkcaD1Codes(j+31)%item=510           ! Land albedo
UkcaD1Codes(j+31)%len_dim1=row_length
UkcaD1Codes(j+31)%len_dim2=rows
IF (.NOT. l_rad_step_prog .OR. .NOT. l_ctile) THEN
  ! Required only on radiation TS, not available if coastal tiling off
  UkcaD1Codes(j+31)%required = .FALSE.
END IF
UkcaD1Codes(j+32)%item=132           ! Seawater DMS emiss
UkcaD1Codes(j+32)%len_dim1=row_length
UkcaD1Codes(j+32)%len_dim2=rows
IF (l_sulpc_dms .OR. (.NOT. l_ukca_chem_aero))                        &
                                  UkcaD1Codes(j+32)%required = .FALSE.
                               ! Only req if CLASSIC DMS OFF and are 
                               ! using aerosol chemistry

! Aerosol prognostics from CLASSIC to calculate surface area
! available for tropospheric heterogeneous chemistry in UKCA.
! They are required only if they are being modelled and at
! the same time tropospheric heterogeneous chemistry on
! CLASSIC aerosols is used.
!
UkcaD1Codes(j+33)%item=108           ! fresh soot mmr
UkcaD1Codes(j+33)%len_dim1=row_length
UkcaD1Codes(j+33)%len_dim2=rows
UkcaD1Codes(j+33)%len_dim3=tr_levels
!
UkcaD1Codes(j+34)%item=109           ! aged soot mmr
UkcaD1Codes(j+34)%len_dim1=row_length
UkcaD1Codes(j+34)%len_dim2=rows
UkcaD1Codes(j+34)%len_dim3=tr_levels
!
IF (.NOT. l_soot .OR. .NOT. l_ukca_classic_hetchem) THEN
  UkcaD1Codes(j+33:j+34)%required=.FALSE.
END IF
!
UkcaD1Codes(j+35)%item=114           ! fresh OCFF
UkcaD1Codes(j+35)%len_dim1=row_length
UkcaD1Codes(j+35)%len_dim2=rows
UkcaD1Codes(j+35)%len_dim3=tr_levels
!
UkcaD1Codes(j+36)%item=115           ! aged OCFF
UkcaD1Codes(j+36)%len_dim1=row_length
UkcaD1Codes(j+36)%len_dim2=rows
UkcaD1Codes(j+36)%len_dim3=tr_levels
!
IF (.NOT. l_ocff .OR. .NOT. l_ukca_classic_hetchem) THEN
  UkcaD1Codes(j+35:j+36)%required=.FALSE.
END IF
!
UkcaD1Codes(j+37)%item=252           ! CO2 MMR
UkcaD1Codes(j+37)%len_dim1=row_length
UkcaD1Codes(j+37)%len_dim2=rows
UkcaD1Codes(j+37)%len_dim3=tr_levels
!
IF (.NOT. (l_co2_interactive .AND. ANY(speci(:) == 'CO2       '))) THEN
  UkcaD1Codes(j+37)%required=.FALSE.
END IF
!
UkcaD1Codes(j+38)%item=351           ! biogenic
UkcaD1Codes(j+38)%len_dim1=row_length
UkcaD1Codes(j+38)%len_dim2=rows
UkcaD1Codes(j+38)%len_dim3=tr_levels
!
IF (.NOT. l_use_biogenic .OR. .NOT. l_ukca_classic_hetchem) THEN
  UkcaD1Codes(j+38)%required=.FALSE.
END IF

UkcaD1Codes(j+39)%item=359           ! SO4 Accumulation Mode climatology
UkcaD1Codes(j+39)%len_dim1=row_length
UkcaD1Codes(j+39)%len_dim2=rows
UkcaD1Codes(j+39)%len_dim3=tr_levels
IF (.NOT. l_use_arclsulp) UkcaD1Codes(j+39)%required=.FALSE.
UkcaD1Codes(j+40)%item=360           ! SO4 Aitken Mode climatology
UkcaD1Codes(j+40)%len_dim1=row_length
UkcaD1Codes(j+40)%len_dim2=rows
UkcaD1Codes(j+40)%len_dim3=tr_levels
IF (.NOT. l_use_arclsulp) UkcaD1Codes(j+40)%required=.FALSE.

! Diagnostics from section zero
j = j + n_in_progs

UkcaD1Codes(j+1:j+n_in_diags0)%section    = 0
UkcaD1Codes(j+1:j+n_in_diags0)%prognostic = .FALSE.
UkcaD1Codes(j+1:j+n_in_diags0)%required   = .TRUE.
UkcaD1Codes(j+1)%item=406           ! Exner Press on theta levels
UkcaD1Codes(j+1)%len_dim1=row_length
UkcaD1Codes(j+1)%len_dim2=rows
UkcaD1Codes(j+1)%len_dim3=model_levels
UkcaD1Codes(j+2)%item=407           ! P on Rho Levels
UkcaD1Codes(j+2)%len_dim1=row_length
UkcaD1Codes(j+2)%len_dim2=rows
UkcaD1Codes(j+2)%len_dim3=model_levels
UkcaD1Codes(j+3)%item=408           ! P on Theta Levels
UkcaD1Codes(j+3)%len_dim1=row_length
UkcaD1Codes(j+3)%len_dim2=rows
UkcaD1Codes(j+3)%len_dim3=model_levels
UkcaD1Codes(j+4)%item=409           ! Surface Pressure
UkcaD1Codes(j+4)%len_dim1=row_length
UkcaD1Codes(j+4)%len_dim2=rows

! Diagnostic items from section 1 (SW radiation)
j = j + n_in_diags0

UkcaD1Codes(j+1:j+n_in_diags1)%section    = 1
UkcaD1Codes(j+1:j+n_in_diags1)%prognostic = .FALSE.   ! always needed
UkcaD1Codes(j+1:j+n_in_diags1)%required   = .TRUE.    ! default option
UkcaD1Codes(j+1)%item=201           ! Net downward surface SW flux
UkcaD1Codes(j+1)%len_dim1=row_length
UkcaD1Codes(j+1)%len_dim2=rows
UkcaD1Codes(j+2)%item=235           ! Total downward surface SW flux
UkcaD1Codes(j+2)%len_dim1=row_length
UkcaD1Codes(j+2)%len_dim2=rows
!
! Sea-salt aerosol number is only required if two conditions are met at the
! same time: heterogeneous chemistry on CLASSIC aerosols is used and
! sea-salt aerosol is being modelled.
IF ( (.NOT. l_use_seasalt_autoconv .AND. .NOT. l_use_seasalt_sulpc   .AND. &
      .NOT. l_use_seasalt_indirect .AND. .NOT. l_use_seasalt_direct) .OR.  &
     (.NOT. l_ukca_classic_hetchem) ) THEN
  UkcaD1Codes(j+3:j+4)%required = .FALSE.
END IF
UkcaD1Codes(j+3)%item = 247    ! Film mode sea-salt aerosol number
UkcaD1Codes(j+4)%item = 248    ! Jet  mode sea-salt aerosol number
UkcaD1Codes(j+3:j+4)%len_dim1 = row_length
UkcaD1Codes(j+3:j+4)%len_dim2 = rows
UkcaD1Codes(j+3:j+4)%len_dim3 = model_levels
!
! Diagnostic items from section 2 (LW radiation)
j = j + n_in_diags1
UkcaD1Codes(j+1:j+n_in_diags2)%section    = 2
UkcaD1Codes(j+1:j+n_in_diags2)%prognostic = .FALSE.
UkcaD1Codes(j+1)%item=284                     ! Sulphate optical depth
UkcaD1Codes(j+1)%len_dim1=row_length
UkcaD1Codes(j+1)%len_dim2=rows
UkcaD1Codes(j+1)%len_dim3=sw_spectrum(1)%basic%n_band
! sw wave band pseudo levels
WRITE(umMessage,'(a,i3)') 'n_band: ', sw_spectrum(1)%basic%n_band
CALL umPrint(umMessage,src='ukca_setd1defs')
! Do not actually use sulphate AOD. If needed have to add checks so
! only used on LW timesteps. Also possibly issues about the
! number of pseudo levels
UkcaD1Codes(j+1)%required=.FALSE.

! Diagnostic variables in section 3 (Boundary Layer)
j = j + n_in_diags2

UkcaD1Codes(j+1:j+n_in_diags3)%section    = 3
UkcaD1Codes(j+1:j+n_in_diags3)%prognostic = .FALSE.
UkcaD1Codes(j+1:j+n_in_diags3)%required   = .TRUE.
UkcaD1Codes(j+1)%item=217           ! Surface heat flux
UkcaD1Codes(j+1)%len_dim1=row_length
UkcaD1Codes(j+1)%len_dim2=rows
UkcaD1Codes(j+2)%item=462            ! Stomatal conductance
UkcaD1Codes(j+2)%len_dim1=row_length
UkcaD1Codes(j+2)%len_dim2=rows
UkcaD1Codes(j+2)%len_dim3=npft
UkcaD1Codes(j+3)%item=465            ! Surface friction velocity
UkcaD1Codes(j+3)%len_dim1=row_length
UkcaD1Codes(j+3)%len_dim2=rows
UkcaD1Codes(j+4)%item=60            ! rhokh_mix
UkcaD1Codes(j+4)%len_dim1=row_length
UkcaD1Codes(j+4)%len_dim2=rows
Ukcad1codes(j+4)%len_dim3=bl_levels
UkcaD1Codes(j+5)%item=64            ! dtrdz_charney_grid
UkcaD1Codes(j+5)%len_dim1=row_length
UkcaD1Codes(j+5)%len_dim2=rows
Ukcad1codes(j+5)%len_dim3=bl_levels
UkcaD1Codes(j+6)%item=65            ! kent
UkcaD1Codes(j+6)%len_dim1=row_length
UkcaD1Codes(j+6)%len_dim2=rows
UkcaD1Codes(j+7)%item=66            ! we_lim
UkcaD1Codes(j+7)%len_dim1=row_length
UkcaD1Codes(j+7)%len_dim2=rows
Ukcad1codes(j+7)%len_dim3=npft
UkcaD1Codes(j+8)%item=67            ! t_frac
UkcaD1Codes(j+8)%len_dim1=row_length
UkcaD1Codes(j+8)%len_dim2=rows
Ukcad1codes(j+8)%len_dim3=npft
UkcaD1Codes(j+9)%item=68            ! zrzi
UkcaD1Codes(j+9)%len_dim1=row_length
UkcaD1Codes(j+9)%len_dim2=rows
Ukcad1codes(j+9)%len_dim3=npft
UkcaD1Codes(j+10)%item=69            ! kent_dsc
UkcaD1Codes(j+10)%len_dim1=row_length
UkcaD1Codes(j+10)%len_dim2=rows
UkcaD1Codes(j+11)%item=70            ! we_lim_dsc
UkcaD1Codes(j+11)%len_dim1=row_length
UkcaD1Codes(j+11)%len_dim2=rows
Ukcad1codes(j+11)%len_dim3=npft
UkcaD1Codes(j+12)%item=71            ! t_frac_dsc
UkcaD1Codes(j+12)%len_dim1=row_length
UkcaD1Codes(j+12)%len_dim2=rows
Ukcad1codes(j+12)%len_dim3=npft
UkcaD1Codes(j+13)%item=72            ! zrzi_dsc
UkcaD1Codes(j+13)%len_dim1=row_length
UkcaD1Codes(j+13)%len_dim2=rows
Ukcad1codes(j+13)%len_dim3=npft
UkcaD1Codes(j+14)%item=73            ! zhsc
UkcaD1Codes(j+14)%len_dim1=row_length
UkcaD1Codes(j+14)%len_dim2=rows
UkcaD1Codes(j+15)%item=230           ! 10 m wind speed on C grid
UkcaD1Codes(j+15)%len_dim1=row_length
UkcaD1Codes(j+15)%len_dim2=rows
UkcaD1Codes(j+16)%item=473           ! Turbulent Kinetic Energy
UkcaD1Codes(j+16)%len_dim1=row_length
UkcaD1Codes(j+16)%len_dim2=rows
UkcaD1Codes(j+16)%len_dim3=bl_levels
IF (.NOT. l_ukca_arg_act) UkcaD1Codes(j+16)%required=.FALSE.

! Diagnostic variables in section 4 (LS Precipitation)
! We currently do not use 4.227 so set required to .FALSE.
j = j + n_in_diags3

UkcaD1Codes(j+1:j+n_in_diags4)%section    = 4
UkcaD1Codes(j+1:j+n_in_diags4)%prognostic = .FALSE.
UkcaD1Codes(j+1:j+n_in_diags4)%required   = .TRUE.
UkcaD1Codes(j+1)%item=205           ! Cloud Liquid Water after LS
UkcaD1Codes(j+1)%len_dim1=row_length
UkcaD1Codes(j+1)%len_dim2=rows
UkcaD1Codes(j+1)%len_dim3=model_levels
UkcaD1Codes(j+2)%item=222           ! Rainfall rate out of model levs
UkcaD1Codes(j+2)%len_dim1=row_length
UkcaD1Codes(j+2)%len_dim2=rows
UkcaD1Codes(j+2)%len_dim3=model_levels
UkcaD1Codes(j+3)%item=223           ! Snowfall rate out of model levs
UkcaD1Codes(j+3)%len_dim1=row_length
UkcaD1Codes(j+3)%len_dim2=rows
UkcaD1Codes(j+3)%len_dim3=model_levels
UkcaD1Codes(j+4)%item=227           ! Rain fraction (3C ppn only)
UkcaD1Codes(j+4)%len_dim1=row_length
UkcaD1Codes(j+4)%len_dim2=rows
UkcaD1Codes(j+4)%len_dim3=model_levels
UkcaD1Codes(j+4)%required=.FALSE.
UkcaD1Codes(j+5)%item=247           ! Riming rate of ice crystals
UkcaD1Codes(j+5)%len_dim1=row_length
UkcaD1Codes(j+5)%len_dim2=rows
UkcaD1Codes(j+5)%len_dim3=model_levels
IF (.NOT. l_ukca_mode) UkcaD1Codes(j+5)%required=.FALSE.
UkcaD1Codes(j+6)%item=248           ! Riming rate of ice aggregates
UkcaD1Codes(j+6)%len_dim1=row_length
UkcaD1Codes(j+6)%len_dim2=rows
UkcaD1Codes(j+6)%len_dim3=model_levels
IF (.NOT. l_ukca_mode) UkcaD1Codes(j+6)%required=.FALSE.
UkcaD1Codes(j+7)%item=253           ! Melting rate of ice crystals
UkcaD1Codes(j+7)%len_dim1=row_length
UkcaD1Codes(j+7)%len_dim2=rows
UkcaD1Codes(j+7)%len_dim3=model_levels
UkcaD1Codes(j+8)%item=254           ! Melting rate of snow
UkcaD1Codes(j+8)%len_dim1=row_length
UkcaD1Codes(j+8)%len_dim2=rows
UkcaD1Codes(j+8)%len_dim3=model_levels
UkcaD1Codes(j+9)%item=257           ! Rain autoconversion rate
UkcaD1Codes(j+9)%len_dim1=row_length
UkcaD1Codes(j+9)%len_dim2=rows
UkcaD1Codes(j+9)%len_dim3=model_levels
UkcaD1Codes(j+10)%item=258           ! Rain accretion rate
UkcaD1Codes(j+10)%len_dim1=row_length
UkcaD1Codes(j+10)%len_dim2=rows
UkcaD1Codes(j+10)%len_dim3=model_levels

! Diagnostic variables in section 5 (Convection)
j = j + n_in_diags4

UkcaD1Codes(j+1:j+n_in_diags5)%section    = 5
UkcaD1Codes(j+1:j+n_in_diags5)%prognostic = .FALSE.
UkcaD1Codes(j+1:j+n_in_diags5)%required   = .TRUE.
UkcaD1Codes(j+1)%item=227           ! 3D Convective rainfall rate
UkcaD1Codes(j+1)%len_dim1=row_length
UkcaD1Codes(j+1)%len_dim2=rows
UkcaD1Codes(j+1)%len_dim3=model_levels
UkcaD1Codes(j+2)%item=228           ! 3D Convective snowfall rate
UkcaD1Codes(j+2)%len_dim1=row_length
UkcaD1Codes(j+2)%len_dim2=rows
UkcaD1Codes(j+2)%len_dim3=model_levels
! REMOVED FROM PROGNOSTICS (00014 & 00015) AND NOW TAKE DIAGNOSTICS 05218 & 05219
UkcaD1Codes(j+3)%item=218             ! Conv cloud base level
UkcaD1Codes(j+3)%len_dim1=row_length
UkcaD1Codes(j+3)%len_dim2=rows
UkcaD1Codes(j+4)%item=219             ! Conv cloud top level
UkcaD1Codes(j+4)%len_dim1=row_length
UkcaD1Codes(j+4)%len_dim2=rows


! Diagnostic variables in section 8 (Hydrology)
j = j + n_in_diags5

UkcaD1Codes(j+1)%section    = 8
UkcaD1Codes(j+1)%prognostic = .FALSE.
UkcaD1Codes(j+1)%required   = L_ukca_qch4inter
UkcaD1Codes(j+1)%item       = 242        ! CH4 wetland flux
UkcaD1Codes(j+1)%len_dim1   = row_length
UkcaD1Codes(j+1)%len_dim2   = rows

! Diagnostic variables in section 15 (Processed Climate Diagnostics)
j = j + n_in_diags8

UkcaD1Codes(j+1)%section    = 15
UkcaD1Codes(j+1)%prognostic = .FALSE.
UkcaD1Codes(j+1)%required   = .TRUE.
UkcaD1Codes(j+1)%item       = 218           ! PV on theta levels
UkcaD1Codes(j+1)%len_dim1   = row_length
UkcaD1Codes(j+1)%len_dim2   = rows
UkcaD1Codes(j+1)%len_dim3   = model_levels

! Diagnostic variables in section 30 (Processed dynamics diagnostics)
j = j + n_in_diags15

! Take tropopause height here. Needed only for volcanic SO2 emissions into
! the stratosphere. (Always required as in call to emission_ctl)
UkcaD1Codes(j+1)%section = 30
UkcaD1Codes(j+1)%prognostic = .FALSE.
UkcaD1Codes(j+1)%required   = .TRUE.
UkcaD1Codes(j+1)%item=453           ! Tropopause height
UkcaD1Codes(j+1)%len_dim1=row_length
UkcaD1Codes(j+1)%len_dim2=rows

! Diagnostic variables for stratospheric fluxes in section 38
j = j + n_in_diags30
istrat_first = j+1
istrat_last  = j+nmax_strat_fluxdiags

n_strat_fluxdiags = 0
DO i=1,nmax_strat_fluxdiags
  DO idiag=1,ndiag     ! search for stash requests
    IF (modl_b(idiag) == submodel_for_sm(atmos_im) .AND.         &
        isec_b(idiag) == stashcode_glomap_sec .AND.              &
        item_b(idiag) == item1_stratflux+i-1) THEN
      n_strat_fluxdiags = n_strat_fluxdiags + 1
      UkcaD1Codes(j+i)%section    = UKCA_sect
      UkcaD1Codes(j+i)%item       = item_b(idiag)
      UkcaD1Codes(j+i)%len_dim1   = row_length
      UkcaD1Codes(j+i)%len_dim2   = rows
      UkcaD1Codes(j+i)%len_dim3   = model_levels
      UkcaD1Codes(j+i)%required   = .FALSE.
      UkcaD1Codes(j+i)%prognostic = .FALSE.
      IF (.NOT. (L_ukca_stratflux)) L_ukca_stratflux = .TRUE.
      EXIT
    END IF
  END DO
END DO

! Diagnostic variables for UKCA-MODE in section 38
j = j + nmax_strat_fluxdiags
imode_first = j+1
imode_last  = j+nmax_mode_diags

IF (L_ukca_mode) THEN
  n_mode_diags = 0
  DO i=1,nmax_mode_diags
    DO idiag=1,ndiag     ! Search for stash requests. Exclude CMIP6 items 
                         ! 485-545 because ukca_mode_diags handles them.
      IF (modl_b(idiag) == submodel_for_sm(atmos_im) .AND.        &
        isec_b(idiag) == stashcode_glomap_sec .AND.               &
        item_b(idiag) == item1_mode_diags+i-1 .AND.               &
        (item_b(idiag) < 485 .OR. item_b(idiag) > 545)  ) THEN
        n_mode_diags = n_mode_diags + 1
        UkcaD1Codes(j+i)%section    = stashcode_glomap_sec
        UkcaD1Codes(j+i)%item       = item_b(idiag)
        UkcaD1Codes(j+i)%len_dim1   = row_length
        UkcaD1Codes(j+i)%len_dim2   = rows
        UkcaD1Codes(j+i)%len_dim3   = model_levels
        UkcaD1Codes(j+i)%required   = .FALSE.
        UkcaD1Codes(j+i)%prognostic = .FALSE.
        IF (.NOT. (L_ukca_mode_diags)) L_ukca_mode_diags = .TRUE.
        EXIT
      END IF
    END DO
  END DO
END IF

! None of the fields apart from tracers are required if this  
! run does not use chemistry (e.g. Age-air only)
IF ( .NOT. l_ukca_chem )                             &
  UkcaD1Codes(n_use_tracers+1:Nukca_D1items)%required = .FALSE.

IF ( printstatus >= prstatus_diag ) THEN
  WRITE(umMessage,'(A2)') ' '
  CALL umPrint(umMessage,src='ukca_setd1defs')
  WRITE(umMessage,'(A20,I5)') 'n_use_tracers: ', n_use_tracers
  CALL umPrint(umMessage,src='ukca_setd1defs')
  WRITE(umMessage,'(A20,I5)') 'n_use_emissions: ', n_use_emissions
  CALL umPrint(umMessage,src='ukca_setd1defs')
  WRITE(umMessage,'(A,I5)') 'Total no of items required = ',Nukca_D1items
  CALL umPrint(umMessage,src='ukca_setd1defs')
  WRITE(umMessage,'(A2)') ' '
  CALL umPrint(umMessage,src='ukca_setd1defs')
  WRITE(umMessage,'(A)') ' UKCA: UkcaD1Codes required items:'
  CALL umPrint(umMessage,src='ukca_setd1defs')
  WRITE(umMessage,'(A)') 'section,item,prognostic,required,dim1,dim2,dim3'
  CALL umPrint(umMessage,src='ukca_setd1defs')
  DO i=1,Nukca_D1items
    IF (UkcaD1Codes(i)%required) THEN
      WRITE(umMessage,'(3I5,2L7,I5,2I8)') i,                        &
                    UkcaD1Codes(i)%section,                         &
                    UkcaD1Codes(i)%item,                            &
                    UkcaD1Codes(i)%prognostic,                      &
                    UkcaD1Codes(i)%required,                        &
                    UkcaD1Codes(i)%len_dim1,                        &
                    UkcaD1Codes(i)%len_dim2,                        &
                    UkcaD1Codes(i)%len_dim3
      CALL umPrint(umMessage,src='ukca_setd1defs')
    END IF
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_setd1defs
END MODULE ukca_setd1defs_mod

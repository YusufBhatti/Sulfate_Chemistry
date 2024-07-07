! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module to hold all UKCA variables in RUN_UKCA
!          namelist
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
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
MODULE ukca_option_mod

USE missing_data_mod,      ONLY: rmdi, imdi
USE ukca_photo_scheme_mod, ONLY: i_ukca_nophot
USE ukca_tracer_stash,     ONLY: a_max_ukcavars

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE filenamelength_mod, ONLY: filenamelength
USE errormessagelength_mod, ONLY: errormessagelength


IMPLICIT NONE


PRIVATE :: a_max_ukcavars

! Declarations for UKCA sub-model
! -----------------------------------------------------------------------------
! Namelist items

LOGICAL :: l_ukca           =.FALSE. ! True when UKCA is switched on
LOGICAL :: l_ukca_aie1      =.FALSE. ! True when 1st aerosol ind effect required
LOGICAL :: l_ukca_aie2      =.FALSE. ! True when 2nd aerosol ind effect required
LOGICAL :: l_ukca_chem_plev =.FALSE. ! True when section 34 on pressure levels
                                     ! required 
LOGICAL :: l_ukca_asad_plev =.FALSE. ! True when section 50 on pressure levels
                                     ! required

! Main chemistry namelist inputs:
INTEGER :: i_ukca_chem         = 0      ! chemistry scheme to use
LOGICAL :: l_ukca_chem_aero    =.FALSE. ! add aerosol precursors to chemistry
LOGICAL :: l_ukca_trophet      =.FALSE. ! T for tropospheric heterogeneous chem
LOGICAL :: l_ukca_mode         =.FALSE. ! True for UKCA-MODE aerosol scheme
LOGICAL :: l_ukca_dust         =.FALSE. ! True for UKCA-mode dust aerosol 
LOGICAL :: l_ukca_qch4inter    =.FALSE. ! True for interact wetland CH4 ems
LOGICAL :: l_ukca_het_psc      =.FALSE. ! True for Het/PSC chemistry
LOGICAL :: l_ukca_limit_nat    =.FALSE. ! True for limiting NatPSC formation
                                        ! below specified height
LOGICAL :: l_ukca_sa_clim      =.FALSE. ! True to use SPARC surface area density
LOGICAL :: l_ukca_h2o_feedback =.FALSE. ! True for H2O feedback from chem
LOGICAL :: l_ukca_rado3        =.FALSE. ! T when using UKCA O3 in radiation
LOGICAL :: l_ukca_radch4       =.FALSE. ! T when using UKCA CH4 in radiation
LOGICAL :: l_ukca_radn2o       =.FALSE. ! T when using UKCA N2O in radiation
LOGICAL :: l_ukca_radf11       =.FALSE. ! T when using UKCA CFC-11 in radn
LOGICAL :: l_ukca_radf12       =.FALSE. ! T when using UKCA CFC-12 in radn
LOGICAL :: l_ukca_radf113      =.FALSE. ! T when using UKCA CFC-113 in radn
LOGICAL :: l_ukca_radf22       =.FALSE. ! T when using UKCA HCFC-22 in radn
LOGICAL :: l_ukca_radaer       =.FALSE. ! Radiative effects of UKCA aerosols
LOGICAL :: l_ukca_radaer_sustrat =.FALSE. ! Use H2SO4 for stratospheric sulphate
LOGICAL :: l_ukca_intdd        =.FALSE. ! T when using interact dry deposition
LOGICAL :: l_ukca_prescribech4 =.FALSE. ! T when prescribing surface ch4
LOGICAL :: l_ukca_set_trace_gases =.FALSE. ! T to use UM values for fCO2 etc
LOGICAL :: l_ukca_use_background_aerosol =.FALSE. ! use bg aerosol climatology

! T to pass columns to ASAD rather than theta_field
LOGICAL :: l_ukca_asad_columns =.FALSE. 

! T to use LOG(p) to distribute lightning NOx in the vertical
LOGICAL :: l_ukca_linox_scaling =.FALSE. 

! T to enable additional print statements to debug asad chemistry solver
LOGICAL :: l_ukca_debug_asad = .FALSE.
INTEGER :: chem_timestep = imdi         ! Chemical timestep in seconds for N-R
                                        ! and Offline oxidant schemes
INTEGER :: dts0 = 300                   ! Default Backward Euler timestep
INTEGER :: nit  = 8                     ! Number of iterations of BE Solver

INTEGER :: i_ukca_photol = i_ukca_nophot ! Photolysis scheme to use

! Directory pathname for 2d photolysis rates
CHARACTER (LEN=filenamelength) :: phot2d_dir  = 'phot2d dir is unset'
! Dir pathname for 2d upper boundary data
CHARACTER (LEN=filenamelength) :: strat2d_dir = 'strat2d dir is unset'

INTEGER :: fastjx_numwl = imdi        ! No. of wavelengths to use (8, 12, 18)
INTEGER :: fastjx_mode  = imdi        ! 1 = use just 2D above prescutoff,
                                      ! 2 = merge, 3 = just fastjx)
REAL :: fastjx_prescutoff             ! Press for 2D stratospheric photolysis
CHARACTER(LEN=filenamelength) :: jvspec_file ='jvspec file is unset' 
                                      ! FastJX spectral file
CHARACTER(LEN=filenamelength) :: jvscat_file ='jvscat file is unset' 
                                      ! FastJX scatter file
CHARACTER(LEN=filenamelength) :: jvsolar_file ='jvsolar file is unset' 
                                      ! FastJX solar file
CHARACTER(LEN=filenamelength) :: jvspec_dir  ='jvspec dir is unset'  
                                      ! Dir for jvspec file

! Dir for stratospheric aerosol file
CHARACTER(LEN=filenamelength) :: dir_strat_aer  = 'dir_strat_aer is unset'
! File for stratospheric aerosol file
CHARACTER(LEN=filenamelength) :: file_strat_aer = 'file_strat_aer is unset'

! Switch to choose scheme to use for interactive sea-air exchange of DMS
INTEGER :: i_ukca_dms_flux = imdi
                                 ! 1=LissMerl; 2=Wannin, 3=Nightingale

! UKCA_MODE control features:
LOGICAL :: l_ukca_primsu    =.FALSE. ! T for primary sulphate aerosol emissions
LOGICAL :: l_ukca_primss    =.FALSE. ! T for primary sea-salt aerosol emissions
LOGICAL :: l_ukca_primbcoc  =.FALSE. ! T for primary BC/OC aerosol emissions
LOGICAL :: l_ukca_prim_moc  =.FALSE. ! T for primary marine OC aerosol emissions
LOGICAL :: l_ukca_primdu    =.FALSE. ! T for primary dust aerosol emissions
LOGICAL :: l_ukca_use_2dtop =.FALSE. ! T for using 2-D top boundary files
LOGICAL :: l_bcoc_ff        =.FALSE. ! T for primary fossil fuel BC/OC emiss.
LOGICAL :: l_bcoc_bf        =.FALSE. ! T for primary biofuel BC/OC emissions
LOGICAL :: l_bcoc_bm        =.FALSE. ! T for primary biomass BC/OC emissions
LOGICAL :: l_ukca_scale_biom_aer_ems = .FALSE. ! Apply scaling factor to
                                     ! biomass burning BC/OC aerosol emissions
LOGICAL :: l_ukca_scale_seadms_ems = .FALSE. ! Apply scaling to marine DMS 
                                     ! emissions.
LOGICAL :: l_mode_bhn_on    =.TRUE.  ! T for binary sulphate nucleation
LOGICAL :: l_mode_bln_on    =.TRUE.  ! T for BL sulphate nucleation
LOGICAL :: l_ukca_arg_act   =.FALSE. ! T when using AR&G aerosol activation
LOGICAL :: l_ukca_sfix      =.FALSE. ! T for diagnosing UKCA CCN at
                                     ! fixed supersaturation
! These are switches for the UKCA-iBVOC coupling (in ukca_emission_ctl) 
LOGICAL :: l_ukca_ibvoc     =.FALSE. ! True for interactive bVOC emissions 

LOGICAL :: l_ukca_scale_soa_yield = .FALSE. ! Apply scaling factor to SOA 
                                            ! production from monoterpene
INTEGER :: i_mode_setup     = imdi     ! Defines MODE aerosol scheme
INTEGER :: i_mode_nzts      = imdi     ! No. of substeps for nucleation/
                                       ! sedimentation
INTEGER :: i_mode_bln_param_method = 1 ! 1=activ; 2=kinetc; 3=PNAS/Metzer
                                       ! maps to IBLN in GLOMAP

! Not included in namelist at present:
LOGICAL :: l_ukca_plume_scav   =.FALSE. ! use plume scavenging for aerosol
                                        ! tracers
LOGICAL :: l_ukca_conserve_h = .FALSE.  ! Include hydrogen conservation when
                                        ! l_ukca_h2o_feedback is true
INTEGER :: i_mode_ss_scheme = 1         ! Defines sea-salt emission scheme.
INTEGER :: i_mode_nucscav   = 3         ! Choice of nucl. scavenging co-effs:
                                        ! 1=original, 2=ECHAM5-HAM 
                                        ! 3=as(1) but no scav of modes 6&7

REAL :: mode_parfrac         = rmdi ! Fraction of SO2 emissions as aerosol(%)
REAL :: mode_aitsol_cvscav   = rmdi ! Plume scavenging fraction for AITSOL
REAL :: mode_activation_dryr = rmdi ! Activation dry radius in nm
REAL :: mode_incld_so2_rfrac = rmdi 
! fraction of in-cloud oxidised SO2 removed by precipitation 
REAL :: biom_aer_ems_scaling = rmdi ! Biomass-burning aerosol emissions scaling
REAL :: soa_yield_scaling    = rmdi ! Monoterpene SOA yield scaling factor
REAL :: ukca_MeBrMMR         = rmdi ! UKCA trace gas mixing value
REAL :: ukca_MeClMMR         = rmdi ! UKCA trace gas mixing value
REAL :: ukca_CH2Br2MMR       = rmdi ! UKCA trace gas mixing value
REAL :: ukca_H2MMR           = rmdi ! UKCA trace gas mixing value
REAL :: ukca_N2MMR           = rmdi ! UKCA trace gas mixing value
REAL :: ukca_CFC115MMR       = rmdi ! UKCA trace gas mixing value
REAL :: ukca_CCl4MMR         = rmdi ! UKCA trace gas mixing value
REAL :: ukca_MeCCl3MMR       = rmdi ! UKCA trace gas mixing value
REAL :: ukca_HCFC141bMMR     = rmdi ! UKCA trace gas mixing value
REAL :: ukca_HCFC142bMMR     = rmdi ! UKCA trace gas mixing value
REAL :: ukca_H1211MMR        = rmdi ! UKCA trace gas mixing value
REAL :: ukca_H1202MMR        = rmdi ! UKCA trace gas mixing value
REAL :: ukca_H1301MMR        = rmdi ! UKCA trace gas mixing value
REAL :: ukca_H2402MMR        = rmdi ! UKCA trace gas mixing value
REAL :: ukca_COSMMR          = rmdi ! UKCA trace gas mixing value

! Variables for new UKCA emission system (based on NetCDF input files)
INTEGER, PARAMETER :: nr_cdf_files      = 100     ! Max nr of NetCDF files
                                                  ! allowed in name list
! path of emiss files
CHARACTER(LEN=filenamelength) :: ukca_em_dir = 'ukca_em_dir is unset' 
CHARACTER(LEN=filenamelength)  :: ukca_em_files(nr_cdf_files) = (/            &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset', 'ukca_em_files is unset', 'ukca_em_files is unset', &
'ukca_em_files is unset'                                                      &
/)
! Names of emission files

! Information on netCDF files for offline oxidants
INTEGER, PARAMETER :: max_offline_files = 5            ! Max no of NetCDF files
CHARACTER(LEN=filenamelength) :: ukca_offline_dir = 'ukca_em_dir is unset' 
                                                       ! directory
CHARACTER(LEN=filenamelength)  :: ukca_offline_files(max_offline_files) = (/  &
 'ukca_offline_files is unset', 'ukca_offline_files is unset',                &
 'ukca_offline_files is unset', 'ukca_offline_files is unset',                &
 'ukca_offline_files is unset'/)
! Names of offline oxidants files

! control of lower boundary condition scenario
INTEGER ::            i_ukca_scenario = imdi           ! 0=UM; 1=WMOA1; 2=RCP
CHARACTER(LEN=filenamelength) :: ukca_RCPdir     = 'ukca_RCPdir is unset'  
                                                       ! file location
CHARACTER(LEN=filenamelength) :: ukca_RCPfile    = 'ukca_RCPfile is unset' 
                                                       ! file name

! UKCA LBC inputs =1 if tr in lbc file
INTEGER ::  tc_lbc_ukca(a_max_ukcavars) = 0

! Options to use with ENDGAME
INTEGER :: i_ukca_conserve_method = 0
                   ! Use separate conservation method for UKCA ?
                   ! 0 = Default tracer conservation
                   ! 1 = Priestley (old) kept for continuity
                   ! 2 = optimised Priestley (recommended default)
                   ! 3 = Conservation Off - only applied to UKCA
                   !      -- was 'conserve_ukca_tracers?'

INTEGER :: i_ukca_hiorder_scheme = imdi
                   ! Use different scheme for High order interpolation?
                   ! Use the same codes as Moisture/tracers
                   ! only active for conserve_meth 1 & 2

LOGICAL :: L_ukca_src_in_conservation = .TRUE.
                   ! physics2 sources in conservation ?
                   ! only active for conserve_meth 1 & 2

! Option values for i_ukca_conserve_method
INTEGER, PARAMETER :: ukca_conserve_um  = 0   ! Original/ ADAS scheme
INTEGER, PARAMETER :: priestley_old     = 1
INTEGER, PARAMETER :: priestley_optimal = 2   ! recommended
INTEGER, PARAMETER :: ukca_no_conserve  = 3   ! Do not conserve tracers

LOGICAL :: l_ukca_ageair =   .FALSE.   ! Allows user to include the
                                       ! Age-of-air tracer on its own or
                                       ! with any chemistry scheme

! Allows the use of heterogeneous chemistry on aerosol surfaces
! from CLASSIC within UKCA.
LOGICAL :: l_ukca_classic_hetchem = .FALSE.

! change resistance based dry deposition scheme to apply deposition 
! losses only in level 1
LOGICAL :: l_ukca_ddep_lev1 = .FALSE.

! RADAER lookup tables and optical properties namelists.
CHARACTER(LEN=filenamelength) :: ukcaaclw = 'ukcaaclw is unset' 
                               !  Aitken + Insol acc mode (LW)
CHARACTER(LEN=filenamelength) :: ukcaacsw = 'ukcaacsw is unset' 
                               !  Aitken + Insol acc mode (SW)
CHARACTER(LEN=filenamelength) :: ukcaanlw = 'ukcaanlw is unset' 
                               !  Soluble accum mode (LW)
CHARACTER(LEN=filenamelength) :: ukcaansw = 'ukcaansw is unset' 
                               !  Soluble accum mode (SW)
CHARACTER(LEN=filenamelength) :: ukcacrlw = 'ukcacrlw is unset' 
                               !  Coarse mode (LW)
CHARACTER(LEN=filenamelength) :: ukcacrsw = 'ukcacrsw is unset' 
                               !  Coarse mode (SW)
CHARACTER(LEN=filenamelength) :: ukcaprec = 'ukcaprec is unset' 
                               !  Precomputed values

CHARACTER(LEN=filenamelength) :: ukcasto3 = 'ukcasto3 is unset'
                               ! UKCA Standard temp and O3 file
CHARACTER(LEN=filenamelength) :: ukcastrd = 'ukcastrd is unset'
                               ! UKCA Photolysis table

REAL  :: lightnox_scale_fac = rmdi   ! Lightning NOX ems scale factor

LOGICAL :: l_ukca_so2ems_expvolc = .FALSE.  ! If True, SO2 emissions from 
                ! specific explosive volcanic eruptions are included for
                ! StratTrop + GLOMAP configuration
INTEGER :: i_ukca_solcyc = 0  ! Use solar cycle in photolysis
REAL :: seadms_ems_scaling = rmdi    ! Marine DMS emission scaling factor

! options for quasi-Newton (Broyden) Method to reduce number of iterations
! in asad_spimpmjp
LOGICAL :: l_ukca_quasinewton       = .FALSE. 
         ! F=do not perform, T=perform
INTEGER :: i_ukca_quasinewton_start = imdi 
         ! iter to start quasi-Newton step (>=2,<=50 2 recommended)
INTEGER :: i_ukca_quasinewton_end   = imdi 
         ! iter to stop quasi-Newton step (>=2,<=50 3 recommended) 

INTEGER :: i_ukca_sad_months = imdi  ! Used for length of SAD ancil file
INTEGER :: i_ukca_sad_start_year = imdi ! Used to set start year of SAD file 

! Parameters to control how the near-surface values in Age-of-air tracer
! are reset to zero - based on model-level-number (1), 
!   or height-above-ground (2)
INTEGER :: i_ageair_reset_method = imdi 
INTEGER :: max_ageair_reset_level = imdi  ! Max level to which to reset
REAL    :: max_ageair_reset_height = rmdi ! Max height (m) to which to reset

INTEGER, PARAMETER :: reset_age_by_level = 1     
INTEGER, PARAMETER :: reset_age_by_height = 2
     
! Define the RUN_UKCA namelist

NAMELIST/run_ukca/ l_ukca, l_ukca_aie1, l_ukca_aie2,              &
         i_ukca_chem, l_ukca_chem_aero, l_ukca_ageair,            &
         i_ukca_photol,                                           &
         l_ukca_mode,                                             &
         l_ukca_dust,                                             &
         l_ukca_qch4inter,                                        &
         l_ukca_het_psc, l_ukca_sa_clim,                          &
         l_ukca_h2o_feedback,                                     &
         l_ukca_rado3, l_ukca_radch4, l_ukca_radn2o,              &
         l_ukca_radf11, l_ukca_radf12, l_ukca_radf113,            &
         l_ukca_radf22, l_ukca_radaer, l_ukca_radaer_sustrat,     &
         l_ukca_intdd, l_ukca_trophet, l_ukca_prescribech4,       &
         l_ukca_set_trace_gases, l_ukca_use_background_aerosol,   &
         l_ukca_asad_columns,                                     &
         l_ukca_primsu, l_ukca_primss,                            &
         l_ukca_primbcoc, l_ukca_prim_moc,                        &
         l_ukca_primdu, l_ukca_use_2dtop,                         &
         l_bcoc_ff, l_bcoc_bf, l_bcoc_bm, l_mode_bhn_on,          &
         l_mode_bln_on, l_ukca_arg_act,                           &
         l_ukca_sfix, i_mode_setup, i_mode_nzts,                  &
         i_mode_bln_param_method, mode_parfrac,                   &
         mode_aitsol_cvscav, mode_activation_dryr,                & 
         mode_incld_so2_rfrac,                                    &
         l_ukca_scale_biom_aer_ems, biom_aer_ems_scaling,         &
         l_ukca_scale_soa_yield, soa_yield_scaling,               &
         chem_timestep, dts0, nit,                                &
         jvspec_dir, jvspec_file, jvscat_file, jvsolar_file,      &
         phot2d_dir,strat2d_dir, fastjx_numwl, fastjx_mode,       &
         fastjx_prescutoff, dir_strat_aer, file_strat_aer,        &
         i_ukca_scenario, ukca_RCPdir, ukca_RCPfile,              &
         ukca_MeBrmmr, ukca_MeClmmr, ukca_CH2Br2mmr, ukca_H2mmr,  &
         ukca_N2mmr, ukca_CFC115mmr, ukca_CCl4mmr,                &
         ukca_MeCCl3mmr, ukca_HCFC141bmmr, ukca_HCFC142bmmr,      &
         ukca_H1211mmr, ukca_H1202mmr, ukca_H1301mmr,             &
         ukca_H2402mmr, ukca_COSmmr,                              &
         ukca_em_dir, ukca_em_files, tc_lbc_ukca,                 &
         i_ukca_conserve_method, i_ukca_hiorder_scheme,           &
         L_ukca_src_in_conservation,                              &
         i_ukca_dms_flux,                                         &
         ukca_offline_dir, ukca_offline_files, l_ukca_ibvoc,      &
         l_ukca_chem_plev, l_ukca_asad_plev,                      & 
         l_ukca_classic_hetchem, l_ukca_ddep_lev1,                &
         ukcaaclw, ukcaacsw, ukcaanlw, ukcaansw, ukcacrlw,        &
         ukcacrsw, ukcaprec, lightnox_scale_fac,                  &
         L_ukca_so2ems_expvolc, i_ukca_solcyc,                    &
         l_ukca_scale_seadms_ems, seadms_ems_scaling,             &
         l_ukca_quasinewton, i_ukca_quasinewton_start,            &
         i_ukca_quasinewton_end,                                  &
         i_ageair_reset_method, max_ageair_reset_level,           &
         max_ageair_reset_height,                                 &
         i_ukca_sad_months, i_ukca_sad_start_year,                &
         l_ukca_limit_nat, l_ukca_linox_scaling,                  &
         l_ukca_debug_asad

! -----------------------------------------------------------------------------
! These are set in ukca_setup_chem_mod after the namelist is read

LOGICAL :: l_ukca_chem      =.FALSE. ! True when UKCA chemistry is on
LOGICAL :: l_ukca_trop      =.FALSE. ! True for tropospheric chemistry (B-E)
LOGICAL :: l_ukca_raq       =.FALSE. ! True for regional air quality chem (B-E)
LOGICAL :: l_ukca_raqaero   =.FALSE. ! True for regional air quality chem (B-E)
                                     !         with aerosols
LOGICAL :: l_ukca_offline_be=.FALSE. ! True for offline oxidants chem. (B-E)
LOGICAL :: l_ukca_tropisop  =.FALSE. ! True for trop chemistry + isoprene
LOGICAL :: l_ukca_strat     =.FALSE. ! True for strat+reduced trop chemistry
LOGICAL :: l_ukca_strattrop =.FALSE. ! True for std strat+trop chemistry
LOGICAL :: l_ukca_offline   =.FALSE. ! True for offline oxidants chemistry N-R
LOGICAL :: l_ukca_achem     =.FALSE. ! add aerosol chemistry to scheme (NR)
LOGICAL :: l_ukca_aerchem   =.FALSE. ! True for trop+aerosol chemistry (B-E)
LOGICAL :: l_ukca_advh2o    =.FALSE. ! True when H2O treated as tracer by ASAD
LOGICAL :: l_ukca_nr_aqchem =.FALSE. ! True when aqueous chem required for N-R

! These schemes are not yet included but logicals used in code so needed here
LOGICAL :: l_ukca_stratcfc  =.FALSE. ! True for extended strat chemistry

! File  and directory for reference sulphur aerosol file
! Not currently used as l_use_stratclim in ukca_fastjx is FALSE
CHARACTER(LEN=filenamelength) :: dir_reff_sulp  = 'dir_reff_sulp is unset'
CHARACTER(LEN=filenamelength) :: file_reff_sulp = 'file_reff_sulp is unset'

! Tracers and chemistry integers:

INTEGER :: ukca_int_method  = 0      ! Defines chemical integration method
INTEGER :: jpctr  = 0                ! No. of chemical tracers
INTEGER :: jpspec = 0                ! No. of chemical species
INTEGER :: jpbk   = 0                ! No. of bimolecular reactions
INTEGER :: jptk   = 0                ! No. of termolecular reactions
INTEGER :: jppj   = 0                ! No. of photolytic reactions
INTEGER :: jphk   = 0                ! No. of heterogeneous reactions
INTEGER :: jpnr   = 0                ! jpbk + jptk + jppj + jphk
INTEGER :: jpdd   = 0                ! No. of dry deposited species
INTEGER :: jpdw   = 0                ! No. of wet deposited species

! Logical array controlling tracers - set up in primary based
! on calls to tstmsk
LOGICAL :: tr_ukca_a (0:a_max_ukcavars) = .FALSE.

! Controls whether UKCA tracers are conserved with same option as 
! CLASSIC tracers. Avoids duplication when same (Priestley/ ADAS) scheme 
! and hi-order interpolation method are being used.
! This affects whether UKCA tracers are lumped with CLASSIC ones
! in the original or Priestley scheme.
LOGICAL :: l_conserve_ukca_with_tr

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_OPTION_MOD'

CONTAINS

SUBROUTINE check_run_ukca()

! Description:
!   Subroutine to apply logic checks based on the
!   options selected in the run_ukca namelist.
!   Note that some chemistry scheme specific checks
!   are done in ukca_setup_chem

USE ukca_photo_scheme_mod, ONLY: i_ukca_nophot, i_ukca_phot2d, i_ukca_fastjx
USE ukca_chem_schemes_mod, ONLY: i_ukca_chem_off, i_ukca_chem_offline, &
                                 i_ukca_chem_offline_be,               &
                                 i_ukca_chem_tropisop,                 &
                                 i_ukca_chem_strattrop,                &
                                 i_ukca_chem_strat,                    &
                                 i_ukca_chem_offline
USE gen_phys_inputs_mod,   ONLY: l_mr_physics, l_use_methox
USE umPrintMgr,            ONLY: PrStatus_Normal,PrintStatus, &
                                 umPrint, prnt_writers, outputAll

USE sl_input_mod,          ONLY: l_priestley_correct_tracers,          &
                                 tr_priestley_opt, tracer_sl,          &
                                 high_order_scheme, l_conserve_tracers
USE ereport_mod,           ONLY: ereport
USE mphys_inputs_mod,      ONLY: l_mcr_arcl, l_autoconv_murk

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER   :: RoutineName = 'CHECK_RUN_UKCA'
CHARACTER (LEN=errormessagelength)            :: cmessage   ! Error message
INTEGER                        :: errcode    ! Variable passed to ereport

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

errcode = 0                   ! Initialise 
! Check photolysis switches
IF ( i_ukca_photol /= i_ukca_nophot .AND. &
     i_ukca_photol /= i_ukca_phot2d .AND. &
     i_ukca_photol /= i_ukca_fastjx ) THEN
  cmessage='Unknown photolysis scheme.'
  errcode=ABS(i_ukca_photol)
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! GLOMAP-mode related switches
IF ( l_ukca_mode .AND. .NOT. l_ukca ) THEN
  cmessage='Cannot use GLOMAP-mode aerosols without UKCA'
  errcode=1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF ( l_ukca_mode ) THEN
 ! Set Plume scavenging of Aerosol tracers ON by default
  l_ukca_plume_scav = .TRUE.
  IF ( PrintStatus > PrStatus_Normal )                             &
    CALL umPrint('Plume Scavenging used for GLOMAP-mode aerosols', &
               src='ukca_option_mod')
END IF

IF ( l_ukca_dust .AND. .NOT. l_ukca_mode ) THEN
  cmessage='Cannot use dust without GLOMAP-mode aerosols'
  errcode=5
  CALL ereport(RoutineName,errcode,cmessage)
END IF  

! Direct aerosol effects
IF ( l_ukca_radaer .AND. .NOT. l_ukca_mode ) THEN
  cmessage='Cannot use RADAER without GLOMAP-mode aerosols'
  errcode=2
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Indirect aerosol effects
IF ( (l_ukca_aie1 .OR. l_ukca_aie2) .AND. .NOT. l_ukca_mode ) THEN
  cmessage='Cannot use AIE without GLOMAP-mode aerosols'
  errcode=2
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Aerosol indirect effect
IF (l_ukca_aie2 .AND. l_autoconv_murk) THEN
  cmessage='Cannot set both l_ukca_aie2 and l_autoconv_murk to .true.'
  errcode=34
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Aerosol indirect effect
IF (l_ukca_aie2 .AND. l_mcr_arcl) THEN
  cmessage='Cannot set both l_ukca_aie2 and l_mcr_arcl to .true.'
  errcode=34
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Mixing ratio physics is not suitable for UKCA
IF (l_mr_physics) THEN
  ! q is mixing ratio, not yet appropriate for UKCA
  cmessage = ' UKCA cannot be run with H2O as a mixing ratio'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Methane oxidation - l_use_methox should not be on
! if l_ukca_h2o_feedback also on

IF (l_ukca_h2o_feedback .AND. l_use_methox) THEN
  cmessage='Cannot use parameterised CH4 oxidation and ' // &
    'water vapour feedback from UKCA'
  errcode=2
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Check UKCA conservation options if ENDGame

! The original or ADAS scheme cannot be applied currently if the 
! Priestley scheme has been selected for other tracers.
IF ( l_priestley_correct_tracers .AND.                           &
   (i_ukca_conserve_method == ukca_conserve_um) ) THEN
  errcode = 3
  cmessage = 'Mismatch in Tracer conservation options. ' // &
    ' Select Priestley scheme for UKCA - conserve_method = 1,2'
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Determine whether UKCA tracers can be corrected at the same time  
! as other (CLASSIC, murk, idealised) tracers.  
! This can happen under three conditions -  
!   i. If conservation is off
!  ii. If default (ADAS) scheme is selected for UKCA and other tracers 
! iii. If hi-order scheme and Priestley options are same for both 
!
! The logical actually controls if UKCA is grouped with the other
! tracers in a super-array, during correction/ conservation.
 
IF (.NOT. l_conserve_tracers) THEN
  l_conserve_ukca_with_tr = .TRUE.
ELSE IF ( L_priestley_correct_tracers ) THEN 
  ! check if related options requested are identical 
  l_conserve_ukca_with_tr =                                   & 
     ( i_ukca_conserve_method == tr_priestley_opt )

  IF ( i_ukca_hiorder_scheme > 0 .AND.                        & 
       i_ukca_hiorder_scheme /= high_order_scheme(tracer_sl) )&
   l_conserve_ukca_with_tr = .FALSE.
 
ELSE 
  ! check if a non-default (Priestley) scheme is selected for UKCA only 
  l_conserve_ukca_with_tr =                                    & 
  ( i_ukca_conserve_method == ukca_conserve_um ) 
END IF
WRITE(cmessage,'(A,L1)')'Conserve UKCA with other tracers? ',&
     l_conserve_ukca_with_tr
CALL umPrint(cmessage,src='ukca_option_mod') 


! biom_aer_ems_scaling - valid range 0. - 30.
IF ( (biom_aer_ems_scaling < 0.0 .OR. biom_aer_ems_scaling > 30.0) & 
   .AND. l_ukca_scale_biom_aer_ems ) THEN
  cmessage='biom_aer_ems_scaling should be between 0 - 30'
  errcode = 3
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! seadms_ems_scaling - valid range 0. - 10.
IF ( (seadms_ems_scaling < 0.0 .OR. seadms_ems_scaling > 10.0) & 
   .AND. l_ukca_scale_seadms_ems ) THEN
  cmessage='seadms_ems_scaling should be between 0 - 10'
  errcode = 3
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! soa_yield_scaling - valid range 0. - 5.
IF ( (soa_yield_scaling < 0.0 .OR. soa_yield_scaling > 5.0) & 
   .AND. l_ukca_scale_soa_yield ) THEN
  cmessage='soa_yield_scaling should be between 0 - 5'
  errcode = 3
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF ( (mode_activation_dryr < 20.0 .OR. mode_activation_dryr > 100.0) .AND.   &
     l_ukca_mode) THEN
  cmessage=' mode_activation_dryr should be between 20.0 and 100.0'
  errcode = 4
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF ( (mode_incld_so2_rfrac < 0.0 .OR. mode_incld_so2_rfrac > 1.0) .AND.   &
     l_ukca_mode) THEN
  cmessage=' mode_incld_so2_rfrac should be between 0.0 and 1.0'
  errcode = 5
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Check value of Lightning NOX scaling factor.
! This is not relevant to the Age-of-air only or Offline schemes
IF ( i_ukca_chem /= i_ukca_chem_off      .AND.          &
     i_ukca_chem /= i_ukca_chem_offline  .AND.          &
     i_ukca_chem /= i_ukca_chem_offline_be ) THEN
  IF ( lightnox_scale_fac < 0.0 .OR. lightnox_scale_fac > 5.0 )   THEN
    cmessage=' Light-NOx scale factor should be between 0.0 and 5.0'
    errcode = 6
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF

! Check the values of parameters that control resetting of age-of-air tracer
IF ( l_ukca_ageair ) THEN
  SELECT CASE(i_ageair_reset_method)
  CASE (reset_age_by_level)
    IF ( max_ageair_reset_level < 1 ) THEN
       errcode = 7
       WRITE(cmessage,'(A,I0)')'Inconsistent Age-of-air reset level: ',   &
           max_ageair_reset_level
    END IF
  CASE (reset_age_by_height)
    IF  ( max_ageair_reset_height < 0.0 .OR.   &
          max_ageair_reset_height > 20000.0 ) THEN
       errcode = 8
       WRITE(cmessage,'(A,F16.4)')'Inconsistent Age-of-air reset height: ',&
           max_ageair_reset_height
    END IF
  CASE DEFAULT
    errcode = 9
    WRITE(cmessage,'(A,I0)') ' Inconsistent Age-of-air reset method: ',    &
           i_ageair_reset_method
  END SELECT
  IF ( errcode > 0 ) CALL ereport(RoutineName,errcode,cmessage)

END IF     ! l_ukca_ageair

! Check settings for quasi-Newton method 
!  - only do for Newton-Raphson schemes
IF ( i_ukca_chem == i_ukca_chem_strat      .OR.                       &
     i_ukca_chem == i_ukca_chem_strattrop  .OR.                       &
     i_ukca_chem == i_ukca_chem_tropisop   .OR.                       &
     i_ukca_chem == i_ukca_chem_offline ) THEN
! check settings in quasi-Newton step to ensure they are sensible
   IF (l_ukca_quasinewton) THEN
      IF (i_ukca_quasinewton_end < i_ukca_quasinewton_start) THEN
         WRITE(cmessage,'(A,A,I4,I4)')                                &
              ' i_ukca_quasinewton_start must be less than or equal ',&
              ' to i_ukca_quasinewton_end', i_ukca_quasinewton_start, &
              i_ukca_quasinewton_end
         errcode = 10         
         CALL ereport(RoutineName,errcode,cmessage)
      END IF
      IF ((i_ukca_quasinewton_start < 2) .OR.                         &
           (i_ukca_quasinewton_start > 50)) THEN
         WRITE(cmessage,'(A,I4)')                                     &
              ' i_ukca_quasinewton_start must be between 2 & 50',     &
              i_ukca_quasinewton_start
         errcode = 11         
         CALL ereport(RoutineName,errcode,cmessage)
      END IF
      IF ((i_ukca_quasinewton_end < 2) .OR.                           &
           (i_ukca_quasinewton_end > 50)) THEN
         WRITE(cmessage,'(A,I4)')                                     &
              ' i_ukca_quasinewton_end must be between 2 & 50',       &
              i_ukca_quasinewton_end
         errcode = 12         
         CALL ereport(RoutineName,errcode,cmessage)
      END IF
   END IF
END IF

IF ( .NOT. (i_ukca_chem == i_ukca_chem_strat      .OR.                &
      i_ukca_chem == i_ukca_chem_strattrop  .OR.                      &
      i_ukca_chem == i_ukca_chem_tropisop   .OR.                      &
      i_ukca_chem == i_ukca_chem_offline )) THEN
! check settings in quasi-Newton step to ensure they are sensible
   IF (l_ukca_asad_columns) THEN
    errcode = 10
    WRITE(cmessage,'(A,L1)')                                          &
         ' Column-call can only be for Newton-Raphson schemes: ',     &
         l_ukca_asad_columns
   END IF
END IF

IF ( .NOT. (i_ukca_chem == i_ukca_chem_strat      .OR.                &
      i_ukca_chem == i_ukca_chem_strattrop  .OR.                      &
      i_ukca_chem == i_ukca_chem_tropisop   .OR.                      &
      i_ukca_chem == i_ukca_chem_offline )) THEN
   IF (l_ukca_debug_asad) THEN
    errcode = 10
    WRITE(cmessage,'(A,L1)')                                          &
         ' ASAD debugging can only be for Newton-Raphson schemes: ',     &
         l_ukca_debug_asad
   END IF
END IF

! Check settings for l_ukca_debug_asad
! - only do for Newton-Raphson schemes
IF ( i_ukca_chem == i_ukca_chem_strat      .OR.                       &
     i_ukca_chem == i_ukca_chem_strattrop  .OR.                       &
     i_ukca_chem == i_ukca_chem_tropisop   .OR.                       &
     i_ukca_chem == i_ukca_chem_offline ) THEN
! Check to ensure that the ASAD solver iteration counter
! only outputs when writing from all MPI tasks
   IF (l_ukca_debug_asad) THEN
      IF (prnt_writers /= outputAll) THEN
         WRITE(cmessage,'(A,L1)')                                           &
              ' prnt_writers must be set to 1 when l_ukca_debug_asad = ', &
              l_ukca_debug_asad
         errcode = 13
         CALL ereport(routinename,errcode,cmessage)
      END IF
   END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE check_run_ukca

SUBROUTINE print_nlist_run_ukca()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
INTEGER :: i ! loop counter
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_UKCA'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_ukca', &
    src='ukca_option_mod')

WRITE(lineBuffer,'(A33,L1)')' l_ukca = ',l_ukca
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_aie1 = ',l_ukca_aie1
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_aie2 = ',l_ukca_aie2
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I6)')' i_ukca_chem = ',i_ukca_chem
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_chem_aero = ',l_ukca_chem_aero
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I6)')' i_ukca_photol = ',i_ukca_photol
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_mode = ',l_ukca_mode
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_dust = ',l_ukca_dust
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_qch4inter = ',l_ukca_qch4inter
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_het_psc = ',l_ukca_het_psc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_limit_nat = ',l_ukca_limit_nat
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_sa_clim = ',l_ukca_sa_clim
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_h2o_feedback = ',l_ukca_h2o_feedback
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_rado3 = ',l_ukca_rado3
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_radch4 = ',l_ukca_radch4
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_radn2o = ',l_ukca_radn2o
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_radf11 = ',l_ukca_radf11
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_radf12 = ',l_ukca_radf12
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_radf113 = ',l_ukca_radf113
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_radf22 = ',l_ukca_radf22
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_radaer = ',l_ukca_radaer
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_radaer_sustrat = ',l_ukca_radaer_sustrat
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_intdd = ',l_ukca_intdd
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_trophet = ',l_ukca_trophet
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_prescribech4 = ',l_ukca_prescribech4
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_set_trace_gases = ',        &
      l_ukca_set_trace_gases
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_use_background_aerosol = ', &
      l_ukca_use_background_aerosol
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_asad_columns = ', &
      l_ukca_asad_columns
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')' l_ukca_linox_scaling = ', &
      l_ukca_linox_scaling
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_debug_asad = ', &
      l_ukca_debug_asad
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_primsu = ',l_ukca_primsu
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_primss = ',l_ukca_primss
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_primbcoc = ',l_ukca_primbcoc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_prim_moc = ',l_ukca_prim_moc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_primdu = ',l_ukca_primdu
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_use_2dtop = ',l_ukca_use_2dtop
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_bcoc_ff = ',l_bcoc_ff
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_bcoc_bf = ',l_bcoc_bf
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_bcoc_bm = ',l_bcoc_bm
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_mode_bhn_on = ',l_mode_bhn_on
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_mode_bln_on = ',l_mode_bln_on
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_arg_act = ',l_ukca_arg_act
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_sfix = ',l_ukca_sfix
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I6)')' i_mode_setup = ',i_mode_setup
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I6)')' i_mode_nzts = ',i_mode_nzts
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I6)')' i_mode_bln_param_method = ',i_mode_bln_param_method
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' mode_parfrac = ',mode_parfrac
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,L1)')' l_ukca_scale_biom_aer_ems = ', &
     l_ukca_scale_biom_aer_ems
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' biom_aer_ems_scaling = ', &
     biom_aer_ems_scaling
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,L1)')' l_ukca_scale_soa_yield = ', &
     l_ukca_scale_soa_yield
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' soa_yield_scaling = ', &
     soa_yield_scaling
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' mode_incld_so2_rfrac = ', &
     mode_incld_so2_rfrac
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' mode_activation_dryr = ', &
     mode_activation_dryr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I6)')' chem_timestep = ',chem_timestep
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I6)')' dts0 = ',dts0
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I6)')' nit = ',nit
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A22,A)')' jvspec_dir = ',jvspec_dir
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A22,A)')' jvspec_file = ',jvspec_file
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A22,A)')' jvscat_file = ',jvscat_file
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A22,A)')' jvsolar_file = ',jvsolar_file
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A22,A)')' phot2d_dir = ',phot2d_dir
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A22,A)')' strat2d_dir = ',strat2d_dir
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I0)')' fastjx_numwl = ',fastjx_numwl
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I0)')' fastjx_mode = ',fastjx_mode
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' fastjx_prescutoff = ',fastjx_prescutoff
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A22,A)')' dir_strat_aer = ',dir_strat_aer
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A22,A)')' file_strat_aer = ',file_strat_aer
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_MeBrmmr = ',ukca_MeBrmmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_MeClmmr = ',ukca_MeClmmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_CH2Br2mmr = ',ukca_CH2Br2mmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_H2mmr = ',ukca_H2mmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_N2mmr = ',ukca_N2mmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_CFC115mmr = ',ukca_CFC115mmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_CCl4mmr = ',ukca_CCl4mmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_MeCCl3mmr = ',ukca_MeCCl3mmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_HCFC141bmmr = ',ukca_HCFC141bmmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_HCFC142bmmr = ',ukca_HCFC142bmmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_H1211mmr = ',ukca_H1211mmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_H1202mmr = ',ukca_H1202mmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_H1301mmr = ',ukca_H1301mmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_H2402mmr = ',ukca_H2402mmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,E15.6)')' ukca_COSmmr = ',ukca_COSmmr
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I6)')' i_ukca_scenario = ',i_ukca_scenario
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A22,A)')' ukca_RCPdir = ',ukca_RCPdir
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A22,A)')' ukca_RCPfile = ',ukca_RCPfile
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I6)')' i_ukca_conserve_method = ',i_ukca_conserve_method
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I6)')' i_ukca_hiorder_scheme = ',i_ukca_hiorder_scheme
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A33,L1)')' l_ukca_src_in_conservation = ',       &
                               l_ukca_src_in_conservation
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A22,A)')' ukca_em_dir = ',ukca_em_dir
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I3,A,A)')' ukca_em_files(',1,') = ',ukca_em_files(1)
CALL umPrint(lineBuffer,src='ukca_option_mod')
DO i=2, nr_cdf_files
  IF (ukca_em_files(i) /= 'ukca_em_files is unset') THEN
    WRITE(lineBuffer,'(A,I3,A,A)')' ukca_em_files(',i,') = ',ukca_em_files(i)
    CALL umPrint(lineBuffer,src='ukca_option_mod')
  END IF
END DO
WRITE(lineBuffer,'(A25,A)')' ukca_offline_dir = ',ukca_offline_dir
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I1,A,A)')' ukca_offline_files(',1,') = ',  &
                                ukca_offline_files(1)
CALL umPrint(lineBuffer,src='ukca_option_mod')
DO i=2, max_offline_files
  IF (ukca_offline_files(i) /= 'ukca_offline_files is unset') THEN
    WRITE(lineBuffer,'(A,I1,A,A)')' ukca_offline_files(',i,') = ',  &
                                    ukca_offline_files(i)
    CALL umPrint(lineBuffer,src='ukca_option_mod')
  END IF
END DO
WRITE(lineBuffer,'(A)') 'tc_lbc_ukca'
CALL umPrint(lineBuffer,src='ukca_option_mod')
DO i=1,a_max_ukcavars
  WRITE(lineBuffer,'(I4,1X,I1)') i, tc_lbc_ukca(i)
  CALL umPrint(lineBuffer,src='ukca_option_mod')
END DO
WRITE(lineBuffer,'(A,I6)')' i_ukca_dms_flux = ',i_ukca_dms_flux
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A16,L1)')' l_ukca_ibvoc = ',l_ukca_ibvoc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')' L_ukca_ageair = ',L_ukca_ageair
CALL umPrint(lineBuffer,src='ukca_option_mod')

WRITE(lineBuffer,'(A20,L1)')' l_ukca_chem_plev = ',l_ukca_chem_plev
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A20,L1)')' l_ukca_asad_plev = ',l_ukca_asad_plev
CALL umPrint(lineBuffer,src='ukca_option_mod')

WRITE(lineBuffer,'(A25,L1)') 'l_ukca_classic_hetchem = ', & 
                              l_ukca_classic_hetchem
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A20,L1)') 'l_ukca_ddep_lev1 = ', & 
                              l_ukca_ddep_lev1
CALL umPrint(lineBuffer,src='ukca_option_mod')

WRITE(linebuffer,"(A,A)")' ukcaaclw = ', TRIM(ukcaaclw)
CALL umprint(linebuffer,src='ukca_option_mod')
WRITE(linebuffer,"(A,A)")' ukcaacsw = ', TRIM(ukcaacsw)
CALL umprint(linebuffer,src='ukca_option_mod')
WRITE(linebuffer,"(A,A)")' ukcaanlw = ', TRIM(ukcaanlw)
CALL umprint(linebuffer,src='ukca_option_mod')
WRITE(linebuffer,"(A,A)")' ukcaansw = ', TRIM(ukcaansw)
CALL umprint(linebuffer,src='ukca_option_mod')
WRITE(linebuffer,"(A,A)")' ukcacrlw = ', TRIM(ukcacrlw)
CALL umprint(linebuffer,src='ukca_option_mod')
WRITE(linebuffer,"(A,A)")' ukcacrsw = ', TRIM(ukcacrsw)
CALL umprint(linebuffer,src='ukca_option_mod')
WRITE(linebuffer,"(A,A)")' ukcaprec = ', TRIM(ukcaprec)
CALL umprint(linebuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' lightnox_scale_fac = ',lightnox_scale_fac
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')' L_ukca_so2ems_expvolc= ',L_ukca_so2ems_expvolc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I0)')' i_ukca_solcyc= ',i_ukca_solcyc
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A30,L1)')' l_ukca_scale_seadms_ems = ', &
     l_ukca_scale_seadms_ems
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' seadms_ems_scaling = ',seadms_ems_scaling
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,L1)')' l_ukca_quasinewton = ',l_ukca_quasinewton
CALL umprint(linebuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I0)')' i_ukca_quasinewton_start = ',             &
                                              i_ukca_quasinewton_start
CALL umprint(linebuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I0)')' i_ukca_quasinewton_end = ',               &
                                              i_ukca_quasinewton_end
CALL umprint(linebuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I0)')' i_ageair_reset_method = ',i_ageair_reset_method
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,I0)')' max_ageair_reset_level = ',max_ageair_reset_level
CALL umPrint(lineBuffer,src='ukca_option_mod')
WRITE(lineBuffer,'(A,F16.4)')' max_ageair_reset_height= ', &
                                                    max_ageair_reset_height
CALL umPrint(lineBuffer,src='ukca_option_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='ukca_option_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_ukca

#if !defined(LFRIC)
SUBROUTINE read_nml_run_ukca(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER (LEN=errormessagelength) :: iomessage
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_UKCA'

! set number of each type of variable in my_namelist type
! No of variable types in namelist
INTEGER, PARAMETER :: no_of_types = 4               
! No of integer variables in namelist
INTEGER, PARAMETER :: n_int = 21 + a_max_ukcavars   
! No of real variables in namelist
INTEGER, PARAMETER :: n_real = 25                   
! No of logical variables in namelist
INTEGER, PARAMETER :: n_log = 53
! No of string variables in namelist
INTEGER, PARAMETER :: n_chars = 10 * filenamelength                   &
                        + filenamelength * (1+ nr_cdf_files)          &
                        + filenamelength * (1+ max_offline_files)     &
                        + filenamelength * 7  ! RADAER namelists  

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_ukca_chem
  INTEGER :: i_ukca_photol
  INTEGER :: i_mode_setup
  INTEGER :: i_mode_nzts
  INTEGER :: i_mode_bln_param_method
  INTEGER :: chem_timestep
  INTEGER :: dts0
  INTEGER :: nit
  INTEGER :: fastjx_numwl
  INTEGER :: fastjx_mode
  INTEGER :: i_ukca_scenario
  INTEGER :: tc_lbc_ukca(a_max_ukcavars)
  INTEGER :: i_ukca_conserve_method
  INTEGER :: i_ukca_hiorder_scheme
  INTEGER :: i_ukca_dms_flux
  INTEGER :: i_ukca_quasinewton_start
  INTEGER :: i_ukca_quasinewton_end
  INTEGER :: i_ageair_reset_method
  INTEGER :: max_ageair_reset_level
  INTEGER :: i_ukca_sad_months
  INTEGER :: i_ukca_sad_start_year
  INTEGER :: i_ukca_solcyc
  REAL :: mode_parfrac
  REAL :: mode_aitsol_cvscav
  REAL :: mode_activation_dryr
  REAL :: mode_incld_so2_rfrac
  REAL :: fastjx_prescutoff
  REAL :: ukca_MeBrmmr
  REAL :: ukca_MeClmmr
  REAL :: ukca_CH2Br2mmr
  REAL :: ukca_H2mmr
  REAL :: ukca_N2mmr
  REAL :: ukca_CFC115mmr
  REAL :: ukca_CCl4mmr
  REAL :: ukca_MeCCl3mmr
  REAL :: ukca_HCFC141bmmr
  REAL :: ukca_HCFC142bmmr
  REAL :: ukca_H1211mmr
  REAL :: ukca_H1202mmr
  REAL :: ukca_H1301mmr
  REAL :: ukca_H2402mmr
  REAL :: ukca_COSmmr
  REAL :: biom_aer_ems_scaling
  REAL :: soa_yield_scaling
  REAL :: lightnox_scale_fac
  REAL :: seadms_ems_scaling
  REAL :: max_ageair_reset_height
  LOGICAL :: l_ukca
  LOGICAL :: l_ukca_aie1
  LOGICAL :: l_ukca_aie2
  LOGICAL :: l_ukca_chem_aero
  LOGICAL :: l_ukca_mode
  LOGICAL :: l_ukca_dust
  LOGICAL :: l_ukca_qch4inter
  LOGICAL :: l_ukca_het_psc
  LOGICAL :: l_ukca_sa_clim
  LOGICAL :: l_ukca_h2o_feedback
  LOGICAL :: l_ukca_rado3
  LOGICAL :: l_ukca_radch4
  LOGICAL :: l_ukca_radn2o
  LOGICAL :: l_ukca_radf11
  LOGICAL :: l_ukca_radf12
  LOGICAL :: l_ukca_radf113
  LOGICAL :: l_ukca_radf22
  LOGICAL :: l_ukca_radaer
  LOGICAL :: l_ukca_radaer_sustrat
  LOGICAL :: l_ukca_intdd
  LOGICAL :: l_ukca_trophet
  LOGICAL :: l_ukca_prescribech4
  LOGICAL :: l_ukca_set_trace_gases
  LOGICAL :: l_ukca_use_background_aerosol
  LOGICAL :: l_ukca_asad_columns
  LOGICAL :: l_ukca_primsu
  LOGICAL :: l_ukca_primss
  LOGICAL :: l_ukca_primbcoc
  LOGICAL :: l_ukca_prim_moc
  LOGICAL :: l_ukca_primdu
  LOGICAL :: l_ukca_use_2dtop
  LOGICAL :: l_bcoc_ff
  LOGICAL :: l_bcoc_bf
  LOGICAL :: l_bcoc_bm
  LOGICAL :: l_mode_bhn_on
  LOGICAL :: l_mode_bln_on
  LOGICAL :: l_ukca_arg_act
  LOGICAL :: l_ukca_sfix
  LOGICAL :: L_ukca_src_in_conservation
  LOGICAL :: L_ukca_ibvoc
  LOGICAL :: l_ukca_scale_soa_yield
  LOGICAL :: l_ukca_scale_biom_aer_ems
  LOGICAL :: l_ukca_scale_seadms_ems
  LOGICAL :: l_ukca_ageair
  LOGICAL :: l_ukca_chem_plev
  LOGICAL :: l_ukca_asad_plev
  LOGICAL :: l_ukca_classic_hetchem
  LOGICAL :: l_ukca_ddep_lev1
  LOGICAL :: l_ukca_so2ems_expvolc
  LOGICAL :: l_ukca_quasinewton
  LOGICAL :: l_ukca_limit_nat  
  LOGICAL :: l_ukca_linox_scaling
  LOGICAL :: l_ukca_debug_asad
  CHARACTER (LEN=filenamelength) :: jvspec_dir
  CHARACTER (LEN=filenamelength) :: jvspec_file
  CHARACTER (LEN=filenamelength) :: jvscat_file
  CHARACTER (LEN=filenamelength) :: jvsolar_file
  CHARACTER (LEN=filenamelength) :: phot2d_dir
  CHARACTER (LEN=filenamelength) :: strat2d_dir
  CHARACTER (LEN=filenamelength) :: dir_strat_aer
  CHARACTER (LEN=filenamelength) :: file_strat_aer
  CHARACTER (LEN=filenamelength) :: ukca_RCPdir
  CHARACTER (LEN=filenamelength) :: ukca_RCPfile
  CHARACTER (LEN=filenamelength) :: ukca_em_dir
  CHARACTER (LEN=filenamelength) :: ukca_em_files(nr_cdf_files)
  CHARACTER (LEN=filenamelength) :: ukca_offline_dir
  CHARACTER (LEN=filenamelength) :: ukca_offline_files(max_offline_files)
  CHARACTER (LEN=filenamelength) :: ukcaaclw 
  CHARACTER (LEN=filenamelength) :: ukcaacsw 
  CHARACTER (LEN=filenamelength) :: ukcaanlw
  CHARACTER (LEN=filenamelength) :: ukcaansw
  CHARACTER (LEN=filenamelength) :: ukcacrlw
  CHARACTER (LEN=filenamelength) :: ukcacrsw
  CHARACTER (LEN=filenamelength) :: ukcaprec
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,        &
                    n_real_in=n_real, n_log_in=n_log, n_chars_in=n_chars)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=run_ukca, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_UKCA", iomessage)

  my_nml % i_ukca_chem     = i_ukca_chem
  my_nml % i_ukca_photol   = i_ukca_photol
  my_nml % i_mode_setup    = i_mode_setup
  my_nml % i_mode_nzts     = i_mode_nzts
  my_nml % i_mode_bln_param_method = i_mode_bln_param_method
  my_nml % chem_timestep   = chem_timestep
  my_nml % dts0            = dts0
  my_nml % nit             = nit
  my_nml % fastjx_numwl    = fastjx_numwl
  my_nml % fastjx_mode     = fastjx_mode
  my_nml % i_ukca_scenario = i_ukca_scenario
  my_nml % tc_lbc_ukca     = tc_lbc_ukca
  my_nml % i_ukca_conserve_method = i_ukca_conserve_method
  my_nml % i_ukca_hiorder_scheme = i_ukca_hiorder_scheme
  my_nml % i_ukca_dms_flux = i_ukca_dms_flux
  my_nml % i_ukca_quasinewton_start = i_ukca_quasinewton_start
  my_nml % i_ukca_quasinewton_end = i_ukca_quasinewton_end
  my_nml % i_ageair_reset_method  = i_ageair_reset_method
  my_nml % max_ageair_reset_level  = max_ageair_reset_level
  my_nml % i_ukca_sad_months  = i_ukca_sad_months
  my_nml % i_ukca_sad_start_year  = i_ukca_sad_start_year
  my_nml % i_ukca_solcyc       = i_ukca_solcyc
  ! end of integers
  my_nml % mode_parfrac       = mode_parfrac
  my_nml % mode_aitsol_cvscav = mode_aitsol_cvscav
  my_nml % mode_activation_dryr = mode_activation_dryr
  my_nml % mode_incld_so2_rfrac = mode_incld_so2_rfrac
  my_nml % fastjx_prescutoff  = fastjx_prescutoff
  my_nml % ukca_MeBrmmr       = ukca_MeBrmmr
  my_nml % ukca_MeClmmr       = ukca_MeClmmr
  my_nml % ukca_CH2Br2mmr     = ukca_CH2Br2mmr
  my_nml % ukca_H2mmr         = ukca_H2mmr
  my_nml % ukca_N2mmr         = ukca_N2mmr
  my_nml % ukca_CFC115mmr     = ukca_CFC115mmr
  my_nml % ukca_CCl4mmr       = ukca_CCl4mmr
  my_nml % ukca_MeCCl3mmr     = ukca_MeCCl3mmr
  my_nml % ukca_HCFC141bmmr   = ukca_HCFC141bmmr
  my_nml % ukca_HCFC142bmmr   = ukca_HCFC142bmmr
  my_nml % ukca_H1211mmr      = ukca_H1211mmr
  my_nml % ukca_H1202mmr      = ukca_H1202mmr
  my_nml % ukca_H1301mmr      = ukca_H1301mmr
  my_nml % ukca_H2402mmr      = ukca_H2402mmr
  my_nml % ukca_COSmmr        = ukca_COSmmr
  my_nml % biom_aer_ems_scaling = biom_aer_ems_scaling
  my_nml % soa_yield_scaling  = soa_yield_scaling
  my_nml % lightnox_scale_fac = lightnox_scale_fac
  my_nml % seadms_ems_scaling = seadms_ems_scaling
  my_nml % max_ageair_reset_height  = max_ageair_reset_height
  ! end of reals
  my_nml % l_ukca              = l_ukca
  my_nml % l_ukca_aie1         = l_ukca_aie1
  my_nml % l_ukca_aie2         = l_ukca_aie2
  my_nml % l_ukca_chem_aero    = l_ukca_chem_aero
  my_nml % l_ukca_mode         = l_ukca_mode
  my_nml % l_ukca_dust         = l_ukca_dust
  my_nml % l_ukca_qch4inter    = l_ukca_qch4inter
  my_nml % l_ukca_het_psc      = l_ukca_het_psc
  my_nml % l_ukca_sa_clim      = l_ukca_sa_clim
  my_nml % l_ukca_h2o_feedback = l_ukca_h2o_feedback
  my_nml % l_ukca_rado3        = l_ukca_rado3
  my_nml % l_ukca_radch4       = l_ukca_radch4
  my_nml % l_ukca_radn2o       = l_ukca_radn2o
  my_nml % l_ukca_radf11       = l_ukca_radf11
  my_nml % l_ukca_radf12       = l_ukca_radf12
  my_nml % l_ukca_radf113      = l_ukca_radf113
  my_nml % l_ukca_radf22       = l_ukca_radf22
  my_nml % l_ukca_radaer       = l_ukca_radaer
  my_nml % l_ukca_radaer_sustrat = l_ukca_radaer_sustrat
  my_nml % l_ukca_intdd        = l_ukca_intdd
  my_nml % l_ukca_trophet      = l_ukca_trophet
  my_nml % l_ukca_prescribech4 = l_ukca_prescribech4
  my_nml % l_ukca_set_trace_gases = l_ukca_set_trace_gases
  my_nml % l_ukca_use_background_aerosol = l_ukca_use_background_aerosol
  my_nml % l_ukca_asad_columns = l_ukca_asad_columns
  my_nml % l_ukca_primsu       = l_ukca_primsu
  my_nml % l_ukca_primss       = l_ukca_primss
  my_nml % l_ukca_primbcoc     = l_ukca_primbcoc
  my_nml % l_ukca_prim_moc     = l_ukca_prim_moc
  my_nml % l_ukca_primdu       = l_ukca_primdu
  my_nml % l_ukca_use_2dtop    = l_ukca_use_2dtop
  my_nml % l_bcoc_ff           = l_bcoc_ff
  my_nml % l_bcoc_bf           = l_bcoc_bf
  my_nml % l_bcoc_bm           = l_bcoc_bm
  my_nml % l_mode_bhn_on       = l_mode_bhn_on
  my_nml % l_mode_bln_on       = l_mode_bln_on
  my_nml % l_ukca_arg_act      = l_ukca_arg_act
  my_nml % l_ukca_sfix         = l_ukca_sfix
  my_nml % L_ukca_src_in_conservation = L_ukca_src_in_conservation
  my_nml % l_ukca_ibvoc        = l_ukca_ibvoc
  my_nml % l_ukca_scale_biom_aer_ems = l_ukca_scale_biom_aer_ems
  my_nml % l_ukca_scale_seadms_ems = l_ukca_scale_seadms_ems
  my_nml % l_ukca_scale_soa_yield = l_ukca_scale_soa_yield
  my_nml % l_ukca_ageair       = l_ukca_ageair
  my_nml % l_ukca_chem_plev    = l_ukca_chem_plev
  my_nml % l_ukca_asad_plev    = l_ukca_asad_plev
  my_nml % l_ukca_classic_hetchem = l_ukca_classic_hetchem
  my_nml % l_ukca_ddep_lev1 = l_ukca_ddep_lev1
  my_nml % l_ukca_so2ems_expvolc = l_ukca_so2ems_expvolc
  my_nml % l_ukca_quasinewton = l_ukca_quasinewton
  my_nml % l_ukca_limit_nat = l_ukca_limit_nat
  my_nml % l_ukca_linox_scaling = l_ukca_linox_scaling
  my_nml % l_ukca_debug_asad = l_ukca_debug_asad
  ! end of logicals
  my_nml % jvspec_dir     = jvspec_dir
  my_nml % jvspec_file    = jvspec_file
  my_nml % jvscat_file    = jvscat_file
  my_nml % jvsolar_file   = jvsolar_file
  my_nml % phot2d_dir     = phot2d_dir
  my_nml % strat2d_dir    = strat2d_dir
  my_nml % dir_strat_aer  = dir_strat_aer
  my_nml % file_strat_aer = file_strat_aer
  my_nml % ukca_RCPdir    = ukca_RCPdir
  my_nml % ukca_RCPfile   = ukca_RCPfile
  my_nml % ukca_em_dir    = ukca_em_dir
  my_nml % ukca_em_files  = ukca_em_files
  my_nml % ukca_offline_dir = ukca_offline_dir
  my_nml % ukca_offline_files = ukca_offline_files
  my_nml % ukcaaclw = ukcaaclw 
  my_nml % ukcaacsw = ukcaacsw 
  my_nml % ukcaanlw = ukcaanlw 
  my_nml % ukcaansw = ukcaansw 
  my_nml % ukcacrlw = ukcacrlw 
  my_nml % ukcacrsw = ukcacrsw 
  my_nml % ukcaprec = ukcaprec 
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  i_ukca_chem     = my_nml % i_ukca_chem
  i_ukca_photol   = my_nml % i_ukca_photol
  i_mode_setup    = my_nml % i_mode_setup
  i_mode_nzts     = my_nml % i_mode_nzts
  i_mode_bln_param_method = my_nml % i_mode_bln_param_method
  chem_timestep   = my_nml % chem_timestep
  dts0            = my_nml % dts0
  nit             = my_nml % nit
  fastjx_numwl    = my_nml % fastjx_numwl
  fastjx_mode     = my_nml % fastjx_mode
  i_ukca_scenario = my_nml % i_ukca_scenario
  tc_lbc_ukca     = my_nml % tc_lbc_ukca
  i_ukca_conserve_method = my_nml % i_ukca_conserve_method
  i_ukca_hiorder_scheme = my_nml % i_ukca_hiorder_scheme
  i_ukca_dms_flux = my_nml % i_ukca_dms_flux
  i_ukca_quasinewton_start = my_nml % i_ukca_quasinewton_start
  i_ukca_quasinewton_end = my_nml % i_ukca_quasinewton_end
  i_ageair_reset_method  = my_nml % i_ageair_reset_method
  max_ageair_reset_level  = my_nml % max_ageair_reset_level
  i_ukca_sad_months  = my_nml % i_ukca_sad_months
  i_ukca_sad_start_year  = my_nml % i_ukca_sad_start_year
  i_ukca_solcyc       = my_nml % i_ukca_solcyc
  ! end of integers
  mode_parfrac       = my_nml % mode_parfrac
  mode_aitsol_cvscav = my_nml % mode_aitsol_cvscav
  mode_activation_dryr = my_nml % mode_activation_dryr
  mode_incld_so2_rfrac = my_nml % mode_incld_so2_rfrac
  fastjx_prescutoff  = my_nml % fastjx_prescutoff
  ukca_MeBrmmr       = my_nml % ukca_MeBrmmr
  ukca_MeClmmr       = my_nml % ukca_MeClmmr
  ukca_CH2Br2mmr     = my_nml % ukca_CH2Br2mmr
  ukca_H2mmr         = my_nml % ukca_H2mmr
  ukca_N2mmr         = my_nml % ukca_N2mmr
  ukca_CFC115mmr     = my_nml % ukca_CFC115mmr
  ukca_CCl4mmr       = my_nml % ukca_CCl4mmr
  ukca_MeCCl3mmr     = my_nml % ukca_MeCCl3mmr
  ukca_HCFC141bmmr   = my_nml % ukca_HCFC141bmmr
  ukca_HCFC142bmmr   = my_nml % ukca_HCFC142bmmr
  ukca_H1211mmr      = my_nml % ukca_H1211mmr
  ukca_H1202mmr      = my_nml % ukca_H1202mmr
  ukca_H1301mmr      = my_nml % ukca_H1301mmr
  ukca_H2402mmr      = my_nml % ukca_H2402mmr
  ukca_COSmmr        = my_nml % ukca_COSmmr
  biom_aer_ems_scaling = my_nml % biom_aer_ems_scaling
  soa_yield_scaling   = my_nml % soa_yield_scaling
  lightnox_scale_fac  = my_nml % lightnox_scale_fac
  seadms_ems_scaling  = my_nml % seadms_ems_scaling
  max_ageair_reset_height  = my_nml % max_ageair_reset_height
  ! end of reals
  l_ukca              = my_nml % l_ukca
  l_ukca_aie1         = my_nml % l_ukca_aie1
  l_ukca_aie2         = my_nml % l_ukca_aie2
  l_ukca_chem_aero    = my_nml % l_ukca_chem_aero
  l_ukca_mode         = my_nml % l_ukca_mode
  l_ukca_dust         = my_nml % l_ukca_dust
  l_ukca_qch4inter    = my_nml % l_ukca_qch4inter
  l_ukca_het_psc      = my_nml % l_ukca_het_psc
  l_ukca_sa_clim      = my_nml % l_ukca_sa_clim
  l_ukca_h2o_feedback = my_nml % l_ukca_h2o_feedback
  l_ukca_rado3        = my_nml % l_ukca_rado3
  l_ukca_radch4       = my_nml % l_ukca_radch4
  l_ukca_radn2o       = my_nml % l_ukca_radn2o
  l_ukca_radf11       = my_nml % l_ukca_radf11
  l_ukca_radf12       = my_nml % l_ukca_radf12
  l_ukca_radf113      = my_nml % l_ukca_radf113
  l_ukca_radf22       = my_nml % l_ukca_radf22
  l_ukca_radaer       = my_nml % l_ukca_radaer
  l_ukca_radaer_sustrat = my_nml % l_ukca_radaer_sustrat
  l_ukca_intdd        = my_nml % l_ukca_intdd
  l_ukca_trophet      = my_nml % l_ukca_trophet
  l_ukca_prescribech4 = my_nml % l_ukca_prescribech4
  l_ukca_set_trace_gases = my_nml % l_ukca_set_trace_gases
  l_ukca_use_background_aerosol= my_nml % l_ukca_use_background_aerosol
  l_ukca_asad_columns = my_nml % l_ukca_asad_columns
  l_ukca_primsu       = my_nml % l_ukca_primsu
  l_ukca_primss       = my_nml % l_ukca_primss
  l_ukca_primbcoc     = my_nml % l_ukca_primbcoc
  l_ukca_prim_moc     = my_nml % l_ukca_prim_moc
  l_ukca_primdu       = my_nml % l_ukca_primdu
  l_ukca_use_2dtop    = my_nml % l_ukca_use_2dtop
  l_bcoc_ff           = my_nml % l_bcoc_ff
  l_bcoc_bf           = my_nml % l_bcoc_bf
  l_bcoc_bm           = my_nml % l_bcoc_bm
  l_mode_bhn_on       = my_nml % l_mode_bhn_on
  l_mode_bln_on       = my_nml % l_mode_bln_on
  l_ukca_arg_act      = my_nml % l_ukca_arg_act
  l_ukca_sfix         = my_nml % l_ukca_sfix
  L_ukca_src_in_conservation = my_nml % L_ukca_src_in_conservation
  l_ukca_ibvoc        = my_nml % l_ukca_ibvoc
  l_ukca_scale_biom_aer_ems = my_nml % l_ukca_scale_biom_aer_ems
  l_ukca_scale_seadms_ems = my_nml % l_ukca_scale_seadms_ems
  l_ukca_scale_soa_yield    = my_nml % l_ukca_scale_soa_yield
  l_ukca_ageair       = my_nml % l_ukca_ageair 
  l_ukca_chem_plev    = my_nml % l_ukca_chem_plev
  l_ukca_asad_plev    = my_nml % l_ukca_asad_plev
  l_ukca_classic_hetchem = my_nml % l_ukca_classic_hetchem
  l_ukca_ddep_lev1 = my_nml % l_ukca_ddep_lev1
  l_ukca_so2ems_expvolc  = my_nml % l_ukca_so2ems_expvolc
  l_ukca_quasinewton = my_nml % l_ukca_quasinewton
  l_ukca_limit_nat = my_nml % l_ukca_limit_nat
  l_ukca_linox_scaling = my_nml % l_ukca_linox_scaling
  l_ukca_debug_asad = my_nml % l_ukca_debug_asad
  ! end of logicals
  jvspec_dir     = my_nml % jvspec_dir
  jvspec_file    = my_nml % jvspec_file
  jvscat_file    = my_nml % jvscat_file
  jvsolar_file   = my_nml % jvsolar_file
  phot2d_dir     = my_nml % phot2d_dir
  strat2d_dir    = my_nml % strat2d_dir
  dir_strat_aer  = my_nml % dir_strat_aer
  file_strat_aer = my_nml % file_strat_aer
  ukca_RCPdir    = my_nml % ukca_RCPdir
  ukca_RCPfile   = my_nml % ukca_RCPfile
  ukca_em_dir    = my_nml % ukca_em_dir
  ukca_em_files  = my_nml % ukca_em_files
  ukca_offline_dir   = my_nml % ukca_offline_dir
  ukca_offline_files = my_nml % ukca_offline_files
  ukcaaclw       = my_nml % ukcaaclw 
  ukcaacsw       = my_nml % ukcaacsw 
  ukcaanlw       = my_nml % ukcaanlw 
  ukcaansw       = my_nml % ukcaansw 
  ukcacrlw       = my_nml % ukcacrlw 
  ukcacrsw       = my_nml % ukcacrsw 
  ukcaprec       = my_nml % ukcaprec 
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_ukca
#endif

END MODULE ukca_option_mod

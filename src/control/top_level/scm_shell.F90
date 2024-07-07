#if !defined(LFRIC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! scm_shell is the main calling program for the Single Column Model.
! It sets up the vertical level information read in from the namelist and
! passes the information down to Scm_Main which then performs the run.
!
! Program scm_shell
!=====================================================================
!                     SCM
!           Single Column Unified Model
!                  Master Deck
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
!=====================================================================
PROGRAM scm_shell

USE nlstcall_mod,         ONLY: nlstcall, lcal360
USE planet_constants_mod, ONLY: planet_constants, set_planet_constants, &
                                print_nlist_planet_constants
USE model_domain_mod,     ONLY: model_domain, print_nlist_model_domain, &
                                check_nml_model_domain

! Dependencies for interfacing with C routines
USE io_dependencies

USE dynamics_input_mod
USE dynamics_testing_mod

! SCM Modules
!---------------------------------------------------------------------------
USE scm_utils
USE scm_cntl_mod
USE s_main_force, ONLY: netcdf_file, obs_t0_prf, source
USE s_nc_obs, ONLY: read_nc_obs
USE global_SCMop, ONLY: incdf
USE netcdf
USE init_scm_misc_mod, ONLY: init_scm_misc

! Physics modules
!---------------------------------------------------------------------------
! Convection modules
!---------------------------------------------------------------------------
USE cv_run_mod                         ! Access to all variables
USE cv_param_mod                       ! Access to all variables

! Module for RUN_Precip namelist
USE mphys_inputs_mod, ONLY: RUN_precip, print_nlist_run_precip,   &
                            check_run_precip

USE electric_inputs_mod, ONLY: l_use_electric

! Module for RUN_Cloud namelist
USE cloud_inputs_mod, ONLY: rhcrit, falliceshear_method,          &
 forced_cu,                                                       &
 i_pc2_conv_coupling, i_pc2_erosion_method, i_eacf,               &
 l_ensure_min_in_cloud_qcf,                                       &
 l_fixbug_pc2_mixph, l_micro_eros,                                &
 dbsdtbs_turb_0, i_cld_area, l_add_cca_to_mcica,                  &
 starticeTKelvin, alliceTdegC, cff_spread_rate, ice_width,        &
 i_cld_vn, l_pc2_lbc, i_rhcpt, RUN_cloud, check_run_cloud,        &
 l_od_cld_filter, tau_thresh, l_sharpen_cbh_diags,                &
 l_ceil_cld_filter, print_nlist_run_cloud
USE pc2_constants_mod, ONLY: i_cld_pc2
USE jules_surface_mod, ONLY: l_flake_model, cor_mo_iter

USE jules_vegetation_mod, ONLY: can_rad_mod

! Module for RUN_Murk namelist
USE murk_inputs_mod, ONLY: RUN_murk, print_nlist_run_murk

! Module for RUN_Rivers namelist
USE river_inputs_mod, ONLY: RUN_Rivers, print_nlist_run_rivers

! Module for RUN_Eng_Corr namelist
USE eng_corr_inputs_mod, ONLY: RUN_Eng_Corr, print_nlist_run_eng_corr

! Module for gen_phys_inputs namelist
USE gen_phys_inputs_mod, ONLY: gen_phys_inputs, print_nlist_gen_phys_inputs

! Boundary layer modules
USE bl_option_mod, ONLY: run_bl, tke_levels, shcu_levels, off,        &
      check_run_bl, i_bl_vn, i_bl_vn_1a, nl_bl_levels, non_local_bl,  &
      blending_option, alpha_cd, print_nlist_run_bl

! segments modules
USE tuning_segments_mod, ONLY: tuning_segments, check_tuning_segments, &
                               print_nlist_tuning_segments

! module for UKCA options
USE ukca_option_mod, ONLY: run_ukca,                            &
       l_ukca, l_ukca_aie1, l_ukca_aie2,                        &
       i_ukca_chem,                                             &
       L_ukca_mode, L_ukca_dust, L_ukca_ageair,                 &
       L_ukca_qch4inter,                                        &
       L_ukca_het_psc, L_ukca_sa_clim,                          &
       L_ukca_h2o_feedback,                                     &
       L_ukca_rado3, L_ukca_radch4, L_ukca_radn2o,              &
       L_ukca_radf11, L_ukca_radf12, L_ukca_radf113,            &
       L_ukca_radf22,                                           &
       L_ukca_intdd, L_ukca_trophet, L_ukca_prescribech4,       &
       L_ukca_set_trace_gases, L_ukca_use_background_aerosol,   &
       L_ukca_primsu, L_ukca_primss,                            &
       L_ukca_primbcoc, L_ukca_primdu, L_ukca_use_2dtop,        &
       L_bcoc_ff, L_bcoc_bf, L_bcoc_bm, L_mode_bhn_on,          &
       L_mode_bln_on, L_ukca_arg_act,                           &
       L_ukca_sfix, i_mode_setup, i_mode_nzts,                  &
       i_mode_bln_param_method, mode_parfrac, dts0, nit,        &
       jvspec_dir, jvspec_file, jvscat_file, phot2d_dir,        &
       strat2d_dir, fastjx_numwl, fastjx_mode,                  &
       fastjx_prescutoff, dir_strat_aer, file_strat_aer,        &
       ukca_MeBrmmr, ukca_MeClmmr, ukca_CH2Br2mmr, ukca_H2mmr,  &
       ukca_N2mmr, ukca_CFC115mmr, ukca_CCl4mmr,                &
       ukca_MeCCl3mmr, ukca_HCFC141bmmr, ukca_HCFC142bmmr,      &
       ukca_H1211mmr, ukca_H1202mmr, ukca_H1301mmr,             &
       ukca_H2402mmr, ukca_COSmmr,                              &
       ukca_em_dir, ukca_em_files, print_nlist_run_ukca

USE glomap_clim_option_mod, ONLY: &
    run_glomap_aeroclim,          &
    l_glomap_mode_clim,           &
    l_glomap_clim_arg_act,        &
    l_glomap_clim_aie1,           &
    l_glomap_clim_aie2,           &
    l_glomap_clim_radaer,         &
    l_glomap_clim_radaer_sustrat, &
    i_glomap_clim_setup,          &
    print_nlist_run_glomap_clim

! GWD modules
USE g_wave_input_mod

! JULES
USE switches, ONLY: l_360
USE jules_surface_mod, ONLY: l_aggregate, Limit_ObukhovL,    &
      cor_mo_iter

USE jules_sea_seaice_mod, ONLY: nice, nice_use, buddy_sea

! Modules for reading JULES namelists:
USE read_jules_namelists_mod, ONLY:                                         &
    read_jules_nvegparm,   read_jules_pftparm, read_jules_triffid,          &
    read_jules_elevate,    read_jules_urban_switches,                       &
    read_jules_vegetation, read_jules_urban2t_param,                        &
    read_jules_hydrology,  read_jules_radiation,  read_jules_sea_seaice,    &
    read_jules_snow,       read_jules_surface_types,                        &
    read_jules_soil,       read_jules_soil_biogeochem,                      &
    read_jules_surface
USE surf_couple_allocate_mod, ONLY: surf_couple_allocate
USE jules_surface_types_mod, ONLY: npft, nnvg
USE calc_ntiles_mod,               ONLY: calc_ntiles

! Radiation modules
USE rad_input_mod
USE sw_rad_input_mod
USE lw_rad_input_mod

USE missing_data_mod, ONLY: imdi

USE science_fixes_mod

! UM Modules
!---------------------------------------------------------------------------
USE atm_fields_bounds_mod, ONLY: atm_fields_bounds_init
USE conversions_mod, ONLY: isec_per_day
USE run_aerosol_mod
USE turb_diff_mod, ONLY: run_diffusion, l_subfilter_vert,            &
                         print_nlist_run_diffusion,check_run_diffusion
USE ereport_mod, ONLY: ereport
USE Control_Max_Sizes
USE um_parvars
USE um_parcore, ONLY: mype
USE filenamelength_mod, ONLY: filenamelength
USE UM_Config, ONLY: &
     appInit,         &
     appTerminate,    &
     exe_scm
USE submodel_mod, ONLY: n_internal_model
USE umPrintMgr, ONLY: umPrint, umMessage, newline,                   &
                      PrintStatus, PrStatus_Normal 
USE nlsizes_namelist_mod, ONLY:                                      &
    bl_levels, cloud_levels, model_levels, ntiles,                   &
    ozone_levels, row_length, rows, sm_levels, st_levels, tr_levels, &
    tr_ukca, tr_vars
USE errormessagelength_mod, ONLY: errormessagelength
USE get_env_var_mod, ONLY: get_env_var

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

USE cv_set_dependent_switches_mod, ONLY: cv_set_dependent_switches
IMPLICIT NONE

! Local Variables

INTEGER, PARAMETER :: &
  ntrop       = 20    &! Max number of levels in the troposphere
                       ! STATS forcing
, sec_day     = isec_per_day  & !  =86400 for earth
, nsprog      = 10     ! No. of single level prognostics

INTEGER, PARAMETER :: &
  co2_dim_len = 1     &! Length of a CO2 field row.
, co2_dim_row = 1      ! Number of CO2 field rows.

REAL ::               &
  dummy

INTEGER ::            &
  length              &
, Istatus             &
, Icode               &
, first_blank         &
, level               &
, j,i,k

LOGICAL :: l_ts_log ! Option for timestep information

LOGICAL :: l_print_namelist = .FALSE.  ! print out namelist entries.

CHARACTER (LEN=filenamelength) ::      &
  dname               &! Directory where the sizes file is kept
, filename            &! Sizes filename to read in basic model dimensions
, vert_lev             ! Vertical level file

CHARACTER(LEN=100) ::  dummy_env
CHARACTER(LEN=50000) :: lineBuffer

CHARACTER(LEN=errormessagelength) :: cmessage
                                          ! Error message if ErrorStatus > 0
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SCM_SHELL'

!---------------------------------------------------------------------
! Variables read in which refer to the SCM forcing namelist
!---------------------------------------------------------------------
INTEGER :: land_points      = 0    ! Default number of land points
INTEGER :: nfor             = imdi ! Number terms for observational
                                   ! forcing, requires setting in
                                   ! namelist
INTEGER :: model_levels_nml = imdi ! Number of model levels
                                   ! specified in supplied
                                   ! namelist. Must be set in namelist.

! Variables for NC_obs namelist and to read NetCDF_file
LOGICAL        :: l_netcdf_obs = .FALSE.

INTEGER(incdf) :: STATUS
INTEGER(incdf) :: ncid
INTEGER(incdf) :: time_dimid, id_lsm
INTEGER(incdf) :: nt_in

REAL, ALLOCATABLE :: rdum1d(:)

! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

NAMELIST/vert_res/ vert_lev
NAMELIST/cntlscm/                                                    &
        nfor, model_levels_nml, l_ts_log, land_points, l_netcdf_obs

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

CALL appInit(exe_scm)

!===================================================================
IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
  WRITE(lineBuffer,'(A)') '******************** ' // RoutineName //           &
            ': Atmosphere run-time constants *******************'
  CALL umPrint(lineBuffer,src=RoutineName)
  l_print_namelist = .TRUE.
END IF
!===================================================================
! Default value of CAN_RAD_MOD which is defined in JULES_VEGETATION
!===================================================================
can_rad_mod = 1

!=====================================================================
! Assign variables normally for MPP configs
!=====================================================================
halo_i = 0   ! Halo in i
halo_j = 0   ! Halo in j
offx   = 0   ! Small halo in i.
offy   = 0   ! Small halo in j.

!=====================================================================
! Initialise STASH variables and various bits to zero
!=====================================================================
CALL init_scm_misc

!=====================================================================
!     First read in JOBDIR env var
!=====================================================================

CALL get_env_var('JOBDIR', dname, allow_missing=.TRUE., length=length)

IF (length  < 0) THEN
  icode    = -506
  cmessage =                                                          newline//&
    'Environment variable $JOBDIR not set. Using current directory'
  CALL ereport(routinename, icode, cmessage)
  dname = '.'
END IF

first_blank = LEN_TRIM(dname)+1
dname       = TRIM(dname)


!=====================================================================
!     Now read in the name of the file containing the vertical level
!     information
!=====================================================================
filename = dname(1:first_blank-1)//'/SCM_SET'

OPEN(10, FILE=Filename, IOSTAT=IstatuS, ACTION='READ', STATUS='old', &
     IOMSG=iomessage)
IF (Istatus  /=  0) THEN
  Icode = 500
  WRITE(cmessage,'(A)')                                               newline//&
    "Error opening " // TRIM(ADJUSTL(Filename)) //                    newline//&
    "IoMsg: "//TRIM(iomessage)

  CALL ereport (RoutineName, Icode, cmessage)
END IF

READ(10,vert_res)
READ(10,scm_cntl)

CLOSE(10)

!=====================================================================
!     Read SCM runtime namelist CNTLSCM
!=====================================================================

OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), ACTION='READ', IOSTAT=istatus, &
     STATUS='old', IOMSG=iomessage)

IF (istatus /= 0) THEN

  WRITE(umMessage,'(A)')                                                       &
    "Error opening " // TRIM(ADJUSTL(scm_nml)) //                     newline//&
    "IoMsg: "//TRIM(iomessage) //                                     newline//&
    "Checking for namelist.scm in current directory:"
  CALL umPrint(umMessage,src='scm_shell')

  scm_nml = './namelist.scm'
  OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus, ACTION='READ',  &
       STATUS='old', IOMSG=iomessage)

  IF (istatus /= 0) THEN
    icode = 500
    WRITE(cmessage,'(A)')                                             newline//&
      "Error opening namelist.scm:"//                                 newline//&
      "IoMsg: "//TRIM(iomessage)
    CALL ereport (routinename, icode, cmessage)
  END IF

END IF

READ(10,cntlscm)

CLOSE(10)

! If using external netcdf file, read the nc_obs namelist
! in scm_nml for information the netcdf file.
IF (l_netcdf_obs) THEN
  CALL read_nc_obs
END IF

!=====================================================================
!     Read in namelists from the SHARED file
!=====================================================================
filename = dname(1:first_blank-1)//'/SHARED'

OPEN(10, FILE=Filename, IOSTAT=Istatus, STATUS='old', ACTION='READ', &
     IOMSG=iomessage)
IF (Istatus  /=  0) THEN
  Icode = 500
  WRITE(cmessage,'(A)')                                               newline//&
    "Error opening " // TRIM(ADJUSTL(Filename)) //                    newline//&
    "IoMsg: "//TRIM(iomessage)
  CALL ereport (RoutineName, Icode, cmessage)
END IF

READ(10,nlstcall)
l_360 = lcal360        ! Set the JULES module value from UM module value

READ(10,temp_fixes)
IF (l_print_namelist) CALL print_nlist_temp_fixes()
CALL warn_temp_fixes()

! Model Domain
READ (10, model_domain)
IF (l_print_namelist) CALL print_nlist_model_domain()
CALL check_nml_model_domain()

! Planet constants
READ (10, planet_constants)
CALL set_planet_constants()
! Print needs to be after set otherwise g etc look unset for Earth simulations.
IF (l_print_namelist) CALL print_nlist_planet_constants()

READ(10,run_glomap_aeroclim)
IF (l_print_namelist) CALL print_nlist_run_glomap_clim()

READ(10,run_ukca)
IF (l_print_namelist) CALL print_nlist_RUN_UKCA() 

READ(10,run_gwd)
IF (l_print_namelist) CALL print_nlist_RUN_GWD()

! Murk not yet used by SCM - read in for consistency with UM
READ(10,RUN_Murk)
IF (l_print_namelist) CALL print_nlist_RUN_Murk()

READ(10,RUN_Convection)
IF (l_print_namelist) CALL print_nlist_RUN_Convection() 
! Check convection namelist and stop model if input values not allowed.
CALL check_run_convection()

READ(10,run_bl)
IF (l_print_namelist) CALL print_nlist_RUN_BL()
CALL check_run_bl()

! River routing not yet used by SCM - read in for consistency with UM
READ(10,RUN_Rivers)
IF (l_print_namelist) CALL print_nlist_RUN_RIVERS()

READ(10,RUN_Precip)
IF (l_print_namelist) CALL print_nlist_RUN_Precip()
CALL check_run_precip()

READ(10,RUN_Radiation)
IF (l_print_namelist) CALL print_nlist_RUN_Radiation()
CALL check_run_radiation()

READ(10,RUN_Cloud)
IF (l_print_namelist) CALL print_nlist_RUN_Cloud()
CALL check_run_cloud()

READ(10,RUN_Aerosol)
IF (l_print_namelist) CALL print_nlist_RUN_Aerosol()

! Eng_Corr not yet used by SCM - read in for consistency with UM
READ(10,RUN_Eng_Corr)
IF (l_print_namelist) CALL print_nlist_RUN_Eng_Corr()

! Unit number is passed as argument
CALL read_jules_surface_types( 10 )
CALL read_jules_surface( 10 )
CALL read_jules_radiation( 10 )
CALL read_jules_hydrology( 10 )
CALL read_jules_sea_seaice( 10 )
nice = 1              ! No. of sea ice categories
nice_use  = 1         ! No. of sea ice categories used fully in
                      ! surface exchange
CALL read_jules_soil( 10 )
CALL read_jules_vegetation( 10 )
CALL read_jules_soil_biogeochem( 10 )
CALL read_jules_snow( 10 )
CALL read_jules_urban_switches( 10 )

CLOSE(10)

!=====================================================================
! Read in model dimensions and sections from SIZES file
!=====================================================================

filename = dname(1:first_blank-1)//'/SIZES'

! DEPENDS ON: read_um_nml
CALL read_um_nml ( filename )

!=====================================================================
! Grid definitions for NewDynamics/EndGame
!=====================================================================
CALL atm_fields_bounds_init                                             &
   ( offx, offy, halo_i, halo_j, row_length, rows, rows                 &
   , tr_levels, bl_levels, ozone_levels )

!=====================================================================
!     Calculate the value to use for NTILES
!=====================================================================

ntiles=9
CALL calc_ntiles(l_aggregate,npft,nnvg,ntiles)

!=====================================================================
!     Read in the namelists from the CNTLATM file
!=====================================================================
filename = dname(1:first_blank-1)//'/CNTLATM'

OPEN(10, FILE=Filename, IOSTAT=Istatus, STATUS='old', ACTION='READ',    &
     IOMSG=iomessage)
IF (Istatus  /=  0) THEN
  Icode = 500
  WRITE(cmessage,'(A)')                                               newline//&
    "Error opening " // TRIM(ADJUSTL(Filename)) //                    newline//&
    "IoMsg: "//TRIM(iomessage)
  CALL ereport (RoutineName, Icode, cmessage)
END IF

IF ( i_bl_vn == i_bl_vn_1a ) THEN
  ! A negative "tke_levels" should be equal to bl_levels.
  IF (tke_levels < 0) THEN
    tke_levels = bl_levels
  END IF
  ! A negative "shcu_levels" should be equal to tke_levels.
  IF (shcu_levels < 0) THEN
    shcu_levels = tke_levels
  END IF
END IF
IF (ANY(rhcrit(1:cloud_levels) == rmdi)) THEN
  WRITE (cmessage,'(A62)') 'Unset values of RHcrit detected. '      &
                         //'Please check size in namelist.'
  icode = 98
  CALL ereport(RoutineName, icode, cmessage)
END IF

READ(10,tuning_segments)
IF (l_print_namelist) CALL print_nlist_tuning_segments()
CALL check_tuning_segments()
READ(10,gen_phys_inputs)
IF (l_print_namelist) CALL print_nlist_gen_phys_inputs()

READ(10,RUN_Dyn)
IF (l_print_namelist) CALL print_nlist_RUN_Dyn()
CALL check_run_dyn()

READ(10,RUN_DynTest)
IF (l_print_namelist) CALL print_nlist_RUN_Dyntest()
CALL check_run_dyntest()

READ(10,RUN_Diffusion)
IF (l_print_namelist) CALL print_nlist_RUN_Diffusion()
CALL check_run_diffusion()

! set BL parameters dependent on diffusion options
IF (l_subfilter_vert .AND. blending_option == off) THEN
  non_local_bl = off
END IF

! Set electric scheme to off for the moment as this isn't
! set up to work with the SCM.
l_use_electric = .FALSE.

! Read in Namelists R2SWCLNL and R2LWCLNL and transfer data to
! the structure SW_CONTROL and LW_CONTROL. This part of the
! code uses modules to pass arguments around.

! Set radiation aerosol switches
IF (cusack_aero==2 .OR.  cusack_aero==3    ) l_climat_aerosol    = .TRUE.
IF (cusack_aero==3 .AND. cusack_aero_hgt==1) l_clim_aero_hgt     = .TRUE.
IF (cusack_aero==2)                          l_HadGEM1_clim_aero = .TRUE.
IF (cusack_aero_hgt==2) aero_bl_levels = bl_levels

! Options for the shortwave
CALL sw_input

! Options for the longwave
CALL lw_input


! Read the JULES namelists
! Unit number is passed as argument
CALL read_jules_nvegparm( 10 )

CALL read_jules_pftparm( 10 )

CALL read_jules_triffid( 10 )

CALL read_jules_elevate( 10 )

CALL read_jules_urban2t_param( 10 )

! Initialise after everything has been read
CALL surf_couple_allocate( land_points, ntiles, sm_levels,   &
                                nice ,nice_use )

CLOSE(10)

n_internal_model = 1

! Set other dependent convective switches valid for whole run
CALL cv_set_dependent_switches

!---------------------------------------------------------------------
! Set dimensions held in scm module
!---------------------------------------------------------------------
nobs    = nfor
nbl_lv  = bl_levels
ntile   = ntiles
o3_lv   = ozone_levels
ntr_lv  = tr_levels
ntr_var = tr_vars
nlnd    = land_points
rw      = rows
rw_lng  = row_length
obs_t0  = obs_t0_prf
st_lv   = st_levels
sm_lv   = sm_levels

nml_nmod_lv = model_levels_nml

IF (nfor == imdi) nc_obsfor = .TRUE.

!-------------------------------------------------------
! Error capture
!-------------------------------------------------------

!-----------------------------------------------------------------------------
! Check settings
!-----------------------------------------------------------------------------

IF ( l_use_orog_corr  .OR. l_use_grad_corr .OR.                                &
     l_use_skyview    .OR. l_orog_unfilt    ) THEN

  WRITE(umMessage,'(A)')                                                       &
    ' ===================================================='//         newline//&
    ' | Invalid control options for SCM. Altering the    |'//         newline//&
    ' | following namelist variables:                    |'//         newline//&
    ' |                                                  |'
  CALL umPrint(umMessage,src='scm_shell')

  IF (l_use_orog_corr) THEN
    WRITE(umMessage,'(A)')                                                     &
    ' |  l_use_orog_corr = .FALSE.                       |'
    CALL umPrint(umMessage,src='scm_shell')
    l_use_orog_corr=.FALSE.
  END IF
  IF (l_use_grad_corr) THEN
    WRITE(umMessage,'(A)')                                                     &
    ' |  l_use_grad_corr = .FALSE.                       |'
    CALL umPrint(umMessage,src='scm_shell')
    l_use_grad_corr=.FALSE.
  END IF
  IF (l_use_skyview) THEN
    WRITE(umMessage,'(A)')                                                     &
    ' |  l_use_skyview = .FALSE.                         |'
    CALL umPrint(umMessage,src='scm_shell')
    l_use_skyview=.FALSE.
  END IF
  IF (l_orog_unfilt) THEN
    WRITE(umMessage,'(A)')                                                     &
    ' |  l_orog_unfilt = .FALSE.                         |'
    CALL umPrint(umMessage,src='scm_shell')
    l_orog_unfilt=.FALSE.
  END IF

  WRITE(umMessage,'(A)')                                                       &
    ' |                                                  |'//         newline//&
    ' ===================================================='
  CALL umPrint(umMessage,src='scm_shell')


END IF ! Test for invalid namelist options

IF (l_netcdf_obs) THEN

  STATUS = 0
  IF (source == 2) THEN
    STATUS = nf90_open(TRIM(netcdf_file), nf90_noWrite, ncid)

    IF (STATUS /= nf90_noerr) THEN
      WRITE(umMessage,'(A)') 'Error opening netcdf file'
      CALL umPrint(umMessage,src='scm_shell')
    END IF
    STATUS = nf90_inq_dimid (ncid, 'time', time_dimid)
    STATUS = nf90_inq_dimid (ncid, 'time', time_dimid)
    STATUS = nf90_inquire_dimension (ncid, time_dimid, LEN=nt_in)
    nfor   = nt_in - obs_t0 + 1
    nobs   = nfor

    ALLOCATE( rdum1d (nobs))

    STATUS = nf90_inq_varid (ncid, 'lsm',  id_lsm)
    STATUS = nf90_get_var   (ncid, id_lsm, rdum1d)

    IF (rdum1d(obs_t0) > 0.5) THEN
      land_points = 1
    ELSE
      land_points = 0
    END IF

    nlnd = land_points

    DEALLOCATE(rdum1d)

    STATUS = nf90_close(ncid)
  END IF

ELSE

  IF (model_levels_nml == imdi) THEN
    icode=501
    WRITE(umMessage,'(A)')                                                     &
      newline                                                                //&
      "============================================="//               newline//&
      " Number of model levels (model_levels_nml) in"//               newline//&
      " SCM namelist (&CNTLSCM) has not been set.   "//               newline//&
      "============================================="//               newline//&
      newline                                                                //&
      "Run ABORTED"// newline


    CALL umPrint(umMessage,src='scm_shell')

    WRITE(cmessage,'(A)')                                             newline//&
        'Variable MODEL_LEVELS_NML has not been set'
    CALL ereport(routinename, icode, cmessage)
  END IF

  IF (land_points > row_length*rows) THEN
    icode=502
    WRITE(umMessage,'(A)')                                                     &
      newline                                                                //&
      "============================================="//               newline//&
      " Specified number of land points greater than"//               newline//&
      " row_length*rows.                            "//               newline//&
      "============================================="//               newline//&
      newline                                                                //&
      "Run ABORTED" // newline

    CALL umPrint(umMessage,src='scm_shell')

    WRITE(cmessage,'(A)')                                             newline//&
      'Too many land points specified'
    CALL ereport(routinename, icode, cmessage)
  END IF

END IF ! l_netcdf_obs

IF (buddy_sea /= off) THEN
  WRITE(umMessage,'(A)')                                                       &
    ' ===================================================='//         newline//&
    ' | Invalid surface options for SCM                  |'//         newline//&
    ' | Altering jules_sea_seaice namelist variable      |'//         newline//&
    ' |                                                  |'//         newline//&
    ' |  Buddy_sea                                       |'//         newline//&
    ' |                                                  |'//         newline//&
    ' ===================================================='
  CALL umPrint(umMessage,src='scm_shell')

  buddy_sea = off

END IF ! Test for invalid switches

!---------------------------------------------------------------
! Check that number of Cloud levels do not exceed the
! number of model levels
!---------------------------------------------------------------

IF (cloud_levels > model_levels) THEN
  WRITE(umMessage,'(A)')                                                       &
    ' ===================================================='//         newline//&
    ' | Warning : Cloud_levels > model_levels            |'//         newline//&
    ' | Setting : Cloud_levels = model_levels            |'//         newline//&
    ' ===================================================='
  CALL umPrint(umMessage,src='scm_shell')
  cloud_levels = model_levels
END IF

IF ((tr_levels /= model_levels) .AND. &
    (tr_levels /= 0)            .AND. &
    (tr_levels /= 1)) THEN

  ! This test is enforced because at present, tracers
  ! from different areas of the UM are treated inconsistently.
  ! This causes issues for tr_levels set to anything other
  ! that model_levels
  icode    = -99
  cmessage =                                                          newline//&
    'Setting tracer levels = model levels'
  CALL ereport(routinename, icode, cmessage)

  tr_levels = model_levels
  ntr_lv    = tr_levels

END IF


!-------------------------------------------------------
! Check that if l_CCRAD is selected certain other
! switches are not set so that they conflict.
!-------------------------------------------------------
IF (l_ccrad) THEN

  IF (.NOT. l_3d_cca) THEN
    icode    = 100
    cmessage =                                                        newline//&
      'CCRad only available with 3D '                                        //&
      'cloud field (l_3d_cca = .TRUE.)'
    CALL ereport(routinename, icode, cmessage)
  END IF

  IF (l_fix_udfactor) THEN
    icode    = 101
    cmessage =                                                        newline//&
      'l_ccrad and l_fix_udfactor '                                          //&
      'should not be both set to .TRUE.'
    CALL ereport(routinename, icode, cmessage)
  END IF

  IF (l_pc2_diag_sh) THEN
    icode    = 102
    cmessage =                                                        newline//&
      'l_ccrad and l_pc2_diag_sh '                                           //&
      'should not be both set to .TRUE.'
    CALL ereport(routinename, icode, cmessage)
  END IF

END IF      ! l_ccrad

IF (l_anvil) THEN
  IF (.NOT. l_3d_cca) THEN
    icode    = 104
    cmessage =                                                        newline//&
      'Anvil scheme requires 3D CCA '                                        //&
      'field (l_3d_cca = .TRUE.)'
    CALL ereport(routinename, icode, cmessage)
  END IF
END IF

IF (i_cld_vn == i_cld_pc2) THEN
  ! So that PC2 does NOT affect Section 5 diagnostics
  l_dcpl_cld4pc2 = .TRUE.
END IF

!-------------------------------------------------------
! End error capture
!-------------------------------------------------------


!=====================================================================
!     Call the main part of the model
!=====================================================================
! DEPENDS ON: scm_main
CALL scm_main                                                               &
  ( vert_lev, nfor, l_ts_log, ntrop, sec_day                                &
  , land_points, nsprog, co2_dim_len, co2_dim_row                           &
  , cloud_levels, tr_levels                                                 &
  , tr_vars, tr_ukca, st_levels, sm_levels, bl_levels                       &
  , ozone_levels, ntiles, nice, nice_use                                    &
  ,  l_netcdf_obs )

! Shut down everything
CALL appTerminate()

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END PROGRAM scm_shell
#endif

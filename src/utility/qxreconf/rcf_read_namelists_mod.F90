! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Top level for reading in reconfiguration namelists.

MODULE Rcf_read_namelists_Mod

IMPLICIT NONE

!  Subroutine Rcf_Read_Namelists - read the rcf namelists.
!
! Description:
!   Read the namelists, assigning and freeing Fortran units as
!   required.
!
! Method:
!   The namelists are provided in various files -
!   including RECONA, SIZES, SHARED, VERTLEVS and HORIZGRID
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READ_NAMELISTS_MOD'

CONTAINS

SUBROUTINE rcf_read_namelists ( )

USE Rcf_readnl_nsubmodl_Mod,          ONLY:  Rcf_readnl_nsubmodl
USE Rcf_readnl_rundyn_Mod,            ONLY:  Rcf_readnl_rundyn
USE Rcf_readnl_runukca_Mod,           ONLY:  Rcf_readnl_runukca
USE Rcf_readnl_runglomapaeroclim_Mod, ONLY:  Rcf_readnl_runglomapaeroclim
USE Rcf_readnl_rundust_Mod,           ONLY:  Rcf_readnl_rundust
USE Rcf_readnl_rungwd_Mod,            ONLY:  Rcf_readnl_rungwd
USE Rcf_readnl_runmurk_Mod,           ONLY:  Rcf_readnl_runmurk
USE Rcf_readnl_runbl_Mod,             ONLY:  Rcf_readnl_runbl
USE Rcf_readnl_runrivers_Mod,         ONLY:  Rcf_readnl_runrivers
USE Rcf_readnl_runprecip_Mod,         ONLY:  Rcf_readnl_runprecip
USE Rcf_readnl_runelectric_Mod,       ONLY:  Rcf_readnl_runelectric
USE Rcf_readnl_lamconfig_Mod,         ONLY:  Rcf_readnl_lamconfig
USE Rcf_readnl_runozone_Mod,          ONLY:  Rcf_readnl_runozone
USE Rcf_readnl_runcloud_Mod,          ONLY:  Rcf_readnl_runcloud
USE Rcf_readnl_ancilcta_Mod,          ONLY:  Rcf_readnl_ancilcta
USE Rcf_readnl_vertical_Mod,          ONLY:  Rcf_readnl_vertical,    &
                                             print_nlist_vertical,   &
                                             check_nml_vertical
USE Rcf_readnl_horizont_Mod,          ONLY:  Rcf_readnl_horizont,         &
                                             print_nlist_horizont,        &
                                             check_nml_horizont
USE Rcf_readnl_headers_Mod,           ONLY:  Rcf_readnl_headers
USE Rcf_readnl_items_Mod,             ONLY:  Rcf_readnl_items
USE Rcf_readnl_trans_Mod,             ONLY:  Rcf_readnl_trans
USE Rcf_readnl_runstochastic_Mod,     ONLY:  Rcf_readnl_runstochastic

USE Rcf_assign_vars_nlsizes_Mod,      ONLY:  Rcf_assign_vars_nlsizes
USE rcf_assign_vars_recon_vertical_mod, ONLY: rcf_assign_vars_recon_vertical
USE Rcf_assign_vars_recon_horizontal_mod, ONLY:                           &
                                             Rcf_assign_vars_recon_horizontal

! Module for reading JULES namelists:
USE surf_couple_read_namelists_mod,   ONLY:  surf_couple_read_namelists

USE rcf_nlist_recon_technical_mod, ONLY: &
    read_nml_recon_technical,            &
    rcf_assign_vars_recon_technical,     &
    rcf_check_recon_technical,           &
    print_nlist_recon_technical,         &
    l_trans
USE rcf_readnl_recon_science_Mod,  ONLY: rcf_readnl_recon_science
USE Rcf_readnl_runaerosol_Mod,     ONLY: Rcf_readnl_runaerosol
USE Rcf_readnl_runradiation_Mod,   ONLY: Rcf_readnl_runradiation
USE Rcf_readnl_runconvection_Mod,  ONLY: Rcf_readnl_runconvection
USE Rcf_readnl_runfreetracers_Mod, ONLY: Rcf_readnl_runfreetracers
USE rcf_readnl_encorr_mod,         ONLY: rcf_readnl_run_eng_corr
USE rcf_readnl_model_domain_mod,   ONLY: rcf_readnl_model_domain
USE rcf_readnl_coupling_mod,       ONLY: rcf_readnl_coupling
USE rcf_readnl_nlstcall_mod,       ONLY: rcf_readnl_nlstcall
USE rcf_readnl_carbon_options_mod, ONLY: rcf_readnl_carbon_options
USE rcf_readnl_gen_phys_inputs_mod, ONLY: rcf_readnl_gen_phys_inputs
USE rcf_readnl_recon_idealised_mod, ONLY: rcf_readnl_recon_idealised
USE idealise_run_mod,              ONLY: &
    read_nml_idealised,                  &
    print_nlist_idealised,               &
    check_nlist_idealised
USE easyaerosol_option_mod, ONLY: read_nml_easyaerosol, &
                                  print_nlist_easyaerosol 

USE vertnamelist_mod, ONLY:  read_nml_vertlevs,                          &
                             check_nml_vertlevs,                         &
                             print_nlist_vertlevs

USE vrhoriz_grid_mod, ONLY    : read_nml_horizgrid,                       &
                                print_nlist_horizgrid,                    &
                                check_nml_horizgrid

USE planet_constants_mod, ONLY: read_nml_planet_constants, &
                                set_planet_constants,      &
                                print_nlist_planet_constants
USE Rcf_readnl_rundyntest_Mod,    ONLY:  Rcf_readnl_rundyntest
USE dynamics_testing_mod,         ONLY: problem_number
USE problem_mod,                  ONLY: standard

USE science_fixes_mod
USE check_iostat_mod
USE missing_data_mod, ONLY: imdi

USE calc_ntiles_mod, ONLY:                                              &
  calc_ntiles

USE nlsizes_namelist_mod, ONLY:                                          &
    vert_lev,                                                            &
    var_grid

USE jules_surface_mod, ONLY:                                             &
  l_aggregate

USE jules_surface_types_mod, ONLY: npft, nnvg
USE land_tile_ids, ONLY: set_tile_id_arrays

USE nlsizes_namelist_mod, ONLY: &
    read_nml_nlsizes,           &
    print_nlist_nlsizes,        &
    check_nlsizes,              &
    ntiles

USE Ereport_Mod, ONLY:&
    Ereport

USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    newline,                &
    PrintStatus,            &
    PrStatus_Oper

USE filenamelength_mod, ONLY: &
    filenamelength

USE file_manager, ONLY:  &
    assign_file_unit,    &
    release_file_unit

USE nlcfiles_namelist_mod, ONLY: read_nml_nlcfiles, print_nlist_nlcfiles

USE get_env_var_mod, ONLY: get_env_var
USE errormessagelength_mod, ONLY: errormessagelength

USE um_parcore, ONLY: mype

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

USE model_domain_mod, ONLY: l_regular

IMPLICIT NONE
! Arguments

! Local Variables/Paramters
CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'RCF_READ_NAMELISTS'
CHARACTER (LEN=errormessagelength) :: Cmessage
CHARACTER (LEN=errormessagelength) :: iomessage
CHARACTER (LEN=filenamelength)     :: FileName
INTEGER                            :: ErrorStatus
INTEGER                            :: nft
INTEGER                            :: nft_sizes
INTEGER                            :: nft_shared
INTEGER                            :: nft_vertlevs  ! File unit - VERTLEVS
INTEGER                            :: nft_horizgrid  ! File unit - HORIZGRID
INTEGER                            :: iostatus
LOGICAL                            :: l_print_namelists = .FALSE.
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  l_print_namelists = .TRUE.
END IF
! ---------------------------------------------------
! Read namelists common to the UM and reconfiguration
! Part 1 - miscellaneous namelists in file SHARED
! ---------------------------------------------------
! First, find the namelist filename from Env Vars
CALL get_env_var( 'SHARED_NLIST', FileName )

FileName = TRIM( FileName )

! Open the file containing namelists from SHARED file
IF ( mype == 0 ) THEN
  CALL assign_file_unit(FileName, nft_shared, handler="fortran")
  OPEN( UNIT=nft_shared, FILE=FileName, ACTION='READ',   &
                    IOSTAT=iostatus, IOMSG=iomessage )
  IF ( iostatus /= 0 ) THEN
    cmessage = 'Error opening Shared namelists file:'    &
               //FileName // ' :' // newline //  TRIM(iomessage)
    CALL Ereport( RoutineName, iostatus, Cmessage )
  END IF
END IF

IF ( PrintStatus >= PrStatus_Oper ) THEN
  WRITE(umMessage,'(''Shared namelists file: '',a)') FileName
  CALL umPrint(umMessage,src='rcf_read_namelists_mod')
END IF

! Read in namelist containing name of output dump.
CALL read_nml_nlcfiles( nft_shared )
IF ( l_print_namelists ) CALL print_nlist_nlcfiles()

CALL rcf_readnl_nlstcall (nft_shared)

! ---------------------------------------------------
! read in any temporary UM fixes and warn if required
! ---------------------------------------------------

CALL read_nml_temp_fixes(nft_shared)
CALL warn_temp_fixes()

! ---------------------------------------------------
CALL rcf_readnl_carbon_options(nft_shared)
CALL rcf_readnl_coupling(nft_shared)
CALL rcf_readnl_model_domain (nft_shared)

! Planet constants
CALL read_nml_planet_constants(nft_shared)
CALL set_planet_constants()
! Print needs to be after set otherwise g etc look unset for Earth simulations.
IF (l_print_namelists)  CALL print_nlist_planet_constants

! Read the JULES namelists
CALL surf_couple_read_namelists("RECON", nft_shared, nft_shared)

! Rewind the unit when we are done with the JULES namelists
IF ( mype == 0 ) REWIND(nft_shared)

! ---------------------------------------------------
! Read namelists common to the UM and reconfiguration
! Part 2 - model SIZES namelists
! ---------------------------------------------------

! First, find the namelist filename from Env Vars
CALL get_env_var( 'SIZES_NLIST', FileName )

FileName = TRIM( FileName )

! Open the file containing namelists from SIZES file
IF ( mype == 0 ) THEN
  CALL assign_file_unit( FileName, nft_sizes, handler="fortran")
  OPEN( UNIT=nft_sizes, FILE=FileName, ACTION='READ',               &
                        IOSTAT=iostatus, IOMSG=iomessage )
  IF ( iostatus /= 0 ) THEN
    cmessage = 'Error opening Sizes namelists file:'//FileName //  &
               ' :' // newline // TRIM(iomessage)
    CALL Ereport( RoutineName, iostatus, Cmessage )
  END IF
END IF

IF ( PrintStatus >= PrStatus_Oper ) THEN
  WRITE(umMessage,'(''Sizes namelists file: '',a)') FileName
  CALL umPrint(umMessage,src='rcf_read_namelists_mod')
END IF

CALL rcf_readnl_nsubmodl (nft_sizes)

CALL read_nml_nlsizes(nft_sizes)
IF (l_print_namelists)  CALL print_nlist_nlsizes()
CALL check_nlsizes()
CALL rcf_assign_vars_nlsizes()

! ------------------------------------------------
! Begin reading reconfiguration-specific namelists
! ------------------------------------------------

! First, find the namelist filename from Env Vars
CALL get_env_var( 'RCF_NAMELIST', FileName )

FileName = TRIM( FileName )

! Open the file containing namelists from RECONA file
IF ( mype == 0 ) THEN
  CALL assign_file_unit( Filename, nft, handler="fortran")
  OPEN( UNIT=nft, FILE=FileName, ACTION='READ',                    &
                  IOSTAT=iostatus, IOMSG=iomessage )
  IF ( iostatus /= 0 ) THEN
    cmessage = 'Error opening Namelist file:'//FileName //  &
               ' :' // newline //  TRIM(iomessage)
    CALL Ereport( RoutineName, iostatus, Cmessage )
  END IF
END IF

IF ( PrintStatus >= PrStatus_Oper ) THEN
  WRITE(umMessage,'(''Namelist file: '',a)') FileName
  CALL umPrint(umMessage,src='rcf_read_namelists_mod')
END IF

! Open the file containing vertical levels namelist
! This must be done after NLSIZES as that's where the filename is specified.
CALL assign_file_unit( vert_lev, nft_vertlevs, handler="fortran" )
OPEN( UNIT=nft_vertlevs, FILE=vert_lev, ACTION='READ',             &
                         IOSTAT=iostatus, IOMSG=iomessage )
IF ( iostatus /= 0 ) THEN
  WRITE(umMessage,'(''Vertical Levels file: '',a)') vert_lev
  CALL umPrint(umMessage,src='rcf_readnl_vertical_Mod')
  cmessage = 'Error opening Vertical Levels namelist file:'                &
                   //TRIM(vert_lev) // ' :' // newline //  TRIM(iomessage)
  CALL Ereport( RoutineName, iostatus, Cmessage )
END IF

IF ( PrintStatus >= PrStatus_Oper ) THEN
  WRITE(umMessage,'(''Vertical Levels file: '',a)') TRIM(vert_lev)
  CALL umPrint(umMessage,src='rcf_read_namelists_mod')
END IF

IF (.NOT. l_regular) THEN
  IF ( mype == 0 )  THEN
    CALL assign_file_unit( var_grid, nft_horizgrid, handler="fortran" )
    OPEN(UNIT=nft_horizgrid, FILE=var_grid, ACTION='READ',         &
                            IOSTAT=iostatus, IOMSG=iomessage)
    IF ( iostatus /= 0 ) THEN
      WRITE(umMessage,'(2A)') 'Horizontal grid file: ', TRIM(var_grid)
      CALL umPrint(umMessage,src='rcf_readnl_horizont_mod')
      cmessage = 'Error opening Horizontal grid namelist file:'   &
                   //TRIM(var_grid) // ' :' // newline // TRIM(iomessage)
      CALL Ereport( RoutineName, iostatus, Cmessage )
    END IF
  END IF

  ! Write out namelist for diagnostic
  IF ( PrintStatus >= PrStatus_Oper ) THEN
    IF ( mype == 0 ) THEN
      WRITE(umMessage,'(2A)') 'horizontal grid file: ',TRIM(var_grid)
      CALL umPrint(umMessage,src='rcf_readnl_horizont_mod')
    END IF
  END IF
END IF

! RECON must be read after NLSIZES due to using model_levels
CALL read_nml_recon_technical (nft)
IF (l_print_namelists)  CALL print_nlist_recon_technical()
CALL rcf_check_recon_technical()
CALL rcf_assign_vars_recon_technical()

CALL rcf_readnl_recon_science  (nft)

! -----------------------------------------
! Resume reading SHARED and SIZES namelists
! -----------------------------------------
CALL calc_ntiles(l_aggregate,npft,nnvg,ntiles)
CALL set_tile_id_arrays ( )

CALL rcf_readnl_rundust  (nft_shared)
CALL rcf_readnl_runglomapaeroclim (nft_shared)
CALL rcf_readnl_runukca (nft_shared)
CALL rcf_readnl_rungwd (nft_shared)
CALL rcf_readnl_runmurk (nft_shared)
CALL rcf_readnl_runconvection (nft_shared)
CALL rcf_readnl_runbl (nft_shared)
CALL rcf_readnl_runrivers (nft_shared)
CALL rcf_readnl_runprecip (nft_shared)
CALL rcf_readnl_runradiation (nft_shared)
CALL rcf_readnl_runcloud (nft_shared)
CALL rcf_readnl_runaerosol (nft_shared)
CALL rcf_readnl_lamconfig (nft_shared)
CALL rcf_readnl_runozone (nft_shared)
CALL rcf_readnl_runfreetracers (nft_shared)
CALL rcf_readnl_run_eng_corr (nft_shared)
CALL rcf_readnl_gen_phys_inputs (nft_shared)
CALL rcf_readnl_rundyn (nft_shared)
CALL rcf_readnl_rundyntest(nft_shared)
CALL read_nml_easyaerosol (nft_shared)
IF (l_print_namelists) CALL print_nlist_easyaerosol()
CALL rcf_readnl_runstochastic (nft_shared)
CALL rcf_readnl_runelectric (nft_shared)

! For vertical, VERTICAL and VERTLEVS are read in
CALL rcf_readnl_vertical (nft)
IF (l_print_namelists)  CALL print_nlist_vertical()
CALL check_nml_vertical()
! Read VERTLEVS Namelist
CALL read_nml_vertlevs(nft_vertlevs)
IF (l_print_namelists)  CALL print_nlist_vertlevs()
CALL check_nml_vertlevs()
CALL rcf_assign_vars_recon_vertical()
! For horizontal, HORIZONT and HORIZGRID are read in
CALL rcf_readnl_horizont (nft)
IF (l_print_namelists) CALL print_nlist_horizont ()
CALL check_nml_horizont()
IF (.NOT. l_regular) THEN
  CALL read_nml_horizgrid (nft_horizgrid)
  IF (l_print_namelists) CALL print_nlist_horizgrid ()
  CALL check_nml_horizgrid()
END IF
CALL rcf_assign_vars_recon_horizontal()
CALL rcf_readnl_headers  (nft)

IF (l_trans) THEN
  CALL rcf_readnl_trans  (nft)
END IF

CALL rcf_readnl_ancilcta (nft_shared)
CALL rcf_readnl_items    (nft_shared)

! ---------------------------------------------------
! Read namelists common to the UM and reconfiguration
! Part 3 - idealised namelists
! ---------------------------------------------------
IF (problem_number /= standard) THEN
  CALL read_nml_idealised()
  IF (l_print_namelists) CALL print_nlist_idealised()
  CALL check_nlist_idealised()

  CALL rcf_readnl_recon_idealised (nft)
END IF

IF ( mype == 0 ) THEN
  CLOSE( UNIT=nft )
  CLOSE( UNIT=nft_sizes )
  CLOSE( UNIT=nft_shared )
  CLOSE( UNIT=nft_vertlevs )

  CALL release_file_unit( nft, handler="fortran" )
  CALL release_file_unit( nft_sizes, handler="fortran" )
  CALL release_file_unit( nft_shared, handler="fortran" )
  CALL release_file_unit( nft_vertlevs, handler="fortran" )
  IF (.NOT. l_regular) THEN
    CLOSE( UNIT=nft_horizgrid )
    CALL release_file_unit( nft_horizgrid, handler="fortran" )
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE  Rcf_read_namelists

END MODULE Rcf_read_namelists_Mod

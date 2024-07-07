! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   A module containing tunable parameters and dependent logicals
!   used in the aerosol scheme
!
MODULE run_aerosol_mod

!
! Description:
!   This module contains declarations for tunable parameters
!   used to set the emission heights/levels of manmade aerosols
!   and dependent Classic aerosol logicals
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: 17 CLASSIC Aerosols
!
! Code description:
!   Language: Fortran 90
!   This code is written to UM programming standards UMDP 3 vn8.2
!
USE missing_data_mod,      ONLY: imdi
USE ereport_mod,           ONLY: ereport
USE umPrintMgr
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!
! Control logicals for Classic aerosol options
! ============================================
!
! Sulphur Cycle scheme
LOGICAL :: l_sulpc_so2 = .FALSE.  ! S Cycle: SO2 MMR included
! with emission options
LOGICAL :: l_so2_surfem = .FALSE. ! SO2 Surface Emissions
LOGICAL :: l_so2_hilem = .FALSE.  ! SO2 High Level Emissions
INTEGER :: SO2_high_level = imdi   ! for chimney SO2 emission
LOGICAL :: l_so2_natem = .FALSE.  ! SO2 Natural Emissions
! Dimethyl sulphide
LOGICAL :: l_sulpc_dms = .FALSE.  ! S Cycle: DMS MMR included
LOGICAL :: l_dms_em = .FALSE.     ! DMS Emissions
LOGICAL :: l_dms_em_inter = .FALSE. ! Interactive DMS Emissions
! Switch to choose scheme to use for interactive sea-air exchange of DMS
INTEGER :: i_dms_flux = imdi
! Ozone included for oxidation of DMS and SO2
LOGICAL :: l_sulpc_ozone = .FALSE. ! S Cycle: Ozone oxidation
! Online oxidants from UKCA
LOGICAL :: l_sulpc_online_oxidants = .FALSE. ! S Cycle : Online oxidants
! Depleted oxidants are passed back to UKCA (.FALSE. is deemed unsafe)
LOGICAL :: l_sulpc_2_way_coupling = .TRUE.  ! S Cycle : Depleted oxidants
! SO2+O3 reaction NOT buffered by NH3
LOGICAL :: l_sulpc_so2_o3_nonbuffered = .FALSE. ! S Cycle: SO2+O3 reaction
! Ammonia
LOGICAL :: l_sulpc_nh3 = .FALSE.  ! S Cycle: NH3 tracer included
LOGICAL :: l_nh3_em = .FALSE.     ! S Cycle: NH3 emissions

! Soot scheme
LOGICAL :: l_soot = .FALSE.  ! Soot scheme
! with emission options
LOGICAL :: l_soot_surem = .FALSE. ! Soot Surface Emissions
LOGICAL :: l_soot_hilem = .FALSE.  ! Soot High Level Emissions

! Biomass scheme
LOGICAL :: l_biomass = .FALSE.  ! Biomass scheme
! with emission options
LOGICAL :: l_bmass_surem = .FALSE. ! Biomass Surface Emissions
LOGICAL :: l_bmass_hilem = .FALSE.  ! Biomass High Level Emissions

! OCFF (Organic Carbon from Fossil Fuels) scheme
LOGICAL :: l_ocff = .FALSE.  ! OCFF scheme
! with emission options
LOGICAL :: l_ocff_surem = .FALSE.  ! OCFF Surface Emissions
LOGICAL :: l_ocff_hilem = .FALSE.  ! OCFF High Level Emissions

! Ammonium nitrate aerosol
LOGICAL :: l_nitrate = .FALSE.     ! Ammonium nitrate aerosol

! Additional Constraints

! Use sulphate aerosol no. in S-cycle
LOGICAL :: l_use_sulphate_sulpc = .FALSE.

! Use biomass aerosol no. in S-cycle
LOGICAL :: l_use_bmass_sulpc = .FALSE.

! Use organic carbon fossil fuel aerosol no. in S-cycle
LOGICAL :: l_use_ocff_sulpc = .FALSE.

! Use ammonium nitrate aerosol no. in S-cycle
LOGICAL :: l_use_nitrate_sulpc = .FALSE.

! Use sea-salt aerosol no. in S-cycle
LOGICAL :: l_use_seasalt_sulpc = .FALSE.

! Use sea-salt aerosol no. in PM diagnostics (unconditional choice)
LOGICAL :: l_use_seasalt_pm = .FALSE.

! Use diurnal and weekly emission cycles for primary emissions in CLASSIC
LOGICAL :: l_temporal_emi = .FALSE.

! Use variable injection heights for biomass burning emissions
LOGICAL :: l_bmass_hilem_variable = .FALSE.

!
! Control logicals for Classic aerosol LBC options
! ================================================
!
! Sulphur Cycle scheme
LOGICAL :: l_so2_lbc = .FALSE.  ! S Cycle LBCs included
! Dimethyl sulphide
LOGICAL :: l_dms_lbc = .FALSE.  ! S Cycle: DMS LBC included
! Ammonia
LOGICAL :: l_nh3_lbc = .FALSE.  ! S Cycle: NH3 LBC included
! Soot scheme
LOGICAL :: l_soot_lbc = .FALSE.  ! Soot scheme LBCs included
! Biomass scheme
LOGICAL :: l_bmass_lbc = .FALSE.  ! Biomass scheme LBCs included
! OCFF (Organic Carbon from Fossil Fuels) scheme
LOGICAL :: l_ocff_lbc = .FALSE.  ! OCFF scheme LBCs included
! Ammonium nitrate aerosol
LOGICAL :: l_nitr_lbc = .FALSE.     ! Nitrate scheme LBCs included

!
! Dependent logicals for Classic aerosol options
! ==============================================
!
LOGICAL :: l_so2 = .FALSE.
LOGICAL :: l_so4_aitken = .FALSE.
LOGICAL :: l_so4_accu = .FALSE.
LOGICAL :: l_so4_diss = .FALSE.
LOGICAL :: l_dms = .FALSE.
LOGICAL :: l_nh3 = .FALSE.

LOGICAL :: l_soot_new = .FALSE.
LOGICAL :: l_soot_agd = .FALSE.
LOGICAL :: l_soot_cld = .FALSE.

LOGICAL :: l_bmass_new = .FALSE.
LOGICAL :: l_bmass_agd = .FALSE.
LOGICAL :: l_bmass_cld = .FALSE.

LOGICAL :: l_ocff_new = .FALSE.
LOGICAL :: l_ocff_agd = .FALSE.
LOGICAL :: l_ocff_cld = .FALSE.

LOGICAL :: l_nitr_acc  = .FALSE.
LOGICAL :: l_nitr_diss = .FALSE.

!
! Dependent logicals for Classic aerosol LBC options
! ==================================================
!
LOGICAL :: l_so4_aitken_lbc = .FALSE.
LOGICAL :: l_so4_accu_lbc = .FALSE.
LOGICAL :: l_so4_diss_lbc = .FALSE.

LOGICAL :: l_soot_new_lbc = .FALSE.
LOGICAL :: l_soot_agd_lbc = .FALSE.
LOGICAL :: l_soot_cld_lbc = .FALSE.

LOGICAL :: l_bmass_new_lbc = .FALSE.
LOGICAL :: l_bmass_agd_lbc = .FALSE.
LOGICAL :: l_bmass_cld_lbc = .FALSE.

LOGICAL :: l_ocff_new_lbc = .FALSE.
LOGICAL :: l_ocff_agd_lbc = .FALSE.
LOGICAL :: l_ocff_cld_lbc = .FALSE.

LOGICAL :: l_nitr_acc_lbc  = .FALSE.
LOGICAL :: l_nitr_diss_lbc = .FALSE.

!
! Heights, i.e. levels, used for manmade aerosol emissions
! ========================================================
!
      ! Aerosol Modelling - model level heights...

INTEGER :: bmass_high_level_1 = imdi ! Lowest and highest
INTEGER :: bmass_high_level_2 = imdi ! for biomass emissions.
INTEGER :: soot_high_level = imdi  ! for chimney soot emission
INTEGER :: ocff_high_level = imdi  ! for chimney OCFF emission

!
! Fixed values no longer in namelist
! ==================================
!
      ! improved nonhydrostatic weights in tracer1 calculation
LOGICAL, PARAMETER :: L_tracer1_non_hydro = .TRUE. ! not in namelist

! For Aero_Ctl (Sulph cycle or Soot) No.of times chem called per timestep
INTEGER, PARAMETER :: call_chem_freq = 1 ! Fixed; no longer in namelist

! =======================================================================
!
NAMELIST/RUN_Aerosol/   &
  l_sulpc_so2, l_so2_surfem, l_so2_hilem, l_so2_natem,  &
  SO2_high_level,  &
  l_sulpc_dms, l_dms_em, l_dms_em_inter, i_dms_flux,  &
  l_sulpc_ozone, l_sulpc_online_oxidants, l_sulpc_2_way_coupling,  &
  l_sulpc_so2_o3_nonbuffered, l_sulpc_nh3, l_nh3_em,  &
  l_soot, l_soot_surem, l_soot_hilem, soot_high_level,  &
  l_biomass, l_bmass_surem, l_bmass_hilem, bmass_high_level_1,  &
  bmass_high_level_2,  &
  l_ocff, l_ocff_surem, l_ocff_hilem, ocff_high_level, l_nitrate,  &
  l_use_sulphate_sulpc, l_use_nitrate_sulpc, l_use_seasalt_sulpc,  &
  l_use_seasalt_pm, l_temporal_emi, l_use_bmass_sulpc, l_use_ocff_sulpc,  &
  l_so2_lbc, l_dms_lbc, l_nh3_lbc, l_soot_lbc, l_bmass_lbc,  &
  l_ocff_lbc, l_nitr_lbc, l_bmass_hilem_variable

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RUN_AEROSOL_MOD'

CONTAINS
!
! Internal subroutine to check that entries to RUN_aerosol are consistent
! =======================================================================
!
SUBROUTINE run_aerosol_check( )

USE UM_ParCore, ONLY: mype
USE umPrintMgr, ONLY: umPrint, umMessage
USE ukca_option_mod, ONLY: i_ukca_chem, l_ukca_classic_hetchem
USE ukca_chem_schemes_mod
USE rad_input_mod, ONLY: l_use_sulpc_direct, l_use_sulpc_indirect_lw, &
                         l_use_sulpc_indirect_sw, l_use_seasalt_indirect,&
                         l_use_nitrate_direct, l_use_nitrate_indirect, &
                         l_use_soot_direct, l_use_soot_indirect, &
                         l_use_bmass_direct, l_use_bmass_indirect, &
                         l_use_ocff_direct, l_use_ocff_indirect
!
USE nlsizes_namelist_mod, ONLY: &
    model_levels

IMPLICIT NONE
!
!
CHARACTER (LEN=errormessagelength)           :: cmessage
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'RUN_AEROSOL_CHECK'
INTEGER :: icode
LOGICAL :: l_ukca_chem
REAL(KIND=jprb) :: zhook_handle
!

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check there is a valid UKCA scheme active
!
l_ukca_chem = (i_ukca_chem==i_ukca_chem_trop)      .OR.  &
              (i_ukca_chem==i_ukca_chem_raq)       .OR.  &
              (i_ukca_chem==i_ukca_chem_tropisop)  .OR.  &
              (i_ukca_chem==i_ukca_chem_strattrop) .OR.  &
              (i_ukca_chem==i_ukca_chem_strat)

! Check that all user-defined schemes and emissions heights are sensible
!
! Sulphur scheme
!
icode = 0

IF (l_ukca_classic_hetchem .AND. .NOT. l_sulpc_so2) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)')'Error -RUN_aerosol- Cannot run ' // &
                          'l_ukca_classic_hetchem without sulpc_so2'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_sulpc_dms .AND. .NOT. l_sulpc_so2) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)')'Error -RUN_aerosol- sulpc_dms but not sulpc_so2'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
ELSE
  l_dms = l_sulpc_dms
END IF

IF (l_dms_lbc .AND. .NOT. l_dms) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error -RUN_aerosol- dms_lbc but not dms'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_dms_em .AND. .NOT. l_sulpc_dms) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error -RUN_aerosol- dms_em but not sulpc_dms'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_dms_em_inter .AND. .NOT. l_dms_em) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)')'Error -RUN_aerosol- dms_em_inter but not dms_em'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (i_dms_flux/=imdi .AND. .NOT. l_dms_em_inter) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)')'Error-RUN_aerosol-i_dms_flux but not dms_em_inter'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_sulpc_ozone .AND. .NOT. l_sulpc_so2) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)')'Error-RUN_aerosol-sulpc_ozone but not sulpc_so2'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_sulpc_nh3 .AND. .NOT. l_sulpc_so2) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error -RUN_aerosol- sulpc_nh3 but not sulpc_so2'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
ELSE
  l_nh3 = l_sulpc_nh3
END IF

IF (l_nh3_lbc .AND. .NOT. l_nh3) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error - RUN_aerosol - nh3_lbc but not nh3'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_nh3_em .AND. .NOT. l_sulpc_nh3) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error - RUN_aerosol - nh3_em but not sulpc_nh3'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (.NOT. l_sulpc_ozone .OR. .NOT. l_sulpc_nh3) THEN
  IF (l_sulpc_so2_o3_nonbuffered) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A)') 'Error -RUN_aerosol- sulpc_so2_o3_nonbuffered'
      CALL umPrint(umMessage,src='run_aerosol_check')
      WRITE(umMessage,'(A)')'not allowed unless sulpc_ozone/sulpc_nh3 both true'
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  END IF
END IF

IF (l_so2_surfem .AND. .NOT. l_sulpc_so2) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)')'Error -RUN_aerosol- so2_surfem but not sulpc_so2'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_so2_hilem .AND. .NOT. l_sulpc_so2) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error -RUN_aerosol- so2_hilem but not sulpc_so2'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF
IF (l_so2_hilem) THEN
  IF (SO2_high_level == imdi) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A)')'Error-RUN_aerosol-SO2 emission but level not set'
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  ELSE IF (SO2_high_level < 1 .OR. SO2_high_level > model_levels) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A,2I6)') 'Error -RUN_aerosol- Query SO2_high_level',&
                          SO2_high_level,model_levels
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  END IF
END IF


IF (l_so2_natem .AND. .NOT. l_sulpc_so2) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error -RUN_aerosol- so2_natem but not sulpc_so2'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_sulpc_so2) THEN
  ! Dependent logicals normally all needed if l_sulpc_so2 is true
  l_so2 = .TRUE.
  l_so4_aitken = .TRUE.
  l_so4_accu = .TRUE.
  l_so4_diss = .TRUE.
END IF

IF (l_sulpc_so2 .AND. l_so2_lbc) THEN
  ! Dependent logicals normally all needed if l_so2_lbc also true
  l_so4_aitken_lbc = .TRUE.
  l_so4_accu_lbc = .TRUE.
  l_so4_diss_lbc = .TRUE.
END IF

IF (l_sulpc_online_oxidants .AND. (.NOT. l_ukca_chem .OR. .NOT. l_sulpc_so2)) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error-RUN_aerosol-sulpc_online_oxidants'
    CALL umPrint(umMessage,src='run_aerosol_check')
    WRITE(umMessage,'(A)') ' but need both ukca_chem and sulpc_so2'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_use_sulphate_sulpc .AND. .NOT. l_sulpc_so2) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)')'Error-RUN_aerosol-use_sulphate_sulpc but not sulpc_so2'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

!
! SOOT scheme
!
IF (l_soot) THEN
  ! Dependent logicals normally all needed if l_soot is true
  l_soot_new = .TRUE.
  l_soot_agd = .TRUE.
  l_soot_cld = .TRUE.
END IF
IF (l_soot .AND. l_soot_lbc) THEN
  ! Dependent logicals normally all needed if l_soot_lbc also true
  l_soot_new_lbc = .TRUE.
  l_soot_agd_lbc = .TRUE.
  l_soot_cld_lbc = .TRUE.
END IF
!
IF (l_soot_surem .AND. .NOT. l_soot) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error - RUN_aerosol - soot_surem but not soot'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_soot_hilem .AND. .NOT. l_soot) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error - RUN_aerosol - soot_hilem but not soot'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (L_Soot_hilem) THEN
  IF (Soot_high_level == imdi) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A)')'Error-RUN_aerosol-Soot emission but level not set'
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  ELSE IF (Soot_high_level < 1 .OR. Soot_high_level > model_levels) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A,2I6)') 'Error -RUN_aerosol- Query Soot_high_level',&
                          Soot_high_level,model_levels
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  END IF
END IF

!
! BIOMASS scheme
!
IF (l_biomass) THEN
  ! Dependent logicals normally all needed if l_biomass is true
  l_bmass_new = .TRUE.
  l_bmass_agd = .TRUE.
  l_bmass_cld = .TRUE.
END IF
IF (l_biomass .AND. l_bmass_lbc) THEN
  ! Dependent logicals normally all needed if l_bmass_lbc also true
  l_bmass_new_lbc = .TRUE.
  l_bmass_agd_lbc = .TRUE.
  l_bmass_cld_lbc = .TRUE.
END IF
!
IF (l_bmass_surem .AND. .NOT. l_biomass) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error -RUN_aerosol- bmass_surem but not biomass'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_bmass_hilem .AND. .NOT. l_biomass) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error -RUN_aerosol- bmass_hilem but not biomass'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

! If injection heights for biomass burning emissions are fixed then check that
! the first model layer to inject them has been set correctly by the user.
IF (L_bmass_hilem .AND. .NOT. (l_bmass_hilem_variable)) THEN
  IF (bmass_high_level_1 == imdi) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A)')'Error-RUN_aerosol-bmass emission level 1 not set'
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  ELSE IF (bmass_high_level_1<1 .OR. bmass_high_level_1>model_levels) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A,A,2I6)') 'Error - RUN_aerosol - ', &
           'Query bmass_high_level_1',bmass_high_level_1,model_levels
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  END IF
END IF
!
! If injection heights for biomass burning emissions are fixed then check that
! the last model layer to inject them has been set correctly by the user.
IF (L_bmass_hilem .AND. .NOT. (l_bmass_hilem_variable)) THEN
  IF (bmass_high_level_2 == imdi) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A)')'Error-RUN_aerosol-bmass emission level 2 not set'
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  ELSE IF (bmass_high_level_2<1 .OR. bmass_high_level_2>model_levels) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A,A,2I6)') 'Error - RUN_aerosol - ', &
           'Query bmass_high_level_2',bmass_high_level_2,model_levels
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  END IF
END IF

IF (l_use_bmass_sulpc .AND. (.NOT. l_biomass .OR. .NOT. l_use_sulphate_sulpc)) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error - RUN_aerosol - use_bmass_sulpc '
    CALL umPrint(umMessage,src='run_aerosol_check')
    WRITE(umMessage,'(A)') 'but not biomass/use_sulphate_sulpc'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

!
! Organic Carbon Fossil Fuel scheme
!
IF (l_ocff) THEN
  ! Dependent logicals normally all needed if l_ocff is true
  l_ocff_new = .TRUE.
  l_ocff_agd = .TRUE.
  l_ocff_cld = .TRUE.
END IF
IF (l_ocff .AND. l_ocff_lbc) THEN
  ! Dependent logicals normally all needed if l_ocff_lbc also true
  l_ocff_new_lbc = .TRUE.
  l_ocff_agd_lbc = .TRUE.
  l_ocff_cld_lbc = .TRUE.
END IF

IF (l_ocff_surem .AND. .NOT. l_ocff) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error - RUN_aerosol - ocff_surem but not ocff'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_ocff_hilem .AND. .NOT. l_ocff) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error - RUN_aerosol - ocff_hilem but not ocff'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (L_ocff_hilem) THEN
  IF (ocff_high_level == imdi) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A)')'Error-RUN_aerosol-ocff emission but level not set'
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  ELSE IF (ocff_high_level < 1 .OR. ocff_high_level > model_levels) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A,2I6)')'Error-RUN_aerosol-Query ocff_high_level',&
                          ocff_high_level,model_levels
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  END IF
END IF

IF (l_use_ocff_sulpc .AND. (.NOT. l_ocff .OR. .NOT. l_use_sulphate_sulpc)) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error - RUN_aerosol - use_ocff_sulpc '
    CALL umPrint(umMessage,src='run_aerosol_check')
    WRITE(umMessage,'(A)') 'but not ocff/use_sulphate_sulpc'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

!
! Ammonium nitrate aerosol scheme
!
IF (l_nitrate .AND. (.NOT. l_nh3 .OR. .NOT. l_ukca_chem)) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)')'Error -RUN_aerosol- nitrate but not nh3/ukca_chem'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF

IF (l_nitrate) THEN
  ! Dependent logicals normally all needed if l_nitrate is true
  l_nitr_acc  = .TRUE.
  l_nitr_diss = .TRUE.
END IF
IF (l_nitrate .AND. l_nitr_lbc) THEN
  ! Dependent logicals normally all needed if l_nitr_lbc also true
  l_nitr_acc_lbc  = .TRUE.
  l_nitr_diss_lbc = .TRUE.
END IF

IF (.NOT. l_nitrate .OR. .NOT. l_use_sulphate_sulpc) THEN
  IF (l_use_nitrate_sulpc) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,'(A)')'Error-RUN_aerosol-use_nitrate_sulpc not allowed'
      CALL umPrint(umMessage,src='run_aerosol_check')
      WRITE(umMessage,'(A)') 'unless nitrate/use_sulphate_sulpc both true'
      CALL umPrint(umMessage,src='run_aerosol_check')
    END IF
    icode = icode + 1
  END IF
END IF

!
! Sea Salt scheme
!
IF (l_use_seasalt_sulpc .AND. .NOT. l_use_sulphate_sulpc) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') 'Error - RUN_aerosol - use_seasalt_sulpc '
    CALL umPrint(umMessage,src='run_aerosol_check')
    WRITE(umMessage,'(A)') 'but not use_sulphate_sulpc'
    CALL umPrint(umMessage,src='run_aerosol_check')
  END IF
  icode = icode + 1
END IF


IF (icode /= 0) THEN
  cmessage='*** Error(s) in levels or logicals set in RUN_Aerosol namelist'
  CALL ereport(RoutineName, icode, cmessage)
END IF

! Also check validity of aerosol radiative effect choices here,
! as this needs to be done after run_aerosol has been read.
IF (.NOT. l_sulpc_so2) THEN
  l_use_sulpc_direct = .FALSE.
  l_use_sulpc_indirect_lw = .FALSE.
  l_use_sulpc_indirect_sw = .FALSE.
END IF
IF (.NOT. l_use_sulpc_indirect_lw .OR. .NOT. l_use_sulpc_indirect_sw) THEN
  l_use_seasalt_indirect = .FALSE.
END IF
IF (.NOT. l_nitrate) THEN
  l_use_nitrate_direct = .FALSE.
  l_use_nitrate_indirect = .FALSE.
ELSE IF (.NOT. l_use_sulpc_indirect_lw .OR. .NOT. l_use_sulpc_indirect_sw) THEN
  l_use_nitrate_indirect = .FALSE.
END IF
IF (.NOT. l_soot) THEN
  l_use_soot_direct = .FALSE.
END IF
l_use_soot_indirect = .FALSE.  ! Always false
IF (.NOT. l_biomass) THEN
  l_use_bmass_direct = .FALSE.
  l_use_bmass_indirect = .FALSE.
ELSE IF (.NOT. l_use_sulpc_indirect_lw .OR. .NOT. l_use_sulpc_indirect_sw) THEN
  l_use_bmass_indirect = .FALSE.
END IF
IF (.NOT. l_ocff) THEN
  l_use_ocff_direct = .FALSE.
  l_use_ocff_indirect = .FALSE.
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE run_aerosol_check

SUBROUTINE print_nlist_run_aerosol()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_AEROSOL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_aerosol', &
    src='run_aerosol_mod')

WRITE(lineBuffer,*)' l_sulpc_so2 = ',l_sulpc_so2
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_so2_surfem = ',l_so2_surfem
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_so2_hilem = ',l_so2_hilem
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_so2_natem = ',l_so2_natem
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_so2_lbc = ',l_so2_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_sulpc_dms = ',l_sulpc_dms
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_dms_em = ',l_dms_em
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_dms_em_inter = ',l_dms_em_inter
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' i_dms_flux = ',i_dms_flux
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_sulpc_ozone = ',l_sulpc_ozone
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_sulpc_online_oxidants = ',l_sulpc_online_oxidants
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_sulpc_nh3 = ',l_sulpc_nh3
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_nh3_em = ',l_nh3_em
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_sulpc_so2_o3_nonbuffered=',l_sulpc_so2_o3_nonbuffered
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_soot = ',l_soot
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_soot_surem = ',l_soot_surem
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_soot_hilem = ',l_soot_hilem
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_biomass = ',l_biomass
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_bmass_surem = ',l_bmass_surem
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_bmass_hilem = ',l_bmass_hilem
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_ocff = ',l_ocff
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_ocff_surem = ',l_ocff_surem
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_ocff_hilem = ',l_ocff_hilem
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_nitrate = ',l_nitrate
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_use_sulphate_sulpc = ',l_use_sulphate_sulpc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_use_nitrate_sulpc = ',l_use_nitrate_sulpc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_use_seasalt_sulpc = ',l_use_seasalt_sulpc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_use_seasalt_pm = ',l_use_seasalt_pm
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,'(A18, L7)')' l_temporal_emi = ',l_temporal_emi
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_use_bmass_sulpc = ',l_use_bmass_sulpc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_use_ocff_sulpc = ',l_use_ocff_sulpc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' so2_high_level = ',so2_high_level
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' soot_high_level = ',soot_high_level
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' bMass_high_level_1 = ',bMass_high_level_1
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' bMass_high_level_2 = ',bMass_high_level_2
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' ocff_high_level = ',ocff_high_level
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_so2_lbc = ',l_so2_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_dms_lbc = ',l_dms_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_so4_aitken_lbc = ',l_so4_aitken_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_so4_accu_lbc = ',l_so4_accu_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_so4_diss_lbc = ',l_so4_diss_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_nh3_lbc = ',l_nh3_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_soot_new_lbc = ',l_soot_new_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_soot_agd_lbc = ',l_soot_agd_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_soot_cld_lbc = ',l_soot_cld_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_bmass_new_lbc = ',l_bmass_new_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_bmass_agd_lbc = ',l_bmass_agd_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_bmass_cld_lbc = ',l_bmass_cld_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_ocff_new_lbc = ',l_ocff_new_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_ocff_agd_lbc = ',l_ocff_agd_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_ocff_cld_lbc = ',l_ocff_cld_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_nitr_acc_lbc = ',l_nitr_acc_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,*)' l_nitr_diss_lbc = ',l_nitr_diss_lbc
CALL umPrint(lineBuffer,src='run_aerosol_mod')
WRITE(lineBuffer,'(A26, L1)') ' l_bmass_hilem_variable = ', &
                                l_bmass_hilem_variable
CALL umPrint(lineBuffer,src='run_aerosol_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='run_aerosol_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_aerosol

#if !defined(LFRIC)
SUBROUTINE read_nml_run_aerosol(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_AEROSOL'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 6
INTEGER, PARAMETER :: n_log = 38


TYPE my_namelist
  SEQUENCE
  INTEGER :: SO2_high_level
  INTEGER :: i_dms_flux
  INTEGER :: soot_high_level
  INTEGER :: bmass_high_level_1
  INTEGER :: bmass_high_level_2
  INTEGER :: ocff_high_level
  LOGICAL :: l_sulpc_so2
  LOGICAL :: l_so2_surfem
  LOGICAL :: l_so2_hilem
  LOGICAL :: l_so2_natem
  LOGICAL :: l_sulpc_dms
  LOGICAL :: l_dms_em
  LOGICAL :: l_dms_em_inter
  LOGICAL :: l_sulpc_ozone
  LOGICAL :: l_sulpc_online_oxidants
  LOGICAL :: l_sulpc_2_way_coupling
  LOGICAL :: l_sulpc_so2_o3_nonbuffered
  LOGICAL :: l_sulpc_nh3
  LOGICAL :: l_nh3_em
  LOGICAL :: l_soot
  LOGICAL :: l_soot_surem
  LOGICAL :: l_soot_hilem
  LOGICAL :: l_biomass
  LOGICAL :: l_bmass_surem
  LOGICAL :: l_bmass_hilem
  LOGICAL :: l_ocff
  LOGICAL :: l_ocff_surem
  LOGICAL :: l_ocff_hilem
  LOGICAL :: l_nitrate
  LOGICAL :: l_use_sulphate_sulpc
  LOGICAL :: l_use_nitrate_sulpc
  LOGICAL :: l_use_seasalt_sulpc
  LOGICAL :: l_use_seasalt_pm
  LOGICAL :: l_temporal_emi
  LOGICAL :: l_use_bmass_sulpc
  LOGICAL :: l_use_ocff_sulpc
  LOGICAL :: l_so2_lbc
  LOGICAL :: l_dms_lbc
  LOGICAL :: l_nh3_lbc
  LOGICAL :: l_soot_lbc
  LOGICAL :: l_bmass_lbc
  LOGICAL :: l_ocff_lbc
  LOGICAL :: l_nitr_lbc
  LOGICAL :: l_bmass_hilem_variable
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,    &
                    n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Aerosol, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Aerosol", iomessage)

  my_nml % SO2_high_level     = SO2_high_level
  my_nml % i_dms_flux         = i_dms_flux
  my_nml % soot_high_level    = soot_high_level
  my_nml % bmass_high_level_1 = bmass_high_level_1
  my_nml % bmass_high_level_2 = bmass_high_level_2
  my_nml % ocff_high_level    = ocff_high_level
  ! end of integers
  my_nml % l_sulpc_so2    = l_sulpc_so2
  my_nml % l_so2_surfem   = l_so2_surfem
  my_nml % l_so2_hilem    = l_so2_hilem
  my_nml % l_so2_natem    = l_so2_natem
  my_nml % l_sulpc_dms    = l_sulpc_dms
  my_nml % l_dms_em       = l_dms_em
  my_nml % l_dms_em_inter = l_dms_em_inter
  my_nml % l_sulpc_ozone  = l_sulpc_ozone
  my_nml % l_sulpc_online_oxidants    = l_sulpc_online_oxidants
  my_nml % l_sulpc_2_way_coupling     = l_sulpc_2_way_coupling
  my_nml % l_sulpc_so2_o3_nonbuffered = l_sulpc_so2_o3_nonbuffered
  my_nml % l_sulpc_nh3    = l_sulpc_nh3
  my_nml % l_nh3_em       = l_nh3_em
  my_nml % l_soot         = l_soot
  my_nml % l_soot_surem   = l_soot_surem
  my_nml % l_soot_hilem   = l_soot_hilem
  my_nml % l_biomass      = l_biomass
  my_nml % l_bmass_surem  = l_bmass_surem
  my_nml % l_bmass_hilem  = l_bmass_hilem
  my_nml % l_ocff         = l_ocff
  my_nml % l_ocff_surem   = l_ocff_surem
  my_nml % l_ocff_hilem   = l_ocff_hilem
  my_nml % l_nitrate      = l_nitrate
  my_nml % l_use_sulphate_sulpc = l_use_sulphate_sulpc
  my_nml % l_use_nitrate_sulpc  = l_use_nitrate_sulpc
  my_nml % l_use_seasalt_sulpc  = l_use_seasalt_sulpc
  my_nml % l_use_seasalt_pm     = l_use_seasalt_pm
  my_nml % l_temporal_emi        = l_temporal_emi
  my_nml % l_use_bmass_sulpc    = l_use_bmass_sulpc
  my_nml % l_use_ocff_sulpc     = l_use_ocff_sulpc
  my_nml % l_so2_lbc      = l_so2_lbc
  my_nml % l_dms_lbc      = l_dms_lbc
  my_nml % l_nh3_lbc      = l_nh3_lbc
  my_nml % l_soot_lbc     = l_soot_lbc
  my_nml % l_bmass_lbc    = l_bmass_lbc
  my_nml % l_ocff_lbc     = l_ocff_lbc
  my_nml % l_nitr_lbc     = l_nitr_lbc
  my_nml % l_bmass_hilem_variable  = l_bmass_hilem_variable

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  SO2_high_level     = my_nml % SO2_high_level
  i_dms_flux         = my_nml % i_dms_flux
  soot_high_level    = my_nml % soot_high_level
  bmass_high_level_1 = my_nml % bmass_high_level_1
  bmass_high_level_2 = my_nml % bmass_high_level_2
  ocff_high_level    = my_nml % ocff_high_level
  ! end of integers
  l_sulpc_so2    = my_nml % l_sulpc_so2
  l_so2_surfem   = my_nml % l_so2_surfem
  l_so2_hilem    = my_nml % l_so2_hilem
  l_so2_natem    = my_nml % l_so2_natem
  l_sulpc_dms    = my_nml % l_sulpc_dms
  l_dms_em       = my_nml % l_dms_em
  l_dms_em_inter = my_nml % l_dms_em_inter
  l_sulpc_ozone  = my_nml % l_sulpc_ozone
  l_sulpc_online_oxidants    = my_nml % l_sulpc_online_oxidants
  l_sulpc_2_way_coupling     = my_nml %  l_sulpc_2_way_coupling
  l_sulpc_so2_o3_nonbuffered = my_nml % l_sulpc_so2_o3_nonbuffered
  l_sulpc_nh3    = my_nml % l_sulpc_nh3
  l_nh3_em       = my_nml % l_nh3_em
  l_soot         = my_nml % l_soot
  l_soot_surem   = my_nml % l_soot_surem
  l_soot_hilem   = my_nml % l_soot_hilem
  l_biomass      = my_nml % l_biomass
  l_bmass_surem  = my_nml % l_bmass_surem
  l_bmass_hilem  = my_nml % l_bmass_hilem
  l_ocff         = my_nml % l_ocff
  l_ocff_surem   = my_nml % l_ocff_surem
  l_ocff_hilem   = my_nml % l_ocff_hilem
  l_nitrate      = my_nml % l_nitrate
  l_use_sulphate_sulpc = my_nml % l_use_sulphate_sulpc
  l_use_nitrate_sulpc  = my_nml % l_use_nitrate_sulpc
  l_use_seasalt_sulpc  = my_nml % l_use_seasalt_sulpc
  l_use_seasalt_pm     = my_nml % l_use_seasalt_pm
  l_temporal_emi       = my_nml % l_temporal_emi
  l_use_bmass_sulpc    = my_nml % l_use_bmass_sulpc
  l_use_ocff_sulpc     = my_nml % l_use_ocff_sulpc
  l_so2_lbc      = my_nml % l_so2_lbc
  l_dms_lbc      = my_nml % l_dms_lbc
  l_nh3_lbc      = my_nml % l_nh3_lbc
  l_soot_lbc     = my_nml % l_soot_lbc
  l_bmass_lbc    = my_nml % l_bmass_lbc
  l_ocff_lbc     = my_nml % l_ocff_lbc
  l_nitr_lbc     = my_nml % l_nitr_lbc
  l_bmass_hilem_variable =  my_nml % l_bmass_hilem_variable

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_aerosol
#endif

END MODULE run_aerosol_mod

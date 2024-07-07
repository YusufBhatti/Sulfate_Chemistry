! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  List of stashcode magic numbers

MODULE pws_vertlev_choices_mod

! Description:
!   Section 20 diags have some hard wired choices on what levels to employ
!   on a diag by diag basis, based upon the number of levels in the model.
!   In future the model should be set up to derive these rather than have
!   a magic number list.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS diagnostics
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

! List of model levels that Fieldcalc diags will work for
INTEGER, PARAMETER :: Model_Levels_38 = 38
INTEGER, PARAMETER :: Model_Levels_50 = 50
INTEGER, PARAMETER :: Model_Levels_80 = 80
!!! Add this for AQUM 63 level 40km top
INTEGER, PARAMETER :: Model_Levels_63 = 63
INTEGER, PARAMETER :: Model_Levels_70 = 70
INTEGER, PARAMETER :: Model_Levels_85 = 85
REAL, PARAMETER :: FortyKM         = 40000.0
REAL, PARAMETER :: SixtythreeKM    = 63000.0
REAL, PARAMETER :: EightyKM        = 80000.0
REAL, PARAMETER :: EightyfiveKM    = 85000.0
REAL, PARAMETER :: ThirtyEightKM   = 38500.0

! Variables for model level ranges required to derive diagnostics.

INTEGER :: TP_NumLevs_Gl = imdi   ! TropHeight Fields:
INTEGER :: TP_ZeroLev_Gl = imdi   !
INTEGER :: IT_NumLevs_Gl = imdi   ! Isotherm Fields:
INTEGER :: IT_ZeroLev_Gl = imdi   !
INTEGER :: CT_NumLevs_Gl = imdi   ! Contrail Fields:
INTEGER :: CT_ZeroLev_Gl = imdi   !
INTEGER :: MX_NumLevs_Gl = imdi   ! MaxWind Fields:
INTEGER :: MX_ZeroLev_Gl = imdi   !
INTEGER :: CA_NumLevs_Gl = imdi   ! CAT Fields:
INTEGER :: CA_ZeroLev_Gl = imdi   !
INTEGER :: MW_NumLevs_Gl = imdi   ! MtnStress Fields:
INTEGER :: MW_ZeroLev_Gl = imdi   !

INTEGER :: WT_NumLevs_Gl = imdi   ! WAFC CAT turb:
INTEGER :: WT_ZeroLev_Gl = imdi   !
INTEGER :: IC_NumLevs_Gl = imdi   ! Icing Fields (original alg):
INTEGER :: IC_ZeroLev_Gl = imdi   !
INTEGER :: LI_NumLevs_Gl = imdi   ! Icing Fields:
INTEGER :: LI_ZeroLev_Gl = imdi   !
INTEGER :: ICT_NumLevs_Gl = imdi  ! In (Layer) cloud turb Fields:
INTEGER :: ICT_ZeroLev_Gl = imdi  !
INTEGER :: CB_NumLevs_Gl = imdi   ! CB Fields:
INTEGER :: CB_ZeroLev_Gl = imdi   !
INTEGER :: DT_NumLevs_Gl = imdi   ! Dust concs
INTEGER :: DT_ZeroLev_Gl = imdi   ! 
 
CONTAINS



SUBROUTINE pws_vertlev_choices_init(model_lid)

USE nlsizes_namelist_mod, ONLY: model_levels
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

REAL, INTENT(IN) :: model_lid

INTEGER :: errorstatus

  !---------------------------------------------------------------
  ! 1.4.1 Initialise model level ranges now that Model_Levels is known
  !---------------------------------------------------------------

  SELECT CASE ( model_levels )

  CASE ( Model_Levels_38 )

    IF (model_lid /= FortyKM) THEN
      ErrorStatus = -30
      CALL ereport('pws_vertlev_choices_init', ErrorStatus, & !error- no output
                 "Section 20 diags may be incorrect as many expect "//&
                 "model lid to be at 40km for 38 level model sets." )
    END IF

    TP_NumLevs_Gl = 27    ! TropHeight Fields:
    TP_ZeroLev_Gl = 5     !   Model Levels 6-32
    IT_NumLevs_Gl = 32    ! Isotherm Fields:
    IT_ZeroLev_Gl = 0     !   Model Levels 1-32
    CT_NumLevs_Gl = 30    ! Contrail Fields:
    CT_ZeroLev_Gl = 0     !   Model Levels 1-30
    MX_NumLevs_Gl = 32    ! MaxWind Fields:
    MX_ZeroLev_Gl = 0     !   Model Levels 1-32
    CA_NumLevs_Gl = 24    ! CAT Fields:
    CA_ZeroLev_Gl = 5     !   Model Levels 6-29
    MW_NumLevs_Gl = 15    ! MtnStress Fields:
    MW_ZeroLev_Gl = 13    !   Model Levels 14-28

    WT_NumLevs_Gl = 32    ! WAFC CAT turb ...
    WT_ZeroLev_Gl = 0     !   Model Levels 1-32
    IC_NumLevs_Gl = 24    ! Icing Fields (original alg):
    IC_ZeroLev_Gl = 0     !   Model Levels 1-24
    LI_NumLevs_Gl = 33    ! Icing Fields:
    LI_ZeroLev_Gl = 0     !   Model Levels 1-33
    ICT_NumLevs_Gl= 24    ! In (Layer) cloud turb Fields:
    ICT_ZeroLev_Gl= 0     !   Model Levels 1-24
    CB_NumLevs_Gl = 32    ! CB Fields:
    CB_ZeroLev_Gl = 0     !   Model Levels 1-32

    DT_NumLevs_Gl = 10    ! Dust concs
    DT_ZeroLev_Gl = 0     !   Model Levels 1-10


  CASE ( Model_Levels_50 )

    IF (model_lid /= SixtythreeKM) THEN
      ErrorStatus = -30
      CALL ereport('pws_vertlev_choices_init', ErrorStatus, & !error- no output
                 "Section 20 diags may be incorrect as many expect "//&
                 "model lid to be at 63km for 50 level model sets." )
    END IF

    TP_NumLevs_Gl = 28    ! TropHeight Fields:
    TP_ZeroLev_Gl = 5     !   Model Levels 6-33
    IT_NumLevs_Gl = 33    ! Isotherm Fields:
    IT_ZeroLev_Gl = 0     !   Model Levels 1-33
    CT_NumLevs_Gl = 30    ! Contrail Fields:
    CT_ZeroLev_Gl = 0     !   Model Levels 1-30
    MX_NumLevs_Gl = 33    ! MaxWind Fields:
    MX_ZeroLev_Gl = 0     !   Model Levels 1-33
    CA_NumLevs_Gl = 24    ! CAT Fields:
    CA_ZeroLev_Gl = 5     !   Model Levels 6-29
    MW_NumLevs_Gl = 15    ! MtnStress Fields:
    MW_ZeroLev_Gl = 13    !   Model Levels 14-28

    WT_NumLevs_Gl = 33    ! WAFC CAT turb ...
    WT_ZeroLev_Gl = 0     !   Model Levels 1-33
    IC_NumLevs_Gl = 24    ! Icing Fields (original alg):
    IC_ZeroLev_Gl = 0     !   Model Levels 1-24
    LI_NumLevs_Gl = 33    ! Icing Fields:
    LI_ZeroLev_Gl = 0     !   Model Levels 1-33
    ICT_NumLevs_Gl= 24    ! In (Layer) cloud turb Fields:
    ICT_ZeroLev_Gl= 0     !   Model Levels 1-24
    CB_NumLevs_Gl = 33    ! CB Fields:
    CB_ZeroLev_Gl = 0     !   Model Levels 1-33

    DT_NumLevs_Gl = 10    ! Dust concs
    DT_ZeroLev_Gl = 0     !   Model Levels 1-10


  CASE ( Model_Levels_70 )

    IF (model_lid /= EightyKM) THEN
      ErrorStatus = -30
      CALL ereport('pws_vertlev_choices_init', ErrorStatus, & !error- no output
                 "Section 20 diags may be incorrect as many expect "//&
                 "model lid to be at 80km for 70 level model sets." )
    END IF

    TP_NumLevs_Gl = 46    ! TropHeight Fields:
    TP_ZeroLev_Gl = 8     !   Model Levels 9-54
    IT_NumLevs_Gl = 54    ! Isotherm Fields:
    IT_ZeroLev_Gl = 0     !   Model Levels 1-54
    CT_NumLevs_Gl = 50    ! Contrail Fields:
    CT_ZeroLev_Gl = 0     !   Model Levels 1-50
    MX_NumLevs_Gl = 54    ! MaxWind Fields:
    MX_ZeroLev_Gl = 0     !   Model Levels 1-54
    CA_NumLevs_Gl = 41    ! CAT Fields:
    CA_ZeroLev_Gl = 8     !   Model Levels 9-49
    MW_NumLevs_Gl = 25    ! MtnStress Fields:
    MW_ZeroLev_Gl = 22    !   Model Levels 23-47

    WT_NumLevs_Gl = 54    ! WAFC CAT turb ...
    WT_ZeroLev_Gl = 0     !   Model Levels 1-54
    IC_NumLevs_Gl = 41    ! Icing Fields:
    IC_ZeroLev_Gl = 0     !   Model Levels 1-41
    LI_NumLevs_Gl = 54    ! Icing Fields:
    LI_ZeroLev_Gl = 0     !   Model Levels 1-54
    ICT_NumLevs_Gl= 41    ! In (Layer) cloud turb Fields:
    ICT_ZeroLev_Gl= 0     !   Model Levels 1-41
    CB_NumLevs_Gl = 54    ! CB Fields:
    CB_ZeroLev_Gl = 0     !   Model Levels 1-54

    DT_NumLevs_Gl = 16    ! Dust concs
    DT_ZeroLev_Gl = 0     !   Model Levels 1-16

  CASE ( Model_Levels_63 )

    IF (model_lid /= FortyKM) THEN
      ErrorStatus = -30
      CALL ereport('pws_vertlev_choices_init', ErrorStatus, & !error- no output
                 "Section 20 diags may be incorrect as many expect "//&
                 "model lid to be at 40km for 63 level model sets." )
    END IF

    TP_NumLevs_Gl = 46    ! TropHeight Fields:
    TP_ZeroLev_Gl = 8     !   Model Levels 9-54
    IT_NumLevs_Gl = 54    ! Isotherm Fields:
    IT_ZeroLev_Gl = 0     !   Model Levels 1-54
    CT_NumLevs_Gl = 50    ! Contrail Fields:
    CT_ZeroLev_Gl = 0     !   Model Levels 1-50
    MX_NumLevs_Gl = 54    ! MaxWind Fields:
    MX_ZeroLev_Gl = 0     !   Model Levels 1-54
    CA_NumLevs_Gl = 41    ! CAT Fields:
    CA_ZeroLev_Gl = 8     !   Model Levels 9-49
    MW_NumLevs_Gl = 25    ! MtnStress Fields:
    MW_ZeroLev_Gl = 22    !   Model Levels 23-47

    WT_NumLevs_Gl = 54    ! WAFC CAT turb ...
    WT_ZeroLev_Gl = 0     !   Model Levels 1-54
    IC_NumLevs_Gl = 41    ! Icing Fields:
    IC_ZeroLev_Gl = 0     !   Model Levels 1-41
    LI_NumLevs_Gl = 54    ! Icing Fields:
    LI_ZeroLev_Gl = 0     !   Model Levels 1-54
    ICT_NumLevs_Gl= 41    ! In (Layer) cloud turb Fields:
    ICT_ZeroLev_Gl= 0     !   Model Levels 1-41
    CB_NumLevs_Gl = 54    ! CB Fields:
    CB_ZeroLev_Gl = 0     !   Model Levels 1-54

    DT_NumLevs_Gl = 16    ! Dust concs
    DT_ZeroLev_Gl = 0     !   Model Levels 1-16

  CASE ( Model_Levels_80 )
! This case was added solely to output the updraft helicity diagnostic for the
! SingV/RA1T configurations. As such, only MX_NumLevs_Gl has been set up 
! correctly and the other levels in the fields below may be incorrect. 

    IF (model_lid /= ThirtyEightKM) THEN
      ErrorStatus = -30
      CALL EReport('pws_vertlev_choices_init',ErrorStatus, & !error - no output
                 "Section 20 diags may be incorrect as many expect "//&
                 "model lid to be at 38.5km for 80 level model sets." )
    END IF

    TP_NumLevs_Gl = 46    ! TropHeight Fields:
    TP_ZeroLev_Gl = 8     !   Model Levels 9-54
    IT_NumLevs_Gl = 54    ! Isotherm Fields:
    IT_ZeroLev_Gl = 0     !   Model Levels 1-54
    CT_NumLevs_Gl = 50    ! Contrail Fields:
    CT_ZeroLev_Gl = 0     !   Model Levels 1-50
    MX_NumLevs_Gl = 63    ! MaxWind Fields:
    MX_ZeroLev_Gl = 0     !   Model Levels 1-63
    CA_NumLevs_Gl = 41    ! CAT Fields:
    CA_ZeroLev_Gl = 8     !   Model Levels 9-49
    MW_NumLevs_Gl = 25    ! MtnStress Fields:
    MW_ZeroLev_Gl = 22    !   Model Levels 23-47

    WT_NumLevs_Gl = 54    ! WAFC CAT turb ...
    WT_ZeroLev_Gl = 0     !   Model Levels 1-54
    IC_NumLevs_Gl = 41    ! Icing Fields:
    IC_ZeroLev_Gl = 0     !   Model Levels 1-41
    LI_NumLevs_Gl = 54    ! Icing Fields:
    LI_ZeroLev_Gl = 0     !   Model Levels 1-54
    ICT_NumLevs_Gl= 41    ! In (Layer) cloud turb Fields:
    ICT_ZeroLev_Gl= 0     !   Model Levels 1-41
    CB_NumLevs_Gl = 54    ! CB Fields:
    CB_ZeroLev_Gl = 0     !   Model Levels 1-54

    DT_NumLevs_Gl = 16    ! Dust concs
    DT_ZeroLev_Gl = 0     !   Model Levels 1-16

  CASE ( Model_Levels_85 )

    IF (model_lid /= EightyfiveKM) THEN
      ErrorStatus = -30
      CALL EReport('pws_vertlev_choices_init', ErrorStatus, & !error- no output
                 "Section 20 diags may be incorrect as many expect "//&
                 "model lid to be at 85km for 85 level model sets." )
    END IF

    TP_NumLevs_Gl = 46    ! TropHeight Fields:
    TP_ZeroLev_Gl = 8     !   Model Levels 9-54
    IT_NumLevs_Gl = 54    ! Isotherm Fields:
    IT_ZeroLev_Gl = 0     !   Model Levels 1-54
    CT_NumLevs_Gl = 50    ! Contrail Fields:
    CT_ZeroLev_Gl = 0     !   Model Levels 1-50
    MX_NumLevs_Gl = 54    ! MaxWind Fields:
    MX_ZeroLev_Gl = 0     !   Model Levels 1-54
    CA_NumLevs_Gl = 41    ! CAT Fields:
    CA_ZeroLev_Gl = 8     !   Model Levels 9-49
    MW_NumLevs_Gl = 25    ! MtnStress Fields:
    MW_ZeroLev_Gl = 22    !   Model Levels 23-47
    WT_NumLevs_Gl = 54    ! WAFC CAT turb ...

    WT_ZeroLev_Gl = 0     !   Model Levels 1-54
    IC_NumLevs_Gl = 41    ! Icing Fields:
    IC_ZeroLev_Gl = 0     !   Model Levels 1-41
    LI_NumLevs_Gl = 54    ! Icing Fields:
    LI_ZeroLev_Gl = 0     !   Model Levels 1-54
    ICT_NumLevs_Gl= 41    ! In (Layer) cloud turb Fields:
    ICT_ZeroLev_Gl= 0     !   Model Levels 1-41
    CB_NumLevs_Gl = 54    ! CB Fields:
    CB_ZeroLev_Gl = 0     !   Model Levels 1-54

    DT_NumLevs_Gl = 16    ! Dust concs
    DT_ZeroLev_Gl = 0     !   Model Levels 1-16


  CASE DEFAULT

     ErrorStatus = 33
    CALL ereport('pws_vertlev_choices_init', ErrorStatus, & ! error - no output
                 "Number of levels not catered "//&
                 "for but trying to continue." )
  END SELECT


END  SUBROUTINE pws_vertlev_choices_init

END MODULE pws_vertlev_choices_mod


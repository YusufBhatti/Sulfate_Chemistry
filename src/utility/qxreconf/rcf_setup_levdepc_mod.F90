! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  sets up the output dump level dependent constants

MODULE Rcf_Setup_LevDepC_Mod

!  Subroutine Rcf_Setup_LevDepC - sets up the output dump
!                                 level dependent contants.
!
! Description:
!   The level dependent constants for the output dump are constructed.
!
! Method:
!   The Level dependent constants are constructed from the output
!   grid details - see UMDP F3 for full  details.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SETUP_LEVDEPC_MOD'

CONTAINS

SUBROUTINE Rcf_Setup_LevDepC( Hdr_Out, Hdr_In, Grid )

USE Rcf_Grid_Type_Mod, ONLY: &
    Grid_Type

USE Rcf_UMhead_Mod, ONLY: &
    UM_Header_Type

USE Rcf_Generate_Heights_Mod, ONLY: &
    height_gen_smooth,               &
    height_gen_linear

USE Rcf_Headaddress_Mod, ONLY:            &
    LDC_EtaTheta,            LDC_EtaRho,   &
    LDC_RHCrit,              SoilDepths,   &
    LDC_ZseaTheta,           LDC_CkTheta,  &
    LDC_ZseaRho,             LDC_CkRho

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( UM_Header_Type ), INTENT(IN)    :: Hdr_In
TYPE( UM_Header_Type ), INTENT(INOUT) :: Hdr_Out
TYPE( Grid_Type ), INTENT(IN)         :: Grid      ! Output Grid

! Local Variables
INTEGER                               :: i        ! Looper
INTEGER                               :: j        ! Looper
INTEGER                               :: Len2Min  ! Minimum length

CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SETUP_LEVDEPC'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Eta Theta Values
DO i = 0, Grid % model_levels
  Hdr_Out % LevDepC(i+1,LDC_EtaTheta) = Grid % eta_theta_levels( i )
END DO

! Eta Rho Values
DO i = 1, Grid % model_levels
  Hdr_Out % LevDepC(i,LDC_EtaRho) = Grid % eta_rho_levels( i )
END DO

! RHCrit values
DO i = 1, Grid % model_levels
  Hdr_Out % LevDepC(i,LDC_RHCrit) = Grid % rhcrit(i)
END DO

! Soil Depths
DO i = 1, Grid % sm_levels
  Hdr_Out % LevDepC(i,SoilDepths) = Grid % soil_depths( i )
END DO

! Only need level dependent constants 5-8 if we are aren't
! using the original height generation method (or converting from ECMWF)
IF ( Grid % Height_Gen_Method == height_gen_smooth ) THEN

  ! Zsea values for Theta
  DO i = 0, Grid % model_levels
    Hdr_Out % LevDepC(i+1,LDC_ZseaTheta) = Grid % eta_theta_levels(i)&
                                         * Grid % z_top_of_model
  END DO

  ! Ck values for Theta
  DO i = 0, grid % first_constant_r_rho_level - 1
    Hdr_Out % LevDepC(i+1,LDC_CkTheta) =(1.0 -                       &
            grid % eta_theta_levels( i )/                            &
            grid % eta_rho_levels(grid % first_constant_r_rho_level))&
             ** 2
  END DO

  DO i = grid % first_constant_r_rho_level , grid % model_levels
    Hdr_Out % LevDepC(i+1,LDC_CkTheta) = 0.0
  END DO

  ! Zsea values for Rho
  DO i = 1, Grid % model_levels
    Hdr_Out % LevDepC(i,LDC_ZseaRho) = Grid % eta_rho_levels( i ) *  &
                                       Grid % z_top_of_model
  END DO

  ! Ck values for Rho
  DO i = 1, grid % first_constant_r_rho_level
    Hdr_Out % LevDepC(i,LDC_CkRho) = (1.0 -                          &
            grid % eta_rho_levels( i ) /                             &
            grid % eta_rho_levels(grid % first_constant_r_rho_level))&
            ** 2
  END DO

  DO i = grid % first_constant_r_rho_level + 1, grid % model_levels
    Hdr_Out % LevDepC(i,LDC_CkRho) = 0.0
  END DO

ELSE IF ( Grid % Height_Gen_Method == height_gen_linear ) THEN

  ! Zsea values for Theta
  DO i = 0, grid % model_levels
    Hdr_Out % LevDepC(i+1,LDC_ZseaTheta) = Grid % eta_theta_levels(i) &
                                         * Grid % z_top_of_model
  END DO

  ! Ck values for Theta
  DO i = 0, grid % model_levels
    Hdr_Out % LevDepC(i+1,LDC_CkTheta) = 1.0 - grid % eta_theta_levels(i)
  END DO

  ! Zsea values for Rho
  DO i = 1, grid % model_levels
    Hdr_Out % LevDepC(i,LDC_ZseaRho) = Grid % eta_rho_levels(i) *    &
                                       Grid % z_top_of_model
  END DO

  ! Ck values for Rho
  DO i = 1, grid % model_levels
    Hdr_Out % LevDepC(i,LDC_CkRho) = 1.0 - grid % eta_rho_levels(i)
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Setup_LevDepC

END MODULE Rcf_Setup_LevDepC_Mod

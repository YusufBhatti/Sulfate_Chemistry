! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_assign_vars_recon_vertical_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                  &
                          ModuleName='RCF_ASSIGN_VARS_RECON_VERTICAL_MOD'

! Subroutine rcf_assign_vars_recon_vertical
!
! Description:
!   Having read the vertical namelists, allocating space and assigning values
!   to the relevant bits of the Output Grid defined type
!
! Method:
!   Use the values read from the namelists directly from the modules
!   they're in and allocate space within the output_grid defined type
!   to store those values.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!

CONTAINS

SUBROUTINE rcf_assign_vars_recon_vertical ( )

USE Rcf_Grid_Type_Mod, ONLY: Output_Grid

USE Rcf_Generate_Heights_Mod, ONLY:                                      &
                             height_gen_smooth,                          &
                             height_gen_original

USE jules_soil_mod, ONLY:    dzsoil_io

USE Ereport_Mod, ONLY:       Ereport

USE vertnamelist_mod, ONLY:                                              &
                             first_constant_r_rho_level,                 &
                             z_top_of_model,                             &
                             eta_theta,                                  &
                             eta_rho

USE Atmos_Max_Sizes, ONLY:   model_levels_max

USE errormessagelength_mod, ONLY:                                        &
                             errormessagelength

USE yomhook,   ONLY:         lhook,                                      &
                             dr_hook

USE parkind1,  ONLY:         jprb,                                       &
                             jpim

IMPLICIT NONE

! default parameter, which has become fixed, but can be changed if required.
INTEGER, PARAMETER                 :: height_method = height_gen_smooth
INTEGER                            :: ErrorStatus

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (LEN=errormessagelength) :: Cmessage
CHARACTER (LEN=*), PARAMETER       ::                                    &
                           RoutineName = 'RCF_ASSIGN_VARS_RECON_VERTICAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ---------------------------------------------------
! Some sanity checks before we allocate space
! ---------------------------------------------------
IF ( height_method /= height_gen_original .AND.   &
     height_method /= height_gen_smooth ) THEN
  ErrorStatus = 90
  Cmessage = 'Unrecognised height generation method'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

IF ( Output_Grid % model_levels < 1 ) THEN
  Cmessage = 'Model levels for output grid is < 1!'
  ErrorStatus = 20
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
ELSE IF ( Output_Grid % model_levels > model_levels_max ) THEN
  WRITE(Cmessage,'(A,I0,A)' ) 'Model Levels for output grid is > ' //    &
       ' model_levels_max (', model_levels_max, ') - Is this correct?'
  ErrorStatus = 30
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! ---------------------------------------------------
! Allocate space and fill up relevant parts of Output_Grid
! ---------------------------------------------------
ALLOCATE (Output_Grid % eta_theta_levels                                 &
                       ( 0 : Output_Grid % model_levels), STAT=ErrorStatus )
IF ( ErrorStatus /= 0 ) THEN
  Cmessage = 'Unable to allocate memory for Output_Grid % eta_theta_levels'
  ErrorStatus = 50
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF
ALLOCATE ( Output_Grid % eta_rho_levels( Output_Grid % model_levels ),   &
           STAT=ErrorStatus )
IF ( ErrorStatus /= 0 ) THEN
  Cmessage = 'Unable to allocate memory for Output_Grid % eta_rho_levels'
  ErrorStatus = 60
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF
ALLOCATE ( Output_Grid % soil_depths( Output_Grid % sm_levels ),         &
           STAT=ErrorStatus )
IF ( ErrorStatus /= 0 ) THEN
  Cmessage = 'Unable to allocate memory for Output_Grid % soil_depths'
  ErrorStatus = 70
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! Readjusting EG levels from the 1:n+1 namelists that are read in on
! to 0:n EG expects
Output_Grid % eta_theta_levels( 0 : Output_Grid % model_levels ) =       &
                     eta_theta( 1 : Output_Grid % model_levels + 1)
Output_Grid % eta_rho_levels( 1 : Output_Grid % model_levels ) =         &
                     eta_rho( 1 : Output_Grid % model_levels )
Output_Grid % soil_depths( 1 : Output_Grid % sm_levels ) =               &
                dzsoil_io( 1 : Output_Grid % sm_levels )
Output_Grid % height_gen_method = height_method
Output_Grid % z_top_of_model = z_top_of_model
Output_Grid % first_constant_r_rho_level = first_constant_r_rho_level

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE  rcf_assign_vars_recon_vertical

END MODULE rcf_assign_vars_recon_vertical_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Sets the number of land points when specified manually by the items namelist.
!
MODULE Rcf_Set_LSM_Land_Points_mod

IMPLICIT NONE

! Description:
!   If the items namelist is being used to set the number of land points
!   by some fixed method (i.e. not read from ancillary or taken from the input
!   dump) then this routine calculates the number of land points expected.
!   Note: this is not where the land sea mask is read/calculated for 
!   putting into the output dump - this is done later in Rcf_Setup_LSM_Out.
!
! Method:
!   Check which method is being used to set the LSM and set the number of land
!   points to all or none as appropriate.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SET_LSM_LAND_POINTS_MOD'

CONTAINS

! Subroutine interface:
SUBROUTINE Rcf_Set_LSM_Land_Points(Output_Grid)

USE Ereport_Mod, ONLY:                 &
    Ereport

USE errormessagelength_mod, ONLY:      &
    errormessagelength

USE nlsizes_namelist_mod, ONLY:        &
    land_field

USE Rcf_Grid_Type_Mod, ONLY:           &
    grid_type

USE Rcf_Lsm_Mod, ONLY:                 &
    glob_land_out,                     &
    lsm_fixed_value,                   &
    lsm_source

USE items_nml_mod, ONLY :              &
    set_to_const,                      &
    set_to_zero

USE umPrintMgr, ONLY:                  &
    newline

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (grid_type), TARGET, INTENT(INOUT)   :: Output_Grid

! Local Data
INTEGER                      :: glob_land_points    ! land_points on all PEs
INTEGER                      :: errorstatus
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_SET_LSM_LAND_POINTS'
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (lsm_source == set_to_zero) THEN

  ! All sea
  glob_land_points = 0

ELSE IF (lsm_source == set_to_const) THEN

  IF (ABS(lsm_fixed_value - 0.0) < EPSILON(1.0)) THEN

    ! All sea
    glob_land_points = 0

  ELSE IF (ABS(lsm_fixed_value - 1.0) < EPSILON(1.0)) THEN

    ! All land
    glob_land_points = output_grid % glob_p_field

  ELSE

    errorstatus = 10
    WRITE(cmessage,'(A,F16.4)')                                       &
      'Error setting constant value for land sea mask.' // newline // &
      'Only 0.0 (sea) or 1.0 (land) are permitted. Input value is: ', &
      lsm_fixed_value
   CALL ereport(RoutineName, errorstatus, cmessage)

  END IF

END IF

! Override current values in output grid and namelist:
land_field = glob_land_points
glob_land_out = glob_land_points
Output_Grid % glob_land_field = glob_land_points

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Set_LSM_Land_Points
END MODULE Rcf_Set_LSM_Land_Points_mod

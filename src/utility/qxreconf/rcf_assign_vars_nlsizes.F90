! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE rcf_assign_vars_nlsizes_mod

IMPLICIT NONE

! Description:
!  Assigns variables to the various grid levels/sizes to be used,
!  based on the contents of the nlsizes namelist
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!  Language: FORTRAN 2003
!  This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_ASSIGN_VARS_NLSIZES_MOD'

CONTAINS

SUBROUTINE rcf_assign_vars_nlsizes()

USE nlsizes_namelist_mod, ONLY:                                          &
    st_levels,                                                           &
    sm_levels,                                                           &
    global_rows,                                                         &
    global_row_length,                                                   &
    river_rows,                                                          &
    river_row_length,                                                    &
    model_levels,                                                        &
    ozone_levels,                                                        &
    bl_levels,                                                           &
    land_field,                                                          &
    tpps_ozone_levels,                                                   &
    tr_ukca,                                                             &
    a_len2_rowdepc,                                                      &
    a_len2_coldepc ,                                                     &
    a_len_inthd,                                                         &
    a_len_realhd,                                                        &
    a_len2_levdepc,                                                      &
    a_len2_flddepc,                                                      &
    a_len_extcnst                                                             

USE lbc_mod, ONLY:                                                       &
    rimwidtha

USE rimtypes, ONLY:                                                      &
    rima_type_norm

USE rcf_lsm_mod, ONLY:                                                   &
    glob_land_out

USE rcf_grid_type_mod, ONLY:                                             &
    output_grid

USE missing_data_mod, ONLY:                                              &
    imdi

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_ASSIGN_VARS_NLSIZES'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
tr_ukca           = 0
river_rows        = 180   ! 1 degree default
river_row_length  = 360   ! 1 degree default
rimwidtha(rima_type_norm) = 0
a_len2_rowdepc    = imdi
a_len2_coldepc    = imdi
a_len_inthd       = imdi
a_len_realhd      = imdi
a_len2_levdepc    = imdi
a_len2_flddepc    = imdi
a_len_extcnst     = imdi

! Set the number of soil moisture levels to be the
! same as the number of soil temperature levels
! Always true for 8A hydrology which is the only
! option
sm_levels=st_levels
output_grid % sm_levels = sm_levels

! Set output variables as required
output_grid % glob_p_rows       = global_rows
output_grid % glob_p_row_length = global_row_length
output_grid % glob_r_rows       = river_rows
output_grid % glob_r_row_length = river_row_length
output_grid % model_levels      = model_levels
output_grid % ozone_levels      = ozone_levels
output_grid % st_levels         = st_levels
output_grid % bl_levels         = bl_levels

! Default value - not required in the dump
output_grid % cloud_levels      = imdi

! Set some land/sea maks values
glob_land_out                   = land_field
output_grid % glob_land_field   = land_field

! Set values for tropopause-based ozone.
! In future developments, obtain directly from namelist.
tpps_ozone_levels               = ozone_levels

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_assign_vars_nlsizes

END MODULE rcf_assign_vars_nlsizes_mod

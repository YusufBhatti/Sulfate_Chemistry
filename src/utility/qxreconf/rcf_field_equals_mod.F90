! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Copy one field into another - *IGNORING THE DATA*

MODULE Rcf_Field_Equals_Mod

!  Subroutine Rcf_Field_Equals - copys a field description into another
!
! Description:
!   Field size descriptions etc are copied from one field to another.
!   No data is copied or data pointers set.
!
! Method:
!   Stashmaster pointers are pointed to!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

INTERFACE rcf_field_equals
MODULE PROCEDURE rcf_field_equals_field, rcf_field_equals_grid
END INTERFACE

PRIVATE :: lhook, dr_hook, jpim, jprb

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_FIELD_EQUALS_MOD'

CONTAINS

SUBROUTINE Rcf_Field_Equals_Field( field_out, field_in )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

IMPLICIT NONE

! Arguments
TYPE (field_type), INTENT(IN) :: field_in
TYPE (field_type), INTENT(OUT) :: field_out
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_FIELD_EQUALS_FIELD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
field_out % levels          = field_in % levels
field_out % rows            = field_in % rows
field_out % row_len         = field_in % row_len
field_out % level_size      = field_in % level_size
field_out % bottom_level    = field_in % bottom_level
field_out % top_level       = field_in % top_level
field_out % glob_rows       = field_in % glob_rows
field_out % glob_row_len    = field_in % glob_row_len
field_out % glob_level_size = field_in % glob_level_size
field_out % dump_pos        = field_in % dump_pos
field_out % interp          = field_in % interp
field_out % stashmaster    => field_in % stashmaster

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Field_Equals_Field

SUBROUTINE Rcf_Field_Equals_Grid( field_out, grid_in )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_no_op

USE missing_data_mod, ONLY: &
    imdi

IMPLICIT NONE

! Arguments
TYPE (grid_type),  INTENT(IN) :: grid_in
TYPE (field_type), INTENT(OUT) :: field_out
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_FIELD_EQUALS_GRID'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Assume 1 level for now.  Can be changed outside this call.
field_out % levels          = 1
field_out % bottom_level    = -1
field_out % top_level       = -1
field_out % rows            = grid_in % loc_p_rows
field_out % row_len         = grid_in % loc_p_row_length
field_out % level_size      = grid_in % loc_p_field
field_out % glob_rows       = grid_in % glob_p_rows
field_out % glob_row_len    = grid_in % glob_p_row_length
field_out % glob_level_size = grid_in % glob_p_field
field_out % dump_pos        = imdi
field_out % interp          = interp_no_op
NULLIFY(field_out % stashmaster)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Field_Equals_Grid

END MODULE Rcf_Field_Equals_Mod

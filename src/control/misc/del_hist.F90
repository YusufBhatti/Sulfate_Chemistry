! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE del_hist_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DEL_HIST_MOD'

CONTAINS
!
!  Subroutine: DEL_HIST -------------------------------------------
!
!  Purpose: delete a history file -called if problems writing
!           out partial sums.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc
!------------------------------------------------------------------

SUBROUTINE del_hist()

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE filenamelength_mod, ONLY:                                    &
    filenamelength
USE MPPIO_file_utils, ONLY:                                      &
    file_delete

USE get_env_var_mod, ONLY: get_env_var

IMPLICIT NONE

INTEGER                       :: icode  ! error from get_file
CHARACTER(LEN=filenamelength) :: filename

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEL_HIST'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL get_env_var("HISTORY_TEMP",filename)
CALL file_delete(TRIM(filename))
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE del_hist
END MODULE del_hist_mod

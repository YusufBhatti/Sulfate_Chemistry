! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Reads model vertical namelist
MODULE read_lev_info_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READ_LEV_INFO_MOD'

CONTAINS

SUBROUTINE read_lev_info(file_lev)

! Central UM - vertical model namelist
USE vertnamelist_mod, ONLY:                                              &
    first_constant_r_rho_level, z_top_of_model                           &
    ,eta_theta, eta_rho , vertlevs

USE filenamelength_mod, ONLY:                                           &
  filenamelength

USE file_manager, ONLY: assign_file_unit, release_file_unit

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Description:
! Subroutine to read in level info - need to know number of levels
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
!-------------------------------------------------------------------------------

CHARACTER(LEN=filenamelength), INTENT(IN) ::  file_lev

!--------------------------------------------------------------------------
! Local variables

INTEGER :: error           ! 0 - no error in this routine
INTEGER :: STATUS          ! 0 - no error in this routine
INTEGER :: lev_unit

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_LEV_INFO'

!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!--------------------------------------------------------------------------
! Namelist for model levels - altered slighly to read as UM
!--------------------------------------------------------------------------

CALL assign_file_unit(TRIM(ADJUSTL(file_lev)), lev_unit, handler="fortran")

OPEN(UNIT=lev_unit,FILE=TRIM(ADJUSTL(file_lev)),ACTION='READ',IOSTAT=STATUS)

eta_theta(:)=0.0
eta_rho(:)  =0.0

READ (UNIT=lev_unit, NML=vertlevs, IOSTAT=Error)
! No UM routine for writing this namelist out after it has been read in.
! Can no longer use the Write(6) method
! WRITE(6,VERTLEVS)  ! write out namelist info

CLOSE(lev_unit)          ! close file

CALL release_file_unit(lev_unit, handler="fortran")
!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!--------------------------------------------------------------------------
RETURN
END SUBROUTINE read_lev_info

END MODULE read_lev_info_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads the coupling_control namelist
!
MODULE rcf_readnl_coupling_mod

IMPLICIT NONE

! Description:
!  rcf_readnl_coupling reads in the coupling_control namelist.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 programming standards.
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_COUPLING_MOD'

CONTAINS

SUBROUTINE rcf_readnl_coupling (nft)
USE coupling_control_mod, ONLY: read_nml_coupling_control, &
                               check_nml_coupling_control, &
                               print_nlist_coupling_control
USE umPrintMgr,  ONLY: PrintStatus, PrStatus_Oper
USE UM_parcore,  ONLY: mype
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)  :: nft  ! unit number
! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_READNL_COUPLING'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Read the namelist
CALL read_nml_coupling_control(nft)
IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_coupling_control()
END IF
CALL check_nml_coupling_control()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_readnl_coupling

END MODULE rcf_readnl_coupling_mod

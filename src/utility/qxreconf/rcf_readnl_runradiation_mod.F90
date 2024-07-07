! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads the RUN_Radiation namelist
!
MODULE rcf_readnl_runradiation_mod

IMPLICIT NONE

! Description:
!  rcf_readnl_runradiation reads in the RUN_Radiation namelist.
!
! Method:
!  Data read in via unit number passed as argument.
!  Only one element of this namelist is required by the
!  reconfiguration: i_rad_topography
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v8.5 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_RUNRADIATION_MOD'

CONTAINS

SUBROUTINE rcf_readnl_runradiation (nft)

USE rad_input_mod, ONLY:                                                 &
  RUN_Radiation, check_run_radiation, print_nlist_run_radiation,         &
  read_nml_run_radiation


USE umPrintMgr, ONLY:                                                    &
  PrintStatus, PrStatus_Oper

USE UM_parcore, ONLY:                                                    &
  mype

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
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_READNL_RUNRADIATION'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Defaults are set in the namelist module
! Read the namelist
CALL read_nml_run_radiation(nft)
IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_run_radiation()
END IF
CALL check_run_radiation()



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_readnl_runradiation

END MODULE rcf_readnl_runradiation_mod

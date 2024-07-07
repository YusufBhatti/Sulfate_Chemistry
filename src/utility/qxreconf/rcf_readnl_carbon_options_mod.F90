! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Read the CARBON_OPTIONS namelist

MODULE rcf_readnl_carbon_options_mod

IMPLICIT NONE

!  Subroutine Rcf_Readnl_Carbon_options - Read the CARBON_OPTIONS namelist
!
! Description:
!   Read the CARBON_OPTIONS namelist of Carbon control variables.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_READNL_CARBON_OPTIONS_MOD'

CONTAINS

SUBROUTINE rcf_readnl_carbon_options (nft)

USE umPrintMgr,   ONLY: PrintStatus, PrStatus_Oper
USE UM_ParCore,   ONLY: mype
USE carbon_options_mod, ONLY: check_carbon_options, print_nlist_carbon_options,&
     read_nml_carbon_options
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

INTEGER :: nft
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_READNL_CARBON_OPTIONS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL read_nml_carbon_options(nft)
IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_carbon_options()
END IF
CALL check_carbon_options()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE  Rcf_readnl_carbon_options
END MODULE Rcf_readnl_carbon_options_Mod

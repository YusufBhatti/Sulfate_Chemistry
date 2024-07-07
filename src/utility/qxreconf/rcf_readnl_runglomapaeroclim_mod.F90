! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Read the RUN_GLOMAP_AEROCLIM namelist

MODULE rcf_readnl_runglomapaeroclim_mod
IMPLICIT NONE

!  Subroutine Rcf_readnl_runglomapaeroclim - Read RUN_GLOMAP_AEROCLIM
!
! Description:
!   Read the RUN_GLOMAP_AEROCLIM namelist
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration_NetCDF
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_READNL_RUNGLOMAPAEROCLIM_MOD'

CONTAINS

SUBROUTINE rcf_readnl_runglomapaeroclim (nft)

USE glomap_clim_option_mod,      ONLY: &
    print_nlist_run_glomap_clim,       &
    read_nml_run_glomap_clim,          &
    check_glomap_clim_options

USE UM_ParCore,                  ONLY: mype

USE umPrintMgr,                  ONLY: &
    PrintStatus,                       &
    PrStatus_Oper

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

INTEGER :: nft
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_READNL_RUNGLOMAPAEROCLIM'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL read_nml_run_glomap_clim(nft)

IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_run_glomap_clim()
END IF

CALL check_glomap_clim_options()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE  rcf_readnl_runglomapaeroclim
END MODULE rcf_readnl_runglomapaeroclim_mod

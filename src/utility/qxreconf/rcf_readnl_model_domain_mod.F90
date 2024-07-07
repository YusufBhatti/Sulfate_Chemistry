! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Read the model_domain namelist

MODULE rcf_readnl_model_domain_mod

IMPLICIT NONE

!  Subroutine rcf_readnl_model_domain - Read the model_domain namelist
!
! Description:
!   Read the model_domain namelist of Atmosphere model control variables.
!
! Method:
!   Reads variables into the rcf_model_domain_mod module.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_MODEL_DOMAIN_MOD'

CONTAINS

SUBROUTINE rcf_readnl_model_domain(nft)

USE umPrintMgr,       ONLY: PrintStatus, PrStatus_Oper
USE UM_ParCore,       ONLY: mype
USE model_domain_mod, ONLY: read_nml_model_domain, print_nlist_model_domain, &
                            check_nml_model_domain
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim
IMPLICIT NONE

INTEGER :: nft
CHARACTER (LEN=50000) :: locMessage
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_READNL_MODEL_DOMAIN'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL read_nml_model_domain(nft)
IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_model_domain()
END IF
CALL check_nml_model_domain()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_readnl_model_domain
END MODULE rcf_readnl_model_domain_mod

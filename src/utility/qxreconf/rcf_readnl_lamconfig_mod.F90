! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads the lam_config namelist
!
MODULE rcf_readnl_lamconfig_mod

IMPLICIT NONE

! Description:
!  rcf_readnl_lamconfig reads in the lam_config namelist.
!
! Method:
!  Data read in via unit number passed as argument.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!  Language: FORTRAN 95.
!  This runcode is written to UMDP3 v8.6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_LAMCONFIG_MOD'

CONTAINS

SUBROUTINE rcf_readnl_lamconfig (nft)

USE lam_config_inputs_mod, ONLY: print_nlist_lam_config, read_nml_lam_config, &
        frstlona, frstlata, polelona, polelata, check_nml_lam_config

USE umPrintMgr, ONLY:                                                    &
     umPrint,                                                            &
     umMessage,                                                          &
     PrintStatus, PrStatus_Oper
USE UM_parcore, ONLY: mype

USE model_domain_mod, ONLY: model_type, mt_global

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
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_READNL_LAMCONFIG'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Set values for global jobs
IF ( model_type == mt_global) THEN
  frstlona =   0.0
  frstlata = -90.0
  polelona =   0.0
  polelata =  90.0
END IF

! Read the namelist
CALL read_nml_lam_config(nft)

IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_lam_config()
END IF
CALL check_nml_lam_config()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_readnl_lamconfig
END MODULE rcf_readnl_lamconfig_mod

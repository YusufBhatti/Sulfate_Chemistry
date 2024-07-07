! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads the Gen_phys_inputs namelist
!
MODULE rcf_readnl_gen_phys_inputs_mod

IMPLICIT NONE

! Description:
!  rcf_readnl_gen_phys_inputs reads in the Gen_phys_inputs namelist.
!
! Method:
!  Data read in via unit number passed as argument.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 programming standards.
CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_READNL_GEN_PHYS_INPUTS_MOD'

CONTAINS

SUBROUTINE rcf_readnl_gen_phys_inputs (nft)

USE gen_phys_inputs_mod, ONLY:                  &
  gen_phys_inputs, print_nlist_gen_phys_inputs, &
  read_nml_gen_phys_inputs

USE umPrintMgr, ONLY:                           &
  PrintStatus, PrStatus_Oper

USE UM_parcore, ONLY:                           &
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
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_READNL_GEN_PHYS_INPUTS'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Defaults are set in the namelist module
! Read the namelist
CALL read_nml_gen_phys_inputs(nft)

IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_gen_phys_inputs()
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_readnl_gen_phys_inputs
END MODULE rcf_readnl_gen_phys_inputs_mod

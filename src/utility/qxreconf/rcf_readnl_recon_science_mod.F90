! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_readnl_recon_science_Mod

!  Subroutine Rcf_Readnl_Recon_science - read the RECON_SCIENCE namelist
!
! Description:
!   Read the recon_science namelist, print out the values and check for
!   consistency.
!
! Method:
!   Variables read into the Rcf_Recon_science_Mod module.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

PRIVATE
PUBLIC :: rcf_readnl_recon_science

CHARACTER(LEN=*), PARAMETER :: ModuleName='RCF_READNL_RECON_SCIENCE_MOD'

CONTAINS

SUBROUTINE rcf_readnl_recon_science( nft )

USE rcf_nlist_recon_science_mod, ONLY: &
    w_zero_end,                        &
    w_zero_start,                      &
    print_nlist_recon_science,         &
    check_nml_recon_science,           &
    read_nml_recon_science

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

USE umPrintMgr, ONLY:       &
    PrintStatus,            &
    PrStatus_Oper

USE UM_ParCore, ONLY: &
    mype

USE nlsizes_namelist_mod, ONLY: &
    model_levels

IMPLICIT NONE

! Arguments
INTEGER                      :: nft     ! Unit number

! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_READNL_RECON_SCIENCE'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Read namelist
CALL read_nml_recon_science(nft)
IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_recon_science()
END IF

! If either w_zero_start or w_zero_end are -1 then set to
! model_levels (from same namelist) as default value.
IF (w_zero_start == -1) THEN
  w_zero_start = model_levels
END IF
IF (w_zero_end == -1) THEN
  w_zero_end = model_levels
END IF
IF (w_zero_start == 0) THEN  !  W at surface is always set to zero.
  w_zero_start = 1
END IF

CALL check_nml_recon_science()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_readnl_recon_science

END MODULE rcf_readnl_recon_science_mod

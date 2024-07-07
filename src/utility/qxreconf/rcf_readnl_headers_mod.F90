! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Read the HEADERS namelist

MODULE rcf_readnl_headers_mod

IMPLICIT NONE

!  Subroutine Rcf_Readnl_Headers - Read the HEADERS namelist
!
! Description:
!   Read the HEADERS namelist for header overrides.
!
! Method:
!   Variables read into Rcf_Headers_Mod module.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_HEADERS_MOD'

CONTAINS

SUBROUTINE rcf_readnl_headers( nft )

USE umPrintMgr, ONLY:      &
    umPrint,                &
    PrintStatus,            &
    PrStatus_Oper

USE UM_ParCore, ONLY: &
    mype

USE rcf_headers_mod, ONLY: &
    rcf_override_headers,   &
    print_nlist_headers,    &
    read_nml_headers,       &
    check_nml_headers

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT (IN)    :: nft
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'RCF_READNL_HEADERS'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Read Namelist
CALL read_nml_headers(nft)

! Write out namelist for checking/logging
IF (PrintStatus >= PrStatus_Oper) THEN
  IF (mype == 0) THEN
    CALL print_nlist_headers()
  END IF
END IF

! Perform checks on values read in namelist
CALL check_nml_headers()

CALL rcf_override_headers()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE rcf_readnl_headers

END MODULE rcf_readnl_headers_mod

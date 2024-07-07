! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Read the NLSTCALL namelist

MODULE rcf_readnl_nlstcall_mod

IMPLICIT NONE

!  Subroutine Rcf_Readnl_Nlstcall - Read the NLSTCALL namelist
!
! Description:
!   Read the NLSTCALL namelist of Atmosphere model control variables.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_NLSTCALL_MOD'

CONTAINS

SUBROUTINE rcf_readnl_nlstcall (nft)

USE umPrintMgr,   ONLY: PrintStatus, PrStatus_Oper
USE ereport_mod,  ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE UM_ParCore,   ONLY: mype
USE nlstcall_mod, ONLY: nlstcall, print_nlist_nlstcall,    &
     read_nml_nlstcall, ltimer
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

LOGICAL :: old_ltimer
INTEGER :: nft
INTEGER :: errorstatus
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'RCF_READNL_NLSTCALL'
CHARACTER (LEN=errormessagelength) :: cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
old_ltimer = ltimer
CALL read_nml_nlstcall(nft)

! LTimer can be overwritten by this namelist. We'll let the 
! RCF_TIMER version take precedence and tell users.
IF (ltimer .NEQV. old_ltimer) THEN
  ltimer = old_ltimer

  IF (ltimer) THEN
    WRITE(cmessage, '(A)')  'Using RCF_TIMER rather than ltimer.' // &
                               ' Timing ON'
  ELSE
    WRITE(cmessage, '(A)')  'Using RCF_TIMER rather than ltimer.' // &
                               ' Timing OFF'
  END IF
  errorstatus = -100     ! just a warning...
  CALL ereport(routinename, errorstatus, cmessage)
END IF


IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_nlstcall()
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE  Rcf_readnl_nlstcall
END MODULE Rcf_readnl_nlstcall_Mod

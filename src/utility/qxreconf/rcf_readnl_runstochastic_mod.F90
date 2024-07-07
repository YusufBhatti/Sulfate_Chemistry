! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads the RUN_STOCHASTIC namelist
!
MODULE rcf_readnl_runstochastic_mod

IMPLICIT NONE

! Description:
!  rcf_readnl_runstochastic reads in the RUN_STOCHASTIC namelist.
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

CONTAINS

SUBROUTINE rcf_readnl_runstochastic (nft)

USE stochastic_physics_run_mod, ONLY:                                    &
  print_nlist_run_stochastic, read_nml_run_stochastic

USE umPrintMgr, ONLY:                                                    &
  PrintStatus, PrStatus_Oper

USE UM_parcore, ONLY:                                                    &
  mype

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)  :: nft  ! unit number

! Local variables
INTEGER              :: ErrorStatus

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'rcf_readnl_runstochastic'

! Defaults are set in the namelist module
! Read the namelist
CALL read_nml_run_stochastic(nft)

IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_run_stochastic()
END IF

RETURN
END SUBROUTINE rcf_readnl_runstochastic

END MODULE rcf_readnl_runstochastic_mod

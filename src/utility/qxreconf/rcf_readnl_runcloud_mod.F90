! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads the RUN_Cloud namelist
!
MODULE rcf_readnl_runcloud_mod

IMPLICIT NONE

! Description:
!  rcf_readnl_runcloud reads in the RUN_Cloud namelist.
!
! Method:
!  Data read in via unit number passed as argument.
!  Only one element of this namelist is required by the reconfiguration: rhcrit.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v8.5 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_RUNCLOUD_MOD'

CONTAINS

SUBROUTINE rcf_readnl_runcloud (nft)

USE cloud_inputs_mod, ONLY:                                              &
  RUN_Cloud,                                                             &
  rhcrit, print_nlist_run_cloud, read_nml_run_cloud

USE rcf_grid_type_mod, ONLY:                                             &
  output_grid

USE atmos_max_sizes, ONLY:                                               &
  model_levels_max

USE umPrintMgr, ONLY:                                                    &
    umPrint,                                                             &
    umMessage,                                                           &
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
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_READNL_RUNCLOUD'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Defaults are set in the namelist module
! Read the namelist
CALL read_nml_run_cloud(nft)

IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_RUN_Cloud()
END IF

! Allocate space and fill up relevant part of output_grid
ALLOCATE ( output_grid % rhcrit( output_grid % model_levels ) )

output_grid % rhcrit( 1 : output_grid % model_levels ) =                 &
              rhcrit( 1 : output_grid % model_levels )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_readnl_runcloud

END MODULE rcf_readnl_runcloud_mod

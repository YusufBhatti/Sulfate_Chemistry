! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reads the NSUBMODL namelist

MODULE Rcf_Readnl_Nsubmodl_Mod

!  Subroutine Rcf_Readnl_Nsubmodl -
!
! Description:
!  Rcf_Readnl_Nsubmodl initialises the model with information specifying
!  internal model and submodel partitions for the run,
! These are hard-wired as we only have ATMOS
!
! Method:
!   Data read into the Submodel_Mod module.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_NSUBMODL_MOD'

CONTAINS

SUBROUTINE Rcf_Readnl_Nsubmodl( nft )

USE submodel_mod, ONLY:                                                        &
    internal_model_list, atmos_sm, submodel_for_sm, n_internal_for_sm,         &
    submodel_for_im, submodel_partition_index,                                 &
    submodel_partition_list, atmos_im

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Oper

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

INTEGER      ::  Nft     ! unit number for namelist

! Local constants
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_READNL_NSUBMODL'

! Local scalars:
INTEGER            ::  &
     s,                & ! submodel loop
     i,                & ! internal model loop
     sm,               & ! submodel identifier
     im,               & ! internal model identifier
     sm_prev,          & ! previous submodel identifier
     im_prev             ! previous internal model identifier

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Only Atmosphere Model catered for so no need to read in NSUBMODL.
internal_model_list(1) = atmos_im
submodel_for_im(1)     = atmos_sm

im = atmos_im
sm = atmos_sm
submodel_partition_list(1) = sm
submodel_for_sm(im) = 1
submodel_partition_index(1)=sm
n_internal_for_sm(sm)=1

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Readnl_Nsubmodl
END MODULE Rcf_Readnl_Nsubmodl_Mod

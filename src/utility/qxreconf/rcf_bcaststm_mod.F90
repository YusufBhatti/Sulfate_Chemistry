! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  broadcasts a STASHmaster record from 1 to many PEs

MODULE Rcf_bcastSTM_Mod
IMPLICIT NONE

!  Subroutine Rcf_bcastSTM - broadcasts a STM record
!
! Description:
!   Used to broadcast a STASHmaster record to all PEs for PE source.
!   An assumption is made that the STASHmaster records are
!   contiguous in memory (ie sequenced).
!
! Method:
!   Two calls to GCOM are used - 1 for integer and 1 for character
!   data. Sizes need to change if STASHmaster format is changed
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_BCASTSTM_MOD'

CONTAINS

SUBROUTINE Rcf_BcastSTM( STASHmaster , source )

USE Rcf_Ppx_Info_Mod, ONLY: &
    Ppxref_pack_profs,   &
    STM_IntDataLen, &
    STM_CharDataLen, &
    STM_record_type

USE UM_ParCore, ONLY:  &
    nproc

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER                :: source      ! PE that holds source record
TYPE (STM_record_type) :: STASHmaster ! The record to broadcast

! Local variables
INTEGER                :: msg         ! message tag
INTEGER                :: info        ! dummy
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_BCASTSTM'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!------------------------------------------------------------------
! Need 2 broadcasts to do whole thing!
!------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

msg = 8001
CALL gc_ibcast( msg, STM_IntDataLen, source, nproc, info, &
                STASHmaster % model )

msg = 8002
CALL gc_cbcast( msg, STM_CharDataLen, source, nproc, info, &
                STASHmaster % NAME )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_BcastSTM

END MODULE Rcf_BcastSTM_Mod

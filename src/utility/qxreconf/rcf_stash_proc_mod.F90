! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Interface to STASH processing

MODULE Rcf_Stash_Proc_Mod

IMPLICIT NONE
!  Subroutine Rcf_Stash_Proc - STASH processing
!
! Description:
! This routine is an interface to most of the stash processing -
! reading in the  stashmaster files and setting up the addressing
! of the output dump
!
! Method:
!   Does some initialisation, then calls rcf_getppx, submodl and
!   rcf_address.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_STASH_PROC_MOD'

CONTAINS


SUBROUTINE Rcf_Stash_Proc( EnvDir )

USE Rcf_Getppx_Mod, ONLY:     &
    Rcf_Getppx

USE Rcf_NRecon_Mod, ONLY:    &
    DumpProgLevs,             &
    PrimDataLen

USE Rcf_PPX_Info_Mod, ONLY: &
    STM_Record

USE h_vers_mod, ONLY: h_vers

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments
CHARACTER (LEN=*), INTENT(IN)   :: EnvDir ! Env. Variable for STASHmaster
                                          ! directory

! Local Data
CHARACTER (LEN=*), PARAMETER    :: RoutineName = 'RCF_STASH_PROC'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialisation of data length arrays
h_vers( :, : )  = 0

DumpProgLevs(:) = 0
PrimDataLen (:) = 0

NULLIFY( STM_Record )

! Read STASHmaster files into look-up arrays PPXI, PPXC
CALL Rcf_Getppx( 'STASHmaster_A', EnvDir )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Stash_Proc

END MODULE Rcf_Stash_Proc_Mod

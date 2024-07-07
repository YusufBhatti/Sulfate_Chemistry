! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initialisation of STASH related variables

MODULE Rcf_Stash_Init_Mod

!  Subroutine Rcf_Stash_Init - initialisation of STASH variables
!
! Description:
!   Reads headers for stashmasters and calls rcf_stash_proc to
!   read stashmaster and set up output dump addressing.
!
! Method:
!   Calls rcf_hdppxrf and rcf_stash_proc
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_STASH_INIT_MOD'

CONTAINS

SUBROUTINE Rcf_Stash_Init( )

USE filenamelength_mod, ONLY: &
    filenamelength

USE UM_ParCore, ONLY: &
    mype

USE Rcf_hdppxrf_mod, ONLY: &
    Rcf_hdppxrf

USE Rcf_Ppx_Info_Mod, ONLY: &
    ppxrecs

USE Rcf_Stash_Proc_mod, ONLY: &
    Rcf_Stash_Proc

USE submodel_mod, ONLY: internal_model_list, n_internal_model, atmos_im

USE umPrintMgr, ONLY:      &
    umPrint,  &
    umMessage,&
    PrintStatus,            &
    PrStatus_Normal

USE get_env_var_mod, ONLY: get_env_var

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Local variables
INTEGER                      :: i
INTEGER                      :: length    ! length of returned string
! STASHmaster
CHARACTER (LEN=11)   :: stashmaster = 'STASHMASTER'
! Temporary file to try reading env to.
CHARACTER (LEN=filenamelength) :: temp_name
! Routine name
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_STASH_INIT'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------
! Read in STASHmaster headers
!-------------------------------------------------------------------
DO i = 1, N_Internal_Model
  IF ( Internal_Model_List(i) == atmos_im ) THEN
    CALL Rcf_hdppxrf( 'STASHmaster_A',stashmaster )
  END IF
END DO

IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
  WRITE(umMessage,*) 'Rcf_Stash_Init : Total No of STASHmaster records ', &
               ppxRecs
  CALL umPrint(umMessage,src='rcf_stash_init_mod')
END IF

!---------------------------------------------------------------------
! Read in the STASHmasters, and do the output dump addressing
!---------------------------------------------------------------------
CALL Rcf_Stash_Proc( stashmaster )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Stash_Init

END MODULE Rcf_Stash_Init_Mod

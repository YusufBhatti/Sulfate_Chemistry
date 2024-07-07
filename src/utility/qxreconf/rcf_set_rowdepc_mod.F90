! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  sets up the output dump row dependent constants

MODULE Rcf_Setup_RowDepC_Mod

!  Subroutine Rcf_Setup_RowDepC - sets up the output dump
!                                 row dependent contants.
!
! Description:
!   The row dependent constants for the output dump are constructed.
!
! Method:
!   The row dependent constants are constructed from the output
!   grid details - see UMDP F3 for full  details.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SETUP_ROWDEPC_MOD'

CONTAINS
SUBROUTINE Rcf_Setup_RowDepC( Hdr_Out, Grid )

USE Rcf_Grid_Type_Mod, ONLY: &
    Grid_Type

USE Rcf_UMhead_Mod, ONLY: &
    UM_Header_Type

USE Rcf_Headaddress_Mod, ONLY:  &
    RDC_Phi_input_p, RDC_Phi_input_v

USE vrhoriz_grid_mod, ONLY: &
    Phi_input_p, Phi_input_v

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( UM_Header_Type ), INTENT(INOUT) :: Hdr_Out
TYPE( Grid_Type ), INTENT(IN)         :: Grid      ! Output Grid

! Local Variables
INTEGER                          :: i        ! Looper

CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SETUP_ROWDEPC'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i = 1,  Hdr_Out % Len1RowDepC
  Hdr_Out % RowDepC( i, RDC_Phi_input_p) =  Phi_input_p(i)
  Hdr_Out % RowDepC( i, RDC_Phi_input_v) =  Phi_input_v(i)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Setup_RowDepC
END MODULE Rcf_Setup_RowDepC_Mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  sets up the output dump column dependent constants

MODULE Rcf_Setup_ColDepC_Mod

!  Subroutine Rcf_Setup_ColDepC - sets up the output dump
!                                 column dependent contants.
!
! Description:
!   The column dependent constants for the output dump are constructed.
!
! Method:
!   The column dependent constants are constructed from the output
!   grid details - see UMDP F3 for full  details.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SETUP_COLDEPC_MOD'

CONTAINS
SUBROUTINE Rcf_Setup_ColDepC( Hdr_Out, Grid )

USE Rcf_Grid_Type_Mod, ONLY: &
    Grid_Type

USE Rcf_UMhead_Mod, ONLY: &
    UM_Header_Type

USE Rcf_Headaddress_Mod, ONLY:  &
    CDC_Lambda_input_p, CDC_Lambda_input_u

USE vrhoriz_grid_mod, ONLY: &
    Lambda_input_p, Lambda_input_u

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
INTEGER                         :: i        ! Looper
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_SETUP_COLDEPC'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i = 1,  Hdr_Out % Len1ColDepC
  Hdr_Out % ColDepC( i, CDC_Lambda_input_p) =  Lambda_input_p(i)
  Hdr_Out % ColDepC( i, CDC_Lambda_input_u) =  Lambda_input_u(i)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Setup_ColDepC
END MODULE Rcf_Setup_ColDepC_Mod

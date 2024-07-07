! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reconfiguration Ancillary Processing

MODULE Rcf_Ancil_Mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

!  Subroutine Rcf_Ancil  - Reconfiguration Ancillary Processing
!
! Description:
!    Controls ancillary processing in the reconfiguration
!
! Method:
!    Reserves an unit number for the ancillary files and
!    calls Rcf_Ancil_Atmos for atmosphere ancillary processing.
!    Ocean ancillary processing to be added later.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_ANCIL_MOD'

CONTAINS

SUBROUTINE Rcf_Ancil ( Hdr_In, Hdr_Out,              &
                       Fields_In, Field_Count_In,    &
                       Fields_Out, Field_Count_Out, data_source )

USE rcf_ancil_atmos_mod, ONLY: &
    rcf_ancil_atmos

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE file_manager, ONLY: &
    assign_file_unit, &
    release_file_unit

USE Ancil_mod, ONLY: &
    AncF_UnitNo

USE nlstcall_mod, ONLY:   &
    LTimer

USE rcf_data_source_mod, ONLY: &
    data_source_type

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (um_header_type), INTENT(IN)   :: Hdr_In
TYPE (um_header_type), INTENT(IN)   :: Hdr_Out
TYPE (field_type), POINTER          :: Fields_In (:)
TYPE (field_type), POINTER          :: Fields_Out (:)
INTEGER, INTENT(IN)                 :: Field_Count_In
INTEGER, INTENT(IN)                 :: Field_Count_Out
TYPE (data_source_type), POINTER    :: data_source(:)

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_ANCIL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (LTimer) CALL Timer( 'Rcf_Ancil', 3)
!-------------------------------------------------------------------
! Get a unit number for Ancillary files
!-------------------------------------------------------------------
CALL assign_file_unit ( "rcf_ancillary_files_unit", &
                         Ancf_UnitNo, handler="portio" )

!-------------------------------------------------------------------
! Do processing for Atmosphere Ancillaries
!-------------------------------------------------------------------
CALL Rcf_Ancil_Atmos ( Hdr_In, Hdr_Out,                &
                         Fields_In, Field_Count_In,      &
                         Fields_Out, Field_Count_Out, data_source )

!-------------------------------------------------------------------
! Free the unit number for Ancillary files
!-------------------------------------------------------------------
CALL release_file_unit ( Ancf_UnitNo, handler="portio" )

IF (LTimer) CALL Timer( 'Rcf_Ancil', 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Ancil

END MODULE Rcf_Ancil_Mod

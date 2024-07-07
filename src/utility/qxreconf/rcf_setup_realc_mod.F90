! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sets up the real constants in the output header

MODULE Rcf_setup_RealC_Mod

!  Subroutine Rcf_Setup_RealC - initialises real constants in header.
!
! Description:
!   Sets the real constants in the output dump header.
!
! Method:
!    Uses namelist information to set header.
!    UMDP F3 defines the real constants.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SETUP_REALC_MOD'

CONTAINS

SUBROUTINE Rcf_Setup_RealC( Hdr_In, Hdr_Out )

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Oper

USE Rcf_UMhead_Mod, ONLY: &
    Um_Header_type

USE rcf_headers_mod, ONLY: &
    RelHd

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE lam_config_inputs_mod, ONLY: &
    frstlona,  frstlata,         &
    polelona,  polelata,         &
    delta_lon, delta_lat

USE Rcf_HeadAddress_Mod, ONLY:                    &
    RC_PoleLong,                  RC_LongSpacing,  &
    RC_LatSpacing,                RC_PoleLat,      &
    RC_FirstLong,                 RC_FirstLat,     &
    RC_ModelTop,                  rc_swldeg,       &
    RC_WEdgeDeg

USE missing_data_mod, ONLY: rmdi 

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (Um_Header_Type), INTENT(IN) :: Hdr_In
TYPE (Um_Header_Type), TARGET       :: Hdr_Out

! Local vars
REAL, POINTER                  :: RealC(:)
INTEGER                        :: i

CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SETUP_REALC'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------------------------------------
! Clean start - RMDI for everything
!---------------------------------------------------------------
RealC => Hdr_Out % RealC
RealC(:) = rmdi

!--------------------------------------------------------------
! Set values we have numbers for
!--------------------------------------------------------------

RealC( RC_LongSpacing ) = delta_lon
RealC( RC_LatSpacing  ) = delta_lat
RealC( RC_FirstLat    ) = frstlata
RealC( RC_FirstLong   ) = frstlona
RealC( RC_PoleLat     ) = polelata
RealC( RC_PoleLong    ) = polelona

! Submodel specifics
RealC( RC_ModelTop    ) = Output_Grid % z_top_of_model

!--------------------------------------------------------------
! Overrides from module
!--------------------------------------------------------------
DO i = 1, Hdr_Out % LenRealC
  IF ( RelHd(i) /= rmdi ) THEN
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,*) 'RealC(',i,') has been reset from ', RealC(i), &
                  ' to ', RelHd(i)
      CALL umPrint(umMessage,src='rcf_setup_realc_mod')
    END IF

    RealC(i) = RelHd(i)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Setup_RealC
END MODULE Rcf_Setup_RealC_Mod

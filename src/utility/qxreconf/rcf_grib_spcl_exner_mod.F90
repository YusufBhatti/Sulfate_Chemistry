! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculate Exner from pressure levels

MODULE Rcf_Grib_Spcl_Exner_Mod

! SUBROUTINE Rcf_Grib_Spcl_Exner
!
! Description: Calculate Exner from pressure level.
!
! Method: Using (P/P0)^k
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_SPCL_EXNER_MOD'

CONTAINS
SUBROUTINE Rcf_Grib_Spcl_Exner(Current,FpData)

USE Rcf_GRIB_Block_Params_Mod, ONLY: &
    Grib_Record,     LenArrayMax,        &
    p_Lvl_Type,      p_Lvl_Desc_1,       &
    Tb3_Pressure

USE EReport_Mod, ONLY:     &
    EReport

USE planet_constants_mod, ONLY: p_zero, kappa

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Declarations:

! Subroutine arguments

!< Scalar arguments with intent(In):>

!< Array  arguments with intent(InOut):>
TYPE (Grib_Record),POINTER          :: Current
REAL, INTENT(INOUT)                 :: FpData(LenArrayMax)


! Local constants
CHARACTER (LEN=*), PARAMETER     :: RoutineName='RCF_GRIB_SPCL_EXNER'

! Local variables

INTEGER                          :: Pressure
REAL                             :: Exner
CHARACTER (LEN=errormessagelength)  :: Cmessage   ! used for EReport
INTEGER                          :: ErrorStatus   ! used for EReport
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!=======================================================================
!  Routine Code Start :
!=======================================================================

!=======================================================================
!  Data on a Pressure level
!=======================================================================

IF (Current % Block_1 (p_Lvl_Type) == Tb3_Pressure ) THEN

  ! Get pressure value for level
  Pressure = Current % Block_1 (p_Lvl_Desc_1) * 100 ! convert from HPa

  ! Calculate Exner for Level
  FpData(:) = ( Pressure / P_zero ) ** kappa

  !=======================================================================
  !  Data _NOT_ on a Pressure level
  !=======================================================================

ELSE

  Cmessage = 'Cant handle level type other than Pressure'
  ErrorStatus = 10
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF
!=======================================================================
!  Code End :
!=======================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Grib_Spcl_Exner
END MODULE Rcf_Grib_Spcl_Exner_Mod

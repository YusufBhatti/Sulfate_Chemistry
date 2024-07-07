! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Converts between P and exner (on rho levels)

MODULE Rcf_Exner_P_Convs_Mod

!  Subroutine Rcf_Conv_Exner_P - converts exner to P
!  Subroutine Rcf_Conv_P_Exner - converts P to exner
!
! Description:
!   Performs conversions between exner and P on rho levels.
!
! Method:
!   Data is stored in the *original* field data pointer!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE submodel_mod, ONLY: atmos_im

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_EXNER_P_CONVS_MOD'

CONTAINS

! *****************************************************************
! This routine converts Exner to P
! *****************************************************************
SUBROUTINE Rcf_Conv_Exner_P( exner_field )

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE um_stashcode_mod, ONLY: &
    stashcode_exner,           &
    stashcode_p,               &
    stashcode_prog_sec

USE Rcf_Exppx_Mod, ONLY: &
    Rcf_Exppx

USE planet_constants_mod, ONLY: pref, recip_kappa

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Arguments
TYPE( field_type), INTENT(INOUT)    :: exner_field


!  Local variables
INTEGER                        :: i
INTEGER                        :: k
INTEGER                        :: ErrorStatus
CHARACTER (LEN=*), PARAMETER   :: RoutineName = 'RCF_CONV_EXNER_P'
CHARACTER (LEN=errormessagelength)       :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------
! Make sure field actually is exner...
!-----------------------------------------------------------------
IF ( exner_field % stashmaster % section /= stashcode_prog_sec .OR.  &
     exner_field % stashmaster % item    /= stashcode_exner ) THEN
  ErrorStatus = 10
  Cmessage = 'Stashcode does not match expected (exner) data'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!------------------------------------------------------------------
! Convert exner to P
!------------------------------------------------------------------
DO k = 1, exner_field % levels
  DO i = 1, exner_field % level_size
    exner_field % DATA(i,k) =                                       &
                  (exner_field % DATA(i,k) ** recip_kappa)* pref
  END DO
END DO

!------------------------------------------------------------------
! Set field stashmaster to P
!------------------------------------------------------------------
exner_field % stashmaster => Rcf_Exppx( atmos_im, stashcode_prog_sec, &
                                        stashcode_p )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Conv_Exner_P


! ******************************************************************
! Routine to convert P to exner
! ******************************************************************

SUBROUTINE Rcf_Conv_P_Exner( exner_field )

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE um_stashcode_mod, ONLY: &
    stashcode_exner,           &
    stashcode_p,               &
    stashcode_prog_sec

USE Rcf_Exppx_Mod, ONLY: &
    Rcf_Exppx


USE planet_constants_mod, ONLY: kappa, recip_p_zero

USE errormessagelength_mod, ONLY: errormessagelength

USE ereport_mod, ONLY: ereport
IMPLICIT NONE

! Arguments
TYPE( field_type), INTENT(INOUT)    :: exner_field

!  Local variables
INTEGER                        :: i
INTEGER                        :: k
INTEGER                        :: ErrorStatus
CHARACTER (LEN=*), PARAMETER   :: RoutineName = 'RCF_CONV_P_EXNER'
CHARACTER (LEN=errormessagelength)  :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------
! Make sure field actually is P...
!-----------------------------------------------------------------
IF ( exner_field % stashmaster % section /= stashcode_prog_sec .OR.   &
     exner_field % stashmaster % item    /= stashcode_p ) THEN
  ErrorStatus = 10
  Cmessage = 'Stashcode does not match expected (pressure) data'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!------------------------------------------------------------------
! Convert P to Exner
!------------------------------------------------------------------
DO k = 1, exner_field % levels
  DO i = 1, exner_field % level_size
    exner_field % DATA(i,k) =                                       &
                  (exner_field % DATA(i,k) * recip_p_zero) ** kappa
  END DO
END DO

!------------------------------------------------------------------
! Set field stashmaster to exner
!------------------------------------------------------------------
exner_field % stashmaster => Rcf_Exppx(atmos_im, stashcode_prog_sec,  &
                                       stashcode_exner)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Conv_P_Exner

END MODULE Rcf_Exner_P_Convs_Mod

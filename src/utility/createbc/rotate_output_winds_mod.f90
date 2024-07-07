! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE rotate_output_winds_mod
!
! Description:
! Adjust the LBC output winds vectors to account for a rotated pole.
! 
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ROTATE_OUTPUT_WINDS_MOD'

CONTAINS

SUBROUTINE rotate_output_winds(wind_rotation_coeff, output_u_field, output_v_field)

USE wind_rotation_coeff_mod,   ONLY: wind_rotation_coeff_type
USE field_mod,                 ONLY: field_type, ASSIGNMENT(=)
USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength
IMPLICIT NONE
! Arguments
TYPE(wind_rotation_coeff_type), INTENT(IN) :: wind_rotation_coeff
TYPE(field_type), INTENT(INOUT) :: output_u_field
TYPE(field_type), INTENT(INOUT) :: output_v_field
! Local variables
INTEGER :: icode, i, j, counter, k, num_levels
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'ROTATE_OUTPUT_WINDS'
REAL :: orig_u_wind, orig_v_wind ! Temp values for winds

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(wind_rotation_coeff%output_wind_coeff1)) THEN
  icode = 10
  cmessage = "Output wind rotation coefficient 1 array is not allocated."
  CALL ereport(routinename, icode, cmessage)
END IF
IF (.NOT. ALLOCATED(wind_rotation_coeff%output_wind_coeff2)) THEN
  icode = 20
  cmessage = "Output wind rotation coefficient 2 array is not allocated."
  CALL ereport(routinename, icode, cmessage)
END IF
IF (SIZE(output_u_field%lbc_rdata,1) /= SIZE(output_v_field%lbc_rdata,1)) THEN
  icode = 30
  WRITE(cmessage, '(A,I0,A,I0)') "Field sizes do not match. single level u field size = ",           &
       SIZE(output_u_field%lbc_rdata,1), " single level v field size = ",                            &
       SIZE(output_v_field%lbc_rdata,1) 
  CALL ereport(routinename, icode, cmessage)
END IF

! Need to check single level is same size as wind coefficents
IF (SIZE(output_u_field%lbc_rdata,1) /= SIZE(wind_rotation_coeff%output_wind_coeff1)) THEN
  icode = 40
  WRITE(cmessage, '(A,I0,A,I0)') "Wind field data array and wind coefficient array " //              &
       "sizes do not match. u field size = ",SIZE(output_u_field%rdata,1),                           &
       " wind coefficient array size = ", SIZE(wind_rotation_coeff%output_wind_coeff1)
  CALL ereport(routinename, icode, cmessage)
END IF

IF (output_u_field%get_num_levels() /= output_v_field%get_num_levels()) THEN
  icode = 50
  WRITE(cmessage, '(A,I0,A,I0)') "Number of levels of U and V fields do not match. Num U levels = ", &
       output_u_field%get_num_levels(), " num V levels = ", output_v_field%get_num_levels()
  CALL ereport(routinename, icode, cmessage)
END IF

! Avoid using the function procedure inside the PARALLEL region to avoid
! potential issues with race conditions (ifort), compiler bugs (gfortran), etc.
num_levels = output_u_field%get_num_levels()

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k, i, orig_u_wind, orig_v_wind) &
!$OMP SHARED(num_levels, output_u_field, output_v_field, wind_rotation_coeff)
DO k = 1, num_levels
  DO i = 1, SIZE(output_u_field%lbc_rdata,1)
    orig_u_wind = output_u_field%lbc_rdata(i,k)
    orig_v_wind = output_v_field%lbc_rdata(i,k)
    output_u_field%lbc_rdata(i,k) = (orig_u_wind * wind_rotation_coeff%output_wind_coeff1(i)) - &
                                    (orig_v_wind * wind_rotation_coeff%output_wind_coeff2(i)) 
    output_v_field%lbc_rdata(i,k) = (orig_v_wind * wind_rotation_coeff%output_wind_coeff1(i)) + &
                                    (orig_u_wind * wind_rotation_coeff%output_wind_coeff2(i))
  END DO
END DO
!$OMP END PARALLEL DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE rotate_output_winds

END MODULE rotate_output_winds_mod

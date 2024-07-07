! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE unrotate_input_winds_mod
! Description:
!  Apply the wind rotation coefficients to the U and V wind.
!  U and V winds must be on the same grid. This should transform
!  the wind components from a rotated grid to a standard polar
!  axis.
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

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UNROTATE_INPUT_WINDS_MOD'

CONTAINS

SUBROUTINE unrotate_input_winds(input_u_field, input_v_field, wind_rotation_coeff)

USE field_mod,                 ONLY: field_type, ASSIGNMENT(=)
USE wind_rotation_coeff_mod,   ONLY: wind_rotation_coeff_type
USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength

IMPLICIT NONE

TYPE(field_type), INTENT(INOUT) :: input_u_field
TYPE(field_type), INTENT(INOUT) :: input_v_field
TYPE(wind_rotation_coeff_type), INTENT(IN) :: wind_rotation_coeff
! Local variables
INTEGER :: icode, i, j, counter, k
INTEGER :: num_levels, num_rows, num_cols
REAL :: orig_rotated_u_wind, orig_rotated_v_wind ! Temp values for winds
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'UNROTATE_INPUT_WINDS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL timer( 'unrotate_input_winds', 5)

IF (.NOT. ALLOCATED(wind_rotation_coeff%input_wind_coeff1)) THEN
  icode = 10
  cmessage = "Input wind rotation coefficient 1 array is not allocated."
  CALL ereport(routinename, icode, cmessage)
END IF
IF (.NOT. ALLOCATED(wind_rotation_coeff%input_wind_coeff2)) THEN
  icode = 20
  cmessage = "Input wind rotation coefficient 2 array is not allocated."
  CALL ereport(routinename, icode, cmessage)
END IF

IF (SIZE(input_u_field%rdata) /= SIZE(input_v_field%rdata)) THEN
  icode = 30
  WRITE(cmessage, '(A,I0,A,I0)') "Field sizes do not match. u field size = ",    &
       SIZE(input_u_field%rdata), " v field size = ",  SIZE(input_v_field%rdata) 
  CALL ereport(routinename, icode, cmessage)
END IF

! Need to check single level size x * y
IF (SIZE(input_u_field%rdata,1)*SIZE(input_u_field%rdata,2) & 
     /= SIZE(wind_rotation_coeff%input_wind_coeff1)) THEN
  icode = 40
  WRITE(cmessage, '(A,I0,A,I0)') "Wind field data array and wind coefficient array " // &
       "sizes do not match. u field size = ",                                           &
       SIZE(input_u_field%rdata,1) *  SIZE(input_u_field%rdata,2),                      &
       " wind coefficient array size = ", SIZE(wind_rotation_coeff%input_wind_coeff1)
  CALL ereport(routinename, icode, cmessage)
END IF

IF (input_u_field%get_num_levels() /= input_v_field%get_num_levels()) THEN
  icode = 50
  WRITE(cmessage, '(A,I0,A,I0)') "Number of levels of U and V fields do not match. Num U levels = ", &
       input_u_field%get_num_levels(), " num V levels = ", input_v_field%get_num_levels()
  CALL ereport(routinename, icode, cmessage)
END IF

! Avoid using these function procedures inside the PARALLEL region to avoid a
! gfortran compiler bug (tested up to gfortran v6.3.0), potential race
! conditions with ifort and other possible compiler issues.
num_levels = input_u_field%get_num_levels()
num_rows   = input_u_field%get_num_rows()
num_cols   = input_u_field%get_num_cols()

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k, i, j, counter, orig_rotated_u_wind, orig_rotated_v_wind) &
!$OMP SHARED(input_u_field, input_v_field, wind_rotation_coeff, num_levels, num_rows, num_cols)
DO k = 1, num_levels
  counter = 0
  DO j = 1, num_rows
    DO i = 1, num_cols
      counter = counter + 1
      orig_rotated_u_wind = input_u_field%rdata(i,j,k)
      orig_rotated_v_wind = input_v_field%rdata(i,j,k)
      input_u_field%rdata(i,j,k) = (orig_rotated_u_wind *                             &
                                    wind_rotation_coeff%input_wind_coeff1(counter)) + &
                                   (orig_rotated_v_wind *                             &
                                    wind_rotation_coeff%input_wind_coeff2(counter)) 
      input_v_field%rdata(i,j,k) = (orig_rotated_v_wind *                             &
                                    wind_rotation_coeff%input_wind_coeff1(counter)) - &
                                   (orig_rotated_u_wind *                             &
                                    wind_rotation_coeff%input_wind_coeff2(counter))
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
CALL timer( 'unrotate_input_winds', 6)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE unrotate_input_winds
END MODULE unrotate_input_winds_mod

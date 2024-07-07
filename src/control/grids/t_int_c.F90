! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE T_INT_C-------------------------------------------------
!
! Purpose:
!   Carries out linear interpolation in time between two fields at
!   times T1 and T2. If the missing data indicator is present at one
!   of the times, the value at the other time is used. The interpolation
!   is controlled by a field ZI. A prescribed value is inserted where ZI=0.
!   If ZI changes between 0 and non-zero in the period T1 - T2, then
!   the field is linearly interpolated between its value at the time when
!   ZI is non-zero and the prescibed value at the time when ZI becomes zero.
!   The fractional time at which ZI changes between 0 and non-zero in
!   period T1 - T2 must be provided as input.
!
!
! Programming standard: UMDP3 v10.3
!
! System component: S191
!
! System task: S1
!
! Documentation:  The interpolation formulae are described in
!                 unified model on-line documentation paper S1.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids
MODULE t_int_c_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'T_INT_C_MOD'

CONTAINS

SUBROUTINE t_int_c ( data_t1, t1, data_t2, t2, data_t3, t3, points             &
                   , frac_time, zi_t1, pres_value )


USE missing_data_mod, ONLY: rmdi
USE t_int_mod,        ONLY: t_int

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: points           ! Number of points to be processed

REAL, INTENT(IN)  :: data_t1(points)    ! Data at T1
REAL, INTENT(IN)  :: data_t2(points)    ! Data at T2
REAL, INTENT(IN)  :: zi_t1(points)      ! Value of controlling fieled at T1
REAL, INTENT(IN)  :: pres_value(points) ! Prescribed value of Data when ZI=0
REAL, INTENT(IN)  :: frac_time(points)  ! Fractional time at which ZI changes
                                        ! between zero and non-zero in this
                                        ! time range

REAL, INTENT(IN)  :: t1    ! Time of first data field
REAL, INTENT(IN)  :: t2    ! Time of second data field
REAL, INTENT(IN)  :: t3    ! Time at which new field is required T1<=T3<=T2

REAL, INTENT(OUT) :: data_t3(points)  ! Data at T3

! Local arrays:---------------------------------------------------------
REAL ::  int_time(points)

! Local variables:------------------------------------------------------
REAL :: alpha         ! Fractional time
REAL :: epsi          ! Rounding error
REAL :: alpha_plus    ! Add rounding error to alpha
REAL :: alpha_minus   ! Alpha minus rounding error

INTEGER :: i          ! Loop index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'T_INT_C'

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set rounding error
epsi = 1.0e-6

CALL t_int(data_t1,t1,data_t2,t2,data_t3,t3,points)

alpha = (t3-t1)/(t2-t1)
alpha_plus  = alpha + epsi
alpha_minus = alpha - epsi

DO i=1, points

  IF (frac_time(i) /= rmdi) THEN
    IF (zi_t1(i) == 0.0) THEN
      IF (alpha_minus < frac_time(i)) THEN
        data_t3(i) = pres_value(i)
      ELSE
        int_time(i) = (alpha-frac_time(i))/(1.0-frac_time(i))
        data_t3(i)  = pres_value(i)*(1.0-int_time(i))                          &
                    + data_t2(i)*int_time(i)
      END IF
    ELSE
      IF (alpha_plus > frac_time(i)) THEN
        data_t3(i)  = pres_value(i)
      ELSE
        int_time(i) = (frac_time(i)-alpha)/(frac_time(i))
        data_t3(i)  = pres_value(i)*(1.0-int_time(i))                          &
                    + data_t1(i)*int_time(i)
      END IF
    END IF
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE t_int_c

END MODULE t_int_c_mod

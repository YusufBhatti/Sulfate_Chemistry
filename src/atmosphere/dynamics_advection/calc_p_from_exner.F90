! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Calc_Exner_at_theta

!
! Subroutine Calc_P_from_Exner
!

SUBROUTINE Calc_P_from_Exner(                                     &
                      p,                                          &
                      row_length, rows, levels,             &
                      off_x, off_y,                               &
                      exner,L_include_halos)

! Purpose:
!          Calculates pressure from exner pressure.
!
! Method:
!          p=p_zero*exner**(1/kappa)
!
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE atm_fields_bounds_mod, ONLY: tdims_s

USE planet_constants_mod, ONLY:  recip_kappa, p_zero

USE vectlib_mod, ONLY:                                                 &
      powr_v

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER ::                                                        &
  row_length                                                      &
                   ! number of point on a row.
, rows                                                            &
                   ! number of rows.
, levels                                                    &
                   ! number of model levels.
, off_x                                                           &
             ! Size of small halo in i
, off_y      ! Size of small halo in j.
LOGICAL ::                                                        &
  L_include_halos  ! If .TRUE. then include halo regions
                   ! when performing the calculations

REAL ::                                                           &
  exner(1-off_x:row_length+off_x,                                 &
        1-off_y:rows+off_y, levels)

! Arguments with Intent OUT. ie: Output variables.
REAL ::                                                          &
 p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     levels)

! Local Variables.

INTEGER ::                                                        &
  i, j, k                                                         &
               ! Loop indices
,  i_start,i_end                                                  &
                   ! Loop bounds
,  j_start,j_end                                                  &
                   ! Loop bounds
,  vector_length   !


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_P_FROM_EXNER'


! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Calculate p from exner.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (L_include_halos) THEN
  i_start=1-off_x
  i_end=row_length+off_x
  j_start=1-Off_y
  j_end=rows+Off_y
ELSE
  i_start=1
  i_end=row_length
  j_start=1
  j_end=rows
END IF

vector_length = i_end -  i_start +1

DO k = 1, levels
  DO j = j_start, j_end
    CALL powr_v(vector_length,exner(i_start,j,k),recip_kappa,     &
                p(i_start,j,k))
    DO i = i_start, i_end
      p(i,j,k) =   p_zero * p(i,j,k)
    END DO
  END DO
END DO

! End of routine
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Calc_P_from_Exner

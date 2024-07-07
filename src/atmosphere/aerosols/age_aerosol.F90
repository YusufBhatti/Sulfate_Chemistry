! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Convert a proportion of fresh aerosol to aged aerosol
MODULE age_aerosol_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='AGE_AEROSOL_MOD'

CONTAINS

SUBROUTINE age_aerosol(                                                 &
  row_length, rows, off_x, off_y,                                       &
  model_levels, aerosol_in,                                             &
  delta_aerosol, flag  )

! Purpose:
!   To convert a proportion of the fresh aerosol to aged aerosol. This
!   conversion takes place as an exponential decay with an e-folding
!   time which defines the rate.
! bmas is mmr of fresh smoke
! ocff is mmr of ocff
! soot is mmr of soot
!
!   Called by Aero_ctl
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
!
! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP3 v8 programming standards
!
! Documentation: UMDP20

USE timestep_mod, ONLY: timestep

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with intent IN:

INTEGER ::  row_length         ! no. of pts along a row
INTEGER ::  rows               ! no. of rows
INTEGER ::  off_x              ! size of small halo in i
INTEGER ::  off_y              ! size of small halo in j.
INTEGER ::  model_levels       ! no. of model levels

! Arguments with intent IN:
REAL  :: aerosol_in(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)

CHARACTER(LEN=4) :: flag

! Arguments with intent OUT:
!aerosol increment due to ageing
REAL  :: delta_aerosol(row_length, rows, model_levels)

! Local variables:

INTEGER            ::  i,j,k  ! loop counters
REAL               ::  rate
REAL               ::  a   ! Local workspace, equal to 1.0-exp(-rate*timestep)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='AGE_AEROSOL'

! Cycle through all points in the field, calculating the amount
! of aerosol converted on each point on this timestep.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (flag == 'bmas') THEN
! conversion rate. Equivalent to 1/(6.0 hours expressed in seconds)
  rate = 4.6296e-5
ELSE IF (flag == 'ocff') THEN
! conversion rate. Equivalent to 1/(1.0 day expressed in seconds)
  rate = 1.157e-5
ELSE IF (flag == 'soot') THEN
! conversion rate equivalent to 1/(1.0 days expressed in seconds)
  rate = 1.157e-5
END IF

a = 1.0 - EXP(-rate*timestep)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i, j, k)           &
!$OMP SHARED(model_levels, rows, row_length, aerosol_in, delta_aerosol, a)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      delta_aerosol(i,j,k) = a * aerosol_in(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE age_aerosol

END MODULE age_aerosol_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Energy Correction
MODULE flux_diag_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FLUX_DIAG_MOD'
CONTAINS

SUBROUTINE flux_diag (flux,area,row_length, rows,off_x,off_y,     &
                      conv_fac,flux_sum, tstep)

!    SUBROUTINE FLUX_DIAG----------------------------------------------
!
!    PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
!              TO SUM FLUXES INTO THE ATMOSPHERE
!              Note global summation now done elsewhere at end of
!              energy correction period to save CPU.
!
!    NOT SUITABLE FOR SINGLE COLUMN MODEL USE
!
!    PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!
!    DOCUMENTATION :
!
!
USE planet_constants_mod, ONLY: planet_radius

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!----------------------------------------------------------------------
! VECTOR LENGTHS
!----------------------------------------------------------------------
!
INTEGER :: row_length, rows                                          &
,     off_x,off_y

!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------

REAL :: flux(row_length, rows) ! IN FLUX TO BE SUMMED
!
REAL :: area(1-off_x:row_length+off_x,                               &
          1-off_y:rows+off_y) ! IN AREA OF GRID BOX
!
REAL :: tstep               ! IN TIMESTEP
!
REAL :: conv_fac            ! IN CONVERSION FACTOR TO TRANSLATE
                         !    FLUX INTO ENERGY UNITS
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE IN AND OUT
!----------------------------------------------------------------------
!

REAL :: flux_sum(row_length,rows)  ! INOUT sum of fluxes into
                                ! the atmosphere.
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
!
REAL :: factor        ! local factor

INTEGER :: i,j ! LOOP COUNTER

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FLUX_DIAG'

!
!
!----------------------------------------------------------------------
! EXTERNAL SUBROUTINE CALLS  -  NONE
! ---------------------------------------------------------------------
!  Flux_sum holds the summation of fluxes into the atmosphere as a full
!l global field weighted by area. The fluxes are summed over the period
!  A_energysteps. At the end of this period the global sum is formed.
!  Global summation is costly on a MPP machine and is therefore only
!  done at the end of the period instead of every time fluxes are added.
!  Note this appraoch means the array flux_sum must be held in memory
!  adding to the memory overhead.
!
!  Add additional fluxes to current flux_sum.
! ---------------------------------------------------------------------
!  factor to multiply flux

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
factor = planet_radius*planet_radius*conv_fac*tstep

DO j=1,rows
  DO i=1,row_length
    flux_sum(i,j) = flux_sum(i,j)+                                &
                           factor*area(i,j)*flux(i,j)
  END DO
END DO
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE flux_diag
END MODULE flux_diag_mod

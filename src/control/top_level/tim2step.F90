! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: TIM2STEP -------------------------------------------------
!
!    Purpose: Converts from an integer number of elapsed whole days and
!             seconds since the model basis time to elapsed timesteps,
!             given a general definition of the submodel timestep.
!             Forms a service routine for model date/time and internal
!             clock purposes, written for 32-bit portability.
!
!    Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!    Logical components covered: S620
!
!    Project task: S62
!
!    External documentation: On-line UM document C0 - The top-level
!                            control system
!
!    -------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

SUBROUTINE tim2step(elapsed_days,elapsed_secs,                    &
                    steps_in_period,secs_in_period,               &
                    elapsed_steps)
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!
INTEGER ::                                                        &
     elapsed_days,                                                &
                             ! IN  - elapsed days since ref time
     elapsed_secs,                                                &
                             ! IN  - elapsed secs in part of day
!                                  !       or days since ref time
           steps_in_period,                                             &
                                   ! IN  - steps in period defining T/S
           secs_in_period,                                              &
                                   ! IN  - secs  in period defining T/S
           elapsed_steps           ! OUT - elapsed steps since ref time
! Local Parameters
INTEGER :: secs_in_day        ! no. of seconds in 1 day
PARAMETER(secs_in_day=3600*24)
! Local Scalars
INTEGER :: factor             ! ratio of integers is exact factor

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TIM2STEP'
! ----------------------------------------------------------------------
!  1. Perform integer arithmetic to compute elapsed steps from elapsed
!     days/seconds and timestep definition.  The subroutine assumes
!     that SECS_IN_PERIOD is a whole number multiple of a day, or a
!     whole number divisor of a day, but does not check explicitly.
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (secs_in_period >= secs_in_day) THEN
  factor       = secs_in_period/secs_in_day ! no. days in period
  elapsed_steps =                                                &
     steps_in_period*(elapsed_days/factor) +                     &
     (((elapsed_days-(elapsed_days/factor)*factor)*secs_in_day + &
       elapsed_secs)*steps_in_period)/secs_in_period
ELSE          ! period is less than 1 day
  factor       = secs_in_day/secs_in_period ! no. periods in day
  elapsed_steps = elapsed_days*steps_in_period*factor +          &
                 elapsed_secs/(secs_in_period/steps_in_period)
END IF
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE tim2step

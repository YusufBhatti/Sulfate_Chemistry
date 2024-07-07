! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: TIME_DF --------------------------------------------------
!
!    Purpose: Subroutine to obtain a new model time in days and seconds
!             from some reference point, given an increment in days and
!             seconds.  Note that the seconds and days increments are
!             treated independently so that -ve increments or seconds
!             increments larger than the no of seconds in a day are
!             handled correctly.
!             Forms a service routine for model date/time and internal
!             clock purposes, written for 32-bit portability.
!
!    Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!    Logical components covered:
!
!    Project task:
!
!    External documentation: On-line UM document C0 - The top-level
!                            control system
!
!    -------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level
SUBROUTINE time_df(days1,secs1,del_days,del_secs,days2,secs2)
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!
INTEGER ::                                                        &
     days1,secs1,                                                 &
                             ! IN  - days/seconds (input time)
     del_days,del_secs,                                           &
                             ! IN  - days/seconds increments
     days2,secs2             ! OUT - days/seconds (output time)
!
INTEGER ::                                                        &
     secs_in_day            ! No of seconds in a day

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TIME_DF'
PARAMETER                                                         &
    (secs_in_day=24*60*60)
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
days2 = days1 + del_days + del_secs/secs_in_day
secs2 = secs1 + MOD(del_secs,secs_in_day)
!
IF (secs2 <  0) THEN
  secs2 = secs2 + secs_in_day
  days2 = days2 - 1
END IF
!
IF (secs2 >= secs_in_day) THEN
  secs2 = secs2 - secs_in_day
  days2 = days2 + 1
END IF
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE time_df

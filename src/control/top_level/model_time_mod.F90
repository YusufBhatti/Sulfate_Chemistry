! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Derived model time/step information including start/end
!              step numbers and frequencies (in steps) of interface field
!              generation, boundary field updating, ancillary field
!              updating; and assimilation start/end times.
!              NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
!              Also contains current time/date information, current
!              step number (echoed in history file) and steps-per-group.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!             This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90

MODULE model_time_mod

USE submodel_mod, ONLY: internal_id_max
USE control_max_sizes, ONLY: max_n_intf_a

IMPLICIT NONE 

INTEGER :: i_year               ! Current model time (years)
INTEGER :: i_month              ! Current model time (months)
INTEGER :: i_day                ! Current model time (days)
INTEGER :: i_hour               ! Current model time (hours)
INTEGER :: i_minute             ! Current model time (minutes)
INTEGER :: i_second             ! Current model time (seconds)
INTEGER :: i_day_number         ! Current model time (day no)
INTEGER :: previous_time(7)     ! Model time at previous step
INTEGER :: iau_dtresetstep      ! Data time reset step for IAU run

INTEGER :: basis_time_days  ! Integral no of days to basis time
INTEGER :: basis_time_secs  ! No of seconds-in-day at basis time

LOGICAL :: l_c360dy ! Use 360 day calendar

REAL    :: forecast_hrs         ! Hours since Data Time (ie T+nn)
REAL    :: data_minus_basis_hrs ! Data time - basis time (hours)

INTEGER :: stepim(internal_id_max)  ! Step no since basis time
INTEGER :: groupim(internal_id_max) ! Number of steps per group

! Finish step number this run
INTEGER :: target_end_stepim(internal_id_max)

! Timestep length in secs
REAL    :: secs_per_stepim(internal_id_max) 

! Frequency of interface field generation in steps
INTEGER :: interface_stepsim(max_n_intf_a,internal_id_max)

! Start steps for interface field generation
INTEGER :: interface_fstepim(max_n_intf_a,internal_id_max)

! End steps for interface field generation
INTEGER :: interface_lstepim(max_n_intf_a,internal_id_max)

! Frequency of  updating boundary fields in steps
INTEGER :: boundary_stepsim(internal_id_max)

! No of steps from boundary data prior to basis time to model basis time
INTEGER :: bndary_offsetim(internal_id_max)

! Lowest frequency for updating of ancillary fields in steps
INTEGER :: ancillary_stepsim(internal_id_max)

! Start steps for assimilation
INTEGER :: assim_firststepim(internal_id_max)

! Number of assimilation steps to analysis
INTEGER :: assim_stepsim(internal_id_max)

! Number of assimilation steps after analysis
INTEGER :: assim_extrastepsim(internal_id_max)

END MODULE model_time_mod

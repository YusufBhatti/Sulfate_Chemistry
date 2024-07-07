#if defined(C97_3A) || defined(RECON) || defined(VOMEXT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE timer_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TIMER_MOD'

CONTAINS
!
!   SUBROUTINE TIMER ------------------------------------------------
!
!                      Purpose:
!   Allows the recording of time spent in any section of the program
!   Two types of timings are supported:
!   non-inclusive : if a timed section of code (1) contains another
!                   timed section of code (2), then the timing for
!                   section (1) will not include the time spent in
!                   section (2). This is the normal use for the timer
!                   routine in the UM up to vn3.4
!   inclusive     : allows the user to measure the time taken between
!                   any two points in the code, irrespective of any
!                   other calls to the timer routine within the timed
!                   section
!
!   NB: Non-inclusive timers DO INCLUDE any inclusive timer sections
!       contained within them. If this section of code should not be
!       included, then also time it with a non-inclusive timer
!
!   Timer now also records the time spent in itself
!   Parameters:
!   section_name - 30 byte character string containing name of
!                  timer reference
!
!   action:
!    1 -> first call to timer (timer initialisation)
!    2 -> last call to timer (prints out the collected data)
!    3 -> non-inclusive start timer
!    4 -> non-inclusive end timer
!    5 -> inclusive start timer
!    6 -> inclusive end timer
!
!   Timer should be called with action=1 before the first executable
!   statement, and with action=2 after the last executable statement.
!
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Timer

SUBROUTINE timer(section_name,action_arg)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE UM_ParVars
USE UM_ParCore, ONLY: nproc
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE get_cpu_time_mod, ONLY: get_cpu_time

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
USE timer_output_mod, ONLY: timer_output
 
IMPLICIT NONE


! Arguments:
CHARACTER(LEN=*)    ::  section_name  ! reference name for timed
                                      ! section
INTEGER, INTENT(IN) ::  action_arg    ! what action to take


! Local variables:

INTEGER :: action             ! Local copy of action_arg
INTEGER :: section_ref        ! reference of the current section
INTEGER :: info               ! gcom return code
INTEGER :: i                  ! looper

REAL :: cpu_time_into_timer        ! cpu time at which timer routine entered
REAL :: wallclock_time_into_timer  ! wallclock time at which timer routine
                                   ! entered

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'TIMER'

INTEGER, PARAMETER :: max_timers = 300 ! maximum number of timings to be handled

! Magic numbers (action types):
INTEGER, PARAMETER :: first_call_to_timer  = 1
INTEGER, PARAMETER :: last_call_to_timer   = 2
INTEGER, PARAMETER :: non_incl_start_timer = 3
INTEGER, PARAMETER :: non_incl_end_timer   = 4
INTEGER, PARAMETER :: incl_start_timer     = 5
INTEGER, PARAMETER :: incl_end_timer       = 6
INTEGER, PARAMETER :: intermediate_output  = 7

! Saved/static data
! ni prefix = non-inclusive timings
! in prefix = inclusive timings

! names of timer references
CHARACTER(LEN=30), SAVE ::             &
                 ni_timer_name(max_timers) = '                             '
CHARACTER(LEN=30), SAVE ::             &
                 in_timer_name(max_timers) = '                             '

! number of times that a section of code has been timed
INTEGER, SAVE :: ni_number_of_times_timed(max_timers) = 0
INTEGER, SAVE :: in_number_of_times_timed(max_timers) = 0

! the reference of the timer stopped when a new one is started
INTEGER, SAVE :: old_timer(max_timers) = 0

! for non-inclusive timer - current section of code being timed
INTEGER, SAVE :: current_timer = 0

! number of timers currently known about
INTEGER, SAVE :: ni_number_of_timers = 0
INTEGER, SAVE :: in_number_of_timers = 0

! total amount of cpu time spent in a section of code
REAL, SAVE :: ni_cpu_time_elapsed(max_timers) = 0.0
REAL, SAVE :: in_cpu_time_elapsed(max_timers) = 0.0

! total amount of wallclock time
REAL, SAVE :: ni_wallclock_time_elapsed(max_timers) = 0.0
REAL, SAVE :: in_wallclock_time_elapsed(max_timers) = 0.0

! cpu time of starting
REAL, SAVE :: ni_cpu_time_started = 0.0
REAL, SAVE :: in_cpu_time_started(max_timers) = 0.0

! wallclock time of starting
REAL, SAVE :: ni_wallclock_time_started = 0.0
REAL, SAVE :: in_wallclock_time_started(max_timers) = 0.0

! is a particular timer running?
LOGICAL, SAVE :: in_timer_running(max_timers) = .FALSE.

! is the timer enabled
LOGICAL, SAVE :: timer_on = .FALSE.


! External timing functions

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
action = action_arg   ! allows local modifications of the argument

IF (( action  ==  first_call_to_timer) .OR.                                    &
    ( action  ==  non_incl_start_timer) .OR.                                   &
    ( action  ==  incl_start_timer)) THEN
  CALL gc_gsync(nproc,info)
END IF

IF (action  >   100) action=action-100
! The following line is useful for general debugging purposes
! It prints out the name of every timed routine on entry and exit
!         WRITE(umMessage,'(A,I2)') section_name,' action= ',action
!         CALL umPrint(umMessage,src='timer-timer3a')

! start up the timer timer

cpu_time_into_timer = get_cpu_time()
wallclock_time_into_timer = get_wallclock_time()
in_number_of_times_timed(1) = in_number_of_times_timed(1) + 1

! check the length of the section_name

IF (LEN(section_name)  >   30) THEN
  WRITE(umMessage,'(A)') 'TIMER has detected a non-fatal ERROR'
  CALL umPrint(umMessage,src='timer-timer3a')
  WRITE(umMessage,'(3A)') 'Section name ',section_name,' is too long.'
  CALL umPrint(umMessage,src='timer-timer3a')
  WRITE(umMessage,'(A)') 'Maximum of 30 characters is allowed'
  CALL umPrint(umMessage,src='timer-timer3a')
  WRITE(umMessage,'(2A)') section_name,' will be truncated to 30 chars.'
  CALL umPrint(umMessage,src='timer-timer3a')
END IF


! diagnose what action to take:

IF (action  ==  first_call_to_timer) THEN

  ! First call to timer - do initialisation

  DO i=1,max_timers
    ni_timer_name(i)            = '                              '
    in_timer_name(i)            = '                              '

    ni_number_of_times_timed(i) = 0
    in_number_of_times_timed(i) = 0

    ni_cpu_time_elapsed(i)      = 0.0
    in_cpu_time_elapsed(i)      = 0.0
    ni_wallclock_time_elapsed(i)= 0.0
    in_wallclock_time_elapsed(i)= 0.0

    in_timer_running(i)=.FALSE.
  END DO

  timer_on = .TRUE.

  current_timer = 1
  ni_number_of_timers = 1
  in_number_of_timers = 1
  in_timer_name(1) = 'TIMER'

  ! and start the timer running

  ni_cpu_time_started = get_cpu_time()
  ni_wallclock_time_started = get_wallclock_time()
  ni_number_of_times_timed(current_timer) = 1
  old_timer(current_timer) = 0
  ni_timer_name(current_timer) = section_name

  ! ----------------------------------------------------------------------

ELSE IF (timer_on .AND.                                                        &
  ( (action  ==  last_call_to_timer) .OR.                                      &
    (action  ==  intermediate_output) ) ) THEN

  ! Last call to timer - or intermediate output required, so
  ! print out table of results

  IF (action  ==  last_call_to_timer) THEN
    ! the only active timer should be no.1

    IF (current_timer  /=  1) THEN
      WRITE(umMessage,'(A)') 'TIMER has detected an ERROR'
      CALL umPrint(umMessage,src='timer-timer3a')
      WRITE(umMessage,'(2A)') 'Attempted to print results without switching ', &
               'off all running non-inclusive timers.'
      CALL umPrint(umMessage,src='timer-timer3a')
      WRITE(umMessage,'(A)') '** TIMER SWITCHED OFF FOR REST OF RUN **'
      CALL umPrint(umMessage,src='timer-timer3a')
      timer_on = .FALSE.
      GO TO 9999
    END IF

    ! Make sure there are no inclusive timers still running

    section_ref = 0
    DO i=1,in_number_of_timers
      IF (in_timer_running(i)) section_ref = i
    END DO

    IF (section_ref  /= 0) THEN
      WRITE(umMessage,'(A)') 'TIMER has detected an ERROR'
      CALL umPrint(umMessage,src='timer-timer3a')
      WRITE(umMessage,'(2A)') 'Attempted to print results without switching ', &
               'off all running inclusive timers.'
      CALL umPrint(umMessage,src='timer-timer3a')
      WRITE(umMessage,'(A)') '** TIMER SWITCHED OFF FOR REST OF RUN **'
      CALL umPrint(umMessage,src='timer-timer3a')
      timer_on = .FALSE.
      GO TO 9999
    END IF

    ! Just to make sure that timer isn't called again
    timer_on = .FALSE.

    ! and switch off the top level non-inclusive timer

    ni_cpu_time_elapsed(current_timer) =                                       &
        ni_cpu_time_elapsed(current_timer) +                                   &
        get_cpu_time() - ni_cpu_time_started

    ni_wallclock_time_elapsed(current_timer) =                                 &
        ni_wallclock_time_elapsed(current_timer) +                             &
        get_wallclock_time() - ni_wallclock_time_started

  END IF ! If this is the final call to timer

  CALL timer_output(                                                           &
    in_number_of_timers, ni_number_of_timers,                                  &
    in_cpu_time_elapsed, ni_cpu_time_elapsed,                                  &
    in_wallclock_time_elapsed, ni_wallclock_time_elapsed,                      &
    in_number_of_times_timed, ni_number_of_times_timed,                        &
    in_timer_name, ni_timer_name,                                              &
    action,section_name)

  ! ----------------------------------------------------------------------

ELSE IF (timer_on .AND.                                                        &
  (action  ==  non_incl_start_timer) ) THEN

  ! Start a non-inclusive timer running

  ! Switch off the current timer
  ni_cpu_time_elapsed(current_timer) =                                         &
    ni_cpu_time_elapsed(current_timer) +  get_cpu_time() - ni_cpu_time_started

  ni_wallclock_time_elapsed(current_timer) =                                   &
    ni_wallclock_time_elapsed(current_timer) +                                 &
    get_wallclock_time() - ni_wallclock_time_started

  ! See if we're already keeping records for this section

  section_ref = 0
  DO i=1,ni_number_of_timers
    IF (ni_timer_name(i)  ==  section_name) section_ref = i
  END DO

  ! Check to make sure that there is no timer already running for
  ! this section
  ! (NB an inclusive timer with the same reference name is allowed
  !  to run simultaneously with this one)

  IF (section_ref  ==  current_timer) THEN
    WRITE(umMessage,'(A)') 'TIMER has detected an ERROR'
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(3A)') 'Simultaneous non-inclusive timers attempted ',    &
               'for section ',section_name
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(A)') '** TIMER SWITCHED OFF FOR REST OF RUN **'
    CALL umPrint(umMessage,src='timer-timer3a')
    timer_on = .FALSE.
    GO TO 9999
  END IF

  ! calculate the section reference for the new timer

  IF (section_ref  ==  0) THEN
    ! this is a new section
    section_ref = ni_number_of_timers+1
    ni_timer_name(section_ref) = section_name
    ni_number_of_timers = section_ref
  END IF

  ! check that max_timers isn't exceeded:
  IF (ni_number_of_timers  >   max_timers) THEN
    WRITE(umMessage,'(A)') 'TIMER has detected an ERROR'
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(A,I4,2A)') 'More than ',max_timers,' non-inclusive ',    &
               'timers is not allowed.'
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(A)') '** TIMER SWITCHED OFF FOR REST OF RUN **'
    CALL umPrint(umMessage,src='timer-timer3a')
    timer_on = .FALSE.
    GO TO 9999
  END IF

  ! set up old_timer so that when this new timer is stopped, the
  ! current timer (that we've just stopped) can be restarted

  old_timer(section_ref)=current_timer

  ! now start up the new timer

  current_timer = section_ref
  ni_number_of_times_timed(current_timer) =                                    &
  ni_number_of_times_timed(current_timer) + 1
  ni_cpu_time_started = get_cpu_time()
  ni_wallclock_time_started = get_wallclock_time()

  ! ----------------------------------------------------------------------

ELSE IF (timer_on .AND.                                                        &
  (action  ==  non_incl_end_timer) ) THEN

  ! Stop a non-inclusive timer

  ! Make sure that we're being asked to end a timer that's actually
  ! running.

  IF (ni_timer_name(current_timer)  /=  section_name) THEN
    WRITE(umMessage,'(A)') 'TIMER has detected an ERROR'
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(3A)') 'Attempted to stop a non-active ',                 &
               'non-inclusive timer ',section_name
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(A)') '** TIMER SWITCHED OFF FOR REST OF RUN **'
    CALL umPrint(umMessage,src='timer-timer3a')
    timer_on = .FALSE.
    GO TO 9999
  END IF

  ! OK - so stop this timer:

  ni_cpu_time_elapsed(current_timer) =                                         &
    ni_cpu_time_elapsed(current_timer) + get_cpu_time() - ni_cpu_time_started

  ni_wallclock_time_elapsed(current_timer) =                                   &
    ni_wallclock_time_elapsed(current_timer) +                                 &
    get_wallclock_time() - ni_wallclock_time_started

  ! and now restart the old timer (ie. the one that was in
  ! operation at the time this one was started)

  IF (old_timer(current_timer)  ==  0) THEN
    ! this means I have just stopped the top level timer - there
    ! are no more to stop. This is an error - I should do this
    ! by calling the timer with action=2
    WRITE(umMessage,'(A)') 'TIMER has detected an ERROR'
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(A)') 'The top-level timer has been stopped'
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(A)') '** TIMER SWITCHED OFF FOR REST OF RUN **'
    CALL umPrint(umMessage,src='timer-timer3a')
    timer_on = .FALSE.
    GO TO 9999
  END IF

  current_timer=old_timer(current_timer)
  ni_cpu_time_started=get_cpu_time()
  ni_wallclock_time_started=get_wallclock_time()

  ! ----------------------------------------------------------------------

ELSE IF (timer_on .AND.                                                        &
  (action  ==  incl_start_timer) ) THEN

  ! Start an inclusive timer running

  ! See if we're already keeping records for this section

  section_ref = 0
  DO i=1,in_number_of_timers
    IF (in_timer_name(i)  ==  section_name) section_ref = i
  END DO

  ! and calculate the section reference

  IF (section_ref  ==  0) THEN
    !       this is a new one
    section_ref = in_number_of_timers + 1
    in_timer_name(section_ref) = section_name
    in_number_of_timers = section_ref
  END IF

  ! check that max_timers isn't exceeded:

  IF (in_number_of_timers  >   max_timers) THEN
    WRITE(umMessage,'(A)') 'TIMER has detected an ERROR'
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(A,I4,2A)') 'More than ',max_timers,' inclusive ',        &
               'timers is not allowed.'
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(A)') '** TIMER SWITCHED OFF FOR REST OF RUN **'
    CALL umPrint(umMessage,src='timer-timer3a')
    timer_on = .FALSE.
    GO TO 9999
  END IF

  ! Check to make sure that there is no timer already running for
  ! this section
  ! (NB a non-inclusive timer with the same reference name is allowed
  !  to run simultaneously with this one)

  IF (in_timer_running(section_ref)) THEN
    WRITE(umMessage,'(A)') 'TIMER has detected an ERROR'
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(2A)') 'Inclusive timer already running for ',            &
               section_name
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(A)') '** TIMER SWITCHED OFF FOR REST OF RUN **'
    CALL umPrint(umMessage,src='timer-timer3a')
    timer_on = .FALSE.
    GO TO 9999
  END IF

  ! so now we can start the timer for this section
  in_number_of_times_timed(section_ref) =                                      &
  in_number_of_times_timed(section_ref) + 1
  in_timer_running(section_ref) = .TRUE.
  in_cpu_time_started(section_ref) = get_cpu_time()
  in_wallclock_time_started(section_ref) = get_wallclock_time()

  ! ----------------------------------------------------------------------

ELSE IF (timer_on .AND.                                                        &
  (action  ==  incl_end_timer) ) THEN

  ! Stop an inclusive timer

  ! Find out what the reference number of this timer is

  section_ref = 0
  DO i=1,in_number_of_timers
    IF (in_timer_name(i)  ==  section_name) section_ref = i
  END DO

  IF (section_ref  ==  0) THEN
    WRITE(umMessage,'(A)') 'TIMER has detected an ERROR'
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(3A)') 'Attempting to stop a non-existent ',              &
               'inclusive timer ',section_name
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(A)') '** TIMER SWITCHED OFF FOR REST OF RUN **'
    CALL umPrint(umMessage,src='timer-timer3a')
    timer_on = .FALSE.
    GO TO 9999
  END IF

  ! Make sure this timer is actually running at the moment

  IF (.NOT. in_timer_running(section_ref)) THEN
    WRITE(umMessage,'(A)') 'TIMER has detected an ERROR'
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(3A)') 'Attempting to stop a non-running ',               &
             'inclusive timer ',section_name
    CALL umPrint(umMessage,src='timer-timer3a')
    WRITE(umMessage,'(A)') '** TIMER SWITCHED OFF FOR REST OF RUN **'
    CALL umPrint(umMessage,src='timer-timer3a')
    timer_on = .FALSE.
    GO TO 9999
  END IF

  ! now we can stop it
  in_cpu_time_elapsed(section_ref) =                                           &
   in_cpu_time_elapsed(section_ref) +                                          &
   get_cpu_time() - in_cpu_time_started(section_ref)

  in_wallclock_time_elapsed(section_ref) =                                     &
   in_wallclock_time_elapsed(section_ref) +                                    &
   get_wallclock_time() - in_wallclock_time_started(section_ref)

  in_timer_running(section_ref) = .FALSE.

  ! ----------------------------------------------------------------------

ELSE IF (timer_on) THEN

  WRITE(umMessage,'(A)') 'TIMER has detected an ERROR'
  CALL umPrint(umMessage,src='timer-timer3a')
  WRITE(umMessage,'(A,I2,3A)') 'incorrect action= ',action,' supplied by ',    &
             'section ',section_name
  CALL umPrint(umMessage,src='timer-timer3a')
  WRITE(umMessage,'(A)') 'Non-fatal error - TIMER will continue'
  CALL umPrint(umMessage,src='timer-timer3a')

END IF

9999 CONTINUE

! stop the timer timer
in_cpu_time_elapsed(1) = in_cpu_time_elapsed(1) +                            &
  get_cpu_time() - cpu_time_into_timer
in_wallclock_time_elapsed(1) = in_wallclock_time_elapsed(1) +                &
  get_wallclock_time() - wallclock_time_into_timer

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE timer

END MODULE timer_mod

#endif

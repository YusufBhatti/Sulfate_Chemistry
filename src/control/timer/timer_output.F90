! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE timer_output_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TIMER_OUTPUT_MOD'

CONTAINS
!
!   SUBROUTINE TIMER_OUTPUT -----------------------------------------
!
!                      Purpose:
!
!   Reports timing information from timer.

! ********************************************************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Timer

SUBROUTINE timer_output(                                                       &
  in_number_of_timers, ni_number_of_timers,                                    &
  in_cpu_time_elapsed, ni_cpu_time_elapsed,                                    &
  in_wallclock_time_elapsed, ni_wallclock_time_elapsed,                        &
  in_number_of_times_timed, ni_number_of_times_timed,                          &
  in_timer_name, ni_timer_name,                                                &
  action,message)

USE IOS_Common, ONLY:                                                          &
    L_IOS_Active,                                                              &
    io_servers,                                                                &
    l_io_server,                                                               &
    global_procs,                                                              &
    IOS_Server_Groups,                                                         &
    IOS_tasks_per_server,                                                      &
    model_procs
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE nlstcall_mod, ONLY: ltimer, ltimers_user, lstashdumptimer
USE coupling_control_mod, ONLY: l_oasis_timers

!$ USE omp_lib                 ! Note OpenMP sentinel

USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParCore, ONLY: mype, nproc, nproc_max
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength

#if !defined (RECON_SERIAL)
USE mpl, ONLY: mpl_real
#endif

IMPLICIT NONE


! Arguments

INTEGER ::                                                                     &
  in_number_of_timers                                                          &
                       ! IN number of inclusive timers
, ni_number_of_timers                                                          &
                       ! IN number of non-inclusive timers
, in_number_of_times_timed(in_number_of_timers)                                &
!                            ! IN number of times timed - inclusive
, ni_number_of_times_timed(ni_number_of_timers)                                &
!                            ! IN number of times timed - non-incl.
, action  ! final output or intermediate

REAL ::                                                                        &
  in_cpu_time_elapsed(in_number_of_timers)                                     &
!                            ! IN elapsed inclusive CPU time
, ni_cpu_time_elapsed(ni_number_of_timers)                                     &
!                            ! IN elapsed non-inclusive CPU time
, in_wallclock_time_elapsed(in_number_of_timers)                               &
!                            ! IN elapsed inclusive wallclock time
, ni_wallclock_time_elapsed(ni_number_of_timers)

CHARACTER(LEN=*) ::                                                            &
  in_timer_name(in_number_of_timers)                                           &
!                            ! IN name of timed section - inclusive
, ni_timer_name(ni_number_of_timers)
!                            ! IN name of timed section - non-incl.

CHARACTER(LEN=*) ::                                                            &
   message                   ! IN message to print


! Local variables

INTEGER, PARAMETER :: max_timers = 300

INTEGER, PARAMETER :: last_call_to_timer = 2
INTEGER, PARAMETER ::intermediate_output = 7

INTEGER, PARAMETER :: idx_wallclock = 1        ! Index for wallclock time
INTEGER, PARAMETER :: idx_cpu       = 2        ! Index for cpu time

INTEGER ::                                                                     &
  number_of_timers                                                             &
, local_number_of_times_timed(max_timers)

REAL ::                                                                        &
  local_cpu_time_elapsed(max_timers)                                           &
, local_wallclock_time_elapsed(max_timers)

INTEGER :: timer_name_length ! Length of timer name
INTEGER :: pe_count_digits   ! Number of digits in PE count

! Lengths of various spacings and column widths used for formatting
CHARACTER(LEN=5) :: name_col_space
CHARACTER(LEN=1) :: n_col_space
CHARACTER(LEN=1) :: pe_width
CHARACTER(LEN=1) :: pe_space
CHARACTER(LEN=1) :: pe_col_space

CHARACTER(LEN=LEN(in_timer_name(1))) :: local_timer_name(max_timers)

! Variables required for using intermediate timers
! They record the values on the last call to this routine
INTEGER :: last_in_number_of_times_timed(max_timers) = 0 
INTEGER :: last_ni_number_of_times_timed(max_timers) = 0

REAL :: last_in_cpu_time_elapsed(max_timers)       = 0.0
REAL :: last_ni_cpu_time_elapsed(max_timers)       = 0.0
REAL :: last_in_wallclock_time_elapsed(max_timers) = 0.0
REAL :: last_ni_wallclock_time_elapsed(max_timers) = 0.0

LOGICAL :: first_intermediate_timer_call = .TRUE.


INTEGER :: sortwork_int    ! work variable for sort
REAL :: sortwork_real   ! work variable for sort
CHARACTER(LEN=LEN(in_timer_name(1))) :: sortwork_char ! work variable for sort

REAL :: total_cpu_time,                                                        &
                             ! total cpu time spent in program
     total_wallclock_time,                                                     &
                             ! total wallclock time spent in
!                                  ! program
     average_cpu_elapsed,                                                      &
                             ! average cpu elapsed time
     average_wallclock_elapsed,                                                &
                                ! average wallclock elapsed time
     percent_of_cpu_total,                                                     &
                             ! % of cpu time spent in a section
     percent_of_wallclock_total,                                               &
                            ! % of wallclock time spent in a
!                                 ! section
     speed_up               ! speed_up=cpu/wallclock

! These are the declarations for MPP timer

INTEGER :: info,                                                               &
  wallclock_max_pe(max_timers),wallclock_min_pe(max_timers),                   &
  cpu_max_pe(max_timers),cpu_min_pe(max_timers)

REAL :: wallclock_mean(max_timers),cpu_mean(max_timers),                       &
     wallclock_median(max_timers),cpu_median(max_timers),                      &
     wallclock_sd(max_timers),cpu_sd(max_timers),                              &
     wallclock_max(max_timers),wallclock_min(max_timers),                      &
     cpu_max(max_timers),cpu_min(max_timers),                                  &
     cpu_total(max_timers)

INTEGER ::                                                                     &
  summ_n_timers                                                                &
                   ! number of routines for ni summary
, routine_id          ! routine id on this processor

INTEGER :: my_comm ! communicator

REAL :: wallclock_times(0:nproc_max-1)   ! wallclock time from each proc
REAL :: cpu_times(0:nproc_max-1)         ! cpu time from each proc
REAL :: my_time(2)                       ! temporary for comms
REAL :: my_times(2,0:nproc_max-1)        ! temporary for comms
REAL :: total_cpu                        ! total cpu time
REAL :: max_wall                         ! maxumum wallclock times

! names of sections
CHARACTER(LEN=LEN(in_timer_name(1))) :: summ_section(max_timers) 

! Variables for loops etc.
INTEGER :: i,j,k,timer_kind

! to get environment variable OMP_NUM_THREADS
CHARACTER(LEN=8) :: env_num_threads 

! Internal error message
CHARACTER(LEN=errormessagelength) :: cmessage    

! to store value read from env_num_threads
INTEGER :: atm_num_threads, omp_max_threads  
INTEGER :: icode
CHARACTER(LEN=*) :: RoutineName
PARAMETER (RoutineName = 'TIMER_OUTPUT')

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Check to see if this is an intermediate output, and the first
! time it has been called
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm,info)

IF (L_IO_Server) THEN
  ! IOS makes no use of timer.
  WRITE(umMessage,'(A,A)')'timer_output: This IO Server is not ',              &
      'reporting timing information'
  CALL umPrint(umMessage,src='timer_output')
ELSE
#if !defined(UTILIO) 
! If ltimer is true then we will use the timer information from all
! processors to get maxima and minima and do some statistical analysis.
! ni variables are non-inclusive values, in variable are inclusive. If
! l_oasis_timers or lstashdumptimer is true then the same process will be used 
! for the timings in the OASIS / IO routines. These have been separated out 
! to allow for unbiased load balancing and performance measurements.
  IF(ltimers_user .OR. ltimer .OR. l_oasis_timers .OR. lstashdumptimer) THEN
#endif
    IF ((action  ==  intermediate_output) .AND.                                &
        (first_intermediate_timer_call  ) ) THEN
      ! Copy the arguments into the last_* arrays
      first_intermediate_timer_call=.FALSE.
  
      DO i=1,in_number_of_timers
        last_in_number_of_times_timed(i)  = in_number_of_times_timed(i)
        last_in_cpu_time_elapsed(i)       = in_cpu_time_elapsed(i)
        last_in_wallclock_time_elapsed(i) = in_wallclock_time_elapsed(i)
      END DO
  
      DO i=1,ni_number_of_timers
        last_ni_number_of_times_timed(i)  = ni_number_of_times_timed(i)
        last_ni_cpu_time_elapsed(i)       = ni_cpu_time_elapsed(i)
        last_ni_wallclock_time_elapsed(i) = ni_wallclock_time_elapsed(i)
      END DO
  
      GO TO 9999  ! jump to end - no output on first call
    END IF

    ! Calculate amount of space to leave after heading in routine name column
    timer_name_length = LEN(in_timer_name(1))
    name_col_space = ""                                  ! 7 is the length
    WRITE(name_col_space, "(I0)") timer_name_length - 7  ! of "ROUTINE"

    ! Calculate the number of spaces to leave in the index/N column (this
    ! is just the width of the digits of the largest timer number)
    DO i = 1,3
      IF (MAX(ni_number_of_timers, in_number_of_timers) / 10**i == 0) EXIT
    END DO
    WRITE(n_col_space, "(I0)") i

    ! Find out the space width required for the PE number columns
    ! assume a maximum width of 9 since the characters below are LEN 1
    ! and a 9-digit PE count is suitably far off enough for now!
    DO i = 1,9
      IF (nproc / 10**i == 0) EXIT
    END DO
    pe_count_digits = i
    ! The space to leave before the pe column header should be one less than
    ! the width of the pe count itself, except for when the pe count width is
    ! less than 3, where it should always be 1
    WRITE(pe_col_space, "(I0)") MAX(pe_count_digits -2, 0) + 1
    ! The space to leave before the pe entry should always be 1 except in the
    ! case of single digit pe counts, where it should be 2
    WRITE(pe_space, "(I0)") MAX(pe_count_digits, 2)/pe_count_digits
    ! Finally the width of the pe number itself, which is just the value
    WRITE(pe_width, "(I0)") pe_count_digits
  
    WRITE(umMessage,'(A)')
    CALL umPrint(umMessage,src='timer_output')
    WRITE(umMessage,'(A)') '******************************************'
    CALL umPrint(umMessage,src='timer_output')
    WRITE(umMessage,'(A)')
    CALL umPrint(umMessage,src='timer_output')
  
    DO timer_kind=1,2  ! 1 is non-inclusive and 2 is inclusive
      ! Copy arguments into local arrays
      IF (action  ==  last_call_to_timer) THEN
        WRITE(umMessage,'(A)') 'END OF RUN - TIMER OUTPUT'
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A)') 'Timer information is for whole run'
        CALL umPrint(umMessage,src='timer_output')
        IF (timer_kind  ==  1) THEN  ! non-inclusive timer
          number_of_timers=ni_number_of_timers
          DO i=1,number_of_timers
            local_timer_name(i)             = ni_timer_name(i)
            local_cpu_time_elapsed(i)       = ni_cpu_time_elapsed(i)
            local_wallclock_time_elapsed(i) = ni_wallclock_time_elapsed(i)
            local_number_of_times_timed(i)  = ni_number_of_times_timed(i)
          END DO
        ELSE ! timer_kind  ==  2 - inclusive timer
          number_of_timers=in_number_of_timers
          DO i=1,number_of_timers
            local_timer_name(i)             = in_timer_name(i)
            local_cpu_time_elapsed(i)       = in_cpu_time_elapsed(i)
            local_wallclock_time_elapsed(i) = in_wallclock_time_elapsed(i)
            local_number_of_times_timed(i)  = in_number_of_times_timed(i)
          END DO
        END IF ! which timer kind this was
      ELSE  ! this is an intermediate output call
        WRITE(umMessage,'(A,A)') 'INTERMEDIATE TIMER OUTPUT :',message
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A,A)') 'Timer information is only for code ',        &
                   'executed since last intermediate timer output.'
        CALL umPrint(umMessage,src='timer_output')
        IF (timer_kind  ==  1) THEN  ! non-inclusive timer
          number_of_timers=ni_number_of_timers
          DO i=1,number_of_timers
            local_timer_name(i)=ni_timer_name(i)
            local_cpu_time_elapsed(i)=ni_cpu_time_elapsed(i)-                  &
                                      last_ni_cpu_time_elapsed(i)
            local_wallclock_time_elapsed(i)=                                   &
              ni_wallclock_time_elapsed(i) - last_ni_wallclock_time_elapsed(i)
            local_number_of_times_timed(i)=                                    &
              ni_number_of_times_timed(i) - last_ni_number_of_times_timed(i)
          END DO
        ELSE ! timer kind  ==  2 - inclusive timer
          number_of_timers=in_number_of_timers
          DO i=1,number_of_timers
            local_timer_name(i)=in_timer_name(i)
            local_cpu_time_elapsed(i)=in_cpu_time_elapsed(i)-                  &
                                      last_in_cpu_time_elapsed(i)
            local_wallclock_time_elapsed(i)=                                   &
              in_wallclock_time_elapsed(i) - last_in_wallclock_time_elapsed(i)
            local_number_of_times_timed(i)=                                    &
              in_number_of_times_timed(i) - last_in_number_of_times_timed(i)
          END DO
        END IF  ! what timer type
      END IF  ! what action to perform
  
      ! Do work for non-inclusive timers
  
      ! Calculate the total time in the program (based on non-inclusive
      ! timers)
      IF (timer_kind  ==  1) THEN
        total_cpu_time = 0.0
        total_wallclock_time = 0.0
        DO i=1,number_of_timers
          total_cpu_time = total_cpu_time + local_cpu_time_elapsed(i)
          total_wallclock_time =                                               &
            total_wallclock_time + local_wallclock_time_elapsed(i)
        END DO
  
        WRITE(umMessage,'(A,I0,A,F12.2)') 'PE ', mype,                         &
                   ' Elapsed CPU Time:       ', total_cpu_time
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A,I0,A,F12.2)') 'PE ', mype,                         &
                   ' Elapsed Wallclock Time: ', total_wallclock_time
        CALL umPrint(umMessage,src='timer_output')
  
        ! Calculate the total cpu time over all processors and the
        ! maximum elapsed time - so allowing a speedup to be caclulated
  
        total_cpu=total_cpu_time
        max_wall=total_wallclock_time
  
        CALL gc_rsum(1,nproc,info,total_cpu)
        CALL gc_rmax(1,nproc,info,max_wall)
  
        max_wall=MAX(max_wall,0.000001)
        WRITE(umMessage,'(A)')
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A,F12.2)') 'Total Elapsed CPU Time:         ',       &
          total_cpu
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A,F12.2)') 'Maximum Elapsed Wallclock Time: ',       &
          max_wall
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A,F12.2)')  'Speedup:                        ',      &
          total_cpu/max_wall
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A)') '--------------------------------------------'
        CALL umPrint(umMessage,src='timer_output')
  
      END IF
  
      ! Sort subroutines into time order (based on wallclock time)
  
      DO i=1,number_of_timers-1
        DO j=(i+1),number_of_timers
          IF (local_wallclock_time_elapsed(j)  >                               &
              local_wallclock_time_elapsed(i)) THEN
  
            !             Swap the two entries
            sortwork_real = local_cpu_time_elapsed(i)
            local_cpu_time_elapsed(i) = local_cpu_time_elapsed(j)
            local_cpu_time_elapsed(j) = sortwork_real
  
            sortwork_real = local_wallclock_time_elapsed(i)
            local_wallclock_time_elapsed(i) = local_wallclock_time_elapsed(j)
            local_wallclock_time_elapsed(j) = sortwork_real
  
            sortwork_int = local_number_of_times_timed(i)
            local_number_of_times_timed(i) = local_number_of_times_timed(j)
            local_number_of_times_timed(j) = sortwork_int
  
            sortwork_char = local_timer_name(i)
            local_timer_name(i) = local_timer_name(j)
            local_timer_name(j) = sortwork_char
  
          END IF
        END DO
      END DO

      IF (timer_kind  ==  1) THEN

        WRITE(umMessage,'(20X,A,I0)') 'Non-Inclusive Timer Summary for PE ',mype
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A,'//n_col_space//'X'//                              &
                         ',A,'//TRIM(name_col_space)//'X'//                    &
                         ',8X,A,7X,A,7X,A,6X,A,7X,A,3X,A,2X,A,2X,A)')          &
          'N','ROUTINE','CALLS','TOT CPU','AVERAGE','TOT WALL','AVERAGE',      &
          '% CPU','% WALL','SPEED-UP'
        CALL umPrint(umMessage,src='timer_output')
      ELSE

        WRITE(umMessage,'(20X,A,i0)') 'Inclusive Timer Summary for PE ',mype
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A,'//n_col_space//'X'//                              &
                         ',A,'//TRIM(name_col_space)//'X'//                    &
                         ',8X,A,7X,A,7X,A,6X,A,7X,A,2X,A)')                    &
          'N','ROUTINE','CALLS','TOT CPU','AVERAGE','TOT WALL',                &
          'AVERAGE','SPEED-UP'
        CALL umPrint(umMessage,src='timer_output')
      END IF
  
      DO i=1,number_of_timers
        IF (local_number_of_times_timed(i)  /=  0) THEN
          average_cpu_elapsed =  local_cpu_time_elapsed(i)/                    &
                                 local_number_of_times_timed(i)
          average_wallclock_elapsed = local_wallclock_time_elapsed(i)/         &
                                      local_number_of_times_timed(i)
        ELSE
          average_cpu_elapsed = 0.0
          average_wallclock_elapsed = 0.0
        END IF
  
        IF (local_wallclock_time_elapsed(i)  >   0) THEN
          speed_up=local_cpu_time_elapsed(i) / local_wallclock_time_elapsed(i)
        ELSE
          speed_up=1.0
        END IF
  
        IF (timer_kind  ==  1) THEN  ! non-inclusive timer has some
          !                                      ! extra output
  
          percent_of_cpu_total = 100.0*local_cpu_time_elapsed(i)/              &
                                 total_cpu_time
          percent_of_wallclock_total =                                         &
            100.0*local_wallclock_time_elapsed(i) / total_wallclock_time
  
          WRITE(umMessage,'(I'//n_col_space//'.'//n_col_space//                &
                          ',1X,A,1X,I12,4(2X,F12.2),2(2X,F6.2),4X,F6.2)')      &
                      i,local_timer_name(i),                                   &
                      local_number_of_times_timed(i),                          &
                      local_cpu_time_elapsed(i),average_cpu_elapsed,           &
                      local_wallclock_time_elapsed(i),                         &
                      average_wallclock_elapsed,                               &
                      percent_of_cpu_total,                                    &
                      percent_of_wallclock_total,speed_up
          CALL umPrint(umMessage,src='timer_output')
  
        ELSE ! inclusive timer has slightly less to output
  
          WRITE(umMessage,'(I'//n_col_space//'.'//n_col_space//                &
                          ',1X,A,1X,I12,4(2X,F12.2),4X,F6.2)')                 &
                      i,local_timer_name(i),                                   &
                      local_number_of_times_timed(i),                          &
                      local_cpu_time_elapsed(i),average_cpu_elapsed,           &
                      local_wallclock_time_elapsed(i),                         &
                      average_wallclock_elapsed,speed_up
          CALL umPrint(umMessage,src='timer_output')
  
        END IF
  
      END DO
  
  
#if !defined(UTILIO) && !defined(UTILHIST)
  
      ! We only want to process the statistics of the timings for real mpp
      ! jobs where nproc > 1.  Any utilities where UTILIO and UTILHIST
      ! are defined for the cpp use nproc = 1 as parameter in parvars.h.
  
      ! And now to assemble an overall timing assesment on PE0
      ! Each PE sends it total wallclock and cpu time spent in each routine
      ! to PE0, which calculates the average, s.d., max and min, and
      ! sorts on the basis of the average wallclock time
  
  
      ! We'll use the list of routines that PE0 already has as the master
      ! list.
  
      IF (mype  ==  0) THEN
  
        CALL model_setup_timer_output()

        summ_n_timers=number_of_timers
        DO i=1,summ_n_timers
          summ_section(i)=local_timer_name(i)
        END DO
      END IF
  
      ! tell everyone else how many routines to do summary on - and which
      ! routines they are
      CALL gc_ibcast(3213,1,0,nproc,info,summ_n_timers)
      CALL gc_cbcast(3214,20*summ_n_timers,0,nproc,info,summ_section)
  
  
      DO i=1,summ_n_timers
  
        ! which section_ref is this for me?
  
        routine_id=0
        DO j=1,number_of_timers
          IF (local_timer_name(j)  ==  summ_section(i))                        &
            routine_id=j
        END DO
  
        IF (routine_id  >   0) THEN
          my_time(idx_wallclock) = local_wallclock_time_elapsed(routine_id)
          my_time(idx_cpu)       = local_cpu_time_elapsed(routine_id)
        ELSE
          my_time(idx_wallclock) = 0.0
          my_time(idx_cpu)       = 0.0
        END IF
  
#if defined(RECON_SERIAL)
        wallclock_times = my_time(idx_wallclock)
        cpu_times       = my_time(idx_cpu)
#else
        ! send my information to PE 0. Copy in and out of temporary space
        ! to use just one communication
        CALL mpl_gather(my_time,  2,       mpl_real,                           &
                        my_times, 2,       mpl_real,                           &
                        0,        my_comm, info)
        wallclock_times = my_times(idx_wallclock,:)
        cpu_times       = my_times(idx_cpu,:)
#endif
  
        ! Use a barrier to ensure that multiple gathers (when there are 
        ! multiple timers) don't interfere with each other.
        CALL gc_gsync(nproc,info)
  
        IF (mype  ==  0) THEN
          ! collect all the information - and start calculating the statistics
          wallclock_mean(i)=0.0
          cpu_total(i)=0.0
          wallclock_max(i)=-1.0e30
          wallclock_min(i)=1.0e30
          cpu_max(i)=-1.0e30
          cpu_min(i)=1.0e30
  
          DO j=0,nproc-1
  
            wallclock_mean(i)=wallclock_mean(i)+wallclock_times(j)
            cpu_total(i)=cpu_total(i)+cpu_times(j)
  
            IF (wallclock_times(j) >  wallclock_max(i)) THEN
              wallclock_max(i)=wallclock_times(j)
              wallclock_max_pe(i)=j
            END IF
            IF (wallclock_times(j) <  wallclock_min(i)) THEN
              wallclock_min(i)=wallclock_times(j)
              wallclock_min_pe(i)=j
            END IF
            IF (cpu_times(j) >  cpu_max(i)) THEN
              cpu_max(i)=cpu_times(j)
              cpu_max_pe(i)=j
            END IF
            IF (cpu_times(j) <  cpu_min(i)) THEN
              cpu_min(i)=cpu_times(j)
              cpu_min_pe(i)=j
            END IF
  
          END DO ! loop over processors
  
          ! and calculate the statistics
          ! first calculate the means
          wallclock_mean(i)=wallclock_mean(i)/nproc
          cpu_mean(i)=cpu_total(i)/nproc
          ! To stop a divide by zero later:
          IF (wallclock_mean(i)  ==  0.0) wallclock_mean(i)=1.0e-20
          IF (cpu_mean(i)  ==  0.0) cpu_mean(i)=1.0e-20
          ! and now the standard deviation
          wallclock_sd(i)=0.0
          cpu_sd(i)=0.0
          DO j=0,nproc-1
            wallclock_sd(i)=wallclock_sd(i)+                                   &
              (wallclock_times(j)-wallclock_mean(i))*                          &
              (wallclock_times(j)-wallclock_mean(i))
            cpu_sd(i)=cpu_sd(i)+(cpu_times(j)-cpu_mean(i))*                    &
                          (cpu_times(j)-cpu_mean(i))
          END DO
          wallclock_sd(i)=SQRT(wallclock_sd(i)/nproc)
          cpu_sd(i)=SQRT(cpu_sd(i)/nproc)
  
          ! Calculate the median
          DO j=0,nproc-2
            DO k=j+1,nproc-1
              IF (wallclock_times(k)  >   wallclock_times(j)) THEN
                sortwork_real=wallclock_times(j)
                wallclock_times(j)=wallclock_times(k)
                wallclock_times(k)=sortwork_real
              END IF
              IF (cpu_times(k)  >   cpu_times(j)) THEN
                sortwork_real=cpu_times(j)
                cpu_times(j)=cpu_times(k)
                cpu_times(k)=sortwork_real
              END IF
            END DO
          END DO
  
          IF (MOD(nproc,2)  ==  0) THEN
            wallclock_median(i)=(wallclock_times((nproc/2)-1)+                 &
                                 wallclock_times(nproc/2))*0.5
            cpu_median(i)=(cpu_times((nproc/2)-1)+                             &
                           cpu_times(nproc/2))*0.5
          ELSE
            wallclock_median(i)=wallclock_times(nproc/2)
            cpu_median(i)=cpu_times(nproc/2)
          END IF
  
        END IF ! am I PE 0?
  
      END DO ! loop over sections
  
      ! Sort and output the information on PE 0
  
      IF (mype  ==  0) THEN
  
        DO i=1,summ_n_timers-1
          DO j=(i+1),summ_n_timers
            IF (wallclock_max(j)  >   wallclock_max(i)) THEN
  
              ! Swap the entries I and J
  
              sortwork_char=summ_section(i)
              summ_section(i)=summ_section(j)
              summ_section(j)=sortwork_char
  
              sortwork_real=wallclock_mean(i)
              wallclock_mean(i)=wallclock_mean(j)
              wallclock_mean(j)=sortwork_real
  
              sortwork_real=wallclock_median(i)
              wallclock_median(i)=wallclock_median(j)
              wallclock_median(j)=sortwork_real
  
              sortwork_real=wallclock_sd(i)
              wallclock_sd(i)=wallclock_sd(j)
              wallclock_sd(j)=sortwork_real
  
              sortwork_real=wallclock_max(i)
              wallclock_max(i)=wallclock_max(j)
              wallclock_max(j)=sortwork_real
  
              sortwork_real=wallclock_min(i)
              wallclock_min(i)=wallclock_min(j)
              wallclock_min(j)=sortwork_real
  
              sortwork_int=wallclock_min_pe(i)
              wallclock_min_pe(i)=wallclock_min_pe(j)
              wallclock_min_pe(j)=sortwork_int
  
              sortwork_int=wallclock_max_pe(i)
              wallclock_max_pe(i)=wallclock_max_pe(j)
              wallclock_max_pe(j)=sortwork_int
  
              sortwork_real=cpu_mean(i)
              cpu_mean(i)=cpu_mean(j)
              cpu_mean(j)=sortwork_real
  
              sortwork_real=cpu_median(i)
              cpu_median(i)=cpu_median(j)
              cpu_median(j)=sortwork_real
  
              sortwork_real=cpu_sd(i)
              cpu_sd(i)=cpu_sd(j)
              cpu_sd(j)=sortwork_real
  
              sortwork_real=cpu_max(i)
              cpu_max(i)=cpu_max(j)
              cpu_max(j)=sortwork_real
  
              sortwork_real=cpu_min(i)
              cpu_min(i)=cpu_min(j)
              cpu_min(j)=sortwork_real
  
              sortwork_int=cpu_min_pe(i)
              cpu_min_pe(i)=cpu_min_pe(j)
              cpu_min_pe(j)=sortwork_int
  
              sortwork_int=cpu_max_pe(i)
              cpu_max_pe(i)=cpu_max_pe(j)
              cpu_max_pe(j)=sortwork_int
  
            END IF
          END DO
        END DO
  
        ! and write out the information
        WRITE(umMessage,'(A)')
        CALL umPrint(umMessage,src='timer_output')
        IF (timer_kind  ==  1) THEN
          WRITE(umMessage,'(A)') 'MPP : Non Inclusive timer summary'
          CALL umPrint(umMessage,src='timer_output')
        ELSE
          WRITE(umMessage,'(A)') 'MPP : Inclusive timer summary'
          CALL umPrint(umMessage,src='timer_output')
        END IF
  
        WRITE(umMessage,'(A)')
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A)')  'WALLCLOCK  TIMES'
        CALL umPrint(umMessage,src='timer_output')
  
        WRITE(umMessage,'(A,'//n_col_space//'X'//                              &
                        ',A,'//TRIM(name_col_space)//'X'//                     &
                        ',9X,A,7X,A,8X,A,3X,A'//                               &
                        ',2(10X,A,'//pe_col_space//'X,A))')                    &
          'N','ROUTINE','MEAN','MEDIAN','SD','% of mean','MAX','(PE)',         &
          'MIN','(PE)'
        CALL umPrint(umMessage,src='timer_output')

        DO i=1,summ_n_timers
  
          WRITE(umMessage,                                                     &
          '(I'//n_col_space//'.'//n_col_space//                                &
          ',1X,A,2(1X,F12.2),1X,F9.2,5X,F6.2,A'//                              &
          ',2(1X,F12.2,'//pe_space//'X,A,I'//pe_width//',A) )')                &
                     i,summ_section(i),                                        &
                     wallclock_mean(i),wallclock_median(i),                    &
                     wallclock_sd(i),                                          &
                     (wallclock_sd(i)/wallclock_mean(i))*100.0,'%',            &
                     wallclock_max(i),'(',wallclock_max_pe(i),')',             &
                     wallclock_min(i),'(',wallclock_min_pe(i),')'
          CALL umPrint(umMessage,src='timer_output')
        END DO
  
        WRITE(umMessage,'(A)')
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A)')  'CPU TIMES (sorted by wallclock times)'
        CALL umPrint(umMessage,src='timer_output')
  
        WRITE(umMessage,'(A,'//n_col_space//'X'//                              &
                        ',A,'//TRIM(name_col_space)//'X'//                     &
                        ',9X,A,7X,A,8X,A,3X,A'//                               &
                        ',2(10X,A,'//pe_col_space//'X,A))')                    &
          'N','ROUTINE','MEAN','MEDIAN','SD','% of mean','MAX','(PE)',         &
          'MIN','(PE)'
        CALL umPrint(umMessage,src='timer_output')
  
        DO i=1,summ_n_timers
          WRITE(umMessage,                                                     &
          '(I'//n_col_space//'.'//n_col_space//                                &
          ',1X,A,2(1X,F12.2),1X,F9.2,5X,F6.2,A'//                              &
          ',2(1X,F12.2,'//pe_space//'X,A,I'//pe_width//',A) )')                &
                     i,summ_section(i),                                        &
                     cpu_mean(i),cpu_median(i),                                &
                     cpu_sd(i),                                                &
                     (cpu_sd(i)/cpu_mean(i))*100.0,'%',                        &
                     cpu_max(i),'(',cpu_max_pe(i),')',                         &
                     cpu_min(i),'(',cpu_min_pe(i),')'
          CALL umPrint(umMessage,src='timer_output')
        END DO
  
      END IF
      WRITE(umMessage,'(A)')
      CALL umPrint(umMessage,src='timer_output')
  
#endif
  
    END DO ! loop over timer kind
  
    ! Finally copy the timer info into the last_* arrays so that the
    ! intermediate timer can calculate the timings since this point
  
    DO i=1,in_number_of_timers
      last_in_number_of_times_timed(i) = in_number_of_times_timed(i)
      last_in_cpu_time_elapsed(i)      = in_cpu_time_elapsed(i)
      last_in_wallclock_time_elapsed(i)= in_wallclock_time_elapsed(i)
    END DO
  
    DO i=1,ni_number_of_timers
      last_ni_number_of_times_timed(i) = ni_number_of_times_timed(i)
      last_ni_cpu_time_elapsed(i)      = ni_cpu_time_elapsed(i)
      last_ni_wallclock_time_elapsed(i)= ni_wallclock_time_elapsed(i)
    END DO

    9999 CONTINUE

#if !defined(UTILIO) 
  ELSE    ! ltimer=false
! If ltimer is false then we present the timing information for the 
! whole run on pe0 only. This involves no communication.

    IF (action  ==  last_call_to_timer) THEN
      IF (mype  ==  0) THEN

        CALL model_setup_timer_output()

        WRITE(umMessage,'(A)')
        CALL umPrint(umMessage,src='timer_output')

        WRITE(umMessage,'(A)') 'END OF RUN - TIMER OUTPUT'
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A)') 'Timer information is for whole run'
        CALL umPrint(umMessage,src='timer_output')

        total_cpu_time = 0.0
        total_wallclock_time = 0.0
        DO i=1,ni_number_of_timers
          total_cpu_time = total_cpu_time + ni_cpu_time_elapsed(i)   
          total_wallclock_time =                                             &
          total_wallclock_time + ni_wallclock_time_elapsed(i)    
        END DO

        WRITE(umMessage,'(A,I6,A,F12.2)') 'PE ', mype,                       &
               '  Elapsed CPU Time: ', total_cpu_time
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A,I6,A,F12.2)') 'PE ', mype,                       &
               '  Elapsed Wallclock Time: ', total_wallclock_time
        CALL umPrint(umMessage,src='timer_output')
        WRITE(umMessage,'(A)')
        CALL umPrint(umMessage,src='timer_output')

      END IF ! mype==0
    END IF   ! last timer action
  END IF !ltimer
#endif
END IF ! L_IO_Server

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

CONTAINS

SUBROUTINE model_setup_timer_output()
! Small internal subroutine which contains code that is called
! twice in the main routine. Describes model resources used.

USE get_env_var_mod, ONLY: get_env_var

IMPLICIT NONE

INTEGER :: length  ! Length of returned string

WRITE(umMessage,'(A)')
CALL umPrint(umMessage,src='timer_output')
CALL umPrint( 'MPP Timing information : ',src='timer_output')
WRITE(umMessage,'(I0,A,I0,A,I0)')nproc,                                &
     ' processors in atmosphere configuration ',                       &
     nproc_x,' x ',nproc_y
CALL umPrint(umMessage,src='timer_output')
IF (L_IOS_Active()) THEN
  WRITE(umMessage,'(I0,A,I0,A,I0,A)')SIZE(io_servers),                 &
       ' IO servers in configuration '                                 &
       ,IOS_Server_Groups,' x ',IOS_tasks_per_server,                  &
       ' are not included in stats'
  CALL umPrint(umMessage,src='timer_output')
ELSE IF (global_procs /= model_procs) THEN
  WRITE(umMessage,'(I0,A)')global_procs-model_procs,                   &
       ' processes are not included in stats (misconfigured IO servers?)'
  CALL umPrint(umMessage,src='timer_output')
END IF
! write out number of threads, warn user if thread number has been
! changed during run of um
CALL get_env_var('OMP_NUM_THREADS',env_num_threads,allow_missing=.TRUE., &
                 length=length)

! If OMP_NUM_THREADS not found, check for Intel's KMP_NUM_THREADS instead:
IF (length < 0) THEN
  CALL get_env_var('KMP_NUM_THREADS',env_num_threads,allow_missing=.TRUE.)
END IF

IF (length > 0) THEN
  
  READ(env_num_threads,'(I4)') atm_num_threads
  omp_max_threads = atm_num_threads
  
  ! Only want OpenMP section executing if OpenMP is compiled in,
  ! so protect by sentinal
!$  omp_max_threads = omp_get_max_threads()
  
  IF ( omp_max_threads /= atm_num_threads) THEN
!$  cmessage = 'Environment variable '//                              &
!$  'OMP_NUM_THREADS does not match current number of threads'
!$  ICODE = -100

!$  CALL ereport(routinename,ICODE,cmessage)

!$  WRITE(umMessage,'(A,I0)') 'Number of threads : ', omp_max_threads
!$  CALL umPrint(umMessage,src='timer_output')
!$  WRITE(umMessage,'(2A,I0)') 'Environment variable ',               &
!$       'OMP_NUM_THREADS : ', atm_num_threads
!$  CALL umPrint(umMessage,src='timer_output')
  ELSE
!$  WRITE(umMessage,'(A,I0)') 'Number of OMP threads : ',             &
!$        atm_num_threads
!$  CALL umPrint(umMessage,src='timer_output')
  END IF
END IF

END SUBROUTINE model_setup_timer_output

END SUBROUTINE timer_output

END MODULE timer_output_mod

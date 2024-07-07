! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc.
!
! Description
!   This module seeks out segment sizes which give the smallest elapsed time per
!   gridpoint. The algorithm is based on the Metropolis algorithm -- simulated
!   annealing.

MODULE autotune_mod

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

USE UM_ParCore, ONLY: mype

USE umPrintMgr, ONLY:   &
    umPrint,            &
    umMessage

USE tuning_segments_mod, ONLY:                        &
  l_autotune_verbose,                                 &
  autotune_trial_radius, autotune_granularity,        &
  autotune_cooling_fraction, autotune_cooling_steps,  &
  autotune_min_seg_size, autotune_max_seg_size,       &
  autotune_init_time_coeff, autotune_num_leaders

USE get_wallclock_time_mod, ONLY: get_wallclock_time

USE mpl

IMPLICIT NONE
PRIVATE

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER :: ModuleName='AUTOTUNE_MOD'

!DrHook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL, PARAMETER :: high_value = HUGE(1.0)
REAL, PARAMETER :: low_value  = TINY(1.0)

!-------------------------------------------------------------------------------
! Module types
!-------------------------------------------------------------------------------

TYPE :: autotune_type
  PRIVATE

  !Region name / tag
  CHARACTER(LEN=80) :: name
  CHARACTER(LEN=10) :: tag

  !Segment size handling
  INTEGER  :: min_seg_size
  INTEGER  :: max_seg_size
  INTEGER  :: granularity
  INTEGER  :: start_size
  INTEGER  :: trial_size
  INTEGER  :: current_size
  INTEGER  :: trial_radius

  !Elapsed time handling
  REAL     :: current_timing
  REAL     :: trial_timing
  REAL     :: start_time
  REAL     :: stop_time
  INTEGER  :: total_points

  !Acceptance ratio handling
  REAL     :: acc_ratio
  INTEGER  :: num_trials
  INTEGER  :: num_accepted
  REAL     :: uphill_acc_ratio
  INTEGER  :: num_uphill_trials
  INTEGER  :: num_uphill_accepted

  !Fictitious temperature handling
  REAL     :: temperature
  REAL     :: initial_temperature
  REAL     :: init_time_coeff
  REAL     :: temperature_fraction
  LOGICAL  :: temperature_initialised
  REAL     :: cooling_coeff
  REAL     :: cooling_fraction
  INTEGER  :: cooling_steps
  INTEGER  :: calls_per_step

  !Leader board handling
  INTEGER  :: consecutive_count
  INTEGER  :: num_leaders
  INTEGER, ALLOCATABLE  :: leader_size(:)
  INTEGER, ALLOCATABLE  :: leader_count(:)

  !Random seed
  INTEGER, ALLOCATABLE  :: rand_seed(:)

END TYPE

!-------------------------------------------------------------------------------
! Public entities
!-------------------------------------------------------------------------------

PUBLIC :: autotune_type
PUBLIC :: autotune_init
PUBLIC :: autotune_start_region
PUBLIC :: autotune_stop_region
PUBLIC :: autotune_advance
PUBLIC :: autotune_get_trial_size
PUBLIC :: autotune_report

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
CONTAINS

!-------------------------------------------------------------------------------
! SYNOPSIS
!   autotune_init(this, name, tag, start_size, calls_per_step)
!
! DESCRIPTION
!   Constructor for an autotuning state object.
!
! ARGUMENTS
!   this           -- The autotune state object to initialise.
!   name           -- A string unique to a particular segmented region.
!   tag            -- A short string unique to a segment region.
!   start_size     -- The initial setting of the segment size.
!   calls_per_step -- Optional. The number of calls made to the segmented region
!                     within a single model timestep.
!
! NOTES
!   Divided into four parts:
!     Part 1: Variables taken from the namelist.
!     Part 2: Variables set via constructor arguments.
!     Part 3: Miscellaneous initialisations.
!     Part 4: Initialisations requiring computation.
!-------------------------------------------------------------------------------

SUBROUTINE autotune_init(this, name, tag, start_size, calls_per_step)
IMPLICIT NONE

TYPE(autotune_type), INTENT(OUT) :: this
CHARACTER(LEN=*),    INTENT(IN)  :: name
CHARACTER(LEN=*),    INTENT(IN)  :: tag
INTEGER,             INTENT(IN)  :: start_size
INTEGER, OPTIONAL,   INTENT(IN)  :: calls_per_step

INTEGER              :: cooling_steps
INTEGER              :: seed_size
INTEGER, ALLOCATABLE :: local_rand_seed(:)

CHARACTER(LEN=*), PARAMETER :: RoutineName='AUTOTUNE_INIT'
REAL(KIND=jprb) :: zhook_handle

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!Part 1: Namelist variables

this%trial_radius     = autotune_trial_radius
this%init_time_coeff  = autotune_init_time_coeff
this%cooling_fraction = autotune_cooling_fraction
this%cooling_steps    = autotune_cooling_steps
this%num_leaders      = autotune_num_leaders
this%min_seg_size     = autotune_min_seg_size
this%max_seg_size     = autotune_max_seg_size
this%granularity      = autotune_granularity

!Part 2: Constructor arguments

this%name         = TRIM(name)
this%tag          = TRIM(tag)
this%start_size   = start_size

IF(PRESENT(calls_per_step)) THEN
  this%calls_per_step = calls_per_step
ELSE
  this%calls_per_step = 1
END IF

!Part 3: Miscellaneous initialisations

this%num_trials              = 0
this%num_accepted            = 0
this%num_uphill_trials       = 0
this%num_uphill_accepted     = 0
this%consecutive_count       = 0
this%acc_ratio               = 0.0
this%uphill_acc_ratio        = 0.0

this%total_points            = 0
this%trial_timing            = high_value
this%trial_size              = this%start_size
this%current_timing          = high_value
this%current_size            = this%start_size

this%temperature             = low_value
this%temperature_initialised = .FALSE.
this%initial_temperature     = 0.0
this%temperature_fraction    = 0.0

!Leader board.
ALLOCATE(this%leader_size (1:this%num_leaders))
ALLOCATE(this%leader_count(1:this%num_leaders))
this%leader_size(:)   = 0
this%leader_count(:)  = 0

!Random seed. Use a temporary allocatable in the call
!to RANDOM_SEED to avoid an internal compiler error with ifort 15.
CALL RANDOM_SEED(SIZE=seed_size)
ALLOCATE(this%rand_seed(1:seed_size))
ALLOCATE(local_rand_seed(1:seed_size))
CALL RANDOM_SEED(GET=local_rand_seed)
this%rand_seed(:) = local_rand_seed(:)
DEALLOCATE(local_rand_seed)

!Part 4: Initialisations requiring computation.

!Cooling coefficient
cooling_steps = this%cooling_steps * this%calls_per_step
this%cooling_coeff = this%cooling_fraction    &
                           ** (1.0 / REAL(cooling_steps))

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE autotune_init

!-------------------------------------------------------------------------------
! SYNOPSIS
!   CALL autotune_start_region(this)
!
! DESCRIPTION
!   Logs the start time for a timed region, local to each MPI task.
!   Also re-initialises the (local) stop time.
!
! ARGUMENTS
!   this -- The autotune state object.
!-------------------------------------------------------------------------------

SUBROUTINE autotune_start_region(this)
IMPLICIT NONE

TYPE(autotune_type), INTENT(INOUT) :: this

CHARACTER(LEN=*), PARAMETER :: RoutineName='AUTOTUNE_START_REGION'
REAL(KIND=jprb) :: zhook_handle

!Omit DrHook: small routine, and slightly reduces impact of DrHook overheads
!interfering with measured elapsed times.

this%start_time   = get_wallclock_time()
this%stop_time    = 0.0
this%trial_timing = 0.0

END SUBROUTINE autotune_start_region

!-------------------------------------------------------------------------------
! SYNOPSIS
!   CALL autotune_stop_region(this, total_points)
!
! DESCRIPTION
!   Logs the stop time for a timed region, local to each MPI task. Computes the
!   local elapsed time for the local MPI task, and then performs a reduction to
!   obtain the total elapsed time across all tasks. For this reason, the routine
!   must be called on all MPI tasks to prevent a model hang.
!
! ARGUMENTS
!   this         -- The autotune state object.
!   total_points -- The total number of gridpoints on which work was done.
!
! NOTES
!   This routine logs the "stop" time after the first DrHook call, and therefore
!   would include DrHook overheads from this routine. However, if DrHook was
!   switched on, then there would be other DrHook overheads from other sources.
!   In the context of an autotuning run, DrHook would be off. The benefit of
!   DrHook in this module leans more towards debugging/tracebacks rather than
!   profiling. For this reason, we keep all executable code inside the callipers
!   for consistency.
!-------------------------------------------------------------------------------

SUBROUTINE autotune_stop_region(this, local_points)
IMPLICIT NONE

TYPE(autotune_type), INTENT(INOUT) :: this
INTEGER,             INTENT(IN)    :: local_points

REAL    :: local_timing
INTEGER :: my_comm
INTEGER :: ierr

CHARACTER(LEN=*), PARAMETER :: RoutineName='AUTOTUNE_STOP_REGION'
REAL(KIND=jprb) :: zhook_handle

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

this%stop_time   = get_wallclock_time()
local_timing     = this%stop_time - this%start_time
this%trial_timing = 0.0
this%total_points = 0

CALL gc_get_communicator(my_comm, ierr)
CALL mpl_reduce(local_timing, this%trial_timing,   &
                1, MPL_REAL, MPL_SUM, 0, my_comm, ierr)
CALL mpl_reduce(local_points, this%total_points,   &
                1, MPL_INTEGER, MPL_SUM, 0, my_comm, ierr)

IF (this%total_points > 0) THEN
  this%trial_timing = this%trial_timing   &
        / REAL(this%total_points)
END IF

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE autotune_stop_region

!-------------------------------------------------------------------------------
! SYNOPSIS
!   CALL autotune_advance(this)
!
! DESCRIPTION
!   Performs a Monte Carlo trial, using information (e.g. timings) stored in the
!   state object. Performs acceptance ratio calculations and housekeeping.
!   Calls other routines to ensure that the leader board is updated. Applies
!   cooling to the fictitious temperature.
!
!   Initialisation of the fictitious temperature requires an initial timing for
!   the segmented region. Hence the fictitious temperature is initialised on the
!   first pass through this routine, and not in the constructor; a valid elapsed
!   time must be available.
!
! ARGUMENTS
!   this -- The autotune state object.
!-------------------------------------------------------------------------------

SUBROUTINE autotune_advance(this)
IMPLICIT NONE

TYPE(autotune_type), INTENT(INOUT) :: this

REAL :: exponential_comparator
REAL :: rand_num
REAL :: real_trial_radius

CHARACTER(LEN=*), PARAMETER :: RoutineName='AUTOTUNE_ADVANCE'
REAL(KIND=jprb) :: zhook_handle

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!If the fictitious temperature is not yet initialised, determine an initial
!fictitious temperature which would give some user-defined proportional
!increase in execution time a ballpark 50% chance of success.
IF (.NOT. this%temperature_initialised) THEN
  this%initial_temperature =                &
    ((this%init_time_coeff-1.0)*this%trial_timing) / LOG(2.0)
  this%temperature = this%initial_temperature
  this%temperature_initialised = .TRUE.
END IF

!Metropolis procedure proper:

exponential_comparator = EXP( -(this%trial_timing-this%current_timing) &
                              / this%temperature )
rand_num = autotune_get_random_number(this)

!Increment total number of trials, regardless.
this%num_trials = this%num_trials + 1

!Timing is smaller than before; always accept.
IF(this%trial_timing <= this%current_timing) THEN

  !Accept
  this%num_accepted = this%num_accepted + 1
  CALL autotune_accept(this)

!Timing is greater than before; accept with some probability.
ELSE IF ( rand_num < exponential_comparator ) THEN

  !Accept
  this%num_accepted        = this%num_accepted        + 1
  this%num_uphill_accepted = this%num_uphill_accepted + 1
  this%num_uphill_trials   = this%num_uphill_trials   + 1
  CALL autotune_accept(this)

!Uphill move not accepted. Increment uphill trial counter only.
ELSE

  this%num_uphill_trials = this%num_uphill_trials + 1

END IF

!Keep a record of the number of times this segment size has occurred
!concurrently. If not in the leader board yet, it might be in the future.
this%consecutive_count = this%consecutive_count + 1

!Update the leader board.
CALL autotune_update_leaders(this)

!Acceptance ratios. There may not have been any uphill trials yet; need to
!prevent potential division by zero.
this%acc_ratio = REAL(this%num_accepted)  &
               / REAL(this%num_trials)

this%uphill_acc_ratio = REAL(this%num_uphill_accepted) &
                      / REAL(MAX(this%num_uphill_trials,1))

!Apply cooling
this%temperature = this%temperature * this%cooling_coeff
this%temperature_fraction = this%temperature / this%initial_temperature

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE autotune_advance

!-------------------------------------------------------------------------------
! SYNOPSIS
!   CALL autotune_accept(this)
!
! DESCRIPTION
!   Performs necessary actions for acceptance of a segment size change. Note
!   that counters needed for acceptance ratios are not in scope, because they
!   depend on the precise path to acceptance.
!
!   Specific functions:
!     * The trial size and its associated timing become the current size and
!       timing.
!     * If the trial size is not the same as the current size, then reset
!       consecutive_count which records the number of consecutive times that the
!       current segment size has been the current segment size. Catch the case
!       where a lower timing has occurred for the same segment size (noise); in
!       this case, do not reset the consecutive_counter.
!
! ARGUMENTS
!   this -- The autotune state object.
!-------------------------------------------------------------------------------

SUBROUTINE autotune_accept(this)
  IMPLICIT NONE

  !Arguments
  TYPE(autotune_type), INTENT(INOUT) :: this

  !Omit DrHook: small routine.

  !Zero if trial segment size is a new candidate.
  IF (this%trial_size /= this%current_size) THEN
    this%consecutive_count = 0
  END IF

  !Accept
  this%current_size   = this%trial_size
  this%current_timing = this%trial_timing

END SUBROUTINE autotune_accept

!-------------------------------------------------------------------------------
! SYNOPSIS
!   CALL autotune_get_trial_size(this)

! DESCRIPTION
!   Changes the segment size in either direction based on a random number.
!   Segment sizes are conditioned to have the granularity specified via the
!   namelist, and are constrained to exist between a specified (non-zero)
!   minimum and maximum.  The min/max, granularity, and the trial radius
!   (maximum move in either direction) all enter via the state object.
!
!   The segment size (currently) needs to be the same on all MPI tasks. This
!   routine performs the necessary broadcast. Corollary: this routine must be
!   called on all MPI tasks to prevent a model hang.
!
!   Under one condition, this routine does not adjust the segment size and
!   simply returns the original starting value: when the segment size is not
!   finite. In this latter case, it may be that the segment sizes have been
!   controlled via the number of segments, rather than the segment size. Since
!   this module acts on segment sizes, not segment numbers, this case is ignored
!   and the segment size is not adjusted.
!
! ARGUMENTS
!   this -- The autotune state object.
!
! RETURNS
!   autotune_get_trial_size -- The new segment size, adjusted or not.
!
! SIDE EFFECTS
!   This is not a pure function: the trial size from the state object is updated
!   on each call, and has the same value as that returned by the function.
!-------------------------------------------------------------------------------

FUNCTION autotune_get_trial_size(this)
IMPLICIT NONE

INTEGER                            :: autotune_get_trial_size
TYPE(autotune_type), INTENT(INOUT) :: this

REAL     :: rand_num
INTEGER  :: trial_size_diff
INTEGER  :: new_size

REAL     :: real_trial_radius
INTEGER  :: new_trial_radius
INTEGER  :: my_comm
INTEGER  :: ierr

CHARACTER(LEN=*), PARAMETER :: RoutineName='AUTOTUNE_GET_TRIAL_SIZE'
REAL(KIND=jprb) :: zhook_handle

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (this%start_size > 0) THEN

  !Proceed and generate a trial segment size
  real_trial_radius = REAL(this%trial_radius)
  rand_num = autotune_get_random_number(this)
  trial_size_diff = NINT((2.0*rand_num-1.0)*real_trial_radius)

  new_size = this%current_size + trial_size_diff

  new_size = MAX(new_size, this%min_seg_size, this%granularity)
  new_size = MIN(new_size, this%max_seg_size)

  new_size = new_size - MOD(new_size, this%granularity)

  CALL gc_get_communicator(my_comm, ierr)
  CALL mpl_bcast(new_size, 1, MPL_INTEGER, 0, my_comm, ierr)

  this%trial_size = new_size
  autotune_get_trial_size = new_size

ELSE

  this%trial_size = this%start_size
  autotune_get_trial_size = this%start_size

END IF

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION autotune_get_trial_size

!-------------------------------------------------------------------------------
! SYNOPSIS
!   CALL autotune_report(state, quiet)
!
! DESCRIPTION
!   Writes autotuning information to stdout. The output can be quietened via the
!   `quiet` argument; this way, we can change the semantics of what `quiet`
!   means at a later date. For the moment, it means silence.
!
!   The report verbosity is controlled by l_autotune_verbose. If false (not
!   verbose), then the region name, current segment size and leader board are
!   written regardless.
!
! ARGUMENTS
!   state -- The autotune state object.
!   quiet -- Optional. If present and true, output is suppressed.
!-------------------------------------------------------------------------------

SUBROUTINE autotune_report(this, quiet)
IMPLICIT NONE

TYPE(autotune_type), INTENT(IN) :: this
LOGICAL, OPTIONAL,   INTENT(IN) :: quiet

LOGICAL :: quiet_now
INTEGER :: i

CHARACTER(LEN=*), PARAMETER :: format_head    = '(a,">",1x,a)'
CHARACTER(LEN=*), PARAMETER :: format_name    = '(a,">",1x,a25,1x,a)'
CHARACTER(LEN=*), PARAMETER :: format_f       = '(a,">",1x,a25,1x,f10.5)'
CHARACTER(LEN=*), PARAMETER :: format_es      = '(a,">",1x,a25,1x,es14.5)'
CHARACTER(LEN=*), PARAMETER :: format_integer = '(a,">",1x,a25,1x,i0)'
CHARACTER(LEN=*), PARAMETER :: format_integer_dual &
                                              = '(a,">",1x,a25,2(1x,i0))'
CHARACTER(LEN=*), PARAMETER :: format_integer_es &
                                       = '(a,">",1x,a25,1x,i0,1x,es14.4)'

!Leader board formatting
CHARACTER(LEN=*), PARAMETER :: format_leader        &
                               = '(a,">",1x,i6,2x,i12,2x,i0)'
CHARACTER(LEN=*), PARAMETER :: format_leader_head   &
                               = '(a,">",1x,a6,2x,a12,2x,a5)'
CHARACTER(LEN=*), PARAMETER :: format_leader_dashes &
                               = '(a,">",1x,6("-"),2x,12("-"),2x,5("-"))'

CHARACTER(LEN=*), PARAMETER :: RoutineName='AUTOTUNE_REPORT'
REAL(KIND=jprb) :: zhook_handle

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

quiet_now = .FALSE.
IF(PRESENT( quiet)) quiet_now = quiet

IF (.NOT. quiet_now) THEN

  WRITE(umMessage, format_head) TRIM(this%tag), &
    '==> Segment size auto-tuning <=='
  CALL umPrint(umMessage, src='autotune_mod', pe=mype)

  WRITE(umMessage, format_name) TRIM(this%tag),     &
    'Region name:', TRIM(this%name)
  CALL umPrint(umMessage, src='autotune_mod', pe=mype)

  IF (l_autotune_verbose) THEN

    WRITE(umMessage, format_integer_dual) TRIM(this%tag),  &
      'Min/Max segment size:',                             &
       this%min_seg_size, this%max_seg_size
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)

    WRITE(umMessage, format_integer) TRIM(this%tag),  &
      'Trial radius:', this%trial_radius
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)

    WRITE(umMessage, format_integer) TRIM(this%tag),  &
      'Granularity:', this%granularity
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)

    WRITE(umMessage, format_es) TRIM(this%tag),       &
      'Cooling fraction:', this%cooling_fraction
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)

    WRITE(umMessage, format_integer) TRIM(this%tag),  &
      'Cooling steps:', this%cooling_steps
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)

    WRITE(umMessage, format_integer) TRIM(this%tag),  &
      'Calls per step:', this%calls_per_step
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)

    WRITE(umMessage, format_f) TRIM(this%tag),        &
      'Overall acceptance ratio:', this%acc_ratio
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)

    WRITE(umMessage, format_f) TRIM(this%tag),        &
      'Uphill acceptance ratio:', this%uphill_acc_ratio
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)

    WRITE(umMessage, format_es) TRIM(this%tag),       &
      'Temperature fraction:', this%temperature_fraction
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)

    WRITE(umMessage, format_es) TRIM(this%tag),       &
      'Current timing:', this%current_timing
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)

    WRITE(umMessage, format_integer) TRIM(this%tag),  &
      'Total num. points:', this%total_points
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)
  END IF

  WRITE(umMessage, format_integer_es) TRIM(this%tag), &
    'Trial size/timing:', this%trial_size, this%trial_timing
  CALL umPrint(umMessage, src='autotune_mod', pe=mype)

  WRITE(umMessage, format_integer) TRIM(this%tag),  &
    'Current segment size:', this%current_size
  CALL umPrint(umMessage, src='autotune_mod', pe=mype)

  !Write the leader board
  WRITE(umMessage, format_leader_dashes) TRIM(this%tag)
  CALL umPrint(umMessage, src='autotune_mod', pe=mype)
  WRITE(umMessage, format_leader_head) TRIM(this%tag), &
       'Leader', 'Segment size', 'Count'
  CALL umPrint(umMessage, src='autotune_mod', pe=mype)
  WRITE(umMessage, format_leader_dashes) TRIM(this%tag)
  CALL umPrint(umMessage, src='autotune_mod', pe=mype)
  WRITE(umMessage, format_leader_dashes) TRIM(this%tag)

  DO i = 1, this%num_leaders
    WRITE(umMessage, format_leader) TRIM(this%tag), &
         i,  this%leader_size(i), this%leader_count(i)
    CALL umPrint(umMessage, src='autotune_mod', pe=mype)
  END DO

  WRITE(umMessage, format_leader_dashes) TRIM(this%tag)
  CALL umPrint(umMessage, src='autotune_mod', pe=mype)

END IF !quiet_now

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE autotune_report

!-------------------------------------------------------------------------------
! SYNOPSIS
!   CALL autotune_update_leaders(this)
!
! DESCRIPTION
!   Maintains the leader board of best segment sizes. Sorts the entries such
!   that the most frequently visited segment size is at the top. This routine
!   does not write the leader board to the output file.
!
!   If the current segment size already has an entry in the table, then its
!   existing count is increased by 1.  If it does not, then this is a new entry
!   and the count is taken from consecutive_count.
!
! ARGUMENTS
!   this -- The autotune state object.
!-------------------------------------------------------------------------------

SUBROUTINE autotune_update_leaders(this)
IMPLICIT NONE

TYPE(autotune_type), INTENT(INOUT) :: this

INTEGER :: i
INTEGER :: pos

INTEGER :: tmp_count
INTEGER :: tmp_size
LOGICAL :: converged

CHARACTER(LEN=*), PARAMETER :: RoutineName='AUTOTUNE_UPDATE_LEADERS'
REAL(KIND=jprb) :: zhook_handle

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!Already in the table?
pos = 0
DO i = 1, this%num_leaders
  IF (this%leader_size(i) == this%current_size) THEN
    pos = i
  END IF
END DO

!Enter into the table
IF(pos>0) THEN
  this%leader_count(pos) = this%leader_count(pos) + 1
ELSE

  !Insert at the bottom, but only if the consecutive count exceeds the count
  !for the segment size currently at the bottom.
  IF (this%consecutive_count > this%leader_count(this%num_leaders)) THEN
    this%leader_count(this%num_leaders) = this%consecutive_count
    this%leader_size (this%num_leaders) = this%current_size
  END IF

END IF

!Sort table
converged = .FALSE.
DO WHILE (.NOT. converged)
  converged = .TRUE.

  DO i = 1, this%num_leaders-1
    IF ( this%leader_count(i+1) > this%leader_count(i) ) THEN
      converged = .FALSE.

      tmp_count = this%leader_count(i)
      tmp_size  = this%leader_size (i)

      this%leader_count(i) = this%leader_count(i+1)
      this%leader_size (i) = this%leader_size (i+1)

      this%leader_count(i+1) = tmp_count
      this%leader_size (i+1) = tmp_size
    END IF
  END DO
END DO

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE autotune_update_leaders

!-------------------------------------------------------------------------------
! SYNOPSIS
!   random_number = autotune_get_random_number(this)
!
! DESCRIPTION
!   Automatic tuning can interfere with model evolution, due to calls to
!   Fortran's RANDOM_NUMBER: the generator state is advanced more frequently,
!   leading to different random numbers in subsequent calls, in other parts of
!   the model.
!
!   This function wraps RANDOM_NUMBER in such a way as to transparently maintain
!   a generator state for this particular autotuning state, and leaves the state
!   unaltered as far as the rest of the code is concerned.
!
! ARGUMENTS
!   this -- The autotune state object.
!
! RETURNS
!   autotune_get_random_number -- The random number from RANDOM_NUMBER call.
!-------------------------------------------------------------------------------

FUNCTION autotune_get_random_number(this)
IMPLICIT NONE

TYPE(autotune_type), INTENT(INOUT) :: this

REAL    :: autotune_get_random_number
REAL    :: rand_num
INTEGER :: seed_size
INTEGER, ALLOCATABLE :: tmp_seed(:)
INTEGER, ALLOCATABLE :: local_rand_seed(:)

CHARACTER(LEN=*), PARAMETER :: RoutineName='AUTOTUNE_GET_RANDOM_NUMBER'
REAL(KIND=jprb) :: zhook_handle

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!Use local_rand_seed to avoid an internal compiler error with ifort 15.
CALL RANDOM_SEED(SIZE=seed_size)
ALLOCATE(tmp_seed(1:seed_size))
ALLOCATE(local_rand_seed(1:seed_size))
local_rand_seed(:) = this%rand_seed(:)

CALL RANDOM_SEED(GET=tmp_seed)
CALL RANDOM_SEED(PUT=local_rand_seed)

CALL RANDOM_NUMBER(rand_num)

CALL RANDOM_SEED(GET=local_rand_seed)
CALL RANDOM_SEED(PUT=tmp_seed)

this%rand_seed(:) = local_rand_seed(:)
DEALLOCATE(local_rand_seed)
DEALLOCATE(tmp_seed)

autotune_get_random_number = rand_num

IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION autotune_get_random_number

END MODULE autotune_mod


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------------
MODULE get_wallclock_time_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GET_WALLCLOCK_TIME_MOD'

CONTAINS
!  Gets the elapsed time from the system

! Function Interface:
REAL FUNCTION get_wallclock_time()

!$ USE omp_lib

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Description:
!   A system function is used to return the numbers of seconds which 
!   have elapsed. The exact routine used under the hood depends on 
!   cpp flags. These should be matched to the facilities available on
!   the platform used. Options include
!   
!  cpp                        routine used
!  ---                        ------------
!  IBM_XL_FORTRAN             rtc          - only available on AIX
!  MPL_WTIME                  mpl_wtime    - only available with MPI/MPL/GCOM
!  (default)                  system_clock - the Fortan intrinsic

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Timer

! Code Description:
!   Language: FORTRAN 90


!- End of header

! Local variables
INTEGER :: tid, tmax          ! threading variables

INTEGER, ALLOCATABLE, SAVE :: start_count(:)
INTEGER, ALLOCATABLE, SAVE :: old_count(:)
REAL, ALLOCATABLE, SAVE    :: rollover(:)
REAL, ALLOCATABLE, SAVE    :: oneover_count_rate(:)

INTEGER       :: my_count, count_rate, count_max, elapsed_count

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_WALLCLOCK_TIME'

#if defined(IBM_XL_FORTRAN)
! IBM Fortran offers a better clock so we don't use the intrinsic.
REAL(8) :: rtc

#elif defined(MPL_WTIME)
! MPI can give a good timer
REAL, EXTERNAL :: mpl_wtime

#else
! Intrinsic procedures called:
INTRINSIC SYSTEM_CLOCK

#endif


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if defined(IBM_XL_FORTRAN)
! On IBM we can use the high-resolution routine rtc which returns
! the number of seconds since the initial value of the machine's
! real-time clock. Nothing mentioned in manual about roll-over so
! are ignoring this aspect.
!
! Of course the 8 byte long float returned by this call is immediately cast
! down to the default float size which the function returns. At the moment
! we use "-r8" to make that 8 bytes but it need not be.

get_wallclock_time = rtc()

#elif defined(MPL_WTIME)
! mpl_wtime has a size cast appropriately by GCOM
! It returns a number seconds - no considerations for roll-over.
get_wallclock_time = mpl_wtime()

#else

tid = 0
!$ tid = omp_get_thread_num()
!$OMP CRITICAL(gwt_alloc)
IF (.NOT. ALLOCATED(start_count)) THEN
  tmax = 1
!$ tmax = omp_get_max_threads()
  ALLOCATE(start_count(0:tmax-1))
  start_count(:) = -1
  ALLOCATE(old_count(0:tmax-1))
  old_count(:) = 0
  ALLOCATE(rollover(0:tmax-1))
  rollover(:) = 0.0
  ALLOCATE(oneover_count_rate(0:tmax-1))
  oneover_count_rate(:) = 0.0
END IF
!$OMP END CRITICAL(gwt_alloc)

CALL SYSTEM_CLOCK(COUNT=my_count, count_rate=count_rate, &
                  count_max=count_max)

IF ((old_count(tid)  <   start_count(tid)) .AND.                       &
    ((my_count <   old_count(tid)) .OR.                                &
     (my_count >   start_count(tid)))) THEN
  IF (count_rate /= 0) THEN
    rollover(tid)=rollover(tid)+(REAL(count_max)/REAL(count_rate))
  END IF
END IF

IF (start_count(tid)  ==  -1) THEN
  start_count(tid) = my_count
  IF (count_rate /= 0) THEN
    oneover_count_rate(tid)=1.0/REAL(count_rate)
  ELSE
    oneover_count_rate(tid) = 1.0
  END IF
END IF

elapsed_count = my_count - start_count(tid)

IF (elapsed_count  <   0) elapsed_count = elapsed_count + count_max

get_wallclock_time = rollover(tid)+                                    &
                    (REAL(elapsed_count)*oneover_count_rate(tid))
old_count(tid) = my_count
#endif

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION get_wallclock_time
END MODULE get_wallclock_time_mod

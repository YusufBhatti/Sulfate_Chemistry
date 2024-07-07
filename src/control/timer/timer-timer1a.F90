#if defined(C97_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE timer_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TIMER_MOD'

CONTAINS

SUBROUTINE timer(sub,i_arg)
!   ..................................................................
!   SUBROUTINE TIMER
!   ----------------
!   A PROGRAM TO RECORD THE FLOW THROUGH AND TIME SPENT IN EACH
!   SUBROUTINE OF A MULTI SUBROUTINED PROGRAM.
!   CALLS TO TIMER MUST BE MADE BEFORE THE FIRST EXECUTABLE STATEMENT
!   OF EACH SUBROUTINE AND BEFORE EACH RETURN POINT.

!   PARAMETERS:
!   SUB - 8 BYTE CHARACTER STRING CONTAINING NAME OF CALLING
!         PROGRAM WHICH HAS BEEN LEFT ADJUSTED

!   I  -  I=1 - FIRST CALL FROM MAIN (OR THE HIGHEST LEVEL PROGRAM)
!         I=2 - LAST CALL FROM MAIN ( OR THE HIGHEST LEVEL PROGRAM)
!         I=3 - FIRST CALL FROM LOWER SUBROUTINE
!         I=4 CALL BEFORE RETURN STATEMENTS IN LOWER SUBROUTINE

!   ---NOTE :THE CRAY FACILITY PERFTRACE IS MORE APPROPRIATE FOR SINGLE
!   TASKED TIMING DIAGNOSTICS, BUT A BREAKDOWN OF ELAPSE TIME BY
!   SUBROUTINE OF MULTITASKED JOBS IS UNAVAILABLE WITHOUT THIS HOMEMADE
!   TIMER ROUTINE

!   REVISED TO CONFORM WITH UM SOFTWARE STANDARDS, THIS VERSION RUNS
!   ON BOTH THE CRAY & HP WORKSTATIONS. IF THERE ARE > 200 DIFFERENT
!   SUBROUTINE CALLS TO TIMER, FINAL TABLE REPLACED BY WARNING MESSAGE.
!   ..................................................................
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Timer
!  ---------------------------------------------------------------------



USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE get_cpu_time_mod, ONLY: get_cpu_time 
USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
IMPLICIT NONE


INTEGER, PARAMETER :: klimit=200
CHARACTER(LEN=8) :: sub,subname(klimit),returnam(klimit),swork
INTEGER :: i_arg    ! Argument to routine
INTEGER :: nenter(klimit),i,l,j,k,nsubs,istop,ncalls,iwork
REAL :: elapse(klimit),tote,elpend,elpstart,cpustart,tot,at,cpuend
REAL :: totelap,avelap,pcent,rwork,totlap,speedup,p,time(klimit)
REAL :: uwork
LOGICAL :: Found ! Used in DO loops, replacing GO TO
SAVE subname,returnam,nenter,elapse,time
SAVE k,j,nsubs,elpstart,istop,cpustart


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TIMER'
!      -----------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
i = i_arg   ! allow changes to input argument

IF (i  >   100) i=i-100
IF (i == 1) THEN

  !      First call to timer from the main program
  !      Set up initial values of variables

  k         = 1
  j         = 0
  istop     = 0
  nsubs     = 1
  subname(1)= sub

  DO l=1,klimit
    elapse(l) = 0.0
    time(l) = 0.0
    nenter(l) = 0

  END DO

  nenter(1)=1
  elpstart = get_wallclock_time()
  cpustart = get_cpu_time()

  !      -----------------------------------------------------------------
ELSE IF ((i == 2) .AND. (istop == 0)) THEN

  !      Last call to timer from main program
  !      Print out table of results

  cpuend = get_cpu_time()
  elpend = get_wallclock_time()
  elapse(1)=elapse(1)+(elpend-elpstart)
  time(1)=time(1)+cpuend-cpustart

  !        Calculate total time in program
  tote=0.0
  tot=0.0
  DO k=1,nsubs
    tote=tote+elapse(k)
    tot=tot+time(k)
  END DO

  !        Sort subroutines into time order

  DO k=1,(nsubs-1)

    DO j=(k+1),nsubs

      IF (time(j) >  time(k)) THEN

        !              Swap the values:
        rwork=time(k)
        time(k)=time(j)
        time(j)=rwork
        uwork=elapse(k)
        elapse(k)=elapse(j)
        elapse(j)=uwork
        iwork=nenter(k)
        nenter(k)=nenter(j)
        nenter(j)=iwork
        swork=subname(k)
        subname(k)=subname(j)
        subname(j)=swork

      END IF

    END DO

  END DO


  !        Output timing information
  WRITE(umMessage,'("1",20X," FLOW TRACE SUMMARY")')
  CALL umPrint(umMessage,src='timer-timer1a')
  CALL umPrint('',src='timer-timer1a')
  WRITE(umMessage,'(4x,a,6x,a,6x,a,3x,a,2x,a,4x,a,4x,a,6x,a,1x,a)')                    &
   'ROUTINE','CPU','%','CALLS','AVERAGE','ELAPSE','%','AVERAGE','CPU'
  CALL umPrint(umMessage,src='timer-timer1a')
  WRITE(umMessage,'(17x,a,4x,a,9x,a,4x,a,3x,a,4x,a,2x,a)')                             &
   'TIME','CPU','CPUTIME','TIME','ELAPSE','ELAPSE','SPEEDUP'
  CALL umPrint(umMessage,src='timer-timer1a')

  DO k=1,nsubs

    sub=subname(k)
    totlap=time(k)
    totelap=elapse(k)
    ncalls=nenter(k)
    avelap=totelap/ncalls
    p=100.0*totlap/tot
    pcent=100.0*totelap/tote
    at=time(k)/nenter(k)
    IF (avelap == 0.0) THEN
      speedup=0.0
    ELSE
      speedup=at/avelap
    END IF
    CALL umPrint('',src='timer-timer1a')
    WRITE(umMessage,'(T1,I3,T5,A8,T13,F10.4,T25,F5.2,T30,'// &
          'I5,T35,F10.4,T45,F10.4,T57,F5.2,T62,F10.4,T74,F5.2)')   &
          k,sub,totlap,p,nenter(k),at,totelap,pcent,avelap,speedup
    CALL umPrint(umMessage,src='timer-timer1a')

  END DO

  speedup=tot/tote
  CALL umPrint('',src='timer-timer1a')
  WRITE(umMessage,'(T3,''**TOTAL'',T12,F11.4,T44,F11.4,T74,F5.2)') &
      tot,tote,speedup
  CALL umPrint(umMessage,src='timer-timer1a')

  !      -----------------------------------------------------------------
ELSE IF ((i == 3) .AND. (istop == 0)) THEN

  !      First call in subroutine

  !        Switch off timer
  cpuend = get_cpu_time()
  elpend = get_wallclock_time()
  elapse(k)=elapse(k)+(elpend-elpstart)
  time(k)=time(k)+cpuend-cpustart

  !        Save name of calling subroutine
  j=j+1
  returnam(j)=subname(k)

  !        Check subroutine name
  Found=.FALSE.
  DO k=1,nsubs
    IF (subname(k) == sub) THEN
      Found=.TRUE.
      EXIT
    END IF
  END DO

  !        New subroutine entered
  IF (.NOT. Found) THEN
    nsubs=nsubs+1
    IF (nsubs  <=  klimit) THEN
      subname(nsubs)=sub
      k=nsubs
    ELSE
      WRITE(umMessage,'(A,i4,A)')'WARNING: More than', &
          klimit,' different subroutine calls to TIMER'
      CALL umPrint(umMessage,src='timer-timer1a')
      istop=1
      GO TO 9999
    END IF
  END IF

  !        Start timer for subroutine
  nenter(k)=nenter(k)+1
  elpstart = get_wallclock_time()
  cpustart = get_cpu_time()
  !      -----------------------------------------------------------------
ELSE IF ((i == 4) .AND. (istop == 0)) THEN

  !      Return from subroutine

  !        Stop timer
  cpuend = get_cpu_time()
  elpend = get_wallclock_time()
  elapse(k)=elapse(k)+(elpend-elpstart)
  time(k)=time(k)+cpuend-cpustart

  !        Find name of calling program
  Found=.FALSE.
  DO k=1,nsubs
    IF (subname(k) == returnam(j)) THEN
      Found=.TRUE.
      EXIT
    END IF
  END DO

  IF (.NOT. Found ) THEN
    WRITE(umMessage,'(3X,A,1X,A8,1X,A,1X,A8)') &
        'Calling prog:-', &
        returnam(j), &
        'not found, now in', &
        subname(j+1)
    CALL umPrint(umMessage,src='timer-timer1a')
    CALL umPrint('TIMER being DISABLED for the rest of this run', &
        src='timer-timer1a')
    istop=1
    GO TO 9999
  END IF

  !        Start timer for calling program
  elpstart = get_wallclock_time()
  cpustart = get_cpu_time()
  j=j-1

  !      -----------------------------------------------------------------
ELSE IF ((i <  1) .OR. (i >  6)) THEN

  !      If I<1 or I>6then there is an error. If 4<I<=6 then this call
  !      to TIMER is ignored. These values of I are recognised by the
  !      TIMER3A version.

  WRITE(umMessage,'(3x,A,1X,A8)') &
      'Illegal call to TIMER by subroutine',sub
  CALL umPrint(umMessage,src='timer-timer1a')

END IF

9999    CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE timer

END MODULE timer_mod

#endif

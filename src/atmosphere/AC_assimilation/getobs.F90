! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    MODULE GETOBS_MOD
!
!    SUBROUTINE GETOBS
!
!    Purpose : Extract appropriate observations from the OBS array   
!              (passed from RDOBS via argument list)
!              GETOBS called from AC gets list of obs in time window
!              GETOB2 called from AC2 gets lat, long, time and
!                     model obs type for an obs type.
!              GETOB3 called from VERTANL gets data specific to a type
!                     (eg data values and assoc error ratios)
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE getobs_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GETOBS_MOD'

CONTAINS

SUBROUTINE getobs (kact,assm_time,inobs,obs_no,lenobt,            &
                   obs,obs_flag,                                  &
                   timestep_no,dg_only,                           &
                   inobsdim,icode,cmessage)
!
!     --------------------------------------------------------
!   CALLED BY SUBROUTINE AC
!   THIS OBTAINS ALL OBSERVATIONS TOGETHER WITH THEIR POSITION
!   OF TYPE LACT(KACT) WHICH ARE WITHIN THE INSERTION PERIOD RELATIVE
!   TO THE ASSIMILATION TIME
!   (for scatwinds and sat120 a subset may be chosen)
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype, nproc
USE umPrintMgr
USE comobs_mod, ONLY: nobtypmx, obs_info
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER :: Iproc
!-----------------------------------------------------------------------
INTEGER :: kact      ! IN Index for this obs type (pos in LACT etc)
INTEGER :: inobs           ! IN  Total no of obs
INTEGER :: inobsdim        ! IN  Total no of obs (for dimensioning)
INTEGER :: obs_no(inobsdim)   ! IN  Pointers to obs to be assimilated
INTEGER :: obs_flag(inobsdim) ! IN  Observation flags
INTEGER :: timestep_no     ! IN  Timestep number
INTEGER :: lenobt          ! OUT No of obs to be assimilated

REAL :: obs(inobsdim,*) ! IN  Observation data for this type
REAL :: assm_time       ! IN  Assm time relative to start

LOGICAL :: dg_only         ! IN  Diagnostic only switch
INTEGER :: icode           ! IN  Return code
CHARACTER(LEN=errormessagelength) :: cmessage  ! IN  Error message
!-----------------------------------------------------------------------
!     Local variables
INTEGER :: inb       !  No of obs before insertion period
INTEGER :: ina       !  No of obs after insertion period
INTEGER :: inf       !  No of obs that have been flagged
INTEGER :: INT       !  No of obs thinned out (skipped)
INTEGER :: job       !  Loop counter over no of obs
INTEGER :: istart    !  First obs
INTEGER :: iskip     !  Skip ISKIP obs in loop over obs
INTEGER :: ip_time   !  Pointer to times in observation data
INTEGER :: obthin    !  Thinning factor for this ob type
INTEGER :: istat     !  for calls to GCOM routines
REAL :: tdiff     !  Time difference between obs and assm time
REAL :: tgetobb   !  Insertion period before obs time
REAL :: tgetoba   !  Insertion period after obs time

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GETOBS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ltimer_ac) CALL timer('GETOBS  ',3)

inb=0
ina=0
inf=0
INT=0

! SET UP INDEX VECTOR FOR GATHERING OBS
!  ALLOWING FOR OPTION TO THIN OBS :likely use GLOSS/SAT120/ERS-1 DATA
!  BUT NOT ON DIAGNOSTIC ONLY ITERATION
!     Get OBTHIN for this observation type
obthin = def_obthin( type_index(kact) )
IF (dg_only) THEN
  istart = 1
  iskip  = 1
ELSE
  IF (obthin >  1 .AND. inobs >  obthin) THEN
    istart = 1 + MOD(timestep_no,obthin)
    iskip  = obthin
  ELSE
    istart = 1
    iskip  = 1
  END IF
END IF

!     Get TGETOBB and TGETOBA for this observation type
tgetobb = def_tgetobb( type_index(kact) )
tgetoba = def_tgetoba( type_index(kact) )

!     Pointer to observation time for this obs type
ip_time = 3

IF (inobs >  0) THEN


  DO job = istart,inobs,iskip

    !     Get time difference between obs time and current assm time
    tdiff = obs(job,ip_time)-assm_time

    !     Use obs if : time diff <  tgetobb-0.5
    !         and if : time diff > -(tgetoba-0.5)
    !         and if : it has not been flagged
    !   -0.5 minutes helps make results same on different machines

    IF (tdiff >= tgetobb-0.5) THEN    !  Obs yet to be assimilated
      inb=inb+1
    ELSE IF (tdiff <= -tgetoba+0.5) THEN !  Obs has been assimilated
      ina=ina+1
    ELSE IF (obs_flag(job) /= 0) THEN  !  Obs has been flagged
      inf=inf+1
    ELSE                              !  Assimilate this observation
      lenobt=lenobt+1
      obs_no(lenobt) = obs_info % obs_no_st(kact)+job
    END IF
    !
  END DO

END IF  !INOBS >  0

INT=inobs-(lenobt+inb+ina+inf)
!
IF (ldiagac) THEN
  CountA(1)=lenobt
  CountA(2)=inb
  CountA(3)=ina
  CountA(4)=inf
  CountA(5)=INT
  IF (mype == 0) THEN
    DO job=1,5
      CountC(job)=CountA(job)
    END DO
  END IF


  DO iproc=1,nproc-1
    IF (mype == 0) THEN
      ! PE0 receives
      CALL gc_irecv(iproc,5,iproc,istat,countb,counta)
      DO job=1,5
        CountC(job)=CountC(job)+CountB(job)
      END DO
    ELSE IF (mype == iproc) THEN
      ! other PEs send
      CALL gc_isend(iproc,5,0,istat,countb,counta)
    END IF
  END DO


  WRITE(umMessage,'(A,I4,A,I6,A,I6,A,I6,A,I6,A,I6,A)')                  &
  ' TYPE',lact(kact),'  LENOBT',countc(1),'  OMITTED WERE :',    &
  countc(2),' (BEFORE) ',countc(3),' (AFTER) ',                  &
  countc(4),' (FLAGGED) ',countc(5),' (THINNED) '
  CALL umPrint(umMessage,src='getobs')
END IF
!
IF (ltimer_ac) CALL timer('GETOBS  ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE getobs
!
END MODULE getobs_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Get time records (indices as well as values) from a NetCDF
!   file that correspond to the closest times
!   before and after the current model_time.
!
! Method:
! 1) Call time2sec to get current model time in hours and seconds,
!    and then calculate the target time in the model in fractional
!    hours (fhr_now).
! 2) Loop over time records in the file, calling NETCDF_GET_TIME_REC
!    to get time for the previous and next record.
!    When fhr_now is between both records then get the indices and
!    current times (both as a date and in hours since base time) exit loop.
!    l_notime is used to identify when there are no more time records
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE glomap_clim_get_netcdffile_rec_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName =                          &
                                          'GLOMAP_CLIM_GET_NETCDFFILE_REC_MOD'

CONTAINS

SUBROUTINE glomap_clim_get_netcdffile_rec ( row_length, rows, fid, nfield,    &
                                            iupdn_type, model_time, l_first,  &
                                            irec, crec, fhr_now, ftrec,       &
                                            l_exact_match )

USE conversions_mod,           ONLY: &
    rsec_per_day

USE ereport_mod,               ONLY: &
    ereport

USE errormessagelength_mod,    ONLY: &
    errormessagelength

USE glomap_clim_netcdf_io_mod, ONLY: &
    netcdf_get_time_rec,             &
    iupdn_single,                    &
    iupdn_serial,                    &
    iupdn_cyclic,                    &
    cyclic_year_netcdf

USE nlstcall_mod,              ONLY: &
    lcal360

USE parkind1,                  ONLY: &
    jpim,                            &
    jprb

USE umPrintMgr,                ONLY: &
    umPrint,                         &
    umMessage,                       &
    PrintStatus,                     &
    PrStatus_Normal,                 &
    PrStatus_Diag

USE yomhook,                   ONLY: &
    lhook,                           &
    dr_hook

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: row_length     ! model horiz dimensions
INTEGER, INTENT(IN) :: rows           ! model horiz dimensions
INTEGER, INTENT(IN) :: fid            ! NetCDF file ID (already open)
INTEGER, INTENT(IN) :: nfield         ! nr of NetCDF field
INTEGER, INTENT(IN) :: iupdn_type     ! iupdn_serial=1, iupdn_cyclic=2
INTEGER, INTENT(IN) :: model_time(6)  ! Current model time:
                                      ! yr, mon, day, hr, min, sec

LOGICAL, INTENT(IN) :: l_first

INTEGER, INTENT(OUT):: irec(2)        ! Time record indices: prev & next
INTEGER, INTENT(OUT):: crec(2,6)      ! Readable time of records
                                      ! 1st dim is for prev & next
                                      ! 2nd dim for yr, mn, day, hr, mon, sec

LOGICAL, INTENT(OUT):: l_exact_match  ! T if no time interpolation required

REAL,    INTENT(OUT):: fhr_now        ! Current model time in fract hours
                                      ! since base time (0:00 h on 1st Jan)
REAL,    INTENT(OUT):: ftrec(2)       ! Prev/next time in hours since
                                      ! the same base time

! Local variables
INTEGER             :: daynow, secnow ! current model time in days and secs
                                      ! elapsed since 0:00 on 1st Jan

INTEGER             :: tday1, tday2, tsec1, tsec2, fdaynum ! prev/next times

INTEGER             :: model_yr                 ! current model year
INTEGER             :: i_time, i_time1, i_time2 ! time rec counters for NCDF 
                                                !  files

REAL, PARAMETER     :: fsec_to_day  = 1.0/rsec_per_day ! Conversion factor
REAL, PARAMETER     :: fsec_to_hour = 1.0/3600.0       ! Conversion factor

LOGICAL             :: l_notime       ! True if end-of-file
LOGICAL             :: l_new_cdf_yr   ! True if periodic data and need to read
                                      ! Dec of year i and Jan of year i+1

LOGICAL             :: l_wrap_yr      ! True when model time 16 Dec to 15 Jan

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

INTEGER                           :: errorstatus
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*),    PARAMETER  :: RoutineName='GLOMAP_CLIM_GET_NETCDFFILE_REC'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Start with negative indices for prev/next time records
! until they are found in the TIME_REC loop
irec(:) = -1

! Time in NetCDF files is in days since the base time (0:00 h on 1st Jan
! of a given year). The code below gets current model days/secs since
! the same base time so that we can interpolate.

SELECT CASE (iupdn_type)
CASE (iupdn_single, iupdn_serial) ! Single Time or Time series
  model_yr = model_time(1)
CASE (iupdn_cyclic) ! Periodic time
  ! When file is periodic then ignore year. Set to cyclic_year_netcdf to
  ! ensure consistency with netcdf_get_time_rec
  model_yr = cyclic_year_netcdf
CASE DEFAULT
  errorstatus = 54999
  WRITE (cmessage,'(A,I0)') 'Single time (0), Time series (1) and ' //        &
                             'Periodic (2) are the only supported time ' //   &
                             'record options. Value received: ', iupdn_type
  CALL ereport( RoutineName, errorstatus, cmessage )
END SELECT

! Get current time in days and seconds (daynow & secnow)
! elapsed since base time (0:00 h on 1st Jan)
! DEPENDS ON: time2sec
CALL time2sec ( model_yr, model_time(2), model_time(3), model_time(4),        &
                model_time(5), model_time(6), 0, 0, daynow, secnow, lcal360 )

! Calculate real hours from ref time
! Make target time to which we will interpolate
! equal to mid of update period and convert to real
fhr_now = REAL( daynow * 24 + INT (secnow * fsec_to_hour) )

! If extra diags messg then print time info only for first NetCDF field
IF (PrintStatus >= PrStatus_Diag .AND. nfield == 1) THEN
  WRITE (umMessage,'(A,A,I4,5(1x,I2),1x,I6,1x,I6,1x,F12.3)')                  &
         RoutineName, ' : Get data field for ',                               &
         model_time(1:6), daynow, secnow, fhr_now
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Deal with single time separately

IF ( iupdn_type  == iupdn_single ) THEN
  ! Get previous time (from a register in NetCDF file) in days and secs:
  ! tday1, tsec1
  i_time    = 1
  irec(1)   = 1
  irec(2)   = 1
  CALL netcdf_get_time_rec ( fid, i_time, iupdn_type,                         &
                             l_first, l_new_cdf_yr, tday1, tsec1, l_notime)
  
  ! Single time in file - get the first time record
  ftrec(1) = REAL(tday1) * 24.0 + REAL(tsec1) * fsec_to_hour
  ftrec(2) = REAL(tday1) * 24.0 + REAL(tsec1) * fsec_to_hour
  
  ! Get time yr-mon-day-hr-min-sec and store them in crec (1,1:6) for
  ! prev time record (given as tday1-tsec1). Useful for debugging.
  ! DEPENDS ON: sec2time
  CALL sec2time (tday1, tsec1, 0, 0, crec(1,1), crec(1,2), crec(1,3),         &
                 crec(1,4), crec(1,5), crec(1,6), fdaynum, lcal360 )
  
  crec  (2,:) = crec  (1,:)     ! record times: yr, mon, day, hr, min, sec
  l_exact_match = .TRUE.        ! no time interpolation required
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Time is either serial or cyclic

! Loop over time records in file until we find two consecutive records
! before and after the model time.
i_time       = 1        ! Initialise time rec counter
l_notime     = .FALSE.  ! Initially assume that time records found
l_wrap_yr    = .FALSE.  ! Initially assume model time not 16 Dec to 15 Jan
l_new_cdf_yr = .FALSE.  ! Initially assume that no need to change year
                        !   in the NetCDF files

time_rec: DO    ! loop over time records
  
  ! Get next time (from a register in NetCDF file) in days and secs:
  ! tday2, tsec2
  i_time2    = i_time
  CALL netcdf_get_time_rec ( fid, i_time2, iupdn_type,                        &
                             l_first, l_new_cdf_yr, tday2, tsec2, l_notime)
  
  ! The line below checks that the file has not run out of time records.
  IF ( l_notime ) THEN
    cmessage = 'Invalid time record 2, End-of-file reached ?'
    CALL ereport ( RoutineName, i_time2, cmessage )
  END IF
  
  ! Convert record times to fractional hours to compare with fhr_now
  ftrec(2) = REAL(tday2) * 24.0 + REAL(tsec2) * fsec_to_hour ! next record
  
  ! Check if reference time is greater than model time
  IF ( fhr_now <= ftrec(2) ) THEN
    IF ( ( iupdn_type == iupdn_cyclic ) .AND. ( i_time == 1 ) ) THEN
      ! Situation where model time between 1st and 15th January
      i_time    = 13
      l_wrap_yr = .TRUE.
      EXIT time_rec
    ELSE
      ! Found i_time 2
      EXIT time_rec
    END IF
  END IF
  
  ! increment time rec counter => next rec in file
  i_time = i_time + 1
  
  ! Check if we have exceeded 12th month for cyclic time
  ! This is situation where model time between 16th and 31st December
  IF ( ( iupdn_type == iupdn_cyclic ) .AND. ( i_time == 13 ) ) THEN
    l_wrap_yr = .TRUE.
    EXIT time_rec
  END IF
  
END DO  time_rec       ! loop over time records

! Get previous time (from a register in NetCDF file) in days and secs:
! tday1, tsec1
i_time1         = i_time - 1
CALL netcdf_get_time_rec ( fid, i_time1, iupdn_type,                          &
                           l_first, l_new_cdf_yr, tday1, tsec1, l_notime)

! The line below checks that the file has not run out of time records.
IF ( l_notime ) THEN
  cmessage = 'Invalid time record 1, End-of-file reached ?'
  CALL ereport ( RoutineName, i_time1, cmessage )
END IF

! Convert record times to fractional hours to compare with fhr_now
ftrec(1) = REAL(tday1) * 24.0 + REAL(tsec1) * fsec_to_hour ! prev record

! Check if current time is smaller than first time series record
IF ( ( iupdn_type == iupdn_serial ) .AND. ( ftrec(1) > fhr_now ) ) THEN
  cmessage = 'Invalid time record 1, time series begins after current time'
  CALL ereport ( RoutineName, i_time, cmessage )
END IF

! Deal with situation where model time is between mid December and mid January
IF (l_wrap_yr) THEN
  ! Get next time (from a register in NetCDF file) in days and secs:
  ! tday2, tsec2
  i_time2 = 1
  ! We need to change year in the NetCDF files
  l_new_cdf_yr   = .TRUE.
  CALL netcdf_get_time_rec ( fid, i_time2, iupdn_type,                        &
                             l_first, l_new_cdf_yr, tday2, tsec2, l_notime)
  
  ! The line below checks that the file has not run out of time records.
  IF ( l_notime ) THEN
    cmessage = 'Invalid time record 2, End-of-file reached ?'
    CALL ereport ( RoutineName, i_time2, cmessage )
  END IF
  
  ! Convert record times to fractional hours to compare with fhr_now
  ftrec(2) = REAL(tday2) * 24.0 + REAL(tsec2) * fsec_to_hour ! next record
  
  ! Recalculate fhr_now if model time between 1st January and 15th January
  IF ( model_time(2) == 1 ) THEN
    
    ! Recalculate the model time 'fhr_now' as done before
    ! but indicating that year is i+1 so that interpolations can be made.
    ! DEPENDS ON: time2sec
    CALL time2sec ( model_yr + 1,  model_time(2), model_time(3),              &
                    model_time(4), model_time(5), model_time(6), 0, 0,        &
                    daynow, secnow, lcal360 )
    
    fhr_now = REAL( daynow * 24 + INT(secnow * fsec_to_hour) )
    
    IF (PrintStatus >= PrStatus_Diag .AND. nfield == 1) THEN
      WRITE (umMessage,'(A,A,I4,5(1x,I2),1x,I6,1x,I6,1x,F12.3,1x,A)')         &
                       RoutineName, ' : Get climatology for ',                &
                       model_time(1:6), daynow, secnow, fhr_now,              &
                       '(after shifting model time one year forward)'
      CALL umPrint(umMessage,src=RoutineName)
    END IF
  END IF
END IF

irec(1) = i_time1             ! Prev time record for interpolation
irec(2) = i_time2             ! Next time record for interpolation

! Get time yr-mon-day-hr-min-sec and store them in crec (1,1:6) for
! prev time record (given as tday1-tsec1). Useful for debugging.
! DEPENDS ON: sec2time
CALL sec2time ( tday1, tsec1, 0, 0, crec(1,1), crec(1,2), crec(1,3),          &
                crec(1,4), crec(1,5), crec(1,6), fdaynum, lcal360 )

! Similarly, get time yr-mon-day-hr-min-sec for next time record
! DEPENDS ON: sec2time
CALL sec2time ( tday2, tsec2, 0, 0, crec(2,1), crec(2,2), crec(2,3),          &
                crec(2,4), crec(2,5), crec(2,6), fdaynum, lcal360 )

! Assume that not exact match, i.e. interpolation required
l_exact_match = .FALSE.

! If exact time match, no interpolation required. Therefore
! set both recs to same value
IF ( ABS (ftrec(1) - fhr_now) < EPSILON(fhr_now) ) THEN
  irec  (2)     = irec  (1)     ! record indices
  crec  (2,:)   = crec  (1,:)   ! record times: yr, mon, day, hr, min, sec
  ftrec (2)     = ftrec (1)     ! hours since base time
  l_exact_match = .TRUE.
END IF

IF ( ABS (ftrec(2) - fhr_now) < EPSILON(fhr_now) ) THEN
  irec  (1)     = irec  (2)
  crec  (1,:)   = crec  (2,:)
  ftrec (1)     = ftrec (2)
  l_exact_match = .TRUE.
END IF

! If extra diag mssg then print time info for required NetCDF fields
IF (PrintStatus >= PrStatus_Diag .AND. nfield == 1) THEN
  IF ( l_exact_match ) THEN
    WRITE (umMessage, '(A,A,1X,I5,1X,F12.3)')                                 &
                      RoutineName,                                            &
                      ' : Time found for this reg. in NetCDF file: ',         &
                      irec(1), ftrec(1)
    CALL umPrint(umMessage,src=RoutineName)
  ELSE
    WRITE (umMessage, '(A,A,1X,I5,1X,I5,1X,F12.3,1X,F12.3)')                  &
                      RoutineName,                                            &
                      ' : Time found between these regs. in NetCDF file:',    &
                      irec(1), irec(2), ftrec(1), ftrec(2)
    CALL umPrint(umMessage,src=RoutineName)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
END SUBROUTINE glomap_clim_get_netcdffile_rec

END MODULE glomap_clim_get_netcdffile_rec_mod

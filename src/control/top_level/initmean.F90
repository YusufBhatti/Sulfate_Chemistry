! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE initmean (                                             &
   isubmodl,icode,cmessage)

! Purpose:
!   To set up and check the input parameters for the
!   meaning subroutines
!
! Code description:
!   Language: Fortran 95.
!   This code is written to UMDP3 standards.
!
! External documentation:
!   On-line UM document C5 - Control of means calculations
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE um_parparams
USE control_max_sizes
USE lookup_addresses
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim,   &
    dumpfreqim, meanfreqim, mean_reftimeim
USE submodel_mod, ONLY: atmos_sm
USE nlstcall_mod, ONLY: model_basis_time, lcal360
USE history, ONLY: h_stepim, mean_offsetim, offset_dumpsim,     &
    mean_numberim
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lookup,                            &
    len_dumphist, len_fixhd, mpp_len1_lookup

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Arguments
INTEGER,           INTENT(IN) :: isubmodl ! Submodel identifier
INTEGER,           INTENT(OUT):: icode    ! Return code; success=0, error> 0
CHARACTER(LEN=errormessagelength), INTENT(OUT):: cmessage 
                                          ! Error message if ICODE > 0


!
!      Local variables and arrays

INTEGER  :: ifind,i           ! Loop counts
INTEGER  :: nmeans            ! No. of means chosen (fixed)
INTEGER  :: mean_reftime_days ! Reference time for period means (D)
INTEGER  :: mean_reftime_secs ! Reference time for period means (s)
INTEGER  :: mean_start_days   ! Start time for model run (days)
INTEGER  :: mean_start_secs   ! Start time for model run (s)
INTEGER  :: mean_offset_days  ! Offset from mean ref. time (days)
INTEGER  :: mean_offset_secs  ! Offset from mean ref. time (s)
INTEGER  :: mean_offset_steps ! Offset from mean ref. time (steps)
INTEGER  :: mean_freq_dumps   ! Mean frequency in dumps (updated)
INTEGER  :: iyear             ! Local variables for year,
INTEGER  :: imonth            ! month etc
INTEGER  :: iday
INTEGER  :: ihour
INTEGER  :: iminute
INTEGER  :: isecond

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INITMEAN'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
cmessage=' '

icode=0

! Define mode of use for means program
IF (isubmodl == atmos_sm) THEN
  CALL umPrint('',src='initmean')
  CALL umPrint('',src='initmean')
  CALL umPrint( 'INITMEAN: ***** Called in ATMOSPHERIC mode *****', &
      src='initmean')
END IF

! Check input parameters
nmeans=0

DO ifind = 1,4
  IF (isubmodl == atmos_sm) THEN
    IF (meanfreqim(ifind,atmos_sm) >  0) THEN
      nmeans = nmeans+1
    END IF
    IF (meanfreqim(ifind,atmos_sm) == 1 .OR.                      &
       meanfreqim(ifind,atmos_sm) <  0) THEN
      icode = 1
      cmessage='INITMEAN: Invalid atmos mean frequency'
      WRITE(umMessage,'(A,I3,I3,A,I5)') 'INITMEAN: MEANFREQ(',ifind,atmos_sm, &
         ') set to ', meanfreqim(ifind,atmos_sm)
      CALL umPrint(umMessage,src='initmean')
      GO TO 9999
    END IF
  END IF
END DO ! IFIND
IF (isubmodl == atmos_sm) THEN
  mean_numberim(isubmodl) = nmeans
  IF (mean_numberim(isubmodl) == 0) THEN
    CALL umPrint( 'INITMEAN: No means requested',src='initmean')
    GO TO 9999
  END IF
END IF

!      If means are to be created:
!      Establish whether an offset exists between the
!      reference time for means creation and the start
!      of the integration
!      N.B. In the case of a restart, this check is ignored

IF (isubmodl == atmos_sm) THEN

  iyear   = mean_reftimeim(1,isubmodl)
  imonth  = mean_reftimeim(2,isubmodl)
  iday    = mean_reftimeim(3,isubmodl)
  ihour   = mean_reftimeim(4,isubmodl)
  iminute = mean_reftimeim(5,isubmodl)
  isecond = mean_reftimeim(6,isubmodl)

  IF (iyear   == 0 .AND.                                           &
     imonth  == 0 .AND.                                           &
     iday    == 0 .AND.                                           &
     ihour   == 0 .AND.                                           &
     iminute == 0 .AND.                                           &
     isecond == 0) THEN

    mean_offsetim(isubmodl) = mean_numberim(isubmodl)
    offset_dumpsim(isubmodl) = 0
    CALL umPrint('INITMEAN: No offset specified for means creation', &
        src='initmean')

  ELSE

    ! DEPENDS ON: time2sec
    CALL time2sec(iyear,imonth,iday,ihour,iminute,isecond   &
       ,0,0,mean_reftime_days,mean_reftime_secs,            &
       lcal360)

    iyear   = model_basis_time(1)
    imonth  = model_basis_time(2)
    iday    = model_basis_time(3)
    ihour   = model_basis_time(4)
    iminute = model_basis_time(5)
    isecond = model_basis_time(6)

    ! DEPENDS ON: time2sec
    CALL time2sec(iyear,imonth,iday,ihour,iminute,isecond   &
       ,0,0,mean_start_days,mean_start_secs,                &
       lcal360)

    ! DEPENDS ON: tim2step
    CALL tim2step(mean_start_days-mean_reftime_days,             &
       mean_start_secs-mean_reftime_secs,                        &
       steps_per_periodim(isubmodl),secs_per_periodim(isubmodl), &
       mean_offset_steps)
    offset_dumpsim(isubmodl) = mean_offset_steps/dumpfreqim(isubmodl)
    CALL umPrint('INITMEAN: Offset set up for means creation',src='initmean')
    WRITE(umMessage,'(A,I7)')' OFFSET_DUMPSim(ISUBMODL)=', &
        offset_dumpsim(isubmodl)
    CALL umPrint(umMessage,src='initmean')
    IF (h_stepim(isubmodl) == 0) THEN ! firststep
      mean_offsetim(isubmodl) = 0
      DO i = 1,mean_numberim(isubmodl) ! meannum
        IF (i == 1) THEN ! i1
          IF (lcal360) THEN ! 360 day year
            mean_freq_dumps = meanfreqim(i,isubmodl)
          ELSE
            ! gregorian calendar
            IF (iday == 1 .AND. ihour == 0 .AND. iminute == 0 .AND. &
               isecond == 0) THEN
              mean_offsetim(isubmodl) = 1
            ELSE
              mean_offsetim(isubmodl) = 0
            END IF
          END IF ! 360 day year
        ELSE
          mean_freq_dumps = mean_freq_dumps*meanfreqim(i,isubmodl)
        END IF  ! i1
        ! The concept of a user defined mean offset only makes sense for
        ! meanings done as multiples of dump frequency. For the gregorian
        ! calendar the meaning must be done at the end of a month, season,
        ! or year with no exceptions.
        IF (lcal360) THEN
          IF (MOD(offset_dumpsim(isubmodl),mean_freq_dumps) == 0) THEN
            mean_offsetim(isubmodl) = mean_offsetim(isubmodl)+1
          END IF
        END IF
      END DO  ! meannum
    END IF  ! firststep
  END IF
END IF

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE initmean


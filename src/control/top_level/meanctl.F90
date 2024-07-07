! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!      Subroutine: MEANCTL-----------------------------------------
!
!      Purpose: To accumulate partial sums and create time-meaned data
!
!      Programming standard: UM Doc Paper 3
!
!      External documentation: UMDP C5 - Control of means calculations
!
!      ------------------------------------------------------------
!      Interface and arguments:

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE meanctl (                                              &
  ind_im,meanlev,icode,cmessage)
!
USE filenamelength_mod,       ONLY: &
    filenamelength
USE yomhook,                  ONLY: lhook, dr_hook
USE parkind1,                 ONLY: jprb, jpim
USE io_configuration_mod,     ONLY: io_field_padding
USE io, ONLY: setpos
USE io_constants, ONLY: ioOpenReadWrite, ioNoDelete
USE model_file,               ONLY: model_file_close, model_file_open
USE atm_fields_bounds_mod
USE UM_ParVars
USE decomp_params, ONLY: decomp_standard_atmos, decomp_standard_ocean, &
                         decomp_nowrap_ocean
USE Control_Max_Sizes
USE Decomp_DB
USE lookup_addresses
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim,    &
    dumpfreqim, meanfreqim, mean_reftimeim, pp_len2_meanim, &
    psum_filename_base
USE dump_headers_mod, ONLY: a_fixhd, a_lookup
USE submodel_mod, ONLY: submodel_for_sm, atmos_im
USE stash_array_mod, ONLY: totitems, stlist
USE stparam_mod, ONLY: st_macrotag, st_d1pos, s_modl, st_freq_code,&
                       st_dump_level_output_length
USE nlstcall_mod, ONLY: lmean, lclimrealyr, lcal360
USE conversions_mod, ONLY: isec_per_day
USE history, ONLY: mean_offsetim, offset_dumpsim,     &
                   mean_numberim, run_meanctl_indicim

USE file_manager, ONLY: assign_file_unit, release_file_unit
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY:                                       &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,   &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,    &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_data,   &
    a_len_extcnst, a_len_inthd, a_len_realhd, intf_len2_coldepc,      &
    intf_len2_levdepc, intf_len2_rowdepc, intf_lookupsa, land_field,  &
    len1_lookup, len_dumphist, len_fixhd,                             &
    len_tot, max_intf_model_levels, max_lbcrow_length, max_lbcrows,   &
    model_levels, mpp_len1_lookup, n_cca_lev, n_intf_a, n_obj_d1_max, &
    pp_len_inthd, pp_len_realhd, sm_levels, st_levels,                &
    theta_field_size, tpps_ozone_levels, tr_lbc_ukca, tr_lbc_vars,    &
    tr_ukca, tr_vars
USE d1_array_mod, ONLY: d1_imodl, d1_no_levels, d1_lookup_ptr, d1,    &
                        d1_addr, no_obj_d1
USE model_time_mod, ONLY: &
    i_day, i_hour, i_minute, i_month, i_second, i_year, stepim
USE errormessagelength_mod, ONLY: errormessagelength
USE setperlen_mod, ONLY: setperlen

USE Packing_Codes_Mod, ONLY: PC_Cray32_Packing

IMPLICIT NONE

!   Arguments


INTEGER ::                                                    &
   ind_im,                                                    &
                              ! IN Internal model indicator
   meanlev,                                                   &
                              ! INOUT  Mean level indicator
   icode             ! OUT return code; successful=0, error> 0


CHARACTER(LEN=errormessagelength) ::                          &
   cmessage          ! OUT Error message if ICODE > 0


!      Local variables and arrays


INTEGER ::                                                        &
   head_len                                                       &
                              ! Length of each record in header
   ,head_size              ! Size of header on disk

PARAMETER(                                                        &
   head_len=3                                                     &
   )

INTEGER ::                                                        &
   head_out(head_len,totitems) ! Header info for output ps file

INTEGER ::                                                        &
   ifind,                                                     &
                              ! Loop counter
   step_dumps,                                                &
                              ! Timestep (in multiples of restart
                              !           dump frequency)
   nmeans,                                                    &
                              ! No. of means chosen (fixed, unless
                              ! there is a mean offset)
   residu,                                                    &
                              ! Reference for partial sum data
   means_total,                                               &
                              ! Number of mean dumps this timestep
   nmvals(4),                                                 &
                              ! Absolute meaning periods (in
                              ! multiples of restart dump frequency)
   ps_flag(4),                                                &
                              ! Flag for partial sum updating
   periodlen,                                                 &
                              ! Length in days of current period N
   periodlendm,                                               &
                              ! Length in dumps of period 1 or
                              ! in days of current period N>1
   dumps_per_day                                              &
                              ! Number of restart dumps per day
   ,orig_decomp                                               &
                              ! Used to check for change in
   ,new_decomp                                                &
                              ! decomposition
   ,tag                                                       &
                              ! Stash tag
   ,ptd1                                                      &
                              ! Pointer to D1_ADDR information
   ,modnum                                                    &
                              ! Pointer to D1_ADDR submodel
   ,ptl                                                       &
                              ! Pointer to LOOKUP information
   ,global_length                                             &
                              ! Length of global field
   ,citems              ! Counter for no of objects to mean

INTEGER ::                                                        &
   index_read,                                                &
                              ! Specific index no for reading
   index_write          ! Specific index no for writing
INTEGER ::                                                        &
   ft_read,                                                   &
                              ! Unit number for partial sum read
   ft_write,                                                  &
                              ! Unit number for partial sum write
   ft_to_close
                              ! Unit number of ps file which is
                              ! soon to be closed

INTEGER :: internal_model
INTEGER :: d1_addr_submodel_id  ! Submodel number in D1_ADDR

LOGICAL ::                                                        &
   lmeaninc             ! increment MEANS_TOTAL or not

INTEGER ::                                                        &
   ie                                                              &
                              ! loop counter over items
   , maxsize                                                       &
                              ! maximum dump output length
   ,lmaxsize                                                       &
                              ! Maximum output length per level
   ,levsize                ! Size of each level
INTEGER :: dump_freq
INTEGER :: mean_freq1,mean_freq2,mean_freq3,mean_freq4
INTEGER :: stashoutputtime,mean1freqts

INTEGER :: i_dump_period_days ! Dump period in days
INTEGER :: i_dump_period_secs ! Dump period in seconds
INTEGER :: i_sec_of_today     ! Number of seconds in day of current timestep

! Partial sum files  
CHARACTER(LEN=filenamelength) :: psname_read      ! Partial sum read filename
CHARACTER(LEN=filenamelength) :: psname_write     ! Partial sum write filename
CHARACTER(LEN=filenamelength) :: psname_to_close  ! Partial sum file temporary
! When toggling between partial sum pairs, the letters "a" and "b" are used in
! the filenames, this array maps the letters to indices
CHARACTER(LEN=1), PARAMETER   :: psum_letter(2) = ["a","b"]

INTEGER :: copy_len                  ! length of temp. copy
REAL, ALLOCATABLE :: d1_save(:)      ! temporary save space

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MEANCTL'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

cmessage=' '

!      Define mode of use for means program

IF (ind_im == 1) THEN
  CALL umPrint('MEANCTL: ***** Called in ATMOSPHERIC mode *****', &
      src='meanctl')
END IF
!
! ----------------------------------------------------------------------
!      Find out which mean datasets need to be created
!      on this timestep (if any) and set MEANS_TOTAL accordingly
! ----------------------------------------------------------------------
!

icode=0
means_total=0

!      Initially check validity of call to subroutine

IF (dumpfreqim(ind_im) == 0) THEN! Is mean dump production off?
  icode=1
  cmessage='MEANCTL: Invalid call to subroutine'
  WRITE(umMessage,'(A,I3,A,I5)') 'meanctl: dumpfreq(',ind_im,')= ', &
     dumpfreqim(ind_im)
  CALL umPrint(umMessage,src='meanctl')
  GO TO 9999
ELSE IF (MOD(stepim(ind_im),dumpfreqim(ind_im)) /= 0) THEN
  icode=2        ! This is not a dumping timestep
  cmessage='MEANCTL: Incorrect timestep to call subroutine'
  CALL umPrint( 'MEANCTL: STEP is not a multiple of DUMPFREQ',src='meanctl')
  WRITE(umMessage,'(A,I3,A,I5,A,I3,A,I5)')          &
     '         step(',   ind_im,            &
     ')= ',              stepim(ind_im),    &
     ' dumpfreq(',       ind_im,            &
     ')= ',              dumpfreqim(ind_im)
  CALL umPrint(umMessage,src='meanctl')
  GO TO 9999
ELSE
  step_dumps = (stepim(ind_im)/dumpfreqim(ind_im))+             &
     offset_dumpsim(ind_im)
END IF

!      Pick up number of means chosen from history file (MEAN_NUMBERim)
!      or number determined by the offset from the reference time whilst
!      the staggered start of means production unwinds  (MEAN_OFFSETim)

DO ifind=1,mean_numberim(ind_im)
  IF (ifind == 1)nmvals(ifind) = meanfreqim(ifind,ind_im)
  IF (ifind >  1)nmvals(ifind) = meanfreqim(ifind,ind_im)*       &
     nmvals(ifind-1)
END DO

IF (lclimrealyr) THEN
  dumps_per_day=(24*3600*steps_per_periodim(ind_im))/             &
     (dumpfreqim(ind_im)*secs_per_periodim(ind_im))
END IF

IF (mean_offsetim(ind_im) == mean_numberim(ind_im)) THEN
  nmeans=mean_numberim(ind_im)
ELSE
  ! there is a non-zero offset so increment mean_offset
  ! when run is partway into the current mean period.
  IF (lclimrealyr) THEN
    ! Turn on monthly/seasonal/annual meaning once past start of
    ! first full month/season/year.
    ! initmean would have turned it on already if run started at
    ! start of a month/season/year
    i_dump_period_days = dumpfreqim(atmos_im)/steps_per_periodim(atmos_im)
    i_dump_period_secs = (isec_per_day*dumpfreqim(atmos_im))/ &
       steps_per_periodim(atmos_im) - isec_per_day*i_dump_period_days
    i_sec_of_today = i_hour * 3600 + i_minute * 60 +i_second

    IF (i_sec_of_today == i_dump_period_secs .AND.          &
       (i_day - i_dump_period_days) == 1) THEN
      IF (mean_offsetim(ind_im)  ==  0) THEN

        mean_offsetim(ind_im) = 1 ! months
      END IF

      IF (mean_offsetim(ind_im)  ==  1) THEN
        IF (MOD((i_month-mean_reftimeim(2,ind_im)),3) == 0) THEN
          mean_offsetim(ind_im) = 2 ! seasons
        END IF
      END IF

      IF (mean_offsetim(ind_im)  ==  2) THEN
        IF (i_month == mean_reftimeim(2,ind_im)) THEN
          mean_offsetim(ind_im) = 3 ! years
        END IF  ! of test on ifind, i_day, etc.
      END IF
    END IF
  ELSE  ! for 360d year, use array of no. of mean dumps
    DO ifind=mean_offsetim(ind_im)+1,mean_numberim(ind_im)
      IF (step_dumps <  0) THEN ! mean ref time not reached yet
        residu=1-nmvals(ifind)
      ELSE                     ! mean ref time has been passed
        residu=1
      END IF
      IF (MOD(step_dumps,nmvals(ifind)) == residu) THEN
        mean_offsetim(ind_im)=mean_offsetim(ind_im)+1
      END IF
    END DO
  END IF
  nmeans=mean_offsetim(ind_im)
  WRITE(umMessage,'(A,I3,A,I5)')' mean_offset(',ind_im,')=', &
      mean_offsetim(ind_im)
  CALL umPrint(umMessage,src='meanctl')
END IF  ! end of setting NMEANS
!
!       If no processing is required (because of staggered
!       start in means production) then skip to end of subroutine
!
IF (nmeans == 0) THEN
  icode=-1
  cmessage='MEANCTL: No accumulation/meaning done this step'
  CALL umPrint('MEANCTL: No accm/meaning due to staggered start', &
      src='meanctl')
  GO TO 9999
END IF

!
!       Output message whilst stagger unwinds
!
DO ifind=1,mean_numberim(ind_im)
  IF (nmeans <  ifind) THEN
    WRITE(umMessage,'(A,I1,A,A)')'meanctl: period_',ifind, &
        ' mean not activated',   &
        ' because of staggered start in means production'
    CALL umPrint(umMessage,src='meanctl')
  END IF
END DO

DO ifind=1,nmeans
  ps_flag(ifind)=0
  lmeaninc=.FALSE. ! Initialise to avoid false positive results

  ! Find out if the end of any meaning period has been reached. If so,
  ! increment MEANS_TOTAL
  IF (lclimrealyr) THEN
    IF ((i_day == 1) .AND. (i_hour == 0) .AND.         &
       (i_minute == 0) .AND. (i_second == 0)) THEN
      ! Climate meaning periods for Gregorian calendar must always
      ! be aligned with start of the month
      IF ((ifind == 1) .AND.                       &
         (stepim(ind_im) > steps_per_periodim(ind_im))) THEN
        lmeaninc=.TRUE.                               ! monthly
      ELSE IF ((ifind == 2) .AND.                  &
         ! seasonal
         (stepim(ind_im) > steps_per_periodim(ind_im)) .AND.   &
         MOD((i_month-mean_reftimeim(2,ind_im)),3) == 0) THEN
        lmeaninc=.TRUE.
      ELSE IF ((ifind == 3) .AND.                  &
         ! annual
         (stepim(ind_im) > steps_per_periodim(ind_im)) .AND.   &
         MOD((i_month-mean_reftimeim(2,ind_im)),12) == 0) THEN
        lmeaninc=.TRUE.
      END IF
    END IF
  ELSE
    ! 360-day calendar
    IF (MOD(step_dumps,nmvals(ifind)) == 0) THEN
      lmeaninc=.TRUE.
    END IF
  END IF  !
  IF (lmeaninc) THEN
    means_total=means_total+1
    WRITE(umMessage,'(A,I1,A)') 'meanctl: period_',ifind,    &
       ' mean to be created this timestep'
    CALL umPrint(umMessage,src='meanctl')
  END IF  ! end of IF test on lmeaninc

  ! Find out if run is one period(N-1) into the period(N), in which case
  ! set PS_FLAG(IFIND)=1 because there will not already be a file for
  ! that partial sum on the disk, so ACUMPS must get it from D1.

  IF (lclimrealyr) THEN
    IF (ifind  ==  1) THEN
      ! Convert dump periods into seconds
      i_dump_period_secs = (isec_per_day*dumpfreqim(atmos_im))/ &
                                 steps_per_periodim(atmos_im)
      ! How many seconds are we into the current month
      i_sec_of_today = (i_day - 1) * isec_per_day + i_hour * 3600 + &
         i_minute * 60 + i_second
      ! If equal, we are one mean dump period into the period 1 mean period
      IF (i_sec_of_today == i_dump_period_secs) THEN
        ps_flag(1)=1                              ! months
      END IF
    ELSE IF (ifind  ==  2 .AND.                             &
       ! First mean dump of a period 2 period is at the start of the month
       (i_day  ==  1) .AND. (i_hour  ==  0) .AND.           &
       (i_minute == 0) .AND. (i_second == 0) .AND.          &
       MOD((i_month-1-mean_reftimeim(2,ind_im)),3) == 0)    &
       THEN
      ps_flag(2)=1                                   ! seasons
    ELSE IF (ifind  ==  3 .AND.                             &
       ! First mean dump of a period 3 period is at the start of the month
       (i_day  ==  1) .AND. (i_hour  ==  0) .AND.           &
       (i_minute == 0) .AND. (i_second == 0) .AND.          &
       MOD(((i_month-3)-mean_reftimeim(2,ind_im)),12) == 0) &
       THEN
      ps_flag(3)=1                                   ! years
    END IF
  ELSE  ! for 360d year, check using array of no. of mean dumps
    IF (step_dumps <  0) THEN
      IF (ifind == 1)residu=1-nmvals(ifind)
      IF (ifind /= 1)residu=nmvals(ifind-1)-nmvals(ifind)
    ELSE
      IF (ifind == 1)residu=1
      IF (ifind /= 1)residu=nmvals(ifind-1)
    END IF
    IF (MOD(step_dumps,nmvals(ifind)) == residu) THEN
      ps_flag(ifind)=1
    END IF
  END IF  ! end of test on lclimrealyr and setting of PS_FLAG


END DO       ! end of loop over IFIND from 1 to NMEANS

! check STASH items output time against restart dump period
dump_freq=dumpfreqim(ind_im)     ! Unit: timestep
mean_freq1=meanfreqim(1,ind_im)
mean_freq2=meanfreqim(2,ind_im)
mean_freq3=meanfreqim(3,ind_im)
mean_freq4=meanfreqim(4,ind_im)
! Calculating mean_freq1 in unit of timestep
mean1freqts=mean_freq1*dump_freq
! Create header for partial sum file
IF (ind_im == 1) THEN
  ! Set pointer to submodel info in D1_ADDR
  modnum=submodel_for_sm(ind_im)
  citems=0
  DO ie=1,totitems
    tag=stlist(st_macrotag,ie)/1000
    ! Get pointer to element in D1_ADDR array
    ptd1=stlist(st_d1pos,ie)
    IF (tag /= 0 .AND.                                              &
       stlist(s_modl,ie) == d1_addr(d1_imodl,ptd1,modnum)) THEN
      stashoutputtime=stlist(st_freq_code,ie)  ! Unit: timestep
      IF (stashoutputtime >= mean1freqts) THEN
        CALL umPrint('ERROR: STASH output time frequency greater', &
            src='meanctl')
        CALL umPrint('       than climate mean period frequency!', &
        src='meanctl')
        WRITE(umMessage,'(A,I5,A,A,I5,A)')                           &
           'stash output frequency=',stashoutputtime,'(ts)'  &
           ,'  period 1 period=',mean1freqts,'(ts)'
        CALL umPrint(umMessage,src='meanctl')
        icode=4
        cmessage='MEANCTL: ERROR - STASHoutput > CMEAN period'
        GO TO 9999
      END IF
      ! Tagged for meaning and submodel information matches.
      ! Counter for no of variables that will be processed
      citems=citems+1
      ptl=d1_addr(d1_lookup_ptr,ptd1,modnum)
      global_length=stlist(st_dump_level_output_length,ie)
      global_length=global_length+1
      head_out(1,citems)=global_length
      head_out(2,citems)=1
      head_out(3,citems)=(((global_length+io_field_padding)       &
         /io_field_padding)*io_field_padding)
      IF (MOD((a_lookup(lbpack,ptl)),10)  ==  PC_Cray32_Packing) THEN
        IF (a_lookup(data_type,ptl)  ==  1) THEN
          ! Data to be packed so reduce data length accordingly
          global_length=(global_length+1)/2
          head_out(1,citems)=global_length
          head_out(2,citems)=2
          head_out(3,citems)=(((global_length+io_field_padding)   &
             /io_field_padding)*io_field_padding)
        END IF
      END IF
    END IF
  END DO
END IF

! Find the largest output field size to dimension I/O buffer array
! and the total size of storage required to set file size
lmaxsize=1
head_size=citems*head_len+io_field_padding+2
head_size=(head_size/io_field_padding)*io_field_padding
citems=0

DO ie=1,totitems
  tag=stlist(st_macrotag,ie)/1000
  ! Get pointer to element in D1_ADDR array
  ptd1=stlist(st_d1pos,ie)
  IF (tag /= 0 .AND.                                             &
     stlist(s_modl,ie) == d1_addr(d1_imodl,ptd1,modnum)) THEN
    ! Tagged for meaning and submodel information matches.
    citems=citems+1
    lmaxsize=MAX(lmaxsize,                                      &
       stlist(st_dump_level_output_length,ie))
  END IF
END DO
! Adding 1 extra space allowing for checksum
lmaxsize=lmaxsize+1

lmaxsize=((lmaxsize+io_field_padding)/io_field_padding)             &
   *io_field_padding

! Dump-packed fields will be rounded up to io_field_padding. Therefore,
! when unpacked, they will be rounded up to 2*io_field_padding. Add
! more space to account for this:
lmaxsize = lmaxsize + io_field_padding

!
!  *********************************************************************
!                   LOGICAL SUB-PROCESS C51
!      Start of default process: updating period_1 partial sum data
!  *********************************************************************
!
! Generate the first partial sum file names
index_read   =     run_meanctl_indicim(1,ind_im)
index_write  = 3 - run_meanctl_indicim(1,ind_im)
WRITE(psname_read, "(A,I1.1,A)") &
    TRIM(psum_filename_base), 1, psum_letter(index_read)
WRITE(psname_write, "(A,I1.1,A)") &
    TRIM(psum_filename_base), 1, psum_letter(index_write)

!      Temporary check on unit numbers

IF (ps_flag(1) /= 1) THEN
  WRITE(umMessage,'(A)') 'Period_1 data read from: '//TRIM(psname_read)
  CALL umPrint(umMessage,src='meanctl')
END IF
WRITE(umMessage,'(A)') 'Period_1 data written to: '//TRIM(psname_write)
CALL umPrint(umMessage,src='meanctl')

!
!                            STEP 1
!      Update or create period_1 partial sum data and write out
!      to period_1 partial sum dump
!

orig_decomp=current_decomp_type
new_decomp=orig_decomp

IF (ind_im==1 .AND. orig_decomp/=decomp_standard_atmos) THEN
  new_decomp = decomp_standard_atmos
ELSE IF (ind_im==2 .AND. orig_decomp/=decomp_nowrap_ocean) THEN
  new_decomp = decomp_nowrap_ocean
END IF
IF (new_decomp  /=  orig_decomp) THEN
  CALL change_decomposition(new_decomp,icode)
  IF (icode /= 0) THEN
    CALL umPrint('ERROR : MEANCTL',src='meanctl')
    WRITE(umMessage,'(A,I5)')'Failed to change decomposition to ',new_decomp
    CALL umPrint(umMessage,src='meanctl')
    cmessage='MEANCTL1: Failed to change decomposition'
    GO TO 9999
  END IF
END IF
IF (means_total >  0) THEN
  !
  !                  STEP 1
  !      Copy instantaneous dump to temporary memory space
  !
  !  Strictly, only climate mean tagged diagnostics need to be saved.
  IF (ind_im == 1) THEN
    copy_len = a_len_data+(model_levels+1)*theta_field_size
    ALLOCATE( d1_save(copy_len) )
    d1_save(1:copy_len) = d1(1:copy_len)
  END IF


END IF

!--preset the file lengths prior to the open
CALL assign_file_unit(psname_read,  ft_read,  handler="portio")
CALL assign_file_unit(psname_write, ft_write, handler="portio")

!      Open input and output partial sum files (preassigned names)

CALL model_file_open(ft_read, psname_read,                      &
                     read_write=ioOpenReadWrite, error=icode)
IF (icode /= 0) GO TO 9999

CALL model_file_open(ft_write, psname_write,                    &
                     read_write=ioOpenReadWrite, error=icode)
IF (icode /= 0) GO TO 9999

IF (ind_im == 1) THEN
  d1_addr_submodel_id = submodel_for_sm(1)
  ! DEPENDS ON: acumps
  CALL acumps(                                                  &
     no_obj_d1(d1_addr_submodel_id),                            &
     d1_addr(1,1,d1_addr_submodel_id),                          &
     a_len_data,d1,                                             &
     lmaxsize,means_total,                                      &
     ps_flag(1),ft_read,ft_write,lclimrealyr,meanlev,           &
     i_month,i_year,                                            &
     head_out,head_len,head_size,                               &
     stepim(ind_im),citems,a_fixhd(12),                         &
     icode,cmessage)
END IF

!      Check return code from ACUMPS

IF (icode /= 0) THEN
  GO TO 9999
END IF

!      Close input and output partial sum files


CALL model_file_close(ft_read, psname_read,                     &
                      delete=ioNoDelete, error=icode)

CALL model_file_close(ft_write, psname_write,                   &
                      delete=ioNoDelete, error=icode)

CALL release_file_unit(ft_read,  handler="portio")
CALL release_file_unit(ft_write, handler="portio")

!      Update RUN_MEANCTL_INDICim for period_1 data

run_meanctl_indicim(1,ind_im)=3-run_meanctl_indicim(1,ind_im)


!  *********************************************************************
!      End of default process: updating period_1 partial sum data
!  *********************************************************************

IF (means_total >  0) THEN

  !  *********************************************************************
  !                   LOGICAL SUB-PROCESS C52
  !      Start of means processing and updating of subsequent
  !      partial sum dump
  !  *********************************************************************

  DO meanlev=1,means_total

    !
    !                              STEP 2
    !      Generate period_N time-meaned data and store in main data block
    !
    !       If real-period meaning selected, find length of current period
    IF (lclimrealyr) THEN
      CALL setperlen (meanlev,i_month,i_year,periodlen)
      IF (meanlev == 1) THEN ! divisor only needs to be in terms
        ! of restart dumps for Period_1
        periodlendm=periodlen*dumps_per_day
      ELSE
        periodlendm=periodlen
      END IF
    END IF
    IF (ind_im == 1) THEN

      !      Open input partial sum file (preassigned or calculated name)

      IF (meanlev == 1) THEN
        index_read=run_meanctl_indicim(meanlev,ind_im)
      END IF

      ! Generate the name
      WRITE(psname_read, "(A,I1.1,A)") &
          TRIM(psum_filename_base), meanlev, psum_letter(index_read)

      WRITE(umMessage,'(a,I1,a)') 'period_',meanlev, &
          ' data read: '//TRIM(psname_read)
      CALL umPrint(umMessage,src='meanctl')

      CALL assign_file_unit(psname_read, ft_read, handler="portio")

      CALL model_file_open(ft_read, psname_read,                    &
                           read_write=ioOpenReadWrite, error=icode)
      IF (icode /= 0) GO TO 9999

      IF (lclimrealyr) THEN  ! Real-period meaning selected
        ! DEPENDS ON: meanps
        CALL meanps(                                                &
           no_obj_d1(d1_addr_submodel_id),                          &
           d1_addr(1,1,d1_addr_submodel_id),                        &
           a_len_data,d1,                                           &
        periodlendm                                                 &
           )
      ELSE            ! 360d year meaning selected
        ! DEPENDS ON: meanps
        CALL meanps(                                              &
           no_obj_d1(d1_addr_submodel_id),                        &
           d1_addr(1,1,d1_addr_submodel_id),                      &
           a_len_data,d1,                                         &
        meanfreqim(meanlev,ind_im)                                &
           )

      END IF            ! end of test on lclimrealyr
    END IF            ! end of test on IND_IM == 1

    !      Check return code from MEANPS

    IF (icode /= 0) THEN
      WRITE(umMessage,'(A,I3)') 'MEANCTL: Error at mean period ',meanlev
      CALL umPrint(umMessage,src='meanctl')
      GO TO 9999
    END IF

    psname_to_close=psname_read

    CALL setpos(ft_read,0,icode)
    ft_to_close=ft_read
    !
    !                           STEP 3.1
    !      Calculate mean diagnostics and extract PPfields from mean data
    !
    IF (ind_im == 1) THEN
      ! Find the largest output field size to dimension I/O buffer array
      maxsize=1

      DO ie=1,totitems
        tag=stlist(st_macrotag,ie)/1000
        IF (MOD(tag/(2**(meanlev-1)),2)  ==  1) THEN
          maxsize=MAX(maxsize,                                    &
             stlist(st_dump_level_output_length,ie))
        END IF
      END DO
      ! Extract mean diagnostics
      ! DEPENDS ON: meandiag
      CALL meandiag (                                             &
      ind_im,meanlev,pp_len2_meanim(meanlev,ind_im),step_dumps,  &
         nmvals(meanlev),                                       &
         maxsize,                                               &
         icode,cmessage)
    END IF
    IF (icode >  0) GO TO 9999
    !
    !                           STEP 4
    !      Check to see if period_N+1 partial sum data needs to updated
    !      or created.
    !      If so, proceed and write out to period_N+1 partial sum dump
    !

    IF (meanlev /= nmeans) THEN

      index_read   =     run_meanctl_indicim(meanlev+1,ind_im)
      index_write  = 3 - run_meanctl_indicim(meanlev+1,ind_im)

      IF (ind_im == 1) THEN

        !      Open input and output partial sum files (calculated names)
        WRITE(psname_read, "(A,I1.1,A)") &
            TRIM(psum_filename_base), meanlev + 1, psum_letter(index_read)

        CALL assign_file_unit(psname_read, ft_read, handler="portio")

        CALL model_file_open(ft_read, psname_read,     &
                             read_write=ioOpenReadWrite, error=icode)
        IF (icode /= 0) GO TO 9999

        WRITE(psname_write, "(A,I1.1,A)") &
            TRIM(psum_filename_base), meanlev + 1, psum_letter(index_write)

        CALL assign_file_unit(psname_write, ft_write, handler="portio")

        CALL model_file_open(ft_write,psname_write,    &
                             read_write=ioOpenReadWrite, error=icode)
        IF (icode /= 0) GO TO 9999

        d1_addr_submodel_id = submodel_for_sm(1)
        ! DEPENDS ON: acumps
        CALL acumps(                                                       &
           no_obj_d1(d1_addr_submodel_id),                                 &
           d1_addr(1,1,d1_addr_submodel_id),                               &
           a_len_data,d1,                                                  &
           lmaxsize,means_total,                                           &
           ps_flag(meanlev+1),ft_read,ft_write,                            &
           lclimrealyr,meanlev,i_month,i_year,                             &
           head_out,head_len,head_size,                                    &
           stepim(ind_im),citems,a_fixhd(12),                              &
        icode,cmessage)
      END IF

      !      Check return code from ACUMPS

      IF (icode /= 0) THEN
        GO TO 9999
      END IF


      !      Update RUN_MEANCTL_INDICim for period_N+1 data

      run_meanctl_indicim(meanlev+1,ind_im)=                      &
         3-run_meanctl_indicim(meanlev+1,ind_im)

    END IF

    !      Decide disposition of period_N+1 partial sum dumps
    !      NB: for restartability it is NOT safe to delete partial
    !      sum files.

    IF (meanlev /= nmeans) THEN

      CALL setpos(ft_read,0,icode)

      CALL model_file_close(ft_read, psname_read, delete=ioNoDelete, &
                            error=icode)
      CALL release_file_unit(ft_read, handler="portio")

      CALL setpos(ft_write,0,icode)

      CALL model_file_close(ft_write, psname_write, delete=ioNoDelete, &
                            error=icode)
      CALL release_file_unit(ft_write, handler="portio")

    END IF

    !      Decide disposition of remaining period_N partial sum dump
    !      NB: for restartability it is NOT safe to delete partial sum
    !      files.

    CALL model_file_close(ft_to_close, psname_to_close, delete=ioNoDelete, &
                          error=icode)
    CALL release_file_unit(ft_to_close, handler="portio")


  END DO ! Loop over meanlev
  !
  !                     STEP 6
  !      Read back instantaneous dump from temporary memory space
  !
  IF (ind_im == 1) THEN

    d1( 1:copy_len ) = d1_save( 1:copy_len )
    DEALLOCATE( d1_save )

  END IF
END IF
!
!  *********************************************************************
!      End of means processing and updating of subsequent
!      partial sum dumps
!  *********************************************************************
!

orig_decomp=current_decomp_type
new_decomp=orig_decomp
IF (ind_im == 1 .AND. orig_decomp /= decomp_standard_atmos) THEN
  new_decomp=decomp_standard_atmos
ELSE IF (ind_im==2 .AND. orig_decomp/=decomp_standard_ocean) THEN
  new_decomp=decomp_standard_ocean
END IF
IF (new_decomp  /=  orig_decomp) THEN
  CALL change_decomposition(new_decomp,icode)
  IF (icode /= 0) THEN
    CALL umPrint('ERROR : MEANCTL',src='meanctl')
    WRITE(umMessage,'(A,I5)')'Failed to change decomposition to ',new_decomp
    CALL umPrint(umMessage,src='meanctl')
    cmessage='MEANCTL1: Failed to change decomposition'
  END IF
END IF

9999 CONTINUE

!      Reset MEANLEV to zero

meanlev=0
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE meanctl

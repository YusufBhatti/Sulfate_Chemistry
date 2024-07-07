! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Routine: SETTSCTL -------------------------------------------------
!
!    Purpose: Sets timestep loop control switches and STASHflags.
!             Note STEP on entry is the values at the N+1 (ie. updated)
!             timelevel; control switches are therefore set using
!             A_STEP/O_STEP, whereas physical switches are set using
!             A_STEP-1/O_STEP-1 to ensure correct synchronisation.
!             Note also that step on entry is N (i.e. not updated) when
!             called from INITIAL.
!
!
!    Programming standard: UM Doc Paper 3, version 8.2
!
!    Project task: C0
!
!    External documentation: On-line UM document C0 - The top-level
!                            control system
!
!     ------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Top Level

SUBROUTINE settsctl (                                             &
       internal_model,linitial,meanlev,icode,cmessage)

USE rad_input_mod,     ONLY: l_radiation
USE set_rad_steps_mod, ONLY: set_l_rad_step

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IAU_mod, ONLY: l_iau
USE UM_ParParams
USE Control_Max_Sizes
USE eng_corr_inputs_mod, ONLY: a_energysteps, lenergy, lflux_reset

USE lookup_addresses
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim,     &
     dumpfreqim, dumptimesim, dumptimes_len1,                     &
     meanfreqim, meanfreq_len1,                                   &
     greg_dump_freq, greg_monthly_dump, end_of_run_dump
USE conversions_mod, ONLY: isec_per_day
USE dump_headers_mod, ONLY: a_fixhd, a_lookup
USE submodel_mod, ONLY: n_internal_model, atmos_im
USE stash_array_mod, ONLY:                                          &
    sf_calc, totitems, nitems, stlist, nsttims, sttabl, nsects, sf
USE stparam_mod, ONLY: st_model_code, st_sect_no_code, st_item_code,&
                  st_proc_no_code, st_freq_code, st_start_time_code,&
                 st_end_time_code, st_infinite_time, st_end_of_list,&
                 st_replace_code, st_input_code 

USE nlstcall_mod, ONLY: model_basis_time, &
                         ancil_reftime, &
                         lpp, &
                         lnc, &
                         ldump, &
                         lmean, &
                         lprint, &
                         lexit, &
                         lancillary, &
                         lboundary, &
                         lassimilation, &
                         lclimrealyr, &
                         l_fastrun, &
                         lcal360

USE file_manager, ONLY: init_file_loop, um_file_type

USE history, ONLY: model_data_time, offset_dumpsim

USE nlsizes_namelist_mod, ONLY:                                             &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,         &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,          &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst,      &
    a_len_inthd, a_len_realhd, a_prog_lookup, intf_len2_coldepc,            &
    intf_len2_levdepc, intf_len2_rowdepc, intf_lookupsa, len1_lookup,       &
    len_dumphist, len_fixhd,                                                &
    max_intf_model_levels, max_lbcrow_length, max_lbcrows, mpp_len1_lookup, &
    n_intf_a, pp_len_inthd, pp_len_realhd

USE model_time_mod, ONLY:                                                    &
    ancillary_stepsim, assim_extrastepsim, assim_firststepim, assim_stepsim, &
    basis_time_days, basis_time_secs, bndary_offsetim, boundary_stepsim,     &
    i_day, i_hour, i_minute, i_month, i_second, i_year, iau_dtresetstep,     &
    stepim, target_end_stepim
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport
USE horiz_grid_mod, ONLY: cartesian_grid

IMPLICIT NONE

INTEGER :: internal_model   ! IN  - internal model identifier
LOGICAL :: linitial         ! IN  - true in called from INITIAL
INTEGER :: meanlev          ! OUT - Mean level indicator
INTEGER :: icode            ! Out - Return code
CHARACTER(LEN=errormessagelength) :: cmessage  ! Out - Return error message

! ----------------------------------------------------------------------

!  Local variables

INTEGER :: i,nftunit        ! Loop counters
INTEGER :: itime            ! Loop index over selected dump times
INTEGER :: step             ! A_STEP or O_STEP for atmos/ocean
INTEGER :: IS,ii,il,im,ntab,it,ie  ! STASH variables
INTEGER :: modl             ! Int model no, read from STASH list arra
INTEGER :: step_assim       ! model step rel. to assim start

INTEGER :: jintf            ! Interface area index
INTEGER :: ancil_ref_days,ancil_ref_secs
INTEGER :: ancil_offset_steps ! offset of ref. from basis time
INTEGER :: secs_per_step    ! seconds per timestep
INTEGER :: months_in        ! Number of months into forecast
INTEGER :: secs_per_period
INTEGER :: steps_per_period
INTEGER :: dumpfreq
INTEGER :: offset_dumps
INTEGER :: exitfreq
INTEGER :: target_end_step
INTEGER :: ancillary_steps
INTEGER :: boundary_steps
INTEGER :: bndary_offset
INTEGER :: dumptimes(dumptimes_len1)
INTEGER :: meanfreq (meanfreq_len1)
INTEGER :: elapsed_days
INTEGER :: elapsed_secs
INTEGER :: elapsed_steps
INTEGER :: step_compare

TYPE(um_file_type), POINTER :: um_file

LOGICAL :: iau_resetdt ! If .TRUE., reset data time.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETTSCTL'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

icode=0

! Initialise internal_model
internal_model = 1 ! Atmosphere

!  1. Set timestep loop top-level control switches
!
!  1.0  Initialise control switches which are shared in coupled runs
!

step=                        stepim(internal_model)
secs_per_period=  secs_per_periodim(internal_model)
steps_per_period=steps_per_periodim(internal_model)
secs_per_step=secs_per_period/steps_per_period
dumpfreq=                dumpfreqim(internal_model)
offset_dumps=        offset_dumpsim(internal_model)
target_end_step=  target_end_stepim(internal_model)
ancillary_steps=  ancillary_stepsim(internal_model)
boundary_steps =   boundary_stepsim(internal_model)
bndary_offset  =   bndary_offsetim(internal_model)

exitfreq=dumpfreq

DO i=1,dumptimes_len1
  dumptimes(i)=       dumptimesim(i,internal_model)
END DO ! i
DO i=1,meanfreq_len1
  meanfreq(i)=         meanfreqim(i,internal_model)
END DO ! i
meanlev=0

lassimilation=.FALSE.
ldump=        .FALSE.
lexit=        .FALSE.
lmean=        .FALSE.
lprint=       .FALSE.
lancillary=   .FALSE.
lboundary=    .FALSE.
!
!  1.1  Set up PPfile and NCfile switches for the timestep
!
lpp =.FALSE.
lnc =.FALSE.

DO i = 1, 2
  NULLIFY(um_file)
  SELECT CASE(i)
    CASE(1)
      um_file => init_file_loop(handler="portio")
    CASE(2)
      um_file => init_file_loop(handler="netcdf")
  END SELECT
  DO WHILE (ASSOCIATED(um_file))
    IF (um_file % meta % is_output_file) THEN

      ! Ensure we fallback to switched off 
      um_file % meta % initialise = .FALSE.

      step_compare = step
      ! If init_file_type == 'c' compare init_start_step with previous timestep
      IF (ASSOCIATED(um_file % pp_meta)) THEN
        IF (um_file % pp_meta % init_file_type == 'c') THEN
          step_compare = step - 1
        END IF
      END IF

      ! Allow mix of real-month & regular reinit. periods on
      ! different units:
      IF (um_file % meta % init_steps <  0) THEN ! Real-month reinit.

        ! Select files to be reinitialised on this timestep:
        ! First initialisation may be on any timestep, not just 0
        IF (step_compare == um_file % meta % init_start_step) THEN

          um_file % meta % initialise =.TRUE.

        ELSE IF (.NOT. lcal360                                                &
                 .AND. (TRIM(um_file % id) /= "ppvar")                        &
                 .AND. (i_day     ==  1)                                      &
                 .AND. (i_hour    ==  secs_per_step/3600           )          &
                 .AND. (i_minute  ==  (MOD(secs_per_step,3600))/60 )          &
                 .AND. (i_second  ==   MOD(secs_per_step,60)       )          &
                 .AND. (step /= 1)) THEN

          ! Gregorian calendar:
          ! Reinitialisation can only be selected when day=1 and
          ! hours, minutes and seconds equals exactly one time-step
          ! This code relies on the rule that the time-step length
          ! divides into a day.

          ! Additionally initialisation only takes place on 1-month,
          ! 3-month (season) or 12-month boundaries. So calculate months
          ! since start of job, and request reinitialisation at the
          ! appropriate time according to which option selected.
          months_in = i_month - model_basis_time(2) +           &
             12 * (i_year - model_basis_time(1))

          IF (months_in >= um_file % meta % init_start_step .AND.    &
             (um_file % meta % init_end_step <  0 .OR.               &
              months_in <= um_file % meta % init_end_step)) THEN

            IF (um_file % meta % init_steps == -1) THEN
              um_file % meta % initialise = .TRUE. ! Months
            ELSE IF (um_file % meta % init_steps == -3 .AND.         &
                     MOD(months_in -                                 &
                      um_file % meta % init_start_step,3)  == 0) THEN
              um_file % meta % initialise = .TRUE. ! Seasons
            ELSE IF (um_file % meta % init_steps == -12 .AND.        &
                     MOD(months_in -                                 &
                      um_file % meta % init_start_step,12) == 0) THEN
              um_file % meta % initialise = .TRUE. ! Years
            END IF ! Of um_file % meta % init_steps = reinit period

          END IF ! Of months_in within reinitialisation limits

        END IF  ! Of lcal360 and nftunit etc.

      ELSE IF (um_file % meta % init_steps >  0) THEN ! Regular reinit.

        ! Select files to be reinitialised on this timestep:
        ! First initialisation may be on any timestep, not just 0
        IF (step_compare == um_file % meta % init_start_step) THEN

          um_file % meta % initialise = .TRUE.

          ! Deal with subsequent reinitialisation after end of period
        ELSE IF ((step-1) > um_file % meta % init_start_step .AND.   &
                (um_file % meta % init_end_step <  0 .OR.            &
                (step-1) <= um_file % meta % init_end_step)) THEN

          IF (MOD((step-1) - um_file % meta % init_start_step,       &
              um_file % meta % init_steps) == 0)                     &
              um_file % meta % initialise = .TRUE.

          ! Do not reinitialise files on step 1 if they were
          ! initialised at step 0:
          IF (step == 1 .AND. um_file % meta % init_start_step == 0) &
            um_file % meta % initialise = .FALSE.

          ! Sub model id must be atmos
          um_file % meta % initialise =                  &
              um_file % meta % initialise                &
              .AND. internal_model == atmos_im

        END IF

      ELSE  ! for files not reinitialised, ie. init_steps == 0

        !           Initialise at step 0
        um_file % meta % initialise = step == 0 .AND.  &
        (um_file % meta % init_steps == 0       .OR.   &
        um_file % meta % init_start_step == 0)

      END IF  ! Of init_steps lt, gt or == 0, ie. reinit. type
    ELSE ! Of um_file % meta % is_output_file
      um_file % meta % initialise = .FALSE.
    END IF

    SELECT CASE(i)
      CASE(1)
        lpp = lpp .OR. um_file % meta % initialise
      CASE(2)
        lnc = lnc .OR. um_file % meta % initialise
    END SELECT

    ! Increment pointer to next file
    um_file => um_file % next
  END DO  
END DO  

!
!  1.2   Set switches for general internal models.
!        For coupled models dump related switches can only be set when
!        the last internal model in a submodel has completed its group
!        of timesteps. For coupled models the only safe restart point
!        is at the completion of all groups within a model timestep.

IF (n_internal_model == 1 .OR. MOD(step,steps_per_period) == 0) THEN
  ! if not coupled model, or last step in group
  ! ldump   : Write-up dump on this timestep
  gregdmp: IF (lclimrealyr) THEN
    ! Calculate the elapsed time since time zero to ensure dumping consistency
    ! between NRUNs and CRUNs
    ! DEPENDS ON: time2sec
    CALL time2sec(i_year, i_month, i_day, i_hour, i_minute, i_second, &
       0, 0, elapsed_days, elapsed_secs, lcal360)
    elapsed_steps = (elapsed_days*isec_per_day+elapsed_secs) / secs_per_step
    IF (greg_dump_freq > 0) THEN
      ldump = (MOD(elapsed_steps,greg_dump_freq) == 0)
    END IF
    IF (greg_monthly_dump .AND. i_day == 1 .AND. i_hour == 0 &
       .AND. i_minute == 0 .AND. i_second == 0) THEN
      ldump = .TRUE.
    END IF
    IF (greg_dump_freq <= 0 .AND. .NOT. greg_monthly_dump) THEN
      DO itime=1,dumptimes_len1
        ldump=ldump .OR. (step == dumptimes(itime))
      END DO
    END IF
    exitdump: IF (stepim(1) == target_end_stepim(1) .AND. end_of_run_dump) THEN
      ! make sure there is a dump at the end of the run to ensure
      ! restartability
      ldump = .TRUE.
    END IF exitdump
  ELSE
    IF (dumpfreq >  0) THEN
      ldump=       (MOD(step,dumpfreq)    == 0)
    ELSE
      ldump=.FALSE.
      DO itime=1,dumptimes_len1
        ldump=ldump .OR. (step == dumptimes(itime))
      END DO
    END IF
  END IF gregdmp

  !  LMEAN   : Perform climate-meaning from dumps on this timestep
  lmean = .FALSE.
  IF (dumpfreq >  0 .AND. meanfreq(1) >  0) THEN
    lmean=     (MOD(step,dumpfreq)       == 0)
  END IF

  !  LEXIT   : Check for exit condition on this timestep
  IF (exitfreq >  0) THEN ! Implies climate meaning

    lexit=  ( (MOD(step,exitfreq)    == 0)  .OR.                  &
              (step  >=  target_end_step) )

  ELSE                    ! No climate meaning

    lexit=    (step  >=  target_end_step)

  END IF

  !  lancillary: Update ancillary fields on this timestep

  !   Convert ancillary reference time to days & secs
  ! DEPENDS ON: time2sec
  CALL time2sec(ancil_reftime(1),ancil_reftime(2),                &
                ancil_reftime(3),ancil_reftime(4),                &
                ancil_reftime(5),ancil_reftime(6),                &
                0,0,ancil_ref_days,ancil_ref_secs,lcal360)

  !   Compute offset in timesteps of basis time from ancillary ref.time
  ! DEPENDS ON: tim2step
  CALL tim2step(basis_time_days-ancil_ref_days,                   &
                basis_time_secs-ancil_ref_secs,                   &
                steps_per_period,secs_per_period,ancil_offset_steps)

  IF (ancillary_steps >  0 .AND. step >  0)                         &
    lancillary=(MOD(step+ancil_offset_steps,ancillary_steps) == 0)

END IF     ! Test for non-coupled or coupled + last step in group

!  lboundary    : Update boundary fields on this timestep

IF (boundary_steps >  0)                                          &
    lboundary= (MOD(step+bndary_offset,boundary_steps)  == 0)

!  1.2.1  Set switches for atmosphere timestep
!
IF (internal_model == atmos_im) THEN
  !  1.2.2 Set switches for all cases

  ! Energy correction switches
  ! Set on the assumption that the energy correction is evaluated at the
  ! end of a timestep
  IF (a_energysteps >  0) THEN
    ! true if this is the step on which to do energy calculation
    lenergy =  (MOD(step,a_energysteps) == 0)
    ! True if this is the step after last energy calculation
    lflux_reset = (MOD(step-1,a_energysteps) == 0)
  END IF


  !  lassimilation: Perform data assimilation on this timestep

  lassimilation= (step >  assim_firststepim(atmos_im) .AND.       &
   (step-assim_firststepim(atmos_im)  <=                          &
         assim_stepsim(atmos_im)+assim_extrastepsim(atmos_im)))
  lassimilation=lassimilation .AND. l_iau
  !
  !       Reset Data Time fields in dump header at new analysis time
  !       ( also in history block )
  !       ( This now includes resetting prognostic field LOOKUP headers )
  !       NB: done even if assimilation suppressed by l_fastrun
  !       Set MEANLEV to -1  and ldump to TRUE to force dump of analysis
  !       - otherwise set MEANLEV to 0 (ie. instantaneous)
  !
  iau_resetdt = .FALSE.
  IF (l_iau .OR. l_fastrun) THEN
    IF (step  ==  iau_dtresetstep) iau_resetdt = .TRUE.
  END IF
  IF ( (step == assim_firststepim(atmos_im)+assim_stepsim(atmos_im) &
        .AND. lassimilation) .OR. iau_resetdt) THEN
    DO i=21,27
      a_fixhd(i)=a_fixhd(i+7)
    END DO
    DO i=1,6
      model_data_time(i)=a_fixhd(20+i)
    END DO
    DO i=1,a_prog_lookup
      a_lookup(lbyrd ,i)=a_fixhd(21)
      a_lookup(lbmond,i)=a_fixhd(22)
      a_lookup(lbdatd,i)=a_fixhd(23)
      a_lookup(lbhrd ,i)=a_fixhd(24)
      a_lookup(lbmind,i)=a_fixhd(25)
      a_lookup(lbsecd,i)=a_fixhd(26)
    END DO
    meanlev=-1
  ELSE
    meanlev=0
  END IF
  !
  !       Suppress assimilation using l_fastrun if necessary
  lassimilation=lassimilation .AND. .NOT. l_fastrun

  ! The final assimilation flag has been set so check it's safe to continue:
  IF (lassimilation .AND. cartesian_grid) THEN
    icode = 10
    cmessage = 'Data assimilation does not support Cartesian grids'
    CALL ereport(RoutineName, icode, cmessage)
  END IF

  ! If using radiation, set whether this is a radiation timestep.
  ! If not using radiation, do nothing; leave the l_rad_step flags in
  ! their initialised states (false, set in set_rad_steps_mod).
  IF (l_radiation)  CALL set_l_rad_step(step)

END IF ! Test for atmos_im
!
! ----------------------------------------------------------------------
!  2. Set STASHflags to activate diagnostics this timestep
!
!     IS is section number
!     IE is item number within section
!     IL is item number within STASHlist
!     II is counter within given section/item sublist for repeated items
!     IM is cumulative active item number within section

!   Clear all STASHflags

DO IS=0,nsects
  DO im=0,nitems
    sf     (im,IS) = .FALSE.
    sf_calc(im,IS) = .FALSE.
  END DO
END DO

!   Loop over all items in STASHlist, enabling flags for diagnostics
!   which are active this step --
!     note that atmosphere/ocean diagnostics must be discriminated

DO il =1,totitems
  modl=stlist(st_model_code  ,il)
  IS  =stlist(st_sect_no_code,il)
  im  =stlist(st_item_code   ,il)
  ! Skip diagnostics which are not active
  IF (stlist(st_proc_no_code,il) == 0) GO TO 200
  ! Skip diagnostics which don't correspond to the submodel this timestep
  IF (internal_model /= modl) GO TO 200

  !     STASHflag is off by default.
  !     But reset ...

  !       ... if required every timestep between start and end

  IF (stlist(st_freq_code,il) == 1) THEN
    IF (step >= stlist(st_start_time_code,il) .AND.                &
       (step <= stlist(st_end_time_code,il) .OR.                   &
        stlist(st_end_time_code,il) == st_infinite_time))         &
    THEN
      sf(im,IS)=.TRUE.
      sf(0 ,IS)=.TRUE.
    END IF

    !       ... if required at specified times and this is one of them

  ELSE IF (stlist(st_freq_code,il) <  0) THEN
    ntab=-stlist(st_freq_code,il)
    DO it=1,nsttims
      IF (sttabl(it,ntab) == st_end_of_list) EXIT
      IF (step == sttabl(it,ntab)) THEN
        sf(im,IS)=.TRUE.
        sf(0 ,IS)=.TRUE.
      END IF
    END DO

    !       ... if required every N timesteps and this is one of them

  ELSE IF (stlist(st_freq_code,il) >  0) THEN
    IF (MOD((step-stlist(st_start_time_code,il)),                &
              stlist(st_freq_code,il)) == 0 .AND.                  &
         step >= stlist(st_start_time_code,il) .AND.               &
        (step <= stlist(st_end_time_code,il) .OR.                  &
         stlist(st_end_time_code,il) == st_infinite_time))        &
    THEN
      sf(im,IS)=.TRUE.
      sf(0 ,IS)=.TRUE.
    END IF
  END IF

  ! We now want just a subset of the above where we only want to calculate
  IF (sf(im,IS)) THEN
    ! See whether we are performing a replace.
    IF ( stlist(st_proc_no_code,il) == st_replace_code ) THEN
      ! We only want to know about replaces which are not dependent on other
      ! records.  If positive this STASH request requires a calculation.
      IF ( stlist(st_input_code,il) > 0 ) THEN
        sf_calc(im,IS)  = .TRUE.
        sf_calc(0 ,IS)  = .TRUE.
      END IF
      ! Everything not a replace we just let through.
    ELSE
      sf_calc(im,IS)  = .TRUE.
      sf_calc(0 ,IS)  = .TRUE.
    END IF

  END IF

  !     Next item

  200 CONTINUE
END DO
! ----------------------------------------------------------------------
!  3. If errors occurred in setting STASHflags, set error return code
IF (icode == 1) cmessage='SETTSCTL: STASH code error'
9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
! ----------------------------------------------------------------------
END SUBROUTINE settsctl

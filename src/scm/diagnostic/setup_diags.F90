! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Initialise the SCM output diagnostic system

SUBROUTINE setup_diags                                                        &
  ( row_length, rows, bl_levels, sm_levels                                    &
  , st_levels, land_points, ntiles, nsmax, n_vis_thresh                       &
  , cloud_levels, total_nsteps  &
  , timestep, full_daysteps, a_sw_radstep_prog, a_sw_radstep_diag             &
  , ntrad1, daycount, stepcount, nscmdpkgs, l_scmdiags, scmop )

USE conversions_mod, ONLY: rsec_per_day

USE scm_cntl_mod, ONLY:                                                     &
    scm_nml, scmop_type, i64, maxndomprof, maxnstreams, main_diag_switch    &
  , netcdf_chunksize, diags, listlength, inot_written, lsname, DEFAULT      &
  , strm_switch, strm_unit, strm_format, strm_filename, strm_dumpstep       &
  , strm_rejectlist, strm_acceptlist                                        &
  , strm_heed_hardwired, strm_heed_acceptlist, strm_heed_rejectlist         &
  , l_scmdiag_gen,   l_scmdiag_rad, l_scmdiag_bl,  l_scmdiag_surf           &
  , l_scmdiag_land,  l_scmdiag_sea, l_scmdiag_lsp, l_scmdiag_conv           &
  , l_scmdiag_lscld, l_scmdiag_pc2, l_scmdiag_gwd, l_scmdiag_forc           &
  , l_scmdiag_incs,  l_scmdiag_mlsnow, l_scmdiag_convss

USE nlsizes_namelist_mod, ONLY: model_levels

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE s_scmop_mod, ONLY: DoNotWrite                                           &
  , t_inst, t_avg, t_max, t_min, t_acc, t_div, t_mult, t_acc_div            &
  , t_acc_mult, t_const, only_radsteps                                      &
  , d_sl, d_soilt, d_bl, d_wet, d_all, d_soilm, d_tile, d_vis, d_point      &
  , d_allxtra, d_land, d_cloud, d_mlsnow, default_streams                   &
  , SCMDiag_gen, SCMDiag_rad, SCMDiag_bl,   SCMDiag_surf,  SCMDiag_land     &
  , SCMDiag_sea, SCMDiag_lsp, SCMDiag_conv, SCMDiag_lscld, SCMDiag_pc2      &
  , SCMDiag_forc, SCMDiag_incs, SCMDiag_gwd, SCMDiag_MLSnow, SCMDiag_convss
USE scmoutput_mod, ONLY: scmoutput

USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE


! Description:
!   Perform all pre-timestepping initialisation for the output
!   diagnostic system, i.e. initialise SCMop. In particular, the
!   DIAGS namelist is read in, and the output streams and domain
!   profiles are set up. The output files are -not- initialised
!   here since this requires information which is not available
!   until the end of at least one timestep (see dump_streams_init).
!   All arguments are INTENT In, apart from main_diag_switch,
!   the l_SCMdiags logicals and SCMop which are INTENT Out.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran90

TYPE(SCMop_type) :: SCMop ! Out The derived-type structure
                          ! containing all the diagnostic
                          ! information

INTEGER ::          &
  row_length        &
, rows              &
, bl_levels         &
, sm_levels         &
, st_levels         &
, land_points       &
, ntiles            &
, nsmax             &
, n_vis_thresh      &
, cloud_levels      &
, total_nsteps      &
, full_daysteps     &
, ntrad1            &
, a_sw_radstep_prog &
, a_sw_radstep_diag


REAL :: timestep

! SCMop has pointers which will point to stepcount and daycount
! so that it can know which step it's on
INTEGER, TARGET ::  &
  daycount          &
, stepcount

INTEGER ::          &
  i,j               &! Counters.
, pos               &! Current position in string.
, last_comma        &! Position of last comma.
, n_words            ! No. of comma separated words found.

CHARACTER (LEN=5) :: suffix        ! Holds frmt-dependent filename suffix
CHARACTER (LEN=listlength) :: list ! Temporarily holds either an accept or
                                   ! reject list

CHARACTER (LEN=lsname) :: NAME     ! Temporarily holds a (short) name of a
                                   ! diagnostic
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER (LEN=*), PARAMETER ::  RoutineName = 'SETUP_DIAGS'
! The ASCII value is used to specify a backslash to avoid the problem of
! backslash having implicit meaning for some compilers.
CHARACTER (LEN=1), PARAMETER :: Backslash = ACHAR(92)

INTEGER :: istatus                 ! Error code

! Some values which will be used to scale diagnostics
REAL :: ntimestepsperday,oneKsecday

INTEGER :: nSCMdpkgs             ! No of SCM diagnostics packages
                                 ! (set in Scm_Main)
LOGICAL :: l_SCMdiags(nSCMdpkgs) ! Logicals array for SCM
                                 ! diagnostics packages

! Dr Hook
!=============================================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

SCMop%first_pass     =  .TRUE.       ! Diagnostic list not yet finalised
SCMop%on             =  .FALSE.      ! Timestepping not yet started
SCMop%daycount       => daycount
SCMop%stepcount      => stepcount
SCMop%full_daysteps  =  full_daysteps
SCMop%ntrad          =  a_sw_radstep_diag
SCMop%ntrad1         =  ntrad1
SCMop%maxnentries    =  0            ! No memory allocated for diagnostics yet
SCMop%nentries       =  0            ! No diagnostics requested
                                     ! to be stored yet
SCMop%n_output       =  0            ! No diagnostics requested
                                     ! to be output yet
SCMop%nSCMoutput     =  0            ! No calls to SCMoutput yet
SCMop%num_substeps   =  1            ! Expected number of sub-steps
SCMop%substep_number =  0            ! No substepping started yet

! Initialise one element of the domain profiles to flag
! undefined profiles
DO i=1, maxndomprof
  SCMop%d_lev1(i) = -1
END DO

! Allocate memory for streams
ALLOCATE(SCMop%strm(maxnstreams))

! Allocate some token space to the diag_mem array (ensures the
! SIZE(diag_mem) operation in SCMoutput doesn't fall over the
! first time)
ALLOCATE(SCMop%diag_mem(1,2))


! Modify default values of some namelist variables...

! Diagnostic system on by default (cannot be set by data
! statement because it is intent out).
main_diag_switch = 1

! Default NetCDF file chunk size
#if defined(NEC)
netcdf_chunksize = 32000000
#else
netcdf_chunksize = 8192
#endif

! l_SCMdiags logicals (intent Out) are also set to default values
l_SCMdiag_gen   = .TRUE.
l_SCMdiag_rad   = .FALSE.
l_SCMdiag_bl    = .FALSE.
l_SCMdiag_surf  = .FALSE.
l_SCMdiag_land  = .FALSE.
l_SCMdiag_sea   = .FALSE.
l_SCMdiag_lsp   = .FALSE.
l_SCMdiag_conv  = .FALSE.
l_SCMdiag_lscld = .FALSE.
l_SCMdiag_PC2   = .FALSE.
l_SCMdiag_gwd   = .FALSE.
l_SCMdiag_MLSnow= .FALSE.
l_SCMdiag_forc  = .FALSE.
l_SCMdiag_incs  = .FALSE.
l_SCMdiag_convss= .FALSE.

! Stream 1 will be on by default
strm_switch(1) = 1

! The default output units for each stream (37 upwards)...
DO i=1, maxnstreams
  strm_unit(i) = 36+i
END DO

! Read the DIAGS namelist scm_nml
OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus, STATUS='old',       &
     ACTION='READ', IOMSG=iomessage)
IF (istatus == 0) THEN
  READ(10, diags, IOSTAT=istatus, IOMSG=iomessage)
  CLOSE(10)
  IF (istatus /= 0) THEN
    WRITE(umMessage,*)'------------------------------------------------'
    CALL umPrint(umMessage,src='setup_diags')
    WRITE(umMessage,*)'SETUP_DIAGS WARNING: DIAGS namelist not found in'
    CALL umPrint(umMessage,src='setup_diags')
    WRITE(umMessage,*)'namelist file. Diagnostic system will operate on'
    CALL umPrint(umMessage,src='setup_diags')
    WRITE(umMessage,*)'defaults.'
    CALL umPrint(umMessage,src='setup_diags')
    WRITE(umMessage,*)'Error: ',TRIM(iomessage)
    CALL umPrint(umMessage,src='setup_diags')
    WRITE(umMessage,*)'------------------------------------------------'
    CALL umPrint(umMessage,src='setup_diags')
  END IF
ELSE
  WRITE(umMessage,*) "-------------------------------------------------"
  CALL umPrint(umMessage,src='setup_diags')
  WRITE(umMessage,*) " Error opening " // TRIM(ADJUSTL(scm_nml))
  CALL umPrint(umMessage,src='setup_diags')
  WRITE(umMessage,*) " IO status = ", istatus
  CALL umPrint(umMessage,src='setup_diags')
  WRITE(umMessage,*) " Error: ", TRIM(iomessage)
  CALL umPrint(umMessage,src='setup_diags')
  WRITE(umMessage,*) " Using default diagnostic settings"
  CALL umPrint(umMessage,src='setup_diags')
  WRITE(umMessage,*) "-------------------------------------------------"
  CALL umPrint(umMessage,src='setup_diags')
END IF

! No point doing any more if the diagnostic system has been
! switched off
IF (main_diag_switch == 0) THEN
  IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! NetCDF I/O is done using C routines which don't utilise Fortran
! unit numbers. Set to -1 so that this is reflected in the stream
! summary print-out.
DO i=1, maxnstreams
  IF (strm_format(i) == 4) THEN
    strm_unit(i)=-1
  END IF
END DO

! Set up the output streams...
DO i=1, maxnstreams
  SCMop%strm(i)%switch  = strm_switch(i)     ! Stream on/off
  SCMop%strm(i)%op_unit = strm_unit(i)       ! Output unit
  SCMop%strm(i)%FORMAT  = strm_format(i)     ! Output format

  ! Output filename...
  IF (strm_filename(i) /= DEFAULT) THEN
    SCMop%strm(i)%filename = strm_filename(i)
  ELSE

    ! The filename suffix will be format dependent
    suffix='.out'     ! Unrecognised format
    IF (SCMop%strm(i)%FORMAT == 0) THEN
      suffix = '.dat'   ! Old Wave format
    ELSE IF (SCMop%strm(i)%FORMAT == 1) THEN
      suffix = '.txt'   ! Easy to read format
    ELSE IF (SCMop%strm(i)%FORMAT == 2) THEN
      suffix = '.fssi'  ! Format for FSSI
    ELSE IF (SCMop%strm(i)%FORMAT == 3) THEN
      suffix = '.dat'   ! New Wave format
    ELSE IF (SCMop%strm(i)%FORMAT == 4) THEN
      suffix = '.nc'    ! NetCDF format
    END IF
    WRITE(SCMop%strm(i)%filename,'(A,I2.2,A)')'stream',i,suffix
  END IF

  ! If the dumping period is positive, take this as a number
  ! of timesteps. If it is negative take it as a number of
  ! seconds. If it is zero, replace it with the number of
  ! timesteps in the whole run (i.e. one dump on the
  ! last step).
  IF (strm_dumpstep(i) == 0) strm_dumpstep(i)=total_nsteps

  IF (strm_dumpstep(i) >  0) THEN

    ! The dumping period has been specified in timesteps
    SCMop%strm(i)%dump_step = strm_dumpstep(i)
  ELSE IF (strm_dumpstep(i) <  0) THEN

    ! The dumping period has been specified in seconds, convert
    ! to nearest non-zero number of timesteps
    SCMop%strm(i)%dump_step = MAX(NINT(-strm_dumpstep(i)/timestep),1)
  END IF

  ! Accept diagnostics sent to stream by respective SCMoutput
  ! parameter? (0=no)
  SCMop%strm(i)%heed_hardwired = strm_heed_hardwired(i)

  ! Accept diagnostics sent to stream by namelist? (0=no)
  SCMop%strm(i)%heed_acceptlist = strm_heed_acceptlist(i)

  ! Reject diagnostics sent to stream by any method? (0=no)
  SCMop%strm(i)%heed_rejectlist = strm_heed_rejectlist(i)

  ! No. of diagnostics sent to this stream (none until at least
  ! first call to SCMoutput)
  SCMop%strm(i)%n_output = 0

  ! We now want to extract the individual diagnostic names from
  ! the comma-separated lists in strm_acceptlist(i) and
  ! strm_rejectlist(i). We will scan through both strings
  ! twice, first to count the number of names so the correct
  ! space can be allocated (j=1 and j=3) and then again to put
  ! the names into the allocated arrays (j=2 and j=4).
  DO j=1, 4

    IF (j <= 2) THEN
      ! First two loops work on the acceptlist
      list = strm_acceptlist(i)
    ELSE
      ! Next two are indepedent of the first and work
      ! on the reject list
      list = strm_rejectlist(i)
    END IF

    ! We're going to scan through the list of
    ! comma-separated words one character at a time
    ! starting from the beginning.
    ! Make sure there's a trailing comma on the list
    list = TRIM(list)//','
    pos = 1               ! Current position
    last_comma = 0        ! Pos'n of last comma found
    n_words = 0           ! No. of words found so far

    ! Loop through the characters of the string
    DO WHILE (pos <= LEN_TRIM(list))

      IF (list(pos:pos) == Backslash) THEN
        ! This character is a backslash - it escapes the
        ! next character causing us to jump forward one.
        pos = pos+1

      ELSE IF (list(pos:pos) == ' ' .AND.                        &
               last_comma == pos-1) THEN
        ! This is a space trailing a comma - ignore by
        ! "moving" the last comma forward one
        last_comma = pos
      ELSE IF (list(pos:pos) == ',') THEN

        ! We've found a (non-escaped) comma. Extract the name
        ! enclosed by this comma and the last.
        NAME = list(last_comma+1:pos-1)

        ! Update the position of the last comma found
        last_comma = pos

        ! Count the number of names we've found
        n_words = n_words+1

        ! On the second pass for each list, fill the allocated
        ! arrays with the names.
        IF (j == 2) SCMop%strm(i)%accept_list(n_words) = NAME
        IF (j == 4) SCMop%strm(i)%reject_list(n_words) = NAME
      END IF

        ! Update the position
      pos = pos+1

    END DO ! over characters in a list

    ! After the first pass for each list, allocate the arrays
    ! into which the names will be put on the second pass.
    IF (j == 1) ALLOCATE(SCMop%strm(i)%accept_list(n_words))
    IF (j == 3) ALLOCATE(SCMop%strm(i)%reject_list(n_words))

  END DO ! j=1,4

   ! FSSI-format output files require a dump for the first
   ! timestep as well as at the end of every dumping period. The
   ! way this has been coded in routine dump_streams means that
   ! if the dumping period is equal to one then this additional
   ! dump will not be produced and the file will be one line
   ! shorter than expected. I am assured that the SSFM will never
   ! be run with a dumping period of one, but it's worth printing
   ! a warning in case it ever is.
  IF (SCMop%strm(i)%FORMAT == 2 .AND.                             &
       SCMop%strm(i)%dump_step == 1) THEN
    WRITE(umMessage,*)'---'
    CALL umPrint(umMessage,src='setup_diags')
    WRITE(umMessage,*)'setup_diags WARNING: a FSSI-format output file'
    CALL umPrint(umMessage,src='setup_diags')
    WRITE(umMessage,*)'has been requested with a dumping period of'
    CALL umPrint(umMessage,src='setup_diags')
    WRITE(umMessage,*)'one. This file will not have an additional dump'
    CALL umPrint(umMessage,src='setup_diags')
    WRITE(umMessage,*)'for the first timestep and so may be one line'
    CALL umPrint(umMessage,src='setup_diags')
    WRITE(umMessage,*)'shorter than expected.'
    CALL umPrint(umMessage,src='setup_diags')
    WRITE(umMessage,*)'---'
    CALL umPrint(umMessage,src='setup_diags')
  END IF

END DO                     ! i=1,maxnstreams


! Set up the domain profiles
! DEPENDS ON: define_domprof
CALL define_domprof(d_all,'all_levels', &
     1,row_length,                      &! Every column
     1,rows,                            &! Every row
     1,model_levels,SCMop)               ! Every level

! DEPENDS ON: define_domprof
CALL define_domprof(d_sl,'single_level',&
     1,row_length,                      &! Every column
     1,rows,                            &! Every row
     1,1,SCMop)                          ! Only the lowest level

! DEPENDS ON: define_domprof
! From vn10.2 Wet levels is no longer different to model_levels.
! This domain profile may want to be removed in future.
CALL define_domprof(d_wet,'wet_levels',            &
     1,row_length,                                 &
     1,rows,                                       &
     1,model_levels,SCMop)

! DEPENDS ON: define_domprof
CALL define_domprof(d_bl,'bl_levels',              &
     1,row_length,                                 &
     1,rows,                                       &
     1,bl_levels,SCMop)

! DEPENDS ON: define_domprof
CALL define_domprof(d_soilm,'soil_moist_levels',   &
     1,row_length,                                 &
     1,rows,                                       &
     1,sm_levels,SCMop)

! DEPENDS ON: define_domprof
CALL define_domprof(d_soilt,'soil_temp_levels',    &
     1,row_length,                                 &
     1,rows,                                       &
     1,st_levels,SCMop)

! DEPENDS ON: define_domprof
IF ( land_points > 0 ) THEN
  CALL define_domprof(d_land,'land_points',        &
       1,land_points,                              &! This is a fudge
       1,1,                                        &
       1,1,SCMop)
END IF

! DEPENDS ON: define_domprof
IF ( land_points > 0 .AND. nsmax > 0 ) THEN
  CALL define_domprof(d_mlsnow,'MLSnow_points',    &
       1,land_points,                              &
       1,ntiles,                                   &
       1,nsmax,SCMop)
END IF

! DEPENDS ON: define_domprof
CALL define_domprof(d_allxtra,'all_levs_plus1',    &
     1,row_length,                                 &
     1,rows,                                       &
     1,model_levels+1,SCMop)

IF ( land_points > 0 ) THEN
  ! DEPENDS ON: define_domprof
  CALL define_domprof(d_tile,'tile_types',         &
       1,land_points,                              &
       1,1,                                        &
       1,ntiles,SCMop)
END IF

! DEPENDS ON: define_domprof
CALL define_domprof(d_point,'single_point',1,1,1,1,1,1,SCMop)

! DEPENDS ON: define_domprof
CALL define_domprof(d_vis,'vis_thresholds',        &
     1,row_length,                                 &
     1,rows,                                       &
     1,n_vis_thresh,SCMop)

! DEPENDS ON: define_domprof
CALL define_domprof(d_cloud,'cloud_levels',        &
     1,row_length,                                 &
     1,rows,                                       &
     1,cloud_levels,SCMop)

! Some diagnostics are just constants which will be combined
! with other, variable diagnostics. May as well just call
! SCMoutput for them once then, rather than every timestep.

! Note: The first input to scmoutput must be an array, not scalar,
! even if its size is 1, because Fortran doesn't allow passing scalar
! arguments to array dummy arguments.  To pass in the following constants,
! they must be wrapped in the array constructor [].

CALL scmoutput([rsec_per_day], 'sec_day'                           &
  , 'No. of seconds per day','-'                                   &
  , t_const, d_point, DoNotWrite(default_streams), '', RoutineName)

oneKsecday=1000.0*rsec_per_day

CALL scmoutput([oneKsecday], 'oneKsecday'                          &
  , '1000.*rsec_day','-'                                           &
  , t_const, d_point, DoNotWrite(default_streams), '', RoutineName)

ntimestepsperday = rsec_per_day/timestep

CALL scmoutput([ntimestepsperday], 'ntspday'                       &
  , 'No. of timesteps per day', '-'                                &
  , t_const, d_point, DoNotWrite(default_streams), '', RoutineName)

! Place details from read in diagnostics packages logicals
! into logicals array
! Note: general is default true, the rest default false
l_SCMdiags(SCMdiag_gen)   = l_SCMdiag_gen
l_SCMdiags(SCMdiag_rad)   = l_SCMdiag_rad
l_SCMdiags(SCMdiag_bl)    = l_SCMdiag_bl
l_SCMdiags(SCMdiag_surf)  = l_SCMdiag_surf
l_SCMdiags(SCMdiag_land)  = l_SCMdiag_land
l_SCMdiags(SCMdiag_sea)   = l_SCMdiag_sea
l_SCMdiags(SCMdiag_lsp)   = l_SCMdiag_lsp
l_SCMdiags(SCMdiag_conv)  = l_SCMdiag_conv
l_SCMdiags(SCMdiag_lscld) = l_SCMdiag_lscld
l_SCMdiags(SCMdiag_PC2)   = l_SCMdiag_PC2
l_SCMdiags(SCMdiag_gwd)   = l_SCMdiag_gwd
l_SCMdiags(SCMdiag_MLSnow)= l_SCMdiag_MLSnow
l_SCMdiags(SCMdiag_forc)  = l_SCMdiag_forc
l_SCMdiags(SCMdiag_incs)  = l_SCMdiag_incs
l_SCMdiags(SCMdiag_convss)= l_SCMdiag_convss

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE setup_diags

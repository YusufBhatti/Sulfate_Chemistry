! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Subroutine: AC
!
!    Purpose : Main control routine for AC Scheme
!     This makes the analysis correction to model fields at the start
!   of each timestep of model integration.  The corrections are
!   calculated at analysis grid points,then interpolated to the model
!   grid before being added to the model fields. The corrections are
!   weighted sums of increments (ob-model) obtained at analysis grid
!   points. The weights depend on the distance between observation and
!   grid point, Observational error and the timeliness of the
!   observation.  Observations are analysed sequentially by observation
!   type in the order given by the array lact in ac_control_mod. A
!   vertical profile of observational increments and errors at
!   observation points are obtained before the weighted sums are
!   calculated. Horizontal filtering of increments on the model grid may
!   be done,and hydrostatic and geostrophic increments are calculated on
!   the model grid if appropriate.mean and rms diagnostics may be
!   calculated.
!
!    This section of code is now only used for the assimilation of
!    rainfall data from MOPS.
!
!
!    Programming standard: Unified Model Documentation Paper No. 3
!
!    Project Task : P3
!
!    Arguments and declarations:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE ac_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='AC_MOD'

CONTAINS

SUBROUTINE ac (                                                   &
  row_length, p_rows, bl_levels,                                  &
  p_field,                                                        &
  timestep_no, timestep,                                          &
  exner, pstar, pressure,                                         &
  theta, rh, qcl, qcf,                                            &
  conv_cld, ls_rain, ls_snow, conv_rain, conv_snow,               &
  layer_cloud, d_theta_dt_conv,d_theta_dt_ls, rhcrit,             &
  obs_flag,obs,                                                   &
  stindex,                                                        &
  stlist, len_stlist, si, sf,                                     &
  stashwork, stash_levels,                                        &
  num_stash_levels, stash_pseudo_levels, num_stash_pseudo,        &
  lambda_p,phi_p,                                                 &
  icode, cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr
USE UM_ParVars
USE UM_ParCore, ONLY: mype, nproc
USE comobs_mod, ONLY: nobtypmx, obs_info
USE ac2_mod, ONLY: ac2
USE getobs_mod, ONLY: getobs
USE hmrtorh_mod, ONLY: hmrtorh
USE rdobs_mod, ONLY: rdobs
USE ac_control_mod
USE nlsizes_namelist_mod, ONLY: model_levels
USE num_obs_mod, ONLY: ac_num_obs_max, ac_tot_obs_size_max

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Imported arguments (INTENT=IN):
INTEGER ::     bl_levels               ! Bdy layer levels
INTEGER ::     row_length              ! Number of points on row
INTEGER ::     p_rows                  ! Number of rows (pstar)
INTEGER ::     p_field                 ! Number of points in
                                       ! mass field
INTEGER ::     timestep_no             ! Current model timestep
REAL ::        conv_cld(p_field, model_levels)  ! conv cld amount
REAL ::        ls_rain(p_field)        ! large scale rain rate
REAL ::        ls_snow(p_field)        ! large scale snow rate
REAL ::        conv_rain(p_field)      ! convective rain rate
REAL ::        conv_snow(p_field)      ! convective snow rate
                                       ! above rates diagnostic
REAL ::        timestep                ! Timestep in seconds
REAL ::        rhcrit(model_levels)    ! Critical rh array
                                       ! for cloud
REAL ::        d_theta_dt_conv(p_field,model_levels)
                                ! convective latent heating rate
REAL ::        d_theta_dt_ls(p_field,model_levels)
                                ! large scale latent heating rate
! Stash variables:
REAL ::        stashwork(*)

INTEGER ::     len_stlist              ! Length of STLIST
INTEGER ::     num_stash_levels        ! No of Stash levels lists
INTEGER ::     num_stash_pseudo        ! No of Stash pseudo lists
INTEGER ::     stindex(2, *)
INTEGER ::     stlist(len_stlist, *)
INTEGER ::     si(*)                   ! Stash Index
INTEGER ::     stash_levels(num_stash_levels +1, *)
                                       ! Levels lists
INTEGER ::     stash_pseudo_levels(num_stash_pseudo +1, *)
                                       ! Pseudo lists

LOGICAL ::     sf(*)                   ! Stash Flags

! Import/export arguments (INTENT=INOUT):
REAL ::        exner(p_field, model_levels)        ! exner on theta levels
REAL ::        pressure(p_field, model_levels)     ! p  on theta levels
REAL ::        layer_cloud (p_field, model_levels) ! as name says
REAL ::        pstar(p_field)                      ! Prognostic variable pstar
REAL ::        theta(p_field, model_levels)        ! Prognostic variable theta
REAL ::        rh(p_field, model_levels)           ! Prognostic variable hmr
REAL ::        qcl(p_field, model_levels)          ! Prognostic variable qcl
REAL ::        qcf(p_field, model_levels)          ! Prognostic variable qcf

! Exported arguments (INTENT=OUT):
INTEGER ::     icode                   ! Non zero for failure

CHARACTER(LEN=errormessagelength) :: cmessage           ! Reason for failure

! These are used in variable resolution runs
REAL :: lambda_p(1-halo_i:row_length+halo_i)
REAL :: phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)

!  For load balancing
INTEGER :: lenob_total

INTEGER :: obs_flag(ac_num_obs_max)    ! Observation flags
INTEGER :: obs_no(ac_num_obs_max)      ! Observation numbers
REAL  ::   obs(ac_tot_obs_size_max)    ! Observation data
! Local variables:
INTEGER ::     tndv                    ! Total no of data values
                                       ! for obs in assimilation
                                       ! time window
INTEGER ::     tnobs                   ! Total no of obs in
                                       ! assimilation time window
INTEGER ::     navt                    ! Analysis variable type
                                       ! 1/2/3/4/5/6
                                       ! p*/theta/winds/rh/precip
                                       ! /tracer
INTEGER ::     wklen                   ! Length of vertical
                                       ! dimension of array for
                                       ! derived increments made
                                       ! by HYDRST, GEOSTR &
                                       ! WINDBAL.
INTEGER ::     nptobt                  ! Offset to first obs of
                                       ! each type in group
INTEGER ::     no_wt_levs              ! No of weight levels
                                       ! for group of obs types
INTEGER ::     no_anal_levs            ! No of analysis levels
                                       ! for group of obs types
INTEGER ::     no_anal_var             ! No of variables being
                                       ! analysed (2 for winds)
INTEGER ::     lenmg                   ! Length of model grid
INTEGER ::     lenag                   ! Length of analysis grid
INTEGER ::     lenob                   ! No of obs in group
INTEGER ::     lenobt                  ! No of obs for obs type
INTEGER ::     itotal                  ! Total no of iterations so
                                       ! far in whole assimilation
INTEGER ::     idiag                   ! Diagnostic control for
                                       ! this timestep
INTEGER ::     total_no_iters          ! Total no of iterations
                                       ! to be done this timestep
                                       ! (including diagnostic
                                       ! only iterations)
INTEGER ::     iter_no                 ! Iteration number
                                       ! excluding diagnostic
                                       ! only iterations
INTEGER ::     kiter                   ! Iteration no used in
                                       ! diagnostic output.
INTEGER ::     iactf                   ! First obs type in group
INTEGER ::     iactl                   ! Last  obs type in group
INTEGER ::     inobs                   ! No of obs for this type
INTEGER ::     inobsdim                !  "  " (for dimensioning)
INTEGER ::     inc_type                ! Pointer into MODEL_INCR
                                       ! if FI used
INTEGER ::     itnobs(nobtypmx)        ! No of obs in each group
                                       ! assimilated on timestep
INTEGER ::     i             !? tracer
INTEGER ::     ipt_tracer    !? tracer ! pointer to one tracer
                                       ! within full array TRACERS
INTEGER ::     mode_hanal              ! FI or HORINF for this
                                       ! group
INTEGER ::     jiter                   ! Loop conter for iteration
INTEGER ::     j
INTEGER ::     jgroup,jjgroup          ! Loop counter over groups
INTEGER ::     jact,jjact              ! Loop counter over obs
                                       ! types in group
INTEGER ::     igroup                  ! Group index in group
                                       ! dependent arrays
INTEGER ::     hmrmode                 ! Mode for action in
                                       ! HMRTORH
                                       ! 1 = Convert RH to HMR
                                       ! 2 = Convert HMR to RH
INTEGER ::     vmode                  ! Mode for action in
                                      ! VISTOLV and AEROTOLA
                                      ! 1 = Convert VIS to LOGVIS
                                      ! 2 = Convert LOGVIS to VIS
INTEGER ::     istat                  ! status for GCOM

REAL ::        assm_time               ! Assimilation time
                                       ! Relative to start

LOGICAL ::     lwind                   ! Indicator if group
                                       ! has wind obs types
LOGICAL ::     dg_this                 ! ) Switches to control
LOGICAL ::     dg_between              ! ) whether diagnostics
LOGICAL ::     dg_end                  ! ) required this timestep
                                       ! ) , between iterations
                                       ! ) and at end of timestep
LOGICAL ::     dg_only                 ! Indicator if diagnostic
                                       ! only iteration
LOGICAL ::     old_ldiagac             ! Original value of LDIAGAC
                                       ! stored during timestep
LOGICAL ::     l_ac_timestep           ! Switch to control if
                                       ! iteration done this
                                       ! timestep

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='AC'

! Subroutine calls:

! Extra bits:
DATA itotal/0/
SAVE itotal

!- End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! If no valid ACOB files available then jump straight out
IF (no_obs_files == 0) THEN
  WRITE(umMessage,*) '***WARNING: NO VALID ACOB FILES TO ASSIMILATE'
  CALL umPrint(umMessage,src='ac-ac1a')
  WRITE(umMessage,*) 'TIMESTEP_NO TIMESTEP ',timestep_no,timestep
  CALL umPrint(umMessage,src='ac-ac1a')
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

!
!  1.0 Call RDOBS to read in observations either from cache file or
!      from observation files depending on timestep number.
CALL rdobs (no_obs_files, timestep_no, timestep, obs, obs_flag,   &
            tndv, tnobs,                                          &
            lambda_p,phi_p,p_rows,row_length,                     &
            icode, cmessage)

IF (mype == 0 .AND. printstatus >= prstatus_normal) THEN
  WRITE(umMessage,*)'RDOBS done, TIMESTEP_NO = ',timestep_no
  CALL umPrint(umMessage,src='ac-ac1a')
END IF

! Check for errors (drop out if one has occured):
IF (icode  >   0) GO TO 999

!  2.0 Set up Diagnostic Iteration Control from MACDIAG(TIMESTEP_NO)
!     MACDIAG = 0   no diagnostics this timestep
!     MACDIAG = 32  to get diagnostics on each iteration this timestep
!     MACDIAG = +64 for extra diagnostic only iteration between
!                   iterations
!     MACDIAG = + 8 for extra diagnostic only iteration at end
!                   of timestep

idiag = macdiag(1 + MOD(timestep_no -1, modeacp))

! Remember LDIAGAC for re-use at end of timestep
old_ldiagac = ldiagac

! Diagnostics required on each iteration this timestep ?
dg_this    = MOD(idiag/32, 2)  ==  1

! Extra diagnostic only iteration required between iterations ?
dg_between = dg_this .AND. MOD(idiag/64, 2)  ==  1

! Extra diagnostic only iteration required at end of timestep ?
dg_end     = dg_this .AND. MOD(idiag/8, 2)  ==  1

! Any diagnostics required this timestep ?
ldiagac = ldiagac .AND. (dg_this .OR. dg_between .OR. dg_end)

!  3.0 Set up iteration control for ac
!     NO_ITERATIONS : no of iterations for each group.
!     INTERVAL_ITER : interval (in timesteps) between iterations
!                   : for each group.

! Determine total no of iterations.
total_no_iters = 0

DO jgroup = 1, n_groups
  igroup = group_index(jgroup)

  IF (MOD(timestep_no -1, def_interval_iter(igroup))  ==  0) THEN
    ! Do iteration(s) from group JGROUP this timestep.
    total_no_iters = MAX(total_no_iters,def_no_iterations(igroup))

  END IF
END DO

IF (dg_between) total_no_iters = total_no_iters*2
IF (dg_end)     total_no_iters = total_no_iters+1

!  4.0 Start loop iterating ac.
iter_no = 0   !  Iteration No excluding Diagnostic only iters.

! Loop over total no of iterations.
DO jiter =1, total_no_iters
  ! Diagnostics only this iteration
  dg_only = (dg_between .AND. MOD(jiter,2)  ==  1) .OR.           &
            (dg_end     .AND. jiter  ==  total_no_iters)

  IF (.NOT. dg_only) THEN
    itotal = itotal +1
    iter_no = iter_no +1

  END IF

  ! KITER is used for printed output.
  kiter = iter_no

  IF (dg_only) kiter = 0

  DO j = 1, nobtypmx
    ! ITNOBS holds a count of obs used this step of each type
    itnobs(j) = 0

  END DO
  navt =0

  !       4.1 Loop over groups of ob types
  DO jgroup = 1, n_groups
    jjgroup=jgroup
    igroup     = group_index(jgroup)
    l_ac_timestep = iter_no  <=  def_no_iterations(igroup) .AND.  &
                MOD(timestep_no-1,def_interval_iter(igroup)) == 0

    IF (dg_only .OR. (.NOT. dg_only .AND. l_ac_timestep) ) THEN

      ! First obs type in group
      iactf = group_first(jgroup)

      ! Last obs type in group
      iactl = group_last (jgroup)

      !           4.1.1 Define analysis variable indiciator for this
      !                   group (NAVT)
                 ! NAVT=4 for humidity (cloud), 5 for precip rate
                 ! defined by (1st OB type in group)/100
      navt  = lact(iactf) / 100
      lwind = navt  ==  3

      IF (ldiagac .AND. mype == 0) THEN
        WRITE(umMessage,'(A,I3,A,I3,A,I2,A,(20I4))')                 &
              ' AC STEP =', timestep_no, ',', kiter,              &
              ' STARTING GROUP NO ', jgroup,                      &
              ' , OBS TYPES ', (lact(j), j = iactf, iactl)
        CALL umPrint(umMessage,src='ac-ac1a')
      END IF

      ! Convert Model MIXING RATIO to RELATIVE HUMIDITY
      ! DO for all variables (ensures RH preserved when theta changed)
      ! (except for humidity with 2A cloud microphysics in use and
      !  doing a group of MOPS cloud data
      !  -assumes obs type 406 will always be in group on its own)


      IF (navt /= 4 ) THEN
        hmrmode = 1
        CALL hmrtorh (hmrmode, exner, pressure, theta, rh,        &
               p_field, icode, cmessage)
      END IF

       ! Check for error - drop out if bad.
      IF (icode >  0) GO TO 999

      ! Get Horizontal analysis mode : FI or HORINF?
      mode_hanal = def_mode_hanal(igroup)

      ! Get No of Weight/Analysis Levels for this group.
      no_anal_levs = def_no_anal_levs(igroup)
      no_wt_levs   = def_no_wt_levs  (igroup)

      ! Set up work area for derived theta incrs from lhn

      wklen = 1

      IF (navt == 5 .AND. l_lhn) wklen = model_levels

      no_anal_var = 1

      ! Set size of model grid to that used in lower
      ! level routines
      lenmg = p_field

      !           4.1.2  Decide on analysis grid for this group of types.
      ! parameters for FI method
      lenag       = 1   ! ANAL_INCR not used in FI method
      inc_type = no_anal_var +1


      !           4.1.3 Get a list of obs relevant to this time (getobs)
      !                 for each type in the current group of types
                 ! ASSM_TIME is time since start in minutes (for GETOBS)
      assm_time = REAL(timestep_no -1) * timestep / 60.0

      nptobt = 0

      ! Loop over observation types in group
      DO jact = iactf, iactl
        jjact=jact
         ! LENOBT is no of obs. for type to be used on
         ! this timestep and is determined in GETOBS
        lenobt = 0

        ! INOBS is total no of observations for type
        inobs = obs_info % nobs(jact)
        inobsdim=MAX(inobs,1)

         ! Check on available workspace in OBS_NO
        IF (nptobt+inobs  <=  ac_num_obs_max) THEN
          CALL getobs (jjact, assm_time, inobs, obs_no(nptobt+1),&
                  lenobt, obs(obs_info % mdispobt(jact)+1),      &
                  obs_flag(obs_info % obs_no_st(jact)+1),        &
                  timestep_no, dg_only,inobsdim,icode, cmessage)

          IF (ldiagac .AND. mype == 0) THEN
            WRITE(umMessage,*)'AC af GETOBS - LENOBT ',lenobt
            CALL umPrint(umMessage,src='ac-ac1a')
          END IF
                     ! Check for error - drop out if bad
          IF (icode >  0) GO TO 999

        ELSE
          icode = 2
          WRITE(umMessage,'(A,I0,I0,I0)')                    &
              'TEST FAILED - NPTOBT,INOBS,AC_NUM_OBS_MAX ',  &
              nptobt,inobs,ac_num_obs_max
          CALL umPrint(umMessage,src='ac-ac1a')

                     ! Drop out of assimilation
          GO TO 990

        END IF

        lenact(jact) = lenobt
        nptobt       = nptobt + lenobt

      END DO   ! End of loop over obs types in group (JACT)

      ! LENOB is now the no of observations in this group
      ! to be used on this timestep
      lenob = nptobt

      !
      !   HOW MANY OBS (lenob_total) for all pe's?
      !
      lenob_total = lenob
      CALL gc_isum(1, nproc, istat, lenob_total)

      IF (ldiagac .AND. mype == 0) THEN
        WRITE(umMessage,*)'b4 AC2 - lenob_total,TNDV,LENOB ',    lenob_total,   &
                                                         tndv, lenob
        CALL umPrint(umMessage,src='ac-ac1a')
      END IF
      IF (lenob_total >  0) THEN!  Skip if no obs for this group
        ! increment ob counter array
        itnobs(jgroup) = itnobs(jgroup) + lenob_total


        ! Lower level routine of ac to begin here now that
        ! all the array dimensions are known - lenob, lenag,
        ! no_anal_levs
        IF (lwind) THEN
          WRITE(umMessage,*) 'FATAL ERROR LWIND=T NO LONGER SUPPORTED'
          CALL umPrint(umMessage,src='ac-ac1a')
          cmessage='FATAL ERROR LWIND=T NO LONGER SUPPORTED'
          icode=101
          IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, &
                                  zhook_handle)
          RETURN
        END IF
        CALL ac2(bl_levels, row_length, p_rows,                   &
            p_field, timestep_no, iter_no,                        &
            timestep, obs, tndv, exner, pstar,                    &
            theta, rh, qcl, qcf,                                  &
            conv_cld, ls_rain, ls_snow, conv_rain, conv_snow,     &
               d_theta_dt_conv,d_theta_dt_ls,                     &
               layer_cloud,pressure,                              &
            rhcrit,                                               &
            obs_no, lenob, no_anal_levs, no_wt_levs, no_anal_var, &
            lenag, lenmg, wklen, inc_type, navt, jjgroup, lwind,  &
            iactf, iactl, dg_only, stindex, stlist, len_stlist,   &
            si, sf, stashwork, stash_levels, num_stash_levels,    &
            stash_pseudo_levels, num_stash_pseudo,                &
            lambda_p,phi_p,                                       &
            icode, cmessage)

        ! Check for error - drop out if bad.
        IF (icode >  0) GO TO 999

      END IF

      ! Convert RELATIVE HUMIDITY back to MIXING RATIO
      ! (except for MOPS with 2A cld microphys)

      IF (navt /= 4 ) THEN
        hmrmode = 2
        CALL hmrtorh (hmrmode, exner, pressure, theta, rh,        &
               p_field, icode, cmessage)
      END IF

       ! Check for error - drop out if bad
      IF (icode  >   0) GO TO 999


      IF (ldiagac .AND. mype == 0) THEN
        WRITE(umMessage,'(A,I3,A,I3,A,I4)') ' AC STEP',                  &
            timestep_no, ',', kiter, ' END OF GROUP NO', jgroup
        CALL umPrint(umMessage,src='ac-ac1a')

      END IF
    END IF
  END DO   ! End of loop over groups (JGROUP)
END DO   ! End of loop over total no of iterations (JITER)

IF (mype == 0) THEN
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='ac-ac1a')
  WRITE(umMessage,'(A,I3)')  ' End of AC for time step:',timestep_no
  CALL umPrint(umMessage,src='ac-ac1a')
  WRITE(umMessage,'(A,(10I8))') ' Group No   ',(j,j=1,n_groups)
  CALL umPrint(umMessage,src='ac-ac1a')
  WRITE(umMessage,'(A,(10I8))') ' No of obs  ',(itnobs(j),j=1,n_groups)
  CALL umPrint(umMessage,src='ac-ac1a')
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='ac-ac1a')

END IF

ldiagac = old_ldiagac

! Fix to allow vectorization of JACT loop by removing character
! operations.
990  CONTINUE
IF (icode  ==  2) THEN
  icode    = 1
  cmessage = 'AC : Insufficient space in array OBS_NO.'//         &
                                             'Increase AC_NUM_OBS_MAX'

END IF

999  CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ac

!    Arguments & declarations:
END MODULE ac_mod

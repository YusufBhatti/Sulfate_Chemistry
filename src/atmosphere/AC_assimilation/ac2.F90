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
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!    Arguments and declarations:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE ac2_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='AC2_MOD'

CONTAINS

SUBROUTINE ac2(bl_levels, row_length, p_rows, p_field,            &
               timestep_no, iter_no, timestep,                    &
               obs, tndv,                                         &
               exner, pstar, theta, rh, qcl, qcf,                 &
               conv_cld, ls_rain, ls_snow, conv_rain, conv_snow,  &
               d_theta_dt_conv,d_theta_dt_ls,                     &
               layer_cloud,pressure,                              &
               rhcrit,                                            &
               obs_no, lenob, no_anal_levs, no_wt_levs,           &
               no_anal_var,                                       &
               lenag, lenmg, wklen, inc_type,                     &
               navt, jgroup, lwind, iactf, iactl, dg_only,        &
               stindex, stlist, len_stlist, si, sf,               &
               stashwork, stash_levels, num_stash_levels,         &
               stash_pseudo_levels, num_stash_pseudo,             &
               lambda_p,phi_p,                                    &
               icode, cmessage)

USE umPrintMgr
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore,  ONLY: mype
USE Field_Types, ONLY: fld_type_p
USE comobs_mod, ONLY: nobtypmx, obs_info
USE ac_stash_mod, ONLY: ac_stash
USE addinc_mod, ONLY: addinc
USE diago_mod, ONLY: diago
USE diagopr_mod, ONLY: diagopr
USE fi_mod, ONLY: fi
USE getob2_mod, ONLY: getob2
USE hintcf_mod, ONLY: hintcf
USE lhn_inc_mod, ONLY: lhn_inc
USE lhn_inc_1A_mod, ONLY: lhn_inc_1A
USE mmspt_mod, ONLY: mmspt
USE relaxc_mod, ONLY: relaxc
USE set_relax_cf_mod, ONLY: set_relax_cf
USE vertanl_mod, ONLY: vertanl
USE nlsizes_namelist_mod, ONLY: model_levels

USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Import arguments (INTENT=IN):
INTEGER ::     bl_levels               ! Bdy layer levels
INTEGER ::     row_length              ! Number of points on row
INTEGER ::     p_rows                  ! Number of rows (pstar)
INTEGER ::     p_field                 ! Number of points in
                                       ! mass field
INTEGER ::     tndv
INTEGER ::     timestep_no             ! Current model timestep
INTEGER ::     iter_no                 ! Iteration number
                                       ! excluding diagnostic
INTEGER ::     lenmg                   ! Length of model grid
INTEGER ::     lenag                   ! Length of analysis grid
INTEGER ::     wklen                   ! Length of 2nd dimension
                                       ! of array for derived incs
                                       ! produced by HTDRST,
                                       ! GEOSTR & WINDBAL.
INTEGER ::     inc_type                ! Index for data types in
                                       ! MODEL_INCR
INTEGER ::     lenob                   ! No of obs in group
INTEGER :: i,iproc
INTEGER ::     no_anal_levs            ! No of analysis levels
                                       ! for group of obs types
INTEGER ::     no_anal_var             ! No of variables being
                                       ! analysed (2 for winds)
INTEGER ::     no_wt_levs              ! No of weight levels
                                       ! for group of obs types
INTEGER ::     navt                    ! Analysis variable type
                                       ! 1/2/3/4/5 =
                                       ! p*/theta/winds/rh/precip
INTEGER ::     jgroup                  ! Loop counter over groups
INTEGER ::     iactf                   ! First obs type in group
INTEGER ::     iactl                   ! Last  obs type in group
INTEGER ::     obs_no(lenob+1)

LOGICAL ::     lwind                   ! Indicator if group
                                       ! has wind obs types
LOGICAL ::     dg_only                 ! Indicator if diagnostic
                                       ! only iteration
REAL ::        timestep                ! Timestep in seconds
REAL ::        obs(tndv+1)             ! Observation data (set up
                                       ! in RDOBS)
REAL ::        conv_cld(p_field, model_levels) ! conv cld amount
REAL ::        ls_rain(p_field)        ! large scale rain rate
REAL ::        ls_snow(p_field)        ! large scale snow rate
REAL ::        conv_rain(p_field)      ! convective rain rate
REAL ::        conv_snow(p_field)      ! convective snow rate
                                       ! above rates diagnostic
REAL ::        d_theta_dt_conv(p_field,model_levels)
                                ! convective latent heating rate
REAL ::        d_theta_dt_ls(p_field,model_levels)
                                ! large scale latent heating rate

REAL ::        rhcrit(model_levels)    ! Critical RH array
                                       ! for cloud
! These are used in variable resolution runs
REAL :: lambda_p(1-halo_i:row_length+halo_i)
REAL :: phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)

! Import / export arguments (INTENT=INOUT):
REAL ::        exner(p_field, model_levels)       ! exner on theta levels
REAL ::        layer_cloud(p_field, model_levels) ! as name says
REAL ::        pressure (p_field, model_levels)   ! p on theta levels
REAL ::        pstar(p_field)                     ! Prognostic variable P*
REAL ::        theta(p_field, model_levels)       ! Prognostic variable theta
REAL ::        rh(p_field, model_levels)          ! Prognostic variable HMR
! for 2A cloud microphysics, but otherwise contains rh at this stage
REAL ::        qcl(p_field, model_levels)         ! Prognostic variable QCL
REAL ::        qcf(p_field, model_levels)         ! Prognostic variable QCF

! Export arguments (INTENT=OUT):
INTEGER ::     icode                   ! Error code (0 = OK)

CHARACTER(LEN=errormessagelength) :: cmessage                ! Error message

! Stash variables
REAL ::        stashwork(*)

INTEGER ::     len_stlist
INTEGER ::     num_stash_levels
INTEGER ::     num_stash_pseudo
INTEGER ::     stindex(2, *)
INTEGER ::     stlist(len_stlist, *)
INTEGER ::     si(*)
INTEGER ::     stash_levels(num_stash_levels +1, *)
INTEGER ::     stash_pseudo_levels(num_stash_pseudo +1, *)

LOGICAL ::     sf(*)

!
!  Single iteration of analysis correction data assimilation
!  for observations or group of observations of variable type
!  defined by NAVT

! Local (dynamic) arrays:
INTEGER ::     model_obs_type(lenob+1)
                                   ! Model observation type No
INTEGER ::     ne_ag_pt(lenob+1)     ! Nearest analysis grid
                                      ! point to obs.
INTEGER ::     nppt(lenob+1, 4)       ! horizontal model ->
                                      ! obs interp pointers

REAL ::        obs_lat (lenob+1)      ! Observation latitude
REAL ::        obs_long(lenob+1)      ! "longitude
REAL ::        obs_time(lenob+1)      ! "time
REAL ::        obs_incr(lenob+1, no_anal_levs, no_anal_var)!
                                      ! Obs - Model increments
                                      ! at obs points.
REAL ::        normf(lenob+1, no_anal_levs)
                                         ! Normalization factor
REAL ::        anal_incr_local(lenmg*no_anal_levs*no_anal_var) !
                                       ! Accumulated increments
                                       ! on analysis grid
REAL ::        model_incr(lenmg, no_anal_levs, inc_type) !
                                       ! Accumulated increments
                                       ! on model grid
REAL ::        drincs(p_field, wklen)  ! Array to hold derived
                                       ! increment fields
REAL ::        cfpt(lenob+1, 4)
                                       ! obs interp coeffs
REAL ::        relax_cf(p_rows)        ! Relaxation coeffs for
                                       ! current group

LOGICAL ::     lmissd(lenob+1, no_anal_levs)
                                       ! missing data in obs.
REAL ::        pstgs(p_field)          ! PSTAR at U points
INTEGER :: dum3
REAL :: dum2(1,1), dum1(1)
INTEGER :: bc,bcount

! Local scalars:
INTEGER ::     j                       ! Miscelaneous loop cntr
INTEGER ::     jact,jjact              ! Loop counter over obs
                                       ! types in group
INTEGER ::     janl                    ! Loop counter over
                                       ! analysis levels
INTEGER ::     jlev,jjlev              ! Loop counter over model
                                       ! levels
INTEGER ::     jl,jjl
INTEGER ::     jk,jjk
INTEGER ::     jpoint                  ! Loop counter over points

INTEGER ::     ipass                   ! 1 if weights pass
                                       ! 2 if increment pass
INTEGER ::     iag                     ! Offset for level in
                                       ! ANAL_INCR array.
INTEGER ::     ijk                     ! Level number in HYDRST
                                       ! & GEOSTR loops
INTEGER ::     inobs                   ! No of obs for this type
INTEGER ::     nptobt                  ! Offset to first obs of
                                       ! each type in group
INTEGER ::     lenobt                  ! No of obs for obs type
INTEGER ::     ndv                     ! No of data values
INTEGER ::     n_rows                  ! No of rows
INTEGER ::     no_anal_incr_levs       ! No of analysis increment
                                       ! levels.
INTEGER ::     mode_hanal              ! FI or HORINF for this
                                       ! group
INTEGER ::     stash_item_no           ! Stash Item Number
LOGICAL ::     surface_wind            ! True if surface wind data

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='AC2'

!- End of header section.


!  Subroutine structure:
!  1.0 Initialize:
! Mode of Horizontal Analysis for this group
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
mode_hanal = def_mode_hanal(group_index(jgroup))

!  2.0 Weights calculation starts here

IF (ldiagac .AND. mype == 0) THEN
  CALL umPrint(' ',src='ac2')
  WRITE(umMessage,'(A,(10I4))')                                          &
      ' START OF WTS ANALYSIS PHASE FOR TYPES ', &
      (lact(j), j = iactf, iactl)
  CALL umPrint(umMessage,src='ac2')
END IF

!  2.1 Setup for interpolation model -> obs & vertical analysis
!      Get common format data for all obs in group       (getob2)
! Vectorisation over all observation types in group
! Loop over obs types in group
IF (lenob /= 0) THEN
  nptobt = 0

  DO jact = iactf, iactl
    jjact=jact
    inobs  = obs_info % nobs(jact)
    lenobt = lenact(jact)

    IF (inobs >  0 .AND. lenobt >  0) THEN
      CALL getob2 (jjact, obs(obs_info % mdispobt(jact) +1), inobs, &
                  obs_lat(nptobt +1), obs_long(nptobt +1),          &
                  obs_time(nptobt +1), model_obs_type(nptobt +1),   &
                  obs_no(nptobt +1), lenobt, icode, cmessage)

      ! Check for error - exit routine if occured
      IF (icode  >   0) GO TO 999

    END IF

    ! Move pointer to next obs type in group
    nptobt = nptobt + lenobt

  END DO

  !  2.2  Horizontal interpolation coefficients model -> obs  (hintcf)
  CALL hintcf(lwind, lenob, obs_lat, obs_long, obs_no, row_length,  &
             p_rows, cfpt(1,1), cfpt(1,2), cfpt(1,3), cfpt(1,4),    &
             nppt(1,1), nppt(1,2), nppt(1,3), nppt(1,4), icode,     &
             lambda_p,phi_p,cmessage)

  ! Check for error - exit routine if occured
  IF (icode  >   0) GO TO 999

  !  2.3 Do vertical analysis
  !      Make increment & error factor vectors for each type in this
  !      group. The details of processing method depend on ob type
  !      so vectorization is over LENOBT obs in one type
  !      rather than LENOB obs in the group. NPTOBT gives the
  !      increment to point to each section in turn of those
  !      vectors which go over all obs in group.

  ! Loop over obs types in group
END IF ! lenob /= 0
nptobt = 0

DO jact = iactf, iactl
  jjact=jact
  inobs  = obs_info % nobs(jact)
  lenobt = lenact(jact)
  !
  !
  !
  IF (inobs == 0 .OR. lenobt == 0) THEN
    IF (ldiagac) THEN
      IF (lact(jact) == 506) THEN
        CALL diagopr (1,dum1,dum2,lenobt,lenob,ndv,                     &

                      lmissd,dum3,nptobt,no_anal_levs)
      ELSE
        bcount=2
        IF (lact(jact) == 101 .OR.                                          &
           (lact(jact) >= 202 .AND. lact(jact) <= 204) .OR.                  &
           (lact(jact) >= 302 .AND. lact(jact) <= 306) .OR.                  &
           (lact(jact) >= 402 .AND. lact(jact) <= 404) .OR.                  &
           lact(jact) == 406 .OR.                                          &
           lact(jact) == 901)bcount=1
        !
        ! BARRIERS FOR DIAGO
        !
        DO bc=1,bcount
          DO jlev=1,glsize(3,fld_type_p)
            DO i=0,8
              s_stat(jlev,i)=0.0
            END DO
          END DO
          IF (lact(jact) == 406) THEN
            ! dummy diago call
            CALL diago ('MULTI-LEVEL', lact(jact),6,obs_incr,normf,          &
                         obs_lat, obs_long, lmissd, lenob, lenobt, nptobt,   &
                         no_anal_levs, no_anal_var)
          ELSE
            ! dummy diago call
            CALL diago ('VAN?', lact(jact), 3+bc-bcount,                     &
                         obs_incr, normf,                                    &
                         obs_lat, obs_long, lmissd, lenob, lenobt, nptobt,   &
                         no_anal_levs, no_anal_var)
          END IF ! LACT(JACT) == 406
        END DO
      END IF ! LACT(JACT) == 506
    END IF ! LDIAGAC
  ELSE !inobs == 0.or.lenobt == 0

    IF (inobs  >   0 .AND. lenobt  >   0) THEN
      ndv   = obs_info % ndatav(jact) - obs_info % ndvhdr
      ipass = 1

      CALL vertanl (pstgs,                                          &
            jjact, ipass, lenobt, ndv, obs_no(nptobt +1),           &
            model_obs_type(nptobt +1), obs_lat(nptobt +1),          &
            obs_long(nptobt +1), obs(obs_info % mdispobt(jact) +1), &
            inobs, exner, pstar, theta, rh, qcl, qcf, conv_cld,     &
            layer_cloud,pressure,                                   &
            ls_rain, ls_snow, conv_rain, conv_snow,                 &
            rhcrit,                                                 &
            p_field, model_incr, lenmg, no_wt_levs,                 &
            cfpt(nptobt+1, 1), cfpt(nptobt+1, 2),                   &
            cfpt(nptobt+1, 3), cfpt(nptobt+1, 4),                   &
            nppt(nptobt+1, 1), nppt(nptobt+1, 2),                   &
            nppt(nptobt+1, 3), nppt(nptobt+1, 4),                   &
            obs_incr, normf, lmissd,                                &
            bl_levels, row_length, p_rows,                          &
            lenob, nptobt, no_anal_levs, no_anal_var,               &
            icode, cmessage)

      ! Check for error - exit routine if occured
      IF (icode  >   0) GO TO 999

    END IF

    ! Move pointer to next obs type in group
    nptobt = nptobt + lenobt

  END IF !inobs == 0.or.lenobt == 0
END DO

! Vectorisation elsewhere over all obs in group
nptobt = 0
lenobt = lenob

IF (ldiagac) THEN
  IF (iactf  /=  iactl) THEN
    ! The group has >1 type in it, so group stats are worthwhile.
    ! Print mean & rms stats for all obs in group

    IF (lldac(2)) THEN
      CALL diago ('AC', lact(iactf), 1, obs_incr, normf,           &
                 obs_lat, obs_long, lmissd, lenob, lenobt, nptobt, &
                 no_anal_levs, no_anal_var)

    END IF
  END IF
END IF


!  2.4 Do horizontal analysis (for weights)
ipass = 1

IF (mode_hanal  ==  2) THEN   !  Use FI method
  IF (lenob /= 0) THEN

    CALL fi (ipass, timestep_no, timestep, navt, jgroup,            &
            no_anal_var, lenob, no_anal_levs, no_wt_levs, 1, lenmg, &
            row_length, p_rows, obs_no, obs_lat, obs_long, obs_time,&
            normf, obs_incr, cfpt, nppt, inc_type, model_incr,      &
            pstar, p_field, icode, cmessage)
  ELSE
    DO jlev=1,no_anal_levs
      DO j=1,lenmg
        model_incr(j,jlev,1)=0.0
      END DO
    END DO
  END IF


   ! Check for error - exit from routine if occured
  IF (icode  >   0) GO TO 999

END IF ! End of FI horizontal analysis step

IF (ldiagac) THEN
  CALL umPrint('END OF WTS ANALYSIS PHASE',src='ac2',pe=0)
END IF

! Save all levels of wts on model grid

DO jlev = 1, no_wt_levs   !  Loop over no of weight levels.
  jjlev=jlev

   ! Save weights field

  stash_item_no = 200+navt
  IF (sf(stash_item_no)) THEN
    CALL ac_stash (stash_item_no,jjlev, jgroup,                   &
      n_groups, timestep_no, stindex, stlist, len_stlist, si, sf, &
      stashwork, stash_levels, num_stash_levels,                  &
      stash_pseudo_levels, num_stash_pseudo, model_incr(1,jlev,1),&
      lenmg, 'Weights Field', icode, cmessage)

  END IF

  ! Get statistics on sum of weights
  IF (ldiagac .AND. lldac(4)) THEN
    IF (navt  ==  4) THEN
      CALL mmspt (model_incr(1,jlev,1), jjlev, 0,                 &
                          'SUM OF WTS - RH ', row_length, p_rows, &
                           phi_p)

    ELSE IF (navt  ==  5) THEN
      CALL mmspt (model_incr(1,jlev,1), jjlev, 0,                 &
                          'SUM OF WTS - PR ', row_length, p_rows, &
                           phi_p)
    ELSE      ! NAVT values other than 4 or 5 excluded
      icode = 1
      cmessage = 'AC : Should not reach here (3)'
      GO TO 999   ! RETURN

    END IF
  END IF
END DO   ! End of loop over no of weight levels.

!  2.6 Calculate ob normalisation factors
! Vectorization here is over types within groups
nptobt = 0

DO jact = iactf, iactl
  jjact=jact
  inobs  = obs_info % nobs(jact)
  lenobt = lenact(jact)
  !
  !
  !
  IF (inobs == 0 .OR. lenobt == 0) THEN
    IF (ldiagac) THEN
      IF (lact(jact) /= 406 .AND. lact(jact) /= 506) THEN
        bcount=2
        IF (lact(jact) == 101 .OR.                                          &
           (lact(jact) >= 202 .AND. lact(jact) <= 204) .OR.                  &
           (lact(jact) >= 302 .AND. lact(jact) <= 306) .OR.                  &
           (lact(jact) >= 402 .AND. lact(jact) <= 404) .OR.                  &
           lact(jact) == 901)bcount=1
        !
        ! BARRIERS FOR DIAGO
        !
        DO bc=1,bcount
          DO jlev=1,glsize(3,fld_type_p)
            DO i=0,8
              s_stat(jlev,i)=0.0
            END DO
          END DO
          ! dummy diago call
          CALL diago ('VAN?', lact(jact), 5+bc-bcount,                     &
                       obs_incr, normf,                                    &
                       obs_lat, obs_long, lmissd, lenob, lenobt, nptobt,   &
                       no_anal_levs, no_anal_var)
        END DO
      END IF ! LACT(JACT) /= 406 (and 506)
    END IF ! LDIAGAC
  ELSE ! lenob == 0

    IF (inobs >  0 .AND. lenobt >  0) THEN
      ndv   = obs_info % ndatav(jact) - obs_info % ndvhdr
      ipass = 2

      CALL vertanl (pstgs,                                          &
             jjact, ipass, lenobt, ndv, obs_no(nptobt +1),          &
             model_obs_type(nptobt +1), obs_lat(nptobt +1),         &
             obs_long(nptobt+1), obs(obs_info % mdispobt(jact) +1), &
             inobs, exner, pstar, theta, rh, qcl, qcf, conv_cld,    &
             layer_cloud,pressure,                                  &
             ls_rain, ls_snow, conv_rain, conv_snow,                &
             rhcrit,                                                &
             p_field, model_incr, lenmg, no_anal_levs,              &
             cfpt(nptobt+1, 1), cfpt(nptobt+1, 2),                  &
             cfpt(nptobt+1, 3), cfpt(nptobt+1, 4),                  &
             nppt(nptobt+1, 1), nppt(nptobt+1, 2),                  &
             nppt(nptobt+1, 3), nppt(nptobt+1, 4),                  &
             obs_incr, normf, lmissd,                               &
             bl_levels, row_length, p_rows,                         &
             lenob, nptobt, no_anal_levs, no_anal_var,              &
             icode, cmessage)

      ! Check for errors:
      IF (icode  >   0) GO TO 999

    END IF

    nptobt = nptobt + lenobt
  END IF ! lenob == 0

END DO

! Vectorisation elsewhere over all obs in group
nptobt = 0
lenobt = lenob

!  3.0 Increments calculation starts here

IF (ldiagac .AND. mype == 0) THEN
  WRITE(umMessage,'(A,(10I4))')                                          &
      ' START OF INCS ANALYSIS PHASE FOR TYPES  ',                       &
      (lact(j),j=iactf,iactl)
  CALL umPrint(umMessage,src='ac2')

END IF

!  3.1 Horizontal analysis step (increments)
ipass = 2

IF (mode_hanal  ==  2) THEN   ! Use FI method

  IF (lenob /= 0) THEN
    CALL fi (ipass, timestep_no, timestep, navt, jgroup,            &
           no_anal_var, lenob, no_anal_levs, no_wt_levs, 1, lenmg,  &
           row_length, p_rows, obs_no, obs_lat, obs_long, obs_time, &
           normf, obs_incr, cfpt, nppt, inc_type, model_incr,       &
           pstar,                                                   &
           p_field, icode,                                          &
           cmessage)
  ELSE
    DO jlev=1,no_anal_levs
      DO j=1,lenmg
        model_incr(j,jlev,1)=0.0
      END DO
    END DO
  END IF

   ! Check for errors:
  IF (icode  >   0) GO TO 999

END IF ! End of FI horizontal analysis step

IF (ldiagac) THEN
  CALL umPrint('END OF INCS ANALYSIS PHASE',src='ac2',pe=0)
END IF

!  3.3 Set up Relaxation Coefficients for this group
n_rows = p_rows

CALL set_relax_cf (jgroup, n_rows, relax_cf, lwind, timestep,     &
     timestep_no, iter_no, icode, cmessage)

! Check for errors:
IF (icode  >   0) GO TO 999

!  3.4 Interpolation from analysis grid to model grid
!      Loop over levels to interpolate incrs from analysis grid to
!      model grid. Add incrs to model.
DO jlev = 1, no_anal_levs
  jjlev=jlev

  !      3.4.1 Scale increment by relaxation coefficent

  CALL relaxc (model_incr(1,jlev,1), lenmg, row_length,         &
    p_rows, relax_cf, icode, cmessage)

  ! Check for errors:
  IF (icode  >   0) GO TO 999


  !      3.4.2 Save increment fields into output file

  stash_item_no = 210+navt
  IF (sf(stash_item_no)) THEN
    CALL ac_stash (stash_item_no, jjlev, jgroup,                &
      n_groups, timestep_no, stindex, stlist, len_stlist, si,   &
      sf, stashwork, stash_levels, num_stash_levels,            &
      stash_pseudo_levels, num_stash_pseudo,                    &
      model_incr(1,jlev,1), lenmg, 'Increment Field', icode,    &
      cmessage)

  END IF


  !      3.4.10 Humidity increments
  IF (navt  ==  4) THEN

    ! Add RH increments to RH fields
    CALL addinc (rh(1,jlev), model_incr(1,jlev,1), lenmg,         &
      row_length, p_rows, navt, icode, cmessage)

    ! Check for errors:
    IF (icode  >   0) GO TO 999

    ! Statistics of increments on model grid
    IF (ldiagac .AND. lldac(4)) THEN
      CALL mmspt (model_incr(1,jlev,1), jjlev, 0,                 &
        'RH INCREMENTS   ', row_length, p_rows,                   &
        phi_p)

    END IF
  END IF   ! NAVT = 4

  !      3.4.11 latent heat nudging
  IF (navt  ==  5) THEN
    !  Stash LS latent heating theta increments if required (rates in K/s)
    stash_item_no = 271

    IF (sf(stash_item_no)) THEN
      DO jl = 1 , model_levels
        jjl = jl
        CALL ac_stash (stash_item_no, jjl,                      &
          jgroup, n_groups, timestep_no,                        &
          stindex, stlist, len_stlist, si, sf, stashwork,       &
          stash_levels, num_stash_levels,                       &
          stash_pseudo_levels, num_stash_pseudo,                &
          d_theta_dt_ls(1,jl),p_field,                          &
          'LS latent heating rates',                            &
          icode, cmessage)
      END DO    ! JL
    END IF

    IF (l_lhn) THEN
      IF(l_lhn_1a) THEN
        ! Call the 1A versions of the LHN routines which use a different
        ! decomposition strategy.
        CALL lhn_inc_1a(d_theta_dt_conv,d_theta_dt_ls,                    &
                     ls_rain,ls_snow,conv_rain,                           &
                     conv_snow,p_field,model_incr,timestep,               &
                     p_rows,row_length,                                   &
                     lhn_range, phi_p,lambda_p,                           &
                     drincs,lenmg,icode,cmessage)
      ELSE        
        CALL lhn_inc(d_theta_dt_conv,d_theta_dt_ls,                       &
                     ls_rain,ls_snow,conv_rain,                           &
                     conv_snow,p_field,model_incr,timestep,               &
                     p_rows,row_length,                                   &
                     lhn_range, phi_p,lambda_p,                           &
                     drincs,lenmg,icode,cmessage)
      END IF
      !C Add on increments
      DO j=1,model_levels
        CALL addinc(theta(1,j), drincs(1,j), p_field, row_length, &
                    p_rows, navt,icode,cmessage)
      END DO     !J

      !  Stash LHN Theta Increments if required (amounts in K)
      stash_item_no = 272

      IF (sf(stash_item_no)) THEN
        DO jl = 1 , model_levels
          jjl = jl
          CALL ac_stash (stash_item_no, jjl,                      &
            jgroup, n_groups, timestep_no,                        &
            stindex, stlist, len_stlist, si, sf, stashwork,       &
            stash_levels, num_stash_levels,                       &
            stash_pseudo_levels, num_stash_pseudo,                &
            drincs(1,jl),p_field,'LHN theta increments',          &
            icode, cmessage)
        END DO    ! JL
      END IF

      IF (ldiagac .AND. lldac(4) ) THEN
        DO jl = 1 , model_levels
          jjl = jl
          CALL mmspt(drincs(1,jl), jjl, 0,                          &
                      'LHN Theta incrs ', row_length, p_rows,       &
                       phi_p)
        END DO     ! JL
      END IF
    END IF       ! L_LHN

  END IF    ! NAVT = 5

END DO   ! JLEV


!  5.0 Exit from AC2
999  CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ac2
END MODULE ac2_mod

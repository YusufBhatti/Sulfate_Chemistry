! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE FI -------------------------------------------------
!
!    Purpose :
!
!       Perform the horizontal spreading of observational increments.
!       There are three stages:
!       1. calculate the horizontal scale to be used
!       2. spread the observation weights to get normalization factors
!       3. spread the normalized increments.
!       Spreading is performed using the adjoint of the model->ob interp
!       followed by a filter on the model grid using RF.
!
!       This is an alternative to the method using HORINF.
!
!   NB: This routine is now only used for MOPS rain and since
!   NPASS_RF=0 FILT=F so recursive filtering routines are no longer
!   called. For this reason the code has not been upgraded to
!   handle variable resolution grids in this mode.
!
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Logical components covered:
!
!    Project Task : P3
!
!    Documentation: UM Doc Paper 30 section 6.4
!
!

!    Arguments:---------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE fi_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FI_MOD'

CONTAINS

SUBROUTINE fi(ipass,timestep_no,timestep,navt,jgroup,             &
  no_anal_var,lenob,no_anal_levs,no_wt_levs,first_level,          &
  lenmg,row_length,p_rows,obs_no,                                 &
  obs_lat,obs_long,obs_time,normf,obs_incr,cfpt,nppt,             &
  ip_scale,model_incr,pstar,                                      &
  p_field,                                                        &
  icode,cmessage)

USE planet_constants_mod, ONLY: planet_radius

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx
USE rfcsl_mod, ONLY: rfcsl
USE rfvsg_mod, ONLY: rfvsg
USE rfvsl_mod, ONLY: rfvsl
USE rfvvg_mod, ONLY: rfvvg
USE rfvvl_mod, ONLY: rfvvl
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

INTEGER :: ipass,timestep_no
INTEGER :: jgroup,first_level
INTEGER :: navt,no_anal_var,no_anal_levs,no_wt_levs,ip_scale
INTEGER :: lenob,lenmg,row_length,p_rows
REAL ::   timestep
INTEGER :: icode  ! 0 if OK
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER :: obs_no(lenob),nppt(lenob+1,4)
REAL ::                                                           &
     obs_lat(lenob),obs_time(lenob),cfpt(lenob+1,4),              &
     obs_long(lenob),normf(lenob+1,no_anal_levs),                 &
     obs_incr(lenob+1,no_anal_levs,no_anal_var),                  &
     model_incr(lenmg,no_anal_levs,ip_scale),                     &
     pstar(row_length*p_rows)

INTEGER ::                                                        &
 p_field   !IN PTS ON U GRID;P GRID


!-INTENT=IN--------------------------------------------------------
!     IPASS : mode
!            1- do stages 1. calculate scale
!                       & 2. spread NORMF
!            0- do stage  3. spread INCR
!     TIMESTEP_NO: current model timestep (TIMESTEP_NO=1,.....N)
!     NAVT: observation variable type
!     JGROUP : ob group code
!     NO_ANAL_LEVS: number of levels analysed ( >=  NO_WT_LEVS)
!     FIRST_LEVEL : model level corresponding to analysis level 1
! N.B. current AC code assumes this is 1; this could be relaxed later
!     NO_WT_LEVS: number of levels for which NORMF is needed
!     LENOB: number of observations in group
!     LENMG : NUMBER OF POINTS IN FIELD ON model GRID
!     ROW_LENGTH : length of model grid row
!     P_ROWS : no of row on model p grid (assumed > U_ROWS)
!     TIMESTEP   : is model timestep in secs
!     OBS_NO : list of ob numbers for this group within ACOBS file
!     OBS_LAT  :  co-lat of obs
!     OBS_LONG : longitude of obs
!     OBS_TIME : time of obs rel to assm start
!     NORMF    : normalization factors for obs (from routine VERTANL)
!     OBS_INCR : increment (ob-model) at ob pt
!     CFPT : model->ob interpolation coeffs
!     NPPT : model->ob interpolation pointers
!     PSTAR : surface pressure (used in DIVFILT for mass weighting)
!     IP_SCALE : 3rd index to MODEL_INCR, pointing to space for scale
!     IP_SCALE=2 when doing a scalar, =3 when doing a vector analysis.
!-INTENT=OUT---------------------------------------------------------
!    stage 1:
!     MODEL_INCR(,,IP_SCALE)     : horizontal scale for stages 2 & 3
!     ( MODEL_INCR(,,1:2) )      : is used as work space in stage 1
!    stage 2:
!     MODEL_INCR(,,1)            : data density field D for norm.factor
!                                : (only first NO_WT_LEVS levels done)
!    stage 3: (second call to FI)
!     MODEL_INCR(,,1:NO_ANAL_VAR): increments
!-INTENT=IN   (second call to FI)
!     MODEL_INCR(,,IP_SCALE)     : horizontal scale
! ---------------------------------------------------------------------

!    Workspace usage:-------------------------------------------------
REAL :: cos_lat(p_rows) ! sin of co-latitude of model grid rows
REAL :: s_incr(lenob,no_anal_levs,2) ! incs at obs to be scattered
REAL :: tf(lenob) ! Time Factor for each ob  as in HORINF
REAL :: s(lenob)  ! Scale for each ob        as in HORINF
REAL :: wrk1(lenob)
INTEGER :: iwk1(lenob)

!----------------------------------------------------------------------
!    Define local variables
INTEGER :: jrow,job,jobt,j,jlev,jc ! loop indices
INTEGER :: i,itype
INTEGER :: irows ! number of rows on model grid for current variable
LOGICAL :: iluv  ! .true. if winds
LOGICAL :: filt  ! .true. if filtering on model grid is reqd.
INTEGER :: m_grid ! 0=Ltd Area, 1=pts at pole, 2=staggd from pole
REAL :: timeb,timea,radinf
REAL :: cscale_start,cscale_obtime,cscale_end
INTEGER :: first_obs,last_obs,first_type,last_type
REAL :: z,zz,zintegral,zioa
REAL :: fi_var_factor  ! to change SCALE to match tuning to HORINF
REAL, ALLOCATABLE :: delta_lat(:)
REAL, ALLOCATABLE :: delta_long(:)

CHARACTER(LEN=10) :: label

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FI'

!----------------------------------------------------------------------

!  **     PRELIMINARIES
!         -------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ltimer_ac) CALL timer('FI      ',3)

icode=0
filt = npass_rf >  0
IF (filt) THEN
  zioa=(25.1327-18.2655/npass_rf)/(dlatmg*dlongmg)
  !     ZIOA is a constant in the formula for correlation area/box area
  !     ZIOA=v (6.5.15) / (dlat*dlong)
  !     For NPASS_RF =2, v=16.; NPASS_RF=infinity, v=8*PI.
  !     v for other NPASS_RF is got from an approximate interpolation.
  !     the COS_LAT(JROW) term in the grid box area is done in the loops
  !
  !     HORINF truncated the influence area at RADINF corr scales
  !     and uses an EXP(-r/SCALE) correlation function for wind -
  !     FI_VAR_FACTOR compensates for this so FI gives similar
  !     results to HORINF with the default values for scale etc..
  fi_var_factor=def_fi_var_factor(group_index(jgroup))
  fi_var_factor=SQRT(fi_var_factor)
ELSE
  zioa=1.0
  fi_var_factor=1.0
END IF  !  FILT

IF (filt) THEN ! these only required if FILT=T so using recursive filters
  iluv=navt == 3
  IF (iluv) THEN ! model U grid
    z=row1mguv
    irows=p_rows-1
    m_grid=2
  ELSE          ! model P grid
    z=row1mgth
    irows=p_rows
    m_grid=1
  END IF
  !     INITIALIZE COS_LAT
  DO      j=1,irows
    cos_lat(j)=SIN(z)
    z=z+dlatmg
  END DO ! J
END IF
m_grid=0              ! no pole in ltd area version

! *** define time-factors and scales for each observation
!     The following code for section 2.4 was copied from HORINF (2.6)
!     (preliminary version 1/6/92) and updated 16/2/93
!     so that FI behaves as similarly as possible to HORINF for testing
!
!L2.4 TF which is time factor R(dT) in 3.18 and 3.21
!L##  AND S(DT) BEING THE CORRELATION SCALE AS IN EQUATION 3.19

!     for synoptic insertion mode treat all obs times as though at
!     the nearest synoptic hour I.E. T(OB) = T+0
IF (lsyn) THEN
  DO 2400 job=1,lenob
    obs_time(job) = obtime_nom
    !      value of OBTIME_NOM is relative to start of assimilation
  2400   END DO ! JOB
END IF

first_obs = 1

!     Get first and last types for group being processed.
first_type = group_first(jgroup)
last_type  = group_last (jgroup)

DO jobt = first_type,last_type   !  Loop over types in group

  IF (lenact(jobt) >  0) THEN   !  Any obs for this type

    last_obs = first_obs + lenact(jobt) - 1

    !         Get analysis parameters for this type
    itype = type_index(jobt)

    timeb         = def_timeb(itype)
    timea         = def_timea(itype)
    radinf        = def_radinf(itype)
    cscale_start  = def_cscale_start(itype)
    cscale_obtime = def_cscale_obtime(itype)
    cscale_end    = def_cscale_end(itype)

    DO job = first_obs,last_obs
      !       WRK1= time relative obs time (+ve before and -ve after)
      wrk1(job) = obs_time(job) - (timestep_no-1)*timestep/60.0
      tf(job)=0.0
      s(job)=1.0
      !       note that these initialised values will not change for an ob
      !       with relative time outside the insertion period (as required)
    END DO ! JOB

    DO job = first_obs,last_obs
      !     calculate ramp functions between start of insertion & ob time
      !     data on edge of insertion period not fetched by GETOBS, so may
      !     be discounted here too.
      IF (wrk1(job) >  0.0 .AND. wrk1(job) <= timeb) THEN
        IF (mrampfn == 2) THEN
          ! quadratic function
          tf(job) =                                                        &
          wrk1(job)*wrk1(job)*(timef_start-timef_obtime)/(timeb*timeb)     &
           + timef_obtime
        ELSE IF (mrampfn == 1) THEN
          ! linear function
          tf(job) =                                                        &
          wrk1(job)*(timef_start-timef_obtime)/timeb + timef_obtime
        END IF

        s(job) =                                                         &
        wrk1(job) * ((cscale_start-cscale_obtime)/timeb) + cscale_obtime
      END IF
    END DO ! JOB

    DO job = first_obs,last_obs
      !     calculate ramp functions between ob time and end of insertion
      IF (wrk1(job) <= 0.0 .AND. wrk1(job) >= -timea) THEN
        IF (mrampfn == 2) THEN
          ! quadratic function
          tf(job) =                                                        &
          wrk1(job)*wrk1(job)*(timef_end-timef_obtime)/(timea*timea)       &
           + timef_obtime
        ELSE IF (mrampfn == 1) THEN
          ! linear function
          tf(job) =                                                        &
          wrk1(job) * (timef_obtime-timef_end)/timea + timef_obtime
        END IF

        s(job) =                                                         &
        wrk1(job) * ((cscale_obtime-cscale_end)/timea) + cscale_obtime
      END IF
    END DO ! JOB

    !     WRK1 can be reused.

    !     IWK1,WRK1 CAN BE REUSED.

    IF (MOD(mhcorfn/16,2) == 1) THEN
      DO job = first_obs,last_obs
        !         reduce time factor where scale is large so that the ramp
        !         function specified holds for the area integrated weight
        wrk1(job) = cscale_obtime/s(job)
        wrk1(job) = wrk1(job)*wrk1(job)
        tf(job)   = tf(job)*wrk1(job)
      END DO ! JOB
    END IF

  END IF   !  Jump here if no obs for this type.

  first_obs = first_obs + lenact(jobt)
END DO    !  End of JOBT loop over obs types in group.

IF (ldiagac) THEN
  IF (lldac(1) .AND. ndac >  0 .AND. ipass == 1) THEN
  END IF
END IF


IF (ipass == 1) THEN! do stages 1 & 2
  !  **                  <<<<    STAGE 1: SCALE    >>>>
  !  ** 1.  calculate scale as smoothed average of that of local obs
  !         using UM Doc 30 (6.4.17)
  !     1.1 Initialize MODEL_INCR
  DO      jc=1,2
    DO      jlev=1,no_anal_levs
      DO      j=1,lenmg
        model_incr(j,jlev,jc)=0.0
      END DO ! J
    END DO ! JLEV
  END DO ! JC

  !     1.2 put values for each ob,level into S_INCR
  DO      jlev=1,no_anal_levs
    DO      job=1,lenob
      !       see remark in UM Doc P30 5.1 about the multiple uses of NORMF
      !       NORMF here is Q1m of (5.1.2) i.e. 1/E**2
      !       TF is used as in (6.1.4)
      s_incr(job,jlev,1)=normf(job,jlev)*tf(job)
      s_incr(job,jlev,2)=normf(job,jlev)*tf(job)*                     &
        s(job)*fi_scale_factor(first_level-1+jlev) ! in metres
      !       FI allows a different scale for each level
      !       To get results comparable with HORINF use FI_SCALE_FACTOR=1.
    END DO ! JOB
  END DO ! JLEV

  IF (filt) THEN
    !     1.3 add S_INCR into MODEL_INCR using transpose of model->obs intrp

    DO     job=1,lenob
      !       the 4 passes through the loop over JC can be concurrent as
      !       long as NPPT(J,1) NPPT(J,2) NPPT(J,3) & NPPT(J,4) are different.
      !       The version of HINTCF introduced with FI ensures that they are.
      DO      jc=1,4
        !       DO      J=1,2               ! the loops on J & JLEV should
        !       DO      JLEV=1,NO_ANAL_LEVS !collapse into a single vector loop
        j=1                     ! collapse loops by hand

        DO      jlev=1,no_anal_levs*2   ! collapse loops by hand
          model_incr(nppt(job,jc),jlev,j)=                               &
          model_incr(nppt(job,jc),jlev,j)+s_incr(job,jlev,j)*cfpt(job,jc)
        END DO ! JLEV
        !       ENDDO ! J                       ! collapse loops by hand
      END DO ! JC
    END DO

    !     1.4 divide by grid-box areas, to convert from transpose to adjoint
    !         (The grid box area is Px in (6.4.15) ).
    !         multiply by area integral of correlation function
    zintegral=zioa*(fi_scale/planet_radius)**2 ! convert to radians
    i=0
    DO jrow=1,irows
      zz=zintegral/cos_lat(jrow)
      DO      jc=1,2
        DO      jlev=1,no_anal_levs
          z=zz*fi_scale_factor(first_level-1+jlev)**2
          DO      j=1,row_length
            model_incr(i+j,jlev,jc)=model_incr(i+j,jlev,jc)*z
          END DO ! J
        END DO ! JLEV
      END DO ! JC
      i=i+row_length
    END DO ! JROW

    !     1.5 filter weighted scales, and weights
    IF (ltimer_ac) CALL timer('RFCSloop',3)
    ! Allocate arrays to store row/col deltas.  This is needed since rfcsl
    ! routine was updated for variable resolution while this routine is not going
    ! to be (see comment in header)

    ALLOCATE(delta_lat(irows))
    ALLOCATE(delta_long(row_length))
    delta_lat(:) = dlatmg
    delta_long(:)= dlongmg

    !     collapse JC & JLEV loops so that multitasking can be over both
    !     DO      JC=1,2
    !     DO      JLEV=1,NO_ANAL_LEVS
    DO      j=   0,no_anal_levs*2-1 ! J loop replaces JC & JLEV loops
      jc=MOD(j,2)+1                   ! J loop replaces JC & JLEV loops
      jlev=j/2+1                      ! J loop replaces JC & JLEV loops
      z=fi_scale*fi_scale_factor(first_level-1+jlev)/planet_radius

      !       These calls can be done concurrently
      CALL rfcsl(model_incr(1,jlev,jc),irows,row_length,m_grid,       &
       z,cos_lat,delta_lat,delta_long,z,npass_rf)
      !     ENDDO ! JLEV
      !     ENDDO ! JC
    END DO ! J                       ! J loop replaces JC & JLEV loops
    DEALLOCATE(delta_lat)
    DEALLOCATE(delta_long)

    IF (ltimer_ac) CALL timer('RFCSloop',4)

    !     1.6 divide to give SCALE, using (6.4.15). metre=>radians using A.
    DO      jlev=1,no_anal_levs
      DO      j=1,lenmg
        model_incr(j,jlev,ip_scale)= (model_incr(j,jlev,2) +            &
        0.1*fi_scale*fi_scale_factor(first_level-1+jlev))/              &
           ((model_incr(j,jlev,1)+0.1)*planet_radius*fi_var_factor)
      END DO ! J
    END DO ! JLEV


  END IF  !  FILT

  !  **                  <<<<    STAGE 2: NORMF    >>>>
  !  ** 2.  calculate data-density and hence Q to be used in VERTANL
  !         using UM Doc 30 (6.4.15)
  !     2.1 Initialize MODEL_INCR
  DO      jlev=1,no_wt_levs
    DO      j=1,lenmg
      model_incr(j,jlev,1)=0.0
    END DO ! J
  END DO ! JLEV

  !     2.2 put values for each ob,level into S_INCR
  !     This need not be done; the values are already there from step 1.2
  !     DO JLEV=1,NO_WT_LEVS
  !     DO JOB=1,LENOB
  !       S_INCR(JOB,JLEV,1)=NORMF(JOB,JLEV)*TF(JOB)
  !     ENDDO ! JOB
  !     ENDDO ! JLEV

  !     2.3 add S_INCR into MODEL_INCR using transpose of model->obs intrp

  DO  job=1,lenob
    !       the 4 passes through the loop over JC could be concurrent as
    !       long as NPPT(J,1) NPPT(J,2) NPPT(J,3) & NPPT(J,4) are different.
    !       The HINTCF introducted with FI makes sure that they differ.
    DO      jc=1,4
      DO      jlev=1,no_wt_levs
        model_incr(nppt(job,jc),jlev,1)=                               &
        model_incr(nppt(job,jc),jlev,1)+s_incr(job,jlev,1)*cfpt(job,jc)
      END DO ! JLEV
    END DO ! JC
  END DO

  IF (filt) THEN
    !     2.4 divide by grid-box areas, to convert from transpose to adjoint
    !         (The grid box area is Px in (6.4.15) ).
    !         multiply by area integral of correlation function
    !         Do not add weight given to background; this is done in VERTANL
    i=0
    DO jrow=1,irows
      z=zioa/cos_lat(jrow)
      DO      jlev=1,no_wt_levs
        DO      j=1,row_length
          model_incr(i+j,jlev,1)=model_incr(i+j,jlev,1)*z*              &
                                 model_incr(i+j,jlev,ip_scale)**2
        END DO ! J
      END DO ! JLEV
      i=i+row_length
    END DO ! JROW

    !     2.5 filter weights to get effective data density
    IF (ltimer_ac) CALL timer('RFVSloop',3)
    DO      jlev=1,no_wt_levs
      !       These calls can be done concurrently
      CALL rfvsl(model_incr(1,jlev,1),irows,row_length,m_grid,0.0,     &
       cos_lat,dlatmg,dlongmg,model_incr(1,jlev,ip_scale),npass_rf)
    END DO ! JLEV
    IF (ltimer_ac) CALL timer('RFVSloop',4)

  END IF  !  FILT

ELSE ! (IPASS == 2)THEN ! do stage  3

  !  **                  <<<<    STAGE 3: INCR     >>>>
  !  ** 3.  calculate INCR on model grid
  !         using UM Doc 30 (6.4.16)
  !     3.1 Initialize MODEL_INCR
  DO      jc=1,no_anal_var
    DO      jlev=1,no_anal_levs
      DO      j=1,lenmg
        model_incr(j,jlev,jc)=0.0
      END DO ! J
    END DO ! JLEV
  END DO ! JC

  !     3.2 put values for each ob,level into S_INCR

  DO   jc=1,no_anal_var
    DO      jlev=1,no_anal_levs
      DO      job=1,lenob
        !       NORMF (from VERTANL) is now Q2m of (6.1.3) i.e. (1/E**2) * Q
        !       TF is used as in (6.1.1)
        s_incr(job,jlev,jc)=obs_incr(job,jlev,jc)*tf(job)**2*           &
                               normf(job,jlev)
      END DO ! JLEV
    END DO ! JOB
  END DO ! JOB


  !     3.3 add S_INCR into MODEL_INCR using transpose of model->obs intrp
  DO job=1,lenob
    !       the 4 passes through the loop over JC could be concurrent as
    !       long as NPPT(J,1) NPPT(J,2) NPPT(J,3) & NPPT(J,4) are different.
    DO      jc=1,4
      !       DO      J=1,NO_ANAL_VAR     ! the loops on J & JLEV should
      !       DO      JLEV=1,NO_ANAL_LEVS !collapse into a single vector loop
      !       the JLEV*J loop is too short to multitask
      j=1                              ! collapse explicitly
      DO      jlev=1,no_anal_levs*no_anal_var  ! collapse explicitly
        model_incr(nppt(job,jc),jlev,j)=                               &
        model_incr(nppt(job,jc),jlev,j)+s_incr(job,jlev,j)*cfpt(job,jc)
      END DO ! JLEV  & collapsed J
      !       ENDDO ! J                                ! collapse explicitly
    END DO ! JC
  END DO


  IF (filt) THEN
    !     3.4 divide by grid-box areas, to convert from transpose to adjoint
    !         (The grid box area is Px in (6.4.15) ).
    !         multiply by area integral of correlation function
    !             which is about 16*SCALE**2 for a 2-pass filter.
    i=0
    DO jrow=1,irows
      z=zioa/cos_lat(jrow)
      DO      jc=1,no_anal_var ! JC is never equal to IP_SCALE
        DO      jlev=1,no_anal_levs
          DO      j=1,row_length
            model_incr(i+j,jlev,jc)=model_incr(i+j,jlev,jc)*z*            &
                                    model_incr(i+j,jlev,ip_scale)**2
          END DO ! J
        END DO ! JLEV
      END DO ! JC
      i=i+row_length
    END DO ! JROW
    IF (ldiagac) THEN
      IF (lldac(1) .AND. ndac >  0) THEN
        DO jlev=1,no_anal_levs
          DO job=1,lenob ! interpolate incr at ob position
            wrk1(job)=model_incr(nppt(job,1 ),jlev,1)*cfpt(job, 1)        &
                     +model_incr(nppt(job,2 ),jlev,1)*cfpt(job, 2)        &
                     +model_incr(nppt(job,3 ),jlev,1)*cfpt(job, 3)        &
                     +model_incr(nppt(job,4 ),jlev,1)*cfpt(job, 4)
          END DO ! JOB
          IF (iluv) THEN
            DO job=1,lenob ! interpolate V incr at ob position
              wrk1(job)=model_incr(nppt(job,1 ),jlev,2)*cfpt(job, 1)        &
                       +model_incr(nppt(job,2 ),jlev,2)*cfpt(job, 2)        &
                       +model_incr(nppt(job,3 ),jlev,2)*cfpt(job, 3)        &
                       +model_incr(nppt(job,4 ),jlev,2)*cfpt(job, 4)
            END DO ! JOB
          END IF ! ILUV
        END DO ! JLEV
      END IF
    END IF

    !     3.5 filter increments
    IF (iluv) THEN ! filter winds as a vector field
      IF (ltimer_ac) CALL timer('RFVVloop',3)
      DO jlev=1,no_anal_levs
        !       These calls can be done concurrently
        IF (model_type == mt_global) THEN
          CALL rfvvg(model_incr(1,jlev,1),model_incr(1,jlev,2),           &
           irows,row_length,m_grid,cos_lat,row1mguv,                      &
           dlatmg,dlongmg,model_incr(1,jlev,ip_scale),npass_rf)
        ELSE
          CALL rfvvl(model_incr(1,jlev,1),model_incr(1,jlev,2),           &
           irows,row_length,m_grid,cos_lat,row1mguv,                      &
           dlatmg,dlongmg,model_incr(1,jlev,ip_scale),npass_rf)
        END IF
      END DO ! JLEV
      IF (ltimer_ac) CALL timer('RFVVloop',4)
    ELSE          ! filter scalar increment field
      IF (ltimer_ac) CALL timer('RFVSloop',3)
      DO jlev=1,no_anal_levs
        !       These calls can be done concurrently
        IF (model_type == mt_global) THEN
          CALL rfvsg(model_incr(1,jlev,1),irows,row_length,m_grid,        &
           cos_lat,dlatmg,dlongmg,model_incr(1,jlev,ip_scale),npass_rf)
        ELSE
          CALL rfvsl(model_incr(1,jlev,1),irows,row_length,m_grid,0.0,     &
           cos_lat,dlatmg,dlongmg,model_incr(1,jlev,ip_scale),npass_rf)
        END IF
      END DO ! JLEV
      IF (ltimer_ac) CALL timer('RFVSloop',4)
    END IF ! ILUV

  END IF ! FILT

  IF (ldiagac) THEN
    IF (lldac(1) .AND. ndac >  0) THEN
      DO jlev=1,no_anal_levs
        DO job=1,lenob ! interpolate incr at ob position
          wrk1(job)=model_incr(nppt(job,1 ),jlev,1)*cfpt(job, 1)        &
                   +model_incr(nppt(job,2 ),jlev,1)*cfpt(job, 2)        &
                   +model_incr(nppt(job,3 ),jlev,1)*cfpt(job, 3)        &
                   +model_incr(nppt(job,4 ),jlev,1)*cfpt(job, 4)
        END DO ! JOB
        IF (iluv) THEN
          DO job=1,lenob ! interpolate V incr at ob position
            wrk1(job)=model_incr(nppt(job,1 ),jlev,2)*cfpt(job, 1)        &
                     +model_incr(nppt(job,2 ),jlev,2)*cfpt(job, 2)        &
                     +model_incr(nppt(job,3 ),jlev,2)*cfpt(job, 3)        &
                     +model_incr(nppt(job,4 ),jlev,2)*cfpt(job, 4)
          END DO ! JOB
        END IF ! ILUV
      END DO ! JLEV
    END IF
  END IF


END IF ! IPASS

IF (ltimer_ac) CALL timer('FI      ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fi
END MODULE fi_mod

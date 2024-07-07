! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINES VANMOPS_MIXED_PHASE and CC_TO_RHTOT-------------------
!
!  Purpose : Performs vertical analysis of MOPS cloud data
!            for (2A) mixed phase cloud microphysics scheme
!     When IPASS=1,
!     model field is interpolated to ob locations and increments
!     calculated. Preliminary weight normalisation factors are also
!     calculated.
!
!     When IPASS=2,
!     data density is interpolated to ob locations and
!     final weight normalisation factors are calculated.
!
!  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!  Project Task : P3
!
!  Arguments:---------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE vanmops_mixed_phase_mod
USE umPrintMgr
 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VANMOPS_MIXED_PHASE_MOD'

CONTAINS

SUBROUTINE vanmops_mixed_phase(                                  &
                   kact,ipass,theta,exner,conv_cld,              &
                   cf,pressure,                                  &
                   q,qcl,qcf,lenfld,nlvfld,rhcrit,               &
                   obdata,cf1pt,cf2pt,cf3pt,cf4pt,               &
                   np1pt,np2pt,np3pt,np4pt,                      &
                   qinc,normf,obs_lat,obs_long,                  &
                   obs_no,lmissd,nptobt,lenobt,ndv,lenob,        &
                   no_anal_levs,no_anal_var,                     &
                   bl_levels,icode,cmessage)

USE planet_constants_mod, ONLY: cp
USE water_constants_mod, ONLY: lc
USE conversions_mod, ONLY: zerodegc
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE comobs_mod, ONLY: nobtypmx, obs_info
USE cc_to_rhtot_mod, ONLY: cc_to_rhtot
USE diago_mod, ONLY: diago
USE hintmo_mod, ONLY: hintmo
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength
USE nlsizes_namelist_mod, ONLY: model_levels
!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    qsat_wat_new     => qsat_wat,                       &
                    l_new_qsat_acassim !Currently defaults to FALSE

IMPLICIT NONE

INTEGER :: lenfld,nlvfld,lenobt,lenob,ndv,no_anal_levs,           &
           bl_levels,no_anal_var
REAL ::                                                           &
               theta(lenfld,model_levels),                        &
               exner(lenfld,model_levels),                        &
               pressure (lenfld,model_levels),                    &
               conv_cld(lenfld,nlvfld),                           &
               cf(lenfld,nlvfld),                                 &
               q   (lenfld,nlvfld),                               &
               qcl (lenfld,nlvfld),                               &
               qcf (lenfld,nlvfld),                               &
               rhcrit  (nlvfld),                                  &
               obdata  (lenobt,ndv),                              &
               qinc(lenob+1,no_anal_levs),                        &
               normf(lenob+1,no_anal_levs),                       &
               obs_lat(lenobt),obs_long(lenobt),                  &
               cf1pt(lenobt), cf2pt(lenobt),                      &
               cf3pt(lenobt), cf4pt(lenobt)
LOGICAL :: lmissd(lenob+1,no_anal_levs)
INTEGER :: obs_no(lenobt)
INTEGER :: np1pt(lenobt),np2pt(lenobt),np3pt(lenobt)
INTEGER :: np4pt(lenobt)

INTEGER :: kact,ipass,nptobt
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage

!-INTENT=IN--------------------------------------------------------
!     KACT      - Observation type
!     IPASS     - Calculate incrs and preliminary norm factors for
!                 IPASS=1 and final norm factors for IPASS=2.
!     LENOBT    - no of obs in type
!     LENOB     - no of obs in group
!     NPTOBT    - pointer to present ob type within group
!     NDV       - no of data values in observation
!     NO_ANAL_LEVS - no of analysis levels
!     NO_ANAL_VAR  - no of analysis variables
!     P_LEVELS     - no of model levels
!     EXNER     - diagnostic variable exner pressure
!     CONV_CLD  - convective cloud amount on model levels
!     CF         - layer cloud amount (liq+ice) on model levels
!     PRESSURE   - p on theta levels
!     Q         - model specific humidity      (IPASS=1)
!               - data density on model grid   (IPASS=2)
!     QCF       - model cloud ice
!     QCL       - model cld liq water
!     THETA     - model potential temperature
!     LENFLD    - length of model field
!     NLVFLD    - no of levels in Q field
!     RHCRIT    - critical rh values for cloud formation (fraction)
!     OBDATA    - observed values
!     OBS_LAT   - ob co-lats
!     OBS_LONG  - ob longitudes
!     OBS_NO    - observation numbers
!     NP1PT,NP2PT,NP3PT,NP4PT,
!               - pointers to nearest model grid points for each ob
!     CF1PT,CF2PT,CF3PT,CF4PT,
!               - weighting of model points in bilinear interpolation
!
!-INTENT=INOUT-----------------------------------------------------
!     QINC        -  ob-model humidity increments
!     NORMF       -  normalisation factors
!-INTENT=OUT-----------------------------------------------------
!     LMISSD          -  logical to indicate missing data
!     ICODE,CMESSAGE  - error code and message
!
!----------------------------------------------------------------------
!     Workspace usage
!-----------------------------------------------------------------------
!     DYNAMIC ALLOCATION ON CRAY
REAL ::        tl(lenfld,nlvfld)

REAL ::        cld_inc (lenob+1,no_anal_levs)
REAL ::        cld_norm(lenob+1,no_anal_levs)

REAL ::        cf_lyrob(lenobt), cf_convb(lenobt)
REAL ::        qb(lenobt), rhtot_ob(lenobt)
REAL ::        cfb(lenobt), tlb(lenobt), qclb(lenobt)
REAL ::        qcfb(lenobt), pb(lenobt)
REAL ::        qsat_iceb(lenobt), qsat_watb(lenobt)

!     TL         - model liquid water temperature
!     CLD_INC    - obs - model cloud fraction (layer+convective)
!     CLD_NORM   - normalisation array needed for input to DIAGO call
!     CF_LYROB   - target layer cloud fraction derived from cloud ob
!     CF_CONVB   - model convective cloud fraction interp'd to obs pts
!     QB         - model humidity (or obs density) at obs point
!     RHTOT_OB   - target rhtot derived from CF_LYROB
!     CFB        - CF at obs points
!     TLB        - TL at obs points
!     QCLB,QCFB  - QCL,QCF at obs points
!     PB         - model level pressure at obs pt
!     QSAT_ICEB  - QSAT wrt ice at obs pt
!     QSAT_WATB  - QSAT wrt water at obs pt
!
!
!----------------------------------------------------------------------
!     Define local variables
!----------------------------------------------------------------------
REAL :: tqsatice,tqsatwat,dtqsat
INTEGER :: ktype,nobipt,noblpt,neropt,jk,job,j,error

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='VANMOPS_MIXED_PHASE'

!     KTYPE    -   Observation type
!     NOBIPT   -   Pointer to data values within OBDATA array
!     NOBLPT   -   Pointer to obs levels within OBDATA array
!     NEROPT   -   Pointer to obs errors within OBDATA array
!     JK       -   Loop counter in loops over levels
!     JOB      -         "      in loops over obs
!     J        -         "      in loops over points
!     ERROR    -   error code
!     TQSATICE -   new cloud likely to be ice below this temp
!                  (see THOMO in C_LSPMIC)
!     TQSATWAT -   new cloud likely to be liquid above this temp
!                  (see TNUC in C_LSPMIC)
!     DTQSAT   -   difference between two previous temperatures
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ltimer_ac) CALL timer('VANMMP  ',3)

!
!  ** 1. PRELIMINARIES
!        -------------

ktype = lact (kact)

IF (ktype == 406) THEN
  nobipt = 1
END IF
neropt = obs_info % nerlev1(kact) - obs_info % ndvhdr

!     Define temps used in choosing qsat value when creating cloud
!     (These temps are the same as the homogeneous and heterogeneous
!      nucleation thresholds in the mixed phase precip scheme,
!      but need not be, hence separate definition)
tqsatice = zerodegc - 40.0
tqsatwat = zerodegc - 10.0
dtqsat   = tqsatwat - tqsatice

IF (ktype == 406) THEN
  !     =================
  IF (ipass == 1) THEN
    !  ** 2.1 Current model cloud fraction is consistent with
    !         latest humidity, since ls_cld called just before AC.
    !         But need to get current TL for use in qsat calcns.


    !     2.1.2 calculate TL
    DO jk=1,nlvfld
      DO j=1,lenfld
        tl(j,jk) = theta(j,jk)*exner(j,jk)- lc * qcl(j,jk) /cp
      END DO
    END DO

    !  ** 2.2  HORIZONTAL INTERP OF FIELDS ON MODEL GRID


  END IF ! IPASS == 1

  ! begin loop over levels
  DO jk=1,no_anal_levs

    !       PUT Q OR OBS DENSITY INTERPOLATED TO OB POSITIONS IN QB
    CALL hintmo( q(1,jk),cf1pt,cf2pt,cf3pt,cf4pt,                   &
                 np1pt,np2pt,np3pt,np4pt,                           &
                 lenfld,1,lenobt,qb,icode,cmessage)
    IF (icode >  0) GO TO 9999

    !  ** 2.3  CALCULATION OF INCRS AND WEIGHT FACTORS FOR OBS INC

    !     2.3.1  FIRST SET UP BIT ARRAY FOR NO OBSERVATIONAL DATA
    DO job=1,lenobt
      lmissd(nptobt+job,jk) = obdata(job,nobipt+jk-1) == obs_info % missd &
                         .OR. obdata(job,neropt+jk-1) == obs_info % missd
    END DO   !JOB

    IF (ipass == 1) THEN

      !      2.3.2 put conv_cld interp'd to ob locations in CF_CONVB
      CALL hintmo( conv_cld(1,jk),cf1pt,cf2pt,cf3pt,cf4pt,             &
                   np1pt,np2pt,np3pt,np4pt,                            &
                   lenfld,1,lenobt,cf_convb,icode,cmessage)
      IF (icode >  0) GO TO 9999

      !      2.3.3 calculate target layer cloud cover in CF_LYROB
      !            based on total cloud = conv_cld + ls_cld *(1-conv_cld)
      !            as assumed in radiation scheme.
      DO job=1,lenobt
        IF (cf_convb(job) /=  1.0) THEN
          cf_lyrob(job)=(obdata(job,nobipt+jk-1)-cf_convb(job))/      &
                                               (1.0-cf_convb(job))
        ELSE
          !         set target layer cloud of zero where conv cover = 1
          cf_lyrob(job)=0.0
        END IF
        !         set target layer cloud of zero where conv cover > MOPS
        !         (will apply also where MOPS data is missing data value<0,
        !          though obs cloud not used then)
        cf_lyrob(job)=MAX(0.0,cf_lyrob(job) )
      END DO

      !      2.3.4 convert target layer cloud to rhtot
      CALL cc_to_rhtot (cf_lyrob,lenobt,rhcrit(jk),rhtot_ob,jk,       &
                                                 bl_levels)

      !      2.3.5 interp model b'ground fields to ob pts for incr calcns
      !            calculate derived quantities

      !        2.3.5.1 put TL    interp'd to ob pts in TLB
      CALL hintmo( tl(1,jk),cf1pt,cf2pt,cf3pt,cf4pt,                &
                   np1pt,np2pt,np3pt,np4pt,                         &
                   lenfld,1,lenobt,tlb,icode,cmessage)
      IF (icode >  0) GO TO 9999

      !        2.3.5.2 put QCL   interp'd to ob pts in QCLB
      CALL hintmo( qcl(1,jk),cf1pt,cf2pt,cf3pt,cf4pt,               &
                   np1pt,np2pt,np3pt,np4pt,                         &
                   lenfld,1,lenobt,qclb,icode,cmessage)
      IF (icode >  0) GO TO 9999

      !        2.3.5.3 put QCF   interp'd to ob pts in QCFB
      CALL hintmo( qcf(1,jk),cf1pt,cf2pt,cf3pt,cf4pt,               &
                   np1pt,np2pt,np3pt,np4pt,                         &
                   lenfld,1,lenobt,qcfb,icode,cmessage)
      IF (icode >  0) GO TO 9999

      !        2.3.5.4 put TOTAL LAYER CLOUD FRACTION interp'd to
      !                                                   ob pts in CFB
      CALL hintmo( cf(1,jk),cf1pt,cf2pt,cf3pt,cf4pt,                &
                   np1pt,np2pt,np3pt,np4pt,                         &
                   lenfld,1,lenobt,cfb,icode,cmessage)
      IF (icode >  0) GO TO 9999

      !        2.3.5.5 put MODEL LEVEL PRESSURE interp'd to ob pts in PB
      CALL hintmo( pressure(1,jk),cf1pt,cf2pt,cf3pt,cf4pt,              &
                       np1pt,np2pt,np3pt,np4pt,                         &
                       lenfld,1,lenobt,pb,icode,cmessage)
      IF (icode >  0) GO TO 9999

      !        2.3.5.6 get QSAT(TL,P) wrt ice in QSAT_ICEB
      !        2.3.5.7 get QSAT(TL,P) wrt water in QSAT_WATB
      IF ( l_new_qsat_acassim ) THEN
        CALL qsat_new (qsat_iceb, tlb, pb, lenobt)
        CALL qsat_wat_new (qsat_watb, tlb, pb, lenobt)
      ELSE
        ! DEPENDS ON: qsat
        CALL qsat (qsat_iceb, tlb, pb, lenobt)
        ! DEPENDS ON: qsat_wat
        CALL qsat_wat (qsat_watb, tlb, pb, lenobt)
      END IF

      !        2.3.5.8 update RHTOT_OB where ice present
      !                covers eqns 12 & 13 in working paper
      DO job=1,lenobt
        IF (qcfb(job)  >   0.0) THEN
          rhtot_ob(job) = rhtot_ob(job) *                          &
           (qclb(job) + qcfb(job)*                                 &
                           qsat_iceb(job)/qsat_watb(job) )         &
          /( qclb(job) + qcfb(job) )
        END IF
      END DO

      !      2.3.6 calculate increments to Q and normalisation factors
      DO job=1,lenobt
        IF (.NOT. lmissd(nptobt+job,jk)) THEN
          !             Data is valid
          !             Calculate incrs to Q
          IF ( cfb(job) == cf_lyrob(job) ) THEN
            !               model and ob cloud equal, so no humidity incr
            qinc(nptobt+job,jk) = 0.0
          ELSE IF (cfb(job) >  0.0) THEN
            !               eqns 7 & 8 in working paper
            qinc(nptobt+job,jk) = rhtot_ob(job)*qsat_watb(job) -    &
                                  ( qb(job)+qclb(job)+qcfb(job) )
          ELSE IF (cfb(job) == 0.0 .AND. cf_lyrob(job) >  0.0) THEN
            !               need to create cloud in model,
            !               refer to section 2.5 of working paper
            IF (tlb(job) >  tqsatwat) THEN
              qinc(nptobt+job,jk) = rhtot_ob(job)                   &
                                   *qsat_watb(job) - qb(job)
            ELSE IF (tlb(job) <  tqsatice) THEN
              qinc(nptobt+job,jk) = rhtot_ob(job)                   &
                                   *qsat_iceb(job) - qb(job)
            ELSE
              qinc(nptobt+job,jk)  = rhtot_ob(job) *                &
                (qsat_iceb(job)*(1-(tlb(job)-tqsatice)/dtqsat) +    &
                 qsat_watb(job)*(tlb(job)-tqsatice)/dtqsat )        &
                - qb(job)
            END IF
          END IF
          !             check humidity incr reasonable, or reset to zero
          !             (humidity and cloud incrs must have same sign)
          IF ( qinc(nptobt+job,jk)*                                  &
                        (cf_lyrob(job)-cfb(job)) <  0.0 ) THEN
            qinc(nptobt+job,jk) = 0.0
          END IF
          !             Calculate normalisation factor
          normf(nptobt+job,jk) = 1.0 /                              &
              ( obdata(job,neropt+jk-1)*obdata(job,neropt+jk-1) )
        ELSE
          !             Missing data
          qinc(nptobt+job,jk)  = 0.0
          normf(nptobt+job,jk) = 0.0
        END IF ! test on missing data
      END DO  !  JOB

      !     2.3.7  store cloud increments (oktas) for diagnostics
      IF (ldiagac) THEN
        DO job=1,lenobt
          IF (.NOT. lmissd(nptobt+job,jk)) THEN
            cld_inc (job,jk) = 8*(obdata(job,nobipt+jk-1)-cf_convb(job) &
                                 -cfb(job)*(1.0 - cf_convb(job) )  )
            cld_norm(job,jk) = 1.0
          ELSE
            cld_inc (job,jk) = 0.0
            cld_norm(job,jk) = 0.0
          END IF
        END DO

      END IF


    ELSE IF (ipass == 2) THEN
      DO job=1,lenobt
        IF (.NOT. lmissd(nptobt+job,jk)) THEN
          normf(nptobt+job,jk) = normf(nptobt+job,jk) /             &
                                ( qb(job) + 1.0 )
        ELSE
          normf(nptobt+job,jk) = 0.0
        END IF
      END DO !  JOB
    END IF ! IPASS


  END DO   ! JK

  !  ** 3.  Diagnostics
  !         ===========
  !     measure fit of model to MOPS cloud cover
  IF (ldiagac .AND. ipass == 1 ) THEN
    CALL umPrint(' ',src='vanmops_mixed_phase',pe=0)
    WRITE(umMessage,'('' DIAGO called from VANMOPS_MIXED_PHASE'//&
        ' - Obs Type'',I5,T50,''No of obs '',I6)') ktype,lenobt
    CALL umPrint(umMessage,src='vanmops_mixed_phase',pe=0)
    CALL diago ('MULTI-LEVEL',ktype,6,                            &
        cld_inc,cld_norm,obs_lat,obs_long,lmissd,         &
        lenobt,lenobt,0,no_anal_levs,no_anal_var)
  END IF

END IF  ! KTYPE
!     =====

IF (ltimer_ac) CALL timer('VANMMP ',4)

9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE vanmops_mixed_phase



END MODULE vanmops_mixed_phase_mod

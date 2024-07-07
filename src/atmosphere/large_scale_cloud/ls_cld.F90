! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Large-scale Cloud Scheme.
! Subroutine Interface:
SUBROUTINE ls_cld(                                                      &
!      Pressure related fields
 p_theta_levels, rhcrit,                                                &
!      Array dimensions
 levels, bl_levels,                                                     &
 rhc_row_length,rhc_rows,                                               &
!      From convection diagnosis (only used if A05_4A)
 ntml, cumulus, l_mixing_ratio,                                         &
!      Prognostic Fields
 t, cf, q, qcf, qcl,                                                    &
!      Liquid and frozen ice cloud fractions
 cfl, cff,                                                              &
 error)

USE cv_run_mod,            ONLY: l_param_conv
USE vectlib_mod,           ONLY: powr_v
USE conversions_mod,       ONLY: pi
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: pdims, tdims
USE cloud_inputs_mod,      ONLY: cloud_fraction_method,               &
 ice_fraction_method, i_eacf, overlap_ice_liquid, ctt_weight,         &
 t_weight, qsat_fixed, sub_cld, smith_orig, cloud_top_temp,           &
 min_liq_overlap, all_clouds, not_mixph

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_wat_new     => qsat_wat,                       &
                    qsat_wat_mix_new => qsat_wat_mix,                   &
                    l_new_qsat_lsc !Currently defaults to FALSE

USE c_cldsgs_mod, ONLY: qcfmin
USE missing_data_mod, ONLY: rmdi
IMPLICIT NONE
!
! Purpose:
!   This subroutine calculates liquid and ice cloud fractional cover
!   for use with the enhanced precipitation microphysics scheme.

! Method:
!   Statistical cloud scheme separates input moisture into specific
!   humidity and cloud liquid water. Temperature calculated from liquid
!   water temperature. Cloud fractions calculated from statistical
!   relation between cloud fraction and cloud liquid/ice water content.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP No. 29


!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
INTEGER ::                                                            &
                      !, INTENT(IN)
 levels,                                                              &
!       No. of levels being processed.
   bl_levels,                                                           &
!       No. of boundary layer levels
   rhc_row_length,rhc_rows

REAL ::                                                               &
                      !, INTENT(IN)
 qcf(           tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,levels),                    &
!       Cloud ice content at processed levels (kg water per kg air).
   p_theta_levels(pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end,levels),                    &
!       pressure at all points (Pa).
   rhcrit(rhc_row_length,rhc_rows,levels)
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.

INTEGER ::                                                            &
 ntml(          tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end)
!       IN Height of diagnosed BL top

LOGICAL ::                                                            &
 l_mixing_ratio
!       IN true if using mixing ratios

LOGICAL ::                                                            &
 cumulus(       tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end)
!       IN Logical indicator of convection

REAL ::                                                               &
                      !, INTENT(INOUT)
 q(             tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,levels),                    &
!       On input : Total water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                  (kg water per kg air).
   t(             tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,levels)
!       On input : Liquid/frozen water temperature (TL) (K).
!       On output: Temperature at processed levels (K).

REAL ::                                                               &
                      !, INTENT(OUT)
 cf(            tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,levels),                    &
!       Cloud fraction at processed levels (decimal fraction).
   qcl(           tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,levels),                    &
!       Cloud liquid water content at processed levels (kg per kg air).
   cfl(           tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,levels),                    &
!       Liquid cloud fraction at processed levels (decimal fraction).
   cff(           tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,levels)
!       Frozen cloud fraction at processed levels (decimal fraction).

!     Error Status:
INTEGER :: error     !, INTENT(OUT)  0 if OK; 1 if bad arguments.

!  Local parameters and other physical constants------------------------
REAL :: rootwo       ! Sqrt(2.)
REAL :: subgrid      ! Subgrid parameter in ice cloud calculation

!  Local scalars--------------------------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
REAL ::                                                               &
 phiqcf,                                                              &
                      ! Arc-cosine term in Cloud ice fraction calc.
 cosqcf,                                                              &
                      ! Cosine term in Cloud ice fraction calc.
 overlap_max,                                                         &
                      ! Maximum possible overlap
 overlap_min,                                                         &
                      ! Minimum possible overlap
 overlap_random,                                                      &
                      ! Random overlap
 temp0,                                                               &
 temp1,                                                               &
 temp2,                                                               &
                      ! Temporaries for combinations of the
 qn_imp,                                                              &
 qn_adj
!                       ! overlap parameters

!  (b) Others.
INTEGER :: k,i,j       ! Loop counters: K - vertical level index.
!                                        I,J - horizontal field indices.

INTEGER :: qc_points,                                                 &
                      ! No. points with non-zero cloud
        multrhc   ! Zero if (rhc_row_length*rhc_rows) le 1, else 1

!  Local dynamic arrays-------------------------------------------------
!    6 blocks of real workspace are required.
REAL ::                                                               &
 qcfrbs(        tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end),                           &
!       qCF / bs
   qsl(           tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end),                           &
   qsl_out(       tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end),                           &
!       Saturated specific humidity for temp TL or T.
   qsl_ctt(       tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end),                           &
   qsl_ctt_out(   tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end),                           &
!       Saturated specific humidity wrt liquid at cloud top temperature
   qn(            tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end),                           &
!       Cloud water normalised with BS.
   grid_qc(       tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,levels),                    &
!       Gridbox mean saturation excess at processed levels
!        (kg per kg air). Set to RMDI when cloud is absent.
   bs(            tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,levels),                    &
!       Maximum moisture fluctuation /6*sigma at processed levels
!        (kg per kg air). Set to RMDI when cloud is absent.
   ctt(           tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,levels)
!       Ice cloud top temperature (K) - as coded it is really TL

REAL, ALLOCATABLE :: cfl_max(:,:,:)
!      Maximum value of liquid cloud in a column

LOGICAL ::                                                            &
 lqc(           tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end)
!       True for points with non-zero cloud
INTEGER ::                                                            &
 idx(tdims%j_len*tdims%i_len,2),                                    &
!       Index for points with non-zero cloud
   llwic(         tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end)
!       Last Level With Ice Cloud
REAL :: rhcritx
!       Scalar copy of RHCRIT(I,J,K)

!       Variables for cache-blocking
INTEGER            :: jj             !Block index

!       jblock is not a parameter. Value needs to be flexible at runtime.
INTEGER            :: jblock


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LS_CLD'

!- End of Header

!Set value of blocking factor, jblock. May need to be changed on other
! platforms or for OMP threads > 4
jblock = 4

! ----------------------------------------------------------------------
!  Check input arguments for potential over-writing problems.
! ----------------------------------------------------------------------
error=0

IF ( (rhc_row_length * rhc_rows)  >   1) THEN
  multrhc = 1
ELSE
  multrhc = 0
END IF

! ==Main Block==--------------------------------------------------------
! Subroutine structure :
! Loop round levels to be processed.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Initialize cloud-top-temperature and last-level-with-ice-cloud arrays

! Levels_do1:
!$OMP  PARALLEL DEFAULT(NONE)                                                 &
!$OMP& SHARED(llwic, ctt, levels, qcl, cfl, grid_qc, bs, p_theta_levels,      &
!$OMP& l_mixing_ratio, qcf, ice_fraction_method, ctt_weight, i_eacf,          &
!$OMP& rhc_row_length, rhc_rows, bl_levels, cloud_fraction_method, cf,        &
!$OMP& overlap_ice_liquid, cff, t_weight, sub_cld, q, t, jblock, cfl_max,     &
!$OMP& qsat_fixed, multrhc, cumulus, ntml, rhcrit, l_param_conv, tdims,       &
!$OMP& l_new_qsat_lsc)                                                        &
!$OMP& PRIVATE(k, j, i, rhcritx, qc_points, rootwo, subgrid, qsl, qsl_ctt,    &
!$OMP& phiqcf, cosqcf, qn_imp, qn_adj, overlap_max, overlap_min,              &
!$OMP& overlap_random, temp0, temp1, temp2, qn, lqc, idx, qcfrbs, jj,         &
!$OMP& qsl_out, qsl_ctt_out)

!Cache-blocking applied to loop over j. It is used here to allow
!OpenMP parallelism to be over j (looping over k is order-dependent),
!but still have j and k in the correct order for good cache use.
!If jblock=rows=(tdims%j_len),
!then the inner "do j" loop does the most work, hence
!the loops over j and k are still the correct way around.
IF (ice_fraction_method  ==  cloud_top_temp) THEN
!$OMP DO SCHEDULE(DYNAMIC)
  DO jj = tdims%j_start, tdims%j_end, jblock
    DO j = jj, MIN((jj+jblock)-1,tdims%j_len)
      DO i = tdims%i_start, tdims%i_end
        llwic(i,j)=0
      END DO
    END DO

    DO k = levels, 1, -1
      DO j = jj, MIN((jj+jblock)-1,tdims%j_len)
        DO i = tdims%i_start, tdims%i_end

          IF (llwic(i,j)  /=  k+1) THEN
            ctt(i,j,k)=t(i,j,k)
          ELSE
            ctt(i,j,k)=ctt(i,j,k+1)
          END IF
          IF (qcf(i,j,k)  >   qcfmin) THEN
            llwic(i,j)=k
          END IF
        END DO
      END DO
    END DO
  END DO
!$OMP END DO
END IF

!$OMP DO SCHEDULE(DYNAMIC)
DO k = 1, levels

  ! ----------------------------------------------------------------------
  ! 1. Calculate QSAT at liquid/ice water temperature, TL, and initialize
  !    cloud water, sub-grid distribution and fraction arrays.
  !    This requires a preliminary calculation of the pressure.
  !    NB: On entry to the subroutine 'T' is TL and 'Q' is QW.
  ! ----------------------------------------------------------------------
  !
  !CDIR COLLAPSE
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      qcl(i,j,k) = 0.0
      cfl(i,j,k) = 0.0
      grid_qc(i,j,k) = rmdi
      bs(i,j,k) = rmdi
    END DO ! i
  END DO ! j

  IF ( l_new_qsat_lsc ) THEN
    IF ( l_mixing_ratio ) THEN
      CALL qsat_wat_mix_new(qsl,t(:,:,k),p_theta_levels(:,:,k),               &
            tdims%i_len,tdims%j_len)
    ELSE
      CALL qsat_wat_new(qsl,t(:,:,k),p_theta_levels(:,:,k),                   &
            tdims%i_len,tdims%j_len)
    END IF
  ELSE
    ! DEPENDS ON: qsat_wat_mix
    CALL qsat_wat_mix(qsl,t(1,1,k),p_theta_levels(1,1,k),             &
          tdims%i_len*tdims%j_len,l_mixing_ratio)
  END IF

  DO j = tdims%j_start, tdims%j_end

    DO i = tdims%i_start, tdims%i_end
      IF (multrhc==1) THEN
        rhcritx = rhcrit(i,j,k)
      ELSE
        rhcritx = rhcrit(1,1,k)
      END IF

      ! Omit CUMULUS points below (and including) NTML+1

      IF ( .NOT. l_param_conv .OR. (l_param_conv .AND.                      &
         (.NOT. cumulus(i,j) .OR. ( cumulus(i,j)                            &
         .AND. (k  >   ntml(i,j)+1) ))) ) THEN

        ! Rhcrit_if:
        IF (rhcritx  <   1.0) THEN
          ! -----------------------------------------------------------------
          ! 2. Calculate the quantity QN = QC/BS = (QW/QSL-1)/(1-RHcrit)
          !    if RHcrit is less than 1
          ! -----------------------------------------------------------------

          qn(i,j) = (q(i,j,k) / qsl(i,j) - 1.0) /                    &
                    (1.0 - rhcritx)

          ! -----------------------------------------------------------------
          ! 3. Set logical variable for cloud, LQC, for the case RHcrit < 1;
          !    where QN > -1, i.e. qW/qSAT(TL,P) > RHcrit, there is cloud
          ! -----------------------------------------------------------------

          lqc(i,j) = (qn(i,j)  >   -1.0)
        ELSE
          ! -----------------------------------------------------------------
          ! 2.a Calculate QN = QW - QSL if RHcrit equals 1
          ! -----------------------------------------------------------------

          qn(i,j) = q(i,j,k) - qsl(i,j)

          ! -----------------------------------------------------------------
          ! 3.a Set logical variable for cloud, LQC, for the case RHcrit = 1;
          !     where QN > 0, i.e. qW > qSAT(TL,P), there is cloud
          ! -----------------------------------------------------------------

          lqc(i,j) = (qn(i,j)  >   0.0)
        END IF ! Rhcrit_if
      ELSE IF (l_param_conv) THEN
        lqc(i,j) = .FALSE.
      END IF  ! Test on CUMULUS and NTML for A05_4A only
    END DO ! i
  END DO ! j

  ! ----------------------------------------------------------------------
  ! 4. Form index of points where non-zero liquid cloud fraction
  ! ----------------------------------------------------------------------

  qc_points=0

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF (lqc(i,j)) THEN
        qc_points = qc_points + 1
        idx(qc_points,1) = i
        idx(qc_points,2) = j
      END IF
    END DO ! i
  END DO ! j

  ! ----------------------------------------------------------------------
  ! 5. Call LS_CLD_C to calculate cloud water content, specific humidity,
  !                  water cloud fraction and determine temperature.
  ! ----------------------------------------------------------------------
  ! Qc_points_if:
  IF (qc_points  >   0) THEN
    ! DEPENDS ON: ls_cld_c
    CALL ls_cld_c(p_theta_levels(1,1,k),rhcrit(1,1,k),qsl,qn,       &
                  q(1,1,k),t(1,1,k),                                &
                  qcl(1,1,k),cfl(1,1,k),grid_qc(1,1,k),bs(1,1,k),   &
                  idx,qc_points,rhc_row_length,rhc_rows,            &
                  bl_levels,k, l_mixing_ratio)
  END IF ! Qc_points_if

END DO !k loop
!$OMP END DO
!$OMP BARRIER
  ! ----------------------------------------------------------------------
  ! 6. Calculate cloud fractions for ice clouds.
  !    THIS IS STILL HIGHLY EXPERIMENTAL.
  !    Begin by calculating Qsat_wat(T,P), at Temp. T, for estimate of bs.
  ! ----------------------------------------------------------------------
!$OMP MASTER
IF (ice_fraction_method == min_liq_overlap) THEN
  ALLOCATE(cfl_max(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,levels))

  k = levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      cfl_max(i,j,k) = cfl(i,j,k)
    END DO
  END DO
  DO k = levels-1, 1, -1
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (cfl(i,j,k) > cfl_max(i,j,k+1)) THEN
          ! Cloud frac has increased over peak, so use this
          cfl_max(i,j,k) = cfl(i,j,k)
        ELSE
          ! Cloud frac has decreased, so keep the same
          cfl_max(i,j,k) = cfl_max(i,j,k+1)
        END IF
      END DO
    END DO
  END DO

ELSE
  ALLOCATE(cfl_max(1,1,1))
END IF
!$OMP END MASTER
!$OMP BARRIER

rootwo = SQRT(2.0)

!$OMP DO SCHEDULE(DYNAMIC)
DO k = 1, levels

  IF ( l_new_qsat_lsc ) THEN
    IF ( l_mixing_ratio ) THEN
      CALL qsat_wat_mix_new(qsl,t(:,:,k),p_theta_levels(:,:,k),               &
            tdims%i_len,tdims%j_len)
    ELSE
      CALL qsat_wat_new(qsl,t(:,:,k),p_theta_levels(:,:,k),                   &
            tdims%i_len,tdims%j_len)
    END IF
  ELSE
    ! DEPENDS ON: qsat_wat_mix
    CALL qsat_wat_mix(qsl,t(1,1,k),p_theta_levels(1,1,k),                     &
          tdims%i_len*tdims%j_len,l_mixing_ratio)
  END IF

  IF (ice_fraction_method  ==  cloud_top_temp) THEN
    ! Use cloud top temperature and a fixed qsat to give QCFRBS

    IF ( l_new_qsat_lsc ) THEN
      IF ( l_mixing_ratio ) THEN
        CALL qsat_wat_mix_new(qsl_ctt,ctt(:,:,k),p_theta_levels(:,:,k),       &
                              tdims%i_len,tdims%j_len)
      ELSE
        CALL qsat_wat_new(qsl_ctt,ctt(:,:,k),p_theta_levels(:,:,k),           &
                              tdims%i_len,tdims%j_len)
      END IF
    ELSE
      ! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsl_ctt,ctt(1,1,k),p_theta_levels(1,1,k),             &
          tdims%i_len*tdims%j_len,l_mixing_ratio)
    END IF

    CALL powr_v(                                                    &
       tdims%i_len*tdims%j_len,                                     &
       qsl_ctt,ctt_weight,qsl_ctt_out )
    qsl_ctt=qsl_ctt_out
    CALL powr_v(                                                    &
       tdims%i_len*tdims%j_len,                                     &
       qsl,t_weight,qsl_out )
    qsl=qsl_out

    subgrid = sub_cld ** (1.0-t_weight)                             &
              / qsat_fixed ** (1.0-t_weight-ctt_weight)
  END IF ! ice_fraction_method eq 2

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF ( multrhc== 1) THEN
        rhcritx = rhcrit(i,j,k)
      ELSE
        rhcritx = rhcrit(1,1,k)
      END IF
      ! ----------------------------------------------------------------------
      ! 6a Calculate qCF/bs.
      ! ----------------------------------------------------------------------
      ! Rhcrit_if2:
      IF (rhcritx  <   1.0) THEN

        IF (ice_fraction_method  ==  smith_orig) THEN
          qcfrbs(i,j)=  qcf(i,j,k) / ((1.0-rhcritx) * qsl(i,j))
        ELSE IF (ice_fraction_method  ==  cloud_top_temp) THEN
          qcfrbs(i,j) = subgrid * qcf(i,j,k) / ((1.0-rhcritx)        &
                   * qsl_ctt(i,j)*qsl(i,j))
        ELSE IF (ice_fraction_method == min_liq_overlap) THEN
          ! method described in appendix of Abel et al (2017, JAS)
          qcfrbs(i,j)=  MAX(1.0 - cfl_max(i,j,k), 0.05) * qcf(i,j,k) /  &
                        ((1.0-rhcritx) * qsl(i,j))
        ELSE
          ! No ice cloud fraction method defined
        END IF ! ice_fraction_method

        ! ----------------------------------------------------------------------
        ! 6b Calculate frozen cloud fraction from frozen cloud water content.
        ! ----------------------------------------------------------------------
        IF (qcfrbs(i,j)  <=  0.0) THEN
          cff(i,j,k) = 0.0
        ELSE IF (0.0<qcfrbs(i,j) .AND. (6.0*qcfrbs(i,j)) <= 1.0) THEN
          cff(i,j,k) = 0.5 * ((6.0 * qcfrbs(i,j))**(2.0/3.0))
        ELSE IF (1.0<(6.0*qcfrbs(i,j)) .AND. qcfrbs(i,j) < 1.0) THEN
          phiqcf = ACOS(rootwo * 0.75 * (1.0 - qcfrbs(i,j)))
          cosqcf = COS((phiqcf + (4.0 * pi)) / 3.0)
          cff(i,j,k) = 1.0 - (4.0 * cosqcf * cosqcf)
        ELSE IF (qcfrbs(i,j)  >=  1.0) THEN
          cff(i,j,k) = 1.0
        END IF
        IF (i_eacf == all_clouds .OR.                                   &
             (i_eacf == not_mixph .AND. cfl(i,j,k) < 0.05) ) THEN
          ! Empirically adjusted cloud fraction
          ! Back out QN
          IF (0.0< qcfrbs(i,j) .AND. (6.0*qcfrbs(i,j)) <=  1.0) THEN
            qn_imp=SQRT(2.0*cff(i,j,k))-1.0
          ELSE IF (1.0<(6.0*qcfrbs(i,j)) .AND. qcfrbs(i,j)<1.0) THEN
            qn_imp=1.0-SQRT((1.0-cff(i,j,k))*2.0)
          ELSE
            qn_imp = 1.0
          END IF

          ! Modify QN with EACF relationship
          IF (k >  bl_levels) THEN
            qn_adj=(qn_imp+0.0955)/(1.0-0.0955)
          ELSE
            qn_adj=(qn_imp+0.184)/(1.0-0.184)
          END IF

          ! Recalculate ice cloud fraction with modified QN
          IF (qcfrbs(i,j)  <=  0.0) THEN
            cff(i,j,k) = 0.0
          ELSE IF (qn_adj  <=  0.0) THEN
            cff(i,j,k) = 0.5 * (1.0 + qn_adj) * (1.0 + qn_adj)
          ELSE IF (qn_adj  <   1.0) THEN
            cff(i,j,k) = 1.0 - 0.5 * (1.0-qn_adj) * (1.0-qn_adj)
          ELSE
            cff(i,j,k) = 1.0
          END IF

        END IF  ! i_eacf


      ELSE ! RHcrit = 1, set cloud fraction to 1 or 0

        IF (qcf(i,j,k)  >   0.0) THEN
          cff(i,j,k) = 1.0
        ELSE
          cff(i,j,k) = 0.0
        END IF

      END IF
    END DO ! i
  END DO ! j

  ! ----------------------------------------------------------------------
  ! 6c Calculate combined cloud fraction.
  ! ----------------------------------------------------------------------

  IF (cloud_fraction_method  ==  1) THEN

    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        !             Use minimum overlap condition
        cf(i,j,k) = MIN(cfl(i,j,k)+cff(i,j,k), 1.0)
      END DO
    END DO

  ELSE IF (cloud_fraction_method  ==  2) THEN

      !Hand-unrolled loop
    DO j = tdims%j_start, tdims%j_len -                             &
              MOD(tdims%j_len,2), 2

      DO i = tdims%i_start, tdims%i_end
        ! Calculate possible overlaps between ice and liquid in THIS layer
        overlap_max=MIN(cfl(i,j,k),cff(i,j,k))
        overlap_min=MAX(cfl(i,j,k)+cff(i,j,k)-1.0,0.0)
        overlap_random=cfl(i,j,k)*cff(i,j,k)
        ! Now use the specified degree of overlap to calculate the total
        ! cloud fraction (= cfice + cfliq - overlap)
        temp0=overlap_random
        temp1=0.5*(overlap_max-overlap_min)
        temp2=0.5*(overlap_max+overlap_min)-overlap_random
        cf(i,j,k)=cfl(i,j,k)+cff(i,j,k)                             &
                -(temp0+temp1*overlap_ice_liquid                    &
                +temp2*overlap_ice_liquid*overlap_ice_liquid)
        ! Check that the overlap wasnt negative
        cf(i,j,k)=MIN(cf(i,j,k),cfl(i,j,k)+cff(i,j,k))

      END DO

      DO i = tdims%i_start, tdims%i_end
        ! Calculate possible overlaps between ice and liquid in THIS layer
        overlap_max=MIN(cfl(i,j+1,k),cff(i,j+1,k))
        overlap_min=MAX(cfl(i,j+1,k)+cff(i,j+1,k)-1.0,0.0)
        overlap_random=cfl(i,j+1,k)*cff(i,j+1,k)
        ! Now use the specified degree of overlap to calculate the total
        ! cloud fraction (= cfice + cfliq - overlap)
        temp0=overlap_random
        temp1=0.5*(overlap_max-overlap_min)
        temp2=0.5*(overlap_max+overlap_min)-overlap_random
        cf(i,j+1,k)=cfl(i,j+1,k)+cff(i,j+1,k)                       &
                -(temp0+temp1*overlap_ice_liquid                    &
                +temp2*overlap_ice_liquid*overlap_ice_liquid)
        ! Check that the overlap wasnt negative
        cf(i,j+1,k)=MIN(cf(i,j+1,k),cfl(i,j+1,k)+cff(i,j+1,k))

      END DO
    END DO

      !Post-conditioning
    IF (MOD(tdims%j_len,2) == 1) THEN
      DO i = tdims%i_start, tdims%i_end
        ! Calculate possible overlaps between ice and liquid in THIS layer
        overlap_max=MIN( cfl(i,tdims%j_len,k),                      &
                         cff(i,tdims%j_len,k) )
        overlap_min=MAX( cfl(i,tdims%j_len,k) +                     &
                         cff(i,tdims%j_len,k)-1.0,                  &
                         0.0 )
        overlap_random=  cfl(i,tdims%j_len,k) *                     &
                         cff(i,tdims%j_len,k)
        ! Now use the specified degree of overlap to calculate the total
        ! cloud fraction (= cfice + cfliq - overlap)
        temp0=overlap_random
        temp1=0.5*(overlap_max-overlap_min)
        temp2=0.5*(overlap_max+overlap_min)-overlap_random
        cf(i,tdims%j_len,k)=                                        &
          cfl(i,tdims%j_len,k) +                                    &
          cff(i,tdims%j_len,k) -                                    &
          (temp0+temp1*overlap_ice_liquid +                         &
          temp2*overlap_ice_liquid*overlap_ice_liquid)
        ! Check that the overlap wasnt negative
        cf(i,tdims%j_len,k)=                                        &
          MIN(cf(i,tdims%j_len,k),                                  &
             cfl(i,tdims%j_len,k) +                                 &
             cff(i,tdims%j_len,k))

      END DO
    END IF !Post-conditioning

    ! CFF + CFL >= 1 implies CF = 1.
    ! Deal with this case separately to avoid roundoff issues:
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (cfl(i,j,k)+cff(i,j,k) >= 1.0) cf(i,j,k) = 1.0
      END DO
    END DO

  ELSE
    ! No total cloud fraction method defined


  END IF ! cloud_fraction_method

END DO ! Levels_do
!$OMP END DO

!$OMP END PARALLEL

DEALLOCATE(cfl_max)

9999 CONTINUE ! Error exit
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ls_cld

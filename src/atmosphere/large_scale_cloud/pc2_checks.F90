! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! PC2 Cloud Scheme: Checking cloud parameters
MODULE pc2_checks_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'PC2_CHECKS_MOD'
CONTAINS

SUBROUTINE pc2_checks(                                                  &
!   Pressure related fields
 p_theta_levels,                                                        &
!   Prognostic Fields
 t, cf, cfl, cff, q, qcl, qcf,                                          &
!   Logical control
 l_mixing_ratio)

USE mphys_ice_mod,         ONLY: thomo
USE planet_constants_mod,  ONLY: lcrcp, lfrcp, lsrcp, r, repsilon
USE conversions_mod,       ONLY: zerodegc
USE water_constants_mod,   ONLY: lc
USE pc2_constants_mod,     ONLY: cloud_rounding_tol,                  &
                                 one_over_qcf0,                       &
                                 min_in_cloud_qcf,                    &
                                 one_over_min_in_cloud_qcf,           &
                                 condensate_limit,                    &
                                 wcgrow

USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: pdims, tdims
USE cloud_inputs_mod,      ONLY: i_pc2_checks_cld_frac_method,                 &
                                 l_ensure_min_in_cloud_qcf

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_wat_new     => qsat_wat,                       &
                    qsat_wat_mix_new => qsat_wat_mix,                   &
                    l_new_qsat_lsc !Currently defaults to FALSE

IMPLICIT NONE

! Purpose:
!   This subroutine checks that cloud fractions, liquid water and
!   vapour contents take on physically reasonable values.
!
! Method:
!   Apply checks sequentially to the input values. It is more important
!   to ensure that liquid does not go negative than to ensure that
!   saturation deficit is correct.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation

!  Subroutine Arguments:------------------------------------------------

REAL ::                                                               &
                      !, INTENT(IN)
 p_theta_levels(tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                1:tdims%k_end)
!    pressure at all points (Pa)

LOGICAL ::                                                            &
                      !, INTENT(IN)
 l_mixing_ratio
!    Use mixing ratio formulation

REAL ::                                                               &
                      !, INTENT(INOUT)
 t(  tdims%i_start:tdims%i_end,                                       &
     tdims%j_start:tdims%j_end,                                       &
                 1:tdims%k_end),                                      &
!    Temperature (K)
   cf( tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,                                       &
                   1:tdims%k_end),                                      &
!    Total cloud fraction (no units)
   cfl(tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,                                       &
                   1:tdims%k_end),                                      &
!    Liquid cloud fraction (no units)
   cff(tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,                                       &
                   1:tdims%k_end),                                      &
!    Ice cloud fraction (no units)
   q(  tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,                                       &
                   1:tdims%k_end),                                      &
!    Vapour content (kg water per kg air)
   qcl(tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,                                       &
                   1:tdims%k_end),                                      &
!    Liquid content (kg water per kg air)
   qcf(tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,                                       &
                   1:tdims%k_end)
!    Ice content (kg water per kg air)

!  External functions:

!  Local parameters and other physical constants------------------------
!
! Options for i_pc2_checks_cld_frac_method
INTEGER, PARAMETER :: original = 0
! Original. No change to liquid cloud fraction (CFL) is made when QCL is
! increased when qv>qsat.
INTEGER, PARAMETER :: force_cfl_cf_unity = 1
! Force CFL (and CF) to unity when creating some extra QCL when qv>qsat.
INTEGER, PARAMETER :: set_cfl_via_gi = 3
! Set CFL based on QCL value using Gultepe and Isaac (gi).

!  Local scalars--------------------------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
REAL ::                                                               &
 alpha,                                                               &
!      Rate of change of saturation specific humidity with
!      temperature calculated at dry-bulb temperature
!      (kg kg-1 K-1)
   al,                                                                  &
!      1 / (1 + alpha L/cp)  (no units)
   sd,                                                                  &
!      Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)
   cfl_old
!      temp store for old cfl value

!  (b) Others.
INTEGER :: k,i,j       ! Loop counters: K - vertical level index
!                           I,J - horizontal position index

!  Local dynamic arrays-------------------------------------------------
!    1 block of real workspace is required.
REAL ::                                                               &
 qsl_t(tdims%i_start:tdims%i_end,                                     &
       tdims%j_start:tdims%j_end)
!       Saturated specific humidity for dry bulb temperature T

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL   (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PC2_CHECKS'

!- End of Header

! ==Main Block==--------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Loop round levels to be processed
! Levels_do1:

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k, al,&
!$OMP& alpha, sd, cfl_old, qsl_t)
DO k = 1,tdims%k_end

  ! ----------------------------------------------------------------------
  ! 1. Calculate Saturated Specific Humidity with respect to liquid water
  !    for dry bulb temperatures.
  ! ----------------------------------------------------------------------

  IF ( l_new_qsat_lsc ) THEN
    IF ( l_mixing_ratio ) THEN
      CALL qsat_wat_mix_new(qsl_t,t(:,:,k),p_theta_levels(:,:,k),             &
            tdims%i_len,tdims%j_len)
    ELSE
      CALL qsat_wat_new(qsl_t,t(:,:,k),p_theta_levels(:,:,k),                 &
            tdims%i_len,tdims%j_len)
    END IF
  ELSE
    ! DEPENDS ON: qsat_wat_mix
    CALL qsat_wat_mix(qsl_t,t(tdims%i_start,tdims%j_start,k),                 &
          p_theta_levels(tdims%i_start,tdims%j_start,k),                      &
          tdims%i_len*tdims%j_len,                                            &
          l_mixing_ratio)
  END IF

  ! Rows_do1:
  DO j = tdims%j_start, tdims%j_end
    ! Row_length_do1:
    DO i = tdims%i_start, tdims%i_end

      !----------------------------------------------------------------------
      ! 2. Calculate the saturation deficit.
      ! ----------------------------------------------------------------------

      ! Need to estimate the rate of change of saturated specific humidity
      ! with respect to temperature (alpha) first, then use this to calculate
      ! factor aL.
      alpha=repsilon*lc*qsl_t(i,j)/(r*t(i,j,k)**2)
      al=1.0/(1.0+lcrcp*alpha)

      ! Calculate the saturation deficit SD

      sd=al*(qsl_t(i,j)-q(i,j,k))

      ! ----------------------------------------------------------------------
      !  3. Checks are applied here for liquid cloud
      ! ----------------------------------------------------------------------

      ! Earlier versions checked whether saturation deficit is zero (or less
      ! than zero). If so, then the liquid cloud fraction was forced to one.
      ! This check has been suspended for numerical reasons.
      !            IF (SD  <=  0.0 .OR. CFL(i,j,k)  >   1.0) THEN

      ! Instead, check simply whether input values of liquid cloud fraction
      ! are, or are between, zero and one. If not, adjust them to zero or one.
      ! Adjust the total cloud fractions accordingly.

      IF (cfl(i,j,k) > (1.0 - cloud_rounding_tol)) THEN
        cfl(i,j,k)=1.0
        cf(i,j,k) =1.0
      END IF

      ! Check also whether the liquid water content is less than zero, and
      ! set liquid cloud fraction to zero if it is.

      IF (qcl(i,j,k) < condensate_limit .OR.                          &
          cfl(i,j,k) < cloud_rounding_tol) THEN
        cfl(i,j,k)=0.0
        cf(i,j,k) =cff(i,j,k)
      END IF

      ! Check whether the saturation deficit is less than zero. If it is
      ! then condense some liquid to bring the saturation deficit to zero.
      ! Adjust the temperature for the latent heating.

      IF (sd < 0.0) THEN
        q(i,j,k)   = q(i,j,k)   + sd
        qcl(i,j,k) = qcl(i,j,k) - sd
        t(i,j,k)   = t(i,j,k)   - sd * lcrcp

        ! In original PC2 (i_pc2_checks_cld_frac_method=original)
        ! there is no change to the cloud fraction
        ! as a result of the above QCL change.
        IF (i_pc2_checks_cld_frac_method /= original) THEN
          ! Various options for adjusting the cloud fraction.
          IF (i_pc2_checks_cld_frac_method == force_cfl_cf_unity) THEN
            ! Force CFL to 1 and hence also set CF to 1.
            cfl(i,j,k) = 1.0
            cf(i,j,k)  = 1.0
          ELSE
            ! Increase the cloud fraction in order to
            ! keep the in-cloud condensate amount the same
            ! providing there is a well-defined in cloud value.
            cfl_old = cfl(i,j,k)
            IF (cfl(i,j,k) > 0.0 .AND. (qcl(i,j,k)+sd) > 0.0) THEN
               ! There is a well defined in cloud value
              cfl(i,j,k) = MAX(0.0,MIN(1.0,                          &
                           qcl(i,j,k)*cfl(i,j,k)/(qcl(i,j,k)+sd) ))
            ELSE

              IF (i_pc2_checks_cld_frac_method >= set_cfl_via_gi) THEN
                ! Initiate cloud using Gultepe and Isaac formulation.
                ! NB 5.57 is valid for  10 km gridbox.
                !    4.37 is valid for 100 km gridbox.
                ! For now just use a parameter value valid for dx=10km
                ! as variation with gridbox size is not huge and the actual
                ! value of qcl is probably not perfect anyway, 
                ! since just found from excess over 100% RH.
                ! Ensure that cloud cover is at least 0.5 though, since
                ! we are dealing with cases with RH>100%.
                cfl(i,j,k) = MAX(MIN(5.57*((qcl(i,j,k)*1000.0)**0.775),1.0),0.5)
              ELSE
                 ! grow cloud at a fixed amount
                cfl(i,j,k) = MAX(0.0,MIN(1.0,                          &
                             qcl(i,j,k)/wcgrow ))
              END IF

            END IF
            ! update total cloud fraction
            cf(i,j,k) = MAX(MIN(cf(i,j,k)+cfl(i,j,k)-cfl_old,1.0),0.0)
          END IF
        END IF ! i_pc2_checks_cld_frac_method /= original

      END IF ! sd < 0.0

      ! Check whether the saturation deficit
      ! is greater than zero but the liquid cloud fraction is one. If it is
      ! then evaporate some liquid (provided there is enough
      ! liquid) to bring the saturation deficit to zero. Adjust the
      ! temperature for the latent heating.

      IF (sd > 0.0 .AND. cfl(i,j,k) == 1.0 .AND. qcl(i,j,k) > sd) THEN
        q(i,j,k)   = q(i,j,k)   + sd
        qcl(i,j,k) = qcl(i,j,k) - sd
        t(i,j,k)   = t(i,j,k)   - sd * lcrcp
      END IF

      ! Check whether the liquid content is less than zero, or whether it is
      ! greater than zero but the liquid cloud fraction is zero. If so then
      ! condense or evaporate liquid to bring the liquid water to zero. Adjust
      ! the temperature for latent heating.

      IF (qcl(i,j,k) < condensate_limit .OR.                          &
         (qcl(i,j,k) >  0.0 .AND. cfl(i,j,k) == 0.0) ) THEN
        q(i,j,k)   = q(i,j,k) + qcl(i,j,k)
        t(i,j,k)   = t(i,j,k) - qcl(i,j,k) * lcrcp
        qcl(i,j,k) = 0.0
      END IF

      ! ----------------------------------------------------------------------
      !  4. Check that ice content and ice cloud fraction are sensible.
      ! ----------------------------------------------------------------------

      ! Check whether ice content is zero (or less than zero). If so then
      ! force the ice cloud fraction to zero. Also check whether input values
      ! of ice cloud fraction are, or are between, zero and one. If not,
      ! adjust them to zero or one. Adjust the total cloud fractions
      ! accordingly.

      IF (cff(i,j,k) > (1.0 - cloud_rounding_tol)) THEN
        cff(i,j,k) = 1.0
        cf(i,j,k)  = 1.0
      END IF

      IF (qcf(i,j,k) < condensate_limit .OR.                          &
          cff(i,j,k) < cloud_rounding_tol) THEN
        cff(i,j,k) = 0.0
        cf(i,j,k)  = cfl(i,j,k)
      END IF

      ! If ice content is negative then condense some vapour to remove the
      ! negative part. Adjust the temperature for the latent heat.

      IF (qcf(i,j,k) < condensate_limit) THEN
        q(i,j,k)   = q(i,j,k) + qcf(i,j,k)
        t(i,j,k)   = t(i,j,k) - qcf(i,j,k) * lsrcp
        qcf(i,j,k) = 0.0
      END IF

      ! If ice content is positive but ice cloud fraction negative, create
      ! some ice cloud fraction.

      IF (qcf(i,j,k) > 0.0 .AND. cff(i,j,k) == 0.0) THEN
        cff(i,j,k) = qcf(i,j,k) * one_over_qcf0
        !             CF is not adjusted here but is checked below
      END IF

      IF (l_ensure_min_in_cloud_qcf) THEN
        IF (cff(i,j,k) > 0.0) THEN
          IF (qcf(i,j,k)/cff(i,j,k) < min_in_cloud_qcf) THEN
            cff(i,j,k)=qcf(i,j,k)*one_over_min_in_cloud_qcf
          END IF
        END IF
      END IF

      ! ----------------------------------------------------------------------
      !  5. Check that total cloud fraction is sensible.
      ! ----------------------------------------------------------------------

      ! Total cloud fraction must be bounded by
      ! i) The maximum of the ice and liquid cloud fractions (maximum overlap)

      IF (cf(i,j,k) < MAX(cfl(i,j,k),cff(i,j,k))) THEN
        cf(i,j,k) = MAX(cfl(i,j,k),cff(i,j,k))
      END IF

      ! ii) The sum of the ice and liquid cloud fractions or one, whichever
      ! is the lower (minimum overlap)

      IF (cf(i,j,k) > MIN( (cfl(i,j,k)+cff(i,j,k)),1.0 ) ) THEN
        cf(i,j,k) = MIN( (cfl(i,j,k)+cff(i,j,k)),1.0 )
      END IF

      ! ----------------------------------------------------------------------
      !  6. Homogeneous nucleation.
      ! ----------------------------------------------------------------------

      IF (t(i,j,k) < (zerodegc+thomo)) THEN
        ! Turn all liquid to ice
        IF (qcl(i,j,k) > 0.0) THEN
          cff(i,j,k) = cf(i,j,k)
          cfl(i,j,k) = 0.0
        END IF
        qcf(i,j,k) = qcf(i,j,k) + qcl(i,j,k)
        t(i,j,k)   = t(i,j,k)   + qcl(i,j,k) * lfrcp
        qcl(i,j,k) = 0.0
      END IF

    END DO !i

  END DO !j

END DO !k
!$OMP END PARALLEL DO

! End of the subroutine

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pc2_checks
END MODULE pc2_checks_mod

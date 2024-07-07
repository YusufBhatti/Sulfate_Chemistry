! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  PC2 Cloud Scheme: Initiation

MODULE pc2_initiate_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'PC2_INITIATE_MOD'
CONTAINS

SUBROUTINE pc2_initiate(                                                &
!      Pressure related fields
 p_theta_levels, ccb, cumulus, rhcrit,                                  &
!      Array dimensions
 nlevels,                                                               &
 rhc_row_length,rhc_rows,zlcl_mixed,r_theta_levels,                     &
!      Prognostic Fields
   t, cf, cfl, cff, q, qcl, rhts,                                       &
!      Logical control
   l_mixing_ratio)

USE conversions_mod,       ONLY: zerodegc
USE water_constants_mod,   ONLY: lc
USE planet_constants_mod,  ONLY: lcrcp, r, repsilon
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: pdims, tdims, pdims_l
USE pc2_constants_mod,     ONLY: init_iterations, rhcrit_tol,         &
                                 pdf_power, rhcpt_tke_based,          &
                                 cloud_pc2_tol
USE cv_run_mod,            ONLY: l_param_conv
USE cloud_inputs_mod,      ONLY: i_rhcpt

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_wat_new     => qsat_wat,                       &
                    qsat_wat_mix_new => qsat_wat_mix,                   &
                    l_new_qsat_lsc !Currently defaults to FALSE

USE pc2_total_cf_mod, ONLY: pc2_total_cf
IMPLICIT NONE

! Purpose:
!   Initiate liquid and total cloud fraction and liquid water content

! Method:
!   Uses the method proposed in Annex C of the PC2 cloud scheme project
!   report, which considers a Smith-like probability density
!   distribution whose width is given by a prescribed value of RHcrit.
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
INTEGER ::                                                            &
                      !, INTENT(IN)
 nlevels,                                                             &
!       No. of levels being processed.
    rhc_row_length,rhc_rows,                                            &
!       Dimensions of the RHCRIT variable.
   ccb(           tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end)
!       Convective cloud base

LOGICAL ::                                                            &
                      !, INTENT(IN)
 cumulus(       tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end),                           &
!       Is this a boundary layer cumulus point
   l_mixing_ratio  ! Use mixing ratio formulation

! height of lcl in a well-mixed BL (types 3 or 4), 0 otherwise
REAL :: zlcl_mixed(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

REAL ::                                                               &
                      !, INTENT(IN)
 p_theta_levels(pdims%i_start:pdims%i_end,                            &
                pdims%j_start:pdims%j_end,                            &
                nlevels),                                             &
!       Pressure at all points (Pa)
   r_theta_levels(pdims_l%i_start:pdims_l%i_end,                        &
                  pdims_l%j_start:pdims_l%j_end,                        &
                  0:nlevels),                                           &
!       Pressure at all points (Pa)
   rhcrit(rhc_row_length,rhc_rows,nlevels),                             &
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.
   cff(           tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  nlevels)
!       Ice cloud fraction (no units)

REAL ::                                                               &
                      !, INTENT(INOUT)
 t(             tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                nlevels),                                             &
!       Temperature (K)
   cf(            tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  nlevels),                                             &
!       Total cloud fraction (no units)
   cfl(           tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  nlevels),                                             &
!       Liquid cloud fraction (no units)
   q(             tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  nlevels),                                             &
!       Vapour content (kg water per kg air)
   qcl(           tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  nlevels),                                             &
!       Liquid content (kg water per kg air)
   rhts(          tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  nlevels)
!       Variable carrying initial RHT wrt TL from start of timestep

!  External functions:

!  Local scalars--------------------------------------------------------

!  (a)  Scalars effectively expanded to workspace by the Cray (using
!       vector registers).
REAL ::                                                               &
 alpha,                                                               &
!       Rate of change of saturation specific humidity with
!       temperature calculated at dry-bulb temperature
!       (kg kg-1 K-1)
   al,                                                                  &
!       1 / (1 + alpha L/cp)  (no units)
   bs,                                                                  &
!       Width of distribution (kg kg-1)
   c_1,                                                                 &
!       New liquid cloud fraction
   deltal,                                                              &
!       Change in liquid (kg kg-1)
   descent_factor,                                                      &
!       Rate at which to relax to new liquid water content
   l_bs,                                                                &
!       New liquid water content divided by PDF width (no units)
   l_out,                                                               &
!       New liquid water content (kg kg-1)
   q_out,                                                               &
!       New vapour content (kg kg-1)
   qc,                                                                  &
!       aL (q + l - qsat(TL) )  (kg kg-1)
   qn,                                                                  &
!       Normalized value of QC for the probability density function.
   rht,                                                                 &
!       Total relative humidity (liquid+vapour)/qsat (kg kg-1)
   rh0,                                                                 &
!       Equivalent critical relative humidity
   sd
!       Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)

!  (b)  Others.
INTEGER :: k,i,j,l,                                                   &
!       Loop counters: K   - vertical level index
!                      I,J - horizontal position index
!                      L   - counter for iterations
          irhi,irhj,                                                    &
!                      Indices for RHcrit array
          multrhc,                                                      &
!                      Zero if (rhc_row_length*rhc_rows) le 1, else 1
          npti, npt
!                      Number of points to iterate over

!  Local dynamic arrays-------------------------------------------------
REAL ::                                                                 &
 qsl_tl(    tdims%i_len*                                                &
            tdims%j_len ),                                              &
!       Saturated specific humidity for liquid temperature TL
   tl_c(    tdims%i_len*                                                &
            tdims%j_len ),                                              &
!       Liquid temperature (= T - L/cp QCL)  (K)
   p_c(     tdims%i_len*                                                &
            tdims%j_len ),                                              &
   cf_c(    tdims%i_len*tdims%j_len),                                   &
!       Total cloud fraction on compressed points (no units)
   cfl_c(   tdims%i_len*tdims%j_len),                                   &
!       Liquid cloud fraction on compressed points (no units)
   cff_c(   tdims%i_len*tdims%j_len),                                   &
!       Ice cloud fraction on compressed points (no units)
   deltacl_c(                                                           &
            tdims%i_len*tdims%j_len),                                   &
!       Change in liquid cloud fraction (no units)
   deltacf_c(                                                           &
            tdims%i_len*tdims%j_len),                                   &
!       Change in ice cloud fraction (no units)
   p_theta_levels_c(                                                    &
            tdims%i_len*tdims%j_len),                                   &
!       Pressure on compressed points (Pa)
   q_c(     tdims%i_len*tdims%j_len),                                   &
!       Vapour content on compressed points (kg kg-1)
   qcl_c(   tdims%i_len*tdims%j_len),                                   &
!       Liquid water content on compressed points (kg kg-1)
   qn_c(    tdims%i_len*tdims%j_len),                                   &
!      QN on compressed points (no units)
   qsl_t_c( tdims%i_len*tdims%j_len),                                   &
!       Saturated specific humidity for dry-bulb temperature T on
!       compressed points (kg kg-1)
   qsl_tl_c(                                                            &
            tdims%i_len*tdims%j_len),                                   &
!       Saturated specific humidity for liquid temperature TL on
!       compressed points (kg kg-1)
   rh0_c(   tdims%i_len*tdims%j_len),                                   &
!       Equivalent critical relative humidity on compressed points
!       (no units)
   t_c(     tdims%i_len*tdims%j_len)
!       Temperature on compressed points (K)

INTEGER ::                                                            &
 ni(      tdims%i_len*tdims%j_len),                                   &
 nj(      tdims%i_len*tdims%j_len),                                   &
!       Compressed point counters
   ind_i(   tdims%i_len*tdims%j_len),                                 &
   ind_j(   tdims%i_len*tdims%j_len)
!       Compressed point counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL   (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PC2_INITIATE'

LOGICAL, PARAMETER :: l_simplify_pc2_init_logic=.FALSE.

!- End of Header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up a flag to state whether RHcrit is a single parameter or defined
! on all points.

IF (rhc_row_length*rhc_rows > 1) THEN
  multrhc=1
ELSE
  multrhc=0
END IF

! ==Main Block==--------------------------------------------------------

! Loop round levels to be processed
! Levels_do1:

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(npt,        &
!$OMP& npti, ind_i, ind_j, irhj, irhi, rht, rh0, qn, c_1, ni, nj, t_c,  &
!$OMP& qn_c, rh0_c, qsl_tl_c, cf_c, cfl_c, cff_c, p_theta_levels_c,     &
!$OMP& qcl_c, q_c, deltacl_c, deltacf_c, p_c, qsl_t_c, l_out, l_bs, al, &
!$OMP& deltal, descent_factor, q_out, bs, i, j, k, l, qsl_tl, tl_c,     &
!$OMP& alpha)
DO k = 1, nlevels

  IF (l_simplify_pc2_init_logic) THEN

    ! Copy points into compressed arrays
    npt = 0
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF ( ( cfl(i,j,k) < cloud_pc2_tol                               &
              .AND. (.NOT. l_param_conv .OR. .NOT. cumulus(i,j))        &
              .AND. (r_theta_levels(i,j,k)-r_theta_levels(i,j,0))       &
                    > zlcl_mixed(i,j) )                                 &
            .OR. cfl(i,j,k) > 1.0 - cloud_pc2_tol ) THEN
          npt = npt+1
          ind_i(npt) = i
          ind_j(npt) = j
          tl_c(npt)=t(i,j,k)-lcrcp*qcl(i,j,k)
          p_c(npt) = p_theta_levels(i,j,k)
        END IF
      END DO
    END DO

  ELSE

    ! Copy points into compressed arrays
    npt = 0
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF ( ( (cfl(i,j,k)  ==  0.0                                       &
           .OR. (cfl(i,j,k)  <   0.05 .AND. t(i,j,k)  <   zerodegc))      &
           .AND. (.NOT. cumulus(i,j))                                     &
           .AND. (r_theta_levels(i,j,k)-r_theta_levels(i,j,0))            &
                    > zlcl_mixed(i,j) ) .OR.                              &
                      cfl(i,j,k) == 1.0 ) THEN
          npt = npt+1
          ind_i(npt) = i
          ind_j(npt) = j
          tl_c(npt)=t(i,j,k)-lcrcp*qcl(i,j,k)
          p_c(npt) = p_theta_levels(i,j,k)
        END IF
      END DO
    END DO

  END IF ! (l_simplify_pc2_init_logic)

  ! ----------------------------------------------------------------------
  ! 2. Calculate Saturated Specific Humidity with respect to liquid water
  !    for liquid temperatures.
  ! ----------------------------------------------------------------------
  IF (npt > 0) THEN
    IF ( l_new_qsat_lsc ) THEN
      IF ( l_mixing_ratio ) THEN
        CALL qsat_wat_mix_new(qsl_tl, tl_c, p_c, npt)
      ELSE
        CALL qsat_wat_new(qsl_tl, tl_c, p_c, npt)
      END IF
    ELSE
      ! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsl_tl, tl_c, p_c, npt, l_mixing_ratio)
    END IF
  END IF
  ! Set number of points to perform the compressed calculation over
  ! to zero

  npti=0.0

  DO i = 1, npt

    ! ----------------------------------------------------------------------
    ! 3. Calculate equivalent critical relative humidity values (RH0) for
    !    initiation.
    ! ----------------------------------------------------------------------

    ! Set up index pointers to critical relative humidity value

    irhi = (multrhc * (ind_i(i) - 1)) + 1
    irhj = (multrhc * (ind_j(i) - 1)) + 1

    ! Calculate relative total humidity with respect to the liquid
    ! temperature

    rht=(q(ind_i(i),ind_j(i),k)+qcl(ind_i(i),ind_j(i),k))           &
        /qsl_tl(i)

    ! Initialize the equivalent critical relative humidity to a dummy
    ! negative value

    rh0 = -1.0

    IF (l_simplify_pc2_init_logic) THEN

      ! Do we need to force some initiation from zero liquid cloud amount?
      ! If so, set RH0 equal the critical relative humidity. 

      IF (      cfl(ind_i(i),ind_j(i),k) < cloud_pc2_tol                       &
          .AND. (.NOT. l_param_conv .OR. .NOT. cumulus(ind_i(i),ind_j(i)) )    &
          .AND. (r_theta_levels(ind_i(i),ind_j(i),k)                           &
                -r_theta_levels(ind_i(i),ind_j(i),0))                          &
                > zlcl_mixed(ind_i(i),ind_j(i))                                &
          .AND. rht > rhts(ind_i(i),ind_j(i),k)                                &
          .AND. rht > (rhcrit(irhi,irhj,k)+rhcrit_tol)  ) THEN
        rh0 = rhcrit(irhi,irhj,k)
      END IF
 
     ! Do we need to force some initiation from total liquid cloud cover?
     ! If so, set RH0 equal the critical relative humidity.
 
     IF (      cfl(ind_i(i),ind_j(i),k) > 1.0 - cloud_pc2_tol                  &
         .AND. rht < rhts(ind_i(i),ind_j(i),k)                                 &
         .AND. rht < (2.0-rhcrit(irhi,irhj,k) - rhcrit_tol) ) THEN
       rh0 = rhcrit(irhi,irhj,k)
     END IF

    ELSE

    ! Do we need to force some initiation from zero liquid cloud amount?
    ! If so, set RH0 to equal the critical relative humidity

      IF ( (cfl(ind_i(i),ind_j(i),k)  ==  0.0                         &
         .OR. (cfl(ind_i(i),ind_j(i),k)  <   0.05                     &
                  .AND. t(ind_i(i),ind_j(i),k) <   zerodegc))         &
         .AND. (.NOT. cumulus(ind_i(i),ind_j(i)))                     &
         .AND. (r_theta_levels(ind_i(i),ind_j(i),k)                   &
               -r_theta_levels(ind_i(i),ind_j(i),0))                  &
               > zlcl_mixed(ind_i(i),ind_j(i))                        &
         .AND.  rht  >   rhts(ind_i(i),ind_j(i),k)                    &
         .AND. rht  >   (rhcrit(irhi,irhj,k)+rhcrit_tol)  ) THEN
        rh0 = rhcrit(irhi,irhj,k)
      END IF

      ! Do we need to force some initiation from total liquid cloud cover?
      ! If so, set RH0 to equal the critical relative humidity

      IF ( cfl(ind_i(i),ind_j(i),k) == 1.0                            &
      .AND. rht < rhts(ind_i(i),ind_j(i),k)                           &
          .AND. rht < (2.0-rhcrit(irhi,irhj,k)                        &
          - rhcrit_tol) ) THEN
        rh0 = rhcrit(irhi,irhj,k)
      END IF

    END IF

    ! ----------------------------------------------------------------------
    ! 4. Calculate the cloud fraction to initiate to
    ! ----------------------------------------------------------------------

    IF (rh0 > 0.0) THEN
      npti=npti+1

      ! Is the liquid cloud cover coming down from one or up from zero?
      ! If it is coming down then reverse the value of RHT and work with
      ! saturation deficit taking the place of liquid water content.

      IF (cfl(ind_i(i),ind_j(i),k)  >   0.5) THEN
        rht=2.0-rht
      END IF

      ! Calculate the new liquid cloud fraction. Start by calculating QN then
      ! use the Smith scheme relationships to convert this to a cloud fraction

      ! when using the TKE based RHcrit parametrization, force the scheme
      ! to always use a symmetric triangular PDF, otherwise use the original
      ! PC2 method

      qn=(rht-1.0)/(1.0-rh0)
      IF (qn <= -1.0) THEN
        c_1 = 0.0
      ELSE IF (qn < 0.0 .AND. i_rhcpt == rhcpt_tke_based) THEN
        c_1 = 0.5*(1.0+qn)**2
      ELSE IF (qn < 0.0 .AND. i_rhcpt /= rhcpt_tke_based) THEN
        c_1 = 0.5*(1.0+qn)**(pdf_power+1.0)
      ELSE IF (qn < 1.0 .AND. i_rhcpt == rhcpt_tke_based) THEN
        c_1 = 1.0 - 0.5*(1.0-qn)**2
      ELSE IF (qn < 1.0 .AND. i_rhcpt /= rhcpt_tke_based) THEN
        c_1 = 1.0 - 0.5*(1.0-qn)**(pdf_power+1.0)
      ELSE
        c_1 = 1.0
      END IF

      ! Reverse the new cloud fraction back to its correct value if the cloud
      ! fraction is being decreased from a high value.

      IF (cfl(ind_i(i),ind_j(i),k) > 0.5) THEN
        c_1=1.0-c_1
      END IF

      ! Calculate change in total cloud fraction. This depends upon the sign
      ! of the change of liquid cloud fraction. Change in ice cloud fraction
      ! is zero.

      deltacl_c(npti) = c_1 - cfl(ind_i(i),ind_j(i),k)
      deltacf_c(npti) = 0.0

      ! ----------------------------------------------------------------------
      ! 5. Calculate the liquid water content to initiate to.
      ! ----------------------------------------------------------------------

      ! We need to iterate to obtain the liquid content because the
      ! rate of change of gradient of the saturation specific humidity is
      ! specified as a function of the temperature, not the liquid
      ! temperature. We will only iterate over the necessary points where
      ! the initiation works out that a change in cloud fraction is required.

      ! Gather variables

      ni              (npti) = ind_i(i)
      nj              (npti) = ind_j(i)
      t_c             (npti) = t(ind_i(i),ind_j(i),k)
      qn_c            (npti) = qn
      rh0_c           (npti) = rh0
      qsl_tl_c        (npti) = qsl_tl(i)
      cf_c            (npti) = cf(ind_i(i),ind_j(i),k)
      cfl_c           (npti) = cfl(ind_i(i),ind_j(i),k)
      cff_c           (npti) = cff(ind_i(i),ind_j(i),k)
      qcl_c           (npti) = qcl(ind_i(i),ind_j(i),k)
      q_c             (npti) = q(ind_i(i),ind_j(i),k)
      p_theta_levels_c(npti) = p_c(i)

      ! End if for RH0 gt 0
    END IF

  END DO !i

  ! Calculate change in total cloud fraction.

  IF (npti > 0) THEN
    CALL pc2_total_cf(                                                &
          npti,cfl_c,cff_c,deltacl_c,deltacf_c,cf_c)
  END IF

  ! Iterations_do1:
  DO l = 1, init_iterations

    ! Calculate saturated specific humidity with respect to the temperature
    ! (not the liquid temperature).

    IF ( l_new_qsat_lsc ) THEN
      IF ( l_mixing_ratio ) THEN
        CALL qsat_wat_mix_new(qsl_t_c,t_c,p_theta_levels_c,npti)
      ELSE
        CALL qsat_wat_new(qsl_t_c,t_c,p_theta_levels_c,npti)
      END IF
    ELSE
      ! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsl_t_c,t_c,p_theta_levels_c,npti,l_mixing_ratio)
    END IF

    ! Now loop over the required points

    ! Points_do1
    DO i = 1, npti
      alpha=repsilon*lc*qsl_t_c(i)/(r*t_c(i)**2)
      al=1.0/(1.0+lcrcp*alpha)
      bs=al*(1.0-rh0_c(i))*qsl_tl_c(i)

      ! when using the TKE based RHcrit parametrization, force the scheme
      ! to always use a symmetric triangular PDF, otherwise use the original
      ! PC2 method

      IF (qn_c(i) <= -1.0) THEN
        l_bs = 0.0
      ELSE IF (qn_c(i) < 0.0 .AND. i_rhcpt == rhcpt_tke_based) THEN
        l_bs = 0.5/3.0*(1.0+qn_c(i))**3
      ELSE IF (qn_c(i) < 0.0 .AND. i_rhcpt /= rhcpt_tke_based) THEN
        l_bs = 0.5/(pdf_power+2.0)*(1.0+qn_c(i))**(pdf_power+2.0)
      ELSE IF (qn_c(i) < 1.0 .AND. i_rhcpt == rhcpt_tke_based) THEN
        l_bs = (qn_c(i)+0.5/3.0                                       &
               *(1.0-qn_c(i))**3)
      ELSE IF (qn_c(i) < 1.0 .AND. i_rhcpt /= rhcpt_tke_based) THEN
        l_bs = (qn_c(i)+0.5/(pdf_power+2.0)                           &
               *(1.0-qn_c(i))**(pdf_power+2.0))
      ELSE
        l_bs = qn_c(i)
      END IF

      l_out = l_bs * bs

      ! Are we working with saturation deficit instead of liquid water?
      ! If so, convert back to be the correct way around.

      IF (cfl_c(i) > 0.5) THEN
        q_out = qsl_t_c(i) - l_out / al
        l_out = qcl_c(i) + q_c(i) - q_out
      END IF


      ! Calculate amount of condensation. To ensure convergence move
      ! slowly towards the calculated liquid water content L_OUT

      descent_factor = al
      deltal         = descent_factor*(l_out-qcl_c(i))
      q_c(i)         = q_c(i)   - deltal
      qcl_c(i)       = qcl_c(i) + deltal
      t_c(i)         = t_c(i)   + deltal*lcrcp

    END DO ! Points_do1

  END DO ! Iterations_do1

  ! Now scatter back values which have been changed

  DO i = 1, npti

    q  (ni(i),nj(i),k) = q_c(i)
    qcl(ni(i),nj(i),k) = qcl_c(i)
    t  (ni(i),nj(i),k) = t_c(i)
    cf (ni(i),nj(i),k) = cf_c(i)
    cfl(ni(i),nj(i),k) = cfl_c(i) + deltacl_c(i)

  END DO !i

END DO ! Levels_do1:
!$OMP END PARALLEL DO

! End of the subroutine

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pc2_initiate
END MODULE pc2_initiate_mod

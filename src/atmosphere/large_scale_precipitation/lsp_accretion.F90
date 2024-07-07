! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Accretion of droplets by raindrops
! Subroutine Interface:
MODULE lsp_accretion_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_ACCRETION_MOD'

CONTAINS

SUBROUTINE lsp_accretion(                                                     &
  points, timestep,                                                           &
                                          ! Number of points and tstep
  qcl, qrain,                                                                 &
                                          ! Water contents
!    &, area_liq, area_mix, area_ice      ! Partition information
  cfliq,                                                                      &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  rainfrac, rain_liq, rain_mix,                                               &
                                          ! Rain fraction information
!    &, rain_ice, rain_clear
!    &, cf, cfl                           ! Current cloud fractions for
!                                         ! updating
  rho, corr, dhir,                                                            &
                                          ! Parametrization information
  ptransfer,                                                                  &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                               &
                                          ! 1/(timestep*iterations)
!    &, cftransfer,cfltransfer            ! Cloud transfer diagnostics
!    &, rftransfer                        ! Rain fraction transfer
  r_theta_levels_c, fv_cos_theta_latitude_c,                                  &
  f_arr1, f_arr2, f_arr3                                                      &
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
USE lsprec_mod, ONLY: cx, constp, acc_pref, acc_qc, acc_qr, qcfmin,           &
                      f_cons, fsd_eff_lam, fsd_eff_phi, rad_mcica_sigma,      &
                      two_d_fsd_factor, zero, half, one

  ! Microphysics Modules- logicals and integers
USE mphys_constants_mod, ONLY: l_inhomog
USE mphys_inputs_mod,    ONLY: l_rainfall_as, l_warm_new, c_r_correl,         &
                               l_mcr_qrain, l_fsd_generator

! Use in KIND for large scale precip, used for compressed variables passed down
! from here
USE um_types,            ONLY: real_lsprec

! Dr Hook Modules
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

! FSD parameters module- logicals and integers
USE fsd_parameters_mod,  ONLY: ip_fsd_constant
! Constant FSD value and ratio between 1D and 2D FSD
USE rad_input_mod,       ONLY: i_fsd

IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of accretion of cloud
!   droplets by rain drops

! Method:
!   Solve implicitly the microphysical transfer equation for a
!   specified distribution of raindrops sweeping out a volume
!   of air uniformally inhabited by cloud water droplets.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.


! Subroutine Arguments

INTEGER, INTENT(IN) ::                                                        &
  points
                        ! Number of points to calculate

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  timestep,                                                                   &
                        ! Timestep / s
  one_over_tsi,                                                               &
                        ! 1/(timestep*iterations)
  rain_liq(points),                                                           &
                        ! Overlap fraction of rain and liquid
  rain_mix(points),                                                           &
                        ! Overlap fraction of rain and mixed phase
  rainfrac(points),                                                           &
                        ! Rain fraction
  cfliq(points),                                                              &
                        ! Liquid cloud fraction at start of timestep
!    &, cf(points)        ! Current cloud fraction
!    &, cfl(points)       ! Current liquid cloud fraction
    rho(points),                                                              &
                          ! Air density / kg m-3
    corr(points),                                                             &
                          ! Fall velocity correction factor (no units)
    r_theta_levels_c(points),                                                 &
                          ! Distance from centre of Earth and ...
    fv_cos_theta_latitude_c(points),                                          &
                          ! ... grid info for working out gridbox size.
    dhir(points),                                                             &
                          ! Depth of layer / timestep  / m s-1
    f_arr1(points),                                                           &
    f_arr2(points),                                                           &
    f_arr3(points)

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  qcl(points),                                                                &
                        ! Liquid water content / kg kg-1
  qrain(points),                                                              &
                        ! Rain mixing ratio / kg kg-1
!    &, rain_ice(points)  ! Overlap fraction of rain and ice
!    &, rain_clear(points)! Overlap fraction of rain and clear sky
    ptransfer(points) ! Mass rimed in this timestep / kg kg-1
!    &, rftransfer(points)! Transfer of rain fraction this timestep
!    &, cftransfer(points) ! Cumulative cloud fraction increment
!    &, cfltransfer(points)! Cumulative liquid cloud fraction increment

! Local Variables

INTEGER ::                                                                    &
  i                 ! Loop counter for points

REAL (KIND=real_lsprec) ::                                                    &
  dpr(points),                                                                &
                        ! Transfer of mixing ratio  / kg kg-1
  temp1(points),                                                              &
                        ! Temporary in rate calculations
  lambda(points),                                                             &
                        ! Temporary lambda for Abel/Shipway
  lambda_h1(points),                                                          &
                        ! Temporary lambda + h1r for Abel/Shipway
  lambda_h2(points),                                                          &
                        ! Temporary lambda + h2r for Abel/Shipway
  temp7(points),                                                              &
                        ! Rain free liquid cloud fraction (no units)
  qclnew(points)    ! Updated value of liquid water / kg kg-1

REAL (KIND=real_lsprec) ::                                                    &
  fsd_qc(points),                                                             &
          ! fractional standard dev of cloud water content
  fsd_qr(points),                                                             &
          ! fractional standard dev of rain water content
  bias,                                                                       &
          ! accretion rate bias
  x_in_km ! grid-box size in km

! Local compression variable
INTEGER ::                                                            &
  npts,                                                               &
                        ! Number of points to compute
  c,                                                                  &
                        ! Compressed point index 
  ix(points)
                        ! Original index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_ACCRETION'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! When working in 32-bit, a Cray compiler bug on power breaks PROC
! comparability.  Conditionally using NOVECTOR makes this go
! away. This directive will be valid until the end of the program unit
! or when a VECTOR directive (without any extra string) is used.
! Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAY_FORTRAN) && (CRAY_FORTRAN <8004000)
!DIR$ NOVECTOR
#endif

! In all the following loops, the 'i' variable is used for the
! original index while the 'c' variable is used for the compressed
! index.


IF (l_warm_new) THEN

  ! Identify the points that need to be calculated.
  npts = 0
  DO i = 1, points
    
    IF (rainfrac(i) >  zero .AND. qrain(i) >  zero                           &
        .AND. cfliq(i) >  zero) THEN
      
      npts = npts + 1
      ix(npts) = i
      
    END IF ! qcf(i) >  m0 etc.
    
  END DO

  !---------------------------------------------
  ! Use new accretion scheme
  !---------------------------------------------
  !DIR$ IVDEP
  !DIR$ VECTOR ALWAYS
  DO c = 1, npts

    i = ix(c)


    IF (l_inhomog) THEN

      ! calculate bias E factor based on Boutle et al 2012 QJ
      x_in_km = 0.001_real_lsprec*SQRT (r_theta_levels_c(i) * fsd_eff_lam   &
                             * r_theta_levels_c(i) * fsd_eff_phi            &
                             * fv_cos_theta_latitude_c(i)     )

      IF (l_fsd_generator) THEN ! Use same FSD param as cld generator
        IF (i_fsd == ip_fsd_constant) THEN
          fsd_qc(i) = rad_mcica_sigma
        ELSE
          IF ( cfliq(i) < one ) THEN
            fsd_qc(i) = (f_arr2(i)-f_arr3(i)*cfliq(i))                      &
              *(((x_in_km*cfliq(i))**0.333_real_lsprec)                     &
              *((f_cons(1)*x_in_km*cfliq(i))**f_cons(2)+one)                &
             **(f_cons(3)))
          ELSE
            fsd_qc(i) = f_arr1(i) * ((x_in_km**0.333_real_lsprec)           &
              *((f_cons(1)*x_in_km)**f_cons(2)+one)**(f_cons(3)))
          END IF
        END IF
      ELSE ! Use FSD param from Boutle et al 2012 QJ
        IF ( cfliq(i) < one ) THEN
          fsd_qc(i)=(0.45_real_lsprec-0.25_real_lsprec*cfliq(i))            &
                  *(((x_in_km*cfliq(i))**0.333_real_lsprec)                 &
                 *((0.06_real_lsprec*x_in_km*cfliq(i))**1.5_real_lsprec+one)&
                 **(-0.17_real_lsprec))
        ELSE
          fsd_qc(i)=0.11_real_lsprec*(((x_in_km*cfliq(i))**0.333_real_lsprec)&
                   *((0.06_real_lsprec*x_in_km*cfliq(i))                    &
                   **1.5_real_lsprec+one)**(-0.17_real_lsprec))
        END IF
      END IF

      fsd_qr(i) = (1.1_real_lsprec-0.8_real_lsprec*rainfrac(i))             &
               *(((x_in_km*rainfrac(i))**0.333_real_lsprec)                 &
               *((0.11_real_lsprec*x_in_km*rainfrac(i))                     &
               **1.14_real_lsprec+one)**(-0.22_real_lsprec))

      fsd_qc(i) = fsd_qc(i)*two_d_fsd_factor
      fsd_qr(i) = fsd_qr(i)*two_d_fsd_factor

    END IF

  END DO

  !DIR$ IVDEP
  !DIR$ VECTOR ALWAYS
  DO c = 1, npts

    i = ix(c)


    IF (l_inhomog) THEN


      bias = ((one+fsd_qc(i)**2)**(-half*acc_qc))*                          &
             ((one+fsd_qc(i)**2)**(half*acc_qc**2))*                        &
             ((one+fsd_qr(i)**2)**(-half*acc_qr))*                          &
             ((one+fsd_qr(i)**2)**(half*acc_qr**2))*                        &
             EXP(c_r_correl*acc_qc*acc_qr*                                  &
             SQRT(LOG(one+fsd_qc(i)**2)*                                    &
             LOG(one+fsd_qr(i)**2)))

    ELSE ! no inhomog param

      bias = one

    END IF


    dpr(i) = acc_pref * bias *                                              &
             ( (qcl(i)/cfliq(i))**acc_qc) *                                 &
             ( (qrain(i)/rainfrac(i))**acc_qr )                             &
             * timestep * (rain_liq(i)+rain_mix(i))
    dpr(i) = MAX(MIN(dpr(i),                                                &
             qcl(i)*(rain_liq(i)+rain_mix(i))/cfliq(i)-qcfmin),zero)

    ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

    !------------------------------------------------
    ! Adjust liquid and rain contents
    !------------------------------------------------
    qrain(i) = qrain(i) + dpr(i)
    qcl(i)   = qcl(i)   - dpr(i)

  END DO 

ELSE ! original accretion

  ! Identify the points that need to be calculated.
  npts = 0
  DO i = 1, points

    IF (rainfrac(i) >  zero .AND. qrain(i) >  zero                            &
        .AND. cfliq(i) >  zero) THEN
    
      npts = npts + 1
      ix(npts) = i
      
    END IF ! qcf(i) >  m0 etc.
    
  END DO

  !DIR$ IVDEP
  !DIR$ VECTOR ALWAYS
  DO c = 1, npts

    i = ix(c)


    IF (l_mcr_qrain) THEN
            ! Use the prognostic rain formulation
      temp1(i) = qrain(i)*rho(i)*constp(50)/(rainfrac(i))
    ELSE
            ! Use the diagnostic rain formulation
      temp1(i) = qrain(i)*rho(i)*dhir(i)                                    &
               /(rainfrac(i)*constp(42)*corr(i))

    END IF  ! l_mcr_qrain

        !-----------------------------------------------
        ! Calculate the new local value of qcl
        !-----------------------------------------------
    IF (l_rainfall_as) THEN

      IF (l_mcr_qrain) THEN

           !Prognostic Abel and Shipway version

           !First need to calculate lambda:
           ! lambda = (1 / (rain fraction temp1)) to power of cx(52)

        lambda(i) = ( one / (rainfrac(i)*temp1(i)) )**cx(52)

        lambda_h1(i) = lambda(i)+cx(56)
        lambda_h2(i) = lambda(i)+cx(57)

        qclnew(i)= qcl(i) / (                                               &
                   cfliq(i)+(cfliq(i)*timestep*corr(i)*  (                  &
                   (constp(51) * (lambda(i)**cx(46))) /                     &
                    (lambda_h1(i)**cx(61))     +                            &
                   (constp(52) * (lambda(i)**cx(46))) /                     &
                    (lambda_h2(i)**cx(62))     )  )                         &
                           )


      ELSE

           !Diagnostic Abel and Shipway version
           !First need lambda, which for this case is
           !lambda = (1 /(rain fraction temp1)) to power of cx(48)

        lambda(i) = ( one / (rainfrac(i)*temp1(i)) )**cx(48)

        lambda_h1(i) = lambda(i)+cx(56)
        lambda_h2(i) = lambda(i)+cx(57)

           !Now we have lambda, should be the same as above.
           !Keeping prognostic and diagnostic separate for now so
           !I can change one without affecting the other.

        qclnew(i)= qcl(i) / (                                               &
                  cfliq(i)+(cfliq(i)*timestep*corr(i)*  (                   &
                  (constp(51) * (lambda(i)**cx(46))) /                      &
                   (lambda_h1(i)**cx(61))     +                             &
                  (constp(52) * (lambda(i)**cx(46))) /                      &
                   (lambda_h2(i)**cx(62))     )  )                          &
                          )


      END IF ! l_mcr_qrain

    ELSE

           ! Use Sachinananda and Zrnic (1986) rain velocity formula

      IF (l_mcr_qrain) THEN
        qclnew(i)=qcl(i)/((cfliq(i)+cfliq(i)*constp(49)*corr(i)             &
                                 *timestep*temp1(i)**cx(53)))
      ELSE
        qclnew(i)=qcl(i)/((cfliq(i)+cfliq(i)*constp(49)*corr(i)             &
                                 *timestep*temp1(i)**cx(50)))
      END IF ! l_mcr_qrain

    END IF ! L_rainfall_as

        !-----------------------------------------------
        ! Convert qclnew to a gridbox mean
        !-----------------------------------------------
        ! temp7 is the rain free region of liquid cloud
    temp7(i) = MAX(cfliq(i)-rain_liq(i)-rain_mix(i),zero)
    qclnew(i) = qclnew(i)*(rain_liq(i)+rain_mix(i))                         &
               +qcl(i)/cfliq(i)*temp7(i)

        !-----------------------------------------------
        ! Calculate transfer
        !-----------------------------------------------
    dpr(i) = qcl(i) - qclnew(i)

    ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

      !------------------------------------------------
      ! Adjust liquid and rain contents
      !------------------------------------------------
    qrain(i) = qrain(i) + dpr(i)
    qcl(i)   = qcl(i)   - dpr(i)

      !------------------------------------------------
      ! Update cloud fractions
      !------------------------------------------------
    !         These are commented out since there is currently no
    !         cloud fraction update associated with the accretion term.
    !          cf_transfer_rate(i)  = 0.0 / (timestep*iterations)
    !          cfl_transfer_rate(i) = 0.0 / (timestep*iterations)


    !            cftransfer(i)  = cftransfer(i)  + cf_transfer_rate(i)
    !            cfltransfer(i) = cfltransfer(i) + cfl_transfer_rate(i)

    !
    !            cf(i)  = cf(i)  +cf_transfer_rate(i) *timestep*iterations
    !            cfl(i) = cfl(i) +cfl_transfer_rate(i)*timestep*iterations
    !

              !-----------------------------------------------
              ! Update rain fractions
              !-----------------------------------------------
    !         These are commented out since there is currently no
    !         rain fraction update associated with the accretion term.
    !          rf_final(i) = rainfrac(i)
    !          rftransfer(i) = rftransfer(i) + (rf_final(i) - rainfrac(i))

    !
    !          rainfrac(i)= rf_final(i)
    !          rain_liq(i)= min(area_liq(i),rainfrac(i))
    !          rain_mix(i)= min(area_mix(i),rainfrac(i)-rain_liq(i))
    !          rain_ice(i)=
    !     &         min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
    !          rain_clear(i)=
    !     &         rainfrac(i) - rain_liq(i)-rain_mix(i)-rain_ice(i)
    !

  END DO  

END IF ! l_warm_new

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_accretion
END MODULE lsp_accretion_mod

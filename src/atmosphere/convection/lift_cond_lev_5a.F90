! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Lifting condensation level calculation
!
MODULE lift_cond_lev_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'LIFT_COND_LEV_MOD'
CONTAINS

SUBROUTINE lift_cond_lev (npnts, nlev, k_plume,                          &
                          pstar, q, t,                                   &
                          p_theta_lev, exner_rho, z_rho,                 &
                          T_lcl, p_lcl, z_lcl, qsat_lcl )

USE planet_constants_mod, ONLY: kappa, pref, repsilon, recip_kappa

USE cv_diag_param_mod, ONLY:                                             &
    a_bolton, b_bolton, c_bolton, d_bolton

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    qsat_mix_new     => qsat_mix,                       &
                    l_new_qsat_conv !Currently defaults to FALSE

USE gen_phys_inputs_mod, ONLY: l_mr_physics

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   This routine calculates the lifting condensation level (LCL) temperature,
!   pressure and height.
!
!  Is designed to work on compressed arrays for just a selected set of points.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                & ! Number of points
 ,nlev                   ! Number of model levels for calculations

INTEGER, INTENT(IN) :: &
  k_plume(npnts)         ! Starting model level for plume ascent

REAL, INTENT(IN) ::       &
  pstar(npnts)            & ! Surface pressure (Pa)
 ,q(npnts,nlev)           & ! water vapour on model levels (kg/kg)
 ,t(npnts,nlev)           & ! Temperature on model levels (K)
 ,p_theta_lev(npnts,nlev) & ! Pressure on theta levels (Pa)
 ,exner_rho(npnts,nlev)   & ! Exner Pressure on rho levels
 ,z_rho(npnts,nlev)         ! Hieght of rho levels  (m)

REAL, INTENT(OUT) ::      &
  T_lcl(npnts)            & ! Temperature of LCL  (K)
 ,p_lcl(npnts)            & ! Pressure of LCL  (Pa)
 ,z_lcl(npnts)            & ! Height of LCL  (m)
 ,qsat_lcl(npnts)           ! qsaturation at zlcl (i.e. cloud base) kg/kg

!-------------------------------------------------------------------------------
! Local variables

INTEGER ::               &
  i,k                      ! loop counter

REAL ::                  &
  exner_lcl              & ! Exner pressure at LCL
 ,exner_surf               ! Exner pressure at surface

REAL ::                  &
  factor                 & ! factor used in interpolation
 ,vap_press                ! vapour pressure

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LIFT_COND_LEV'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! Calculate temperature and pressure of lifting condensation level
!     using approximations from Bolton (1980)
!-------------------------------------------------------------------------------
!
!   vapour pressure e ~ qp/epsilon       q specific humidity
!   vapour pressure e ~ qp/(epsilon+q)   q mixing ratio
!-------------------------------------------------------------------------------

IF (l_mr_physics) THEN   ! expression for mixing ratio

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(i, vap_press)                      &
!$OMP& SHARED(npnts,k_plume,pstar,q,t,p_theta_lev,t_lcl,p_lcl,              &
!$OMP&        recip_kappa,repsilon)                                         & 
!$OMP& SCHEDULE(STATIC)  
  DO i=1, npnts

    vap_press = 0.01*q(i,k_plume(i)) * p_theta_lev(i,k_plume(i))            &
                                      / (repsilon+q(i,k_plume(i)) )
    IF (vap_press  >   0.0) THEN
      T_lcl(i) = a_bolton + b_bolton/                                       &
                              (c_bolton*LOG(t(i,k_plume(i)))                &
                                         - LOG(vap_press) - d_bolton )

      p_lcl(i) = p_theta_lev(i,k_plume(i)) *                                &
                     ( T_lcl(i) / t(i,k_plume(i)) )**recip_kappa

    ELSE
      ! If no moisture present, set LCL to the top of the atmosphere
      ! so that the diagnosis parcel ascent will never reach it
      T_lcl(i) = 0.0
      p_lcl(i) = 0.0
    END IF

  END DO     ! i loop
!$OMP END PARALLEL DO

ELSE          ! expression for specific humidity

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(i, vap_press)                      &
!$OMP& SHARED(npnts,k_plume,pstar,q,t,p_theta_lev,t_lcl,p_lcl,              &
!$OMP&        recip_kappa,repsilon)                                         &
!$OMP& SCHEDULE(STATIC)  
  DO i=1, npnts
    vap_press = q(i,k_plume(i)) *                                           &
                       p_theta_lev(i,k_plume(i)) / ( 100.0*repsilon )
    IF (vap_press  >   0.0) THEN
      T_lcl(i) = a_bolton + b_bolton/                                       &
                         (c_bolton*LOG(t(i,k_plume(i)))                     &
                                         - LOG(vap_press) - d_bolton )
      p_lcl(i) = p_theta_lev(i,k_plume(i)) *                                &
                        ( T_lcl(i) / t(i,k_plume(i)) )**recip_kappa
    ELSE
      ! If no moisture present, set LCL to the top of the atmosphere
      ! so that the diagnosis parcel ascent will never reach it
      T_lcl(i) = 0.0
      p_lcl(i) = 0.0
    END IF

  END DO
!$OMP END PARALLEL DO

END IF ! test on l_mr_physics

! work out qsat at LCL
IF ( l_new_qsat_conv ) THEN
  IF ( l_mr_physics ) THEN
    CALL qsat_mix_new(qsat_lcl,T_lcl,p_lcl,npnts)
  ELSE
    CALL qsat_new(qsat_lcl,T_lcl,p_lcl,npnts)
  END IF
ELSE
  ! DEPENDS ON: qsat_mix
  CALL qsat_mix(qsat_lcl,T_lcl,p_lcl,npnts,l_mr_physics)
END IF

!-----------------------------------------------------------------------
! Accurate calculation of height of LCL using exner_lcl rather than p_lcl
! as UM interpolates exner linearly in height but not pressure.
!-----------------------------------------------------------------------

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(i, k, factor, exner_lcl, exner_surf) &
!$OMP& SHARED(npnts,nlev,exner_rho,z_rho,z_lcl,p_lcl,pstar,pref,kappa)        &
!$OMP& SCHEDULE(STATIC)  
! The main index of the loop is "i" so that it can be
! parallelised using OpenMP. 
DO i=1, npnts
  exner_lcl  = (p_lcl(i)/pref)**kappa
  exner_surf = (pstar(i)/pref)**kappa

  k = 1

  IF ( exner_lcl >= exner_surf) THEN
    z_lcl(i) = 0.0           ! at or below surface
  ELSE IF (exner_lcl < exner_surf                                    &
                   .AND. exner_lcl > exner_rho(i,k)) THEN
    factor= (exner_rho(i,k) - exner_lcl)/                            &
                      (exner_rho(i,k) - exner_surf)
    z_lcl(i) = (1.0-factor)*z_rho(i,k)
  END IF

  DO k=2,nlev
    IF (exner_lcl >= exner_rho(i,k)                                    &
                        .AND. exner_lcl < exner_rho(i,k-1) ) THEN
      factor= (exner_rho(i,k) - exner_lcl)/                            &
                        (exner_rho(i,k) - exner_rho(i,k-1))
      z_lcl(i) = (1.0-factor)*z_rho(i,k)+factor*z_rho(i,k-1)
    END IF
  END DO         ! level loop

  ! If LCL pressure is lower than the pressure at the top of the
  ! model, set z_lcl to the model-top
  IF ( exner_lcl < exner_rho(i,nlev) ) THEN
    z_lcl(i) = z_rho(i,nlev)
  END IF

! Check z_lcl not less than a minimum value
! Note z_lcl currently used by diurnal cycle diagnosis and
! also by Deep turbulence scheme. The diurnal cycle diagnosis will be
! happy with a value of zero but this will not work for the deep
! turbulence scheme.
! Set z_lcl lowest model depth - fix may imply some model resolution dependence.

  z_lcl(i) =MAX(z_lcl(i),z_rho(i,2))
END DO
!$OMP END PARALLEL DO

!-------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lift_cond_lev
END MODULE lift_cond_lev_mod

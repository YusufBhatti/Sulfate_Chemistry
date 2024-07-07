! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!   New Shallow convection scheme - Calculation of turbulent fluxes
!
MODULE shconv_turb_fluxes_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SHCONV_TURB_FLUXES_MOD'
CONTAINS

SUBROUTINE shconv_turb_fluxes(n_sh,ntra, max_cldlev, nlev         &
,                      max_cldtrlev, trlev                        &
,                      imethod_precip ,l_tracer                   &
,                      kterm_thv,kterm_tracer                     &
,                      fng_hfunc, gql_func                        &
,                      wql_cld2, mf_cld, q_cld,qsat_cld           &
,                      dqsatdt_cld,theta_cld, exner_thcld         &
,                      rho_cld, precip_product_cld,wup_h_cld      &
,                      wthetav_surf,zcld,wstar_up,wtracer_cb      &
,                      wtheta, wq, wthetav                        &
,                      wthetal, wqt, wtracer)

!
! Purpose:
!   Calculate turbulent fluxes w'theta' and w'q' on in cloud levels
!
!
!   Called by SHALLOW_CONV5A
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
USE water_constants_mod, ONLY: lc
USE planet_constants_mod, ONLY: cp, c_virtual, rv

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
INTEGER, INTENT(IN) ::                                            &
  n_sh                                                            &
                ! No. of shallow convection points
, ntra                                                            &
                ! No. of tracers
, max_cldlev                                                      &
                ! Maximum number of convective cloud levels
, nlev                                                            &
                ! Maximum number of convective cloud levels
, max_cldtrlev                                                    &
                ! Maximum number of convective cloud levels
, trlev                                                           &
                ! Maximum number of tracer levels
, imethod_precip  ! swicth for precip production 1 for on

LOGICAL, INTENT(IN) ::                                            &
  l_tracer      ! true - tracers present

REAL, INTENT(IN) ::                                               &
  kterm_thv(n_sh,nlev)                                            &
                          ! thetav gradient (kterm)
, fng_hfunc(n_sh,nlev)                                            &
                          ! Fng on rho levels
, gql_func(n_sh,nlev)                                             &
                          ! gql on rho levels
, wql_cld2(n_sh,nlev)                                             &
                          ! wql on uv cloud lev (similarity func)
, mf_cld(n_sh,nlev)                                               &
                          ! mass_flux on theta levels in cloud
, q_cld(n_sh,nlev)                                                &
                          ! q on theta levels in cloud
, qsat_cld(n_sh,nlev)                                             &
                          ! qsat on theta levels in cloud
, dqsatdt_cld(n_sh,nlev)                                          &
                          !dqsat/dT on theta levels in cloud
, theta_cld(n_sh,nlev)                                            &
                          ! theta on cloud levels
, exner_thcld(n_sh,nlev)                                          &
                          ! exner on th cloud levels
, rho_cld(n_sh,nlev)                                              &
                          ! density on rho levels
, precip_product_cld(n_sh,nlev)                                   &
                                ! precip production rho levels
, wup_h_cld (n_sh,nlev)                                           &
                          ! wup on rho levels
, kterm_tracer(n_sh,nlev,ntra)                                    &
                                ! tracer gradient (kterm)
!                                   dimensioned on nlev not trlev
      , wthetav_surf(n_sh)                                              &
                                ! wthetav flux at surface
      , zcld(n_sh)                                                      &
                                ! cloud depth (m)
      , wstar_up(n_sh)                                                  &
                                ! in cloud convective velocity scale
      , wtracer_cb(n_sh,ntra)   ! w'tracer' at cloud base (kg/kg m/s)
!
! Arguments with intent INOUT:
!
!     None
!
! Arguments with intent OUT:
!
REAL, INTENT(OUT) ::                                              &
 wtheta(n_sh,nlev)                                                &
                    ! wtheta flux on cloud levels  (K m/s)
,wq(n_sh,nlev)                                                    &
                    ! wq flux on cloud levels      (kg/kg m/s)
,wthetav(n_sh,nlev)                                               &
                    ! wthetav flux on cloud levels  (K m/s)
,wthetal(n_sh,nlev)                                               &
                    ! wthetal flux on cloud levels  (K m/s)
,wqt(n_sh,nlev)                                                   &
                    ! wqt flux on cloud levels      (kg/kg m/s)
,wtracer(n_sh,trlev,ntra) ! wtracer flux on cloud levels(kg/kg m/s)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

INTEGER :: i,k,ktra          ! loop counters

REAL ::                                                           &
  wtheta_vl_cb(n_sh)         ! wtheta_vl flux at cloud base.

! temporary stores
REAL ::                                                           &
  fngterm2,  A_term                                               &
, div_term
REAL ::                                                           &
 lc_o_cp      ! Lc/Cp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SHCONV_TURB_FLUXES'

!-----------------------------------------------------------------------
! 1.0 Initialise arrays and variables
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Constants evaluated once at start of routine to save CPU

lc_o_cp = lc/cp

!-----------------------------------------------------------------------
! 1.0 Initialise arrays - set all functions to zero
!-----------------------------------------------------------------------
! Initialise dhdz
DO k = 1,max_cldlev
  DO i = 1,n_sh
    wq(i,k) = 0.0
    wqt(i,k) = 0.0
    wtheta(i,k) = 0.0
    wthetav(i,k) = 0.0
    wthetal(i,k) = 0.0
  END DO
END DO


!-----------------------------------------------------------------------
! 2.0 Level loop to calculate in cloud values of fluxes on uv levels
!-----------------------------------------------------------------------
!  Numbering of levels
!  level 1 is cloud base for uv - in cloud levels
!
!-----------------------------------------------------------------------
DO i = 1,n_sh

  wtheta_vl_cb(i) = -0.16*wthetav_surf(i)

END DO


DO k = 1,max_cldlev
  DO i = 1,n_sh

    ! thetav flux
    fngterm2 = wtheta_vl_cb(i)*fng_hfunc(i,k)

    IF (imethod_precip == 1) THEN
      ! additional term for precip added to original w'thetav'

      wthetav(i,k) = -kterm_thv(i,k) + fngterm2                    &
            + (lc_o_cp/exner_thcld(i,k)-c_virtual*theta_cld(i,k))  &
     *zcld(i)*(wup_h_cld(i,k)/wstar_up(i))*precip_product_cld(i,k) &
                     /rho_cld(i,k)
    ELSE
      wthetav(i,k) = -kterm_thv(i,k) + fngterm2
    END IF


    ! calculations depend on interpolation for q and qsat
    ! preferred method
    IF (gql_func(i,k) /= 0.0) THEN
      A_term = (1.0+lc_o_cp*dqsatdt_cld(i,k))*wql_cld2(i,k)      &
      /gql_func(i,k) - mf_cld(i,k)*(q_cld(i,k) - qsat_cld(i,k))
    ELSE  ! not in cloud
      a_term=0.0
    END IF

    ! Should I use mf on th or uv levels?

    div_term=1.0+c_virtual*exner_thcld(i,k)*theta_cld(i,k)      &
                                            *dqsatdt_cld(i,k)

    wthetal(i,k)=(wthetav(i,k)-c_virtual*theta_cld(i,k)*a_term)&
   /div_term
    wqt(i,k) =( A_term + wthetav(i,k)*exner_thcld(i,k)         &
                   *dqsatdt_cld(i,k))/div_term


  END DO
END DO

! ----------------------------------------------------------------------
! 3.0 Tracers
! ----------------------------------------------------------------------

IF (l_tracer) THEN

  ! initialise
  DO ktra = 1,ntra
    DO k = 1,max_cldtrlev
      DO i = 1,n_sh
        wtracer(i,k,ktra) = 0.0
      END DO
    END DO
  END DO

  DO ktra = 1,ntra
    DO k = 1,max_cldtrlev
      DO i = 1,n_sh

        ! non-gradient term
        fngterm2 = wtracer_cb(i,ktra)*fng_hfunc(i,k)

        wtracer(i,k,ktra) = -kterm_tracer(i,k,ktra) + fngterm2

      END DO
    END DO
  END DO
END IF

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE shconv_turb_fluxes
END MODULE shconv_turb_fluxes_mod

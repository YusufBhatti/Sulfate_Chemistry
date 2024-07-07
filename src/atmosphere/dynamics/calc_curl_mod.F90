! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

MODULE calc_curl_mod

USE model_domain_mod,      ONLY: mt_global, model_type
USE um_parparams,          ONLY: pnorth, psouth
USE um_parvars,            ONLY: at_extremity

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_CURL_MOD'

CONTAINS
SUBROUTINE calc_curl_vort(u,v,w,curl_u,curl_v,curl_w)

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod, ONLY: udims, vdims, wdims,                          &
                                 udims_s, vdims_s, wdims_s,                    &
                                 pdims, pdims_s

USE horiz_grid_mod,        ONLY: xi1_u, xi2_v, xi1_p, xi2_p,                   &
                                 intw_rho2w, intw_p2u, intw_p2v,               &
                                 Csxi2_p, Csxi2_v

USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels,             &
                                 xi3_at_theta => r_theta_levels,               &
                                 xi3_at_rho   => r_rho_levels,                 &
                                 xi3_at_u     => r_at_u,                       &
                                 xi3_at_v     => r_at_v

USE metric_terms_mod,      ONLY: deta_xi3_u, deta_xi3_v, deta_xi3,             &
                                 deta_xi3_theta,                               &
                                 h1_xi1_u, h2_xi2_v, h3_p_eta

USE planet_constants_mod,   ONLY: planet_radius

USE curl_at_poles_mod,     ONLY: curl_at_poles

IMPLICIT NONE
!
! Description:
!   Calculates curl of a wind field.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ONLY coded for spherical geometry and no orography !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Method: ENDGame formulation version 1.01,
!         section 7.2.
  
! Subroutine arguments

! Input
REAL, INTENT(IN) :: u(udims_s%i_start:udims_s%i_end,                           &
                      udims_s%j_start:udims_s%j_end,                           &
                      udims_s%k_start:udims_s%k_end)

REAL, INTENT(IN) :: v(vdims_s%i_start:vdims_s%i_end,                           &
                      vdims_s%j_start:vdims_s%j_end,                           &
                      vdims_s%k_start:vdims_s%k_end)

REAL, INTENT(IN) :: w(wdims_s%i_start:wdims_s%i_end,                           &
                      wdims_s%j_start:wdims_s%j_end,                           &
                      wdims_s%k_start:wdims_s%k_end)

! Output
REAL, INTENT(OUT) :: curl_u(pdims_s%i_start:pdims_s%i_end,                     &
                            vdims_s%j_start:vdims_s%j_end,                     &
                            wdims_s%k_start:wdims_s%k_end)

REAL, INTENT(OUT) :: curl_v(udims_s%i_start:udims_s%i_end,                     &
                            pdims_s%j_start:pdims_s%j_end,                     &
                            wdims_s%k_start:wdims_s%k_end)

REAL, INTENT(OUT) :: curl_w(udims_s%i_start:udims_s%i_end,                     &
                            vdims_s%j_start:vdims_s%j_end,                     &
                            pdims_s%k_start:pdims_s%k_end)

! Local
INTEGER :: i, j, k                             ! loop counters
INTEGER :: j_start, j_end, pole, j_pole
REAL    :: height_domain, recip_height_domain
REAL    :: rdxi1_u(udims_s%i_start:udims_s%i_end),                             &
           rdxi1_p(pdims_s%i_start:pdims_s%i_end),                             &
           rdxi2_v(vdims_s%j_start:vdims_s%j_end),                             &
           rdxi2_p(pdims_s%j_start:pdims_s%j_end),                             &
           rdz_w(wdims%k_start:wdims%k_end),                                   &
           rdz_rho(pdims%k_start:pdims%k_end)

REAL    :: curl_xi3_pole(udims%k_start:udims%k_end)

REAL    :: work(wdims_s%i_start:wdims_s%i_end,                                 &
                wdims_s%j_start:wdims_s%j_end,                                 &
                wdims_s%k_start:wdims_s%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_CURL_VORT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Grid spacing
DO i=pdims_s%i_start, pdims_s%i_end
  rdxi1_p(i) = 1.0/(xi1_u(i)-xi1_u(i-1))
END DO

DO i=udims_s%i_start, udims_s%i_end
  rdxi1_u(i) = 1.0/(xi1_p(i+1)-xi1_p(i))
END DO

DO j=vdims_s%j_start, vdims_s%j_end
  rdxi2_v(j) = 1.0/(xi2_p(j+1)-xi2_p(j))
END DO

DO j=pdims_s%j_start, pdims_s%j_end
  rdxi2_p(j) = 1.0/(xi2_v(j)-xi2_v(j-1))
END DO

height_domain=xi3_at_theta(1,1,wdims%k_end)-planet_radius
recip_height_domain=1.0/height_domain

rdz_w(0) = recip_height_domain/(eta_rho_levels(1) - eta_theta_levels(0))
DO k=1, wdims%k_end-1
  rdz_w(k) = recip_height_domain/(eta_rho_levels(k+1) - eta_rho_levels(k))
END DO
rdz_w(wdims%k_end) = recip_height_domain/(eta_theta_levels(wdims%k_end) -      &
                            eta_rho_levels(wdims%k_end))

DO k=pdims%k_start,pdims%k_end
  rdz_rho(k) = recip_height_domain/(eta_theta_levels(k)-eta_theta_levels(k-1))
END DO


! ----------------------------------------------------------------------------
! U-component
! ----------------------------------------------------------------------------

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k)                                                          &
!$OMP& SHARED(wdims,vdims,pdims,curl_u,w,rdxi2_v,xi3_at_rho,v,rdz_w,           &
!$OMP&        xi3_at_theta)
DO k=wdims%k_start+1, wdims%k_end-1
  DO j=vdims%j_start, vdims%j_end
    DO i=pdims%i_start, pdims%i_end
      curl_u(i,j,k)=((w(i,j+1,k)-w(i,j,k))*rdxi2_v(j) -                        &
                     (xi3_at_rho(i,j,k+1)*v(i,j,k+1) -                         &
                      xi3_at_rho(i,j,k)*v(i,j,k))*rdz_w(k))/xi3_at_theta(i,j,k)

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

DO j=vdims_s%j_start, vdims_s%j_end
  DO i=pdims_s%i_start, pdims_s%i_end
    curl_u(i,j,wdims%k_start)=0.0
    curl_u(i,j,wdims%k_end)=0.0
  END DO
END DO

! ----------------------------------------------------------------------------
! V-component
! ----------------------------------------------------------------------------
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k)                                                          &
!$OMP& SHARED(wdims,pdims,udims,curl_v,xi3_at_rho,u,rdz_w,w,rdxi1_u,           &
!$OMP&        Csxi2_p,xi3_at_theta)
DO k=wdims%k_start+1, wdims%k_end-1
  DO j=pdims%j_start, pdims%j_end
    DO i=udims%i_start, udims%i_end
      curl_v(i,j,k)=((xi3_at_rho(i,j,k+1)*u(i,j,k+1) -                         &
                      xi3_at_rho(i,j,k)*u(i,j,k))*rdz_w(k) -                   &
                     (w(i+1,j,k)-w(i,j,k))*rdxi1_u(i)/Csxi2_p(j)) /            &
                     xi3_at_theta(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

DO j=pdims_s%j_start, pdims_s%j_end
  DO i=udims_s%i_start, udims_s%i_end
    curl_v(i,j,wdims%k_start)=0.0
    curl_v(i,j,wdims%k_end)=0.0
  END DO
END DO

! ----------------------------------------------------------------------------
! W-component
! ----------------------------------------------------------------------------
pole=0
j_start=vdims%j_start
j_end=vdims%j_end
IF (model_type == mt_global) THEN
  IF (at_extremity(psouth)) THEN
    pole=-1
    j_pole=j_start
    j_start=j_pole+1
    CALL curl_at_poles(pole,u,curl_xi3_pole)
    DO k=udims%k_start,udims%k_end
      DO i=udims%i_start, udims%i_end
        curl_w(i,j_pole,k)=curl_xi3_pole(k)
      END DO
    END DO
  END IF
  IF (at_extremity(pnorth)) THEN
    pole=1
    j_pole=j_end
    j_end=j_pole-1
    CALL curl_at_poles(pole,u,curl_xi3_pole)
    DO k=udims%k_start,udims%k_end
      DO i=udims%i_start, udims%i_end
        curl_w(i,j_pole,k)=curl_xi3_pole(k)
      END DO
    END DO
  END IF
END IF
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k)                                                          &
!$OMP& SHARED(pdims,j_start,j_end,udims,curl_w,v,rdxi1_u,Csxi2_p,u,rdxi2_v,    &
!$OMP&        xi3_at_rho,Csxi2_v)
DO k=pdims%k_start, pdims%k_end
  DO j=j_start, j_end
    DO i=udims%i_start, udims%i_end
      curl_w(i,j,k)=((v(i+1,j,k)-v(i,j,k))*rdxi1_u(i) -                        &
                     (Csxi2_p(j+1)*u(i,j+1,k)-Csxi2_p(j)*u(i,j,k))*rdxi2_v(j))/&
                    (xi3_at_rho(i,j,k)*Csxi2_v(j))
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! Need calls to swap bounds for vorticity
! The first and third components also require polar values

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_curl_vort

! =============================================================================
SUBROUTINE calc_curl_velo(u,v,w,curl_u,curl_v,curl_w)

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod, ONLY: udims, vdims, wdims,                          &
                                 udims_s, vdims_s, wdims_s,                    &
                                 pdims, pdims_s

USE horiz_grid_mod,        ONLY: xi1_u, xi2_v, xi1_p, xi2_p,                   &
                                 intw_rho2w, intw_p2u, intw_p2v,               &
                                 Csxi2_p, Csxi2_v

USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels,             &
                                 xi3_at_theta => r_theta_levels,               &
                                 xi3_at_rho   => r_rho_levels,                 &
                                 xi3_at_u     => r_at_u,                       &
                                 xi3_at_v     => r_at_v

USE metric_terms_mod,      ONLY: deta_xi3_u, deta_xi3_v, deta_xi3,             &
                                 deta_xi3_theta,                               &
                                 h1_xi1_u, h2_xi2_v, h3_p_eta

USE planet_constants_mod,   ONLY: planet_radius

IMPLICIT NONE
!
! Description:
!   Calculates curl of a wind field.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ONLY coded for spherical geometry and no orography !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Method: ENDGame formulation version 1.01,
!         section 7.2.

! Subroutine arguments

! Input
REAL, INTENT(IN) :: u(pdims_s%i_start:pdims_s%i_end,                           &
                       vdims_s%j_start:vdims_s%j_end,                          &
                       wdims_s%k_start:wdims_s%k_end)

REAL, INTENT(IN) :: v(udims_s%i_start:udims_s%i_end,                           &
                       pdims_s%j_start:pdims_s%j_end,                          &
                       wdims_s%k_start:wdims_s%k_end)

REAL, INTENT(IN) :: w(udims_s%i_start:udims_s%i_end,                           &
                       vdims_s%j_start:vdims_s%j_end,                          &
                       pdims_s%k_start:pdims_s%k_end)


! Output
REAL, INTENT(OUT) :: curl_u(udims%i_start:udims%i_end,                         &
                            udims%j_start:udims%j_end,                         &
                            udims%k_start:udims%k_end)

REAL, INTENT(OUT) :: curl_v(vdims%i_start:vdims%i_end,                         &
                            vdims%j_start:vdims%j_end,                         &
                            vdims%k_start:vdims%k_end)

REAL, INTENT(OUT) :: curl_w(wdims%i_start:wdims%i_end,                         &
                            wdims%j_start:wdims%j_end,                         &
                            wdims%k_start:wdims%k_end)



! Local
INTEGER :: i, j, k, j_start, j_end, pole, j_pole
REAL    :: height_domain, recip_height_domain
REAL    :: rdxi1_u(udims%i_start:udims%i_end),                                 &
           rdxi1_p(pdims%i_start:pdims%i_end),                                 &
           rdxi2_v(vdims%j_start:vdims%j_end),                                 &
           rdxi2_p(pdims%j_start:pdims%j_end),                                 &
           rdz_w(wdims%k_start:wdims%k_end),                                   &
           rdz_rho(pdims%k_start:pdims%k_end)

REAL    :: work(wdims_s%i_start:wdims_s%i_end,                                 &
                wdims_s%j_start:wdims_s%j_end,                                 &
                wdims_s%k_start:wdims_s%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_CURL_VELO'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Grid spacing
DO i=pdims%i_start, pdims%i_end
  rdxi1_p(i) = 1.0/(xi1_u(i)-xi1_u(i-1))
END DO

DO i=udims%i_start, udims%i_end
  rdxi1_u(i) = 1.0/(xi1_p(i+1)-xi1_p(i))
END DO

DO j=vdims%j_start, vdims%j_end
  rdxi2_v(j) = 1.0/(xi2_p(j+1)-xi2_p(j))
END DO

DO j=pdims%j_start, pdims%j_end
  rdxi2_p(j) = 1.0/(xi2_v(j)-xi2_v(j-1))
END DO


height_domain=xi3_at_theta(1,1,wdims%k_end)-planet_radius
recip_height_domain=1.0/height_domain

rdz_w(0) = recip_height_domain/(eta_rho_levels(1) - eta_theta_levels(0))
DO k=1, wdims%k_end-1
  rdz_w(k) = recip_height_domain/(eta_rho_levels(k+1) - eta_rho_levels(k))
END DO
rdz_w(wdims%k_end) = recip_height_domain/(eta_theta_levels(wdims%k_end) -      &
                            eta_rho_levels(wdims%k_end))

DO k=pdims%k_start,pdims%k_end
  rdz_rho(k) = recip_height_domain/(eta_theta_levels(k)-eta_theta_levels(k-1))
END DO

! ----------------------------------------------------------------------------
! U-component
! ----------------------------------------------------------------------------

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k)                                                          &
!$OMP& SHARED(udims,curl_u,w,rdxi2_p,xi3_at_theta,v,rdz_rho,xi3_at_rho)
DO k=udims%k_start, udims%k_end
  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end
      curl_u(i,j,k)=((w(i,j,k)-w(i,j-1,k))*rdxi2_p(j) -                        &
                     (xi3_at_theta(i,j,k)*v(i,j,k) -                           &
                      xi3_at_theta(i,j,k-1)*v(i,j,k-1))*rdz_rho(k)) /          &
                    xi3_at_rho(i,j,k)

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------------
! V-component
! ----------------------------------------------------------------------------
j_start=vdims%j_start
j_end=vdims%j_end
IF (model_type == mt_global) THEN
  IF (at_extremity(psouth)) THEN
    j_pole=j_start
    j_start=j_pole+1
    DO k=vdims%k_start, vdims%k_end
      DO i=vdims%i_start, vdims%i_end
        curl_v(i,j_pole,k)=0.0
      END DO
    END DO
  END IF
  IF (at_extremity(pnorth)) THEN
    j_pole=j_end
    j_end=j_pole-1
    DO k=vdims%k_start, vdims%k_end
      DO i=vdims%i_start, vdims%i_end
        curl_v(i,j_pole,k)=0.0
      END DO
    END DO
  END IF
END IF

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k)                                                          &
!$OMP& SHARED(vdims,j_start,j_end,curl_v,xi3_at_theta,u,rdz_rho,w,rdxi1_p,     &
!$OMP&        Csxi2_v,xi3_at_rho)
DO k=vdims%k_start, vdims%k_end
  DO j=j_start, j_end
    DO i=vdims%i_start, vdims%i_end
      curl_v(i,j,k)=((xi3_at_theta(i,j,k)*u(i,j,k) -                           &
                      xi3_at_theta(i,j,k-1)*u(i,j,k-1))*rdz_rho(k) -           &
                     (w(i,j,k)-w(i-1,j,k))*rdxi1_p(i)/Csxi2_v(j)) /            &
                     xi3_at_rho(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------------
! W-component
! ----------------------------------------------------------------------------

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k)                                                          &
!$OMP& SHARED(wdims,curl_w,v,rdxi1_p,Csxi2_v,u,rdxi2_p,xi3_at_theta,Csxi2_p)
DO k=wdims%k_start+1, wdims%k_end-1
  DO j=wdims%j_start, wdims%j_end
    DO i=wdims%i_start, wdims%i_end
      curl_w(i,j,k)=((v(i,j,k)-v(i-1,j,k))*rdxi1_p(i) -                        &
                     (Csxi2_v(j)*u(i,j,k)-Csxi2_v(j-1)*u(i,j-1,k))*rdxi2_p(j))/&
                    (xi3_at_theta(i,j,k)*Csxi2_p(j))
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

DO j=wdims%j_start, wdims%j_end
  DO i=wdims%i_start, wdims%i_end
    curl_w(i,j,wdims%k_start)=0.0
    curl_w(i,j,wdims%k_end)=0.0
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_curl_velo

END MODULE calc_curl_mod

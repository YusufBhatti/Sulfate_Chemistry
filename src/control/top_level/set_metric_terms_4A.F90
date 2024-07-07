! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
MODULE set_metric_terms_4A_mod
IMPLICIT NONE

REAL,    PARAMETER  :: tny=1.0e-20

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SET_METRIC_TERMS_4A_MOD'

CONTAINS


! Subroutine to determine phi iteratively for  use with oblate
! spherioldal coords (ENDGame Documentation Appendix B)
!
FUNCTION eg_phi(xi2,xi3,a,eps)
USE ereport_mod, ONLY: ereport
IMPLICIT NONE
REAL,    PARAMETER  :: tol    = 1.0e-10
INTEGER, PARAMETER  :: it_max = 100
REAL                :: eg_phi
REAL,    INTENT(IN) :: xi2, xi3, a, eps

REAL                :: err, phi, cxi2, sxi2, z
REAL                :: tmp, ep2
INTEGER             :: it, err_code

phi = ABS(SIN(xi2))
cxi2 = COS(xi2)**2
sxi2 = SIN(xi2)**2
z    = xi3/a
ep2  = eps**2

IF ( ep2 > 1.0e-8 ) THEN
  DO it = 1, it_max
    tmp = phi*z*SQRT( (1.0-ep2)/(1.0-ep2*phi**2) )
    tmp = tmp**(2.0*ep2)

    tmp = SQRT( 1.0 - cxi2/( cxi2 + sxi2*tmp/(1.0-ep2) ) )
    err = ABS(tmp-phi)
    phi = tmp
    IF ( err < tol ) EXIT
  END DO

  IF ( err > tol ) THEN
    err_code = 1
    CALL Ereport('eg_set_metric_terms',err_code,                      &
    'Covergence Failure calculating phi in setup')
  END IF
  eg_phi = ASIN(phi)*SIGN(1.0,xi2)
ELSE
  eg_phi = xi2
END IF

END FUNCTION eg_phi

#if !defined(RECON)

SUBROUTINE set_metric_terms_4A(l_shallow)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE planet_constants_mod, ONLY: planet_radius
USE dynamics_input_mod,   ONLY: eccentricity
USE level_heights_mod,    ONLY: xi3_at_theta=>r_theta_levels,         &
                                xi3_at_rho=>r_rho_levels,             &
                                xi3_at_u=>r_at_u,                     &
                                xi3_at_v=>r_at_v,                     &
                                xi3_at_u_w=>r_at_u_w,                 &
                                xi3_at_v_w=>r_at_v_w,                 &
                                eta_rho_levels, eta_theta_levels
USE horiz_grid_mod

USE UM_ParVars
USE UM_ParParams, ONLY: pnorth, psouth

USE ereport_mod, ONLY: ereport

USE atm_fields_bounds_mod

USE metric_terms_mod

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE
!
! Description: calculates metric terms for an SOS/Cartesian coordinate
!              system and trigonometric functions.
!
!
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

LOGICAL, INTENT(IN) :: l_shallow

! Local variables

INTEGER :: i, j, k, model_levels, err_code

REAL :: gi, gj, a, phi, eps, tmp
REAL :: rdxi1, rdxi2, rdeta

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_METRIC_TERMS_4A'

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL init_metric_terms()

!------------------------------------------------------------------------!
!  Compute deta(xi3) quantities
!------------------------------------------------------------------------!

model_levels = pdims%k_end

DO k=1, model_levels
  rdeta = 1.0/( eta_theta_levels(k) - eta_theta_levels(k-1) )
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      deta_xi3(i,j,k) = ( xi3_at_theta(i,j,k) -                   &
                                    xi3_at_theta(i,j,k-1))        &
                       *rdeta
    END DO
  END DO
END DO

k = 0
rdeta = 1.0/(eta_rho_levels(k+1)-eta_theta_levels(k))
DO j=pdims%j_start, pdims%j_end
  DO i=pdims%i_start, pdims%i_end
    deta_xi3_theta(i,j,k) =                                     &
             (xi3_at_rho(i,j,k+1)-xi3_at_theta(i,j,k))*rdeta
  END DO
END DO

DO k=1, model_levels-1
  rdeta = 1.0/( eta_rho_levels(k+1) - eta_rho_levels(k) )
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      deta_xi3_theta(i,j,k) = ( xi3_at_rho(i,j,k+1) -             &
                                       xi3_at_rho(i,j,k) )*rdeta
    END DO
  END DO
END DO

k = model_levels
rdeta = 1.0/(eta_theta_levels(k)-eta_rho_levels(k))
DO j=pdims%j_start, pdims%j_end
  DO i=pdims%i_start, pdims%i_end
    deta_xi3_theta(i,j,k) =                                     &
             (xi3_at_theta(i,j,k)-xi3_at_rho(i,j,k))*rdeta
  END DO
END DO

! Interpolate full level height field on u-points

DO k=1, model_levels
  rdeta = 1.0/(eta_theta_levels(k)-eta_theta_levels(k-1))
  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end
      deta_xi3_u(i,j,k) = (xi3_at_u_w(i,j,k) -                    &
                                  xi3_at_u_w(i,j,k-1))*rdeta
    END DO
  END DO
END DO

! Interpolate full level height field on v-points

DO k=1, model_levels
  rdeta = 1.0/(eta_theta_levels(k)-eta_theta_levels(k-1))
  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end
      deta_xi3_v(i,j,k) = (xi3_at_v_w(i,j,k) -                    &
                                  xi3_at_v_w(i,j,k-1))*rdeta
    END DO
  END DO
END DO

DO k=0, model_levels
  !CDIR NOUNROLL
  DO j=pdims%j_start, pdims%j_end
    rdxi2 = 1.0/(xi2_v(j)-xi2_v(j-1))
    DO i=pdims%i_start, pdims%i_end
      rdxi1           = 1.0/(xi1_u(i)-xi1_u(i-1))
      dxi1_xi3(i,j,k) = (xi3_at_u_w(i,j,k)-xi3_at_u_w(i-1,j,k))   &
                        *rdxi1
      dxi2_xi3(i,j,k) = (xi3_at_v_w(i,j,k)-xi3_at_v_w(i,j-1,k))   &
                        *rdxi2
    END DO
  END DO
END DO


eps = eccentricity                 ! Hardwired for spherical coords
a   = planet_radius

IF ( eps /= 0.0 ) THEN
  err_code = 1
  CALL Ereport('eg_set_metric_terms',err_code,                   &
  'Code will not work correctly for nonzero Eccentricity')
END IF

IF ( cartesian_grid ) THEN

  ! set metric terms.
  k=0
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      h1_p_eta(i,j,k) = 1.0
      h2_p_eta(i,j,k) = 1.0
      h3_p_eta(i,j,k) = 1.0
    END DO
  END DO

  DO k= pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        h1_p(i,j,k)     = 1.0
        h2_p(i,j,k)     = 1.0
        h3_p(i,j,k)     = 1.0
        h1_p_eta(i,j,k) = 1.0
        h2_p_eta(i,j,k) = 1.0
        h3_p_eta(i,j,k) = 1.0
      END DO
    END DO
  END DO

  DO k=udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        h1_xi1_u(i,j,k) = 1.0
        h2_xi1_u(i,j,k) = 1.0
        h3_xi1_u(i,j,k) = 1.0
      END DO
    END DO
  END DO

  DO k=vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        h1_xi2_v(i,j,k) = 1.0
        h2_xi2_v(i,j,k) = 1.0
        h3_xi2_v(i,j,k) = 1.0
      END DO
    END DO
  END DO

  DO i = udims%i_start, udims%i_end
    csxi1_u(i) = 1.0
    snxi1_u(i) = 0.0
  END DO

  DO i = pdims%i_start, pdims%i_end
    csxi1_p(i) = 1.0
    snxi1_p(i) = 0.0
  END DO

  DO j = vdims_l%j_start, vdims_l%j_end
    csxi2_v(j) = 1.0
    snxi2_v(j) = 0.0
  END DO

  DO j = pdims_l%j_start, pdims_l%j_end
    csxi2_p(j) = 1.0
    snxi2_p(j) = 0.0
  END DO

ELSE  ! .not. cartesian_grid
  ! Precalculate sines and cosines

  DO i = udims%i_start, udims%i_end
    csxi1_u(i) = COS(xi1_u(i))
    snxi1_u(i) = SIN(xi1_u(i))
  END DO

  DO i = pdims%i_start, pdims%i_end
    csxi1_p(i) = COS(xi1_p(i))
    snxi1_p(i) = SIN(xi1_p(i))
  END DO

  DO j = vdims_l%j_start, vdims_l%j_end
    csxi2_v(j) = COS(xi2_v(j))
    snxi2_v(j) = SIN(xi2_v(j))
  END DO

  DO j = pdims_l%j_start, pdims_l%j_end
    csxi2_p(j) = COS(xi2_p(j))
    snxi2_p(j) = SIN(xi2_p(j))
  END DO

  ! set metric terms.
  IF ( l_shallow ) THEN
    k=0
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        h1_p_eta(i,j,k) = planet_radius * csxi2_p(j)
        h2_p_eta(i,j,k) = planet_radius
        h3_p_eta(i,j,k) = 1.0
      END DO
    END DO

    DO k= pdims%k_start, pdims%k_end
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          h1_p(i,j,k)     = planet_radius * csxi2_p(j)
          h2_p(i,j,k)     = planet_radius
          h3_p(i,j,k)     = 1.0

          h1_p_eta(i,j,k) = planet_radius * csxi2_p(j)
          h2_p_eta(i,j,k) = planet_radius
          h3_p_eta(i,j,k) = 1.0
        END DO
      END DO
    END DO

    DO k= udims%k_start, udims%k_end
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          h1_xi1_u(i,j,k) = planet_radius * csxi2_p(j)
          h2_xi1_u(i,j,k) = planet_radius
          h3_xi1_u(i,j,k) = 1.0
        END DO
      END DO
    END DO

    DO k= vdims%k_start, vdims%k_end
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          h1_xi2_v(i,j,k) = planet_radius * csxi2_v(j)
          h2_xi2_v(i,j,k) = planet_radius
          h3_xi2_v(i,j,k) = 1.0
        END DO
      END DO
    END DO

  ELSE   ! not L_shallow
    k=0
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        phi = eg_phi(xi2_p(j),xi3_at_theta(i,j,k),a,eps)
        tmp = SQRT( 1.0 - (eps*SIN(phi))**2 )

        h1_p_eta(i,j,k) = xi3_at_theta(i,j,k)*COS(phi)/tmp
        h2_p_eta(i,j,k) = xi3_at_theta(i,j,k)*SIN(2.0*phi+tny)          &
                          /(tmp*SIN(2.0*xi2_p(j)+tny))
        h3_p_eta(i,j,k) = tmp
        phi_at_eta(i,j,k) = phi
      END DO
    END DO

    DO k= pdims%k_start, pdims%k_end
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          phi = eg_phi(xi2_p(j),xi3_at_rho(i,j,k),a,eps)
          tmp = SQRT( 1.0 - (eps*SIN(phi))**2 )

          h1_p(i,j,k)     = xi3_at_rho(i,j,k)*COS(phi)/tmp
          h2_p(i,j,k)     = xi3_at_rho(i,j,k)*SIN(2.0*phi+tny)            &
                            /(tmp*SIN(2.0*xi2_p(j)+tny))
          h3_p(i,j,k)     = tmp
          phi_at_p(i,j,k) = phi

          phi = eg_phi(xi2_p(j),xi3_at_theta(i,j,k),a,eps)
          tmp = SQRT( 1.0 - (eps*SIN(phi))**2 )

          h1_p_eta(i,j,k) = xi3_at_theta(i,j,k)*COS(phi)/tmp
          h2_p_eta(i,j,k) = xi3_at_theta(i,j,k)*SIN(2.0*phi+tny)          &
                            /(tmp*SIN(2.0*xi2_p(j)+tny))
          h3_p_eta(i,j,k) = tmp
          phi_at_eta(i,j,k) = phi
        END DO
      END DO
    END DO

    DO k= udims%k_start, udims%k_end
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          phi = eg_phi(xi2_p(j),xi3_at_u(i,j,k),a,eps)
          tmp = SQRT( 1.0 - (eps*SIN(phi))**2 )

          h1_xi1_u(i,j,k) = xi3_at_u(i,j,k)*COS(phi)/tmp
          h2_xi1_u(i,j,k) = xi3_at_u(i,j,k)*SIN(2.0*phi+tny)              &
                            /(tmp*SIN(2.0*xi2_p(j)+tny))
          h3_xi1_u(i,j,k) = tmp
          phi_at_u(i,j,k) = phi
        END DO
      END DO
    END DO

    DO k= vdims%k_start, vdims%k_end
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          phi = eg_phi(xi2_v(j),xi3_at_v(i,j,k),a,eps)
          tmp = SQRT( 1.0 - (eps*SIN(phi))**2 )

          h1_xi2_v(i,j,k) = xi3_at_v(i,j,k)*COS(phi)/tmp
          h2_xi2_v(i,j,k) = xi3_at_v(i,j,k)*SIN(2.0*phi+tny)              &
                            /(tmp*SIN(2.0*xi2_v(j)+tny))
          h3_xi2_v(i,j,k) = tmp
          phi_at_v(i,j,k) = phi
        END DO
      END DO
    END DO
  END IF  ! End of If Shallow

  IF (model_type == mt_global) THEN
    IF ( at_extremity(PSouth) ) h1_xi2_v(:,vdims%j_start,:) = 0.0
    IF ( at_extremity(PNorth) ) h1_xi2_v(:,vdims%j_end,:) = 0.0
  END IF

END IF   ! End of If Cartesian

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE set_metric_terms_4A

#endif

END MODULE set_metric_terms_4A_mod

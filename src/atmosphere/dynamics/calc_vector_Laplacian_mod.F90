! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE calc_vector_Laplacian_mod

IMPLICIT NONE

  ! Description: 
  !  Module contain routines/calls to implement explicit viscosity
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: DYNAMICS
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_VECTOR_LAPLACIAN_MOD'

CONTAINS
SUBROUTINE calc_vector_Laplacian(u, v, w, etadot, visco_u, visco_v, visco_w)

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod, ONLY: udims, vdims, wdims,                          &
                                 udims_s, vdims_s, wdims_s,                    &
                                 pdims, pdims_s, array_dims
USE Field_Types,           ONLY: fld_type_u, fld_type_v, fld_type_p
USE dynamics_input_mod,    ONLY: l_viscosity, horiz_viscosity,                 &
                                 vert_viscosity
USE calc_div_mod,          ONLY: calc_div
USE calc_grad_mod,         ONLY: calc_grad
USE calc_curl_mod,         ONLY: calc_curl_vort, calc_curl_velo
USE eg_v_at_poles_mod,     ONLY: eg_v_at_poles
USE um_parvars,            ONLY: at_extremity
USE um_parparams,          ONLY: pnorth, psouth
USE model_domain_mod,      ONLY: mt_global, model_type
USE halo_exchange,         ONLY: swap_bounds
USE mpp_conf_mod,          ONLY: swap_field_is_vector, swap_field_is_scalar

IMPLICIT NONE
!
! Description:
!   Calculates divergence of a vector field, whose components are
!   distributed on a C-grid precisely as the wind field.
!
! Method: ENDGame formulation version 1.01,
!         section 7.2.
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

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

REAL, INTENT(IN) :: etadot(wdims_s%i_start:wdims_s%i_end,                      &
                           wdims_s%j_start:wdims_s%j_end,                      &
                           wdims_s%k_start:wdims_s%k_end)

! Output: Note that limits match those of ustar, vstar, wstar in
!         SISL_Init_uvw
REAL, INTENT(OUT) :: visco_u(udims%i_start:udims%i_end,                        &
                             udims%j_start:udims%j_end,                        &
                             udims%k_start:udims%k_end)

REAL, INTENT(OUT) :: visco_v(vdims%i_start:vdims%i_end,                        &
                             vdims%j_start:vdims%j_end,                        &
                             vdims%k_start:vdims%k_end)

REAL, INTENT(OUT) :: visco_w(wdims%i_start:wdims%i_end,                        &
                             wdims%j_start:wdims%j_end,                        &
                             wdims%k_start:wdims%k_end)

INTEGER :: i, j, k

REAL :: div(pdims_s%i_start:pdims_s%i_end,                                     &
            pdims_s%j_start:pdims_s%j_end,                                     &
            pdims_s%k_start:pdims_s%k_end)

REAL :: div_surf(pdims_s%i_start:pdims_s%i_end,                                &
                 pdims_s%j_start:pdims_s%j_end)

REAL :: grad_u(udims%i_start:udims%i_end,                                      &
               udims%j_start:udims%j_end,                                      &
               udims%k_start:udims%k_end)

REAL :: grad_v(vdims%i_start:vdims%i_end,                                      &
               vdims%j_start:vdims%j_end,                                      &
               vdims%k_start:vdims%k_end)

REAL :: grad_w(wdims%i_start:wdims%i_end,                                      &
               wdims%j_start:wdims%j_end,                                      &
               wdims%k_start:wdims%k_end)

REAL :: curl_u(pdims_s%i_start:pdims_s%i_end,                                  &
               vdims_s%j_start:vdims_s%j_end,                                  &
               wdims_s%k_start:wdims_s%k_end)

REAL :: curl_v(udims_s%i_start:udims_s%i_end,                                  &
               pdims_s%j_start:pdims_s%j_end,                                  &
               wdims_s%k_start:wdims_s%k_end)

REAL :: curl_w(udims_s%i_start:udims_s%i_end,                                  &
               vdims_s%j_start:vdims_s%j_end,                                  &
               pdims_s%k_start:pdims_s%k_end)

TYPE(array_dims) :: dims

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_VECTOR_LAPLACIAN'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

visco_u=0.0
visco_v=0.0
visco_w=0.0

CALL calc_div(u,v,etadot,pdims_s,div)
div_surf(:,:)=div(:,:,1)

CALL calc_grad(div,div_surf,grad_u,grad_v,grad_w)


CALL calc_curl_vort(u,v,w,curl_u,curl_v,curl_w)

CALL swap_bounds(curl_u, pdims%i_len, vdims%j_len, wdims%k_len,                &
                 pdims_s%halo_i, vdims_s%halo_j,                               &
                 fld_type_v, swap_field_is_vector)

CALL swap_bounds(curl_v, udims%i_len, pdims%j_len, wdims%k_len,                &
                 udims_s%halo_i, pdims_s%halo_j,                               &
                 fld_type_u, swap_field_is_vector)

CALL swap_bounds(curl_w, udims%i_len, vdims%j_len, pdims%k_len,                &
                 udims_s%halo_i, vdims_s%halo_j,                               &
                 fld_type_p, swap_field_is_scalar)

CALL calc_curl_velo(curl_u,curl_v,curl_w,visco_u,visco_v,visco_w)

DO k=udims%k_start, udims%k_end
  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end
      visco_u(i,j,k)=grad_u(i,j,k)-visco_u(i,j,k)
    END DO
  END DO
END DO

DO k=vdims%k_start, vdims%k_end
  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end
      visco_v(i,j,k)=grad_v(i,j,k)-visco_v(i,j,k)
    END DO
  END DO
END DO



DO k=wdims%k_start, wdims%k_end
  DO j=wdims%j_start, wdims%j_end
    DO i=wdims%i_start, wdims%i_end
      visco_w(i,j,k)=grad_w(i,j,k)-visco_w(i,j,k)
    END DO
  END DO
END DO

IF (model_type == mt_global) THEN
  IF (at_extremity(psouth)) THEN
    CALL eg_v_at_poles(visco_u,visco_v,1.0,                                    &
                       udims%j_start,vdims%j_start,udims,vdims)
  END IF
  IF (at_extremity(pnorth)) THEN
    CALL eg_v_at_poles(visco_u,visco_v,-1.0,                                   &
                       udims%j_end,vdims%j_end,udims,vdims)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_vector_Laplacian
END MODULE calc_vector_Laplacian_mod

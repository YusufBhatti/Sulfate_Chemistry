! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Subroutine Interface:
MODULE eg_set_adv_winds_mod
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SET_ADV_WINDS_MOD'

CONTAINS
SUBROUTINE eg_set_adv_winds(u,v,w,u_adv,v_adv,w_adv,                  &
                           row_length,rows,n_rows,model_levels,       &
                           halo_i, halo_j, l_shallow)


USE metric_terms_mod,    ONLY: h1_xi1_u, h2_xi2_v
USE atm_fields_bounds_mod
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar,swap_field_is_vector
USE planet_constants_mod, ONLY: planet_radius
USE level_heights_mod,   ONLY: r_at_u,r_at_v
USE field_types
USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

LOGICAL, INTENT(IN) :: l_shallow

INTEGER, INTENT(IN) :: row_length,rows,n_rows,model_levels,           &
                       halo_i, halo_j

REAL, INTENT(IN) :: u(udims_s%i_start:udims_s%i_end,                  &
                      udims_s%j_start:udims_s%j_end,                  &
                      udims_s%k_start:udims_s%k_end),                 &
                    v(vdims_s%i_start:vdims_s%i_end,                  &
                      vdims_s%j_start:vdims_s%j_end,                  &
                      vdims_s%k_start:vdims_s%k_end),                 &
                    w(wdims_s%i_start:wdims_s%i_end,                  &
                      wdims_s%j_start:wdims_s%j_end,                  &
                      wdims_s%k_start:wdims_s%k_end)

REAL, INTENT(INOUT) :: u_adv(udims_l%i_start:udims_l%i_end,           &
                             udims_l%j_start:udims_l%j_end,           &
                             udims_l%k_start:udims_l%k_end),          &
                       v_adv(vdims_l%i_start:vdims_l%i_end,           &
                             vdims_l%j_start:vdims_l%j_end,           &
                             vdims_l%k_start:vdims_l%k_end),          &
                       w_adv(wdims_l%i_start:wdims_l%i_end,           &
                             wdims_l%j_start:wdims_l%j_end,           &
                             wdims_l%k_start:wdims_l%k_end)

INTEGER :: i,j,k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SET_ADV_WINDS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                             &
!$OMP          SHARED(model_type, udims, vdims, model_levels, u_adv,    &
!$OMP                 u, h1_xi1_u, v_adv, v, h2_xi2_v, planet_radius,   &
!$OMP                 r_at_u, r_at_v, w_adv, w, pdims, l_shallow)
IF (model_type /= mt_global) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        u_adv(i,j,k) = u(i,j,k)/h1_xi1_u(i,j,k)
      END DO
    END DO

    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        v_adv(i,j,k) = v(i,j,k)/h2_xi2_v(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
ELSE
  IF ( l_shallow ) THEN

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          u_adv(i,j,k) = u(i,j,k)/planet_radius
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          v_adv(i,j,k) = v(i,j,k)/planet_radius
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

  ELSE

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          u_adv(i,j,k) = u(i,j,k)/r_at_u(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          v_adv(i,j,k) = v(i,j,k)/r_at_v(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

  END IF
END IF

!$OMP DO SCHEDULE(STATIC)
DO k = 0, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      w_adv(i,j,k) = w(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

CALL swap_bounds(u_adv,                                                &
                 udims_l%i_len - 2*udims_l%halo_i,                     &
                 udims_l%j_len - 2*udims_l%halo_j,                     &
                 udims_l%k_len,                                        &
                 udims_l%halo_i, udims_l%halo_j,                       &
                 fld_type_u,swap_field_is_vector)
CALL swap_bounds(v_adv,                                                &
                 vdims_l%i_len - 2*vdims_l%halo_i,                     &
                 vdims_l%j_len - 2*vdims_l%halo_j,                     &
                 vdims_l%k_len,                                        &
                 vdims_l%halo_i, vdims_l%halo_j,                       &
                 fld_type_v,swap_field_is_vector)
CALL swap_bounds(w_adv,                                                &
                 wdims_l%i_len - 2*wdims_l%halo_i,                     &
                 wdims_l%j_len - 2*wdims_l%halo_j,                     &
                 wdims_l%k_len,                                        &
                 wdims_l%halo_i, wdims_l%halo_j,                       &
                 fld_type_p,swap_field_is_scalar)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_set_adv_winds
END MODULE eg_set_adv_winds_mod


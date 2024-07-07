! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE eg_calc_ax_mod


! Description: This subroutine multiplies the input field
!              by the Constant part of the Helmholtz operator
!
! Method: ENDGame formulation version 3.02
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

PRIVATE

PUBLIC :: eg_Calc_Ax
INTERFACE eg_Calc_Ax
MODULE PROCEDURE eg_Calc_Ax_dp_dp,                    &
     eg_Calc_Ax_sp_dp,                                &
     eg_Calc_Ax_dp_sp,                                &
     eg_Calc_Ax_sp_sp
END INTERFACE


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_CALC_AX_MOD'

CONTAINS

SUBROUTINE eg_Calc_Ax_dp_dp(Ax,x,dims,dims_s)

USE atm_fields_bounds_mod
USE Field_Types
USE helmholtz_const_matrix_mod
USE um_types, ONLY: real64
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

USE eg_helmholtz_mod, ONLY: HM_p, HM_theta, HM_pp
USE ref_pro_mod
USE horiz_grid_mod, ONLY: intw_w2rho
USE atm_fields_bounds_mod

USE pressure_grad_mod
USE div_pu_mod
USE dynamics_input_mod, ONLY: l_inc_solver

USE um_parvars,          ONLY: at_extremity
USE um_ParParams,        ONLY: Pnorth, Psouth, Pwest, Peast
USE model_domain_mod,    ONLY: model_type, mt_lam

IMPLICIT NONE

TYPE (array_dims) dims, dims_s

! x is the input vector and Ax is the result of multiplying
! x by the matrix A

REAL(KIND=real64), INTENT(INOUT) :: Ax(dims_s%i_start:dims_s%i_end,      &
                                       dims_s%j_start:dims_s%j_end,      &
                                       dims_s%k_start:dims_s%k_end),     &
                                     x(dims_s%i_start:dims_s%i_end,      &
                                       dims_s%j_start:dims_s%j_end,      &
                                       dims_s%k_start:dims_s%k_end)

REAL ::    u(udims_s%i_start:udims_s%i_end,                              &
             udims_s%j_start:udims_s%j_end,                              &
             udims_s%k_start:udims_s%k_end)

REAL ::    v(vdims_s%i_start:vdims_s%i_end,                              &
             vdims_s%j_start:vdims_s%j_end,                              &
             vdims_s%k_start:vdims_s%k_end)

REAL ::    w(wdims%i_start:wdims%i_end,                                  &
             wdims%j_start:wdims%j_end,                                  &
             wdims%k_start:wdims%k_end)

REAL :: etad(wdims%i_start:wdims%i_end,                                  &
             wdims%j_start:wdims%j_end,                                  &
             wdims%k_start:wdims%k_end)

REAL :: hm_tp

INTEGER                       :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CALC_AX_DP_DP'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( l_inc_solver ) THEN

  CALL swap_bounds(x, dims%i_len, dims%j_len, dims%k_len,                &
                   dims_s%halo_i, dims_s%halo_j,                         &
                   fld_type_p, swap_field_is_scalar,                     &
                   do_corners_arg=.FALSE.)

! Calculate pressure gradient as implied increments to u,v,w (etadot)
  CALL pressure_grad(x, u, v, w, etad)

! Take divergence of result
  CALL div_Pu(Ax, u,v,etad)

! Now add remaining terms to get full Helmholtz result

  DO k = dims%k_start,dims%k_end
    DO j = dims%j_start, dims%j_end
      DO i = dims%i_start, dims%i_end
        hm_tp     = rho_ref_pro(i,j,k)/thetav_ref_eta(i,j,k)
        hm_tp     = hm_tp*( intw_w2rho(k,1)*HM_theta(i,j,k)*w(i,j,k)       &
                          + intw_w2rho(k,2)*HM_theta(i,j,k-1)*w(i,j,k-1) )

        Ax(i,j,k) = (HM_pp*rho_ref_pro(i,j,k)                              &
                     /exner_ref_pro(i,j,k))*x(i,j,k) - Ax(i,j,k) - hm_tp

        Ax(i,j,k) = Ax(i,j,k)*Hlm_Lp(i,j,k)
 
      END DO
    END DO
  END DO

  IF ( model_type == mt_lam ) THEN
    i = dims%i_end - 1         ! one less p-point in LAM's
    IF ( at_extremity(PWest) ) THEN
      Ax(1,:,:) = x(1,:,:)
    END IF
    IF ( at_extremity(PEast) ) THEN
      Ax(i,:,:)   = x(i,:,:)
      Ax(i+1,:,:) = x(i+1,:,:)
    END IF

    j = dims%j_end
    IF ( at_extremity(PSouth) ) THEN
      Ax(:,1,:) = x(:,1,:)
    END IF
    IF ( at_extremity(PNorth) ) THEN
     Ax(:,j,:) = x(:,j,:)
    END IF
  END IF

ELSE

#include "eg_calc_ax.h"

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_calc_ax_dp_dp

!!!!!!!!!!!!!!!! END of subroutine double double   !!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE eg_Calc_Ax_sp_dp(Ax,x,dims,dims_s)

USE atm_fields_bounds_mod
USE Field_Types
USE helmholtz_const_matrix_mod
USE um_types, ONLY: real32, real64
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

IMPLICIT NONE

TYPE (array_dims) dims, dims_s

! x is the input vector and Ax is the result of multiplying
! x by the matrix A

REAL (KIND=real32),     INTENT(INOUT) ::                               &
     Ax(dims_s%i_start:dims_s%i_end,                                   &
        dims_s%j_start:dims_s%j_end,                                   &
        dims_s%k_start:dims_s%k_end)
REAL (KIND=real64), INTENT(INOUT) ::                                   &
      x(dims_s%i_start:dims_s%i_end,                                   &
        dims_s%j_start:dims_s%j_end,                                   &
        dims_s%k_start:dims_s%k_end)


INTEGER                       :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CALC_AX_SP_DP'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#include "eg_calc_ax.h"

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_calc_ax_sp_dp

!!!!!!!!!!!!!!!!!!!!!!!!! End single double !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE eg_Calc_Ax_dp_sp(Ax,x,dims,dims_s)

USE atm_fields_bounds_mod
USE Field_Types
USE helmholtz_const_matrix_mod
USE um_types, ONLY: real32, real64
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
IMPLICIT NONE

TYPE (array_dims) dims, dims_s

! x is the input vector and Ax is the result of multiplying
! x by the matrix A

REAL(KIND=real64),     INTENT(INOUT) ::                               &
     Ax(dims_s%i_start:dims_s%i_end,                                  &
        dims_s%j_start:dims_s%j_end,                                  &
        dims_s%k_start:dims_s%k_end)
REAL (KIND=real32), INTENT(INOUT) ::                                  &
     x(dims_s%i_start:dims_s%i_end,                                   &
        dims_s%j_start:dims_s%j_end,                                  &
        dims_s%k_start:dims_s%k_end)


INTEGER                       :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CALC_AX_DP_SP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#include "eg_calc_ax.h"

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_calc_ax_dp_sp

!!!!!!!!!!!!!!!!!!! End of Double Single !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE eg_Calc_Ax_sp_sp(Ax,x,dims,dims_s)

USE atm_fields_bounds_mod
USE Field_Types
USE helmholtz_const_matrix_mod
USE um_types, ONLY: real32
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
IMPLICIT NONE

TYPE (array_dims) dims, dims_s

! x is the input vector and Ax is the result of multiplying
! x by the matrix A

REAL (KIND=real32),     INTENT(INOUT) ::                              &
     Ax(dims_s%i_start:dims_s%i_end,                                  &
        dims_s%j_start:dims_s%j_end,                                  &
        dims_s%k_start:dims_s%k_end),                                 &
      x(dims_s%i_start:dims_s%i_end,                                  &
        dims_s%j_start:dims_s%j_end,                                  &
        dims_s%k_start:dims_s%k_end)


INTEGER                       :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CALC_AX_SP_SP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#include "eg_calc_ax.h"

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_calc_ax_sp_sp


END MODULE eg_calc_ax_mod

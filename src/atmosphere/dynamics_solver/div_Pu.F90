! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE div_Pu_mod
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIV_PU_MOD'

CONTAINS
SUBROUTINE div_Pu(div,u, v, w)
USE atm_fields_bounds_mod
USE metric_terms_mod
USE eg_helmholtz_mod
USE horiz_grid_mod
USE ref_pro_mod
USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels,           &
                             xi3_at_theta=>r_theta_levels

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
!
! Description: Calculate the divergence of the input
!              vector (u,v,w) (NOTE w=etadot)
!
!
! Method: Appendix F, ENDGame formulation version 1.01
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

REAL,  INTENT(OUT) :: div(pdims_s%i_start:pdims_s%i_end,                       &
                          pdims_s%j_start:pdims_s%j_end,                       &
                          pdims_s%k_start:pdims_s%k_end)

REAL,  INTENT(IN)  :: u(udims_s%i_start:udims_s%i_end,                         &
                        udims_s%j_start:udims_s%j_end,                         &
                        udims_s%k_start:udims_s%k_end)

REAL,  INTENT(IN)  :: v(vdims_s%i_start:vdims_s%i_end,                         &
                        vdims_s%j_start:vdims_s%j_end,                         &
                        vdims_s%k_start:vdims_s%k_end)

REAL,  INTENT(IN)  :: w(wdims%i_start:wdims%i_end,                             &
                        wdims%j_start:wdims%j_end,                             &
                        wdims%k_start:wdims%k_end)

REAL               :: r_dxi1, r_dxi2, r_deta

INTEGER            :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DIV_PU'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,                &
!$OMP                    r_deta,r_dxi2,r_dxi1)                                 &
!$OMP&            SHARED(pdims, div, u,v,w,                                    &
!$OMP&                   eta_theta_levels,xi2_v,xi1_u,                         &
!$OMP&                   hm_rhox, hm_rhoy, hm_rhoz, hm_vol)
DO k = pdims%k_start, pdims%k_end
   r_deta = 1.0/(eta_theta_levels(k)-eta_theta_levels(k-1))
   DO j = pdims%j_start, pdims%j_end
      r_dxi2 = 1.0/(xi2_v(j) - xi2_v(j-1))
      DO i = pdims%i_start, pdims%i_end
         r_dxi1     = 1.0/(xi1_u(i) - xi1_u(i-1))
         div(i,j,k) = ( ( HM_rhox(i,j,k)*u(i,j,k)                              &
                         -HM_rhox(i-1,j,k)*u(i-1,j,k) )*r_dxi1                 &
                       +( HM_rhoy(i,j,k)*v(i,j,k)                              &
                         -HM_rhoy(i,j-1,k)*v(i,j-1,k) )*r_dxi2                 &
                       +( HM_rhoz(i,j,k)*w(i,j,k) -                            &
                          HM_rhoz(i,j,k-1)*w(i,j,k-1) )*r_deta                 &
                      )*HM_vol(i,j,k)
      END DO
   END DO
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE div_Pu
END MODULE div_Pu_mod

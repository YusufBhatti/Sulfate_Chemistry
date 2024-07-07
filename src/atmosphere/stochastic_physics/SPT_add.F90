! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:  SPT_add gets increments from SPT and adds them in both
!             ENDGAME cycles for fast physics.
!
!     Code Owner: Please refer to the UM file CodeOwners.txt
!     This file belongs in section: Stochastic Physics
MODULE SPT_add_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SPT_ADD_MOD'

CONTAINS

SUBROUTINE SPT_add( theta_star,q_star)

! Bounds of arrays
USE atm_fields_bounds_mod, ONLY:                                  &
    tdims, tdims_s,                                               &
    pdims, pdims_s,                                               &
    udims, udims_s,                                               &
    vdims, vdims_s
! Load r_u variables and spt increment arrays
USE fields_rhs_mod,       ONLY:                                   &
    r_u_p2,r_v_p2,theta_spt,q_spt,r_u_spt,r_v_spt

! DrHook modules
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER :: i,j,k
      ! Integer for DO loops

REAL ::                                                                 &
     theta_star(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                tdims%k_start:tdims%k_end)                              &
      ! theta after atmos_physics2


 ,   q_star(tdims%i_start:tdims%i_end,                              &
              tdims%j_start:tdims%j_end,                            &
              tdims%k_start:tdims%k_end)
      ! q after atmos_physics2

! DrHook stuff
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SPT_ADD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


DO k=tdims%k_start,tdims%k_end
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      theta_star(i,j,k)= theta_star(i,j,k) + theta_spt(i,j,k)
    END DO
  END DO
END DO

DO k=tdims%k_start,tdims%k_end
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      q_star(i,j,k)= q_star(i,j,k) + q_spt(i,j,k)
    END DO
  END DO
END DO

DO k=udims%k_start,udims%k_end
  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end
      r_u_p2(i,j,k) = r_u_p2(i,j,k) + r_u_spt(i,j,k)
    END DO
  END DO
END DO

DO k=vdims%k_start,vdims%k_end
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end
      r_v_p2(i,j,k) = r_v_p2(i,j,k) + r_v_spt(i,j,k)
    END DO
  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE SPT_add

END MODULE SPT_add_mod

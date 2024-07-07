! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE save_conv_diags_mod

IMPLICIT NONE
!---------------------------------------------------------------------------
! Description:
!   Saves the convective increments etc (which are defined locally in the
!   convection control routine) into diagnostic arrays in module
!   cv_diagnostic_array_mod, so that they can be output later.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection
!
! code description:
!   language: fortran 90
!   This code is written to umdp3 programming standards version 10.8
!---------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SAVE_CONV_DIAGS_MOD'

CONTAINS

SUBROUTINE save_conv_diags( theta_inc, q_inc, qcl_inc, qcf_inc,        &
                            bulk_cf_inc, cf_liquid_inc, cf_frozen_inc, &
                            exner_theta_levels, l_calc_dxek )

USE atm_fields_bounds_mod, ONLY: tdims, tdims_s

USE cv_stash_flg_mod, ONLY:                                            &
    l_theta_incr_conv,                                                 &
    l_t_incr_conv, l_q_incr_conv, l_qcl_incr_conv, l_qcf_incr_conv,    &
    l_cfl_incr_conv, l_cff_incr_conv, l_bcf_incr_conv,                 &
    l_T_conv_only, l_q_conv_only, l_qcl_incr_cinh, l_qcf_incr_cinh,    &
    l_cfl_incr_cinh, l_cff_incr_cinh, l_bcf_incr_cinh

USE cv_diagnostic_array_mod, ONLY:                                     &
    theta_incr_diag_conv,                                              &
    T_incr_diag_conv, q_incr_diag_conv,                                &
    qcl_incr_diag_conv, qcf_incr_diag_conv,                            &
    cf_liquid_incr_diag_conv, cf_frozen_incr_diag_conv,                &
    bulk_cf_incr_diag_conv,                                            &
    T_incr_conv_only, q_incr_conv_only,                                &
    qcl_incr_inhom_diag, qcf_incr_inhom_diag,                          &
    cf_liquid_incr_inhom_diag, cf_frozen_incr_inhom_diag,              &
    bulk_cf_incr_inhom_diag

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

! Convection increments to thermodynamic, moisture and cloud fields.
! These are the fields after convection itself but before PC2 erosion.
REAL, INTENT(IN) :: theta_inc     ( tdims%i_start:tdims%i_end,         &
                                    tdims%j_start:tdims%j_end,         &
                                    1:tdims%k_end )
REAL, INTENT(IN) :: q_inc         ( tdims%i_start:tdims%i_end,         &
                                    tdims%j_start:tdims%j_end,         &
                                    1:tdims%k_end )
REAL, INTENT(IN) :: qcl_inc       ( tdims%i_start:tdims%i_end,         &
                                    tdims%j_start:tdims%j_end,         &
                                    1:tdims%k_end )
REAL, INTENT(IN) :: qcf_inc       ( tdims%i_start:tdims%i_end,         &
                                    tdims%j_start:tdims%j_end,         &
                                    1:tdims%k_end )
REAL, INTENT(IN) :: bulk_cf_inc   ( tdims%i_start:tdims%i_end,         &
                                    tdims%j_start:tdims%j_end,         &
                                    1:tdims%k_end )
REAL, INTENT(IN) :: cf_liquid_inc ( tdims%i_start:tdims%i_end,         &
                                    tdims%j_start:tdims%j_end,         &
                                    1:tdims%k_end )
REAL, INTENT(IN) :: cf_frozen_inc ( tdims%i_start:tdims%i_end,         &
                                    tdims%j_start:tdims%j_end,         &
                                    1:tdims%k_end )

! Exner, for computing T increment from theta increment
REAL, INTENT(IN) :: exner_theta_levels( tdims_s%i_start:tdims_s%i_end, &
                                        tdims_s%j_start:tdims_s%j_end, &
                                        tdims_s%k_start:tdims_s%k_end )

! Flag for whether convection is updating the prognostic cloud variables
LOGICAL, INTENT(IN) :: l_calc_dxek

! Loop counters
INTEGER :: i, j, k

CHARACTER(LEN=*), PARAMETER :: RoutineName='SAVE_CONV_DIAGS'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i, j, k)                         &
!$OMP& SHARED( tdims, theta_inc, q_inc, qcl_inc, qcf_inc,              &
!$OMP&         bulk_cf_inc, cf_liquid_inc, cf_frozen_inc,              &
!$OMP&         exner_theta_levels, l_calc_dxek,                        &
!$OMP&         l_theta_incr_conv, theta_incr_diag_conv,                &
!$OMP&         l_t_incr_conv, t_incr_diag_conv,                        &
!$OMP&         l_q_incr_conv, q_incr_diag_conv,                        &
!$OMP&         l_qcl_incr_conv, qcl_incr_diag_conv,                    &
!$OMP&         l_qcf_incr_conv, qcf_incr_diag_conv,                    &
!$OMP&         l_cfl_incr_conv, cf_liquid_incr_diag_conv,              &
!$OMP&         l_cff_incr_conv, cf_frozen_incr_diag_conv,              &
!$OMP&         l_bcf_incr_conv, bulk_cf_incr_diag_conv,                &
!$OMP&         l_T_conv_only, T_incr_conv_only,                        &
!$OMP&         l_q_conv_only, q_incr_conv_only,                        &
!$OMP&         l_qcl_incr_cinh, qcl_incr_inhom_diag,                   &
!$OMP&         l_qcf_incr_cinh, qcf_incr_inhom_diag,                   &
!$OMP&         l_cfl_incr_cinh, cf_liquid_incr_inhom_diag,             &
!$OMP&         l_cff_incr_cinh, cf_frozen_incr_inhom_diag,             &
!$OMP&         l_bcf_incr_cinh, bulk_cf_incr_inhom_diag )


! Total convective increment diagnostics.
! These get updated with increments from PC2 erosion
! in the call to pc2_from_conv_ctl.

IF (l_theta_incr_conv) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        theta_incr_diag_conv(i,j,k) = theta_inc(i,j,k)
      END DO ! i
    END DO ! j
  END DO  ! k
!$OMP END DO NOWAIT
END IF

IF (l_T_incr_conv) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        T_incr_diag_conv(i,j,k) = theta_inc(i,j,k)                     &
                                * exner_theta_levels(i,j,k)
      END DO ! i
    END DO ! j
  END DO  ! k
!$OMP END DO NOWAIT
END IF

IF (l_q_incr_conv) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q_incr_diag_conv(i,j,k) = q_inc(i,j,k)
      END DO ! i
    END DO ! j
  END DO  ! k
!$OMP END DO NOWAIT
END IF

IF (l_qcl_incr_conv) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qcl_incr_diag_conv(i,j,k) = qcl_inc(i,j,k)
      END DO ! i
    END DO ! j
  END DO  ! k
!$OMP END DO NOWAIT
END IF

IF (l_qcf_incr_conv) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qcf_incr_diag_conv(i,j,k) = qcf_inc(i,j,k)
      END DO ! i
    END DO ! j
  END DO  ! k
!$OMP END DO NOWAIT
END IF

IF ( l_calc_dxek ) THEN

  IF (l_cfl_incr_conv) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          cf_liquid_incr_diag_conv(i,j,k) = cf_liquid_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO  ! k
!$OMP END DO NOWAIT
  END IF

  IF (l_cff_incr_conv) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          cf_frozen_incr_diag_conv(i,j,k) = cf_frozen_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO  ! k
!$OMP END DO NOWAIT
  END IF

  IF (l_bcf_incr_conv) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          bulk_cf_incr_diag_conv(i,j,k) = bulk_cf_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO  ! k
!$OMP END DO NOWAIT
  END IF

END IF  ! ( l_calc_dxek )


! Increments from convection without the contribution from PC2 erosion.
! These get saved with their current values and are not updated by
! the call to pc2_from_conv_ctl.

IF (l_T_conv_only) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        T_incr_conv_only(i,j,k) = theta_inc(i,j,k)                     &
                                * exner_theta_levels(i,j,k)
      END DO ! i
    END DO ! j
  END DO  ! k
!$OMP END DO NOWAIT
END IF

IF (l_q_conv_only) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q_incr_conv_only(i,j,k) = q_inc(i,j,k)
      END DO ! i
    END DO ! j
  END DO  ! k
!$OMP END DO NOWAIT
END IF

IF (l_qcl_incr_cinh) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qcl_incr_inhom_diag(i,j,k) = qcl_inc(i,j,k)
      END DO ! i
    END DO ! j
  END DO  ! k
!$OMP END DO NOWAIT
END IF

IF (l_qcf_incr_cinh) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qcf_incr_inhom_diag(i,j,k) = qcf_inc(i,j,k)
      END DO ! i
    END DO ! j
  END DO  ! k
!$OMP END DO NOWAIT
END IF

IF ( l_calc_dxek ) THEN

  IF (l_cfl_incr_cinh) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          cf_liquid_incr_inhom_diag(i,j,k) = cf_liquid_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO  ! k
!$OMP END DO NOWAIT
  END IF

  IF (l_cff_incr_cinh) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          cf_frozen_incr_inhom_diag(i,j,k) = cf_frozen_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO  ! k
!$OMP END DO NOWAIT
  END IF

  IF (l_bcf_incr_cinh) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          bulk_cf_incr_inhom_diag(i,j,k) = bulk_cf_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO  ! k
!$OMP END DO NOWAIT
  END IF

END IF  ! ( l_calc_dxek )

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE save_conv_diags

END MODULE save_conv_diags_mod

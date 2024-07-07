! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_q_to_mix_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_Q_TO_MIX_MOD'

CONTAINS

SUBROUTINE eg_q_to_mix                                                &
                  (tdims,mixdims,                                     &
                   q, qcl, qcf,                                       &
                   qcf2, qrain, qgraup,                               &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,             &
                   mix_v, mix_cl, mix_cf,                             &
                   mix_cf2, mix_rain, mix_graup , swap_in             &
                   )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE Field_Types

USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE atm_fields_bounds_mod, ONLY: array_dims

IMPLICIT NONE
!
! Description:
!          Convert from specific humidities to mixing ratios
!
!
! Method:
!          See ENDGame Formulation version 3.03 Section 7.10
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.6.

! Subroutine arguments

TYPE (array_dims) , INTENT(IN) :: tdims, mixdims

! Arguments with Intent IN. ie: Input

REAL, INTENT(IN) :: q      (tdims%i_start:tdims%i_end,                &
                            tdims%j_start:tdims%j_end,                &
                            tdims%k_start:tdims%k_end)

REAL, INTENT(IN) :: qcl    (tdims%i_start:tdims%i_end,                &
                            tdims%j_start:tdims%j_end,                &
                            tdims%k_start:tdims%k_end)

REAL, INTENT(IN) :: qcf    (tdims%i_start:tdims%i_end,                &
                            tdims%j_start:tdims%j_end,                &
                            tdims%k_start:tdims%k_end)

REAL, INTENT(IN) :: qcf2   (tdims%i_start:tdims%i_end,                &
                            tdims%j_start:tdims%j_end,                &
                            tdims%k_start:tdims%k_end)

REAL, INTENT(IN) :: qrain  (tdims%i_start:tdims%i_end,                &
                            tdims%j_start:tdims%j_end,                &
                            tdims%k_start:tdims%k_end)

REAL, INTENT(IN) :: qgraup (tdims%i_start:tdims%i_end,                &
                            tdims%j_start:tdims%j_end,                &
                            tdims%k_start:tdims%k_end)

LOGICAL, INTENT(IN) :: l_mcr_qcf2     ! true if second cloud ice active
LOGICAL, INTENT(IN) :: l_mcr_qrain    ! true if rain active
LOGICAL, INTENT(IN) :: l_mcr_qgraup   ! true if graupel active

! Arguments with Intent INOUT.

! Note that the mixing ratio variables need to be INTENT(INOUT) as 
! initialisation happens before this routine is called and some of the
! variables (e.g. cf2) can be used without logical protection in the 
! solver. INTENT(OUT) can result in random memory taking the place of
! the initialisation, and thus failures or corruption.

REAL, INTENT(INOUT) :: mix_v     (mixdims%i_start:mixdims%i_end,        &
                                mixdims%j_start:mixdims%j_end,        &
                                mixdims%k_start:mixdims%k_end)

REAL, INTENT(INOUT) :: mix_cl    (mixdims%i_start:mixdims%i_end,        &
                                mixdims%j_start:mixdims%j_end,        &
                                mixdims%k_start:mixdims%k_end)

REAL, INTENT(INOUT) :: mix_cf    (mixdims%i_start:mixdims%i_end,        &
                                mixdims%j_start:mixdims%j_end,        &
                                mixdims%k_start:mixdims%k_end)

REAL, INTENT(INOUT) :: mix_cf2   (mixdims%i_start:mixdims%i_end,        &
                                mixdims%j_start:mixdims%j_end,        &
                                mixdims%k_start:mixdims%k_end)

REAL, INTENT(INOUT) :: mix_rain  (mixdims%i_start:mixdims%i_end,        &
                                mixdims%j_start:mixdims%j_end,        &
                                mixdims%k_start:mixdims%k_end)

REAL, INTENT(INOUT) :: mix_graup (mixdims%i_start:mixdims%i_end,        &
                                mixdims%j_start:mixdims%j_end,        &
                                mixdims%k_start:mixdims%k_end)

! Optional Arguments

LOGICAL, INTENT(IN), OPTIONAL :: swap_in

! local variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0 !DrHook
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1 !DrHook
REAL(KIND=jprb)               :: zhook_handle  !DrHook

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_Q_TO_MIX'

REAL :: moist_conv

INTEGER :: i, j, k

LOGICAL :: swap

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(swap_in)) THEN
  swap = swap_in
ELSE
  swap = .TRUE.
END IF

! ----------------------------------------------------------------------
! Section 1. convert q, qcl,qcf to mix_v, mix_cl,mix_cf
! ----------------------------------------------------------------------

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                      &
!$OMP& PRIVATE(k,j,i,moist_conv)                                       &
!$OMP& SHARED(mix_v,mix_cl,mix_cf,mix_cf2,mix_rain,                    &
!$OMP&        mix_graup,q,qcl,qcf,qcf2,qrain,qgraup,mixdims,           &
!$OMP&        l_mcr_qcf2,l_mcr_qrain,l_mcr_qgraup)
DO k=mixdims%k_start,mixdims%k_end
  DO j= mixdims%j_start + mixdims%halo_j, mixdims%j_end - mixdims%halo_j
    DO i= mixdims%i_start + mixdims%halo_i, mixdims%i_end - mixdims%halo_i

      moist_conv = 1 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k)

      IF (l_mcr_qcf2) THEN
        moist_conv = moist_conv - qcf2(i,j,k)
      END IF

      IF (l_mcr_qrain) THEN
        moist_conv = moist_conv - qrain(i,j,k)
      END IF

      IF (l_mcr_qgraup) THEN
        moist_conv = moist_conv - qgraup(i,j,k)
      END IF

      moist_conv = 1.0 / moist_conv

      mix_v (i,j,k) = q  (i,j,k) * moist_conv
      mix_cl(i,j,k) = qcl(i,j,k) * moist_conv
      mix_cf(i,j,k) = qcf(i,j,k) * moist_conv

      IF (l_mcr_qcf2) THEN
        mix_cf2(i,j,k) = qcf2(i,j,k) * moist_conv
      END IF

      IF (l_mcr_qrain) THEN
        mix_rain(i,j,k) = qrain(i,j,k) * moist_conv
      END IF

      IF (l_mcr_qgraup) THEN
        mix_graup(i,j,k) = qgraup(i,j,k) * moist_conv
      END IF

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (swap) THEN
  CALL swap_bounds(mix_v,                                                &
                   mixdims%i_end - mixdims%i_start+1 - 2*mixdims%halo_i, &
                   mixdims%j_end - mixdims%j_start+1 - 2*mixdims%halo_j, &
                   mixdims%k_end - mixdims%k_start + 1,                  &
                   mixdims%halo_i, mixdims%halo_j,                       &
                   fld_type_p, swap_field_is_scalar)
  CALL swap_bounds(mix_cl,                                               &
                   mixdims%i_end - mixdims%i_start+1 - 2*mixdims%halo_i, &
                   mixdims%j_end - mixdims%j_start+1 - 2*mixdims%halo_j, &
                   mixdims%k_end - mixdims%k_start + 1,                  &
                   mixdims%halo_i, mixdims%halo_j,                       &
                   fld_type_p, swap_field_is_scalar)
  CALL swap_bounds(mix_cf,                                               &
                   mixdims%i_end - mixdims%i_start+1 - 2*mixdims%halo_i, &
                   mixdims%j_end - mixdims%j_start+1 - 2*mixdims%halo_j, &
                   mixdims%k_end - mixdims%k_start + 1,                  &
                   mixdims%halo_i, mixdims%halo_j,                       &
                   fld_type_p, swap_field_is_scalar)
  IF (l_mcr_qcf2  )                                                      &
    CALL swap_bounds(mix_cf2,                                            &
                   mixdims%i_end - mixdims%i_start+1 - 2*mixdims%halo_i, &
                   mixdims%j_end - mixdims%j_start+1 - 2*mixdims%halo_j, &
                   mixdims%k_end - mixdims%k_start + 1,                  &
                   mixdims%halo_i, mixdims%halo_j,                       &
                   fld_type_p, swap_field_is_scalar)
  IF (l_mcr_qrain )                                                      &
    CALL swap_bounds(mix_rain,                                           &
                   mixdims%i_end - mixdims%i_start+1 - 2*mixdims%halo_i, &
                   mixdims%j_end - mixdims%j_start+1 - 2*mixdims%halo_j, &
                   mixdims%k_end - mixdims%k_start + 1,                  &
                   mixdims%halo_i, mixdims%halo_j,                       &
                   fld_type_p, swap_field_is_scalar)
  IF (l_mcr_qgraup)                                                      &
    CALL swap_bounds(mix_graup,                                          &
                   mixdims%i_end - mixdims%i_start+1 - 2*mixdims%halo_i, &
                   mixdims%j_end - mixdims%j_start+1 - 2*mixdims%halo_j, &
                   mixdims%k_end - mixdims%k_start + 1,                  &
                   mixdims%halo_i, mixdims%halo_j,                       &
                   fld_type_p, swap_field_is_scalar)

END IF

! end of routine

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE eg_q_to_mix

END MODULE eg_q_to_mix_mod

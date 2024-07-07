! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_mix_to_q_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_MIX_TO_Q_MOD'

CONTAINS

SUBROUTINE eg_mix_to_q                                                &
                  (tdims,mixdims,                                     &
                   mix_v, mix_cl, mix_cf,                             &
                   mix_cf2, mix_rain, mix_graup,                      &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,             &
                   q, qcl, qcf,                                       &
                   qcf2,qrain,qgraup                                  &
                   )


USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod, ONLY: array_dims

IMPLICIT NONE
!
! Description:
!
!          Convert from mixing ratios to specific humidities
!
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.6.

! Subroutine arguments

! Arguments with Intent IN. ie: Input

TYPE (array_dims) , INTENT(IN) :: tdims,mixdims

REAL, INTENT(IN) :: mix_v     (mixdims%i_start:mixdims%i_end,         &
                               mixdims%j_start:mixdims%j_end,         &
                               mixdims%k_start:mixdims%k_end)

REAL, INTENT(IN) :: mix_cl    (mixdims%i_start:mixdims%i_end,         &
                               mixdims%j_start:mixdims%j_end,         &
                               mixdims%k_start:mixdims%k_end)

REAL, INTENT(IN) :: mix_cf    (mixdims%i_start:mixdims%i_end,         &
                               mixdims%j_start:mixdims%j_end,         &
                               mixdims%k_start:mixdims%k_end)

REAL, INTENT(IN) :: mix_cf2   (mixdims%i_start:mixdims%i_end,         &
                               mixdims%j_start:mixdims%j_end,         &
                               mixdims%k_start:mixdims%k_end)

REAL, INTENT(IN) :: mix_rain  (mixdims%i_start:mixdims%i_end,         &
                               mixdims%j_start:mixdims%j_end,         &
                               mixdims%k_start:mixdims%k_end)

REAL, INTENT(IN) :: mix_graup (mixdims%i_start:mixdims%i_end,         &
                               mixdims%j_start:mixdims%j_end,         &
                               mixdims%k_start:mixdims%k_end)

LOGICAL, INTENT(IN) :: l_mcr_qcf2   ! true if second cloud ice active
LOGICAL, INTENT(IN) :: l_mcr_qrain  ! true if rain active
LOGICAL, INTENT(IN) :: l_mcr_qgraup ! true if graupel active

! Arguments with Intent OUT. ie: Output

REAL, INTENT(OUT) :: q      (tdims%i_start:tdims%i_end,               &
                             tdims%j_start:tdims%j_end,               &
                             tdims%k_start:tdims%k_end)

REAL, INTENT(OUT) :: qcl    (tdims%i_start:tdims%i_end,               &
                             tdims%j_start:tdims%j_end,               &
                             tdims%k_start:tdims%k_end)

REAL, INTENT(OUT) :: qcf    (tdims%i_start:tdims%i_end,               &
                             tdims%j_start:tdims%j_end,               &
                             tdims%k_start:tdims%k_end)

REAL, INTENT(OUT) :: qcf2   (tdims%i_start:tdims%i_end,               &
                             tdims%j_start:tdims%j_end,               &
                             tdims%k_start:tdims%k_end)

REAL, INTENT(OUT) :: qrain  (tdims%i_start:tdims%i_end,               &
                             tdims%j_start:tdims%j_end,               &
                             tdims%k_start:tdims%k_end)

REAL, INTENT(OUT) :: qgraup (tdims%i_start:tdims%i_end,               &
                             tdims%j_start:tdims%j_end,               &
                             tdims%k_start:tdims%k_end)

! local variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0 !DrHook
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1 !DrHook
REAL(KIND=jprb)               :: zhook_handle  !DrHook

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_MIX_TO_Q'

REAL :: moist_conv

INTEGER :: i, j, k

! ---------------------------------------------------------------------
! Section 1. convert mix_v, mix_cl,mix_cf to q, qcl,qcf
! ---------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                     &
!$OMP& PRIVATE(k,j,i,moist_conv)                                      &
!$OMP& SHARED(mix_v,mix_cl,mix_cf,mix_cf2,mix_rain,mix_graup,         &
!$OMP&        q,qcl,qcf,qcf2,qrain,qgraup,mixdims,                    &
!$OMP&        l_mcr_qcf2,l_mcr_qrain,l_mcr_qgraup)
DO k=mixdims%k_start,mixdims%k_end
  DO j= mixdims%j_start, mixdims%j_end
    DO i= mixdims%i_start, mixdims%i_end


      moist_conv = 1.0 + mix_v(i,j,k) + mix_cl(i,j,k) + mix_cf(i,j,k)

      IF (l_mcr_qcf2) THEN
        moist_conv = moist_conv + mix_cf2(i,j,k)
      END IF

      IF (l_mcr_qrain) THEN
        moist_conv = moist_conv + mix_rain(i,j,k)
      END IF

      IF (l_mcr_qgraup) THEN
        moist_conv = moist_conv + mix_graup(i,j,k)
      END IF

      moist_conv = 1.0 / moist_conv
      q  (i,j,k) = mix_v (i,j,k) * moist_conv
      qcl(i,j,k) = mix_cl(i,j,k) * moist_conv
      qcf(i,j,k) = mix_cf(i,j,k) * moist_conv

      IF (l_mcr_qcf2) THEN
        qcf2(i,j,k)  = mix_cf2(i,j,k) * moist_conv
      END IF

      IF (l_mcr_qrain) THEN
        qrain(i,j,k)  = mix_rain(i,j,k) * moist_conv
      END IF

      IF (l_mcr_qgraup) THEN
        qgraup(i,j,k)  = mix_graup(i,j,k) * moist_conv
      END IF

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! end of routine

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_mix_to_q
END MODULE eg_mix_to_q_mod

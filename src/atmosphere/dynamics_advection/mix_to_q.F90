! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine mix_to_q

SUBROUTINE mix_to_q                                               &
                  (row_length, rows, halo_i, halo_j,              &
                   mix_v, mix_cl, mix_cf,                         &
                   mix_cf2, mix_rain, mix_graup,                  &
                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
                   q, qcl, qcf                                    &
                  ,qcf2,qrain,qgraup                              &
                   )

! Purpose:
!          Convert from mixing ratios to specific humidities
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE nlsizes_namelist_mod, ONLY: model_levels

USE yomhook, ONLY: lhook, dr_hook

USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER ::                                                        &
  row_length                                                      &
                  ! number of points on a row
, rows                                                            &
                  ! number of rows of data
, halo_i                                                          &
, halo_j

LOGICAL ::                                                        &
  L_mcr_qcf2                                                      &
                    ! true if second cloud ice active
, L_mcr_qrain                                                     &
                    ! true if rain active
, L_mcr_qgraup      ! true if graupel active

! Arguments with Intent IN. ie: Input

REAL ::                                                           &
  mix_v  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
          model_levels)                                           &
, mix_cl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
          model_levels)                                           &
, mix_cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
          model_levels)                                           &
, mix_cf2   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
             model_levels)                                        &
, mix_rain  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
             model_levels)                                        &
, mix_graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
             model_levels)

! Arguments with Intent OUT. ie: Output

REAL ::                                                           &
  q    (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
        model_levels)                                             &
, qcl  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
        model_levels)                                             &
, qcf  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
        model_levels)                                             &
, qcf2    (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
           model_levels)                                          &
, qrain   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
           model_levels)                                          &
, qgraup  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
           model_levels)

! local variables
REAL :: moist   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)

INTEGER ::                                                        &
 i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MIX_TO_Q'


! ----------------------------------------------------------------------
! Section 1. convert mix_v, mix_cl,mix_cf to q, qcl,qcf
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k,moist)
DO k = 1, model_levels

  DO j = 1, rows
    DO i = 1, row_length
      moist(i,j) =1.0+ mix_v (i,j,k)+mix_cl(i,j,k)+mix_cf(i,j,k)
    END DO
  END DO

  IF (L_mcr_qcf2) THEN

    DO j = 1, rows
      DO i = 1, row_length
        moist(i,j) = moist(i,j)+mix_cf2(i,j,k)
      END DO
    END DO

  END IF
  IF (L_mcr_qrain) THEN

    DO j = 1, rows
      DO i = 1, row_length
        moist(i,j)= moist(i,j)+mix_rain(i,j,k)
      END DO
    END DO

  END IF
  IF (L_mcr_qgraup) THEN

    DO j = 1, rows
      DO i = 1, row_length
        moist(i,j)= moist(i,j)+mix_graup(i,j,k)
      END DO
    END DO

  END IF

  DO j = 1, rows
    DO i = 1, row_length
      moist(i,j)= 1.0/moist(i,j)
    END DO
  END DO

  DO j = 1, rows
    DO i = 1, row_length
      q  (i,j,k) = mix_v (i,j,k) * moist(i,j)
      qcl(i,j,k) = mix_cl(i,j,k) * moist(i,j)
      qcf(i,j,k) = mix_cf(i,j,k) * moist(i,j)
    END DO
  END DO

  IF (L_mcr_qcf2) THEN

    DO j = 1, rows
      DO i = 1, row_length
        qcf2(i,j,k)  = mix_cf2(i,j,k) * moist(i,j)
      END DO
    END DO

  END IF
  IF (L_mcr_qrain) THEN

    DO j = 1, rows
      DO i = 1, row_length
        qrain(i,j,k)  = mix_rain(i,j,k) * moist(i,j)
      END DO
    END DO

  END IF
  IF (L_mcr_qgraup) THEN

    DO j = 1, rows
      DO i = 1, row_length
        qgraup(i,j,k)  = mix_graup(i,j,k) * moist(i,j)
      END DO
    END DO

  END IF

END DO
!$OMP END PARALLEL DO

! end of routine

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE mix_to_q

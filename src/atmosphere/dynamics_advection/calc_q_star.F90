! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine calc_q_star

SUBROUTINE calc_q_star                                            &
                  (row_length, rows,                              &
                   halo_i, halo_j, offx, offy,                    &
                   mix_v, mix_cl, mix_cf,                         &
                   mix_cf2, mix_rain, mix_graup,                  &
                   mix_v_star, mix_cl_star, mix_cf_star,          &
                   mix_cf2_star, mix_rain_star, mix_graup_star,   &
                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
                   q, qcl, qcf,                                   &
                   qcf2, qrain, qgraup,                           &
                   q_star, qcl_star, qcf_star                     &
                  ,qcf2_star, qrain_star, qgraup_star             &
                   )

! Purpose:
!          calculate specific humidity increments
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: Fortran 90 + extensions
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
, halo_j                                                          &
, offx                                                            &
, offy

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

REAL ::                                                           &
  mix_v_star  (row_length,rows, model_levels)                     &
, mix_cl_star (row_length,rows, model_levels)                     &
, mix_cf_star (row_length,rows, model_levels)                     &
, mix_cf2_star   (row_length,rows, model_levels)                  &
, mix_rain_star  (row_length,rows, model_levels)                  &
, mix_graup_star (row_length,rows, model_levels)

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

! Arguments with Intent OUT. ie: Output

REAL ::                                                           &
  q_star     (row_length, rows, model_levels)                     &
, qcl_star   (row_length, rows, model_levels)                     &
, qcf_star   (row_length, rows, model_levels)                     &
, qcf2_star  (row_length, rows, model_levels)                     &
, qrain_star (row_length, rows, model_levels)                     &
, qgraup_star(row_length, rows, model_levels)

! local variables
REAL ::                                                           &
 conv,sum_mix,sum_mix_star                                        &
, moist(row_length, rows, model_levels)                           &
, moist_star(row_length, rows, model_levels)

INTEGER ::                                                        &
 i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_Q_STAR'


! ----------------------------------------------------------------------
! Section 1. calculate q_star, qcl_star, qcf_star
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,conv)

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      moist(i,j,k) = 1 + mix_v(i,j,k) +                           &
                         mix_cl(i,j,k) + mix_cf(i,j,k)
      moist_star(i,j,k) =  mix_v_star(i,j,k) +                    &
                         mix_cl_star(i,j,k) + mix_cf_star(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO
IF (L_mcr_qcf2) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        moist(i,j,k)      = moist(i,j,k) + mix_cf2(i,j,k)
        moist_star(i,j,k) =                                       &
                       moist_star(i,j,k) + mix_cf2_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

END IF
IF (L_mcr_qrain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        moist(i,j,k)      = moist(i,j,k) + mix_rain(i,j,k)
        moist_star(i,j,k) =                                       &
                       moist_star(i,j,k) + mix_rain_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
END IF
IF (L_mcr_qgraup) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        moist(i,j,k)      = moist(i,j,k) + mix_graup(i,j,k)
        moist_star(i,j,k) =                                       &
                       moist_star(i,j,k) + mix_graup_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
END IF

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      conv= 1.0/(  moist(i,j,k) *                                &
                 (moist(i,j,k) + moist_star(i,j,k)) )
      q_star(i,j,k)  =                                          &
             ( mix_v_star(i,j,k)*moist(i,j,k) -                 &
               mix_v(i,j,k)*moist_star(i,j,k) )  * conv
      qcl_star(i,j,k)  =                                        &
             ( mix_cl_star(i,j,k)*moist(i,j,k) -                &
               mix_cl(i,j,k)*moist_star(i,j,k) )  * conv
      qcf_star(i,j,k)  =                                        &
             ( mix_cf_star(i,j,k)*moist(i,j,k) -                &
               mix_cf(i,j,k)*moist_star(i,j,k) )  * conv

    END DO
  END DO
END DO
!$OMP END DO

IF (L_mcr_qcf2) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        conv= 1.0/(  moist(i,j,k) *                                &
                   (moist(i,j,k) + moist_star(i,j,k)) )
        qcf2_star(i,j,k)  =                                       &
               ( mix_cf2_star(i,j,k)*moist(i,j,k) -               &
                 mix_cf2(i,j,k)*moist_star(i,j,k) )  * conv
      END DO
    END DO
  END DO
!$OMP END DO
END IF
IF (L_mcr_qrain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        conv= 1.0/(  moist(i,j,k) *                                &
                   (moist(i,j,k) + moist_star(i,j,k)) )
        qrain_star(i,j,k)  =                                      &
               ( mix_rain_star(i,j,k)*moist(i,j,k) -              &
                 mix_rain(i,j,k)*moist_star(i,j,k) ) * conv
      END DO
    END DO
  END DO
!$OMP END DO
END IF
IF (L_mcr_qgraup) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        conv= 1.0/(  moist(i,j,k) *                                &
                   (moist(i,j,k) + moist_star(i,j,k)) )
        qgraup_star(i,j,k)  =                                     &
               ( mix_graup_star(i,j,k)*moist(i,j,k) -             &
                 mix_graup(i,j,k)*moist_star(i,j,k) ) * conv
      END DO
    END DO
  END DO
!$OMP END DO
END IF

! end of routine

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_q_star


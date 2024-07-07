! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to deal with the calculation of increment diagnostics for advection
! correction.
! Declares variables and allocatable arrays which are used in the calculation
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to the UM programming standards in UMDP3

MODULE adv_correct_incs_mod

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain

USE atm_fields_bounds_mod, ONLY: pdims_s, tdims_s, tdims

USE stash_array_mod, ONLY: sf

IMPLICIT NONE

! STASH section for diagnostics
INTEGER, PARAMETER :: sect = 12

! Variables for increment diagnostics
REAL, ALLOCATABLE::                                                     &
    advcor_inc_t(:,:,:), advcor_inc_m_v(:,:,:), advcor_inc_m_cl(:,:,:), &
    advcor_inc_m_cf(:,:,:), advcor_inc_m_r(:,:,:),                      &
    advcor_inc_m_gr(:,:,:), advcor_inc_m_cf2(:,:,:),                    &
    advcor_inc_q(:,:,:), advcor_inc_qcl(:,:,:), advcor_inc_qcf(:,:,:),  &
    advcor_inc_qrain(:,:,:), advcor_inc_qgraup(:,:,:),                  &
    advcor_inc_qcf2(:,:,:)

! Loop counters
INTEGER :: i, j, k

CONTAINS

SUBROUTINE adv_correct_incs_init(exner, thetav,                         &
                                 m_v, m_cl, m_cf, m_r, m_gr, m_cf2)

  USE planet_constants_mod, ONLY: recip_epsilon
  USE horiz_grid_mod, ONLY: intw_rho2w

  IMPLICIT NONE

  ! Subroutine arguments

  REAL, INTENT(IN) :: exner  (pdims_s%i_start:pdims_s%i_end,                 &
                              pdims_s%j_start:pdims_s%j_end,                 &
                              pdims_s%k_start:pdims_s%k_end+1)
  REAL, INTENT(IN) :: thetav (tdims_s%i_start:tdims_s%i_end,                 &
                              tdims_s%j_start:tdims_s%j_end,                 &
                              tdims_s%k_start:tdims_s%k_end)
  ! Water Vapour Mixing Ratio
  REAL, INTENT(IN) :: m_v    (tdims_s%i_start:tdims_s%i_end,                 &
                              tdims_s%j_start:tdims_s%j_end,                 &
                              tdims_s%k_start:tdims_s%k_end)
  ! Cloud Liquid Water Mixing Ratio
  REAL, INTENT(IN) :: m_cl   (tdims_s%i_start:tdims_s%i_end,                 &
                              tdims_s%j_start:tdims_s%j_end,                 &
                              tdims_s%k_start:tdims_s%k_end)
  ! Cloud Ice Mixing Ratio
  REAL, INTENT(IN) :: m_cf   (tdims_s%i_start:tdims_s%i_end,                 &
                              tdims_s%j_start:tdims_s%j_end,                 &
                              tdims_s%k_start:tdims_s%k_end)
  ! Rain Mixing Ratio
  REAL, INTENT(IN) :: m_r    (tdims_s%i_start:tdims_s%i_end,                 &
                              tdims_s%j_start:tdims_s%j_end,                 &
                              tdims_s%k_start:tdims_s%k_end)
  ! Graupel Mixing Ratio
  REAL, INTENT(IN) :: m_gr   (tdims_s%i_start:tdims_s%i_end,                 &
                              tdims_s%j_start:tdims_s%j_end,                 &
                              tdims_s%k_start:tdims_s%k_end)
  ! Cloud Ice (2nd Species) Mixing Ratio
  REAL, INTENT(IN) :: m_cf2  (tdims_s%i_start:tdims_s%i_end,                 &
                              tdims_s%j_start:tdims_s%j_end,                 &
                              tdims_s%k_start:tdims_s%k_end)

  ! Local Variables

  REAL :: exner_theta

  REAL :: moist_conv(tdims%i_start:tdims%i_end, &
                     tdims%j_start:tdims%j_end, &
                                 1:tdims%k_end)

  !---------------------------------------------------------------------------
  ! Allocate and initialise arrays for advection correction increments
  !---------------------------------------------------------------------------

  ! Calculate moisture conversion
  IF ( sf(382,sect) .OR. sf(383,sect) .OR. sf(384,sect) .OR.            &
       (l_mcr_qrain .AND. sf(389,sect)) .OR.                            &
       (l_mcr_qgraup .AND. sf(390,sect)) .OR.                           &
       (l_mcr_qcf2 .AND. sf(391,sect)) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,moist_conv,m_v,m_cl,m_cf,m_r,m_gr,m_cf2,             &
!$OMP        l_mcr_qgraup,l_mcr_qrain,l_mcr_qcf2)
    DO k=            1,tdims%k_end
      DO j=tdims%j_start,tdims%j_end
        DO i=tdims%i_start,tdims%i_end
          moist_conv(i,j,k) = 1.0 + m_v(i,j,k) + m_cl(i,j,k) + m_cf(i,j,k)
          IF (l_mcr_qcf2) THEN
            moist_conv(i,j,k) = moist_conv(i,j,k) + m_cf2(i,j,k)
          END IF
          IF (l_mcr_qrain) THEN
            moist_conv(i,j,k) = moist_conv(i,j,k) + m_r(i,j,k)
          END IF
          IF (l_mcr_qgraup) THEN
            moist_conv(i,j,k) = moist_conv(i,j,k) + m_gr(i,j,k)
          END IF
          moist_conv(i,j,k) = 1.0 / moist_conv(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Temperature
  IF ( sf(381,sect) ) THEN
    ALLOCATE ( advcor_inc_t(tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                                        1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i,exner_theta)                                        &
!$OMP SHARED(tdims,advcor_inc_t,intw_rho2w,exner,thetav,m_v,recip_epsilon)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          exner_theta = (intw_rho2w(k,2)*exner(i,j,k) +                 &
                          intw_rho2w(k,1)*exner(i,j,k+1))
          advcor_inc_t(i,j,k) = exner_theta * thetav(i,j,k)             &
                              / (1.0 + m_v(i,j,k)*recip_epsilon)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_t( 1,1,1) )
  END IF

  ! Specific Humidity
  IF ( sf(382,sect) ) THEN
    ALLOCATE ( advcor_inc_q(tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                                        1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_q,m_v,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_q(i,j,k) = m_v(i,j,k) * moist_conv(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_q( 1,1,1) )
  END IF

  ! Cloud Liquid Water
  IF ( sf(383,sect) ) THEN
    ALLOCATE ( advcor_inc_qcl(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_qcl,m_cl,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_qcl(i,j,k) = m_cl(i,j,k) * moist_conv(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_qcl( 1,1,1) )
  END IF

  ! Cloud Ice
  IF ( sf(384,sect) ) THEN
    ALLOCATE ( advcor_inc_qcf(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_qcf,m_cf,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_qcf(i,j,k) = m_cf(i,j,k) * moist_conv(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_qcf( 1,1,1) )
  END IF

  ! Prognostic Rain
  IF ( l_mcr_qrain .AND. sf(389,sect) ) THEN
    ALLOCATE ( advcor_inc_qrain(tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                            1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_qrain,m_r,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_qrain(i,j,k) = m_r(i,j,k) * moist_conv(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_qrain( 1,1,1) )
  END IF

  ! Prognostic Graupel
  IF ( l_mcr_qgraup .AND. sf(390,sect) ) THEN
    ALLOCATE ( advcor_inc_qgraup(tdims%i_start:tdims%i_end,             &
                                 tdims%j_start:tdims%j_end,             &
                                             1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_qgraup,m_gr,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_qgraup(i,j,k) = m_gr(i,j,k) * moist_conv(i,j,k)
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE ( advcor_inc_qgraup( 1,1,1) )
  END IF

  ! Prognostic Cloud Ice - 2nd species
  IF ( l_mcr_qcf2 .AND. sf(391,sect) ) THEN
    ALLOCATE ( advcor_inc_qcf2(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,               &
                                           1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_qcf2,m_cf2,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_qcf2(i,j,k) = m_cf2(i,j,k) * moist_conv(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_qcf2( 1,1,1) )
  END IF

  ! Water Vapour Mixing Ratio
  IF ( sf(395,sect) ) THEN
    ALLOCATE ( advcor_inc_m_v(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_v,m_v)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_v(i,j,k) = m_v(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_m_v( 1,1,1) )
  END IF

  ! Cloud Liquid Water Mixing Ratio
  IF ( sf(396,sect) ) THEN
    ALLOCATE ( advcor_inc_m_cl(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,               &
                                           1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_cl,m_cl)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_cl(i,j,k) = m_cl(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_m_cl( 1,1,1) )
  END IF

  ! Cloud Ice Mixing Ratio
  IF ( sf(397,sect) ) THEN
    ALLOCATE ( advcor_inc_m_cf(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,               &
                                           1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_cf,m_cf)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_cf(i,j,k) = m_cf(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_m_cf( 1,1,1) )
  END IF

  ! Rain Mixing Ratio
  IF (l_mcr_qrain .AND. sf(398,sect) ) THEN
    ALLOCATE ( advcor_inc_m_r(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_r,m_r)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_r(i,j,k) = m_r(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_m_r( 1,1,1) )
  END IF

  ! Graupel Mixing Ratio
  IF (l_mcr_qgraup .AND. sf(399,sect) ) THEN
    ALLOCATE ( advcor_inc_m_gr(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,               &
                                           1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_gr,m_gr)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_gr(i,j,k) = m_gr(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_m_gr( 1,1,1) )
  END IF

  ! Cloud Ice (2nd species) Mixing Ratio
  IF (l_mcr_qcf2 .AND. sf(400,sect) ) THEN
    ALLOCATE ( advcor_inc_m_cf2(tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                            1:tdims%k_end))

!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_cf2,m_cf2)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_cf2(i,j,k) = m_cf2(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( advcor_inc_m_cf2( 1,1,1) )
  END IF

END SUBROUTINE adv_correct_incs_init

SUBROUTINE adv_correct_incs_calc(exner, thetav,                         &
                                 m_v, m_cl, m_cf, m_r, m_gr, m_cf2)

  USE planet_constants_mod, ONLY: recip_epsilon
  USE horiz_grid_mod, ONLY: intw_rho2w

  IMPLICIT NONE

  ! Subroutine arguments

  REAL, INTENT(IN) :: exner  (pdims_s%i_start:pdims_s%i_end,            &
                              pdims_s%j_start:pdims_s%j_end,            &
                              pdims_s%k_start:pdims_s%k_end+1)
  REAL, INTENT(IN) :: thetav (tdims_s%i_start:tdims_s%i_end,            &
                              tdims_s%j_start:tdims_s%j_end,            &
                              tdims_s%k_start:tdims_s%k_end)
  ! Water Vapour Mixing Ratio
  REAL, INTENT(IN) :: m_v    (tdims_s%i_start:tdims_s%i_end,            &
                              tdims_s%j_start:tdims_s%j_end,            &
                              tdims_s%k_start:tdims_s%k_end)
  ! Cloud Liquid Water Mixing Ratio
  REAL, INTENT(IN) :: m_cl   (tdims_s%i_start:tdims_s%i_end,            &
                              tdims_s%j_start:tdims_s%j_end,            &
                              tdims_s%k_start:tdims_s%k_end)
  ! Cloud Ice Mixing Ratio
  REAL, INTENT(IN) :: m_cf   (tdims_s%i_start:tdims_s%i_end,            &
                              tdims_s%j_start:tdims_s%j_end,            &
                              tdims_s%k_start:tdims_s%k_end)
  ! Rain Mixing Ratio
  REAL, INTENT(IN) :: m_r    (tdims_s%i_start:tdims_s%i_end,            &
                              tdims_s%j_start:tdims_s%j_end,            &
                              tdims_s%k_start:tdims_s%k_end)
  ! Graupel Mixing Ratio
  REAL, INTENT(IN) :: m_gr   (tdims_s%i_start:tdims_s%i_end,            &
                              tdims_s%j_start:tdims_s%j_end,            &
                              tdims_s%k_start:tdims_s%k_end)
  ! Cloud Ice (2nd Species) Mixing Ratio
  REAL, INTENT(IN) :: m_cf2  (tdims_s%i_start:tdims_s%i_end,            &
                              tdims_s%j_start:tdims_s%j_end,            &
                              tdims_s%k_start:tdims_s%k_end)

  ! Local Variables

  REAL :: temp, exner_theta

  REAL :: moist_conv(tdims%i_start:tdims%i_end, &
                     tdims%j_start:tdims%j_end, &
                                 1:tdims%k_end)

  !---------------------------------------------------------------------------
  ! Calculate moisture correction increments
  !---------------------------------------------------------------------------

  ! Calculate moisture conversion
  IF ( sf(382,sect) .OR. sf(383,sect) .OR. sf(384,sect) .OR.            &
       (l_mcr_qrain .AND. sf(389,sect)) .OR.                            &
       (l_mcr_qgraup .AND. sf(390,sect)) .OR.                           &
       (l_mcr_qcf2 .AND. sf(391,sect)) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,moist_conv,m_v,m_cl,m_cf,m_r,m_gr,m_cf2,             &
!$OMP        l_mcr_qgraup,l_mcr_qrain,l_mcr_qcf2)
    DO k=            1,tdims%k_end
      DO j=tdims%j_start,tdims%j_end
        DO i=tdims%i_start,tdims%i_end
          moist_conv(i,j,k) = 1.0 + m_v(i,j,k) + m_cl(i,j,k) + m_cf(i,j,k)
          IF (l_mcr_qcf2) THEN
            moist_conv(i,j,k) = moist_conv(i,j,k) + m_cf2(i,j,k)
          END IF
          IF (l_mcr_qrain) THEN
            moist_conv(i,j,k) = moist_conv(i,j,k) + m_r(i,j,k)
          END IF
          IF (l_mcr_qgraup) THEN
            moist_conv(i,j,k) = moist_conv(i,j,k) + m_gr(i,j,k)
          END IF
          moist_conv(i,j,k) = 1.0 / moist_conv(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Temperature
  IF ( sf(381,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i,exner_theta,temp)                                   &
!$OMP SHARED(tdims,advcor_inc_t,intw_rho2w,exner,thetav,m_v,recip_epsilon)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          exner_theta = (intw_rho2w(k,2)*exner(i,j,k) +                 &
                         intw_rho2w(k,1)*exner(i,j,k+1))
          temp = exner_theta * thetav(i,j,k)          &
                  / (1.0 + m_v(i,j,k)*recip_epsilon)
          advcor_inc_t(i,j,k) = temp - advcor_inc_t(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Specific Humidity
  IF ( sf(382,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_q,m_v,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_q(i,j,k) = m_v(i,j,k) * moist_conv(i,j,k)          &
                              - advcor_inc_q(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Cloud Liquid Water
  IF ( sf(383,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_qcl,m_cl,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_qcl(i,j,k) = m_cl(i,j,k) * moist_conv(i,j,k)       &
                                - advcor_inc_qcl(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Cloud Ice
  IF ( sf(384,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_qcf,m_cf,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_qcf(i,j,k) = m_cf(i,j,k) * moist_conv(i,j,k)       &
                                - advcor_inc_qcf(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Prognostic Rain
  IF (l_mcr_qrain .AND. sf(389,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_qrain,m_r,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_qrain(i,j,k) = m_r(i,j,k) * moist_conv(i,j,k)      &
                                  - advcor_inc_qrain(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Prognostic Graupel
  IF (l_mcr_qgraup .AND. sf(390,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_qgraup,m_gr,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_qgraup(i,j,k) = m_gr(i,j,k) * moist_conv(i,j,k)    &
                                    - advcor_inc_qgraup(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Prognostic Cloud Ice - 2nd species
  IF (l_mcr_qcf2 .AND. sf(391,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_qcf2,m_cf2,moist_conv)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_qcf2(i,j,k) = m_cf2(i,j,k) * moist_conv(i,j,k)     &
                                  - advcor_inc_qcf2(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Water Vapour Mixing Ratio
  IF ( sf(395,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_v,m_v)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_v(i,j,k) = m_v(i,j,k) - advcor_inc_m_v(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Cloud Liquid Water Mixing Ratio
  IF ( sf(396,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_cl,m_cl)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_cl(i,j,k) = m_cl(i,j,k) - advcor_inc_m_cl(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Cloud Ice Mixing Ratio
  IF ( sf(397,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_cf,m_cf)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_cf(i,j,k) = m_cf(i,j,k) - advcor_inc_m_cf(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Rain Mixing Ratio
  IF (l_mcr_qrain .AND. sf(398,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_r,m_r)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_r(i,j,k) = m_r(i,j,k) - advcor_inc_m_r(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Graupel Mixing Ratio
  IF (l_mcr_qgraup .AND. sf(399,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_gr,m_gr)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_gr(i,j,k) = m_gr(i,j,k) - advcor_inc_m_gr(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Cloud Ice (2nd Species) Mixing Ratio
  IF (l_mcr_qcf2 .AND. sf(400,sect) ) THEN
!$OMP PARALLEL DO                                                       &
!$OMP SCHEDULE(STATIC)                                                  &
!$OMP DEFAULT(NONE)                                                     &
!$OMP PRIVATE(k,j,i)                                                    &
!$OMP SHARED(tdims,advcor_inc_m_cf2,m_cf2)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          advcor_inc_m_cf2(i,j,k) = m_cf2(i,j,k) - advcor_inc_m_cf2(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

END SUBROUTINE adv_correct_incs_calc

SUBROUTINE adv_correct_incs_dealloc()
  IMPLICIT NONE
  IF (ALLOCATED(advcor_inc_t)) DEALLOCATE ( advcor_inc_t )
  IF (ALLOCATED(advcor_inc_q)) DEALLOCATE ( advcor_inc_q )
  IF (ALLOCATED(advcor_inc_qcl)) DEALLOCATE ( advcor_inc_qcl )
  IF (ALLOCATED(advcor_inc_qcf)) DEALLOCATE ( advcor_inc_qcf )
  IF (ALLOCATED(advcor_inc_qrain)) DEALLOCATE ( advcor_inc_qrain )
  IF (ALLOCATED(advcor_inc_qgraup)) DEALLOCATE ( advcor_inc_qgraup )
  IF (ALLOCATED(advcor_inc_qcf2)) DEALLOCATE ( advcor_inc_qcf2 )
  IF (ALLOCATED(advcor_inc_m_v)) DEALLOCATE ( advcor_inc_m_v )
  IF (ALLOCATED(advcor_inc_m_cl)) DEALLOCATE ( advcor_inc_m_cl )
  IF (ALLOCATED(advcor_inc_m_cf)) DEALLOCATE ( advcor_inc_m_cf )
  IF (ALLOCATED(advcor_inc_m_r)) DEALLOCATE ( advcor_inc_m_r )
  IF (ALLOCATED(advcor_inc_m_gr)) DEALLOCATE ( advcor_inc_m_gr )
  IF (ALLOCATED(advcor_inc_m_cf2)) DEALLOCATE ( advcor_inc_m_cf2 )
END SUBROUTINE adv_correct_incs_dealloc

END MODULE adv_correct_incs_mod

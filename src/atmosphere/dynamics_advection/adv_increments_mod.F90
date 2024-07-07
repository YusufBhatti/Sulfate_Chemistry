! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to deal with the calculation of advection increment diagnostics
! Declares variables and allocatable arrays which are used in the calculation
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to the UM programming standards in UMDP3

MODULE adv_increments_mod

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
USE cloud_inputs_mod, ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2

USE atm_fields_bounds_mod, ONLY: udims_s, vdims_s, tdims, udims, vdims

USE stash_array_mod, ONLY: sf

IMPLICIT NONE

! STASH section for diagnostics
INTEGER, PARAMETER :: sect = 12

! Variables for increment diagnostics
REAL, ALLOCATABLE::                                                          &
    adv_inc_t(:,:,:), adv_inc_u(:,:,:), adv_inc_v(:,:,:), adv_inc_q(:,:,:),  &
    adv_inc_qcl(:,:,:), adv_inc_qcf(:,:,:), adv_inc_qrain(:,:,:),            &
    adv_inc_qgraup(:,:,:), adv_inc_qcf2(:,:,:), adv_inc_m_v(:,:,:),          &
    adv_inc_m_cl(:,:,:), adv_inc_m_cf(:,:,:), adv_inc_m_r(:,:,:),            &
    adv_inc_m_gr(:,:,:), adv_inc_m_cf2(:,:,:), adv_inc_cf(:,:,:),            &
    adv_inc_cfl(:,:,:), adv_inc_cff(:,:,:)

! Loop counters
INTEGER :: i, j, k

CONTAINS

SUBROUTINE adv_incs_init(r_u, r_v, theta_star,                               &
    q_star, qcl_star, qcf_star, qrain_star, qgraup_star, qcf2_star,          &
    m_star, mcl_star, mcf_star, mrain_star, mgraup_star, mcf2_star,          &
    cf_star, cfl_star, cff_star)

  IMPLICIT NONE

  ! Subroutine arguments

  REAL, TARGET     :: r_u        (udims_s%i_start:udims_s%i_end,             &
                                  udims_s%j_start:udims_s%j_end,             &
                                  udims_s%k_start:udims_s%k_end)
  REAL, TARGET     :: r_v        (vdims_s%i_start:vdims_s%i_end,             &
                                  vdims_s%j_start:vdims_s%j_end,             &
                                  vdims_s%k_start:vdims_s%k_end)

  REAL, INTENT(IN) :: theta_star (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)

  ! Water Vapour
  REAL, INTENT(IN) :: q_star     (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  ! Cloud Liquid Water 
  REAL, INTENT(IN) :: qcl_star   (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  ! Cloud Ice
  REAL, INTENT(IN) :: qcf_star   (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  ! Rain
  REAL, INTENT(IN) :: qrain_star (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  ! Graupel
  REAL, INTENT(IN) :: qgraup_star(tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  ! Cloud Ice (2nd Species)
  REAL, INTENT(IN) :: qcf2_star  (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)

  ! Water Vapour Mixing Ratio
  REAL, INTENT(IN) :: m_star     (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  ! Cloud Liquid Water Mixing Ratio
  REAL, INTENT(IN) :: mcl_star   (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  ! Cloud Ice Mixing Ratio
  REAL, INTENT(IN) :: mcf_star   (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  ! Rain Mixing Ratio
  REAL, INTENT(IN) :: mrain_star (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  ! Graupel Mixing Ratio
  REAL, INTENT(IN) :: mgraup_star(tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  ! Cloud Ice (2nd Species) Mixing Ratio
  REAL, INTENT(IN) :: mcf2_star  (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)

  REAL, INTENT(IN) :: cf_star    (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  REAL, INTENT(IN) :: cfl_star   (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)
  REAL, INTENT(IN) :: cff_star   (tdims%i_start:tdims%i_end,                 &
                                  tdims%j_start:tdims%j_end,                 &
                                  tdims%k_start:tdims%k_end)

  !---------------------------------------------------------------------------
  ! Allocate and initialise arrays for advection increments
  !---------------------------------------------------------------------------

  ! Temperature
  ! Note that we are holding theta, conversion to temperature will be done
  !  in calculation of increment.
  IF ( sf(181,sect) ) THEN
    ALLOCATE ( adv_inc_t(tdims%i_start:tdims%i_end,             &
                         tdims%j_start:tdims%j_end,             &
                                     1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_t,theta_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_t(i,j,k) = theta_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_t( 1,1,1) )
  END IF

  ! Prognostic water vapour
  IF (sf(182,sect)) THEN
    ALLOCATE ( adv_inc_q  (tdims%i_start:tdims%i_end,           &
                           tdims%j_start:tdims%j_end,           &
                                       1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_q,q_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_q(i,j,k) = q_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_q(1,1,1) )
  END IF

  ! Prognostic cloud liquid water
  IF (sf(183,sect) .OR. sf(170,sect) .OR. sf(171,sect)) THEN
    ALLOCATE ( adv_inc_qcl(tdims%i_start:tdims%i_end,           &
                           tdims%j_start:tdims%j_end,           &
                                       1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_qcl,qcl_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_qcl(i,j,k) = qcl_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_qcl(1,1,1) )
  END IF

  ! Prognostic cloud ice
  IF (sf(184,sect) .OR. sf(172,sect) .OR. sf(173,sect)) THEN
    ALLOCATE ( adv_inc_qcf(tdims%i_start:tdims%i_end,           &
                           tdims%j_start:tdims%j_end,           &
                                       1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_qcf,qcf_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_qcf(i,j,k) = qcf_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_qcf(1,1,1) )
  END IF

  ! U wind
  IF ( sf(185,sect) ) THEN
    ALLOCATE ( adv_inc_u(udims%i_start:udims%i_end,               &
                         udims%j_start:udims%j_end,               &
                         udims%k_start:udims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(udims,adv_inc_u,r_u)
    DO k= udims%k_start, udims%k_end
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          adv_inc_u(i,j,k) = r_u(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_u( 1,1,1) )
  END IF

  ! V wind
  IF ( sf(186,sect) ) THEN
    ALLOCATE ( adv_inc_v(vdims%i_start:vdims%i_end,               &
                         vdims%j_start:vdims%j_end,               &
                         vdims%k_start:vdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(vdims,adv_inc_v,r_v)
    DO k= vdims%k_start, vdims%k_end
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          adv_inc_v(i,j,k) = r_v(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_v( 1,1,1) )
  END IF

  ! Extra prognostic rain
  IF (l_mcr_qrain .AND. sf(189,sect) ) THEN
    ALLOCATE ( adv_inc_qrain(tdims%i_start:tdims%i_end,         &
                             tdims%j_start:tdims%j_end,         &
                                         1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_qrain,qrain_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_qrain(i,j,k) = qrain_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_qrain(1,1,1) )
  END IF

  ! Extra prognostic graupel
  IF (l_mcr_qgraup .AND. sf(190,sect) ) THEN
    ALLOCATE ( adv_inc_qgraup(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                                          1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_qgraup,qgraup_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_qgraup(i,j,k) = qgraup_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_qgraup(1,1,1) )
  END IF

  ! Extra prognostic cloud ice
  IF (l_mcr_qcf2 .AND. sf(191,sect) ) THEN
    ALLOCATE ( adv_inc_qcf2(tdims%i_start:tdims%i_end,          &
                            tdims%j_start:tdims%j_end,          &
                                        1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_qcf2,qcf2_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_qcf2(i,j,k) = qcf2_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_qcf2(1,1,1) )
  END IF

  ! If PC2 and prognostic bulk cloud fraction
  IF (i_cld_vn == i_cld_pc2 .AND. sf(192,sect) ) THEN
    ALLOCATE ( adv_inc_cf(tdims%i_start:tdims%i_end,            &
                          tdims%j_start:tdims%j_end,            &
                                      1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_cf,cf_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_cf(i,j,k) = cf_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_cf(1,1,1) )
  END IF

  ! If PC2 and prognostic liquid cloud fraction
  IF (i_cld_vn == i_cld_pc2 .AND.                               &
      (sf(193,sect) .OR. sf(176,sect) .OR. sf(177,sect))) THEN
    ALLOCATE ( adv_inc_cfl(tdims%i_start:tdims%i_end,           &
                           tdims%j_start:tdims%j_end,           &
                                       1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_cfl,cfl_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_cfl(i,j,k) = cfl_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_cfl(1,1,1) )
  END IF

  ! If PC2 and prognostic frozen cloud fraction
  IF (i_cld_vn == i_cld_pc2 .AND.                               &
      (sf(194,sect) .OR. sf(178,sect) .OR. sf(179,sect))) THEN
    ALLOCATE ( adv_inc_cff(tdims%i_start:tdims%i_end,           &
                           tdims%j_start:tdims%j_end,           &
                                       1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_cff,cff_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_cff(i,j,k) = cff_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_cff(1,1,1) )
  END IF

  ! Prognostic mixing ratio water vapour
  IF (sf(195,sect)) THEN
    ALLOCATE ( adv_inc_m_v (tdims%i_start:tdims%i_end,          &
                            tdims%j_start:tdims%j_end,          &
                                        1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_m_v,m_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_m_v(i,j,k)  = m_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_m_v(1,1,1) )
  END IF

  ! Prognostic mixing ratio cloud liquid water
  IF (sf(196,sect)) THEN
    ALLOCATE ( adv_inc_m_cl(tdims%i_start:tdims%i_end,          &
                            tdims%j_start:tdims%j_end,          &
                                        1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_m_cl,mcl_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_m_cl(i,j,k) = mcl_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_m_cl(1,1,1) )
  END IF

  ! Prognostic mixing ratio cloud ice
  IF (sf(197,sect)) THEN
    ALLOCATE ( adv_inc_m_cf(tdims%i_start:tdims%i_end,          &
                            tdims%j_start:tdims%j_end,          &
                                        1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_m_cf,mcf_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_m_cf(i,j,k) = mcf_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_m_cf(1,1,1) )
  END IF

  ! Extra prognostic mixing ratio rain
  IF (l_mcr_qrain .AND. sf(198,sect) ) THEN
    ALLOCATE ( adv_inc_m_r(tdims%i_start:tdims%i_end,           &
                           tdims%j_start:tdims%j_end,           &
                                       1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_m_r,mrain_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_m_r(i,j,k) = mrain_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_m_r(1,1,1) )
  END IF

  ! Extra prognostic mixing ratio graupel
  IF (l_mcr_qgraup .AND. sf(199,sect) ) THEN
    ALLOCATE ( adv_inc_m_gr(tdims%i_start:tdims%i_end,          &
                            tdims%j_start:tdims%j_end,          &
                                        1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_m_gr,mgraup_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_m_gr(i,j,k) = mgraup_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_m_gr(1,1,1) )
  END IF

  ! Extra prognostic mixing ratio cloud ice
  IF (l_mcr_qcf2 .AND. sf(200,sect) ) THEN
    ALLOCATE ( adv_inc_m_cf2(tdims%i_start:tdims%i_end,         &
                             tdims%j_start:tdims%j_end,         &
                                         1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_m_cf2,mcf2_star)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          adv_inc_m_cf2(i,j,k) = mcf2_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( adv_inc_m_cf2(1,1,1) )
  END IF

END SUBROUTINE adv_incs_init

SUBROUTINE adv_incs_dealloc()
  IMPLICIT NONE
  IF (ALLOCATED(adv_inc_t)) DEALLOCATE (adv_inc_t)
  IF (ALLOCATED(adv_inc_u)) DEALLOCATE (adv_inc_u)
  IF (ALLOCATED(adv_inc_v)) DEALLOCATE (adv_inc_v)
  IF (ALLOCATED(adv_inc_q)) DEALLOCATE (adv_inc_q)
  IF (ALLOCATED(adv_inc_qcl)) DEALLOCATE (adv_inc_qcl)
  IF (ALLOCATED(adv_inc_qcf)) DEALLOCATE (adv_inc_qcf)
  IF (ALLOCATED(adv_inc_qrain)) DEALLOCATE (adv_inc_qrain)
  IF (ALLOCATED(adv_inc_qgraup)) DEALLOCATE (adv_inc_qgraup)
  IF (ALLOCATED(adv_inc_qcf2)) DEALLOCATE (adv_inc_qcf2)
  IF (ALLOCATED(adv_inc_m_v)) DEALLOCATE (adv_inc_m_v)
  IF (ALLOCATED(adv_inc_m_cl)) DEALLOCATE (adv_inc_m_cl)
  IF (ALLOCATED(adv_inc_m_cf)) DEALLOCATE (adv_inc_m_cf)
  IF (ALLOCATED(adv_inc_m_r)) DEALLOCATE (adv_inc_m_r)
  IF (ALLOCATED(adv_inc_m_gr)) DEALLOCATE (adv_inc_m_gr)
  IF (ALLOCATED(adv_inc_m_cf2)) DEALLOCATE (adv_inc_m_cf2)
  IF (ALLOCATED(adv_inc_cf)) DEALLOCATE (adv_inc_cf)
  IF (ALLOCATED(adv_inc_cfl)) DEALLOCATE (adv_inc_cfl)
  IF (ALLOCATED(adv_inc_cff)) DEALLOCATE (adv_inc_cff)
END SUBROUTINE adv_incs_dealloc

END MODULE adv_increments_mod

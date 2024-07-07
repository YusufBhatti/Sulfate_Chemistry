! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to deal with the calculation of EOT increment diagnostics
! Declares variables and allocatable arrays which are used in the calculation
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Climate Diagnostics
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to the UM programming standards in UMDP3

MODULE eot_increments_mod

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
USE cloud_inputs_mod, ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2

USE atm_fields_mod, ONLY: u, v, w, q, qcf, qcl, qcf2, qrain, qgraup ,     &
    exner_theta_levels, theta, cf_bulk, cf_liquid, cf_frozen, m_v, m_cl,  &
    m_cf, m_cf2, m_gr, m_r, thetav, dryrho

USE atm_fields_bounds_mod, ONLY: pdims, tdims, udims, vdims, wdims

USE stash_array_mod, ONLY: sf

IMPLICIT NONE

! STASH section for diagnostics
INTEGER, PARAMETER :: sect = 30

! Variables for increment diagnostics
REAL, ALLOCATABLE::                                                           &
    eot_inc_t(:,:,:), eot_inc_rho(:,:,:),                                     &
    eot_inc_theta(:,:,:), eot_inc_thetav(:,:,:), eot_inc_dryrho(:,:,:),       &
    eot_inc_u(:,:,:), eot_inc_v(:,:,:), eot_inc_w(:,:,:),                     &
    eot_inc_q(:,:,:), eot_inc_qcl(:,:,:), eot_inc_qcf(:,:,:),                 &
    eot_inc_qrain(:,:,:), eot_inc_qgraup(:,:,:), eot_inc_qcf2(:,:,:),         &
    eot_inc_m_v(:,:,:), eot_inc_m_cl(:,:,:), eot_inc_m_cf(:,:,:),             &
    eot_inc_m_r(:,:,:), eot_inc_m_gr(:,:,:), eot_inc_m_cf2(:,:,:),            &
    eot_inc_cf(:,:,:), eot_inc_cfl(:,:,:), eot_inc_cff(:,:,:)

! Loop counters
INTEGER :: i, j, k

CONTAINS

SUBROUTINE eot_incs_init()

IMPLICIT NONE

!---------------------------------------------------------------------------
! Allocate and initialise arrays for total timestep increments
! NOTE it is very important this is always done before any process changes
!      the fields at the beginning of the timestep.
!---------------------------------------------------------------------------

! Temperature
IF ( sf(181,sect) ) THEN
  ALLOCATE ( eot_inc_t(tdims%i_start:tdims%i_end,               &
                       tdims%j_start:tdims%j_end,               &
                                   1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_t,theta,exner_theta_levels)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_t(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_t( 1,1,1) )
END IF

! Prognostic water vapour
IF (sf(182,sect)) THEN
  ALLOCATE ( eot_inc_q  (tdims%i_start:tdims%i_end,           &
                         tdims%j_start:tdims%j_end,           &
                                     1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_q,q)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_q(i,j,k) = q(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_q(1,1,1) )
END IF

! Prognostic cloud liquid water
IF (sf(183,sect)) THEN
  ALLOCATE ( eot_inc_qcl(tdims%i_start:tdims%i_end,           &
                         tdims%j_start:tdims%j_end,           &
                                     1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_qcl,qcl)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_qcl(i,j,k) = qcl(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_qcl(1,1,1) )
END IF

! Prognostic cloud ice
IF (sf(184,sect)) THEN
  ALLOCATE ( eot_inc_qcf(tdims%i_start:tdims%i_end,           &
                         tdims%j_start:tdims%j_end,           &
                                     1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_qcf,qcf)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_qcf(i,j,k) = qcf(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_qcf(1,1,1) )
END IF

! U or U**2
IF ( sf(175,sect) .OR. sf(185,sect) ) THEN
  ALLOCATE ( eot_inc_u(udims%i_start:udims%i_end,               &
                       udims%j_start:udims%j_end,               &
                       udims%k_start:udims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(udims,eot_inc_u,u)
  DO k= udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        eot_inc_u(i,j,k) = u(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_u( 1,1,1) )
END IF

! V or V**2
IF ( sf(176,sect) .OR. sf(186,sect) ) THEN
  ALLOCATE ( eot_inc_v(vdims%i_start:vdims%i_end,               &
                       vdims%j_start:vdims%j_end,               &
                       vdims%k_start:vdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(vdims,eot_inc_v,v)
  DO k= vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        eot_inc_v(i,j,k) = v(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_v( 1,1,1) )
END IF

! W or W**2
IF ( sf(177,sect) .OR. sf(187,sect) ) THEN
  ALLOCATE ( eot_inc_w(wdims%i_start:wdims%i_end,               &
                       wdims%j_start:wdims%j_end,               &
                                   1:wdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(wdims,eot_inc_w,w)
  DO k= 1, wdims%k_end
    DO j = wdims%j_start, wdims%j_end
      DO i = wdims%i_start, wdims%i_end
        eot_inc_w(i,j,k) = w(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_w( 1,1,1) )
END IF

! Rho increments can only be allocated here, not initialised as
! this is done later in atm_step_4A.F90
IF ( sf(188,sect) ) THEN
  ALLOCATE ( eot_inc_rho(pdims%i_start:pdims%i_end,             &
                         pdims%j_start:pdims%j_end,             &
                         pdims%k_start:pdims%k_end) )
ELSE
  ALLOCATE ( eot_inc_rho( 1,1,1) )
END IF

! Extra prognostic rain
IF (l_mcr_qrain .AND. sf(189,sect) ) THEN
  ALLOCATE ( eot_inc_qrain(tdims%i_start:tdims%i_end,         &
                           tdims%j_start:tdims%j_end,         &
                                       1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_qrain,qrain)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_qrain(i,j,k) = qrain(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_qrain(1,1,1) )
END IF

! Extra prognostic graupel
IF (l_mcr_qgraup .AND. sf(190,sect) ) THEN
  ALLOCATE ( eot_inc_qgraup(tdims%i_start:tdims%i_end,        &
                            tdims%j_start:tdims%j_end,        &
                                        1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_qgraup,qgraup)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_qgraup(i,j,k) = qgraup(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_qgraup(1,1,1) )
END IF

! Extra prognostic cloud ice
IF (l_mcr_qcf2 .AND. sf(191,sect) ) THEN
  ALLOCATE ( eot_inc_qcf2(tdims%i_start:tdims%i_end,          &
                          tdims%j_start:tdims%j_end,          &
                                      1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_qcf2,qcf2)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_qcf2(i,j,k) = qcf2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_qcf2(1,1,1) )
END IF

! If PC2 and prognostic bulk cloud fraction
IF (i_cld_vn == i_cld_pc2 .AND. sf(192,sect) ) THEN
  ALLOCATE ( eot_inc_cf(tdims%i_start:tdims%i_end,            &
                        tdims%j_start:tdims%j_end,            &
                                    1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_cf,cf_bulk)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_cf(i,j,k) = cf_bulk(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_cf(1,1,1) )
END IF

! If PC2 and prognostic liquid cloud fraction
IF (i_cld_vn == i_cld_pc2 .AND. sf(193,sect) ) THEN
  ALLOCATE ( eot_inc_cfl(tdims%i_start:tdims%i_end,           &
                         tdims%j_start:tdims%j_end,           &
                                     1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_cfl,cf_liquid)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_cfl(i,j,k) = cf_liquid(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_cfl(1,1,1) )
END IF

! If PC2 and prognostic frozen cloud fraction
IF (i_cld_vn == i_cld_pc2 .AND. sf(194,sect) ) THEN
  ALLOCATE ( eot_inc_cff(tdims%i_start:tdims%i_end,           &
                         tdims%j_start:tdims%j_end,           &
                                     1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_cff,cf_frozen)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_cff(i,j,k) = cf_frozen(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_cff(1,1,1) )
END IF

! Prognostic mixing ratio water vapour
IF (sf(195,sect)) THEN
  ALLOCATE ( eot_inc_m_v (tdims%i_start:tdims%i_end,          &
                          tdims%j_start:tdims%j_end,          &
                                      1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_v,m_v)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_v(i,j,k)  = m_v(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_m_v(1,1,1) )
END IF

! Prognostic mixing ratio cloud liquid water
IF (sf(196,sect)) THEN
  ALLOCATE ( eot_inc_m_cl(tdims%i_start:tdims%i_end,          &
                          tdims%j_start:tdims%j_end,          &
                                      1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_cl,m_cl)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_cl(i,j,k) = m_cl(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_m_cl(1,1,1) )
END IF

! Prognostic mixing ratio cloud ice
IF (sf(197,sect)) THEN
  ALLOCATE ( eot_inc_m_cf(tdims%i_start:tdims%i_end,          &
                          tdims%j_start:tdims%j_end,          &
                                      1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_cf,m_cf)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_cf(i,j,k) = m_cf(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_m_cf(1,1,1) )
END IF

! Extra prognostic mixing ratio rain
IF (l_mcr_qrain .AND. sf(198,sect) ) THEN
  ALLOCATE ( eot_inc_m_r(tdims%i_start:tdims%i_end,           &
                         tdims%j_start:tdims%j_end,           &
                                     1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_r,m_r)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_r(i,j,k) = m_r(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_m_r(1,1,1) )
END IF

! Extra prognostic mixing ratio graupel
IF (l_mcr_qgraup .AND. sf(199,sect) ) THEN
  ALLOCATE ( eot_inc_m_gr(tdims%i_start:tdims%i_end,          &
                          tdims%j_start:tdims%j_end,          &
                                      1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_gr,m_gr)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_gr(i,j,k) = m_gr(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_m_gr(1,1,1) )
END IF

! Extra prognostic mixing ratio cloud ice
IF (l_mcr_qcf2 .AND. sf(200,sect) ) THEN
  ALLOCATE ( eot_inc_m_cf2(tdims%i_start:tdims%i_end,         &
                           tdims%j_start:tdims%j_end,         &
                                       1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_cf2,m_cf2)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_cf2(i,j,k) = m_cf2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_m_cf2(1,1,1) )
END IF

! Theta
IF ( sf(901,sect) ) THEN
  ALLOCATE ( eot_inc_theta(tdims%i_start:tdims%i_end,               &
                           tdims%j_start:tdims%j_end,               &
                                       1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_theta,theta)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_theta(i,j,k) = theta(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_theta( 1,1,1) )
END IF

! ThetaV
IF ( sf(902,sect) ) THEN
  ALLOCATE ( eot_inc_thetav(tdims%i_start:tdims%i_end,               &
                            tdims%j_start:tdims%j_end,               &
                                        1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_thetav,thetav)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_thetav(i,j,k) = thetav(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_thetav( 1,1,1) )
END IF

! Dry Rho
IF ( sf(903,sect) ) THEN
  ALLOCATE ( eot_inc_dryrho(pdims%i_start:pdims%i_end,               &
                            pdims%j_start:pdims%j_end,               &
                            pdims%k_start:pdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(pdims,eot_inc_dryrho,dryrho)
  DO k= pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        eot_inc_dryrho(i,j,k) = dryrho(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( eot_inc_dryrho( 1,1,1) )
END IF

END SUBROUTINE eot_incs_init

SUBROUTINE eot_incs_calc()

IMPLICIT NONE

!---------------------------------------------------------------------------
! Calculation of total increments
! NOTE - this must be after all processes which alter model prognostics
!---------------------------------------------------------------------------

! Temperature from prognostic theta
IF (sf(181,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_t,theta,exner_theta_levels)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_t(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)      &
              - eot_inc_t(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Prognostic water vapour
IF (sf(182,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_q,q)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_q(i,j,k)   = q(i,j,k) - eot_inc_q(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Prognostic cloud liquid water
IF (sf(183,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_qcl,qcl)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_qcl(i,j,k) = qcl(i,j,k) - eot_inc_qcl(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Prognostic cloud ice
IF (sf(184,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_qcf,qcf)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_qcf(i,j,k) = qcf(i,j,k) - eot_inc_qcf(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! U or U**2
IF (sf(175,sect) .OR. sf(185,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(udims,eot_inc_u,u)
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        eot_inc_u(i,j,k) = u(i,j,k) - eot_inc_u(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! V or V**2
IF (sf(176,sect) .OR. sf(186,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(vdims,eot_inc_v,v)
  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        eot_inc_v(i,j,k) = v(i,j,k) - eot_inc_v(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! W or W**2
IF (sf(177,sect) .OR. sf(187,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(wdims,eot_inc_w,w)
  DO k =             1, wdims%k_end
    DO j = wdims%j_start, wdims%j_end
      DO i = wdims%i_start, wdims%i_end
        eot_inc_w(i,j,k) = w(i,j,k) - eot_inc_w(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Extra prognostic rain
IF ( l_mcr_qrain .AND. sf(189,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_qrain,qrain)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_qrain(i,j,k) = qrain(i,j,k) - eot_inc_qrain(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Extra prognostic graupel
IF ( l_mcr_qgraup .AND. sf(190,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_qgraup,qgraup)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_qgraup(i,j,k) = qgraup(i,j,k) - eot_inc_qgraup(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Extra prognostic cloud ice
IF ( l_mcr_qcf2 .AND. sf(191,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_qcf2,qcf2)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_qcf2(i,j,k) = qcf2(i,j,k) - eot_inc_qcf2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! If PC2 and prognostic cf
IF ( i_cld_vn == i_cld_pc2 .AND. sf(192,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_cf,cf_bulk)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_cf(i,j,k) = cf_bulk(i,j,k) - eot_inc_cf(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! If PC2 and prognostic cfl
IF ( i_cld_vn == i_cld_pc2 .AND. sf(193,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_cfl,cf_liquid)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_cfl(i,j,k) = cf_liquid(i,j,k) - eot_inc_cfl(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! If PC2 and prognostic cff
IF ( i_cld_vn == i_cld_pc2 .AND. sf(194,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_cff,cf_frozen)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_cff(i,j,k) = cf_frozen(i,j,k) - eot_inc_cff(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Prognostic mixing ratio water vapour
IF (sf(195,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_v,m_v)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_v(i,j,k)  = m_v(i,j,k) - eot_inc_m_v(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Prognostic mixing ratio cloud water
IF (sf(196,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_cl,m_cl)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_cl(i,j,k) = m_cl(i,j,k) - eot_inc_m_cl(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Prognostic mixing ratio cloud ice
IF (sf(197,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_cf,m_cf)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_cf(i,j,k) = m_cf(i,j,k) - eot_inc_m_cf(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Extra prognostic mixing ratio rain
IF ( l_mcr_qrain .AND. sf(198,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_r,m_r)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_r(i,j,k) = m_r(i,j,k) - eot_inc_m_r(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Extra prognostic mixing ratio graupel
IF ( l_mcr_qgraup .AND. sf(199,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_gr,m_gr)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_gr(i,j,k) = m_gr(i,j,k) - eot_inc_m_gr(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Extra prognostic mixing ratio cloud ice
IF ( l_mcr_qcf2 .AND. sf(200,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_m_cf2,m_cf2)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_m_cf2(i,j,k) = m_cf2(i,j,k) - eot_inc_m_cf2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Theta
IF (sf(901,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_theta,theta)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_theta(i,j,k) = theta(i,j,k) - eot_inc_theta(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! ThetaV
IF (sf(902,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,eot_inc_thetav,thetav)
  DO k= 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        eot_inc_thetav(i,j,k) = thetav(i,j,k) - eot_inc_thetav(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Dry Rho
IF (sf(903,sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(pdims,eot_inc_dryrho,dryrho)
  DO k= pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        eot_inc_dryrho(i,j,k) = dryrho(i,j,k) - eot_inc_dryrho(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

END SUBROUTINE eot_incs_calc

SUBROUTINE eot_incs_dealloc()
IMPLICIT NONE
IF (ALLOCATED(eot_inc_dryrho)) DEALLOCATE (eot_inc_dryrho)
IF (ALLOCATED(eot_inc_thetav)) DEALLOCATE (eot_inc_thetav)
IF (ALLOCATED(eot_inc_theta)) DEALLOCATE (eot_inc_theta)
IF (ALLOCATED(eot_inc_m_cf2)) DEALLOCATE (eot_inc_m_cf2)
IF (ALLOCATED(eot_inc_m_gr)) DEALLOCATE (eot_inc_m_gr)
IF (ALLOCATED(eot_inc_m_r)) DEALLOCATE (eot_inc_m_r)
IF (ALLOCATED(eot_inc_m_cf)) DEALLOCATE (eot_inc_m_cf)
IF (ALLOCATED(eot_inc_m_cl)) DEALLOCATE (eot_inc_m_cl)
IF (ALLOCATED(eot_inc_m_v)) DEALLOCATE (eot_inc_m_v)
IF (ALLOCATED(eot_inc_cff)) DEALLOCATE (eot_inc_cff)
IF (ALLOCATED(eot_inc_cfl)) DEALLOCATE (eot_inc_cfl)
IF (ALLOCATED(eot_inc_cf)) DEALLOCATE (eot_inc_cf)
IF (ALLOCATED(eot_inc_qcf2)) DEALLOCATE (eot_inc_qcf2)
IF (ALLOCATED(eot_inc_qgraup)) DEALLOCATE (eot_inc_qgraup)
IF (ALLOCATED(eot_inc_qrain)) DEALLOCATE (eot_inc_qrain)
IF (ALLOCATED(eot_inc_rho)) DEALLOCATE (eot_inc_rho)
IF (ALLOCATED(eot_inc_w)) DEALLOCATE (eot_inc_w)
IF (ALLOCATED(eot_inc_v)) DEALLOCATE (eot_inc_v)
IF (ALLOCATED(eot_inc_u)) DEALLOCATE (eot_inc_u)
IF (ALLOCATED(eot_inc_qcf)) DEALLOCATE (eot_inc_qcf)
IF (ALLOCATED(eot_inc_qcl)) DEALLOCATE (eot_inc_qcl)
IF (ALLOCATED(eot_inc_q)) DEALLOCATE (eot_inc_q)
IF (ALLOCATED(eot_inc_t)) DEALLOCATE (eot_inc_t)
END SUBROUTINE eot_incs_dealloc

END MODULE eot_increments_mod

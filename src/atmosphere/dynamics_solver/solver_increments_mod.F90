! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to deal with the calculation of solver increment diagnostics
! Declares variables and allocatable arrays which are used in the calculation
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Dynamics Solver
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to the UM programming standards in UMDP3

MODULE solver_increments_mod

USE atm_fields_bounds_mod, ONLY: pdims_s, tdims_s, udims_s, vdims_s, wdims_s, &
                                 pdims, tdims, udims, vdims, wdims

USE stash_array_mod, ONLY: sf

IMPLICIT NONE

! STASH section for diagnostics
INTEGER, PARAMETER :: sect = 10

! Variables for increment diagnostics
REAL, ALLOCATABLE:: solv_inc_t(:,:,:), solv_inc_u(:,:,:),                     &
                    solv_inc_v(:,:,:), solv_inc_w(:,:,:)

! Loop counters
INTEGER :: i, j, k

CONTAINS

SUBROUTINE solv_incs_init(exner_theta, theta, r_u, r_v, r_w, u, v, w)

  IMPLICIT NONE

  ! Subroutine arguments

  REAL, INTENT(IN) :: exner_theta (tdims_s%i_start:tdims_s%i_end,             &
                                   tdims_s%j_start:tdims_s%j_end,             &
                                   tdims_s%k_start:tdims_s%k_end)
  REAL, INTENT(IN) :: theta       (tdims%i_start:tdims%i_end,                 &
                                   tdims%j_start:tdims%j_end,                 &
                                   tdims%k_start:tdims%k_end)

  ! Latest wind increments prior to solver
  REAL, INTENT(IN) :: r_u         (udims_s%i_start:udims_s%i_end,             &
                                   udims_s%j_start:udims_s%j_end,             &
                                   udims_s%k_start:udims_s%k_end)
  REAL, INTENT(IN) :: r_v         (vdims_s%i_start:vdims_s%i_end,             &
                                   vdims_s%j_start:vdims_s%j_end,             &
                                   vdims_s%k_start:vdims_s%k_end)
  REAL, INTENT(IN) :: r_w         (wdims%i_start:wdims%i_end,                 &
                                   wdims%j_start:wdims%j_end,                 &
                                   wdims%k_start:wdims%k_end)

  ! Beginning of timestep winds
  REAL, INTENT(IN) :: u           (udims_s%i_start:udims_s%i_end,             &
                                   udims_s%j_start:udims_s%j_end,             &
                                   udims_s%k_start:udims_s%k_end)
  REAL, INTENT(IN) :: v           (vdims_s%i_start:vdims_s%i_end,             &
                                   vdims_s%j_start:vdims_s%j_end,             &
                                   vdims_s%k_start:vdims_s%k_end)
  REAL, INTENT(IN) :: w           (wdims_s%i_start:wdims_s%i_end,             &
                                   wdims_s%j_start:wdims_s%j_end,             &
                                   wdims_s%k_start:wdims_s%k_end)

  !---------------------------------------------------------------------------
  ! Allocate and initialise arrays for solver increments
  !---------------------------------------------------------------------------

  ! Temperature
  IF ( sf(181,sect) ) THEN
    ALLOCATE ( solv_inc_t(tdims%i_start:tdims%i_end,             &
                          tdims%j_start:tdims%j_end,             &
                                      1:tdims%k_end))

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,solv_inc_t,theta,exner_theta)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          solv_inc_t(i,j,k) = exner_theta(i,j,k) * theta(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( solv_inc_t( 1,1,1) )
  END IF

  ! U Wind
  IF ( sf(185,sect) ) THEN
    ALLOCATE ( solv_inc_u(udims%i_start:udims%i_end,             &
                          udims%j_start:udims%j_end,             &
                          udims%k_start:udims%k_end))

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(udims,solv_inc_u,u,r_u)
    DO k= udims%k_start, udims%k_end
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          solv_inc_u(i,j,k) = u(i,j,k) + r_u(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( solv_inc_u( 1,1,1) )
  END IF

  ! V wind
  IF ( sf(186,sect) ) THEN
    ALLOCATE ( solv_inc_v(vdims%i_start:vdims%i_end,             &
                          vdims%j_start:vdims%j_end,             &
                          vdims%k_start:vdims%k_end))

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(vdims,solv_inc_v,v,r_v)
    DO k= vdims%k_start, vdims%k_end
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          solv_inc_v(i,j,k) = v(i,j,k) + r_v(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( solv_inc_v( 1,1,1) )
  END IF

  ! W wind
  IF ( sf(187,sect) ) THEN
    ALLOCATE ( solv_inc_w(wdims%i_start:wdims%i_end,             &
                          wdims%j_start:wdims%j_end,             &
                                      1:wdims%k_end))

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(wdims,solv_inc_w,w,r_w)
    DO k=              1, wdims%k_end
      DO j = wdims%j_start, wdims%j_end
        DO i = wdims%i_start, wdims%i_end
          solv_inc_w(i,j,k) = w(i,j,k) + r_w(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE ( solv_inc_w( 1,1,1) )
  END IF

END SUBROUTINE solv_incs_init

SUBROUTINE solv_incs_calc(exner, thetav, m_v, u, v, w)

  USE planet_constants_mod, ONLY: recip_epsilon
  USE horiz_grid_mod, ONLY: intw_rho2w

  IMPLICIT NONE

  ! Subroutine arguments

  REAL, INTENT(IN) :: exner      (pdims_s%i_start:pdims_s%i_end,             &
                                  pdims_s%j_start:pdims_s%j_end,             &
                                  pdims_s%k_start:pdims_s%k_end+1)
  REAL, INTENT(IN) :: thetav     (tdims_s%i_start:tdims_s%i_end,             &
                                  tdims_s%j_start:tdims_s%j_end,             &
                                  tdims_s%k_start:tdims_s%k_end)
  REAL, INTENT(IN) :: m_v        (tdims_s%i_start:tdims_s%i_end,             &
                                  tdims_s%j_start:tdims_s%j_end,             &
                                  tdims_s%k_start:tdims_s%k_end)

  REAL, INTENT(IN) :: u          (udims_s%i_start:udims_s%i_end,             &
                                  udims_s%j_start:udims_s%j_end,             &
                                  udims_s%k_start:udims_s%k_end)
  REAL, INTENT(IN) :: v          (vdims_s%i_start:vdims_s%i_end,             &
                                  vdims_s%j_start:vdims_s%j_end,             &
                                  vdims_s%k_start:vdims_s%k_end)
  REAL, INTENT(IN) :: w          (wdims_s%i_start:wdims_s%i_end,             &
                                  wdims_s%j_start:wdims_s%j_end,             &
                                  wdims_s%k_start:wdims_s%k_end)

  REAL :: temp, exner_theta

  !---------------------------------------------------------------------------
  ! Calculate solver increments
  !---------------------------------------------------------------------------

  ! Temperature
  IF ( sf(181,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i,temp, exner_theta)                                        &
!$OMP SHARED(tdims,solv_inc_t,thetav,m_v,recip_epsilon,intw_rho2w,exner)
    DO k=              1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          exner_theta = (intw_rho2w(k,2)*exner(i,j,k)+      &
                          intw_rho2w(k,1)*exner(i,j,k+1))
          temp = exner_theta * thetav(i,j,k) / &
                  (1.0 + m_v(i,j,k)*recip_epsilon)
          solv_inc_t(i,j,k) = temp - solv_inc_t(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! U Wind
  IF ( sf(185,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(udims,solv_inc_u,u)
    DO k= udims%k_start, udims%k_end
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          solv_inc_u(i,j,k) = u(i,j,k) - solv_inc_u(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! V wind
  IF ( sf(186,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(vdims,solv_inc_v,v)
    DO k= vdims%k_start, vdims%k_end
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          solv_inc_v(i,j,k) = v(i,j,k) - solv_inc_v(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! W wind
  IF ( sf(187,sect) ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(wdims,solv_inc_w,w)
    DO k=              1, wdims%k_end
      DO j = wdims%j_start, wdims%j_end
        DO i = wdims%i_start, wdims%i_end
          solv_inc_w(i,j,k) = w(i,j,k) - solv_inc_w(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

END SUBROUTINE solv_incs_calc

SUBROUTINE solv_incs_dealloc()
  IMPLICIT NONE
  IF (ALLOCATED(solv_inc_u)) DEALLOCATE ( solv_inc_u )
  IF (ALLOCATED(solv_inc_v)) DEALLOCATE ( solv_inc_v )
  IF (ALLOCATED(solv_inc_w)) DEALLOCATE ( solv_inc_w )
  IF (ALLOCATED(solv_inc_t)) DEALLOCATE ( solv_inc_t )
END SUBROUTINE solv_incs_dealloc

END MODULE solver_increments_mod

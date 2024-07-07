! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to deal with the calculation of diffusion increment diagnostics
! Declares variables and allocatable arrays which are used in the calculation
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to the UM programming standards in UMDP3

MODULE diff_increments_mod

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
USE cloud_inputs_mod, ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2

USE atm_fields_bounds_mod, ONLY: tdims_s, udims_s, vdims_s, &
                                 tdims, udims, vdims, wdims
USE nlsizes_namelist_mod, ONLY: row_length, rows

USE stash_array_mod, ONLY: sf

IMPLICIT NONE

! STASH section for diagnostics
INTEGER, PARAMETER :: sect = 13

! Variables for increment diagnostics
REAL, ALLOCATABLE::                                                          &
    diff_inc_u(:,:,:), diff_inc_v(:,:,:), diff_inc_w(:,:,:),                 &
    diff_inc_t(:,:,:), diff_inc_q(:,:,:), w_local_mask(:,:)

! Loop counters
INTEGER :: i, j, k

CONTAINS

SUBROUTINE diff_incs_init(s_u, s_v, s_w, s_theta, s_q)

  USE atm_step_local, ONLY: cycleno
  USE dynamics_input_mod, ONLY: numcycles
  
  IMPLICIT NONE

  ! Subroutine arguments

  REAL, TARGET     :: s_u        (udims_s%i_start:udims_s%i_end,             &
                                  udims_s%j_start:udims_s%j_end,             &
                                  udims_s%k_start:udims_s%k_end)
  REAL, TARGET     :: s_v        (vdims_s%i_start:vdims_s%i_end,             &
                                  vdims_s%j_start:vdims_s%j_end,             &
                                  vdims_s%k_start:vdims_s%k_end)
  REAL, TARGET     :: s_w        (wdims  %i_start:wdims  %i_end,             &
                                  wdims  %j_start:wdims  %j_end,             &
                                  wdims  %k_start:wdims  %k_end)

  REAL, INTENT(IN) :: s_theta    (tdims_s%i_start:tdims_s%i_end,             &
                                  tdims_s%j_start:tdims_s%j_end,             &
                                  tdims_s%k_start:tdims_s%k_end)

  ! Water Vapour
  REAL, INTENT(IN) :: s_q        (tdims  %i_start:tdims  %i_end,             &
                                  tdims  %j_start:tdims  %j_end,             &
                                  tdims  %k_start:tdims  %k_end)

  !---------------------------------------------------------------------------
  ! Allocate and initialise arrays for advection increments
  !---------------------------------------------------------------------------

  ! Apply diagnostics only at last cycle.
  IF ( cycleno == numcycles ) THEN

    ! Temperature
    ! Note that we are holding theta, conversion to temperature will be done
    !  in calculation of increment.
    IF ( sf(181,sect) ) THEN
      ALLOCATE ( diff_inc_t(tdims%i_start:tdims%i_end,            &
                            tdims%j_start:tdims%j_end,            &
                                        1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,diff_inc_t,s_theta)
      DO k=              1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            diff_inc_t(i,j,k) = s_theta(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    ELSE
      ALLOCATE ( diff_inc_t( 1,1,1) )
    END IF

    ! Prognostic water vapour
    IF (sf(182,sect)) THEN
      ALLOCATE ( diff_inc_q  (tdims%i_start:tdims%i_end,          &
                              tdims%j_start:tdims%j_end,          &
                                          1:tdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,diff_inc_q,s_q)
      DO k=              1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            diff_inc_q(i,j,k) = s_q(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    ELSE
      ALLOCATE ( diff_inc_q(1,1,1) )
    END IF

    ! U wind
    IF ( sf(185,sect) ) THEN
      ALLOCATE ( diff_inc_u(udims%i_start:udims%i_end,              &
                            udims%j_start:udims%j_end,              &
                            udims%k_start:udims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(udims,diff_inc_u,s_u)
      DO k= udims%k_start, udims%k_end
        DO j = udims%j_start, udims%j_end
          DO i = udims%i_start, udims%i_end
            diff_inc_u(i,j,k) = s_u(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    ELSE
      ALLOCATE ( diff_inc_u( 1,1,1) )
    END IF

    ! V wind
    IF ( sf(186,sect) ) THEN
      ALLOCATE ( diff_inc_v(vdims%i_start:vdims%i_end,              &
                            vdims%j_start:vdims%j_end,              &
                            vdims%k_start:vdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(vdims,diff_inc_v,s_v)
      DO k= vdims%k_start, vdims%k_end
        DO j = vdims%j_start, vdims%j_end
          DO i = vdims%i_start, vdims%i_end
            diff_inc_v(i,j,k) = s_v(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    ELSE
      ALLOCATE ( diff_inc_v( 1,1,1) )
    END IF

    ! W wind
    IF ( sf(187,sect) ) THEN
      ALLOCATE ( diff_inc_w(wdims%i_start:wdims%i_end,              &
                            wdims%j_start:wdims%j_end,              &
                                        1:wdims%k_end))
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(wdims,diff_inc_w,s_w)
      DO k=            1, wdims%k_end
        DO j = wdims%j_start, wdims%j_end
          DO i = wdims%i_start, wdims%i_end
            diff_inc_w(i,j,k) = s_w(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    ELSE
      ALLOCATE ( diff_inc_w( 1,1,1) )
    END IF

  END IF  ! cycleno == numcycles

  IF ( cycleno == 1 ) THEN

    IF ( sf(201,sect) ) THEN
      ALLOCATE ( w_local_mask(row_length,rows) )
    ELSE
      ALLOCATE ( w_local_mask(1,1) )
    END IF

  END IF  ! cycleno == 1

END SUBROUTINE diff_incs_init

SUBROUTINE diff_incs_dealloc()
  IMPLICIT NONE
  IF (ALLOCATED(diff_inc_t)) DEALLOCATE (diff_inc_t)
  IF (ALLOCATED(diff_inc_q)) DEALLOCATE (diff_inc_q)
  IF (ALLOCATED(diff_inc_u)) DEALLOCATE (diff_inc_u)
  IF (ALLOCATED(diff_inc_v)) DEALLOCATE (diff_inc_v)
  IF (ALLOCATED(diff_inc_w)) DEALLOCATE (diff_inc_w)
  IF (ALLOCATED(w_local_mask)) DEALLOCATE (w_local_mask)
END SUBROUTINE diff_incs_dealloc

END MODULE diff_increments_mod

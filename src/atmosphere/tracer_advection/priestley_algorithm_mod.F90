! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
MODULE priestley_algorithm_mod

USE global_2d_sums_mod, ONLY: global_2d_sums

IMPLICIT NONE

! Description:
!  This module contains routines that apply conservation -using
!  the Priestley scheme- to moisture and tracer fields.
!
!  There are two versions 'original' and 'optimised' of the
!   Priestley algorithm, applied depending on user selection.
!
! Method: ENDGame formulation v4.00 (Feb 2014)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Tracer Advection
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PRIESTLEY_ALGORITHM_MOD'

CONTAINS

!==================================================================
SUBROUTINE Priestley_algorithm(rhon, rhoh, rhol, source, voln, vol,     &
                               its, itf, jts, jtf, kts, ktf,            &
                               irs, irf, jrs, jrf, krs, krf,            &
                               nt                                    )


!----------------------------------------------------------------------------
! This routine imposes the following mass conservation constraint:
!
!  SUM{ rhoh * vol } = SUM{ rhon * voln  + source * vol }  (1)
!
!  The left-hand side of (1) represents the mass at the present time
!  and the right-hand side represents the mass at previous time plus
!  the mass due sources at present time
!
!  The total mass inside this routine is defined as SUM{rho * vol}.
!  This routine should work with any meaning of "rho" and "vol" as long
!  as they are set appropriately outside this routine, while the
!  true mass (which is the quantity to be conserved) remains SUM{rho * vol}
!---------------------------------------------------------------------------

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER, INTENT(IN)  :: its, itf, jts, jtf, kts, ktf    ! Dimension info
INTEGER, INTENT(IN)  :: irs, irf, jrs, jrf, krs, krf, nt

REAL, INTENT(IN)     :: rhon(its:itf,jts:jtf,kts:ktf,nt)
REAL, INTENT(INOUT)  :: rhoh(its:itf,jts:jtf,kts:ktf,nt)
REAL, INTENT(IN)     :: rhol(its:itf,jts:jtf,kts:ktf,nt)
REAL, INTENT(IN)     :: source(its:itf,jts:jtf,kts:ktf,nt)

REAL, INTENT(IN)     :: voln(irs:irf,jrs:jrf,krs:krf)
REAL, INTENT(IN)     :: vol (irs:irf,jrs:jrf,krs:krf)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRIESTLEY_ALGORITHM'

! Local variables

REAL    :: w(its:itf,jts:jtf,kts:ktf,nt), du(its:itf,jts:jtf,kts:ktf,nt)
INTEGER :: flag(its:itf,jts:jtf,kts:ktf,nt)

REAL    :: eps = 1.0e-16    ! Threshold value
REAL    :: a(nt), b(nt), c(nt), c_star(nt)
REAL    :: wa, surplus, temp
REAL    :: local_sum(irs:irf,jrs:jrf,nt,3)
REAL    :: local_sum2(irs:irf,jrs:jrf,2)
REAL    :: mz_sum(nt,3)
REAL    :: mz_sum2(2)
INTEGER :: i,j,k,kk

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!====================================================
! compute some global sums and setup defaults
!====================================================

w(:,:,:,:) = 1.0
mz_sum(:,:) = 0.0

DO kk = 1, nt
  local_sum(:,:,kk,:) = 0.0
  DO k = krs, krf
    DO j = jrs, jrf
      DO i = irs, irf
        local_sum(i,j,kk,1) = local_sum(i,j,kk,1)                     &
                               + rhon(i,j,k,kk) * voln(i,j,k)         &
                               + source(i,j,k,kk) *  vol(i,j,k)
        local_sum(i,j,kk,2) = local_sum(i,j,kk,2)                     &
                               + (rhoh(i,j,k,kk)-rhol(i,j,k,kk))      &
                               * vol(i,j,k)
        local_sum(i,j,kk,3) = local_sum(i,j,kk,3)                     &
                              + rhol(i,j,k,kk) * vol(i,j,k)
      END DO
    END DO
  END DO
END DO

CALL global_2d_sums(local_sum,irf-irs+1,jrf-jrs+1,0,0,3*nt,mz_sum)

c(:) = mz_sum(:,1)
a(:) = mz_sum(:,2)
b(:) = c(:) - mz_sum(:,3)

DO kk = 1, nt
  IF ( a(kk) > b(kk) ) THEN

    c_star(kk) = b(kk)
    DO k = krs, krf
      DO j = jrs, jrf
        DO i = irs, irf
          du(i,j,k,kk) = vol(i,j,k) *(rhoh(i,j,k,kk)-rhol(i,j,k,kk))
        END DO
      END DO
    END DO

  ELSE

    c_star(kk) = -b(kk)
    DO k = krs, krf
      DO j = jrs, jrf
        DO i = irs, irf
          du(i,j,k,kk) = - vol(i,j,k) * (rhoh(i,j,k,kk)-rhol(i,j,k,kk))
        END DO
      END DO
    END DO

  END IF
END DO

DO kk = 1, nt
  DO k = krs, krf
    DO j = jrs, jrf
      DO i = irs, irf

        IF ( du(i,j,k,kk) <= 0.0 ) THEN
          flag(i,j,k,kk) = 1
        ELSE
          w(i,j,k,kk)    = 0.0
          flag(i,j,k,kk) = 0
        END IF
      END DO
    END DO
  END DO
END DO

!====================================================
!  perform a Priestly algorithm on each field kk=1:nt
!====================================================

DO kk = 1, nt

  local_sum2 = 0.0
  DO k = krs, krf
    DO j = jrs, jrf
      DO i = irs, irf
        IF ( flag(i,j,k,kk) == 1 ) THEN
          local_sum2(i,j,1) = local_sum2(i,j,1) + du(i,j,k,kk)
        ELSE
          local_sum2(i,j,2) = local_sum2(i,j,2) + du(i,j,k,kk)
        END IF
      END DO
    END DO
  END DO

  CALL global_2d_sums(local_sum2,irf-irs+1,jrf-jrs+1,0,0,2,mz_sum2)

  surplus = c_star(kk) - mz_sum2(1)
  wa = 1.0
  IF ( ABS(mz_sum2(2)) > eps ) THEN
    wa = surplus / mz_sum2(2)
  END IF

  IF ( wa*(1.0-wa) > 0.0 ) THEN
    DO k = krs, krf
      DO j = jrs, jrf
        DO i = irs, irf
          IF ( flag(i,j,k,kk) == 0 ) THEN
            w(i,j,k,kk) = wa
          END IF
        END DO
      END DO
    END DO

    DO k = krs, krf
      DO j = jrs, jrf
        DO i = irs, irf
          rhoh(i,j,k,kk) =   w(i,j,k,kk)  * rhoh(i,j,k,kk) +           &
                            ( 1.0 - w(i,j,k,kk)) * rhol(i,j,k,kk)
        END DO
      END DO
    END DO

  ELSE

    temp = mz_sum(kk,2) + mz_sum(kk,3)
    IF ( ABS(temp) > eps ) THEN
      temp = mz_sum(kk,1) / temp
      rhoh(:,:,:,kk) = rhoh(:,:,:,kk) * temp
    END IF

  END IF    ! IF wa*(1-Wa) etc


END DO ! number of fields nt

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)


RETURN

END SUBROUTINE priestley_algorithm

!==================================================================

SUBROUTINE priestley_algorithm2 (rhon, rhoh, rhol, source, voln, vol,   &
                                 its, itf, jts, jtf, kts, ktf,          &
                                 irs, irf, jrs, jrf, krs, krf,          &
                                 nt                                   )

!----------------------------------------------------------------------------
! This routine imposes the following mass conservation constraint:
!
!  SUM{ rhoh * vol } = SUM{ rhon * voln  + source * vol }  (1)
!
!  The left-hand side of (1) represents the mass at the present time
!  and the right-hand side represents the mass at previous time plus
!  the mass due sources at present time
!
!  The total mass inside this routine is defined as SUM{rho * vol}.
!  This routine should work with any meaning of "rho" and "vol" as long
!  as they are set appropriately outside this routine, while the
!  true mass (which is the quantity to be conserved) remains SUM{rho * vol}
!---------------------------------------------------------------------------

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER, INTENT(IN)  :: its, itf, jts, jtf, kts, ktf    ! Dimension info
INTEGER, INTENT(IN)  :: irs, irf, jrs, jrf, krs, krf, nt

REAL, INTENT(IN)     :: rhon(its:itf,jts:jtf,kts:ktf,nt)
REAL, INTENT(INOUT)  :: rhoh(its:itf,jts:jtf,kts:ktf,nt)
REAL, INTENT(IN)     :: rhol(its:itf,jts:jtf,kts:ktf,nt)
REAL, INTENT(IN)     :: source(its:itf,jts:jtf,kts:ktf,nt)

REAL, INTENT(IN)     :: voln(irs:irf,jrs:jrf,krs:krf)
REAL, INTENT(IN)     :: vol (irs:irf,jrs:jrf,krs:krf)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRIESTLEY_ALGORITHM2'

! Local variables

REAL    :: du(its:itf,jts:jtf,kts:ktf,nt)
REAL    :: local_sum(irs:irf,jrs:jrf,nt,4)

REAL    :: eps = 1.0e-16    ! Threshold value
REAL    :: w1, w2
REAL    :: s1, s2, s3, deno, temp
REAL    :: global_sum(nt,4)
INTEGER :: i,j,k,kk

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!====================================================
! compute some global sums and setup defaults
!====================================================

global_sum(:,:) = 0.0

!$OMP PARALLEL IF(nt > 1)                                                &
!$OMP DEFAULT(NONE)                                                      &
!$OMP PRIVATE(kk, k, j, i)                                               &
!$OMP SHARED(nt, krs, krf, jrs, jrf, irs, irf, du, vol, rhoh, rhol,      &
!$OMP        local_sum, rhon, voln, source, global_sum)

DO kk = 1, nt ! nt can range from 1 to 100+
  DO k = krs, krf
!$OMP DO SCHEDULE(STATIC)
    DO j = jrs, jrf
      DO i = irs, irf
        du(i,j,k,kk) = vol(i,j,k) *(rhoh(i,j,k,kk)-rhol(i,j,k,kk))
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO


!$OMP DO SCHEDULE(STATIC)
  DO j = jrs, jrf
    DO i = irs, irf
      local_sum(i,j,kk,:) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT

  DO k = krs, krf
!$OMP DO SCHEDULE(STATIC)
    DO j = jrs, jrf
      DO i = irs, irf   
        local_sum(i,j,kk,1) = local_sum(i,j,kk,1)                      &
                                + rhon(i,j,k,kk)   * voln(i,j,k)       &
                                + source(i,j,k,kk) *  vol(i,j,k)

        local_sum(i,j,kk,2) = local_sum(i,j,kk,2)                      &
                                + rhol(i,j,k,kk) *  vol(i,j,k)

        IF (du(i,j,k,kk) < 0.0 ) THEN
          local_sum(i,j,kk,3) = local_sum(i,j,kk,3) - du(i,j,k,kk)
        ELSE
          local_sum(i,j,kk,4) = local_sum(i,j,kk,4) + du(i,j,k,kk)
        END IF

      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
END DO    ! Loop kk

!$OMP END PARALLEL

CALL global_2d_sums(local_sum, irf-irs+1, jrf-jrs+1, 0,0,4*nt,global_sum)

!$OMP PARALLEL IF(nt > 1)                                                &
!$OMP DEFAULT(NONE)                                                      &
!$OMP PRIVATE(kk, k, j, i, deno, s1, s2, s3, w1, w2)                     &
!$OMP SHARED(nt, global_sum, eps, rhoh, rhol, krs,                       &
!$OMP        krf, jrs, jrf, irs, irf, du)
 

DO kk = 1, nt ! nt can range from 1 to 100+

  s1 = global_sum(kk,4)
  s2 = global_sum(kk,3)
  s3 = global_sum(kk,1) - global_sum(kk,2)

  deno = MAX(s1*s1 + s2*s2, eps)
  w1 = ( s2*s2 + s1*s3 + s2*s1 )/ deno
  w2 = ( s1*s1 - s2*s3 + s2*s1 )/ deno


  IF ( ABS(w1) < eps .AND. ABS(w2) < eps ) THEN ! case of w1=w2=0

    s1 = global_sum(kk,1)/MAX(global_sum(kk,2),eps)
!$OMP DO SCHEDULE(STATIC)
    DO k = krs, krf
      DO j = jrs, jrf
        DO i = irs, irf
          rhoh(i,j,k,kk) = rhoh(i,j,k,kk) * s1
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO SCHEDULE(STATIC)
    DO k = krs, krf
      DO j = jrs, jrf
        DO i = irs, irf
          IF (du(i,j,k,kk) < 0.0 ) THEN
            rhoh(i,j,k,kk) = w2*rhoh(i,j,k,kk) + (1.0-w2)*rhol(i,j,k,kk)
          ELSE
            rhoh(i,j,k,kk) = w1*rhoh(i,j,k,kk) + (1.0-w1)*rhol(i,j,k,kk)
          END IF
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

  END IF     ! ABS(w1) < eps

END DO      ! Loop kk
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE priestley_algorithm2

!================================================================== 

SUBROUTINE priestley_algorithm3(dims, dims_s, nfields,                  &
                                volume_n, volume_np1,                   &
                                field_n, field_np1_lo, field_np1_hi,    &
                                source)
   
USE atm_fields_bounds_mod, ONLY: array_dims
   
IMPLICIT NONE
!
! Description: Algorithm for conserving the global integral of a
!              field. Based on Priestley 1993, Monthly Weather 
!              Review, 121, 621--629.
!
! Method:
!
! This algorithm attempts to ensure that the global integral of a field
! at time t_(n+1) is equal to its global integral at time t_n.
! Two approximations are required for the field at time t_(n+1). One
! solution is obtained by a high-order method, and the other by a
! low-order method. The difference between these two solutions provides
! a measure of the local truncation error in the high-order scheme.
!
! Errors are largest where the field is least smooth and it is here
! that conservation errors are likely to be significant, too. 
! The algorithm corrects these errors by blending the high-order and
! low-order solutions. See UMDP16 for further details.
!
! Subroutine arguments
! Input array dimensions
TYPE(array_dims), INTENT(IN) :: dims_s, dims

! Number of fields to be corrected
INTEGER, INTENT(IN)          :: nfields

! Fields at time t_n
REAL, INTENT(IN)    :: field_n(dims_s%i_start:dims_s%i_end,             &
                               dims_s%j_start:dims_s%j_end,             &
                               dims_s%k_start:dims_s%k_end, nfields)

! Fields at time t_(n+1) = t_n + dt:
! - low-order solution
REAL, INTENT(IN)    :: field_np1_lo(dims_s%i_start:dims_s%i_end,        &
                                    dims_s%j_start:dims_s%j_end,        &
                                    dims_s%k_start:dims_s%k_end,        &
                                    nfields)

! - high-order solution, corrected by this subroutine
REAL, INTENT(INOUT) :: field_np1_hi(dims_s%i_start:dims_s%i_end,        &
                                    dims_s%j_start:dims_s%j_end,        &
                                    dims_s%k_start:dims_s%k_end,        &
                                    nfields)

! Cell volumes (or field integration weights)
! - volumes at time t_n
REAL, INTENT(IN)    :: volume_n(dims%i_start:dims%i_end,                &
                                dims%j_start:dims%j_end,                &
                                dims%k_start:dims%k_end)
! - volumes at time t_(n+1)
REAL, INTENT(IN)    :: volume_np1(dims%i_start:dims%i_end,              &
                                  dims%j_start:dims%j_end,              &
                                  dims%k_start:dims%k_end)

! Optional source fields
REAL, INTENT(IN), OPTIONAL :: source(dims_s%i_start:dims_s%i_end,       &
                                     dims_s%j_start:dims_s%j_end,       &
                                     dims_s%k_start:dims_s%k_end,       &
                                     nfields)

! Local parameters
! Number of summations and their meanings:
INTEGER, PARAMETER :: nsums       = 4
INTEGER, PARAMETER :: i_time_n    = 1 ! integral of field at time t_n
INTEGER, PARAMETER :: i_low_order = 2 ! integral of low order solution
INTEGER, PARAMETER :: i_neg_diff  = 3 ! integral where hi - lo < 0
INTEGER, PARAMETER :: i_pos_diff  = 4 ! integral where hi - lo > 0

REAL, PARAMETER    :: eps = 1.0e-16   ! threshold for pathological case
                                      ! in which algorithm breaks down

! Local variables
INTEGER :: i, j, k, ifield            ! loop counters

! Integrated difference between high and low order solutions
REAL    :: diff_hi_lo(dims%i_start:dims%i_end,                          &
                      dims%j_start:dims%j_end,                          &
                      dims%k_start:dims%k_end, nfields)

! Summations required by algorithm on this processor
REAL    :: local_sum(dims%i_start:dims%i_end,                           &
                     dims%j_start:dims%j_end, nfields, nsums)

! Summations across all processors
REAL    :: global_sum(nfields, nsums)     

! Weights for blending high- and low-order solutions:
REAL    :: w_pos(nfields) ! weight to use where hi > lo
REAL    :: w_neg(nfields) ! weight to use where hi < lo

! Working variables for the algorithm
REAL    :: s1, s2, s3, deno

!====================================================
! compute some global sums and setup defaults
!====================================================         

! Weighted ifference between high and low order solutions
DO ifield = 1, nfields
  DO k = dims%k_start, dims%k_end
    DO j = dims%j_start, dims%j_end
      DO i = dims%i_start, dims%i_end
        diff_hi_lo(i,j,k,ifield) = (field_np1_hi(i,j,k,ifield) -        &
                                    field_np1_lo(i,j,k,ifield))         &
                                   * volume_np1(i,j,k)
      END DO                   
    END DO
  END DO  
END DO

! Compute sums
global_sum(:,:)  = 0.0
local_sum(:,:,:,:) = 0.0
DO ifield = 1, nfields
  DO k = dims%k_start, dims%k_end
    DO j = dims%j_start, dims%j_end
      DO i = dims%i_start, dims%i_end

! Integral of field at time t_n
        local_sum(i, j, ifield ,i_time_n) =                             &
          local_sum(i, j, ifield, i_time_n)                             &
            + field_n(i,j,k,ifield) * volume_n(i,j,k)

! Integral of low-order solution at time t_(n+1)
        local_sum(i, j, ifield, i_low_order) =                          &
          local_sum(i, j, ifield, i_low_order)                          &
            + field_np1_lo(i,j,k,ifield) * volume_np1(i,j,k)

        IF (diff_hi_lo(i,j,k,ifield) < 0.0 ) THEN
! Absolute integral of difference hi - lo where this is negative        
          local_sum(i, j, ifield, i_neg_diff) =                         &
            local_sum(i, j, ifield, i_neg_diff)                         &
              - diff_hi_lo(i,j,k,ifield)
        ELSE
! Absolute integral of difference hi - lo where this is positive
          local_sum(i, j, ifield, i_pos_diff) =                         &
            local_sum(i, j, ifield, i_pos_diff)                         &
              + diff_hi_lo(i,j,k,ifield)
        END IF 
                
      END DO
    END DO
  END DO
END DO

! Add source, if present
IF (PRESENT(source)) THEN
  DO ifield = 1, nfields
    DO k = dims%k_start, dims%k_end
      DO j = dims%j_start, dims%j_end
        DO i = dims%i_start, dims%i_end
          local_sum(i, j, ifield, i_time_n) =                           &
            local_sum(i, j, ifield, i_time_n)                           &
              + source(i,j,k,ifield) * volume_np1(i,j,k)
        END DO
      END DO
    END DO
  END DO
END IF       

! Sum across all processors
CALL global_2d_sums(local_sum, dims%i_len, dims%j_len, dims%halo_i,     &
                    dims%halo_j, nfields*nsums, global_sum) 

! Compute the weights for blending high- and low-order solutions          
DO ifield = 1, nfields
   
  s3 = global_sum(ifield, i_time_n) - global_sum(ifield, i_low_order)
  s2 = global_sum(ifield, i_neg_diff)
  s1 = global_sum(ifield, i_pos_diff)

  deno = MAX(s1*s1 + s2*s2, eps) 
  w_pos(ifield) = ( s2*s2 + s1*s3 + s2*s1 )/ deno
  w_neg(ifield) = ( s1*s1 - s2*s3 + s2*s1 )/ deno
         
END DO              

! Blend high- and low-order solutions for each field
DO ifield = 1, nfields
 
 ! Trap the pathological case, w_pos = w_neg = 0
  IF ( ABS(w_pos(ifield)) < eps .AND. ABS(w_neg(ifield)) < eps ) THEN
    s1 = global_sum(ifield, i_time_n) /                                 &
         MAX(global_sum(ifield, i_low_order), eps) 
    field_np1_hi(:,:,:,ifield) = field_np1_hi(:,:,:,ifield) * s1   
  ELSE 
    DO k = dims%k_start, dims%k_end
      DO j = dims%j_start, dims%j_end
        DO i = dims%i_start, dims%i_end
! Blend high- and low-order solutions, depending on the sign of 
! their difference
          IF (diff_hi_lo(i,j,k,ifield) < 0.0 ) THEN
            field_np1_hi(i,j,k,ifield) =                                &
                      w_neg(ifield) * field_np1_hi(i,j,k,ifield)        &
              + (1.0-w_neg(ifield) )* field_np1_lo(i,j,k,ifield) 
          ELSE
            field_np1_hi(i,j,k,ifield) =                                &
                      w_pos(ifield) * field_np1_hi(i,j,k,ifield)        &
              + (1.0-w_pos(ifield)) * field_np1_lo(i,j,k,ifield)  
          END IF                    
        END DO
      END DO                       
    END DO
  END IF      ! ABS(w_pos) < eps 

END DO       ! Loop over fields
         
RETURN      

END SUBROUTINE PRIESTLEY_ALGORITHM3

! ==================================================================
! Imposes a minimum value f > fmin (e.g. for removing negative values)
! while maintaining conservation
!
! (i.e., it imposes the following conditions):
!
!      f2 = max(f,fmin)               (1)
! and
!   SUM{ f2 * vol } = SUM{ f * vol }  (2)
!
! at exit  f --> f2
!
! ==================================================================

SUBROUTINE impose_minima_with_conservation(f, fmin, vol,                &
                                       its, itf, jts, jtf, kts, ktf,    &
                                       irs, irf, jrs, jrf, krs, krf, nt)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER, INTENT(IN)  :: its, itf, jts, jtf, kts, ktf    ! Dimension info
INTEGER, INTENT(IN)  :: irs, irf, jrs, jrf, krs, krf, nt

REAL, INTENT(INOUT)  :: f(its:itf,jts:jtf,kts:ktf,nt)
REAL, INTENT(IN)     :: vol(irs:irf,jrs:jrf,krs:krf)
REAL, INTENT(IN)     :: fmin(nt)                ! Targeted min value

!    locals

REAL    :: f2(its:itf,jts:jtf,kts:ktf,nt)
REAL    :: global_sum(nt,2), local_sum(irs:irf,jrs:jrf,nt,2), &
           local_sum2(irs:irf,jrs:jrf,1)
REAL    :: sum_vol(1), deno, a, b
INTEGER :: i,j,k,kk

REAL, PARAMETER :: eps = 1.0e-16   ! Threshold value

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IMPOSE_MINIMA_WITH_CONSERVATION'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

sum_vol(1) = 0.0
global_sum(:,:) = 0.0

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP PRIVATE(kk, k, j, i)                                             &
!$OMP SHARED(nt, f, f2, fmin, local_sum,                               &
!$OMP        local_sum2, jrs, jrf, irs, irf, krs, krf, vol)            

!$OMP DO SCHEDULE(STATIC)
DO kk = 1, nt ! nt can range from 1 to 100+
  DO k = krs, krf
    DO j = jrs, jrf
      DO i = irs, irf    
        f2(i,j,k,kk) = MAX( f(i,j,k,kk), fmin(kk) )
      END DO
    END DO
  END DO
  local_sum(:,:,kk,:)  = 0.0
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO kk = 1, nt ! nt can range from 1 to 100+
  DO k = krs, krf
    DO j = jrs, jrf
      DO i = irs, irf
        local_sum(i,j,kk,1) = local_sum(i,j,kk,1)                &
                               + f(i,j,k,kk)* vol(i,j,k)
        local_sum(i,j,kk,2) = local_sum(i,j,kk,2)                &
                               + f2(i,j,k,kk) * vol(i,j,k)
      END DO
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO j = jrs, jrf
  DO i = irs, irf 
    local_sum2(i,j,1) = 0.0
  END DO
END DO
!$OMP END DO NOWAIT

DO k = krs, krf
!$OMP DO SCHEDULE(STATIC)
  DO j = jrs, jrf
    DO i = irs, irf 
      local_sum2(i,j,1) = local_sum2(i,j,1) + vol(i,j,k)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO

!$OMP END PARALLEL

CALL global_2d_sums(local_sum, irf-irs+1, jrf-jrs+1, 0, 0,        &
                    2*nt, global_sum)

CALL global_2d_sums(local_sum2,irf-irs+1,jrf-jrs+1,0,0,1,sum_vol)

DO kk = 1, nt ! nt can range from 1 to 100+
  deno = global_sum(kk,2) - fmin(kk)*sum_vol(1)

  IF (ABS(deno) > eps ) THEN

    a = ( global_sum(kk,1)-fmin(kk)*sum_vol(1) )/deno
    b = fmin(kk)*( 1.0 - a )

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP PRIVATE(k, j, i)                                                 &
!$OMP SHARED( krs, krf, jrs, jrf, kk, irs, irf, f, f2, a, b)  
    DO k = krs, krf
      DO j = jrs, jrf
        DO i = irs, irf
          f(i,j,k,kk) = a * f2(i,j,k,kk) + b
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE impose_minima_with_conservation

END MODULE priestley_algorithm_mod


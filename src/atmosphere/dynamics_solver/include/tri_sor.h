! -- begin tri_sor.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver
INTEGER,  INTENT(IN)    :: NoIts

#if defined(C_DP_HLM)
INTEGER, PARAMETER :: real_eg_hlm_kind=real64
#else
INTEGER, PARAMETER :: real_eg_hlm_kind=real32
#endif

! Input vector x and output vector Px (i.e. x multiplied by the
! Pre/Post-conditioner).

REAL(KIND=x_kind),     INTENT(INOUT) ::                       &
                           x(dims_s%i_start:dims_s%i_end,     &
                             dims_s%j_start:dims_s%j_end,     &
                             dims_s%k_start:dims_s%k_end)

REAL(KIND=Px_kind),    INTENT(OUT)   ::                       &
                          Px(dims_s%i_start:dims_s%i_end,     &
                             dims_s%j_start:dims_s%j_end,     &
                             dims_s%k_start:dims_s%k_end)

REAL(KIND=real_eg_hlm_kind),    INTENT(IN)   ::               &
            ltHlm_Ln   (dims%i_start:pdims%i_end,             &
                        dims%j_start:dims%j_end,              &
                        dims%k_start:dims%k_end),             &
            ltHlm_Ls   (dims%i_start:dims%i_end,              &
                        dims%j_start:dims%j_end,              &
                        dims%k_start:dims%k_end),             &
            ltHlm_Le   (dims%i_start:dims%i_end,              &
                        dims%j_start:dims%j_end,              &
                        dims%k_start:dims%k_end),             &
            ltHlm_Lw   (dims%i_start:dims%i_end,              &
                        dims%j_start:dims%j_end,              &
                        dims%k_start:dims%k_end),             &
            ltHlm_Lu   (dims%i_start:dims%i_end,              &
                        dims%j_start:dims%j_end,              &
                        dims%k_start:dims%k_end),             &
            ltHlm_Ld   (dims%i_start:dims%i_end,              &
                        dims%j_start:dims%j_end,              &
                        dims%k_start:dims%k_end),             &
            ltHu_k     (dims%i_start:dims%i_end,              &
                        dims%k_start:dims%k_end,              &
                        dims%j_start:dims%j_end),             &
            ltHd_k     (dims%i_start:dims%i_end,              &
                        dims%k_start:dims%k_end,              &
                        dims%j_start:dims%j_end)

!     Local variables

REAL(KIND=local_kind) :: rlx
INTEGER :: io, jo, ioff, i, j, k, it, isweep
INTEGER :: i1, i2, j1
INTEGER :: i_start, i_end, j_start, j_end

REAL(KIND=local_kind) :: tmp(dims_s%i_start:dims_s%i_end,     &
               dims_s%k_start:dims_s%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

io = MOD(datastart(1) + datastart(2), 2)
rlx = 1.5

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k)   &
!$OMP SHARED( dims_s, Px ) 
DO k = dims_s%k_start,dims_s%k_end
  Px(:,:,k) = 0.0
END DO
!$OMP END PARALLEL DO

i_start = dims%i_start 
j_start = dims%j_start 
i_end   = dims%i_end   
j_end   = dims%j_end   
    
DO it = 1, NoIts
   
  ! Build RHS using the red/black ordering

  DO isweep = 0, 1

    jo      = io + isweep + 1

    ! Now solve tridiagonal matrix

    ! Forward sweep

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)          &
!$OMP PRIVATE(i,j,k,ioff, tmp, i1, i2, j1)                &
!$OMP SHARED(j_start, j_end, jo, it, isweep, i_start,     &
!$OMP i_end, x, dims_s, lthlm_ld, lthd_k, px,             &
!$OMP lthlm_lw, lthlm_le, lthlm_ls, lthlm_ln, lthlm_lu,   &
!$OMP lthu_k, rlx)
    DO j = j_start, j_end

      ! Calculate offsets and indexing
      ioff = MOD(jo+j,2)
      i1 = -(i_start+ioff)
      i2 = (i_end - (i_start+ioff))/2 + (i_start+ioff)
      IF (ioff == 0) THEN
        j1 = 0
      ELSE
        j1 = (i_end - i_start)/2
      END IF

      IF ( it == 1 .AND. isweep == 0 ) THEN
        k = 1
        DO i = i_start+ioff, i2
          tmp(i,k) = x(2*i+i1,j,k)
        END DO

        DO k = 2, dims_s%k_end-1
          DO i = i_start+ioff, i2
            tmp(i,k) = x(2*i+i1,j,k)

            tmp(i,k) = (tmp(i,k) - ltHlm_Ld(i+j1,j,k)*tmp(i,k-1))    &
                       *ltHd_k(i+j1,k,j)
          END DO
        END DO

        k = dims_s%k_end
        DO i = i_start+ioff, i2

          tmp(i,k) = x(2*i+i1,j,k)
          tmp(i,k) = (tmp(i,k) - ltHlm_Ld(i+j1,j,k)*tmp(i,k-1))      &
                       *ltHd_k(i+j1,k,j)
        END DO
      ELSE
        k = 1
        DO i = i_start+ioff, i2
          tmp(i,k) = x(2*i+i1,j,k) - Px(2*i+i1,j,k)                  &
                      - ltHlm_Lw(i+j1,j,k)*Px(2*i+i1-1,j,k)          &
                      - ltHlm_Le(i+j1,j,k)*Px(2*i+i1+1,j,k)          &
                      - ltHlm_Ls(i+j1,j,k)*Px(2*i+i1,j-1,k)          &
                      - ltHlm_Ln(i+j1,j,k)*Px(2*i+i1,j+1,k)          &
                      - ltHlm_Lu(i+j1,j,k)*Px(2*i+i1,j,k+1)
        END DO

        DO k = 2,dims_s%k_end-1
          DO i = i_start+ioff, i2
            tmp(i,k) = x(2*i+i1,j,k) -  Px(2*i+i1,j,k)               &
                        - ltHlm_Lw(i+j1,j,k)*Px(2*i+i1-1,j,k)        &
                        - ltHlm_Le(i+j1,j,k)*Px(2*i+i1+1,j,k)        &
                        - ltHlm_Ls(i+j1,j,k)*Px(2*i+i1,j-1,k)        &
                        - ltHlm_Ln(i+j1,j,k)*Px(2*i+i1,j+1,k)        &
                        - ltHlm_Ld(i+j1,j,k)*Px(2*i+i1,j,k-1)        &
                        - ltHlm_Lu(i+j1,j,k)*Px(2*i+i1,j,k+1)

            tmp(i,k) = (tmp(i,k) - ltHlm_Ld(i+j1,j,k)*tmp(i,k-1))    &
                       *ltHd_k(i+j1,k,j)
          END DO
        END DO

        k = dims_s%k_end
        DO i = i_start+ioff, i2
          tmp(i,k) = x(2*i+i1,j,k) - Px(2*i+i1,j,k)                  &
                   - ltHlm_Lw(i+j1,j,k)*Px(2*i+i1-1,j,k)             &
                   - ltHlm_Le(i+j1,j,k)*Px(2*i+i1+1,j,k)             &
                   - ltHlm_Ls(i+j1,j,k)*Px(2*i+i1,j-1,k)             &
                   - ltHlm_Ln(i+j1,j,k)*Px(2*i+i1,j+1,k)             &
                   - ltHlm_Ld(i+j1,j,k)*Px(2*i+i1,j,k-1)

          tmp(i,k) = (tmp(i,k) - ltHlm_Ld(i+j1,j,k)*tmp(i,k-1))      &
                       *ltHd_k(i+j1,k,j)
        END DO
      END IF

      ! Backward sweep

      k = dims_s%k_end
      DO i = i_start+ioff, i2
        Px(2*i+i1,j,k) = Px(2*i+i1,j,k) + rlx*tmp(i,k)
      END DO

      DO k = dims_s%k_end-1, 1, -1
        DO i = i_start+ioff, i2
          tmp(i,k) = tmp(i,k) - ltHu_k(i+j1,k,j)*tmp(i,k+1)
          Px(2*i+i1,j,k)  = Px(2*i+i1,j,k)  + rlx*tmp(i,k)
        END DO
      END DO

    END DO
!$OMP END PARALLEL DO
    IF ( it < NoIts .OR. isweep == 0 ) THEN

      ! isweep == 0, red   elements
      ! isweep == 1, black elements
      CALL swap_bounds_RB(Px,                                  &
           dims_s%i_len - 2*dims_s%halo_i,                     &
           dims_s%j_len - 2*dims_s%halo_j,                     &
           dims_s%k_len, (isweep==0))


    END IF
  END DO

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
! -- end tri_sor.h --

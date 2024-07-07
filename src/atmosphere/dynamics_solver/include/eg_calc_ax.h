! -- begin eg_calc_ax.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver

CALL swap_bounds(x,                                                         &
           dims_s%i_len - 2*dims_s%halo_i,                                  &
           dims_s%j_len - 2*dims_s%halo_j,                                  &
           dims_s%k_len,                                                    &
           dims_s%halo_i, dims_s%halo_j, fld_type_p, swap_field_is_scalar,  &
           do_corners_arg=.FALSE.)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)             &
!$OMP&            SHARED(dims, Ax, x, Hlm_lw, Hlm_le, Hlm_ls,               &
!$OMP&                   Hlm_ln, Hlm_lu, Hlm_ld, dims_s)
DO k = dims_s%k_start,dims_s%k_end
  DO j = dims%j_start, dims%j_end
    DO i = dims%i_start, dims%i_end
      Ax(i,j,k) =  x(i,j,k) +                                               &
                      Hlm_Lw(i,j,k)*x(i-1,j,k)+Hlm_Le(i,j,k)*x(i+1,j,k)     &
                     +Hlm_Ls(i,j,k)*x(i,j-1,k)+Hlm_Ln(i,j,k)*x(i,j+1,k)

    END DO
  END DO

  IF ( k == 1 ) THEN
    DO j = dims%j_start, dims%j_end
      DO i = dims%i_start, dims%i_end
        Ax(i,j,k) = Ax(i,j,k) + Hlm_Lu(i,j,k)*x(i,j,k+1)
      END DO
    END DO
  ELSE IF ( k == dims_s%k_end) THEN
    DO j = dims%j_start, dims%j_end
      DO i = dims%i_start, dims%i_end
        Ax(i,j,k) = Ax(i,j,k) + Hlm_Ld(i,j,k)*x(i,j,k-1)
      END DO
    END DO
  ELSE
    DO j = dims%j_start, dims%j_end
      DO i = dims%i_start, dims%i_end
        Ax(i,j,k) = Ax(i,j,k) +                                         &
                        Hlm_Ld(i,j,k)*x(i,j,k-1)+Hlm_Lu(i,j,k)*x(i,j,k+1)
      END DO
    END DO
  END IF

END DO
!$OMP END PARALLEL DO

! -- end eg_calc_ax.h --

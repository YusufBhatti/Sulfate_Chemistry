! -- begin gmres1.h -
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver

row_length = pdims%i_end
rows       = pdims%j_end
levels     = pdims%k_end

DO k = 1, 3     
  DO j = 1, rows
    DO i = 1, row_length
      k_sum(i,j,k) = 0.0
    END DO
  END DO
END DO

len1 = 4
len2 = ( rows - 1 )/len1 + 1

! The following code is the straightforward implementation of the
! OpenMP code following (for reference)
!      DO k = 1, levels
!        DO j = 1, rows
!            DO i = 1, row_length
!              k_sum(i,j,1) = k_sum(i,j,1) + t(i,j,k)**2
!              k_sum(i,j,2) = k_sum(i,j,2) + t(i,j,k)*s(i,j,k)
!              k_sum(i,j,3) = k_sum(i,j,3) + s(i,j,k)**2
!            END DO
!         END DO
!      END DO
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                      &
!$OMP&            PRIVATE(jj,j1,j2,i,j,k)                             &
!$OMP&            SHARED(len1, len2, rows, levels, row_length,        &
!$OMP&            k_sum, t, s)
DO jj = 1, len1

  j1 = 1 + (jj-1)*len2
  j2 = MIN( j1+len2-1, rows)

  DO k = 1, levels-1, 2

    DO j = j1, j2
      DO i = 1, row_length
        k_sum(i,j,1) = k_sum(i,j,1)+t(i,j,k)**2 + t(i,j,k+1)**2
        k_sum(i,j,2) = k_sum(i,j,2)+t(i,j,k)*s(i,j,k)           &
                                   +t(i,j,k+1)*s(i,j,k+1)
        k_sum(i,j,3) = k_sum(i,j,3)+s(i,j,k)**2 + s(i,j,k+1)**2
      END DO
    END DO
  END DO

  IF ( MOD(levels,2) == 1 ) THEN
    k = levels
    DO j = j1, j2
      DO i = 1, row_length
        k_sum(i,j,1) = k_sum(i,j,1) + t(i,j,k)**2
        k_sum(i,j,2) = k_sum(i,j,2) + t(i,j,k)*s(i,j,k)
        k_sum(i,j,3) = k_sum(i,j,3) + s(i,j,k)**2
      END DO
    END DO
  END IF
END DO
!$OMP END PARALLEL DO

! -- end gmres1.h -

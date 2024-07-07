! -- begin eg_inner_prod.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver


row_length = dims%i_end
rows       = dims%j_end
levels     = dims%k_end

DO j = 1, rows
  DO i = 1, row_length
    k_sum(i,j) = 0.0
  END DO
END DO

len1 = 4
len2 = ( rows - 1 )/len1 + 1

! The following code is the simple algorithm that is OpenMPed below and
! is given for reference.
!      DO k = 1, levels
!        DO j = 1, rows
!           DO i = 1, row_length
!               k_sum(i,j) = k_sum(i,j) + x(i,j,k)*y(i,j,k)
!            END DO
!         END DO
!      END DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jj,j1,j2,i,j,k) DEFAULT(NONE) &
!$OMP SHARED( len1,len2,rows,levels,row_length,k_sum,x,y)
DO jj = 1, len1
  j1 = 1 + (jj-1)*len2
  j2 = MIN( j1+len2-1, rows)
  DO k = 1, levels-1, 2
    DO j = j1, j2
      DO i = 1, row_length
        k_sum(i,j) = k_sum(i,j) + x(i,j,k)*y(i,j,k)        &
                                + x(i,j,k+1)*y(i,j,k+1)
      END DO
    END DO
  END DO
  IF ( MOD(levels,2) == 1 ) THEN
    k = levels
    DO j = j1, j2
      DO i = 1, row_length
        k_sum(i,j) = k_sum(i,j) + x(i,j,k)*y(i,j,k)
      END DO
    END DO
  END IF
END DO
!$OMP END PARALLEL DO

! -- end eg_inner_prod.h --

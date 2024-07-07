! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE gmres1_coef_mod

! Description: Function to calculate the inner product of two fields
!              for use by the Linear Solvers
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
IMPLICIT NONE

PRIVATE

PUBLIC :: gmres1_coef
INTERFACE gmres1_coef
MODULE PROCEDURE gmres1_coef_s,                               &
     gmres1_coef_d
END INTERFACE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GMRES1_COEF_MOD'

CONTAINS

FUNCTION gmres1_coef_d(s,t,err,dims,dims_s)

USE global_2d_sums_mod,   ONLY: global_2d_sums
USE atm_fields_bounds_mod
USE um_types,             ONLY: real32, real64

IMPLICIT NONE

TYPE (array_dims) dims, dims_s

REAL(KIND=real64),  INTENT(IN)   ::                             &
     s(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end),                          &
     t(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end)

REAL(KIND=real64)             :: gmres1_coef_d
REAL(KIND=real64),INTENT(OUT) :: err
REAL(KIND=real64)             :: temp(3)

INTEGER                       :: i, j, k
INTEGER                       :: rows, row_length, levels
REAL(KIND=real64)             :: k_sum(dims%i_end,dims%j_end,3)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GMRES1_COEF_D'

INTEGER :: j1,j2,jj,len1,len2

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#include "gmres1.h"

CALL global_2d_sums(k_sum, row_length, rows, 0, 0, 3, temp)

gmres1_coef_d = temp(2)/temp(1)

err         = temp(3) + gmres1_coef_d**2 *temp(1)                 &
                      - 2.0*gmres1_coef_d*temp(2)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION gmres1_coef_d

!!!!!!!!!!!! End of double

FUNCTION gmres1_coef_s(s,t,err,dims,dims_s)

USE global_2d_sums_mod,   ONLY: global_2d_sums_sp, global_sum_reprod_dd
USE atm_fields_bounds_mod
USE um_types,             ONLY: real32, real64

IMPLICIT NONE

TYPE (array_dims) dims, dims_s

REAL(KIND=real32),  INTENT(IN)   ::                             &
     s(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end),                          &
     t(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end)

REAL(KIND=real32)             :: gmres1_coef_s
REAL(KIND=real32),    INTENT(OUT) :: err
REAL(KIND=real32)                        :: temp(3)

INTEGER                       :: i, j, k
INTEGER                       :: rows, row_length, levels
REAL(KIND=real32)             :: k_sum(dims%i_end,dims%j_end,3)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GMRES1_COEF_S'

INTEGER :: j1,j2,jj,len1,len2

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#include "gmres1.h"

! now have SP version
CALL global_2d_sums_sp(k_sum, row_length, rows, 0, 0, 3, temp, &
                       global_sum_method_override=global_sum_reprod_dd)

gmres1_coef_s = temp(2)/temp(1)

err         = temp(3) + gmres1_coef_s**2 *temp(1)                 &
                      - 2.0*gmres1_coef_s*temp(2)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION gmres1_coef_s

END MODULE gmres1_coef_mod

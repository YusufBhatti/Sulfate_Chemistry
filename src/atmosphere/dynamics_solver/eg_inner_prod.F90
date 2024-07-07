! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_inner_prod_mod

USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
IMPLICIT NONE

PRIVATE

PUBLIC :: eg_inner_prod
INTERFACE eg_inner_prod
MODULE PROCEDURE eg_inner_prod_dp_dp,                           &
     eg_inner_prod_dp_sp,                                       &
     eg_inner_prod_sp_dp,                                       &
     eg_inner_prod_sp_sp
END INTERFACE


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_INNER_PROD_MOD'

CONTAINS

FUNCTION eg_inner_prod_dp_dp(x,y,dims,dims_s)

USE global_2d_sums_mod,   ONLY: global_2d_sums
USE atm_fields_bounds_mod
USE um_types, ONLY: real64

IMPLICIT NONE

!
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
!   This code is written to UM programming standards
!

! Array dimensions

TYPE (array_dims) dims,dims_s

! Input vectors

REAL(KIND=real64), INTENT(IN)   ::                              &
     x(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end),                          &
     y(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end)

REAL(KIND=real64)             :: eg_inner_prod_dp_dp

INTEGER                       :: i, j, k
REAL(KIND=real64)             :: temp(1)
REAL(KIND=real64)             :: k_sum(dims%i_end,dims%j_end)
INTEGER                       :: rows, row_length, levels

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_INNER_PROD_DP_DP'

INTEGER :: j1,j2,jj,len1,len2

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#include "eg_inner_prod.h"

CALL global_2d_sums(k_sum,dims%i_end,dims%j_end , 0, 0, 1, temp)

eg_inner_prod_dp_dp = temp(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION eg_inner_prod_dp_dp


!!!!!!!!!!!!!!!!!  End  of function dp_dp !!!!!!!!!!!!!!


FUNCTION eg_inner_prod_dp_sp(x,y,dims,dims_s)

USE global_2d_sums_mod,   ONLY: global_2d_sums_sp, global_sum_reprod_dd
USE atm_fields_bounds_mod
USE um_types, ONLY: real32,real64

IMPLICIT NONE

!
! Description: Function to calculate the inner product of two fields
!              for use by the Linear Solvers
!
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards.
!

! Array dimensions

TYPE (array_dims) dims,dims_s

! Input vectors

REAL(KIND=real64), INTENT(IN)   ::                              &
     x(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end)
REAL (KIND=real32) , INTENT(IN) ::                              &
     y(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end)

REAL (KIND=real32)            :: eg_inner_prod_dp_sp

INTEGER                       :: i, j, k
REAL (KIND=real32)            :: temp(1)
REAL(KIND=real32)             :: k_sum(dims%i_end,dims%j_end)
INTEGER                       :: rows, row_length, levels

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_INNER_PROD_DP_SP'

INTEGER :: j1,j2,jj,len1,len2

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#include "eg_inner_prod.h"

! now have SP version
CALL global_2d_sums_sp(k_sum,dims%i_end,dims%j_end , 0, 0, 1, temp, &
                       global_sum_method_override=global_sum_reprod_dd)

eg_inner_prod_dp_sp = temp(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION eg_inner_prod_dp_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! End of DP_SP !!!!


FUNCTION eg_inner_prod_sp_dp(x,y,dims,dims_s)

USE global_2d_sums_mod,   ONLY: global_2d_sums_sp, global_sum_reprod_dd
USE atm_fields_bounds_mod
USE um_types, ONLY: real32,real64

IMPLICIT NONE

!
! Description: Function to calculate the inner product of two fields
!              for use by the Linear Solvers
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards.
!

! Array dimensions

TYPE (array_dims) dims,dims_s

! Input vectors

REAL (KIND=real32), INTENT(IN)   ::                             &
     x(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end)
REAL (KIND=real64), INTENT(IN) ::                               &
     y(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end)

REAL(KIND=real32)             :: eg_inner_prod_sp_dp

INTEGER                       :: i, j, k
REAL(KIND=real32)             :: temp(1)
REAL (KIND=real32)            :: k_sum(dims%i_end,dims%j_end)
INTEGER                       :: rows, row_length, levels

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_INNER_PROD_SP_DP'

INTEGER :: j1,j2,jj,len1,len2

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#include "eg_inner_prod.h"

! now have sp version
CALL global_2d_sums_sp(k_sum,dims%i_end,dims%j_end , 0, 0, 1, temp, &
                       global_sum_method_override=global_sum_reprod_dd)

eg_inner_prod_sp_dp = temp(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION eg_inner_prod_sp_dp

!!!!!!!!!!!!!!  single double


FUNCTION eg_inner_prod_sp_sp(x,y,dims,dims_s)

USE global_2d_sums_mod,   ONLY: global_2d_sums_sp, global_sum_reprod_dd
USE atm_fields_bounds_mod
USE um_types

IMPLICIT NONE

!
! Description: Function to calculate the inner product of two fields
!              for use by the Linear Solvers
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards.
!

! Array dimensions

TYPE (array_dims) dims,dims_s

! Input vectors

REAL (KIND=real32), INTENT(IN)   ::                                    &
     x(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end)
REAL (KIND=real32), INTENT(IN) ::                              &
     y(dims_s%i_start:dims_s%i_end,                           &
       dims_s%j_start:dims_s%j_end,                           &
       dims_s%k_start:dims_s%k_end)

REAL (KIND=real32)            :: eg_inner_prod_sp_sp

INTEGER                       :: i, j, k
REAL(KIND=real32)             :: temp(1)
REAL (KIND=real32)            :: k_sum(dims%i_end,dims%j_end)
INTEGER                       :: rows, row_length, levels

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_INNER_PROD_SP_SP'

INTEGER :: j1,j2,jj,len1,len2

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#include "eg_inner_prod.h"

!now have SP version
CALL global_2d_sums_sp(k_sum,dims%i_end,dims%j_end , 0, 0, 1, temp, &
                       global_sum_method_override=global_sum_reprod_dd)

eg_inner_prod_sp_sp = temp(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION eg_inner_prod_sp_sp

END MODULE eg_inner_prod_mod

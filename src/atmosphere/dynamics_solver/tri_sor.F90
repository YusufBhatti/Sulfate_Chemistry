! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE tri_sor_mod

! Description: SOR type preconditioner using the
!              tridiagonal part of the matrix
!              + red-black ordering to improve parallelization
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver
!
!   Language: Fortran 90.

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod
USE Field_Types
USE helmholtz_const_matrix_mod
USE UM_ParVars, ONLY: datastart
USE halo_exchange, ONLY: swap_bounds_RB
USE um_types, ONLY: real32, real64

IMPLICIT NONE


PRIVATE

PUBLIC :: Tri_sor
INTERFACE Tri_sor
MODULE PROCEDURE Tri_sor_dp_dp,                                &
     Tri_sor_dp_sp,                                            &
     Tri_sor_sp_dp,                                            &
     Tri_sor_sp_sp
END INTERFACE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TRI_SOR_MOD'

CONTAINS

SUBROUTINE Tri_sor_dp_dp(Px,x,NoIts,dims,dims_s,             &
ltHu_k,ltHd_k,ltHlm_Lw,ltHlm_Le,ltHlm_Ls,ltHlm_Ln,ltHlm_Ld,ltHlm_Lu)

IMPLICIT NONE

TYPE (array_dims) dims,dims_s

INTEGER, PARAMETER :: x_kind     = real64
INTEGER, PARAMETER :: Px_kind    = real64
INTEGER, PARAMETER :: local_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='TRI_SOR_DP_DP'

#include "tri_sor.h"

END SUBROUTINE Tri_sor_dp_dp

!!!!!!!!!!!!!!!!! END OF Double Double !!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Tri_sor_dp_sp(Px,x,NoIts,dims, dims_s,             &
ltHu_k,ltHd_k,ltHlm_Lw,ltHlm_Le,ltHlm_Ls,ltHlm_Ln,ltHlm_Ld,ltHlm_Lu)

IMPLICIT NONE

TYPE (array_dims) dims,dims_s

INTEGER, PARAMETER :: x_kind     = real32
INTEGER, PARAMETER :: Px_kind    = real64
INTEGER, PARAMETER :: local_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='TRI_SOR_DP_SP'

#include "tri_sor.h"

END SUBROUTINE Tri_sor_dp_sp

!!!!!!!!!!!! END OF Double Single !!!!!!!!!!!!!!!!!!

SUBROUTINE Tri_sor_sp_dp(Px,x,NoIts,dims,dims_s,             &
ltHu_k,ltHd_k,ltHlm_Lw,ltHlm_Le,ltHlm_Ls,ltHlm_Ln,ltHlm_Ld,ltHlm_Lu)

IMPLICIT NONE

TYPE (array_dims) dims,dims_s

INTEGER, PARAMETER :: x_kind     = real64
INTEGER, PARAMETER :: Px_kind    = real32
INTEGER, PARAMETER :: local_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='TRI_SOR_SP_DP'

#include "tri_sor.h"

END SUBROUTINE Tri_sor_sp_dp

!!!!!!!!!!!! END OF Single Double !!!!!!!!!!!!!!!!!!

SUBROUTINE Tri_sor_sp_sp(Px,x,NoIts,dims, dims_s,             &
ltHu_k,ltHd_k,ltHlm_Lw,ltHlm_Le,ltHlm_Ls,ltHlm_Ln,ltHlm_Ld,ltHlm_Lu)

IMPLICIT NONE

TYPE (array_dims) dims,dims_s

INTEGER, PARAMETER :: x_kind     = real32
INTEGER, PARAMETER :: Px_kind    = real32
INTEGER, PARAMETER :: local_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='TRI_SOR_SP_SP'

#include "tri_sor.h"

END SUBROUTINE Tri_sor_sp_sp

!!!!!!!!!!!! END OF Single Single !!!!!!!!!!!!!!!!!!

END MODULE tri_sor_mod

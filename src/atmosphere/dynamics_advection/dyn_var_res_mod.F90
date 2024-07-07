! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Constants necessary for variable resolution advection scheme.

MODULE dyn_var_res_Mod

! Description:
! This module is used to hold variable resolution values set in atmos_init.
!
! Method:
! The arrays are calculated in routine control/top_level/
! set_var_res called from atmos_init,
! and used in variable resolution dynamics_advection routines.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

! Arguments

! VarRes Array co-ordinates in radians

IMPLICIT NONE

REAL, ALLOCATABLE :: glambda_p ( : )
REAL, ALLOCATABLE :: glambda_u ( : )
REAL, ALLOCATABLE :: lambda_p ( : )
REAL, ALLOCATABLE :: lambda_u ( : )
REAL, ALLOCATABLE :: phi_p ( : , : )
REAL, ALLOCATABLE :: phi_v( : , : )
REAL, ALLOCATABLE :: gdlambda_p ( : )
REAL, ALLOCATABLE :: gdlambda_u ( : )
REAL, ALLOCATABLE :: dphi_p ( : , : )
REAL, ALLOCATABLE :: dphi_v( : , : )
REAL, ALLOCATABLE :: grecip_dlamp(:)
REAL, ALLOCATABLE :: grecip_dlamu(:)
REAL, ALLOCATABLE :: recip_dphip( : , : )
REAL, ALLOCATABLE :: recip_dphiv( : , : )
REAL, ALLOCATABLE :: wt_lambda_p(:)
REAL, ALLOCATABLE :: wt_lambda_u(:)
REAL, ALLOCATABLE :: wt_phi_p( : , : )
REAL, ALLOCATABLE :: wt_phi_v( : , : )
REAL, ALLOCATABLE :: lambda_p_rm(:)
REAL, ALLOCATABLE :: lambda_p_rp(:)
REAL, ALLOCATABLE :: lambda_u_rm(:)
REAL, ALLOCATABLE :: lambda_u_rp(:)
REAL, ALLOCATABLE :: phi_p_rm( : , : )
REAL, ALLOCATABLE :: phi_p_rp( : , : )
REAL, ALLOCATABLE :: phi_v_rm( : , : )
REAL, ALLOCATABLE :: phi_v_rp( : , : )
REAL, ALLOCATABLE :: Recip_lambda_p_m ( : )
REAL, ALLOCATABLE :: Recip_lambda_p_0 ( : )
REAL, ALLOCATABLE :: Recip_lambda_p_p ( : )
REAL, ALLOCATABLE :: Recip_lambda_p_p2( : )
REAL, ALLOCATABLE :: Recip_lambda_u_m ( : )
REAL, ALLOCATABLE :: Recip_lambda_u_0 ( : )
REAL, ALLOCATABLE :: Recip_lambda_u_p ( : )
REAL, ALLOCATABLE :: Recip_lambda_u_p2( : )
REAL, ALLOCATABLE :: Recip_phi_p_m ( : , : )
REAL, ALLOCATABLE :: Recip_phi_p_0 ( : , : )
REAL, ALLOCATABLE :: Recip_phi_p_p ( : , : )
REAL, ALLOCATABLE :: Recip_phi_p_p2( : , : )
REAL, ALLOCATABLE :: Recip_phi_v_m ( : , : )
REAL, ALLOCATABLE :: Recip_phi_v_0 ( : , : )
REAL, ALLOCATABLE :: Recip_phi_v_p ( : , : )
REAL, ALLOCATABLE :: Recip_phi_v_p2( : , : )

END MODULE dyn_var_res_Mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Constants necessary for East-West coefficients in diffusion scheme.

MODULE diff_coeff_Mod

! Description:
! This module is used to hold diffusion coefficients set in atmos_init.
!
! Method:
! The diff_coeff_u/v arrays are calculated in routine control/top_level/
! set_diff called from atmos_init,
! and used in diffusion_and_filtering routines.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

! Arguments

! East-West diffusion coefficient Arrays

REAL, ALLOCATABLE  ::  diff_coeff_u (:,:)
REAL, ALLOCATABLE  ::  diff_coeff_v (:,:)

END MODULE diff_coeff_Mod

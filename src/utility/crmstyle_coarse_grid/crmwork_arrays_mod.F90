! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Holds work arrays for reading in full data from fields files

MODULE crmwork_arrays_mod

USE word_sizes_mod, ONLY: iwp, wp

IMPLICIT NONE
SAVE

! Description:
!   Module containing work arrays, use for interpolating in height toegther
!  with information on heights, orography and cos latitude.
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CRMstyle_coarse_grid
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3  programming standards.
!
! Declarations:

! Note these fields are only available on PE 0

! Full 64 bit arrays
REAL, ALLOCATABLE  ::            &
  xcoslat_full(:,:)                ! cos(latitude full grid)

REAL(wp), ALLOCATABLE  ::        &
  prec_full(:,:)                   ! total surface precipitation full grid
                                   ! (kg/m2/s or mm/s) 

REAL(wp), ALLOCATABLE  ::        &
  orog_full(:,:)                   ! orography (m)

REAL, ALLOCATABLE  ::            &
  r_theta_levels_local(:,:,:)    & ! r_theta local area
 ,r_rho_levels_local(:,:,:)        ! r_rho local area

REAL(wp), ALLOCATABLE  ::        &
  r_theta_sea2(:)                  ! radius of sea points 2nd copy
REAL(wp), ALLOCATABLE  ::        &
  r_theta_sea(:)                 & ! radius of sea points  
 ,r_rho_sea(:)                     ! radius of sea points rho levels 


! Interpolation in vertical
INTEGER(iwp), ALLOCATABLE  ::  &
  uv_km_level(:,:,:)           & ! level number below required fixed height
 ,th_km_level_full(:,:,:)      & ! level number below required fixed height
 ,th_km_level(:,:,:)             ! level number below required fixed height


REAL(wp), ALLOCATABLE  ::    &
  uv_weight(:,:,:)           & ! weights for interpolation from uv levels
 ,th_weight(:,:,:)           & ! weights for interpolation from theta levels
 ,th_weight_full(:,:,:)        ! weights for interpolation from theta levels

LOGICAL, ALLOCATABLE  :: &
  mask(:,:,:)              ! true if above surface

REAL(wp), ALLOCATABLE ::  &
  h_theta_sea(:)          &  ! Required output heights
 ,h_rho_sea(:)            &  ! Required rho heights over sea
 ,dz_theta_sea(:)         &  ! Thickness of output theta levels
 ,z_theta_full(:,:,:)     &  ! heights from surface of earth theta levels
 ,z_theta(:,:,:)          &  ! heights from surface of earth theta levels
 ,z_rho(:,:,:)               ! heights from surface of earth rho levels

INTEGER, ALLOCATABLE ::      &
  index_col(:)               & ! column on coarse grid
 ,index_row(:)                 ! row on coarse grid

END MODULE crmwork_arrays_mod

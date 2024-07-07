! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Holds hiresolution data arrays on fixed level set

MODULE hires_data_mod

USE word_sizes_mod, ONLY: wp, iwp

IMPLICIT NONE
SAVE

! Description:
!   Module containing high resolution data on fixed levels
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:CRMstyle_coarse_grid
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3  programming standards.
!
! Declarations:

INTEGER(iwp), ALLOCATABLE  ::      &
  iclass_col(:,:)                    ! column classification based on surface 
                                     ! precipitation 
                                     ! 0 - non-precipitating
                                     ! 1 - convective
                                     ! 2 - stratiform

REAL(wp), ALLOCATABLE  ::          &
  orog(:,:)                        & ! orography (m)
 ,landsea(:,:)                     & ! Land sea mask
 ,precip(:,:)                      & ! total instantaneous precip (kg/m2/s)
 ,rain(:,:)                        & ! instantaneous rain (kg/m2/s)
 ,snow(:,:)                        & ! instantaneous snow (kg/m2/s)
 ,zh(:,:)                          & ! boundary layer depth (m)
 ,lh(:,:)                          & ! latent heat
 ,sh(:,:)                          & ! sensible heat
 ,pstar(:,:)                       & ! surface pressure (Pa)
 ,tstar(:,:)                         ! surface temperature (K)

REAL(wp), ALLOCATABLE  ::              &
  theta(:,:,:)                     & ! potential temperature (K)
 ,thetav(:,:,:)                    & ! virtual potential temperature (K)
 ,t(:,:,:)                         & ! temperature (K)
 ,density(:,:,:)                   & ! density of air kg/m3
 ,rh(:,:,:)                        & ! relative humdity
 ,q(:,:,:)                         & ! water vapour (kg/kg)
 ,qcl(:,:,:)                       & ! cloud liquid water (kg/kg)
 ,qcf(:,:,:)                       & ! cloud ice water (kg/kg)
 ,qrain(:,:,:)                     & ! rain (kg/kg)
 ,qgraup(:,:,:)                    & ! graupel (kg/kg)
 ,p_theta_lev(:,:,:)               & ! pressure on theta levels (Pa)
 ,p_rho_lev(:,:,:)                 & ! pressure on rho levels (Pa)
 ,exner(:,:,:)                     & ! Exner on theta levels
 ,u(:,:,:)                         & ! U wind (m/s)
 ,v(:,:,:)                         & ! V wind (m/s)
 ,w(:,:,:)                           ! W wind (m/s)

REAL(wp), ALLOCATABLE  ::              &
  dpdx(:,:,:)                      & ! pressure gradient dp/dx
 ,dpdy(:,:,:)                        ! pressure gradient dp/dy

REAL(wp), ALLOCATABLE  :: &
  dt1(:,:,:)          & ! Temperature tendency from SW (K/timestep)
 ,dt2(:,:,:)          & ! Temperature tendency from LW (K/timestep)
 ,dt4(:,:,:)          & ! Temperature tendency from microphysics (K/timestep)
 ,dt9(:,:,:)          & ! Temperature tendency from BL & cld (K/timestep)
 ,dt12(:,:,:)         & ! Temperature tendency from advection (K/timestep)
 ,dt30(:,:,:)           ! Temperature tendency total (K/timestep)

REAL(wp), ALLOCATABLE  :: &
  dq4(:,:,:)          & ! water vapour tendency from microphysics (K/timestep)
 ,dq9(:,:,:)          & ! water vapour tendency from BL & cld (K/timestep)
 ,dq12(:,:,:)         & ! water vapour tendency from advection (K/timestep)
 ,dq30(:,:,:)           ! water vapour tendency total (K/timestep)

REAL(wp), ALLOCATABLE  :: &
  dqcl4(:,:,:)        & ! cloud water tendency from microphysics (K/timestep)
 ,dqcl9(:,:,:)        & ! cloud water tendency from BL & cld (K/timestep)
 ,dqcl12(:,:,:)       & ! cloud water tendency from advection (K/timestep)
 ,dqcl30(:,:,:)         ! cloud water tendency total (K/timestep)

REAL(wp), ALLOCATABLE  :: &
  dqcf4(:,:,:)        & ! cloud ice tendency from microphysics (K/timestep)
 ,dqcf3(:,:,:)        & ! cloud ice tendency from BL  (K/timestep)
 ,dqcf12(:,:,:)       & ! cloud ice tendency from advection (K/timestep)
 ,dqcf30(:,:,:)         ! cloud ice tendency total (K/timestep)

REAL(wp), ALLOCATABLE  :: &
  dqrain30(:,:,:)         & ! rain tendency total (kg/kg/timestep) 
 ,dqgr30(:,:,:)           & ! graupel tendency total (kg/kg/timestep) 
 ,drho(:,:,:)               ! density tendency total (kg/m3/timestep) 

END MODULE hires_data_mod

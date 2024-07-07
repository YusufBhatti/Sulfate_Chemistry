! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Holds arrays for coarse grids

MODULE crmstyle_sample_arrays_mod

! Description:
!  Arrays to hold CRM divisions

USE word_sizes_mod, ONLY: iwp,wp       ! Allows use of 4 byte words to reduce
                                        ! memory

IMPLICIT NONE
SAVE

! Description:
! Arrays for coarse grids
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utilty - crmstyle_coarse_grid
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

!----------------------------------------------------------------------
! Series of possible grids - only nres used
!----------------------------------------------------------------------
! Single level fields

REAL(wp), ALLOCATABLE  ::  &
  all_zh(:,:)          & ! mean zh
 ,all_sh(:,:)          & ! mean sh
 ,all_lh(:,:)          & ! mean lh
 ,all_pstar(:,:)       & ! mean pstar
 ,all_tstar(:,:)       & ! mean tstar
 ,all_rain(:,:)        & ! mean rain
 ,all_snow(:,:)        & ! mean snow
 ,all_precip(:,:)      & ! mean precip
 ,all_orog(:,:)        & ! mean orog
 ,all_land(:,:)          ! mean land frac

REAL(wp), ALLOCATABLE  ::  &
  all_sd_zh(:,:)            & ! sd zh
 ,all_sd_sh(:,:)            & ! sd sh
 ,all_sd_lh(:,:)            & ! sd lh
 ,all_sd_pstar(:,:)         & ! sd pstar
 ,all_sd_tstar(:,:)         & ! sd tstar
 ,all_sd_rain(:,:)          & ! sd rain
 ,all_sd_snow(:,:)          & ! sd snow
 ,all_sd_precip(:,:)        & ! sd precip
 ,all_sd_orog(:,:)          & ! sd orog
 ,all_sd_land(:,:)            ! sd land frac

! Column classification means
REAL(wp), ALLOCATABLE  ::  &
  fract_conv(:,:)          & ! Fraction of region classed as convective
 ,fract_strat(:,:)         & ! Fraction of region classed as stratiform
 ,prec_conv(:,:)           & ! prec classed as convective
 ,prec_strat(:,:)            ! prec classed as stratiform

!----------------------------------------------------------------------
! CAPE and CIN from mean profiles
!----------------------------------------------------------------------
REAL(wp), ALLOCATABLE  ::  &
  cape(:,:)                & ! Convectively available potential energy (J/kg) 
 ,cin(:,:)                 & ! Convective Inhibition (J/kg)
 ,zneutral(:,:)            & ! height of neutral buoyancy (m)
 ,zfree(:,:)               & ! Level of free convection (m)
 ,zlcl(:,:)                  ! Lifting condensation level above sea level(m)

!----------------------------------------------------------------------
! PCAPE and changes in PCAPE 
!----------------------------------------------------------------------
REAL(wp), ALLOCATABLE  ::  &
  bcu_pcape(:,:)           & ! Convectively available potential energy 
 ,bcu_dpcapedt(:,:)        & ! Rate of change of PCAPE 
 ,bcu_dpcapedt_bl(:,:)     & ! Rate of change of PCAPE in BL
 ,bcw_pcape(:,:)           & ! Convectively available potential energy 
 ,bcw_dpcapedt(:,:)        & ! Rate of change of PCAPE 
 ,bcw_dpcapedt_bl(:,:)     & ! Total rate of change of Tv in BL
 ,ppd_dpcapedt_bl(:,:)       ! Rate of change of Tv in DD part of BL

! 3d fields used by CAPE, CIN and PCAPE calculations
REAL(wp), ALLOCATABLE  ::  &
  all_exner(:,:,:)         & ! mean exner
 ,all_tv(:,:,:)              ! mean virtual temperature (K)

!----------------------------------------------------------------------
! dilute CAPE, rate of change in dilute CAPE over convective cloud
!----------------------------------------------------------------------
REAL(wp), ALLOCATABLE  ::  &
  bcu_dilcape(:,:)         & ! Dilute convectively available potential energy 
 ,bcu_dcapedt(:,:)         & ! Rate of change of CAPE over BCU
 ,bcw_dilcape(:,:)         & ! Dilute convectively available potential energy 
 ,bcw_dcapedt(:,:)           ! Rate of change of CAPE over BCW

!----------------------------------------------------------------------
! Grid 1
REAL(wp), ALLOCATABLE  ::     &
 p_theta_hydro(:,:,:)          ! Hydrostatic pressure for mean profile

REAL(wp), ALLOCATABLE  ::     &
  all_th(:,:,:)       & ! mean theta (K)
 ,all_thv(:,:,:)      & ! mean thetav (K)
 ,all_t(:,:,:)        & ! mean temperature (K)
 ,all_q(:,:,:)        & ! mean q (kg/kg)
 ,all_qcl(:,:,:)      & ! mean qcl (kg/kg)
 ,all_qcf(:,:,:)      & ! mean qcf (kg/kg)
 ,all_qrain(:,:,:)    & ! mean rain (kg/kg)
 ,all_qgraup(:,:,:)   & ! mean qraupel (kg/kg)
 ,all_u(:,:,:)        & ! u wind (m/s)
 ,all_v(:,:,:)        & ! v wind (m/s)
 ,all_w(:,:,:)        & ! vertical velocity (m/s)
 ,all_ptheta(:,:,:)   & ! pressure theta levels (Pa)
 ,all_a(:,:,:)        & ! fraction of area above surface
 ,all_rh(:,:,:)       & ! relative humidity
 ,all_rho(:,:,:)        ! density

! Mean tendencies
REAL(wp), ALLOCATABLE  ::     &
  all_dt1(:,:,:)      & ! dT from SW radiation K/timestep
 ,all_dt2(:,:,:)      & ! dT from LW radiation K/timestep
 ,all_dt4(:,:,:)      & ! dT from microphysics K/timestep
 ,all_dt9(:,:,:)      & ! dT from BL & cloud K/timestep
 ,all_dt12(:,:,:)     & ! dT from advection K/timestep
 ,all_dt30(:,:,:)     & ! dT from total K/timestep
 ,all_dq4(:,:,:)      & ! dq from microphysics K/timestep
 ,all_dq9(:,:,:)      & ! dq from BL & cloud K/timestep
 ,all_dq12(:,:,:)     & ! dq from advection K/timestep
 ,all_dq30(:,:,:)     & ! dq from total K/timestep
 ,all_dqcl4(:,:,:)    & ! dqcl from microphysics K/timestep
 ,all_dqcl9(:,:,:)    & ! dqcl from BL & cloud K/timestep
 ,all_dqcl12(:,:,:)   & ! dqcl from advection K/timestep
 ,all_dqcl30(:,:,:)   & ! dqcl from total K/timestep
 ,all_dqcf4(:,:,:)    & ! dqcf from microphysics K/timestep
 ,all_dqcf3(:,:,:)    & ! dqcf from BL & cloud K/timestep
 ,all_dqcf12(:,:,:)   & ! dqcf from advection K/timestep
 ,all_dqcf30(:,:,:)   & ! dqcf from total kg/kg/timestep
 ,all_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,all_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,all_drho(:,:,:)       ! drho from total kg/m3/timestep

REAL(wp), ALLOCATABLE  ::     &
  all_sd_th(:,:,:)       & ! sd theta (K)
 ,all_sd_thv(:,:,:)      & ! sd thetav (K)
 ,all_sd_t(:,:,:)        & ! sd t (K)
 ,all_sd_q(:,:,:)        & ! sd q (kg/kg)
 ,all_sd_qcl(:,:,:)      & ! sd qcl (kg/kg)
 ,all_sd_qcf(:,:,:)      & ! sd qcf (kg/kg)
 ,all_sd_qrain(:,:,:)    & ! sd rain (kg/kg)
 ,all_sd_qgraup(:,:,:)   & ! sd qraupel (kg/kg)
 ,all_sd_u(:,:,:)        & ! u wind (m/s)
 ,all_sd_v(:,:,:)        & ! v wind (m/s)
 ,all_sd_w(:,:,:)        & ! vertical velocity (m/s)
 ,all_sd_ptheta(:,:,:)   & ! pressure theta levels (Pa)
 ,all_sd_rh(:,:,:)       & ! relative humidity
 ,all_sd_rho(:,:,:)        ! density

! sd tendencies
REAL(wp), ALLOCATABLE  ::     &
  all_sd_dt1(:,:,:)      & ! dT from SW radiation K/timestep
 ,all_sd_dt2(:,:,:)      & ! dT from LW radiation K/timestep
 ,all_sd_dt4(:,:,:)      & ! dT from microphysics K/timestep
 ,all_sd_dt9(:,:,:)      & ! dT from BL & cloud K/timestep
 ,all_sd_dt12(:,:,:)     & ! dT from advection K/timestep
 ,all_sd_dt30(:,:,:)     & ! dT from total K/timestep
 ,all_sd_dq4(:,:,:)      & ! dq from microphysics K/timestep
 ,all_sd_dq9(:,:,:)      & ! dq from BL & cloud K/timestep
 ,all_sd_dq12(:,:,:)     & ! dq from advection K/timestep
 ,all_sd_dq30(:,:,:)     & ! dq from total K/timestep
 ,all_sd_dqcl4(:,:,:)    & ! dqcl from microphysics K/timestep
 ,all_sd_dqcl9(:,:,:)    & ! dqcl from BL & cloud K/timestep
 ,all_sd_dqcl12(:,:,:)   & ! dqcl from advection K/timestep
 ,all_sd_dqcl30(:,:,:)   & ! dqcl from total K/timestep
 ,all_sd_dqcf4(:,:,:)    & ! dqcf from microphysics K/timestep
 ,all_sd_dqcf3(:,:,:)    & ! dqcf from BL & cloud K/timestep
 ,all_sd_dqcf12(:,:,:)   & ! dqcf from advection K/timestep
 ,all_sd_dqcf30(:,:,:)   & ! dqcf from total K/timestep
 ,all_sd_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,all_sd_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,all_sd_drho(:,:,:)       ! drho from total kg/m3/timestep

! Mean products
REAL(wp), ALLOCATABLE  ::     &
  all_thw(:,:,:)      & ! theta'w'
 ,all_thvw(:,:,:)     & ! thetav'w'
 ,all_uth(:,:,:)      & ! u'theta'
 ,all_uthv(:,:,:)     & ! u'thetav'
 ,all_vth(:,:,:)      & ! v'theta'
 ,all_vthv(:,:,:)     & ! v'thetav'
 ,all_uq(:,:,:)       & ! u'q'
 ,all_vq(:,:,:)       & ! v'q'
 ,all_wp(:,:,:)       & ! w'p'/density
 ,all_wp_hydro(:,:,:) & ! w'p'/density
 ,all_qw(:,:,:)       & ! q'w'
 ,all_qclw(:,:,:)     & ! qcl'w'
 ,all_qcfw(:,:,:)     & ! qcf'w'
 ,all_qrainw(:,:,:)   & ! qrain'w'
 ,all_qgraupw(:,:,:)  & ! qgraup'w'
 ,all_uw(:,:,:)       & ! u'w'
 ,all_vw(:,:,:)       & ! v'w'
 ,all_ww(:,:,:)       & ! w'w'
 ,all_w3(:,:,:)       & ! w'w'w'
 ,all_vv(:,:,:)       & ! v'v'
 ,all_uu(:,:,:)       & ! u'u'
 ,all_uv(:,:,:)       & ! u'v
 ,all_dpx(:,:,:)      & !
 ,all_dpy(:,:,:)

! All cloudy points

REAL(wp), ALLOCATABLE ::    &
  acc_h_w(:,:,:)        &
 ,acc_h_th(:,:,:)       &
 ,acc_h_thv(:,:,:)      &
 ,acc_h_rho(:,:,:)      &
 ,acc_h_u(:,:,:)        &
 ,acc_h_v(:,:,:)        &
 ,acc_h_rh(:,:,:)       &
 ,acc_h_dt1(:,:,:)      &
 ,acc_h_dt2(:,:,:)      &
 ,acc_h_dt4(:,:,:)      &
 ,acc_h_dt9(:,:,:)      &
 ,acc_h_dt12(:,:,:)     &
 ,acc_h_dq4(:,:,:)      &
 ,acc_h_dq9(:,:,:)      &
 ,acc_h_dq12(:,:,:)     &
 ,acc_h_dqcl4(:,:,:)    &
 ,acc_h_dqcl9(:,:,:)    &
 ,acc_h_dqcl12(:,:,:)   &
 ,acc_h_dqcf4(:,:,:)    &
 ,acc_h_dqcf3(:,:,:)    &
 ,acc_h_dqcf12(:,:,:)   &
 ,acc_h_dt30(:,:,:)     &
 ,acc_h_dq30(:,:,:)     &
 ,acc_h_dqcl30(:,:,:)   &
 ,acc_h_dqcf30(:,:,:)   &
 ,acc_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,acc_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,acc_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,acc_h_q(:,:,:)        &
 ,acc_h_qcl(:,:,:)      &
 ,acc_h_qcf(:,:,:)      &
 ,acc_h_qrain(:,:,:)    &
 ,acc_h_qgraup(:,:,:)   &
 ,acc_h_a(:,:,:)        & ! fraction of columns in divison
 ,acc_h_thw(:,:,:)      & ! theta'w'
 ,acc_h_thvw(:,:,:)     & ! thetav'w'
 ,acc_h_qw(:,:,:)       & ! q'w'
 ,acc_h_qclw(:,:,:)     & ! qcl'w'
 ,acc_h_qcfw(:,:,:)     & ! qcf'w'
 ,acc_h_qrainw(:,:,:)   & ! qrain'w'
 ,acc_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,acc_h_uw(:,:,:)       & ! u'w'
 ,acc_h_vw(:,:,:)       & ! v'w'
 ,acc_h_ww(:,:,:)       & ! w'w'
 ,acc_h_vv(:,:,:)       & ! v'v'
 ,acc_h_uu(:,:,:)       & ! u'u'
 ,acc_h_uv(:,:,:)       & ! u'v'
 ,acc_h_w3(:,:,:)       & ! w'w'w'
 ,acc_h_dpx(:,:,:)      & ! dp/dx
 ,acc_h_dpy(:,:,:)      & ! dp/dy
 ,acc_h_uth(:,:,:)      & ! u'theta'
 ,acc_h_uthv(:,:,:)     & ! u'thetav'
 ,acc_h_vth(:,:,:)      & ! v'theta'
 ,acc_h_vthv(:,:,:)     & ! v'thetav'
 ,acc_h_uq(:,:,:)       & ! u'q'
 ,acc_h_vq(:,:,:)       & ! v'q'
 ,acc_h_wp(:,:,:)         ! w'p'/density

! All cloudy updraughts

REAL(wp), ALLOCATABLE ::    &
  acu_h_w(:,:,:)        &
 ,acu_h_th(:,:,:)       &
 ,acu_h_thv(:,:,:)      &
 ,acu_h_rho(:,:,:)      &
 ,acu_h_u(:,:,:)        &
 ,acu_h_v(:,:,:)        &
 ,acu_h_rh(:,:,:)       &
 ,acu_h_dt1(:,:,:)      &
 ,acu_h_dt2(:,:,:)      &
 ,acu_h_dt4(:,:,:)      &
 ,acu_h_dt9(:,:,:)      &
 ,acu_h_dt12(:,:,:)     &
 ,acu_h_dq4(:,:,:)      &
 ,acu_h_dq9(:,:,:)      &
 ,acu_h_dq12(:,:,:)     &
 ,acu_h_dqcl4(:,:,:)    &
 ,acu_h_dqcl9(:,:,:)    &
 ,acu_h_dqcl12(:,:,:)   &
 ,acu_h_dqcf4(:,:,:)    &
 ,acu_h_dqcf3(:,:,:)    &
 ,acu_h_dqcf12(:,:,:)   &
 ,acu_h_dt30(:,:,:)     &
 ,acu_h_dq30(:,:,:)     &
 ,acu_h_dqcl30(:,:,:)   &
 ,acu_h_dqcf30(:,:,:)   &
 ,acu_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,acu_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,acu_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,acu_h_q(:,:,:)        &
 ,acu_h_qcl(:,:,:)      &
 ,acu_h_qcf(:,:,:)      &
 ,acu_h_qrain(:,:,:)    &
 ,acu_h_qgraup(:,:,:)   &
 ,acu_h_a(:,:,:)        & ! fraction of columns in divison
 ,acu_h_thw(:,:,:)      & ! theta'w'
 ,acu_h_thvw(:,:,:)     & ! thetav'w'
 ,acu_h_qw(:,:,:)       & ! q'w'
 ,acu_h_qclw(:,:,:)     & ! qcl'w'
 ,acu_h_qcfw(:,:,:)     & ! qcf'w'
 ,acu_h_qrainw(:,:,:)   & ! qrain'w'
 ,acu_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,acu_h_uw(:,:,:)       & ! u'w'
 ,acu_h_vw(:,:,:)       & ! v'w'
 ,acu_h_ww(:,:,:)       & ! w'w'
 ,acu_h_vv(:,:,:)       & ! v'v'
 ,acu_h_uu(:,:,:)       & ! u'u'
 ,acu_h_uv(:,:,:)       & ! u'v'
 ,acu_h_w3(:,:,:)       & ! w'w'w'
 ,acu_h_dpx(:,:,:)      & ! dp/dx
 ,acu_h_dpy(:,:,:)      & ! dp/dy
 ,acu_h_uth(:,:,:)      & ! u'theta'
 ,acu_h_uthv(:,:,:)     & ! u'thetav'
 ,acu_h_vth(:,:,:)      & ! v'theta'
 ,acu_h_vthv(:,:,:)     & ! v'thetav'
 ,acu_h_uq(:,:,:)       & ! u'q'
 ,acu_h_vq(:,:,:)       & ! v'q'
 ,acu_h_wp(:,:,:)         ! w'p'/density

! Buoyant cloudy updraughts  - relative to mean w

REAL(wp), ALLOCATABLE ::    &
  bcu_h_w(:,:,:)        &
 ,bcu_h_th(:,:,:)       &
 ,bcu_h_thv(:,:,:)      &
 ,bcu_h_rho(:,:,:)      &
 ,bcu_h_u(:,:,:)        &
 ,bcu_h_v(:,:,:)        &
 ,bcu_h_rh(:,:,:)       &
 ,bcu_h_dt1(:,:,:)      &
 ,bcu_h_dt2(:,:,:)      &
 ,bcu_h_dt4(:,:,:)      &
 ,bcu_h_dt9(:,:,:)      &
 ,bcu_h_dt12(:,:,:)     &
 ,bcu_h_dq4(:,:,:)      &
 ,bcu_h_dq9(:,:,:)      &
 ,bcu_h_dq12(:,:,:)     &
 ,bcu_h_dqcl4(:,:,:)    &
 ,bcu_h_dqcl9(:,:,:)    &
 ,bcu_h_dqcl12(:,:,:)   &
 ,bcu_h_dqcf4(:,:,:)    &
 ,bcu_h_dqcf3(:,:,:)    &
 ,bcu_h_dqcf12(:,:,:)   &
 ,bcu_h_dt30(:,:,:)     &
 ,bcu_h_dq30(:,:,:)     &
 ,bcu_h_dqcl30(:,:,:)   &
 ,bcu_h_dqcf30(:,:,:)   &
 ,bcu_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,bcu_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,bcu_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,bcu_h_q(:,:,:)        &
 ,bcu_h_qcl(:,:,:)      &
 ,bcu_h_qcf(:,:,:)      &
 ,bcu_h_qrain(:,:,:)    &
 ,bcu_h_qgraup(:,:,:)   &
 ,bcu_h_a(:,:,:)        & ! fraction of columns in divison
 ,bcu_h_thw(:,:,:)      & ! theta'w'
 ,bcu_h_thvw(:,:,:)     & ! thetav'w'
 ,bcu_h_qw(:,:,:)       & ! q'w'
 ,bcu_h_qclw(:,:,:)     & ! qcl'w'
 ,bcu_h_qcfw(:,:,:)     & ! qcf'w'
 ,bcu_h_qrainw(:,:,:)   & ! qrain'w'
 ,bcu_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,bcu_h_uw(:,:,:)       & ! u'w'
 ,bcu_h_vw(:,:,:)       & ! v'w'
 ,bcu_h_ww(:,:,:)       & ! w'w'
 ,bcu_h_vv(:,:,:)       & ! v'v'
 ,bcu_h_uu(:,:,:)       & ! u'u'
 ,bcu_h_uv(:,:,:)       & ! u'v'
 ,bcu_h_w3(:,:,:)       & ! w'w'w'
 ,bcu_h_dpx(:,:,:)      & ! dp/dx
 ,bcu_h_dpy(:,:,:)      & ! dp/dy
 ,bcu_h_uth(:,:,:)      & ! u'theta'
 ,bcu_h_uthv(:,:,:)     & ! u'thetav'
 ,bcu_h_vth(:,:,:)      & ! v'theta'
 ,bcu_h_vthv(:,:,:)     & ! v'thetav'
 ,bcu_h_uq(:,:,:)       & ! u'q'
 ,bcu_h_vq(:,:,:)       & ! v'q'
 ,bcu_h_wp(:,:,:)         ! w'p'/density

! Strong buoyant updraughts (>1)  - relative to mean w

REAL(wp), ALLOCATABLE ::    &
  wg1_h_w(:,:,:)        &
 ,wg1_h_th(:,:,:)       &
 ,wg1_h_thv(:,:,:)      &
 ,wg1_h_rho(:,:,:)      &
 ,wg1_h_u(:,:,:)        &
 ,wg1_h_v(:,:,:)        &
 ,wg1_h_rh(:,:,:)       &
 ,wg1_h_dt1(:,:,:)      &
 ,wg1_h_dt2(:,:,:)      &
 ,wg1_h_dt4(:,:,:)      &
 ,wg1_h_dt9(:,:,:)      &
 ,wg1_h_dt12(:,:,:)     &
 ,wg1_h_dq4(:,:,:)      &
 ,wg1_h_dq9(:,:,:)      &
 ,wg1_h_dq12(:,:,:)     &
 ,wg1_h_dqcl4(:,:,:)    &
 ,wg1_h_dqcl9(:,:,:)    &
 ,wg1_h_dqcl12(:,:,:)   &
 ,wg1_h_dqcf4(:,:,:)    &
 ,wg1_h_dqcf3(:,:,:)    &
 ,wg1_h_dqcf12(:,:,:)   &
 ,wg1_h_dt30(:,:,:)     &
 ,wg1_h_dq30(:,:,:)     &
 ,wg1_h_dqcl30(:,:,:)   &
 ,wg1_h_dqcf30(:,:,:)   &
 ,wg1_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,wg1_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,wg1_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,wg1_h_q(:,:,:)        &
 ,wg1_h_qcl(:,:,:)      &
 ,wg1_h_qcf(:,:,:)      &
 ,wg1_h_qrain(:,:,:)    &
 ,wg1_h_qgraup(:,:,:)   &
 ,wg1_h_a(:,:,:)        & ! fraction of columns in divison
 ,wg1_h_thw(:,:,:)      & ! theta'w'
 ,wg1_h_thvw(:,:,:)     & ! thetav'w'
 ,wg1_h_qw(:,:,:)       & ! q'w'
 ,wg1_h_qclw(:,:,:)     & ! qcl'w'
 ,wg1_h_qcfw(:,:,:)     & ! qcf'w'
 ,wg1_h_qrainw(:,:,:)   & ! qrain'w'
 ,wg1_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,wg1_h_uw(:,:,:)       & ! u'w'
 ,wg1_h_vw(:,:,:)       & ! v'w'
 ,wg1_h_ww(:,:,:)       & ! w'w'
 ,wg1_h_vv(:,:,:)       & ! v'v'
 ,wg1_h_uu(:,:,:)       & ! u'u'
 ,wg1_h_uv(:,:,:)       & ! u'v'
 ,wg1_h_w3(:,:,:)       & ! w'w'w'
 ,wg1_h_dpx(:,:,:)      & ! dp/dx
 ,wg1_h_dpy(:,:,:)      & ! dp/dy
 ,wg1_h_uth(:,:,:)      & ! u'theta'
 ,wg1_h_uthv(:,:,:)     & ! u'thetav'
 ,wg1_h_vth(:,:,:)      & ! v'theta'
 ,wg1_h_vthv(:,:,:)     & ! v'thetav'
 ,wg1_h_uq(:,:,:)       & ! u'q'
 ,wg1_h_vq(:,:,:)       & ! v'q'
 ,wg1_h_wp(:,:,:)         ! w'p'/density

! precipitation downdraughts

REAL(wp), ALLOCATABLE ::    &
  ppd_h_w(:,:,:)        &
 ,ppd_h_th(:,:,:)       &
 ,ppd_h_thv(:,:,:)      &
 ,ppd_h_rho(:,:,:)      &
 ,ppd_h_u(:,:,:)        &
 ,ppd_h_v(:,:,:)        &
 ,ppd_h_rh(:,:,:)       &
 ,ppd_h_dt1(:,:,:)      &
 ,ppd_h_dt2(:,:,:)      &
 ,ppd_h_dt4(:,:,:)      &
 ,ppd_h_dt9(:,:,:)      &
 ,ppd_h_dt12(:,:,:)     &
 ,ppd_h_dq4(:,:,:)      &
 ,ppd_h_dq9(:,:,:)      &
 ,ppd_h_dq12(:,:,:)     &
 ,ppd_h_dqcl4(:,:,:)    &
 ,ppd_h_dqcl9(:,:,:)    &
 ,ppd_h_dqcl12(:,:,:)   &
 ,ppd_h_dqcf4(:,:,:)    &
 ,ppd_h_dqcf3(:,:,:)    &
 ,ppd_h_dqcf12(:,:,:)   &
 ,ppd_h_dt30(:,:,:)     &
 ,ppd_h_dq30(:,:,:)     &
 ,ppd_h_dqcl30(:,:,:)   &
 ,ppd_h_dqcf30(:,:,:)   &
 ,ppd_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,ppd_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,ppd_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,ppd_h_q(:,:,:)        &
 ,ppd_h_qcl(:,:,:)      &
 ,ppd_h_qcf(:,:,:)      &
 ,ppd_h_qrain(:,:,:)    &
 ,ppd_h_qgraup(:,:,:)   &
 ,ppd_h_a(:,:,:)        & ! fraction of columns in divison
 ,ppd_h_thw(:,:,:)      & ! theta'w'
 ,ppd_h_thvw(:,:,:)     & ! thetav'w'
 ,ppd_h_qw(:,:,:)       & ! q'w'
 ,ppd_h_qclw(:,:,:)     & ! qcl'w'
 ,ppd_h_qcfw(:,:,:)     & ! qcf'w'
 ,ppd_h_qrainw(:,:,:)   & ! qrain'w'
 ,ppd_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,ppd_h_uw(:,:,:)       & ! u'w'
 ,ppd_h_vw(:,:,:)       & ! v'w'
 ,ppd_h_ww(:,:,:)       & ! w'w'
 ,ppd_h_vv(:,:,:)       & ! v'v'
 ,ppd_h_uu(:,:,:)       & ! u'u'
 ,ppd_h_uv(:,:,:)       & ! u'v'
 ,ppd_h_w3(:,:,:)       & ! w'w'w'
 ,ppd_h_dpx(:,:,:)      & ! dp/dx
 ,ppd_h_dpy(:,:,:)      & ! dp/dy
 ,ppd_h_uth(:,:,:)      & ! u'theta'
 ,ppd_h_uthv(:,:,:)     & ! u'thetav'
 ,ppd_h_vth(:,:,:)      & ! v'theta'
 ,ppd_h_vthv(:,:,:)     & ! v'thetav'
 ,ppd_h_uq(:,:,:)       & ! u'q'
 ,ppd_h_vq(:,:,:)       & ! v'q'
 ,ppd_h_wp(:,:,:)         ! w'p'/density


! Negatively buoyant precipitation downdraughts

REAL(wp), ALLOCATABLE ::    &
  nbd_h_w(:,:,:)        &
 ,nbd_h_th(:,:,:)       &
 ,nbd_h_thv(:,:,:)      &
 ,nbd_h_rho(:,:,:)      &
 ,nbd_h_u(:,:,:)        &
 ,nbd_h_v(:,:,:)        &
 ,nbd_h_rh(:,:,:)       &
 ,nbd_h_dt1(:,:,:)      &
 ,nbd_h_dt2(:,:,:)      &
 ,nbd_h_dt4(:,:,:)      &
 ,nbd_h_dt9(:,:,:)      &
 ,nbd_h_dt12(:,:,:)     &
 ,nbd_h_dq4(:,:,:)      &
 ,nbd_h_dq9(:,:,:)      &
 ,nbd_h_dq12(:,:,:)     &
 ,nbd_h_dqcl4(:,:,:)    &
 ,nbd_h_dqcl9(:,:,:)    &
 ,nbd_h_dqcl12(:,:,:)   &
 ,nbd_h_dqcf4(:,:,:)    &
 ,nbd_h_dqcf3(:,:,:)    &
 ,nbd_h_dqcf12(:,:,:)   &
 ,nbd_h_dt30(:,:,:)     &
 ,nbd_h_dq30(:,:,:)     &
 ,nbd_h_dqcl30(:,:,:)   &
 ,nbd_h_dqcf30(:,:,:)   &
 ,nbd_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,nbd_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,nbd_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,nbd_h_q(:,:,:)        &
 ,nbd_h_qcl(:,:,:)      &
 ,nbd_h_qcf(:,:,:)      &
 ,nbd_h_qrain(:,:,:)    &
 ,nbd_h_qgraup(:,:,:)   &
 ,nbd_h_a(:,:,:)        & ! fraction of columns in divison
 ,nbd_h_thw(:,:,:)      & ! theta'w'
 ,nbd_h_thvw(:,:,:)     & ! thetav'w'
 ,nbd_h_qw(:,:,:)       & ! q'w'
 ,nbd_h_qclw(:,:,:)     & ! qcl'w'
 ,nbd_h_qcfw(:,:,:)     & ! qcf'w'
 ,nbd_h_qrainw(:,:,:)   & ! qrain'w'
 ,nbd_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,nbd_h_uw(:,:,:)       & ! u'w'
 ,nbd_h_vw(:,:,:)       & ! v'w'
 ,nbd_h_ww(:,:,:)       & ! w'w'
 ,nbd_h_vv(:,:,:)       & ! v'v'
 ,nbd_h_uu(:,:,:)       & ! u'u'
 ,nbd_h_uv(:,:,:)       & ! u'v'
 ,nbd_h_w3(:,:,:)       & ! w'w'w'
 ,nbd_h_dpx(:,:,:)      & ! dp/dx
 ,nbd_h_dpy(:,:,:)      & ! dp/dy
 ,nbd_h_uth(:,:,:)      & ! u'theta'
 ,nbd_h_uthv(:,:,:)     & ! u'thetav'
 ,nbd_h_vth(:,:,:)      & ! v'theta'
 ,nbd_h_vthv(:,:,:)     & ! v'thetav'
 ,nbd_h_uq(:,:,:)       & ! u'q'
 ,nbd_h_vq(:,:,:)       & ! v'q'
 ,nbd_h_wp(:,:,:)         ! w'p'/density

!------------------------------------------------------------------------------
! Negatively buoyant precipitation with ice downdraughts
!------------------------------------------------------------------------------

REAL(wp), ALLOCATABLE ::    &
  nid_h_w(:,:,:)        &
 ,nid_h_th(:,:,:)       &
 ,nid_h_thv(:,:,:)      &
 ,nid_h_rho(:,:,:)      &
 ,nid_h_u(:,:,:)        &
 ,nid_h_v(:,:,:)        &
 ,nid_h_rh(:,:,:)       &
 ,nid_h_dt1(:,:,:)      &
 ,nid_h_dt2(:,:,:)      &
 ,nid_h_dt4(:,:,:)      &
 ,nid_h_dt9(:,:,:)      &
 ,nid_h_dt12(:,:,:)     &
 ,nid_h_dq4(:,:,:)      &
 ,nid_h_dq9(:,:,:)      &
 ,nid_h_dq12(:,:,:)     &
 ,nid_h_dqcl4(:,:,:)    &
 ,nid_h_dqcl9(:,:,:)    &
 ,nid_h_dqcl12(:,:,:)   &
 ,nid_h_dqcf4(:,:,:)    &
 ,nid_h_dqcf3(:,:,:)    &
 ,nid_h_dqcf12(:,:,:)   &
 ,nid_h_dt30(:,:,:)     &
 ,nid_h_dq30(:,:,:)     &
 ,nid_h_dqcl30(:,:,:)   &
 ,nid_h_dqcf30(:,:,:)   &
 ,nid_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,nid_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,nid_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,nid_h_q(:,:,:)        &
 ,nid_h_qcl(:,:,:)      &
 ,nid_h_qcf(:,:,:)      &
 ,nid_h_qrain(:,:,:)    &
 ,nid_h_qgraup(:,:,:)   &
 ,nid_h_a(:,:,:)        & ! fraction of columns in divison
 ,nid_h_thw(:,:,:)      & ! theta'w'
 ,nid_h_thvw(:,:,:)     & ! thetav'w'
 ,nid_h_qw(:,:,:)       & ! q'w'
 ,nid_h_qclw(:,:,:)     & ! qcl'w'
 ,nid_h_qcfw(:,:,:)     & ! qcf'w'
 ,nid_h_qrainw(:,:,:)   & ! qrain'w'
 ,nid_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,nid_h_uw(:,:,:)       & ! u'w'
 ,nid_h_vw(:,:,:)       & ! v'w'
 ,nid_h_ww(:,:,:)       & ! w'w'
 ,nid_h_vv(:,:,:)       & ! v'v'
 ,nid_h_uu(:,:,:)       & ! u'u'
 ,nid_h_uv(:,:,:)       & ! u'v'
 ,nid_h_w3(:,:,:)       & ! w'w'w'
 ,nid_h_dpx(:,:,:)      & ! dp/dx
 ,nid_h_dpy(:,:,:)      & ! dp/dy
 ,nid_h_uth(:,:,:)      & ! u'theta'
 ,nid_h_uthv(:,:,:)     & ! u'thetav'
 ,nid_h_vth(:,:,:)      & ! v'theta'
 ,nid_h_vthv(:,:,:)     & ! v'thetav'
 ,nid_h_uq(:,:,:)       & ! u'q'
 ,nid_h_vq(:,:,:)       & ! v'q'
 ,nid_h_wp(:,:,:)         ! w'p'/density

!------------------------------------------------------------------------------
! ADU - Dry upward air  (where upward in not relative)
!------------------------------------------------------------------------------

REAL(wp), ALLOCATABLE ::    &
  adu_h_w(:,:,:)        &
 ,adu_h_th(:,:,:)       &
 ,adu_h_thv(:,:,:)      &
 ,adu_h_rho(:,:,:)      &
 ,adu_h_u(:,:,:)        &
 ,adu_h_v(:,:,:)        &
 ,adu_h_rh(:,:,:)       &
 ,adu_h_dt1(:,:,:)      &
 ,adu_h_dt2(:,:,:)      &
 ,adu_h_dt4(:,:,:)      &
 ,adu_h_dt9(:,:,:)      &
 ,adu_h_dt12(:,:,:)     &
 ,adu_h_dq4(:,:,:)      &
 ,adu_h_dq9(:,:,:)      &
 ,adu_h_dq12(:,:,:)     &
 ,adu_h_dqcl4(:,:,:)    &
 ,adu_h_dqcl9(:,:,:)    &
 ,adu_h_dqcl12(:,:,:)   &
 ,adu_h_dqcf4(:,:,:)    &
 ,adu_h_dqcf3(:,:,:)    &
 ,adu_h_dqcf12(:,:,:)   &
 ,adu_h_dt30(:,:,:)     &
 ,adu_h_dq30(:,:,:)     &
 ,adu_h_dqcl30(:,:,:)   &
 ,adu_h_dqcf30(:,:,:)   &
 ,adu_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,adu_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,adu_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,adu_h_q(:,:,:)        &
 ,adu_h_qcl(:,:,:)      &
 ,adu_h_qcf(:,:,:)      &
 ,adu_h_qrain(:,:,:)    &
 ,adu_h_qgraup(:,:,:)   &
 ,adu_h_a(:,:,:)        & ! fraction of columns in divison
 ,adu_h_thw(:,:,:)      & ! theta'w'
 ,adu_h_thvw(:,:,:)     & ! thetav'w'
 ,adu_h_qw(:,:,:)       & ! q'w'
 ,adu_h_qclw(:,:,:)     & ! qcl'w'
 ,adu_h_qcfw(:,:,:)     & ! qcf'w'
 ,adu_h_qrainw(:,:,:)   & ! qrain'w'
 ,adu_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,adu_h_uw(:,:,:)       & ! u'w'
 ,adu_h_vw(:,:,:)       & ! v'w'
 ,adu_h_ww(:,:,:)       & ! w'w'
 ,adu_h_vv(:,:,:)       & ! v'v'
 ,adu_h_uu(:,:,:)       & ! u'u'
 ,adu_h_uv(:,:,:)       & ! u'v'
 ,adu_h_w3(:,:,:)       & ! w'w'w'
 ,adu_h_dpx(:,:,:)      & ! dp/dx
 ,adu_h_dpy(:,:,:)      & ! dp/dy
 ,adu_h_uth(:,:,:)      & ! u'theta'
 ,adu_h_uthv(:,:,:)     & ! u'thetav'
 ,adu_h_vth(:,:,:)      & ! v'theta'
 ,adu_h_vthv(:,:,:)     & ! v'thetav'
 ,adu_h_uq(:,:,:)       & ! u'q'
 ,adu_h_vq(:,:,:)       & ! v'q'
 ,adu_h_wp(:,:,:)         ! w'p'/density
!------------------------------------------------------------------------------
! ACW - All cloudy upwards air (where upward in not relative)
!------------------------------------------------------------------------------

REAL(wp), ALLOCATABLE ::    &
  acw_h_w(:,:,:)        &
 ,acw_h_th(:,:,:)       &
 ,acw_h_thv(:,:,:)      &
 ,acw_h_rho(:,:,:)      &
 ,acw_h_u(:,:,:)        &
 ,acw_h_v(:,:,:)        &
 ,acw_h_rh(:,:,:)       &
 ,acw_h_dt1(:,:,:)      &
 ,acw_h_dt2(:,:,:)      &
 ,acw_h_dt4(:,:,:)      &
 ,acw_h_dt9(:,:,:)      &
 ,acw_h_dt12(:,:,:)     &
 ,acw_h_dq4(:,:,:)      &
 ,acw_h_dq9(:,:,:)      &
 ,acw_h_dq12(:,:,:)     &
 ,acw_h_dqcl4(:,:,:)    &
 ,acw_h_dqcl9(:,:,:)    &
 ,acw_h_dqcl12(:,:,:)   &
 ,acw_h_dqcf4(:,:,:)    &
 ,acw_h_dqcf3(:,:,:)    &
 ,acw_h_dqcf12(:,:,:)   &
 ,acw_h_dt30(:,:,:)     &
 ,acw_h_dq30(:,:,:)     &
 ,acw_h_dqcl30(:,:,:)   &
 ,acw_h_dqcf30(:,:,:)   &
 ,acw_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,acw_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,acw_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,acw_h_q(:,:,:)        &
 ,acw_h_qcl(:,:,:)      &
 ,acw_h_qcf(:,:,:)      &
 ,acw_h_qrain(:,:,:)    &
 ,acw_h_qgraup(:,:,:)   &
 ,acw_h_a(:,:,:)        & ! fraction of columns in divison
 ,acw_h_thw(:,:,:)      & ! theta'w'
 ,acw_h_thvw(:,:,:)     & ! thetav'w'
 ,acw_h_qw(:,:,:)       & ! q'w'
 ,acw_h_qclw(:,:,:)     & ! qcl'w'
 ,acw_h_qcfw(:,:,:)     & ! qcf'w'
 ,acw_h_qrainw(:,:,:)   & ! qrain'w'
 ,acw_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,acw_h_uw(:,:,:)       & ! u'w'
 ,acw_h_vw(:,:,:)       & ! v'w'
 ,acw_h_ww(:,:,:)       & ! w'w'
 ,acw_h_vv(:,:,:)       & ! v'v'
 ,acw_h_uu(:,:,:)       & ! u'u'
 ,acw_h_uv(:,:,:)       & ! u'v'
 ,acw_h_w3(:,:,:)       & ! w'w'w'
 ,acw_h_dpx(:,:,:)      & ! dp/dx
 ,acw_h_dpy(:,:,:)      & ! dp/dy
 ,acw_h_uth(:,:,:)      & ! u'theta'
 ,acw_h_uthv(:,:,:)     & ! u'thetav'
 ,acw_h_vth(:,:,:)      & ! v'theta'
 ,acw_h_vthv(:,:,:)     & ! v'thetav'
 ,acw_h_uq(:,:,:)       & ! u'q'
 ,acw_h_vq(:,:,:)       & ! v'q'
 ,acw_h_wp(:,:,:)         ! w'p'/density
!------------------------------------------------------------------------------
! BCW - All buoyant cloudy upwards air (where upward in not relative)
!------------------------------------------------------------------------------

REAL(wp), ALLOCATABLE ::    &
  bcw_h_w(:,:,:)        &
 ,bcw_h_th(:,:,:)       &
 ,bcw_h_thv(:,:,:)      &
 ,bcw_h_rho(:,:,:)      &
 ,bcw_h_u(:,:,:)        &
 ,bcw_h_v(:,:,:)        &
 ,bcw_h_rh(:,:,:)       &
 ,bcw_h_dt1(:,:,:)      &
 ,bcw_h_dt2(:,:,:)      &
 ,bcw_h_dt4(:,:,:)      &
 ,bcw_h_dt9(:,:,:)      &
 ,bcw_h_dt12(:,:,:)     &
 ,bcw_h_dq4(:,:,:)      &
 ,bcw_h_dq9(:,:,:)      &
 ,bcw_h_dq12(:,:,:)     &
 ,bcw_h_dqcl4(:,:,:)    &
 ,bcw_h_dqcl9(:,:,:)    &
 ,bcw_h_dqcl12(:,:,:)   &
 ,bcw_h_dqcf4(:,:,:)    &
 ,bcw_h_dqcf3(:,:,:)    &
 ,bcw_h_dqcf12(:,:,:)   &
 ,bcw_h_dt30(:,:,:)     &
 ,bcw_h_dq30(:,:,:)     &
 ,bcw_h_dqcl30(:,:,:)   &
 ,bcw_h_dqcf30(:,:,:)   &
 ,bcw_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,bcw_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,bcw_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,bcw_h_q(:,:,:)        &
 ,bcw_h_qcl(:,:,:)      &
 ,bcw_h_qcf(:,:,:)      &
 ,bcw_h_qrain(:,:,:)    &
 ,bcw_h_qgraup(:,:,:)   &
 ,bcw_h_a(:,:,:)        & ! fraction of columns in divison
 ,bcw_h_thw(:,:,:)      & ! theta'w'
 ,bcw_h_thvw(:,:,:)     & ! thetav'w'
 ,bcw_h_qw(:,:,:)       & ! q'w'
 ,bcw_h_qclw(:,:,:)     & ! qcl'w'
 ,bcw_h_qcfw(:,:,:)     & ! qcf'w'
 ,bcw_h_qrainw(:,:,:)   & ! qrain'w'
 ,bcw_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,bcw_h_uw(:,:,:)       & ! u'w'
 ,bcw_h_vw(:,:,:)       & ! v'w'
 ,bcw_h_ww(:,:,:)       & ! w'w'
 ,bcw_h_vv(:,:,:)       & ! v'v'
 ,bcw_h_uu(:,:,:)       & ! u'u'
 ,bcw_h_uv(:,:,:)       & ! u'v'
 ,bcw_h_w3(:,:,:)       & ! w'w'w'
 ,bcw_h_dpx(:,:,:)      & ! dp/dx
 ,bcw_h_dpy(:,:,:)      & ! dp/dy
 ,bcw_h_uth(:,:,:)      & ! u'theta'
 ,bcw_h_uthv(:,:,:)     & ! u'thetav'
 ,bcw_h_vth(:,:,:)      & ! v'theta'
 ,bcw_h_vthv(:,:,:)     & ! v'thetav'
 ,bcw_h_uq(:,:,:)       & ! u'q'
 ,bcw_h_vq(:,:,:)       & ! v'q'
 ,bcw_h_wp(:,:,:)         ! w'p'/density


!------------------------------------------------------------------------------
! UCU - unstable cloudy updraughts (-dln(thesat)/dz > 0.0, w > 0.0)
!------------------------------------------------------------------------------

REAL(wp), ALLOCATABLE ::    &
  ucu_h_w(:,:,:)        &
 ,ucu_h_th(:,:,:)       &
 ,ucu_h_thv(:,:,:)      &
 ,ucu_h_rho(:,:,:)      &
 ,ucu_h_u(:,:,:)        &
 ,ucu_h_v(:,:,:)        &
 ,ucu_h_rh(:,:,:)       &
 ,ucu_h_dt1(:,:,:)      &
 ,ucu_h_dt2(:,:,:)      &
 ,ucu_h_dt4(:,:,:)      &
 ,ucu_h_dt9(:,:,:)      &
 ,ucu_h_dt12(:,:,:)     &
 ,ucu_h_dq4(:,:,:)      &
 ,ucu_h_dq9(:,:,:)      &
 ,ucu_h_dq12(:,:,:)     &
 ,ucu_h_dqcl4(:,:,:)    &
 ,ucu_h_dqcl9(:,:,:)    &
 ,ucu_h_dqcl12(:,:,:)   &
 ,ucu_h_dqcf4(:,:,:)    &
 ,ucu_h_dqcf3(:,:,:)    &
 ,ucu_h_dqcf12(:,:,:)   &
 ,ucu_h_dt30(:,:,:)     &
 ,ucu_h_dq30(:,:,:)     &
 ,ucu_h_dqcl30(:,:,:)   &
 ,ucu_h_dqcf30(:,:,:)   &
 ,ucu_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,ucu_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,ucu_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,ucu_h_q(:,:,:)        &
 ,ucu_h_qcl(:,:,:)      &
 ,ucu_h_qcf(:,:,:)      &
 ,ucu_h_qrain(:,:,:)    &
 ,ucu_h_qgraup(:,:,:)   &
 ,ucu_h_a(:,:,:)        & ! fraction of columns in divison
 ,ucu_h_thw(:,:,:)      & ! theta'w'
 ,ucu_h_thvw(:,:,:)     & ! thetav'w'
 ,ucu_h_qw(:,:,:)       & ! q'w'
 ,ucu_h_qclw(:,:,:)     & ! qcl'w'
 ,ucu_h_qcfw(:,:,:)     & ! qcf'w'
 ,ucu_h_qrainw(:,:,:)   & ! qrain'w'
 ,ucu_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,ucu_h_uw(:,:,:)       & ! u'w'
 ,ucu_h_vw(:,:,:)       & ! v'w'
 ,ucu_h_ww(:,:,:)       & ! w'w'
 ,ucu_h_vv(:,:,:)       & ! v'v'
 ,ucu_h_uu(:,:,:)       & ! u'u'
 ,ucu_h_uv(:,:,:)       & ! u'v'
 ,ucu_h_w3(:,:,:)       & ! w'w'w'
 ,ucu_h_dpx(:,:,:)      & ! dp/dx
 ,ucu_h_dpy(:,:,:)      & ! dp/dy
 ,ucu_h_uth(:,:,:)      & ! u'theta'
 ,ucu_h_uthv(:,:,:)     & ! u'thetav'
 ,ucu_h_vth(:,:,:)      & ! v'theta'
 ,ucu_h_vthv(:,:,:)     & ! v'thetav'
 ,ucu_h_uq(:,:,:)       & ! u'q'
 ,ucu_h_vq(:,:,:)       & ! v'q'
 ,ucu_h_wp(:,:,:)         ! w'p'/density

!------------------------------------------------------------------------------
! PPW - precipitation downdraughts (w<0.0)
!------------------------------------------------------------------------------

REAL(wp), ALLOCATABLE ::    &
  ppw_h_w(:,:,:)        &
 ,ppw_h_th(:,:,:)       &
 ,ppw_h_thv(:,:,:)      &
 ,ppw_h_rho(:,:,:)      &
 ,ppw_h_u(:,:,:)        &
 ,ppw_h_v(:,:,:)        &
 ,ppw_h_rh(:,:,:)       &
 ,ppw_h_dt1(:,:,:)      &
 ,ppw_h_dt2(:,:,:)      &
 ,ppw_h_dt4(:,:,:)      &
 ,ppw_h_dt9(:,:,:)      &
 ,ppw_h_dt12(:,:,:)     &
 ,ppw_h_dq4(:,:,:)      &
 ,ppw_h_dq9(:,:,:)      &
 ,ppw_h_dq12(:,:,:)     &
 ,ppw_h_dqcl4(:,:,:)    &
 ,ppw_h_dqcl9(:,:,:)    &
 ,ppw_h_dqcl12(:,:,:)   &
 ,ppw_h_dqcf4(:,:,:)    &
 ,ppw_h_dqcf3(:,:,:)    &
 ,ppw_h_dqcf12(:,:,:)   &
 ,ppw_h_dt30(:,:,:)     &
 ,ppw_h_dq30(:,:,:)     &
 ,ppw_h_dqcl30(:,:,:)   &
 ,ppw_h_dqcf30(:,:,:)   &
 ,ppw_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,ppw_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,ppw_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,ppw_h_q(:,:,:)        &
 ,ppw_h_qcl(:,:,:)      &
 ,ppw_h_qcf(:,:,:)      &
 ,ppw_h_qrain(:,:,:)    &
 ,ppw_h_qgraup(:,:,:)   &
 ,ppw_h_a(:,:,:)        & ! fraction of columns in divison
 ,ppw_h_thw(:,:,:)      & ! theta'w'
 ,ppw_h_thvw(:,:,:)     & ! thetav'w'
 ,ppw_h_qw(:,:,:)       & ! q'w'
 ,ppw_h_qclw(:,:,:)     & ! qcl'w'
 ,ppw_h_qcfw(:,:,:)     & ! qcf'w'
 ,ppw_h_qrainw(:,:,:)   & ! qrain'w'
 ,ppw_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,ppw_h_uw(:,:,:)       & ! u'w'
 ,ppw_h_vw(:,:,:)       & ! v'w'
 ,ppw_h_ww(:,:,:)       & ! w'w'
 ,ppw_h_vv(:,:,:)       & ! v'v'
 ,ppw_h_uu(:,:,:)       & ! u'u'
 ,ppw_h_uv(:,:,:)       & ! u'v'
 ,ppw_h_w3(:,:,:)       & ! w'w'w'
 ,ppw_h_dpx(:,:,:)      & ! dp/dx
 ,ppw_h_dpy(:,:,:)      & ! dp/dy
 ,ppw_h_uth(:,:,:)      & ! u'theta'
 ,ppw_h_uthv(:,:,:)     & ! u'thetav'
 ,ppw_h_vth(:,:,:)      & ! v'theta'
 ,ppw_h_vthv(:,:,:)     & ! v'thetav'
 ,ppw_h_uq(:,:,:)       & ! u'q'
 ,ppw_h_vq(:,:,:)       & ! v'q'
 ,ppw_h_wp(:,:,:)         ! w'p'/density

!------------------------------------------------------------------------------
! NBW - negatively buoyant precipitation downdraughts (w<0.0)
!------------------------------------------------------------------------------

REAL(wp), ALLOCATABLE ::    &
  nbw_h_w(:,:,:)        &
 ,nbw_h_th(:,:,:)       &
 ,nbw_h_thv(:,:,:)      &
 ,nbw_h_rho(:,:,:)      &
 ,nbw_h_u(:,:,:)        &
 ,nbw_h_v(:,:,:)        &
 ,nbw_h_rh(:,:,:)       &
 ,nbw_h_dt1(:,:,:)      &
 ,nbw_h_dt2(:,:,:)      &
 ,nbw_h_dt4(:,:,:)      &
 ,nbw_h_dt9(:,:,:)      &
 ,nbw_h_dt12(:,:,:)     &
 ,nbw_h_dq4(:,:,:)      &
 ,nbw_h_dq9(:,:,:)      &
 ,nbw_h_dq12(:,:,:)     &
 ,nbw_h_dqcl4(:,:,:)    &
 ,nbw_h_dqcl9(:,:,:)    &
 ,nbw_h_dqcl12(:,:,:)   &
 ,nbw_h_dqcf4(:,:,:)    &
 ,nbw_h_dqcf3(:,:,:)    &
 ,nbw_h_dqcf12(:,:,:)   &
 ,nbw_h_dt30(:,:,:)     &
 ,nbw_h_dq30(:,:,:)     &
 ,nbw_h_dqcl30(:,:,:)   &
 ,nbw_h_dqcf30(:,:,:)   &
 ,nbw_h_dqrain30(:,:,:) & ! dqrain from total kg/kg/timestep
 ,nbw_h_dqgr30(:,:,:)   & ! dqgr from total kg/kg/timestep
 ,nbw_h_drho(:,:,:)     & ! drho from total kg/m3/timestep
 ,nbw_h_q(:,:,:)        &
 ,nbw_h_qcl(:,:,:)      &
 ,nbw_h_qcf(:,:,:)      &
 ,nbw_h_qrain(:,:,:)    &
 ,nbw_h_qgraup(:,:,:)   &
 ,nbw_h_a(:,:,:)        & ! fraction of columns in divison
 ,nbw_h_thw(:,:,:)      & ! theta'w'
 ,nbw_h_thvw(:,:,:)     & ! thetav'w'
 ,nbw_h_qw(:,:,:)       & ! q'w'
 ,nbw_h_qclw(:,:,:)     & ! qcl'w'
 ,nbw_h_qcfw(:,:,:)     & ! qcf'w'
 ,nbw_h_qrainw(:,:,:)   & ! qrain'w'
 ,nbw_h_qgraupw(:,:,:)  & ! qgraup'w'
 ,nbw_h_uw(:,:,:)       & ! u'w'
 ,nbw_h_vw(:,:,:)       & ! v'w'
 ,nbw_h_ww(:,:,:)       & ! w'w'
 ,nbw_h_vv(:,:,:)       & ! v'v'
 ,nbw_h_uu(:,:,:)       & ! u'u'
 ,nbw_h_uv(:,:,:)       & ! u'v'
 ,nbw_h_w3(:,:,:)       & ! w'w'w'
 ,nbw_h_dpx(:,:,:)      & ! dp/dx
 ,nbw_h_dpy(:,:,:)      & ! dp/dy
 ,nbw_h_uth(:,:,:)      & ! u'theta'
 ,nbw_h_uthv(:,:,:)     & ! u'thetav'
 ,nbw_h_vth(:,:,:)      & ! v'theta'
 ,nbw_h_vthv(:,:,:)     & ! v'thetav'
 ,nbw_h_uq(:,:,:)       & ! u'q'
 ,nbw_h_vq(:,:,:)       & ! v'q'
 ,nbw_h_wp(:,:,:)         ! w'p'/density
!------------------------------------------------------------------------------
! Full size arrays for Field minus mean of field
!------------------------------------------------------------------------------
REAL(wp), ALLOCATABLE ::       &
  th_h_prime(:,:,:)        & !
 ,thv_h_prime(:,:,:)       & !
 ,q_h_prime(:,:,:)         & !
 ,qcl_h_prime(:,:,:)       & !
 ,qcf_h_prime(:,:,:)       & !
 ,qrain_h_prime(:,:,:)     & !
 ,qgraup_h_prime(:,:,:)    & !
 ,w_h_prime(:,:,:)         & !
 ,u_h_prime(:,:,:)         & !
 ,v_h_prime(:,:,:)         & !
 ,p_h_prime(:,:,:)


!------------------------------------------------------------------------------
! plume counting
!------------------------------------------------------------------------------
REAL(wp), ALLOCATABLE ::       &
  n_plume(:,:,:)               & ! number of plumes
 ,plume_size(:,:,:)            & ! mean plume size
 ,plume_diam_pdf(:,:,:,:)        ! pdf of plume diameters


!------------------------------------------------------------------------------
! Mask of buoyant points
!------------------------------------------------------------------------------
LOGICAL, ALLOCATABLE ::    &
  bcu_mask(:,:,:)            ! T if buoyant otherwise F

REAL(wp), ALLOCATABLE ::   &
  bcu_mask_w(:,:,:)          ! w value

INTEGER, ALLOCATABLE ::    &
  nn_mask(:,:,:)             !  number of buoyant points in each level

END MODULE crmstyle_sample_arrays_mod


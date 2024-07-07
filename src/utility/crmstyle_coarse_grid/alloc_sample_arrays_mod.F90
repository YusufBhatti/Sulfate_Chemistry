! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
! ALLOCATE arrays for crm sampling
!
MODULE alloc_sample_arrays_mod

! ------------------------------------------------------------------------------
! Description:
!   ALLOCATE and initialise arrays required for CRM sampling.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CRMstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE

PUBLIC :: alloc_sample_arrays

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOC_SAMPLE_ARRAYS_MOD'

CONTAINS

SUBROUTINE alloc_sample_arrays(local_num_x,local_num_y,mlevs,row_length,rows,&
                               bins_diam)

USE crmstyle_sample_arrays_mod
USE crmstyle_cntl_mod, ONLY: l_sect30, l_qgraup, l_cape, l_pcape,           &
  l_bcu, l_wg1, l_acc, l_acu, l_ppd, l_nbd, l_nid, l_adu, l_acw, l_bcw,     &
  l_ucu, l_ppw, l_nbw

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::  &
  local_num_x           & ! Coarse grid columns
 ,local_num_y           & ! Coarse grid rows
 ,mlevs                 & ! Number of levels
 ,row_length            & ! Number of columns
 ,rows                  & ! Number of rows
 ,bins_diam               ! Number of bins for diameter pdf

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOC_SAMPLE_ARRAYS'

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

! single level
ALLOCATE (all_zh(local_num_x,local_num_y) )
ALLOCATE (all_lh(local_num_x,local_num_y) )
ALLOCATE (all_sh(local_num_x,local_num_y) )
ALLOCATE (all_pstar(local_num_x,local_num_y) )
ALLOCATE (all_tstar(local_num_x,local_num_y) )
ALLOCATE (all_rain(local_num_x,local_num_y) )
ALLOCATE (all_snow(local_num_x,local_num_y) )
ALLOCATE (all_precip(local_num_x,local_num_y) )
ALLOCATE (all_orog(local_num_x,local_num_y) )
ALLOCATE (all_land(local_num_x,local_num_y) )
! single level
ALLOCATE (all_sd_zh(local_num_x,local_num_y) )
ALLOCATE (all_sd_lh(local_num_x,local_num_y) )
ALLOCATE (all_sd_sh(local_num_x,local_num_y) )
ALLOCATE (all_sd_pstar(local_num_x,local_num_y) )
ALLOCATE (all_sd_tstar(local_num_x,local_num_y) )
ALLOCATE (all_sd_rain(local_num_x,local_num_y) )
ALLOCATE (all_sd_snow(local_num_x,local_num_y) )
ALLOCATE (all_sd_precip(local_num_x,local_num_y) )
ALLOCATE (all_sd_orog(local_num_x,local_num_y) )
ALLOCATE (all_sd_land(local_num_x,local_num_y) )


! LCL level needed by CAPE & PCAPE
IF (l_cape .OR. l_pcape) THEN
  ALLOCATE (zlcl(local_num_x,local_num_y) )
ELSE
  ALLOCATE (zlcl(1,1) )
END IF

! CAPE & CIN
IF (l_cape) THEN
  ALLOCATE (cape(local_num_x,local_num_y) )
  ALLOCATE (cin(local_num_x,local_num_y) )
  ALLOCATE (zfree(local_num_x,local_num_y) )
  ALLOCATE (zneutral(local_num_x,local_num_y) )
ELSE  ! unit arrays so that can be passed in arguements
  ALLOCATE (cape(1,1) )
  ALLOCATE (cin(1,1) )
  ALLOCATE (zfree(1,1) )
  ALLOCATE (zneutral(1,1) )
END IF
! PCAPE 
IF (l_pcape .AND. l_bcu) THEN
  ALLOCATE (bcu_pcape(local_num_x,local_num_y) )
  ALLOCATE (bcu_dpcapedt(local_num_x,local_num_y) )
  ALLOCATE (bcu_dpcapedt_bl(local_num_x,local_num_y) )
  ALLOCATE (bcu_dilcape(local_num_x,local_num_y) )
  ALLOCATE (bcu_dcapedt(local_num_x,local_num_y) )
ELSE  ! unit arrays so that can be passed in arguements
  ALLOCATE (bcu_pcape(1,1) )
  ALLOCATE (bcu_dpcapedt(1,1) )
  ALLOCATE (bcu_dpcapedt_bl(1,1) )
  ALLOCATE (bcu_dilcape(1,1) )
  ALLOCATE (bcu_dcapedt(1,1) )
END IF
IF (l_pcape .AND. l_bcu .AND. l_ppd) THEN
  ALLOCATE (ppd_dpcapedt_bl(local_num_x,local_num_y) )
ELSE  ! unit arrays so that can be passed in arguements
  ALLOCATE (ppd_dpcapedt_bl(1,1) )
END IF

IF (l_pcape .AND. l_bcw) THEN
  ALLOCATE (bcw_pcape(local_num_x,local_num_y) )
  ALLOCATE (bcw_dpcapedt(local_num_x,local_num_y) )
  ALLOCATE (bcw_dpcapedt_bl(local_num_x,local_num_y) )
  ALLOCATE (bcw_dilcape(local_num_x,local_num_y) )
  ALLOCATE (bcw_dcapedt(local_num_x,local_num_y) )
ELSE  ! unit arrays so that can be passed in arguements
  ALLOCATE (bcw_pcape(1,1) )
  ALLOCATE (bcw_dpcapedt(1,1) )
  ALLOCATE (bcw_dpcapedt_bl(1,1) )
  ALLOCATE (bcw_dilcape(1,1) )
  ALLOCATE (bcw_dcapedt(1,1) )
END IF

IF (l_cape  .OR. l_pcape) THEN
  ALLOCATE (all_exner(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_tv(local_num_x,local_num_y,mlevs) )
END IF


ALLOCATE (all_th(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_t(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_thv(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_q(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_qcl(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_qcf(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_qrain(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_qgraup(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_u(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_v(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_w(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_ptheta(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_a(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_rh(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_rho(local_num_x,local_num_y,mlevs) )

! Mean tendencies

ALLOCATE (all_dt1(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dt2(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dt4(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dt9(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dt12(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dq4(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dq9(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dq12(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dqcl4(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dqcl9(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dqcl12(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dqcf4(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dqcf3(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dqcf12(local_num_x,local_num_y,mlevs) )

ALLOCATE (all_sd_th(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_thv(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_t(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_q(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_qcl(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_qcf(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_qrain(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_qgraup(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_u(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_v(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_w(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_ptheta(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_rh(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_rho(local_num_x,local_num_y,mlevs) )

! Section 30 
IF (l_sect30) THEN
  ALLOCATE (all_dt30(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_dq30(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_dqcl30(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_dqcf30(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_dqrain30(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_drho(local_num_x,local_num_y,mlevs) )
ELSE  ! unit arrays so that can be passed in arguements
  ALLOCATE (all_dt30(1,1,1) )
  ALLOCATE (all_dq30(1,1,1) )
  ALLOCATE (all_dqcl30(1,1,1) )
  ALLOCATE (all_dqcf30(1,1,1) )
  ALLOCATE (all_dqrain30(1,1,1) )
  ALLOCATE (all_drho(1,1,1) )
END IF

! s.d tendencies

ALLOCATE (all_sd_dt1(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dt2(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dt4(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dt9(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dt12(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dq4(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dq9(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dq12(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dqcl4(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dqcl9(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dqcl12(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dqcf4(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dqcf3(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_sd_dqcf12(local_num_x,local_num_y,mlevs) )


! Section 30 
IF (l_sect30) THEN
  ALLOCATE (all_sd_dt30(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_sd_dq30(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_sd_dqcl30(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_sd_dqcf30(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_sd_dqrain30(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_sd_drho(local_num_x,local_num_y,mlevs) )
ELSE  ! unit arrays so that can be passed in arguements
  ALLOCATE (all_sd_dt30(1,1,1) )
  ALLOCATE (all_sd_dq30(1,1,1) )
  ALLOCATE (all_sd_dqcl30(1,1,1) )
  ALLOCATE (all_sd_dqcf30(1,1,1) )
  ALLOCATE (all_sd_dqrain30(1,1,1) )
  ALLOCATE (all_sd_drho(1,1,1) )
END IF

IF (l_sect30 .AND. l_qgraup) THEN
  ALLOCATE (all_dqgr30(local_num_x,local_num_y,mlevs) )
  ALLOCATE (all_sd_dqgr30(local_num_x,local_num_y,mlevs) )
ELSE  ! unit arrays so that can be passed in arguements
  ALLOCATE (all_dqgr30(1,1,1) )
  ALLOCATE (all_sd_dqgr30(1,1,1) )
END IF

! Mean products
ALLOCATE (all_thw(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_thvw(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_qw(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_qclw(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_qcfw(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_qrainw(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_qgraupw(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_uw(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_vw(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_ww(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_w3(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_vv(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_uu(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_uv(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dpx(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_dpy(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_uth(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_vth(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_uthv(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_vthv(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_uq(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_vq(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_wp(local_num_x,local_num_y,mlevs) )
ALLOCATE (all_wp_hydro(local_num_x,local_num_y,mlevs) )

ALLOCATE (p_theta_hydro(local_num_x,local_num_y,mlevs) )
! --------------------------------------------------------------------------
! Partitions
! --------------------------------------------------------------------------
! All cloudy points

IF (l_acc) THEN
  ALLOCATE (acc_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acc_h_wp(local_num_x,local_num_y,mlevs) )

  IF (l_sect30) THEN
    ALLOCATE (acc_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acc_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acc_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acc_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acc_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acc_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (acc_h_dt30(1,1,1) )
    ALLOCATE (acc_h_dq30(1,1,1) )
    ALLOCATE (acc_h_dqcl30(1,1,1) )
    ALLOCATE (acc_h_dqcf30(1,1,1) )
    ALLOCATE (acc_h_dqrain30(1,1,1) )
    ALLOCATE (acc_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (acc_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (acc_h_dqgr30(1,1,1) )
  END IF
END IF

! All cloudy updraughts
IF (l_acu) THEN
  ALLOCATE (acu_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acu_h_wp(local_num_x,local_num_y,mlevs) )

  IF (l_sect30) THEN
    ALLOCATE (acu_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acu_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acu_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acu_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acu_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acu_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (acu_h_dt30(1,1,1) )
    ALLOCATE (acu_h_dq30(1,1,1) )
    ALLOCATE (acu_h_dqcl30(1,1,1) )
    ALLOCATE (acu_h_dqcf30(1,1,1) )
    ALLOCATE (acu_h_dqrain30(1,1,1) )
    ALLOCATE (acu_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (acu_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (acu_h_dqgr30(1,1,1) )
  END IF
END IF

! Buoyant cloudy updraughts  - relative to mean w

IF (l_bcu) THEN
  ALLOCATE (bcu_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcu_h_wp(local_num_x,local_num_y,mlevs) )

  IF (l_sect30) THEN
    ALLOCATE (bcu_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (bcu_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (bcu_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (bcu_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (bcu_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (bcu_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (bcu_h_dt30(1,1,1) )
    ALLOCATE (bcu_h_dq30(1,1,1) )
    ALLOCATE (bcu_h_dqcl30(1,1,1) )
    ALLOCATE (bcu_h_dqcf30(1,1,1) )
    ALLOCATE (bcu_h_dqrain30(1,1,1) )
    ALLOCATE (bcu_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (bcu_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (bcu_h_dqgr30(1,1,1) )
  END IF
END IF

! Strong buoyant updraughts (>1)  - relative to mean w
IF (l_wg1) THEN
  ALLOCATE (wg1_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (wg1_h_wp(local_num_x,local_num_y,mlevs) )

  IF (l_sect30) THEN
    ALLOCATE (wg1_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (wg1_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (wg1_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (wg1_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (wg1_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (wg1_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (wg1_h_dt30(1,1,1) )
    ALLOCATE (wg1_h_dq30(1,1,1) )
    ALLOCATE (wg1_h_dqcl30(1,1,1) )
    ALLOCATE (wg1_h_dqcf30(1,1,1) )
    ALLOCATE (wg1_h_dqrain30(1,1,1) )
    ALLOCATE (wg1_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (wg1_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (wg1_h_dqgr30(1,1,1) )
  END IF
END IF
! precipitation downdraughts
IF (l_ppd) THEN
  ALLOCATE (  ppd_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppd_h_wp(local_num_x,local_num_y,mlevs) )

  IF (l_sect30) THEN
    ALLOCATE (ppd_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ppd_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ppd_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ppd_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ppd_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ppd_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (ppd_h_dt30(1,1,1) )
    ALLOCATE (ppd_h_dq30(1,1,1) )
    ALLOCATE (ppd_h_dqcl30(1,1,1) )
    ALLOCATE (ppd_h_dqcf30(1,1,1) )
    ALLOCATE (ppd_h_dqrain30(1,1,1) )
    ALLOCATE (ppd_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (ppd_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (ppd_h_dqgr30(1,1,1) )
  END IF
END IF


! Negatively buoyant precipitation downdraughts
IF (l_nbd) THEN
  ALLOCATE (  nbd_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbd_h_wp(local_num_x,local_num_y,mlevs) )
  IF (l_sect30) THEN
    ALLOCATE (nbd_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nbd_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nbd_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nbd_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nbd_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nbd_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (nbd_h_dt30(1,1,1) )
    ALLOCATE (nbd_h_dq30(1,1,1) )
    ALLOCATE (nbd_h_dqcl30(1,1,1) )
    ALLOCATE (nbd_h_dqcf30(1,1,1) )
    ALLOCATE (nbd_h_dqrain30(1,1,1) )
    ALLOCATE (nbd_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (nbd_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (nbd_h_dqgr30(1,1,1) )
  END IF
END IF

!------------------------------------------------------------------------------
! Negatively buoyant precipitation with ice downdraughts
!------------------------------------------------------------------------------
IF (l_nid) THEN
  ALLOCATE (nid_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nid_h_wp(local_num_x,local_num_y,mlevs) )
  IF (l_sect30) THEN
    ALLOCATE (nid_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nid_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nid_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nid_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nid_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nid_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (nid_h_dt30(1,1,1) )
    ALLOCATE (nid_h_dq30(1,1,1) )
    ALLOCATE (nid_h_dqcl30(1,1,1) )
    ALLOCATE (nid_h_dqcf30(1,1,1) )
    ALLOCATE (nid_h_dqrain30(1,1,1) )
    ALLOCATE (nid_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (nid_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (nid_h_dqgr30(1,1,1) )
  END IF
END IF

!------------------------------------------------------------------------------
! ADU - Dry upward air  (where upward in not relative)
!------------------------------------------------------------------------------
IF (l_adu) THEN
  ALLOCATE (adu_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (adu_h_wp(local_num_x,local_num_y,mlevs) )
  IF (l_sect30) THEN
    ALLOCATE (adu_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (adu_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (adu_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (adu_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (adu_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (adu_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (adu_h_dt30(1,1,1) )
    ALLOCATE (adu_h_dq30(1,1,1) )
    ALLOCATE (adu_h_dqcl30(1,1,1) )
    ALLOCATE (adu_h_dqcf30(1,1,1) )
    ALLOCATE (adu_h_dqrain30(1,1,1) )
    ALLOCATE (adu_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (adu_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (adu_h_dqgr30(1,1,1) )
  END IF
END IF
!------------------------------------------------------------------------------
! ACW - All cloudy upwards air (where upward in not relative)
!------------------------------------------------------------------------------
IF (l_acw) THEN
  ALLOCATE (acw_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (acw_h_wp(local_num_x,local_num_y,mlevs) )
  IF (l_sect30) THEN
    ALLOCATE (acw_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acw_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acw_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acw_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acw_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (acw_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (acw_h_dt30(1,1,1) )
    ALLOCATE (acw_h_dq30(1,1,1) )
    ALLOCATE (acw_h_dqcl30(1,1,1) )
    ALLOCATE (acw_h_dqcf30(1,1,1) )
    ALLOCATE (acw_h_dqrain30(1,1,1) )
    ALLOCATE (acw_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (acw_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (acw_h_dqgr30(1,1,1) )
  END IF
END IF
!------------------------------------------------------------------------------
!BCW - All buoyant cloudy upwards air (where upward in not relative)
!------------------------------------------------------------------------------
IF (l_bcw) THEN
  ALLOCATE (bcw_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (bcw_h_wp(local_num_x,local_num_y,mlevs) )
  IF (l_sect30) THEN
    ALLOCATE (bcw_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (bcw_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (bcw_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (bcw_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (bcw_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (bcw_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (bcw_h_dt30(1,1,1) )
    ALLOCATE (bcw_h_dq30(1,1,1) )
    ALLOCATE (bcw_h_dqcl30(1,1,1) )
    ALLOCATE (bcw_h_dqcf30(1,1,1) )
    ALLOCATE (bcw_h_dqrain30(1,1,1) )
    ALLOCATE (bcw_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (bcw_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (bcw_h_dqgr30(1,1,1) )
  END IF
END IF

!------------------------------------------------------------------------------
! UCU - unstable cloudy updraughts (-dln(thesat)/dz > 0.0, w > 0.0)
!------------------------------------------------------------------------------
IF (l_ucu) THEN
  ALLOCATE (ucu_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ucu_h_wp(local_num_x,local_num_y,mlevs) )
  IF (l_sect30) THEN
    ALLOCATE (ucu_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ucu_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ucu_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ucu_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ucu_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ucu_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (ucu_h_dt30(1,1,1) )
    ALLOCATE (ucu_h_dq30(1,1,1) )
    ALLOCATE (ucu_h_dqcl30(1,1,1) )
    ALLOCATE (ucu_h_dqcf30(1,1,1) )
    ALLOCATE (ucu_h_dqrain30(1,1,1) )
    ALLOCATE (ucu_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (ucu_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (ucu_h_dqgr30(1,1,1) )
  END IF
END IF

!------------------------------------------------------------------------------
! PPW - precipitating downdraughts (w<0)
!------------------------------------------------------------------------------
IF (l_ppw) THEN
  ALLOCATE (ppw_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (ppw_h_wp(local_num_x,local_num_y,mlevs) )
  IF (l_sect30) THEN
    ALLOCATE (ppw_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ppw_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ppw_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ppw_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ppw_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (ppw_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (ppw_h_dt30(1,1,1) )
    ALLOCATE (ppw_h_dq30(1,1,1) )
    ALLOCATE (ppw_h_dqcl30(1,1,1) )
    ALLOCATE (ppw_h_dqcf30(1,1,1) )
    ALLOCATE (ppw_h_dqrain30(1,1,1) )
    ALLOCATE (ppw_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (ppw_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (ppw_h_dqgr30(1,1,1) )
  END IF
END IF
!------------------------------------------------------------------------------
! NBW - negatively buoyant precipitating downdraughts (w<0)
!------------------------------------------------------------------------------
IF (l_nbw) THEN
  ALLOCATE (nbw_h_w(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_th(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_thv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_rho(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_u(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_v(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_rh(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dt1(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dt2(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dt4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dt9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dt12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dq4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dq9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dq12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dqcl4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dqcl9(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dqcl12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dqcf4(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dqcf3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dqcf12(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_q(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_qcl(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_qcf(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_qrain(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_qgraup(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_a(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_thw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_thvw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_qw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_qclw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_qcfw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_qrainw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_qgraupw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_uw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_vw(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_ww(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_vv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_uu(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_uv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_w3(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dpx(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_dpy(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_uth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_vth(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_uthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_vthv(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_uq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_vq(local_num_x,local_num_y,mlevs) )
  ALLOCATE (nbw_h_wp(local_num_x,local_num_y,mlevs) )
  IF (l_sect30) THEN
    ALLOCATE (nbw_h_dt30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nbw_h_dq30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nbw_h_dqcl30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nbw_h_dqcf30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nbw_h_dqrain30(local_num_x,local_num_y,mlevs) )
    ALLOCATE (nbw_h_drho(local_num_x,local_num_y,mlevs) )
  ELSE
    ALLOCATE (nbw_h_dt30(1,1,1) )
    ALLOCATE (nbw_h_dq30(1,1,1) )
    ALLOCATE (nbw_h_dqcl30(1,1,1) )
    ALLOCATE (nbw_h_dqcf30(1,1,1) )
    ALLOCATE (nbw_h_dqrain30(1,1,1) )
    ALLOCATE (nbw_h_drho(1,1,1) )
  END IF
  IF (l_sect30 .AND. l_qgraup) THEN
    ALLOCATE (nbw_h_dqgr30(local_num_x,local_num_y,mlevs) )
  ELSE  ! unit arrays so that can be passed in arguements
    ALLOCATE (nbw_h_dqgr30(1,1,1) )
  END IF
END IF
!------------------------------------------------------------------------------
! Full size arrays for Field minus mean of field
!------------------------------------------------------------------------------

ALLOCATE (th_h_prime(row_length,rows,mlevs) )
ALLOCATE (thv_h_prime(row_length,rows,mlevs) )
ALLOCATE (q_h_prime(row_length,rows,mlevs) )
ALLOCATE (qcl_h_prime(row_length,rows,mlevs) )
ALLOCATE (qcf_h_prime(row_length,rows,mlevs) )
ALLOCATE (qrain_h_prime(row_length,rows,mlevs) )
ALLOCATE (qgraup_h_prime(row_length,rows,mlevs) )
ALLOCATE (w_h_prime(row_length,rows,mlevs) )
ALLOCATE (u_h_prime(row_length,rows,mlevs) )
ALLOCATE (v_h_prime(row_length,rows,mlevs) )
ALLOCATE (p_h_prime(row_length,rows,mlevs) )

!------------------------------------------------------------------------------
! Plume arrays
!------------------------------------------------------------------------------

ALLOCATE (n_plume(local_num_x,local_num_y,mlevs) )
ALLOCATE (plume_size(local_num_x,local_num_y,mlevs) )
ALLOCATE (plume_diam_pdf(bins_diam,local_num_x,local_num_y,mlevs) )

!------------------------------------------------------------------------------
! Mask of buoyant points
!------------------------------------------------------------------------------
ALLOCATE (bcu_mask(row_length,rows,mlevs) )
ALLOCATE (bcu_mask_w(row_length,rows,mlevs) )

ALLOCATE (  nn_mask(local_num_x,local_num_y,mlevs) )

!------------------------------------------------------------------------------
! Column classification arrays based on surface precipitation
!------------------------------------------------------------------------------

ALLOCATE (fract_conv(local_num_x,local_num_y) )
ALLOCATE (fract_strat(local_num_x,local_num_y) )
ALLOCATE (prec_conv(local_num_x,local_num_y) )
ALLOCATE (prec_strat(local_num_x,local_num_y) )

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE alloc_sample_arrays
END MODULE alloc_sample_arrays_mod

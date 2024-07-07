! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  Zero arrays for partition

MODULE crmstyle_zero_arrays_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_ZERO_ARRAYS_MOD'

CONTAINS

! ------------------------------------------------------------------------------
! Description:
!  zeros arrays for partition
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------


SUBROUTINE crmstyle_zero_arrays(crm_factor,                               &  
  crm_h_w, crm_h_u, crm_h_v, crm_h_th, crm_h_thv, crm_h_rho,              &
  crm_h_rh, crm_h_a, crm_h_dpx, crm_h_dpy,                                &
  crm_h_q,   crm_h_qcl, crm_h_qcf, crm_h_qrain, crm_h_qgraup,             &
  crm_h_dt1, crm_h_dt2, crm_h_dt4, crm_h_dt9, crm_h_dt12, crm_h_dt30,     &
  crm_h_dq4, crm_h_dq9, crm_h_dq12, crm_h_dq30,                           &
  crm_h_dqcl4, crm_h_dqcl9, crm_h_dqcl12, crm_h_dqcl30,                   &
  crm_h_dqcf4, crm_h_dqcf3, crm_h_dqcf12, crm_h_dqcf30,                   &
  crm_h_dqrain30, crm_h_dqgr30, crm_h_drho, crm_h_thw, crm_h_thvw,        &
  crm_h_qw, crm_h_qclw, crm_h_qcfw, crm_h_qrainw, crm_h_qgraupw,          &
  crm_h_uth, crm_h_vth, crm_h_uthv, crm_h_vthv, crm_h_uq, crm_h_vq,       &
  crm_h_wp,                                                               &
  crm_h_uw, crm_h_vw, crm_h_ww, crm_h_w3, crm_h_vv, crm_h_uu, crm_h_uv )


USE crmstyle_cntl_mod, ONLY:                                      &
   mlevs, new_res, l_qgraup, l_sect30

USE crmstyle_grid_info_mod, ONLY:                                 &
   local_new_x, local_new_y

USE word_sizes_mod, ONLY: iwp,wp   ! Allows use of 4 byte words to reduce
                                   ! memory

USE umPrintMgr                              ! for writing output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE


REAL, INTENT(INOUT) ::                          &
  crm_factor(local_new_x,local_new_y,mlevs)

REAL(wp), INTENT(INOUT) ::                      &
  crm_h_w(local_new_x,local_new_y,mlevs)        &
 ,crm_h_u(local_new_x,local_new_y,mlevs)        &
 ,crm_h_v(local_new_x,local_new_y,mlevs)        &
 ,crm_h_th(local_new_x,local_new_y,mlevs)       &
 ,crm_h_thv(local_new_x,local_new_y,mlevs)      &
 ,crm_h_rho(local_new_x,local_new_y,mlevs)      &
 ,crm_h_rh(local_new_x,local_new_y,mlevs)       &
 ,crm_h_a(local_new_x,local_new_y,mlevs)        &
 ,crm_h_dpx(local_new_x,local_new_y,mlevs)      &
 ,crm_h_dpy(local_new_x,local_new_y,mlevs)      &
 ,crm_h_q(local_new_x,local_new_y,mlevs)        &
 ,crm_h_qcl(local_new_x,local_new_y,mlevs)      &
 ,crm_h_qcf(local_new_x,local_new_y,mlevs)      &
 ,crm_h_qrain(local_new_x,local_new_y,mlevs)    &
 ,crm_h_qgraup(local_new_x,local_new_y,mlevs)   &
 ,crm_h_dt1(local_new_x,local_new_y,mlevs)      &
 ,crm_h_dt2(local_new_x,local_new_y,mlevs)      &
 ,crm_h_dt4(local_new_x,local_new_y,mlevs)      &
 ,crm_h_dt9(local_new_x,local_new_y,mlevs)      &
 ,crm_h_dt12(local_new_x,local_new_y,mlevs)     &
 ,crm_h_dt30(local_new_x,local_new_y,mlevs)     &
 ,crm_h_dq4(local_new_x,local_new_y,mlevs)      &
 ,crm_h_dq9(local_new_x,local_new_y,mlevs)      &
 ,crm_h_dq12(local_new_x,local_new_y,mlevs)     &
 ,crm_h_dq30(local_new_x,local_new_y,mlevs)     &
 ,crm_h_dqcl4(local_new_x,local_new_y,mlevs)    &
 ,crm_h_dqcl9(local_new_x,local_new_y,mlevs)    &
 ,crm_h_dqcl12(local_new_x,local_new_y,mlevs)   &
 ,crm_h_dqcl30(local_new_x,local_new_y,mlevs)   &
 ,crm_h_dqcf4(local_new_x,local_new_y,mlevs)    &
 ,crm_h_dqcf3(local_new_x,local_new_y,mlevs)    &
 ,crm_h_dqcf12(local_new_x,local_new_y,mlevs)   &
 ,crm_h_dqcf30(local_new_x,local_new_y,mlevs)   &
 ,crm_h_dqrain30(local_new_x,local_new_y,mlevs) &
 ,crm_h_dqgr30(local_new_x,local_new_y,mlevs)   &
 ,crm_h_drho(local_new_x,local_new_y,mlevs)     &
 ,crm_h_thw(local_new_x,local_new_y,mlevs)      &
 ,crm_h_thvw(local_new_x,local_new_y,mlevs)     &
 ,crm_h_qw(local_new_x,local_new_y,mlevs)       &
 ,crm_h_qclw(local_new_x,local_new_y,mlevs)     &
 ,crm_h_qcfw(local_new_x,local_new_y,mlevs)     &
 ,crm_h_qrainw(local_new_x,local_new_y,mlevs)   &
 ,crm_h_qgraupw(local_new_x,local_new_y,mlevs)  &
 ,crm_h_uth(local_new_x,local_new_y,mlevs)      &
 ,crm_h_vth(local_new_x,local_new_y,mlevs)      &
 ,crm_h_uthv(local_new_x,local_new_y,mlevs)     &
 ,crm_h_vthv(local_new_x,local_new_y,mlevs)     &
 ,crm_h_uq(local_new_x,local_new_y,mlevs)       &
 ,crm_h_vq(local_new_x,local_new_y,mlevs)       &
 ,crm_h_wp(local_new_x,local_new_y,mlevs)       &
 ,crm_h_uw(local_new_x,local_new_y,mlevs)       &
 ,crm_h_vw(local_new_x,local_new_y,mlevs)       &
 ,crm_h_ww(local_new_x,local_new_y,mlevs)       &
 ,crm_h_w3(local_new_x,local_new_y,mlevs)       &
 ,crm_h_uu(local_new_x,local_new_y,mlevs)       &
 ,crm_h_vv(local_new_x,local_new_y,mlevs)       &
 ,crm_h_uv(local_new_x,local_new_y,mlevs)

!---------------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------------

INTEGER :: i,j,k            ! loop counters


CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_ZERO_ARRAYS"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Zero arrays
!-------------------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(k, i, j) DEFAULT(SHARED)

DO k=1,mlevs
  ! loop over coarse grid
  DO j = 1,local_new_y
    DO i = 1, local_new_x
      crm_factor(i,j,k)=0.0 
      crm_h_w(i,j,k)=0.0
      crm_h_rho(i,j,k)=0.0
      crm_h_u(i,j,k)=0.0
      crm_h_v(i,j,k)=0.0
      crm_h_th(i,j,k)=0.0
      crm_h_thv(i,j,k)=0.0
      crm_h_rh(i,j,k)=0.0
      crm_h_dt1(i,j,k)=0.0
      crm_h_dt2(i,j,k)=0.0
      crm_h_dt4(i,j,k)=0.0
      crm_h_dt9(i,j,k)=0.0
      crm_h_dt12(i,j,k)=0.0
      crm_h_dq4(i,j,k)=0.0
      crm_h_dq9(i,j,k)=0.0
      crm_h_dq12(i,j,k)=0.0
      crm_h_dqcl4(i,j,k)=0.0
      crm_h_dqcl9(i,j,k)=0.0
      crm_h_dqcl12(i,j,k)=0.0
      crm_h_dqcf4(i,j,k)=0.0
      crm_h_dqcf3(i,j,k)=0.0
      crm_h_dqcf12(i,j,k)=0.0
      crm_h_q(i,j,k)=0.0
      crm_h_qcl(i,j,k)=0.0
      crm_h_qcf(i,j,k)=0.0
      crm_h_qrain(i,j,k)=0.0
      crm_h_a(i,j,k)=0.0
      crm_h_thw(i,j,k) = 0.0
      crm_h_thvw(i,j,k) = 0.0
      crm_h_qw(i,j,k) = 0.0
      crm_h_qclw(i,j,k) = 0.0
      crm_h_qcfw(i,j,k) = 0.0
      crm_h_qrainw(i,j,k) = 0.0
      crm_h_uw(i,j,k) = 0.0
      crm_h_vw(i,j,k) = 0.0
      crm_h_uu(i,j,k) = 0.0
      crm_h_vv(i,j,k) = 0.0
      crm_h_uv(i,j,k) = 0.0
      crm_h_ww(i,j,k) = 0.0
      crm_h_w3(i,j,k) = 0.0
      crm_h_dpx(i,j,k)=0.0
      crm_h_dpy(i,j,k)=0.0
      crm_h_uth(i,j,k) = 0.0
      crm_h_vth(i,j,k) = 0.0
      crm_h_uthv(i,j,k) = 0.0
      crm_h_vthv(i,j,k) = 0.0
      crm_h_uq(i,j,k) = 0.0
      crm_h_vq(i,j,k) = 0.0
      crm_h_wp(i,j,k) = 0.0
    END DO
  END DO
  IF (l_qgraup) THEN
    DO j=1,local_new_y
      DO i=1,local_new_x
        crm_h_qgraup(i,j,k)  = 0.0
        crm_h_qgraupw(i,j,k) = 0.0
      END DO
    END DO
  END IF
  IF (l_sect30) THEN
    DO j=1,local_new_y
      DO i=1,local_new_x
        crm_h_dt30(i,j,k)=0.0
        crm_h_dq30(i,j,k)=0.0
        crm_h_dqcl30(i,j,k)=0.0
        crm_h_dqcf30(i,j,k)=0.0
        crm_h_dqrain30(i,j,k)=0.0
        crm_h_drho(i,j,k)=0.0
      END DO
    END DO  
  END IF
  IF (l_sect30 .AND. l_qgraup ) THEN
    DO j=1,local_new_y
      DO i=1,local_new_x
        crm_h_dqgr30(i,j,k)=0.0
      END DO
    END DO  
  END IF

END DO   ! loop over k
!$OMP END PARALLEL DO

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_zero_arrays

END MODULE crmstyle_zero_arrays_mod

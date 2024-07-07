! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  scale fields and set missing data where no values

MODULE crmstyle_scale_rmdi_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_SCALE_RMDI_MOD'

CONTAINS

! ------------------------------------------------------------------------------
! Description:
!  Calculates mean value of partition or sets missing data if no value
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------


SUBROUTINE crmstyle_scale_rmdi(all_factor, crm_factor,                    &  
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

USE missing_data_mod, ONLY: rmdi

USE umPrintMgr                              ! for writing output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

REAL, INTENT(IN) ::                          &
  all_factor(local_new_x,local_new_y,mlevs)  & ! Sum of points above surface
 ,crm_factor(local_new_x,local_new_y,mlevs)    ! Sum of points for partition

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

INTEGER :: i,j,k,ii,jj           ! loop counters

REAL ::                                        &
 factor, ftot


CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_SCALE_RMDI"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

! Number of points in each coarse grid
ftot = REAL(new_res(1)*new_res(1))

!-------------------------------------------------------------------------------
! Create means and set missing data
!-------------------------------------------------------------------------------


! divide by total number of points

!$OMP PARALLEL DO PRIVATE(k, ii, jj, factor) DEFAULT(SHARED)

DO k=1,mlevs
  ! loop over coarse grid
  DO jj = 1,local_new_y
    DO ii = 1, local_new_x

      IF (all_factor(ii,jj,k) > 0.0) THEN
        ! Field remains zero if above surface but no points 
        IF (crm_factor(ii,jj,k) > 0.0) THEN
          factor = 1.0/crm_factor(ii,jj,k)
          crm_h_a(ii,jj,k) = crm_factor(ii,jj,k)/ftot
          crm_h_rho(ii,jj,k) = crm_h_rho(ii,jj,k)/factor
          crm_h_w(ii,jj,k) = crm_h_w(ii,jj,k)*factor
          crm_h_u(ii,jj,k) = crm_h_u(ii,jj,k)*factor
          crm_h_v(ii,jj,k) = crm_h_v(ii,jj,k)*factor
          crm_h_th(ii,jj,k) = crm_h_th(ii,jj,k)*factor
          crm_h_thv(ii,jj,k) = crm_h_thv(ii,jj,k)*factor
          crm_h_rh(ii,jj,k) = crm_h_rh(ii,jj,k)*factor
          crm_h_dpx(ii,jj,k) = crm_h_dpx(ii,jj,k)*factor
          crm_h_dpy(ii,jj,k) = crm_h_dpy(ii,jj,k)*factor
          crm_h_q(ii,jj,k) = crm_h_q(ii,jj,k)*factor
          crm_h_qcl(ii,jj,k) = crm_h_qcl(ii,jj,k)*factor
          crm_h_qcf(ii,jj,k) = crm_h_qcf(ii,jj,k)*factor
          crm_h_qrain(ii,jj,k) = crm_h_qrain(ii,jj,k)*factor
          crm_h_thw(ii,jj,k) = crm_h_thw(ii,jj,k)*factor
          crm_h_thvw(ii,jj,k) = crm_h_thvw(ii,jj,k)*factor
          crm_h_qw(ii,jj,k) = crm_h_qw(ii,jj,k)*factor
          crm_h_qclw(ii,jj,k) = crm_h_qclw(ii,jj,k)*factor
          crm_h_qcfw(ii,jj,k) = crm_h_qcfw(ii,jj,k)*factor
          crm_h_qrainw(ii,jj,k) = crm_h_qrainw(ii,jj,k)*factor
          crm_h_ww(ii,jj,k) = crm_h_ww(ii,jj,k)*factor
          crm_h_w3(ii,jj,k) = crm_h_w3(ii,jj,k)*factor
          crm_h_uw(ii,jj,k) = crm_h_uw(ii,jj,k)*factor
          crm_h_vw(ii,jj,k) = crm_h_vw(ii,jj,k)*factor
          crm_h_uu(ii,jj,k) = crm_h_uu(ii,jj,k)*factor
          crm_h_vv(ii,jj,k) = crm_h_vv(ii,jj,k)*factor
          crm_h_uv(ii,jj,k) = crm_h_uv(ii,jj,k)*factor
          crm_h_dt1(ii,jj,k) = crm_h_dt1(ii,jj,k)*factor
          crm_h_dt2(ii,jj,k) = crm_h_dt2(ii,jj,k)*factor
          crm_h_dt4(ii,jj,k) = crm_h_dt4(ii,jj,k)*factor
          crm_h_dt9(ii,jj,k) = crm_h_dt9(ii,jj,k)*factor
          crm_h_dt12(ii,jj,k) = crm_h_dt12(ii,jj,k)*factor
          crm_h_dq4(ii,jj,k) = crm_h_dq4(ii,jj,k)*factor
          crm_h_dq9(ii,jj,k) = crm_h_dq9(ii,jj,k)*factor
          crm_h_dq12(ii,jj,k) = crm_h_dq12(ii,jj,k)*factor
          crm_h_dqcl4(ii,jj,k) = crm_h_dqcl4(ii,jj,k)*factor
          crm_h_dqcl9(ii,jj,k) = crm_h_dqcl9(ii,jj,k)*factor
          crm_h_dqcl12(ii,jj,k) = crm_h_dqcl12(ii,jj,k)*factor
          crm_h_dqcf4(ii,jj,k) = crm_h_dqcf4(ii,jj,k)*factor
          crm_h_dqcf3(ii,jj,k) = crm_h_dqcf3(ii,jj,k)*factor
          crm_h_dqcf12(ii,jj,k) = crm_h_dqcf12(ii,jj,k)*factor
          crm_h_uth(ii,jj,k)   = crm_h_uth(ii,jj,k)*factor
          crm_h_uthv(ii,jj,k)  = crm_h_uthv(ii,jj,k)*factor
          crm_h_vth(ii,jj,k)   = crm_h_vth(ii,jj,k)*factor
          crm_h_vthv(ii,jj,k)  = crm_h_vthv(ii,jj,k)*factor
          crm_h_uq(ii,jj,k)    = crm_h_uq(ii,jj,k)*factor
          crm_h_vq(ii,jj,k)    = crm_h_vq(ii,jj,k)*factor
          crm_h_wp(ii,jj,k)    = crm_h_wp(ii,jj,k)*factor
          IF (l_qgraup) THEN
            crm_h_qgraup(ii,jj,k)  = crm_h_qgraup(ii,jj,k)*factor
            crm_h_qgraupw(ii,jj,k) = crm_h_qgraupw(ii,jj,k)*factor
          END IF
          IF (l_sect30) THEN
            crm_h_dt30(ii,jj,k) = crm_h_dt30(ii,jj,k)*factor
            crm_h_dq30(ii,jj,k) = crm_h_dq30(ii,jj,k)*factor
            crm_h_dqcl30(ii,jj,k) = crm_h_dqcl30(ii,jj,k)*factor
            crm_h_dqcf30(ii,jj,k) = crm_h_dqcf30(ii,jj,k)*factor
            crm_h_dqrain30(ii,jj,k) = crm_h_dqrain30(ii,jj,k)*factor
            crm_h_drho(ii,jj,k) = crm_h_drho(ii,jj,k)*factor
          END IF
          IF (l_sect30 .AND. l_qgraup ) THEN
            crm_h_dqgr30(ii,jj,k) = crm_h_dqgr30(ii,jj,k)*factor 
          END IF
        END IF    

      ELSE      ! all points below surface
        crm_h_a(ii,jj,k) = rmdi
        crm_h_rho(ii,jj,k) = rmdi
        crm_h_w(ii,jj,k) =  rmdi
        crm_h_u(ii,jj,k) =  rmdi
        crm_h_v(ii,jj,k) =  rmdi
        crm_h_th(ii,jj,k) =  rmdi
        crm_h_thv(ii,jj,k) = rmdi
        crm_h_rh(ii,jj,k) =  rmdi
        crm_h_dpx(ii,jj,k) = rmdi
        crm_h_dpy(ii,jj,k) = rmdi 
        crm_h_q(ii,jj,k) =  rmdi
        crm_h_qcl(ii,jj,k) =  rmdi
        crm_h_qcf(ii,jj,k) =  rmdi
        crm_h_qrain(ii,jj,k) = rmdi
        crm_h_thw(ii,jj,k) =  rmdi
        crm_h_thvw(ii,jj,k) =  rmdi
        crm_h_qw(ii,jj,k) =  rmdi
        crm_h_qclw(ii,jj,k) =   rmdi
        crm_h_qcfw(ii,jj,k) =   rmdi
        crm_h_qrainw(ii,jj,k) =   rmdi
        crm_h_ww(ii,jj,k) =   rmdi
        crm_h_w3(ii,jj,k) =   rmdi
        crm_h_uw(ii,jj,k) =   rmdi
        crm_h_vw(ii,jj,k) =   rmdi
        crm_h_uu(ii,jj,k) =  rmdi
        crm_h_vv(ii,jj,k) =  rmdi
        crm_h_uv(ii,jj,k) =  rmdi
        crm_h_dt1(ii,jj,k) = rmdi
        crm_h_dt2(ii,jj,k) = rmdi
        crm_h_dt4(ii,jj,k) = rmdi
        crm_h_dt9(ii,jj,k) = rmdi
        crm_h_dt12(ii,jj,k) = rmdi
        crm_h_dq4(ii,jj,k) = rmdi
        crm_h_dq9(ii,jj,k) = rmdi
        crm_h_dq12(ii,jj,k) = rmdi
        crm_h_dqcl4(ii,jj,k) = rmdi
        crm_h_dqcl9(ii,jj,k) = rmdi
        crm_h_dqcl12(ii,jj,k) = rmdi
        crm_h_dqcf4(ii,jj,k) = rmdi
        crm_h_dqcf3(ii,jj,k) = rmdi
        crm_h_dqcf12(ii,jj,k) = rmdi
        crm_h_uth(ii,jj,k)   = rmdi
        crm_h_uthv(ii,jj,k)  = rmdi
        crm_h_vth(ii,jj,k)   = rmdi
        crm_h_vthv(ii,jj,k)  = rmdi
        crm_h_uq(ii,jj,k)    = rmdi
        crm_h_vq(ii,jj,k)    = rmdi
        crm_h_wp(ii,jj,k)    = rmdi

        IF (l_qgraup) THEN
          crm_h_qgraup(ii,jj,k) = rmdi
          crm_h_qgraupw(ii,jj,k) = rmdi
        END IF
        IF (l_sect30) THEN
          crm_h_dt30(ii,jj,k) = rmdi
          crm_h_dq30(ii,jj,k) = rmdi
          crm_h_dqcl30(ii,jj,k) = rmdi
          crm_h_dqcf30(ii,jj,k) = rmdi
          crm_h_dqrain30(ii,jj,k) = rmdi
          crm_h_drho(ii,jj,k) = rmdi
        END IF
        IF (l_sect30 .AND. l_qgraup ) THEN
          crm_h_dqgr30(ii,jj,k) = rmdi
        END IF
      END IF

    END DO   ! loop over ii
  END DO     ! loop over jj
END DO   ! loop over k
!$OMP END PARALLEL DO

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_scale_rmdi

END MODULE crmstyle_scale_rmdi_mod

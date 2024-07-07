! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  Calculate means and s.d plus deviation from mean

MODULE crmstyle_all_means_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_ALL_MEANS_MOD'

CONTAINS

! ------------------------------------------------------------------------------
! Description:
!  Calculate means and s.d plus deviation from mean for all 3d fields
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CRMstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------


SUBROUTINE crmstyle_all_means(mask)

USE crmwork_arrays_mod, ONLY:                                             &
   index_col, index_row

USE crmstyle_cntl_mod, ONLY:                                              &
   mlevs, new_res, l_qgraup, l_sect30, l_cape, l_pcape

USE crmstyle_grid_info_mod, ONLY:                                        &
  local_new_x, local_new_y, local_row_len, local_rows

USE hires_data_mod, ONLY:                                                 &
  orog, landsea, precip, density, snow, rain, zh, lh, sh, pstar, tstar,   &
  theta, thetav, q, qcl, qcf, qrain, qgraup, p_theta_lev, p_rho_lev,      &
  exner, u, v, w, dpdx, dpdy, t, rh,                                      &
  dt1, dt2, dt4, dt9, dt12, dt30,                                         &
  dq4, dq9, dq12, dq30, dqcl4, dqcl9, dqcl12, dqcl30,                     &
  dqcf4, dqcf3, dqcf12, dqcf30

USE crmstyle_sample_arrays_mod, ONLY:                                       &
  all_zh, all_sh, all_lh,  all_pstar, all_tstar, all_rain,                  &
  all_snow, all_precip, all_orog,  all_land,                                &
  all_sd_zh, all_sd_sh, all_sd_lh,  all_sd_pstar, all_sd_tstar, all_sd_rain,&
  all_sd_snow, all_sd_precip, all_sd_orog,  all_sd_land,                    &
  all_w, all_u,all_v, all_th, all_thv, all_rho,                             &
  all_q, all_qcl, all_qcf, all_qrain, all_qgraup,                           &
  all_ptheta, all_rh, all_t,                                                &
  all_dt1, all_dt2, all_dt4, all_dt9, all_dt12,all_dt30,                    &
  all_dq4, all_dq9, all_dq12, all_dq30,                                     &
  all_dqcl4, all_dqcl9, all_dqcl12, all_dqcl30,                             &
  all_dqcf4, all_dqcf3, all_dqcf12, all_dqcf30,                             &
  all_dqrain30, all_dqgr30, all_drho,                                       &
  all_sd_w, all_sd_u,all_sd_v, all_sd_th, all_sd_t, all_sd_thv, all_sd_rho, &
  all_sd_q, all_sd_qcl, all_sd_qcf, all_sd_qrain, all_sd_qgraup,            &
  all_sd_ptheta, all_sd_rh,                                                 &
  all_sd_dt1, all_sd_dt2, all_sd_dt4, all_sd_dt9, all_sd_dt12,all_sd_dt30,  &
  all_sd_dq4, all_sd_dq9, all_sd_dq12, all_sd_dq30,                         &
  all_sd_dqcl4, all_sd_dqcl9, all_sd_dqcl12, all_sd_dqcl30,                 &
  all_sd_dqcf4, all_sd_dqcf3, all_sd_dqcf12, all_sd_dqcf30,                 &
  all_sd_dqrain30, all_sd_dqgr30, all_sd_drho,                              &
  u_h_prime, v_h_prime, w_h_prime, th_h_prime, thv_h_prime, p_h_prime,      &
  q_h_prime, qcl_h_prime, qcf_h_prime, qrain_h_prime, qgraup_h_prime,       &
  all_exner, all_tv

USE planet_constants_mod, ONLY: kappa, c_virtual, pref

USE missing_data_mod, ONLY: rmdi
USE word_sizes_mod, ONLY: iwp,wp       ! Allows use of 4 byte words

! Subroutines
USE mean_f3d_mask_mod
USE mean_f3d_large_mod
USE mean_f2d_mod


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

LOGICAL, INTENT(IN) ::                    &
  mask(local_row_len,local_rows,mlevs)      ! mask true if above surface

!---------------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------------

INTEGER :: i,j,k        ! loop counters

REAL(wp) ::                                    &
  work_h_prime(local_row_len,local_rows,mlevs)   ! work array for case where f'
                                                 !  not requried later

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CRMSTYLE_ALL_MEANS'

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! mean single fields
!-------------------------------------------------------------------------------
CALL mean_f2d(local_row_len,local_rows,local_new_y,local_new_x,new_res(1), &
              index_row,index_col,                                         &
                  zh, all_zh, all_sd_zh)
CALL mean_f2d(local_row_len,local_rows,local_new_y,local_new_x,new_res(1), &
              index_row,index_col,                                         &
                  lh, all_lh, all_sd_lh)
CALL mean_f2d(local_row_len,local_rows,local_new_y,local_new_x,new_res(1), &
              index_row,index_col,                                         &
                  sh, all_sh, all_sd_sh)
CALL mean_f2d(local_row_len,local_rows,local_new_y,local_new_x,new_res(1), &
              index_row,index_col,                                         &
                  pstar, all_pstar, all_sd_pstar)
CALL mean_f2d(local_row_len,local_rows,local_new_y,local_new_x,new_res(1), &
              index_row,index_col,                                         &
                  tstar, all_tstar, all_sd_tstar)
CALL mean_f2d(local_row_len,local_rows,local_new_y,local_new_x,new_res(1), &
              index_row,index_col,                                         &
                  rain, all_rain, all_sd_rain)
CALL mean_f2d(local_row_len,local_rows,local_new_y,local_new_x,new_res(1), &
              index_row,index_col,                                         &
                  snow, all_snow, all_sd_snow)
CALL mean_f2d(local_row_len,local_rows,local_new_y,local_new_x,new_res(1), &
              index_row,index_col,                                         &
                  precip, all_precip, all_sd_precip)
CALL mean_f2d(local_row_len,local_rows,local_new_y,local_new_x,new_res(1), &
              index_row,index_col,                                         &
                  orog, all_orog, all_sd_orog)
CALL mean_f2d(local_row_len,local_rows,local_new_y,local_new_x,new_res(1), &
              index_row,index_col,                                         &
                  landsea, all_land, all_sd_land)

!-------------------------------------------------------------------------------
! mean each field
!-------------------------------------------------------------------------------

CALL mean_f3d_large(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, theta,                  &
                         th_h_prime, all_th, all_sd_th)


CALL mean_f3d_large(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, thetav,                 &
                         thv_h_prime, all_thv, all_sd_thv)


CALL mean_f3d_large(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, t,                      &
                         work_h_prime, all_t, all_sd_t)

CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, w,                      &
                         w_h_prime, all_w, all_sd_w)

CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, u,                      &
                         u_h_prime, all_u, all_sd_u)

CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, v,                      &
                         v_h_prime, all_v, all_sd_v)

CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, q,                      &
                         q_h_prime, all_q, all_sd_q)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, qcl,                    &
                         qcl_h_prime, all_qcl, all_sd_qcl)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, qcf,                    &
                         qcf_h_prime, all_qcf, all_sd_qcf)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, qrain,                  &
                         qrain_h_prime, all_qrain, all_sd_qrain)

IF (l_qgraup) THEN
  CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                           index_row,index_col,mask, qgraup,                 &
                           qgraup_h_prime, all_qgraup, all_sd_qgraup)

END IF

CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, density,                &
                         work_h_prime, all_rho, all_sd_rho)

CALL mean_f3d_large(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, p_theta_lev,             &
                         p_h_prime, all_ptheta, all_sd_ptheta)

CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, rh,                     &
                         work_h_prime, all_rh, all_sd_rh)

! Tendencies
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dt1,                    &
                         work_h_prime, all_dt1, all_sd_dt1)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dt2,                    &
                         work_h_prime, all_dt2, all_sd_dt2)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dt4,                    &
                         work_h_prime, all_dt4, all_sd_dt4)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dt9,                    &
                         work_h_prime, all_dt9, all_sd_dt9)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dt12,                   &
                         work_h_prime, all_dt12, all_sd_dt12)

CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dq4,                    &
                         work_h_prime, all_dq4, all_sd_dq4)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dq9,                    &
                         work_h_prime, all_dq9, all_sd_dq9)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dq12,                   &
                         work_h_prime, all_dq12, all_sd_dq12)

CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dqcl4,                  &
                         work_h_prime, all_dqcl4, all_sd_dqcl4)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dqcl9,                  &
                         work_h_prime, all_dqcl9, all_sd_dqcl9)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dqcl12,                 &
                         work_h_prime, all_dqcl12, all_sd_dqcl12)

CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dqcf4,                  &
                         work_h_prime, all_dqcf4, all_sd_dqcf4)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dqcf3,                  &
                         work_h_prime, all_dqcf3, all_sd_dqcf3)
CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                         index_row,index_col,mask, dqcf12,                 &
                         work_h_prime, all_dqcf12, all_sd_dqcf12)

IF (l_sect30) THEN
  CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                           index_row,index_col,mask, dt30,                   &
                           work_h_prime, all_dt30, all_sd_dt30)
  CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                           index_row,index_col,mask, dq30,                   &
                           work_h_prime, all_dq30, all_sd_dq30)
  CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                           index_row,index_col,mask, dqcl30,                 &
                           work_h_prime, all_dqcl30, all_sd_dqcl30)
  CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                           index_row,index_col,mask, dqcf30,                 &
                           work_h_prime, all_dqcf30, all_sd_dqcf30)
  CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                           index_row,index_col,mask, dqcf30,                 &
                           work_h_prime, all_dqrain30, all_sd_dqrain30)
  CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                           index_row,index_col,mask, dqcf30,                 &
                           work_h_prime, all_drho, all_sd_drho)
END IF

IF (l_sect30 .AND. l_qgraup) THEN
  CALL mean_f3d_mask(local_row_len,local_rows,mlevs,local_new_y,local_new_x, &
                           index_row,index_col,mask, dqcf30,                 &
                           work_h_prime, all_dqgr30, all_sd_dqgr30)
END IF

!------------------------------------------------------------------------------
! Calculate additional mean fields for use later from other mean fields
!------------------------------------------------------------------------------

IF (l_cape .OR. l_pcape) THEN

! Assume zero if below surface
  DO k=1,mlevs
    DO j=1,local_new_y
      DO i=1,local_new_x
        all_exner(i,j,k) = all_t(i,j,k)/all_th(i,j,k)
        all_tv(i,j,k) = all_t(i,j,k) *(1.0 + c_virtual*all_q(i,j,k))
      END DO
    END DO
  END DO

END IF

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_all_means

END MODULE crmstyle_all_means_mod

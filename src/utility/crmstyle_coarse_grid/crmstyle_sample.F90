! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  Calculate lots of different partitions for convective analysis

MODULE crmstyle_sample_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_SAMPLE_MOD'

CONTAINS

! ------------------------------------------------------------------------------
! Description:
!  Calculate lots of different partitions for convective analysis.
! all - all points
! acc - all cloudy points
! acu - all cloudy updraughts (w'>0)
! bcu - buoyant cloudy updraughts (w'>0)
! wg1 - strong (w'>1) buoyant cloudy updraughts
! ppd - precipitating downdraughts (rain+graupel) (w' < 0.0)
! nbd - Negatively buoyant precipitating downdraughts (w' < 0.0)
! nid - Negatively buoyant precipitating +ice downdraughts relative
! adu - upward dry air (no cloud) upward not relative
! acw - upward cloudy air upward not relative (w>0.0)
! bcw - upward buoyant cloudy air upward not relative  (w>0.0)
! ucu - unstable cloudy updraughts (-dln(thesat)/dz > 0.0, w > 0.0)
! ppw - precipitating downdraughts (rain+graupel) (w < 0.0)
! nbw - Negatively buoyant precipitating downdraughts (w < 0.0)
! 
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------


SUBROUTINE crmstyle_sample(mask)

USE crmwork_arrays_mod, ONLY:                      &
   index_col, index_row

USE crmstyle_cntl_mod, ONLY:                                             &
   mlevs, new_res, l_qgraup, l_sect30, iprint,                           &
   l_bcu, l_wg1, l_acc, l_acu, l_ppd, l_nbd, l_nid, l_adu, l_acw, l_bcw, &
   l_ucu, l_ppw, l_nbw, nbins_diam, nbins_size, nbins_fract, nbins_dxdy

USE crmstyle_grid_info_mod, ONLY:                                 &
   local_new_x, local_new_y, local_row_len, local_rows

USE hires_data_mod, ONLY:                                              &
  orog, landsea, precip,                                               &
  theta, thetav, q, qcl, qcf, qrain, qgraup, p_theta_lev, p_rho_lev,   &
  exner, u, v, w, dpdx, dpdy, t, rh, density,                          &
  dt1, dt2, dt4, dt9, dt12, dt30,                                      &
  dq4, dq9, dq12, dq30, dqcl4, dqcl9, dqcl12, dqcl30,                  &
  dqcf4, dqcf3, dqcf12, dqcf30, dqrain30, dqgr30, drho

USE crmstyle_sample_arrays_mod, ONLY:                                     &
  all_w, all_u,all_v, all_th, all_thv, all_rho,                           &
  all_q, all_qcl, all_qcf, all_qrain, all_qgraup,                         &
  all_ptheta, all_rh, all_a, all_dpx, all_dpy,                            &
  all_dt1, all_dt2, all_dt4, all_dt9, all_dt12,all_dt30,                  &
  all_dq4, all_dq9, all_dq12, all_dq30,                                   &
  all_dqcl4, all_dqcl9, all_dqcl12, all_dqcl30,                           &
  all_dqcf4, all_dqcf3, all_dqcf12, all_dqcf30,                           &
  all_dqrain30, all_dqgr30, all_drho,                                     &
  all_thw, all_thvw,all_wp_hydro, p_theta_hydro,                          &
  all_qw, all_qclw, all_qcfw, all_qrainw, all_qgraupw,                    &
  all_uw, all_vw, all_ww, all_w3, all_vv, all_uu, all_uv,                 &
  all_uth, all_vth, all_uthv, all_vthv, all_uq, all_vq, all_wp,           &
  acc_h_w, acc_h_u, acc_h_v, acc_h_th, acc_h_thv, acc_h_rho,              &
  acc_h_rh, acc_h_a, acc_h_dpx, acc_h_dpy,                                &
  acc_h_q,   acc_h_qcl, acc_h_qcf, acc_h_qrain, acc_h_qgraup,             &
  acc_h_dt1, acc_h_dt2, acc_h_dt4, acc_h_dt9, acc_h_dt12, acc_h_dt30,     &
  acc_h_dq4, acc_h_dq9, acc_h_dq12, acc_h_dq30,                           &
  acc_h_dqcl4, acc_h_dqcl9, acc_h_dqcl12, acc_h_dqcl30,                   &
  acc_h_dqcf4, acc_h_dqcf3, acc_h_dqcf12, acc_h_dqcf30,                   &
  acc_h_dqrain30, acc_h_dqgr30, acc_h_drho,                               &
  acc_h_thw, acc_h_thvw, acc_h_wp,                                        &
  acc_h_qw, acc_h_qclw, acc_h_qcfw, acc_h_qrainw, acc_h_qgraupw,          &
  acc_h_uw, acc_h_vw, acc_h_ww, acc_h_w3, acc_h_vv, acc_h_uu, acc_h_uv,   &
  acc_h_uth, acc_h_vth, acc_h_uthv, acc_h_vthv, acc_h_uq, acc_h_vq,       &
  acu_h_w, acu_h_u, acu_h_v, acu_h_th, acu_h_thv, acu_h_rho,              &
  acu_h_rh, acu_h_a, acu_h_dpx, acu_h_dpy,                                &
  acu_h_q,   acu_h_qcl, acu_h_qcf, acu_h_qrain, acu_h_qgraup,             &
  acu_h_dt1, acu_h_dt2, acu_h_dt4, acu_h_dt9, acu_h_dt12, acu_h_dt30,     &
  acu_h_dq4, acu_h_dq9, acu_h_dq12, acu_h_dq30,                           &
  acu_h_dqcl4, acu_h_dqcl9, acu_h_dqcl12, acu_h_dqcl30,                   &
  acu_h_dqcf4, acu_h_dqcf3, acu_h_dqcf12, acu_h_dqcf30,                   &
  acu_h_dqrain30, acu_h_dqgr30, acu_h_drho,                               &
  acu_h_thw, acu_h_thvw, acu_h_wp,                                        &
  acu_h_qw, acu_h_qclw, acu_h_qcfw, acu_h_qrainw, acu_h_qgraupw,          &
  acu_h_uw, acu_h_vw, acu_h_ww, acu_h_w3, acu_h_vv, acu_h_uu, acu_h_uv,   &
  acu_h_uth, acu_h_vth, acu_h_uthv, acu_h_vthv, acu_h_uq, acu_h_vq,       &
  bcu_h_w, bcu_h_u, bcu_h_v, bcu_h_th, bcu_h_thv, bcu_h_rho,              &
  bcu_h_rh, bcu_h_a, bcu_h_dpx, bcu_h_dpy,                                &
  bcu_h_q,   bcu_h_qcl, bcu_h_qcf, bcu_h_qrain, bcu_h_qgraup,             &
  bcu_h_dt1, bcu_h_dt2, bcu_h_dt4, bcu_h_dt9, bcu_h_dt12, bcu_h_dt30,     &
  bcu_h_dq4, bcu_h_dq9, bcu_h_dq12, bcu_h_dq30,                           &
  bcu_h_dqcl4, bcu_h_dqcl9, bcu_h_dqcl12, bcu_h_dqcl30,                   &
  bcu_h_dqcf4, bcu_h_dqcf3, bcu_h_dqcf12, bcu_h_dqcf30,                   &
  bcu_h_dqrain30, bcu_h_dqgr30, bcu_h_drho,                               &
  bcu_h_thw, bcu_h_thvw, bcu_h_wp,                                        &
  bcu_h_qw, bcu_h_qclw, bcu_h_qcfw, bcu_h_qrainw, bcu_h_qgraupw,          &
  bcu_h_uw, bcu_h_vw, bcu_h_ww, bcu_h_w3, bcu_h_vv, bcu_h_uu, bcu_h_uv,   &
  bcu_h_uth, bcu_h_vth, bcu_h_uthv, bcu_h_vthv, bcu_h_uq, bcu_h_vq,       &
  wg1_h_w, wg1_h_u, wg1_h_v, wg1_h_th, wg1_h_thv, wg1_h_rho,              &
  wg1_h_rh, wg1_h_a, wg1_h_dpx, wg1_h_dpy,                                &
  wg1_h_q,   wg1_h_qcl, wg1_h_qcf, wg1_h_qrain, wg1_h_qgraup,             &
  wg1_h_dt1, wg1_h_dt2, wg1_h_dt4, wg1_h_dt9, wg1_h_dt12, wg1_h_dt30,     &
  wg1_h_dq4, wg1_h_dq9, wg1_h_dq12, wg1_h_dq30,                           &
  wg1_h_dqcl4, wg1_h_dqcl9, wg1_h_dqcl12, wg1_h_dqcl30,                   &
  wg1_h_dqcf4, wg1_h_dqcf3, wg1_h_dqcf12, wg1_h_dqcf30,                   &
  wg1_h_dqrain30, wg1_h_dqgr30, wg1_h_drho,                               &
  wg1_h_thw, wg1_h_thvw, wg1_h_wp,                                        &
  wg1_h_qw, wg1_h_qclw, wg1_h_qcfw, wg1_h_qrainw, wg1_h_qgraupw,          &
  wg1_h_uw, wg1_h_vw, wg1_h_ww, wg1_h_w3, wg1_h_vv, wg1_h_uu, wg1_h_uv,   &
  wg1_h_uth, wg1_h_vth, wg1_h_uthv, wg1_h_vthv, wg1_h_uq, wg1_h_vq,       &
  ppd_h_w, ppd_h_u, ppd_h_v, ppd_h_th, ppd_h_thv, ppd_h_rho,              &
  ppd_h_rh, ppd_h_a, ppd_h_dpx, ppd_h_dpy,                                &
  ppd_h_q,   ppd_h_qcl, ppd_h_qcf, ppd_h_qrain, ppd_h_qgraup,             &
  ppd_h_dt1, ppd_h_dt2, ppd_h_dt4, ppd_h_dt9, ppd_h_dt12, ppd_h_dt30,     &
  ppd_h_dq4, ppd_h_dq9, ppd_h_dq12, ppd_h_dq30,                           &
  ppd_h_dqcl4, ppd_h_dqcl9, ppd_h_dqcl12, ppd_h_dqcl30,                   &
  ppd_h_dqcf4, ppd_h_dqcf3, ppd_h_dqcf12, ppd_h_dqcf30,                   &
  ppd_h_dqrain30, ppd_h_dqgr30, ppd_h_drho,                               &
  ppd_h_thw, ppd_h_thvw, ppd_h_wp,                                        &
  ppd_h_qw, ppd_h_qclw, ppd_h_qcfw, ppd_h_qrainw, ppd_h_qgraupw,          &
  ppd_h_uw, ppd_h_vw, ppd_h_ww, ppd_h_w3, ppd_h_vv, ppd_h_uu, ppd_h_uv,   &
  ppd_h_uth, ppd_h_vth, ppd_h_uthv, ppd_h_vthv, ppd_h_uq, ppd_h_vq,       &
  nbd_h_w, nbd_h_u, nbd_h_v, nbd_h_th, nbd_h_thv, nbd_h_rho,              &
  nbd_h_rh, nbd_h_a, nbd_h_dpx, nbd_h_dpy,                                &
  nbd_h_q,   nbd_h_qcl, nbd_h_qcf, nbd_h_qrain, nbd_h_qgraup,             &
  nbd_h_dt1, nbd_h_dt2, nbd_h_dt4, nbd_h_dt9, nbd_h_dt12, nbd_h_dt30,     &
  nbd_h_dq4, nbd_h_dq9, nbd_h_dq12, nbd_h_dq30,                           &
  nbd_h_dqcl4, nbd_h_dqcl9, nbd_h_dqcl12, nbd_h_dqcl30,                   &
  nbd_h_dqcf4, nbd_h_dqcf3, nbd_h_dqcf12, nbd_h_dqcf30,                   &
  nbd_h_dqrain30, nbd_h_dqgr30, nbd_h_drho,                               &
  nbd_h_thw, nbd_h_thvw, nbd_h_wp,                                        &
  nbd_h_qw, nbd_h_qclw, nbd_h_qcfw, nbd_h_qrainw, nbd_h_qgraupw,          &
  nbd_h_uw, nbd_h_vw, nbd_h_ww, nbd_h_w3, nbd_h_vv, nbd_h_uu, nbd_h_uv,   &
  nbd_h_uth, nbd_h_vth, nbd_h_uthv, nbd_h_vthv, nbd_h_uq, nbd_h_vq,       &
  nid_h_w, nid_h_u, nid_h_v, nid_h_th, nid_h_thv, nid_h_rho,              &
  nid_h_rh, nid_h_a, nid_h_dpx, nid_h_dpy,                                &
  nid_h_q,   nid_h_qcl, nid_h_qcf, nid_h_qrain, nid_h_qgraup,             &
  nid_h_dt1, nid_h_dt2, nid_h_dt4, nid_h_dt9, nid_h_dt12, nid_h_dt30,     &
  nid_h_dq4, nid_h_dq9, nid_h_dq12, nid_h_dq30,                           &
  nid_h_dqcl4, nid_h_dqcl9, nid_h_dqcl12, nid_h_dqcl30,                   &
  nid_h_dqcf4, nid_h_dqcf3, nid_h_dqcf12, nid_h_dqcf30,                   &
  nid_h_dqrain30, nid_h_dqgr30, nid_h_drho,                               &
  nid_h_thw, nid_h_thvw, nid_h_wp,                                        &
  nid_h_qw, nid_h_qclw, nid_h_qcfw, nid_h_qrainw, nid_h_qgraupw,          &
  nid_h_uw, nid_h_vw, nid_h_ww, nid_h_w3, nid_h_vv, nid_h_uu, nid_h_uv,   &
  nid_h_uth, nid_h_vth, nid_h_uthv, nid_h_vthv, nid_h_uq, nid_h_vq,       &
  adu_h_w, adu_h_u, adu_h_v, adu_h_th, adu_h_thv, adu_h_rho,              &
  adu_h_rh, adu_h_a, adu_h_dpx, adu_h_dpy,                                &
  adu_h_q,   adu_h_qcl, adu_h_qcf, adu_h_qrain, adu_h_qgraup,             &
  adu_h_dt1, adu_h_dt2, adu_h_dt4, adu_h_dt9, adu_h_dt12, adu_h_dt30,     &
  adu_h_dq4, adu_h_dq9, adu_h_dq12, adu_h_dq30,                           &
  adu_h_dqcl4, adu_h_dqcl9, adu_h_dqcl12, adu_h_dqcl30,                   &
  adu_h_dqcf4, adu_h_dqcf3, adu_h_dqcf12, adu_h_dqcf30,                   &
  adu_h_dqrain30, adu_h_dqgr30, adu_h_drho,                               &
  adu_h_thw, adu_h_thvw, adu_h_wp,                                        &
  adu_h_qw, adu_h_qclw, adu_h_qcfw, adu_h_qrainw, adu_h_qgraupw,          &
  adu_h_uw, adu_h_vw, adu_h_ww, adu_h_w3, adu_h_vv, adu_h_uu, adu_h_uv,   &
  adu_h_uth, adu_h_vth, adu_h_uthv, adu_h_vthv, adu_h_uq, adu_h_vq,       &
  acw_h_w, acw_h_u, acw_h_v, acw_h_th, acw_h_thv, acw_h_rho,              &
  acw_h_rh, acw_h_a, acw_h_dpx, acw_h_dpy,                                &
  acw_h_q,   acw_h_qcl, acw_h_qcf, acw_h_qrain, acw_h_qgraup,             &
  acw_h_dt1, acw_h_dt2, acw_h_dt4, acw_h_dt9, acw_h_dt12, acw_h_dt30,     &
  acw_h_dq4, acw_h_dq9, acw_h_dq12, acw_h_dq30,                           &
  acw_h_dqcl4, acw_h_dqcl9, acw_h_dqcl12, acw_h_dqcl30,                   &
  acw_h_dqcf4, acw_h_dqcf3, acw_h_dqcf12, acw_h_dqcf30,                   &
  acw_h_dqrain30, acw_h_dqgr30, acw_h_drho,                               &
  acw_h_thw, acw_h_thvw, acw_h_wp,                                        &
  acw_h_qw, acw_h_qclw, acw_h_qcfw, acw_h_qrainw, acw_h_qgraupw,          &
  acw_h_uw, acw_h_vw, acw_h_ww, acw_h_w3, acw_h_vv, acw_h_uu, acw_h_uv,   &
  acw_h_uth, acw_h_vth, acw_h_uthv, acw_h_vthv, acw_h_uq, acw_h_vq,       &
  bcw_h_w, bcw_h_u, bcw_h_v, bcw_h_th, bcw_h_thv, bcw_h_rho,              &
  bcw_h_rh, bcw_h_a, bcw_h_dpx, bcw_h_dpy,                                &
  bcw_h_q,   bcw_h_qcl, bcw_h_qcf, bcw_h_qrain, bcw_h_qgraup,             &
  bcw_h_dt1, bcw_h_dt2, bcw_h_dt4, bcw_h_dt9, bcw_h_dt12, bcw_h_dt30,     &
  bcw_h_dq4, bcw_h_dq9, bcw_h_dq12, bcw_h_dq30,                           &
  bcw_h_dqcl4, bcw_h_dqcl9, bcw_h_dqcl12, bcw_h_dqcl30,                   &
  bcw_h_dqcf4, bcw_h_dqcf3, bcw_h_dqcf12, bcw_h_dqcf30,                   &
  bcw_h_dqrain30, bcw_h_dqgr30, bcw_h_drho,                               &
  bcw_h_thw, bcw_h_thvw, bcw_h_wp,                                        &
  bcw_h_qw, bcw_h_qclw, bcw_h_qcfw, bcw_h_qrainw, bcw_h_qgraupw,          &
  bcw_h_uw, bcw_h_vw, bcw_h_ww, bcw_h_w3, bcw_h_vv, bcw_h_uu, bcw_h_uv,   &
  bcw_h_uth, bcw_h_vth, bcw_h_uthv, bcw_h_vthv, bcw_h_uq, bcw_h_vq,       &
  ucu_h_w, ucu_h_u, ucu_h_v, ucu_h_th, ucu_h_thv, ucu_h_rho,              &
  ucu_h_rh, ucu_h_a, ucu_h_dpx, ucu_h_dpy,                                &
  ucu_h_q,   ucu_h_qcl, ucu_h_qcf, ucu_h_qrain, ucu_h_qgraup,             &
  ucu_h_dt1, ucu_h_dt2, ucu_h_dt4, ucu_h_dt9, ucu_h_dt12, ucu_h_dt30,     &
  ucu_h_dq4, ucu_h_dq9, ucu_h_dq12, ucu_h_dq30,                           &
  ucu_h_dqcl4, ucu_h_dqcl9, ucu_h_dqcl12, ucu_h_dqcl30,                   &
  ucu_h_dqcf4, ucu_h_dqcf3, ucu_h_dqcf12, ucu_h_dqcf30,                   &
  ucu_h_dqrain30, ucu_h_dqgr30, ucu_h_drho,                               &
  ucu_h_thw, ucu_h_thvw, ucu_h_wp,                                        &
  ucu_h_qw, ucu_h_qclw, ucu_h_qcfw, ucu_h_qrainw, ucu_h_qgraupw,          &
  ucu_h_uw, ucu_h_vw, ucu_h_ww, ucu_h_w3, ucu_h_vv, ucu_h_uu, ucu_h_uv,   &
  ucu_h_uth, ucu_h_vth, ucu_h_uthv, ucu_h_vthv, ucu_h_uq, ucu_h_vq,       &
  ppw_h_w, ppw_h_u, ppw_h_v, ppw_h_th, ppw_h_thv, ppw_h_rho,              &
  ppw_h_rh, ppw_h_a, ppw_h_dpx, ppw_h_dpy,                                &
  ppw_h_q,   ppw_h_qcl, ppw_h_qcf, ppw_h_qrain, ppw_h_qgraup,             &
  ppw_h_dt1, ppw_h_dt2, ppw_h_dt4, ppw_h_dt9, ppw_h_dt12, ppw_h_dt30,     &
  ppw_h_dq4, ppw_h_dq9, ppw_h_dq12, ppw_h_dq30,                           &
  ppw_h_dqcl4, ppw_h_dqcl9, ppw_h_dqcl12, ppw_h_dqcl30,                   &
  ppw_h_dqcf4, ppw_h_dqcf3, ppw_h_dqcf12, ppw_h_dqcf30,                   &
  ppw_h_dqrain30, ppw_h_dqgr30, ppw_h_drho,                               &
  ppw_h_thw, ppw_h_thvw, ppw_h_wp,                                        &
  ppw_h_qw, ppw_h_qclw, ppw_h_qcfw, ppw_h_qrainw, ppw_h_qgraupw,          &
  ppw_h_uw, ppw_h_vw, ppw_h_ww, ppw_h_w3, ppw_h_vv, ppw_h_uu, ppw_h_uv,   &
  ppw_h_uth, ppw_h_vth, ppw_h_uthv, ppw_h_vthv, ppw_h_uq, ppw_h_vq,       &
  nbw_h_w, nbw_h_u, nbw_h_v, nbw_h_th, nbw_h_thv, nbw_h_rho,              &
  nbw_h_rh, nbw_h_a, nbw_h_dpx, nbw_h_dpy,                                &
  nbw_h_q,   nbw_h_qcl, nbw_h_qcf, nbw_h_qrain, nbw_h_qgraup,             &
  nbw_h_dt1, nbw_h_dt2, nbw_h_dt4, nbw_h_dt9, nbw_h_dt12, nbw_h_dt30,     &
  nbw_h_dq4, nbw_h_dq9, nbw_h_dq12, nbw_h_dq30,                           &
  nbw_h_dqcl4, nbw_h_dqcl9, nbw_h_dqcl12, nbw_h_dqcl30,                   &
  nbw_h_dqcf4, nbw_h_dqcf3, nbw_h_dqcf12, nbw_h_dqcf30,                   &
  nbw_h_dqrain30, nbw_h_dqgr30, nbw_h_drho,                               &
  nbw_h_thw, nbw_h_thvw, nbw_h_wp,                                        &
  nbw_h_qw, nbw_h_qclw, nbw_h_qcfw, nbw_h_qrainw, nbw_h_qgraupw,          &
  nbw_h_uw, nbw_h_vw, nbw_h_ww, nbw_h_w3, nbw_h_vv, nbw_h_uu, nbw_h_uv,   &
  nbw_h_uth, nbw_h_vth, nbw_h_uthv, nbw_h_vthv, nbw_h_uq, nbw_h_vq,       &
  u_h_prime, v_h_prime, w_h_prime, th_h_prime, thv_h_prime, p_h_prime,    &
  q_h_prime, qcl_h_prime, qcf_h_prime, qrain_h_prime, qgraup_h_prime,     &
  bcu_mask, bcu_mask_w, nn_mask,                                          &
  n_plume, plume_size, plume_diam_pdf

USE crmwork_arrays_mod, ONLY: h_theta_sea

USE word_sizes_mod, ONLY: iwp,wp   ! Allows use of 4 byte words to reduce
                                    ! memory

USE planet_constants_mod, ONLY: r, lcrcp, lsrcp
USE water_constants_mod, ONLY: tm

USE missing_data_mod, ONLY: rmdi

! subroutines
USE count_plumes_mod, ONLY: count_plumes
USE crmstyle_scale_rmdi_mod, ONLY: crmstyle_scale_rmdi
USE crmstyle_zero_arrays_mod, ONLY: crmstyle_zero_arrays

USE umPrintMgr                              ! for writing output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE


LOGICAL, INTENT(IN) ::                 &
  mask(local_row_len,local_rows,mlevs)   ! mask true if above surface

!---------------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------------

INTEGER :: i,j,k,ii,jj,ij          ! loop counters
INTEGER ::    &
  ix1,ix2     & ! start and stop for x gridpoints
 ,iy1,iy2       ! start and stop for y gridpoints

INTEGER ::    &
  nups          ! number of plumes

REAL ::                                        &
  acc_factor(local_new_x,local_new_y,mlevs)    &
 ,bcu_factor(local_new_x,local_new_y,mlevs)    &
 ,wg1_factor(local_new_x,local_new_y,mlevs)    &
 ,acu_factor(local_new_x,local_new_y,mlevs)    &
 ,ppd_factor(local_new_x,local_new_y,mlevs)    &
 ,nbd_factor(local_new_x,local_new_y,mlevs)    &
 ,nid_factor(local_new_x,local_new_y,mlevs)    &
 ,all_factor(local_new_x,local_new_y,mlevs)    &
 ,adu_factor(local_new_x,local_new_y,mlevs)    &
 ,acw_factor(local_new_x,local_new_y,mlevs)    &
 ,bcw_factor(local_new_x,local_new_y,mlevs)    &
 ,ucu_factor(local_new_x,local_new_y,mlevs)    &
 ,ppw_factor(local_new_x,local_new_y,mlevs)    &
 ,nbw_factor(local_new_x,local_new_y,mlevs)    &
 ,factor, ftot

LOGICAL ::                            &
  mask_temp(new_res(1),new_res(1))   ! copy of bcu mask for count_plume

REAL(wp) ::               &
  diam_pdf(nbins_diam)    &  ! pdf of diameters
 ,size_pdf(nbins_size)    &  ! pdf of sizes
 ,fract_pdf(nbins_fract)  &  ! pdf of fractions
 ,dxdy_pdf(nbins_dxdy)       ! pdf of directions

REAL(wp) ::                                  &
  qs_over_t(local_row_len,local_rows,mlevs)  & ! qsat / T    (kg/kgK)
 ,dthesdz(local_row_len,local_rows,mlevs)      ! -d(ln(thesat))/dz

REAL ::                          &
  qw_work                        & ! q'w'    for a level
 ,qclw_work                      & ! qcl'w'
 ,qcfw_work                      & ! qcf'w'
 ,qrainw_work                    & ! qrain'w'
 ,qgraupw_work                   & ! qgruap'w'
 ,thw_work                       & ! theta'w'
 ,thvw_work                      & ! thetav'w'
 ,uth_work                       & ! theta'u'
 ,uthv_work                      & ! thetav'u'
 ,vth_work                       & ! theta'v'
 ,vthv_work                      & ! thetav'v'
 ,uq_work                        & ! u'q'
 ,vq_work                        & ! v'q'
 ,ww_work                        & ! w'w'
 ,w3_work                        & ! w'w'w'
 ,wp_work                        & ! w'p'/density 
 ,wp_work2                       & ! w'p'/density 
 ,uu_work                        & ! u'u'
 ,vv_work                        & ! v'v'
 ,uv_work                        & ! v'v'
 ,uw_work                        & ! u'w'
 ,vw_work                        & ! v'w'
 ,dpx_work                       & ! rho*dp/dx
 ,dpy_work                         ! rho*dp/dy

REAL ::                      &
  precip_water               & ! total precipitation for downdraught tests
 ,qw_test                    & ! total condensate
 ,lrcp                       & ! lc/cp or (lc+lf)/cp
 ,rdz                        & ! 1/difference height between theta levels (/m)
 ,dqstdz                     & ! d(qsat/T)/dz     (kg/kgKm)
 ,dthdz                        ! 1/th * dTH/dz     (/Km)

REAL, PARAMETER :: qw_crit = 1.0e-5  ! critical cloud water/ice

CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_SAMPLE"

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
! Initialise arrays
!-------------------------------------------------------------------------------



!$OMP PARALLEL DO PRIVATE(k, i, j) DEFAULT(SHARED)

DO k=1,mlevs

  DO j=1,local_new_y
    DO i=1,local_new_x
      nn_mask(i,j,k) = 0           ! number of buoyant points in coarse area

      all_factor(i,j,k) =0.0
      all_rho(i,j,k) = 0.0
      all_thw(i,j,k) = 0.0
      all_thvw(i,j,k) = 0.0
      all_qw(i,j,k) = 0.0
      all_qclw(i,j,k) = 0.0
      all_qcfw(i,j,k) = 0.0
      all_qrainw(i,j,k) = 0.0
      all_uw(i,j,k) = 0.0
      all_vw(i,j,k) = 0.0
      all_ww(i,j,k) = 0.0
      all_w3(i,j,k) = 0.0
      all_uu(i,j,k) = 0.0
      all_vv(i,j,k) = 0.0
      all_uv(i,j,k) = 0.0
      all_dpx(i,j,k) = 0.0
      all_dpy(i,j,k) = 0.0
      all_uth(i,j,k) = 0.0
      all_vth(i,j,k) = 0.0
      all_uthv(i,j,k) = 0.0
      all_vthv(i,j,k) = 0.0
      all_uq(i,j,k) = 0.0
      all_vq(i,j,k) = 0.0
      all_wp(i,j,k) = 0.0
      all_wp_hydro(i,j,k) = 0.0
    END DO
  END DO
  IF (l_qgraup) THEN
    DO j=1,local_new_y
      DO i=1,local_new_x
        all_qgraup(i,j,k) = 0.0
        all_qgraupw(i,j,k) = 0.0
      END DO
    END DO
  END IF

  ! Full local grid arrays
  DO j=1,local_rows
    DO i=1,local_row_len
      bcu_mask(i,j,k) = .FALSE.    ! buoyant mask initialise
      bcu_mask_w(i,j,k) = 0.0      ! buoyant mask initialise
! calculated in read_ff_input
!      ! Calculate density for original grid
!      IF (mask(i,j,k)) THEN
!        density(i,j,k) = p_theta_lev(i,j,k)/(r*t(i,j,k))
!      END IF
    END DO
  END DO

END DO
!$OMP END PARALLEL DO

!-------------------------------------------------------------------------
! Zero arrays
!-------------------------------------------------------------------------
IF (l_acc) THEN
  CALL crmstyle_zero_arrays(acc_factor, acc_h_w, acc_h_u, acc_h_v,        &
  acc_h_th, acc_h_thv, acc_h_rho, acc_h_rh, acc_h_a, acc_h_dpx, acc_h_dpy,&
  acc_h_q,   acc_h_qcl, acc_h_qcf, acc_h_qrain, acc_h_qgraup,             &
  acc_h_dt1, acc_h_dt2, acc_h_dt4, acc_h_dt9, acc_h_dt12, acc_h_dt30,     &
  acc_h_dq4, acc_h_dq9, acc_h_dq12, acc_h_dq30,                           &
  acc_h_dqcl4, acc_h_dqcl9, acc_h_dqcl12, acc_h_dqcl30,                   &
  acc_h_dqcf4, acc_h_dqcf3, acc_h_dqcf12, acc_h_dqcf30,                   &
  acc_h_dqrain30, acc_h_dqgr30, acc_h_drho, acc_h_thw, acc_h_thvw,        &
  acc_h_qw, acc_h_qclw, acc_h_qcfw, acc_h_qrainw, acc_h_qgraupw,          &
  acc_h_uth, acc_h_vth, acc_h_uthv, acc_h_vthv, acc_h_uq, acc_h_vq,       &
  acc_h_wp,                                                               &
  acc_h_uw, acc_h_vw, acc_h_ww, acc_h_w3, acc_h_vv, acc_h_uu, acc_h_uv )
END IF

IF (l_acu) THEN
  CALL crmstyle_zero_arrays(acu_factor, acu_h_w, acu_h_u, acu_h_v,        &
  acu_h_th, acu_h_thv, acu_h_rho, acu_h_rh, acu_h_a, acu_h_dpx, acu_h_dpy,&
  acu_h_q,   acu_h_qcl, acu_h_qcf, acu_h_qrain, acu_h_qgraup,             &
  acu_h_dt1, acu_h_dt2, acu_h_dt4, acu_h_dt9, acu_h_dt12, acu_h_dt30,     &
  acu_h_dq4, acu_h_dq9, acu_h_dq12, acu_h_dq30,                           &
  acu_h_dqcl4, acu_h_dqcl9, acu_h_dqcl12, acu_h_dqcl30,                   &
  acu_h_dqcf4, acu_h_dqcf3, acu_h_dqcf12, acu_h_dqcf30,                   &
  acu_h_dqrain30, acu_h_dqgr30, acu_h_drho, acu_h_thw, acu_h_thvw,        &
  acu_h_qw, acu_h_qclw, acu_h_qcfw, acu_h_qrainw, acu_h_qgraupw,          &
  acu_h_uth, acu_h_vth, acu_h_uthv, acu_h_vthv, acu_h_uq, acu_h_vq,       &
  acu_h_wp,                                                               &
  acu_h_uw, acu_h_vw, acu_h_ww, acu_h_w3, acu_h_vv, acu_h_uu, acu_h_uv )
END IF

IF (l_bcu) THEN
  CALL crmstyle_zero_arrays(bcu_factor, bcu_h_w, bcu_h_u, bcu_h_v,        &
  bcu_h_th, bcu_h_thv, bcu_h_rho, bcu_h_rh, bcu_h_a, bcu_h_dpx, bcu_h_dpy,&
  bcu_h_q,   bcu_h_qcl, bcu_h_qcf, bcu_h_qrain, bcu_h_qgraup,             &
  bcu_h_dt1, bcu_h_dt2, bcu_h_dt4, bcu_h_dt9, bcu_h_dt12, bcu_h_dt30,     &
  bcu_h_dq4, bcu_h_dq9, bcu_h_dq12, bcu_h_dq30,                           &
  bcu_h_dqcl4, bcu_h_dqcl9, bcu_h_dqcl12, bcu_h_dqcl30,                   &
  bcu_h_dqcf4, bcu_h_dqcf3, bcu_h_dqcf12, bcu_h_dqcf30,                   &
  bcu_h_dqrain30, bcu_h_dqgr30, bcu_h_drho, bcu_h_thw, bcu_h_thvw,        &
  bcu_h_qw, bcu_h_qclw, bcu_h_qcfw, bcu_h_qrainw, bcu_h_qgraupw,          &
  bcu_h_uth, bcu_h_vth, bcu_h_uthv, bcu_h_vthv, bcu_h_uq, bcu_h_vq,       &
  bcu_h_wp,                                                               &
  bcu_h_uw, bcu_h_vw, bcu_h_ww, bcu_h_w3, bcu_h_vv, bcu_h_uu, bcu_h_uv )
END IF

IF (l_wg1) THEN
  CALL crmstyle_zero_arrays(wg1_factor, wg1_h_w, wg1_h_u, wg1_h_v,        &
  wg1_h_th, wg1_h_thv, wg1_h_rho, wg1_h_rh, wg1_h_a, wg1_h_dpx, wg1_h_dpy,&
  wg1_h_q,   wg1_h_qcl, wg1_h_qcf, wg1_h_qrain, wg1_h_qgraup,             &
  wg1_h_dt1, wg1_h_dt2, wg1_h_dt4, wg1_h_dt9, wg1_h_dt12, wg1_h_dt30,     &
  wg1_h_dq4, wg1_h_dq9, wg1_h_dq12, wg1_h_dq30,                           &
  wg1_h_dqcl4, wg1_h_dqcl9, wg1_h_dqcl12, wg1_h_dqcl30,                   &
  wg1_h_dqcf4, wg1_h_dqcf3, wg1_h_dqcf12, wg1_h_dqcf30,                   &
  wg1_h_dqrain30, wg1_h_dqgr30, wg1_h_drho, wg1_h_thw, wg1_h_thvw,        &
  wg1_h_qw, wg1_h_qclw, wg1_h_qcfw, wg1_h_qrainw, wg1_h_qgraupw,          &
  wg1_h_uth, wg1_h_vth, wg1_h_uthv, wg1_h_vthv, wg1_h_uq, wg1_h_vq,       &
  wg1_h_wp,                                                               &
  wg1_h_uw, wg1_h_vw, wg1_h_ww, wg1_h_w3, wg1_h_vv, wg1_h_uu, wg1_h_uv )
END IF

IF (l_ppd) THEN
  CALL crmstyle_zero_arrays(ppd_factor, ppd_h_w, ppd_h_u, ppd_h_v,        &
  ppd_h_th, ppd_h_thv, ppd_h_rho, ppd_h_rh, ppd_h_a, ppd_h_dpx, ppd_h_dpy,&
  ppd_h_q,   ppd_h_qcl, ppd_h_qcf, ppd_h_qrain, ppd_h_qgraup,             &
  ppd_h_dt1, ppd_h_dt2, ppd_h_dt4, ppd_h_dt9, ppd_h_dt12, ppd_h_dt30,     &
  ppd_h_dq4, ppd_h_dq9, ppd_h_dq12, ppd_h_dq30,                           &
  ppd_h_dqcl4, ppd_h_dqcl9, ppd_h_dqcl12, ppd_h_dqcl30,                   &
  ppd_h_dqcf4, ppd_h_dqcf3, ppd_h_dqcf12, ppd_h_dqcf30,                   &
  ppd_h_dqrain30, ppd_h_dqgr30, ppd_h_drho, ppd_h_thw, ppd_h_thvw,        &
  ppd_h_qw, ppd_h_qclw, ppd_h_qcfw, ppd_h_qrainw, ppd_h_qgraupw,          &
  ppd_h_uth, ppd_h_vth, ppd_h_uthv, ppd_h_vthv, ppd_h_uq, ppd_h_vq,       &
  ppd_h_wp,                                                               &
  ppd_h_uw, ppd_h_vw, ppd_h_ww, ppd_h_w3, ppd_h_vv, ppd_h_uu, ppd_h_uv )
END IF

IF (l_nbd) THEN
  CALL crmstyle_zero_arrays(nbd_factor, nbd_h_w, nbd_h_u, nbd_h_v,        &
  nbd_h_th, nbd_h_thv, nbd_h_rho, nbd_h_rh, nbd_h_a, nbd_h_dpx, nbd_h_dpy,&
  nbd_h_q,   nbd_h_qcl, nbd_h_qcf, nbd_h_qrain, nbd_h_qgraup,             &
  nbd_h_dt1, nbd_h_dt2, nbd_h_dt4, nbd_h_dt9, nbd_h_dt12, nbd_h_dt30,     &
  nbd_h_dq4, nbd_h_dq9, nbd_h_dq12, nbd_h_dq30,                           &
  nbd_h_dqcl4, nbd_h_dqcl9, nbd_h_dqcl12, nbd_h_dqcl30,                   &
  nbd_h_dqcf4, nbd_h_dqcf3, nbd_h_dqcf12, nbd_h_dqcf30,                   &
  nbd_h_dqrain30, nbd_h_dqgr30, nbd_h_drho, nbd_h_thw, nbd_h_thvw,        &
  nbd_h_qw, nbd_h_qclw, nbd_h_qcfw, nbd_h_qrainw, nbd_h_qgraupw,          &
  nbd_h_uth, nbd_h_vth, nbd_h_uthv, nbd_h_vthv, nbd_h_uq, nbd_h_vq,       &
  nbd_h_wp,                                                               &
  nbd_h_uw, nbd_h_vw, nbd_h_ww, nbd_h_w3, nbd_h_vv, nbd_h_uu, nbd_h_uv )
END IF

IF (l_nid) THEN
  CALL crmstyle_zero_arrays(nid_factor, nid_h_w, nid_h_u, nid_h_v,        &
  nid_h_th, nid_h_thv, nid_h_rho, nid_h_rh, nid_h_a, nid_h_dpx, nid_h_dpy,&
  nid_h_q,   nid_h_qcl, nid_h_qcf, nid_h_qrain, nid_h_qgraup,             &
  nid_h_dt1, nid_h_dt2, nid_h_dt4, nid_h_dt9, nid_h_dt12, nid_h_dt30,     &
  nid_h_dq4, nid_h_dq9, nid_h_dq12, nid_h_dq30,                           &
  nid_h_dqcl4, nid_h_dqcl9, nid_h_dqcl12, nid_h_dqcl30,                   &
  nid_h_dqcf4, nid_h_dqcf3, nid_h_dqcf12, nid_h_dqcf30,                   &
  nid_h_dqrain30, nid_h_dqgr30, nid_h_drho, nid_h_thw, nid_h_thvw,        &
  nid_h_qw, nid_h_qclw, nid_h_qcfw, nid_h_qrainw, nid_h_qgraupw,          &
  nid_h_uth, nid_h_vth, nid_h_uthv, nid_h_vthv, nid_h_uq, nid_h_vq,       &
  nid_h_wp,                                                               &
  nid_h_uw, nid_h_vw, nid_h_ww, nid_h_w3, nid_h_vv, nid_h_uu, nid_h_uv )
END IF

IF (l_adu) THEN
  CALL crmstyle_zero_arrays(adu_factor, adu_h_w, adu_h_u, adu_h_v,        &
  adu_h_th, adu_h_thv, adu_h_rho, adu_h_rh, adu_h_a, adu_h_dpx, adu_h_dpy,&
  adu_h_q,   adu_h_qcl, adu_h_qcf, adu_h_qrain, adu_h_qgraup,             &
  adu_h_dt1, adu_h_dt2, adu_h_dt4, adu_h_dt9, adu_h_dt12, adu_h_dt30,     &
  adu_h_dq4, adu_h_dq9, adu_h_dq12, adu_h_dq30,                           &
  adu_h_dqcl4, adu_h_dqcl9, adu_h_dqcl12, adu_h_dqcl30,                   &
  adu_h_dqcf4, adu_h_dqcf3, adu_h_dqcf12, adu_h_dqcf30,                   &
  adu_h_dqrain30, adu_h_dqgr30, adu_h_drho, adu_h_thw, adu_h_thvw,        &
  adu_h_qw, adu_h_qclw, adu_h_qcfw, adu_h_qrainw, adu_h_qgraupw,          &
  adu_h_uth, adu_h_vth, adu_h_uthv, adu_h_vthv, adu_h_uq, adu_h_vq,       &
  adu_h_wp,                                                               &
  adu_h_uw, adu_h_vw, adu_h_ww, adu_h_w3, adu_h_vv, adu_h_uu, adu_h_uv )
END IF

IF (l_acw) THEN
  CALL crmstyle_zero_arrays(acw_factor, acw_h_w, acw_h_u, acw_h_v,        &
  acw_h_th, acw_h_thv, acw_h_rho, acw_h_rh, acw_h_a, acw_h_dpx, acw_h_dpy,&
  acw_h_q,   acw_h_qcl, acw_h_qcf, acw_h_qrain, acw_h_qgraup,             &
  acw_h_dt1, acw_h_dt2, acw_h_dt4, acw_h_dt9, acw_h_dt12, acw_h_dt30,     &
  acw_h_dq4, acw_h_dq9, acw_h_dq12, acw_h_dq30,                           &
  acw_h_dqcl4, acw_h_dqcl9, acw_h_dqcl12, acw_h_dqcl30,                   &
  acw_h_dqcf4, acw_h_dqcf3, acw_h_dqcf12, acw_h_dqcf30,                   &
  acw_h_dqrain30, acw_h_dqgr30, acw_h_drho, acw_h_thw, acw_h_thvw,        &
  acw_h_qw, acw_h_qclw, acw_h_qcfw, acw_h_qrainw, acw_h_qgraupw,          &
  acw_h_uth, acw_h_vth, acw_h_uthv, acw_h_vthv, acw_h_uq, acw_h_vq,       &
  acw_h_wp,                                                               &
  acw_h_uw, acw_h_vw, acw_h_ww, acw_h_w3, acw_h_vv, acw_h_uu, acw_h_uv )
END IF

IF (l_bcw) THEN
  CALL crmstyle_zero_arrays(bcw_factor, bcw_h_w, bcw_h_u, bcw_h_v,        &
  bcw_h_th, bcw_h_thv, bcw_h_rho, bcw_h_rh, bcw_h_a, bcw_h_dpx, bcw_h_dpy,&
  bcw_h_q,   bcw_h_qcl, bcw_h_qcf, bcw_h_qrain, bcw_h_qgraup,             &
  bcw_h_dt1, bcw_h_dt2, bcw_h_dt4, bcw_h_dt9, bcw_h_dt12, bcw_h_dt30,     &
  bcw_h_dq4, bcw_h_dq9, bcw_h_dq12, bcw_h_dq30,                           &
  bcw_h_dqcl4, bcw_h_dqcl9, bcw_h_dqcl12, bcw_h_dqcl30,                   &
  bcw_h_dqcf4, bcw_h_dqcf3, bcw_h_dqcf12, bcw_h_dqcf30,                   &
  bcw_h_dqrain30, bcw_h_dqgr30, bcw_h_drho, bcw_h_thw, bcw_h_thvw,        &
  bcw_h_qw, bcw_h_qclw, bcw_h_qcfw, bcw_h_qrainw, bcw_h_qgraupw,          &
  bcw_h_uth, bcw_h_vth, bcw_h_uthv, bcw_h_vthv, bcw_h_uq, bcw_h_vq,       &
  bcw_h_wp,                                                               &
  bcw_h_uw, bcw_h_vw, bcw_h_ww, bcw_h_w3, bcw_h_vv, bcw_h_uu, bcw_h_uv )
END IF

IF (l_ucu) THEN
  CALL crmstyle_zero_arrays(ucu_factor, ucu_h_w, ucu_h_u, ucu_h_v,        &
  ucu_h_th, ucu_h_thv, ucu_h_rho, ucu_h_rh, ucu_h_a, ucu_h_dpx, ucu_h_dpy,&
  ucu_h_q,   ucu_h_qcl, ucu_h_qcf, ucu_h_qrain, ucu_h_qgraup,             &
  ucu_h_dt1, ucu_h_dt2, ucu_h_dt4, ucu_h_dt9, ucu_h_dt12, ucu_h_dt30,     &
  ucu_h_dq4, ucu_h_dq9, ucu_h_dq12, ucu_h_dq30,                           &
  ucu_h_dqcl4, ucu_h_dqcl9, ucu_h_dqcl12, ucu_h_dqcl30,                   &
  ucu_h_dqcf4, ucu_h_dqcf3, ucu_h_dqcf12, ucu_h_dqcf30,                   &
  ucu_h_dqrain30, ucu_h_dqgr30, ucu_h_drho, ucu_h_thw, ucu_h_thvw,        &
  ucu_h_qw, ucu_h_qclw, ucu_h_qcfw, ucu_h_qrainw, ucu_h_qgraupw,          &
  ucu_h_uth, ucu_h_vth, ucu_h_uthv, ucu_h_vthv, ucu_h_uq, ucu_h_vq,       &
  ucu_h_wp,                                                               &
  ucu_h_uw, ucu_h_vw, ucu_h_ww, ucu_h_w3, ucu_h_vv, ucu_h_uu, ucu_h_uv )
END IF

IF (l_ppw) THEN
  CALL crmstyle_zero_arrays(ppw_factor, ppw_h_w, ppw_h_u, ppw_h_v,        &
  ppw_h_th, ppw_h_thv, ppw_h_rho, ppw_h_rh, ppw_h_a, ppw_h_dpx, ppw_h_dpy,&
  ppw_h_q,   ppw_h_qcl, ppw_h_qcf, ppw_h_qrain, ppw_h_qgraup,             &
  ppw_h_dt1, ppw_h_dt2, ppw_h_dt4, ppw_h_dt9, ppw_h_dt12, ppw_h_dt30,     &
  ppw_h_dq4, ppw_h_dq9, ppw_h_dq12, ppw_h_dq30,                           &
  ppw_h_dqcl4, ppw_h_dqcl9, ppw_h_dqcl12, ppw_h_dqcl30,                   &
  ppw_h_dqcf4, ppw_h_dqcf3, ppw_h_dqcf12, ppw_h_dqcf30,                   &
  ppw_h_dqrain30, ppw_h_dqgr30, ppw_h_drho, ppw_h_thw, ppw_h_thvw,        &
  ppw_h_qw, ppw_h_qclw, ppw_h_qcfw, ppw_h_qrainw, ppw_h_qgraupw,          &
  ppw_h_uth, ppw_h_vth, ppw_h_uthv, ppw_h_vthv, ppw_h_uq, ppw_h_vq,       &
  ppw_h_wp,                                                               &
  ppw_h_uw, ppw_h_vw, ppw_h_ww, ppw_h_w3, ppw_h_vv, ppw_h_uu, ppw_h_uv )
END IF

IF (l_nbw) THEN
  CALL crmstyle_zero_arrays(nbw_factor, nbw_h_w, nbw_h_u, nbw_h_v,        &
  nbw_h_th, nbw_h_thv, nbw_h_rho, nbw_h_rh, nbw_h_a, nbw_h_dpx, nbw_h_dpy,&
  nbw_h_q,   nbw_h_qcl, nbw_h_qcf, nbw_h_qrain, nbw_h_qgraup,             &
  nbw_h_dt1, nbw_h_dt2, nbw_h_dt4, nbw_h_dt9, nbw_h_dt12, nbw_h_dt30,     &
  nbw_h_dq4, nbw_h_dq9, nbw_h_dq12, nbw_h_dq30,                           &
  nbw_h_dqcl4, nbw_h_dqcl9, nbw_h_dqcl12, nbw_h_dqcl30,                   &
  nbw_h_dqcf4, nbw_h_dqcf3, nbw_h_dqcf12, nbw_h_dqcf30,                   &
  nbw_h_dqrain30, nbw_h_dqgr30, nbw_h_drho, nbw_h_thw, nbw_h_thvw,        &
  nbw_h_qw, nbw_h_qclw, nbw_h_qcfw, nbw_h_qrainw, nbw_h_qgraupw,          &
  nbw_h_uth, nbw_h_vth, nbw_h_uthv, nbw_h_vthv, nbw_h_uq, nbw_h_vq,       &
  nbw_h_wp,                                                               &
  nbw_h_uw, nbw_h_vw, nbw_h_ww, nbw_h_w3, nbw_h_vv, nbw_h_uu, nbw_h_uv )
END IF
!-------------------------------------------------------------------------

! Calculation of dthes/dz 
IF (l_ucu) THEN

  DO k=1,mlevs-1 
    rdz = 1./(h_theta_sea(k+1) - h_theta_sea(k))
    ! Full local grid arrays
    DO j=1,local_rows
      DO i=1,local_row_len
        IF (mask(i,j,k)) THEN
          ! Calculate qsat/T using RH
          qs_over_t(i,j,k) =  100.0 *q(i,j,k)/(rh(i,j,k)*t(i,j,k))
          dthdz = (theta(i,j,k+1) -theta(i,j,k))*rdz/theta(i,j,k)
          IF (t(i,j,k) < tm) THEN
            lrcp = lsrcp
          ELSE 
            lrcp = lcrcp 
          END IF
          dqstdz = lrcp*( qs_over_t(i,j,k+1) - qs_over_t(i,j,k) )*rdz 
          dthesdz(i,j,k) = -dthdz -dqstdz
        END IF    ! test on mask
      END DO
    END DO
  END DO
  ! Top level set to zero  
  DO j=1,local_rows
    DO i=1,local_row_len
      dthesdz(i,j,k) = 0.0
    END DO
  END DO
END IF

! -------------------------------------------------------------------
!  Main work creating selective sampling
! -------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(k, i, j, ii, jj, precip_water,                   &
!$OMP& qw_work, qclw_work, qcfw_work, qrainw_work, qgraupw_work, thw_work, &
!$OMP& thvw_work, ww_work, w3_work, uw_work, vw_work, uu_work,             &
!$OMP& vv_work, uv_work, dpx_work, dpy_work, qw_test, uq_work, vq_work,    &
!$OMP& uth_work, vth_work, uthv_work, vthv_work, wp_work,                  &
!$OMP& wp_work2) DEFAULT(SHARED)

DO k=1,mlevs

  ! Loop over fine grid but need to assign to a coarse grid ii, jj

  DO j=1,local_rows

    jj = index_row(j)

    DO i=1,local_row_len

      ii = index_col(i)

      IF (mask(i,j,k)) THEN    ! above surface of model

        all_factor(ii,jj,k) = all_factor(ii,jj,k)+1.0

        qw_work     = q_h_prime(i,j,k)    *w_h_prime(i,j,k)
        qclw_work   = qcl_h_prime(i,j,k)  *w_h_prime(i,j,k)
        qcfw_work   = qcf_h_prime(i,j,k)  *w_h_prime(i,j,k)
        qrainw_work = qrain_h_prime(i,j,k)*w_h_prime(i,j,k)
        IF (l_qgraup) THEN
          qgraupw_work = qgraup_h_prime(i,j,k)*w_h_prime(i,j,k)
        END IF
        thw_work    = th_h_prime(i,j,k)*w_h_prime(i,j,k)
        thvw_work   = thv_h_prime(i,j,k)*w_h_prime(i,j,k)
        ww_work     = w_h_prime(i,j,k)*w_h_prime(i,j,k)
        w3_work     = ww_work*w_h_prime(i,j,k)
        uw_work     = u_h_prime(i,j,k)*w_h_prime(i,j,k)
        vw_work     = v_h_prime(i,j,k)*w_h_prime(i,j,k)
        uu_work     = u_h_prime(i,j,k)*u_h_prime(i,j,k)
        vv_work     = v_h_prime(i,j,k)*v_h_prime(i,j,k)
        uv_work     = u_h_prime(i,j,k)*v_h_prime(i,j,k)
        dpx_work    = dpdx(i,j,k)/density(i,j,k)
        dpy_work    = dpdy(i,j,k)/density(i,j,k)
        uq_work     = q_h_prime(i,j,k) *u_h_prime(i,j,k)
        vq_work     = q_h_prime(i,j,k) *v_h_prime(i,j,k)
        uth_work    = th_h_prime(i,j,k) *u_h_prime(i,j,k)
        uthv_work   = thv_h_prime(i,j,k)*u_h_prime(i,j,k)
        vth_work    = th_h_prime(i,j,k) *v_h_prime(i,j,k)
        vthv_work   = thv_h_prime(i,j,k)*v_h_prime(i,j,k)
        wp_work2    = p_h_prime(i,j,k)*w_h_prime(i,j,k)/density(i,j,k)
        wp_work     = (p_theta_lev(i,j,k)- p_theta_hydro(ii,jj,k)) &
                              *w_h_prime(i,j,k)/density(i,j,k)

        ! products for a level

        all_thw(ii,jj,k)   = all_thw(ii,jj,k)   + thw_work
        all_thvw(ii,jj,k)  = all_thvw(ii,jj,k)  + thvw_work
        all_qw(ii,jj,k)    = all_qw(ii,jj,k)    + qw_work
        all_qclw(ii,jj,k)  = all_qclw(ii,jj,k)  + qclw_work
        all_qcfw(ii,jj,k)  = all_qcfw(ii,jj,k)  + qcfw_work
        all_qrainw(ii,jj,k)= all_qrainw(ii,jj,k)+ qrainw_work
        all_ww(ii,jj,k)    = all_ww(ii,jj,k)    + ww_work
        all_w3(ii,jj,k)    = all_w3(ii,jj,k)    + w3_work
        all_uw(ii,jj,k)    = all_uw(ii,jj,k)    + uw_work
        all_vw(ii,jj,k)    = all_vw(ii,jj,k)    + vw_work
        all_uu(ii,jj,k)    = all_uu(ii,jj,k)    + uu_work
        all_vv(ii,jj,k)    = all_vv(ii,jj,k)    + vv_work
        all_uv(ii,jj,k)    = all_uv(ii,jj,k)    + uv_work
        all_dpx(ii,jj,k) = all_dpx(ii,jj,k) + dpx_work
        all_dpy(ii,jj,k) = all_dpy(ii,jj,k) + dpy_work
        all_rho(ii,jj,k) = all_rho(ii,jj,k) + density(i,j,k)
        all_uth(ii,jj,k)   = all_uth(ii,jj,k)   + uth_work
        all_uthv(ii,jj,k)  = all_uthv(ii,jj,k)  + uthv_work
        all_vth(ii,jj,k)   = all_vth(ii,jj,k)   + vth_work
        all_vthv(ii,jj,k)  = all_vthv(ii,jj,k)  + vthv_work
        all_uq(ii,jj,k)    = all_uq(ii,jj,k)   + uq_work
        all_vq(ii,jj,k)    = all_vq(ii,jj,k)   + vq_work
        all_wp(ii,jj,k)    = all_wp(ii,jj,k)   + wp_work2
        all_wp_hydro(ii,jj,k) = all_wp_hydro(ii,jj,k) + wp_work

        IF (l_qgraup) THEN
          all_qgraupw(ii,jj,k)= all_qgraupw(ii,jj,k)+ qgraupw_work
        END IF

        ! No cloud condensate or precip present but upwards
        IF (l_qgraup) THEN
          qw_test = qcf(i,j,k)+qcl(i,j,k)+qrain(i,j,k)+qgraup(i,j,k)
        ELSE
          qw_test = qcf(i,j,k)+qcl(i,j,k)+qrain(i,j,k)
        END IF

        !------------------------------------------------------------
        ! Partitions
        !------------------------------------------------------------
        IF (l_adu .AND. qw_test < qw_crit .AND. (w(i,j,k) > 0.0) ) THEN
          adu_factor(ii,jj,k) = adu_factor(ii,jj,k)+1.0
          adu_h_rho(ii,jj,k) = adu_h_rho(ii,jj,k) + density(i,j,k)
          adu_h_w(ii,jj,k) = adu_h_w(ii,jj,k)+ w(i,j,k)
          adu_h_u(ii,jj,k) = adu_h_u(ii,jj,k)+ u(i,j,k)
          adu_h_v(ii,jj,k) = adu_h_v(ii,jj,k)+ v(i,j,k)
          adu_h_th(ii,jj,k)  = adu_h_th(ii,jj,k)  +th_h_prime(i,j,k)
          adu_h_thv(ii,jj,k) = adu_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
          adu_h_rh(ii,jj,k)    = adu_h_rh(ii,jj,k) +rh(i,j,k)
          adu_h_q(ii,jj,k)     = adu_h_q(ii,jj,k)   +q(i,j,k)
          adu_h_qcl(ii,jj,k)   = adu_h_qcl(ii,jj,k) +qcl(i,j,k)
          adu_h_qcf(ii,jj,k)   = adu_h_qcf(ii,jj,k) +qcf(i,j,k)
          adu_h_qrain(ii,jj,k) = adu_h_qrain(ii,jj,k) +qrain(i,j,k)
          adu_h_dt1(ii,jj,k) = adu_h_dt1(ii,jj,k) +dt1(i,j,k)
          adu_h_dt2(ii,jj,k) = adu_h_dt2(ii,jj,k) +dt2(i,j,k)
          adu_h_dt4(ii,jj,k) = adu_h_dt4(ii,jj,k) +dt4(i,j,k)
          adu_h_dt9(ii,jj,k) = adu_h_dt9(ii,jj,k) +dt9(i,j,k)
          adu_h_dt12(ii,jj,k) = adu_h_dt12(ii,jj,k) +dt12(i,j,k)
          adu_h_dq4(ii,jj,k) = adu_h_dq4(ii,jj,k) +dq4(i,j,k)
          adu_h_dq9(ii,jj,k) = adu_h_dq9(ii,jj,k) +dq9(i,j,k)
          adu_h_dq12(ii,jj,k) = adu_h_dq12(ii,jj,k) +dq12(i,j,k)
          adu_h_dqcl4(ii,jj,k) = adu_h_dqcl4(ii,jj,k)   +dqcl4(i,j,k)
          adu_h_dqcl9(ii,jj,k) = adu_h_dqcl9(ii,jj,k)   +dqcl9(i,j,k)
          adu_h_dqcl12(ii,jj,k) = adu_h_dqcl12(ii,jj,k) +dqcl12(i,j,k)
          adu_h_dqcf4(ii,jj,k) = adu_h_dqcf4(ii,jj,k)   +dqcf4(i,j,k)
          adu_h_dqcf3(ii,jj,k) = adu_h_dqcf3(ii,jj,k)   +dqcf3(i,j,k)
          adu_h_dqcf12(ii,jj,k) = adu_h_dqcf12(ii,jj,k) +dqcf12(i,j,k)
          adu_h_thw(ii,jj,k)   = adu_h_thw(ii,jj,k)   + thw_work
          adu_h_thvw(ii,jj,k)  = adu_h_thvw(ii,jj,k)  + thvw_work
          adu_h_qw(ii,jj,k)    = adu_h_qw(ii,jj,k)    + qw_work
          adu_h_qclw(ii,jj,k)  = adu_h_qclw(ii,jj,k)  + qclw_work
          adu_h_qcfw(ii,jj,k)  = adu_h_qcfw(ii,jj,k)  + qcfw_work
          adu_h_qrainw(ii,jj,k)= adu_h_qrainw(ii,jj,k)+ qrainw_work
          adu_h_ww(ii,jj,k)    = adu_h_ww(ii,jj,k)    + ww_work
          adu_h_w3(ii,jj,k)    = adu_h_w3(ii,jj,k)    + w3_work
          adu_h_uw(ii,jj,k)    = adu_h_uw(ii,jj,k)    + uw_work
          adu_h_vw(ii,jj,k)    = adu_h_vw(ii,jj,k)    + vw_work
          adu_h_uu(ii,jj,k)    = adu_h_uu(ii,jj,k)    + uu_work
          adu_h_vv(ii,jj,k)    = adu_h_vv(ii,jj,k)    + vv_work
          adu_h_uv(ii,jj,k)    = adu_h_uv(ii,jj,k)    + uv_work
          adu_h_dpx(ii,jj,k)   = adu_h_dpx(ii,jj,k)   + dpx_work
          adu_h_dpy(ii,jj,k)   = adu_h_dpy(ii,jj,k)   + dpy_work
          adu_h_uth(ii,jj,k)   = adu_h_uth(ii,jj,k)   + uth_work
          adu_h_uthv(ii,jj,k)  = adu_h_uthv(ii,jj,k)  + uthv_work
          adu_h_vth(ii,jj,k)   = adu_h_vth(ii,jj,k)   + vth_work
          adu_h_vthv(ii,jj,k)  = adu_h_vthv(ii,jj,k)  + vthv_work
          adu_h_uq(ii,jj,k)    = adu_h_uq(ii,jj,k)   + uq_work
          adu_h_vq(ii,jj,k)    = adu_h_vq(ii,jj,k)   + vq_work
          adu_h_wp(ii,jj,k)    = adu_h_wp(ii,jj,k)   + wp_work
          IF (l_qgraup) THEN
            adu_h_qgraup(ii,jj,k)= adu_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
            adu_h_qgraupw(ii,jj,k)= adu_h_qgraupw(ii,jj,k)+ qgraupw_work
          END IF
          IF (l_sect30) THEN
            adu_h_dt30(ii,jj,k) = adu_h_dt30(ii,jj,k) +dt30(i,j,k)
            adu_h_dq30(ii,jj,k) = adu_h_dq30(ii,jj,k) +dq30(i,j,k)
            adu_h_dqcl30(ii,jj,k) = adu_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
            adu_h_dqcf30(ii,jj,k) = adu_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
            adu_h_dqrain30(ii,jj,k) = adu_h_dqrain30(ii,jj,k) +dqrain30(i,j,k)
            adu_h_drho(ii,jj,k) = adu_h_drho(ii,jj,k) +drho(i,j,k)
          END IF
          IF (l_sect30 .AND. l_qgraup ) THEN
            adu_h_dqgr30(ii,jj,k) = adu_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
          END IF
        END IF
        !------------------------------------------------------------
        ! Upward, Cloudy, and conditionally Unstable
        !------------------------------------------------------------
        IF (l_ucu .AND. qw_test > qw_crit .AND. (w(i,j,k) > 0.0)   &
              .AND. dthesdz(i,j,k) > 0.0 ) THEN
          ucu_factor(ii,jj,k) = ucu_factor(ii,jj,k)+1.0
          ucu_h_rho(ii,jj,k) = ucu_h_rho(ii,jj,k) + density(i,j,k)
          ucu_h_w(ii,jj,k) = ucu_h_w(ii,jj,k)+ w(i,j,k)
          ucu_h_u(ii,jj,k) = ucu_h_u(ii,jj,k)+ u(i,j,k)
          ucu_h_v(ii,jj,k) = ucu_h_v(ii,jj,k)+ v(i,j,k)
          ucu_h_th(ii,jj,k)  = ucu_h_th(ii,jj,k)  +th_h_prime(i,j,k)
          ucu_h_thv(ii,jj,k) = ucu_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
          ucu_h_rh(ii,jj,k)    = ucu_h_rh(ii,jj,k) +rh(i,j,k)
          ucu_h_q(ii,jj,k)     = ucu_h_q(ii,jj,k)   +q(i,j,k)
          ucu_h_qcl(ii,jj,k)   = ucu_h_qcl(ii,jj,k) +qcl(i,j,k)
          ucu_h_qcf(ii,jj,k)   = ucu_h_qcf(ii,jj,k) +qcf(i,j,k)
          ucu_h_qrain(ii,jj,k) = ucu_h_qrain(ii,jj,k) +qrain(i,j,k)
          ucu_h_dt1(ii,jj,k) = ucu_h_dt1(ii,jj,k) +dt1(i,j,k)
          ucu_h_dt2(ii,jj,k) = ucu_h_dt2(ii,jj,k) +dt2(i,j,k)
          ucu_h_dt4(ii,jj,k) = ucu_h_dt4(ii,jj,k) +dt4(i,j,k)
          ucu_h_dt9(ii,jj,k) = ucu_h_dt9(ii,jj,k) +dt9(i,j,k)
          ucu_h_dt12(ii,jj,k) = ucu_h_dt12(ii,jj,k) +dt12(i,j,k)
          ucu_h_dq4(ii,jj,k) = ucu_h_dq4(ii,jj,k) +dq4(i,j,k)
          ucu_h_dq9(ii,jj,k) = ucu_h_dq9(ii,jj,k) +dq9(i,j,k)
          ucu_h_dq12(ii,jj,k) = ucu_h_dq12(ii,jj,k) +dq12(i,j,k)
          ucu_h_dqcl4(ii,jj,k) = ucu_h_dqcl4(ii,jj,k)   +dqcl4(i,j,k)
          ucu_h_dqcl9(ii,jj,k) = ucu_h_dqcl9(ii,jj,k)   +dqcl9(i,j,k)
          ucu_h_dqcl12(ii,jj,k) = ucu_h_dqcl12(ii,jj,k) +dqcl12(i,j,k)
          ucu_h_dqcf4(ii,jj,k) = ucu_h_dqcf4(ii,jj,k)   +dqcf4(i,j,k)
          ucu_h_dqcf3(ii,jj,k) = ucu_h_dqcf3(ii,jj,k)   +dqcf3(i,j,k)
          ucu_h_dqcf12(ii,jj,k) = ucu_h_dqcf12(ii,jj,k) +dqcf12(i,j,k)
          ucu_h_thw(ii,jj,k)   = ucu_h_thw(ii,jj,k)   + thw_work
          ucu_h_thvw(ii,jj,k)  = ucu_h_thvw(ii,jj,k)  + thvw_work
          ucu_h_qw(ii,jj,k)    = ucu_h_qw(ii,jj,k)    + qw_work
          ucu_h_qclw(ii,jj,k)  = ucu_h_qclw(ii,jj,k)  + qclw_work
          ucu_h_qcfw(ii,jj,k)  = ucu_h_qcfw(ii,jj,k)  + qcfw_work
          ucu_h_qrainw(ii,jj,k)= ucu_h_qrainw(ii,jj,k)+ qrainw_work
          ucu_h_ww(ii,jj,k)    = ucu_h_ww(ii,jj,k)    + ww_work
          ucu_h_w3(ii,jj,k)    = ucu_h_w3(ii,jj,k)    + w3_work
          ucu_h_uw(ii,jj,k)    = ucu_h_uw(ii,jj,k)    + uw_work
          ucu_h_vw(ii,jj,k)    = ucu_h_vw(ii,jj,k)    + vw_work
          ucu_h_uu(ii,jj,k)    = ucu_h_uu(ii,jj,k)    + uu_work
          ucu_h_vv(ii,jj,k)    = ucu_h_vv(ii,jj,k)    + vv_work
          ucu_h_uv(ii,jj,k)    = ucu_h_uv(ii,jj,k)    + uv_work
          ucu_h_dpx(ii,jj,k)   = ucu_h_dpx(ii,jj,k)   + dpx_work
          ucu_h_dpy(ii,jj,k)   = ucu_h_dpy(ii,jj,k)   + dpy_work
          ucu_h_uth(ii,jj,k)   = ucu_h_uth(ii,jj,k)   + uth_work
          ucu_h_uthv(ii,jj,k)  = ucu_h_uthv(ii,jj,k)  + uthv_work
          ucu_h_vth(ii,jj,k)   = ucu_h_vth(ii,jj,k)   + vth_work
          ucu_h_vthv(ii,jj,k)  = ucu_h_vthv(ii,jj,k)  + vthv_work
          ucu_h_uq(ii,jj,k)    = ucu_h_uq(ii,jj,k)   + uq_work
          ucu_h_vq(ii,jj,k)    = ucu_h_vq(ii,jj,k)   + vq_work
          ucu_h_wp(ii,jj,k)    = ucu_h_wp(ii,jj,k)   + wp_work
          IF (l_qgraup) THEN
            ucu_h_qgraup(ii,jj,k)= ucu_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
            ucu_h_qgraupw(ii,jj,k)= ucu_h_qgraupw(ii,jj,k)+ qgraupw_work
          END IF
          IF (l_sect30) THEN
            ucu_h_dt30(ii,jj,k) = ucu_h_dt30(ii,jj,k) +dt30(i,j,k)
            ucu_h_dq30(ii,jj,k) = ucu_h_dq30(ii,jj,k) +dq30(i,j,k)
            ucu_h_dqcl30(ii,jj,k) = ucu_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
            ucu_h_dqcf30(ii,jj,k) = ucu_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
            ucu_h_dqrain30(ii,jj,k) = ucu_h_dqrain30(ii,jj,k) +dqrain30(i,j,k)
            ucu_h_drho(ii,jj,k) = ucu_h_drho(ii,jj,k) +drho(i,j,k)
          END IF
          IF (l_sect30 .AND. l_qgraup ) THEN
            ucu_h_dqgr30(ii,jj,k) = ucu_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
          END IF
        END IF


        ! Cloud condensate present

        IF (qcf(i,j,k)+qcl(i,j,k) > qw_crit) THEN

          IF (l_acc) THEN  
            acc_factor(ii,jj,k) = acc_factor(ii,jj,k)+1.0
            acc_h_rho(ii,jj,k) = acc_h_rho(ii,jj,k) + density(i,j,k)
            acc_h_w(ii,jj,k) = acc_h_w(ii,jj,k)+ w(i,j,k)
            acc_h_u(ii,jj,k) = acc_h_u(ii,jj,k)+ u(i,j,k)
            acc_h_v(ii,jj,k) = acc_h_v(ii,jj,k)+ v(i,j,k)
            acc_h_th(ii,jj,k)  = acc_h_th(ii,jj,k)  +th_h_prime(i,j,k)
            acc_h_thv(ii,jj,k) = acc_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
            acc_h_rh(ii,jj,k)    = acc_h_rh(ii,jj,k) +rh(i,j,k)
            acc_h_q(ii,jj,k)     = acc_h_q(ii,jj,k)   +q(i,j,k)
            acc_h_qcl(ii,jj,k)   = acc_h_qcl(ii,jj,k) +qcl(i,j,k)
            acc_h_qcf(ii,jj,k)   = acc_h_qcf(ii,jj,k) +qcf(i,j,k)
            acc_h_qrain(ii,jj,k) = acc_h_qrain(ii,jj,k) +qrain(i,j,k)
            acc_h_dt1(ii,jj,k) = acc_h_dt1(ii,jj,k) +dt1(i,j,k)
            acc_h_dt2(ii,jj,k) = acc_h_dt2(ii,jj,k) +dt2(i,j,k)
            acc_h_dt4(ii,jj,k) = acc_h_dt4(ii,jj,k) +dt4(i,j,k)
            acc_h_dt9(ii,jj,k) = acc_h_dt9(ii,jj,k) +dt9(i,j,k)
            acc_h_dt12(ii,jj,k) = acc_h_dt12(ii,jj,k) +dt12(i,j,k)
            acc_h_dq4(ii,jj,k) = acc_h_dq4(ii,jj,k) +dq4(i,j,k)
            acc_h_dq9(ii,jj,k) = acc_h_dq9(ii,jj,k) +dq9(i,j,k)
            acc_h_dq12(ii,jj,k) = acc_h_dq12(ii,jj,k) +dq12(i,j,k)
            acc_h_dqcl4(ii,jj,k) = acc_h_dqcl4(ii,jj,k)   +dqcl4(i,j,k)
            acc_h_dqcl9(ii,jj,k) = acc_h_dqcl9(ii,jj,k)   +dqcl9(i,j,k)
            acc_h_dqcl12(ii,jj,k) = acc_h_dqcl12(ii,jj,k) +dqcl12(i,j,k)
            acc_h_dqcf4(ii,jj,k) = acc_h_dqcf4(ii,jj,k)   +dqcf4(i,j,k)
            acc_h_dqcf3(ii,jj,k) = acc_h_dqcf3(ii,jj,k)   +dqcf3(i,j,k)
            acc_h_dqcf12(ii,jj,k) = acc_h_dqcf12(ii,jj,k) +dqcf12(i,j,k)
            acc_h_thw(ii,jj,k)   = acc_h_thw(ii,jj,k)   + thw_work
            acc_h_thvw(ii,jj,k)  = acc_h_thvw(ii,jj,k)  + thvw_work
            acc_h_qw(ii,jj,k)    = acc_h_qw(ii,jj,k)    + qw_work
            acc_h_qclw(ii,jj,k)  = acc_h_qclw(ii,jj,k)  + qclw_work
            acc_h_qcfw(ii,jj,k)  = acc_h_qcfw(ii,jj,k)  + qcfw_work
            acc_h_qrainw(ii,jj,k)= acc_h_qrainw(ii,jj,k)+ qrainw_work
            acc_h_ww(ii,jj,k)    = acc_h_ww(ii,jj,k)    + ww_work
            acc_h_w3(ii,jj,k)    = acc_h_w3(ii,jj,k)    + w3_work
            acc_h_uw(ii,jj,k)    = acc_h_uw(ii,jj,k)    + uw_work
            acc_h_vw(ii,jj,k)    = acc_h_vw(ii,jj,k)    + vw_work
            acc_h_uu(ii,jj,k)    = acc_h_uu(ii,jj,k)    + uu_work
            acc_h_vv(ii,jj,k)    = acc_h_vv(ii,jj,k)    + vv_work
            acc_h_uv(ii,jj,k)    = acc_h_uv(ii,jj,k)    + uv_work
            acc_h_dpx(ii,jj,k)   = acc_h_dpx(ii,jj,k)   + dpx_work
            acc_h_dpy(ii,jj,k)   = acc_h_dpy(ii,jj,k)   + dpy_work
            acc_h_uth(ii,jj,k)   = acc_h_uth(ii,jj,k)   + uth_work
            acc_h_uthv(ii,jj,k)  = acc_h_uthv(ii,jj,k)  + uthv_work
            acc_h_vth(ii,jj,k)   = acc_h_vth(ii,jj,k)   + vth_work
            acc_h_vthv(ii,jj,k)  = acc_h_vthv(ii,jj,k)  + vthv_work
            acc_h_uq(ii,jj,k)    = acc_h_uq(ii,jj,k)   + uq_work
            acc_h_vq(ii,jj,k)    = acc_h_vq(ii,jj,k)   + vq_work
            acc_h_wp(ii,jj,k)    = acc_h_wp(ii,jj,k)   + wp_work
            IF (l_qgraup) THEN
              acc_h_qgraup(ii,jj,k)= acc_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
              acc_h_qgraupw(ii,jj,k)= acc_h_qgraupw(ii,jj,k)+ qgraupw_work
            END IF
            IF (l_sect30) THEN
              acc_h_dt30(ii,jj,k) = acc_h_dt30(ii,jj,k) +dt30(i,j,k)
              acc_h_dq30(ii,jj,k) = acc_h_dq30(ii,jj,k) +dq30(i,j,k)
              acc_h_dqcl30(ii,jj,k) = acc_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
              acc_h_dqcf30(ii,jj,k) = acc_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
              acc_h_dqrain30(ii,jj,k) = acc_h_dqrain30(ii,jj,k) +dqrain30(i,j,k)
              acc_h_drho(ii,jj,k) = acc_h_drho(ii,jj,k) +drho(i,j,k)
            END IF
            IF (l_sect30 .AND. l_qgraup ) THEN
              acc_h_dqgr30(ii,jj,k) = acc_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
            END IF
          END IF   ! acc

          ! updraughts relative to mean vertical velocity

          IF (w_h_prime(i,j,k) > 0.0) THEN
            IF (l_acu) THEN
              acu_factor(ii,jj,k) = acu_factor(ii,jj,k)+1.0
              acu_h_rho(ii,jj,k) = acu_h_rho(ii,jj,k) + density(i,j,k)
              acu_h_w(ii,jj,k) = acu_h_w(ii,jj,k)+ w(i,j,k)
              acu_h_u(ii,jj,k) = acu_h_u(ii,jj,k)+ u(i,j,k)
              acu_h_v(ii,jj,k) = acu_h_v(ii,jj,k)+ v(i,j,k)
              acu_h_th(ii,jj,k) = acu_h_th(ii,jj,k) +th_h_prime(i,j,k)
              acu_h_thv(ii,jj,k) = acu_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
              acu_h_rh(ii,jj,k) = acu_h_rh(ii,jj,k) +rh(i,j,k)
              acu_h_q(ii,jj,k) = acu_h_q(ii,jj,k) +q(i,j,k)
              acu_h_qcl(ii,jj,k) = acu_h_qcl(ii,jj,k) +qcl(i,j,k)
              acu_h_qcf(ii,jj,k) = acu_h_qcf(ii,jj,k) +qcf(i,j,k)
              acu_h_qrain(ii,jj,k) = acu_h_qrain(ii,jj,k) +qrain(i,j,k)
              acu_h_dt1(ii,jj,k) = acu_h_dt1(ii,jj,k) +dt1(i,j,k)
              acu_h_dt2(ii,jj,k) = acu_h_dt2(ii,jj,k) +dt2(i,j,k)
              acu_h_dt4(ii,jj,k) = acu_h_dt4(ii,jj,k) +dt4(i,j,k)
              acu_h_dt9(ii,jj,k) = acu_h_dt9(ii,jj,k) +dt9(i,j,k)
              acu_h_dt12(ii,jj,k) = acu_h_dt12(ii,jj,k) +dt12(i,j,k)
              acu_h_dq4(ii,jj,k) = acu_h_dq4(ii,jj,k) +dq4(i,j,k)
              acu_h_dq9(ii,jj,k) = acu_h_dq9(ii,jj,k) +dq9(i,j,k)
              acu_h_dq12(ii,jj,k) = acu_h_dq12(ii,jj,k) +dq12(i,j,k)
              acu_h_dqcl4(ii,jj,k) = acu_h_dqcl4(ii,jj,k)   +dqcl4(i,j,k)
              acu_h_dqcl9(ii,jj,k) = acu_h_dqcl9(ii,jj,k)   +dqcl9(i,j,k)
              acu_h_dqcl12(ii,jj,k) = acu_h_dqcl12(ii,jj,k) +dqcl12(i,j,k)
              acu_h_dqcf4(ii,jj,k) = acu_h_dqcf4(ii,jj,k)   +dqcf4(i,j,k)
              acu_h_dqcf3(ii,jj,k) = acu_h_dqcf3(ii,jj,k)   +dqcf3(i,j,k)
              acu_h_dqcf12(ii,jj,k) = acu_h_dqcf12(ii,jj,k) +dqcf12(i,j,k)
              acu_h_thw(ii,jj,k)   = acu_h_thw(ii,jj,k)   + thw_work
              acu_h_thvw(ii,jj,k)  = acu_h_thvw(ii,jj,k)  + thvw_work
              acu_h_qw(ii,jj,k)    = acu_h_qw(ii,jj,k)    + qw_work
              acu_h_qclw(ii,jj,k)  = acu_h_qclw(ii,jj,k)  + qclw_work
              acu_h_qcfw(ii,jj,k)  = acu_h_qcfw(ii,jj,k)  + qcfw_work
              acu_h_qrainw(ii,jj,k)= acu_h_qrainw(ii,jj,k)+ qrainw_work
              acu_h_ww(ii,jj,k)    = acu_h_ww(ii,jj,k)    + ww_work
              acu_h_w3(ii,jj,k)    = acu_h_w3(ii,jj,k)    + w3_work
              acu_h_uw(ii,jj,k)    = acu_h_uw(ii,jj,k)    + uw_work
              acu_h_vw(ii,jj,k)    = acu_h_vw(ii,jj,k)    + vw_work
              acu_h_uu(ii,jj,k)    = acu_h_uu(ii,jj,k)    + uu_work
              acu_h_vv(ii,jj,k)    = acu_h_vv(ii,jj,k)    + vv_work
              acu_h_uv(ii,jj,k)    = acu_h_uv(ii,jj,k)    + uv_work
              acu_h_dpx(ii,jj,k)   = acu_h_dpx(ii,jj,k)   + dpx_work
              acu_h_dpy(ii,jj,k)   = acu_h_dpy(ii,jj,k)   + dpy_work
              acu_h_uth(ii,jj,k)   = acu_h_uth(ii,jj,k)   + uth_work
              acu_h_uthv(ii,jj,k)  = acu_h_uthv(ii,jj,k)  + uthv_work
              acu_h_vth(ii,jj,k)   = acu_h_vth(ii,jj,k)   + vth_work
              acu_h_vthv(ii,jj,k)  = acu_h_vthv(ii,jj,k)  + vthv_work
              acu_h_uq(ii,jj,k)    = acu_h_uq(ii,jj,k)   + uq_work
              acu_h_vq(ii,jj,k)    = acu_h_vq(ii,jj,k)   + vq_work
              acu_h_wp(ii,jj,k)    = acu_h_wp(ii,jj,k)   + wp_work
              IF (l_qgraup) THEN
                acu_h_qgraup(ii,jj,k)= acu_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
                acu_h_qgraupw(ii,jj,k)= acu_h_qgraupw(ii,jj,k)+ qgraupw_work
              END IF
              IF (l_sect30) THEN
                acu_h_dt30(ii,jj,k) = acu_h_dt30(ii,jj,k) +dt30(i,j,k)
                acu_h_dq30(ii,jj,k) = acu_h_dq30(ii,jj,k) +dq30(i,j,k)
                acu_h_dqcl30(ii,jj,k) = acu_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
                acu_h_dqcf30(ii,jj,k) = acu_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
                acu_h_dqrain30(ii,jj,k) = acu_h_dqrain30(ii,jj,k)  &
                                                            +dqrain30(i,j,k)
                acu_h_drho(ii,jj,k) = acu_h_drho(ii,jj,k) +drho(i,j,k)
              END IF
              IF (l_sect30 .AND. l_qgraup ) THEN
                acu_h_dqgr30(ii,jj,k) = acu_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
              END IF
            END IF  ! acu

            ! Buoyant plumes
            IF (l_bcu .AND. thv_h_prime(i,j,k) > 0.0) THEN ! BCU

              bcu_mask(i,j,k) = .TRUE.
              bcu_mask_w(i,j,k) = w(i,j,k)
              nn_mask(ii,jj,k) = nn_mask(ii,jj,k) + 1

              bcu_factor(ii,jj,k) = bcu_factor(ii,jj,k)+1.0
              bcu_h_rho(ii,jj,k) = bcu_h_rho(ii,jj,k) + density(i,j,k)
              bcu_h_w(ii,jj,k) = bcu_h_w(ii,jj,k)+ w(i,j,k)
              bcu_h_u(ii,jj,k) = bcu_h_u(ii,jj,k)+ u(i,j,k)
              bcu_h_v(ii,jj,k) = bcu_h_v(ii,jj,k)+ v(i,j,k)
              bcu_h_th(ii,jj,k)  = bcu_h_th(ii,jj,k) +th_h_prime(i,j,k)
              bcu_h_thv(ii,jj,k) = bcu_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
              bcu_h_rh(ii,jj,k) = bcu_h_rh(ii,jj,k) +rh(i,j,k)
              bcu_h_q(ii,jj,k)   = bcu_h_q(ii,jj,k) +q(i,j,k)
              bcu_h_qcl(ii,jj,k) = bcu_h_qcl(ii,jj,k) +qcl(i,j,k)
              bcu_h_qcf(ii,jj,k) = bcu_h_qcf(ii,jj,k) +qcf(i,j,k)
              bcu_h_qrain(ii,jj,k) = bcu_h_qrain(ii,jj,k) +qrain(i,j,k)
              bcu_h_dt1(ii,jj,k) = bcu_h_dt1(ii,jj,k) +dt1(i,j,k)
              bcu_h_dt2(ii,jj,k) = bcu_h_dt2(ii,jj,k) +dt2(i,j,k)
              bcu_h_dt4(ii,jj,k) = bcu_h_dt4(ii,jj,k) +dt4(i,j,k)
              bcu_h_dt9(ii,jj,k) = bcu_h_dt9(ii,jj,k) +dt9(i,j,k)
              bcu_h_dt12(ii,jj,k) = bcu_h_dt12(ii,jj,k) +dt12(i,j,k)
              bcu_h_dq4(ii,jj,k)  = bcu_h_dq4(ii,jj,k)    + dq4(i,j,k)
              bcu_h_dq9(ii,jj,k)  = bcu_h_dq9(ii,jj,k)    + dq9(i,j,k)
              bcu_h_dq12(ii,jj,k) = bcu_h_dq12(ii,jj,k)   + dq12(i,j,k)
              bcu_h_dqcl4(ii,jj,k)  = bcu_h_dqcl4(ii,jj,k)  + dqcl4(i,j,k)
              bcu_h_dqcl9(ii,jj,k)  = bcu_h_dqcl9(ii,jj,k)  + dqcl9(i,j,k)
              bcu_h_dqcl12(ii,jj,k) = bcu_h_dqcl12(ii,jj,k) + dqcl12(i,j,k)
              bcu_h_dqcf4(ii,jj,k)  = bcu_h_dqcf4(ii,jj,k)  + dqcf4(i,j,k)
              bcu_h_dqcf3(ii,jj,k)  = bcu_h_dqcf3(ii,jj,k)  + dqcf3(i,j,k)
              bcu_h_dqcf12(ii,jj,k) = bcu_h_dqcf12(ii,jj,k) + dqcf12(i,j,k)

              bcu_h_thw(ii,jj,k)   = bcu_h_thw(ii,jj,k)   + thw_work
              bcu_h_thvw(ii,jj,k)  = bcu_h_thvw(ii,jj,k)  + thvw_work
              bcu_h_qw(ii,jj,k)    = bcu_h_qw(ii,jj,k)    + qw_work
              bcu_h_qclw(ii,jj,k)  = bcu_h_qclw(ii,jj,k)  + qclw_work
              bcu_h_qcfw(ii,jj,k)  = bcu_h_qcfw(ii,jj,k)  + qcfw_work
              bcu_h_qrainw(ii,jj,k)= bcu_h_qrainw(ii,jj,k)+ qrainw_work
              bcu_h_ww(ii,jj,k)    = bcu_h_ww(ii,jj,k)    + ww_work
              bcu_h_w3(ii,jj,k)    = bcu_h_w3(ii,jj,k)    + w3_work
              bcu_h_uw(ii,jj,k)    = bcu_h_uw(ii,jj,k)    + uw_work
              bcu_h_vw(ii,jj,k)    = bcu_h_vw(ii,jj,k)    + vw_work
              bcu_h_uu(ii,jj,k)    = bcu_h_uu(ii,jj,k)    + uu_work
              bcu_h_vv(ii,jj,k)    = bcu_h_vv(ii,jj,k)    + vv_work
              bcu_h_uv(ii,jj,k)    = bcu_h_uv(ii,jj,k)    + uv_work
              bcu_h_dpx(ii,jj,k)   = bcu_h_dpx(ii,jj,k)   + dpx_work
              bcu_h_dpy(ii,jj,k)   = bcu_h_dpy(ii,jj,k)   + dpy_work
              bcu_h_uth(ii,jj,k)   = bcu_h_uth(ii,jj,k)   + uth_work
              bcu_h_uthv(ii,jj,k)  = bcu_h_uthv(ii,jj,k)  + uthv_work
              bcu_h_vth(ii,jj,k)   = bcu_h_vth(ii,jj,k)   + vth_work
              bcu_h_vthv(ii,jj,k)  = bcu_h_vthv(ii,jj,k)  + vthv_work
              bcu_h_uq(ii,jj,k)    = bcu_h_uq(ii,jj,k)   + uq_work
              bcu_h_vq(ii,jj,k)    = bcu_h_vq(ii,jj,k)   + vq_work
              bcu_h_wp(ii,jj,k)    = bcu_h_wp(ii,jj,k)   + wp_work
              IF (l_qgraup) THEN
                bcu_h_qgraup(ii,jj,k) = bcu_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
                bcu_h_qgraupw(ii,jj,k)= bcu_h_qgraupw(ii,jj,k)+ qgraupw_work
              END IF
              IF (l_sect30) THEN
                bcu_h_dt30(ii,jj,k) = bcu_h_dt30(ii,jj,k) +dt30(i,j,k)
                bcu_h_dq30(ii,jj,k) = bcu_h_dq30(ii,jj,k) +dq30(i,j,k)
                bcu_h_dqcl30(ii,jj,k) = bcu_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
                bcu_h_dqcf30(ii,jj,k) = bcu_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
                bcu_h_dqrain30(ii,jj,k) = bcu_h_dqrain30(ii,jj,k)    &
                                            +dqrain30(i,j,k)
                bcu_h_drho(ii,jj,k) = bcu_h_drho(ii,jj,k) +drho(i,j,k)
              END IF
              IF (l_sect30 .AND. l_qgraup ) THEN
                bcu_h_dqgr30(ii,jj,k) = bcu_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
              END IF

            END IF     ! bcu
          END IF       ! Updraughts relative

          !-----------------------------------------------------------------
          ! WG1 - Strong updraughts
          !-----------------------------------------------------------------

          IF (l_wg1 .AND. w_h_prime(i,j,k) > 1.0) THEN

            IF (thv_h_prime(i,j,k) > 0.0) THEN

              wg1_factor(ii,jj,k) = wg1_factor(ii,jj,k)+1.0
              wg1_h_rho(ii,jj,k) = wg1_h_rho(ii,jj,k) + density(i,j,k)
              wg1_h_w(ii,jj,k) = wg1_h_w(ii,jj,k)+ w(i,j,k)
              wg1_h_u(ii,jj,k) = wg1_h_u(ii,jj,k)+ u(i,j,k)
              wg1_h_v(ii,jj,k) = wg1_h_v(ii,jj,k)+ v(i,j,k)
              wg1_h_th(ii,jj,k) = wg1_h_th(ii,jj,k) +th_h_prime(i,j,k)
              wg1_h_thv(ii,jj,k) = wg1_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
              wg1_h_rh(ii,jj,k) = wg1_h_rh(ii,jj,k) +rh(i,j,k)
              wg1_h_q(ii,jj,k) = wg1_h_q(ii,jj,k) +q(i,j,k)
              wg1_h_qcl(ii,jj,k) = wg1_h_qcl(ii,jj,k) +qcl(i,j,k)
              wg1_h_qcf(ii,jj,k) = wg1_h_qcf(ii,jj,k) +qcf(i,j,k)
              wg1_h_qrain(ii,jj,k) = wg1_h_qrain(ii,jj,k) +qrain(i,j,k)
              wg1_h_dt1(ii,jj,k)  = wg1_h_dt1(ii,jj,k)  + dt1(i,j,k)
              wg1_h_dt2(ii,jj,k)  = wg1_h_dt2(ii,jj,k)  + dt2(i,j,k)
              wg1_h_dt4(ii,jj,k)  = wg1_h_dt4(ii,jj,k)  + dt4(i,j,k)
              wg1_h_dt9(ii,jj,k)  = wg1_h_dt9(ii,jj,k)  + dt9(i,j,k)
              wg1_h_dt12(ii,jj,k) = wg1_h_dt12(ii,jj,k) + dt12(i,j,k)
              wg1_h_dq4(ii,jj,k)  = wg1_h_dq4(ii,jj,k)  + dq4(i,j,k)
              wg1_h_dq9(ii,jj,k)  = wg1_h_dq9(ii,jj,k)  + dq9(i,j,k)
              wg1_h_dq12(ii,jj,k) = wg1_h_dq12(ii,jj,k) + dq12(i,j,k)
              wg1_h_dqcl4(ii,jj,k)  = wg1_h_dqcl4(ii,jj,k)   +dqcl4(i,j,k)
              wg1_h_dqcl9(ii,jj,k)  = wg1_h_dqcl9(ii,jj,k)   +dqcl9(i,j,k)
              wg1_h_dqcl12(ii,jj,k) = wg1_h_dqcl12(ii,jj,k) +dqcl12(i,j,k)
              wg1_h_dqcf4(ii,jj,k)  = wg1_h_dqcf4(ii,jj,k)   +dqcf4(i,j,k)
              wg1_h_dqcf3(ii,jj,k)  = wg1_h_dqcf3(ii,jj,k)   +dqcf3(i,j,k)
              wg1_h_dqcf12(ii,jj,k) = wg1_h_dqcf12(ii,jj,k) +dqcf12(i,j,k)
              wg1_h_thw(ii,jj,k)   = wg1_h_thw(ii,jj,k)   + thw_work
              wg1_h_thvw(ii,jj,k)  = wg1_h_thvw(ii,jj,k)  + thvw_work
              wg1_h_qw(ii,jj,k)    = wg1_h_qw(ii,jj,k)    + qw_work
              wg1_h_qclw(ii,jj,k)  = wg1_h_qclw(ii,jj,k)  + qclw_work
              wg1_h_qcfw(ii,jj,k)  = wg1_h_qcfw(ii,jj,k)  + qcfw_work
              wg1_h_qrainw(ii,jj,k)= wg1_h_qrainw(ii,jj,k)+ qrainw_work
              wg1_h_ww(ii,jj,k)    = wg1_h_ww(ii,jj,k)    + ww_work
              wg1_h_w3(ii,jj,k)    = wg1_h_w3(ii,jj,k)    + w3_work
              wg1_h_uw(ii,jj,k)    = wg1_h_uw(ii,jj,k)    + uw_work
              wg1_h_vw(ii,jj,k)    = wg1_h_vw(ii,jj,k)    + vw_work
              wg1_h_uu(ii,jj,k)    = wg1_h_uu(ii,jj,k)    + uu_work
              wg1_h_vv(ii,jj,k)    = wg1_h_vv(ii,jj,k)    + vv_work
              wg1_h_uv(ii,jj,k)    = wg1_h_uv(ii,jj,k)    + uv_work
              wg1_h_dpx(ii,jj,k) = wg1_h_dpx(ii,jj,k) + dpx_work
              wg1_h_dpy(ii,jj,k) = wg1_h_dpy(ii,jj,k) + dpy_work
              wg1_h_uth(ii,jj,k)   = wg1_h_uth(ii,jj,k)   + uth_work
              wg1_h_uthv(ii,jj,k)  = wg1_h_uthv(ii,jj,k)  + uthv_work
              wg1_h_vth(ii,jj,k)   = wg1_h_vth(ii,jj,k)   + vth_work
              wg1_h_vthv(ii,jj,k)  = wg1_h_vthv(ii,jj,k)  + vthv_work
              wg1_h_uq(ii,jj,k)    = wg1_h_uq(ii,jj,k)   + uq_work
              wg1_h_vq(ii,jj,k)    = wg1_h_vq(ii,jj,k)   + vq_work
              wg1_h_wp(ii,jj,k)    = wg1_h_wp(ii,jj,k)   + wp_work

              IF (l_qgraup) THEN
                wg1_h_qgraup(ii,jj,k)= wg1_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
                wg1_h_qgraupw(ii,jj,k)= wg1_h_qgraupw(ii,jj,k)+ qgraupw_work
              END IF
              IF (l_sect30) THEN
                wg1_h_dt30(ii,jj,k) = wg1_h_dt30(ii,jj,k) +dt30(i,j,k)
                wg1_h_dq30(ii,jj,k) = wg1_h_dq30(ii,jj,k) +dq30(i,j,k)
                wg1_h_dqcl30(ii,jj,k) = wg1_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
                wg1_h_dqcf30(ii,jj,k) = wg1_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
                wg1_h_dqrain30(ii,jj,k) = wg1_h_dqrain30(ii,jj,k)  &
                                             +dqrain30(i,j,k)
                wg1_h_drho(ii,jj,k) = wg1_h_drho(ii,jj,k) +drho(i,j,k)
              END IF
              IF (l_sect30 .AND. l_qgraup ) THEN
                wg1_h_dqgr30(ii,jj,k) = wg1_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
              END IF
            END IF      !  WG1
          END IF       !   WG1
          !-----------------------wg1 tests --------------------------------

          !--------------------------------------------------------------------
          ! Upwards
          !--------------------------------------------------------------------
          IF (w(i,j,k) > 0.0) THEN
            IF (l_acw) THEN
              acw_factor(ii,jj,k) = acw_factor(ii,jj,k)+1.0
              acw_h_rho(ii,jj,k) = acw_h_rho(ii,jj,k) + density(i,j,k)
              acw_h_w(ii,jj,k) = acw_h_w(ii,jj,k)+ w(i,j,k)
              acw_h_u(ii,jj,k) = acw_h_u(ii,jj,k)+ u(i,j,k)
              acw_h_v(ii,jj,k) = acw_h_v(ii,jj,k)+ v(i,j,k)
              acw_h_th(ii,jj,k) = acw_h_th(ii,jj,k) +th_h_prime(i,j,k)
              acw_h_thv(ii,jj,k) = acw_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
              acw_h_rh(ii,jj,k) = acw_h_rh(ii,jj,k) +rh(i,j,k)
              acw_h_q(ii,jj,k) = acw_h_q(ii,jj,k) +q(i,j,k)
              acw_h_qcl(ii,jj,k) = acw_h_qcl(ii,jj,k) +qcl(i,j,k)
              acw_h_qcf(ii,jj,k) = acw_h_qcf(ii,jj,k) +qcf(i,j,k)
              acw_h_qrain(ii,jj,k) = acw_h_qrain(ii,jj,k) +qrain(i,j,k)
              acw_h_dt1(ii,jj,k) = acw_h_dt1(ii,jj,k) +dt1(i,j,k)
              acw_h_dt2(ii,jj,k) = acw_h_dt2(ii,jj,k) +dt2(i,j,k)
              acw_h_dt4(ii,jj,k) = acw_h_dt4(ii,jj,k) +dt4(i,j,k)
              acw_h_dt9(ii,jj,k) = acw_h_dt9(ii,jj,k) +dt9(i,j,k)
              acw_h_dt12(ii,jj,k) = acw_h_dt12(ii,jj,k) +dt12(i,j,k)
              acw_h_dq4(ii,jj,k) = acw_h_dq4(ii,jj,k) +dq4(i,j,k)
              acw_h_dq9(ii,jj,k) = acw_h_dq9(ii,jj,k) +dq9(i,j,k)
              acw_h_dq12(ii,jj,k) = acw_h_dq12(ii,jj,k) +dq12(i,j,k)
              acw_h_dqcl4(ii,jj,k) = acw_h_dqcl4(ii,jj,k)   +dqcl4(i,j,k)
              acw_h_dqcl9(ii,jj,k) = acw_h_dqcl9(ii,jj,k)   +dqcl9(i,j,k)
              acw_h_dqcl12(ii,jj,k) = acw_h_dqcl12(ii,jj,k) +dqcl12(i,j,k)
              acw_h_dqcf4(ii,jj,k) = acw_h_dqcf4(ii,jj,k)   +dqcf4(i,j,k)
              acw_h_dqcf3(ii,jj,k) = acw_h_dqcf3(ii,jj,k)   +dqcf3(i,j,k)
              acw_h_dqcf12(ii,jj,k) = acw_h_dqcf12(ii,jj,k) +dqcf12(i,j,k)
              acw_h_thw(ii,jj,k)   = acw_h_thw(ii,jj,k)   + thw_work
              acw_h_thvw(ii,jj,k)  = acw_h_thvw(ii,jj,k)  + thvw_work
              acw_h_qw(ii,jj,k)    = acw_h_qw(ii,jj,k)    + qw_work
              acw_h_qclw(ii,jj,k)  = acw_h_qclw(ii,jj,k)  + qclw_work
              acw_h_qcfw(ii,jj,k)  = acw_h_qcfw(ii,jj,k)  + qcfw_work
              acw_h_qrainw(ii,jj,k)= acw_h_qrainw(ii,jj,k)+ qrainw_work
              acw_h_ww(ii,jj,k)    = acw_h_ww(ii,jj,k)    + ww_work
              acw_h_w3(ii,jj,k)    = acw_h_w3(ii,jj,k)    + w3_work
              acw_h_uw(ii,jj,k)    = acw_h_uw(ii,jj,k)    + uw_work
              acw_h_vw(ii,jj,k)    = acw_h_vw(ii,jj,k)    + vw_work
              acw_h_uu(ii,jj,k)    = acw_h_uu(ii,jj,k)    + uu_work
              acw_h_vv(ii,jj,k)    = acw_h_vv(ii,jj,k)    + vv_work
              acw_h_uv(ii,jj,k)    = acw_h_uv(ii,jj,k)    + uv_work
              acw_h_dpx(ii,jj,k)   = acw_h_dpx(ii,jj,k)   + dpx_work
              acw_h_dpy(ii,jj,k)   = acw_h_dpy(ii,jj,k)   + dpy_work
              acw_h_uth(ii,jj,k)   = acw_h_uth(ii,jj,k)   + uth_work
              acw_h_uthv(ii,jj,k)  = acw_h_uthv(ii,jj,k)  + uthv_work
              acw_h_vth(ii,jj,k)   = acw_h_vth(ii,jj,k)   + vth_work
              acw_h_vthv(ii,jj,k)  = acw_h_vthv(ii,jj,k)  + vthv_work
              acw_h_uq(ii,jj,k)    = acw_h_uq(ii,jj,k)   + uq_work
              acw_h_vq(ii,jj,k)    = acw_h_vq(ii,jj,k)   + vq_work
              acw_h_wp(ii,jj,k)    = acw_h_wp(ii,jj,k)   + wp_work
              IF (l_qgraup) THEN
                acw_h_qgraup(ii,jj,k)= acw_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
                acw_h_qgraupw(ii,jj,k)= acw_h_qgraupw(ii,jj,k)+ qgraupw_work
              END IF
              IF (l_sect30) THEN
                acw_h_dt30(ii,jj,k) = acw_h_dt30(ii,jj,k) +dt30(i,j,k)
                acw_h_dq30(ii,jj,k) = acw_h_dq30(ii,jj,k) +dq30(i,j,k)
                acw_h_dqcl30(ii,jj,k) = acw_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
                acw_h_dqcf30(ii,jj,k) = acw_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
                acw_h_dqrain30(ii,jj,k) = acw_h_dqrain30(ii,jj,k) + &
                                            dqrain30(i,j,k)
                acw_h_drho(ii,jj,k) = acw_h_drho(ii,jj,k) +drho(i,j,k)
              END IF
              IF (l_sect30 .AND. l_qgraup ) THEN
                acw_h_dqgr30(ii,jj,k) = acw_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
              END IF
            END IF  ! acw

            !-------------------------------------------------------------
            ! Buoyant plumes
            !-------------------------------------------------------------
            IF (l_bcw .AND. thv_h_prime(i,j,k) > 0.0) THEN ! BCU

              bcw_factor(ii,jj,k) = bcw_factor(ii,jj,k)+1.0
              bcw_h_rho(ii,jj,k) = bcw_h_rho(ii,jj,k) + density(i,j,k)
              bcw_h_w(ii,jj,k) = bcw_h_w(ii,jj,k)+ w(i,j,k)
              bcw_h_u(ii,jj,k) = bcw_h_u(ii,jj,k)+ u(i,j,k)
              bcw_h_v(ii,jj,k) = bcw_h_v(ii,jj,k)+ v(i,j,k)
              bcw_h_th(ii,jj,k)  = bcw_h_th(ii,jj,k) +th_h_prime(i,j,k)
              bcw_h_thv(ii,jj,k) = bcw_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
              bcw_h_rh(ii,jj,k) = bcw_h_rh(ii,jj,k) +rh(i,j,k)
              bcw_h_q(ii,jj,k)   = bcw_h_q(ii,jj,k) +q(i,j,k)
              bcw_h_qcl(ii,jj,k) = bcw_h_qcl(ii,jj,k) +qcl(i,j,k)
              bcw_h_qcf(ii,jj,k) = bcw_h_qcf(ii,jj,k) +qcf(i,j,k)
              bcw_h_qrain(ii,jj,k) = bcw_h_qrain(ii,jj,k) +qrain(i,j,k)
              bcw_h_dt1(ii,jj,k) = bcw_h_dt1(ii,jj,k) +dt1(i,j,k)
              bcw_h_dt2(ii,jj,k) = bcw_h_dt2(ii,jj,k) +dt2(i,j,k)
              bcw_h_dt4(ii,jj,k) = bcw_h_dt4(ii,jj,k) +dt4(i,j,k)
              bcw_h_dt9(ii,jj,k) = bcw_h_dt9(ii,jj,k) +dt9(i,j,k)
              bcw_h_dt12(ii,jj,k) = bcw_h_dt12(ii,jj,k) +dt12(i,j,k)
              bcw_h_dq4(ii,jj,k)  = bcw_h_dq4(ii,jj,k)    + dq4(i,j,k)
              bcw_h_dq9(ii,jj,k)  = bcw_h_dq9(ii,jj,k)    + dq9(i,j,k)
              bcw_h_dq12(ii,jj,k) = bcw_h_dq12(ii,jj,k)   + dq12(i,j,k)
              bcw_h_dqcl4(ii,jj,k)  = bcw_h_dqcl4(ii,jj,k)  + dqcl4(i,j,k)
              bcw_h_dqcl9(ii,jj,k)  = bcw_h_dqcl9(ii,jj,k)  + dqcl9(i,j,k)
              bcw_h_dqcl12(ii,jj,k) = bcw_h_dqcl12(ii,jj,k) + dqcl12(i,j,k)
              bcw_h_dqcf4(ii,jj,k)  = bcw_h_dqcf4(ii,jj,k)  + dqcf4(i,j,k)
              bcw_h_dqcf3(ii,jj,k)  = bcw_h_dqcf3(ii,jj,k)  + dqcf3(i,j,k)
              bcw_h_dqcf12(ii,jj,k) = bcw_h_dqcf12(ii,jj,k) + dqcf12(i,j,k)

              bcw_h_thw(ii,jj,k)   = bcw_h_thw(ii,jj,k)   + thw_work
              bcw_h_thvw(ii,jj,k)  = bcw_h_thvw(ii,jj,k)  + thvw_work
              bcw_h_qw(ii,jj,k)    = bcw_h_qw(ii,jj,k)    + qw_work
              bcw_h_qclw(ii,jj,k)  = bcw_h_qclw(ii,jj,k)  + qclw_work
              bcw_h_qcfw(ii,jj,k)  = bcw_h_qcfw(ii,jj,k)  + qcfw_work
              bcw_h_qrainw(ii,jj,k)= bcw_h_qrainw(ii,jj,k)+ qrainw_work
              bcw_h_ww(ii,jj,k)    = bcw_h_ww(ii,jj,k)    + ww_work
              bcw_h_w3(ii,jj,k)    = bcw_h_w3(ii,jj,k)    + w3_work
              bcw_h_uw(ii,jj,k)    = bcw_h_uw(ii,jj,k)    + uw_work
              bcw_h_vw(ii,jj,k)    = bcw_h_vw(ii,jj,k)    + vw_work
              bcw_h_uu(ii,jj,k)    = bcw_h_uu(ii,jj,k)    + uu_work
              bcw_h_vv(ii,jj,k)    = bcw_h_vv(ii,jj,k)    + vv_work
              bcw_h_uv(ii,jj,k)    = bcw_h_uv(ii,jj,k)    + uv_work
              bcw_h_dpx(ii,jj,k)   = bcw_h_dpx(ii,jj,k)   + dpx_work
              bcw_h_dpy(ii,jj,k)   = bcw_h_dpy(ii,jj,k)   + dpy_work
              bcw_h_uth(ii,jj,k)   = bcw_h_uth(ii,jj,k)   + uth_work
              bcw_h_uthv(ii,jj,k)  = bcw_h_uthv(ii,jj,k)  + uthv_work
              bcw_h_vth(ii,jj,k)   = bcw_h_vth(ii,jj,k)   + vth_work
              bcw_h_vthv(ii,jj,k)  = bcw_h_vthv(ii,jj,k)  + vthv_work
              bcw_h_uq(ii,jj,k)    = bcw_h_uq(ii,jj,k)   + uq_work
              bcw_h_vq(ii,jj,k)    = bcw_h_vq(ii,jj,k)   + vq_work
              bcw_h_wp(ii,jj,k)    = bcw_h_wp(ii,jj,k)   + wp_work
              IF (l_qgraup) THEN
                bcw_h_qgraup(ii,jj,k) = bcw_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
                bcw_h_qgraupw(ii,jj,k)= bcw_h_qgraupw(ii,jj,k)+ qgraupw_work
              END IF
              IF (l_sect30) THEN
                bcw_h_dt30(ii,jj,k) = bcw_h_dt30(ii,jj,k) +dt30(i,j,k)
                bcw_h_dq30(ii,jj,k) = bcw_h_dq30(ii,jj,k) +dq30(i,j,k)
                bcw_h_dqcl30(ii,jj,k) = bcw_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
                bcw_h_dqcf30(ii,jj,k) = bcw_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
                bcw_h_dqrain30(ii,jj,k) = bcw_h_dqrain30(ii,jj,k)     &
                                             +dqrain30(i,j,k)
                bcw_h_drho(ii,jj,k) = bcw_h_drho(ii,jj,k) +drho(i,j,k)
              END IF
              IF (l_sect30 .AND. l_qgraup ) THEN
                bcw_h_dqgr30(ii,jj,k) = bcw_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
              END IF
            END IF     ! bcw

          END IF       ! upwards
          !--------------------------------------------------------------------
        END IF       ! ACC   cloud present
        !--------------------------------------------------------------------

        ! PPD  Should really have snow and graupel but not in UM and qcf gives
        !   too much ?
        IF (l_qgraup) THEN
          precip_water = qrain(i,j,k) + qgraup(i,j,k)
        ELSE
          precip_water = qrain(i,j,k)
        END IF

        IF (precip_water  > qw_crit .AND. w_h_prime(i,j,k) < 0.0) THEN
          IF (l_ppd) THEN
            ppd_factor(ii,jj,k) = ppd_factor(ii,jj,k)+1.0
            ppd_h_rho(ii,jj,k) = ppd_h_rho(ii,jj,k) + density(i,j,k)
            ppd_h_w(ii,jj,k) = ppd_h_w(ii,jj,k)+ w(i,j,k)
            ppd_h_u(ii,jj,k) = ppd_h_u(ii,jj,k)+ u(i,j,k)
            ppd_h_v(ii,jj,k) = ppd_h_v(ii,jj,k)+ v(i,j,k)
            ppd_h_th(ii,jj,k) = ppd_h_th(ii,jj,k) +th_h_prime(i,j,k)
            ppd_h_thv(ii,jj,k) = ppd_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
            ppd_h_rh(ii,jj,k) = ppd_h_rh(ii,jj,k) +rh(i,j,k)
            ppd_h_q(ii,jj,k) = ppd_h_q(ii,jj,k) +q(i,j,k)
            ppd_h_qcl(ii,jj,k) = ppd_h_qcl(ii,jj,k) +qcl(i,j,k)
            ppd_h_qcf(ii,jj,k) = ppd_h_qcf(ii,jj,k) +qcf(i,j,k)
            ppd_h_qrain(ii,jj,k) = ppd_h_qrain(ii,jj,k) +qrain(i,j,k)

            ppd_h_dt1(ii,jj,k) = ppd_h_dt1(ii,jj,k) +dt1(i,j,k)
            ppd_h_dt2(ii,jj,k) = ppd_h_dt2(ii,jj,k) +dt2(i,j,k)
            ppd_h_dt4(ii,jj,k) = ppd_h_dt4(ii,jj,k) +dt4(i,j,k)
            ppd_h_dt9(ii,jj,k) = ppd_h_dt9(ii,jj,k) +dt9(i,j,k)
            ppd_h_dt12(ii,jj,k) = ppd_h_dt12(ii,jj,k) +dt12(i,j,k)
            ppd_h_dq4(ii,jj,k) = ppd_h_dq4(ii,jj,k) +dq4(i,j,k)
            ppd_h_dq9(ii,jj,k) = ppd_h_dq9(ii,jj,k) +dq9(i,j,k)
            ppd_h_dq12(ii,jj,k) = ppd_h_dq12(ii,jj,k) +dq12(i,j,k)
            ppd_h_dqcl4(ii,jj,k) = ppd_h_dqcl4(ii,jj,k)   +dqcl4(i,j,k)
            ppd_h_dqcl9(ii,jj,k) = ppd_h_dqcl9(ii,jj,k)   +dqcl9(i,j,k)
            ppd_h_dqcl12(ii,jj,k) = ppd_h_dqcl12(ii,jj,k) +dqcl12(i,j,k)
            ppd_h_dqcf4(ii,jj,k) = ppd_h_dqcf4(ii,jj,k)   +dqcf4(i,j,k)
            ppd_h_dqcf3(ii,jj,k) = ppd_h_dqcf3(ii,jj,k)   +dqcf3(i,j,k)
            ppd_h_dqcf12(ii,jj,k) = ppd_h_dqcf12(ii,jj,k) +dqcf12(i,j,k)

            ppd_h_thw(ii,jj,k)   = ppd_h_thw(ii,jj,k)   + thw_work
            ppd_h_thvw(ii,jj,k)  = ppd_h_thvw(ii,jj,k)  + thvw_work
            ppd_h_qw(ii,jj,k)    = ppd_h_qw(ii,jj,k)    + qw_work
            ppd_h_qclw(ii,jj,k)  = ppd_h_qclw(ii,jj,k)  + qclw_work
            ppd_h_qcfw(ii,jj,k)  = ppd_h_qcfw(ii,jj,k)  + qcfw_work
            ppd_h_qrainw(ii,jj,k)= ppd_h_qrainw(ii,jj,k)+ qrainw_work
            ppd_h_ww(ii,jj,k)    = ppd_h_ww(ii,jj,k)    + ww_work
            ppd_h_w3(ii,jj,k)    = ppd_h_w3(ii,jj,k)    + w3_work
            ppd_h_uw(ii,jj,k)    = ppd_h_uw(ii,jj,k)    + uw_work
            ppd_h_vw(ii,jj,k)    = ppd_h_vw(ii,jj,k)    + vw_work
            ppd_h_uu(ii,jj,k)    = ppd_h_uu(ii,jj,k)    + uu_work
            ppd_h_vv(ii,jj,k)    = ppd_h_vv(ii,jj,k)    + vv_work
            ppd_h_uv(ii,jj,k)    = ppd_h_uv(ii,jj,k)    + uv_work
            ppd_h_dpx(ii,jj,k) = ppd_h_dpx(ii,jj,k) + dpx_work
            ppd_h_dpy(ii,jj,k) = ppd_h_dpy(ii,jj,k) + dpy_work
            ppd_h_uth(ii,jj,k)   = ppd_h_uth(ii,jj,k)   + uth_work
            ppd_h_uthv(ii,jj,k)  = ppd_h_uthv(ii,jj,k)  + uthv_work
            ppd_h_vth(ii,jj,k)   = ppd_h_vth(ii,jj,k)   + vth_work
            ppd_h_vthv(ii,jj,k)  = ppd_h_vthv(ii,jj,k)  + vthv_work
            ppd_h_uq(ii,jj,k)    = ppd_h_uq(ii,jj,k)   + uq_work
            ppd_h_vq(ii,jj,k)    = ppd_h_vq(ii,jj,k)   + vq_work
            ppd_h_wp(ii,jj,k)    = ppd_h_wp(ii,jj,k)   + wp_work
            IF (l_qgraup) THEN
              ppd_h_qgraup(ii,jj,k)= ppd_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
              ppd_h_qgraupw(ii,jj,k)= ppd_h_qgraupw(ii,jj,k)+ qgraupw_work
            END IF
            IF (l_sect30) THEN
              ppd_h_dt30(ii,jj,k) = ppd_h_dt30(ii,jj,k) +dt30(i,j,k)
              ppd_h_dq30(ii,jj,k) = ppd_h_dq30(ii,jj,k) +dq30(i,j,k)
              ppd_h_dqcl30(ii,jj,k) = ppd_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
              ppd_h_dqcf30(ii,jj,k) = ppd_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
              ppd_h_dqrain30(ii,jj,k) = ppd_h_dqrain30(ii,jj,k) &
                                                          +dqrain30(i,j,k)
              ppd_h_drho(ii,jj,k) = ppd_h_drho(ii,jj,k) +drho(i,j,k)
            END IF
            IF (l_sect30 .AND. l_qgraup ) THEN
              ppd_h_dqgr30(ii,jj,k) = ppd_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
            END IF
          END IF   ! l_ppd

          ! Precipitation downdraughts but also negatively buoyant
          IF (l_nbd .AND. thv_h_prime(i,j,k) < 0.0) THEN
            nbd_factor(ii,jj,k) = nbd_factor(ii,jj,k)+1.0
            nbd_h_rho(ii,jj,k) = nbd_h_rho(ii,jj,k) + density(i,j,k)
            nbd_h_w(ii,jj,k) = nbd_h_w(ii,jj,k)+ w(i,j,k)
            nbd_h_u(ii,jj,k) = nbd_h_u(ii,jj,k)+ u(i,j,k)
            nbd_h_v(ii,jj,k) = nbd_h_v(ii,jj,k)+ v(i,j,k)
            nbd_h_th(ii,jj,k) = nbd_h_th(ii,jj,k) +th_h_prime(i,j,k)
            nbd_h_thv(ii,jj,k) = nbd_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
            nbd_h_rh(ii,jj,k) = nbd_h_rh(ii,jj,k) +rh(i,j,k)
            nbd_h_q(ii,jj,k) = nbd_h_q(ii,jj,k) +q(i,j,k)
            nbd_h_qcl(ii,jj,k) = nbd_h_qcl(ii,jj,k) +qcl(i,j,k)
            nbd_h_qcf(ii,jj,k) = nbd_h_qcf(ii,jj,k) +qcf(i,j,k)
            nbd_h_qrain(ii,jj,k) = nbd_h_qrain(ii,jj,k) +qrain(i,j,k)
            nbd_h_dt1(ii,jj,k) = nbd_h_dt1(ii,jj,k) +dt1(i,j,k)
            nbd_h_dt2(ii,jj,k) = nbd_h_dt2(ii,jj,k) +dt2(i,j,k)
            nbd_h_dt4(ii,jj,k) = nbd_h_dt4(ii,jj,k) +dt4(i,j,k)
            nbd_h_dt9(ii,jj,k) = nbd_h_dt9(ii,jj,k) +dt9(i,j,k)
            nbd_h_dt12(ii,jj,k) = nbd_h_dt12(ii,jj,k) +dt12(i,j,k)
            nbd_h_dq4(ii,jj,k) = nbd_h_dq4(ii,jj,k) +dq4(i,j,k)
            nbd_h_dq9(ii,jj,k) = nbd_h_dq9(ii,jj,k) +dq9(i,j,k)
            nbd_h_dq12(ii,jj,k) = nbd_h_dq12(ii,jj,k) +dq12(i,j,k)
            nbd_h_dqcl4(ii,jj,k) = nbd_h_dqcl4(ii,jj,k)   +dqcl4(i,j,k)
            nbd_h_dqcl9(ii,jj,k) = nbd_h_dqcl9(ii,jj,k)   +dqcl9(i,j,k)
            nbd_h_dqcl12(ii,jj,k) = nbd_h_dqcl12(ii,jj,k) +dqcl12(i,j,k)
            nbd_h_dqcf4(ii,jj,k) = nbd_h_dqcf4(ii,jj,k)   +dqcf4(i,j,k)
            nbd_h_dqcf3(ii,jj,k) = nbd_h_dqcf3(ii,jj,k)   +dqcf3(i,j,k)
            nbd_h_dqcf12(ii,jj,k) = nbd_h_dqcf12(ii,jj,k) +dqcf12(i,j,k)
            nbd_h_thw(ii,jj,k)   = nbd_h_thw(ii,jj,k)   + thw_work
            nbd_h_thvw(ii,jj,k)  = nbd_h_thvw(ii,jj,k)  + thvw_work
            nbd_h_qw(ii,jj,k)    = nbd_h_qw(ii,jj,k)    + qw_work
            nbd_h_qclw(ii,jj,k)  = nbd_h_qclw(ii,jj,k)  + qclw_work
            nbd_h_qcfw(ii,jj,k)  = nbd_h_qcfw(ii,jj,k)  + qcfw_work
            nbd_h_qrainw(ii,jj,k)= nbd_h_qrainw(ii,jj,k)+ qrainw_work
            nbd_h_ww(ii,jj,k)    = nbd_h_ww(ii,jj,k)    + ww_work
            nbd_h_w3(ii,jj,k)    = nbd_h_w3(ii,jj,k)    + w3_work
            nbd_h_uw(ii,jj,k)    = nbd_h_uw(ii,jj,k)    + uw_work
            nbd_h_vw(ii,jj,k)    = nbd_h_vw(ii,jj,k)    + vw_work
            nbd_h_uu(ii,jj,k)    = nbd_h_uu(ii,jj,k)    + uu_work
            nbd_h_vv(ii,jj,k)    = nbd_h_vv(ii,jj,k)    + vv_work
            nbd_h_uv(ii,jj,k)    = nbd_h_uv(ii,jj,k)    + uv_work
            nbd_h_dpx(ii,jj,k) = nbd_h_dpx(ii,jj,k) + dpx_work
            nbd_h_dpy(ii,jj,k) = nbd_h_dpy(ii,jj,k) + dpy_work
            nbd_h_uth(ii,jj,k)   = nbd_h_uth(ii,jj,k)   + uth_work
            nbd_h_uthv(ii,jj,k)  = nbd_h_uthv(ii,jj,k)  + uthv_work
            nbd_h_vth(ii,jj,k)   = nbd_h_vth(ii,jj,k)   + vth_work
            nbd_h_vthv(ii,jj,k)  = nbd_h_vthv(ii,jj,k)  + vthv_work
            nbd_h_uq(ii,jj,k)    = nbd_h_uq(ii,jj,k)   + uq_work
            nbd_h_vq(ii,jj,k)    = nbd_h_vq(ii,jj,k)   + vq_work
            nbd_h_wp(ii,jj,k)    = nbd_h_wp(ii,jj,k)   + wp_work
            IF (l_qgraup) THEN
              nbd_h_qgraup(ii,jj,k)= nbd_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
              nbd_h_qgraupw(ii,jj,k)= nbd_h_qgraupw(ii,jj,k)+ qgraupw_work
            END IF
            IF (l_sect30) THEN
              nbd_h_dt30(ii,jj,k) = nbd_h_dt30(ii,jj,k) +dt30(i,j,k)
              nbd_h_dq30(ii,jj,k) = nbd_h_dq30(ii,jj,k) +dq30(i,j,k)
              nbd_h_dqcl30(ii,jj,k) = nbd_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
              nbd_h_dqcf30(ii,jj,k) = nbd_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
              nbd_h_dqrain30(ii,jj,k) = nbd_h_dqrain30(ii,jj,k)           &
                                                           +dqrain30(i,j,k)
              nbd_h_drho(ii,jj,k) = nbd_h_drho(ii,jj,k) +drho(i,j,k)
            END IF
            IF (l_sect30 .AND. l_qgraup ) THEN
              nbd_h_dqgr30(ii,jj,k) = nbd_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
            END IF
          END IF       ! nbd

        END IF       ! ppd

        ! Dependent on w not w' 
        IF (precip_water  > qw_crit .AND. w(i,j,k) < 0.0) THEN
          IF (l_ppw) THEN
            ppw_factor(ii,jj,k) = ppw_factor(ii,jj,k)+1.0
            ppw_h_rho(ii,jj,k) = ppw_h_rho(ii,jj,k) + density(i,j,k)
            ppw_h_w(ii,jj,k) = ppw_h_w(ii,jj,k)+ w(i,j,k)
            ppw_h_u(ii,jj,k) = ppw_h_u(ii,jj,k)+ u(i,j,k)
            ppw_h_v(ii,jj,k) = ppw_h_v(ii,jj,k)+ v(i,j,k)
            ppw_h_th(ii,jj,k) = ppw_h_th(ii,jj,k) +th_h_prime(i,j,k)
            ppw_h_thv(ii,jj,k) = ppw_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
            ppw_h_rh(ii,jj,k) = ppw_h_rh(ii,jj,k) +rh(i,j,k)
            ppw_h_q(ii,jj,k) = ppw_h_q(ii,jj,k) +q(i,j,k)
            ppw_h_qcl(ii,jj,k) = ppw_h_qcl(ii,jj,k) +qcl(i,j,k)
            ppw_h_qcf(ii,jj,k) = ppw_h_qcf(ii,jj,k) +qcf(i,j,k)
            ppw_h_qrain(ii,jj,k) = ppw_h_qrain(ii,jj,k) +qrain(i,j,k)

            ppw_h_dt1(ii,jj,k) = ppw_h_dt1(ii,jj,k) +dt1(i,j,k)
            ppw_h_dt2(ii,jj,k) = ppw_h_dt2(ii,jj,k) +dt2(i,j,k)
            ppw_h_dt4(ii,jj,k) = ppw_h_dt4(ii,jj,k) +dt4(i,j,k)
            ppw_h_dt9(ii,jj,k) = ppw_h_dt9(ii,jj,k) +dt9(i,j,k)
            ppw_h_dt12(ii,jj,k) = ppw_h_dt12(ii,jj,k) +dt12(i,j,k)
            ppw_h_dq4(ii,jj,k) = ppw_h_dq4(ii,jj,k) +dq4(i,j,k)
            ppw_h_dq9(ii,jj,k) = ppw_h_dq9(ii,jj,k) +dq9(i,j,k)
            ppw_h_dq12(ii,jj,k) = ppw_h_dq12(ii,jj,k) +dq12(i,j,k)
            ppw_h_dqcl4(ii,jj,k) = ppw_h_dqcl4(ii,jj,k)   +dqcl4(i,j,k)
            ppw_h_dqcl9(ii,jj,k) = ppw_h_dqcl9(ii,jj,k)   +dqcl9(i,j,k)
            ppw_h_dqcl12(ii,jj,k) = ppw_h_dqcl12(ii,jj,k) +dqcl12(i,j,k)
            ppw_h_dqcf4(ii,jj,k) = ppw_h_dqcf4(ii,jj,k)   +dqcf4(i,j,k)
            ppw_h_dqcf3(ii,jj,k) = ppw_h_dqcf3(ii,jj,k)   +dqcf3(i,j,k)
            ppw_h_dqcf12(ii,jj,k) = ppw_h_dqcf12(ii,jj,k) +dqcf12(i,j,k)

            ppw_h_thw(ii,jj,k)   = ppw_h_thw(ii,jj,k)   + thw_work
            ppw_h_thvw(ii,jj,k)  = ppw_h_thvw(ii,jj,k)  + thvw_work
            ppw_h_qw(ii,jj,k)    = ppw_h_qw(ii,jj,k)    + qw_work
            ppw_h_qclw(ii,jj,k)  = ppw_h_qclw(ii,jj,k)  + qclw_work
            ppw_h_qcfw(ii,jj,k)  = ppw_h_qcfw(ii,jj,k)  + qcfw_work
            ppw_h_qrainw(ii,jj,k)= ppw_h_qrainw(ii,jj,k)+ qrainw_work
            ppw_h_ww(ii,jj,k)    = ppw_h_ww(ii,jj,k)    + ww_work
            ppw_h_w3(ii,jj,k)    = ppw_h_w3(ii,jj,k)    + w3_work
            ppw_h_uw(ii,jj,k)    = ppw_h_uw(ii,jj,k)    + uw_work
            ppw_h_vw(ii,jj,k)    = ppw_h_vw(ii,jj,k)    + vw_work
            ppw_h_uu(ii,jj,k)    = ppw_h_uu(ii,jj,k)    + uu_work
            ppw_h_vv(ii,jj,k)    = ppw_h_vv(ii,jj,k)    + vv_work
            ppw_h_uv(ii,jj,k)    = ppw_h_uv(ii,jj,k)    + uv_work
            ppw_h_dpx(ii,jj,k) = ppw_h_dpx(ii,jj,k) + dpx_work
            ppw_h_dpy(ii,jj,k) = ppw_h_dpy(ii,jj,k) + dpy_work
            ppw_h_uth(ii,jj,k)   = ppw_h_uth(ii,jj,k)   + uth_work
            ppw_h_uthv(ii,jj,k)  = ppw_h_uthv(ii,jj,k)  + uthv_work
            ppw_h_vth(ii,jj,k)   = ppw_h_vth(ii,jj,k)   + vth_work
            ppw_h_vthv(ii,jj,k)  = ppw_h_vthv(ii,jj,k)  + vthv_work
            ppw_h_uq(ii,jj,k)    = ppw_h_uq(ii,jj,k)   + uq_work
            ppw_h_vq(ii,jj,k)    = ppw_h_vq(ii,jj,k)   + vq_work
            ppw_h_wp(ii,jj,k)    = ppw_h_wp(ii,jj,k)   + wp_work
            IF (l_qgraup) THEN
              ppw_h_qgraup(ii,jj,k)= ppw_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
              ppw_h_qgraupw(ii,jj,k)= ppw_h_qgraupw(ii,jj,k)+ qgraupw_work
            END IF
            IF (l_sect30) THEN
              ppw_h_dt30(ii,jj,k) = ppw_h_dt30(ii,jj,k) +dt30(i,j,k)
              ppw_h_dq30(ii,jj,k) = ppw_h_dq30(ii,jj,k) +dq30(i,j,k)
              ppw_h_dqcl30(ii,jj,k) = ppw_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
              ppw_h_dqcf30(ii,jj,k) = ppw_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
              ppw_h_dqrain30(ii,jj,k) = ppw_h_dqrain30(ii,jj,k)     &
                                                      +dqrain30(i,j,k)
              ppw_h_drho(ii,jj,k) = ppw_h_drho(ii,jj,k) +drho(i,j,k)
            END IF
            IF (l_sect30 .AND. l_qgraup ) THEN
              ppw_h_dqgr30(ii,jj,k) = ppw_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
            END IF
          END IF   ! l_ppw

          ! Precipitation downdraughts but also negatively buoyant
          IF (l_nbw .AND. thv_h_prime(i,j,k) < 0.0) THEN
            nbw_factor(ii,jj,k) = nbw_factor(ii,jj,k)+1.0
            nbw_h_rho(ii,jj,k) = nbw_h_rho(ii,jj,k) + density(i,j,k)
            nbw_h_w(ii,jj,k) = nbw_h_w(ii,jj,k)+ w(i,j,k)
            nbw_h_u(ii,jj,k) = nbw_h_u(ii,jj,k)+ u(i,j,k)
            nbw_h_v(ii,jj,k) = nbw_h_v(ii,jj,k)+ v(i,j,k)
            nbw_h_th(ii,jj,k) = nbw_h_th(ii,jj,k) +th_h_prime(i,j,k)
            nbw_h_thv(ii,jj,k) = nbw_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
            nbw_h_rh(ii,jj,k) = nbw_h_rh(ii,jj,k) +rh(i,j,k)
            nbw_h_q(ii,jj,k) = nbw_h_q(ii,jj,k) +q(i,j,k)
            nbw_h_qcl(ii,jj,k) = nbw_h_qcl(ii,jj,k) +qcl(i,j,k)
            nbw_h_qcf(ii,jj,k) = nbw_h_qcf(ii,jj,k) +qcf(i,j,k)
            nbw_h_qrain(ii,jj,k) = nbw_h_qrain(ii,jj,k) +qrain(i,j,k)
            nbw_h_dt1(ii,jj,k) = nbw_h_dt1(ii,jj,k) +dt1(i,j,k)
            nbw_h_dt2(ii,jj,k) = nbw_h_dt2(ii,jj,k) +dt2(i,j,k)
            nbw_h_dt4(ii,jj,k) = nbw_h_dt4(ii,jj,k) +dt4(i,j,k)
            nbw_h_dt9(ii,jj,k) = nbw_h_dt9(ii,jj,k) +dt9(i,j,k)
            nbw_h_dt12(ii,jj,k) = nbw_h_dt12(ii,jj,k) +dt12(i,j,k)
            nbw_h_dq4(ii,jj,k) = nbw_h_dq4(ii,jj,k) +dq4(i,j,k)
            nbw_h_dq9(ii,jj,k) = nbw_h_dq9(ii,jj,k) +dq9(i,j,k)
            nbw_h_dq12(ii,jj,k) = nbw_h_dq12(ii,jj,k) +dq12(i,j,k)
            nbw_h_dqcl4(ii,jj,k) = nbw_h_dqcl4(ii,jj,k)   +dqcl4(i,j,k)
            nbw_h_dqcl9(ii,jj,k) = nbw_h_dqcl9(ii,jj,k)   +dqcl9(i,j,k)
            nbw_h_dqcl12(ii,jj,k) = nbw_h_dqcl12(ii,jj,k) +dqcl12(i,j,k)
            nbw_h_dqcf4(ii,jj,k) = nbw_h_dqcf4(ii,jj,k)   +dqcf4(i,j,k)
            nbw_h_dqcf3(ii,jj,k) = nbw_h_dqcf3(ii,jj,k)   +dqcf3(i,j,k)
            nbw_h_dqcf12(ii,jj,k) = nbw_h_dqcf12(ii,jj,k) +dqcf12(i,j,k)
            nbw_h_thw(ii,jj,k)   = nbw_h_thw(ii,jj,k)   + thw_work
            nbw_h_thvw(ii,jj,k)  = nbw_h_thvw(ii,jj,k)  + thvw_work
            nbw_h_qw(ii,jj,k)    = nbw_h_qw(ii,jj,k)    + qw_work
            nbw_h_qclw(ii,jj,k)  = nbw_h_qclw(ii,jj,k)  + qclw_work
            nbw_h_qcfw(ii,jj,k)  = nbw_h_qcfw(ii,jj,k)  + qcfw_work
            nbw_h_qrainw(ii,jj,k)= nbw_h_qrainw(ii,jj,k)+ qrainw_work
            nbw_h_ww(ii,jj,k)    = nbw_h_ww(ii,jj,k)    + ww_work
            nbw_h_w3(ii,jj,k)    = nbw_h_w3(ii,jj,k)    + w3_work
            nbw_h_uw(ii,jj,k)    = nbw_h_uw(ii,jj,k)    + uw_work
            nbw_h_vw(ii,jj,k)    = nbw_h_vw(ii,jj,k)    + vw_work
            nbw_h_uu(ii,jj,k)    = nbw_h_uu(ii,jj,k)    + uu_work
            nbw_h_vv(ii,jj,k)    = nbw_h_vv(ii,jj,k)    + vv_work
            nbw_h_uv(ii,jj,k)    = nbw_h_uv(ii,jj,k)    + uv_work
            nbw_h_dpx(ii,jj,k) = nbw_h_dpx(ii,jj,k) + dpx_work
            nbw_h_dpy(ii,jj,k) = nbw_h_dpy(ii,jj,k) + dpy_work
            nbw_h_uth(ii,jj,k)   = nbw_h_uth(ii,jj,k)   + uth_work
            nbw_h_uthv(ii,jj,k)  = nbw_h_uthv(ii,jj,k)  + uthv_work
            nbw_h_vth(ii,jj,k)   = nbw_h_vth(ii,jj,k)   + vth_work
            nbw_h_vthv(ii,jj,k)  = nbw_h_vthv(ii,jj,k)  + vthv_work
            nbw_h_uq(ii,jj,k)    = nbw_h_uq(ii,jj,k)   + uq_work
            nbw_h_vq(ii,jj,k)    = nbw_h_vq(ii,jj,k)   + vq_work
            nbw_h_wp(ii,jj,k)    = nbw_h_wp(ii,jj,k)   + wp_work
            IF (l_qgraup) THEN
              nbw_h_qgraup(ii,jj,k)= nbw_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
              nbw_h_qgraupw(ii,jj,k)= nbw_h_qgraupw(ii,jj,k)+ qgraupw_work
            END IF
            IF (l_sect30) THEN
              nbw_h_dt30(ii,jj,k) = nbw_h_dt30(ii,jj,k) +dt30(i,j,k)
              nbw_h_dq30(ii,jj,k) = nbw_h_dq30(ii,jj,k) +dq30(i,j,k)
              nbw_h_dqcl30(ii,jj,k) = nbw_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
              nbw_h_dqcf30(ii,jj,k) = nbw_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
              nbw_h_dqrain30(ii,jj,k) = nbw_h_dqrain30(ii,jj,k)    &
                                                              +dqrain30(i,j,k)
              nbw_h_drho(ii,jj,k) = nbw_h_drho(ii,jj,k) +drho(i,j,k)
            END IF
            IF (l_sect30 .AND. l_qgraup ) THEN
              nbw_h_dqgr30(ii,jj,k) = nbw_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
            END IF
          END IF       ! nbd

        END IF       ! 

        ! NID  Negatively buoyant downdraughts including ice
        !      Still using w' < 0.0 condition
        IF (l_nid) THEN
          IF (l_qgraup) THEN
            precip_water = qrain(i,j,k) +qcf(i,j,k)+ qgraup(i,j,k)
          ELSE
            precip_water = qrain(i,j,k) +qcf(i,j,k)
          END IF

          IF (precip_water  > qw_crit .AND. w_h_prime(i,j,k) < 0.0) THEN
            ! Precipitation downdraughts but also negatively buoyant
            IF (thv_h_prime(i,j,k) < 0.0) THEN
              nid_factor(ii,jj,k) = nid_factor(ii,jj,k)+1.0
              nid_h_rho(ii,jj,k) = nid_h_rho(ii,jj,k) + density(i,j,k)
              nid_h_w(ii,jj,k) = nid_h_w(ii,jj,k)+ w(i,j,k)
              nid_h_u(ii,jj,k) = nid_h_u(ii,jj,k)+ u(i,j,k)
              nid_h_v(ii,jj,k) = nid_h_v(ii,jj,k)+ v(i,j,k)
              nid_h_th(ii,jj,k) = nid_h_th(ii,jj,k) +th_h_prime(i,j,k)
              nid_h_thv(ii,jj,k) = nid_h_thv(ii,jj,k) +thv_h_prime(i,j,k)
              nid_h_rh(ii,jj,k) = nid_h_rh(ii,jj,k) +rh(i,j,k)
              nid_h_q(ii,jj,k) = nid_h_q(ii,jj,k) +q(i,j,k)
              nid_h_qcl(ii,jj,k) = nid_h_qcl(ii,jj,k) +qcl(i,j,k)
              nid_h_qcf(ii,jj,k) = nid_h_qcf(ii,jj,k) +qcf(i,j,k)
              nid_h_qrain(ii,jj,k) = nid_h_qrain(ii,jj,k) +qrain(i,j,k)
              nid_h_dt1(ii,jj,k) = nid_h_dt1(ii,jj,k) +dt1(i,j,k)
              nid_h_dt2(ii,jj,k) = nid_h_dt2(ii,jj,k) +dt2(i,j,k)
              nid_h_dt4(ii,jj,k) = nid_h_dt4(ii,jj,k) +dt4(i,j,k)
              nid_h_dt9(ii,jj,k) = nid_h_dt9(ii,jj,k) +dt9(i,j,k)
              nid_h_dt12(ii,jj,k) = nid_h_dt12(ii,jj,k) +dt12(i,j,k)
              nid_h_dq4(ii,jj,k) = nid_h_dq4(ii,jj,k) +dq4(i,j,k)
              nid_h_dq9(ii,jj,k) = nid_h_dq9(ii,jj,k) +dq9(i,j,k)
              nid_h_dq12(ii,jj,k) = nid_h_dq12(ii,jj,k) +dq12(i,j,k)
              nid_h_dqcl4(ii,jj,k) = nid_h_dqcl4(ii,jj,k)   +dqcl4(i,j,k)
              nid_h_dqcl9(ii,jj,k) = nid_h_dqcl9(ii,jj,k)   +dqcl9(i,j,k)
              nid_h_dqcl12(ii,jj,k) = nid_h_dqcl12(ii,jj,k) +dqcl12(i,j,k)
              nid_h_dqcf4(ii,jj,k) = nid_h_dqcf4(ii,jj,k)   +dqcf4(i,j,k)
              nid_h_dqcf3(ii,jj,k) = nid_h_dqcf3(ii,jj,k)   +dqcf3(i,j,k)
              nid_h_dqcf12(ii,jj,k) = nid_h_dqcf12(ii,jj,k) +dqcf12(i,j,k)
              nid_h_thw(ii,jj,k)   = nid_h_thw(ii,jj,k)   + thw_work
              nid_h_thvw(ii,jj,k)  = nid_h_thvw(ii,jj,k)  + thvw_work
              nid_h_qw(ii,jj,k)    = nid_h_qw(ii,jj,k)    + qw_work
              nid_h_qclw(ii,jj,k)  = nid_h_qclw(ii,jj,k)  + qclw_work
              nid_h_qcfw(ii,jj,k)  = nid_h_qcfw(ii,jj,k)  + qcfw_work
              nid_h_qrainw(ii,jj,k)= nid_h_qrainw(ii,jj,k)+ qrainw_work
              nid_h_ww(ii,jj,k)    = nid_h_ww(ii,jj,k)    + ww_work
              nid_h_w3(ii,jj,k)    = nid_h_w3(ii,jj,k)    + w3_work
              nid_h_uw(ii,jj,k)    = nid_h_uw(ii,jj,k)    + uw_work
              nid_h_vw(ii,jj,k)    = nid_h_vw(ii,jj,k)    + vw_work
              nid_h_uu(ii,jj,k)    = nid_h_uu(ii,jj,k)    + uu_work
              nid_h_vv(ii,jj,k)    = nid_h_vv(ii,jj,k)    + vv_work
              nid_h_uv(ii,jj,k)    = nid_h_uv(ii,jj,k)    + uv_work
              nid_h_dpx(ii,jj,k)   = nid_h_dpx(ii,jj,k) + dpx_work
              nid_h_dpy(ii,jj,k)   = nid_h_dpy(ii,jj,k) + dpy_work
              nid_h_uth(ii,jj,k)   = nid_h_uth(ii,jj,k)   + uth_work
              nid_h_uthv(ii,jj,k)  = nid_h_uthv(ii,jj,k)  + uthv_work
              nid_h_vth(ii,jj,k)   = nid_h_vth(ii,jj,k)   + vth_work
              nid_h_vthv(ii,jj,k)  = nid_h_vthv(ii,jj,k)  + vthv_work
              nid_h_uq(ii,jj,k)    = nid_h_uq(ii,jj,k)   + uq_work
              nid_h_vq(ii,jj,k)    = nid_h_vq(ii,jj,k)   + vq_work
              nid_h_wp(ii,jj,k)    = nid_h_wp(ii,jj,k)   + wp_work
              IF (l_qgraup) THEN
                nid_h_qgraup(ii,jj,k)= nid_h_qgraup(ii,jj,k)+ qgraup(i,j,k)
                nid_h_qgraupw(ii,jj,k)= nid_h_qgraupw(ii,jj,k)+ qgraupw_work
              END IF
              IF (l_sect30) THEN
                nid_h_dt30(ii,jj,k) = nid_h_dt30(ii,jj,k) +dt30(i,j,k)
                nid_h_dq30(ii,jj,k) = nid_h_dq30(ii,jj,k) +dq30(i,j,k)
                nid_h_dqcl30(ii,jj,k) = nid_h_dqcl30(ii,jj,k) +dqcl30(i,j,k)
                nid_h_dqcf30(ii,jj,k) = nid_h_dqcf30(ii,jj,k) +dqcf30(i,j,k)
                nid_h_dqrain30(ii,jj,k) = nid_h_dqrain30(ii,jj,k)            &
                                                  +dqrain30(i,j,k)
                nid_h_drho(ii,jj,k) = nid_h_drho(ii,jj,k) +drho(i,j,k)
              END IF
              IF (l_sect30 .AND. l_qgraup ) THEN
                nid_h_dqgr30(ii,jj,k) = nid_h_dqgr30(ii,jj,k) +dqgr30(i,j,k)
              END IF
            END IF       ! buoyant

          END IF       ! precip
        END IF       ! nid

      END IF   ! above surface
    END DO
  END DO

END DO  ! k level loop
!$OMP END PARALLEL DO



! divide by total number of points

!$OMP PARALLEL DO PRIVATE(k, ii, jj, factor) DEFAULT(SHARED)

DO k=1,mlevs
  ! loop over coarse grid
  DO jj = 1,local_new_y
    DO ii = 1, local_new_x

      IF (all_factor(ii,jj,k) > 0.0) THEN
        all_a(ii,jj,k)  = all_factor(ii,jj,k)/ftot

        factor = 1.0/all_factor(ii,jj,k)
        all_rho(ii,jj,k)  = all_rho(ii,jj,k)*factor
        all_thw(ii,jj,k)  = all_thw(ii,jj,k)*factor
        all_thvw(ii,jj,k) = all_thvw(ii,jj,k)*factor
        all_qw(ii,jj,k) = all_qw(ii,jj,k)*factor
        all_qclw(ii,jj,k) = all_qclw(ii,jj,k)*factor
        all_qcfw(ii,jj,k) = all_qcfw(ii,jj,k)*factor
        all_qrainw(ii,jj,k) = all_qrainw(ii,jj,k)*factor
        all_ww(ii,jj,k) = all_ww(ii,jj,k)*factor
        all_w3(ii,jj,k) = all_w3(ii,jj,k)*factor
        all_uw(ii,jj,k) = all_uw(ii,jj,k)*factor
        all_vw(ii,jj,k) = all_vw(ii,jj,k)*factor
        all_uu(ii,jj,k) = all_uu(ii,jj,k)*factor
        all_vv(ii,jj,k) = all_vv(ii,jj,k)*factor
        all_uv(ii,jj,k) = all_uv(ii,jj,k)*factor
        all_dpx(ii,jj,k) = all_dpx(ii,jj,k)*factor
        all_dpy(ii,jj,k) = all_dpy(ii,jj,k)*factor
        all_uth(ii,jj,k)   = all_uth(ii,jj,k)*factor
        all_uthv(ii,jj,k)  = all_uthv(ii,jj,k)*factor
        all_vth(ii,jj,k)   = all_vth(ii,jj,k)*factor
        all_vthv(ii,jj,k)  = all_vthv(ii,jj,k)*factor
        all_uq(ii,jj,k)    = all_uq(ii,jj,k)*factor
        all_vq(ii,jj,k)    = all_vq(ii,jj,k)*factor
        all_wp(ii,jj,k)    = all_wp(ii,jj,k)*factor
        all_wp_hydro(ii,jj,k)    = all_wp_hydro(ii,jj,k)*factor
        
        IF (l_qgraup) THEN
          all_qgraupw(ii,jj,k) = all_qgraupw(ii,jj,k)*factor
        END IF

      ELSE      ! all points below surface
        all_a(ii,jj,k)    = rmdi         ! no values
        all_thw(ii,jj,k)  = rmdi
        all_thvw(ii,jj,k) = rmdi
        all_qw(ii,jj,k)   = rmdi
        all_qclw(ii,jj,k) =  rmdi
        all_qcfw(ii,jj,k) =  rmdi
        all_qrainw(ii,jj,k) =  rmdi
        all_ww(ii,jj,k) =  rmdi
        all_w3(ii,jj,k) =  rmdi
        all_uw(ii,jj,k) =  rmdi
        all_vw(ii,jj,k) =  rmdi
        all_uu(ii,jj,k) =  rmdi
        all_vv(ii,jj,k) =  rmdi
        all_uv(ii,jj,k) =  rmdi
        all_rho(ii,jj,k)  = rmdi
        all_dpx(ii,jj,k)  = rmdi         ! no values
        all_dpy(ii,jj,k)  = rmdi
        all_uth(ii,jj,k)   = rmdi
        all_uthv(ii,jj,k)  = rmdi
        all_vth(ii,jj,k)   = rmdi
        all_vthv(ii,jj,k)  = rmdi
        all_uq(ii,jj,k)    = rmdi
        all_vq(ii,jj,k)    = rmdi
        all_wp(ii,jj,k)    = rmdi
        all_wp_hydro(ii,jj,k)    = rmdi
        IF (l_qgraup) THEN
          all_qgraupw(ii,jj,k) = rmdi
        END IF

      END IF

    END DO   ! loop over ii
  END DO     ! loop over jj
END DO   ! loop over k
!$OMP END PARALLEL DO

!------------------------------------------------------------------------------
! Create means and set missing data below surface for partitions
!------------------------------------------------------------------------------
IF (l_acc) THEN
  CALL crmstyle_scale_rmdi(all_factor, acc_factor,                        &  
  acc_h_w, acc_h_u, acc_h_v, acc_h_th, acc_h_thv, acc_h_rho,              &
  acc_h_rh, acc_h_a, acc_h_dpx, acc_h_dpy,                                &
  acc_h_q,   acc_h_qcl, acc_h_qcf, acc_h_qrain, acc_h_qgraup,             &
  acc_h_dt1, acc_h_dt2, acc_h_dt4, acc_h_dt9, acc_h_dt12, acc_h_dt30,     &
  acc_h_dq4, acc_h_dq9, acc_h_dq12, acc_h_dq30,                           &
  acc_h_dqcl4, acc_h_dqcl9, acc_h_dqcl12, acc_h_dqcl30,                   &
  acc_h_dqcf4, acc_h_dqcf3, acc_h_dqcf12, acc_h_dqcf30,                   &
  acc_h_dqrain30, acc_h_dqgr30, acc_h_drho, acc_h_thw, acc_h_thvw,        &
  acc_h_qw, acc_h_qclw, acc_h_qcfw, acc_h_qrainw, acc_h_qgraupw,          &
  acc_h_uth, acc_h_vth, acc_h_uthv, acc_h_vthv, acc_h_uq, acc_h_vq,       &
  acc_h_wp,                                                               &
  acc_h_uw, acc_h_vw, acc_h_ww, acc_h_w3, acc_h_vv, acc_h_uu, acc_h_uv )
END IF

IF (l_acu) THEN
  CALL crmstyle_scale_rmdi(all_factor, acu_factor,                        &  
  acu_h_w, acu_h_u, acu_h_v, acu_h_th, acu_h_thv, acu_h_rho,              &
  acu_h_rh, acu_h_a, acu_h_dpx, acu_h_dpy,                                &
  acu_h_q,   acu_h_qcl, acu_h_qcf, acu_h_qrain, acu_h_qgraup,             &
  acu_h_dt1, acu_h_dt2, acu_h_dt4, acu_h_dt9, acu_h_dt12, acu_h_dt30,     &
  acu_h_dq4, acu_h_dq9, acu_h_dq12, acu_h_dq30,                           &
  acu_h_dqcl4, acu_h_dqcl9, acu_h_dqcl12, acu_h_dqcl30,                   &
  acu_h_dqcf4, acu_h_dqcf3, acu_h_dqcf12, acu_h_dqcf30,                   &
  acu_h_dqrain30, acu_h_dqgr30, acu_h_drho, acu_h_thw, acu_h_thvw,        &
  acu_h_qw, acu_h_qclw, acu_h_qcfw, acu_h_qrainw, acu_h_qgraupw,          &
  acu_h_uth, acu_h_vth, acu_h_uthv, acu_h_vthv, acu_h_uq, acu_h_vq,       &
  acu_h_wp,                                                               &
  acu_h_uw, acu_h_vw, acu_h_ww, acu_h_w3, acu_h_vv, acu_h_uu, acu_h_uv )
END IF

IF (l_bcu) THEN
  CALL crmstyle_scale_rmdi(all_factor, bcu_factor,                        &  
  bcu_h_w, bcu_h_u, bcu_h_v, bcu_h_th, bcu_h_thv, bcu_h_rho,              &
  bcu_h_rh, bcu_h_a, bcu_h_dpx, bcu_h_dpy,                                &
  bcu_h_q,   bcu_h_qcl, bcu_h_qcf, bcu_h_qrain, bcu_h_qgraup,             &
  bcu_h_dt1, bcu_h_dt2, bcu_h_dt4, bcu_h_dt9, bcu_h_dt12, bcu_h_dt30,     &
  bcu_h_dq4, bcu_h_dq9, bcu_h_dq12, bcu_h_dq30,                           &
  bcu_h_dqcl4, bcu_h_dqcl9, bcu_h_dqcl12, bcu_h_dqcl30,                   &
  bcu_h_dqcf4, bcu_h_dqcf3, bcu_h_dqcf12, bcu_h_dqcf30,                   &
  bcu_h_dqrain30, bcu_h_dqgr30, bcu_h_drho, bcu_h_thw, bcu_h_thvw,        &
  bcu_h_qw, bcu_h_qclw, bcu_h_qcfw, bcu_h_qrainw, bcu_h_qgraupw,          &
  bcu_h_uth, bcu_h_vth, bcu_h_uthv, bcu_h_vthv, bcu_h_uq, bcu_h_vq,       &
  bcu_h_wp,                                                               &
  bcu_h_uw, bcu_h_vw, bcu_h_ww, bcu_h_w3, bcu_h_vv, bcu_h_uu, bcu_h_uv )
END IF

IF (l_wg1) THEN
  CALL crmstyle_scale_rmdi(all_factor, wg1_factor,                        &  
  wg1_h_w, wg1_h_u, wg1_h_v, wg1_h_th, wg1_h_thv, wg1_h_rho,              &
  wg1_h_rh, wg1_h_a, wg1_h_dpx, wg1_h_dpy,                                &
  wg1_h_q,   wg1_h_qcl, wg1_h_qcf, wg1_h_qrain, wg1_h_qgraup,             &
  wg1_h_dt1, wg1_h_dt2, wg1_h_dt4, wg1_h_dt9, wg1_h_dt12, wg1_h_dt30,     &
  wg1_h_dq4, wg1_h_dq9, wg1_h_dq12, wg1_h_dq30,                           &
  wg1_h_dqcl4, wg1_h_dqcl9, wg1_h_dqcl12, wg1_h_dqcl30,                   &
  wg1_h_dqcf4, wg1_h_dqcf3, wg1_h_dqcf12, wg1_h_dqcf30,                   &
  wg1_h_dqrain30, wg1_h_dqgr30, wg1_h_drho, wg1_h_thw, wg1_h_thvw,        &
  wg1_h_qw, wg1_h_qclw, wg1_h_qcfw, wg1_h_qrainw, wg1_h_qgraupw,          &
  wg1_h_uth, wg1_h_vth, wg1_h_uthv, wg1_h_vthv, wg1_h_uq, wg1_h_vq,       &
  wg1_h_wp,                                                               &
  wg1_h_uw, wg1_h_vw, wg1_h_ww, wg1_h_w3, wg1_h_vv, wg1_h_uu, wg1_h_uv )
END IF

IF (l_ppd) THEN
  CALL crmstyle_scale_rmdi(all_factor, ppd_factor,                        &  
  ppd_h_w, ppd_h_u, ppd_h_v, ppd_h_th, ppd_h_thv, ppd_h_rho,              &
  ppd_h_rh, ppd_h_a, ppd_h_dpx, ppd_h_dpy,                                &
  ppd_h_q,   ppd_h_qcl, ppd_h_qcf, ppd_h_qrain, ppd_h_qgraup,             &
  ppd_h_dt1, ppd_h_dt2, ppd_h_dt4, ppd_h_dt9, ppd_h_dt12, ppd_h_dt30,     &
  ppd_h_dq4, ppd_h_dq9, ppd_h_dq12, ppd_h_dq30,                           &
  ppd_h_dqcl4, ppd_h_dqcl9, ppd_h_dqcl12, ppd_h_dqcl30,                   &
  ppd_h_dqcf4, ppd_h_dqcf3, ppd_h_dqcf12, ppd_h_dqcf30,                   &
  ppd_h_dqrain30, ppd_h_dqgr30, ppd_h_drho, ppd_h_thw, ppd_h_thvw,        &
  ppd_h_qw, ppd_h_qclw, ppd_h_qcfw, ppd_h_qrainw, ppd_h_qgraupw,          &
  ppd_h_uth, ppd_h_vth, ppd_h_uthv, ppd_h_vthv, ppd_h_uq, ppd_h_vq,       &
  ppd_h_wp,                                                               &
  ppd_h_uw, ppd_h_vw, ppd_h_ww, ppd_h_w3, ppd_h_vv, ppd_h_uu, ppd_h_uv )
END IF

IF (l_nbd) THEN
  CALL crmstyle_scale_rmdi(all_factor, nbd_factor,                        &  
  nbd_h_w, nbd_h_u, nbd_h_v, nbd_h_th, nbd_h_thv, nbd_h_rho,              &
  nbd_h_rh, nbd_h_a, nbd_h_dpx, nbd_h_dpy,                                &
  nbd_h_q,   nbd_h_qcl, nbd_h_qcf, nbd_h_qrain, nbd_h_qgraup,             &
  nbd_h_dt1, nbd_h_dt2, nbd_h_dt4, nbd_h_dt9, nbd_h_dt12, nbd_h_dt30,     &
  nbd_h_dq4, nbd_h_dq9, nbd_h_dq12, nbd_h_dq30,                           &
  nbd_h_dqcl4, nbd_h_dqcl9, nbd_h_dqcl12, nbd_h_dqcl30,                   &
  nbd_h_dqcf4, nbd_h_dqcf3, nbd_h_dqcf12, nbd_h_dqcf30,                   &
  nbd_h_dqrain30, nbd_h_dqgr30, nbd_h_drho, nbd_h_thw, nbd_h_thvw,        &
  nbd_h_qw, nbd_h_qclw, nbd_h_qcfw, nbd_h_qrainw, nbd_h_qgraupw,          &
  nbd_h_uth, nbd_h_vth, nbd_h_uthv, nbd_h_vthv, nbd_h_uq, nbd_h_vq,       &
  nbd_h_wp,                                                               &
  nbd_h_uw, nbd_h_vw, nbd_h_ww, nbd_h_w3, nbd_h_vv, nbd_h_uu, nbd_h_uv )
END IF

IF (l_nid) THEN
  CALL crmstyle_scale_rmdi(all_factor, nid_factor,                        &  
  nid_h_w, nid_h_u, nid_h_v, nid_h_th, nid_h_thv, nid_h_rho,              &
  nid_h_rh, nid_h_a, nid_h_dpx, nid_h_dpy,                                &
  nid_h_q,   nid_h_qcl, nid_h_qcf, nid_h_qrain, nid_h_qgraup,             &
  nid_h_dt1, nid_h_dt2, nid_h_dt4, nid_h_dt9, nid_h_dt12, nid_h_dt30,     &
  nid_h_dq4, nid_h_dq9, nid_h_dq12, nid_h_dq30,                           &
  nid_h_dqcl4, nid_h_dqcl9, nid_h_dqcl12, nid_h_dqcl30,                   &
  nid_h_dqcf4, nid_h_dqcf3, nid_h_dqcf12, nid_h_dqcf30,                   &
  nid_h_dqrain30, nid_h_dqgr30, nid_h_drho, nid_h_thw, nid_h_thvw,        &
  nid_h_qw, nid_h_qclw, nid_h_qcfw, nid_h_qrainw, nid_h_qgraupw,          &
  nid_h_uth, nid_h_vth, nid_h_uthv, nid_h_vthv, nid_h_uq, nid_h_vq,       &
  nid_h_wp,                                                               &
  nid_h_uw, nid_h_vw, nid_h_ww, nid_h_w3, nid_h_vv, nid_h_uu, nid_h_uv )
END IF

IF (l_adu) THEN
  CALL crmstyle_scale_rmdi(all_factor, adu_factor,                        &  
  adu_h_w, adu_h_u, adu_h_v, adu_h_th, adu_h_thv, adu_h_rho,              &
  adu_h_rh, adu_h_a, adu_h_dpx, adu_h_dpy,                                &
  adu_h_q,   adu_h_qcl, adu_h_qcf, adu_h_qrain, adu_h_qgraup,             &
  adu_h_dt1, adu_h_dt2, adu_h_dt4, adu_h_dt9, adu_h_dt12, adu_h_dt30,     &
  adu_h_dq4, adu_h_dq9, adu_h_dq12, adu_h_dq30,                           &
  adu_h_dqcl4, adu_h_dqcl9, adu_h_dqcl12, adu_h_dqcl30,                   &
  adu_h_dqcf4, adu_h_dqcf3, adu_h_dqcf12, adu_h_dqcf30,                   &
  adu_h_dqrain30, adu_h_dqgr30, adu_h_drho, adu_h_thw, adu_h_thvw,        &
  adu_h_qw, adu_h_qclw, adu_h_qcfw, adu_h_qrainw, adu_h_qgraupw,          &
  adu_h_uth, adu_h_vth, adu_h_uthv, adu_h_vthv, adu_h_uq, adu_h_vq,       &
  adu_h_wp,                                                               &
  adu_h_uw, adu_h_vw, adu_h_ww, adu_h_w3, adu_h_vv, adu_h_uu, adu_h_uv )
END IF

IF (l_acw) THEN
  CALL crmstyle_scale_rmdi(all_factor, acw_factor,                        &  
  acw_h_w, acw_h_u, acw_h_v, acw_h_th, acw_h_thv, acw_h_rho,              &
  acw_h_rh, acw_h_a, acw_h_dpx, acw_h_dpy,                                &
  acw_h_q,   acw_h_qcl, acw_h_qcf, acw_h_qrain, acw_h_qgraup,             &
  acw_h_dt1, acw_h_dt2, acw_h_dt4, acw_h_dt9, acw_h_dt12, acw_h_dt30,     &
  acw_h_dq4, acw_h_dq9, acw_h_dq12, acw_h_dq30,                           &
  acw_h_dqcl4, acw_h_dqcl9, acw_h_dqcl12, acw_h_dqcl30,                   &
  acw_h_dqcf4, acw_h_dqcf3, acw_h_dqcf12, acw_h_dqcf30,                   &
  acw_h_dqrain30, acw_h_dqgr30, acw_h_drho, acw_h_thw, acw_h_thvw,        &
  acw_h_qw, acw_h_qclw, acw_h_qcfw, acw_h_qrainw, acw_h_qgraupw,          &
  acw_h_uth, acw_h_vth, acw_h_uthv, acw_h_vthv, acw_h_uq, acw_h_vq,       &
  acw_h_wp,                                                               &
  acw_h_uw, acw_h_vw, acw_h_ww, acw_h_w3, acw_h_vv, acw_h_uu, acw_h_uv )
END IF

IF (l_bcw) THEN
  CALL crmstyle_scale_rmdi(all_factor, bcw_factor,                        &  
  bcw_h_w, bcw_h_u, bcw_h_v, bcw_h_th, bcw_h_thv, bcw_h_rho,              &
  bcw_h_rh, bcw_h_a, bcw_h_dpx, bcw_h_dpy,                                &
  bcw_h_q,   bcw_h_qcl, bcw_h_qcf, bcw_h_qrain, bcw_h_qgraup,             &
  bcw_h_dt1, bcw_h_dt2, bcw_h_dt4, bcw_h_dt9, bcw_h_dt12, bcw_h_dt30,     &
  bcw_h_dq4, bcw_h_dq9, bcw_h_dq12, bcw_h_dq30,                           &
  bcw_h_dqcl4, bcw_h_dqcl9, bcw_h_dqcl12, bcw_h_dqcl30,                   &
  bcw_h_dqcf4, bcw_h_dqcf3, bcw_h_dqcf12, bcw_h_dqcf30,                   &
  bcw_h_dqrain30, bcw_h_dqgr30, bcw_h_drho, bcw_h_thw, bcw_h_thvw,        &
  bcw_h_qw, bcw_h_qclw, bcw_h_qcfw, bcw_h_qrainw, bcw_h_qgraupw,          &
  bcw_h_uth, bcw_h_vth, bcw_h_uthv, bcw_h_vthv, bcw_h_uq, bcw_h_vq,       &
  bcw_h_wp,                                                               &
  bcw_h_uw, bcw_h_vw, bcw_h_ww, bcw_h_w3, bcw_h_vv, bcw_h_uu, bcw_h_uv )
END IF

IF (l_ucu) THEN
  CALL crmstyle_scale_rmdi(all_factor, ucu_factor,                        &  
  ucu_h_w, ucu_h_u, ucu_h_v, ucu_h_th, ucu_h_thv, ucu_h_rho,              &
  ucu_h_rh, ucu_h_a, ucu_h_dpx, ucu_h_dpy,                                &
  ucu_h_q,   ucu_h_qcl, ucu_h_qcf, ucu_h_qrain, ucu_h_qgraup,             &
  ucu_h_dt1, ucu_h_dt2, ucu_h_dt4, ucu_h_dt9, ucu_h_dt12, ucu_h_dt30,     &
  ucu_h_dq4, ucu_h_dq9, ucu_h_dq12, ucu_h_dq30,                           &
  ucu_h_dqcl4, ucu_h_dqcl9, ucu_h_dqcl12, ucu_h_dqcl30,                   &
  ucu_h_dqcf4, ucu_h_dqcf3, ucu_h_dqcf12, ucu_h_dqcf30,                   &
  ucu_h_dqrain30, ucu_h_dqgr30, ucu_h_drho, ucu_h_thw, ucu_h_thvw,        &
  ucu_h_qw, ucu_h_qclw, ucu_h_qcfw, ucu_h_qrainw, ucu_h_qgraupw,          &
  ucu_h_uth, ucu_h_vth, ucu_h_uthv, ucu_h_vthv, ucu_h_uq, ucu_h_vq,       &
  ucu_h_wp,                                                               &
  ucu_h_uw, ucu_h_vw, ucu_h_ww, ucu_h_w3, ucu_h_vv, ucu_h_uu, ucu_h_uv )
END IF

!------------------------------------------------------------------------------
! Add Plume counting - here loop over levels and coarse areas
! Note doing this only makes sense if the number of points in the coarse
! area is large i.e. ~100x100 say 
!------------------------------------------------------------------------------
!
! nareas = local_new_x*local_new_y
! Currently setup to do this for BCu only
!
IF (l_bcu) THEN
  DO k=1, mlevs
    DO j= 1,local_new_y
      iy1 = (j-1)*new_res(1)

      DO i= 1,local_new_x
        ix1 = (i-1)*new_res(1)

        ! Only count if there are buoyant plumes
        IF (nn_mask(i,j,k) > 0) THEN
          DO jj=1,new_res(1)
            DO ii=1,new_res(1)
              mask_temp(ii,jj) = bcu_mask(ix1+ii,iy1+jj,k)
            END DO
          END DO
          IF (iprint == 2) THEN
            WRITE(umMessage,'(A,4I6)') ' Before call to count_plumes ',      &
                            k,i,j,nn_mask(i,j,k)
            CALL umPrint(umMessage,src=RoutineName)
          END IF
          CALL count_plumes(new_res(1),new_res(1),iprint,nn_mask(i,j,k),    &
                          mask_temp, nups,                                  &
                          diam_pdf,size_pdf,fract_pdf,dxdy_pdf)

          n_plume(i,j,k)    = REAL(nups)
          plume_size(i,j,k) = REAL(nn_mask(i,j,k))/n_plume(i,j,k)
          plume_diam_pdf(:,i,j,k) = diam_pdf(:)
        ELSE    ! no plumes set to zero
          n_plume(i,j,k)    = 0.0
          plume_size(i,j,k) = 0.0
        END IF

      END DO
    END DO
  END DO       ! k
END IF




!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_sample

END MODULE crmstyle_sample_mod

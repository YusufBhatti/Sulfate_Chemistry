! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  write out fieldsfiles

MODULE crmstyle_write_ff_out_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_WRITE_FF_OUT_MOD'

CONTAINS

SUBROUTINE crmstyle_write_ff_out(it, ntimes)

USE crmstyle_cntl_mod, ONLY:                                             &
  nx_start,ny_start,num_x,num_y, mlevs, l_class_col, l_sect30,           &
  l_bcu, l_wg1, l_acc, l_acu, l_ppd, l_nbd, l_nid, l_adu, l_acw, l_bcw,  &
  l_ucu, l_ppw, l_nbw, l_cape, l_pcape

USE crmwork_arrays_mod, ONLY:                                            &
  h_theta_sea

USE crmstyle_pp_data_mod, ONLY:                                          &
  iyear,imon,iday,ihour,imin,isec,isyear,ismon,isday,ishour,ismin,issec, &
  itim, bzy,bzx,bdy,bdx, new_bzy, new_bzx, new_bdy, new_bdx, origin_x1,  &
  origin_y1,pseudo_lat,pseudo_lon

USE missing_data_mod, ONLY: rmdi, imdi

USE crmstyle_sample_arrays_mod, ONLY:                                        &
  all_zh, all_sh, all_lh,  all_pstar, all_tstar, all_rain,                   &
  all_snow, all_precip, all_orog,  all_land,                                 &
  all_sd_zh, all_sd_sh, all_sd_lh,  all_sd_pstar, all_sd_tstar, all_sd_rain, &
  all_sd_snow, all_sd_precip, all_sd_orog,  all_sd_land,                     &
  fract_conv, fract_strat, prec_conv, prec_strat,                            &
  cape, cin, zlcl, zneutral, zfree,                                          &
  bcu_pcape, bcu_dpcapedt, bcu_dpcapedt_bl, ppd_dpcapedt_bl,                 &
  bcw_pcape, bcw_dpcapedt, bcw_dpcapedt_bl,                                  &
  bcu_dilcape, bcu_dcapedt, bcw_dilcape, bcw_dcapedt,                        &
  all_w, all_u,all_v, all_th, all_thv, all_rho,                              &
  all_q, all_qcl, all_qcf, all_qrain, all_qgraup,                            &
  all_ptheta, all_rh, all_a, all_dpx, all_dpy,                               &
  all_dt1, all_dt2, all_dt4, all_dt9, all_dt12,all_dt30,                     &
  all_dq4, all_dq9, all_dq12, all_dq30,  all_t ,                             &
  all_dqcl4, all_dqcl9, all_dqcl12, all_dqcl30,                              &
  all_dqcf4, all_dqcf3, all_dqcf12, all_dqcf30,                              &
  all_thw, all_thvw,                                                         &
  all_qw, all_qclw, all_qcfw, all_qrainw, all_qgraupw,                       &
  all_uw, all_vw, all_ww, all_w3, all_vv,                                    &
  all_uth, all_vth, all_uthv, all_vthv, all_uq, all_vq, all_wp,              &
  all_uu, all_uv,                                                            &
  all_sd_w, all_sd_u,all_sd_v, all_sd_th, all_sd_t, all_sd_thv, all_sd_rho,  &
  all_sd_q, all_sd_qcl, all_sd_qcf, all_sd_qrain, all_sd_qgraup,             &
  all_sd_ptheta, all_sd_rh,                                                  &
  all_sd_dt1, all_sd_dt2, all_sd_dt4, all_sd_dt9, all_sd_dt12,all_sd_dt30,   &
  all_sd_dq4, all_sd_dq9, all_sd_dq12, all_sd_dq30,                          &
  all_sd_dqcl4, all_sd_dqcl9, all_sd_dqcl12, all_sd_dqcl30,                  &
  all_sd_dqcf4, all_sd_dqcf3, all_sd_dqcf12, all_sd_dqcf30,                  &
  acc_h_w, acc_h_u, acc_h_v, acc_h_th, acc_h_thv, acc_h_rho,                 &
  acc_h_rh, acc_h_a, acc_h_dpx, acc_h_dpy,                                   &
  acc_h_q,   acc_h_qcl, acc_h_qcf, acc_h_qrain, acc_h_qgraup,                &
  acc_h_dt1, acc_h_dt2, acc_h_dt4, acc_h_dt9, acc_h_dt12, acc_h_dt30,        &
  acc_h_dq4, acc_h_dq9, acc_h_dq12, acc_h_dq30,                              &
  acc_h_dqcl4, acc_h_dqcl9, acc_h_dqcl12, acc_h_dqcl30,                      &
  acc_h_dqcf4, acc_h_dqcf3, acc_h_dqcf12, acc_h_dqcf30,                      &
  acc_h_dqrain30, acc_h_dqgr30, acc_h_drho,                                  &
  acc_h_thw, acc_h_thvw, acc_h_wp,                                           &
  acc_h_qw, acc_h_qclw, acc_h_qcfw, acc_h_qrainw, acc_h_qgraupw,             &
  acc_h_uw, acc_h_vw, acc_h_ww, acc_h_w3, acc_h_vv, acc_h_uu, acc_h_uv,      &
  acc_h_uth, acc_h_vth, acc_h_uthv, acc_h_vthv, acc_h_uq, acc_h_vq,          &
  acu_h_w, acu_h_u, acu_h_v, acu_h_th, acu_h_thv, acu_h_rho,                 &
  acu_h_rh, acu_h_a, acu_h_dpx, acu_h_dpy,                                   &
  acu_h_q,   acu_h_qcl, acu_h_qcf, acu_h_qrain, acu_h_qgraup,                &
  acu_h_dt1, acu_h_dt2, acu_h_dt4, acu_h_dt9, acu_h_dt12, acu_h_dt30,        &
  acu_h_dq4, acu_h_dq9, acu_h_dq12, acu_h_dq30,                              &
  acu_h_dqcl4, acu_h_dqcl9, acu_h_dqcl12, acu_h_dqcl30,                      &
  acu_h_dqcf4, acu_h_dqcf3, acu_h_dqcf12, acu_h_dqcf30,                      &
  acu_h_dqrain30, acu_h_dqgr30, acu_h_drho,                                  &
  acu_h_thw, acu_h_thvw, acu_h_wp,                                           &
  acu_h_qw, acu_h_qclw, acu_h_qcfw, acu_h_qrainw, acu_h_qgraupw,             &
  acu_h_uw, acu_h_vw, acu_h_ww, acu_h_w3, acu_h_vv, acu_h_uu, acu_h_uv,      &
  acu_h_uth, acu_h_vth, acu_h_uthv, acu_h_vthv, acu_h_uq, acu_h_vq,          &
  bcu_h_w, bcu_h_u, bcu_h_v, bcu_h_th, bcu_h_thv, bcu_h_rho,                 &
  bcu_h_rh, bcu_h_a, bcu_h_dpx, bcu_h_dpy,                                   &
  bcu_h_q,   bcu_h_qcl, bcu_h_qcf, bcu_h_qrain, bcu_h_qgraup,                &
  bcu_h_dt1, bcu_h_dt2, bcu_h_dt4, bcu_h_dt9, bcu_h_dt12, bcu_h_dt30,        &
  bcu_h_dq4, bcu_h_dq9, bcu_h_dq12, bcu_h_dq30,                              &
  bcu_h_dqcl4, bcu_h_dqcl9, bcu_h_dqcl12, bcu_h_dqcl30,                      &
  bcu_h_dqcf4, bcu_h_dqcf3, bcu_h_dqcf12, bcu_h_dqcf30,                      &
  bcu_h_dqrain30, bcu_h_dqgr30, bcu_h_drho,                                  &
  bcu_h_thw, bcu_h_thvw, bcu_h_wp,                                           &
  bcu_h_qw, bcu_h_qclw, bcu_h_qcfw, bcu_h_qrainw, bcu_h_qgraupw,             &
  bcu_h_uw, bcu_h_vw, bcu_h_ww, bcu_h_w3, bcu_h_vv, bcu_h_uu, bcu_h_uv,      &
  bcu_h_uth, bcu_h_vth, bcu_h_uthv, bcu_h_vthv, bcu_h_uq, bcu_h_vq,          &
  wg1_h_w, wg1_h_u, wg1_h_v, wg1_h_th, wg1_h_thv, wg1_h_rho,                 &
  wg1_h_rh, wg1_h_a, wg1_h_dpx, wg1_h_dpy,                                   &
  wg1_h_q,   wg1_h_qcl, wg1_h_qcf, wg1_h_qrain, wg1_h_qgraup,                &
  wg1_h_dt1, wg1_h_dt2, wg1_h_dt4, wg1_h_dt9, wg1_h_dt12, wg1_h_dt30,        &
  wg1_h_dq4, wg1_h_dq9, wg1_h_dq12, wg1_h_dq30,                              &
  wg1_h_dqcl4, wg1_h_dqcl9, wg1_h_dqcl12, wg1_h_dqcl30,                      &
  wg1_h_dqcf4, wg1_h_dqcf3, wg1_h_dqcf12, wg1_h_dqcf30,                      &
  wg1_h_dqrain30, wg1_h_dqgr30, wg1_h_drho,                                  &
  wg1_h_thw, wg1_h_thvw, wg1_h_wp,                                           &
  wg1_h_qw, wg1_h_qclw, wg1_h_qcfw, wg1_h_qrainw, wg1_h_qgraupw,             &
  wg1_h_uw, wg1_h_vw, wg1_h_ww, wg1_h_w3, wg1_h_vv, wg1_h_uu, wg1_h_uv,      &
  wg1_h_uth, wg1_h_vth, wg1_h_uthv, wg1_h_vthv, wg1_h_uq, wg1_h_vq,          &
  ppd_h_w, ppd_h_u, ppd_h_v, ppd_h_th, ppd_h_thv, ppd_h_rho,                 &
  ppd_h_rh, ppd_h_a, ppd_h_dpx, ppd_h_dpy,                                   &
  ppd_h_q,   ppd_h_qcl, ppd_h_qcf, ppd_h_qrain, ppd_h_qgraup,                &
  ppd_h_dt1, ppd_h_dt2, ppd_h_dt4, ppd_h_dt9, ppd_h_dt12, ppd_h_dt30,        &
  ppd_h_dq4, ppd_h_dq9, ppd_h_dq12, ppd_h_dq30,                              &
  ppd_h_dqcl4, ppd_h_dqcl9, ppd_h_dqcl12, ppd_h_dqcl30,                      &
  ppd_h_dqcf4, ppd_h_dqcf3, ppd_h_dqcf12, ppd_h_dqcf30,                      &
  ppd_h_dqrain30, ppd_h_dqgr30, ppd_h_drho,                                  &
  ppd_h_thw, ppd_h_thvw, ppd_h_wp,                                           &
  ppd_h_qw, ppd_h_qclw, ppd_h_qcfw, ppd_h_qrainw, ppd_h_qgraupw,             &
  ppd_h_uw, ppd_h_vw, ppd_h_ww, ppd_h_w3, ppd_h_vv, ppd_h_uu, ppd_h_uv,      &
  ppd_h_uth, ppd_h_vth, ppd_h_uthv, ppd_h_vthv, ppd_h_uq, ppd_h_vq,          &
  nbd_h_w, nbd_h_u, nbd_h_v, nbd_h_th, nbd_h_thv, nbd_h_rho,                 &
  nbd_h_rh, nbd_h_a, nbd_h_dpx, nbd_h_dpy,                                   &
  nbd_h_q,   nbd_h_qcl, nbd_h_qcf, nbd_h_qrain, nbd_h_qgraup,                &
  nbd_h_dt1, nbd_h_dt2, nbd_h_dt4, nbd_h_dt9, nbd_h_dt12, nbd_h_dt30,        &
  nbd_h_dq4, nbd_h_dq9, nbd_h_dq12, nbd_h_dq30,                              &
  nbd_h_dqcl4, nbd_h_dqcl9, nbd_h_dqcl12, nbd_h_dqcl30,                      &
  nbd_h_dqcf4, nbd_h_dqcf3, nbd_h_dqcf12, nbd_h_dqcf30,                      &
  nbd_h_dqrain30, nbd_h_dqgr30, nbd_h_drho,                                  &
  nbd_h_thw, nbd_h_thvw, nbd_h_wp,                                           &
  nbd_h_qw, nbd_h_qclw, nbd_h_qcfw, nbd_h_qrainw, nbd_h_qgraupw,             &
  nbd_h_uw, nbd_h_vw, nbd_h_ww, nbd_h_w3, nbd_h_vv, nbd_h_uu, nbd_h_uv,      &
  nbd_h_uth, nbd_h_vth, nbd_h_uthv, nbd_h_vthv, nbd_h_uq, nbd_h_vq,          &
  nid_h_w, nid_h_u, nid_h_v, nid_h_th, nid_h_thv, nid_h_rho,                 &
  nid_h_rh, nid_h_a, nid_h_dpx, nid_h_dpy,                                   &
  nid_h_q,   nid_h_qcl, nid_h_qcf, nid_h_qrain, nid_h_qgraup,                &
  nid_h_dt1, nid_h_dt2, nid_h_dt4, nid_h_dt9, nid_h_dt12, nid_h_dt30,        &
  nid_h_dq4, nid_h_dq9, nid_h_dq12, nid_h_dq30,                              &
  nid_h_dqcl4, nid_h_dqcl9, nid_h_dqcl12, nid_h_dqcl30,                      &
  nid_h_dqcf4, nid_h_dqcf3, nid_h_dqcf12, nid_h_dqcf30,                      &
  nid_h_dqrain30, nid_h_dqgr30, nid_h_drho,                                  &
  nid_h_thw, nid_h_thvw, nid_h_wp,                                           &
  nid_h_qw, nid_h_qclw, nid_h_qcfw, nid_h_qrainw, nid_h_qgraupw,             &
  nid_h_uw, nid_h_vw, nid_h_ww, nid_h_w3, nid_h_vv, nid_h_uu, nid_h_uv,      &
  nid_h_uth, nid_h_vth, nid_h_uthv, nid_h_vthv, nid_h_uq, nid_h_vq,          &
  adu_h_w, adu_h_u, adu_h_v, adu_h_th, adu_h_thv, adu_h_rho,                 &
  adu_h_rh, adu_h_a, adu_h_dpx, adu_h_dpy,                                   &
  adu_h_q,   adu_h_qcl, adu_h_qcf, adu_h_qrain, adu_h_qgraup,                &
  adu_h_dt1, adu_h_dt2, adu_h_dt4, adu_h_dt9, adu_h_dt12, adu_h_dt30,        &
  adu_h_dq4, adu_h_dq9, adu_h_dq12, adu_h_dq30,                              &
  adu_h_dqcl4, adu_h_dqcl9, adu_h_dqcl12, adu_h_dqcl30,                      &
  adu_h_dqcf4, adu_h_dqcf3, adu_h_dqcf12, adu_h_dqcf30,                      &
  adu_h_dqrain30, adu_h_dqgr30, adu_h_drho,                                  &
  adu_h_thw, adu_h_thvw, adu_h_wp,                                           &
  adu_h_qw, adu_h_qclw, adu_h_qcfw, adu_h_qrainw, adu_h_qgraupw,             &
  adu_h_uw, adu_h_vw, adu_h_ww, adu_h_w3, adu_h_vv, adu_h_uu, adu_h_uv,      &
  adu_h_uth, adu_h_vth, adu_h_uthv, adu_h_vthv, adu_h_uq, adu_h_vq,          &
  ucu_h_w, ucu_h_u, ucu_h_v, ucu_h_th, ucu_h_thv, ucu_h_rho,                 &
  ucu_h_rh, ucu_h_a, ucu_h_dpx, ucu_h_dpy,                                   &
  ucu_h_q,   ucu_h_qcl, ucu_h_qcf, ucu_h_qrain, ucu_h_qgraup,                &
  ucu_h_dt1, ucu_h_dt2, ucu_h_dt4, ucu_h_dt9, ucu_h_dt12, ucu_h_dt30,        &
  ucu_h_dq4, ucu_h_dq9, ucu_h_dq12, ucu_h_dq30,                              &
  ucu_h_dqcl4, ucu_h_dqcl9, ucu_h_dqcl12, ucu_h_dqcl30,                      &
  ucu_h_dqcf4, ucu_h_dqcf3, ucu_h_dqcf12, ucu_h_dqcf30,                      &
  ucu_h_dqrain30, ucu_h_dqgr30, ucu_h_drho,                                  &
  ucu_h_thw, ucu_h_thvw,  ucu_h_wp,                                          &
  ucu_h_qw, ucu_h_qclw, ucu_h_qcfw, ucu_h_qrainw, ucu_h_qgraupw,             &
  ucu_h_uw, ucu_h_vw, ucu_h_ww, ucu_h_w3, ucu_h_vv, ucu_h_uu, ucu_h_uv,      &
  ucu_h_uth, ucu_h_vth, ucu_h_uthv, ucu_h_vthv, ucu_h_uq, ucu_h_vq,          &
  acw_h_w, acw_h_u, acw_h_v, acw_h_th, acw_h_thv, acw_h_rho,                 &
  acw_h_rh, acw_h_a, acw_h_dpx, acw_h_dpy,                                   &
  acw_h_q,   acw_h_qcl, acw_h_qcf, acw_h_qrain, acw_h_qgraup,                &
  acw_h_dt1, acw_h_dt2, acw_h_dt4, acw_h_dt9, acw_h_dt12, acw_h_dt30,        &
  acw_h_dq4, acw_h_dq9, acw_h_dq12, acw_h_dq30,                              &
  acw_h_dqcl4, acw_h_dqcl9, acw_h_dqcl12, acw_h_dqcl30,                      &
  acw_h_dqcf4, acw_h_dqcf3, acw_h_dqcf12, acw_h_dqcf30,                      &
  acw_h_dqrain30, acw_h_dqgr30, acw_h_drho,                                  &
  acw_h_thw, acw_h_thvw, acw_h_wp,                                           &
  acw_h_qw, acw_h_qclw, acw_h_qcfw, acw_h_qrainw, acw_h_qgraupw,             &
  acw_h_uw, acw_h_vw, acw_h_ww, acw_h_w3, acw_h_vv, acw_h_uu, acw_h_uv,      &
  acw_h_uth, acw_h_vth, acw_h_uthv, acw_h_vthv, acw_h_uq, acw_h_vq,          &
  bcw_h_w, bcw_h_u, bcw_h_v, bcw_h_th, bcw_h_thv, bcw_h_rho,                 &
  bcw_h_rh, bcw_h_a, bcw_h_dpx, bcw_h_dpy,                                   &
  bcw_h_q,   bcw_h_qcl, bcw_h_qcf, bcw_h_qrain, bcw_h_qgraup,                &
  bcw_h_dt1, bcw_h_dt2, bcw_h_dt4, bcw_h_dt9, bcw_h_dt12, bcw_h_dt30,        &
  bcw_h_dq4, bcw_h_dq9, bcw_h_dq12, bcw_h_dq30,                              &
  bcw_h_dqcl4, bcw_h_dqcl9, bcw_h_dqcl12, bcw_h_dqcl30,                      &
  bcw_h_dqcf4, bcw_h_dqcf3, bcw_h_dqcf12, bcw_h_dqcf30,                      &
  bcw_h_dqrain30, bcw_h_dqgr30, bcw_h_drho,                                  &
  bcw_h_thw, bcw_h_thvw, bcw_h_wp,                                           &
  bcw_h_qw, bcw_h_qclw, bcw_h_qcfw, bcw_h_qrainw, bcw_h_qgraupw,             &
  bcw_h_uw, bcw_h_vw, bcw_h_ww, bcw_h_w3, bcw_h_vv, bcw_h_uu, bcw_h_uv,      &
  bcw_h_uth, bcw_h_vth, bcw_h_uthv, bcw_h_vthv, bcw_h_uq, bcw_h_vq,          &
  ppw_h_w, ppw_h_u, ppw_h_v, ppw_h_th, ppw_h_thv, ppw_h_rho,                 &
  ppw_h_rh, ppw_h_a, ppw_h_dpx, ppw_h_dpy,                                   &
  ppw_h_q,   ppw_h_qcl, ppw_h_qcf, ppw_h_qrain, ppw_h_qgraup,                &
  ppw_h_dt1, ppw_h_dt2, ppw_h_dt4, ppw_h_dt9, ppw_h_dt12, ppw_h_dt30,        &
  ppw_h_dq4, ppw_h_dq9, ppw_h_dq12, ppw_h_dq30,                              &
  ppw_h_dqcl4, ppw_h_dqcl9, ppw_h_dqcl12, ppw_h_dqcl30,                      &
  ppw_h_dqcf4, ppw_h_dqcf3, ppw_h_dqcf12, ppw_h_dqcf30,                      &
  ppw_h_dqrain30, ppw_h_dqgr30, ppw_h_drho,                                  &
  ppw_h_thw, ppw_h_thvw, ppw_h_wp,                                           &
  ppw_h_qw, ppw_h_qclw, ppw_h_qcfw, ppw_h_qrainw, ppw_h_qgraupw,             &
  ppw_h_uw, ppw_h_vw, ppw_h_ww, ppw_h_w3, ppw_h_vv, ppw_h_uu, ppw_h_uv,      &
  ppw_h_uth, ppw_h_vth, ppw_h_uthv, ppw_h_vthv, ppw_h_uq, ppw_h_vq,          &
  nbw_h_w, nbw_h_u, nbw_h_v, nbw_h_th, nbw_h_thv, nbw_h_rho,                 &
  nbw_h_rh, nbw_h_a, nbw_h_dpx, nbw_h_dpy,                                   &
  nbw_h_q,   nbw_h_qcl, nbw_h_qcf, nbw_h_qrain, nbw_h_qgraup,                &
  nbw_h_dt1, nbw_h_dt2, nbw_h_dt4, nbw_h_dt9, nbw_h_dt12, nbw_h_dt30,        &
  nbw_h_dq4, nbw_h_dq9, nbw_h_dq12, nbw_h_dq30,                              &
  nbw_h_dqcl4, nbw_h_dqcl9, nbw_h_dqcl12, nbw_h_dqcl30,                      &
  nbw_h_dqcf4, nbw_h_dqcf3, nbw_h_dqcf12, nbw_h_dqcf30,                      &
  nbw_h_dqrain30, nbw_h_dqgr30, nbw_h_drho,                                  &
  nbw_h_thw, nbw_h_thvw, nbw_h_wp,                                           &
  nbw_h_qw, nbw_h_qclw, nbw_h_qcfw, nbw_h_qrainw, nbw_h_qgraupw,             &
  nbw_h_uw, nbw_h_vw, nbw_h_ww, nbw_h_w3, nbw_h_vv, nbw_h_uu, nbw_h_uv,      &
  nbw_h_uth, nbw_h_vth, nbw_h_uthv, nbw_h_vthv, nbw_h_uq, nbw_h_vq,          &
  bcu_mask_w,                                                                &
  n_plume, plume_size, plume_diam_pdf

USE crmstyle_output_hdr_mod, ONLY:        &
  all_hdr                & ! UM Headers:  all_file
 ,acc_hdr                & ! UM Headers:  acc_file
 ,acu_hdr                & ! UM Headers:  acu_file
 ,bcu_hdr                & ! UM Headers:  bcu_file
 ,wg1_hdr                & ! UM Headers:  wg1_file
 ,ppd_hdr                & ! UM Headers:  ppd_file
 ,nbd_hdr                & ! UM Headers:  nbd_file
 ,nid_hdr                & ! UM Headers:  nid_file
 ,adu_hdr                & ! UM Headers:  adu_file
 ,acw_hdr                & ! UM Headers:  acw_file
 ,bcw_hdr                & ! UM Headers:  bcw_file
 ,ucu_hdr                & ! UM Headers:  ucu_file
 ,ppw_hdr                & ! UM Headers:  ppw_file
 ,nbw_hdr                & ! UM Headers:  nbw_file
 ,single_hdr               ! UM Headers:  single_file

USE io

USE file_manager, ONLY: assign_file_unit

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type,         &
  lenfixhd

USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal,            &
  EndofFile

! subroutines

USE crmstyle_write_all_ff_mod,    ONLY: crmstyle_write_all_ff
USE crmstyle_write_crm_ff_mod,    ONLY: crmstyle_write_crm_ff
USE crmstyle_write_plume_ff_mod,  ONLY: crmstyle_write_plume_ff
USE crmstyle_write_single_ff_mod, ONLY: crmstyle_write_single_ff
USE ereport_mod, ONLY: ereport, ereport_finalise

USE get_env_var_mod, ONLY: get_env_var

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE packing_codes_mod, ONLY: PC_No_Packing

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Read Land sea mask ancillary file
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.5.

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER,INTENT(IN) ::    &
  it                     & ! call number
 ,ntimes                   ! Total number of times


!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, k, ij,ii,jj                  ! loop counters

INTEGER ::               &
  icall_type             & ! type of call
 ,icall_type2              ! type of call
INTEGER ::               &
  ErrorStatus            & ! Error code from operations on file
! ErrorStatus is declared and passed, but never checked. Should it be deleted?
 ,MaxFldsout               ! max output fields


CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_WRITE_FF_OUT"

TYPE(UM_Header_type)    :: example_hdr   ! Header to copy

TYPE(PP_Field_type) ::       &
  ref_modlev_field(mlevs)      ! pp fields on model levels
TYPE(PP_Field_type) ::       &
  ref_single_field(1)          ! pp fields - single level

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

! Set up reference pp header for output fields
! All output fields on same p grid (i.e. no uv grids etc)
! take date from input pp fields

! Setup output field headers for model level fields correctly

! integer header
ref_modlev_field(:) % Hdr % validyear  = iyear    ! take from input fields
ref_modlev_field(:) % Hdr % validmonth = imon
ref_modlev_field(:) % Hdr % validdate  = iday
ref_modlev_field(:) % Hdr % validhour  = ihour
ref_modlev_field(:) % Hdr % validmin   = imin
ref_modlev_field(:) % Hdr % validsec   = isec
ref_modlev_field(:) % Hdr % datayear  = isyear
ref_modlev_field(:) % Hdr % datamonth = ismon
ref_modlev_field(:) % Hdr % datadate  = isday
ref_modlev_field(:) % Hdr % datahour  = ishour
ref_modlev_field(:) % Hdr % datamin   = ismin
ref_modlev_field(:) % Hdr % datasec   = issec

! Now takes value from input so correct calendar info
ref_modlev_field(:) % Hdr % lbtim     = itim   ! forecast & validity time
ref_modlev_field(:) % Hdr % fcrange   = 0      ! hours from input fields

ref_modlev_field(:) % Hdr % LBLRec     = num_y*num_x

ref_modlev_field(:) % Hdr % lbcode     = 101     ! rotated pole
ref_modlev_field(:) % Hdr % lbhem      = 3       ! limited area
ref_modlev_field(:) % Hdr % NumRows    = num_y
ref_modlev_field(:) % Hdr % NumCols    = num_x
ref_modlev_field(:) % Hdr % lbext      = 0       ! no extra data
ref_modlev_field(:) % Hdr % lbpack     = PC_No_Packing       ! unpacked output
ref_modlev_field(:) % Hdr % lbrel      = 3
ref_modlev_field(:) % Hdr % ppcode     = 0       ! none
ref_modlev_field(:) % Hdr % lbcfc      = 0       ! none
ref_modlev_field(:) % Hdr % lbproc     = 0

ref_modlev_field(:) % Hdr % lbvc       = 1       ! height
ref_modlev_field(:) % Hdr % lbrvc      = 0       ! none
ref_modlev_field(:) % Hdr % lbexp      = 0       ! none
ref_modlev_field(:) % Hdr % MO8proj    = 0       ! none
ref_modlev_field(:) % Hdr % MO8type    = 0       ! none
DO k=1,mlevs
  ref_modlev_field(k) % Hdr % MO8Level   = k       ! level number
END DO
DO k=1,4
  ref_modlev_field(:) % Hdr % LBRsvd(k) = 0    ! No values
END DO
ref_modlev_field(:) % Hdr % stcode   = 0       ! need to set for each fiels
ref_modlev_field(:) % Hdr % LBSrce   = 1111    !
ref_modlev_field(:) % Hdr % LBUser1  = 1       ! reals
ref_modlev_field(:) % Hdr % LBUser3  = 0       !
ref_modlev_field(:) % Hdr % LBUser5  = 0       !
ref_modlev_field(:) % Hdr % LBUser6  = 0       !
ref_modlev_field(:) % Hdr % LBUser7  = 1       ! atmosphere model

! real part
DO k=1,mlevs
  ref_modlev_field(k) % Hdr %RLevel   = h_theta_sea(k)  ! height of level
END DO
ref_modlev_field(:) % Hdr % pseudolat  = pseudo_lat   ! from grid
ref_modlev_field(:) % Hdr % pseudolon  = pseudo_lon
ref_modlev_field(:) % Hdr % bgor       = 0.0
ref_modlev_field(:) % Hdr % ZerothLat  = new_bzy   ! Zeroth, so deduct dy
ref_modlev_field(:) % Hdr % LatInt     = new_bdy
ref_modlev_field(:) % Hdr % ZerothLon  = new_bzx
ref_modlev_field(:) % Hdr % LonInt     = new_bdx
ref_modlev_field(:) % Hdr % bmdi       = rmdi
ref_modlev_field(:) % Hdr % bmks       = 1.0

! allocate space for data
DO k=1,mlevs
  ALLOCATE ( ref_modlev_field(k) % RData( num_y*num_x,1 ) )
END DO

!------------------------------------------------------------------------------
! Set up reference pp header for single fields
! take date from input pp fields
!------------------------------------------------------------------------------

! integer header
ref_single_field(1) % Hdr % validyear  = iyear    ! take from input fields
ref_single_field(1) % Hdr % validmonth = imon
ref_single_field(1) % Hdr % validdate  = iday
ref_single_field(1) % Hdr % validhour  = ihour
ref_single_field(1) % Hdr % validmin   = imin
ref_single_field(1) % Hdr % validsec   = isec
ref_single_field(1) % Hdr % datayear  = isyear
ref_single_field(1) % Hdr % datamonth = ismon
ref_single_field(1) % Hdr % datadate  = isday
ref_single_field(1) % Hdr % datahour  = ishour
ref_single_field(1) % Hdr % datamin   = ismin
ref_single_field(1) % Hdr % datasec   = issec
! Now takes value from input so correct calendar info
ref_single_field(1) % Hdr % lbtim     = itim   ! forecast & validity time
ref_single_field(1) % Hdr % fcrange   = 0      ! hours from input fields

ref_single_field(1) % Hdr % LBLRec     = num_y*num_x

ref_single_field(1) % Hdr % lbcode     = 101     ! rotated pole
ref_single_field(1) % Hdr % lbhem      = 3       ! limited area
ref_single_field(1) % Hdr % NumRows    = num_y
ref_single_field(1) % Hdr % NumCols    = num_x
ref_single_field(1) % Hdr % lbext      = 0       ! no extra data
ref_single_field(1) % Hdr % lbpack     = PC_No_Packing       ! unpacked output
ref_single_field(1) % Hdr % lbrel      = 3
ref_single_field(1) % Hdr % ppcode     = 0       ! none
ref_single_field(1) % Hdr % lbcfc      = 0       ! none
ref_single_field(1) % Hdr % lbproc     = 0
ref_single_field(1) % Hdr % lbexp      = 0       ! none
ref_single_field(1) % Hdr % MO8proj    = 0       ! none
ref_single_field(1) % Hdr % MO8type    = 0       ! none

ref_single_field(1) % Hdr % lbvc       = 129     ! surface
ref_single_field(1) % Hdr % lbrvc      = 0       ! none
ref_single_field(1) % Hdr % LBRsvd(:)  = 0       ! No values
ref_single_field(1) % Hdr % stcode     = 0       ! need to set for each fiels
ref_single_field(1) % Hdr % LBSrce     = 1111    !
ref_single_field(1) % Hdr % LBUser1    = 1       !
ref_single_field(1) % Hdr % LBUser3    = 0       !
ref_single_field(1) % Hdr % LBUser5    = 0       !
ref_single_field(1) % Hdr % LBUser6    = 0       !
ref_single_field(1) % Hdr % LBUser7    = 1       ! atmosphere model

ref_single_field(1) % Hdr % pseudolat  = pseudo_lat  ! from input fields
ref_single_field(1) % Hdr % pseudolon  = pseudo_lon
ref_single_field(1) % Hdr % bgor       = 0.0
ref_single_field(1) % Hdr % ZerothLat  = new_bzy   ! Zeroth
ref_single_field(1) % Hdr % LatInt     = new_bdy
ref_single_field(1) % Hdr % ZerothLon  = new_bzx
ref_single_field(1) % Hdr % LonInt     = new_bdx
ref_single_field(1) % Hdr % bmdi       = rmdi
ref_single_field(1) % Hdr % bmks       = 1.0

! allocate space for data

ALLOCATE ( ref_single_field(1) % RData( num_y*num_x,1 ) )

!=============================================================================
! Setup values for Fields File header
!=============================================================================
ALLOCATE (example_hdr %FixHd(LenFixHd))  ! fix header
example_hdr %FixHd(:) = imdi             ! intialise

! Need to be set ?
! Lengths of integer and real headers seem to need to correspond to pp headers
example_hdr % LenIntC  = 46    ! for a FF ? ( 15 for ancillary)
example_hdr % LenRealC = 38    ! for a FF ? (6 ancillary)
example_hdr % Len1LevDepC = mlevs
example_hdr % Len2LevDepC = 8
example_hdr % Len1RowDepC = 0
example_hdr % Len2RowDepC = 0
example_hdr % Len1ColDepC = 0
example_hdr % Len2ColDepC = 0
example_hdr % Len1FldsOfC = 0
example_hdr % Len2FldsOfC = 0
example_hdr % LenExtraC   = 0
example_hdr % LenHistFile = 0
example_hdr % LenCompFldI1 = 0
example_hdr % LenCompFldI2 = 0
example_hdr % LenCompFldI3 = imdi
example_hdr % Len1Lookup = 64
example_hdr % Len2Lookup = 2    ! small lookup table as used as example

! Note need to allocate these arrays or problems in init_pp called by new_hdr
ALLOCATE (example_hdr % RowDepC(1))  !
example_hdr % RowDepC    = 0.0
ALLOCATE (example_hdr % ColDepC(1))  !
example_hdr % ColDepC    = 0.0

!Positions set by a call to init_pp from new_hdr
example_hdr % FixHd(1)   = 20   !  version
example_hdr % FixHd(2)   = 1    ! Atmosphere
example_hdr % FixHd(3)   = 1    ! hybrid vertical coordinate type
example_hdr % FixHd(4)   = 3    ! LAM no wrap

example_hdr % FixHd(5)   = 3     ! fields file (example come from ancillary)
example_hdr % FixHd(6)   = 0     ! no run identifier
example_hdr % FixHd(8)   = 1     ! gregorian calendar
example_hdr % FixHd(9)   = 3     ! Arakawa C grid
example_hdr % FixHd(10)   = 0     ! gregorian calendar
example_hdr % FixHd(12)   = 803     ! UM version 803

example_hdr % FixHd(21)   = iyear    ! First validity time in FF
example_hdr % FixHd(22)   = imon     !
example_hdr % FixHd(23)   = iday     !
example_hdr % FixHd(24)   = ihour    !
example_hdr % FixHd(25)   = imin     !
example_hdr % FixHd(26)   = isec     !
example_hdr % FixHd(27)   = 0        ! Should be day number


example_hdr % FixHd(28)   = iyear    ! last validity time in FF
example_hdr % FixHd(29)   = imon     !
example_hdr % FixHd(30)   = iday     !
example_hdr % FixHd(31)   = ihour    !
example_hdr % FixHd(32)   = imin     !
example_hdr % FixHd(33)   = isec     !
example_hdr % FixHd(34)   = 0        ! Should be day number

! Assume always the same
example_hdr % FixHd(100) = 256+1    ! start of integer constants
example_hdr % FixHd(101) = 46       ! for a FF  ( 15 for ancillary)

example_hdr % FixHd(112) = imdi     !
example_hdr % FixHd(115) = imdi     ! No row-dependent constants
example_hdr % FixHd(116) = imdi     !
example_hdr % FixHd(117) = imdi     !
example_hdr % FixHd(120) = imdi     ! No column-dependent constants
example_hdr % FixHd(121) = imdi     !
example_hdr % FixHd(122) = imdi     !
example_hdr % FixHd(125:149)=  imdi     ! Not using

example_hdr % FixHd(151) = 64       ! 1st dimension lookup table
example_hdr % FixHd(152) = 2        ! 2nd dimension lookup table

ALLOCATE (example_hdr %IntC(example_hdr % LenIntC))
example_hdr %IntC(:)     = imdi        ! initialise
example_hdr % IntC(6)    = num_x       ! Number of points E-W
example_hdr % IntC(7)    = num_y       ! Number of points N-S
example_hdr % IntC(8)    = mlevs
example_hdr % IntC(9)    = mlevs


ALLOCATE (example_hdr %RealC(example_hdr % LenRealC))
example_hdr %RealC(:) = rmdi          ! initialise

example_hdr % RealC(1)   = new_bdx       ! E-W grid spacing in degrees
example_hdr % RealC(2)   = new_bdy       ! N-S grid spacing in degrees
example_hdr % RealC(3)   = origin_y1     ! Latitude of first position
example_hdr % RealC(4)   = origin_x1     ! Longitude of first position
example_hdr % RealC(5)   = pseudo_lat    ! From original grid
example_hdr % RealC(6)   = pseudo_lon    ! From original grid

ALLOCATE (example_hdr %LevDepC(8*mlevs+1))
example_hdr %LevDepC(:) = rmdi          ! initialise
!=============================================================================

!-----------------------------------------------------------------------
! work out type of call

IF (it == 1) THEN
  icall_type = 1    ! first call need to open and set up header
ELSE IF (it == ntimes) THEN
  icall_type = 3    ! Final call need to close file and write header
ELSE
  icall_type = 2    ! No need to open or close file
END IF
!
!-----------------------------------------------------------------------
! Coarse grid means
!-----------------------------------------------------------------------
IF (icall_type == 1) THEN
  CALL get_env_var("ALL_FILE", all_hdr % FileName)
  
  CALL assign_file_unit(all_hdr % FileName, all_hdr % UnitNum, handler="portio")
END IF

IF (l_sect30) THEN ! 12 more fields in all stream 6 more in other streams
  maxfldsout = mlevs*99*ntimes  ! more than enough
ELSE               !  less fields
  maxfldsout = mlevs*87*ntimes  ! more than enough
END IF

CALL  crmstyle_write_all_ff(icall_type, maxfldsout, example_hdr,           &
                            ref_modlev_field,                              &
                            all_hdr, errorstatus )

!-----------------------------------------------------------------------
! Coarse grid means 
!-----------------------------------------------------------------------
IF (icall_type == 1) THEN
  CALL get_env_var("SINGLE_FILE", single_hdr % FileName)

  CALL assign_file_unit(single_hdr % FileName, single_hdr % UnitNum, &
      handler="portio")
END IF

maxfldsout = 100*ntimes  ! more than enough

IF (l_cape .OR. l_cape) THEN
  maxfldsout = 110*ntimes  ! make sure more space
END IF

CALL  crmstyle_write_single_ff(icall_type, maxfldsout, example_hdr,       &
  ref_single_field,                                                       &
  all_zh, all_sh, all_lh,  all_pstar, all_tstar, all_rain,                &
  all_snow, all_precip, all_orog,  all_land,                              &
  all_sd_zh, all_sd_sh, all_sd_lh,  all_sd_pstar,all_sd_tstar,all_sd_rain,&
  all_sd_snow, all_sd_precip, all_sd_orog,  all_sd_land,                  &
  fract_conv, fract_strat, prec_conv, prec_strat, cape, cin, zlcl, zfree, &
  zneutral,  bcu_pcape, bcu_dpcapedt, bcu_dpcapedt_bl,                    &
  bcw_pcape, bcw_dpcapedt, bcw_dpcapedt_bl, ppd_dpcapedt_bl,              &
  bcu_dilcape, bcu_dcapedt, bcw_dilcape, bcw_dcapedt,                     &
  single_hdr, errorstatus     )
  
!-----------------------------------------------------------------------
! All partition files - same set of fields in most. Output now controlled
! by switches so that don't have to have all partitions.
!-----------------------------------------------------------------------

IF (l_sect30) THEN
  maxfldsout = mlevs*59*ntimes  ! more than enough (~57 fields at present)
ELSE
  maxfldsout = mlevs*57*ntimes  ! more than enough (~51 fields at present)
END IF


!-----------------------------------------------------------------------
! ACC file
!-----------------------------------------------------------------------
IF (l_acc) THEN

  IF (icall_type == 1) THEN
    CALL get_env_var("ACC_FILE", acc_hdr % FileName)

    CALL assign_file_unit(acc_hdr % FileName, acc_hdr % UnitNum, &
                          handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  acc_h_w, acc_h_u, acc_h_v, acc_h_th, acc_h_thv, acc_h_rho,                &
  acc_h_rh, acc_h_a,  acc_h_dpx, acc_h_dpy,                                 &
  acc_h_q,   acc_h_qcl, acc_h_qcf, acc_h_qrain, acc_h_qgraup,               &
  acc_h_dt1, acc_h_dt2, acc_h_dt4, acc_h_dt9, acc_h_dt12, acc_h_dt30,       &
  acc_h_dq4, acc_h_dq9, acc_h_dq12, acc_h_dq30,                             &
  acc_h_dqcl4, acc_h_dqcl9, acc_h_dqcl12,  acc_h_dqcl30,                    &
  acc_h_dqcf4, acc_h_dqcf3, acc_h_dqcf12,  acc_h_dqcf30,                    &
  acc_h_dqrain30, acc_h_dqgr30, acc_h_drho,                                 &
  acc_h_thw, acc_h_thvw, acc_h_qw, acc_h_qclw, acc_h_qcfw, acc_h_qrainw,    &
  acc_h_qgraupw, acc_h_uw, acc_h_vw, acc_h_ww, acc_h_w3, acc_h_vv,          &
  acc_h_uu, acc_h_uv,                                                       &
  acc_h_uth, acc_h_vth, acc_h_uthv, acc_h_vthv, acc_h_uq, acc_h_vq,acc_h_wp,&
  acc_hdr, errorstatus )

END IF  
!-----------------------------------------------------------------------
! ACU file
!-----------------------------------------------------------------------
IF (l_acu) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("ACU_FILE", acu_hdr % FileName)

    CALL assign_file_unit(acu_hdr % FileName, acu_hdr % UnitNum, &
                          handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  acu_h_w, acu_h_u, acu_h_v, acu_h_th, acu_h_thv, acu_h_rho,                &
  acu_h_rh, acu_h_a,  acu_h_dpx, acu_h_dpy,                                 &
  acu_h_q,   acu_h_qcl, acu_h_qcf, acu_h_qrain, acu_h_qgraup,               &
  acu_h_dt1, acu_h_dt2, acu_h_dt4, acu_h_dt9, acu_h_dt12, acu_h_dt30,       &
  acu_h_dq4, acu_h_dq9, acu_h_dq12, acu_h_dq30,                             &
  acu_h_dqcl4, acu_h_dqcl9, acu_h_dqcl12,  acu_h_dqcl30,                    &
  acu_h_dqcf4, acu_h_dqcf3, acu_h_dqcf12,  acu_h_dqcf30,                    &
  acu_h_dqrain30, acu_h_dqgr30, acu_h_drho,                                 &
  acu_h_thw, acu_h_thvw, acu_h_qw, acu_h_qclw, acu_h_qcfw, acu_h_qrainw,    &
  acu_h_qgraupw, acu_h_uw, acu_h_vw, acu_h_ww, acu_h_w3, acu_h_vv,          &
  acu_h_uu, acu_h_uv,                                                       &
  acu_h_uth, acu_h_vth, acu_h_uthv, acu_h_vthv, acu_h_uq, acu_h_vq,acu_h_wp,&
  acu_hdr, errorstatus )

END IF

!-----------------------------------------------------------------------
! BCU file
!-----------------------------------------------------------------------
IF (l_bcu) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("BCU_FILE", bcu_hdr % FileName)

    CALL assign_file_unit(bcu_hdr % FileName, bcu_hdr % UnitNum, &
                          handler="portio")
  END IF

  IF (icall_type == 1 .OR. icall_type == 2) THEN
    icall_type2 = icall_type
  ELSE
    icall_type2 = 2     ! don't close file as still to write plume info
  END IF

  CALL  crmstyle_write_crm_ff(icall_type2,maxfldsout, example_hdr,          &
  ref_modlev_field,                                                         &
  bcu_h_w, bcu_h_u, bcu_h_v, bcu_h_th, bcu_h_thv, bcu_h_rho,                &
  bcu_h_rh, bcu_h_a,  bcu_h_dpx, bcu_h_dpy,                                 &
  bcu_h_q,   bcu_h_qcl, bcu_h_qcf, bcu_h_qrain, bcu_h_qgraup,               &
  bcu_h_dt1, bcu_h_dt2, bcu_h_dt4, bcu_h_dt9, bcu_h_dt12, bcu_h_dt30,       &
  bcu_h_dq4, bcu_h_dq9, bcu_h_dq12, bcu_h_dq30,                             &
  bcu_h_dqcl4, bcu_h_dqcl9, bcu_h_dqcl12,  bcu_h_dqcl30,                    &
  bcu_h_dqcf4, bcu_h_dqcf3, bcu_h_dqcf12,  bcu_h_dqcf30,                    &
  bcu_h_dqrain30, bcu_h_dqgr30, bcu_h_drho,                                 &
  bcu_h_thw, bcu_h_thvw, bcu_h_qw, bcu_h_qclw, bcu_h_qcfw, bcu_h_qrainw,    &
  bcu_h_qgraupw, bcu_h_uw, bcu_h_vw, bcu_h_ww, bcu_h_w3, bcu_h_vv,          &
  bcu_h_uu, bcu_h_uv,                                                       &
  bcu_h_uth, bcu_h_vth, bcu_h_uthv, bcu_h_vthv, bcu_h_uq, bcu_h_vq,bcu_h_wp,&
  bcu_hdr, errorstatus )

  ! Also write out plume info

  CALL  crmstyle_write_plume_ff(icall_type,maxfldsout, example_hdr,         &
                              ref_modlev_field,                             &
                              n_plume,plume_size,plume_diam_pdf,            &
                              bcu_hdr, errorstatus )

END IF
!-----------------------------------------------------------------------
! wg1 file
!-----------------------------------------------------------------------
IF (l_wg1) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("WG1_FILE", wg1_hdr % FileName)

    CALL assign_file_unit(wg1_hdr % FileName, wg1_hdr % UnitNum,  &
                          handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  wg1_h_w, wg1_h_u, wg1_h_v, wg1_h_th, wg1_h_thv, wg1_h_rho,                &
  wg1_h_rh, wg1_h_a,  wg1_h_dpx, wg1_h_dpy,                                 &
  wg1_h_q,   wg1_h_qcl, wg1_h_qcf, wg1_h_qrain, wg1_h_qgraup,               &
  wg1_h_dt1, wg1_h_dt2, wg1_h_dt4, wg1_h_dt9, wg1_h_dt12, wg1_h_dt30,       &
  wg1_h_dq4, wg1_h_dq9, wg1_h_dq12, wg1_h_dq30,                             &
  wg1_h_dqcl4, wg1_h_dqcl9, wg1_h_dqcl12,  wg1_h_dqcl30,                    &
  wg1_h_dqcf4, wg1_h_dqcf3, wg1_h_dqcf12,  wg1_h_dqcf30,                    &
  wg1_h_dqrain30, wg1_h_dqgr30, wg1_h_drho,                                 &
  wg1_h_thw, wg1_h_thvw, wg1_h_qw, wg1_h_qclw, wg1_h_qcfw, wg1_h_qrainw,    &
  wg1_h_qgraupw, wg1_h_uw, wg1_h_vw, wg1_h_ww, wg1_h_w3, wg1_h_vv,          &
  wg1_h_uu, wg1_h_uv,                                                       &
  wg1_h_uth, wg1_h_vth, wg1_h_uthv, wg1_h_vthv, wg1_h_uq, wg1_h_vq,wg1_h_wp,&
  wg1_hdr, errorstatus )

END IF
!-----------------------------------------------------------------------
! ppd file
!-----------------------------------------------------------------------
IF (l_ppd) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("PPD_FILE", ppd_hdr % FileName)

    CALL assign_file_unit(ppd_hdr % FileName, ppd_hdr % UnitNum, &
                          handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  ppd_h_w, ppd_h_u, ppd_h_v, ppd_h_th, ppd_h_thv, ppd_h_rho,                &
  ppd_h_rh, ppd_h_a,  ppd_h_dpx, ppd_h_dpy,                                 &
  ppd_h_q,   ppd_h_qcl, ppd_h_qcf, ppd_h_qrain, ppd_h_qgraup,               &
  ppd_h_dt1, ppd_h_dt2, ppd_h_dt4, ppd_h_dt9, ppd_h_dt12, ppd_h_dt30,       &
  ppd_h_dq4, ppd_h_dq9, ppd_h_dq12, ppd_h_dq30,                             &
  ppd_h_dqcl4, ppd_h_dqcl9, ppd_h_dqcl12,  ppd_h_dqcl30,                    &
  ppd_h_dqcf4, ppd_h_dqcf3, ppd_h_dqcf12,  ppd_h_dqcf30,                    &
  ppd_h_dqrain30, ppd_h_dqgr30, ppd_h_drho,                                 &
  ppd_h_thw, ppd_h_thvw, ppd_h_qw, ppd_h_qclw, ppd_h_qcfw, ppd_h_qrainw,    &
  ppd_h_qgraupw, ppd_h_uw, ppd_h_vw, ppd_h_ww, ppd_h_w3, ppd_h_vv,          &
  ppd_h_uu, ppd_h_uv,                                                       &
  ppd_h_uth, ppd_h_vth, ppd_h_uthv, ppd_h_vthv, ppd_h_uq, ppd_h_vq,ppd_h_wp,&
  ppd_hdr, errorstatus )

END IF
!-----------------------------------------------------------------------
! nbd file
!-----------------------------------------------------------------------
IF (l_nbd) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("NBD_FILE", nbd_hdr % FileName)

    CALL assign_file_unit(nbd_hdr % FileName, nbd_hdr % UnitNum, &
                 handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  nbd_h_w, nbd_h_u, nbd_h_v, nbd_h_th, nbd_h_thv, nbd_h_rho,                &
  nbd_h_rh, nbd_h_a,  nbd_h_dpx, nbd_h_dpy,                                 &
  nbd_h_q,   nbd_h_qcl, nbd_h_qcf, nbd_h_qrain, nbd_h_qgraup,               &
  nbd_h_dt1, nbd_h_dt2, nbd_h_dt4, nbd_h_dt9, nbd_h_dt12, nbd_h_dt30,       &
  nbd_h_dq4, nbd_h_dq9, nbd_h_dq12, nbd_h_dq30,                             &
  nbd_h_dqcl4, nbd_h_dqcl9, nbd_h_dqcl12,  nbd_h_dqcl30,                    &
  nbd_h_dqcf4, nbd_h_dqcf3, nbd_h_dqcf12,  nbd_h_dqcf30,                    &
  nbd_h_dqrain30, nbd_h_dqgr30, nbd_h_drho,                                 &
  nbd_h_thw, nbd_h_thvw, nbd_h_qw, nbd_h_qclw, nbd_h_qcfw, nbd_h_qrainw,    &
  nbd_h_qgraupw, nbd_h_uw, nbd_h_vw, nbd_h_ww, nbd_h_w3, nbd_h_vv,          &
  nbd_h_uu, nbd_h_uv,                                                       &
  nbd_h_uth, nbd_h_vth, nbd_h_uthv, nbd_h_vthv, nbd_h_uq, nbd_h_vq,nbd_h_wp,&
  nbd_hdr, errorstatus )

END IF
!-----------------------------------------------------------------------
! nid file 
!-----------------------------------------------------------------------
IF (l_nid) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("NID_FILE", nid_hdr % FileName)

    CALL assign_file_unit(nid_hdr % FileName, nid_hdr % UnitNum,  &
                          handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  nid_h_w, nid_h_u, nid_h_v, nid_h_th, nid_h_thv, nid_h_rho,                &
  nid_h_rh, nid_h_a,  nid_h_dpx, nid_h_dpy,                                 &
  nid_h_q,   nid_h_qcl, nid_h_qcf, nid_h_qrain, nid_h_qgraup,               &
  nid_h_dt1, nid_h_dt2, nid_h_dt4, nid_h_dt9, nid_h_dt12, nid_h_dt30,       &
  nid_h_dq4, nid_h_dq9, nid_h_dq12, nid_h_dq30,                             &
  nid_h_dqcl4, nid_h_dqcl9, nid_h_dqcl12,  nid_h_dqcl30,                    &
  nid_h_dqcf4, nid_h_dqcf3, nid_h_dqcf12,  nid_h_dqcf30,                    &
  nid_h_dqrain30, nid_h_dqgr30, nid_h_drho,                                 &
  nid_h_thw, nid_h_thvw, nid_h_qw, nid_h_qclw, nid_h_qcfw, nid_h_qrainw,    &
  nid_h_qgraupw, nid_h_uw, nid_h_vw, nid_h_ww, nid_h_w3, nid_h_vv,          &
  nid_h_uu, nid_h_uv,                                                       &
  nid_h_uth, nid_h_vth, nid_h_uthv, nid_h_vthv, nid_h_uq, nid_h_vq,nid_h_wp,&
  nid_hdr, errorstatus )

END IF
!-----------------------------------------------------------------------
! adu file
!-----------------------------------------------------------------------
IF (l_adu) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("ADU_FILE", adu_hdr % FileName)

    CALL assign_file_unit(adu_hdr % FileName, adu_hdr % UnitNum,  &
                          handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  adu_h_w, adu_h_u, adu_h_v, adu_h_th, adu_h_thv, adu_h_rho,                &
  adu_h_rh, adu_h_a,  adu_h_dpx, adu_h_dpy,                                 &
  adu_h_q,   adu_h_qcl, adu_h_qcf, adu_h_qrain, adu_h_qgraup,               &
  adu_h_dt1, adu_h_dt2, adu_h_dt4, adu_h_dt9, adu_h_dt12, adu_h_dt30,       &
  adu_h_dq4, adu_h_dq9, adu_h_dq12, adu_h_dq30,                             &
  adu_h_dqcl4, adu_h_dqcl9, adu_h_dqcl12,  adu_h_dqcl30,                    &
  adu_h_dqcf4, adu_h_dqcf3, adu_h_dqcf12,  adu_h_dqcf30,                    &
  adu_h_dqrain30, adu_h_dqgr30, adu_h_drho,                                 &
  adu_h_thw, adu_h_thvw, adu_h_qw, adu_h_qclw, adu_h_qcfw, adu_h_qrainw,    &
  adu_h_qgraupw, adu_h_uw, adu_h_vw, adu_h_ww, adu_h_w3, adu_h_vv,          &
  adu_h_uu, adu_h_uv,                                                       &
  adu_h_uth, adu_h_vth, adu_h_uthv, adu_h_vthv, adu_h_uq, adu_h_vq,adu_h_wp,&
  adu_hdr, errorstatus )

END IF
!-----------------------------------------------------------------------
! acw file
!-----------------------------------------------------------------------
IF (l_acw) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("ACW_FILE", acw_hdr % FileName)

    CALL assign_file_unit(acw_hdr % FileName, acw_hdr % UnitNum,  &
                          handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  acw_h_w, acw_h_u, acw_h_v, acw_h_th, acw_h_thv, acw_h_rho,                &
  acw_h_rh, acw_h_a,  acw_h_dpx, acw_h_dpy,                                 &
  acw_h_q,   acw_h_qcl, acw_h_qcf, acw_h_qrain, acw_h_qgraup,               &
  acw_h_dt1, acw_h_dt2, acw_h_dt4, acw_h_dt9, acw_h_dt12, acw_h_dt30,       &
  acw_h_dq4, acw_h_dq9, acw_h_dq12, acw_h_dq30,                             &
  acw_h_dqcl4, acw_h_dqcl9, acw_h_dqcl12,  acw_h_dqcl30,                    &
  acw_h_dqcf4, acw_h_dqcf3, acw_h_dqcf12,  acw_h_dqcf30,                    &
  acw_h_dqrain30, acw_h_dqgr30, acw_h_drho,                                 &
  acw_h_thw, acw_h_thvw, acw_h_qw, acw_h_qclw, acw_h_qcfw, acw_h_qrainw,    &
  acw_h_qgraupw, acw_h_uw, acw_h_vw, acw_h_ww, acw_h_w3, acw_h_vv,          &
  acw_h_uu, acw_h_uv,                                                       &
  acw_h_uth, acw_h_vth, acw_h_uthv, acw_h_vthv, acw_h_uq, acw_h_vq,acw_h_wp,&
  acw_hdr, errorstatus )

END IF
!-----------------------------------------------------------------------
! bcw file 
!-----------------------------------------------------------------------
IF (l_bcw) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("BCW_FILE", bcw_hdr % FileName)

    CALL assign_file_unit(bcw_hdr % FileName, bcw_hdr % UnitNum, &
                          handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  bcw_h_w, bcw_h_u, bcw_h_v, bcw_h_th, bcw_h_thv, bcw_h_rho,                &
  bcw_h_rh, bcw_h_a,  bcw_h_dpx, bcw_h_dpy,                                 &
  bcw_h_q,   bcw_h_qcl, bcw_h_qcf, bcw_h_qrain, bcw_h_qgraup,               &
  bcw_h_dt1, bcw_h_dt2, bcw_h_dt4, bcw_h_dt9, bcw_h_dt12, bcw_h_dt30,       &
  bcw_h_dq4, bcw_h_dq9, bcw_h_dq12, bcw_h_dq30,                             &
  bcw_h_dqcl4, bcw_h_dqcl9, bcw_h_dqcl12,  bcw_h_dqcl30,                    &
  bcw_h_dqcf4, bcw_h_dqcf3, bcw_h_dqcf12,  bcw_h_dqcf30,                    &
  bcw_h_dqrain30, bcw_h_dqgr30, bcw_h_drho,                                 &
  bcw_h_thw, bcw_h_thvw, bcw_h_qw, bcw_h_qclw, bcw_h_qcfw, bcw_h_qrainw,    &
  bcw_h_qgraupw, bcw_h_uw, bcw_h_vw, bcw_h_ww, bcw_h_w3, bcw_h_vv,          &
  bcw_h_uu, bcw_h_uv,                                                       &
  bcw_h_uth, bcw_h_vth, bcw_h_uthv, bcw_h_vthv, bcw_h_uq, bcw_h_vq,bcw_h_wp,&
  bcw_hdr, errorstatus )

END IF
!-----------------------------------------------------------------------
! ucu file
!-----------------------------------------------------------------------
IF (l_ucu) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("UCU_FILE", ucu_hdr % FileName)

    CALL assign_file_unit(ucu_hdr % FileName, ucu_hdr % UnitNum,  &
                          handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  ucu_h_w, ucu_h_u, ucu_h_v, ucu_h_th, ucu_h_thv, ucu_h_rho,                &
  ucu_h_rh, ucu_h_a,  ucu_h_dpx, ucu_h_dpy,                                 &
  ucu_h_q,   ucu_h_qcl, ucu_h_qcf, ucu_h_qrain, ucu_h_qgraup,               &
  ucu_h_dt1, ucu_h_dt2, ucu_h_dt4, ucu_h_dt9, ucu_h_dt12, ucu_h_dt30,       &
  ucu_h_dq4, ucu_h_dq9, ucu_h_dq12, ucu_h_dq30,                             &
  ucu_h_dqcl4, ucu_h_dqcl9, ucu_h_dqcl12,  ucu_h_dqcl30,                    &
  ucu_h_dqcf4, ucu_h_dqcf3, ucu_h_dqcf12,  ucu_h_dqcf30,                    &
  ucu_h_dqrain30, ucu_h_dqgr30, ucu_h_drho,                                 &
  ucu_h_thw, ucu_h_thvw, ucu_h_qw, ucu_h_qclw, ucu_h_qcfw, ucu_h_qrainw,    &
  ucu_h_qgraupw, ucu_h_uw, ucu_h_vw, ucu_h_ww, ucu_h_w3, ucu_h_vv,          &
  ucu_h_uu, ucu_h_uv,                                                       &
  ucu_h_uth, ucu_h_vth, ucu_h_uthv, ucu_h_vthv, ucu_h_uq, ucu_h_vq,ucu_h_wp,&
  ucu_hdr, errorstatus )

END IF
!-----------------------------------------------------------------------
! ppw file
!-----------------------------------------------------------------------
IF (l_ppw) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("PPW_FILE", ppw_hdr % FileName)

    CALL assign_file_unit(ppw_hdr % FileName, ppw_hdr % UnitNum,  &
                          handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  ppw_h_w, ppw_h_u, ppw_h_v, ppw_h_th, ppw_h_thv, ppw_h_rho,                &
  ppw_h_rh, ppw_h_a,  ppw_h_dpx, ppw_h_dpy,                                 &
  ppw_h_q,   ppw_h_qcl, ppw_h_qcf, ppw_h_qrain, ppw_h_qgraup,               &
  ppw_h_dt1, ppw_h_dt2, ppw_h_dt4, ppw_h_dt9, ppw_h_dt12, ppw_h_dt30,       &
  ppw_h_dq4, ppw_h_dq9, ppw_h_dq12, ppw_h_dq30,                             &
  ppw_h_dqcl4, ppw_h_dqcl9, ppw_h_dqcl12,  ppw_h_dqcl30,                    &
  ppw_h_dqcf4, ppw_h_dqcf3, ppw_h_dqcf12,  ppw_h_dqcf30,                    &
  ppw_h_dqrain30, ppw_h_dqgr30, ppw_h_drho,                                 &
  ppw_h_thw, ppw_h_thvw, ppw_h_qw, ppw_h_qclw, ppw_h_qcfw, ppw_h_qrainw,    &
  ppw_h_qgraupw, ppw_h_uw, ppw_h_vw, ppw_h_ww, ppw_h_w3, ppw_h_vv,          &
  ppw_h_uu, ppw_h_uv,                                                       &
  ppw_h_uth, ppw_h_vth, ppw_h_uthv, ppw_h_vthv, ppw_h_uq, ppw_h_vq,ppw_h_wp,&
  ppw_hdr, errorstatus )

END IF
!-----------------------------------------------------------------------
! nbw file
!-----------------------------------------------------------------------
IF (l_nbw) THEN
  IF (icall_type == 1) THEN
    CALL get_env_var("NBW_FILE", nbw_hdr % FileName)

    CALL assign_file_unit(nbw_hdr % FileName, nbw_hdr % UnitNum,  &
                          handler="portio")
  END IF

  CALL  crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,           &
  ref_modlev_field,                                                         &
  nbw_h_w, nbw_h_u, nbw_h_v, nbw_h_th, nbw_h_thv, nbw_h_rho,                &
  nbw_h_rh, nbw_h_a,  nbw_h_dpx, nbw_h_dpy,                                 &
  nbw_h_q,   nbw_h_qcl, nbw_h_qcf, nbw_h_qrain, nbw_h_qgraup,               &
  nbw_h_dt1, nbw_h_dt2, nbw_h_dt4, nbw_h_dt9, nbw_h_dt12, nbw_h_dt30,       &
  nbw_h_dq4, nbw_h_dq9, nbw_h_dq12, nbw_h_dq30,                             &
  nbw_h_dqcl4, nbw_h_dqcl9, nbw_h_dqcl12,  nbw_h_dqcl30,                    &
  nbw_h_dqcf4, nbw_h_dqcf3, nbw_h_dqcf12,  nbw_h_dqcf30,                    &
  nbw_h_dqrain30, nbw_h_dqgr30, nbw_h_drho,                                 &
  nbw_h_thw, nbw_h_thvw, nbw_h_qw, nbw_h_qclw, nbw_h_qcfw, nbw_h_qrainw,    &
  nbw_h_qgraupw, nbw_h_uw, nbw_h_vw, nbw_h_ww, nbw_h_w3, nbw_h_vv,          &
  nbw_h_uu, nbw_h_uv,                                                       &
  nbw_h_uth, nbw_h_vth, nbw_h_uthv, nbw_h_vthv, nbw_h_uq, nbw_h_vq,nbw_h_wp,&
  nbw_hdr, errorstatus )

END IF
!-----------------------------------------------------------------------
! Cleanup space

DEALLOCATE(example_hdr%LevDepC)
DEALLOCATE(example_hdr%RealC)
DEALLOCATE(example_hdr%IntC)
DEALLOCATE(example_hdr%ColDepC)
DEALLOCATE(example_hdr%RowDepC)
DEALLOCATE(example_hdr%FixHd)

NULLIFY(ref_single_field(1)%RData )

DO k=1,mlevs
  NULLIFY(ref_modlev_field(k)%RData )
END DO

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_write_ff_out

END MODULE crmstyle_write_ff_out_mod

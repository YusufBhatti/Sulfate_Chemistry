! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Dynamics Diagnostics
!
! Subroutine Interface:
MODULE sl_wind_norm_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SL_WIND_NORM_MOD'

CONTAINS
SUBROUTINE sl_wind_norm(                                                &
                        level_start, level_end, rims,                   &
                        etadot, u_np1, v_np1, w_np1, etadot_np1,        &
                        u, v, w, r_u, r_v, r_w, r_u_d, r_v_d, r_w_d,    &
                        L_do_rims, L_do_halos, L_print_pe )

USE parkind1,          ONLY: jpim, jprb       !DrHook
USE yomhook,           ONLY: lhook, dr_hook   !DrHook

USE array_norm_mod, ONLY: array_norm
USE pr_block4_mod, ONLY: dom_no, dom_eo, dom_s, dom_w
USE nlsizes_namelist_mod, ONLY: row_length, rows, n_rows, model_levels, &
                                global_row_length, global_rows
USE UM_parvars, ONLY: offx, offy, halo_i, halo_j
USE um_parcore, ONLY: mype

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr, ONLY: umMessage, umPrint, umPrintFlush

IMPLICIT NONE
!
! Description:
!   Printing of norms for arrays after eg_sl_full_wind
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SL_WIND_NORM'

! Arguments with Intent IN. ie: Input variables.
LOGICAL, INTENT(IN) :: L_do_rims
LOGICAL, INTENT(IN) :: L_do_halos
LOGICAL, INTENT(IN) :: L_print_pe   ! print diagnostics on all pe's if true

! Integer parameters for advection

INTEGER, INTENT(IN) :: level_start, level_end, rims

! Model dimensions

REAL, INTENT(IN) ::                                                   &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,                 &
         model_levels),                                               &
  v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,               &
        model_levels),                                                &
  w(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                 &
         0:model_levels),                                             &
  u_np1(1-offx:row_length+offx,1-offy:rows+offy, model_levels),       &
  v_np1(1-offx:row_length+offx,1-offy:n_rows+offy, model_levels),     &
  w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,                 &
             0:model_levels)

! Timelevel n arrival point quantities

REAL, INTENT(IN) ::                                                  &
  r_u(1-offx:row_length+offx,1-offy:rows+offy,model_levels),         &
  r_v(1-offx:row_length+offx,1-offy:n_rows+offy,model_levels),       &
  r_w(row_length,rows,0:model_levels),                               &
  r_u_d(1-offx:row_length+offx,1-offy:rows+offy,model_levels),       &
  r_v_d(1-offx:row_length+offx,1-offy:n_rows+offy,model_levels),     &
  r_w_d(row_length,rows,0:model_levels)

! Local variables

INTEGER :: k, num
INTEGER :: levels      !  model_levels + 1 
INTEGER :: start_level, end_level
INTEGER :: k0, s0      !  level zero pointers
INTEGER :: rim_n, rim_s, rim_e, rim_w

REAL ::  two_norm(14)

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( L_do_rims ) THEN
  rim_n = 1
  rim_s = 1
  rim_e = 1
  rim_w = 1
ELSE ! L_do_rims
  rim_n = rims
  rim_s = rims
  rim_e = rims
  rim_w = rims
!  set up block norm if dom_w and dom_s /= 1
  IF ( dom_w /= 1 .AND. dom_s /= 1 ) THEN
    rim_s = dom_s
    rim_w = dom_w
    rim_n = global_rows - dom_no
    rim_e = global_row_length - dom_eo
  END IF ! dom_w /= 1 .AND. dom_s /= 1
END IF ! L_do_rims

levels = model_levels + 1
start_level = level_start
end_level = level_end
s0 = 0
IF ( level_start == level_end ) THEN
  IF ( level_start < 1 ) THEN
    start_level = 1
    end_level   = 1
  ELSE IF ( level_start == model_levels) THEN
    start_level = 1
    s0 = 1
    levels = model_levels
  END IF ! level_start == model_levels 
ELSE     !    0 < level_start < model_levels
  s0 = 1
  levels = model_levels
END IF   !   level_start == level_end

IF ( level_start == level_end ) THEN

  DO k =  start_level, end_level
    k0 = k
    IF( level_start < 1 ) k0 = 0
    num = 1
    CALL array_norm(                                                 &
                    etadot(1-offx,1-offy,s0), row_length, rows,      &
                    levels, k, k,                                    &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,         &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                                 &
                    etadot_np1(1-offx,1-offy,s0), row_length, rows,  &
                    levels, k, k,                                    &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,         &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                              &
                    u, row_length, rows,                          &
                    model_levels, k, k,                           &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,  &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                              &
                    u_np1, row_length, rows,                      &
                    model_levels, k, k,                           &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                              &
                    v, row_length, n_rows,                        &
                    model_levels, k, k,                           &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,  &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                              &
                    v_np1, row_length, n_rows,                    &
                    model_levels, k, k,                           &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                               &
                    w(1-halo_i,1-halo_j,s0), row_length, rows,     &
                    levels, k, k,                                  &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                               &
                    w_np1(1-offx,1-offy,s0), row_length, rows,     &
                    levels, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                              &
                    r_u, row_length, rows,                        &
                    model_levels, k, k,                           &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                              &
                    r_u_d, row_length, rows,                      &
                    model_levels, k, k,                           &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                              &
                    r_v, row_length, n_rows,                      &
                    model_levels, k, k,                           &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                              &
                    r_v_d, row_length, n_rows,                    &
                    model_levels, k, k,                           &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                              &
                    r_w(1,1,s0), row_length, rows,                &
                    levels, k, k,                                 &
                    0, 0, 0, 0, L_do_halos, L_do_rims,            &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                              &
                    r_w_d(1,1,s0), row_length, rows,              &
                    levels, k, k,                                 &
                    0, 0, 0, 0, L_do_halos, L_do_rims,            &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )

    IF ( L_print_pe .OR. mype == 0 ) THEN
      num = 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k0, ' etadot Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
            ' Level ', k0, ' etadot_np1 Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k, ' u Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k, ' u_np1 Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k, ' v Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k, ' v_np1 Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k0, ' w Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k0, ' w_np1 Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k, ' r_u Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k, ' r_u_d Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k, ' r_v Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k, ' r_v_d Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k0, ' r_w Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                      &
                ' Level ', k0, ' r_w_d Two_Norm = ' , two_norm(num)
      CALL umPrint(umMessage, src='SL_WIND_NORM')
    END IF ! L_print_pe .OR. mype == 0

  END DO !  k = 1, model_levels

ELSE   !  3d norms from level_start to level_end  

  num = 1
  CALL array_norm(                                                 &
                  etadot(1-offx,1-offy,s0), row_length, rows,      &
                  levels, start_level, level_end,                  &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,         &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                                 &
                  etadot_np1(1-offx,1-offy,s0), row_length, rows,  &
                  levels, start_level, level_end,                  &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,         &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  u, row_length, rows,                          &
                  model_levels, start_level, level_end,         &
                  halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,  &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  u_np1, row_length, rows,                      &
                  model_levels, start_level, level_end,         &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  v, row_length, n_rows,                        &
                  model_levels, start_level, level_end,         &
                  halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,  &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  v_np1, row_length, n_rows,                    &
                  model_levels, start_level, level_end,         &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  w(1-halo_i,1-halo_j,s0), row_length, rows,    &
                  levels, start_level, level_end,               &
                  halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,  &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  w_np1(1-offx,1-offy,s0), row_length, rows,    &
                  levels, start_level, level_end,               &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  r_u, row_length, rows,                        &
                  model_levels, start_level, level_end,         &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  r_u_d, row_length, rows,                      &
                  model_levels, start_level, level_end,         &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  r_v, row_length, n_rows,                      &
                  model_levels, start_level, level_end,         &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  r_v_d, row_length, n_rows,                    &
                  model_levels, start_level, level_end,         &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  r_w(1,1,s0), row_length, rows,                &
                  levels, start_level, level_end,               &
                  0, 0, 0, 0, L_do_halos, L_do_rims,            &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                              &
                  r_w_d(1,1,s0), row_length, rows,              &
                  levels, start_level, level_end,               &
                  0, 0, 0, 0, L_do_halos, L_do_rims,            &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )

  IF ( L_print_pe .OR. mype == 0 ) THEN
    num = 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' etadot Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' etadot_np1 Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' u Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' u_np1 Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' v Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' v_np1 Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' w Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' w_np1 Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' r_u Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' r_u_d Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' r_v Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' r_v_d Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, E23.16)') ' Levels 1 to', level_end,  &
                  ' r_w Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, E23.16)') ' Levels 1 to', level_end,  &
                  ' r_w_d Two_Norm = ' , two_norm(num)
    CALL umPrint(umMessage, src='SL_WIND_NORM')
  END IF ! L_print_pe .OR. mype == 0

END IF ! level_start == end_level

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sl_wind_norm

END MODULE sl_wind_norm_mod

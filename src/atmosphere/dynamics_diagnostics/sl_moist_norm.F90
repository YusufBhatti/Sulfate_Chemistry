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
MODULE sl_moist_norm_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SL_MOIST_NORM_MOD'

CONTAINS
SUBROUTINE sl_moist_norm(                                             &
                         level_start, level_end, rims,                &
                         exner_star,                                  &
                         r_m_v_d, r_m_cl_d, r_m_cf_d,                 &
                         r_m_r_d, r_m_gr_d, r_m_cf2_d,                &
                         cf_star, cfl_star, cff_star,                 &
                         L_do_rims, L_do_halos, L_print_pe )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE mphys_inputs_mod,  ONLY: l_mcr_qrain, l_mcr_qcf2, l_mcr_qgraup
USE cloud_inputs_mod,  ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2
USE array_norm_mod, ONLY: array_norm
USE pr_block4_mod, ONLY: dom_no, dom_eo, dom_s, dom_w
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels,        &
                                 global_row_length, global_rows
USE UM_parvars, ONLY: offx, offy
USE um_parcore, ONLY: mype

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr, ONLY: umPrint, umMessage, umPrintFlush

!
! Description:  Printing of norms for arrays after eg_sl_moist
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SL_MOIST_NORM'

! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_do_rims
LOGICAL :: L_do_halos
LOGICAL :: L_print_pe   ! print diagnostics on all pe's if true

INTEGER, INTENT(IN) :: level_start, level_end, rims

! Timelevel n departure point quantities
REAL, INTENT(IN) ::                                                    &
  exner_star(1-offx:row_length+offx,1-offy:rows+offy, model_levels),   &
  cf_star  (row_length, rows, 0:model_levels),                         &
  cfl_star (row_length, rows, 0:model_levels),                         &
  cff_star (row_length, rows, 0:model_levels),                         &
  r_m_v_d  (row_length, rows, 0:model_levels),                         &
  r_m_cl_d (row_length, rows, 0:model_levels),                         &
  r_m_cf_d (row_length, rows, 0:model_levels),                         &
  r_m_r_d  (row_length, rows, 0:model_levels),                         &
  r_m_gr_d (row_length, rows, 0:model_levels),                         &
  r_m_cf2_d(row_length, rows, 0:model_levels)

! Local variables

INTEGER :: k, num
INTEGER :: levels      !  model_levels + 1 
INTEGER :: start_level, end_level
INTEGER :: k0, s0      !  level zero pointers
INTEGER :: rim_n, rim_s, rim_e, rim_w

! Local array for storing two_norms
REAL ::  two_norm(10)

! ----------------------------------------------------------------------
! 1.0 Start of subroutine code
! ----------------------------------------------------------------------

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
    CALL array_norm(                                               &
                    r_m_v_d(1,1,s0), row_length, rows,             &
                    levels, k, k,                                  &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    r_m_cl_d(1,1,s0), row_length, rows,            &
                    levels, k, k,                                  &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                               &
                    r_m_cf_d(1,1,s0), row_length, rows,            &
                    levels, k, k,                                  &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    IF ( l_mcr_qrain ) THEN
      num = num + 1
      CALL array_norm(                                             &
                    r_m_r_d(1,1,s0), row_length, rows,             &
                    levels, k, k,                                  &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! l_mcr_qrain
    IF ( l_mcr_qcf2 ) THEN
      num = num + 1
      CALL array_norm(                                             &
                    r_m_cf2_d(1,1,s0), row_length, rows,           &
                    levels, k, k,                                  &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! l_cf2_rain
    IF ( l_mcr_qgraup ) THEN
      num = num + 1
      CALL array_norm(                                             &
                    r_m_gr_d(1,1,s0), row_length, rows,            &
                    levels, k, k,                                  &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! l_mcr_qgraup
    IF ( i_cld_vn == i_cld_pc2 ) THEN
      num = num + 1
      CALL array_norm(                                             &
                      cf_star(1,1,s0), row_length, rows,           &
                      levels, k, k,                                &
                      0, 0, 0, 0, L_do_halos, L_do_rims,           &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num) )
      num = num + 1
      CALL array_norm(                                             &
                      cfl_star(1,1,s0), row_length, rows,          &
                      levels, k, k,                                &
                      0, 0, 0, 0, L_do_halos, L_do_rims,           &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num) )
      num = num + 1
      CALL array_norm(                                             &
                      cff_star(1,1,s0), row_length, rows,          &
                      levels, k, k,                                &
                      0, 0, 0, 0, L_do_halos, L_do_rims,           &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num) )
      num = num + 1
      CALL array_norm(                                             &
                      exner_star, row_length, rows,                &
                      model_levels, k, k,                          &
                      offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! i_cld_vn == i_cld_pc2

    IF ( L_print_pe .OR. mype == 0 ) THEN
      num = 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' r_m_v_d Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_MOIST_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' r_m_cl_d Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_MOIST_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' r_m_cf_d Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_MOIST_NORM')
      IF ( l_mcr_qrain ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
                ' Level ', k0, ' r_m_r_d Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_MOIST_NORM')
      END IF ! l_mcr_qrain
      IF ( l_mcr_qcf2 ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
                ' Level ', k0, ' r_m_cf2_d Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_MOIST_NORM')
      END IF ! l_mcr_qcf2
      IF ( l_mcr_qgraup ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
                ' Level ', k0, ' r_m_gr_d Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_MOIST_NORM')
      END IF ! l_mcr_qgraup
      IF ( i_cld_vn == i_cld_pc2 ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
                ' Level ', k0, ' cf_star Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_MOIST_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
                ' Level ', k0, ' cfl_star Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_MOIST_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
                ' Level ', k0, ' cff_star Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_MOIST_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
        ' Level ', k, ' exner_star Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='SL_MOIST_NORM')
      END IF ! i_cld_vn == i_cld_pc2
    END IF ! L_print_pe .OR. mype == 0

  END DO !  k = 0, model_levels

ELSE   !  3d norms from   level_start to level_end

  num = 1
  CALL array_norm(                                               &
                  r_m_v_d(1,1,s0), row_length, rows,             &
                  levels, start_level, level_end,                &
                  0, 0, 0, 0, L_do_halos, L_do_rims,             &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  r_m_cl_d(1,1,s0), row_length, rows,            &
                  levels, start_level, level_end,                &
                  0, 0, 0, 0, L_do_halos, L_do_rims,             &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                               &
                  r_m_cf_d(1,1,s0), row_length, rows,            &
                  levels, start_level, level_end,                &
                  0, 0, 0, 0, L_do_halos, L_do_rims,             &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  IF ( l_mcr_qrain ) THEN
    num = num + 1
    CALL array_norm(                                             &
                  r_m_r_d(1,1,s0), row_length, rows,             &
                  levels, start_level, level_end,                &
                  0, 0, 0, 0, L_do_halos, L_do_rims,             &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! l_mcr_qrain
  IF ( l_mcr_qcf2 ) THEN
    num = num + 1
    CALL array_norm(                                             &
                  r_m_cf2_d(1,1,s0), row_length, rows,           &
                  levels, start_level, level_end,                &
                  0, 0, 0, 0, L_do_halos, L_do_rims,             &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! l_mcr_qcf2
  IF ( l_mcr_qgraup ) THEN
    num = num + 1
    CALL array_norm(                                             &
                  r_m_gr_d(1,1,s0), row_length, rows,            &
                  levels, start_level, level_end,                &
                  0, 0, 0, 0, L_do_halos, L_do_rims,             &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! l_mcr_qgraup
  IF ( i_cld_vn == i_cld_pc2 ) THEN
    num = num + 1
    CALL array_norm(                                               &
                    cf_star(1,1,s0), row_length, rows,             &
                    levels, start_level, level_end,                &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                               &
                    cfl_star(1,1,s0), row_length, rows,            &
                    levels, start_level, level_end,                &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                               &
                    cff_star(1,1,s0), row_length, rows,            &
                    levels, start_level, level_end,                &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                               &
                    exner_star, row_length, rows,                  &
                    model_levels, start_level, level_end,          &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! i_cld_vn == i_cld_pc2


  IF ( L_print_pe .OR. mype == 0 ) THEN
    num = 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', level_start, ' to', level_end,      &
                  ' r_m_v_d Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='SL_MOIST_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', level_start, ' to', level_end,      &
                  ' r_m_cl_d Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='SL_MOIST_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', level_start, ' to', level_end,      &
                  ' r_m_cf_d Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='SL_MOIST_NORM')
    IF ( l_mcr_qrain ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', level_start, ' to', level_end,        &
                  ' r_m_r_d Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_MOIST_NORM')
    END IF ! l_mcr_qrain
    IF ( l_mcr_qcf2 ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', level_start, ' to', level_end,        &
                  ' r_m_cf2_d Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_MOIST_NORM')
    END IF ! l_mcr_qcf2
    IF ( l_mcr_qgraup ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', level_start, ' to', level_end,        &
                  ' r_m_gr_d Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_MOIST_NORM')
    END IF ! l_mcr_qgraup
    IF ( i_cld_vn == i_cld_pc2 ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', level_start, ' to', level_end,        &
                  ' cf_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_MOIST_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', level_start, ' to', level_end,        &
                  ' cfl_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_MOIST_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', level_start, ' to', level_end,        &
                  ' cff_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_MOIST_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', level_start, ' to', level_end,        &
                  ' exner_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='SL_MOIST_NORM')
    END IF ! i_cld_vn == i_cld_pc2
  END IF ! L_print_pe .OR. mype == 0

END IF ! level_start == level_end

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sl_moist_norm
END MODULE sl_moist_norm_mod

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
MODULE moist_norm_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MOIST_NORM_MOD'

CONTAINS
SUBROUTINE moist_norm(                                                  &
                      level_start, level_end, rims,                     &
                      m_v, m_cl, m_cf, m_r, m_gr, m_cf2,                &
                      L_do_rims, L_do_halos, L_print_pe )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE mphys_inputs_mod,  ONLY: l_mcr_qrain, l_mcr_qcf2, l_mcr_qgraup
USE array_norm_mod, ONLY: array_norm
USE pr_block4_mod, ONLY: dom_no, dom_eo, dom_s, dom_w
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels,        &
                                 global_row_length, global_rows
USE UM_parvars, ONLY: offx, offy
USE um_parcore, ONLY: mype

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr, ONLY: umMessage, umPrint, umPrintFlush

!
! Description:  Printing of moist norms for arrays after eg_sl_helmholtz
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MOIST_NORM'

! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_do_rims
LOGICAL :: L_do_halos
LOGICAL :: L_print_pe   ! print diagnostics on all pe's if true

INTEGER, INTENT(IN) :: level_start, level_end, rims

! Timelevel n departure point quantities
REAL, INTENT(IN) ::                                                     &
  m_v  (1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels),      &
  m_cl (1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels),      &
  m_cf (1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels),      &
  m_r  (1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels),      &
  m_gr (1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels),      &
  m_cf2(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)

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
                    m_v(1-offx,1-offy,s0), row_length, rows,       &
                    levels, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    m_cl(1-offx,1-offy,s0), row_length, rows,      &
                    levels, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                               &
                    m_cf(1-offx,1-offy,s0), row_length, rows,      &
                    levels, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    IF ( l_mcr_qrain ) THEN
      num = num + 1
      CALL array_norm(                                             &
                    m_r(1-offx,1-offy,s0), row_length, rows,       &
                    levels, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! l_mcr_qrain
    IF ( l_mcr_qgraup ) THEN
      num = num + 1
      CALL array_norm(                                             &
                    m_gr(1-offx,1-offy,s0), row_length, rows,      &
                    levels, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! l_mcr_qgraup
    IF ( l_mcr_qcf2 ) THEN
      num = num + 1
      CALL array_norm(                                             &
                    m_cf2(1-offx,1-offy,s0), row_length, rows,     &
                    levels, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! l_mcr_qcf2


    IF ( L_print_pe .OR. mype == 0 ) THEN
      num = 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' water vapour Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='MOIST_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' cloud water Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='MOIST_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' cloud ice Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='MOIST_NORM')
      IF ( l_mcr_qrain ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
                ' Level ', k0, ' rain Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='MOIST_NORM')
      END IF ! l_mcr_qrain
      IF ( l_mcr_qgraup ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
                ' Level ', k0, ' graupel Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='MOIST_NORM')
      END IF ! l_mcr_qgraup
      IF ( l_mcr_qcf2 ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
                ' Level ', k0, ' cf2 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='MOIST_NORM')
      END IF ! l_mcr_qcf2
    END IF ! L_print_pe .OR. mype == 0

  END DO !  start_level, end_level

ELSE   !  3d norms from  start_level to level_end

  num = 1
  CALL array_norm(                                               &
                  m_v(1-offx,1-offy,s0), row_length, rows,        &
                  levels, start_level, level_end,                 &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  m_cl(1-offx,1-offy,s0), row_length, rows,       &
                  levels, start_level, level_end,                 &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                               &
                  m_cf(1-offx,1-offy,s0), row_length, rows,       &
                  levels, start_level, level_end,                 &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  IF ( l_mcr_qrain ) THEN
    num = num + 1
    CALL array_norm(                                             &
                  m_r(1-offx,1-offy,s0), row_length, rows,        &
                  levels, start_level, level_end,                 &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! l_mcr_qrain
  IF ( l_mcr_qgraup ) THEN
    num = num + 1
    CALL array_norm(                                             &
                  m_gr(1-offx,1-offy,s0), row_length, rows,       &
                  levels, start_level, level_end,                 &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! l_mcr_qgraup
  IF ( l_mcr_qcf2 ) THEN
    num = num + 1
    CALL array_norm(                                             &
                  m_cf2(1-offx,1-offy,s0), row_length, rows,      &
                  levels, start_level, level_end,                 &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! l_mcr_qcf2


  IF ( L_print_pe .OR. mype == 0 ) THEN
    num = 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' water vapour Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='MOIST_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' cloud water Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='MOIST_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' cloud ice Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='MOIST_NORM')
    IF ( l_mcr_qrain ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' rain Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='MOIST_NORM')
    END IF ! l_mcr_qrain
    IF ( l_mcr_qgraup ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' graupel Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='MOIST_NORM')
    END IF ! l_mcr_qgraup
    IF ( l_mcr_qcf2 ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,        &
                  ' cf2 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='MOIST_NORM')
    END IF ! l_mcr_qcf2
  END IF ! L_print_pe .OR. mype == 0

END IF ! level_start == level_end

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE moist_norm
END MODULE moist_norm_mod

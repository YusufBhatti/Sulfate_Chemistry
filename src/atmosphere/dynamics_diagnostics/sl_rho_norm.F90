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
MODULE sl_rho_norm_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SL_RHO_NORM_MOD'

CONTAINS
SUBROUTINE sl_rho_norm(                                                 &
                       level_start, level_end, rims,                    &
                       r_rho, r_rho_d,                                  &
                       L_do_rims, L_do_halos, L_print_pe )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE array_norm_mod, ONLY: array_norm
USE pr_block4_mod, ONLY: dom_no, dom_eo, dom_s, dom_w
USE nlsizes_namelist_mod, ONLY: row_length, rows, model_levels,         &
                                global_row_length, global_rows
USE UM_parvars, ONLY: offx, offy
USE um_parcore, ONLY: mype

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr,    ONLY: umMessage, umPrint, umPrintFlush

IMPLICIT NONE
!
! Description:  Printing of norms for arrays after eg_sl_rho
!
! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SL_RHO_NORM'

! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_do_rims
LOGICAL :: L_do_halos
LOGICAL :: L_print_pe   ! print diagnostics on all pe's if true

INTEGER, INTENT(IN) :: level_start, level_end, rims

! Timelevel n arrival point quantities
REAL, INTENT(IN) ::                                                     &
  r_rho(1-offx:row_length+offx,1-offy:rows+offy, model_levels)

! Timelevel n departure point quantities
REAL, INTENT(IN) ::                                                     &
  r_rho_d(1-offx:row_length+offx,1-offy:rows+offy, model_levels)

! Local variables

INTEGER :: k
INTEGER :: rim_n, rim_s, rim_e, rim_w
INTEGER :: start_level, end_level

REAL ::  two_norm1, two_norm2

! 1.0 Start of subroutine code

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

start_level = level_start
end_level   = level_end
IF ( level_start == level_end ) THEN
  IF ( level_start < 1 ) THEN
    start_level = 1
    end_level   = 1
  ELSE IF ( level_start == model_levels) THEN
    start_level = 1
  END IF ! level_start == level_end 
END IF   !   level_start == level_end

IF ( level_start == level_end ) THEN

  DO k =  start_level, end_level

    CALL array_norm(                                               &
                    r_rho, row_length, rows,                       &
                    model_levels, k, k,                            &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm1 )

    CALL array_norm(                                               &
                    r_rho_d, row_length, rows,                     &
                    model_levels, k, k,                            &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm2 )

    IF ( L_print_pe .OR. mype == 0 ) THEN
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
                ' Level ', k, ' r_rho Two_Norm = ' , two_norm1
      CALL umPrint(umMessage, src='SL_RHO_NORM')
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
                ' Level ', k, ' r_rho_d Two_Norm = ' , two_norm2
      CALL umPrint(umMessage, src='SL_RHO_NORM')
    END IF ! L_print_pe .OR. mype == 0

  END DO !  k = 1, model_levels

ELSE   !  3d norms from start_level to level_end  

  CALL array_norm(                                               &
                   r_rho, row_length, rows,                      &
                   model_levels, start_level, level_end,         &
                   offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                   rim_n, rim_s, rim_e, rim_w, two_norm1 )
  CALL array_norm(                                               &
                   r_rho_d, row_length, rows,                    &
                   model_levels, start_level, level_end,         &
                   offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                   rim_n, rim_s, rim_e, rim_w, two_norm2 )

  IF ( L_print_pe .OR. mype == 0 ) THEN
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' r_rho Two_Norm = ' , two_norm1
    CALL umPrint(umMessage, src='SL_RHO_NORM')
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' r_rho_d Two_Norm = ' , two_norm2
    CALL umPrint(umMessage, src='SL_RHO_NORM')
  END IF ! L_print_pe .OR. mype == 0

END IF ! start_level == end_level

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sl_rho_norm
END MODULE sl_rho_norm_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Dynamics Diagnostics
!
MODULE eg_helm_norm_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_HELM_NORM_MOD'

CONTAINS
SUBROUTINE eg_helm_norm(                                                &
                        level_start, level_end, rims,                   &
                        exner_np1, p_star_np1, u_np1, v_np1, w_np1,     &
                        etadot_np1, rho_np1, thetav_np1,                &
                        exner_prime_term,                               &
                        L_do_rims, L_do_halos, L_print_pe )

USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim

USE array_norm_mod, ONLY: array_norm
USE pr_block4_mod, ONLY: dom_no, dom_eo, dom_s, dom_w
USE nlsizes_namelist_mod, ONLY: row_length, rows, n_rows, model_levels, &
                                 global_row_length, global_rows
USE UM_parvars, ONLY: offx, offy
USE um_parcore, ONLY: mype

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr, ONLY: umPrint, umMessage, umPrintFlush

IMPLICIT NONE

!
! Description:  Printing of norms for arrays eg_sl_helmholtz
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'EG_HELM_NORM'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_do_rims
LOGICAL :: L_do_halos
LOGICAL :: L_print_pe   ! print diagnostics on all pe's if true

INTEGER, INTENT(IN) :: level_start, level_end, rims

REAL, INTENT(IN)  ::                                                  &
 exner_np1(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1)

REAL, INTENT(IN)   ::                                                   &
      p_star_np1(1-offx:row_length+offx,1-offy:rows+offy)

REAL, INTENT(IN)   ::                                                   &
      u_np1(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),     &
      v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),   &
      rho_np1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL, INTENT(IN)   ::                                                   &
   thetav_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),  &
        w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),  &
   etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

REAL, INTENT(IN)   ::                                                   &
      exner_prime_term(1-offx:row_length+offx,1-offy:rows+offy,         &
                       model_levels)

! Local variables

INTEGER :: k, num
INTEGER :: levels, levelp
INTEGER :: start_level, end_level, end_press
INTEGER :: k0, s0      !  level zero pointers
INTEGER :: rim_n, rim_s, rim_e, rim_w

! Local array for storing two_norms
REAL ::  two_norm(8)

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
levelp = levels

CALL array_norm(                                                   &
                p_star_np1, row_length, rows, 1, 1, 1,             &
                offx, offy, 0, 0, L_do_halos, L_do_rims,           &
                rim_n, rim_s, rim_e, rim_w, two_norm(1))
IF ( L_print_pe .OR. mype == 0 ) THEN
  WRITE(umMessage, '(A, E23.16)')                                  &
          'Surface       exner_surf_np1 Two_Norm =' , two_norm(1)
  CALL umPrint(umMessage, src='EG_HELM_NORM')
END IF ! L_print_pe .OR. mype == 0

start_level = level_start
end_level = level_end
end_press = end_level + 1
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
                    exner_np1, row_length, rows,                   &
                    levelp, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    exner_prime_term, row_length, rows,            &
                    model_levels, k, k,                            &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    rho_np1, row_length, rows,                     &
                    model_levels, k, k,                            &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    u_np1, row_length, rows,                       &
                    model_levels, k, k,                            &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    v_np1, row_length, n_rows,                     &
                    model_levels, k, k,                            &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    w_np1(1-offx,1-offy,s0), row_length, rows,     &
                    levels, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                 &
                    etadot_np1(1-offx,1-offy,s0), row_length, rows,  &
                    levels, k, k,                                    &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,         &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                 &
                    thetav_np1(1-offx,1-offy,s0), row_length, rows,  &
                    levels, k, k,                                    &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,         &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))

    IF ( L_print_pe .OR. mype == 0 ) THEN
      num = 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k, ' exner_np1 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='EG_HELM_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
         ' Level ', k, ' exner_prime_term Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='EG_HELM_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k, ' dryrho_np1 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='EG_HELM_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k, ' u_np1 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='EG_HELM_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k, ' v_np1 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='EG_HELM_NORM')
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' w_np1 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='EG_HELM_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' etadotv_np1 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='EG_HELM_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' thetav_np1 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='EG_HELM_NORM')
    END IF ! L_print_pe .OR. mype == 0

  END DO !  start_level, end_level

ELSE   !  3d norms from start_level to end_level  

  num = 1
  CALL array_norm(                                                  &
                  exner_np1, row_length, rows,                      &
                  levelp, start_level, end_press,                   &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,          &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  exner_prime_term, row_length, rows,            &
                  model_levels, start_level, level_end,          &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  rho_np1, row_length, rows,                     &
                  model_levels, start_level, level_end,          &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  u_np1, row_length, rows,                       &
                  model_levels, start_level, level_end,          &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  v_np1, row_length, n_rows,                     &
                  model_levels, start_level, level_end,          &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  w_np1(1-offx, 1-offy, s0), row_length, rows,   &
                  levels, start_level, level_end,                &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                                &
                  etadot_np1(1-offx,1-offy,s0), row_length, rows, &
                  levels, start_level, level_end,                 &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                                &
                  thetav_np1(1-offx,1-offy,s0), row_length, rows, &
                  levels, start_level, level_end,                 &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,        &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))

  IF ( L_print_pe .OR. mype == 0 ) THEN
    num = 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', end_press,      &
                  ' exner_np1 Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='EG_HELM_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' exner_prime_term Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='EG_HELM_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' dryrho_np1 Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='EG_HELM_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' u_np1 Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='EG_HELM_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' v_np1 Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='EG_HELM_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' w_np1 Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='EG_HELM_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' etadot_np1 Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='EG_HELM_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' thetav_np1 Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='EG_HELM_NORM')
  END IF ! L_print_pe .OR. mype == 0

END IF !  level_start == level_end

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_helm_norm
END MODULE eg_helm_norm_mod

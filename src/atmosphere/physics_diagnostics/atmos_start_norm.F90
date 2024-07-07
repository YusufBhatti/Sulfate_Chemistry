! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Physics Diagnostics
!
! Subroutine Interface:
MODULE atmos_start_norm_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ATMOS_START_NORM_MOD'

CONTAINS
SUBROUTINE atmos_start_norm(                                            &
                            level_start, level_end, rims,               &
                            exner, rho, u, v, w,                        &
                            theta, q, qcl, qcf, qcf2, qrain, qgraup,    &
                            cf_bulk, cfl, cff,                          &
                            L_do_rims, L_do_halos, L_print_pe )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod,  ONLY: pdims_s, tdims_l, tdims_s,            &
                                  udims_s, vdims_s, wdims_s
USE array_norm_mod, ONLY: array_norm
USE pr_block4_mod, ONLY: dom_no, dom_eo, dom_s, dom_w
USE nlsizes_namelist_mod,  ONLY: row_length, rows,n_rows, model_levels, &
                                 global_row_length, global_rows
USE UM_parvars,  ONLY: offx, offy, halo_i, halo_j

USE mphys_inputs_mod,  ONLY: l_mcr_qrain, l_mcr_qcf2, l_mcr_qgraup
USE cloud_inputs_mod,  ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2
USE um_parcore, ONLY: mype

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr, ONLY: umPrint, umMessage, umPrintFlush

!
! Description:  Printing of norms around atmos_physics1
!
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATMOS_START_NORM'

! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_do_rims
LOGICAL :: L_do_halos
LOGICAL :: L_print_pe   ! print diagnostics on all pe's if true

INTEGER, INTENT(IN) ::  level_start, level_end, rims

! Data arrays

REAL, INTENT (IN) :: exner                                        &
                          (pdims_s%i_start:pdims_s%i_end,         &
                           pdims_s%j_start:pdims_s%j_end,         &
                           pdims_s%k_start:pdims_s%k_end + 1)
REAL, INTENT (IN) :: rho                                          &
                        (pdims_s%i_start:pdims_s%i_end,           &
                         pdims_s%j_start:pdims_s%j_end,           &
                         pdims_s%k_start:pdims_s%k_end)
REAL, INTENT (IN) :: u                                            &
                        (udims_s%i_start:udims_s%i_end,           &
                         udims_s%j_start:udims_s%j_end,           &
                         udims_s%k_start:udims_s%k_end)
REAL, INTENT (IN) :: v                                            &
                        (vdims_s%i_start:vdims_s%i_end,           &
                         vdims_s%j_start:vdims_s%j_end,           &
                         vdims_s%k_start:vdims_s%k_end)
REAL, INTENT(IN)  :: w                                            &
                        (wdims_s%i_start:wdims_s%i_end,           &
                         wdims_s%j_start:wdims_s%j_end,           &
                         wdims_s%k_start:wdims_s%k_end)
REAL, INTENT (IN) :: theta                                        &
                        (tdims_s%i_start:tdims_s%i_end,           &
                         tdims_s%j_start:tdims_s%j_end,           &
                         tdims_s%k_start:tdims_s%k_end)
REAL, INTENT (IN) :: q  (tdims_l%i_start:tdims_l%i_end,           &
                         tdims_l%j_start:tdims_l%j_end,           &
                         tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (IN) :: qcl(tdims_l%i_start:tdims_l%i_end,           &
                         tdims_l%j_start:tdims_l%j_end,           &
                         tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (IN) :: qcf(tdims_l%i_start:tdims_l%i_end,           &
                         tdims_l%j_start:tdims_l%j_end,           &
                         tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (IN) :: qcf2                                         &
                        (tdims_l%i_start:tdims_l%i_end,           &
                         tdims_l%j_start:tdims_l%j_end,           &
                         tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (IN) :: qrain                                        &
                        (tdims_l%i_start:tdims_l%i_end,           &
                         tdims_l%j_start:tdims_l%j_end,           &
                         tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (IN) :: qgraup                                       &
                           (tdims_l%i_start:tdims_l%i_end,        &
                            tdims_l%j_start:tdims_l%j_end,        &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (IN) :: cf_bulk                                      & 
                           (tdims_l%i_start:tdims_l%i_end,        &
                            tdims_l%j_start:tdims_l%j_end,        &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (IN) :: cfl                                          &
                           (tdims_l%i_start:tdims_l%i_end,        &
                            tdims_l%j_start:tdims_l%j_end,        &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (IN) :: cff                                          &
                           (tdims_l%i_start:tdims_l%i_end,        &
                            tdims_l%j_start:tdims_l%j_end,        &
                            tdims_l%k_start:tdims_l%k_end)

! === End of arguments ==============================================

! Local variables

INTEGER :: k, num
INTEGER :: levels, levelp
INTEGER :: start_level, end_level, end_press
INTEGER :: k0, s0      !  level zero pointers
INTEGER :: rim_n, rim_s, rim_e, rim_w

! Local array for storing two_norms  (sized for star and inc fields)
REAL ::  two_norm(15)

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
end_press = level_end + 1
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
                    exner, row_length, rows,                       &
                    levelp, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    rho, row_length, rows,                         &
                    model_levels, k, k,                            &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    u, row_length, rows,                           &
                    model_levels, k, k,                            &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    v, row_length, n_rows,                         &
                    model_levels, k, k,                            &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    w(1-offx,1-offy,s0), row_length, rows,         &
                    levels, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    theta(1-offx,1-offy,s0), row_length, rows,     &
                    levels, k, k,                                  &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    q(1-halo_i,1-halo_j,s0), row_length, rows,     &
                    levels, k, k,                                  &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    qcl(1-halo_i,1-halo_j,s0), row_length, rows,   &
                    levels, k, k,                                  &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    qcf(1-halo_i,1-halo_j,s0), row_length, rows,   &
                    levels, k, k,                                  &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    IF ( L_mcr_qrain ) THEN
      num = num + 1
      CALL array_norm(                                                 &
                      qrain(1-halo_i,1-halo_j,s0), row_length, rows,   &
                      levels, k, k,                                    &
                      halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,     &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF ! L_mcr_qrain
    IF ( L_mcr_qcf2 ) THEN
      num = num + 1
      CALL array_norm(                                             &
                    qcf2(1-halo_i,1-halo_j,s0), row_length, rows,  &
                    levels, k, k,                                  &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! L_mcr_qcf2
    IF ( L_mcr_qgraup ) THEN
      num = num + 1
      CALL array_norm(                                                &
                      qgraup(1-halo_i,1-halo_j,s0), row_length, rows, &
                      levels, k, k,                                   &
                      halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,    &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! L_mcr_qgraup
    IF ( i_cld_vn == i_cld_pc2 ) THEN
      num = num + 1
      CALL array_norm(                                               &
                    cf_bulk(1-halo_i,1-halo_j,s0), row_length, rows, &
                    levels, k, k,                                    &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,     &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
      num = num + 1
      CALL array_norm(                                             &
                    cfl(1-halo_i,1-halo_j,s0), row_length, rows,   &
                    levels, k, k,                                  &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
      num = num + 1
      CALL array_norm(                                             &
                    cff(1-halo_i,1-halo_j,s0), row_length, rows,   &
                    levels, k, k,                                  &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! i_cld_vn == i_cld_pc2

    IF ( L_print_pe .OR. mype == 0 ) THEN
      num = 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k, ' exner Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
         ' Level ', k, ' wet rho_r_sq_n Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k, ' u Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k, ' v Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k0, ' w Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k0, ' theta Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k0, ' q Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k0, ' qcl Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k0, ' qcf Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
      IF ( L_mcr_qrain ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' qrain Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_START_NORM')
      END IF ! L_mcr_qrain
      IF ( L_mcr_qcf2 ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' qcf2 Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_START_NORM')
      END IF ! L_mcr_qcf2
      IF ( L_mcr_qgraup ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' qgraup Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_START_NORM')
      END IF ! L_mcr_qgraup
      IF ( i_cld_vn == i_cld_pc2 ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' cf_bulk Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_START_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' cfl Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_START_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' cff Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_START_NORM')
      END IF ! i_cld_vn == i_cld_pc2

    END IF ! L_print_pe .OR. mype == 0

  END DO !  k = 1, model_levels

ELSE   !  3d norms from start_level to end_level  

  num = 1
  CALL array_norm(                                               &
                  exner, row_length, rows,                       &
                  levelp, start_level, end_press,                &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  rho, row_length, rows,                         &
                  model_levels, start_level, level_end,          &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  u, row_length, rows,                           &
                  model_levels, start_level, level_end,          &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                               &
                  v, row_length, n_rows,                         &
                  model_levels, start_level, level_end,          &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                               &
                  w(1-offx,1-offy,s0), row_length, rows,         &
                  levels, start_level, level_end,                &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                               &
                  theta(1-offx,1-offy,s0), row_length, rows,     &
                  levels, start_level, level_end,                &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  q(1-halo_i,1-halo_j,1), row_length, rows,      &
                  levels, start_level, level_end,                &
                  halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  qcl(1-halo_i,1-halo_j,1), row_length, rows,    &
                  levels, start_level, level_end,                &
                  halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  qcf(1-halo_i,1-halo_j,1), row_length, rows,    &
                  levels, start_level, level_end,                &
                  halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  IF ( L_mcr_qrain ) THEN
    num = num + 1
    CALL array_norm(                                               &
                    qrain(1-halo_i,1-halo_j,1), row_length, rows,  &
                    levels, start_level, level_end,                &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF ! L_mcr_qrain
  IF ( L_mcr_qcf2 ) THEN
    num = num + 1
    CALL array_norm(                                               &
                    qcf2(1-halo_i,1-halo_j,1), row_length, rows,   &
                    levels, start_level, level_end,                &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! L_mcr_qcf2
  IF ( L_mcr_qgraup ) THEN
    num = num + 1
    CALL array_norm(                                               &
                    qgraup(1-halo_i,1-halo_j,1), row_length, rows, &
                    levels, start_level, level_end,                &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! L_mcr_qgraup
  IF ( i_cld_vn == i_cld_pc2 ) THEN
    num = num + 1
    CALL array_norm(                                                &
                    cf_bulk(1-halo_i,1-halo_j,1), row_length, rows, &
                    levels, level_start, level_end,                 &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,    &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                               &
                    cfl(1-halo_i,1-halo_j,1), row_length, rows,    &
                    levels, start_level, level_end,                &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                               &
                    cff(1-halo_i,1-halo_j,1), row_length, rows,    &
                    levels, start_level, level_end,                &
                    halo_i, halo_j, 0, 0, L_do_halos, L_do_rims,   &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! i_cld_vn == i_cld_pc2

  IF ( L_print_pe .OR. mype == 0 ) THEN
    num = 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', end_press,      &
                  ' exner Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_START_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' wet rho_r_sq_n Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_START_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' u Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_START_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' v Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_START_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' w Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_START_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' theta Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_START_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' q Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_START_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' qcl Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_START_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' qcf Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_START_NORM')
    IF ( L_mcr_qrain ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' qrain Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
    END IF ! L_mcr_qrain
    IF ( L_mcr_qcf2 ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' qcf2 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
    END IF ! L_mcr_qcf2
    IF ( L_mcr_qgraup ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' qgraup Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
    END IF ! L_mcr_qgraup
    IF ( i_cld_vn == i_cld_pc2 ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' cf_bulk Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                   ' Levels', start_level, ' to', level_end,     &
                   ' cfl Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' cff Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_START_NORM')
    END IF ! i_cld_vn == i_cld_pc2

  END IF ! L_print_pe .OR. mype == 0

END IF ! start_level == end_level

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE atmos_start_norm
END MODULE atmos_start_norm_mod


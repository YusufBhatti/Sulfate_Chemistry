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
MODULE atmos_phys2_norm_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ATMOS_PHYS2_NORM_MOD'

CONTAINS
SUBROUTINE atmos_phys2_norm(                                            &
                            level_start, level_end, rims,               &
                            theta_star, q_star, qcl_star, qcf_star,     &
                            qrain_star, qcf2_star, qgraup_star,         &
                            cf_star, cfl_star, cff_star, R_u, R_v, R_w, &
                            L_do_rims, L_do_halos, L_print_pe )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE mphys_inputs_mod,  ONLY: l_mcr_qrain, l_mcr_qcf2, l_mcr_qgraup
USE cloud_inputs_mod,  ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2
USE atm_fields_bounds_mod, ONLY: tdims, tdims_s, udims_s, vdims_s, wdims
USE array_norm_mod, ONLY: array_norm
USE pr_block4_mod, ONLY: dom_no, dom_eo, dom_s, dom_w
USE nlsizes_namelist_mod,  ONLY: row_length, rows,n_rows, model_levels, &
                                 global_row_length, global_rows
USE UM_parvars,  ONLY: offx, offy
USE um_parcore, ONLY: mype

USE turb_diff_mod, ONLY: l_flush6
USE umPrintMgr, ONLY: umPrint, umMessage, umPrintFlush

!
! Description:  Printing of norms after atmos_physics2
!
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATMOS_PHYS2_NORM'

! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_do_rims
LOGICAL :: L_do_halos
LOGICAL :: L_print_pe   ! print diagnostics on all pe's if true

INTEGER, INTENT(IN) ::  level_start, level_end, rims

! Data arrays

REAL, INTENT (IN) :: theta_star                                &
                        (tdims%i_start:tdims%i_end,            &
                         tdims%j_start:tdims%j_end,            &
                         tdims%k_start:tdims%k_end)
REAL, INTENT (IN) :: q_star                                    &
                        (tdims%i_start:tdims%i_end,            &
                         tdims%j_start:tdims%j_end,            &
                         tdims%k_start:tdims%k_end)
REAL, INTENT (IN) :: qcl_star                                  &
                        (tdims%i_start:tdims%i_end,            &
                         tdims%j_start:tdims%j_end,            &
                         tdims%k_start:tdims%k_end)
REAL, INTENT (IN) :: qcf_star                                  &
                        (tdims%i_start:tdims%i_end,            &
                         tdims%j_start:tdims%j_end,            &
                         tdims%k_start:tdims%k_end)
REAL, INTENT (IN) :: qrain_star                                &
                        (tdims%i_start:tdims%i_end,            &
                         tdims%j_start:tdims%j_end,            &
                         tdims%k_start:tdims%k_end)
REAL, INTENT (IN) :: qgraup_star                               &
                        (tdims%i_start:tdims%i_end,            &
                         tdims%j_start:tdims%j_end,            &
                         tdims%k_start:tdims%k_end)
REAL, INTENT (IN) :: qcf2_star                                 &
                           (tdims%i_start:tdims%i_end,         &
                            tdims%j_start:tdims%j_end,         &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (IN) :: cf_star                                   &
                        (tdims%i_start:tdims%i_end,            &
                         tdims%j_start:tdims%j_end,            &
                         tdims%k_start:tdims%k_end)
REAL, INTENT (IN) :: cfl_star                                  &
                        (tdims%i_start:tdims%i_end,            &
                         tdims%j_start:tdims%j_end,            &
                         tdims%k_start:tdims%k_end)
REAL, INTENT (IN) :: cff_star                                  &
                           (tdims%i_start:tdims%i_end,         &
                            tdims%j_start:tdims%j_end,         &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (IN) :: R_u                                       &
                        (udims_s%i_start:udims_s%i_end,        &
                         udims_s%j_start:udims_s%j_end,        &
                         udims_s%k_start:udims_s%k_end)
REAL, INTENT (IN) :: R_v                                       &
                        (vdims_s%i_start:vdims_s%i_end,        &
                         vdims_s%j_start:vdims_s%j_end,        &
                         vdims_s%k_start:vdims_s%k_end)
REAL, INTENT (IN) :: R_w                                       &
                        (wdims%i_start:wdims%i_end,            &
                         wdims%j_start:wdims%j_end,            &
                         wdims%k_start:wdims%k_end)

! === End of arguments ==============================================

! Local variables

INTEGER :: k, num
INTEGER :: levels      !  model_levels + 1 
INTEGER :: start_level, end_level
INTEGER :: k0, s0      !  level zero pointers
INTEGER :: rim_n, rim_s, rim_e, rim_w

! Local array for storing two_norms  (sized for star and R_uvw fields)
REAL ::  two_norm(13)

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
    CALL array_norm(                                                &
                    theta_star(1,1,s0), row_length, rows,           &
                    levels, k, k,                                   &
                    0, 0, 0, 0, L_do_halos, L_do_rims,              &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    q_star(1,1,s0), row_length, rows,               &
                    levels, k, k,                                   &
                    0, 0, 0, 0, L_do_halos, L_do_rims,              &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    qcl_star(1,1,s0), row_length, rows,             &
                    levels, k, k,                                   &
                    0, 0, 0, 0, L_do_halos, L_do_rims,              &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                                &
                    qcf_star(1,1,s0), row_length, rows,             &
                    levels, k, k,                                   &
                    0, 0, 0, 0, L_do_halos, L_do_rims,              &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    IF ( L_mcr_qrain ) THEN
      num = num + 1
      CALL array_norm(                                               &
                      qrain_star(1,1,s0), row_length, rows,          &
                      levels, k, k,                                  &
                      0, 0, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num))
    END IF ! L_mcr_qrain
    IF ( L_mcr_qcf2 ) THEN
      num = num + 1
      CALL array_norm(                                               &
                    qcf2_star(1,1,s0), row_length, rows,             &
                    levels, k, k,                                    &
                    0, 0, 0, 0, L_do_halos, L_do_rims,         &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! L_mcr_qcf2
    IF ( L_mcr_qgraup ) THEN
      num = num + 1
      CALL array_norm(                                               &
                      qgraup_star(1,1,s0), row_length, rows,         &
                      levels, k, k,                                  &
                      0, 0, 0, 0, L_do_halos, L_do_rims,       &
                      rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! L_mcr_qgraup
    IF ( i_cld_vn == i_cld_pc2 ) THEN
      num = num + 1
      CALL array_norm(                                               &
                    cf_star(1,1,s0), row_length, rows,               &
                    levels, k, k,                                    &
                    0, 0, 0, 0, L_do_halos, L_do_rims,         &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
      num = num + 1
      CALL array_norm(                                               &
                    cfl_star(1,1,s0), row_length, rows,              &
                    levels, k, k,                                    &
                    0, 0, 0, 0, L_do_halos, L_do_rims,         &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
      num = num + 1
      CALL array_norm(                                               &
                    cff_star(1,1,s0), row_length, rows,              &
                    levels, k, k,                                    &
                    0, 0, 0, 0, L_do_halos, L_do_rims,         &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    END IF ! i_cld_vn == i_cld_pc2
    num = num + 1
    CALL array_norm(                                               &
                    R_u, row_length, rows,                         &
                    model_levels, k, k,                            &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    R_v, row_length, n_rows,                       &
                    model_levels, k, k,                            &
                    offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
    num = num + 1
    CALL array_norm(                                               &
                    R_w(1,1,s0), row_length, rows,                 &
                    levels, k, k,                                  &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))

    IF ( L_print_pe .OR. mype == 0 ) THEN
      num = 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k0, ' theta_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k0, ' q_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k0, ' qcl_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                         &
              ' Level ', k0, ' qcf_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      IF ( L_mcr_qrain ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' qrain_star Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      END IF ! L_mcr_qrain
      IF ( L_mcr_qcf2 ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' qcf2_star Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      END IF ! L_mcr_qcf2
      IF ( L_mcr_qgraup ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' qgraup_star Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      END IF ! L_mcr_qgraup
      IF ( i_cld_vn == i_cld_pc2 ) THEN
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' cf_star Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' cfl_star Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
        num = num + 1
        WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' cff_star Two_Norm =' , two_norm(num)
        CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      END IF ! i_cld_vn == i_cld_pc2
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k, ' r_u_p2 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k, ' r_v_p2 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, E23.16)')                       &
              ' Level ', k0, ' r_w_p2 Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')

    END IF ! L_print_pe .OR. mype == 0

  END DO !  k = start_level, end_level

ELSE   !  3d norms from level_start to level_end

  num = 1
  CALL array_norm(                                               &
                  theta_star(1,1,1), row_length, rows,           &
                  levels, start_level, level_end,                &
                  0, 0, 0, 0, L_do_halos, L_do_rims,             &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  q_star(1, 1, 1), row_length, rows,             &
                  levels, start_level, level_end,                &
                  0, 0, 0, 0, L_do_halos, L_do_rims,             &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  qcl_star(1, 1, 1), row_length, rows,           &
                  levels, start_level, level_end,                &
                  0, 0, 0, 0, L_do_halos, L_do_rims,             &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  num = num + 1
  CALL array_norm(                                               &
                  qcf_star(1, 1, 1), row_length, rows,           &
                  levels, start_level, level_end,                &
                  0, 0, 0, 0, L_do_halos, L_do_rims,             &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num))
  IF ( L_mcr_qrain ) THEN
    num = num + 1
    CALL array_norm(                                               &
                    qrain_star(1, 1, 1), row_length, rows,         &
                    levels, start_level, level_end,                &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num))
  END IF ! L_mcr_qrain
  IF ( L_mcr_qcf2 ) THEN
    num = num + 1
    CALL array_norm(                                               &
                    qcf2_star(1, 1, 1), row_length, rows,          &
                    levels, start_level, level_end,                &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! L_mcr_qcf2
  IF ( L_mcr_qgraup ) THEN
    num = num + 1
    CALL array_norm(                                                &
                    qgraup_star(1, 1, 1), row_length, rows,         &
                    model_levels, start_level, level_end,           &
                    0, 0, 0, 0, L_do_halos, L_do_rims,              &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! L_mcr_qgraup
  IF ( i_cld_vn == i_cld_pc2 ) THEN
    num = num + 1
    CALL array_norm(                                               &
                    cf_star(1, 1, 1), row_length, rows,            &
                    levels, start_level, level_end,                &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                               &
                    cfl_star(1, 1, 1), row_length, rows,           &
                    levels, start_level, level_end,                &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
    num = num + 1
    CALL array_norm(                                               &
                    cff_star(1, 1, 1), row_length, rows,           &
                    levels, start_level, level_end,                &
                    0, 0, 0, 0, L_do_halos, L_do_rims,             &
                    rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  END IF ! i_cld_vn == i_cld_pc2
  num = num + 1
  CALL array_norm(                                               &
                  R_u, row_length, rows,                         &
                  model_levels, start_level, level_end,          &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                               &
                  R_v, row_length, n_rows,                       &
                  model_levels, start_level, level_end,          &
                  offx, offy, 0, 0, L_do_halos, L_do_rims,       &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )
  num = num + 1
  CALL array_norm(                                               &
                  R_w(1,1,1) , row_length, rows,                 &
                  levels, start_level, level_end,                &
                  0, 0, 0, 0, L_do_halos, L_do_rims,             &
                  rim_n, rim_s, rim_e, rim_w, two_norm(num) )

  IF ( L_print_pe .OR. mype == 0 ) THEN
    num = 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' theta_star Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' q_star Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' qcl_star Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' qcf_star Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
    IF ( L_mcr_qrain ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' qrain_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
    END IF ! L_mcr_qrain
    IF ( L_mcr_qcf2 ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' qcf2_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
    END IF ! L_mcr_qcf2
    IF ( L_mcr_qgraup ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' qgraup_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
    END IF ! L_mcr_qgraup
    IF ( i_cld_vn == i_cld_pc2 ) THEN
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' cf_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' cfl_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
      num = num + 1
      WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')              &
                  ' Levels', start_level, ' to', level_end,      &
                  ' cff_star Two_Norm =' , two_norm(num)
      CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
    END IF ! i_cld_vn == i_cld_pc2
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' r_u_p2 Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' r_v_p2 Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')
    num = num + 1
    WRITE(umMessage, '(A, I5, A, I5, A, E23.16)')                &
                  ' Levels', start_level, ' to', level_end,      &
                  ' r_w_p2 Two_Norm =' , two_norm(num)
    CALL umPrint(umMessage, src='ATMOS_PHYS2_NORM')

  END IF ! L_print_pe .OR. mype == 0

END IF ! level_start == level_end

IF ( L_flush6 ) CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE atmos_phys2_norm
END MODULE atmos_phys2_norm_mod


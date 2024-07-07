! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! subroutine print_L2nhelm
SUBROUTINE print_L2nhelm(                                         &
                         HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,      &
                         HM_Czz, HM_Cz, HM_Cxz, HM_Cyz,           &
                         HM_Cxp, HM_Cyp, HM_C2n,                  &
                         HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,      &
                         hm_c3, hm_c4, hm_c5, hm_rhs,             &
                         row_length, rows, model_levels,          &
                         first_constant_r_rho_level_m1,           &
                         start_level, end_level,                  &
                         off_x, off_y, me, n_procx, n_procy,      &
                         proc_row_group, proc_col_group,          &
                         global_row_length, global_rows,          &
                         at_extremity, l_datastart,               &
                         L_do_halos, L_do_rims, rims_to_do,       &
                         L_print_pe )

! Purpose:
!          To calculate and print l2norms of tendencies or star fields
!
! Method:
!          Is described in ;
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.
INTEGER, INTENT(IN) ::                                            &
  row_length                                                      &
                    ! in: no of points per local row
, rows                                                            &
                    ! in: no of local (theta) rows
, model_levels                                                    &
                    ! in: no of model levels
, off_x                                                           &
                     ! standard halo size in east-west
, off_y                                                           &
                     ! standard halo size in North-South
, rims_to_do                                                      &
, me                                                              &
                    ! processor id
, n_procx                                                         &
, n_procy                                                         &
, proc_row_group                                                  &
, proc_col_group                                                  &
, global_row_length                                               &
, global_rows                                                     &
, first_constant_r_rho_level_m1                                   &
                                ! value used to dimension
                                ! arrays, max of (1 and
                                ! first_constant_r_rho_level)
, start_level                                                     &
                    ! start level for norm calculation
, end_level                                                       &
                    ! end level for norm calculation
, l_datastart(3)

REAL, INTENT(IN) ::                                               &
 HM_Cxx1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
,HM_Cxx2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
,HM_Cxy1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
,HM_Cxy2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
,HM_Cyy1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
,HM_Cyy2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
,HM_Cyx1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
,HM_Cyx2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
,HM_Czz (1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
,HM_Cz (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
,HM_C2n (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
         first_constant_r_rho_level_m1)                           &
,hm_c3 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
,hm_c4 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
,hm_c5 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
         first_constant_r_rho_level_m1)                           &
, HM_Cxz (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
         first_constant_r_rho_level_m1)                           &
, HM_Cyz (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
         first_constant_r_rho_level_m1)                           &
, HM_Cxp (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
         first_constant_r_rho_level_m1)                           &
, HM_Cyp (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
         first_constant_r_rho_level_m1)                           &
, hm_rhs(row_length, rows, model_levels)

LOGICAL ::                                                        &
  L_do_halos                                                      &
, L_do_rims                                                       &
, L_print_pe      ! true if  printing on all pe's

LOGICAL ::                                                        &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid

!    local variables
INTEGER :: first_level
INTEGER :: last_level

REAL ::  Two_norm
LOGICAL ::  L_do

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_L2NHELM'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!    set pointers for arrays with levels < first_constant_r_rho_level
L_do = .TRUE.
first_level = start_level
last_level = end_level
IF ( start_level > first_constant_r_rho_level_m1) THEN
  L_do = .FALSE.
ELSE IF ( end_level > first_constant_r_rho_level_m1) THEN
  first_level = start_level
  last_level = first_constant_r_rho_level_m1
END IF ! start_level > first_constant_r_rho_level_m1

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     HM_Cxx1, row_length, rows, model_levels,     &
                     start_level, end_level,                      &
                     off_x, off_y, 1, 0,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'   ***   L2_Norms of solver coefficients  ****  '
  CALL umPrint(umMessage,src='print_l2nhelm')
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_Cxx1 Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     HM_Cxx2, row_length, rows, model_levels,     &
                     start_level, end_level,                      &
                     off_x, off_y, 1, 0,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_Cxx2 Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     HM_Cyy1, row_length, rows, model_levels,     &
                     start_level, end_level,                      &
                     off_x, off_y, 0, 1,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_Cyy1 Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     HM_Cyy2, row_length, rows, model_levels,     &
                     start_level, end_level,                      &
                     off_x, off_y, 0, 1,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_Cyy2 Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     HM_Czz, row_length, rows, model_levels,      &
                     start_level, end_level,                      &
                     off_x, off_y, 0, 0,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_Czz Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     HM_Cz, row_length, rows, model_levels,       &
                     start_level, end_level,                      &
                     off_x, off_y, 0, 0,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_Cz Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

IF (L_do) THEN
  ! DEPENDS ON: two_norm_levels
  CALL Two_Norm_levels(                                             &
                       HM_Cxz, row_length, rows,                    &
                       first_constant_r_rho_level_m1,               &
                       first_level, last_level,                     &
                       off_x, off_y, 1, 0,                          &
                       at_extremity, n_procx, n_procy,              &
                       global_row_length, global_rows,              &
                       proc_col_group, proc_row_group,              &
                       l_datastart, L_do_halos,                     &
                       L_do_rims, rims_to_do, Two_Norm )
  IF ( L_print_pe .OR. me == 0 ) THEN
    WRITE(umMessage,*)'Levels',first_level,' to',last_level,        &
           ' HM_Cxz Two_Norm = ' , Two_Norm
    CALL umPrint(umMessage,src='print_l2nhelm')
  END IF ! L_print_pe .or. me == 0

  ! DEPENDS ON: two_norm_levels
  CALL Two_Norm_levels(                                             &
                       HM_Cyz, row_length, rows,                    &
                       first_constant_r_rho_level_m1,               &
                       first_level, last_level,                     &
                       off_x, off_y, 0, 1,                          &
                       at_extremity, n_procx, n_procy,              &
                       global_row_length, global_rows,              &
                       proc_col_group, proc_row_group,              &
                       l_datastart, L_do_halos,                     &
                       L_do_rims, rims_to_do, Two_Norm )
  IF ( L_print_pe .OR. me == 0 ) THEN
    WRITE(umMessage,*)'Levels',first_level,' to',last_level,        &
           ' HM_Cyz Two_Norm = ' , Two_Norm
    CALL umPrint(umMessage,src='print_l2nhelm')
  END IF ! L_print_pe .or. me == 0

  ! DEPENDS ON: two_norm_levels
  CALL Two_Norm_levels(                                             &
                       HM_Cxp, row_length, rows,                    &
                       first_constant_r_rho_level_m1,               &
                       first_level, last_level,                     &
                       off_x, off_y, 1, 0,                          &
                       at_extremity, n_procx, n_procy,              &
                       global_row_length, global_rows,              &
                       proc_col_group, proc_row_group,              &
                       l_datastart, L_do_halos,                     &
                       L_do_rims, rims_to_do, Two_Norm )
  IF ( L_print_pe .OR. me == 0 ) THEN
    WRITE(umMessage,*)'Levels',first_level,' to',last_level,        &
           ' HM_Cxp Two_Norm = ' , Two_Norm
    CALL umPrint(umMessage,src='print_l2nhelm')
  END IF ! L_print_pe .or. me == 0

  ! DEPENDS ON: two_norm_levels
  CALL Two_Norm_levels(                                             &
                       HM_Cyp, row_length, rows,                    &
                       first_constant_r_rho_level_m1,               &
                       first_level, last_level,                     &
                       off_x, off_y, 0, 1,                          &
                       at_extremity, n_procx, n_procy,              &
                       global_row_length, global_rows,              &
                       proc_col_group, proc_row_group,              &
                       l_datastart, L_do_halos,                     &
                       L_do_rims, rims_to_do, Two_Norm )
  IF ( L_print_pe .OR. me == 0 ) THEN
    WRITE(umMessage,*)'Levels',first_level,' to',last_level,        &
           ' HM_Cyp Two_Norm = ' , Two_Norm
    CALL umPrint(umMessage,src='print_l2nhelm')
  END IF ! L_print_pe .or. me == 0
END IF !  L_do

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     HM_Cxy1, row_length, rows, model_levels,     &
                     start_level, end_level,                      &
                     off_x, off_y, 1, 0,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_Cxy1 Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     HM_Cxy2, row_length, rows, model_levels,     &
                     start_level, end_level,                      &
                     off_x, off_y, 0, 1,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_Cxy2 Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     HM_Cyx1, row_length, rows, model_levels,     &
                     start_level, end_level,                      &
                     off_x, off_y, 0, 1,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_Cyx1 Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     HM_Cyx2, row_length, rows, model_levels,     &
                     start_level, end_level,                      &
                     off_x, off_y, 1, 0,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_Cyx2 Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

IF (L_do) THEN
  ! DEPENDS ON: two_norm_levels
  CALL Two_Norm_levels(                                             &
                       HM_C2n, row_length, rows,                    &
                       first_constant_r_rho_level_m1,               &
                       first_level, last_level,                     &
                       off_x, off_y, 0, 0,                          &
                       at_extremity, n_procx, n_procy,              &
                       global_row_length, global_rows,              &
                       proc_col_group, proc_row_group,              &
                       l_datastart, L_do_halos,                     &
                       L_do_rims, rims_to_do, Two_Norm )
  IF ( L_print_pe .OR. me == 0 ) THEN
    WRITE(umMessage,*)'Levels',first_level,' to',last_level,        &
           ' HM_C2n Two_Norm = ' , Two_Norm
    CALL umPrint(umMessage,src='print_l2nhelm')
  END IF ! L_print_pe .or. me == 0
END IF !  L_do

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     hm_c3, row_length, rows, model_levels,       &
                     start_level, end_level,                      &
                     off_x, off_y, 0, 0,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_C3 Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     hm_c4, row_length, rows, model_levels,       &
                     start_level, end_level,                      &
                     off_x, off_y, 0, 0,                          &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_C4 Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

IF (L_do) THEN
  ! DEPENDS ON: two_norm_levels
  CALL Two_Norm_levels(                                             &
                       hm_c5, row_length, rows,                     &
                       first_constant_r_rho_level_m1,               &
                       first_level, last_level,                     &
                       off_x, off_y, 0, 0,                          &
                       at_extremity, n_procx, n_procy,              &
                       global_row_length, global_rows,              &
                       proc_col_group, proc_row_group,              &
                       l_datastart, L_do_halos,                     &
                       L_do_rims, rims_to_do, Two_Norm )
  IF ( L_print_pe .OR. me == 0 ) THEN
    WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
           ' HM_C5 Two_Norm = ' , Two_Norm
    CALL umPrint(umMessage,src='print_l2nhelm')
  END IF ! L_print_pe .or. me == 0
END IF !  L_do

! DEPENDS ON: two_norm_levels
CALL Two_Norm_levels(                                             &
                     hm_rhs, row_length, rows, model_levels,      &
                     start_level, end_level,                      &
                     0, 0, 0, 0,                                  &
                     at_extremity, n_procx, n_procy,              &
                     global_row_length, global_rows,              &
                     proc_col_group, proc_row_group,              &
                     l_datastart, L_do_halos,                     &
                     L_do_rims, rims_to_do, Two_Norm )
IF ( L_print_pe .OR. me == 0 ) THEN
  WRITE(umMessage,*)'Levels',start_level,' to',end_level,         &
         ' HM_RHS Two_Norm = ' , Two_Norm
  CALL umPrint(umMessage,src='print_l2nhelm')
END IF ! L_print_pe .or. me == 0

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE print_L2nhelm


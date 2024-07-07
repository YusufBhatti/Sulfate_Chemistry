! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Two_Norm_levels

SUBROUTINE Two_Norm_levels(                                       &
                          field, row_length, rows, levels,        &
                          start_level, end_level,                 &
                          halo_i, halo_j, off_u, off_v,           &
                          at_extremity, n_procx, n_procy,         &
                          global_row_length, global_rows,         &
                          gc_proc_col_group, gc_proc_row_group,   &
                          l_datastart, L_do_halos,                &
                          L_do_rims, rims_to_do, Two_Norm )

! Purpose:
!          Calculates the two norm over a range of levels of the
!          input field.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE global_2d_sums_mod, ONLY: global_2d_sums
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE model_domain_mod, ONLY: model_type, mt_global, mt_lam, mt_cyclic_lam

IMPLICIT NONE


! Arguments with Intent IN. ie: Input variables.

LOGICAL ::                                                        &
  at_extremity(4)                                                 &
                   ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid
, L_do_rims                                                       &
                 ! LAMs: include rims, Global do all polar points
, L_do_halos     ! include halos in calculation of norm
!                    NB this means that some points are counted twice

INTEGER ::                                                        &
  row_length                                                      &
                   ! number of point on a row.
, rows                                                            &
                   ! number of rows.
, levels                                                          &
                   ! number of levels in field.
, start_level                                                     &
                   ! start_level for norm calculation
, end_level        ! end_level for norm calculation

INTEGER ::                                                        &
  halo_i                                                          &
, halo_j                                                          &
, off_u                                                           &
                   ! = 1 if field is at u-point
, off_v                                                           &
                   ! = 1 if field is at v-point
, rims_to_do                                                      &
, l_datastart(3)                                                  &
, global_row_length                                               &
, global_rows                                                     &
, gc_proc_col_group                                               &
, gc_proc_row_group                                               &
, n_procx                                                         &
, n_procy

REAL ::                                                           &
  field(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j, levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

! Arguments with Intent OUT. ie: variables Output only

REAL ::                                                           &
  Two_Norm

! Local Variables.

INTEGER ::                                                        &
  i, j, k                                                         &
, i_start                                                         &
, i_stop                                                          &
, j_start                                                         &
, j_stop                                                          &
, rims2                                                           &
, ipoints                                                         &
, istat

REAL ::                                                           &
  points           ! number of points used in norm

! Local arrays for parallel code

REAL ::                                                           &
  two_norm_component( row_length+2*halo_i, rows+2*halo_j )        &
, two_norm_temp( 1 )

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TWO_NORM_LEVELS'




!  External Routines:
! ----------------------------------------------------------------------
! Section 1.   Calculate Norm.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

rims2 = 2* rims_to_do
i_start = 1
i_stop = row_length
j_start = 1
j_stop = rows

IF ( L_do_rims ) THEN

  ipoints = global_row_length * global_rows

  IF ( L_do_halos ) THEN
    i_start = 1 - halo_i
    i_stop = row_length + halo_i
    j_start = 1 - halo_j
    j_stop = rows + halo_j
    ipoints = (global_row_length + 2 * n_procx * halo_i ) *       &
              (global_rows + 2 * n_procy * halo_j )
  END IF ! L_do_halos

  IF ( model_type == mt_global ) THEN
    IF (at_extremity(PNorth)) j_stop = rows - off_v
    ipoints = global_row_length * (global_rows - off_v)
  END IF ! model_type == mt_global

  IF ( model_type == mt_lam ) THEN
    IF (at_extremity(PNorth)) j_stop = rows - off_v
    IF (at_extremity(PEast)) i_stop = row_length - off_u
    ipoints = (global_row_length - off_u) * (global_rows - off_v)
  END IF ! model_type == mt_lam

ELSE  ! do not do halos or rims

  IF ( L_do_halos ) THEN
    WRITE(umMessage,*)' You must do rims if halos are included in L2norms'
    CALL umPrint(umMessage,src='two_norm_levels')
  END IF ! L_do_halos

  IF ( model_type == mt_global ) THEN
    IF (at_extremity(PSouth)) j_start = 2
    IF (at_extremity(PNorth)) j_stop = rows - 1
    ipoints = global_row_length * (global_rows - 2) + 2
  END IF ! model_type == mt_global

  IF ( model_type == mt_lam ) THEN
    IF (at_extremity(PSouth)) j_start = 1 + rims_to_do
    IF (at_extremity(PNorth)) j_stop = rows - rims_to_do
    IF (at_extremity(PEast)) i_stop = row_length - rims_to_do
    IF (at_extremity(PWest)) i_start = 1 + rims_to_do
    ipoints = (global_row_length - rims2) *                       &
              (global_rows  - rims2)

  ELSE IF (model_type == mt_cyclic_lam) THEN
    IF (at_extremity(PSouth)) j_start = 1 + rims_to_do
    IF (at_extremity(PNorth)) j_stop = rows - rims_to_do
    ipoints = global_row_length * (global_rows - rims2)

  END IF ! model_type == mt_lam

END IF ! L_do_rims

points = REAL(ipoints)

two_norm_component = 0.0

IF ( model_type == mt_global .AND. .NOT. L_do_rims) THEN
  ! Global model only calculate norm for one of the polar points
  !  if L_do_rims is .false.

  IF (at_extremity(PSouth) .AND. (l_datastart(1) == 1)) THEN
    DO k = start_level, end_level
      two_norm_component(1,1)= two_norm_component(1,1) +          &
                               field(1,1,k) * field(1,1,k)
    END DO
  END IF
  IF (at_extremity(PNorth) .AND. (l_datastart(1) == 1)) THEN
    DO k = start_level, end_level
      two_norm_component(1,rows)= two_norm_component(1,rows) +    &
                             field(1,rows,k) * field(1,rows,k)
    END DO
  END IF
END IF  !  model_type == mt_global .and. .not. L_do_rims

DO k = start_level, end_level
  DO j = j_start , j_stop
    DO i = i_start, i_stop
      two_norm_component( i + halo_i, j + halo_j ) =              &
          two_norm_component( i + halo_i, j + halo_j ) +          &
                                 field(i,j,k) * field(i,j,k)
    END DO
  END DO
END DO  ! k = start_level, end_level

! note halos included in sum.
CALL global_2d_sums(two_norm_component, row_length+2*halo_i,      &
                    rows + 2*halo_j, 0, 0, 1, two_norm_temp)

! Normalize by dividing by number of points

Two_Norm = SQRT (Two_Norm_temp(1) / points)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Two_Norm_levels

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_ideal_random_perturbation_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_IDEAL_RANDOM_PERTURBATION_MOD'

CONTAINS

! Subroutine rcf_ideal_random_perturbation
!
! Description:
!   Adds a random number field to any prognostic field defined on p-levels,
!   dependent on options specified in the idealised namelist. Could be extended
!   to work with other fields or add uncorrelated perturbations.
!
! Method:
!   Uses RANDOM_NUMBER routine to generate a random number field from
!   0 to 1 for the global field on pe0, then scatters it across all
!   processors. The same random number field is applied on each call (the seed
!   is reset each time), meaning perturbations to different fields are
!   implicitly correlated. An independent scaling factor is applied to the 
!   perturbations for each field.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

SUBROUTINE rcf_ideal_random_perturbation ( hdr_out, fields_out,              &
                                           field_count_out, stashcode,       &
                                           perturb_magnitude, zero)

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type
                                 
USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_Locate_Mod, ONLY: &
    rcf_locate

USE Rcf_Alloc_Field_Mod, ONLY: &
    rcf_alloc_field,           &
    rcf_dealloc_field

USE Rcf_Read_Field_Mod, ONLY: &
    rcf_read_field

USE Rcf_Write_Field_Mod, ONLY: &
    rcf_write_field

USE rcf_scatter_field_mod, ONLY: &
    rcf_scatter_field_real

USE decomp_params, ONLY: &
    decomp_rcf_output

USE um_parvars, ONLY: &
    gc_all_proc_group

USE rcf_nlist_recon_idealised_mod, ONLY: &
    perturb_height,                      &
    perturb_type

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec

USE umPrintMgr, ONLY: &
    newline,          &
    umPrint,          &
    umMessage,        &
    PrintStatus,      &
    PrStatus_Min,     &
    PrStatus_Normal

USE UM_ParCore, ONLY: &
    mype

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE(um_header_type), INTENT(IN) :: hdr_out      ! Output dump header
TYPE(field_type), POINTER :: fields_out(:)       ! Output fields
INTEGER, INTENT(IN)       :: field_count_out     ! No. of output fields
INTEGER, INTENT(IN)       :: stashcode           ! Field to be perturbed
REAL, INTENT(IN)          :: perturb_magnitude   ! Perturbation scale factor
LOGICAL, INTENT(IN), OPTIONAL :: zero            ! Flag to indicate zero-level
                                                 !  fields. Set to TRUE when
                                                 !  field has a zeroth level.

! Local variables
INTEGER :: i
INTEGER :: k
INTEGER :: pos     ! Index for locating a field
INTEGER :: lev_bot ! Bottom level corresponding to perturb_height(1)
INTEGER :: lev_top ! Top level corresponding to perturb_height(2)
INTEGER :: ransize ! Size of the random seed array
INTEGER, ALLOCATABLE :: iseed(:) ! Random seed array

REAL :: r_ref_theta(0:output_grid % model_levels)  ! Heights on theta levels
REAL :: glob_random_field(output_grid % glob_p_field) ! Global random array
REAL :: loc_random_field(output_grid % loc_p_field) ! Local random array

LOGICAL :: zero_lev   ! Set by the optional argument 'zero'.

! Pointer to output field:
TYPE( field_type ), POINTER :: field

CHARACTER (LEN=*),  PARAMETER :: RoutineName='RCF_IDEAL_RANDOM_PERTURBATION'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Process optional argument:
IF (PRESENT(zero)) THEN
  zero_lev = zero
ELSE
  zero_lev = .FALSE.
END IF

! Extract and read a field:
CALL rcf_locate(stashcode_prog_sec, stashcode, fields_out, field_count_out, pos)
field => fields_out(pos)
CALL rcf_alloc_field( field )
CALL rcf_read_field(field, hdr_out, decomp_rcf_output)

! Create array of random numbers on pe0.
IF (mype == 0) THEN
  CALL umPrint('Applying random perturbation to initial data',               &
               src='rcf_ideal_random_perturbation', level=PrStatus_Min)

  ! Allocate and populate the seed array
  CALL RANDOM_SEED(ransize)
  ALLOCATE(iseed(ransize))
  ! Set the seed array to the same values on each pass through the routine.
  DO i = 1, ransize
    iseed(i) = 0
  END DO
  CALL RANDOM_SEED(put=iseed)
  ! Generate random numbers for the global grid
  CALL RANDOM_NUMBER(glob_random_field)
  DEALLOCATE(iseed)

  ! Modify the random number field, if desired.
  ! (Placeholder functionality for later expansion.)
  SELECT CASE(perturb_type)

  CASE(1)
    IF (PrintStatus >= PrStatus_Min) THEN
      WRITE(umMessage,'(A)') 'System-generated random field with no '        &
               // 'additional perturbations selected.'
      CALL umPrint(umMessage, src='rcf_ideal_random_perturbation')
    END IF
    

  CASE DEFAULT
    IF (PrintStatus >= PrStatus_Min) THEN
      WRITE(umMessage,'(A,I0)') 'Unrecognised perturbation type: '           &
               // 'perturb_type = ',perturb_type
      CALL umPrint(umMessage, src='rcf_ideal_random_perturbation')
      WRITE(umMessage,'(A)') 'No additional perturbations applied to '       &
               // 'system-generated random field.'
      CALL umPrint(umMessage, src='rcf_ideal_random_perturbation')
    END IF
  END SELECT

  IF (PrintStatus >= PrStatus_Normal) THEN
    WRITE(umMessage,'(A,F6.3,A,F6.3)') 'Random field to be applied with '    &
               // 'minimum of ',MINVAL(glob_random_field), ' and '           &
               // 'maximum of ',MAXVAL(glob_random_field)
    CALL umPrint(umMessage, src='rcf_ideal_random_perturbation')
  END IF
END IF  ! PE0

! Scatter random number array.
CALL rcf_scatter_field_real(loc_random_field, glob_random_field,             &
                            output_grid % loc_p_row_length,                  &
                            output_grid % loc_p_rows,                        &
                            output_grid % glob_p_row_length,                 &
                            output_grid % glob_p_rows, 0,                    &
                            gc_all_proc_group )

! Perturbations are applied to all levels within the height (not level)
! range specified by perturb_height(1):perturb_height(2).

DO k = 0, output_grid % model_levels
  r_ref_theta(k) = output_grid % eta_theta_levels(k)                         &
                   * output_grid % z_top_of_model
END DO

! If maximum perturbation height is below the first (zero) model level,
! reset it (including a margin for subsequent tests on it).
IF (perturb_height(2) < r_ref_theta(0)) THEN
  perturb_height(2) = r_ref_theta(0) + EPSILON(1.0)
END IF

! Find bottom and top levels corresponding to the given heights.
lev_bot = -1
lev_top = 0
DO k = 0, output_grid % model_levels
  IF (r_ref_theta(k) >= perturb_height(1) .AND.                              &
      r_ref_theta(k) <= perturb_height(2)) THEN
    IF (lev_bot == -1) THEN
      lev_bot = k
    END IF
    lev_top = k
  END IF
END DO

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
  WRITE(umMessage,'(2(A,E13.6,A,I0))') 'Applying perturbations from ',       &
          perturb_height(1), ' m = level ',lev_bot, newline // 'to ',        &
          perturb_height(2), ' m = level ',lev_top
  CALL umPrint(umMessage, src='rcf_ideal_random_perturbation')
END IF

! If the field has a zero level, it will be defined over 1:N+1 levels, not 0:N,
! so modify the perturbation levels accordingly:
IF (zero_lev) THEN
  lev_bot = lev_bot + 1
  lev_top = lev_top + 1
END IF

! --------------------------------
! Apply random perturbations.
! --------------------------------
DO k = lev_bot, lev_top
  DO i = 1, output_grid % loc_p_field
    field % data(i,k) = field % data(i,k)                                    &
                      + 2.0 * perturb_magnitude * (loc_random_field(i) - 0.5)
  END DO
END DO

! Write out and deallocate the modified field:
CALL rcf_write_field   (fields_out(pos), hdr_out, decomp_rcf_output)
CALL rcf_dealloc_field (fields_out(pos))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_ideal_random_perturbation
END MODULE rcf_ideal_random_perturbation_mod

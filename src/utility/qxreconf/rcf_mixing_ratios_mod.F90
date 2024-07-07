! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_mixing_ratios_mod

IMPLICIT NONE

REAL, ALLOCATABLE :: moisture_field_sum(:,:)
! Used to determine how many times the rcf_mixing_ratio
! routine will be called. Can therefore deallocate 
! moisture_field_sum on last call.
INTEGER, SAVE :: total_mixing_ratios_calls = 0
INTEGER, SAVE :: current_mixing_ratio_call = 0

!  Subroutine Rcf_Mixing_Ratios_Mod

! Description:
!   Derive mixing ratios from specific humidities. Note that
!   currently ENDGame requires the optional moisture fields to
!   exist. Even if they only contain 0's.

! Method:
!    Each mixing ratio is calculated by dividing the equivilent
!    standard moisture field by 1 minus sum of all other standard
!    moisture fields.
!
!    mixing_ratio  = moisture_field  /
!                               (1.0 - moisture_field_sum)

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_MIXING_RATIOS_MOD'

CONTAINS

SUBROUTINE rcf_mixing_ratios( fields_out, field_count_out,     &
                              mr_field_stash_code, hdr_out,    &
                              mixing_ratio, data_source)

USE decomp_params,       ONLY: &
    decomp_rcf_output

USE parkind1,  ONLY: &
    jprb,            &
    jpim

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,           &
    rcf_dealloc_field

USE rcf_data_source_Mod, ONLY: &
    data_source_type

USE items_nml_mod, ONLY: &
    field_calcs

USE rcf_field_type_mod,  ONLY: &
    field_type

USE rcf_grid_type_mod,   ONLY: &
    grid_type

USE rcf_locate_mod,      ONLY: &
    rcf_locate

USE rcf_read_field_mod,  ONLY: &
    rcf_read_field

USE rcf_umhead_mod,      ONLY: &
    um_header_type

USE submodel_mod,        ONLY: &
    atmos_im

USE um_stashcode_mod,    ONLY: &
    stashcode_q,               &
    stashcode_qcf,             &
    stashcode_qcl,             &
    stashcode_qcf2,            &
    stashcode_qrain,           &
    stashcode_qgraup,          &
    stashcode_mv,              &
    stashcode_mcl,             &
    stashcode_mcf,             &
    stashcode_mr,              &
    stashcode_mgr,             &
    stashcode_mcf2,            &
    stashcode_prog_sec

USE um_parcore,          ONLY: & 
    mype

USE umPrintMgr,          ONLY: &
    printstatus,               &
    prstatus_normal,           &
    umPrint,                   &
    umMessage

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( um_header_type), INTENT(IN) :: hdr_out
TYPE( field_type ), INTENT(INOUT), TARGET :: mixing_ratio
INTEGER, INTENT(IN)               :: field_count_out
INTEGER, INTENT(IN)               :: mr_field_stash_code
TYPE( data_source_type ), POINTER :: data_source( : )

! Moisture fields
TYPE( field_type ), POINTER       :: temp_field
TYPE( field_type ), POINTER       :: moisture_field

CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_MIXING_RATIOS'

INTEGER                           :: k,i ! loop index
INTEGER                           :: pos
INTEGER                           :: moisture_field_pos
INTEGER                           :: moisture_field_stash_code
INTEGER                           :: stashcode_list(6)

LOGICAL                           :: zero_ok = .FALSE.
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE (umMessage,FMT='(a,i5)') 'Deriving STASH code ',mr_field_stash_code
  CALL umPrint(umMessage,src='rcf_mixing_ratios')
END IF


! On first call calculate how many times this routine will be called. 
! Use this information to deallocate moisture_field_sum at end of last
! call.
current_mixing_ratio_call = current_mixing_ratio_call + 1
IF (current_mixing_ratio_call == 1) THEN
  DO i = 1, field_count_out
    IF (data_source( i ) % Source == Field_Calcs .AND. &
        fields_out( i ) % stashmaster % model == atmos_im .AND. &
        fields_out( i ) % stashmaster % section == stashcode_prog_sec) THEN
      SELECT CASE( fields_out(i) % stashmaster % item )
        CASE ( stashcode_mv,   stashcode_mcl,     &
               stashcode_mcf,  stashcode_mr,      &
               stashcode_mgr,  stashcode_mcf2 )
          total_mixing_ratios_calls = total_mixing_ratios_calls + 1
      END SELECT
    END IF
  END DO
END IF

! Only calculate the moisture field sum once on first pass through routine
IF (current_mixing_ratio_call == 1) THEN
  ! Need to sum up all moisture fields present in output dump
  ALLOCATE (moisture_field_sum(SIZE( mixing_ratio % DATA, 1 ), &
       SIZE( mixing_ratio % DATA, 2 ) ))
  moisture_field_sum(:,:) = 0.0
  ! Need to sum up all moisture fields present in output dump
  stashcode_list = (/stashcode_q,    stashcode_qcf,    stashcode_qcl,  &
                     stashcode_qcf2, stashcode_qrain, stashcode_qgraup/)
  DO i = 1, SIZE(stashcode_list)
    IF (stashcode_list(i) == stashcode_qcf2   .OR.                     &
        stashcode_list(i) == stashcode_qrain  .OR.                     &
        stashcode_list(i) == stashcode_qgraup) THEN
      ! for optional fields, set flag allowing rcf_locate not to find field
      zero_ok = .TRUE.
    ELSE
      ! for required fields, set flag causing rcf_locate to abort if field
      ! not found.
      zero_ok = .FALSE.
    END IF
    CALL rcf_locate(stashcode_prog_sec, stashcode_list(i),             &
                    fields_out, field_count_out, pos,                  &
                    zero_ok_arg=zero_ok)
    ! if rcf_locate didn't find the field, pos will be 0 - do nothing.
    IF (pos /= 0) THEN
      temp_field => fields_out(pos)
      CALL rcf_alloc_field( temp_field )
      CALL rcf_read_field(  temp_field, hdr_out, decomp_rcf_output )
      DO k = 1, mixing_ratio % levels
        moisture_field_sum(:,k) = moisture_field_sum(:,k) + &
                                  temp_field % DATA(:,k)
      END DO
      CALL rcf_dealloc_field( temp_field )
    END IF
  END DO
END IF

! Which mixing ratio are we calculating?
! Need to map standard moisture field to the mixing ratio moisture
! field
SELECT CASE (mr_field_stash_code)
  CASE (stashcode_mv)
    moisture_field_stash_code = stashcode_q
  CASE (stashcode_mcl)
    moisture_field_stash_code = stashcode_qcl
  CASE (stashcode_mcf)
    moisture_field_stash_code = stashcode_qcf
  CASE (stashcode_mr)
    moisture_field_stash_code = stashcode_qrain
  CASE (stashcode_mgr)
    moisture_field_stash_code = stashcode_qgraup
  CASE (stashcode_mcf2)
    moisture_field_stash_code = stashcode_qcf2
END SELECT

CALL rcf_locate(stashcode_prog_sec, moisture_field_stash_code,         &
                fields_out, field_count_out, moisture_field_pos,       &
                zero_ok_arg=.TRUE.)
IF (moisture_field_pos /= 0) THEN
  ! If the standard moisture field is present in output dump
  ! we can calculate the mixing ratio
  moisture_field => fields_out(moisture_field_pos)
  CALL rcf_alloc_field( moisture_field )
  CALL rcf_read_field(  moisture_field, hdr_out, decomp_rcf_output )
  
  DO k = 1, mixing_ratio % levels
    mixing_ratio % DATA(:,k) = moisture_field % DATA(:,k) / &
                               (1.0 - moisture_field_sum(:,k))
  END DO
  CALL rcf_dealloc_field( moisture_field )
ELSE
  ! The requested mixing ratio cannot be calculated as equivalent moisture
  ! field is not present in the output dump.  However ENDGame expects to
  ! the mixing ratio so set data field to zero.
  mixing_ratio % DATA(:,:) = 0.0
END IF

! Tidy up memory on final pass through routine
IF (current_mixing_ratio_call == total_mixing_ratios_calls) THEN
  IF (ALLOCATED(moisture_field_sum)) DEALLOCATE(moisture_field_sum)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_mixing_ratios
END MODULE rcf_mixing_ratios_mod

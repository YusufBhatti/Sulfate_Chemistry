! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE find_fields_to_interp_mod

USE datafile_mod,              ONLY: datafile_type
USE fieldsfile_mod,            ONLY: fieldsfile_type
USE lbc_output_control_mod,    ONLY: max_input_files, lbc_output_control_type,  &
                                     max_fields_per_file, max_stash_items
USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength
USE time_utils_mod,            ONLY: timedate_to_seconds
USE um_stashcode_mod,          ONLY: stashcode_dust1_mmr, stashcode_dust2_mmr,  &
                                     stashcode_dust3_mmr, stashcode_dust4_mmr,  &
                                     stashcode_dust5_mmr, stashcode_dust6_mmr,  &
                                     stashcode_orog
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Description:
!   Module containing class and routines to work out which fields to interpolate
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

LOGICAL, SAVE :: found_orography = .FALSE.
INTEGER, ALLOCATABLE :: fields_to_interp(:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FIND_FIELDS_TO_INTERP_MOD'

CONTAINS

!-------------------------------------------------------------------------------
SUBROUTINE find_fields_to_interp(input_file, lbc_output_control)

USE missing_data_mod, ONLY: imdi
USE transform_stash_order_mod, ONLY: transform_stash_order
USE umPrintMgr, ONLY: umPrint, umMessage, PrStatus_Diag, PrintStatus
IMPLICIT NONE
        
CLASS(datafile_type), INTENT(IN) :: input_file
CLASS(lbc_output_control_type), INTENT(IN) :: lbc_output_control
INTEGER :: fields_list(max_fields_per_file)

INTEGER :: i, j, k, l, m, p ! Loop counters
INTEGER, ALLOCATABLE :: field_nums(:)
INTEGER, ALLOCATABLE :: required_stash(:)
INTEGER :: start_timedate = imdi
INTEGER :: end_timedate = imdi
INTEGER :: data_timedate = imdi
INTEGER :: first_timedate = imdi
LOGICAL :: needed_stash, all_present
INTEGER :: num_dust_bins_in_file = 0
INTEGER :: num_fields_required
LOGICAL, ALLOCATABLE :: used_num(:)

CHARACTER(LEN=15) :: array_format

CHARACTER(LEN=*), PARAMETER :: routinename = 'FIND_FIELDS_TO_INTERP'
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise the list of fields to interpolate to "nothing to do"
fields_list = -99
    
! If the start_time and/or end_time are set in the namelist, convert the
! six-item arrays into an integer number of seconds
IF (.NOT. ALL(lbc_output_control%start_time == imdi)) THEN
  start_timedate = timedate_to_seconds(lbc_output_control%start_time)
END IF
IF (.NOT. ALL(lbc_output_control%end_time == imdi)) THEN
  end_timedate = timedate_to_seconds(lbc_output_control%end_time)
END IF

! If we're dealing with a fields file, count the number of dust bins, otherwise
! assume none
SELECT TYPE(input_file)
CLASS IS (fieldsfile_type)
  num_dust_bins_in_file = input_file%number_of_dust_bins_in_file()
  CLASS DEFAULT
  num_dust_bins_in_file = 0
END SELECT

! Calculate the total number of fields we expect for each STASH time    
IF (lbc_output_control%num_dust_fields > 0) THEN
  num_fields_required = lbc_output_control%num_stash_items + num_dust_bins_in_file
ELSE
  ! No dust required
  num_fields_required = lbc_output_control%num_stash_items
END IF

! If there are no fields to interpolate, generate an error
IF (num_fields_required == 0) THEN
  cmessage = 'No fields specified to interpolate.'
  icode = 40
  CALL ereport(routinename, icode, cmessage)
END IF

IF (lbc_output_control%num_dust_fields > 0 .AND.             &
    num_dust_bins_in_file == 0) THEN
    cmessage = 'Dust requested but not present in input file'
    icode = 41
    CALL ereport(routinename, icode, cmessage)
END IF

IF (input_file%num_fields > max_fields_per_file) THEN
  cmessage = 'Too many fields in input file - please increase parameter ' &
          // 'and recompile'
  icode = 31
  CALL ereport(routinename, icode, cmessage)
END IF

      
ALLOCATE(field_nums(num_fields_required))
ALLOCATE(required_stash(num_fields_required))
ALLOCATE(used_num(input_file%num_fields))
used_num = .FALSE.

! Populate a list of stash codes required - this is equal to the namelist
! variable stash_codes for non-dust fields, and then the available dust
! fields are appended to the list afterwards if required.
required_stash(1:lbc_output_control%num_stash_items) = lbc_output_control%stash_codes(1:lbc_output_control%num_stash_items)

IF (lbc_output_control%num_dust_fields > 0) THEN
  IF (num_dust_bins_in_file == 2) THEN
    required_stash(lbc_output_control%num_stash_items + 1) = stashcode_dust1_mmr
    required_stash(lbc_output_control%num_stash_items + 2) = stashcode_dust2_mmr
  ELSE IF (num_dust_bins_in_file == 6) THEN
    required_stash(lbc_output_control%num_stash_items + 1) = stashcode_dust1_mmr
    required_stash(lbc_output_control%num_stash_items + 2) = stashcode_dust2_mmr
    required_stash(lbc_output_control%num_stash_items + 3) = stashcode_dust3_mmr
    required_stash(lbc_output_control%num_stash_items + 4) = stashcode_dust4_mmr
    required_stash(lbc_output_control%num_stash_items + 5) = stashcode_dust5_mmr
    required_stash(lbc_output_control%num_stash_items + 6) = stashcode_dust6_mmr
  ELSE
    cmessage = 'Invalid number of dust bins in input file'
    icode = 42
    CALL ereport(routinename, icode, cmessage)
  END IF
END IF

! The UM requires the fields to be in a specific order, so we adjust the
! list of STASH codes to be in that order
CALL transform_stash_order(num_fields_required, required_stash)

WRITE(array_format, '(A,I0,A)') '(A,',num_fields_required, '(I0,1x))'

WRITE(umMessage,'(I0,A)') num_fields_required, &
                                          ' fields are required at each time'
CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
WRITE(umMessage,array_format) 'These are ',required_stash
CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')

! Initialise the number of fields to interpolate to be zero
p = 0

! Loop over all fields in input file
next_field: DO j = 1, input_file%num_fields
  field_nums = 0
  needed_stash = .FALSE.

  IF (used_num(j)) THEN
    ! Skip current field if already added to file on previous loop iteration
    CYCLE next_field
  END IF

  WRITE(umMessage,'(A,I0,A,I0,A,A,A)') 'Analysing field number ',j, ' (STASH=',&
               input_file%fields(j)%quantity_ident,'; ',                       &
               input_file%fields(j)%return_validity_time_string(), ')'
  CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')

  ! Test to see if this field matches a STASH code we're interested in
  DO l = 1, num_fields_required
    IF (input_file%fields(j)%has_stashcode(required_stash(l))) THEN
      IF (found_orography .AND. required_stash(l) == stashcode_orog) THEN
        field_nums(l) = -1
      ELSE
        field_nums(l) = j
      END IF
      needed_stash = .TRUE.
    END IF
          
  END DO

  ! If we're not interested in this STASH code, try the next field
  IF (.NOT. needed_stash) THEN
    WRITE(umMessage,'(A)') 'Field with this STASH not required in output file'
    CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
    CYCLE next_field
  END IF

  ! If they're set, check start_time and end_time for this field 
  data_timedate = timedate_to_seconds(input_file%fields(j)%validity_time)
  IF (start_timedate /= imdi) THEN
    IF (data_timedate < start_timedate) THEN
      WRITE(umMessage,'(A)') 'Field validity time before start time'
      CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
      CYCLE next_field
    END IF
  END IF
        
  IF (end_timedate /= imdi) THEN
    IF (data_timedate > end_timedate) THEN
      WRITE(umMessage,'(A)') 'Field validity time after end time'
      CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
      CYCLE next_field
    END IF
  END IF
  
  IF (lbc_output_control%time_interval /= imdi) THEN
    IF (first_timedate /= imdi .AND. MOD(data_timedate - first_timedate,     &
      lbc_output_control%time_interval) /= 0) THEN
      WRITE(umMessage,'(A)') 'Field is not at correct time interval'
      CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
      CYCLE
    END IF
  END IF
      
  ! If we only require a single field (stash_codes in the namelist has a single
  ! entry), set this to be the field to interpolate and go on to the next field
  IF (num_fields_required == 1) THEN
    ! Only a single field at this validity time
    p = p + 1
    fields_list(p) = j
    WRITE(umMessage,'(A)') 'Field is the only field required at this time.'// &
                           ' Added to file.'
    CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
  ELSE
    ! Otherwise look for remaining fields at this validity time
    ! Loop over the remainder of the fields in the file
    WRITE(umMessage,'(A)') 'Searching for other required fields at ' // &
                           'this validity time.'
    CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
    DO k = j+1, input_file%num_fields

      IF (used_num(k)) THEN
        CYCLE
      END IF

      WRITE(umMessage,'(A,I0,A,I0)') 'Field ',k, ' is STASH=',input_file%fields(k)%quantity_ident
      CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')

      DO l = 1, num_fields_required
        ! Test to see if the field matches a STASH code we're interested in,
        ! and that it has the same validity time as the original field, and
        ! that we haven't already added this field to the list of fields to
        ! interpolate
        IF (input_file%fields(k)%has_stashcode(required_stash(l))) THEN
          IF (ALL(input_file%fields(j)%validity_time  == input_file%fields(k)%validity_time)) THEN
            IF (.NOT. used_num(k)) THEN
              WRITE(umMessage,'(A,I0,A)') 'Field ',k, ' has correct STASH and validity time. Saving.'
              CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
              field_nums(l) = k
            END IF
          END IF
        END IF

        ! We only need orography once; if we've already found it, note that.
        IF (found_orography .AND. required_stash(l) == stashcode_orog) THEN
          field_nums(l) = -1
        END IF
      END DO

      ! Test to see if we have a complete set of fields for this time
      all_present = .TRUE.
      WRITE(umMessage,'(A)') 'Analysing saved fields.'
      CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
      DO m = 1, num_fields_required
        IF (field_nums(m) == 0) THEN
          WRITE(umMessage,'(A,I0)') 'Cannot find field with STASH ',required_stash(m)
          CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
          all_present = .FALSE.
        END IF
      END DO

      ! If we have a complete set of fields for this time, add them all to the
      ! list of fields to interpolate.
      IF (all_present) THEN
        DO l = 1, num_fields_required
          IF (required_stash(l) == stashcode_orog .AND. found_orography) THEN
            CYCLE
          END IF
          WRITE(umMessage,'(A,I0,A,I0)') 'Added field ',field_nums(l), ' with STASH ',required_stash(l)
          CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
          p = p + 1
          fields_list(p) = field_nums(l)
          used_num(field_nums(l)) = .TRUE.
        END DO
        found_orography = .TRUE.
        IF (first_timedate == imdi) THEN
          first_timedate = data_timedate
        END IF
        WRITE(umMessage,'(A,I0,A,I0)') 'Done. All required fields found at this validity time.'
        CALL umPrint(umMessage,level=PrStatus_Diag, src='find_fields_to_interp')
        CYCLE next_field
      END IF

    END DO
  END IF

END DO next_field

! Leave deallocated if there are no fields to interpolate in this file
IF (p > 0) THEN
  CALL assign_fields_to_interp(fields_list)
END IF

DEALLOCATE(used_num)
DEALLOCATE(required_stash)
DEALLOCATE(field_nums)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE find_fields_to_interp

!-------------------------------------------------------------------------------

SUBROUTINE assign_fields_to_interp(field_list)
! Given a list of fields, assign them to the array of field numbers to interpolate
IMPLICIT NONE
INTEGER, INTENT(IN) :: field_list(max_fields_per_file)
INTEGER :: num
INTEGER :: i_field
CHARACTER(LEN=*), PARAMETER :: routinename = 'ASSIGN_FIELDS_TO_INTERP'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

num = max_fields_per_file
DO i_field = 1, max_fields_per_file
  IF (field_list(i_field) == -99) THEN
    num = i_field - 1
    EXIT
  END IF
END DO
  
IF (ALLOCATED(fields_to_interp)) THEN
  DEALLOCATE(fields_to_interp)
END IF

ALLOCATE(fields_to_interp(num))
fields_to_interp = -99
fields_to_interp(1:num) = field_list(1:num)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE assign_fields_to_interp

!-------------------------------------------------------------------------------

SUBROUTINE clear_fields_to_interp()
! Clear the list of fields to interpolate in preparation for a new input file
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER :: routinename = 'CLEAR_FIELDS_TO_INTERP'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DEALLOCATE(fields_to_interp)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE clear_fields_to_interp

!-------------------------------------------------------------------------------

FUNCTION found_fields_to_interp() RESULT(r)
! Return true is any fields to interpolate have been set
IMPLICIT NONE
LOGICAL :: r
CHARACTER(LEN=*), PARAMETER :: routinename = 'FOUND_FIELDS_TO_INTERP'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ALLOCATED(fields_to_interp)) THEN
  r = .TRUE.
ELSE
  r= .FALSE.
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION found_fields_to_interp

!-------------------------------------------------------------------------------

END MODULE find_fields_to_interp_mod

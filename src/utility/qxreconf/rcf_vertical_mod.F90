! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Wrapper for vertical interpolation

MODULE Rcf_vertical_Mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

!  Subroutine Rcf_Vertical - wrapper for vertical interpolation
!
! Description:
! This module contains a wrapper subroutine for vertical
! interpolation. Data is left on the current decomposition
! as there should be no spatial dependencies.
!
! Method:
!   If no interp. required, data is copied. Otherwise, loop
!   through levels calling relevant interpolation routine for each.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_VERTICAL_MOD'

CONTAINS

SUBROUTINE Rcf_vertical( field_in, field_out, grid_in, grid_out, &
                         heights_in, heights_out )

USE umPrintMgr, ONLY: umprint, newline, umMessage

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_V_Int_Ctl_Mod, ONLY: &
    v_int_order

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_all,                      &
    interp_v_only,                   &
    interp_h_only,                   &
    interp_copy,                     &
    interp_no_op

USE nlstcall_mod, ONLY:   &
    LTimer

USE vert_interp_mod, ONLY: vert_interp

USE cppxref_mod, ONLY:             &
    ppx_type_real,                  &
    ppx_type_int,                   &
    ppx_type_log

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (field_type), INTENT(INOUT)  :: field_in
TYPE (field_type), INTENT(INOUT)  :: field_out
TYPE (grid_type), INTENT(IN)      :: grid_in
TYPE (grid_type), INTENT(IN)      :: grid_out
REAL                              :: heights_in(field_in % level_size,&
                                       0 : grid_in % model_levels+1)
REAL                              :: heights_out(field_out %level_size,&
                                       0 : grid_out % model_levels+1)

! Local Data
CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_VERTICAL'
CHARACTER (LEN=errormessagelength)  :: Cmessage
INTEGER                           :: ErrorStatus
INTEGER                           :: i
INTEGER                           :: j
INTEGER                           :: k
INTEGER                           :: start_level_in
INTEGER                           :: start_level_out
INTEGER                           :: end_level_in
INTEGER                           :: end_level_out
INTEGER                           :: stash_set_level_size

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (LTimer) CALL Timer( RoutineName, 3 )

! Find bottom and top level for output field.
IF (field_out % bottom_level >= 0) THEN
  start_level_out = field_out % bottom_level
ELSE
  start_level_out = 1
END IF

IF (field_out % top_level >= 0) THEN
  end_level_out = field_out % top_level
ELSE
  end_level_out = field_out % levels
END IF

! Find bottom and top level for input field.
IF (field_in % bottom_level >= 0) THEN
  start_level_in = field_in % bottom_level
ELSE
  start_level_in = 1
END IF

IF (field_in % top_level >= 0) THEN
  end_level_in = field_in % top_level
ELSE
  end_level_in = field_in % levels
END IF


stash_set_level_size = end_level_in - start_level_in + 1
IF (stash_set_level_size /= field_in % levels) THEN
  WRITE(cmessage,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)')                        &
        'Input field information does not match STASHmaster specification!'  &
        // newline // "STASH section ", field_in % stashmaster % section,    &
        " Item no. ", field_in % stashmaster % item, newline //              &
        "STASHmaster defines ", stash_set_level_size,                        &
        " levels : ", start_level_in, " to ", end_level_in , newline //      &
        "field from input file has ", field_in % levels, " levels"
  ErrorStatus = 5
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

stash_set_level_size = end_level_out - start_level_out + 1
IF (stash_set_level_size /= field_out % levels) THEN
  WRITE(cmessage,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)')                        &
        'Output field information does not match STASHmaster specification!' &
        // newline // "STASH section ", field_out % stashmaster % section,   &
        " Item no. ", field_out % stashmaster % item, newline //             &
        "STASHmaster defines ", stash_set_level_size,                        &
        " levels : ", start_level_out, " to ", end_level_out , newline //    &
        "field for output file has ", field_out % levels, " levels"
  ErrorStatus = 6
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF


! sizes should be the same, but will check
IF ( field_in % level_size /= field_out % level_size ) THEN
  WRITE(cmessage,'(A,I0,A,I0,A,I0,A,I0)') &
        "Size of input and output fields does not match !" // newline //      &
        "STASH section ", field_in % stashmaster % section,                   &
        " Item no. ",field_in % stashmaster % item,  newline //               &
        "Level size in input file ", field_in % level_size,                   &
        " does not equal level size in Output file ", field_out % level_size
  ErrorStatus = 7
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! start levels should be either 0 or 1
IF ( ABS(start_level_in - start_level_out) > 1 ) THEN
  WRITE(cmessage,'(A,I0,A,I0,A,I0,A,I0)') &
        "Input and output fields have incompatible start levels !"             &
        // newline // "They should be either 0 or 1." // newline //            &
        "STASH section ", field_in % stashmaster % section,                    &
        " Item no. ",field_in % stashmaster % item,  newline //                &
        "Start level in input file is : ", start_level_in,                     &
        newline // "Start level in Output file is : ", start_level_out
  ErrorStatus = 8
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!-----------------------------------------------------------------
! Is interpolation activated? If not, copy data across is all we
! will do.
!-----------------------------------------------------------------
SELECT CASE( field_in % interp )
CASE ( interp_copy, interp_h_only )         ! Do a copy for vertical

  ! levels should be a supported transition (e.g. ND->EG)
  SELECT CASE( start_level_in - start_level_out )
  CASE ( 1 ) ! ND -> EG
    IF ( field_in % levels /= field_out % levels - 1 ) THEN
      WRITE(cmessage,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
            "Vertical interpolation is off but no. of vertical levels differs" &
            // newline // "by amount not attributable to ND to EG shift"       &
            // newline // "STASH section ", field_in % stashmaster % section,  &
            " Item no. ",field_in % stashmaster % item,  newline //            &
            "For Input file Start level is : ", start_level_in,                &
            " and level count is : ", field_in % levels, newline //            &
            "For Output file Start level is : ", start_level_out,              &
            " and level count is : ", field_out % levels
      ErrorStatus = 12
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

  CASE ( -1 ) ! EG -> ND
    IF ( field_in % levels /= field_out % levels + 1 ) THEN
      WRITE(cmessage,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
            "Vertical interpolation is off but no. of vertical levels differs" &
            // newline // "by amount not attributable to EG to ND shift"       &
            // newline // "STASH section ", field_in % stashmaster % section,  &
            " Item no. ",field_in % stashmaster % item,  newline //            &
            "For Input file Start level is : ", start_level_in,                &
            " and level count is : ", field_in % levels, newline //            &
            "For Output file Start level is : ", start_level_out,              &
            " and level count is : ", field_out % levels
      ErrorStatus = 13
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

  CASE ( 0 ) ! EG -> EG OR ! ND -> ND
    IF ( field_in % levels /= field_out % levels ) THEN
      WRITE(cmessage,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
            "Vertical interpolation is off but no. of vertical levels differs" &
            // newline // "and both fields are of same ND/EG start level"      &
            // newline // "STASH section ", field_in % stashmaster % section,  &
            " Item no. ",field_in % stashmaster % item,  newline //            &
            "For Input file Start level is : ", start_level_in,                &
            " and level count is : ", field_in % levels, newline //            &
            "For Output file Start level is : ", start_level_out,              &
            " and level count is : ", field_out % levels
      ErrorStatus = 14
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

  CASE DEFAULT
    ! Really don't know how you got here - If difference is >1 it should
    ! Have been trapped earlier in a previous check.
    ! Repeating earlier error message just in case......
    WRITE(cmessage,'(A,I0,A,I0,A,I0,A,I0)') &
          "Input and output fields have incompatible start levels !"           &
          // newline // "They should be either 0 or 1." // newline //          &
          "STASH section ", field_in % stashmaster % section,                  &
          " Item no. ",field_in % stashmaster % item,  newline //              &
          "Start level in input file is : ", start_level_in,                   &
          " start level in Output file is : ", start_level_out
    ErrorStatus = 20
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )

  END SELECT

  SELECT CASE( field_in % stashmaster % data_type )
  CASE ( ppx_type_real )
    CALL copy_input_field_real(field_in, field_out, start_level_in,   &
                               start_level_out, end_level_out  )

  CASE ( ppx_type_int )
    IF ( ALLOCATED( field_in %  data ) ) THEN
      CALL copy_input_field_real(field_in, field_out, start_level_in, &
                                 start_level_out, end_level_out  )

    ELSE
      CALL copy_input_field_int(field_in, field_out, start_level_in,  &
                                start_level_out, end_level_out  )

    END IF

  CASE ( ppx_type_log )
    IF ( ALLOCATED( field_in % data ) ) THEN
      CALL copy_input_field_real(field_in, field_out, start_level_in,  &
                                 start_level_out, end_level_out  )

    ELSE
      CALL copy_input_field_log(field_in, field_out, start_level_in,   &
                                start_level_out, end_level_out  )

    END IF

  CASE DEFAULT
    WRITE(cmessage,'(A,I0,A,I0,A,I0)') &
          "Unsupported Data-Type" // newline //                                &
          "STASH section ", field_in % stashmaster % section,                  &
          " Item no. ",field_in % stashmaster % item,  newline //              &
          "Data-type given by STASHMaster is : ",                              &
           field_in % stashmaster % data_type
    ErrorStatus = 20
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )

  END SELECT

CASE ( interp_all, interp_v_only )      ! Do vertical interpolation

  ! Rcf_generate_heights makes sure we are using the correct height
  ! information for the type of level we are on.

      ! Check input level is valid since it is used as index within
      ! heights_in.
  IF ( start_level_in < 0 ) THEN
    Cmessage = 'Input start level is negative.'
    WRITE(cmessage,'(A,I0,A,I0,A,I0)') &
          "Input start level is negative." // newline //                       &
          "STASH section ", field_in % stashmaster % section,                  &
          " Item no. ",field_in % stashmaster % item,  newline //              &
          "Start level is : ", start_level_in
    ErrorStatus=45
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  ! Loop through levels of output grid
  j = 1
  DO i = start_level_out, end_level_out

    CALL vert_interp( field_in % DATA, field_in % level_size,      &
                      field_in % levels, heights_out(1:,i),        &
                      heights_in(:,start_level_in:), v_int_order,  &
                      field_out % DATA(1:,j) )

    j = j + 1

  END DO

  ! If we are using theta level 0 (i.e. STASH bottom level code is 38) then
  ! lets copy theta level 1 to theta level 0.
  ! For vertical wind this is set to 0
  IF ( field_out % bottom_level == 0 ) THEN
    field_out % DATA(:,1) = field_out % DATA(:,2)
  END IF


CASE ( interp_no_op )
  ! Do nothing

END SELECT

IF (LTimer) CALL Timer( RoutineName, 4 )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

CONTAINS

! Small helper routines to copy the required information.  It uses a different
! array inside the field type depending on datatype.
SUBROUTINE copy_input_field_real(field_in, field_out, start_level_in, &
                                 start_level_out, end_level_out  )

IMPLICIT NONE
TYPE (field_type), INTENT(INOUT)  :: field_in
TYPE (field_type), INTENT(INOUT)  :: field_out
INTEGER, INTENT(IN)               :: start_level_in
INTEGER, INTENT(IN)               :: start_level_out
INTEGER, INTENT(IN)               :: end_level_out

INTEGER :: i
INTEGER :: j
INTEGER :: k

j = 1
DO i = start_level_out, end_level_out
  k = start_level_out - start_level_in + j
  k = MAX(k,1)
  k = MIN(k,field_in % levels)
  field_out % DATA(:,j) = field_in % DATA(:,k)
  j = j + 1
END DO

RETURN
END SUBROUTINE copy_input_field_real

SUBROUTINE copy_input_field_int(field_in, field_out, start_level_in, &
                                 start_level_out, end_level_out  )

IMPLICIT NONE
TYPE (field_type), INTENT(INOUT)  :: field_in
TYPE (field_type), INTENT(INOUT)  :: field_out
INTEGER, INTENT(IN)               :: start_level_in
INTEGER, INTENT(IN)               :: start_level_out
INTEGER, INTENT(IN)               :: end_level_out

INTEGER :: i
INTEGER :: j
INTEGER :: k

j = 1
DO i = start_level_out, end_level_out
  k = start_level_out - start_level_in + j
  k = MAX(k,1)
  k = MIN(k,field_in % levels)
  field_out % Data_int(:,j) = field_in % Data_int(:,k)
  j = j + 1
END DO

RETURN
END SUBROUTINE copy_input_field_int

SUBROUTINE copy_input_field_log(field_in, field_out, start_level_in, &
                                 start_level_out, end_level_out  )

IMPLICIT NONE
TYPE (field_type), INTENT(INOUT)  :: field_in
TYPE (field_type), INTENT(INOUT)  :: field_out
INTEGER, INTENT(IN)               :: start_level_in
INTEGER, INTENT(IN)               :: start_level_out
INTEGER, INTENT(IN)               :: end_level_out

INTEGER :: i
INTEGER :: j
INTEGER :: k

j = 1
DO i = start_level_out, end_level_out
  k = start_level_out - start_level_in + j
  k = MAX(k,1)
  k = MIN(k,field_in % levels)
  field_out % Data_log(:,j) = field_in % Data_log(:,k)
  j = j + 1
END DO

RETURN
END SUBROUTINE copy_input_field_log

END SUBROUTINE Rcf_vertical

END MODULE Rcf_vertical_Mod

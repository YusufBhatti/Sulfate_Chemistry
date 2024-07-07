! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Performs Source=8 field initialisation calculations

MODULE rcf_pre_process_calcs_mod

!  Subroutine Rcf_Pre_Process_Calcs - field initialisation calculations.

! Description:
!   Some fields may have hard-coded Source=8 initialisation. These fields
!   are initialised by the calculations in this routine.

!   NOTE: This routine must run before Post-Processing

! Method:
!   Choice of method applied determined by stashcode.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.4 programming standards.

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_PRE_PROCESS_CALCS_MOD'

CONTAINS

SUBROUTINE rcf_pre_process_calcs( fields_in, fields_out,                       &
  field_count_in, field_count_out,                                             &
  data_source, hdr_in, hdr_out )

USE rcf_field_type_mod, ONLY:                                                  &
  field_type

USE rcf_data_source_mod, ONLY:                                                 &
  data_source_type

USE items_nml_mod, ONLY:                                                       &
  field_calcs

USE submodel_mod, ONLY: atmos_im

USE rcf_umhead_mod, ONLY:                                                      &
  um_header_type

USE decomp_params, ONLY:                                                       &
  decomp_rcf_output

USE rcf_alloc_field_mod, ONLY:                                                 &
  rcf_alloc_field,                                                             &
  rcf_dealloc_field

USE rcf_write_field_mod, ONLY:                                                 &
  rcf_write_field

USE nlstcall_mod, ONLY:                                                        &
    LTimer

USE rcf_derv_cloudfrac_mod, ONLY:                                              &
  rcf_derv_cloudfrac

USE um_stashcode_mod, ONLY:                                                    &
  stashcode_prog_sec,                                                          &
  stashcode_liquid_cf, stashcode_frozen_cf,                                    &
  stashcode_bulk_cf

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER,            INTENT(IN)       :: field_count_in
INTEGER,            INTENT(IN)       :: field_count_out
TYPE( field_type ), POINTER          :: fields_in( : )
TYPE( field_type ), POINTER          :: fields_out( : )
TYPE( data_source_type ), POINTER    :: data_source( : )
TYPE( um_header_type ), INTENT(IN)   :: hdr_in
TYPE( um_header_type ), INTENT(IN)   :: hdr_out

! Local variables
INTEGER                           :: i
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_PRE_PROCESS_CALCS'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ltimer) CALL timer( "Pre_Process_Calcs", 3)

!-----------------------------------------------------------------
! Loop around output fields/sources until a source=field_calcs
! is found
!-----------------------------------------------------------------
DO i = 1, field_count_out
  IF ( data_source( i ) % source == field_calcs ) THEN

    CALL rcf_alloc_field( fields_out( i ) )

    SELECT CASE( fields_out( i ) % stashmaster % model )

    CASE ( atmos_im )

      ! -----------
      ! Atmos items
      ! -----------

      SELECT CASE( fields_out( i ) % stashmaster % section )

        ! For the moment only section zero fields can be USEd here.
      CASE ( stashcode_prog_sec )

        SELECT CASE( fields_out(i) % stashmaster % item )

        CASE ( stashcode_liquid_cf, stashcode_frozen_cf,                       &
          stashcode_bulk_cf )
          CALL rcf_derv_cloudfrac(                                             &
            fields_out( i ) % stashmaster % item,                              &
            fields_out, field_count_out, hdr_out,                              &
            fields_out( i ), data_source )
          CALL rcf_write_field( fields_out( i ), hdr_out,                      &
            decomp_rcf_output )

        END SELECT

      END SELECT

    END SELECT                            ! Selection on Internal Model
    !-----------------------------------------------------------------
    ! Clean up
    !-----------------------------------------------------------------
    CALL rcf_dealloc_field( fields_out( i ) )
  END IF

END DO

IF (ltimer) CALL timer( "Pre_Process_Calcs", 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_pre_process_calcs
END MODULE rcf_pre_process_calcs_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reinitialises area, lquid and frozen cloud fractions.

MODULE rcf_derv_cloudfrac_mod

!  Subroutine Rcf_Derv_CloudFrac

! Description:
!   Derive area, liquid and frozen cloud fractions from bulk cf

! Method:
!   A basic copy with no calculation (area) or basic calculation
!   from qcl and qcf (liquid/frozen).

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.4 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DERV_CLOUDFRAC_MOD'

CONTAINS

SUBROUTINE rcf_derv_cloudfrac( stash_item, fields_out,                         &
                               field_count_out, hdr_out,                       &
                               cloudfrac, data_source )

USE rcf_locate_mod, ONLY:                                                     &
  rcf_locate

USE rcf_alloc_field_mod, ONLY:                                                &
  rcf_alloc_field,                                                             &
  rcf_dealloc_field

USE rcf_read_field_mod, ONLY:                                                 &
  rcf_read_field

USE rcf_data_source_mod, ONLY:                                                &
  data_source_type,                                                            &
  already_processed

USE um_stashcode_mod, ONLY:                                                 &
  stashcode_prog_sec,                                                          &
  stashcode_qcl,       stashcode_qcf,                                          &
  stashcode_bulk_cf,   stashcode_area_cf,                                      &
  stashcode_liquid_cf, stashcode_frozen_cf

USE rcf_field_type_mod, ONLY:                                                 &
  field_type

USE rcf_grid_type_mod, ONLY:                                                  &
  output_grid

USE rcf_level_code_mod, ONLY:                                                 &
  rcf_level_code

USE rcf_umhead_mod, ONLY:                                                     &
  um_header_type

USE umPrintMgr, ONLY:                                                         &
    umPrint,                                                                   &
    umMessage,                                                                 &
    printstatus,                                                               &
    prstatus_normal

USE um_parcore, ONLY:                                                         &
  mype

USE decomp_params, ONLY:                                                      &
  decomp_rcf_output

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( um_header_type), INTENT(IN) :: hdr_out
TYPE( field_type ), INTENT(INOUT), TARGET :: cloudfrac
TYPE( data_source_type ), POINTER :: data_source( : )
INTEGER, INTENT(IN)               :: stash_item
INTEGER, INTENT(IN)               :: field_count_out

! Internal variables
TYPE( field_type ), POINTER       ::  area,qcl,qcf

INTEGER                           ::  pos   ! position in array
INTEGER                           ::  i,j,k ! loop index
INTEGER                           ::  start_level, end_level

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_DERV_CLOUDFRAC'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  IF ( stash_item == stashcode_bulk_cf ) THEN
    WRITE(umMessage,*) 'Reinitialising Bulk CF as Area CF'
    CALL umPrint(umMessage,src='rcf_derv_cloudfrac_mod')
  ELSE IF ( stash_item == stashcode_liquid_cf ) THEN
    WRITE(umMessage,*) 'Reinitialising Liquid CF as fraction of Area CF'
    CALL umPrint(umMessage,src='rcf_derv_cloudfrac_mod')
  ELSE IF ( stash_item == stashcode_frozen_cf ) THEN
    WRITE(umMessage,*) 'Reinitialising Frozen CF as fraction of Area CF'
    CALL umPrint(umMessage,src='rcf_derv_cloudfrac_mod')
  END IF
END IF

!----------------------------------------------------------------------
! Find required fields in output dump and read them in where available
!----------------------------------------------------------------------
! Bulk Cloud fraction and QCL/QCF if required; will abort if not found

CALL rcf_locate(stashcode_prog_sec, stashcode_area_cf,                         &
                fields_out, field_count_out, pos)
area => fields_out(pos)
CALL rcf_alloc_field( area )
CALL rcf_read_field( area, hdr_out, decomp_rcf_output )

IF ( stash_item == stashcode_liquid_cf .OR.                                    &
     stash_item == stashcode_frozen_cf ) THEN

  CALL rcf_locate(stashcode_prog_sec, stashcode_qcl,                           &
                  fields_out, field_count_out, pos)
  qcl => fields_out(pos)
  CALL rcf_alloc_field( qcl )
  CALL rcf_read_field( qcl, hdr_out, decomp_rcf_output )

  CALL rcf_locate(stashcode_prog_sec, stashcode_qcf,                           &
                  fields_out, field_count_out, pos)
  qcf => fields_out(pos)
  CALL rcf_alloc_field( qcf )
  CALL rcf_read_field( qcf, hdr_out, decomp_rcf_output )
END IF

! ---------------------------------------------------------------------
! Check levels for anything ENDGame related!
! ---------------------------------------------------------------------
start_level = cloudfrac % bottom_level
end_level   = cloudfrac % top_level

! ---------------------------------------------------------------------
! Calculate cloud fractions
! ---------------------------------------------------------------------

! Initialise cloud fraction to be area cloud fraction
! If ENDGame grid, i.e. includes surface level 0, then need to make sure
!  that copy from area cf (no surface level) is done properly.
!  In this case level 0 is copy of level 1.
IF ( start_level == 0 ) THEN
  cloudfrac % DATA(:,2:end_level+1) = area % DATA(:,:)
  cloudfrac % DATA(:,1) = cloudfrac % DATA(:,2)
ELSE
  cloudfrac % DATA = area % DATA
END IF

! If liquid cloud fraction then partition as fraction of qcl
IF ( stash_item == stashcode_liquid_cf ) THEN
  DO k = 1, qcl % levels
    DO i = 1, qcl % level_size
      IF (qcl % DATA(i,k) > 0.0) THEN
        cloudfrac % DATA(i,k) = cloudfrac % DATA(i,k) *                        &
          qcl % DATA(i,k) / ( qcl % DATA(i,k) + qcf % DATA(i,k) )
      ELSE
        cloudfrac % DATA(i,k) = 0.0
      END IF
    END DO
  END DO
END IF

! If frozen cloud fraction then partition as fraction of qcf
IF ( stash_item == stashcode_frozen_cf ) THEN
  DO k = 1, qcf % levels
    DO i = 1, qcf % level_size
      IF (qcf % DATA(i,k) > 0.0) THEN
        cloudfrac % DATA(i,k) = cloudfrac % DATA(i,k) *                        &
          qcf % DATA(i,k) / ( qcl % DATA(i,k) + qcf % DATA(i,k) )
      ELSE
        cloudfrac % DATA(i,k) = 0.0
      END IF
    END DO
  END DO
END IF

CALL rcf_locate(stashcode_prog_sec, stash_item,                                &
                fields_out, field_count_out, pos)

data_source( pos ) % source = already_processed

!----------------------------------------------------------------------
! Clear up dynamic memory used
!----------------------------------------------------------------------

CALL rcf_dealloc_field( area )
IF ( stash_item == stashcode_liquid_cf .OR.                                    &
     stash_item == stashcode_frozen_cf ) THEN
  CALL rcf_dealloc_field( qcl )
  CALL rcf_dealloc_field( qcf )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_derv_cloudfrac
END MODULE rcf_derv_cloudfrac_mod

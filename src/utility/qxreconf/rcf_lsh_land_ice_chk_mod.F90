! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_lsh_land_ice_chk_mod

IMPLICIT NONE

!  Subroutine rcf_lsh_land_ice_chk
!
! Description:
!   Ensures that land-ice are consistent with LSH fields.
!
! Method:
! Some LSH fields are sensitive to land-ice values so lets make sure we are
! consistent.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_LSH_LAND_ICE_CHK_MOD'

CONTAINS

SUBROUTINE rcf_lsh_land_ice_chk( fields,      &
                                 field_count, &
                                 decomp, hdr )

USE rcf_field_type_mod, ONLY: &
    field_type

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE rcf_write_field_mod, ONLY: &
    rcf_write_field

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE um_stashcode_mod, ONLY: &
    stashcode_frac_surf_type,  &
    stashcode_zw,              &
    stashcode_fsat,            &
    stashcode_fwetl,           &
    stashcode_prog_sec

USE um_parcore, ONLY:         &
    nproc

! Replaces c_topog.h
USE jules_hydrology_mod, ONLY: zw_max

USE umPrintMgr, ONLY:         &
    umPrint,                   &
    umMessage,                 &
    PrNorm

USE jules_surface_types_mod, ONLY: ice

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER        :: fields(:)
TYPE( um_header_type ), INTENT(IN) :: hdr
INTEGER, INTENT(IN)                :: field_count
INTEGER, INTENT(IN)                :: decomp

! Local variables
INTEGER, PARAMETER                 :: stashlist_size = 3
INTEGER, PARAMETER                 :: stashlist(stashlist_size) = &
                                      (/ stashcode_zw,            &
                                         stashcode_fsat,          &
                                         stashcode_fwetl /)

INTEGER                            :: i, j
INTEGER                            :: reset_count
INTEGER                            :: pos
INTEGER                            :: istat
REAL                               :: value_on_land_ice

TYPE( field_type ), POINTER        :: frac_surf
TYPE( field_type ), POINTER        :: field

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_LSH_LAND_ICE_CHK'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-------------------------------------------------------------------
! Locate required parameters in output dump
!-------------------------------------------------------------------

CALL rcf_locate( stashcode_prog_sec, stashcode_frac_surf_type,         &
                 fields, field_count, pos, zero_ok_arg =.TRUE. )

! If we have surface fractions available then lets check land ice level.
IF (pos /= 0 ) THEN
  frac_surf => fields( pos )
  CALL rcf_alloc_field( frac_surf )
  CALL rcf_read_field( frac_surf, hdr, decomp )

  ! Loop over all stashcodes which we want to check.
  DO j = 1, stashlist_size
    reset_count = 0
    CALL rcf_locate( stashcode_prog_sec, stashlist(j),                 &
                     fields, field_count, pos, zero_ok_arg =.TRUE.)
    IF (pos /= 0) THEN
      field => fields(pos)
      CALL rcf_alloc_field(field)
      CALL rcf_read_field(field, hdr, decomp)

      IF (stashlist(j) == stashcode_zw) THEN
        ! Use value from c_topog header.
        value_on_land_ice = zw_max
      ELSE
        ! Currently other stashcodes checks set to zero.
        value_on_land_ice = 0.0
      END IF

      DO i = 1, field % level_size
        ! Land-ice fraction should be 0 or 1 so lets just check for positive.
        IF ( frac_surf % DATA ( i, ice ) > 0.0 .AND. &
             field % DATA ( i, 1 ) /= value_on_land_ice ) THEN
          ! Set new value in field over land ice.
          field % DATA (i, 1) = value_on_land_ice
          reset_count = reset_count + 1
        END IF
      END DO

      CALL gc_isum (1, nproc, istat, reset_count)
      WRITE(umMessage,'(A,I8,A,I5)')                                        &
                     'LSH land-ice check: No of values reset ',             &
                     reset_count, ' for stashcode ', stashlist(j)
      CALL umPrint(umMessage,src='rcf_lsh_land_ice_chk_mod',level=PrNorm,pe=0)
      !---------------------------------------------------------------------
      ! Write out changed field
      !---------------------------------------------------------------------
      IF (reset_count > 0) THEN
        CALL rcf_write_field( field, hdr, decomp )
      END IF
      CALL rcf_dealloc_field( field )
    END IF
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_lsh_land_ice_chk
END MODULE rcf_lsh_land_ice_chk_mod

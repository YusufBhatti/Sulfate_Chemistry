! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Checks that sea ice fraction has a minimum value.

MODULE Rcf_Sea_Ice_Frac_Chk_Mod

!  Subroutine Rcf_Sea_Ice_Frac_Chk
!
! Description:
!   Ensures that sea ice fraction field is reset to zero
!   when it goes below a threshold.
!
! Method:
!   Reset any values below 0.1 to 0.0
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SEA_ICE_FRAC_CHK_MOD'

CONTAINS

SUBROUTINE Rcf_Sea_Ice_Frac_Chk ( fields, field_count, decomp,        &
                                  hdr, data_source)

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Write_Field_Mod, ONLY: &
    Rcf_Write_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE um_stashcode_mod, ONLY: &
    stashcode_icefrac,         &
    stashcode_icethick,        &
    stashcode_prog_sec

USE Rcf_Data_Source_Mod, ONLY: &
    data_source_type

USE items_nml_mod, ONLY: &
    Field_Calcs

USE UM_ParCore, ONLY: &
    nproc,                  &
    mype

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Normal

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
TYPE( data_source_type ), POINTER  :: data_source(:)

! Local variables
INTEGER                            :: pos_IceFrac
INTEGER                            :: pos_IceThick
INTEGER                            :: IceFrac_changed

REAL, PARAMETER                    :: Threshold = 0.1

INTEGER                            :: i
INTEGER                            :: istat
INTEGER                            :: COUNT
TYPE( field_type ), POINTER        :: Ice_Fraction

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_SEA_ICE_FRAC_CHK'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------
! Loacte Ice Fraction field in output dump
!-------------------------------------------------------------------

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_icefrac, &
                  fields, field_count, pos_IceFrac, zero_ok_arg = .TRUE. )

!--------------------------------------------------------------------
! Reset sea ice fraction to 0 if it is below given threshold value.
!--------------------------------------------------------------------

IceFrac_changed = 0

IF (pos_IceFrac /= 0 ) THEN

  Ice_Fraction => fields( pos_IceFrac )
  CALL Rcf_Alloc_Field( Ice_Fraction )
  CALL Rcf_Read_Field( Ice_Fraction, hdr, decomp )

  COUNT = 0
  DO i = 1, Ice_Fraction % level_size

    IF ( Ice_Fraction % DATA (i,1) <  Threshold .AND.                 &
         Ice_Fraction % DATA (i,1) /= 0.0     ) THEN

      Ice_Fraction % DATA (i,1) = 0.0
      IceFrac_changed = 1
      COUNT = COUNT + 1

    END IF

  END DO

  !---------------------------------------------------------------------
  ! Calculate sum of all changes on all PEs
  !---------------------------------------------------------------------
  CALL gc_isum (1, nproc, istat, COUNT)
  IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
    WRITE(umMessage,*) 'Ice Fraction : No of small values reset to zero ', &
                 COUNT
    CALL umPrint(umMessage,src='rcf_sea_ice_frac_chk_mod')
  END IF

  !---------------------------------------------------------------------
  ! Synchronise `changed' flag
  !---------------------------------------------------------------------
  CALL GC_Imax( 1, nproc, istat, IceFrac_changed )

  !-------------------------------------------------------------------
  ! If there have been changes to IceFrac, find IceThick and reset
  ! datasource to recalculate IceThick from new IceFrac
  !-------------------------------------------------------------------
  IF (IceFrac_changed == 1) THEN
    CALL Rcf_Locate ( stashcode_prog_sec, stashcode_iceThick, &
                      fields, field_count, pos_IceThick, zero_ok_arg = .TRUE.)

    IF (pos_IceThick /= 0 ) THEN
      IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
        WRITE(umMessage,*) 'Ice Fraction : Setting Ice Thickness to be',  &
                    ' recalculated'
        CALL umPrint(umMessage,src='rcf_sea_ice_frac_chk_mod')
      END IF
      data_source( pos_IceThick ) % source      = Field_Calcs
    END IF !(pos_IceThick /= 0 )

    !---------------------------------------------------------------------
    ! Write out changed field
    !---------------------------------------------------------------------

    CALL Rcf_Write_Field( Ice_Fraction, hdr, decomp )
  END IF


  CALL Rcf_Dealloc_Field( Ice_Fraction )
END IF !(pos_IceFrac /= 0 ) - found no icefrac, don't do anything.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Sea_Ice_Frac_Chk
END MODULE Rcf_Sea_Ice_Frac_Chk_Mod

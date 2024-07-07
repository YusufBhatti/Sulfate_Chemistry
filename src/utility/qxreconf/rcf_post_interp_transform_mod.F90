! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Performs transforms after a field has been interpolated.

MODULE Rcf_Post_Interp_Transform_Mod

!  Subroutine Rcf_Post_Interp_Transform
!
! Description:
!   Wrapper to perform tranformations/processing on a field after
!   it has been interpolated.
!
! Method:
!   Choice of transform/processing is based on stashcode.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE UM_ParCore, ONLY: &
    mype,                   &
    nproc

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE

PUBLIC :: Rcf_Post_Interp_Transform, field_set_to_min, field_set_to_max

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_POST_INTERP_TRANSFORM_MOD'

CONTAINS


SUBROUTINE Rcf_Post_Interp_Transform( output_field, fields_out, &
                                      field_count_out )

USE Rcf_V_Int_Ctl_Mod, ONLY: &
    v_int_active

USE Rcf_Exner_P_Convs_Mod, ONLY: &
    Rcf_Conv_P_Exner

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_done

USE um_stashcode_mod, ONLY: &
    stashcode_cca,                 stashcode_cc_lwp,        &
    stashcode_w,                   stashcode_w_adv,         &
    stashcode_exner,               stashcode_q,             &
    stashcode_qcf,                 stashcode_qcl,           &
    stashcode_area_cf,             stashcode_bulk_cf,       &
    stashcode_liquid_cf,           stashcode_frozen_cf,     &
    stashcode_mean_canopyw,        stashcode_can_water_tile,&
    stashcode_qcf2,                stashcode_qrain,         &
    stashcode_qgraup,              stashcode_qc,            &
    stashcode_qT,                  stashcode_prog_sec,      &
    stashcode_proc_phys_sec,       stashcode_tracer_sec,    &
    stashcode_ukca_sec,            stashcode_3d_cca,        &
    stashcode_p

USE rcf_nlist_recon_science_mod, ONLY: &
    q_min,                             &
    w_zero_end,                        &
    w_zero_start

USE nlsizes_namelist_mod, ONLY: &
    model_levels

USE Rcf_Exppx_Mod, ONLY: &
    Rcf_Exppx

USE Submodel_Mod, ONLY: &
    atmos_im

IMPLICIT NONE

! Arguments
TYPE( field_type ), INTENT(INOUT) :: output_field
TYPE( field_type ), POINTER         :: fields_out(:)
INTEGER,            INTENT(IN)      :: field_count_out

! Local variables
INTEGER                             :: level
CHARACTER (LEN=*), PARAMETER        :: RoutineName='RCF_POST_INTERP_TRANSFORM'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------------------------------------
! Only do transforms if interpolation is switched on
!---------------------------------------------------------------
IF ( output_field % interp == interp_done ) THEN

  SELECT CASE( output_field % stashmaster % section )

  CASE ( stashcode_prog_sec )

    ! Which fields do we wish to apply transforms to?
    SELECT CASE( output_field % stashmaster % item )

    CASE ( stashcode_exner )
      ! convert interpolated P back to Exner
      IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
        WRITE(umMessage,'(a25)') 'Converting P to exner'
        CALL umPrint(umMessage,src='rcf_post_interp_transform_mod')
      END IF
      ! This would have been changed back to Exner in pre_interp.
      output_field % stashmaster => Rcf_Exppx( atmos_im, stashcode_prog_sec, &
                                               stashcode_p)

      CALL Rcf_Conv_P_Exner( output_field )

    CASE ( stashcode_w, stashcode_w_adv )
      ! Zero first level (surface)
      output_field % DATA( :,1 ) = 0.0
      IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
        WRITE(umMessage,*) ' Setting w to zero, level 0'
        CALL umPrint(umMessage,src='rcf_post_interp_transform_mod')
      END IF

      DO level = w_zero_start+1, w_zero_end+1
        output_field % DATA( :, level ) = 0.0
        IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
          WRITE(umMessage,*) ' Setting w to zero, level ',level-1
          CALL umPrint(umMessage,src='rcf_post_interp_transform_mod')
        END IF
      END DO

      IF (w_zero_end /= model_levels) THEN
        output_field % DATA( :, output_field % levels ) = 0.0
        IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
          WRITE(umMessage,*) ' Setting w to zero, level ',model_levels
          CALL umPrint(umMessage,src='rcf_post_interp_transform_mod')
        END IF
      END IF

    CASE ( stashcode_cca )
      CALL field_set_to_min( "Convective Cloud Amount", &
                              output_field )

      CALL field_set_to_max( "Convective Cloud Amount", &
                              output_field )

    CASE ( stashcode_3d_cca )
      CALL field_set_to_min( "3D Convective Cloud Amount", &
                              output_field )

      CALL field_set_to_max( "3D Convective Cloud Amount", &
                              output_field )

    CASE ( stashcode_cc_lwp )
      CALL field_set_to_min( "Liquid Water Path", &
                              output_field )

    CASE ( stashcode_q )
      CALL field_set_to_min( "Specific Humidity", &
                              output_field,       &
                              field_min=q_min )

    CASE ( stashcode_qcf )
      CALL field_set_to_min( "QCF", &
                              output_field )

    CASE ( stashcode_qcl )
      CALL field_set_to_min( "QCL", &
                              output_field )

    CASE ( stashcode_qcf2)
      CALL field_set_to_min( "QCF2", &
                              output_field )

    CASE ( stashcode_qrain)
      CALL field_set_to_min( "QRain", &
                              output_field )

    CASE ( stashcode_qgraup)
      CALL field_set_to_min( "Qgraup", &
                              output_field )

    CASE ( stashcode_area_cf )
      CALL field_set_to_min( "Area Cloud Fraction", &
                              output_field )

      CALL field_set_to_max( "Area Cloud Fraction", &
                              output_field )

    CASE ( stashcode_liquid_cf )
      CALL field_set_to_min( "Liquid Cloud Fraction", &
                              output_field )

      CALL field_set_to_max( "Liquid Cloud Fraction", &
                              output_field )

    CASE ( stashcode_bulk_cf )
      CALL field_set_to_min( "Bulk Cloud Fraction", &
                              output_field )

      CALL field_set_to_max( "Bulk Cloud Fraction", &
                              output_field )

    CASE ( stashcode_frozen_cf )
      CALL field_set_to_min( "Frozen Cloud Fraction", &
                              output_field )

      CALL field_set_to_max( "Frozen Cloud Fraction", &
                              output_field )

    CASE ( stashcode_mean_canopyw )
      CALL field_set_to_min( "Canopy Water", &
                              output_field )

    CASE ( stashcode_can_water_tile )
      CALL field_set_to_min( "Canopy Water on Tiles", &
                              output_field )

    END SELECT

  CASE ( stashcode_proc_phys_sec )

    ! Which fields do we wish to apply transforms to?
    SELECT CASE( output_field % stashmaster % item )

    CASE ( MOD(stashcode_qc, 1000) )
      CALL field_set_to_min( "Cloud water content (qc)", &
                              output_field )

    CASE ( MOD(stashcode_qT, 1000) )
      CALL field_set_to_min( "Total specific humidity (qT)", &
                              output_field )

    END SELECT

  END SELECT
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Post_Interp_Transform


SUBROUTINE field_set_to_min( field_name, output_field, &
                             field_min)
! Iterates across whole of output_field, resetting any
! values below the minimum to be at the minimum.
! Default minimum is zero.

IMPLICIT NONE

! Arguments
CHARACTER(LEN=*),     INTENT(IN)      :: field_name
TYPE( field_type ), INTENT(INOUT)   :: output_field

!Optional Argument
REAL, OPTIONAL,     INTENT(IN)      :: field_min

! Local variables
INTEGER                               :: i
INTEGER                               :: j
INTEGER                               :: COUNT
INTEGER                               :: istat

REAL                                  :: MIN

CHARACTER (LEN=*), PARAMETER  :: RoutineName='FIELD_SET_TO_MIN'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Default minimum is zero
IF ( PRESENT( field_min ) ) THEN
  MIN = field_min
ELSE
  MIN = 0.0
END IF

! Check for fields too low and reset to minimum.
DO j = 1, output_field % levels
  COUNT = 0
  DO i = 1, output_field % level_size
    IF ( output_field % DATA(i,j) < MIN ) THEN
      output_field % DATA(i,j) = MIN
      COUNT = COUNT + 1
    END IF
  END DO

  CALL gc_isum (1, nproc, istat, COUNT)
  IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
    WRITE(umMessage,'(a,i4,1x,a,i8)') ' Level ',j, &
         TRIM(field_name)//' reset to minimum ',COUNT
    CALL umPrint(umMessage,src='rcf_post_interp_transform_mod')
  END IF

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE field_set_to_min

SUBROUTINE field_set_to_max( field_name, output_field, &
                             field_max)
! Iterates across whole of output_field, resetting any
! values above the maximum to be at the maximum.
! Default maximum is one.

IMPLICIT NONE

! Arguments
CHARACTER(LEN=*),     INTENT(IN)      :: field_name
TYPE( field_type ), INTENT(INOUT)   :: output_field

!Optional Argument
REAL, OPTIONAL,     INTENT(IN)      :: field_max

! Local variables
INTEGER                               :: i
INTEGER                               :: j
INTEGER                               :: COUNT
INTEGER                               :: istat

REAL                                  :: MAX
CHARACTER (LEN=*), PARAMETER  :: RoutineName='FIELD_SET_TO_MAX'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Default maximum is one
IF ( PRESENT( field_max ) ) THEN
  MAX = field_max
ELSE
  MAX = 1.0
END IF

! Check for fields too high and reset to maximum.
DO j = 1, output_field % levels
  COUNT = 0
  DO i = 1, output_field % level_size
    IF ( output_field % DATA(i,j) > MAX ) THEN
      output_field % DATA(i,j) = MAX
      COUNT = COUNT + 1
    END IF
  END DO

  CALL gc_isum (1, nproc, istat, COUNT)
  IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
    WRITE(umMessage,'(a,i4,1x,a,i8)') ' Level ',j, &
    TRIM(field_name)//' reset to maximum ',COUNT
    CALL umPrint(umMessage,src='rcf_post_interp_transform_mod')
  END IF

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE field_set_to_max

END MODULE Rcf_Post_Interp_Transform_Mod

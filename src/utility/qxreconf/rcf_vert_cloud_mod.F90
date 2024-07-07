! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   vertical interpolation of some cloud variables

MODULE Rcf_Vert_Cloud_Mod
IMPLICIT NONE

!  Subroutine Rcf_Vert_Cloud
!
! Description:
!   Performs vertical interpolation for certain cloud variables
!   (convective cloud top and bottom for example)
!
! Method:
!    The variables store a model level - thus the nearest output level
!    to the input level is calculated as the interpolated level.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_VERT_CLOUD_MOD'

CONTAINS

SUBROUTINE Rcf_Vert_Cloud( field_in, field_out, grid_in, grid_out, &
                           heights_in, heights_out )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_v_only,                   &
    interp_h_only,                   &
    interp_all,                      &
    interp_copy

USE UM_ParCore, ONLY:       &
    nproc,                  &
    mype

USE missing_data_mod, ONLY: rmdi

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), INTENT(IN)    :: field_in
TYPE( field_type ), INTENT(INOUT) :: field_out
TYPE( grid_type ),  INTENT(IN)    :: grid_in
TYPE( grid_type ),  INTENT(IN)    :: grid_out
REAL,               INTENT(IN)    :: heights_in(field_in % level_size,&
                                       0 : grid_in % model_levels+1)
REAL,               INTENT(IN)    :: heights_out(field_out %level_size,&
                                       0 : grid_out % model_levels+1)
! Local variables
INTEGER                           :: ErrorStatus
INTEGER                           :: i
INTEGER                           :: j
INTEGER                           :: fixcount
INTEGER                           :: istat
REAL                              :: level_height_in
CHARACTER (LEN=*), PARAMETER      :: RoutineName = 'RCF_VERT_CLOUD'
CHARACTER (LEN=errormessagelength)     :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Initial checks
!------------------------------------------------------------------
IF ( field_out % levels /= 1 .OR. field_in % levels /= 1 ) THEN
  ErrorStatus = 10
  Cmessage = 'Cloud base/top interpolation requires single level field'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!------------------------------------------------------------------
! If no interpolation required, just copy
!------------------------------------------------------------------
SELECT CASE( field_in % interp )
CASE ( interp_h_only, interp_copy )
  IF (field_in % level_size /= field_out % level_size ) THEN
    ErrorStatus = 20
    Cmessage = 'Unable to copy data as initial '//&
        'and final grid sizes are different'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  IF ( field_in % interp == interp_copy) THEN
    field_out % Data_Int(:,:) = field_in % Data_Int(:,:)
  ELSE
    field_out % DATA(:,:) = field_in % DATA(:,:)
  END IF

CASE ( interp_all, interp_v_only )

  !-------------------------------------------------------------------
  ! Find output height nearest to input height
  !-------------------------------------------------------------------
  fixcount = 0
  DO i = 1, field_out % level_size
    ! Do nothing if level is 0
    IF ( NINT(field_in % DATA(i,1)) <= 0 ) THEN
      field_out % DATA(i,1) = 0.0
    ELSE IF ( NINT(field_in % DATA(i,1)) == NINT(rmdi)) THEN
      field_out % DATA(i,1) = rmdi
    ELSE IF ( NINT(field_in % DATA(i,1)) > grid_in % model_levels+1) THEN
      field_out % DATA(i,1) = rmdi
    ELSE
      ! find the input height
      level_height_in = heights_in( i, NINT(field_in % DATA(i,1)) )

      ! set output height to missing data
      field_out % DATA(i,1) = rmdi

      ! First check that height of input is above 1st level of heights_out
      IF ( level_height_in >= heights_out(i,1) ) THEN
        DO j = 1, grid_out % model_levels - 1
          IF ( level_height_in >= heights_out(i, j) .AND.    &
               level_height_in < heights_out(i, j+1) ) THEN
            field_out % DATA(i,1) = REAL(j)
            EXIT
          END IF
        END DO
      ELSE
        ! Input height is less that first level so lets assume its the ground.
        field_out % DATA(i,1) = 0.0
      END IF

    END IF
    IF (field_out % DATA(i,1) == rmdi) THEN
      fixcount = fixcount + 1
      field_out % DATA(i,1) = REAL( grid_out % model_levels )
    END IF
  END DO

  CALL gc_isum(1,nproc, istat, fixcount)
  IF (mype == 0 .AND. fixcount /= 0) THEN
    ErrorStatus = -14
    WRITE (cmessage,'(A,I8)') 'Number of cloud points fixed in ' //           &
                              'rcf_vert_cloud: ', fixcount
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Vert_Cloud
END MODULE Rcf_Vert_Cloud_Mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initialisations for FLake lake model prognostics.

MODULE Rcf_Init_Flake_Mod

IMPLICIT NONE

! Subroutine Rcf_Init_Flake
!
! Description:
!   Initialises FLake model prognostics
!
! Method:
!   Based on land conditions of surface and soil temperatures,
!   and presence of snow.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INIT_FLAKE_MOD'

CONTAINS

SUBROUTINE Rcf_Init_Flake( fields_out,  field_count_out, hdr_out, &
                           flake_field, field_stashcode )


USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec,        &
    stashcode_tstar_tile,      &
    stashcode_snow_tile,       &
    stashcode_soil_temp,       &
    stashcode_flake_depth,     &
    stashcode_flake_t_mean,    &
    stashcode_flake_t_mxl,     &
    stashcode_flake_t_ice,     &
    stashcode_flake_h_mxl,     &
    stashcode_flake_h_ice

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Normal

USE UM_ParCore, ONLY: &
    mype

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_Lsm_Mod, ONLY: &
    local_lsm_out

USE water_constants_mod, ONLY: &
    tm

USE jules_surface_types_mod, ONLY: lake
USE jules_snow_mod,          ONLY: rho_snow_const

USE lake_mod, ONLY:      &
       h_snow_min_flk    &
      ,h_ice_min_flk     &
      ,h_ice_max         &
      ,lake_h_mxl_0      &
      ,lake_h_deltaT     &
      ,lake_h_scale

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
! Arguments
TYPE( field_type ), POINTER          :: fields_out(:)
TYPE( field_type ), INTENT(INOUT)  :: flake_field
TYPE( um_header_type ), INTENT(IN) :: hdr_out

INTEGER, INTENT(IN)                :: field_count_out
INTEGER, INTENT(IN)                :: field_stashcode

! Local variables
TYPE( field_type ), POINTER          :: tstar_tile
TYPE( field_type ), POINTER          :: snow_tile
TYPE( field_type ), POINTER          :: soil_temp
TYPE( field_type ), POINTER          :: lake_depth

REAL                                 :: lake_t_mean &
                                       (flake_field % level_size)
REAL                                 :: lake_t_mxl  &
                                       (flake_field % level_size)
REAL                                 :: lake_t_ice  &
                                       (flake_field % level_size)
REAL                                 :: lake_h_mxl  &
                                       (flake_field % level_size)
REAL                                 :: lake_h_ice  &
                                       (flake_field % level_size)

REAL                                 :: soil_mean_temp     &
                                       ,lake_ice_thickness &
                                       ,lake_snow_depth

INTEGER                              :: pos  ! field position
INTEGER                              :: i    ! looper
INTEGER                              :: j    ! looper
INTEGER                              :: l    ! looper
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_INIT_FLAKE'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,*) 'Initialising FLake variable ',field_stashcode
  CALL umPrint(umMessage,src='rcf_init_flake_mod')
END IF

!----------------------------------------------------------------------
! Find tile T* temperature in output fields and read it in
!----------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_tstar_tile,         &
                 fields_out, field_count_out, pos)
tstar_tile => fields_out(pos)
CALL Rcf_Alloc_Field( tstar_tile )
CALL Rcf_Read_Field( tstar_tile, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Find tile snow in output fields and read it in
!----------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_snow_tile,         &
                 fields_out, field_count_out, pos)
snow_tile => fields_out(pos)
CALL Rcf_Alloc_Field( snow_tile )
CALL Rcf_Read_Field( snow_tile, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Find soil temperature in output fields and read it in
!----------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_soil_temp,         &
                 fields_out, field_count_out, pos)
soil_temp => fields_out(pos)
CALL Rcf_Alloc_Field( soil_temp )
CALL Rcf_Read_Field( soil_temp, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Find lake depth in output fields and read it in
!----------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_flake_depth,        &
                 fields_out, field_count_out, pos)
lake_depth => fields_out(pos)
CALL Rcf_Alloc_Field( lake_depth )
CALL Rcf_Read_Field( lake_depth, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Perform FLake init calculations
!----------------------------------------------------------------------
DO l = 1, flake_field % level_size

  ! calculate mean soil temperature
  !
  soil_mean_temp = 0.0
  DO i = 1, Output_Grid % sm_levels
    soil_mean_temp = soil_mean_temp + soil_temp % DATA(l,i) * Output_Grid % soil_depths(i)
  END DO
  soil_mean_temp= soil_mean_temp / SUM( Output_Grid % soil_depths )

  ! estimate the snow depth
  lake_snow_depth = MAX(0.0, snow_tile % DATA(l, lake ) / rho_snow_const )

  ! 1) empirical estimate of lake ice thickness
  !    based on the soil column temperature
  !
  lake_ice_thickness = 0.0
  IF (soil_mean_temp                    <= tm) THEN
    lake_ice_thickness =   SUM( Output_Grid % soil_depths ) &
        * MIN(1.0,(tm - soil_mean_temp         + 0.05)/lake_h_deltaT)
  ELSE IF (soil_temp % DATA(l, 1)      <= tm) THEN
    lake_ice_thickness = Output_Grid % soil_depths(1) &
        * MIN(1.0,(tm - soil_temp % DATA(l, 1) + 0.05)/lake_h_deltaT)
  ELSE IF (tstar_tile % DATA(l, lake ) <= tm) THEN
    lake_ice_thickness = 2.0 * h_ice_min_flk
  ELSE
    lake_ice_thickness = 0.0
  END IF

  ! 2) empirical estimate of lake ice thickness
  !    based on the snow depth
  lake_ice_thickness = MAX( lake_ice_thickness, &
                            lake_snow_depth / lake_h_scale )


  ! ...and bound by the Mironov & Ritter (2004) maximum
  !    to avoid SQRT(-ve) in FLake_driver
  lake_ice_thickness = MIN( lake_ice_thickness,h_ice_max )

  !----------------------------------------------------------------------
  ! Decide which calculations are necessary
  !----------------------------------------------------------------------
  SELECT CASE( field_stashcode )

  CASE (stashcode_flake_t_mxl,stashcode_flake_t_mean)

    ! set mixed-layer T based on the temperature of the 1st soil level
    lake_t_mxl(l) = MAX(soil_temp % DATA(l,1), tm + 0.10)

    ! set mean T to the mean temperature over the soil column,
    ! BUT constrain to between the mixed-layer temperature and freezing.
    lake_t_mean(l) = MIN(soil_mean_temp,  &
                            lake_t_mxl(l) - 0.05)
    lake_t_mean(l) = MAX(lake_t_mean(l), tm + 0.05)


  CASE (stashcode_flake_h_ice,stashcode_flake_h_mxl)

    ! initialise the ice thickness
    lake_h_ice(l) = lake_ice_thickness

    ! bound the mixed-layer depth by the available unfrozen depth
    lake_h_mxl(l) = MIN(lake_h_mxl_0,  &
                           lake_depth % DATA(l,1)-lake_h_ice(l))
    lake_h_mxl(l) = MAX(lake_h_mxl(l),0.0)


  CASE (stashcode_flake_t_ice)

    ! linear interpolation between T* and Tm
    IF ( lake_ice_thickness + lake_snow_depth > EPSILON(1.0) ) THEN
      lake_t_ice(l) = tstar_tile % DATA(l, lake) +            &
                      ( tm - tstar_tile % DATA(l, lake) ) *   &
                      lake_snow_depth /                       &
                      ( lake_ice_thickness + lake_snow_depth ) 
      !
      ! take the soil temperature instead if it is lower
      lake_t_ice(l) = MIN( lake_t_ice(l), soil_temp % DATA(l, 1) )
    ELSE
      lake_t_ice(l) = tstar_tile % DATA(l, lake)
    END IF

    ! put an upper bound on the ice T
    lake_t_ice(l) = MIN(lake_t_ice(l),tm - 0.05)

  END SELECT


  !----------------------------------------------------------------------
  ! write the correct data to the output field
  !----------------------------------------------------------------------
  SELECT CASE( field_stashcode )

  CASE (stashcode_flake_t_mxl)

    flake_field % DATA(l,1) = lake_t_mxl(l)

  CASE (stashcode_flake_t_mean)

    flake_field % DATA(l,1) = lake_t_mean(l)

  CASE (stashcode_flake_h_ice)

    flake_field % DATA(l,1) = lake_h_ice(l)

  CASE (stashcode_flake_h_mxl)

    flake_field % DATA(l,1) = lake_h_mxl(l)

  CASE (stashcode_flake_t_ice)

    flake_field % DATA(l,1) = lake_t_ice(l)

  END SELECT

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Init_Flake

END MODULE Rcf_Init_Flake_Mod

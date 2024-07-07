! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sets up idealised orography for the reconfiguration

MODULE Rcf_Ideal_Set_Orography_Mod
IMPLICIT NONE

!  Subroutine Rcf_Ideal_Set_Orography - set an idealised output orography.
!
! Description:
! This module sets up output orography for the idealised model,
! choosing one from a range of pre-defined orography types.
!
! Method:
!  Determine which orography type is to be used and populate the
!  field accordingly. The input orography is not used. This must be
!  done here rather than in field_calcs later on because many other
!  fields depend on orography.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_IDEAL_SET_OROGRAPHY_MOD'

CONTAINS

SUBROUTINE Rcf_Ideal_Set_Orography( orography, grid_out, hdr_out )


USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE UM_ParCore, ONLY:       &
    mype

USE UM_ParVars, ONLY:       &
    g_datastart

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Grid_Type_Mod, ONLY: &
    Grid_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_HeadAddress_Mod, ONLY:         &
    RC_LongSpacing,     RC_LatSpacing, &
    RC_FirstLat,        RC_FirstLong

USE ereport_mod, ONLY:      &
    ereport

USE rcf_nlist_recon_idealised_mod, ONLY: &
    h_o,                                 &
    half_width_x,                        &
    half_width_y,                        &
    lambda_fraction,                     &
    phi_fraction,                        &
    surface_type

USE rcf_ideal_baro_geo_p_mod, ONLY: &
    baro_geo_p

USE model_domain_mod, ONLY: &
    model_type,             &
    mt_global,              &
    mt_lam,                 &
    mt_cyclic_lam,          &
    mt_bi_cyclic_lam

USE rcf_ideal_surface_mod, ONLY: &
    surface_zero,                &
    surface_gauss,               &
    surface_schar_ridge,         &
    surface_baroclinic,          &
    surface_dump,                &
    surface_wave_sin,            &
    surface_wave_gauss,          &
    surface_ancil

USE planet_constants_mod, ONLY: &
    g

USE conversions_mod, ONLY: &
    pi_over_180,           &
    pi

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE(field_type),     POINTER           :: orography
TYPE(grid_type),      INTENT(IN)        :: grid_out
TYPE(um_header_type), INTENT(IN)        :: hdr_out

! Local variables
INTEGER :: i
INTEGER :: j
REAL    :: latitude(orography % level_size)
REAL    :: longitude(orography % level_size)
REAL    :: phi_s   ! Southern wave boundary (rad)
REAL    :: phi_n   ! Northern wave boundary (rad)
REAL    :: dphi    ! North-south wave boundary difference (rad)
INTEGER :: icode = 0

! Parameters for use with the bounded sinusoidal and Gaussian waves.
REAL, PARAMETER    :: lat_s = 25.0  ! Southern wave boundary (deg)
REAL, PARAMETER    :: lat_n = 65.0  ! Northern wave boundary (deg)
INTEGER, PARAMETER :: freq  = 2     ! Wave frequency

CHARACTER(LEN=*),   PARAMETER :: RoutineName = 'RCF_IDEAL_SET_OROGRAPHY'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
  CALL umPrint( 'Creating Idealised Orography (stashcode 33) ', &
      src='rcf_ideal_set_orography_mod')
END IF

IF (model_type /= mt_global         .AND.                               &
    model_type /= mt_lam            .AND.                               &
    model_type /= mt_cyclic_lam     .AND.                               &
    model_type /= mt_bi_cyclic_lam) THEN
  icode = model_type
  CALL ereport( RoutineName, icode, 'model_type is not supported'  )
END IF

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
  WRITE(umMessage,'(A,I0)')' Surface_type = ',surface_type
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Calculate latitudes and longitudes in radians:
DO j = 1, grid_out % loc_p_rows
  DO i = 1, grid_out % loc_p_row_length

    latitude(i + (j-1)*grid_out % loc_p_row_length) = &
                    ( hdr_out % RealC(RC_FirstLat) +  &
                         (g_datastart(2,mype)+j-2) *  &
                     hdr_out % RealC(RC_LatSpacing) ) &
                       * pi_over_180

    longitude(i + (j-1)*grid_out % loc_p_row_length) = &
                    ( hdr_out % RealC(RC_FirstLong) +  &
                         (g_datastart(1,mype)+i-2)  *  &
                     hdr_out % RealC(RC_LongSpacing) ) &
                       * pi_over_180

  END DO
END DO

! Rotated grids are not catered for.

! ------------------------
! Define basic orographies
! ------------------------

SELECT CASE(surface_type)

CASE(surface_zero)
! Flat surface

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('Orography: flat (all points 0)',src=RoutineName)
  END IF

  DO i = 1, orography % level_size
    orography % data(i, 1) = 0.0
  END DO


CASE(surface_gauss)
! Gaussian hill surface

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('Orography: Gaussian hill',src=RoutineName)
    WRITE(umMessage, '(3(A,E12.5))') 'height (m): ', h_o,              &
         ', east-west centre (deg): ', lambda_fraction * 360,          &
         ', north-south centre (deg): ',phi_fraction * 180.0 - 90.0
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage, '(2(A,E12.5))') 'half-width x (rad): ', half_width_x, &
                                   ', half-width y (rad): ', half_width_y
    CALL umPrint(umMessage,src=RoutineName)
  END IF
    
  DO i = 1, orography % level_size
    ! Translate latitudes and longitudes:
    longitude(i) = longitude(i) - 2.0 * pi * lambda_fraction
    latitude(i)  = latitude(i) - pi * (phi_fraction - 0.5)

    orography % data(i, 1) = h_o * EXP(-1.0 *                              &
      ( (longitude(i)/half_width_x)**2 + (latitude(i)/half_width_y)**2 ) )
  END DO



CASE(surface_schar_ridge)
! Schar N-S ridge surface

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('Orography: north-south ridge',src=RoutineName)
    WRITE(umMessage, '(3(A,E12.5))') ' height (m): ', h_o,               &
         ', centre (deg): ',lambda_fraction * 360.0,                     &
         ', half-width (rad): ', half_width_x
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  DO i = 1, orography % level_size
    ! Translate longitudes:
    longitude(i) = longitude(i) - 2.0 * pi * lambda_fraction
    orography % data(i, 1) = h_o * COS(pi * longitude(i) / half_width_x)**2  &
                                 * EXP(-(longitude(i)/half_width_y)**2)
  END DO


CASE (surface_baroclinic)
! Baroclinic wave surface
! Surface orography to give lower boundary surface pressure of
! 1000hPa in baroclinic wave test (QJRMS 132, 2943--2975).

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    WRITE(umMessage, '(A)') 'Orography: Jablonowski & Williamson '     &
      // 'baroclinic wave (QJRMS 132, 2943--2975)'
    CALL umPrint(umMessage,src=RoutineName)
  END IF


  DO i = 1, orography % level_size
    orography % data (i, 1) = baro_geo_p(latitude(i), 1.0)/g
  END DO


CASE(surface_wave_sin, surface_wave_gauss)
! Bounded sinusoidal wave surface
! As used in Gerber/Polvani.
!
! Bounded Gaussian wave surface
! As used in Jucker et al. The positive parts of surface_wave_sin.

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    IF (surface_type == surface_wave_sin) THEN
      CALL umPrint('Orography: bounded sinusoidal wave', src=RoutineName)
    ELSE
      CALL umPrint('Orography: bounded Gaussian wave', src=RoutineName)
    END IF
    WRITE(umMessage, '(A,E12.5,2(A,F6.1),A,I0)')        &
      'height (m): ',  h_o,   ', south (deg): ', lat_s, &
    ', north (deg): ', lat_n, ', frequency: ',   freq
    CALL umPrint(umMessage, src=RoutineName)
  END IF

  phi_s = lat_s * pi_over_180
  phi_n = lat_n * pi_over_180
  dphi  = phi_n - phi_s

  DO i = 1, orography % level_size
    IF (latitude(i) > phi_s .AND. latitude(i) < phi_n) THEN
      orography % data(i, 1) = h_o * COS(freq*longitude(i))              &
                             * SIN(pi * (latitude(i) - phi_s) / dphi)**2
    ELSE
      orography % data(i, 1) = 0.0
    END IF
  END DO

  IF (surface_type == surface_wave_gauss) THEN
    DO i = 1, orography % level_size
      orography % data(i, 1) = MAX(0.0, orography % data(i, 1))
    END DO
  END IF


CASE(surface_dump, surface_ancil)
! Surface taken from input dump or ancillary file
! Should already have been handled in rcf_set_orography, not here.

  icode = 100
  CALL ereport(RoutineName, icode, 'Should not be able to reach here!')


CASE DEFAULT

  icode = 200
  WRITE(umMessage, '(A,I0)') 'Unsupported surface_type: ',surface_type
  CALL ereport(RoutineName, icode, umMessage)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Ideal_Set_Orography
END MODULE Rcf_Ideal_Set_Orography_Mod

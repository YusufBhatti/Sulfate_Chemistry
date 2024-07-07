! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine pws_mtn_turb_mod ---------------------------------------
!
!   Purpose: Calculate PWS mountain wave turbulence diags, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE pws_mtn_turb_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_MTN_TURB_MOD'

CONTAINS

SUBROUTINE pws_mtn_turb( icode, cmessage)

! Description/Method:
! This code is based on the fieldcalc mtnstress.F90 routine.
! Diagnostic Mountain Wave Turbulence Predictor produced on uv-points, B-grid.
! Compute 200mb Wind Speed and Direction for use in Advection routine.
! Interpolate orography and model level Gravity Wave Drag stress horizontally.
! Find max stress at each point between specified heights above orography.
! Fieldcalc original used levels 23-47 of 70 level model.
! Height - planet radius - scaled orog gives best signal between 3676-15248m.
! Max stress advected for various thresholds using AdvectMWS_B routine.
! This diagnostic retains a strongly empirical basis.

USE missing_data_mod,      ONLY: rmdi
USE atm_fields_bounds_mod, ONLY: udims, vdims, pdims, pdims_s, vdims_s, &
                                 udims_l, vdims_l
USE conversions_mod,       ONLY: recip_pi_over_180
USE nlsizes_namelist_mod,  ONLY: row_length, rows, n_rows, &
                                 global_row_length
USE UM_ParVars,            ONLY: offx, offy
USE field_types,           ONLY: fld_type_p, fld_type_v
USE halo_exchange,         ONLY: swap_bounds
USE mpp_conf_mod,          ONLY: swap_field_is_scalar
USE atm_fields_real_mod,   ONLY: orography
USE planet_constants_mod,  ONLY: planet_radius
USE level_heights_mod,     ONLY: r_at_v
USE pws_diags_mod,         ONLY: pws_mtn_wave_turb,    &
                                 pws_gwd_stress_lev_u, &
                                 pws_gwd_stress_lev_v, &
                                 pws_wind_ub_200,      &
                                 pws_wind_vb_200,      &
                                 flag_mtn_wave_turb
USE uc_to_ub_mod,          ONLY: uc_to_ub
USE vc_to_vb_mod,          ONLY: vc_to_vb
USE p_to_v_mod,            ONLY: p_to_v
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

IMPLICIT NONE

! Argument list
INTEGER,  INTENT(INOUT) :: icode
CHARACTER (LEN=*) :: Cmessage

! Local Constants:
INTEGER, PARAMETER :: NumContours = 3
REAL,    PARAMETER :: Contours(NumContours) =  &
                       (/0.007,0.065,0.25/) ! The contour thresholds
REAL,    PARAMETER :: MinWSpeed(2)  =   &
                      (/0.0,30.0/) ! Windspeed thresholds for advection
                                   ! here chosen to be 0.0m/s and 30m/s.
REAL,    PARAMETER :: HtLev23of70 = 3767.0 ! Ht(m) of level 23 for 70Lev model
REAL,    PARAMETER :: HtLev47of70 =15248.0 ! Ht(m) of level 47 for 70Lev model
REAL,    PARAMETER :: ScaleOrogLow = 0.5 ! Scale factor for orog at lower end
REAL,    PARAMETER :: ScaleOrogHi = 1.5 ! Scale factor for orog at higher end

! Local variables

INTEGER :: i, j, k
REAL    :: tempStress
REAL    :: MaxStress(udims%i_start:udims%i_end, &
                     vdims%j_start:vdims%j_end)
REAL    :: uB(udims%i_start:udims%i_end, &
              vdims%j_start:vdims%j_end, pdims%k_len)
REAL    :: vB(udims%i_start:udims%i_end, &
              vdims%j_start:vdims%j_end, pdims%k_len)
REAL    :: wind_200_mag(udims_l%i_start:udims_l%i_end, &
                        vdims_l%j_start:vdims_l%j_end)
REAL    :: wind_200_dir(udims_l%i_start:udims_l%i_end, &
                        vdims_l%j_start:vdims_l%j_end)
REAL    :: pws_mtn_wave_turb_halo(udims_l%i_start:udims_l%i_end, &
                                  vdims_l%j_start:vdims_l%j_end)
REAL    :: Stress_mask(udims_l%i_start:udims_l%i_end, &
                       vdims_l%j_start:vdims_l%j_end)
REAL    :: r_at_uv(udims%i_start:udims%i_end, &
                   vdims%j_start:vdims%j_end) 
REAL    :: orog_halo(pdims_s%i_start:pdims_s%i_end, &
                     pdims_s%j_start:pdims_s%j_end) 
REAL    :: orog_vc(vdims%i_start:vdims%i_end, &
                   vdims%j_start:vdims%j_end)
REAL    :: orog_vc_halo(vdims_s%i_start:vdims_s%i_end, &
                        vdims_s%j_start:vdims_s%j_end)
REAL    :: orog_vb(udims%i_start:udims%i_end, &
                        vdims%j_start:vdims%j_end)

REAL    :: CutOff            ! One of the Contours values

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_MTN_TURB'
! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! -----------------------------------
! Compute 200mb Wind Speed and Direction for use in Advection routine.

DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    IF (pws_wind_ub_200(i,j)==0.0 .AND. pws_wind_vb_200(i,j)==0.0) THEN
      wind_200_mag(i,j) = 0.0
      wind_200_dir(i,j) = 0.0
    ELSE 
      wind_200_mag(i,j) = SQRT(pws_wind_ub_200(i,j)*pws_wind_ub_200(i,j)  &
                             + pws_wind_vb_200(i,j)*pws_wind_vb_200(i,j) )
      wind_200_dir(i,j) = ATAN2(pws_wind_ub_200(i,j),pws_wind_vb_200(i,j)) &
                             *recip_pi_over_180
    END IF
  END DO
END DO

! Determine maximum stress between HtLev23of70=3676m and HtLev47of70=15248m.
! Corresponds to model levels 23-47 in 70 level set as used by fieldcalc.
! Doubt this diag really works for wind levels other than close to 200mb.

DO j = pdims%j_start,pdims%j_end
  DO i = pdims%i_start,pdims%i_end
    orog_halo(i,j)=orography(i,j)
  END DO
END DO

CALL swap_bounds(orog_halo,    &
                 pdims%i_len,     &
                 pdims%j_len,  1,  &
                 pdims_s%halo_i, pdims_s%halo_j, fld_type_p,    &
                 swap_field_is_scalar, do_south_arg=.TRUE.)

! Interpolate orography and model level Gravity Wave Drag stress horizontally.

CALL p_to_v(orog_halo,                                          &
                 pdims_s%i_start,pdims_s%i_end,                 &
                 pdims_s%j_start,pdims_s%j_end,                 &
                 vdims%i_start,vdims%i_end,                     &
                 vdims%j_start,vdims%j_end,                     &
                 1,1,orog_vc)

DO j = vdims%j_start,vdims%j_end
  DO i = vdims%i_start,vdims%i_end
    orog_vc_halo(i,j)=orog_vc(i,j)
  END DO
END DO

CALL swap_bounds(orog_vc_halo,    &
                 vdims%i_len,     &
                 vdims%j_len,  1,  &
                 vdims_s%halo_i, vdims_s%halo_j, fld_type_v,     &
                 swap_field_is_scalar, do_west_arg=.TRUE.)

DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    orog_vb(i,j)=0.5*( orog_vc_halo(i+1,j)+orog_vc_halo(i,j) ) 
    MaxStress(i,j) = 0.0
  END DO
END DO

! Find max stress at each point between specified heights above orography.
! Fieldcalc original used levels 23-47 of 70 level model.

CALL  uC_to_uB(pws_gwd_stress_lev_u(:,:,1:pdims%k_len),             &
               row_length,rows,n_rows,                              &
               pdims%k_len,offx,offy,uB)
CALL  vC_to_vB(pws_gwd_stress_lev_v(:,:,1:pdims%k_len),             &
               rows,row_length,n_rows,                              &
               pdims%k_len,offx,offy,global_row_length,vB,uB)

! r_at_v array dimensioned (vdims_l...) in setcona(_4A) 
! so interpolation to r_at_uv is safe.
DO k = 1, pdims%k_len
  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      r_at_uv(i,j)=0.5*( r_at_v(i-1,j,k)+r_at_v(i,j,k) )
    END DO
  END DO

  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
! Height range 3676 to 15248 metre chosen to approximate 
! levels 23 to 47 in standard 70 level 80km lid model set-up.
! Orography scaling is empirical to give some MaxStress signal 
! like in fieldcalc over Andes and Himalayas which is otherwise missing.
      IF ( (r_at_uv(i,j)-planet_radius-                    &
            ScaleOrogLow*orog_vb(i,j) > HtLev23of70) .AND. &
           (r_at_uv(i,j)-planet_radius-                    &
            ScaleOrogHi*orog_vb(i,j) < HtLev47of70) ) THEN
        tempStress = SQRT( uB(i,j,k)**2 + vB(i,j,k)**2 )
        IF ( tempStress > MaxStress(i,j) ) THEN
          MaxStress(i,j) = tempStress
        END IF
      END IF
    END DO
  END DO
END DO

! Initialise output diagnostic
DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    pws_mtn_wave_turb_halo(i,j) = MaxStress(i,j)
  END DO
END DO

CALL swap_bounds( wind_200_mag,                               &
                  udims%i_len,                                &
                  vdims%j_len,  1,                            &
                  udims_l%halo_i, vdims_l%halo_j, fld_type_v, &
                  swap_field_is_scalar)

CALL swap_bounds( wind_200_dir,                               &
                  udims%i_len,                                &
                  vdims%j_len,  1,                            &
                  udims_l%halo_i, vdims_l%halo_j, fld_type_v, &
                  swap_field_is_scalar)

CALL swap_bounds( pws_mtn_wave_turb_halo,                     &
                  udims%i_len,                                &
                  vdims%j_len,  1,                            &
                  udims_l%halo_i, vdims_l%halo_j, fld_type_v, &
                  swap_field_is_scalar)

! Max stress advected for various thresholds using AdvectMWS_B routine.
! This diagnostic retains a strongly empirical basis.

DO k = 1, NumContours     ! loop over the contour threshold values
  !  Determine the arrays of values above the thresholds.
  CutOff = Contours(k)

  stress_mask(:,:)=rmdi
  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      IF ( MaxStress(i,j) > CutOff ) THEN
        stress_mask(i,j)=1.0  ! mask of important stress pts
      END IF
    END DO
  END DO

  CALL AdvectMWS_B( CutOff,                          &
                    MinWSpeed,                       &
                    MaxStress,                       &
                    wind_200_mag, wind_200_dir,      &
                    pws_mtn_wave_turb_halo,          &
                    icode, stress_mask )

END DO

! Copy result into final location
DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    pws_mtn_wave_turb(i,j) = pws_mtn_wave_turb_halo(i,j)
  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_mtn_turb

SUBROUTINE AdvectMWS_B (CutOff,       &  ! in
                        MinWSpeed,    &  ! in
                        MaxStress,    &  ! in
                        WindSpd,      &  ! in
                        WindDir,      &  ! in
                        MWTPred,      &  ! inout
                        ErrorStatus,stress_mask)    ! inout

! Description:
!       The input arrays WindDir and stress_mask are used to calculate
!       the updated output MWTPred which may be modified 
!       downstream of the stress depending upon the windspeed.
!       This is a simplified version of the original fieldcalc algorithm
!       as we need to be MPP compliant.
!       Large halos are used so that we can guarantee that the required
!       pts to update are on processor.
!
! Method:
!       See the technical paper:
!         Turner (1999): "Development of a mountain wave turbulence
!         prediction scheme for civil aviation", Forecasting Research
!         Technical Report No. 265.

USE trignometric_mod, ONLY: sin_v_latitude, FV_cos_theta_latitude
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE atm_fields_bounds_mod, ONLY: udims, udims_l, vdims, vdims_l
USE conversions_mod, ONLY: recip_pi_over_180

USE halo_exchange,          ONLY: swap_bounds
USE mpp_conf_mod, ONLY: swap_field_is_scalar
USE Field_Types
USE UM_ParVars
USE missing_data_mod, ONLY: rmdi

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Subroutine Arguments:
REAL,    INTENT(IN) :: CutOff      ! Threshold value for this call
REAL,    INTENT(IN) :: MinWSpeed(2)! WSpeed needed for advection
REAL,    INTENT(IN) :: MaxStress(udims%i_start:udims%i_end, &
                   vdims%j_start:vdims%j_end)  ! Maximum stress magnitude
REAL,    INTENT(IN) :: WindSpd(udims_l%i_start:udims_l%i_end,   &
                      vdims_l%j_start:vdims_l%j_end)  ! Wind speed field
REAL,    INTENT(IN) :: WindDir(udims_l%i_start:udims_l%i_end,   &
                      vdims_l%j_start:vdims_l%j_end)  ! Wind direction field

REAL,    INTENT(INOUT) :: stress_mask(udims_l%i_start:udims_l%i_end, &
                                      vdims_l%j_start:vdims_l%j_end)  
REAL,    INTENT(INOUT) :: MWTPred(udims_l%i_start:udims_l%i_end,  &
                                  vdims_l%j_start:vdims_l%j_end)

INTEGER, INTENT(INOUT) :: ErrorStatus         ! Return error status

! Local Variables
INTEGER :: i,j,k             ! DO loop variable
INTEGER :: maski, maskj      ! advected flow location
REAL :: latitude             ! Latitude in degrees
REAL :: cosang               ! COS(latitude*deg2r)
REAL :: alpha1               ! The first comparison angle
REAL :: alpha2               ! The second comparison
REAL :: mod_deltalat         ! The absolute value of lat_interval
REAL :: wdir
REAL :: wspd

CHARACTER (LEN=*), PARAMETER :: RoutineName='ADVECTMWS_B'
! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL swap_bounds( stress_mask,     &
                  udims%i_len,     &
                  vdims%j_len,  1, &
                  udims_l%halo_i,  &
                  vdims_l%halo_j,  &
                  fld_type_v,      &
                  swap_field_is_scalar )

DO k=1,2 ! Here we loop over two wind speed thresholds.
         ! To mimic fieldcalc we call this routine with 0.0 and 30.0 m/s
         ! and 'walk' upstream as required using maski,maskj
         ! fieldcalc walked downstream but we need to go the other way
         ! so that we only ever write to non halo points.
         ! In theory we could use any threshold and walk downstream
         ! numerous times as long as the points are available in the halos


  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      maski=i
      maskj=j

! Determine the comparison angles alpha1 and alpha2
! Angles alpha1 and alpha2 equal 27.0 and 63.0 degrees 
! for a square grid on the equator.
      latitude = ASIN(sin_v_latitude(1,maskj))
      cosang = COS(latitude)

! Determine the angles alpha1 and alpha2
      mod_deltalat = delta_phi
      alpha1 = ATAN( delta_lambda * cosang/(2*mod_deltalat))    &
                                  * recip_pi_over_180
                                                         
      alpha2 = ATAN(2*delta_lambda * cosang/mod_deltalat)       &
                                   * recip_pi_over_180 

      wdir = WindDir(maski,maskj) ! wind direction at advected point.
      IF ( wdir < 0.0 ) THEN  ! wdir is the vector (not met.) direction
        wdir = wdir + 360.0
      END IF

      ! Determine the adjacent point in longitude direction upstream
      IF ( (wdir > alpha1) .AND. &
         (wdir < (180.0 - alpha1)) ) THEN  ! NW, W or SW 
        maski=maski-1

      ELSE IF ( (wdir > (180.0 + alpha1)) .AND. &
              (wdir < (360.0 - alpha1)) ) THEN  ! NE, E or SE
        maski=maski+1

      END IF

      ! Determine the adjacent point in latitude direction
      IF ((wdir > (360.0 - alpha2)) .OR.  &
          (wdir < alpha2)) THEN   !NW, N or NE
        maskj=maskj+1

      ELSE IF ((wdir > (180.0 - alpha2)) .AND. &
               (wdir < (180.0 + alpha2))) THEN   !SW, S or SE
        maskj=maskj-1
          
      END IF

      wspd = WindSpd(maski,maskj)  ! wind speed at upstream point.

      IF ( (wspd > MinWSpeed(k)) .AND.   &
             stress_mask(maski,maskj) == 1.0 ) THEN
 
        IF ( MWTPred(i,j) < CutOff ) THEN
           MWTPred(i,j) = MWTPred(maski,maskj)
        END IF

      END IF

    END DO

  END DO

END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE AdvectMWS_B

END MODULE pws_mtn_turb_mod

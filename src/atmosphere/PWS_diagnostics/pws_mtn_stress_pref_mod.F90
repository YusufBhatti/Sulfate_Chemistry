! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
MODULE pws_mtn_stress_pref_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_MTN_STRESS_PREF_MOD'

CONTAINS

SUBROUTINE pws_mtn_stress_pref(pref,  &
                            uwind_plev, vwind_plev, & 
                            gwduB, gwdvB, local_ptheta_b, mwtpred)

! Description
!This calculates
!mountain wave predictor (held in MWTPred) for a certain pressure level 
!It does this by:
!
!1) the pressure levels corresponding to the model levels which UStress and
!   VStress are on are passed into the subroutine (these are in hPa).
!
!2) the pressure level (in hPa) for which the mountain wave stress is to be
!   calculated is also passed in to the subroutine
!
!3) New fields UStress_PRef and VStress_PRef are set up to hold the
!   interpolated values of U and V gravity wave stress at the pressure
!   level of interest.
!
!4) For each gridpoint, the values of UStress and VStress at PRef are found
!   by linear interpolation.
!
!5) For each gridpoint, the magnitude of UStress and VStress is determined
!   and put into the field Stress.
!
!6) Then advect the stress/MWPred field the same way as the Max stress is
!   advected in the original version (this part hasn't been changed).

USE conversions_mod, ONLY: recip_pi_over_180
USE atm_fields_bounds_mod, ONLY: udims, vdims, pdims
USE missing_data_mod, ONLY: rmdi

USE UM_parcore

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine Arguments
REAL , INTENT(IN) :: pref   ! pressure in Pascals.
REAL , INTENT(IN) :: uwind_plev(udims%i_start:udims%i_end, &
                             vdims%j_start:vdims%j_end)   ! U & V wind comps on
REAL , INTENT(IN) :: vwind_plev(udims%i_start:udims%i_end, &
                             vdims%j_start:vdims%j_end)   ! std level
REAL , INTENT(IN) :: gwduB(udims%i_start:udims%i_end, &
                             vdims%j_start:vdims%j_end,pdims%k_len)
REAL , INTENT(IN) :: gwdvB(udims%i_start:udims%i_end, &
                             vdims%j_start:vdims%j_end,pdims%k_len)
REAL , INTENT(IN) :: local_ptheta_b(udims%i_start:udims%i_end, &   ! pascals
                             vdims%j_start:vdims%j_end,pdims%k_len)

REAL , INTENT(OUT) :: mwtpred(udims%i_start:udims%i_end, &
                              vdims%j_start:vdims%j_end)


! Local variables

INTEGER, PARAMETER :: NumContours = 3
REAL,    PARAMETER :: Contours(NumContours) =  &
                          (/0.007,0.065,0.25/) ! The contour thresholds

REAL,    PARAMETER :: params(3) = (/0.0625,42.6666667,1.833333/) 
                                   ! settings fed in previously to fieldcalc
                                   ! to control this calculation. 
REAL,    PARAMETER :: MinWSpeed(2)  =   &
                      (/0.0,30.0/) ! Windspeed thresholds for advection
                                   ! here chosen to be 0.0m/s and 30m/s.


! Local Variables:
INTEGER :: jlwr, jupr, jmid  ! Variables holding upper, lower & mid levels
INTEGER :: i, j, k           ! DO loop variables
REAL :: m_UStress            ! m in y=mx+c for UStress interpolation
REAL :: c_UStress            ! c in y=mx+c for UStress interpolation
REAL :: m_VStress            ! m in y=mx+c for VStress interpolation
REAL :: c_VStress            ! c in y=mx+c for VStress interpolation
REAL :: CutOff               ! One of the Contours values
REAL :: prob                 ! Holds probability

REAL :: WindFPRef(udims%i_start:udims%i_end, &
                  vdims%j_start:vdims%j_end)    ! wind speed at this press level

REAL :: WindDPRef(udims%i_start:udims%i_end, &
                  vdims%j_start:vdims%j_end)    ! wind direction 
                                                ! at this press level

REAL :: UStress_PRef(udims%i_start:udims%i_end, &
                     vdims%j_start:vdims%j_end) ! UStress at PRef

REAL :: VStress_PRef(udims%i_start:udims%i_end, &
                     vdims%j_start:vdims%j_end) ! VStress at PRef

REAL :: Stress(udims%i_start:udims%i_end, &
               vdims%j_start:vdims%j_end)       ! Magnitude of gravity wave 
                                                ! stress at this pressure level

REAL :: Stress_mask(udims%i_start:udims%i_end, &
               vdims%j_start:vdims%j_end) 

! New variables for vectorization

INTEGER :: Lower(udims%i_start:udims%i_end, &
               vdims%j_start:vdims%j_end)
INTEGER :: Upper(udims%i_start:udims%i_end, &
               vdims%j_start:vdims%j_end)


LOGICAL :: following

INTEGER :: errorstatus

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_MTN_STRESS_PREF'
! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

errorstatus=0

DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    IF ( (uwind_plev(i,j) /= rmdi ) .AND. (vwind_plev(i,j) /= rmdi ) ) THEN
      WindFPRef(i,j) = SQRT((uwind_plev(i,j))**2 + (vwind_plev(i,j))**2 )
    ELSE
      WindFPRef(i,j) = rmdi
    END IF
  END DO
END DO

WHERE ( (uwind_plev == 0.0) .AND. (vwind_plev == 0.0) )
  WindDPRef  = 0.0
ELSEWHERE ( (uwind_plev /= rmdi) .AND.(vwind_plev /= rmdi) )
  WindDPRef = ATAN2(uwind_plev , vwind_plev)*recip_pi_over_180
ELSEWHERE
  WindDPRef = rmdi
END WHERE

!--------------------------------------------------------------------------
! Determine model levels above and below this pressure level
! Then determine UStress & VStress at PRef by linear interpolation of
! UStress & VStress from the level above and below. Do this by first
! checking if Ustress is the same above and below PRef (in which case
! m=0) if so, set UStress equal to this value at PRef; similarly for
! VStress. Otherwise perform linear interpolation - first find the
! m and c coefficents for the equation y=mx+c for the Ustress and
! Vstress ready for interpolation.


DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    IF (PRef >= local_ptheta_b(i,j,1)) THEN !all smaller
      Lower(i,j) = 1
      Upper(i,j) = 2
    ELSE IF (PRef < local_ptheta_b(i,j,pdims%k_len)) THEN  
      !all lower
      Lower(i,j) = pdims%k_len-1
      Upper(i,j) = pdims%k_len
    ELSE
      Lower (i,j) = -1 ! One level to find
    END IF
  END DO
END DO

DO k=1 , pdims%k_len-1
  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end

      IF ( ( local_ptheta_b(i,j,k) > Pref  ) &
         .AND. ( local_ptheta_b(i,j,k+1) <= Pref  ) &
         .AND. (Lower(i,j) == -1) ) THEN
        Lower(i,j) = k
        Upper(i,j) = k+1
      END IF

    END DO
  END DO
END DO


DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end

    jlwr=Lower(i,j)
    jupr=Upper(i,j)

    IF (gwduB(i,j,jlwr) == gwduB(i,j,jupr)) THEN
      UStress_PRef(i,j) = gwduB(i,j,jlwr)
    ELSE

      m_UStress =(gwduB(i,j,jlwr)-gwduB(i,j,jupr))  &
               / (local_ptheta_b(i,j,jlwr)- local_ptheta_b(i,j,jupr) )
      c_UStress = gwduB(i,j,jlwr) - &
                  (m_UStress * local_ptheta_b(i,j,jlwr))
      UStress_PRef(i,j) = (m_UStress * PRef )+ c_UStress

    END IF


    IF (gwdvB(i,j,jlwr) == gwdvB(i,j,jupr)) THEN
      VStress_PRef(i,j) = gwdvB(i,j,jlwr)
    ELSE
      m_VStress =(gwdvB(i,j,jlwr)-gwdvB(i,j,jupr)) &
                 /(local_ptheta_b(i,j,jlwr) - local_ptheta_b(i,j,jupr) )
      c_VStress = gwdvB(i,j,jlwr) - &
                  (m_VStress * local_ptheta_b(i,j,jlwr))
      VStress_PRef(i,j) = (m_VStress * PRef )+ c_VStress
    END IF

  END DO
END DO


! Now have 2 fields that hold the interpolated values for Ustress and
! Vstress for this PRef for all grid points. Determine the magnitude at
! each gridpoint of the gravity wave stress.

DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    IF ( (uwind_plev(i,j) /= rmdi ) .AND. (vwind_plev(i,j) /= rmdi ) ) THEN
      stress(i,j) = SQRT((UStress_PRef(i,j))**2 +  &
                              (VStress_PRef(i,j))**2 )
    ELSE
      stress(i,j) = rmdi
    END IF
  END DO
END DO


! Set the MW predictor to the magnitude of the stress at each grid point.
! Then determine whether the stress or wind is strong enough to advect
! the stress and modify the MW predictor if needbe.

MWTPred(:,:) = Stress(:,:)

DO k = 1, NumContours     ! loop over the contour threshold values

  !  Determine the arrays of values above the thresholds.
  CutOff = Contours(k)

  stress_mask(:,:)=rmdi
  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      IF ( Stress(i,j) > CutOff ) THEN
        stress_mask(i,j)=1.0  ! mask of important stress pts
      END IF
    END DO
  END DO  

  ! Determine the downstream points related to input wind thresholds. 
  CALL AdvectMWS_Ext( CutOff,            &
                      MinWSpeed,                       &
                      Stress,                          &
                      WindFPRef,     WindDPRef,        &
                      MWTPred,                         &
                      ErrorStatus,stress_mask)
  
END DO


DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    IF (MWTPred(i,j) < Params(1)) THEN
      prob = 0.0
    ELSE
      prob = Params(2)*MWTPred(i,j) + Params(3)
    END IF
   ! By default these values are worked out from extrapolating the points
   ! 0.0645=>4.5% and 0.25=>15.5% thus when these are extrapolated
   ! we get an equation y=mx+c or y=42.667x+1.8333
   ! this extends up past 100% which needs to stop thus x cannot be
   ! greater than 2.300764
   ! Currently need to scale by 0.1 since it was too strong otherwise.
   ! (This is a worrying inline comment but we keep what was in fieldcalc)
    MWTPred(i,j) = 0.1*MIN(prob,100.0)
 
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE pws_mtn_stress_pref


SUBROUTINE AdvectMWS_Ext (CutOff,       &  ! in
                          MinWSpeed,    &  ! in
                          MaxStress,    &  ! in
                          WindSpd,      &  ! in
                          WindDir,      &  ! in
                          MWTPred,      &  ! inout
                          ErrorStatus,stress_mask)    ! inout

! Description:
!       The input arrays wdir and stress_mask are used to calculate
!       the updated output mwtpred which may be modified 
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

USE trignometric_mod, ONLY: sin_v_latitude
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE atm_fields_bounds_mod, ONLY: udims, vdims, pdims, udims_l, vdims_l, pdims_s
USE conversions_mod, ONLY: recip_pi_over_180

USE halo_exchange,          ONLY: swap_bounds
USE mpp_conf_mod, ONLY: swap_field_is_scalar
USE Field_Types
USE UM_ParVars
USE missing_data_mod, ONLY: rmdi
USE nlsizes_namelist_mod,  ONLY: global_row_length

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
                                 vdims%j_start:vdims%j_end)  
                                                 ! Maximum stress magnitude
REAL,    INTENT(IN) :: WindSpd(udims%i_start:udims%i_end,   &
                               vdims%j_start:vdims%j_end)    
                                                 ! Wind speed field
REAL,    INTENT(IN) :: WindDir(udims%i_start:udims%i_end,   &
                               vdims%j_start:vdims%j_end)    
                                                 ! Wind direction field

REAL,    INTENT(IN) :: stress_mask(udims%i_start:udims%i_end, &
                                 vdims%j_start:vdims%j_end)  

REAL,    INTENT(INOUT) :: MWTPred(udims%i_start:udims%i_end,  &
                                  vdims%j_start:vdims%j_end)

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

REAL :: MWTPred_halo(udims_l%i_start:udims_l%i_end, &
                     vdims_l%j_start:vdims_l%j_end)

CHARACTER (LEN=*), PARAMETER :: RoutineName='ADVECTMWS_EXT'
! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

MWTPred_halo(:,:)=0.0

DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    MWTPred_halo(i,j)=MWTPred(i,j)
  END DO
END DO

CALL swap_bounds( MWTPred_halo,    &
                  udims%i_len,     &
                  vdims%j_len,  1,  &
          udims_l%halo_i, vdims_l%halo_j, fld_type_v,swap_field_is_scalar)


DO k=1,2 ! Here we loop over two wind speed thresholds. 
         ! To mimic fieldcalc we call thsi routine with 0.0 and 30.0 m/s
         ! and 'walk' downstream as required using maski,maskj
         ! In theory we could use any threshold and walk downstream
         ! numerous times as long as the points are available in the halos

  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      wspd = WindSpd(i,j)   ! wind speed at initial stress points
      maski=i
      maskj=j
      IF ( (wspd > MinWSpeed(k)) .AND. stress_mask(i,j) == 1.0 ) THEN
        wdir = WindDir(maski,maskj) ! wind direction at advected point.
        IF ( wdir < 0.0 ) THEN  ! wdir is the vector (not met.) direction
          wdir = wdir + 360.0
        END IF
        ! Determine the comparison angles alpha1 and alpha2
        ! Angles alpha1 and alpha2 equal 27.0 and 63.0 for a square
        ! grid on the equator.
        latitude = ASIN(sin_v_latitude(1,maskj))
        cosang = COS(latitude)

        ! Determine the angles alpha1 and alpha2
        mod_deltalat = delta_phi
        alpha1 = ATAN( delta_lambda * cosang/(2*mod_deltalat))    &
                                          *recip_pi_over_180
                                                         
        alpha2 = ATAN(2* delta_lambda *cosang/mod_deltalat)       &
                                          *recip_pi_over_180 

        ! Determine the adjacent point in longitude direction
        IF ( (wdir > alpha1) .AND. &
           (wdir < (180.0 - alpha1)) ) THEN  ! NW, W or SW
          maski=maski+1
          IF ( mwtpred_halo(maski,maskj) < CutOff ) THEN
            mwtpred_halo(maski,maskj) = MaxStress(i,j)    
          END IF

        ELSE IF ( (wdir > (180.0 + alpha1)) .AND. &
                (wdir < (360.0 - alpha1)) ) THEN  ! NE, E or SE
          maski=maski-1
          IF ( mwtpred_halo(maski,maskj) < CutOff ) THEN
            mwtpred_halo(maski,maskj) = MaxStress(i,j)
          END IF

        END IF

        ! Determine the adjacent point in latitude direction
        IF ( (wdir > (360.0 - alpha2)) .OR.  &
           (wdir < alpha2) ) THEN        !NW, N or NE
          maskj=maskj-1
          IF ( mwtpred_halo(maski,maskj) < CutOff ) THEN
            mwtpred_halo(maski,maskj) = MaxStress(i,j) 
          END IF

        ELSE IF ( (wdir > (180.0 - alpha2)) .AND. &
                (wdir < (180.0 + alpha2)) ) THEN  !SW, S or SE
          maskj=maskj+1
          IF ( mwtpred_halo(maski,maskj) < CutOff ) THEN
            mwtpred_halo(maski,maskj) = MaxStress(i,j)
          END IF

        END IF
 
      END IF

    END DO

  END DO
  CALL swap_bounds( MWTPred_halo,    &
                    udims%i_len,     &
                    vdims%j_len,  1, &
    udims_l%halo_i, vdims_l%halo_j, fld_type_v, swap_field_is_scalar )
END DO



DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    MWTPred(i,j)=mwtpred_halo(i,j)
  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE AdvectMWS_Ext 

END MODULE pws_mtn_stress_pref_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculate contribution of precipitation to extinction for visibility.

! Description:
!   Process fields of precipitation intensity to give scattering coefft
!   in 1/metres. This is added to an input visibility to give a total.
!   Calculated at a single model level or level within surface layer
!   e.g. screen height (1.5m)

! Documentation:
!   Forecasting Research Scientific Paper NO.4
!   Diagnosis of visibility in the UK Met Office Mesoscale Model
!   and the use of a visibility analysis to constrain initial
!   conditions.  SP Ballard, BJ Wright, BW Golding    1992
!     NIMROD diagnostic:
!   Wright, B. J., 1997: Improvements to the Nimrod Visibility
!     Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!   Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!     for Nimrod. Met. Office FR Tech Rep., No. 222.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Service

! Code description:
!   This code is written to UMDP3 standards.

SUBROUTINE vis_precip                                                   &
           (vis_no_precip,                                              &
                                                      !INPUT
            lca,cca,pct,                                                &
                                                      !INPUT
            beta_ls_rain, beta_ls_snow,                                 &
                                                      !INPUT
            beta_c_rain, beta_c_snow,                                   &
                                                      !INPUT
            p_field,points,k1stpt,                                      &
                                                      !INPUT
            vis_overall,vis_lsp,vis_cp,                                 &
                                                      !OUTPUT
            error)                              !OUTPUT

! Modules
USE visbty_constants_mod, ONLY: lnliminalcontrast
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!---------------------------------------------------------------------
! input variables-----------------------------------------------------
!---------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                  &
 p_field,                                                               &
                                      ! IN NO. points in field.
 points,                                                                &
              ! IN Number of gridpoints being processed.
 k1stpt
              ! IN First gridpoint processed within complete field.

REAL, INTENT(IN) ::                                                     &
 vis_no_precip(p_field),                                                &
                                      ! IN Vis outside precip.
 lca(p_field),                                                          &
                                      ! IN Total Layer Cloud.
 cca(p_field),                                                          &
                                      ! IN Convective Cloud.
 beta_ls_rain(p_field),                                                 &
                                      ! IN Scattering in LS Rain.
 beta_ls_snow(p_field),                                                 &
                                      ! IN Scattering in LS Snow.
 beta_c_rain(p_field),                                                  &
                                      ! IN Scattering in Conv Rain
 beta_c_snow(p_field)             ! IN Scattering in Conv Snow

LOGICAL, INTENT(IN) ::                                                  &
 pct                              ! IN T:Cloud amounts are in %
!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
REAL, INTENT(OUT) ::                                                    &
 vis_overall(p_field),                                                  &
                                     ! OUT Visibility overall
 vis_lsp(p_field),                                                      &
                                     ! OUT Visibility in LS Precip.
 vis_cp(p_field)                 ! OUT Visibility in Conv Precip.
INTEGER, INTENT(OUT) :: error    ! OUT Error code
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Local varables:------------------------------------------------------
!  Define local variables ---------------------------------------------
INTEGER :: i       ! Loop counters: I - horizontal field index;

REAL ::                                                                 &
 beta_no_precip,                                                        &
 p_lsp(p_field),                                                        &
 p_cp(p_field)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='VIS_PRECIP'

!---------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
error=0
IF ((k1stpt+points-1) >  p_field) THEN
  error=1
  GO TO 9999
END IF
IF (pct) THEN
  DO i = k1stpt, k1stpt+points-1
    p_cp(i)=cca(i)/100.0
    p_lsp(i)=(1.0-p_cp(i))*lca(i)/100.0
  END DO
ELSE
  DO i = k1stpt, k1stpt+points-1
    p_cp(i)=cca(i)
    p_lsp(i)=(1.0-p_cp(i))*lca(i)
  END DO
END IF

DO i = k1stpt, k1stpt+points-1

  beta_no_precip=-lnliminalcontrast/vis_no_precip(i)

  IF (p_lsp(i)  >   0.0) THEN
    vis_lsp(i) = -lnliminalcontrast /                                   &
      (beta_no_precip +                                                 &
       beta_ls_rain(i) + beta_ls_snow(i))
  ELSE
    vis_lsp(i)=vis_no_precip(i)
  END IF

  IF (p_cp(i)  >   0.0) THEN
    vis_cp(i) = -lnliminalcontrast /                                    &
      (beta_no_precip +                                                 &
       beta_c_rain(i) + beta_c_snow(i))
  ELSE
    vis_cp(i)=vis_no_precip(i)
  END IF

  ! Ensure no rounding problems lead to vis > vis in clear air
  vis_overall(i) = MIN((1.0-p_cp(i)-p_lsp(i))*vis_no_precip(i) +        &
                   p_lsp(i)*vis_lsp(i) +p_cp(i)*vis_cp(i),              &
                   vis_no_precip(i))

END DO

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE vis_precip

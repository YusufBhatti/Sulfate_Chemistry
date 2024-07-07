! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculates visibility probability.

! Description:
!   The visibility probability is similar to the cloud fraction
!   except it records the fraction of a grid box with RH
!   greater than that required for the critical visibility
!   (e.g 1 km), taking into account precipitation.

!   Suitable for single-column use.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Service

! Code description:
!   This code is written to UMDP3 standards.

! Subroutine Interface:
SUBROUTINE calc_vis_prob(                                               &
  pstar,rhcrit,levels,points,pfield,                                    &
                                                  !INPUT
  t,aerosol,l_murk,q,qcl,qcf,vis,nvis,                                  &
                                                  !INPUT
  lca,cca,pct,                                                          &
                                                  !INPUT
  beta_ls_rain, beta_ls_snow,                                           &
                                                  !INPUT
  beta_c_rain, beta_c_snow,                                             &
                                                  !INPUT
  prob_of_vis,                                                          &
                                                  !OUTPUT
  error                                                                 &
                                                  !OUTPUT
 )

! Modules
USE visbty_constants_mod, ONLY: lnliminalcontrast
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                  &
 levels,                                                                &
                         ! IN No. of levels being processed.
 points,                                                                &
                         ! IN No. of gridpoints being processed.
 pfield,                                                                &
                         ! IN No. of points in global field (at one
!                          !    vertical level).
   nvis                ! IN No. of visibility thresholds
REAL, INTENT(IN) ::                                                     &
 pstar(pfield),                                                         &
                         ! IN Surface pressure (Pa).
 rhcrit(levels),                                                        &
                         ! IN Critical relative humidity.  See the
!                          !    the paragraph incorporating eqs P292.11
!                          !    to P292.14; the values need to be tuned
!                          !    for the given set of levels.
   q(pfield,levels),                                                    &
                           ! IN Specific Humidity
!                          !    (kg per kg air).
   qcl(pfield,levels),                                                  &
                           ! Cloud liquid water content at
!                          !     processed levels (kg per kg air).
   qcf(pfield,levels),                                                  &
                           ! Cloud ice content at processed levels
!                          !    (kg per kg air).
   t(pfield,levels),                                                    &
                           ! IN Temperature (K).
   aerosol(pfield,levels),                                              &
                              ! IN Aerosol mixing ratio(ug/kg)
   vis(nvis),                                                           &
                              ! Visibility thresholds
   lca(pfield),                                                         &
                                       ! IN Total Layer Cloud.
   cca(pfield),                                                         &
                                       ! IN Convective Cloud.
   beta_ls_rain(pfield,levels),                                         &
                                       ! IN Scattering in LS Rain.
   beta_ls_snow(pfield,levels),                                         &
                                       ! IN Scattering in LS Snow.
   beta_c_rain(pfield,levels),                                          &
                                       ! IN Scattering in Conv Rain
   beta_c_snow(pfield,levels)      ! IN Scattering in Conv Snow
LOGICAL, INTENT(IN) ::                                                  &
   l_murk               ! IN : Aerosol present

LOGICAL, INTENT(IN) ::                                                  &
 pct                              ! IN T:Cloud amounts are in %

REAL, INTENT(OUT) ::                                                    &
 prob_of_vis(pfield,levels,nvis) ! OUT Vis prob at processed level
!                          !     (decimal fraction).
INTEGER, INTENT(OUT) :: error    ! OUT 0 if OK; 1 if bad arguments.

!---------------------------------------------------------------------
!  Workspace usage----------------------------------------------------
REAL ::                                                                 &
                         ! "Automatic" arrays on Cray.
 vis_threshold(pfield,levels,nvis),                                     &
 prob_of_vis_ls(pfield,levels,nvis),                                    &
 prob_of_vis_c(pfield,levels,nvis),                                     &
 p_lsp(pfield),                                                         &
 p_cp(pfield),                                                          &
 temp ! Temporary

! Parameters
REAL, PARAMETER :: large_distance = 1.0e6

! Local, including SAVE'd, storage------------------------------------
INTEGER :: k,i,j     ! Loop counters: K - vertical level index.
!                       !                I - horizontal field index.
                        !                J - Vis threshold index.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_VIS_PROB'

!-----------------------------------------------------------------------
!  Check input arguments for potential over-writing problems.
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
error=0
IF (points >  pfield) THEN
  error=1
  GO TO 9999
END IF

IF (pct) THEN
  DO i = 1, points
    p_cp(i) =MAX(cca(i)/100.0,0.0)
    p_lsp(i)=MAX((1.0-p_cp(i))*lca(i)/100.0,0.0)
  END DO
ELSE
  DO i = 1, points
    p_cp(i) =MAX(cca(i),0.0)
    p_lsp(i)=MAX((1.0-p_cp(i))*lca(i),0.0)
  END DO
END IF

DO k = 1, levels
  DO i = 1, points
    DO j = 1, nvis
      vis_threshold(i,k,j)=vis(j)
    END DO
  END DO
END DO

! DEPENDS ON: fog_fr
CALL fog_fr(                                                            &
 pstar,rhcrit,levels,pfield,                                            &
 t,aerosol,l_murk,q,qcl,qcf,vis_threshold,prob_of_vis,nvis              &
 )

DO k = 1, levels
  DO i = 1, points
    DO j = 1, nvis
      IF (p_lsp(i)  >   0.0) THEN
        temp = 1.0/vis(j)+                                              &
         (beta_ls_rain(i,k)+beta_ls_snow(i,k))/lnliminalcontrast
        IF (temp  >   0.0) THEN
          vis_threshold(i,k,j)=MAX(1.0/temp , vis(j))
        ELSE
          vis_threshold(i,k,j)=large_distance
        END IF
      ELSE
        vis_threshold(i,k,j)=vis(j)
      END IF
    END DO
  END DO
END DO

! DEPENDS ON: fog_fr
CALL fog_fr(                                                            &
 pstar,rhcrit,levels,pfield,                                            &
 t,aerosol,l_murk,q,qcl,qcf,vis_threshold,prob_of_vis_ls,nvis           &
 )

DO k = 1, levels
  DO i = 1, points
    DO j = 1, nvis
      IF (p_cp(i)  >   0.0) THEN
        temp = 1.0/vis(j)+                                              &
         (beta_c_rain(i,k)+beta_c_snow(i,k))/lnliminalcontrast
        IF (temp  >   0.0) THEN
          vis_threshold(i,k,j)=MAX(1.0/temp , vis(j))
        ELSE
          vis_threshold(i,k,j)=large_distance
        END IF
      ELSE
        vis_threshold(i,k,j)=vis(j)
      END IF
    END DO
  END DO
END DO

! DEPENDS ON: fog_fr
CALL fog_fr(                                                            &
 pstar,rhcrit,levels,pfield,                                            &
 t,aerosol,l_murk,q,qcl,qcf,vis_threshold,prob_of_vis_c,nvis            &
 )

DO k = 1, levels
  DO i = 1, points
    DO j = 1, nvis
      prob_of_vis(i,k,j)=(1.0-p_cp(i)-p_lsp(i))*                        &
                            prob_of_vis(i,k,j) +                        &
                          p_lsp(i)*                                     &
                            prob_of_vis_ls(i,k,j) +                     &
                          p_cp(i)*                                      &
                            prob_of_vis_c(i,k,j)
    END DO
  END DO
END DO

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_vis_prob

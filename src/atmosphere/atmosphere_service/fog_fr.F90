! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculates fraction of grid-box covered with fog.

! Description:
!   Calculates the fraction of a gridsquare with visibility
!   less than threshold, Vis_thresh, given the total water
!   mixing ratio (qt), temperature (T), pressure (p), and the
!   (Murk) aerosol mass mixing ratio (m), assuming a triangular
!   distribution of states about the median, characterised by
!   a critical relative humdity value, RHcrit.
!   NB:  Throughout, levels are counted from the bottom up,
!   i.e. the lowest level under consideration is level 1, the
!   next lowest level 2, and so on.

!   Suitable for single-column use.

! Documentation:
!    Wright, B. J., 1997: Improvements to the Nimrod Visibility
!       Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!    Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!       for Nimrod. Met. Office FR Tech Rep., No. 222.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Service

! Code description:
!   This code is written to UMDP3 standards.

! Subroutine Interface:
SUBROUTINE fog_fr(                                                      &
 p_layer,rhcrit,levels,pfield,                                          &
 t,aerosol,l_murk,q,qcl,qcf,vis,ff,nvis)

! Modules
!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_wat_new     => qsat_wat,                       &
                    l_new_qsat_atm_Serv !Currently defaults to FALSE

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                  &
 levels,                                                                &
                         ! IN No. of levels being processed.
 pfield,                                                                &
                         ! IN No. of points in field (at one level).
 nvis                ! IN No. of visibility thresholds
REAL, INTENT(IN) ::                                                     &
 p_layer(pfield,levels),                                                &
                            ! IN pressure (Pa) at processed levels.
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
   vis(pfield,levels,nvis)  ! Visibility thresholds
LOGICAL, INTENT(IN) ::                                                  &
   l_murk               ! IN : Aerosol present

REAL, INTENT(OUT) ::                                                    &
 ff(pfield,levels,nvis)   ! OUT Vis prob at processed levels
!                          !     (decimal fraction).

! --------------------------------------------------------------------
!    Workspace usage----------------------------------------------------
REAL ::                                                                 &
                         ! "Automatic" arrays on Cray.
 p(pfield),                                                             &
 qt(pfield),                                                            &
                         ! total of cloud water and vapour
 qs(pfield),                                                            &
                         ! Saturated spec humidity for temp T
 qt_thresh(pfield),                                                     &
                         ! modified qt
 bs
!  Local, including SAVE'd, storage------------------------------------
INTEGER :: k,i,j     ! Loop counters: K - vertical level index.
!                       !                I - horizontal field index.
                        !                J - Vis threshold index.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FOG_FR'

!-----------------------------------------------------------------------
!  Subroutine structure :
!  Loop round levels to be processed.
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO k = 1, levels

  !-----------------------------------------------------------------------
  !  1. Calculate Pressure and initialise temporary arrays
  !-----------------------------------------------------------------------

  DO i = 1, pfield
    p(i)=p_layer(i,k)
    qt(i)=q(i,k)+qcl(i,k)
  END DO ! Loop over points

  !-----------------------------------------------------------------------
  !  2.  Calculate total water threshold corresponding to visibility
  !      Since Qs is needed more than once, pre-calculate and pass it
  !-----------------------------------------------------------------------

  IF (l_new_qsat_atm_Serv) THEN
    CALL qsat_wat_new (qs,t(:,k),p,pfield)
  ELSE
    ! DEPENDS ON: qsat_wat
    CALL qsat_wat (qs,t(1,k),p,pfield)
  END IF

  DO j = 1, nvis

    ! DEPENDS ON: vistoqt
    CALL vistoqt( vis(1,k,j), qs, aerosol(1,k), l_murk,                 &
                pfield, qt_thresh )

    !-----------------------------------------------------------------------
    !  3.  Calculate the width of the distribution in total water space, bs:
    !
    !            bs = ( 1 - RHcrit ) * qs(T)
    !
    !-----------------------------------------------------------------------

    DO i = 1, pfield

      bs = (1.0-rhcrit(k)) * qs(i)

      !=======================================================================
      !  4.  Calculate the fraction of states in a triangular
      !      distribution which exceed the total water threshold.
      !=======================================================================

      !-----------------------------------------------------------------------
      !  4.1 If total water threshold value is less than the total water value
      !      minus the width of the distribution, then all of the states have
      !      a total water value exceeding the threshold, so set the
      !      visibility fraction to 1.0
      !-----------------------------------------------------------------------

      IF ( qt_thresh(i)  <=  qt(i)-bs ) THEN

        ff(i,k,j) = 1.0

        !-----------------------------------------------------------------------
        !  4.2 If total water threshold value is greater than the total water
        !      value minus the width of the distribution, but less than the
        !      total water value then the visibility fraction, VF, is given by:
        !
        !                                                     2
        !                              ( qt       - qt + bs  )
        !             VF = 1.0 - 0.5 * (    thresh           )
        !                              ( ------------------- )
        !                              (          bs         )
        !
        !-----------------------------------------------------------------------

      ELSE IF ( qt_thresh(i)  >   qt(i)-bs .AND.                        &
                qt_thresh(i)  <=  qt(i) ) THEN

        ff(i,k,j) = 1.0 - 0.5 *                                         &
             (( qt_thresh(i) - qt(i) + bs )/ bs)**2

        !-----------------------------------------------------------------------
        !  4.3 If total water threshold value is greater than the total water
        !      value, but less than the total water value plus the width of the
        !      distribution, then the visibility fraction, VF, is given by:
        !
        !                                               2
        !                        ( qt + bs - qt        )
        !             VF = 0.5 * (             thresh  )
        !                        ( ------------------- )
        !                        (          bs         )
        !
        !-----------------------------------------------------------------------

      ELSE IF ( qt_thresh(i)  >   qt(i) .AND.                           &
                qt_thresh(i)  <=  qt(i)+bs    ) THEN

        ff(i,k,j)= 0.5 * (( qt(i) + bs - qt_thresh(i))/bs)**2

        !-----------------------------------------------------------------------
        !  4.4 If total water threshold value is greater than the total water
        !      value plus the width of the distribution, then non of the states
        !      have a total water value exceeding the threshold, so set the
        !      visibility fraction to 0.0
        !-----------------------------------------------------------------------

      ELSE

        ff(i,k,j) = 0.0

      END IF

    END DO ! Loop over PFIELD I

  END DO ! Loop over VIS J

END DO ! Loop over levels

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fog_fr

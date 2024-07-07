! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module containing subroutine to calculate the pressure of
!          the 2.0pvu surface and the pressure of the 380K surface.
!          Routine combines them to calculate the pressure of the
!          tropopause and set L_troposphere to .true. for those
!          gridboxes in the troposphere.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_tropopause

USE ukca_constants
USE conversions_mod, ONLY: recip_pi_over_180
USE umPrintMgr
USE ereport_mod, ONLY: ereport
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook
USE vectlib_mod, ONLY:asin_v
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
SAVE
PRIVATE

REAL, PARAMETER :: tpv = 2.0e-6        ! tropopause PV (pvu)
REAL, PARAMETER :: tpt = 380.0         ! tropopause theta (K)

INTEGER,ALLOCATABLE,PUBLIC  :: tropopause_level(:,:)

!     Pressures of the theta, pv, and combined tropopause

REAL,ALLOCATABLE,PUBLIC  :: p_tropopause(:,:)
REAL,ALLOCATABLE,PUBLIC  :: theta_trop(:,:)
REAL,ALLOCATABLE,PUBLIC  :: pv_trop(:,:)

!     Logical set to true for gridpoints within the troposphere

LOGICAL, ALLOCATABLE, PUBLIC :: L_troposphere(:,:,:)

CHARACTER(LEN=errormessagelength) :: cmessage           ! Error message
INTEGER           :: ierr               !   "   code

PUBLIC ukca_calc_tropopause

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_TROPOPAUSE'

CONTAINS

SUBROUTINE ukca_calc_tropopause(                                 &
  row_length, rows, model_levels,                                &
  sin_theta_latitude, theta, pv, pr_boundaries, pr_levels,       &
  p_tropopause, tropopause_level)

!      Description:
!       Subroutine to calculate p_tropopause. This is a weighted
!       average of 2 tropopause definitions. In the extratropics
!       (>= 28 deg), p_tropopause is the pressure (in Pa) of the
!       2.0 PVU surface. In the tropics (<= 13.0 deg), it is the
!       pressure of the 380K isentropic surface. Between the two
!       latitudes, the weighting function is
!       W = A * sech (lat) + B * lat**2 + C and then
!       p_tropopause = W *(PV tropopause) + (1-W) *(380K tropopause)
!
!       This function is from Hoerling, Monthly Weather Review,
!       Vol 121 162-172. "A Global Analysis of STE during Northern
!       Winter."
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

INTEGER, INTENT(OUT) :: tropopause_level(row_length,rows)
REAL, PARAMETER :: lat_etropics = 28.0    ! define extratropics
REAL, PARAMETER :: lat_tropics  = 13.0    ! define tropics
REAL, PARAMETER :: tpv          = 2.0e-6  ! tropopause PV
REAL, PARAMETER :: tpt          = 380.0   ! tropopause theta
REAL, PARAMETER :: fixed_pres   = 40000.0 ! default tropopause

REAL, INTENT(IN):: sin_theta_latitude(row_length,rows) ! sine(latitude)
REAL, INTENT(IN):: theta(row_length,rows,model_levels) ! theta
REAL, INTENT(IN):: pv(row_length,rows,model_levels)    ! pv on model levs
REAL, INTENT(IN):: pr_boundaries(row_length,rows,0:model_levels)
                          ! pressure at layer boundaries
REAL, INTENT(IN):: pr_levels(row_length,rows,1:model_levels)
                          ! pressure at theta levels

REAL, INTENT(OUT) :: p_tropopause(row_length,rows)

!     Local variables

INTEGER                          :: max_trop_level
INTEGER                          :: i,j,l           ! Loop counters
INTEGER                          :: jll,jlu         ! Level indices

REAL                             :: thalph

REAL :: wt (rows)             ! weighting fn
REAL :: wth(model_levels)             ! theta profile
REAL :: wpl(model_levels)             ! pressure profile
REAL :: wpv(model_levels)             ! PV profile
REAL :: theta_latitude(row_length,rows)  ! lat in radians
REAL :: latitude(row_length,rows)        ! lat in degrees

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALC_TROPOPAUSE'


!     Calculate weighting function from Hoerling, 1993 paper
!     W = a sech (lat) + B (lat)**2 + C

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL asin_v(row_length*rows, sin_theta_latitude, theta_latitude)
latitude       = theta_latitude*recip_pi_over_180 ! latitude in deg

DO i = 1,rows
  IF (ABS(latitude(1,i)) >= lat_etropics) THEN
    wt(i) = 1.0                                   ! extratropcs
  ELSE IF (ABS(latitude(1,i)) <= lat_tropics) THEN
    wt(i) = 0.0                                   ! tropics
  ELSE
    wt(i) = 5560.74*2                                            &
          /(EXP(latitude(1,i))+EXP(-1.0*latitude(1,i)))          &
          + 1.67e-3*(latitude(1,i))**2 - 0.307    ! sub-tropics
  END IF
END DO

!     Calculate theta and pv tropopauses

max_trop_level = model_levels-2
DO j = 1,rows
  DO i = 1,row_length
    DO l=1,model_levels
      wth(l) = theta(i,j,l)         ! theta profile
      wpl(l) = pr_levels(i,j,l)     ! pressure profile
      wpv(l) = pv(i,j,l)            ! PV profile
    END DO

    !         Find theta levels which straddle tpt

    DO l = max_trop_level,2,-1
      IF (wth(l) <= tpt .AND. wth(l+1) >= tpt) THEN
        jll = l
        jlu = l+1
        EXIT
      END IF
    END DO

    !         Calculate pressure of theta tropopause

    thalph = (tpt-wth(jll))/(wth(jlu)-wth(jll))
    theta_trop(i,j) = (1.0-thalph)*LOG(wpl(jll))                 &
                    + thalph*LOG(wpl(jlu))
    IF (wpl(jll) < 0.0 .OR. wpl(jlu) < 0.0 .OR. thalph < 0.0     &
      .OR. thalph > 1.0) THEN
      IF (PrintStatus >= PrStatus_Diag) THEN
        DO l=1,model_levels
          WRITE(umMessage,*) i,j,l,theta(i,j,l),                         &
                     pr_levels(i,j,l), pv(i,j,l)
          CALL umPrint(umMessage,src='ukca_tropopause')
        END DO
        WRITE(umMessage,*) jll,jlu,thalph,theta_trop(i,j)
        CALL umPrint(umMessage,src='ukca_tropopause')
      END IF
    END IF
    theta_trop(i,j) = EXP(theta_trop(i,j))

    !         Find theta levels which straddle tpv

    DO l = max_trop_level,2,-1
      IF (ABS(wpv(l)) <= tpv .AND. ABS(wpv(l+1)) >= tpv) THEN
        jll = l
        jlu = l+1
        EXIT
      END IF
    END DO

    !         Calculate pressure of pv tropopause

    thalph = (tpv-ABS(wpv(jll)))/(ABS(wpv(jlu))-ABS(wpv(jll)))
    pv_trop(i,j) = (1.0-thalph)*LOG(wpl(jll))                    &
                 + thalph*LOG(wpl(jlu))
    IF (wpl(jll) < 0.0 .OR. wpl(jlu) < 0.0 .OR. thalph < 0.0     &
      .OR. thalph > 1.0) THEN
      IF (PrintStatus >= PrStatus_Diag) THEN
        WRITE(umMessage,'(A20)') 'UKCA_TROPOPAUSE:'
        CALL umPrint(umMessage,src='ukca_tropopause')
        DO l=1,model_levels
          WRITE(umMessage,'(3I6,3E12.3)') i,j,l,theta(i,j,l),            &
                     pr_levels(i,j,l), pv(i,j,l)
          CALL umPrint(umMessage,src='ukca_tropopause')
        END DO
      END IF
      cmessage = ' Difficulty diagnosing pv tropopause, '//      &
                 'Reverting to default tropopause pressure'
      ierr = -1
      CALL ereport('UKCA_TROPOPAUSE::UKCA_CALC_TROPOPAUSE',      &
                    ierr,cmessage)

      pv_trop(i,j) = LOG(fixed_pres)
    END IF

    pv_trop(i,j) = EXP(pv_trop(i,j))

  END DO
END DO

!     Calculate combined tropopause based on weighting function

L_troposphere = .FALSE.
DO j = 1,rows
  DO i = 1,row_length
    p_tropopause(i,j) = (wt(j)*pv_trop(i,j))                     &
                      + ((1.0-wt(j))*theta_trop(i,j))

    !         Check for whole gridboxes below tropopause and find
    !         the model level which contains the tropopause

    LOOP_trop: DO l=1,model_levels
      IF (pr_boundaries(i,j,l) >= p_tropopause(i,j)) THEN
        L_troposphere(i,j,l)  = .TRUE.
      ELSE
        tropopause_level(i,j) = l
        EXIT LOOP_trop
      END IF
    END DO LOOP_trop

  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calc_tropopause

END MODULE ukca_tropopause

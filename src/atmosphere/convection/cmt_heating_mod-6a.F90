! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the heating due to the dissipation of kinetic energy in
! the convective momentum transport
!
MODULE cmt_heating_mod

IMPLICIT NONE

!
! Description:
!   Calculates the heating due to the dissipation of kinetic energy in
!   the convective momentum transport
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CMT_HEATING_MOD'

CONTAINS

SUBROUTINE cmt_heating(npnts, nlev,                                     &
                       z_theta, z_rho, exner_layer_centres,             &
                       u, v, dubydt, dvbydt,                            &
                       ! Out
                       dthbydt)

USE planet_constants_mod, ONLY: cp

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine arguments

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

INTEGER,INTENT(IN) :: npnts               ! Number of points
INTEGER,INTENT(IN) :: nlev                ! Number of model layers

REAL,INTENT(IN) :: z_theta(npnts,nlev)    ! height of theta levels (m)
REAL,INTENT(IN) :: z_rho(npnts,nlev)      ! height of rho levels (m)
REAL,INTENT(IN) :: exner_layer_centres(npnts,0:nlev)
                                          ! exner pressure at layer centres
REAL,INTENT(IN) :: u(npnts,nlev)          ! Model U field (m/s)
REAL,INTENT(IN) :: v(npnts,nlev)          ! Model V field (m/s)
REAL,INTENT(IN) :: dubydt(npnts,nlev+1)   ! Increments to U due to CMT (m/s2)
REAL,INTENT(IN) :: dvbydt(npnts,nlev+1)   ! Increments to V due to CMT (m/s2)

!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------

REAL,INTENT(INOUT) :: dthbydt(npnts,nlev) ! Increments to potential temp.
                                          ! due to convection (K/s)

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER :: i, k ! loop counters

REAL :: dzkm14        ! Thickness of level k-1/4 (m)
REAL :: dzkp14        ! Thickness of level k+1/4 (m)
REAL :: dzk           ! Thickness of level k (m)
REAL :: uhat          ! U-wind interpolated to theta level (m/s)
REAL :: vhat          ! V-wind interpolated to theta level (m/s)
REAL :: ududt         ! u*du/dt on theta level (m2/s3)
REAL :: vdvdt         ! v*dv/dt on theta level (m2/s3)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CMT_HEATING'

!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO k = 1, nlev-1
  DO i = 1, npnts

    dzkm14 = z_theta(i,k)   -  z_rho(i,k)
    dzkp14 = z_rho(i,k+1)   -  z_theta(i,k)
    dzk    = z_rho(i,k+1)   -  z_rho(i,k)

    uhat   = ( dzkp14 * u(i,k) + dzkm14 * u(i,k+1) ) / dzk
    vhat   = ( dzkp14 * v(i,k) + dzkm14 * v(i,k+1) ) / dzk

    ududt  = uhat * ( dzkp14 * dubydt(i,k) + dzkm14 * dubydt(i,k+1) ) / dzk
    vdvdt  = vhat * ( dzkp14 * dvbydt(i,k) + dzkm14 * dvbydt(i,k+1) ) / dzk

    !      dT/dt on theta_level(k)
    dthbydt(i,k)  = dthbydt(i,k) - (ududt + vdvdt) /                    &
                   (cp * exner_layer_centres(i,k))
  END DO !Loop over points
END DO !Loop over levels

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE cmt_heating
END MODULE cmt_heating_mod

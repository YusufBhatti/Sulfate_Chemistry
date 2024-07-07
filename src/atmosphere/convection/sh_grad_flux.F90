! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!   New Shallow convection scheme - Calculation of turbulent fluxes
!
MODULE sh_grad_flux_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SH_GRAD_FLUX_MOD'
CONTAINS

SUBROUTINE sh_grad_flux(n_sh, ntra, nlev, trlev, maxlev, l_tracer &
,                      r2rho,r2rho_th,dr_across_th                &
,                      wthetav, wthetal, wqt,wtracer              &
,                      dwthetav_dz, dwthetal_dz, dwqt_dz          &
,                      dwtracer_dz)

!
! Purpose:
!   Calculate gradient of the turbulent fluxes
!   w'theta' and w'q' on in cloud levels
!
!
!   Called by SHALLOW_CONV5A
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
INTEGER, INTENT(IN) ::                                            &
  n_sh                                                            &
            ! No. of shallow convection points
, ntra                                                            &
            ! No. of tracer
, maxlev                                                          &
            ! Maximum number of levels where gradient is non zero
, nlev                                                            &
            ! Maximum number of convective cloud levels
, trlev     ! Maximum number of tracer levels

LOGICAL, INTENT(IN) ::                                            &
  l_tracer  ! true - tracers present


REAL, INTENT(IN) ::                                               &
  r2rho(n_sh,nlev)                                                &
                           ! radius**2 density on rho lev (kg/m)
, r2rho_th(n_sh,nlev)                                             &
                           ! radius**2 density on theta lev (kg/m)
, dr_across_th(n_sh,nlev)  !  thickness on theta levels (m)

!
! fluxes  all held on rho levels
!
REAL, INTENT(IN) ::                                               &
  wthetav(n_sh,nlev)                                              &
                          ! w'thetav'
, wthetal(n_sh,nlev)                                              &
                          ! w'theta'
, wqt(n_sh,nlev)                                                  &
                          ! w'qt'
, wtracer(n_sh,trlev,ntra)  ! w'tracer' (kgm/kg/s)
!
! Arguments with intent INOUT:
!
!     None
!
! Arguments with intent OUT:
!
REAL, INTENT(OUT) ::                                              &
 dwthetav_dz(n_sh,nlev)                                           &
                         ! dwthetav/dz  on theta levels
,dwthetal_dz(n_sh,nlev)                                           &
                         ! dwthetal/dz  on theta levels
,dwqt_dz(n_sh,nlev)                                               &
                         ! dwqt/dz      on theta levels
,dwtracer_dz(n_sh,trlev,ntra)
                         ! dwtracer/dz   on theta levels (kg/kg/s)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------


INTEGER :: i,k,ktra         ! loop counters
INTEGER :: max_trlev

REAL ::                                                           &
 rdz                                                              &
           ! 1/(dz*r2*rho)
, r2_kp1                                                          &
             ! r**2 rho at k+1 levels
, r2_k       ! r**2 rho at k levels

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SH_GRAD_FLUX'



! Model constants:
!


!-----------------------------------------------------------------------
! 1.0 Initialise arrays - set all functions to zero
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO k = 1,nlev
  DO i = 1,n_sh
    dwqt_dz(i,k) = 0.0
    dwthetal_dz(i,k) = 0.0
    dwthetav_dz(i,k) = 0.0
  END DO
END DO

!-----------------------------------------------------------------------
! 2.0 Level loop to calculate gradients of fluxes
!-----------------------------------------------------------------------

DO k = 1,maxlev    ! problem if nlev (no rho theta)
  DO i = 1,n_sh

    rdz  =1.0/(dr_across_th(i,k)*r2rho_th(i,k))


    r2_kp1 = r2rho(i,k+1)
    r2_k   = r2rho(i,k)

    dwqt_dz(i,k)     = (r2_kp1*wqt(i,k+1)-r2_k*wqt(i,k))*rdz

    dwthetal_dz(i,k) = (r2_kp1*wthetal(i,k+1)                    &
                                       -r2_k*wthetal(i,k))*rdz

    dwthetav_dz(i,k)= (r2_kp1*wthetav(i,k+1)                     &
                                       -r2_k*wthetav(i,k))*rdz
  END DO
END DO

!-----------------------------------------------------------------------
! 3.0 Tracers
!-----------------------------------------------------------------------


IF (l_tracer) THEN

  ! Check max levels ?
  ! May be problems if tracers on less levels < maxlev ?
  max_trlev = maxlev
  IF (trlev  <=  maxlev) THEN
    max_trlev = trlev-1
  END IF

  DO ktra = 1,ntra
    DO k = 1,trlev
      DO i = 1,n_sh
        dwtracer_dz(i,k,ktra) = 0.0
      END DO
    END DO

    DO k = 1,max_trlev  ! problem if nlev (no rho theta)
      DO i = 1,n_sh

        rdz  =1.0/(dr_across_th(i,k)*r2rho_th(i,k))
        r2_kp1 = r2rho(i,k+1)
        r2_k   = r2rho(i,k)
        dwtracer_dz(i,k,ktra) = (r2_kp1*wtracer(i,k+1,ktra)      &
                                 -r2_k*wtracer(i,k,ktra))*rdz

      END DO
    END DO
  END DO
END IF
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sh_grad_flux
END MODULE sh_grad_flux_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to calculate CMT increments to U and V.
!
MODULE tcs_cmt_incr_6a


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!
! Description:
! This module calculates the cmt increments to U and V winds for
! the tcs warm rain convection schem
!
! Method:
!   <reference to documentation to go here, once available>
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.2 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TCS_CMT_INCR_6A'

CONTAINS

SUBROUTINE calc_cmt_incr(n_xx, nlevs,ntpar_max                       &
   ,                            r2rho, r2rho_th                      &
   ,                            dr_across_rh                         &
   ,                            uw,vw                                &
                              ! output arguements
   ,                            dubydt,dvbydt)

IMPLICIT NONE
!-------------------------------------------------------------------
! Subroutine Arguments
!-------------------------------------------------------------------
!
! Arguments with intent IN:
!
INTEGER, INTENT(IN) ::                                             &
   n_xx                                                            &
   ! total number of congestus convective points
   , nlevs                                                         &
   ! number of model levels
   , ntpar_max
   ! maximum cloud top level

REAL, INTENT(IN) ::                                                &
   r2rho(n_xx,nlevs)                                               &
                            ! r2 rho on rho levels
   , r2rho_th(n_xx,nlevs)                                          &
                            ! r2 rho on theta levels
   , dr_across_rh(n_xx,nlevs)                                      &
                            ! dr across rho levels
   , uw(n_xx,nlevs)                                                &
                            ! U-Component of stress profile (M2)
   , vw(n_xx,nlevs)
! V-Component of stress profile (M2)
!
! Arguments with intent OUT:
!
REAL, INTENT(OUT) ::                                               &
   dubydt(n_xx,nlevs+1)                                            &
                            ! Tendency in U (MS-2)
   , dvbydt(n_xx,nlevs+1)
                            ! Tendency in V (MS-2)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
!
REAL ::                                                            &
   rhodz                    ! r2 rho * dz

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_CMT_INCR'

!-------------------------
! Loop counters
!-------------------------
INTEGER :: i,k

!
!-----------------------------------------------------------------------
! Input to this routine contains stress profile of uw & vw on levels
! from the surface to the top of the inversion.
!-----------------------------------------------------------------------
!
! Lowest uv level surface stress zero
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

k=1
DO i=1,n_xx
  rhodz  = r2rho(i,k)*dr_across_rh(i,k)
  dubydt(i,k) = -uw(i,k)*r2rho_th(i,k)/rhodz
  dvbydt(i,k) = -vw(i,k)*r2rho_th(i,k)/rhodz
END DO

!
! Mixed layer and Cloud layer increments
!
DO k=2,ntpar_max+1
  DO i=1,n_xx

    rhodz  = r2rho(i,k)*dr_across_rh(i,k)
    dubydt(i,k) =                                                 &
       -(r2rho_th(i,k)*uw(i,k)-r2rho_th(i,k-1)*uw(i,k-1))/rhodz
    dvbydt(i,k) =                                                 &
       -(r2rho_th(i,k)*vw(i,k)-r2rho_th(i,k-1)*vw(i,k-1))/rhodz
  END DO
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE calc_cmt_incr

END MODULE tcs_cmt_incr_6a

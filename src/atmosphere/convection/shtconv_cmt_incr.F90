! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!   Shallow convection calculate increments to u and v
!

MODULE shtconv_cmt_incr_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SHTCONV_CMT_INCR_MOD'
CONTAINS

SUBROUTINE shtconv_cmt_incr(n_sh, nlevs,ntpar_max                 &
,                            r2rho, r2rho_th                      &
,                            dr_across_rh                         &
,                            uw,vw                                &
                                  ! output arguements
,                            dubydt,dvbydt)
!
! Purpose:
!  To calculate increments to U and V due to shallow convection.
!
!   Called by Shallow_conv (5A version)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
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
                 ! total number of shallow convcetive points
, nlevs                                                           &
                 ! number of model levels
, ntpar_max      ! maximum cloud top level

REAL, INTENT(IN) ::                                               &
  r2rho(n_sh,nlevs)                                               &
                           ! r2 rho on rho levels
, r2rho_th(n_sh,nlevs)                                            &
                           ! r2 rho on theta levels
, dr_across_rh(n_sh,nlevs)                                        &
                           ! dr across rho levels
, uw(n_sh,nlevs)                                                  &
                     ! U-Component of stress profile (M2
, vw(n_sh,nlevs)     ! V-Component of stress profile (M2
!
! Arguments with intent OUT:
!
REAL, INTENT(OUT) ::                                              &
  dubydt(n_sh,nlevs+1)                                            &
                            ! Tendency in U (MS-2)
, dvbydt(n_sh,nlevs+1)      ! Tendency in V (MS-2)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

INTEGER :: i,k       ! counters
!
REAL ::                                                           &
 rhodz               ! r2 rho * dz

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SHTCONV_CMT_INCR'


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
DO i=1,n_sh
  rhodz  = r2rho(i,k)*dr_across_rh(i,k)
  dubydt(i,k) = -uw(i,k)*r2rho_th(i,k)/rhodz
  dvbydt(i,k) = -vw(i,k)*r2rho_th(i,k)/rhodz
END DO

!
! Mixed layer and Cloud layer increments
!
DO k=2,ntpar_max+1    !
  DO i=1,n_sh

    rhodz  = r2rho(i,k)*dr_across_rh(i,k)
    dubydt(i,k) =                                                 &
        -(r2rho_th(i,k)*uw(i,k)-r2rho_th(i,k-1)*uw(i,k-1))/rhodz
    dvbydt(i,k) =                                                 &
        -(r2rho_th(i,k)*vw(i,k)-r2rho_th(i,k-1)*vw(i,k-1))/rhodz
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE shtconv_cmt_incr
END MODULE shtconv_cmt_incr_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!   Shallow convection scheme - wql & dqsat/dT in cloud
!

MODULE shallow_wql_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SHALLOW_WQL_MOD'
CONTAINS

SUBROUTINE shallow_wql(n_sh, nlev, max_cldlev                     &
,                            ncld_thlev                           &
,                            zcld, wstar_up, root_mb_o_wsc        &
,                            wql_cb, wql_inv                      &
 ! fields on cloud levels
,                           qsat_cld, q_cld, theta_cld            &
,                           t_cld, fql_func                       &
 ! In/Output fields
,                           wql_cld                               &
 ! Output fields
,                       dqsatdt_cld)



!
! Purpose:
!   Shallow convection scheme - routine to calculate the liquid water
!   budget.
!   Make use of the mass flux and qlup calculated earlier.
!
!   Called by SHALLOW_CONV
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
USE planet_constants_mod, ONLY: g, repsilon, rv
USE water_constants_mod, ONLY: lc
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
, nlev                                                            &
                   ! No. of model layers
, max_cldlev       ! maximum noumber of in cloud levels

INTEGER, INTENT(IN) ::                                            &
  ncld_thlev(n_sh) ! number of theta cloud levels

REAL, INTENT(IN)    ::                                            &
  zcld(n_sh)                                                      &
                         ! Depth of cloud layer (m)
, wstar_up(n_sh)                                                  &
                         ! w* convective velocity scale (wsc)
, wql_cb(n_sh)                                                    &
                         ! flux of wql across cloud base
, wql_inv(n_sh)                                                   &
                         ! flux of wql across inversion
, root_mb_o_wsc(n_sh)     ! sqrt of (mb/wsc)


REAL, INTENT(IN) ::                                               &
  theta_cld(n_sh,nlev)                                            &
                        ! theta
, q_cld(n_sh,nlev)                                                &
                        ! q - mixing ratio
, qsat_cld(n_sh,nlev)                                             &
                        ! qsat in cloud
, t_cld(n_sh,nlev)                                                &
                        ! temperature in cloud
, fql_func(n_sh,nlev)   ! fql function


!
! Arguments with intent INOUT:
!
REAL, INTENT(INOUT) ::                                            &
  wql_cld(n_sh,nlev)     ! wql in cloud (includes cloud base)


!
! Arguments with intent OUT:
!
REAL, INTENT(OUT) ::                                              &
  dqsatdt_cld(n_sh,max_cldlev+1)    ! dqsat/dt

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

REAL ::                                                           &
  thetav_cld(n_sh,max_cldlev+1)                                   &
, temp(n_sh)

!
! Loop counters
!

INTEGER :: i,k,itop

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SHALLOW_WQL'

!-----------------------------------------------------------------------
! 1.0  Calculate dqsat/dT  in cloud & store ratio w*/zcld
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!
!   dqsatdt on cloud levels  (expression for mixing ratio)
!
DO k = 1,max_cldlev+1
  DO i = 1,n_sh
    dqsatdt_cld(i,k) = qsat_cld(i,k)*lc                           &
                     /(Rv*t_cld(i,k)*t_cld(i,k))
    thetav_cld(i,k)=theta_cld(i,k)*(1.0+q_cld(i,k)/repsilon)      &
                 /(1.0+q_cld(i,k))
  END DO
END DO
!-----------------------------------------------------------------------
! 2.0 Calculate  in cloud values of wql on uv levels using
!      similarity expression
!-----------------------------------------------------------------------
! values for new water flux parametrisation

k=1
DO i=1,n_sh
  wql_cld(i,k) = wql_cb(i)
  temp(i) =root_mb_o_wsc(i)*(wstar_up(i)**3)/(zcld(i)*g)
END DO

! taking thetav/theta=1.

DO k=2,max_cldlev+1
  DO i=1,n_sh

    wql_cld(i,k) = 0.5*temp(i)*fql_func(i,k)*theta_cld(i,k)     &
                                      /thetav_cld(i,k)

  END DO
END DO

! value at inversion

DO i=1,n_sh
  itop=ncld_thlev(i)+1
  wql_cld(i,itop) = wql_inv(i)
END DO

!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE shallow_wql
END MODULE shallow_wql_mod

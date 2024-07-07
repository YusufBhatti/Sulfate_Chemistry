! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    FUNCTION BDRYV----------------------------------------------------
!
!    Purpose: Special routine to add psuedo source terms to boundary
!             data in limited area.
!    PARAMETERS ARE SPECIFIC TO UK MESOSCALE MODEL
!    Method:  The boundary concentrations are computed using a
!             simple model of transport from sources outside the
!             model. Analysis of the source distribution outside
!             the UK MES shows that it can be well represented by
!             a line source at constant radius from the centre of
!             the model, with a source distribution given by the
!             sum of two Gaussians. Concentrations from these are
!             computed assuming transport using the local windspeed u
!             or 1 m/s, whichever is stronger, over a distance
!             determined from the centroid of the source distribution, x
!             with a linear transformation rate k from emission to
!             aerosol, dry deposition at a rate determined from the
!             dry deposition velocity vd and mean mixed layer depth h.
!             Thus the max concentration is given by
!                 Q/(uh)*k/(k+vd/h)*(1-exp(-k*x/u))
!             The source term is assumed to decrease with level
!             pressure. See forthcoming documentation for details.
!
!    Programming standard: Unified Model Documentation Paper No 3,
!                          Version 7, dated 11/3/93.
!
!
!    Arguments:---------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Tracer Advection
REAL FUNCTION bdryv(                                              &
 wdir                                                             &
,wspeed                                                           &
,press                                                            &
)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

REAL ::                                                           &
 wdir                                                             &
               ! IN Wind direction : Cartesian degrees
,wspeed                                                           &
               ! IN Wind speed m/s
,press         ! IN Pressure

!
!  Local, including SAVE'd, storage------------------------------------
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
REAL ::                                                           &
 zangle1,width1                                                   &
                ! Centre and width of first source Gaussian.
,zangle2,width2                                                   &
                ! Centre and width of second source Gaussian.
,cmax,czer                                                        &
           ! Max concentration and 'background'.
,wdirn,reciproot2pi                                               &
,mixd,travel                                                      &
                ! Average mixed layer depth and travel distance.
,QMAX1                                                            &
                ! Peak height of first source Gaussian.
,qmax2                                                            &
                ! Peak height of second source Gaussian.
,vd                                                               &
                ! Dry deposition velocity.
,k,k1                                                             &
      ! Transformation parameters.
,ph                                                               &
    ! Pressure height scale.
,krat,kt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BDRYV'

PARAMETER(zangle1=178.0)
PARAMETER(width1=5.0)
PARAMETER(zangle2=173.0)
PARAMETER(width2=25.0)
PARAMETER(reciproot2pi=0.3989422803)
PARAMETER(QMAX1=4.7e5,qmax2=9.0e4,mixd=800.0,travel=7.7e5)
PARAMETER(czer=6.0)
PARAMETER(vd=5.0e-3)
PARAMETER(k=3.0e-6, k1=k+vd/mixd, krat=k/k1,kt=-k1*travel)
PARAMETER(ph=3.0e4)
!
!     Max concentration = Q/(uh)*k/(k+vd/h)*(1-exp(-k*x/u))
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
cmax=1.0/MAX(wspeed,1.0)/mixd
cmax=cmax*krat*(1-EXP(kt/MAX(wspeed,1.0)))
wdirn = wdir - zangle1
IF (wdirn  <   -180.0) wdirn=wdirn+360.0
IF (wdirn  >    180.0) wdirn=wdirn-360.0
wdirn=wdirn/width1
bdryv= QMAX1 * EXP(-wdirn*wdirn/2.0)
wdirn = wdir - zangle2
IF (wdirn  <   -180.0) wdirn=wdirn+360.0
IF (wdirn  >    180.0) wdirn=wdirn-360.0
wdirn=wdirn/width2
bdryv= bdryv + qmax2 * EXP(-wdirn*wdirn/2.0)
!
!     Add 'background' value.
!
bdryv= bdryv * cmax + czer
!
!     Reduce concentration with pressure altitude.
!
bdryv= bdryv * EXP(-(1.0e5-press)/ph)
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION bdryv

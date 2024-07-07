! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!   New Shallow convection scheme - calculation of d(w'u')/z
!
MODULE shtconv_grad_stress_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SHTCONV_GRAD_STRESS_MOD'
CONTAINS

SUBROUTINE shtconv_grad_stress(n_sh, nlev,max_cldlev              &
,                       nclev                                     &
,                       timestep, mb, wstar_up                    &
,                       zcld                                      &
,                       mf_cld, w_up, k_func                      &
,                       u,v                                       &
,                      r2rho,r2rho_th,dr_across_rh,dr_across_th   &
                               ! output arguements
,                       uw_cld,vw_cld)
!
! Purpose:
!   Calculate the gradient component of the stress due to shallow
!   convection in the cloud
!
!   Called by shallow_conv.
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
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE tridiag_all_mod, ONLY: tridiag_all
IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!

INTEGER, INTENT(IN) ::                                            &
  nlev                                                            &
           ! No. of model layers
, n_sh                                                            &
           ! No. of shallow convection points
, nclev(n_sh)                                                     &
               ! number of cloud levels
, max_cldlev   ! Maximum number of cloud levels

REAL, INTENT(IN) ::                                               &
  timestep                                                        &
                           ! model timestep (s)
, zcld(n_sh)                                                      &
                           ! cloud layer depth (m) for CMT
, mb(n_sh)                                                        &
                           ! cloud base mass flux
, wstar_up(n_sh)                                                  &
                           ! w* convective velocity scale (m/s)
, w_up(n_sh,nlev)                                                 &
                           ! ensemble vertical velocity (m/s)
, mf_cld(n_sh,nlev)                                               &
                           ! mass flux in (th lev) ensemble (m/s)
, k_func(n_sh,nlev)                                               &
                           ! k function in cld
, u(n_sh,nlev)                                                    &
                           ! U-component of mean wind (MS-1)
, v(n_sh,nlev)             ! V-component of mean wind (MS-1)
REAL, INTENT(IN) ::                                               &
  r2rho(n_sh,nlev)                                                &
, r2rho_th(n_sh,nlev)                                             &
, dr_across_th(n_sh,nlev)                                         &
, dr_across_rh(n_sh,nlev)

!
! Arguments with intent inout:
!
!
! Arguments with intent out:
!
REAL, INTENT(OUT) ::                                              &
  uw_cld(n_sh,nlev)                                               &
                         ! U-component of stress
, vw_cld(n_sh,nlev)      ! V-component of stress

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
!
INTEGER :: i,k    ! counters

REAL ::                                                           &
  visc(n_sh,max_cldlev+1)                                         &
                           ! viscosity profile (M2S-1)
, a(n_sh,max_cldlev+1)                                            &
                           ! tridiagonal matrix elements
, b(n_sh,max_cldlev+1)                                            &
, c(n_sh,max_cldlev+1)                                            &
, ue_tp1(n_sh,max_cldlev+1)                                       &
                             ! after timestep velocity vectors
, ve_tp1(n_sh,max_cldlev+1)                                       &
                             ! for subsequent use
, dz,dz12
INTEGER ::                                                        &
  nclevp1(n_sh)                                                   &
                     ! cloud level for uv cal
, max_cldlev1        ! max cloud levels plus 1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SHTCONV_GRAD_STRESS'


!
!-----------------------------------------------------------------------
! Note 4A code assumed a 2 level inversion now go to 1 level
! ie drop top condition.
!-----------------------------------------------------------------------
!
!  wup and mass_flux input on required levels ?
!
! Calculate the eddy viscosity profile K(z/zcld)
!

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
max_cldlev1 = max_cldlev+1

DO k=1,max_cldlev1  ! from level above cloud base to just
  DO i=1,n_sh

    IF (k <  nclev(i)) THEN         ! in cloud

      visc(i,k)=k_func(i,k)*mf_cld(i,k)                           &
                            *W_up(i,k)*zcld(i)/wstar_up(i)
    ELSE IF (k == nclev(i)) THEN     ! inversion ?

      dz = dr_across_th(i,k)
      visc(i,k)=0.09*mb(i)*dz

    ELSE
      visc(i,k)=0.0
    END IF
  END DO
END DO

!
! Calculate GRADIENT component of stress
!
!
! Use implicit timestepping
!

k=1
DO i=1,n_sh
  dz12 = dr_across_th(i,k)
  dz   = dr_across_rh(i,k)*r2rho(i,k)
  a(i,k)=0.0
  c(i,k)=-visc(i,k)*r2rho_th(i,k)                            &
                          *timestep/(dz*dz12)
  b(i,k)=1.0-a(i,k)-c(i,k)

  nclevp1(i) = nclev(i) + 1
END DO
DO k=2,max_cldlev1
  DO i=1,n_sh
    dz = dr_across_rh(i,k)*r2rho(i,k)
    IF (dz == 0.0) THEN
      WRITE(umMessage,*) ' dz = 0',dr_across_rh(i,k),r2rho(i,k),i,k
      CALL umPrint(umMessage,src='shtconv_grad_stress')
    END IF
    IF (k <= (nclev(i))) THEN
      dz12 = dr_across_th(i,k-1)
      IF (dz12 == 0.0) THEN
        WRITE(umMessage,*) ' dz12 = 0',dr_across_th(i,k-1),i,k
        CALL umPrint(umMessage,src='shtconv_grad_stress')
      END IF
      a(i,k)=-visc(i,k-1)*r2rho_th(i,k-1)                         &
                               *timestep/(dz*dz12)
      dz12 = dr_across_th(i,k)
      c(i,k)=-visc(i,k)*r2rho_th(i,k)                            &
                              *timestep/(dz*dz12)
    ELSE IF (k == (nclev(i)+1)) THEN
      dz12 = dr_across_th(i,k-1)
      a(i,k)=-visc(i,k-1)*r2rho_th(i,k-1)                        &
                              *timestep/(dz*dz12)
      c(i,k)=0.0
    ELSE
      ! elements not required in calculation (zero)
      c(i,k) = 0.0
      a(i,k) = 0.0

    END IF
    b(i,k)=1.0-a(i,k)-c(i,k)

  END DO
END DO

!
! Calculate NEW timestep wind conponents using tridiag
!
CALL TRIDIAG_all(max_cldlev1,n_sh,nclevp1,a,b,c,u,ue_tp1)
CALL TRIDIAG_all(max_cldlev1,n_sh,nclevp1,a,b,c,v,ve_tp1)
!
! Calculate stress profile -Kdu/dz from latest u/v values
!

! Initialise uw,vw
DO k=1,max_cldlev1
  DO i = 1,n_sh
    uw_cld(i,k) = 0.0
    vw_cld(i,k) = 0.0
  END DO
END DO

DO k=1,max_cldlev1
  DO i=1,n_sh
    IF (k <= nclev(i)) THEN
      dz = dr_across_th(i,k)
      uw_cld(i,k)=-visc(i,k)*(ue_tp1(i,k+1)-ue_tp1(i,k))/dz
      vw_cld(i,k)=-visc(i,k)*(ve_tp1(i,k+1)-ve_tp1(i,k))/dz
    END IF
  END DO
END DO

!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE shtconv_grad_stress
END MODULE shtconv_grad_stress_mod

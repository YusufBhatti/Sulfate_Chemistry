! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!   Shallow convection calculate Cloud base stress
!
MODULE shtconv_base_stress_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SHTCONV_BASE_STRESS_MOD'
CONTAINS

SUBROUTINE shtconv_base_stress(n_sh, nlevs, ntml, ntpar,ntpar_max &
,                              timestep                           &
,                              mb,wsc_o_mb                        &
,                              uw0,vw0,du_cb,dv_cb                &
,                              rho_theta,zrho,ztheta              &
,                              flg_uw_shall,flg_vw_shall          &
                                     ! IN/OUT ARGUMENTS
,                              uw,vw                              &
                                    ! OUTPUT ARGUMENTS
,                              uw_shall,vw_shall)

!
! Purpose:
!  To calculate Cloud base stress for shallow cumulus.
!  Also complete calcuation ofstress profile.
!
!   Called by Shallow_conv (5A version)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards


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
              ! Total number of shallow points
, nlevs                                                           &
              ! Number of model levels
, ntml(n_sh)                                                      &
              ! levels of LCL
, ntpar(n_sh)                                                     &
               ! levels of TOP ofcloud layer
,ntpar_max     ! max value of ntpar +1
!
LOGICAL, INTENT(IN) ::                                            &
 flg_uw_shall                                                     &
                  ! STASH FLAGS FOR SHALLOW
,flg_vw_shall     ! CONVECTION STRESS DIAGNOSTIC
!
REAL, INTENT(IN) ::                                               &
  timestep                                                        &
                   ! MODEL timestep (S)
,  mb(n_sh)                                                       &
                   ! Cloud base mass flux (MS-1)
,  wsc_o_mb(n_sh)                                                 &
                   ! Cloud-layer velocity scale (MS-1)
,  uw0(n_sh)                                                      &
                   ! U-component of surface stress (M2S-2)
,  vw0(n_sh)                                                      &
                   ! V-component of surface stress (M2S-2)
,  du_cb(n_sh)                                                    &
                   ! dU across cloud base (m/s)
,  dv_cb(n_sh)                                                    &
                   ! dV across cloud base (m/s)
,  rho_theta(n_sh,nlevs)                                          &
                         ! Density model th levels (kgm-3)
,  zrho(n_sh,nlevs)                                               &
                         ! height of model rho levels (m)
,  ztheta(n_sh,nlevs)    ! height of model theta levels (m)

!
! Arguments with intent INOUT:
!
REAL, INTENT(INOUT) ::                                            &
  uw(n_sh,nlevs)                                                  &
                   ! U-component of STRESS PROFILE (M2S-2)
, vw(n_sh,nlevs)   ! V-component of STRESS PROFILE (M2S-2)
!
! Arguments with intent OUT:
!
REAL, INTENT(OUT) ::                                              &
  uw_shall(n_sh,nlevs)                                            &
                        ! STASH DIAGNOSTIC FOR U-COMP STREss
, vw_shall(n_sh,nlevs)  ! STASH DIAGNOSTIC FOR V-COMP STREss
!
!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

INTEGER :: i,k            ! COUNTERS
!
REAL ::                                                           &
  omg2_jump(n_sh)                                                 &
                          ! jump in Y component of vorticity
, omg1_jump(n_sh)                                                 &
                          ! jump in X-component of vorticity
, zlcl_cmt(n_sh)                                                  &
                          ! lcl for CMT    (m)
, fcmt,dz                                                         &
, a,b,c                                                           &
                          ! coefficients
, t,dz1                                                           &
, z_depth(n_sh)                                                   &
, expadt         ! exp (Adt)

!
REAL ::                                                           &
  beta  = 0.04                                                    &
                   ! Coefficient for CMT calculation
, delta = 2.3                                                     &
                   ! Coefficient for CMT calculation
, GAMMA = 1.63     ! Coefficient for CMT calculation

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SHTCONV_BASE_STRESS'



!-----------------------------------------------------------------------
!
! Calculate jumps in vorticity across cloud base. (This is done by
! assuming that during the time step du and dv vary as exp(-T/TAU).
! Needs to be done to avoid instability around cloud base.)
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO i=1,n_sh

  zlcl_cmt(i) = ztheta(i,ntml(i))   ! cloud base fo CMT
  !
  ! dz (cloud base - theta level below cloud base (zlcl_cmt) )
  !
  dz = zrho(i,ntml(i)+1)-zlcl_cmt(i)
  !
  ! dz across cloud base
  !
  dz1 = ztheta(i,ntml(i)+1) - zlcl_cmt(i)
  !
  ! depth of cloud as defined for uv calculations (different from th &q)
  ! calculations ?
  !
  z_depth(i) = ztheta(i,ntpar(i)) - zlcl_cmt(i)

  !  f(z/zcld) = exp(-fcmt)
  fcmt=beta*wsc_o_mb(i)*(dz1/z_depth(i))
  b=(1.0/zlcl_cmt(i)-(EXP(-fcmt)-1.0)/dz1)
  ! alpha in doc (but extra /dz factor)
  a=zlcl_cmt(i)*mb(i)*b/(delta*dz)
  expadt = EXP(-a*timestep)

  ! beta in documentation

  c = ( b*(1.0-GAMMA/delta) - 1.0/zlcl_cmt(i) )*uw0(i)

  IF (c == 0.0 .AND. du_cb(i) == 0.0) THEN
    omg2_jump=0.0
  ELSE
    t=-LOG( (c*(1.0-expadt)/a+du_cb(i)*(expadt-1.0))/             &
                     ((c-a*du_cb(i))*timestep))/a
    omg2_jump(i)=( c*(1.0-EXP(-a*t))/a + du_cb(i)*EXP(-a*t) )/dz
  END IF

  ! beta in documentation
  c = ( b*(1.0-GAMMA/delta) - 1.0/zlcl_cmt(i) )*vw0(i)

  IF (c == 0.0 .AND. dv_cb(i) == 0.0) THEN
    omg1_jump=0.0
  ELSE
    t=-LOG((c*(1.0-expadt)/a+dv_cb(i)*(expadt-1.0))/              &
                     ((c-a*dv_cb(i))*timestep))/a
    omg1_jump(i)=-(c*(1.0-EXP(-a*t))/a+dv_cb(i)*EXP(-a*t))/dz
  END IF
END DO
!
! Calculate the cloud-base stress components
! Equations 11 & 12 section 5.1
!

DO i=1,n_sh
  uw(i,ntml(i))=zlcl_cmt(i)*(-mb(i)*omg2_jump(i)-                 &
                 GAMMA*uw0(i)/zlcl_cmt(i))/delta+uw0(i)
  vw(i,ntml(i))=zlcl_cmt(i)*(mb(i)*omg1_jump(i)-                  &
                 GAMMA*vw0(i)/zlcl_cmt(i))/delta+vw0(i)
END DO

!
! Calculate non-gradient stress profile in cloud
! Altered numbering of uw arrays
!
! New form of Fcmt   where alpha >1?  (like fng term for thermo?)
! Fcmt(z/zlcd)  = exp(-alpha*(z/zcld))

DO k=1,ntpar_max

  DO i=1,n_sh


    IF (k >= (ntml(i)+1) .AND. k <= (ntpar(i)-1)) THEN
      ! F function

      fcmt=EXP(-1.1*(ztheta(i,k)-zlcl_cmt(i))/z_depth(i))

      ! all cloud levels add non-gradient term to eddy viscosity term

      uw(i,k)=uw(i,k)+uw(i,ntml(i))*fcmt
      vw(i,k)=vw(i,k)+vw(i,ntml(i))*fcmt

    ELSE IF (k <= (ntml(i)-1)) THEN

      ! Which is correct ?
      !            uw(i,k) = uw(i,ntml(i))*ztheta(i,k)/zlcl_cmt(i)
      !            vw(i,k) = vw(i,ntml(i))*ztheta(i,k)/zlcl_cmt(i)
      ! This if assuming uw on rho levels
      uw(i,k) = uw(i,ntml(i))*zrho(i,k)/zlcl_cmt(i)
      vw(i,k) = vw(i,ntml(i))*zrho(i,k)/zlcl_cmt(i)

    END IF

  END DO

END DO

!
! Copy stress to output arrays multiplying by density on theta levels
!
IF (flg_uw_shall) THEN
  DO k=1,ntpar_max+1
    DO i=1,n_sh
      uw_shall(i,k)=uw(i,k)*rho_theta(i,k)
    END DO
  END DO
END IF
IF (flg_vw_shall) THEN
  DO k=1,ntpar_max+1
    DO i=1,n_sh
      vw_shall(i,k)=vw(i,k)*rho_theta(i,k)
    END DO
  END DO
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE shtconv_base_stress
END MODULE shtconv_base_stress_mod

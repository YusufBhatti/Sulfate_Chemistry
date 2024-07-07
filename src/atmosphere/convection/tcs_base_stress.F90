! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing tcs warm rain subroutine for calculating
! cloud base stress
!
MODULE tcs_base_stress


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!
! Description:
!   This routine calculates cloud base stress
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
!

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TCS_BASE_STRESS'

CONTAINS

SUBROUTINE calc_base_stress(n_xx, nlevs, ntml, ntpar,ntpar_max &
   ,                              timestep                     &
   ,                              scales                       &
   ,                              uw0,vw0,du_cb,dv_cb          &
   ,                              rho_theta,zrho,ztheta        &
   ,                              flg_uw_cong,flg_vw_cong      &
                              ! IN/OUT ARGUMENTS
   ,                              uw,vw                        &
                              ! OUTPUT ARGUMENTS
   ,                              uw_cong,vw_cong)


USE tcs_parameters_warm,   ONLY:                             &
   beta_cmt, delta_cmt, gamma_cmt
USE tcs_class_scales,      ONLY:                             &
   scales_conv

IMPLICIT NONE

!--------------------------------------------------------------------
! Subroutine Arguments
!--------------------------------------------------------------------
!
! Arguments with intent IN:
!
INTEGER, INTENT(IN) ::                                            &
   n_xx                                                           &
                            ! Total number of congestus points
   , nlevs                                                        &
                            ! Number of model levels
   , ntml(n_xx)                                                   &
                            ! levels of LCL
   , ntpar(n_xx)                                                  &
                            ! levels of TOP ofcloud layer
   , ntpar_max
! max value of ntpar +1
!
LOGICAL, INTENT(IN) ::                                            &
   flg_uw_cong                                                    &
                            ! STASH FLAGS FOR CONGESTUS
   ,flg_vw_cong
! CONVECTION STRESS DIAGNOSTIC
!
REAL, INTENT(IN) ::                                               &
   timestep         ! MODEL timestep (S)

TYPE(scales_conv), INTENT(IN) :: scales

REAL, INTENT(IN) ::                                               &
   uw0(n_xx)                                                      &
                            ! U-component of surface stress (M2S-2)
   ,  vw0(n_xx)                                                   &
                            ! V-component of surface stress (M2S-2)
   ,  du_cb(n_xx)                                                 &
                            ! dU across cloud base (m/s)
   ,  dv_cb(n_xx)                                                 &
                            ! dV across cloud base (m/s)
   ,  rho_theta(n_xx,nlevs)                                       &
                            ! Density model th levels (kgm-3)
   ,  zrho(n_xx,nlevs)                                            &
                            ! height of model rho levels (m)
   ,  ztheta(n_xx,nlevs)
! height of model theta levels (m)

!
! Arguments with intent INOUT:
!
REAL, INTENT(INOUT) ::                                            &
   uw(n_xx,nlevs)                                                 &
                            ! U-component of STRESS PROFILE (M2S-2)
   , vw(n_xx,nlevs)
! V-component of STRESS PROFILE (M2S-2)
!
! Arguments with intent OUT:
!
REAL, INTENT(OUT) ::                                              &
   uw_cong(n_xx,nlevs)                                            &
                            ! STASH DIAGNOSTIC FOR U-COMP STREss
   , vw_cong(n_xx,nlevs)
! STASH DIAGNOSTIC FOR V-COMP STREss
!
!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

REAL ::                                                           &
   omg2_jump(n_xx)                                                &
                            ! jump in Y component of vorticity
   , omg1_jump(n_xx)                                              &
                            ! jump in X-component of vorticity
   , zlcl_cmt(n_xx)                                               &
                            ! lcl for CMT    (m)
   , fcmt,dz                                                      &
   , a,b,c                                                        &
                            ! coefficients
   , t,dz1                                                        &
   , z_depth(n_xx)                                                &
   , expadt                  ! exp (Adt)

!-------------------------
! Loop counters
!-------------------------
INTEGER :: i,k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_BASE_STRESS'

!-------------------------------------------------------------------
!
! Calculate jumps in vorticity across cloud base. (This is done by
! assuming that during the time step du and dv vary as exp(-T/TAU).
! Needs to be done to avoid instability around cloud base.)
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i=1,n_xx

  zlcl_cmt(i) = ztheta(i,ntml(i))   ! cloud base for CMT
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

  !  f(z/scales%zcld) = exp(-fcmt)
  fcmt=beta_cmt*scales%wsc_o_mb(i)*(dz1/z_depth(i))
  b=(1.0/zlcl_cmt(i)-(EXP(-fcmt)-1.0)/dz1)
  ! alpha in doc (but extra /dz factor)
  a=zlcl_cmt(i)*scales%mb(i)*b/(delta_cmt*dz)
  expadt = EXP(-a*timestep)

  ! beta in documentation

  c = ( b*(1.0-gamma_cmt/delta_cmt) - 1.0/zlcl_cmt(i) )*uw0(i)

  IF (c == 0.0 .AND. du_cb(i) == 0.0) THEN
    omg2_jump=0.0
  ELSE
    t=-LOG( (c*(1.0-expadt)/a+du_cb(i)*(expadt-1.0))/             &
                             ((c-a*du_cb(i))*timestep))/a
    omg2_jump(i)=( c*(1.0-EXP(-a*t))/a + du_cb(i)*EXP(-a*t) )/dz
  END IF

  ! beta in documentation
  c = ( b*(1.0-gamma_cmt/delta_cmt) - 1.0/zlcl_cmt(i) )*vw0(i)

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

DO i=1,n_xx
  uw(i,ntml(i))=zlcl_cmt(i)*(-scales%mb(i)*omg2_jump(i)-          &
                       gamma_cmt*uw0(i)/zlcl_cmt(i))/delta_cmt+uw0(i)
  vw(i,ntml(i))=zlcl_cmt(i)*(scales%mb(i)*omg1_jump(i)-           &
                       gamma_cmt*vw0(i)/zlcl_cmt(i))/delta_cmt+vw0(i)
END DO

!
! Calculate non-gradient stress profile in cloud
! Altered numbering of uw arrays
!
! New form of Fcmt   where alpha >1?  (like fng term for thermo?)
! Fcmt(z/zlcd)  = exp(-alpha*(z/scales%zcld))

DO k=1,ntpar_max

  DO i=1,n_xx


    IF (k >= (ntml(i)+1) .AND. k <= (ntpar(i)-1)) THEN
      ! F function

      fcmt=EXP(-1.1*(ztheta(i,k)-zlcl_cmt(i))/z_depth(i))

      ! all cloud levels add non-gradient term to eddy viscosity term

      uw(i,k)=uw(i,k)+uw(i,ntml(i))*fcmt
      vw(i,k)=vw(i,k)+vw(i,ntml(i))*fcmt

    ELSE IF (k <= (ntml(i)-1)) THEN

      ! This if assuming uw on rho levels
      uw(i,k) = uw(i,ntml(i))*zrho(i,k)/zlcl_cmt(i)
      vw(i,k) = vw(i,ntml(i))*zrho(i,k)/zlcl_cmt(i)

    END IF

  END DO

END DO

!
! Copy stress to output arrays multiplying by density on theta levels
!
IF (flg_uw_cong) THEN
  DO k=1,ntpar_max+1
    DO i=1,n_xx
      uw_cong(i,k)=uw(i,k)*rho_theta(i,k)
    END DO
  END DO
END IF
IF (flg_vw_cong) THEN
  DO k=1,ntpar_max+1
    DO i=1,n_xx
      vw_cong(i,k)=vw(i,k)*rho_theta(i,k)
    END DO
  END DO
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE calc_base_stress

END MODULE tcs_base_stress

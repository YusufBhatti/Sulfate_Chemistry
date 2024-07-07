! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates cloud base stress for shallow convection
!

MODULE shallow_base_stress_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SHALLOW_BASE_STRESS_MOD'
CONTAINS

SUBROUTINE shallow_base_stress(np_field,npnts,n_cumulus,nlevs,        &
                               nterm,                                 &
                               cu_ind,cu_full,nlcl,ntop,mb,wsc,       &
                               zlcl,zcld,uw0,vw0,plcl,ptop,ue,        &
                               ve,phalf,p,rho,timestep,weight_param,  &
                                     ! IN/OUT ARGUMENTS
                               uw,vw,                                 &
                                     ! OUTPUT ARGUMENTS
                               uw_shall,vw_shall)

USE planet_constants_mod, ONLY: g
USE cv_stash_flg_mod,    ONLY: flg_uw_shall, flg_vw_shall
USE model_domain_mod,    ONLY: model_type, mt_single_column

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!------------------------------------------------------------------------
! Description:
!   Calculates cloud base stress for shallow convection
!   (also completes caluclation of stress profile).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-----------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  np_field             & ! Full field length
 ,npnts                & ! Total number of points in segment
 ,n_cumulus            & ! Number of convecting points
 ,nlevs                & ! Number of model levels
 ,nterm                & ! Number of points terminating
 ,cu_ind(nterm)        & ! Indices for terminating points
 ,cu_full(nterm)       & ! Indices of points in output array
 ,nlcl(n_cumulus)      & ! Lifting condensation level
 ,ntop(n_cumulus)        ! Top level of convection

REAL, INTENT(IN)    ::     &
  mb(n_cumulus)            & ! Cloud base mass flux for shallow cu (m/s)
 ,wsc(n_cumulus)           & ! Cloud-layer velocity scale (m/s)
 ,zlcl(npnts)              & ! Height of LCL (m)
 ,zcld(npnts)              & ! cloud layer depth (m)
 ,uw0(npnts)               & ! U component of surface stress (m2/s2)
 ,vw0(npnts)               & ! V component of surface stress (m2/s2)
 ,plcl(n_cumulus)          & ! Pressure at lifting condensation level (Pa)
 ,ptop(n_cumulus)          & ! Pressure at top of cloud layer (Pa)
 ,ue(nlevs,n_cumulus)      & ! U-component of mean wind (m/s)
 ,ve(nlevs,n_cumulus)      & ! V-component of mean wind (m/s)
 ,phalf(nlevs,n_cumulus)   & ! Pressure on model half levels (Pa)
 ,p(nlevs,n_cumulus)       & ! Pressure on model levels (Pa)
 ,rho(nlevs,n_cumulus)     & ! Density, model uv levels (kg/m3/s)
 ,weight_param(n_cumulus)  & ! Grey zone weighting factor
 ,timestep                   ! Model timestep (s)

REAL, INTENT(INOUT) ::     &
  uw(nlevs,n_cumulus)      & ! U-component of stress (m2/s2)
 ,vw(nlevs,n_cumulus)        ! V-component of stress (m2/s2)

REAL, INTENT(OUT) ::       &
  uw_shall(np_field,nlevs) & ! Stash diagnostic for U comp stress
 ,vw_shall(np_field,nlevs)   ! Stash diagnostic for V comp stress

! Local variables

INTEGER ::       &
  i,j,k,n          ! Loop counters


REAL ::                &
  delta_z(nterm)       & ! Layer thickness above LCL (m)
 ,coeff_1(nterm)       & ! Coeffficient
 ,omg2_jump(nterm)     & ! Jump in Y component of vorticity
 ,omg1_jump(nterm)     & ! Jump in X component of vorticity
 ,u_jump(nterm)        & !
 ,v_jump(nterm)        & !
 ,zeta                 & !
 ,dz                   & ! height difference
 ,p_depth              & !
 ,a                    & !
 ,b                    & !
 ,c                    & !
 ,t                    & !
 ,du                   & !
 ,dv                   & !
 ,dz1                  & ! height difference
 ,rho_h                  !


! Parameters (in future consider putting in a module).

REAL, PARAMETER ::  &
  beta=0.04         & !
 ,delta=2.3         & !
 ,GAMMA=1.63

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SHALLOW_BASE_STRESS'

!------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------
! Calculate jumps in vorticity across cloud base (this is done by assuming
! that during the timestep du and dv vary as EXP(-T/TAU). Needs to be done
! to avoid instability around cloud base).

SELECT CASE (model_type)

CASE DEFAULT
  DO i=1, nterm
    j       = cu_ind(i)
    n       = cu_full(i)
    dz      = -(p(nlcl(j),j)-plcl(j))/(g*rho(nlcl(j),j))
    du      = (ue(nlcl(j),j)-ue(nlcl(j)-1,j))
    dv      = (ve(nlcl(j),j)-ve(nlcl(j)-1,j))
    dz1     = -(phalf(nlcl(j)+1,j)-phalf(nlcl(j),j))/(g*rho(nlcl(j)+1,j))
    p_depth = (ptop(j)-plcl(j))
    zeta    = beta*wsc(j)*(phalf(nlcl(j)+1,j)-plcl(j))/(mb(j)*p_depth)
    b       = (1.0/zlcl(n)-(EXP(-zeta)-1.0)/dz1)
    a       = zlcl(n)*mb(j)*b/(delta*dz)
    c       = (b*(1.0-GAMMA/delta)-1.0/zlcl(n))*uw0(n)

    t = -LOG((c*(1.0-EXP(-a*timestep))/a+                          &
          du*(EXP(-a*timestep)-1.0))/                              &
          ((c-a*du)*timestep))/a
    omg2_jump(i) = (c*(1.0-EXP(-a*t))/a+du*EXP(-a*t))/dz

    c = (b*(1.0-GAMMA/delta)-1.0/zlcl(n))*vw0(n)

    t = -LOG((c*(1.0-EXP(-a*timestep))/a+                          &
          dv*(EXP(-a*timestep)-1.0))/                              &
          ((c-a*dv)*timestep))/a
    omg1_jump(i) = -(c*(1.0-EXP(-a*t))/a+dv*EXP(-a*t))/dz
  END DO

CASE (mt_single_column)
  DO i=1, nterm
    j       = cu_ind(i)
    n       = cu_full(i)
    dz      = -(p(nlcl(j),j)-plcl(j))/(g*rho(nlcl(j),j))
    du      = (ue(nlcl(j),j)-ue(nlcl(j)-1,j))
    dv      = (ve(nlcl(j),j)-ve(nlcl(j)-1,j))
    dz1     = -(phalf(nlcl(j)+1,j)-phalf(nlcl(j),j))/(g*rho(nlcl(j)+1,j))
    p_depth = (ptop(j)-plcl(j))
    zeta    = beta*wsc(j)*(phalf(nlcl(j)+1,j)-plcl(j))/(mb(j)*p_depth)
    b       = (1.0/zlcl(n)-(EXP(-zeta)-1.0)/dz1)
    a       = zlcl(n)*mb(j)*b/(delta*dz)
    c       = (b*(1.0-GAMMA/delta)-1.0/zlcl(n))*uw0(n)

    IF (c == 0.0 .AND. (dv == 0.0 .OR. du == 0.0)) THEN
      omg2_jump=0.0
    ELSE
      t = -LOG((c*(1.0-EXP(-a*timestep))/a+                         &
           du*(EXP(-a*timestep)-1.0))/                              &
           ((c-a*du)*timestep))/a
      omg2_jump(i) = (c*(1.0-EXP(-a*t))/a+du*EXP(-a*t))/dz
    END IF
    c = (b*(1.0-GAMMA/delta)-1.0/zlcl(n))*vw0(n)
    IF (c == 0.0 .AND. (dv == 0.0 .OR. du == 0.0)) THEN
      omg1_jump=0.0
    ELSE
      t = -LOG((c*(1.0-EXP(-a*timestep))/a+                         &
           dv*(EXP(-a*timestep)-1.0))/                              &
           ((c-a*dv)*timestep))/a
      omg1_jump(i) = -(c*(1.0-EXP(-a*t))/a+dv*EXP(-a*t))/dz
    END IF

  END DO

END SELECT ! model_type

! Calculate the cloud-base stress components
DO i=1,nterm
  j=cu_ind(i)
  n=cu_full(i)
  uw(nlcl(j),j) = zlcl(n)*( -mb(j)*omg2_jump(i) -                  &
                             GAMMA*uw0(n)/zlcl(n) )/delta + uw0(n)
  vw(nlcl(j),j) = zlcl(n)*( mb(j)*omg1_jump(i) -                   &
                            GAMMA*vw0(n)/zlcl(n)  )/delta + vw0(n)
  uw(nlcl(j),j) = weight_param(j)*uw(nlcl(j),j)
  vw(nlcl(j),j) = weight_param(j)*vw(nlcl(j),j)
END DO

! Calculate non-gradient stress profile

DO i=1,nterm
  j=cu_ind(i)
  n=cu_full(i)
  p_depth=(ptop(j)-plcl(j))
  DO k=nlcl(j)+1,ntop(j)+1
    zeta=beta*wsc(j)*(phalf(k,j)-plcl(j))/(mb(j)*p_depth)
    rho_h=rho(k-1,j)+(rho(k,j)-rho(k-1,j))/(p(k,j)-p(k-1,j))*       &
                                     (phalf(k,j)-p(k-1,j))
    IF (k <  ntop(j)) THEN
      uw(k,j)=rho_h*(uw(k,j)+uw(nlcl(j),j)*EXP(-zeta))
      vw(k,j)=rho_h*(vw(k,j)+vw(nlcl(j),j)*EXP(-zeta))
    ELSE IF (k == ntop(j)) THEN
      uw(k,j)=rho_h*(uw(k,j)+uw(nlcl(j),j)*EXP(-beta*wsc(j)/mb(j)))
      vw(k,j)=rho_h*(vw(k,j)+vw(nlcl(j),j)*EXP(-beta*wsc(j)/mb(j)))
    ELSE IF (k == ntop(j)+1) THEN
      uw(k,j)=rho_h*(uw(k,j)+uw(nlcl(j),j)*                          &
                          EXP(-1.75*beta*wsc(j)/mb(j)))
      vw(k,j)=rho_h*(vw(k,j)+vw(nlcl(j),j)*                          &
                          EXP(-1.75*beta*wsc(j)/mb(j)))
    END IF
  END DO

  ! Weight cloud base stress by rho (omitted from above level loop)
  ! Needs to be done after level loop.

  k=nlcl(j)
  rho_h=rho(k-1,j)+(rho(k,j)-rho(k-1,j))/(p(k,j)-p(k-1,j))*        &
                            (phalf(k,j)-p(k-1,j))
  uw(nlcl(j),j) = rho_h*uw(nlcl(j),j)
  vw(nlcl(j),j) = rho_h*vw(nlcl(j),j)


  IF (flg_uw_shall) THEN
    DO k=nlcl(j),ntop(j)+1
      uw_shall(n,k)=uw(k,j)
    END DO
  END IF
  IF (flg_vw_shall) THEN
    DO k=nlcl(j),ntop(j)+1
      vw_shall(n,k)=vw(k,j)
    END DO
  END IF
END DO   ! nterm

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE shallow_base_stress
END MODULE shallow_base_stress_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Large-scale precipitation scheme. Moist Froude calculation

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

! Subroutine Interface:
MODULE lsp_froude_moist_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_FROUDE_MOIST_MOD'

CONTAINS

SUBROUTINE lsp_froude_moist(levels, points,                                   &
                           u, v, rho, theta,                                  &
                           rho_levels, theta_levels,                          &
                           hmteff, zb)

! Code to calculate moist Brunt-Vaisala frequency squared
! which is then used to calculate the moist Froude number
! and the amount of blocking
! Copied gw_setup.F90  gw_block.F90 from gravity_wave_drag
! and modified

USE planet_constants_mod, ONLY: g
USE mphys_inputs_mod,     ONLY: fcrit

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE c_gwave_mod, ONLY: lambdaz_min, lambdaz_max,                              &
                       nsq_neutral, zav_converge, zav_iterate

IMPLICIT NONE

!----------------------
! Intent(in) variables
!----------------------
INTEGER, INTENT(IN) :: &! INPUT array dimensions
  levels,              &!num of vertical levels
  points                !num of land points

REAL, INTENT(IN) ::           &! ARRAYS
  u(points,levels),           &!zonal wind (on P levels)
  v(points,levels),           &!meridional wind (on P levels)
  rho(points,levels),         &!density
  theta(points,levels),       &!potential temperature
  theta_levels(points,levels), &!height(theta_levels)
  rho_levels(points,levels)    !height (rho_levels)

!----------------------
! Intent(inout) variables
!----------------------
REAL, INTENT(INOUT) ::                                                        &
  zb(points),            & !Depth of flow blocking layer
  hmteff(points)           !Effective mountain height producing ascent

! Local variables
REAL ::                                                                       &
  mt_high(points),    &!sso height (n_sigma*sd_orog)
  nsq(points,levels), &!brunt-vaisala freq squared
  ulow(points),       &!u averaged from z=0.5mt_high to z=mt_high
  vlow(points),       &!v averaged from z=0.5mt_high to z=mt_high
  modu(points),       &!modulus of horizontal wind (i.e. sqrt(u^2+v^2)
  nlow(points),       &!N bulk averaged from z=0.5mt_high to z=mt_high
  rholow(points),     &!rho averaged from z=0.5mt_high to z=mt_high
  zneu(points),       &!depth of near surface neutral layer
  zav(points),        &!depth used to calc fav
  zav1(points),       &!used to calculate zav1
  zav_new(points),    &!used in calculation of zav
  u_n(points),        &!wind speed div by buoyancy freq
  nav(points),        &!bulk averaged n^2 from z=0 to z=zav
  uav(points),        &!u averaged from z=0 to z=zav
  vav(points),        &!v averaged from z=0 to z=zav
  fav(points)          !Froude number for calc zb
                       ! (zb=max(0,mt_high(fcrit-fav))

REAL ::                                                                       &
 wind,        &!wind speed resolved in direction of low-level flow.
 dzt,         &!layer thickness used to calc average u etc.
 dzb,         &!layer thickness used to calc average u etc.
 dzu

INTEGER :: i, k, ii, t
INTEGER :: kbot(points)
INTEGER :: ktop(points)    ! From setup
INTEGER :: ktop_fb(points) ! From block

LOGICAL :: l_sfblock(points),  &   !whether point has a non-zero stress or not
           l_cont(points),                                                    &
           l_cont2(points),                                                   &
           l_cont3(points)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_FROUDE_MOIST'
! ----------------------------------------
! initialise arrays
! ----------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise BV freq profile
DO k = 1, levels
  DO i =1, points
    nsq(i,k)=0.0
  END DO
END DO

! Initialise low-level parameters
DO i=1,points
! Mountain height
  mt_high(i) = hmteff(i)
! Other local variables
  ulow(i)    = 0.0
  vlow(i)    = 0.0
  rholow(i)  = 0.0
  nlow(i)    = 0.0
  ktop(i)    = 2
  ktop_fb(i) = 2
  kbot(i)    = 1
  zav_new(i) = 0.0
  zneu(i)    = 0.0
  uav(i)     = 0.0
  vav(i)     = 0.0
  nav(i)     = 0.0
  l_sfblock(i)  = .TRUE.
  !u(i,1)=0 implies at north/south pole in global model
  IF ((mt_high(i) <= 0.0) .OR. (u(i,1) == 0.0)) THEN
    l_sfblock(i) = .FALSE.
  END IF
  l_cont(i)  = .TRUE.
  l_cont2(i) = .TRUE.
  l_cont3(i) = .TRUE.
END DO


! Code copied from gw_setup
! --------------------------
! Estimate N squared at each level - DRY
  DO k = 2, levels-1
    DO i = 1, points
      IF (l_sfblock(i)) THEN
        nsq(i,k) = 2.0*g*  ( theta(i,k) - theta(i,k-1) )                      &
                     /( ( theta(i,k) + theta(i,k-1) )                         &
                       *( theta_levels(i,k) - theta_levels(i,k-1) ) )
      END IF!(l_sfblock(i)).
    END DO !i = 1, points
  END DO !k= 2, levels-1

  DO i = 1, points
    !Set nsq(1)=nsq(2) as nsq is undefined on level 1
    nsq(i,1)       = nsq(i,2)
    !Set nsq(levels)=nsq(levels-1) as nsq is undefined on top level
    nsq(i,levels)  = nsq(i,levels-1)
  END DO !i = 1, points

! Get bottom and top of averaging layer (top half mountain)
  DO k = 1,levels
    DO i = 1, points
      IF (l_sfblock(i)) THEN
        IF (theta_levels(i,k) <= 0.5*mt_high(i)) THEN
          kbot(i) = k
        END IF
        IF (theta_levels(i,k) < mt_high(i)) THEN
          ktop(i) = k+1
        END IF
      END IF !(l_sfblock(i)).
    END DO !i = 1, points
  END DO !k = 1,levels

! Find low-level averages (top half mountain)
  DO k = 2, levels-1
    DO i = 1, points
      IF (l_sfblock(i)) THEN
        IF ((k > kbot(i)) .AND. (k <= ktop(i))) THEN
          dzt = theta_levels(i,k)   -  0.5 * mt_high(i)
          dzb = theta_levels(i,k)   -  theta_levels(i,k-1)
          IF (k  ==  kbot(i)+1) THEN
            dzu = 0
            dzb = dzt
          ELSE
            dzu = theta_levels(i,k-1) - 0.5 * mt_high(i)
          END IF ! (k == bot(i)+1)
          IF (k == ktop(i)) THEN
            IF (k  ==  kbot(i)+1) THEN
              dzb = dzt
            ELSE
              dzt = mt_high(i) -  0.5 * mt_high(i)
              dzb = mt_high(i) -  theta_levels(i,k-1)
            END IF !(k  ==  kbot(i)+1)
          END IF ! (k == ktop(i))
          !------------------------------------------------------
          ! average u,v and rho from z = 0.5h to current level
          !-----------------------------------------------------
          ulow(i)  = (ulow(i)   * dzu  + u(i,k)  *dzb)/dzt
          vlow(i)  = (vlow(i)   * dzu  + v(i,k)  *dzb)/dzt
          rholow(i)= (rholow(i) * dzu  + rho(i,k)*dzb)/dzt
          ! bulk average N from z = 0.5mt_high tO z = mt_high
          nlow(i)  = (nlow(i)   * dzu  + nsq(i,k)*dzb)/dzt
        END IF ! (k > kbot(i) .AND k <= ktop(i))
      END IF !(l_sfblock(i))
    END DO !i = 1, points
  END DO  ! k = 2, levels-1

  DO i = 1, points
    IF (l_sfblock(i)) THEN

      modu(i)    = SQRT(ulow(i)**2 + vlow(i)**2)

      IF ( modu(i) <= 0.0 ) THEN
        l_sfblock(i) = .FALSE.
        zb(i) = mt_high(i)
        hmteff(i) = 0.0
        fav(i) = -2
      END IF

      IF (nlow(i) > 0.0) THEN
        nlow(i) = SQRT(nlow(i))
      ELSE
        nlow(i)   = 0.0
        l_sfblock(i) = .FALSE.
      END IF !(nlow(i) > 0.)

    END IF !(l_sfblock(i)).
  END DO !i=1,points


! Code from gw_block to get Froude and blocking
! ---------------------------------------------

! ------------------------------------------------
!  calculate zav following Vosper et al (2009)
! ------------------------------------------------
! calculate depth of neutral layer - zneu
  DO k = 1, levels
    DO i = 1, points
      IF (l_sfblock(i)) THEN
        IF (l_cont(i)) THEN
          IF (nsq(i,k) < nsq_neutral) THEN
            zneu(i)    = rho_levels(i,k)
          ELSE
            l_cont(i)  = .FALSE.
          END IF !nsq(i,k) < nsq_neutral
        END IF !l_cont
      END IF !l_sfblock(i)
    END DO!i = 1, points
  END DO !k = 1, levels

! find max(mt_high,zneu) for zav
  DO i = 1, points
    IF (l_sfblock(i)) THEN
      IF (zneu(i) > mt_high(i)) THEN
        zav1(i) = zneu(i)
      ELSE
        zav1(i) = mt_high(i)
      END IF
      zav(i)     = zav1(i)
    END IF!l_sfblock(i)
  END DO !i = 1, points

  DO ii = 1, zav_iterate
    !Reset l_cont3 (logical to test if rho_levels(k) lt zav)
    DO i = 1, points
      l_cont3(i) = .TRUE.
    END DO !i = 1, points
    DO k = 2, levels-1
      DO i = 1, points
        IF (l_sfblock(i)) THEN
          IF (l_cont2(i)) THEN !l_cont2(i) tests if zav is converged
            !-----------------------------------------------------
            ! need an if test for case where zav doesn't
            ! converge in zav_iterate iterations?
            !-----------------------------------------------------
            IF (l_cont3(i)) THEN !l_cont3(i) tests if rho_levels(k) lt zav
              IF (theta_levels(i,k) < zav(i)) THEN
                dzt = theta_levels(i,k)
                IF (k == 2) THEN
                  dzb    = dzt
                ELSE
                  dzb    = theta_levels(i,k) -  theta_levels(i,k-1)
                END IF !(k==2)
              ELSE
                dzt    = zav(i)
                IF (k == 2) THEN
                  dzb    = dzt
                ELSE
                  dzb     = zav(i) - theta_levels(i,k-1)
                  ktop_fb(i) = k
                  l_cont3(i)  = .FALSE.
                END IF !(k==2)
              END IF !theta_levels(i,k)<zav(i)
              !-----------------------------------------------------
              ! average u,v and n from z = 0 to current level
              ! which when k = ktop_fb become z = 0 tO z = zav(i) averages
              ! nav(i) calc is equivalent to bulk average n
              ! i.e. n^2 = sqrt(g/theta0*thetaav-theta0/zav(i))
              !-----------------------------------------------------
              uav(i) = (uav(i)*theta_levels(i,k-1) + u(i,k)*dzb)  / dzt
              vav(i) = (vav(i)*theta_levels(i,k-1) + v(i,k)*dzb)  / dzt
              nav(i) = (nav(i)*theta_levels(i,k-1) + nsq(i,k)*dzb)/ dzt
            END IF    ! l_cont3
          END IF    ! l_cont2
        END IF !l_sfblock
      END DO ! i=1, points
    END DO   ! loop over k levels

    DO i = 1, points
      IF (l_sfblock(i)) THEN
        IF (l_cont2(i)) THEN !l_cont2(i) tests if zav is converged
          ! resolve wind in the direction of the low-level flow
          wind    = (ulow(i)*uav(i) + vlow(i)*vav(i)) /modu(i)
          wind    = ABS(wind)
          IF (nav(i) > nsq_neutral) THEN
            ! first consider stable cases
            u_n(i) = wind/SQRT(nav(i))
          ELSE
            ! now consider neutral (and near neutral) cases
            u_n(i) = wind/SQRT(nsq_neutral)
          END IF
          ! limit u_n to sensible values
          u_n(i) = MAX(lambdaz_min,u_n(i))
          u_n(i) = MIN(lambdaz_max,u_n(i))
          zav_new(i) = zav1(i) + u_n(i)
          ! currently set zav_converge to 0.05
          IF ((zav(i) < zav_new(i)*(1.0+zav_converge)) .AND.                  &
             ( zav(i) > zav_new(i)*(1.0-zav_converge))) THEN
            l_cont2(i)  =  .FALSE.
          END IF ! test if zav(i) is converged
          zav(i)    = zav_new(i)
        END IF  ! l_cont2(i)
      END IF!l_sfblock(i)
    END DO!i = 1, points
  END DO  !ii= 1, zav_iterate

! -----------------------------------------
!  calculate fav and blocked layer depth
! -----------------------------------------
  DO i = 1, points
    l_cont(i)  =  .TRUE.
    IF (l_sfblock(i)) THEN
      ! calculate froude number
      ! prevent div by zero if nav(i) =0.
      IF (nav(i) > 0.0) THEN
        !limit fav with wavelength too
        fav(i) = u_n(i)/mt_high(i)
      ELSE
        fav(i) = -1.0
      END IF
      ! find zb and variable drag coefficient
      IF ((fav(i) > 0.0) .AND. ((fcrit - fav(i)) > 0.0)) THEN
        zb(i) = mt_high(i)*(1 - fav(i)/fcrit)
        hmteff(i) = mt_high(i) - zb(i)
      ELSE
        zb(i) = 0.0
        hmteff(i) = mt_high(i)
      END IF
      ! find level at top of blocked layer
      ktop_fb(i) = 0
    END IF!l_sfblock(i)

  END DO !i = 1, points

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_froude_moist
END MODULE lsp_froude_moist_mod

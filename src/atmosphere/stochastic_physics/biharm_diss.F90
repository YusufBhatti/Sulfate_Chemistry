! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:  Biharmonic numerical dissipation
!     This code provides a routine to compute the biharmonic
!     disipation which is the modulus of the gradient of vorticity
!     times K constant. This way to compute the numerical dissipation
!     produces stronger resolution dependence.
!     Code Owner: Please refer to the UM file CodeOwners.txt
!     This file belongs in section: Stochastic Physics
MODULE biharm_diss_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='BIHARM_DISS_MOD'

CONTAINS

SUBROUTINE biharm_diss(rows, row_length, n_rows                         &
,                      model_levels,nbot,ntop                           &
,                      u, v, w, rho                                     &
,                      delta_lambda, delta_phi                          &
,                      disp                                             &
,                      first_atmstep_call, timestep)

! For swap_bounds in EG
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE Field_Types

! Bounds for arrays
USE atm_fields_bounds_mod, ONLY:                                        &
    pdims, pdims_s,                                                     &
    udims, udims_s,                                                     &
    vdims, vdims_s,                                                     &
    wdims_s

! Model grid trigonometry
USE trignometric_mod, ONLY:                                             &
    cos_theta_latitude, cos_v_latitude

! Model level-height modules
USE level_heights_mod,     ONLY:                                        &
    r_rho_levels,r_theta_levels

USE UM_ParVars
USE UM_ParParams, ONLY: pnorth, psouth

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                  &
  row_length                                                            &
           ! Number of points on a row.
, rows                                                                  &
           ! Number of rows.
, n_rows                                                                &
           ! Number of v rows.
, model_levels                                                          &
           ! Number of model levels.
, nbot, ntop
           ! Bottom and top for model levels
           ! to perform biharmonic dissp

LOGICAL, INTENT (IN) ::                                                 &
  first_atmstep_call
           ! Indicates first time-step

REAL, INTENT (IN) ::                                                    &

  u            (udims_s%i_start:udims_s%i_end,                          &
                udims_s%j_start:udims_s%j_end,                          &
                udims_s%k_start:udims_s%k_end)                          &
           ! main prog. u field
, v            (vdims_s%i_start:vdims_s%i_end,                          &
                vdims_s%j_start:vdims_s%j_end,                          &
                vdims_s%k_start:vdims_s%k_end)                          &
           ! main prog. v field
, w            (wdims_s%i_start:wdims_s%i_end,                          &
                wdims_s%j_start:wdims_s%j_end,                          &
                wdims_s%k_start:wdims_s%k_end)                          &
           ! main prog. w field
, rho          (pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,                              &
                pdims%k_start:pdims%k_end)                              &
           ! Rho for divergence
, delta_lambda                                                          &
, delta_phi                                                             &
, timestep
          ! atmosphere model timestep (in seconds)

 ! Variables with Intent (Out)
REAL, INTENT(OUT) ::                                                    &
 disp         (pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)
          ! Numerical Dissipation (m**2.s**-3)

REAL ::                                                                 &

 vort(vdims_s%i_start:vdims_s%i_end,                                    &
              vdims_s%j_start:vdims_s%j_end,                            &
              vdims_s%k_start:vdims_s%k_end)                            &
           ! 'vertical' component of vorticity
, div(pdims_s%i_start:pdims_s%i_end,                                    &
              pdims_s%j_start:pdims_s%j_end,                            &
              pdims_s%k_start:pdims_s%k_end)
           ! Divergence of winds

INTEGER ::                                                              &
    i,j,k                                                               &
            ! Loop indexes
,   ip1                                                                 &
            ! i plus one
,   im1                                                                 &
            ! i minus one
,   jp1                                                                 &
            ! j plus one
,   jm1                                                                 &
            ! j minus one
,   ismooth
            ! loop index over smoothing iterations
REAL, ALLOCATABLE, SAVE ::                                              &
  dx_v(:,:,:)                                                           &
         ! dx on v points
, dy_u(:,:,:)                                                           &
         ! dy on u points
, dx_u(:,:,:)                                                           &
         ! dx on u points
, dy_v(:,:,:)                                                           &
         ! dy on v points
, latitude_mask(:,:)
           ! Theta-level thickness
           ! Mask-out high latitudes

REAL ::                                                                 &
   biharmonic_x, biharmonic_y, biharmonic_z
   ! Coefficients of the biharmonic matrix for x and y
REAL, SAVE :: amp_K !, amp_z
!  Amplitude of viscosity coefficient K eq to 3/ (128 x deltaT)
! and for z eq alpha / (12 deltaT )

! Variables for polar smoothing
REAL ::                                                                 &
  piover2                                                               &
           ! half of PI
, maxlat                                                                &
           ! Lat threshold for truncation of shear to zero at pole
, mylat
           ! Work value of Lat in Radians

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BIHARM_DISS'

!---------------------------------------------------------------------
! END OF VARIABLES DECLARATION. PROGRAM STARTS HERE
!---------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------
! Declare grid spacing variables
!----------------------------------------------------------------

  ! dx_v on v points
IF (.NOT. ALLOCATED(dx_v)) THEN
  ALLOCATE(dx_v(vdims_s%i_start:vdims_s%i_end,                          &
                vdims_s%j_start:vdims_s%j_end,                          &
                vdims_s%k_start:vdims_s%k_end))
END IF

  ! dy_u on u points
IF (.NOT. ALLOCATED(dy_u)) THEN
  ALLOCATE(dy_u(udims_s%i_start:udims_s%i_end,                          &
                udims_s%j_start:udims_s%j_end,                          &
                udims_s%k_start:udims_s%k_end))
END IF

! dx_u on u points
IF (.NOT. ALLOCATED(dx_u)) THEN
  ALLOCATE(dx_u(udims_s%i_start:udims_s%i_end,                          &
                udims_s%j_start:udims_s%j_end,                          &
                udims_s%k_start:udims_s%k_end))
END IF

! dy_v on v points
IF (.NOT. ALLOCATED(dy_v)) THEN
  ALLOCATE(dy_v(vdims_s%i_start:vdims_s%i_end,                          &
                vdims_s%j_start:vdims_s%j_end,                          &
                vdims_s%k_start:vdims_s%k_end))
END IF

!----------------------------------------------------------------
! Calculate grid spacings
!----------------------------------------------------------------
! Horizontal grid spacings used in calculation of rate of strain
! term delta_lambda and delta_phi are in radians

IF (first_atmstep_call) THEN

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i, j, k)                          &
!$OMP& SHARED(vdims,udims,dx_v,dy_v,dx_u,dy_u,r_rho_levels,delta_phi,   &
!$OMP&    cos_v_latitude,cos_theta_latitude,delta_lambda,at_extremity,  &
!$OMP&    udims_s,vdims_s)

  ! For dx in V points
  ! For dy in V points
!$OMP DO SCHEDULE(STATIC)  
  DO k =  vdims%k_start,vdims%k_end
    DO j = vdims%j_start,vdims%j_end
      DO i = vdims%i_start,vdims%i_end
        dx_v(i,j,k) = r_rho_levels(i,j,k) *                             &
                      delta_lambda * cos_v_latitude(i,j)

        dy_v(i,j,k) = 0.5 * ( r_rho_levels(i,j,k) +                     &
                      r_rho_levels(i,j-1,k) ) * delta_phi
      END DO
    END DO

    ! Set vorticity to zero at Psouth and Pnorth
    IF (at_extremity(psouth)) THEN
      DO i = vdims%i_start,vdims%i_end
        dx_v(i,vdims%j_start,k)=0.0
        dy_v(i,vdims%j_start,k)=0.0
      END DO
    END IF

    IF (at_extremity(pnorth)) THEN
      DO i = vdims%i_start,vdims%i_end
        dx_v(i,vdims%j_end,k)=0.0
        dy_v(i,vdims%j_end,k)=0.0
      END DO
    END IF
  END DO   !k
!$OMP END DO 

  ! swap bounds for dx_v and dy_v
!$OMP MASTER
  CALL swap_bounds(dx_v,                                                 &
                   vdims_s%i_len - 2*vdims_s%halo_i,                     &
                   vdims_s%j_len - 2*vdims_s%halo_j,                     &
                   vdims_s%k_len,                                        &
                   vdims_s%halo_i, vdims_s%halo_j,                       &
                   fld_type_v, swap_field_is_scalar)
  CALL swap_bounds(dy_v,                                                 &
                   vdims_s%i_len - 2*vdims_s%halo_i,                     &
                   vdims_s%j_len - 2*vdims_s%halo_j,                     &
                   vdims_s%k_len,                                        &
                   vdims_s%halo_i, vdims_s%halo_j,                       &
                   fld_type_v, swap_field_is_scalar)
!$OMP END MASTER

!$OMP BARRIER

    ! For dy in U points
    ! For dx in U points
!$OMP DO SCHEDULE(STATIC)  
  DO k =  udims%k_start,udims%k_end
    DO j = udims%j_start,udims%j_end
      DO i = udims%i_start,udims%i_end
        dy_u(i,j,k) = r_rho_levels(i,j,k) * delta_phi

        dx_u(i,j,k) = r_rho_levels(i,j,k) *                             &
                      delta_lambda * cos_theta_latitude(i,j)
      END DO
    END DO
  END DO
!$OMP END DO 

  ! swap bounds for dy_v and dx_u
!$OMP MASTER
  CALL swap_bounds(dx_u,                                                 &
                   udims_s%i_len - 2*udims_s%halo_i,                     &
                   udims_s%j_len - 2*udims_s%halo_j,                     &
                   udims_s%k_len,                                        &
                   udims_s%halo_i, udims_s%halo_j,                       &
                   fld_type_u, swap_field_is_scalar)
  CALL swap_bounds(dy_u,                                                 &
                   udims_s%i_len - 2*udims_s%halo_i,                     &
                   udims_s%j_len - 2*udims_s%halo_j,                     &
                   udims_s%k_len,                                        &
                   udims_s%halo_i, udims_s%halo_j,                       &
                   fld_type_u, swap_field_is_scalar)
!$OMP END MASTER

!$OMP END PARALLEL

  ! +++++++++ Declare K constants for diffusion  ++++++++++

  !  Define amplitude for K 3/128 / timestep
  amp_K=3.0/ (128.0* timestep)

  ! ++++++++ Apply polar filtering +++++
  ! Mask for latitudes > MAXLAT deg
  ! Decreases field with a COSINE-shaped drop off across a half-wave
  !  cycle from 1 at maxlat to 0 at pole
  piover2=ACOS(0.0)
  maxlat=1.46   ! ~84 deg = Arbitrary choice based on energy plots

  IF (.NOT. ALLOCATED(latitude_mask)) THEN
    ALLOCATE(latitude_mask(pdims%i_start:pdims%i_end,                   &
                         pdims%j_start:pdims%j_end))
    ! Default global value of one
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        latitude_mask(i,j) = 1.0
      END DO
    END DO

    ! Apply ramp near poles
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(j,i,mylat)                     &
!$OMP& SHARED(pdims,piover2,maxlat,cos_theta_latitude,latitude_mask)    &
!$OMP& SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      mylat = (piover2 - MAX(maxlat,                                    &
      ACOS(cos_theta_latitude(pdims%i_start,j))))/(piover2 - maxlat)
      mylat = 2.0 * piover2 * mylat     ! convert span over 0 <-> 180 deg
      DO i = pdims%i_start, pdims%i_end
        latitude_mask(i,j) = 0.5 - 0.5 * COS(mylat)
      END DO
    END DO
!$OMP END PARALLEL DO

    ! Set value to zero at pole (ND grid) or last theta-row (EG grid)
    ! to avoid spurious turbulence values close to pole
    IF (at_extremity(psouth)) THEN
      DO i = pdims%i_start, pdims%i_end
        latitude_mask(i, pdims%j_start) = 0.0
      END DO
    END IF
    IF (at_extremity(pnorth)) THEN
      DO i = pdims%i_start, pdims%i_end
        latitude_mask(i, pdims%j_end) = 0.0
      END DO
    END IF
  END IF ! end latitude mask computation

END IF    ! firstimestep


!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i, j, k, biharmonic_x,            &
!$OMP& biharmonic_y,jm1,jp1,im1,ip1)                                    &
!$OMP& SHARED(nbot,ntop,vdims,vort,v,u,delta_lambda,dx_u,dx_v,dy_u,     &
!$OMP&    cos_theta_latitude,delta_phi,cos_v_latitude,r_rho_levels,dy_v,&
!$OMP&    at_extremity,vdims_s,pdims,div,rho,w,r_theta_levels,pdims_s,  &
!$OMP&    amp_k,latitude_mask,disp)

! --------------------------------------------------------------------
! Compute Vorticity field.
! --------------------------------------------------------------------
!   compute vorticity in an Arakawa C staggered grid:
!
!
!         |------ v(i,j) --------|
!         |                      |
!      u(i-1,j)                u(i,j)
!         |                      |
!         |------- v(i,j-1)------|
!
! In order to compute the gradients we need to take +1 instead of -1
! in order to put the vorticity point over the i,j ->
!
!  vort=dv/dx-du/dy=[v(i+1,j)-v(i,j)]/dx-[u(i,j+1)-u(i,j)]/dy
!
!
!                       u(i,j+1)
!                          |
!             v(i,j) --- vort(i,j)--- v(i+1,j)
!                          |
!                       u(i,j)
!
! Vorticity has the same dimension as v since it is over the same
! latitudinal grid, however it longitudinally lies over the u grid.

!$OMP DO SCHEDULE(STATIC)
DO k = nbot,ntop
  DO j = vdims%j_start,vdims%j_end
    DO i = vdims%i_start,vdims%i_end
      vort(i,j,k) = (v(i+1,j,k) - v(i,j,k))/delta_lambda -              &
                    (u(i,j+1,k)*cos_theta_latitude(i,j+1) -             &
                      u(i,j,k)*cos_theta_latitude(i,j))/delta_phi
      vort(i,j,k) =  vort(i,j,k) /                                      &
                     (r_rho_levels(i,j,k)*cos_v_latitude(i,j))
    END DO  ! i
  END DO    ! j

  ! Set vorticity to zero at Psouth and Pnorth
  IF (at_extremity(psouth)) THEN
    DO i = vdims%i_start,vdims%i_end
      vort(i,vdims%j_start,k) = 0.0
    END DO
  END IF

  IF (at_extremity(pnorth)) THEN
    DO i = vdims%i_start,vdims%i_end
      vort(i,vdims%j_end,k) = 0.0
    END DO
  END IF

END DO      ! k
!$OMP END DO

 ! Swap bounds in the vorticity field to get values from other
 ! processors' boxes
!$OMP MASTER
CALL swap_bounds(vort,                                                 &
                 vdims_s%i_len - 2*vdims_s%halo_i,                     &
                 vdims_s%j_len - 2*vdims_s%halo_j,                     &
                 vdims_s%k_len,                                        &
                 vdims_s%halo_i, vdims_s%halo_j,                       &
                 fld_type_v, swap_field_is_scalar)
!$OMP END MASTER

! --------------------------------------------------------------------
! Compute Divergence field.
! --------------------------------------------------------------------
!  div = du/dx + dv/dy = [u(i,j)-u(i-1,j)]/dx + [v(i,j)-v(i,j-1)]/dy =
!       = - 1/rho x d(rho x w)/dz
!
! The vertical staggering of UM is as follows:
!
!                        w(k)    theta k level
!                          |
!                       rho(k)   Rho k level
!                          |
!                        w(k-1)  thete k-1 level
!
!   Rho must be averaged to theta levels so
! d(rho x w ) = rho(k+1/2) w(k) - rho(k-1/2) w(k-1)

!$OMP DO SCHEDULE(STATIC)
  ! +1 and -1 to avoid taking rho(nbot-1) ro rho(ntop+1) which are zero
DO k = nbot+1,ntop-1
  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      div(i,j,k)= ( (rho(i,j,k+1) + rho(i,j,k) )* w(i,j,k) -            &
                    (rho(i,j,k) + rho(i,j,k-1))*w(i,j,k-1) ) * 0.5      &
                   ! Divide by delta_Z for theta levels
                  / (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1) )

      div(i,j,k)= div(i,j,k) /rho(i,j,k) * (-1)

    END DO  ! i
  END DO    ! j
END DO      ! k
!$OMP END DO 

! Set dissipation to zero at ntop and nbot
!$OMP MASTER
DO j = pdims%j_start,pdims%j_end
  DO i = pdims%i_start,pdims%i_end
    div(i,j,ntop)= 0.0
    div(i,j,nbot)= 0.0
  END DO    ! i
END DO      ! j

 ! Swap bounds in the divergece field to get values from other
 ! processors' boxes
CALL swap_bounds(div,                                                  &
                 pdims_s%i_len - 2*pdims_s%halo_i,                     &
                 pdims_s%j_len - 2*pdims_s%halo_j,                     &
                 pdims_s%k_len,                                        &
                 pdims_s%halo_i, pdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
!$OMP END MASTER

!$OMP BARRIER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! Compute Biharmonic diffusion K x |grad(vor) + grad(div) |**2
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!$OMP DO SCHEDULE(STATIC)
DO k = nbot,ntop
  DO j = pdims%j_start,pdims%j_end
    jm1 = j - 1
    jp1 = j + 1
    DO i = pdims%i_start,pdims%i_end
      im1 = i - 1
      ip1 = i + 1

      ! K is the biharmonic viscosity tensor. K = 3/128 x delta_i^4 / delta_t
      ! The tensorial decomposition of Dnum is:
      ! Dnum= Kx {(d(vor)/dx)**2 + (d(div)/dx)**2 ) +
      !       Ky {(d(vor)/dy)**2 + (d(div)/dy)**2 ) +
      !
      ! Derivates are averaged between the boxes and  dx^2 from derivates are
      ! cross out. The final expresion is thus
      ! Dnum_x=[ d(vor)**2 + d(div)**2 ) ] * (dx)**2
      !
      ! d(dvor)/dx is over v points so it does need averaging between
      ! j and j-1 to place it into p points
      ! d(dvor)/dy is over u points so it does need averaging between
      ! i and i-1 to place it into p points.
      ! d(div)/dx is on u points so it does averaging over i and i-1
      ! d(div)/dy is on v points so it does averaging over j and j-1


      biharmonic_x = 0.5*( ((vort(i,j,k)-vort(im1,j,k))*dx_v(i,j,k))**2 &
               + ((vort(i,jm1,k)-vort(im1,jm1,k))*dx_v(i,jm1,k))**2 )   &
               + 0.5 * ( ((div(ip1,j,k)-div(i,j,k))*dx_u(i,j,k))**2     &
               + ((div(i,j,k)-div(im1,j,k))*dx_u(i,j,k))**2 )

      biharmonic_y = 0.5*( ((vort(i,j,k)-vort(i,jm1,k))*dy_u(i,j,k))**2 &
                + ((vort(im1,j,k)-vort(im1,jm1,k))*dy_u(im1,j,k))**2 )  &
                + 0.5 * ( ((div(i,jp1,k)-div(i,j,k))*dy_v(i,j,k))**2    &
                + ((div(i,j,k)-div(i,jm1,k))*dy_v(i,jm1,k))**2 )


      ! Sum both terms for each coordinate
      disp(i,j,k)= (biharmonic_x + biharmonic_y ) * amp_K

      ! apply latitude mask
      disp(i,j,k)=disp(i,j,k) * latitude_mask(i,j)

    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE biharm_diss
END MODULE biharm_diss_mod

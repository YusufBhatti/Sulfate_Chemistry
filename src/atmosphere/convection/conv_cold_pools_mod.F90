! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
! Scheme to model behaviour of convective cold pools.
!   
! Method:
! Based on gravity-current models e.g. Rooney (J.Fluid Mech., 2015).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
MODULE conv_cold_pools_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CONV_COLD_POOLS_MOD'
CONTAINS

! Code Description:
!   Programming standard : UMDP 3
SUBROUTINE conv_cold_pools( global_row_length, global_rows, ddmfx,      &
                            u_steer, v_steer, cca_steer,                &
                            ux_ccp, uy_ccp, um_ccp, g_ccp, h_ccp,       &
                            riso_ccp, rdir_ccp )

USE conversions_mod,        ONLY : pi, pi_over_180
USE planet_constants_mod,   ONLY : g, planet_radius
USE horiz_grid_mod,         ONLY : xi1_u, xi2_v, csxi2_p
USE atm_fields_bounds_mod,  ONLY : pdims, pdims_s
USE timestep_mod,           ONLY : timestep
USE edge_exchange_mod,      ONLY : edge_exchange
USE ereport_mod,            ONLY : ereport
USE errormessagelength_mod, ONLY : errormessagelength
USE yomhook,                ONLY : lhook, dr_hook
USE parkind1,               ONLY : jprb, jpim

IMPLICIT NONE

! IN arguments

REAL, INTENT(IN) ::                                                     &
   ddmfx(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! IN Convective downdraught mass-flux
!                                !    at cloud base (Pa/s)
   cca_steer(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
!                                 ! IN cnv cloud amount cloud base
!                                 ! [0-1]
   u_steer(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
!                                 ! IN U component of wind at cloud base
!                                 ! (m.s-1)
   v_steer(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
!                                 ! IN V component of wind at cloud base
!                                 ! (m.s-1)

INTEGER, INTENT(IN) ::                                                  &
   global_row_length,                                                   &
                                 ! IN number of points on a global row
   global_rows                   ! IN number of global rows

! INOUT arguments

REAL, INTENT (INOUT) ::     ux_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!            x-component of vector sum of c.c.p. front velocities (m/s)
!
                        ,   uy_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!            y-component of vector sum of c.c.p. front velocities (m/s)
!
                        ,   um_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!            scalar sum of c.c.p. front speeds (m/s)
!
                        ,    g_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!            gridbox c.c.p. reduced gravity (m/s^2)
!         (total reduced gravity, was g_new at end of previous timestep)
!
                        ,    h_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!            gridbox c.c.p. height (m)
!         (total height, was h_new at end of previous timestep)
!
                        , riso_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!             remain counter (isotropic) (dimensionless)
!
                        , rdir_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)
!             remain counter (directed) (dimensionless)

! housekeeping

CHARACTER (LEN=errormessagelength) :: cmessage      ! used for ereport
INTEGER                            :: icode         ! used for ereport
INTEGER                            :: n_warn        ! no. of warnings
INTEGER, PARAMETER                 :: max_warn=20
! number of warnings that will trigger a fatal error

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CONV_COLD_POOLS'

! local variables

REAL :: gridbox_x, gridbox_y,                                           &
!       ! gridbox dimensions (m)
        mod_sig_u,                                                      &
!       ! modulus of the vector sum of velocities (m/s)
        direction,                                                      &
!       ! direction that "directed" propagation currently has
        prop_direction,                                                 &
!       ! direction that "directed" propagation may be given
        speed_dir,                                                      &
!       ! speed of directed propagation (m/s)
        speed_np1dir,                                                   &
!       ! projected speed of directed propagation at n+1 (m/s)
        radius_dir,                                                     &
!       ! travel distance of directed propagation inc. "remain" boost(m)
        rad_dt_dir,                                                     &
!       ! travel distance of directed propagation w/o "remain" boost (m)
        g_np1dir,                                                       &
!       ! projected value of directed g' at n+1 (m/s^2)
        h_np1dir,                                                       &
!       ! projected value of directed h at n+1 (m)
        gh_np1dir,                                                      &
!       ! projected value of directed g'h at n+1 (m^2/s^2)
        speed_iso,                                                      &
!       ! speed of isotropic propagation (m/s)
        speed_np1iso,                                                   &
!       ! projected speed of isotropic propagation at n+1 (m/s)
        radius_iso,                                                     &
!       ! travel distance of isotropic propagation inc "remain" boost(m)
        rad_dt_iso,                                                     &
!       ! travel distance of isotropic propagation w/o "remain" boost(m)
        g_np1iso,                                                       &
!       ! projected value of isotropic g' at n+1 (m/s^2)
        h_np1iso,                                                       &
!       ! projected value of isotropic h at n+1 (m)
        gh_np1iso,                                                      &
!       ! projected value of isotropic g'h at n+1 (m^2/s^2)
        ux_unit,                                                        &
!       ! x-component of unit vector
        uy_unit,                                                        &
!       ! y-component of unit vector
        nds,                                                            &
!       ! non-dimensional speed for calculations
        mfa_dd,                                                         &
!       ! downdraught mass flux per unit area at cloud base (kg/m^2/s)
        area_dd
!       ! downdraught area at cloud base (m^2)

! loopers
INTEGER :: i,j,i_ang,iv,jv

! size of halo, used to bound all propagation distances
INTEGER :: halo_size

! ccp local variables
!----------------------

REAL, PARAMETER :: froude = 0.5
! Froude number of gravity-current flow

REAL, PARAMETER :: gamma_l = 2.0
! Proportonality factor from average-current-length model

REAL, PARAMETER :: avg_s = 8.0e3
! representative horizontal gravity-current lengthscale (m)

REAL, PARAMETER :: gh_min = 1.0
! Minimum allowable potential (m^2 s^-2)

REAL, PARAMETER :: cca_min = 1.0e-6
! Minimum convective cloud area to trigger cold pools

REAL, PARAMETER :: speed_min = 1.0e-2
! Minimum equivalent speed to start propagation calculations

REAL, PARAMETER :: spread_dir = 60.0 * pi_over_180
! spread of propagation around the main direction;

REAL, PARAMETER :: rho_a = 1.0
! Representative density for scaling (kg m^-3)

INTEGER, PARAMETER :: n_halo_fields = 5
! number of work fields with haloes to swap at one time

REAL, PARAMETER :: kappa_g = 1.0e4
! g'_DD tuning parameter

REAL, PARAMETER :: kappa_h = 0.5
! h_DD tuning parameter

REAL, PARAMETER :: kappa_b = 1.0e-3
! b_DD tuning parameter

! logicals for edge swapping
LOGICAL, PARAMETER :: l_vector(n_halo_fields) = &
           [ .FALSE.,  .TRUE., .FALSE., .FALSE., .FALSE. ]

LOGICAL, PARAMETER :: l_max_not_add(n_halo_fields) = &
           [ .FALSE., .FALSE., .FALSE.,  .TRUE.,  .TRUE. ]

REAL, PARAMETER :: remain_max = 100.0
! Limit on size of remain counters

INTEGER, PARAMETER :: n_ang = 10
! number of angular subdivisions used to spread propagating currents

LOGICAL :: l_approx_small
! Whether to use the approximation that the propagation distance on
! this timestep is small compared to the total propagation distance.

INTEGER :: i_prop, j_prop, prop_x, prop_y
! x and y propagation in terms of gridbox indices

REAL :: giso   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! isotropic reduced gravity (m/s^2)
      , gdir   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! directed reduced gravity (m/s^2)
      , hiso   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! isotropic cold-pool depth (m)
      , hdir   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! directed cold-pool depth (m)
      , g_dd   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! downdraught reduced gravity scale (m/s^2)
      , h_dd   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! downdraught depth scale (m)
      , b_dd   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! downdraught radial scale (m)
      , pot_dd (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! downdraught cold-pool potential (m^2/s^2)
      , pot_c  (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! collision potential (m^2/s^2)
      , pot_i  (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! isotropic potential (m^2/s^2)
      , phi    (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! fraction of collision potential which is directed
      , delx   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)    &
! gridbox dimension for calculation purposes (m)
      , dely   (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
! gridbox dimension for calculation purposes (m)

!
! quantities with haloes
!

REAL ::  &
  swap_halo_fields(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,&
                   n_halo_fields)

REAL ::  &
        ux_new_a (pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end) &
! x-component of front velocity vector sum at n+1 (m/s)
      , uy_new_a (pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end) &
! y-component of front velocity vector sum at n+1 (m/s)
      , um_new_a (pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end) &
! front speed scalar sum at n+1 (m/s)
      , g_new_a  (pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end) &
! gridbox reduced gravity at n+1 (m/s^2)
      , h_new_a  (pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
! gridbox cold-pool depth at n+1 (m)

REAL ::  &
        ux_new_b (pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end) &
! x-component of front velocity vector sum at n+1 (m/s)
      , uy_new_b (pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end) &
! y-component of front velocity vector sum at n+1 (m/s)
      , um_new_b (pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end) &
! front speed scalar sum at n+1 (m/s)
      , g_new_b  (pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end) &
! gridbox reduced gravity at n+1 (m/s^2)
      , h_new_b  (pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
! gridbox cold-pool depth at n+1 (m)

REAL :: g_new_ab, h_new_ab
! work variables to hold value at a point


LOGICAL :: visited(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
! Array to record whether gridpoints have already
! been visited during directed or isotropic propagation.

REAL :: wind_u, wind_v
! Dummy wind components to be filled in later.

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialisation
!-----------------
icode  = 0
n_warn = 0

l_approx_small = .FALSE.

! dummy values should be zero.
wind_u = 0.0
wind_v = 0.0

! reset arrays to be filled here
!
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    giso(i,j)   = 0.0
    gdir(i,j)   = 0.0
    hiso(i,j)   = 0.0
    hdir(i,j)   = 0.0
    g_dd(i,j)   = 0.0
    h_dd(i,j)   = 0.0
    b_dd(i,j)   = 0.0
    pot_dd(i,j) = 0.0
    pot_c(i,j)  = 0.0
    pot_i(i,j)  = 0.0
    phi(i,j)    = 0.0
  END DO
END DO

DO j = pdims_s%j_start, pdims_s%j_end
  DO i = pdims_s%i_start, pdims_s%i_end
    ux_new_a(i,j) = 0.0
    uy_new_a(i,j) = 0.0
    um_new_a(i,j) = 0.0
    g_new_a(i,j)  = 0.0
    h_new_a(i,j)  = 0.0
    ux_new_b(i,j) = 0.0
    uy_new_b(i,j) = 0.0
    um_new_b(i,j) = 0.0
    g_new_b(i,j)  = 0.0
    h_new_b(i,j)  = 0.0
  END DO
END DO

! set the max. propagation distance in terms of gridboxes
halo_size = MIN(pdims_s%halo_i, pdims_s%halo_j)

!=====================================================

! START initialising and rescaling loop
!----------------------------------------
DO j = pdims%j_start, pdims%j_end

  ! representative y-grid length
  !-------------------------------
  gridbox_y = planet_radius * (xi2_v(j) - xi2_v(j-1))

  DO i = pdims%i_start, pdims%i_end

    ! representative x-grid length
    !-------------------------------
    gridbox_x = planet_radius * (xi1_u(i) - xi1_u(i-1)) * csxi2_p(j)

    delx(i,j) = MAX( gridbox_x, gridbox_y )
    dely(i,j) = delx(i,j)

    ! downdraught area, and mass flux per unit area
    !-----------------------------------------------
    area_dd = cca_steer(i,j) * gridbox_x * gridbox_y
    mfa_dd  = ddmfx(i,j) / g

    ! estimate downdraught-derived g' and h
    !--------------------------------------
    IF (cca_steer(i,j) > cca_min) THEN
      b_dd(i,j) = kappa_b * SQRT(area_dd)
      g_dd(i,j) = kappa_g * (mfa_dd / rho_a)**2.0 / b_dd(i,j)
      h_dd(i,j) = kappa_h * b_dd(i,j)
    END IF

    ! calculate collision potential
    !---------
    pot_c(i,j) = g_ccp(i,j) * h_ccp(i,j)

    ! calculate phi
    !---------
    mod_sig_u = SQRT(ux_ccp(i,j)*ux_ccp(i,j)+uy_ccp(i,j)*uy_ccp(i,j))
    IF ( mod_sig_u > 0.0 ) THEN

      IF ( um_ccp(i,j) <= 0.0 ) THEN ! sanity check
        n_warn = n_warn + 1
        icode=10
        um_ccp(i,j) = EPSILON(1.0)
      END IF

      phi(i,j) = mod_sig_u / um_ccp(i,j)
    ELSE
      phi(i,j) = 0.0
    END IF

    IF ( phi(i,j) > 1.0 + EPSILON(1.0) ) THEN ! sanity check
      n_warn = n_warn + 1
      icode=20
      phi(i,j) = 1.0
    END IF

    ! components of unit vector
    !---------
    IF ( mod_sig_u > 0.0 ) THEN
      ux_unit = ux_ccp(i,j) / mod_sig_u
      uy_unit = uy_ccp(i,j) / mod_sig_u
    ELSE
      ux_unit = 0.0
      uy_unit = 0.0
    END IF

    ! rescale to fraction phi of the Collision Potential
    !---------
    ux_ccp(i,j) = ux_unit * froude * SQRT( phi(i,j) * pot_c(i,j) )
    uy_ccp(i,j) = uy_unit * froude * SQRT( phi(i,j) * pot_c(i,j) )

    ! separate total height into isotropic and directed fractions
    !---------
    hiso(i,j) = ( 1.0 - phi(i,j) ) * h_ccp(i,j)
    hdir(i,j) =         phi(i,j)   * h_ccp(i,j)

    ! initialise isotropic and directed reduced gravity
    !---------
    giso(i,j) = g_ccp(i,j)
    gdir(i,j) = g_ccp(i,j)

    ! Initialise isotropic potential with isotropic collision potl.
    !---------------------------
    pot_i(i,j) = giso(i,j) * hiso(i,j)

  ! END initialising loop
  !-------------------------
  END DO
END DO

!=============================================

! Input from downdraught potential
!-----------------------------------

! START loop over arrays
!-------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    IF ( g_dd(i,j) < 0.0 ) THEN ! sanity check
      n_warn = n_warn + 1
      icode=30
    END IF
    IF ( h_dd(i,j) < 0.0 ) THEN ! sanity check
      n_warn = n_warn + 1
      icode=40
    END IF

    ! initialise downdraught potential
    !---------------------------------
    pot_dd(i,j) = g_dd(i,j) * h_dd(i,j)

    ! update the isotropic properties
    ! using the downdraught properties
    !----------------------------------------
    IF ( pot_dd(i,j) > 0.0 ) THEN
      pot_i(i,j) = pot_i(i,j) + pot_dd(i,j)
      hiso(i,j) = hiso(i,j) + h_dd(i,j)

      IF ( hiso(i,j) <= 0.0 )  THEN ! sanity check
        n_warn = n_warn + 1
        icode=50
        hiso(i,j) = EPSILON(1.0)
      END IF

      giso(i,j) = pot_i(i,j) / hiso(i,j)
    END IF

  ! END  loop
  !-------------
  END DO
END DO

!=============================================

! Directed propagation
!----------------------------

! START "directed" loop over arrays
!-------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    ! initialise the 'visited' array
    !--------------------------------
    DO jv = pdims_s%j_start, pdims_s%j_end
      DO iv = pdims_s%i_start, pdims_s%i_end
        visited(iv,jv) = .FALSE.
      END DO
    END DO

    ! work out the speed, direction,
    ! distance travelled in a timestep
    !------------------------------------------
    speed_dir =   SQRT( ux_ccp(i,j) * ux_ccp(i,j)                    &
                 + uy_ccp(i,j) * uy_ccp(i,j) )

    ! don't work on null points
    !-----------------------------
    IF (speed_dir > speed_min) THEN

      direction = ATAN2( uy_ccp(i,j), ux_ccp(i,j) )
      !
      ! travel distance in one timestep
      !----------------------------------
      rad_dt_dir  = speed_dir * timestep

      ! account (approximately) for being stuck in
      ! the same gridbox over successive timesteps
      !-------------------------------------------
      IF ( rad_dt_dir < delx(i,j) ) THEN
        rdir_ccp(i,j) = MIN( rdir_ccp(i,j) + 1.0, remain_max )
      ELSE
       rdir_ccp(i,j) = 0.0
      END IF
      radius_dir = rad_dt_dir * ( rdir_ccp(i,j) + 1.0 )

      ! which approximation?
      !-----------------------
      IF ( SQRT( phi(i,j) * pot_c(i,j) ) < &
           avg_s / (gamma_l * froude * timestep ) ) THEN
        l_approx_small = .TRUE.
      ELSE
        l_approx_small = .FALSE.
      END IF

      ! new values of g' and h
      !----------------------
      IF (l_approx_small) THEN

        ! base changes on the travel distance for one timestep only!
        !-----------------------------------------------------------
        g_np1dir = gdir(i,j) * (1.0 - gamma_l * rad_dt_dir / avg_s)
        h_np1dir = hdir(i,j) * (1.0 - gamma_l * rad_dt_dir / avg_s)

      ELSE

        ! quasi non-dimensional speed
        !----------------------
        nds = ( ( phi(i,j) * pot_c(i,j) )**0.5) * timestep/delx(i,j)

        IF ( nds <= 0.0 ) THEN ! sanity check
          n_warn = n_warn + 1
          icode=60
          nds = EPSILON(1.0)
        END IF

        g_np1dir = gdir(i,j) *                                       &
       ( (2.0*froude*nds) * (1.0 + 1.0 / (2.0*froude*nds) ) )**(-0.5)
        h_np1dir = hdir(i,j) *                                       &
       ( (2.0*froude*nds) * (1.0 + 1.0 / (2.0*froude*nds) ) )**(-0.5)

      END IF

      ! directed g'h value at n+1
      !----------------------------
      gh_np1dir = g_np1dir * h_np1dir

      IF ( gh_np1dir < 0.0 ) THEN ! sanity check
        n_warn = n_warn + 1
        icode=70
      END IF

      ! work out the velocity components at the new location
      !--------------------------------------------
      speed_np1dir = froude * SQRT( gh_np1dir )

      ! propagate outward
      !---------------------
      DO i_ang = 0, n_ang

        ! i and j 'distances' in terms of array index
        ! with a wind advection component added
        !-----------------------------------------------
        prop_x = NINT( ( radius_dir * COS( direction + spread_dir *  &
                          ( REAL(i_ang) / REAL(n_ang) - 0.5 ) )      &
                        + wind_u * timestep                          &
                      ) / delx(i,j) )
        prop_y = NINT( ( radius_dir * SIN( direction + spread_dir *  &
                          ( REAL(i_ang) / REAL(n_ang) - 0.5 ) )      &
                        + wind_v * timestep                          &
                      ) / dely(i,j) )

        ! limit the absolute distance to the halo size
        !-----------------------------------------------
        prop_x = SIGN(1, prop_x) * MIN( ABS(prop_x), halo_size )
        prop_y = SIGN(1, prop_y) * MIN( ABS(prop_y), halo_size )

        ! indices of arrival point
        !-----------------------------------------------
        i_prop = i + prop_x
        j_prop = j + prop_y

        ! only continue if arrival point is within halo
        !-----------------------------------------------
        IF (       (i_prop >= pdims_s%i_start)                     &
             .AND. (i_prop <= pdims_s%i_end  )                     &
             .AND. (j_prop >= pdims_s%j_start)                     &
             .AND. (j_prop <= pdims_s%j_end  ) ) THEN

          IF ( .NOT. visited( i_prop, j_prop ) ) THEN

            ! If the gridbox isn't exited then the atan
            ! won't make any sense
            !-----------------------------------------------
            IF ((i_prop /= i) .OR. (j_prop /= j)) THEN
              prop_direction = ATAN2( REAL(prop_y), REAL(prop_x) )
            ELSE
              prop_direction = direction
            END IF

            ! add velocity components and modulus into the arrays
            !----------------------------------------------------
            ux_new_a( i_prop, j_prop ) =   ux_new_a( i_prop, j_prop ) &
                                       + speed_np1dir                 &
                                       * COS( prop_direction )
            uy_new_a( i_prop, j_prop ) =   uy_new_a( i_prop, j_prop ) &
                                       + speed_np1dir                 &
                                       * SIN( prop_direction )
            um_new_a( i_prop, j_prop ) =   um_new_a( i_prop, j_prop ) &
                                       + speed_np1dir

            ! update g' and h at the new location
            !--------------------------------------------
            IF ( g_np1dir > g_new_a( i_prop, j_prop ) )            &
             THEN
              g_new_a( i_prop, j_prop ) = g_np1dir
            END IF
            IF ( h_np1dir > h_new_a( i_prop, j_prop ) )            &
             THEN
              h_new_a( i_prop, j_prop ) = h_np1dir
            END IF

            ! mark the point as visited
            !----------------------------
            visited( i_prop, j_prop ) = .TRUE.

          END IF ! not visited
        END IF ! not out-of-bounds
      END DO ! loop over angles
    END IF ! speed_dir > speed_min

! END directed loop over arrays
!-------------------------
  END DO
END DO

!=============================================
! do edges via haloes
!----------------------

DO j = pdims_s%j_start, pdims_s%j_end
  DO i = pdims_s%i_start, pdims_s%i_end
    swap_halo_fields(i,j,1) = ux_new_a(i,j)
    swap_halo_fields(i,j,2) = uy_new_a(i,j)
    swap_halo_fields(i,j,3) = um_new_a(i,j)
    swap_halo_fields(i,j,4) =  g_new_a(i,j)
    swap_halo_fields(i,j,5) =  h_new_a(i,j)
  END DO
END DO

CALL edge_exchange( n_halo_fields,           &
                    l_vector, l_max_not_add, &
                    swap_halo_fields )

DO j = pdims_s%j_start, pdims_s%j_end
  DO i = pdims_s%i_start, pdims_s%i_end
    ux_new_a(i,j) = swap_halo_fields(i,j,1)
    uy_new_a(i,j) = swap_halo_fields(i,j,2)
    um_new_a(i,j) = swap_halo_fields(i,j,3)
    g_new_a(i,j)  = swap_halo_fields(i,j,4)
    h_new_a(i,j)  = swap_halo_fields(i,j,5)
  END DO
END DO

!=============================================

! isotropic propagation
!----------------------------

! START "isotropic" loop over arrays
!-------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    ! initialise the 'visited' array
    !--------------------------------
    DO jv = pdims_s%j_start, pdims_s%j_end
      DO iv = pdims_s%i_start, pdims_s%i_end
        visited(iv,jv) = .FALSE.
      END DO
    END DO

    ! u = Fr (g'h)^1/2
    !-------------------
    speed_iso = froude * SQRT( pot_i(i,j) )

    ! don't work on null points
    !-----------------------------
    IF (speed_iso > speed_min) THEN

      ! physical radius of propagation in a timestep
      !-------------------------------------------------
      radius_iso = timestep * speed_iso

      ! account (approximately) for being stuck in
      ! the same gridbox over successive timesteps
      !-------------------------------------------
      IF ( radius_iso < delx(i,j) ) THEN
        riso_ccp(i,j) = MIN( riso_ccp(i,j) + 1.0, remain_max )
      ELSE
        riso_ccp(i,j) = 0.0
      END IF
      radius_iso = radius_iso * ( riso_ccp(i,j) + 1.0 )

      ! work out new values of isotropic g' and h
      !--------------------------------------------

      ! quasi non-dimensional speed
      !----------------------
      nds = ( ( pot_i(i,j) )**0.5) * timestep/delx(i,j)

      IF ( nds <= 0.0 ) THEN ! sanity check
        n_warn = n_warn + 1
        icode=80
        nds = EPSILON(1.0)
      END IF

      g_np1iso = giso(i,j) *                                         &
       ( (2.0*froude*nds) * (1.0 + 1.0 / (2.0*froude*nds) ) )**(-0.5)
      h_np1iso = hiso(i,j) *                                         &
       ( (2.0*froude*nds) * (1.0 + 1.0 / (2.0*froude*nds) ) )**(-0.5)

      ! isotropic g'h value at n+1
      !---------
      gh_np1iso = g_np1iso * h_np1iso

      IF ( gh_np1iso < 0.0 ) THEN ! sanity check
        n_warn = n_warn + 1
        icode=90
      END IF

      ! work out the radial speed at the new location
      !--------------------------------------------
      speed_np1iso = froude * SQRT( gh_np1iso )

      ! propagate outward
      !---------------------
      DO i_ang = 0, n_ang - 1

        ! i and j 'distances' in terms of array index
        ! with a wind advection component added
        !-----------------------------------------------
        prop_x = NINT( ( radius_iso *                                &
                  COS( 2.0 * pi * REAL(i_ang) / REAL(n_ang) )        &
                        + wind_u * timestep                          &
                      ) / delx(i,j) )
        prop_y = NINT( ( radius_iso *                                &
                  SIN( 2.0 * pi * REAL(i_ang) / REAL(n_ang) )        &
                        + wind_v * timestep                          &
                      ) / dely(i,j) )

        ! limit the absolute distance to the halo size
        !-----------------------------------------------
        prop_x = SIGN(1, prop_x) * MIN( ABS(prop_x), halo_size )
        prop_y = SIGN(1, prop_y) * MIN( ABS(prop_y), halo_size )

        ! indices of arrival point
        !-----------------------------------------------
        i_prop = i + prop_x
        j_prop = j + prop_y

        ! only continue if arrival point is within halo
        !-----------------------------------------------
        IF (       (i_prop >= pdims_s%i_start)                     &
             .AND. (i_prop <= pdims_s%i_end  )                     &
             .AND. (j_prop >= pdims_s%j_start)                     &
             .AND. (j_prop <= pdims_s%j_end  ) ) THEN

          IF ( .NOT. visited( i_prop, j_prop ) ) THEN

            ! add velocity components and modulus into the arrays
            !----------------------------------------------------
            !
            ! For ISOTROPIC propagation, only update the directional
            ! components if the gridbox has been exited.
            !-------------------------------------------------------
            IF ((i_prop /= i) .OR. (j_prop /= j)) THEN
              ux_new_b( i_prop, j_prop ) =   ux_new_b( i_prop, j_prop )  &
               + speed_np1iso *                                          &
                  COS( ATAN2( REAL(prop_y), REAL(prop_x) ) )
              uy_new_b( i_prop, j_prop ) =   uy_new_b( i_prop, j_prop )  &
               + speed_np1iso *                                          &
                  SIN( ATAN2( REAL(prop_y), REAL(prop_x) ) )
            END IF
            um_new_b( i_prop, j_prop ) =   um_new_b( i_prop, j_prop ) &
                                       + speed_np1iso

            ! update g' and h at the new location
            !--------------------------------------------
            IF ( g_np1iso > g_new_b( i_prop, j_prop ) )             &
             THEN
              g_new_b( i_prop, j_prop ) = g_np1iso
            END IF
            IF ( h_np1iso > h_new_b( i_prop, j_prop ) )             &
             THEN
              h_new_b( i_prop, j_prop ) = h_np1iso
            END IF

            ! mark the point as visited
            !----------------------------
            visited( i_prop, j_prop ) = .TRUE.

          END IF ! not visited
        END IF ! not out-of-bounds
      END DO ! loop over angles
    END IF ! speed_iso > speed_min

! END isotropic loop over arrays
!-------------------------
  END DO
END DO

!=============================================
! do edges via haloes
!----------------------

DO j = pdims_s%j_start, pdims_s%j_end
  DO i = pdims_s%i_start, pdims_s%i_end
    swap_halo_fields(i,j,1) = ux_new_b(i,j)
    swap_halo_fields(i,j,2) = uy_new_b(i,j)
    swap_halo_fields(i,j,3) = um_new_b(i,j)
    swap_halo_fields(i,j,4) =  g_new_b(i,j)
    swap_halo_fields(i,j,5) =  h_new_b(i,j)
  END DO
END DO

CALL edge_exchange( n_halo_fields,           &
                    l_vector, l_max_not_add, &
                    swap_halo_fields )

DO j = pdims_s%j_start, pdims_s%j_end
  DO i = pdims_s%i_start, pdims_s%i_end
    ux_new_b(i,j) = swap_halo_fields(i,j,1)
    uy_new_b(i,j) = swap_halo_fields(i,j,2)
    um_new_b(i,j) = swap_halo_fields(i,j,3)
    g_new_b(i,j)  = swap_halo_fields(i,j,4)
    h_new_b(i,j)  = swap_halo_fields(i,j,5)
  END DO
END DO

!=============================================

! Check the number of anomalies
! and issue a fatal error or a warning
! depending on the number.
!---------------------------------------
IF (icode /= 0) THEN
  IF (n_warn > max_warn) THEN
    WRITE(cmessage,FMT='(A)')                    &
     "Anomalous behaviour of code variables, hence stopping."
    CALL ereport(RoutineName, icode, cmessage)
  ELSE
    icode = -icode ! error code sign change
    WRITE(cmessage,FMT='(A)')                    &
     "Some anomalies in code variables have been logged."
    CALL ereport(RoutineName, icode, cmessage)
  END IF
END IF

! update prognostic arrays
! but remove very small values
!----------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    g_new_ab =   MAX( g_new_a(i,j), g_new_b(i,j) )
    h_new_ab =   MAX( h_new_a(i,j), h_new_b(i,j) )

    IF (g_new_ab * h_new_ab >= gh_min) THEN

      ux_ccp(i,j) = ux_new_a(i,j) + ux_new_b(i,j)
      uy_ccp(i,j) = uy_new_a(i,j) + uy_new_b(i,j)
      um_ccp(i,j) = um_new_a(i,j) + um_new_b(i,j)
      g_ccp( i,j) = g_new_ab
      h_ccp( i,j) = h_new_ab

    ELSE

      ux_ccp(i,j) = 0.0
      uy_ccp(i,j) = 0.0
      um_ccp(i,j) = 0.0
      g_ccp( i,j) = 0.0
      h_ccp( i,j) = 0.0

      ! zero remain counters at unoccupied positions
      !-----------------------------------------------
      riso_ccp(i,j) = 0.0
      rdir_ccp(i,j) = 0.0

    END IF

  END DO
END DO

!=====================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE conv_cold_pools
END MODULE conv_cold_pools_mod

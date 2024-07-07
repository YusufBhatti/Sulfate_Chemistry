! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Gravity Wave (Ultra-Simple Spectral Parametrization) Scheme.
! Subroutine Interface:

MODULE gw_ussp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GW_USSP_MOD'

CONTAINS

SUBROUTINE gw_ussp(levels, rows, nrows,                                 &
  off_x, off_y, halo_i, halo_j, row_length,                             &
  global_row_length,n_proc, n_procy, proc_row_group,at_extremity,       &
  r_rho_levels, r_theta_levels, p_layer_boundaries,                     &
  sin_theta_longitude, sin_theta_latitude,                              &
  theta, rho, u, v, totalppn, timestep,                                 &
  r_u, r_v, T_inc,                                                      &
  L_ussp_heating,                                                       &
  gwspec_eflux,gwspec_sflux,gwspec_wflux,gwspec_nflux,                  &
  gwspec_ewacc,gwspec_nsacc,                                            &
  gwspec_eflux_on, gwspec_eflux_p_on, gwspec_sflux_on,                  &
  gwspec_wflux_on, gwspec_wflux_p_on, gwspec_nflux_on,                  &
  gwspec_ewacc_on, gwspec_ewacc_p_on, gwspec_nsacc_on)
!
! purpose:    This subroutine calculates the vertical momentum flux
!             divergence due to gravity waves as parametrised by the
!             Warner and McIntyre Ultra Simple Spectral gravity wave
!             Parametrization adapted for use in the UM.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Gravity Wave Drag
!
! code description:
!   language: fortran 90
!   this code is written to umdp3 programming standards.
!   documentation: Unified Model Documentation Paper 34 (Non-Orog GW)

! Processor addresses (N, E, S, W)
USE um_parparams, ONLY: pnorth, peast, psouth, pwest

USE g_wave_input_mod, ONLY: ussp_launch_factor,wavelstar, L_add_cgw,    &
                            cgw_scale_factor
USE planet_constants_mod, ONLY: two_omega, g, r, cp,                    &
                                planet_radius, recip_a2
USE gw_ussp_params_mod, ONLY: idir

USE tuning_segments_mod, ONLY: ussp_seg_size

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                        &
   udims, vdims, pdims, tdims, pdims,                                   &
   udims_s, vdims_s, pdims_s, tdims_s, pdims_l, tdims_l


! Model level values
USE level_heights_mod, ONLY:                                            &
   eta_theta_levels           ! Eta values of theta levels
USE horiz_grid_mod, ONLY :                                              &
   intw_rho2w                 ! Weight for vertical interpolation
USE vertnamelist_mod, ONLY:                                             &
! Top of model height (metres)
! Lowest level where rho level has constant radius (i.e. not terrain-following)
   z_top_of_model, first_constant_r_rho_level

USE missing_data_mod, ONLY: rmdi
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Field_Types
!$    USE omp_lib
USE model_domain_mod, ONLY: model_type, mt_global
USE segments_mod, ONLY:                                                 &
  segment_type, meta_segment_type,                                      &
  segments_mod_seg_meta, segments_mod_segments
USE umPrintMgr, ONLY: umPrint,ummessage


! ----------------------------------------------------------------------+-------
! Subroutines defined from modules
! ----------------------------------------------------------------------+-------
USE gw_ussp_core_mod, ONLY: gw_ussp_core

USE eg_v_at_poles_mod, ONLY: eg_v_at_poles

USE p_to_t_mod, ONLY: p_to_t
USE p_to_u_mod, ONLY: p_to_u
USE p_to_v_mod, ONLY: p_to_v
USE polar_row_mean_mod, ONLY: polar_row_mean
USE u_to_p_mod, ONLY: u_to_p
USE v_to_p_mod, ONLY: v_to_p
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

! ----------------------------------------------------------------------+-------

IMPLICIT NONE

! ----------------------------------------------------------------------+-------
! Subroutine arguments of GW_USSP
! ----------------------------------------------------------------------+-------
!     Fixed starting theta-level (e.g. for P_layer_boundaries)
INTEGER, PARAMETER :: tkfix0start    = 0
!
!     Fixed starting theta-level (e.g. for STASH diagnostic arrays)
INTEGER, PARAMETER :: tkfix1start    = 1
!
! Number of model levels
INTEGER, INTENT(IN) :: levels
!
! Number of rows for u field
INTEGER, INTENT(IN) :: rows
!
! Number of rows for v field
INTEGER, INTENT(IN) :: nrows
!
! Offset longitude, latitude
INTEGER, INTENT(IN) :: off_x, off_y
!
! Halo in longitude, latitude
INTEGER, INTENT(IN) :: halo_i, halo_j
!
! Number of points on a row of a full global field
INTEGER, INTENT(IN) :: global_row_length
!
! Group id for processors on the same row
INTEGER, INTENT(IN) :: proc_row_group
!
! Total number of processors
INTEGER, INTENT(IN) :: n_proc
!
! Number of processors in latitude
INTEGER, INTENT(IN) :: n_procy
!
! Number of grid points in row
INTEGER, INTENT(IN) :: row_length
! ----------------------------------------------------------------------+-------
!
! Timestep
REAL, INTENT(IN)    :: timestep
!
! Grid point longitudes
REAL, INTENT(IN)    :: sin_theta_longitude(tdims%i_start:tdims%i_end,   &
                                           tdims%j_start:tdims%j_end)
!
! P-GRID Latitudes
REAL, INTENT(IN)    :: sin_theta_latitude(tdims%i_start:tdims%i_end,    &
                                          tdims%j_start:tdims%j_end)
!
! Centre of earth distance to theta-level points
REAL, INTENT(IN)    :: r_theta_levels(tdims_l%i_start:tdims_l%i_end,    &
                                      tdims_l%j_start:tdims_l%j_end,    &
                                                    0:tdims_l%k_end)
!
! Centre of earth distance to rho-level points
REAL, INTENT(IN)    :: r_rho_levels(pdims_l%i_start:pdims_l%i_end,      &
                                    pdims_l%j_start:pdims_l%j_end,      &
                                                    pdims_l%k_end)
!
! Pressure on layer boundaries
REAL, INTENT(IN)    :: p_layer_boundaries(tdims%i_start:tdims%i_end,    &
                                          tdims%j_start:tdims%j_end,    &
                                            tkfix0start:tdims%k_end)
!
! Total Precipitation for CGW source
REAL, INTENT(IN)    :: totalppn(tdims%i_start:tdims%i_end,              &
                                     tdims%j_start:tdims%j_end)
!
! Primary model array for theta
REAL, INTENT(IN)    :: theta(tdims%i_start:tdims%i_end,                 &
                             tdims%j_start:tdims%j_end,                 &
                               tkfix1start:tdims%k_end)
!
! Primary model array for Density x (radius earth)^2.
REAL, INTENT(IN)    :: rho(pdims_s%i_start:pdims_s%i_end,               &
                           pdims_s%j_start:pdims_s%j_end,               &
                           pdims_s%k_start:pdims_s%k_end)
!
! Primary model array for U-wind field
REAL, INTENT(IN)    :: u(udims_s%i_start:udims_s%i_end,                 &
                         udims_s%j_start:udims_s%j_end,                 &
                         udims_s%k_start:udims_s%k_end)
!
! Primary model array for V-wind field
REAL, INTENT(IN)    :: v(vdims_s%i_start:vdims_s%i_end,                 &
                         vdims_s%j_start:vdims_s%j_end,                 &
                         vdims_s%k_start:vdims_s%k_end)
!
! FLux Fp in E azimuth
REAL, INTENT(INOUT) :: gwspec_eflux(udims%i_start:udims%i_end,          &
                                    udims%j_start:udims%j_end,          &
                                      tkfix1start:tdims%k_end)
!
! FLux Fp in S azimuth
REAL, INTENT(INOUT) :: gwspec_sflux(vdims%i_start:vdims%i_end,          &
                                    vdims%j_start:vdims%j_end,          &
                                      tkfix1start:tdims%k_end)
!
! FLux Fp in W azimuth
REAL, INTENT(INOUT) :: gwspec_wflux(udims%i_start:udims%i_end,          &
                                    udims%j_start:udims%j_end,          &
                                      tkfix1start:tdims%k_end)
!
! Fp in N azimuth
REAL, INTENT(INOUT) :: gwspec_nflux(vdims%i_start:vdims%i_end,          &
                                    vdims%j_start:vdims%j_end,          &
                                      tkfix1start:tdims%k_end)
!
! Accel of U wind
REAL, INTENT(INOUT) :: gwspec_ewacc(udims%i_start:udims%i_end,          &
                                    udims%j_start:udims%j_end,          &
                                    udims%k_start:udims%k_end)
!
! Accel of V wind
REAL, INTENT(INOUT) :: gwspec_nsacc(vdims%i_start:vdims%i_end,          &
                                    vdims%j_start:vdims%j_end,          &
                                    vdims%k_start:vdims%k_end)
!
! Temperature increment
REAL, INTENT(INOUT) :: T_inc(tdims%i_start:tdims%i_end,                 &
                             tdims%j_start:tdims%j_end,                 &
                               tkfix1start:tdims%k_end)
!
! U-wind increment diagnostic
REAL, INTENT(INOUT) :: r_u(udims_s%i_start:udims_s%i_end,               &
                           udims_s%j_start:udims_s%j_end,               &
                           udims_s%k_start:udims_s%k_end)
!
! V-wind increment diagnostic
REAL, INTENT(INOUT) :: r_v(vdims_s%i_start:vdims_s%i_end,               &
                           vdims_s%j_start:vdims_s%j_end,               &
                           vdims_s%k_start:vdims_s%k_end)
! ----------------------------------------------------------------------+-------
!
! Switch to calculate heating tendency
LOGICAL, INTENT(IN) :: L_ussp_heating
!
! Switch for diagnostic of Fp in azimuthal direction East
LOGICAL, INTENT(IN) :: gwspec_eflux_on
!
! Switch for diagnostic of Fp in azimuthal direction East interpolated to p-levs
LOGICAL, INTENT(IN) :: gwspec_eflux_p_on
!
! Switch for diagnostic of Fp in azimuthal direction South
LOGICAL, INTENT(IN) :: gwspec_sflux_on
!
! Switch for diagnostic of Fp in azimuthal direction West
LOGICAL, INTENT(IN) :: gwspec_wflux_on
!
! Switch for diagnostic of Fp in azimuthal direction West interpolated to p-levs
LOGICAL, INTENT(IN) :: gwspec_wflux_p_on
!
! Switch for diagnostic of Fp in azimuthal direction North
LOGICAL, INTENT(IN) :: gwspec_nflux_on
!
!Switch for net Eastward acceleration diagnostic
LOGICAL, INTENT(IN) :: gwspec_ewacc_on
!
!Switch for net Eastward acceleration diagnostic interpolated to p-levels
LOGICAL, INTENT(IN) :: gwspec_ewacc_p_on
!
!Switch for net Northward acceleration diagnostic
LOGICAL, INTENT(IN) :: gwspec_nsacc_on
!
! Edge of (global) domain indicator
LOGICAL, INTENT(IN) :: at_extremity(4)
!
! ----------------------------------------------------------------------+-------
! Local parameters
! ----------------------------------------------------------------------+-------
!
INTEGER, PARAMETER :: one = 1     ! Hard coded indexing
!
! Eta (hybrid height) level for model launch
REAL, PARAMETER :: etalaunch         = 0.045
!
! Constant radius height level for model launch (m)
REAL, PARAMETER :: zlaunch           = 3825.0
!
! ----------------------------------------------------------------------+-------
! Security parameters
! ----------------------------------------------------------------------+-------
!
! Minimum allowed value of buoyancy frequency squared
REAL, PARAMETER ::  sqnmin           = 1.0e-4
!
! ----------------------------------------------------------------------+-------
! Local Constants (a) Physical
! ----------------------------------------------------------------------+-------
!
! Cos(phi_j) - azimuthal direction
REAL            :: cosphi(4)
!
! Sin(phi_j) - azimuthal direction
REAL            :: sinphi(4)
!
! Reciprocal of middle atmosphere mean scale height for pressure,
! conventionally assumed to be around 7km: 1 / H = g / (R * 239.145)
REAL :: rscale_h
!
! ----------------------------------------------------------------------+-------
! Local variables (scalars) used in GW_USSP
! Some effectively expanded to workspace (using vector registers)
! ----------------------------------------------------------------------+-------
!
! Longitude index
INTEGER         :: i
!
! Latitude index
INTEGER         :: j
!
! Level index
INTEGER         :: k
!
! Azimuthal direction index
INTEGER         :: jdir
!
! omp block iterator
INTEGER         :: jj
!
! blocking size for omp_block
INTEGER         :: omp_block
!
! Level for gw launch - in principle this could be launchlev(i,j)
INTEGER         :: launchlev
!
! Minimum level in launchlev(i,j)
INTEGER         :: minlaunchlev
!
! Maximum level in launchlev(i,j)
INTEGER         :: maxlaunchlev
INTEGER         :: ii             ! Looper
INTEGER         :: ss             ! Addressing in segment
INTEGER         :: num_parallel   ! Segmentation variable
INTEGER         :: ipar           ! Segmentation variable
INTEGER         :: num_segments   ! Number of segments
!
! Wave-induced force per unit mass due to azimuthal sectors (m s^-2)
REAL            :: g_g
!
! Layer depths for heating calculations
REAL            :: dzb, dzu, dzl
!
! Velocities on theta levels
REAL            :: uhat, vhat
!
! Kinetic energy tendencies on theta levels
REAL            :: ududt, vdvdt
!
! ----------------------------------------------------------------------+-------
! Local variables (dynamic arrays) used in GW_USSP
! ----------------------------------------------------------------------+-------
!
! Constant radius level heights
REAL            :: z_levels(tkfix1start:tdims%k_end)
!
! Pseudomomentum flux integrated over azimuthal sector (kg m^-1 s^-2)
REAL            :: fptot(tdims%i_start:tdims%i_end,                     &
                 tdims%j_start:tdims%j_end,tkfix1start:tdims%k_end,idir)
!
! Zonal component of wave-induced force on density level (m s^-2)
! Note that in our notation G_X equates to DU_DT, the zonal wind tendency
REAL            :: g_x(pdims%i_start:pdims%i_end,                       &
                   pdims%j_start:pdims%j_end,pdims%k_start:pdims%k_end)
!
! Meridional component of wave-induced force on density level (m s^-2)
! Note that in our notation G_Y equates to DV_DT, meridional wind tendency
REAL            :: g_y(pdims%i_start:pdims%i_end,                       &
                   pdims%j_start:pdims%j_end,pdims%k_start:pdims%k_end)
!
! Buoyancy [Brunt Vaisala] frequency on half-levels (rad s^-1)
REAL            :: nbv(tdims%i_start:tdims%i_end,                       &
                      tdims%j_start:tdims%j_end,tkfix1start:tdims%k_end)
!
! Rho on theta levels
REAL            :: rho_th(tdims%i_start:tdims%i_end,                    &
                      tdims%j_start:tdims%j_end,tkfix1start:tdims%k_end)
!
! Component of wind in phi_jdir azimuth direction (m s^-1)
REAL            :: udotk(tdims%i_start:tdims%i_end,                     &
                 tdims%j_start:tdims%j_end,tkfix1start:tdims%k_end,idir)
!
! Zonal wind U interpolated onto Rho grid
REAL            :: uonp(pdims%i_start:pdims%i_end,                      &
                   pdims%j_start:pdims%j_end,pdims%k_start:pdims%k_end)
!
! Meridional wind V interpolated onto Rho grid
REAL            :: vonp(pdims%i_start:pdims%i_end,                      &
                   pdims%j_start:pdims%j_end,pdims%k_start:pdims%k_end)
!
! G_X with halo
REAL            :: g_xp_smallhalo(pdims_s%i_start:pdims_s%i_end,        &
                                  pdims_s%j_start:pdims_s%j_end,        &
                                  pdims_s%k_start:pdims_s%k_end)
!
! G_Y with halo
REAL            :: g_yp_smallhalo(pdims_s%i_start:pdims_s%i_end,        &
                                  pdims_s%j_start:pdims_s%j_end,        &
                                  pdims_s%k_start:pdims_s%k_end)
!
! Rho*(radius_earth)^2 on Theta grid with halo
REAL            :: rhont_smallhalo(tdims_s%i_start:tdims_s%i_end,       &
                                   tdims_s%j_start:tdims_s%j_end,       &
                                       tkfix1start:tdims_s%k_end)
!
! Zonal wind Increment dU / dt on rho grid
REAL            :: uincrt(udims%i_start:udims%i_end,                    &
                   udims%j_start:udims%j_end,udims%k_start:udims%k_end)
!
! Meridional wind Increment  dV / dt on rho grid
REAL            :: vincrt(vdims%i_start:vdims%i_end,                    &
                   vdims%j_start:vdims%j_end,vdims%k_start:vdims%k_end)
!
! Workspace arrays for setting up interpolations
REAL            :: work_TsmallHALO(tdims_s%i_start:tdims_s%i_end,       &
                                   tdims_s%j_start:tdims_s%j_end,       &
                                       tkfix1start:tdims_s%k_end)
!
REAL, ALLOCATABLE :: s_sin_theta_lat(:)
REAL, ALLOCATABLE :: s_totalppn(:)
REAL, ALLOCATABLE :: s_rho_th(:,:)
REAL, ALLOCATABLE :: s_nbv(:,:)
REAL, ALLOCATABLE :: s_udotk(:,:,:)
REAL, ALLOCATABLE :: s_fptot(:,:,:)

! Variables for segmentation
TYPE(segment_type),      ALLOCATABLE  :: segments(:)
TYPE(meta_segment_type)               :: meta_segments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GW_USSP'


!
!  End of Header
!
! ==Main Block==--------------------------------------------------------+-------
DATA cosphi/0,-1,0,1/
DATA sinphi/1,0,-1,0/
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------
! Local Constants (a) Physical
! ----------------------------
rscale_h = g / (r * 239.145)
!
! ----------------------------------------------------------------------+-------
! Find model level for launch
! ----------------------------------------------------------------------+-------
! Geographically fixed value - rethink if variable height sources are used
! WARNING: Not defined relative to tdims%k_start as endgame alters it
launchlev=tkfix1start
DO k=(tkfix1start + 1),tdims%k_end
  IF (eta_theta_levels(k) >= etalaunch .AND.                            &
      eta_theta_levels(k-1) <  etalaunch) THEN
    IF ((eta_theta_levels(k)-etalaunch) <                               &
        (etalaunch-eta_theta_levels(k-1))) THEN
      launchlev=k
    ELSE
      launchlev=k-1
    END IF
  END IF
END DO
!
! ----------------------------------------------------------------------+-------
!     Interpolate : [Vertical]   RHO onto T grid
!             and : [Horizontal] U,V onto P grid
! ----------------------------------------------------------------------+-------
!$OMP PARALLEL DEFAULT(SHARED)                                          &
!$OMP&         PRIVATE( i, j, k, jdir)
CALL p_to_t (row_length, rows, halo_i, halo_j,                          &
             off_x,off_y,levels-1,r_theta_levels,                       &
             r_rho_levels,rho,rhont_smallhalo)
! Need to barrier here to make sure previous routine has fully completed
! before we use the results
!$OMP BARRIER
! Extrapolate topmost level to be a scale height from level below
!$OMP DO SCHEDULE(STATIC)
Rows_do_init: DO j=tdims%j_start,tdims%j_end
  Row_length_do_init: DO i=tdims%i_start,tdims%i_end
    rhont_smallhalo(i,j,tdims%k_end) =                                  &
      rhont_smallhalo(i,j,tdims%k_end-1) *                              &
      EXP(-(r_theta_levels(i,j,tdims%k_end)-                            &
            r_theta_levels(i,j,tdims%k_end-1)) * rscale_h)
  END DO  Row_length_do_init
END DO  Rows_do_init
!$OMP END DO NOWAIT
!
CALL u_to_p(u,                                                    &
                  udims_s%i_start,udims_s%i_end,                  &
                  udims_s%j_start,udims_s%j_end,                  &
                  pdims%i_start,pdims%i_end,                      &
                  pdims%j_start,pdims%j_end,                      &
                  levels, at_extremity,uonp )

CALL v_to_p(v,                                                    &
                  vdims_s%i_start,vdims_s%i_end,                  &
                  vdims_s%j_start,vdims_s%j_end,                  &
                  pdims%i_start,pdims%i_end,                      &
                  pdims%j_start,pdims%j_end,                      &
                  levels, at_extremity,vonp )

! Need to make sure all parallelised calculations in these routines are
! complete before we use the computed values
!$OMP BARRIER

!
! ----------------------------------------------------------------------+-------
!     Initialize local arrays to zero
! ----------------------------------------------------------------------+-------
! ---------------------------------
!     Set winds with haloes to zero
! ---------------------------------

!$OMP DO SCHEDULE(STATIC)
Levels_do1: DO k=pdims_s%k_start,pdims_s%k_end
  Rows_do1: DO j=pdims_s%j_start,pdims_s%j_end
    Row_length_do1: DO i=pdims_s%i_start,pdims_s%i_end
      g_xp_smallhalo(i,j,k) = 0.0
      g_yp_smallhalo(i,j,k) = 0.0
    END DO  Row_length_do1
  END DO  Rows_do1
END DO  Levels_do1
!$OMP END DO NOWAIT
!
!$OMP DO SCHEDULE(STATIC)
Levels_do1a: DO k=pdims%k_start,pdims%k_end
  Rows_do1a: DO j=pdims%j_start,pdims%j_end
    Row_length_do1a: DO i=pdims%i_start,pdims%i_end
! ----------------------------------------------------------------------+-------
!     Zero vertical divergence of pseudomomentum flux
! ----------------------------------------------------------------------+-------
      g_x(i,j,k) = 0.0
      g_y(i,j,k) = 0.0
    END DO  Row_length_do1a
  END DO  Rows_do1a
END DO  Levels_do1a
!$OMP END DO NOWAIT

!
! ----------------------------------------------------------------------+-------
! 1.0   Set variables that are to be defined on all model levels
! ----------------------------------------------------------------------+-------
!

!$OMP DO SCHEDULE(STATIC)
Levels_do2: DO k=(tkfix1start + 1),(tdims%k_end - 1)
! ----------------------------------------------------------------------+-------
! 1.1 Density, buoyancy frequency and altitude for middle levels
! ----------------------------------------------------------------------+-------
  Rows_do2: DO j=tdims%j_start,tdims%j_end
    Row_length_do2: DO i=tdims%i_start,tdims%i_end
      rho_th(i,j,k) = rhont_smallhalo(i,j,k) * recip_a2
!     Buoyancy (Brunt-Vaisala) frequency calculation
      nbv(i,j,k) = ( g*(theta(i,j,k+1) - theta(i,j,k-1)) /              &
                          (theta(i,j,k) *                               &
        (r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k-1)) ) )
      nbv(i,j,k) = MAX(nbv(i,j,k), sqnmin)
      nbv(i,j,k) = SQRT(nbv(i,j,k))
    END DO  Row_length_do2
  END DO  Rows_do2
END DO  Levels_do2
!$OMP END DO

!
! ----------------------------------------------------------------------+-------
!     Density, Buoyancy (Brunt-Vaisala) frequency at top and bottom
! ----------------------------------------------------------------------+-------

!$OMP DO SCHEDULE(STATIC)
Rows_do3: DO j=tdims%j_start,tdims%j_end
  Row_length_do3: DO i=tdims%i_start,tdims%i_end
!   Set rho and theta level value of Z_TH.
    rho_th(i,j,tkfix1start) = rhont_smallhalo(i,j,tkfix1start) *        &
                              recip_a2
    nbv(i,j,tkfix1start)    = nbv(i,j,(tkfix1start + 1))
    rho_th(i,j,tdims%k_end) = rhont_smallhalo(i,j,tdims%k_end) *        &
                              recip_a2
    nbv(i,j,tdims%k_end)    = nbv(i,j,tdims%k_end-1)
  END DO  Row_length_do3
END DO  Rows_do3
!$OMP END DO

!
Levels_do4: DO k=tdims%k_end,(tkfix1start + 1),-1
! ----------------------------------------------------------------------+-------
! 1.2 Set buoyancy frequency constant up to 1km altitude
! ----------------------------------------------------------------------+-------
!$OMP DO SCHEDULE(STATIC)
  Rows_do4: DO j=tdims%j_start,tdims%j_end
    Row_length_do4: DO i=tdims%i_start,tdims%i_end
      IF ( (r_theta_levels(i,j,k) - planet_radius) <  1.0e3)            &
        nbv(i,j,k-1) = nbv(i,j,k)
    END DO  Row_length_do4
  END DO  Rows_do4
!$OMP END DO
END DO  Levels_do4

!
! ----------------------------------------------------------------------+-------
! 1.3 Compute component of wind U in each wave-propagation direction.
!     U is the dot product of (u,v) with k_0 but n.b. UDOTK is half-way
!     between rho levels.
!     Interpolate : [Vertical]   RHO onto T grid
! ----------------------------------------------------------------------+-------

!$OMP DO SCHEDULE(STATIC)
IDir_do1: DO jdir=1,idir
!
  Levels_do5: DO k=tkfix1start,(tdims%k_end - 1)
    Rows_do5: DO j=tdims%j_start,tdims%j_end
      Row_length_do5: DO i=tdims%i_start,tdims%i_end
        udotk(i,j,k,jdir) =                                             &
        0.5*(uonp(i,j,k) + uonp(i,j,k+1))*cosphi(jdir) +                &
        0.5*(vonp(i,j,k) + vonp(i,j,k+1))*sinphi(jdir)
!       Assumes theta levels are half way between rho levels.
!       Below is preferred except that what worked at 8.5 seems not to at 10.1!
!        ( (intw_rho2w(k,2)*uonp(i,j,k)) +                              &
!          (intw_rho2w(k,1)*uonp(i,j,k+1)) )*cosphi(jdir) +             &
!        ( (intw_rho2w(k,2)*vonp(i,j,k)) +                              &
!          (intw_rho2w(k,1)*vonp(i,j,k+1)) )*sinphi(jdir)
      END DO  Row_length_do5
    END DO  Rows_do5
  END DO  Levels_do5
!
! ----------------------------------------------------------------------+-------
! Set wind component for top level, to be equal to that on the top Rho level
! ----------------------------------------------------------------------+-------
  Rows_do5a: DO j=tdims%j_start,tdims%j_end
    Row_length_do5a: DO i=tdims%i_start,tdims%i_end
      udotk(i,j,tdims%k_end,jdir) = uonp(i,j,pdims%k_end)*cosphi(jdir)  &
                                  + vonp(i,j,pdims%k_end)*sinphi(jdir)
    END DO  Row_length_do5a
  END DO  Rows_do5a
!
END DO  IDir_do1
!$OMP END DO

!$OMP END PARALLEL

!
! ----------------------------------------------------------------------+-------
! Core calculations on Physics Grid now within separate subroutine
! ----------------------------------------------------------------------+-------
!
! As presently coded, a single launch level is set but in principle it could be
! variable. Optimization is easier if we define level limits such that:
! 1 <= k < minlaunchlev            : initial total flux = 0.
! minlaunchlev <= k < maxlaunchlev : initial total flux = 0. / launch flux(i,j)
! maxlaunchlev <= k                : initial total flux = launch flux(i,j)
! For now they are simply assigned (also replicated in the core code)
minlaunchlev = launchlev
maxlaunchlev = launchlev
! ----------------------------------------------------------------------+-------
! 3.0 Initialize gravity wave spectrum variables
! ----------------------------------------------------------------------+-------
! 4.0 Calculations carried out for each azimuthal direction
! ----------------------------------------------------------------------+-------
!
! L_add_cgw = F : calculate standard USSP isotropic GW launch flux
! L_add_cgw = T : calculate variable CGW launch flux

! Number of threads is 1 at this point: we are not load-balancing between
! threads, then segmenting here.  Rather, set up segmentation on one thread,
! then distribute the segments among threads.

! Variables needed for segmentation
num_parallel = 1
ipar         = 1
num_segments = -99
!Set up segment meta-information.
CALL segments_mod_seg_meta(meta_segments, ipar, num_parallel,      &
                           rows*row_length, ussp_seg_size, num_segments)
 
! Allocate space for segmentation arrays
ALLOCATE(  segments( meta_segments%num_segments ) )

! Work out starting points and sizes of each segment individually
CALL segments_mod_segments(segments, meta_segments, row_length, rows)

! Main parallelisation and segmentation loop
! For this we allocate space for each segment and then copy necessary data
! in to these variables before calling the ussp_core with them. The result
! is then copied back out. Looping is a little awkward as segments may well
! span partial rows.
!$OMP PARALLEL DEFAULT(SHARED)                                          &
!$OMP PRIVATE(ii, jj, ss, i, k, jdir, s_sin_theta_lat, s_totalppn,      &
!$OMP         s_rho_th, s_nbv, s_udotk, s_fptot)


!$OMP DO SCHEDULE(DYNAMIC)
DO i = 1, meta_segments%num_segments

  ALLOCATE(s_sin_theta_lat(segments(i)%seg_points))
  ALLOCATE(s_totalppn(segments(i)%seg_points))
  ALLOCATE(s_rho_th(segments(i)%seg_points,tkfix1start:tdims%k_end))
  ALLOCATE(s_nbv(segments(i)%seg_points,tkfix1start:tdims%k_end))
  ALLOCATE(s_udotk(segments(i)%seg_points, tkfix1start:tdims%k_end, idir))
  ALLOCATE(s_fptot(segments(i)%seg_points, tkfix1start:tdims%k_end, idir))
 
  ii = segments(i)%first_x
  jj = segments(i)%first_y
  ss = 1
  DO WHILE (ss <= segments(i)%seg_points)
    s_sin_theta_lat(ss) = sin_theta_latitude(ii,jj)
    s_totalppn(ss)      = totalppn(ii,jj)
    ss = ss + 1
    IF (ii == tdims%i_end) THEN
      ii = tdims%i_start
      jj = jj + 1
    ELSE
      ii = ii + 1
    END IF
  END DO
  
  DO k = tkfix1start, tdims%k_end
    ii = segments(i)%first_x
    jj = segments(i)%first_y
    ss = 1
    DO WHILE (ss <= segments(i)%seg_points)
      s_nbv(ss,k) = nbv(ii,jj,k)
        s_rho_th(ss,k) = rho_th(ii,jj,k)
      ss = ss + 1
      IF (ii == tdims%i_end) THEN
        ii = tdims%i_start
        jj = jj + 1
      ELSE
        ii = ii + 1
      END IF
    END DO
  END DO
  
  DO jdir = 1, idir
    DO k = tkfix1start, tdims%k_end
      ii = segments(i)%first_x
      jj = segments(i)%first_y
      ss = 1
      DO WHILE (ss <= segments(i)%seg_points)
        s_udotk(ss,k,jdir)    = udotk(ii,jj,k,jdir)
        ss = ss + 1
        IF (ii == tdims%i_end) THEN
          ii = tdims%i_start
          jj = jj + 1
        ELSE
          ii = ii + 1
        END IF
      END DO
    END DO
  END DO

  CALL gw_ussp_core(one, segments(i)%seg_points, launchlev,               &
         one, segments(i)%seg_points, one, one,                           &
         tkfix1start, tdims%k_end, ussp_launch_factor, wavelstar,         &
         cgw_scale_factor, two_omega, s_sin_theta_lat, s_rho_th, s_nbv,   &
         s_udotk, s_totalppn, s_fptot, L_add_cgw)

  DO jdir = 1, idir
    DO k = tkfix1start, tdims%k_end
      ii = segments(i)%first_x
      jj = segments(i)%first_y
      ss = 1
      DO WHILE (ss <= segments(i)%seg_points)
        fptot(ii,jj,k,jdir) = s_fptot(ss,k,jdir)
        ss = ss + 1
        IF (ii == tdims%i_end) THEN
          ii = tdims%i_start
          jj = jj + 1
        ELSE
          ii = ii + 1
        END IF
      END DO
    END DO
  END DO

  DEALLOCATE(s_fptot)
  DEALLOCATE(s_udotk)
  DEALLOCATE(s_nbv)
  DEALLOCATE(s_rho_th)
  DEALLOCATE(s_totalppn)
  DEALLOCATE(s_sin_theta_lat)
END DO
!$OMP END DO

!$OMP END PARALLEL

DEALLOCATE(segments)
! ----------------------------------------------------------------------+-------
!

!$OMP  PARALLEL DEFAULT(NONE) SHARED(fptot, g_x, g_y, cosphi,           &
!$OMP& p_layer_boundaries, g_xp_smallhalo, g_yp_smallhalo, timestep,    &
!$OMP& sinphi, pdims, minlaunchlev, L_ussp_heating, r_theta_levels,     &
!$OMP& r_rho_levels, uonp, vonp, T_inc, tdims, g, cp)                   &
!$OMP& PRIVATE(g_g, i, j, k, jdir, jj, dzb, dzu, dzl, uhat, vhat,       &
!$OMP& ududt, vdvdt, omp_block)

!
! ----------------------------------------------------------------------+-------
! 5.0   Compute vertical divergence of pseudomomentum flux.
!       Column Integral Converts : [Vertical]   T onto RHO grid
! ----------------------------------------------------------------------+-------

! gives each thread the largest block possible to execute

omp_block = pdims%k_end-(minlaunchlev+1)+1
!$ omp_block = CEILING(REAL(((pdims%k_end - (minlaunchlev+1)) + 1))/    &
!$ omp_get_num_threads())

!$OMP DO SCHEDULE(STATIC)
omp_blocking2: DO jj=minlaunchlev+1, pdims%k_end, omp_block
  Levels_do14: DO k=jj, MIN(jj+omp_block-1,pdims%k_end)
    Rows_do14: DO j=pdims%j_start,pdims%j_end
      Row_length_do14: DO i=pdims%i_start,pdims%i_end
        IDir_do5: DO jdir=1,idir
!       Pseudomomentum flux
          g_g = g * (fptot(i,j,k,jdir) - fptot(i,j,k-1,jdir)) /         &
           (p_layer_boundaries(i,j,k) - p_layer_boundaries(i,j,k-1))
          g_x(i,j,k) = g_x(i,j,k) + g_g* cosphi(jdir)
          g_y(i,j,k) = g_y(i,j,k) + g_g* sinphi(jdir)
        END DO  IDir_do5
      END DO  Row_length_do14
    END DO  Rows_do14
  END DO  Levels_do14
END DO  omp_blocking2
!$OMP END DO

!
! ----------------------------------------------------------------------+-------
! 5.1   Wind and temperature increments from wave dissipation
! ----------------------------------------------------------------------+-------

!$OMP DO SCHEDULE(STATIC)
Levels_do15: DO k=minlaunchlev,pdims%k_end
  Rows_do15: DO j=pdims%j_start,pdims%j_end
    Row_length_do15: DO i=pdims%i_start,pdims%i_end
!     g_xp_smallhalo(i,j,k) = g_x(i,j,k)
!     g_yp_smallhalo(i,j,k) = g_y(i,j,k)
!! ABOVE is preferred but following should reduce bit comparison differences
!! for purposes of initial testing. REMOVE when scheme goes LIVE!
      g_xp_smallhalo(i,j,k) = timestep * g_x(i,j,k)
      g_yp_smallhalo(i,j,k) = timestep * g_y(i,j,k)
    END DO  Row_length_do15
  END DO  Rows_do15
END DO  Levels_do15
!$OMP END DO

!-----------------------------------------------------------------------+-------
! Calculate heating due to gravity wave dissipation
!-----------------------------------------------------------------------+-------
GW_heating: IF ( L_ussp_heating ) THEN

!$OMP DO SCHEDULE(STATIC)
  Levels_do16: DO k = tkfix1start,tdims%k_end-1
    Rows_do16: DO j = tdims%j_start,tdims%j_end
      Row_length_do16: DO i = tdims%i_start,tdims%i_end

        dzb = r_theta_levels(i,j,k) -  r_rho_levels(i,j,k)
        dzu = r_rho_levels(i,j,k+1) -  r_theta_levels(i,j,k)
        dzl = r_rho_levels(i,j,k+1) -  r_rho_levels(i,j,k)

!       u and v on theta_level(k)
        uhat   = dzu * uonp(i,j,k) + dzb * uonp(i,j,k+1)
        vhat   = dzu * vonp(i,j,k) + dzb * vonp(i,j,k+1)

!       u*du and v*dv on theta_level(k)
        ududt  = uhat *( dzu * g_x(i,j,k) + dzb * g_x(i,j,k+1) )
        vdvdt  = vhat *( dzu * g_y(i,j,k) + dzb * g_y(i,j,k+1) )

!       dT/dt on theta_level(k)
        T_inc(i,j,k)=T_inc(i,j,k) - timestep * ( ududt + vdvdt ) /      &
                    ( cp * dzl * dzl )

      END DO  Row_length_do16
    END DO  Rows_do16
  END DO  Levels_do16
!$OMP END DO

END IF GW_heating

!$OMP END PARALLEL

!
! ----------------------------------------------------------------------+-------
! Put U,V increments onto U,V grid after initialising increments.
! Interpolate : [Horizontal] P onto U,V grid
! ----------------------------------------------------------------------+-------
!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP& SHARED( udims, uincrt, vdims, vincrt )                           &
!$OMP& PRIVATE( i, j, k )
!$OMP DO SCHEDULE(STATIC) 
DO k = udims%k_start, udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      uincrt(i,j,k) = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = vdims%k_start, vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      vincrt(i,j,k) = 0.0
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
! ENDGame swaps West - ND swaps East
CALL swap_bounds(g_xp_smallhalo,                                        &
                row_length,rows,levels,                                 &
                off_x,off_y,fld_type_p, swap_field_is_scalar,           &
                do_west_arg=.TRUE., do_east_arg= .FALSE.)
!
! Interpolate : [Horizontal] P onto U grid
!
!$OMP PARALLEL DEFAULT(NONE)                                            &
!$OMP SHARED( g_xp_smallhalo, pdims_s, udims, levels, uincrt )
CALL p_to_u(g_xp_smallhalo,                                             &
            pdims_s%i_start,pdims_s%i_end,pdims_s%j_start,pdims_s%j_end,&
            udims%i_start,  udims%i_end,  udims%j_start,  udims%j_end,  &
            1, levels, uincrt)
!$OMP END PARALLEL
!
! ----------------------------------------------------------------------+-------
! Diagnostic output : Acceleration of Zonal Wind (on rho levels)
! ----------------------------------------------------------------------+-------
! No correction needed for U-increment at poles (initialised to zero)
GWspec_Acc_if1ew: IF (gwspec_ewacc_on .OR. gwspec_ewacc_p_on) THEN
  Levels_do17ew: DO k=udims%k_start,udims%k_end
    Rows_do17ew: DO j=udims%j_start,udims%j_end
      Row_length_do17ew: DO i=udims%i_start,udims%i_end
!       gwspec_ewacc(i,j,k) = uincrt(i,j,k)
!!  ABOVE is preferred but following should reduce bit comparison differences
!!  for purposes of initial testing. REMOVE when scheme goes LIVE!
        gwspec_ewacc(i,j,k) = uincrt(i,j,k) / timestep
      END DO  Row_length_do17ew
    END DO  Rows_do17ew
  END DO  Levels_do17ew
END IF GWspec_Acc_if1ew
! ----------------------------------------------------------------------+-------
!
! ENDGame swaps South - ND swaps North
CALL swap_bounds(g_yp_smallhalo,                                        &
                row_length,rows,levels,                                 &
                off_x,off_y,fld_type_p, swap_field_is_scalar,           &
                do_south_arg=.TRUE., do_north_arg= .FALSE.)
!
! Interpolate : [Horizontal] P onto V grid
!
!$OMP PARALLEL DEFAULT(NONE)                                            &
!$OMP SHARED( g_yp_smallhalo, pdims_s, vdims, levels, vincrt )
CALL p_to_v(g_yp_smallhalo,                                             &
            pdims_s%i_start,pdims_s%i_end,pdims_s%j_start,pdims_s%j_end,&
            vdims%i_start,  vdims%i_end,  vdims%j_start,  vdims%j_end,  &
            1, levels, vincrt)
!$OMP END PARALLEL
!

! Correct V-increment at poles before updating model prognostic
Type_global_if1: IF (model_type == mt_global) THEN
  IF ( at_extremity(psouth) ) THEN
!
    CALL eg_v_at_poles(uincrt, vincrt, 1.0,                             &
                       udims%j_start, vdims%j_start, udims, vdims)
  END IF
!
  IF( at_extremity(pnorth) ) THEN
!
    CALL eg_v_at_poles(uincrt, vincrt, -1.0,                            &
                       udims%j_end, vdims%j_end, udims, vdims)
  END IF
END IF  Type_global_if1

!
! ----------------------------------------------------------------------+-------
! Diagnostic output : Acceleration of Meridional Wind (on rho levels)
! ----------------------------------------------------------------------+-------
GWspec_Acc_if1ns: IF (gwspec_nsacc_on) THEN
  Levels_do18ns: DO k=vdims%k_start,vdims%k_end
    Rows_do18ns: DO j=vdims%j_start,vdims%j_end
      Row_length_do18ns: DO i=vdims%i_start,vdims%i_end
!       gwspec_nsacc(i,j,k) = vincrt(i,j,k)
!!  ABOVE is preferred but following should reduce bit comparison differences
!!  for purposes of initial testing. REMOVE when scheme goes LIVE!
        gwspec_nsacc(i,j,k) = vincrt(i,j,k) / timestep
      END DO  Row_length_do18ns
    END DO  Rows_do18ns
  END DO  Levels_do18ns
END IF  GWspec_Acc_if1ns
! ----------------------------------------------------------------------+-------
! Add increments to wind and temperature
! ----------------------------------------------------------------------+-------
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& SHARED( minlaunchlev, udims, r_u, uincrt, vdims, r_v, vincrt )   &
!$OMP& PRIVATE( i, j, k )
DO k=minlaunchlev,udims%k_end
! NB: assumes that U and V winds will always be on same vertical levels
  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end
!     r_u(i,j,k) = r_u(i,j,k) + (timestep * uincrt(i,j,k))
!!  ABOVE is preferred but following should reduce bit comparison differences
!!  for purposes of initial testing. REMOVE when scheme goes LIVE!
      r_u(i,j,k) = r_u(i,j,k) + (uincrt(i,j,k))
    END DO
  END DO
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end
!     r_v(i,j,k) = r_v(i,j,k) + (timestep * vincrt(i,j,k))
!!  ABOVE is preferred but following should reduce bit comparison differences
!!  for purposes of initial testing. REMOVE when scheme goes LIVE!
      r_v(i,j,k) = r_v(i,j,k) + (vincrt(i,j,k))
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!
! ----------------------------------------------------------------------+-------
! Diagnostic output : Northward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
GWspec_Flux_if1n: IF (gwspec_nflux_on) THEN
  Levels_do20n: DO k=tkfix1start,tdims%k_end
    Rows_do20n: DO j=tdims%j_start,tdims%j_end
      Row_length_do20n: DO i=tdims%i_start,tdims%i_end
        work_TsmallHALO(i,j,k) = fptot(i,j,k,1)
      END DO  Row_length_do20n
    END DO  Rows_do20n
  END DO  Levels_do20n
!
! Usually, ENDGame swaps south - ND swaps north
! Here, ENDGame swaps North and South as interpolation to northernmost global 
!   v-row requires data in small halos on p-points
! ND swaps north only as southernmost v-row can be calculated from p-points
!   without small halos
  CALL swap_bounds(                                                     &
       work_TsmallHALO,                                                 &
       row_length, rows,                                                &
       levels, off_x, off_y, fld_type_p, swap_field_is_scalar,          &
       do_south_arg=.TRUE., do_north_arg= .TRUE.)
!
! Interpolate : [Horizontal] P onto V grid
!
  CALL p_to_v(work_TsmallHALO(tdims_s%i_start,tdims_s%j_start,1),       &
              pdims_s%i_start, pdims_s%i_end,                           &
              pdims_s%j_start, pdims_s%j_end,                           &
              vdims%i_start, vdims%i_end, vdims%j_start, vdims%j_end,   &
              1, levels, gwspec_nflux(vdims%i_start,vdims%j_start,1))
!    

  Type_global_if2n: IF (model_type == mt_global) THEN
    ! set polar rows to common mean value.
    CALL polar_row_mean(gwspec_nflux(vdims%i_start,vdims%j_start,1),  &
                        vdims%i_start,vdims%i_end,                    &
                        vdims%j_start,vdims%j_end,                    &
                        1,levels,                                     &
                        global_row_length,                            &
                        n_proc, n_procy, proc_row_group,              &
                        at_extremity)
  END IF  Type_global_if2n
!
END IF  GWspec_Flux_if1n
!
! ----------------------------------------------------------------------+-------
! Diagnostic output : Westward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
GWspec_Flux_if1w: IF (gwspec_wflux_on .OR. gwspec_wflux_p_on) THEN
  Levels_do20w: DO k=tkfix1start,tdims%k_end
    Rows_do20w: DO j=tdims%j_start,tdims%j_end
      Row_length_do20w: DO i=tdims%i_start,tdims%i_end
        Work_TsmallHALO(i,j,k) = fptot(i,j,k,2)
      END DO  Row_length_do20w
    END DO  Rows_do20w
  END DO  Levels_do20w
!
! ENDGame swaps west - ND swaps east
  CALL swap_bounds(                                                     &
       work_TsmallHALO,                                                 &
       row_length, rows,                                                &
       levels, off_x, off_y, fld_type_p, swap_field_is_scalar,          &
       do_west_arg=.TRUE., do_east_arg= .FALSE.)
!
! Interpolate : [Horizontal] P onto U grid
!
  CALL p_to_u(work_TsmallHALO(tdims_s%i_start,tdims_s%j_start,1),       &
              pdims_s%i_start, pdims_s%i_end,                           &
              pdims_s%j_start, pdims_s%j_end,                           &
              udims%i_start, udims%i_end, udims%j_start,udims%j_end,    &
              1, levels, gwspec_wflux(udims%i_start,udims%j_start,1))
!
END IF  GWspec_Flux_if1w
!
! ----------------------------------------------------------------------+-------
! Diagnostic output : Southward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
GWspec_Flux_if1s: IF (gwspec_sflux_on) THEN
  Levels_do20s: DO k=tkfix1start,tdims%k_end
    Rows_do20s: DO j=tdims%j_start,tdims%j_end
      Row_length_do20s: DO i=tdims%i_start,tdims%i_end
        work_TsmallHALO(i,j,k) = fptot(i,j,k,3)
      END DO  Row_length_do20s
    END DO  Rows_do20s
  END DO  Levels_do20s
!
! Usually, ENDGame swaps south - ND swaps north
! Here, ENDGame swaps North and South as interpolation to northernmost global
!   v-row requires data in small halos on p-points
! ND swaps north only as southernmost v-row can be calculated from p-points
!   without small halos
  CALL swap_bounds(                                                     &
       work_TsmallHALO,                                                 &
       row_length, rows,                                                &
       levels, off_x, off_y, fld_type_p, swap_field_is_scalar,          &
       do_south_arg=.TRUE., do_north_arg=.TRUE.)
!
! Interpolate : [Horizontal] P onto V grid
!
  CALL p_to_v(work_TsmallHALO(tdims_s%i_start,tdims_s%j_start,1),       &
              pdims_s%i_start, pdims_s%i_end,                           &
              pdims_s%j_start, pdims_s%j_end,                           &
              vdims%i_start, vdims%i_end, vdims%j_start, vdims%j_end,   &
              1, levels, gwspec_sflux(vdims%i_start,vdims%j_start,1))
!    
 
  Type_global_if2s: IF (model_type == mt_global) THEN
    ! set polar rows to common mean value.
    CALL polar_row_mean(gwspec_sflux(vdims%i_start,vdims%j_start,1),  &
                        vdims%i_start,vdims%i_end,                    &
                        vdims%j_start,vdims%j_end,                    &
                        1,levels,                                     &
                        global_row_length,                            &
                        n_proc, n_procy, proc_row_group,              &
                        at_extremity)
  END IF  Type_global_if2s
!
END IF  GWspec_Flux_if1s
!
! ----------------------------------------------------------------------+-------
! Diagnostic output : Eastward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
GWspec_Flux_if1e: IF (gwspec_eflux_on .OR. gwspec_eflux_p_on) THEN
  Levels_do20e: DO k=tkfix1start,tdims%k_end
    Rows_do20e: DO j=tdims%j_start,tdims%j_end
      Row_length_do20e: DO i=tdims%i_start,tdims%i_end
        work_TsmallHALO(i,j,k) = fptot(i,j,k,4)
      END DO  Row_length_do20e
    END DO  Rows_do20e
  END DO  Levels_do20e
!
! ENDGame swaps West - ND swaps East
  CALL swap_bounds(                                                     &
       work_TsmallHALO,                                                 &
       row_length, rows,                                                &
       levels, off_x, off_y, fld_type_p, swap_field_is_scalar,          &
       do_west_arg=.TRUE., do_east_arg= .FALSE.)
!
! Interpolate : [Horizontal] P onto U grid
!
  CALL p_to_u(work_TsmallHALO(tdims_s%i_start,tdims_s%j_start, 1),      &
              pdims_s%i_start, pdims_s%i_end,                           &
              pdims_s%j_start, pdims_s%j_end,                           &
              udims%i_start, udims%i_end, udims%j_start, udims%j_end,   &
              1, levels, gwspec_eflux(udims%i_start,udims%j_start,1))
END IF  GWspec_Flux_if1e
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
!
END SUBROUTINE gw_ussp

END MODULE gw_ussp_mod

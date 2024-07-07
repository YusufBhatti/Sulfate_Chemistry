MODULE TWPICE_netcdf

! Description:
!
!
! *********************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *********************************COPYRIGHT*********************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model


  ! Modules
USE netcdf
USE scm_utils
USE s_interp_mod, ONLY: interp1d
USE s_main_force, ONLY: netcdf_file
USE netCDF_arr
USE global_SCMop, ONLY: incdf
USE MCC_data
USE TWPICE_data

USE planet_constants_mod,  ONLY: g
USE conversions_mod,  ONLY: rsec_per_day 

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

IMPLICIT NONE

REAL ::           &
  file_lat        &
, file_long       &
, file_albsoil

REAL, ALLOCATABLE ::   &
  file_t (:,:)         &! Temperature
, file_q (:,:)         &! H2o mixing ratio
, file_u (:,:)         &! Zonal wind
, file_v  (:,:)        &! Meridional wind
, file_omega (:,:)     &! Rate of change of pressure (dp/dt)
, file_divt (:,:)      &! Hor.  adv. of s (dry static temperature)
, file_divq (:,:)      &! Hor.  adv. of q
, file_vertdivt (:,:)  &! Vert. adv. of s (dry static temperature)
, file_vertdivq (:,:)  &! Vert. adv. of q
, file_w   (:,:)       &! Vertical wind         (Calculated)
, file_th  (:,:)        ! Potential temperature (Calculated)

REAL, ALLOCATABLE :: &
  file_time(:)       &! Time levels
, file_year(:)       &
, file_month(:)      &
, file_day(:)        &
, file_hour(:)       &
, file_min (:)       &
, file_z_prf(:)      &! z of obs levels (assumed fixed at initial run time)
, file_layer_thk(:)  &! Layer thickness of obs level (Using mean T of layer)
, file_flx_H (:)     &! Srf. sensible heat flux
, file_flx_E (:)     &! Srf. latent   heat flux
, file_tsskin(:)     &! Srf. Skin temperature
, file_ts_air(:)     &! Srf. Air temperature
, file_p_srf (:)     &! Srf. Mean Pressure
, file_omega_srf(:)  &! Srf. Rate of change of pressure (dp/dt)
, file_qs_srf (:)    &! Srf. Saturation mixing ratio
, file_u_srf  (:)    &! Srf. Zonal wind
, file_v_srf  (:)    &! Srf. Meridional wind
, file_rhs_air(:)    &! Srf. Air Relative humidity
, file_w_srf  (:)    &! Srf. Vertical velocity (Calculated)
, file_th_srf (:)    &! Srf. Potential temperature (Calculated)
, file_q_srf  (:)    &! Srf. Mixing Ratio (Calculated)
, file_p      (:)     ! Pressure levels


INTEGER(incdf) :: &
  STATUS          &! NetCDF Error status
, ncid            &! NetCDF File id
, time_dimid      &! Time  dimension id
, nlev_dimid      &! Nlevs dimension id
, slev_dimid      &! n soil levels dimension id
, xtype           &
, nt_nc           &! Number of time levels in file
, nk_in           &! Number of levels in file
, attnum          &!
, coor_lena       &!
, coor_lenb       &!
, nk_half_in

! Variable Ids for:
INTEGER(incdf) :: &
  time_id         &! Time (days from 2005-12-31, 00:00:00)
, date_id         &!
, year_id         &!
, mon_id          &!
, day_id          &!
, hour_id         &!
, min_id          &!
, sec_id          &!
, lat_id          &! Latitude
, lon_id          &! Longitude
, p_in_id         &! Pressure
, t_id            &! Temperature
, q_id            &! H2o vapour Mixing ratio
, ql_id           &! H2o Liquid Mixing Ratio (kg/kg)
, qi_id           &! H2o Ice    Mixing Ratio (kg/kg)
, u_id            &! Zonal wind
, v_id            &! Meridional wind
, omega_id        &! Rate of change of pressure (dp/dt)
, hdivt_id        &! Hor.  adv. of s (dry static temperature)
, hdivq_id        &! Hor.  adv. of q
, vdivt_id        &! Vert. adv. of s (dry static temperature)
, vdivq_id        &! Vert. adv. of q
, alb_id          &! Albedo
, alb_snow_id     &! Snow Albedo
, flx_H_id        &! Srf. sensible heat flux
, flx_E_id         ! Srf. latent   heat flux

INTEGER(incdf) :: &
  open_sst_id     &!
, orog_id         &!
, t_soil_id       &!
, q_soil_id       &!
, z0m_id          &!
, z0h_id          &!
, tsskin_id       &! Srf. Skin temperature
, ts_air_id       &! Srf. Air temperature
, ps_id           &! Srf. Mean Pressure
, omegas_id       &! Srf. Rate of change of pressure (dp/dt)
, qs_srf_id       &! Srf. Saturation mixing ratio
, u_srf_id        &! Srf. Zonal wind
, v_srf_id        &! Srf. Meridional wind
, land_sea_id     &!
, rhs_air_id       ! Srf. Air Relative humidity

REAL, PRIVATE, ALLOCATABLE :: dummy(:,:)

CONTAINS
!============================================================================

SUBROUTINE alloc_twpice_nc

IMPLICIT NONE


STATUS = nf90_open(TRIM(ADJUSTL(netcdf_file)), nf90_noWrite, ncid)

IF (STATUS /= nf90_noerr) THEN
  WRITE(umMessage,*)'Error opening netcdf file: '// &
      TRIM(ADJUSTL(netcdf_file))
  CALL umPrint(umMessage,src='TWPICE_netCDF')
END IF

! Get NetCDF file Dimension ID numbers
STATUS = nf90_inq_dimid (ncid, 'time', time_dimid)
STATUS = nf90_inq_dimid (ncid, 'lev',  nlev_dimid)


! Get Dimension sizes
STATUS = nf90_inquire_dimension (ncid, time_dimid, LEN=nt_nc)
STATUS = nf90_inquire_dimension (ncid, nlev_dimid, LEN=nk_in)


! Allocate arrays to receive obs data 0 element to hold srf face values
! 1D arrays
ALLOCATE                   &
  ( file_flx_h     (nt_nc) &
  , file_flx_e     (nt_nc) &
  , file_tsskin    (nt_nc) &
  , file_time      (nt_nc) &
  , file_year      (nt_nc) &
  , file_month     (nt_nc) &
  , file_day       (nt_nc) &
  , file_hour      (nt_nc) &
  , file_min       (nt_nc) &
  , file_ts_air    (nt_nc) &
  , file_th_srf    (nt_nc) &
  , file_p_srf     (nt_nc) &
  , file_omega_srf (nt_nc) &
  , file_qs_srf    (nt_nc) &
  , file_q_srf     (nt_nc) &
  , file_rhs_air   (nt_nc) &
  , file_u_srf     (nt_nc) &
  , file_v_srf     (nt_nc) &
  , file_w_srf     (nt_nc) &
  , file_p       (0:nk_in) & ! srf calc'd depends on starting point
  , file_z_prf   (0:nk_in) & ! add srf, i.e. 0m
  , file_layer_thk (nk_in) )

! 2D Arrays
! (extra k level so to add srf value to array for interpolation)
ALLOCATE                           &
  ( file_divT     (nt_nc, 0:nk_in) &! Assume srf value = lowest obs
  , file_vertdivT (nt_nc, 0:nk_in) &! Assume srf value = lowest obs
  , file_divq     (nt_nc, 0:nk_in) &! Assume srf value = lowest obs
  , file_vertdivq (nt_nc, 0:nk_in) &! Assume srf value = lowest obs
  , file_u        (nt_nc, 0:nk_in) &
  , file_v        (nt_nc, 0:nk_in) &
  , file_w        (nt_nc, 0:nk_in) &
  , file_t        (nt_nc, 0:nk_in) &
  , file_th       (nt_nc, 0:nk_in) &
  , file_q        (nt_nc, 0:nk_in) &
  , file_omega    (nt_nc, 0:nk_in) )

STATUS = nf90_close(ncid)

! Initialise Arrays
file_flx_h     = rmdi  ;  file_flx_e  = rmdi  ;  file_tsskin = rmdi
file_time      = rmdi  ;  file_year   = rmdi  ;  file_month  = rmdi
file_day       = rmdi  ;  file_hour   = rmdi  ;  file_min    = rmdi
file_ts_air    = rmdi  ;  file_th_srf = rmdi  ;  file_p_srf  = rmdi
file_omega_srf = rmdi  ;  file_qs_srf = rmdi  ;  file_q_srf  = rmdi
file_rhs_air   = rmdi  ;  file_u_srf  = rmdi  ;  file_v_srf  = rmdi
file_w_srf     = rmdi  ;  file_p      = rmdi  ;  file_z_prf  = rmdi
file_layer_thk = rmdi


file_divt      = rmdi  ;  file_vertdivt = rmdi ; file_divq   = rmdi
file_vertdivq  = rmdi  ;  file_u        = rmdi ; file_v      = rmdi
file_w         = rmdi  ;  file_t        = rmdi ; file_th     = rmdi
file_q         = rmdi  ;  file_omega    = rmdi


END SUBROUTINE alloc_twpice_nc

!============================================================================

SUBROUTINE dealloc_twpice_nc

IMPLICIT NONE


DEALLOCATE(                                                      &
  ! 1D Arrays
  file_flx_H,     file_flx_E,     file_tsskin,    file_time,     &
  file_year,      file_month,     file_day,       file_hour,     &
  file_min,       file_ts_air,    file_th_srf,    file_p_srf,    &
  file_omega_srf, file_qs_srf,    file_q_srf,     file_rhs_air,  &
  file_u_srf,     file_v_srf,     file_w_srf,     file_p,        &
  file_z_prf,     file_layer_thk,                                &

  ! 2D Arrays
  file_divT,      file_vertdivT,  file_divq,      file_vertdivq, &
  file_u,         file_v,         file_w,         file_t,        &
  file_th,        file_q,         file_omega)

END SUBROUTINE dealloc_twpice_nc

!============================================================================

SUBROUTINE read_TWPICE_nc

IMPLICIT NONE


STATUS = nf90_open(TRIM(ADJUSTL(netcdf_file)), nf90_noWrite, ncid)

IF (STATUS /= nf90_noerr) THEN
  WRITE(umMessage,*)'Error opening netcdf file: '// &
      TRIM(ADJUSTL(netcdf_file))
  CALL umPrint(umMessage,src='TWPICE_netCDF')
END IF

! Setup variables
STATUS = nf90_inq_varid (ncid, 'time',      time_id)
STATUS = nf90_inq_varid (ncid, 'year',      year_id)
STATUS = nf90_inq_varid (ncid, 'month',      mon_id)
STATUS = nf90_inq_varid (ncid, 'day',        day_id)
STATUS = nf90_inq_varid (ncid, 'hour',      hour_id)
STATUS = nf90_inq_varid (ncid, 'minute',     min_id)
STATUS = nf90_inq_varid (ncid, 'lat',        lat_id)
STATUS = nf90_inq_varid (ncid, 'lon',        lon_id)

! Forcing variables
STATUS = nf90_inq_varid (ncid, 'lev',       p_in_id)
STATUS = nf90_inq_varid (ncid, 'T',            t_id)
STATUS = nf90_inq_varid (ncid, 'q',            q_id)
STATUS = nf90_inq_varid (ncid, 'u',            u_id)
STATUS = nf90_inq_varid (ncid, 'v',            v_id)
STATUS = nf90_inq_varid (ncid, 'omega',    omega_id)
STATUS = nf90_inq_varid (ncid, 'divT',     hdivt_id)
STATUS = nf90_inq_varid (ncid, 'divq',     hdivq_id)
STATUS = nf90_inq_varid (ncid, 'vertdivT', vdivt_id)
STATUS = nf90_inq_varid (ncid, 'vertdivq', vdivq_id)
STATUS = nf90_inq_varid (ncid, 'alb',        alb_id)
STATUS = nf90_inq_varid (ncid, 'shflx',    flx_H_id)
STATUS = nf90_inq_varid (ncid, 'lhflx',    flx_E_id)
STATUS = nf90_inq_varid (ncid, 'Tskin',   tsskin_id)
STATUS = nf90_inq_varid (ncid, 'Tsair',   ts_air_id)
STATUS = nf90_inq_varid (ncid, 'Ps',          ps_id)
STATUS = nf90_inq_varid (ncid, 'omegas',  omegas_id)
STATUS = nf90_inq_varid (ncid, 'qs',      qs_srf_id)
STATUS = nf90_inq_varid (ncid, 'RHsair', rhs_air_id)
STATUS = nf90_inq_varid (ncid, 'usrf',     u_srf_id)
STATUS = nf90_inq_varid (ncid, 'vsrf',     v_srf_id)


!=======================
! Get obs data from file
!=======================

! Get set-up variables
!=======================
STATUS = nf90_get_var (ncid, lat_id,  file_lat)
STATUS = nf90_get_var (ncid, lon_id,  file_long)
STATUS = nf90_get_var (ncid, alb_id,  file_albsoil)
STATUS = nf90_get_var (ncid, year_id, file_year)
STATUS = nf90_get_var (ncid, mon_id,  file_month)
STATUS = nf90_get_var (ncid, day_id,  file_day)
STATUS = nf90_get_var (ncid, hour_id, file_hour)
STATUS = nf90_get_var (ncid, min_id,  file_min)



! Get surface timeseries of:
!===========================
STATUS = nf90_get_var (ncid, flx_H_id,   file_flx_H)     ! (W/m2)
STATUS = nf90_get_var (ncid, flx_E_id,   file_flx_E)     ! (W/m2)
STATUS = nf90_get_var (ncid, tsskin_id,  file_tsskin)    ! (K)
STATUS = nf90_get_var (ncid, qs_srf_id,  file_qs_srf)    ! (kg/kg)
STATUS = nf90_get_var (ncid, rhs_air_id, file_rhs_air)   ! (%)
STATUS = nf90_get_var (ncid, ts_air_id,  file_ts_air)    ! (K)
STATUS = nf90_get_var (ncid, ps_id,      file_p_srf)     ! (mb)
STATUS = nf90_get_var (ncid, u_srf_id,   file_u_srf)     ! (m/s)
STATUS = nf90_get_var (ncid, v_srf_id,   file_v_srf)     ! (m/s)
STATUS = nf90_get_var (ncid, omegas_id,  file_omega_srf) ! (mb/hr)



! Get time levels
!=========================

STATUS = nf90_get_var (ncid, time_id, file_time) ! Time (days since
                                                 ! 2005-12-31, 00:00:00)


obs_pd_nc = (file_time(2) - file_time(1))*rsec_per_day
yeari_nc  = file_year(obs_t0)
monthi_nc = file_month(obs_t0)
dayi_nc   = file_day(obs_t0)
houri_nc  = file_hour(obs_t0)
mini_nc   = file_min(obs_t0)

albsoil_nc = file_albsoil
lat_nc     = file_lat
long_nc    = file_long
tstari_nc  = file_tsskin(obs_t0)

!      ndayin = INT((  time_in(nt_nc)                            &
!                    - time_in(obs_t0)))

!      nminin = INT((  time_in(nt_nc)                            &
!                    - time_in(obs_t0)                           &
!                    - ndayin)*24.0*60.0)

!      nsecin = INT((((  time_in(nt_nc)                          &
!                      - time_in(obs_t0)                         &
!                      - ndayin)*24.0*60.0)                      &
!                      - nminin)*60.0)


      ! Get pressure co-ord levels
      !===========================
STATUS = nf90_get_var (ncid, p_in_id, file_p(1:nk_in)) ! Pressure(mb)

ALLOCATE( dummy(nk_in, nt_nc) )

! Get atmospheric profiles
!=========================

! Temperature (K)
STATUS = nf90_get_var (ncid, T_id, dummy)
file_t(:,1:nk_in)  = TRANSPOSE(dummy)

! q mixing ratio (g/kg)
STATUS = nf90_get_var (ncid, q_id, dummy)
file_q(:,1:nk_in) = TRANSPOSE(dummy)

! U-wind (m/s)
STATUS = nf90_get_var (ncid, u_id, dummy)
file_u(:,1:nk_in) = TRANSPOSE(dummy)

! V-wind (m/s)
STATUS = nf90_get_var (ncid, v_id, dummy)
file_v(:,1:nk_in) = TRANSPOSE(dummy)

! omega (mb/hr)
STATUS = nf90_get_var (ncid, omega_id, dummy)
file_omega(:,1:nk_in) = TRANSPOSE(dummy)


! Get Forcing tendencies
!=======================

! Horizontal temperature advection
STATUS = nf90_get_var (ncid, hdivt_id, dummy)
file_divT(:,1:nk_in) = TRANSPOSE(dummy)

! Vertical dry static energy(s) advection
STATUS = nf90_get_var (ncid, vdivt_id, dummy)
file_vertdivT(:,1:nk_in) = TRANSPOSE(dummy)

! Horizontal mixing ratio (q) advection
STATUS = nf90_get_var (ncid, hdivq_id, dummy)
file_divq(:,1:nk_in) = TRANSPOSE(dummy)

! Vertical mixing ratio (q) advection
STATUS = nf90_get_var (ncid, vdivq_id, dummy)
file_vertdivq(:,1:nk_in) = TRANSPOSE(dummy)

DEALLOCATE(dummy)

STATUS = nf90_close(ncid)

! Invert Vertical axis, since these were read in order of
! increasing pressure rather than in increasing height.

file_p (1:nk_in)   = file_p(  nk_in:1:-1)
file_t (:,1:nk_in) = file_t(:,nk_in:1:-1)
file_q (:,1:nk_in) = file_q(:,nk_in:1:-1)
file_u (:,1:nk_in) = file_u(:,nk_in:1:-1)
file_v (:,1:nk_in) = file_v(:,nk_in:1:-1)

file_omega    (:,1:nk_in) = file_omega    (:,nk_in:1:-1)
file_divT     (:,1:nk_in) = file_divT     (:,nk_in:1:-1)
file_vertdivT (:,1:nk_in) = file_vertdivT (:,nk_in:1:-1)
file_divq     (:,1:nk_in) = file_divq     (:,nk_in:1:-1)
file_vertdivq (:,1:nk_in) = file_vertdivq (:,nk_in:1:-1)

END SUBROUTINE read_TWPICE_nc

!==============================================================================

!==============================================================================

SUBROUTINE TWPICE_nc_var_to_scm

IMPLICIT NONE


INTEGER :: i,j,k

!=================================
! Add Surface values to the arrays
!=================================
file_p(0)          = file_p_srf(obs_t0)         ! Pressure (Pa)
file_t(:,0)        = file_ts_air                ! Temperature (K)
file_q(:,0)        = (file_rhs_air(:)/100.0) &  ! q mix rat. (g/kg)
                   * file_qs_srf(:)*1000.0
file_u(:,0)        = file_u_srf(:)              ! U-wind (m/s)
file_v(:,0)        = file_v_srf(:)              ! V-wind (m/s)

file_omega(:,0)    = file_omega_srf(:)          ! omega (mb/hr)
file_divT(:,0)     = file_divT(:,1)     ! Hor. temperature adv.
file_vertdivT(:,0) = file_vertdivT(:,1) ! Ver. dry static energy(s) adv.
file_divq(:,0)     = file_divq(:,1)     ! Hor. mixing ratio (q) adv.
file_vertdivq(:,0) = file_vertdivq(:,1) ! Ver. mixing ratio (q) adv.


!===========================================
! Scale/Convert to variables required by SCM
!===========================================

! Pressure (Pa)
file_p(:) = file_p(:) * 100.0

! Mixing ratio: g/kg --> kg/kg
file_q(:,:)        = 1.0e-3 * file_q(:,:)         ! Note: File comments
file_divq(:,:)     = 1.0e-3 * file_divq(:,:)      ! may be wrong, file
file_vertdivq(:,:) = 1.0e-3 * file_vertdivq(:,:)  ! specifies kg/kg but
                                                  ! appear to be in g/kg.
WHERE (file_q == -0.0)        file_q = 0.0        ! Remove -0.0 values
WHERE (file_divq == -0.0)     file_divq = 0.0     ! Remove -0.0 values
WHERE (file_vertdivq == -0.0) file_vertdivq = 0.0 ! Remove -0.0 values

! Mixing ratio ---> Specific humidity
file_q(:,:)        =  file_q(:,:)        / (file_q(:,:)        + 1)
file_divq(:,:)     =  file_divq(:,:)     / (file_divq(:,:)     + 1)
file_vertdivq(:,:) =  file_vertdivq(:,:) / (file_vertdivq(:,:) + 1)

! Omega(mb/hr) --> w-wind(m/s) using -Omega*R*T/P
! Assuming surface pressure does not change
DO j=1, nt_nc
  file_w(j,:) = - r * file_t(j,:)                               &
              * file_omega(j,:)*(1.0/36.0)                      &
              / (file_p(:) * g )
END DO

WHERE (file_w == -0.0) file_w = 0.0 ! Stop creation of -0.0 values

! Temperature  --> Theta
! Assumeing surface pressure does not change
DO j=1, nt_nc
  file_th(j,:) = file_t(j,:)                                    &
                * ((1.0e+5/file_p(:))**(r/cp))
END DO


! Calculate z-profile
!====================

!=================================================================!
! ASSUMING z-heights of forcing data does not vary significantly, !
!          Use initial profile for interpolation on UM grid.      !
!=================================================================!

! Set Surface height to 0m
file_z_prf(0)    = 0.0

DO k=1, nk_in
  file_layer_thk(k) = -(r/g)                                       &
                    * ((  file_t(obs_t0,k)                         &
                      + file_t(obs_t0,k-1))/2.0)                   &
                      * (LOG(  file_p(k)                           &
                             / file_p(k-1)))

  file_z_prf(k)     = file_z_prf(k-1) + file_layer_thk(k)
END DO

END SUBROUTINE TWPICE_nc_var_to_scm

!==============================================================================

SUBROUTINE TWPICE_nc_grid_to_scm

IMPLICIT NONE


INTEGER :: i,j,k,kk

!============================================================
! Ozone
! Ozone is dealt with separately as it is not containing in the
! TWPice netcdf files

!!$      ! Forcing data on UM levels
!!$      ! Initial profiles
!!$      real, intent(inout) :: init_th(:)
!!$      real, intent(inout) :: init_p(:)
!!$      real, intent(inout) :: init_q(:)
!!$      real, intent(inout) :: init_u(:)
!!$      real, intent(inout) :: init_v(:)
!!$      real, intent(inout) :: init_w(:)
!!$      real, intent(inout) :: ozone(:)
!!$
!!$      ! obs forcing
!!$      real, intent(inout) :: t_h_inc(:,:)
!!$      real, intent(inout) :: t_v_inc(:,:)
!!$      real, intent(inout) :: q_h_inc(:,:)
!!$      real, intent(inout) :: q_v_inc(:,:)
!!$      real, intent(inout) :: u_inc(:,:)
!!$      real, intent(inout) :: v_inc(:,:)
!!$      real, intent(inout) :: w_inc(:,:)

DO j=1, rw
  DO i=1, rw_lng
    shflx_nc(i,j,:) = file_flx_H  (obs_t0:nt_nc) ! Sen. heat flux (W/m2)
    lhflx_nc(i,j,:) = file_flx_E  (obs_t0:nt_nc) ! Lat. heat flux (W/m2)
    tstar_nc(i,j,:) = file_tsskin (obs_t0:nt_nc) ! (K)
    CALL interp1d(twp_data_o3_z, twp_data_o3, z_th, o3_nc(i,j,:))
  END DO
END DO

IF (MAXVAL(z_th) > MAXVAL(twp_data_o3_z)) THEN
  CALL alloc_mcc(model_levels)
  CALL interp1d(mcc_trp_z, mcc_trp_o3, z_th, mcc_o3)
  DO k=1, model_levels
    DO j=1, rw
      DO i=1, rw_lng
        IF (o3_nc(i,j,k) == rmdi) o3_nc(i,j,k) = mcc_o3(k)
      END DO
    END DO
  END DO
  CALL dealloc_mcc
END IF

DO j=1, rw
  DO i=1, rw_lng

    CALL interp1d(file_z_prf, file_p(:),             &
                  z_rh,       pi_rh_nc(i,j,:))
    CALL interp1d(file_z_prf, file_q(obs_t0,:),      &
                  z_th,       qi_nc(i,j,:))
    CALL interp1d(file_z_prf, file_th(obs_t0,:),     &
                  z_th,       thi_nc(i,j,:))
    CALL interp1d(file_z_prf, file_u(obs_t0,:),      &
                  z_rh(1:model_levels), ui_nc(i,j,:))
    CALL interp1d(file_z_prf, file_v(obs_t0,:),      &
                  z_rh(1:model_levels), vi_nc(i,j,:))
    CALL interp1d(file_z_prf, file_w(obs_t0,:),      &
                  z_th,       wi_nc(i,j,1:model_levels))
  END DO
END DO

DO kk=obs_t0, nt_nc
  DO j=1, rw
    DO i=1, rw_lng
      ! obs forcing
      CALL interp1d(file_z_prf, file_divt(kk,:),              &
                    z_th, t_h_inc_nc(i,j,kk-obs_t0+1,:))

      CALL interp1d(file_z_prf, file_divq(kk,:),              &
                    z_th, q_h_inc_nc(i,j,kk-obs_t0+1,:))

      CALL interp1d(file_z_prf, file_vertdivt(kk,:),          &
                    z_th, t_v_inc_nc(i,j,kk-obs_t0+1,:))

      CALL interp1d(file_z_prf, file_vertdivq(kk,:),          &
                    z_th, q_v_inc_nc(i,j,kk-obs_t0+1,:))

      CALL interp1d(file_z_prf, file_u(kk,:),                 &
                    z_rh(1:model_levels),                     &
                    u_h_inc_nc(i,j,kk-obs_t0+1,:))

      CALL interp1d(file_z_prf, file_v(kk,:),                 &
                    z_rh(1:model_levels),                     &
                    v_h_inc_nc(i,j,kk-obs_t0+1,:))

      CALL interp1d(file_z_prf, file_w(kk,:),                 &
                    z_th,                                     &
                    w_inc_nc(i,j,kk-obs_t0+1,:))
    END DO
  END DO
END DO


! UM resolution is outside current obs/forcing data
! Use McClatchly routines to fill in missing domain

IF (MAXVAL(z_rh) > MAXVAL(file_z_prf)) THEN

  CALL alloc_mcc(model_levels)

  CALL interp1d(mcc_trp_z, mcc_trp_p,  z_rh, mcc_rh_p)
  CALL interp1d(mcc_trp_z, mcc_trp_th, z_th, mcc_th)
  CALL interp1d(mcc_trp_z, mcc_trp_q,  z_th, mcc_q)

  DO k=1, model_levels
    DO j=1, rw
      DO i=1, rw_lng
        IF (thi_nc (i,j,k) == rmdi) thi_nc(i,j,k) = mcc_th(k)
      END DO
    END DO
  END DO

  DO k=1, model_levels+1
    DO j=1, rw
      DO i=1, rw_lng
        IF (pi_rh_nc(i,j,k) == rmdi) pi_rh_nc(i,j,k) = mcc_rh_p(k)
      END DO
    END DO
  END DO

  DO k=1, model_levels
    DO j=1, rw
      DO i=1, rw_lng
        IF (qi_nc(i,j,k) == rmdi) qi_nc(i,j,k) = mcc_q(k)
      END DO
    END DO
  END DO

  CALL dealloc_mcc


  DO k=1, model_levels
    DO j=1, rw
      DO i=1, rw_lng
        IF (ui_nc(i,j,k) == rmdi) ui_nc(i,j,k) = ui_nc(i,j,k-1)
        IF (vi_nc(i,j,k) == rmdi) vi_nc(i,j,k) = vi_nc(i,j,k-1)
        IF (wi_nc(i,j,k) == rmdi) wi_nc(i,j,k) = wi_nc(i,j,k-1)
      END DO
    END DO
  END DO

  DO k=1, model_levels
    DO kk=1, nobs
      DO j=1, rw
        DO i=1, rw_lng
          IF (t_h_inc_nc(i,j,kk,k) == rmdi)                             &
              t_h_inc_nc(i,j,kk,k) =  t_h_inc_nc(i,j,kk,k-1)

          IF (q_h_inc_nc(i,j,kk,k) == rmdi)                             &
              q_h_inc_nc(i,j,kk,k) =  q_h_inc_nc(i,j,kk,k-1)

          IF (t_v_inc_nc(i,j,kk,k) == rmdi)                             &
              t_v_inc_nc(i,j,kk,k) =  t_v_inc_nc(i,j,kk,k-1)

          IF (q_v_inc_nc(i,j,kk,k) == rmdi)                             &
              q_v_inc_nc(i,j,kk,k) =  q_v_inc_nc(i,j,kk,k-1)

          IF (u_h_inc_nc(i,j,kk,k) == rmdi)                             &
              u_h_inc_nc(i,j,kk,k) =  u_h_inc_nc(i,j,kk,k-1)

          IF (v_h_inc_nc(i,j,kk,k) == rmdi)                             &
              v_h_inc_nc(i,j,kk,k) =  v_h_inc_nc(i,j,kk,k-1)

          IF (w_inc_nc(i,j,kk,k) == rmdi)                               &
              w_inc_nc(i,j,kk,k) =  w_inc_nc(i,j,kk,k-1)
        END DO
      END DO
    END DO
  END DO

END IF

END SUBROUTINE TWPICE_nc_grid_to_scm                                      !

!==============================================================================

END MODULE TWPICE_netCDF

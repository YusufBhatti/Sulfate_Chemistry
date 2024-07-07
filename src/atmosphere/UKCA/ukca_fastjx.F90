! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Main routine for calculating online photolysis rates
!   using fast-jx. Borrows a lot from the fast-j main routine.
!   Further developments that are still required include:
!    (i) use of MODE aerosol optical depths
!    (ii) use of more complete blocking for potential speedup
!    (iii) update to fast-jx 6.6 which will have some small
!    improvements (eg will remove the requirement that jaceto is last)
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
SUBROUTINE ukca_fastjx(                                          &
  dj,                                                            &
  p_layer_boundaries, pstar,                                     &
  p_theta_levels,                                                &
  t_theta_levels,                                                &
  f3u_by_omega,                                                  &
  rho_r2,                                                        &
  z_top_of_model,                                                &
  so4_aitken, so4_accum,                                         &
  qcl, qcf, rel_humid_frac, area_cloud_fraction,                 &
  conv_cloud_lwp, conv_cloud_top, conv_cloud_base,               &
  conv_cloud_amount,                                             &
  surf_albedo,                                                   &
  tropopause_level,                                              &
  ozone,                                                         &
  land_fraction)

USE fastjx_data,    ONLY: fastjx_set_limits,                           &
                          fastjx_allocate_memory,                      &
                          fastjx_deallocate_memory, cmessage,          &
                          jjpnl,nsl,                                   &
                          Blocking,kpcx,jpcl,                          &
                          aer_swband,                                  &
                          tau, daynumber,                              &
                          sza, szafac, sza_2d, szafac_2d, u0,          &
                          sa_block, fl, fl_cyc,                        &
                          rz_3d, rz_all, zzht,                         &
                          pz_3d, pz_all, tz_3d, sa_2d,                 &
                          dm_3d, o3_3d,                                &
                          ods_3d, odw_3d, odi_3d,                      &
                          solcyc_spec
                           
USE atm_fields_bounds_mod,   ONLY:udims
USE ukca_option_mod,         ONLY: jppj,             l_ukca,           &
                                   i_ukca_solcyc
USE level_heights_mod,       ONLY: r_theta_levels,   r_rho_levels,     &
                                   eta_theta_levels
USE trignometric_mod,        ONLY: true_longitude,   sin_theta_latitude
USE spec_sw_lw,              ONLY: sw_spectrum
USE rad_pcf,                 ONLY: ip_aerosol_param_moist,             &
                                   ip_accum_sulphate,                  &
                                   ip_aitken_sulphate
USE UM_ParVars
USE ukca_constants,          ONLY: c_o3, m_s, m_nh42so4
USE planet_constants_mod,    ONLY: g
USE water_constants_mod,     ONLY: tm
USE conversions_mod,         ONLY: pi_over_180
USE umPrintMgr
USE ereport_mod,             ONLY: ereport
USE yomhook,                 ONLY: lhook, dr_hook
USE parkind1,                ONLY: jprb, jpim
USE nlstcall_mod,            ONLY: lcal360
USE cloud_inputs_mod,        ONLY: i_cld_vn
USE pc2_constants_mod,       ONLY: i_cld_pc2
USE cv_run_mod,              ONLY: l_3d_cca
USE fastjx_inphot_mod,       ONLY: fastjx_inphot
USE fastjx_solar2_mod,       ONLY: fastjx_solar2
USE fastjx_photoj_mod,       ONLY: fastjx_photoj
USE submodel_mod,            ONLY: atmos_im
USE nlsizes_namelist_mod,    ONLY:                                    &
    land_field, model_levels, n_cca_lev, ozone_levels, row_length,    &
    rows

USE model_time_mod,          ONLY:                                    &
    i_day_number, i_hour, i_minute, i_month, i_second, i_year,        &
    secs_per_stepim, previous_time
USE ukca_solflux_mod,        ONLY: ukca_solflux

IMPLICIT NONE


! Output photolysis rates
REAL, INTENT(INOUT) :: dj(row_length,rows,model_levels,jppj)
! Pressure on rho levels
REAL, INTENT(IN)    :: p_layer_boundaries(row_length,rows,        &
                                          model_levels+1)
! Pressure on theta levels
REAL, INTENT(IN)    :: p_theta_levels(row_length,rows,model_levels)
! Temperature
REAL, INTENT(IN)    :: t_theta_levels(row_length,rows,model_levels)
! Sulphate aerosol (aitken mode)
REAL, INTENT(IN)    :: f3u_by_omega(row_length,rows)
! = f3_at_u/two_omega. Used as a representation of true latitude.
REAL, INTENT(IN)    :: so4_aitken(row_length,rows,model_levels)
! Sulphate aerosol (accumulation mode)
REAL, INTENT(IN)    :: so4_accum(row_length,rows,model_levels)
! Density * R * R
REAL, INTENT(IN)    :: rho_r2(row_length,rows,model_levels)
! Top of model
REAL, INTENT(IN)    :: z_top_of_model
! liquid water cloud
REAL, INTENT(IN)    :: qcl(row_length,rows,model_levels)
! ice water cloud
REAL, INTENT(IN)    :: qcf(row_length,rows,model_levels)
! Relative humidity fraction
REAL, INTENT(IN)    :: rel_humid_frac(row_length,rows,model_levels)
! cloud area fraction
REAL, INTENT(IN)    :: area_cloud_fraction(row_length,rows,model_levels)
! convective cloud amount
REAL, INTENT(IN)    :: conv_cloud_amount(row_length,rows,model_levels)
! ozone mmr
REAL, INTENT(IN)    :: ozone(row_length,rows,ozone_levels)
! surface pressure
REAL, INTENT(IN)    :: pstar(row_length,rows)
! Convective cloud LWP
REAL, INTENT(IN)    :: conv_cloud_lwp(row_length,rows)
! Surface albedo
REAL, INTENT(IN)    :: surf_albedo(row_length,rows)
! Tropopause level
INTEGER, INTENT(IN) :: tropopause_level(row_length,rows)
! Land fraction
REAL, INTENT(IN)    :: land_fraction(row_length,rows)
! Convective cloud top
INTEGER, INTENT(IN) :: conv_cloud_top(row_length,rows)
! Convective cloud bottom
INTEGER, INTENT(IN) :: conv_cloud_base(row_length,rows)


! Local variables

! Sulphate climatology surface area density
REAL :: so4_sa(row_length,rows,model_levels)
! Sulphate climatology effective radius
REAL :: so4_reff(row_length,rows,model_levels)
! Sulphate climatology mass mixing ratio
REAL :: so4_sa_mmr(row_length,rows,model_levels)

! Parameter to add stratospheric aerosol climatology to OD calculation
! Still not sure of the exact implementation so don't add it to cruntime yet
LOGICAL, PARAMETER  :: L_use_stratclim = .FALSE.

INTEGER :: icode

! Loop variables
INTEGER                       :: i,j,k,l,ia, ix
! Effective diameter
REAL                          :: d_eff
! Sine of true latitude
REAL                          :: sin_true_latitude(row_length, rows)
! Timestep in hours
REAL                          :: timej
! Latitude and Longitude of current point
INTEGER                       :: nslat
INTEGER                       :: nslon
! First time?
LOGICAL                       :: first=.TRUE.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_FASTJX'


! *************************************
! End of Header

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (PrintStatus >= PrStatus_Diag) THEN
  WRITE(umMessage,*) 'UKCA_FASTJX inputs:'
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '1: ',MINVAL(p_layer_boundaries),                    &
                   MAXVAL(p_layer_boundaries),                    &
                   SUM(p_layer_boundaries)/                       &
                   SIZE(p_layer_boundaries)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '2: ',MINVAL(p_theta_levels),MAXVAL(p_theta_levels), &
                   SUM(p_theta_levels)/SIZE(p_theta_levels)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '3: ',MINVAL(t_theta_levels),MAXVAL(t_theta_levels), &
                   SUM(t_theta_levels)/SIZE(t_theta_levels)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '4: ',MINVAL(so4_aitken),MAXVAL(so4_aitken),         &
                   SUM(so4_aitken)/SIZE(so4_aitken)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '5: ',MINVAL(so4_accum),MAXVAL(so4_accum),           &
                   SUM(so4_accum)/SIZE(so4_accum)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '6: ',MINVAL(rho_r2),MAXVAL(rho_r2),                 &
                   SUM(rho_r2)/SIZE(rho_r2)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '7: ',z_top_of_model
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '8: ',MINVAL(qcl),MAXVAL(qcl),SUM(qcl)/SIZE(qcl)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '9: ',MINVAL(qcf),MAXVAL(qcf),SUM(qcf)/SIZE(qcf)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '10:',MINVAL(rel_humid_frac),MAXVAL(rel_humid_frac), &
                   SUM(rel_humid_frac)/SIZE(rel_humid_frac)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '11:',MINVAL(area_cloud_fraction),                   &
                   MAXVAL(area_cloud_fraction),                   &
                   SUM(area_cloud_fraction)/                      &
                   SIZE(area_cloud_fraction)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '12:',MINVAL(conv_cloud_amount),                     &
                   MAXVAL(conv_cloud_amount),                     &
                   SUM(conv_cloud_amount)/SIZE(conv_cloud_amount)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '13:',MINVAL(ozone),MAXVAL(ozone),                   &
                   SUM(ozone)/SIZE(ozone)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '14:',MINVAL(pstar),MAXVAL(pstar),                   &
                   SUM(pstar)/SIZE(pstar)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '15:',MINVAL(conv_cloud_lwp),MAXVAL(conv_cloud_lwp), &
                   SUM(conv_cloud_lwp)/SIZE(conv_cloud_lwp)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '16:',MINVAL(surf_albedo),MAXVAL(surf_albedo),       &
                   SUM(surf_albedo)/SIZE(surf_albedo)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '17:',MINVAL(tropopause_level),                      &
                   MAXVAL(tropopause_level),                      &
                   SUM(tropopause_level)/SIZE(tropopause_level)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '18:',MINVAL(land_fraction),MAXVAL(land_fraction),   &
                   SUM(land_fraction)/SIZE(land_fraction)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '19:',MINVAL(conv_cloud_top),MAXVAL(conv_cloud_top), &
                   SUM(conv_cloud_top)/SIZE(conv_cloud_top)
  CALL umPrint(umMessage,src='ukca_fastjx')
  WRITE(umMessage,*) '20:',MINVAL(conv_cloud_base),                       &
                   MAXVAL(conv_cloud_base),                       &
                   SUM(conv_cloud_base)/SIZE(conv_cloud_base)
  CALL umPrint(umMessage,src='ukca_fastjx')
  icode=0
  CALL umPrintFlush()
END IF

! Initialise
dj(:,:,:,:) = 0.0

! Set Blocking mode:          0) Column-by-column
        !                     1) blocking 1 row
        !                     2) blocking domain
        !                     3) compressed  (not implemented)
        !                     4) load balancing (not implemented)
Blocking%Mode = 2

! Allocate arrays etc.
CALL fastjx_set_limits(row_length,rows,model_levels)
CALL fastjx_allocate_memory

! Read in data from files
IF (first) THEN
  CALL fastjx_inphot
  first=.FALSE.
END IF

! Initialise arrays for levels/units appropriate for fast-j
! Need to update to include aerosols
CALL fastjx_set_arrays

! Set variables concerning model time
! Convert timestep into hours
timej              = secs_per_stepim(atmos_im)/3600.0

  ! Day of the year
daynumber          = i_day_number

  ! Time in hours
tau                = i_hour*1.0+i_minute/60.0+i_second/3600.0      &
                     - timej*0.5

! Calculate solar zenith angles

CALL fastjx_solar2 (timej, sin_true_latitude, true_longitude,   &
                      lcal360)

! Modify solar flux if using solar cycle
IF (i_ukca_solcyc > 0) THEN
  CALL ukca_solflux(previous_time, lcal360,solcyc_spec, fl_cyc)
ELSE
  fl_cyc = 0.
END IF

! Block the data appropriately and call photolysis routines
! Still need to add modes 3 (compressed) and 4 (load balancing)

! Initialise longitude counter to 1
nslon = 1

SELECT CASE (Blocking%Mode)

  ! if blocking point by point
CASE (0)

  ! Loop over rows
  DO j=1,rows

    ! Loop over row, setting longitude counter
    DO i=1,row_length

      ! Set latitude counter to row number
      nslat  = j
      nslon  = i

      ! Loop over row, setting longitude counter
      nsl(1,1)=nslon
      nsl(2,1)=nslat

      ! block using consistent approach
      DO ix=1, kpcx
        sza(ix) = sza_2d(nsl(1,ix),nsl(2,ix))
        szafac(ix) = szafac_2d(nsl(1,ix),nsl(2,ix))

        u0(ix) = sza(ix)*pi_over_180
        u0(ix) = COS(u0(ix))
      END DO

      CALL fastjx_photoj (dj)

    END DO
  END DO

  ! if blocking row by row
CASE (1)

  ! Loop over rows
  DO j=1,rows

    ! Set latitude counter to row number
    nslat  = j

    ! Loop over row, setting longitude counter
    DO i=1,row_length
      nsl(1,i)=nslon+(i-1)
      nsl(2,i)=nslat
    END DO

    ! block using consistent approach
    DO ix=1, kpcx
      sza(ix) = sza_2d(nsl(1,ix),nsl(2,ix))
      szafac(ix) = szafac_2d(nsl(1,ix),nsl(2,ix))

      u0(ix) = sza(ix)*pi_over_180
      u0(ix) = COS(u0(ix))

      sa_block(ix) = sa_2d(nsl(1,ix),nsl(2,ix))
      sa_block(ix) = MIN(1.0, sa_block(ix))
      sa_block(ix) = MAX(0.0, sa_block(ix))

    END DO

    CALL fastjx_photoj (dj)

  END DO

  ! *********************************
  ! If blocking whole domain
CASE (2)

  ! initialise latitude counter to 1
  nslat = 1

  ! Loop over rows
  DO j=1,rows

    ! Loop over columns
    DO i=1,row_length

      ! Calculate positions of longitude and latitude in blocked arrays
      l=i+(j-1)*row_length

      nsl(1,l)=nslon+(i-1)
      nsl(2,l)=nslat+(j-1)
    END DO ! columns
  END DO ! rows

  ! block using consistent approach
  DO ix=1, kpcx
    sza(ix) = sza_2d(nsl(1,ix),nsl(2,ix))
    szafac(ix) = szafac_2d(nsl(1,ix),nsl(2,ix))

    u0(ix) = sza(ix)*pi_over_180
    u0(ix) = COS(u0(ix))

    sa_block(ix) = sa_2d(nsl(1,ix),nsl(2,ix))
    sa_block(ix) = MIN(1.0, sa_block(ix))
    sa_block(ix) = MAX(0.0, sa_block(ix))

  END DO

  CALL fastjx_photoj (dj)


  ! *********************************
  ! Else exit with an error (need to implement blocking types 3 and 4)
CASE DEFAULT

  cmessage='Blocking Mode does not Exist'
  CALL ereport('UKCA_FASTJX', Blocking%Mode,cmessage)

END SELECT

  ! Tidy up at the end
CALL fastjx_deallocate_memory

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

CONTAINS

! ######################################################################
SUBROUTINE fastjx_set_arrays

USE ukca_constants,    ONLY: avogadro, m_air
USE ukca_option_mod,   ONLY: L_ukca_use_background_aerosol

USE ukca_read_aerosol_mod, ONLY: ukca_read_aerosol
USE ukca_read_reff_mod, ONLY: ukca_read_reff
USE ukca_strat_aero_clim_mmr_mod, ONLY: ukca_strat_aero_clim_mmr
IMPLICIT NONE

  ! Loop variables
INTEGER                               :: i,j,k,l, ia

  ! Total mass in column
REAL :: total_mass(1:row_length,1:rows)

  ! Mass column per layer
REAL :: d_mass(1:row_length,1:rows,1:model_levels)

  ! Relative humidity
REAL :: rh(1:row_length,1:rows,1:model_levels)
REAL :: qsat(1:row_length,1:rows,1:model_levels)

  ! Sulphate in accumulation and aitken modes
REAL :: sulph_accu(1:row_length,1:rows,1:model_levels)
REAL :: sulph_aitk(1:row_length,1:rows,1:model_levels)
REAL :: sulphur(1:row_length,1:rows,1:model_levels)

  ! Cloud optical depths
REAL :: odi(1:row_length,1:rows,1:model_levels)
REAL :: odw(1:row_length,1:rows,1:model_levels)
REAL :: ods(1:row_length,1:rows,1:model_levels)

  ! Conversion factor from kg to molecules
REAL                                               :: masfac

  ! Temp variables used to calc aerosol od
INTEGER                                            :: i_humidity
REAL                                               :: delta_humidity
REAL                                               :: f

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_SET_ARRAYS'


  ! *************************************************
  ! EOH
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

  ! Calculate latitude
sin_true_latitude = f3u_by_omega(1:row_length,1:rows)

  ! rz in cm
rz_3d(:,:,1)              =                                     &
    r_Theta_levels(1:row_length,1:rows,0)*100.0
rz_3d(:,:,2:model_levels) =                                     &
    r_rho_levels  (1:row_length,1:rows,2:model_levels)*100.0
rz_3d(:,:,model_levels+1) =                                     &
    r_Theta_levels(1:row_length,1:rows,model_levels)*100.0

  ! calculate pressure at box edges
pz_3d(:,:,1)              = pstar
pz_3d(:,:,2:model_levels) = p_layer_boundaries(:,:,2:model_levels)
pz_3d(:,:,model_levels+1) = 0.0

  ! Calculate mass in box from pressure differences
  ! using hydrostatic approximation
d_mass = (pz_3d(:,:,1:model_levels)                             &
       - pz_3d(:,:,2:model_levels+1))/g

  ! Calculate total mass within convective clouds
total_mass = 0.0
DO k = 1,model_levels
  WHERE ( k <= conv_cloud_top .AND. k >= conv_cloud_base)
    total_mass = total_mass + d_mass(:,:,k)
  END WHERE
END DO


  ! Initialise water vapour to large scale precipitation
odw(:,:,1:model_levels)   = qcl(:,:,1:model_levels)
odi(:,:,1:model_levels)   = qcf(:,:,1:model_levels)

!-----------------------------------------------------------------
! Only add the convective cloud if we are not using the PC2 cloud scheme
IF (i_cld_vn /= i_cld_pc2) THEN

  ! I L-3d_cca (cloud levels?)
  IF (l_3d_cca) THEN
    DO k = 1,n_cca_lev

      WHERE (conv_cloud_top > 0 .AND. t_theta_levels(:,:,k) > tm)
        ! If above freezing point evaluate liquid water
        odw(:,:,k) = odw(:,:,k) + (conv_cloud_lwp(:,:)*           &
                   conv_cloud_amount(:,:,k))/total_mass(:,:)
      ELSEWHERE (conv_cloud_top > 0)
        ! Else evaluate ice water
        odi(:,:,k) = odi(:,:,k) + (conv_cloud_lwp(:,:)*           &
                   conv_cloud_amount(:,:,k))/total_mass(:,:)
      END WHERE

    END DO !  k levels

  ELSE

    ! Else loop over model levels
    DO k=1,model_levels

      WHERE (conv_cloud_top >= k .AND. conv_cloud_base <= k .AND.   &
             t_theta_levels(:,:,k) > tm)
        ! If above freezing point evaluate liquid water
        odw(:,:,k) = odw(:,:,k)+(conv_cloud_lwp(:,:)*               &
                   conv_cloud_amount(:,:,1))/total_mass(:,:)
      ELSEWHERE (conv_cloud_top >= k .AND. conv_cloud_base <= k)
        ! Else evaluate ice water
        odi(:,:,k) = odi(:,:,k)+(conv_cloud_lwp(:,:)*               &
                   conv_cloud_amount(:,:,1))/total_mass(:,:)
      END WHERE

    END DO ! k levels
  END IF ! l_3d_cca

END IF ! i_cld_vn

! Convert mass mixing ratios to column densities.
odw = odw*d_mass
odi = odi*d_mass

!--------------------------------------------------------------
! set effective radii for water drops, different for land&sea
! set effective diameter for ice crystals
! Formulae from John Edwards to convert from column densities to
! optical depths

d_eff=100.0 ! in microns

DO k=1,model_levels

  WHERE (land_fraction > 0.5)
    odw(:,:,k) = odw(:,:,k)*(-8.86964+1.67373e3/6.0)
  ELSEWHERE
    odw(:,:,k) = odw(:,:,k)*(-8.86964+1.67373e3/12.0)
  END WHERE

END DO

odi = odi*(-2.189e-3+3.311e3/d_eff+3.611/d_eff**2)

! Use approach of Briegleb(1992), JGR 97, 7603.
! to account for random overlap of cloud layers
odw(:,:,1:model_levels)    =  odw(:,:,1:model_levels) &
             * (area_cloud_fraction(:,:,1:model_levels))**(1.5)
odi(:,:,1:model_levels)    =  odi(:,:,1:model_levels) &
             * (area_cloud_fraction(:,:,1:model_levels))**(1.5)

! Set optical depth in top layer to be the same as top+1
! could just set to 0. For top level shouldn't make any difference
odw_3d(:,:,1:model_levels)      = odw
odw_3d(:,:,(model_levels+1))    = odw_3d(:,:,(model_levels))
odi_3d(:,:,1:model_levels)      = odi
odi_3d(:,:,(model_levels+1))    = odi_3d(:,:,(model_levels))

! **********************************************************
! Calculate aerosol od columns @600nm
! Note this uses the classic aerosol scheme. Work needs doing to
! incorporate the mode aerosol scheme


! Calculate optical depth from surface aerosol density

! Consistent approach with stratospheric heterogeneous chemistry
IF ( L_use_stratclim) THEN

  CALL ukca_read_aerosol(i_year, i_month, row_length, rows,      &
      model_levels,                                              &
      sin_theta_latitude(1:row_length, 1:rows),                  &
      eta_theta_levels(1:model_levels) * z_top_of_model,         &
      L_ukca_use_background_aerosol,so4_sa)

  CALL ukca_read_reff(i_year, i_month, row_length, rows,         &
      model_levels,                                              &
      sin_theta_latitude(1:row_length, 1:rows),                  &
      eta_theta_levels(1:model_levels) * z_top_of_model,         &
      so4_reff)

  ! Convert surface aerosol density to mmr
  ! Note this routine uses assumptions that cause OD to be
  ! underestimated during Pinatubo period (reff too small)
  CALL ukca_strat_aero_clim_mmr(row_length, rows, model_levels,  &
      r_rho_levels(1:row_length,1:rows,1:model_levels),          &
      rho_r2(1:row_length,1:rows,1:model_levels),                &
      so4_sa(1:row_length,1:rows,1:model_levels),                &
      so4_reff(1:row_length,1:rows,1:model_levels),              &
      so4_sa_mmr(1:row_length,1:rows,1:model_levels))

    ! Add mmr from this process to accumulation mode
  DO i=1, rows
    DO j=1, row_length
      sulph_accu(j,i,(tropopause_level(j,i)):model_levels) =     &
      so4_sa_mmr(j,i,(tropopause_level(j,i)):model_levels)
    END DO
  END DO
ELSE
  ! Use input aerosol concentrations
  sulph_accu(:,:,:) = so4_accum(:,:,:)
  sulph_aitk(:,:,:) = so4_aitken(:,:,:)
END IF

!-----------------------------------------------------------
! Calculate relative humidity

rh = rel_humid_frac !             *1.0E2

! Adopt ukca_fastj approach (see ukca_fast.F90)
! Multiply by molecular weight ratio to convert from mass mixing ratio of
! sulphur atoms to mass mixing ratio of ammonium sulphate.
! Use d_mass to convert to mass of ammonium sulphate.
sulph_accu = sulph_accu*d_mass*m_nh42so4/m_s
sulph_aitk = sulph_aitk*d_mass*m_nh42so4/m_s

! Loop over aerosol types
DO ia=1,sw_spectrum(1)%aerosol%n_aerosol

    ! Select sulfate in aitken and accumulation modes
  IF (sw_spectrum(1)%aerosol%type_aerosol(ia)                     &
        == ip_accum_sulphate .OR.                                 &
      sw_spectrum(1)%aerosol%type_aerosol(ia)                     &
        == ip_aitken_sulphate) THEN

    IF (sw_spectrum(1)%aerosol%i_aerosol_parm(ia) ==              &
           ip_aerosol_param_moist) THEN ! moist aerosol

        ! evaluate size of humidity bins
      delta_humidity = sw_spectrum(1)%aerosol%humidities(2,ia)-   &
                        sw_spectrum(1)%aerosol%humidities(1,ia)

        ! Loop over grid cells calculating optical depths
      DO k=1,model_levels
        DO j = 1, rows
          DO i = 1, row_length

             ! Determine humidity bin of parameterisation
            i_humidity = INT(rh(i,j,k)/delta_humidity)+1

            ! fraction of bin
            f = (rh(i,j,k)-                                       &
              sw_spectrum(1)%aerosol%humidities(i_humidity,ia))/  &
              delta_humidity

             ! Determine optical depths from parameterisation
            ods(i,j,k) =                                          &
            (sw_spectrum(1)%aerosol%ABS(i_humidity,ia,            &
                                               aer_swband)+       &
             sw_spectrum(1)%aerosol%scat(i_humidity,ia,           &
                                               aer_swband))*f+    &
            (sw_spectrum(1)%aerosol%ABS(i_humidity+1,ia,          &
                                               aer_swband)+       &
             sw_spectrum(1)%aerosol%scat(i_humidity+1,ia,         &
                                               aer_swband))*(1-f)
          END DO
        END DO
      END DO

    ELSE ! dry aerosol
      DO k=1,model_levels
        DO j = 1, rows
          DO i = 1, row_length
            ods(i,j,k) =                                          &
              sw_spectrum(1)%aerosol%ABS(1,ia,aer_swband)+        &
              sw_spectrum(1)%aerosol%scat(1,ia,aer_swband)
          END DO
        END DO
      END DO
    END IF

    IF (sw_spectrum(1)%aerosol%type_aerosol(ia)                   &
      == ip_accum_sulphate) sulph_accu = sulph_accu*ods
    IF (sw_spectrum(1)%aerosol%type_aerosol(ia)                   &
      == ip_aitken_sulphate) sulph_aitk = sulph_aitk*ods
  END IF   ! type_aerosol_sw
END DO  ! ia

! Sum the aitkin and accumulation type optical depths
sulphur = sulph_aitk + sulph_accu

! Copy sulphate od to global array
ods_3d(:,:,1:model_levels)  = sulphur(:,:,1:model_levels)
ods_3d(:,:,(model_levels+1))= sulphur(:,:,model_levels)

IF (printstatus >= prstatus_diag) THEN
  DO i=1, model_levels
    WRITE(umMessage,*) 'fastjX converter max ', i,  MAXVAL(so4_sa_mmr(:,:,i))
    CALL umPrint(umMessage,src='ukca_fastjx')
  END DO
END IF

! ********************************************************
! Set other variables here

! Set surface albedo
sa_2d                           = surf_albedo

! Calculate pressure at box edges (include box to TOA)
! NB convert Pa to hPa for fastj routines
pz_all(:,:,1:(model_levels+1))  = p_layer_boundaries(:,:,:)/100.0
pz_all(:,:,(model_levels+2))    = 0.0e0

! Calculate heights of box edges
rz_all(:,:,1:model_levels+1)    = rz_3d(:,:,1:model_levels+1)
rz_all(:,:,model_levels+2)      = rz_all(:,:,model_levels+1)      &
                                  + zzht

! Set temperature to be temperature
! Use top temperature for top+1 level
tz_3d(:,:,1:model_levels)       = t_theta_levels(:,:,:)
tz_3d(:,:,(model_levels+1))     = t_theta_levels(:,:,model_levels)

! Copy air mass from local array
dm_3d(:,:,1:model_levels)       = d_mass(:,:,1:model_levels)
! Evaluate top box from difference with TOA p i.e. 0!
dm_3d(:,:,(model_levels+1))     = pz_3d(:,:,model_levels)/g

! Convert air mass from kg/m2 to molecules/cm2
masfac                          = avogadro/(m_air*1.0e4)
dm_3d(:,:,:)                    = dm_3d(:,:,:)*masfac

! Convert ozone mmr to vmr and write to fj_ozone array
o3_3d(:,:,1:model_levels)       = ozone/c_o3
o3_3d(:,:,(model_levels+1))     = o3_3d(:,:,model_levels)

! Multiply by molecules per m2 to get it in units fastj wants
o3_3d                           = o3_3d*dm_3d

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_set_arrays

!#######################################################################

END SUBROUTINE ukca_fastjx

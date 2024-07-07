! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Description:
!  Main driver routine for chemistry
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_MAIN1.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_chemistry_ctl_col_mod

IMPLICIT NONE 

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName=              &
                                     'UKCA_CHEMISTRY_CTL_COL_MOD'

CONTAINS

SUBROUTINE ukca_chemistry_ctl_col(i_month, i_day_number, i_hour, &
                r_minute, secs_per_step,                         &
                ntracers,                                        &
                sinlat,                                          &
                coslat,                                          &
                true_longitude,                                  &
                pres, temp, q,                                   &
                qcf, qcl, rh,                                    &
                p_layer_boundaries,                              &
                r_theta_levels,                                  &
                z_top_of_model,                                  &
                cos_zenith_angle,                                &
                tracer,                                          &
                all_ntp,                                         &
                t_surf, dzl, z0m, u_s,                           &
                drain, crain,                                    &
                cloud_frac,                                      &
                fastj_dj,                                        &
                volume, mass,                                    &
                land_points, land_index,                         &
                tile_pts, tile_index, tile_frac,                 &
                zbl, surf_hf, seaice_frac, stcon,                &
                soilmc_lp, fland, laift_lp, canhtft_lp,          &
                z0tile_lp, t0tile_lp, canwctile_lp,              &
                have_nat,                                        &
                pv_at_theta,                                     &
                theta,                                           &
                um_ozone3d,                                      &
                uph2so4inaer,                                    &
                delso2_wet_h2o2,                                 &
                delso2_wet_o3,                                   &
                delh2so4_chem,                                   &
                delso2_drydep,                                   &
                delso2_wetdep,                                   &
                so4_sa,                                          &
                nat_psc,                                         &
                atm_ch4_mol,                                     &
                atm_co_mol,                                      &
                atm_n2o_mol,                                     &
                atm_cf2cl2_mol,                                  &
                atm_cfcl3_mol,                                   &
                atm_mebr_mol,                                    &
                atm_h2_mol                                       &
                )

USE conversions_mod, ONLY: pi_over_180, pi
USE jules_surface_types_mod, ONLY: ntype, npft
USE asad_mod
USE asad_chem_flux_diags
USE ukca_d1_defs
USE ukca_cspecies
USE asad_findreaction_mod, ONLY: asad_findreaction
USE UKCA_tropopause
USE UKCA_strat_update
USE ukca_chem_offline,    ONLY: h2o2_offline
USE ukca_constants,       ONLY: c_o3, c_h2o, c_hono2, c_o1d,     &
                                c_co2, avogadro, boltzmann
USE UKCA_dissoc
USE UKCA_phot2d,          ONLY: ukca_photin, ukca_curve,         &
                                ukca_inpr2d, nolev, nlphot,      &
                                ntphot, pjin, pr2d

USE ukca_option_mod,      ONLY: L_ukca_raq, L_ukca_trop,         &
                                L_ukca_trophet, L_ukca_use_2dtop,&
                                L_ukca_chem,                     &
                                L_ukca_het_psc, L_ukca_nr_aqchem,&
                                L_ukca_advh2o, L_ukca_aerchem,   &
                                l_ukca_raqaero,                  &
                                l_ukca_offline, l_ukca_intdd,    &
                                i_ukca_photol, fastjx_prescutoff
USE ukca_photo_scheme_mod, ONLY: i_ukca_phot2d, i_ukca_fastjx
USE ukca_ntp_mod,       ONLY: ntp_type, dim_ntp, name2ntpindex
USE carbon_options_mod, ONLY: l_co2_interactive
USE ukca_um_interf_mod, ONLY: co2_interactive

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE UM_ParVars
USE missing_data_mod,      ONLY: rmdi, imdi

USE nlsizes_namelist_mod, ONLY:                              &
    bl_levels, global_row_length, global_rows, model_levels, &
    row_length, rows, theta_field_size

USE errormessagelength_mod, ONLY: errormessagelength

USE asad_cdrive_mod, ONLY: asad_cdrive
USE ukca_ddepctl_mod, ONLY: ukca_ddepctl
USE ukca_ddeprt_mod, ONLY: ukca_ddeprt
USE ukca_fracdiss_mod, ONLY: ukca_fracdiss
USE ukca_sediment_mod, ONLY: ukca_sediment
USE ukca_stratf_mod, ONLY: ukca_stratf
USE ukca_wdeprt_mod, ONLY: ukca_wdeprt
IMPLICIT NONE


INTEGER, INTENT(IN) :: ntracers          ! no. of tracers
INTEGER, INTENT(IN) :: i_month           ! month
INTEGER, INTENT(IN) :: i_day_number      ! day
INTEGER, INTENT(IN) :: i_hour            ! hour
INTEGER, INTENT(IN) :: uph2so4inaer      ! flag for H2SO4 updating

!       Variables for interactive dry deposition scheme
INTEGER, INTENT(IN) :: land_points
INTEGER, INTENT(IN) :: land_index(land_points)
INTEGER, INTENT(IN) :: tile_pts(ntype)
INTEGER, INTENT(IN) :: tile_index(land_points,ntype)

REAL, INTENT(IN) :: r_minute                           ! minute
REAL, INTENT(IN) :: secs_per_step                      ! time step
REAL, INTENT(IN) :: z_top_of_model                     ! top of model (m)
REAL, INTENT(IN) :: sinlat(row_length, rows)           ! sin(latitude)
REAL, INTENT(IN) :: coslat(row_length, rows)           ! cos(latitude)
REAL, INTENT(IN) :: true_longitude(row_length,rows)    ! longitude
REAL, INTENT(IN) :: pres(row_length,rows,model_levels) ! pressure
REAL, INTENT(IN) :: p_layer_boundaries(row_length,rows,    &
                                       0:model_levels) ! pressure
REAL, INTENT(IN) :: r_theta_levels(row_length,rows,             &
                                   0:model_levels)
REAL, INTENT(IN) :: temp(row_length,rows,model_levels) ! actual temp
REAL, INTENT(IN) :: dzl(row_length, rows, bl_levels)   ! thickness
REAL, INTENT(IN) :: u_s(row_length, rows)              ! ustar
REAL, INTENT(IN) :: z0m(row_length, rows)              ! roughness
REAL, INTENT(IN) :: t_surf(row_length, rows)           ! surface temp
REAL, INTENT(IN) :: drain(row_length,rows,model_levels) ! 3-D LS rain
REAL, INTENT(IN) :: crain(row_length,rows,model_levels) ! 3-D convec
REAL, INTENT(IN) :: cos_zenith_angle(row_length, rows)  ! cosine of ZA
REAL, INTENT(IN) :: volume(row_length,rows,model_levels) ! cell vol.
REAL, INTENT(IN) :: mass(row_length, rows, model_levels) ! cell mass

REAL, INTENT(IN) :: pv_at_theta(row_length, rows, model_levels) ! PV
REAL, INTENT(IN) :: theta(row_length, rows, model_levels)! theta
REAL, INTENT(IN) :: um_ozone3d(row_length, rows, model_levels) ! O3
REAL, INTENT(IN) :: qcf(row_length, rows, model_levels)  ! qcf
REAL, INTENT(IN) :: qcl(row_length, rows, model_levels)  ! qcl
REAL, INTENT(IN) :: rh(row_length, rows, model_levels)   ! RH frac
REAL, INTENT(IN) :: cloud_frac(row_length, rows, model_levels)
REAL, INTENT(IN) :: so4_sa(row_length,rows,model_levels) ! aerosol surface area

!       Variables for interactive dry deposition scheme

REAL, INTENT(IN) :: tile_frac(land_points,ntype)
REAL, INTENT(IN) :: zbl(row_length,rows)
REAL, INTENT(IN) :: surf_hf(row_length,rows)
REAL, INTENT(IN) :: seaice_frac(row_length,rows)
REAL, INTENT(IN) :: stcon(row_length,rows,npft)
REAL, INTENT(IN) :: soilmc_lp(land_points)
REAL, INTENT(IN) :: fland(land_points)
REAL, INTENT(IN) :: laift_lp(land_points,npft)
REAL, INTENT(IN) :: canhtft_lp(land_points,npft)
REAL, INTENT(IN) :: z0tile_lp(land_points,ntype)
REAL, INTENT(IN) :: t0tile_lp(land_points,ntype)
REAL, INTENT(IN) :: canwctile_lp(land_points,ntype)

! Mask to limit formation of Nat below specified height
LOGICAL, INTENT(IN) :: have_nat(row_length, rows, model_levels)

REAL, INTENT(INOUT) :: fastj_dj(row_length,rows,model_levels,   &
                                jppj)
REAL, INTENT(INOUT) :: q(row_length,rows,model_levels)   ! water vapour
REAL, INTENT(INOUT) :: tracer(row_length,rows,                  &
                              model_levels,ntracers)     ! tracer MMR

! Non transported prognostics
TYPE(ntp_type), INTENT(INOUT) :: all_ntp(dim_ntp)

! SO2 increments
REAL, INTENT(INOUT) :: delSO2_wet_H2O2(row_length,rows,           &
                                       model_levels)
REAL, INTENT(INOUT) :: delSO2_wet_O3(row_length,rows,model_levels)
REAL, INTENT(INOUT) :: delh2so4_chem(row_length,rows,model_levels)
REAL, INTENT(INOUT) :: delSO2_drydep(row_length,rows,model_levels)
REAL, INTENT(INOUT) :: delSO2_wetdep(row_length,rows,model_levels)

! Nitric acid trihydrate (kg(nat)/kg(air))
REAL, INTENT(INOUT) :: nat_psc(row_length,rows,model_levels)

! Atmospheric Burden of CH4
REAL, INTENT(INOUT) :: atm_ch4_mol(row_length,rows,model_levels)

! Atmospheric Burden of CO
REAL, INTENT(INOUT) :: atm_co_mol(row_length,rows,model_levels)

! Atmospheric Burden of Nitrous Oxide (N2O)
REAL, INTENT(INOUT) :: atm_n2o_mol(row_length,rows,model_levels)

! Atmospheric Burden of CFC-12
REAL, INTENT(INOUT) :: atm_cf2cl2_mol(row_length,rows,model_levels)

! Atmospheric Burden of CFC-11
REAL, INTENT(INOUT) :: atm_cfcl3_mol(row_length,rows,model_levels)

! Atmospheric Burden of CH3Br
REAL, INTENT(INOUT) :: atm_mebr_mol(row_length,rows,model_levels)

! Atmospheric Burden of H2
REAL, INTENT(INOUT) :: atm_h2_mol(row_length,rows,model_levels)

! Local variables
INTEGER :: nlev_with_ddep(row_length, rows)     ! No levs in bl
INTEGER :: nlev_with_ddep2(model_levels)    ! No levs in bl

INTEGER, SAVE :: first_row
INTEGER, SAVE :: first_column

INTEGER :: i             ! Loop variable
INTEGER :: j             ! loop variable
INTEGER :: js            ! loop variable
INTEGER :: k             ! loop variable
INTEGER :: klevel        ! dummy variable
INTEGER :: l             ! loop variable
INTEGER :: ll            ! loop variable
INTEGER :: m             ! loop variable
INTEGER :: n             ! loop variable
INTEGER :: sp            ! loop variable
INTEGER :: n_pnts        ! no. of pts in 2D passed to CDRIVE

INTEGER           :: ierr                     ! Error code: asad diags routines
INTEGER           :: errcode                  ! Error code: ereport
CHARACTER(LEN=errormessagelength) :: cmessage         ! Error message
CHARACTER(LEN=10)      :: prods(2)                 ! Products
LOGICAL           :: ddmask(model_levels) ! mask

REAL, PARAMETER :: fxb = 23.45 * pi_over_180 ! tropic of capricorn
REAL, PARAMETER :: fxc = 24.0/ pi

!       Pressure level above which top boundary conditions are applied
REAL, PARAMETER :: p_above = 7000.0           ! Pa

REAL :: tgmt               ! GMT time (decimal represent
REAL :: declin             ! declination
REAL :: total_water        ! Total water
REAL :: const              ! constant
REAL, PARAMETER :: limit = 1200.0   ! seconds. For timesteps
                                   ! greater than this we halve
                                   ! the chemical timestep
                                   ! (for IMPACT).

! SO2 increments in molecules/cm^3
REAL :: SO2_wetox_H2O2(model_levels)
REAL :: SO2_wetox_O3(model_levels)
REAL :: SO2_dryox_OH(model_levels)

! array to store H2SO4 when updated in MODE
REAL, ALLOCATABLE :: ystore(:)
REAL :: zftr(model_levels,jpctr)  ! 1-D array of tracers
REAL :: zp  (model_levels)        ! 1-D pressure
REAL :: zt  (model_levels)        ! 1-D temperature
REAL :: zclw(model_levels)        ! 1-D cloud liquid water
REAL :: zfcloud(model_levels)     ! 1-D cloud fraction
REAL :: cdot(model_levels,jpctr)  ! 1-D chem. tendency
REAL :: zq(model_levels)          ! 1-D water vapour vmr
REAL :: co2_1d(model_levels)      ! 1-D CO2 vmr
REAL :: pjinda(rows, ntphot, jppj)    ! PJIN at one level
REAL :: zprt(row_length, rows, jppj)  ! 2-D photolysis rates
! 2-D photolysis rates from strat_photol
REAL :: zprt_strat(row_length, rows, jppj)  
! 3-D (from 2-D) photolysis rates
REAL :: zprt3(row_length, rows, model_levels, jppj)
REAL :: zprt1d(model_levels,jppj) ! 1-D photolysis rates
REAL :: tloc  (row_length, rows)      ! local time
REAL :: daylen(row_length, rows)      ! local daylength
REAL :: cs_hour_ang(row_length, rows) ! cosine hour angle
REAL :: tanlat(row_length, rows)      ! tangens of latitude
REAL :: zdryrt(row_length, rows, jpdd)                ! dry dep rate
REAL :: zdryrt2(model_levels, jpdd)               ! dry dep rate
REAL :: zwetrt(row_length, rows, model_levels, jpdw)  ! wet dep rate
REAL :: zwetrt_theta(theta_field_size, model_levels, jpdw)  ! wet dep rate
REAL :: zwetrt2(model_levels, jpdw)               ! wet dep rat
REAL :: zfrdiss2(model_levels,jpdw,jpeq+1)        ! dissolved fraction
REAL :: zfrdiss(row_length, rows, model_levels, jpdw, jpeq+1)
REAL :: rc_het(model_levels,2)                ! heterog rates for trop chem
REAL :: kp_nh(row_length, rows, model_levels)     ! Dissociation const
REAL :: kp_nh2(model_levels)                  ! Dissociation const
REAL :: ozonecol(row_length, rows, model_levels)  ! for strat chem

! Local stratospheric CH4 loss rate
REAL :: strat_ch4loss_2d(model_levels,model_levels) 
REAL :: pr2dj(nlphot)                             ! 2D photolysis level pres

REAL, SAVE :: first_lat
REAL, SAVE :: dellat
REAL, SAVE :: first_lon
REAL, SAVE :: dellon

LOGICAL, SAVE :: firstcall = .TRUE.


! The calls to ukca_conserve require a logical to be set.
! ukca_conserve calculates and conserves total chlorine, bromine, and
! hydrogen. For these elements closed chemistry should be prescribed.
! Called before chemistry, before_chem, it calculates
! total bromine, chlorine, and hydrogen as 3-D fields. Called afer
! chemistry, after_chem, it rescales the chlorine, bromine
! and hydrogen containing compounds so that total chlorine, bromine
! and hydrogen are conserved under chemistry.
LOGICAL, PARAMETER :: before_chem = .TRUE.
LOGICAL, PARAMETER :: after_chem = .FALSE.

! Variables for heterogeneous chemistry
REAL, ALLOCATABLE :: shno3_3d(:,:,:)
LOGICAL :: stratflag(model_levels)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEMISTRY_CTL_COL'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
n_pnts = model_levels

! dummy variable for compatability with theta_field call
klevel=0

IF (firstcall) THEN

  ! Determine where the domain is in latitude
  first_lat = ASIN(sinlat(1,1))
  dellat = ASIN(sinlat(1,2)) - first_lat
  first_row = INT((first_lat + 0.5*pi) / dellat + 0.00001) + 1

  ! Determine where the domain is in longitude
  first_lon = true_longitude(1,2)
  dellon = true_longitude(2,2) - first_lon
  first_column = INT(first_lon / dellon + 0.00001) + 1

  !         Check whether water vapour is advective tracer. Then,
  !         check whether UM and ASAD advective tracers correspond
  !         to each other.

  IF ((L_ukca_advh2o) .AND. (n_h2o == 0)) THEN
    cmessage='No tracer for advected water vapour'
    errcode = 4
    CALL ereport('UKCA_CHEMISTRY_CTL_COL',errcode,cmessage)
  END IF

  ! Allocate stratospheric photolysis arrays
  IF (L_ukca_strat .OR. L_ukca_stratcfc .OR. L_ukca_strattrop)  &
    CALL ukca_strat_photol_init()

  ! Identify the SO2+OH rate coeff, the products are alternatives depending
  !  on whether the H2SO4 tracer updating is to be done in ASAD or in MODE.
  iso2_oh = 0
  ih2so4_hv = 0
  IF (l_ukca_nr_aqchem) THEN
    IF (l_ukca_offline) THEN
      prods = (/'H2SO4     ','          '/)
    ELSE
      prods = (/'H2SO4     ','HO2       '/)
    END IF
    iso2_oh = asad_findreaction( 'SO2       ', 'OH        ',      &
                           prods, 2, spt, ntrkx, jptk+1, jpspt )

    IF (iso2_oh == 0) THEN   ! check for stratospheric sulphur chemistry
      prods = (/'SO3       ','H2O       '/)
      iso2_oh = asad_findreaction( 'SO2       ', 'OH        ',    &
                           prods, 2, spt, ntrkx, jptk+1, jpspt )
      prods = (/'SO3       ','OH        '/)
      ih2so4_hv = asad_findreaction( 'H2SO4     ', 'PHOTON    ',  &
                           prods, 2, spj, nprkx, jppj+1, jpspj )

    END IF

    IF (iso2_oh == 0 .AND. ih2so4_hv == 0) THEN
      cmessage=' Sulphur chemistry reactions not found'
      WRITE(umMessage,'(A)') cmessage
      CALL umPrint(umMessage,src='ukca_chemistry_ctl_col')
      WRITE(umMessage,'(A,I0,A,I0)') 'iso2_oh: ',iso2_oh,             &
                                     ' ih2so4_hv: ',ih2so4_hv
      CALL umPrint(umMessage,src='ukca_chemistry_ctl_col')
      errcode = 1
      CALL ereport('UKCA_CHEMISTRY_CTL_COL',errcode,cmessage)
    END IF
  END IF   ! L_ukca_achem

END IF  ! of initialization of chemistry subroutine (firstcall)


!       Calculate local time as function of longitude

tgmt = REAL(i_hour) + r_minute/60.0                             &
                    + secs_per_step * 0.5 / 3600.0
!        IF (tgmt < 0.) tgmt = tgmt + 24.

tloc = tgmt + 24.0 * true_longitude/pi/2.0
WHERE (tloc > 24.0) tloc = tloc - 24.0

!       Calculate Declination Angle and Daylength for each row for
!       current day of the year.
!       Ensure COS of HOUR ANGLE does not exceed + or - 1, and set DAY
!       LENGTH of 1st & last rows to be the same as that of adjacent rows
!       to avoid possible problems at the poles (tan(90)=infinity).

daylen = 0.0
DO i=1,rows
  DO j=1,row_length
    IF (ABS(coslat(j,i)) < 1e-10) THEN
      IF (sinlat(j,i) >= 0.0) THEN
        tanlat(j,i) = 1.0e20
      ELSE
        tanlat(j,i) = -1.0e20
      END IF
    ELSE
      tanlat(j,i)=sinlat(j,i)/coslat(j,i)
    END IF
  END DO
END DO

declin = fxb * SIN(pi_over_180*(266.0+i_day_number))
cs_hour_ang = -tanlat * TAN(declin)
WHERE (cs_hour_ang < -1.0) cs_hour_ang = -1.0
WHERE (cs_hour_ang >  1.0) cs_hour_ang =  1.0
daylen = fxc * ACOS(cs_hour_ang)

! Call routine to calculate dry deposition rates.
zdryrt  = 0.0
zdryrt2 = 0.0
IF (ndepd /= 0) THEN

  IF (L_ukca_intdd) THEN           ! Call interactive dry dep

    CALL ukca_ddepctl(row_length, rows, bl_levels,             &
      land_points, land_index, tile_pts, tile_index,           &
      secs_per_step, sinlat, tile_frac, t_surf,                &
      p_layer_boundaries(:,:,0), dzl, zbl, surf_hf, u_s,       &
      rh, stcon, soilmc_lp, fland, seaice_frac, laift_lp,      &
      canhtft_lp, z0tile_lp, t0tile_lp, canwctile_lp,          &
      nlev_with_ddep, zdryrt)

  ELSE                             ! Call prescribed dry dep

     CALL ukca_ddeprt(daylen, tloc, theta_field_size, dzl,     &
                             bl_levels,                        &
                             z0m, u_s, t_surf,                 &
                             sinlat, i_month,                  &
                             1, theta_field_size,              &
                             zdryrt)

  END IF
END IF

!       Call routine to calculate wet deposition rates.

zwetrt  = 0.0
zwetrt2 = 0.0
zwetrt_theta  = 0.0
IF (ndepw /= 0) THEN

   CALL ukca_wdeprt(RESHAPE(drain,(/theta_field_size,model_levels/)), &
        RESHAPE(crain,(/theta_field_size,model_levels/)),             &
        theta_field_size,                                             &
        RESHAPE(temp,(/theta_field_size,model_levels/)),              &
        RESHAPE(sinlat,(/theta_field_size/)), secs_per_step,          &
        1,theta_field_size, zwetrt_theta)
   zwetrt=RESHAPE(zwetrt_theta,(/row_length, rows, model_levels, jpdw/))
END IF

! Calculate dissolved fraction
IF (L_ukca_aerchem .OR. l_ukca_raqaero) THEN
  CALL ukca_fracdiss(row_length, rows,                            &
                   temp, pres, rh, qcl, zfrdiss, kp_nh)
END IF

!       Call routine to read in 2D photolysis rates once per day and
!       interpolate to model latitude and levels.


IF ( (i_ukca_photol == i_ukca_phot2d) .AND.                     &
     ((i_hour == 0 .AND. r_minute == 0.0) .OR. firstcall) ) THEN
  CALL ukca_photin(i_day_number, row_length, rows,              &
                   model_levels, theta_field_size,              &
                   first_row, global_row_length, global_rows,   &
                   RESHAPE(sinlat,(/theta_field_size/)),        &
                   pres, jppj)
ELSE IF (i_ukca_photol == i_ukca_fastjx) THEN
  CALL ukca_inpr2d(pr2d,pr2dj)  ! for stratosphere only
END IF

!       Calculate ozone column for stratospheric photolysis

IF (L_ukca_strat .OR. L_ukca_stratcfc .OR. L_ukca_strattrop)    &
  THEN
  IF (n_o3 > 0) THEN
!   ukca_calc_ozonecol calculates ozone column above gridbox in 
!   molecules/cm2 given ozone concentration in VMR
    CALL ukca_calc_ozonecol(model_levels, rows, row_length,     &
                   z_top_of_model, p_layer_boundaries, pres,    &
                   tracer(:,:,:,n_o3)/c_o3,                     &
                   ozonecol)
  END IF

  ! Calculate total chlorine and total bromine before chemistry

  IF (nn_cl > 0) THEN  
    CALL ukca_conserve(row_length, rows, model_levels, ntracers,&
         tracer, pres, drain, crain, before_chem)
  END IF
END IF    ! L_ukca_strat etc

! if heterogeneous chemistry is selected, allocate solid HNO3 array
IF (L_ukca_het_psc) THEN
  IF (.NOT. ALLOCATED(shno3_3d)) &
       ALLOCATE(shno3_3d(row_length, rows, model_levels))
  shno3_3d = 0.0
END IF

!       Initialize budget variables
strat_ch4loss_2d = 0.0

! Reduce over-prediction of H2O2 using ancillary value.
IF (L_ukca_offline) THEN
   WHERE (tracer(:,:,:,n_h2o2) > h2o2_offline(:,:,:))            &
        tracer(:,:,:,n_h2o2) = h2o2_offline(:,:,:)
END IF

zprt3 = 0.0
DO k=1,model_levels
  ! Calculate photolysis.

  ! reset zprt to zero here - it is used to fill zprt3 array so is in effect
  ! only temporary. This is required as it is NOT being used as we go through
  ! this loop
  zprt = 0.0

  ! Calculate TOMCAT-heritage photolysis if PHOT2D is selected.
  IF (i_ukca_photol == i_ukca_phot2d) THEN

    !         Interpolate 2D photolysis rates as function of longitude
    !         per model step.

    pjinda(:,:,:) = pjin(1:rows,k,:,:)
    CALL ukca_curve(pjinda, RESHAPE(tloc,(/theta_field_size/)),       &
         RESHAPE(daylen,(/theta_field_size/)),theta_field_size,       &
         model_levels, rows, row_length, zprt1d)
    zprt = RESHAPE(zprt1d,(/row_length, rows, jppj/))
  END IF

  ! Calculate SLIMCAT-heritage photolysis if PHOT2D is selected or for
  ! middle-atmosphere chemistry in the mesosphere. STRAT_PHOTOL merges
  ! the SLIMCAT rates with TOMCAT rates as calculated above.
  ! If FASTJX is selected, only photolysis for wavelengths less than 177 nm
  ! is calculated, and the result is added to the FAST-JX rates. This
  ! only matters at heights > 60 km as below the short wavelengths can be
  ! ignored.
  
  IF (L_ukca_strat .OR. L_ukca_strattrop .OR. L_ukca_stratcfc) THEN

    IF (i_ukca_photol == i_ukca_fastjx) THEN

      ! if any pressure point on domain is below cut-off value 
      ! call strat_photol. Be careful here as need rates to bit-compare
      ! when changing the size of the domain, so only copy columns when
      ! required
      zprt_strat = 0.0
      IF (MINVAL(pres(:,:,k)) < fastjx_prescutoff .AND.           &
          fastjx_mode /= 3) THEN
        CALL ukca_strat_photol(pres(:,:,k), temp(:,:,k),          &
             ozonecol(:,:,k), cos_zenith_angle, secs_per_step,    &
             zprt_strat)
      END IF

      ! only overwrite zprt on columns where pressure is below minimum
      DO i=1,rows
         DO j=1,row_length
            IF (pres(j,i,k) < fastjx_prescutoff) THEN
               ! FILL ZPRT DEPENDING ON FASTJX_MODE
               IF (fastjx_mode == 1) THEN
                  ! (fastjx_mode == 1): from strat_photol
                  zprt(j,i,:) = zprt_strat(j,i,:)
               ELSE IF (fastjx_mode == 2) THEN
                  ! (fastjx_mode == 2): from strat_photol & FJX
                  zprt(j,i,:) = zprt_strat(j,i,:) + fastj_dj(j,i,k,:)
               ELSE IF (fastjx_mode == 3) THEN
                  ! (fastjx_mode == 3): from FJX
                  zprt(j,i,:) = fastj_dj(j,i,k,:)
               ELSE
                  ! no other options should be available
                  WRITE(cmessage,'(A,A,I0)') 'Incorrect option for ',   &
                       'fastjx_mode: ',fastjx_mode
                  errcode = 3
                  CALL ereport('UKCA_CHEMISTRY_CTL_COL',errcode,        &
                               cmessage)
               END IF 
            ELSE
               ! only take FJX here
               zprt(j,i,:) = fastj_dj(j,i,k,:)
            END IF
         END DO
      END DO
      
    ELSE
      ! No FAST-JX here.
      CALL ukca_strat_photol(pres(:,:,k), temp(:,:,k),            &
           ozonecol(:,:,k), cos_zenith_angle, secs_per_step,      &
           zprt)
    END IF

  ELSE

    ! tropospheric chemistry selected here. Use FAST-JX rates or previously
    ! calculated 2-D rates.
    ! If using FastJX, over-write zprt with FJX values
    IF (i_ukca_photol == i_ukca_fastjx) THEN
      zprt = fastj_dj(:,:,k,:)
    END IF
 END IF

 ! make 3D zprt array for use when calling ASAD
 zprt3(:,:,k,:) = zprt(:,:,:)
END DO


! Model levels loop
!$OMP PARALLEL DEFAULT(NONE)                                               &
!$OMP PRIVATE(cdot, cmessage, errcode, i, ierr,                            &
!$OMP         j, js, l, rc_het, stratflag,                                 &
!$OMP         ystore, zclw, zdryrt2, zfcloud, zftr,                        &
!$OMP         zp, zprt1d, zq, zt, co2_1d, zwetrt2)                         &
!$OMP SHARED(all_ntp, atm_cf2cl2_mol, atm_cfcl3_mol, atm_ch4_mol,          &
!$OMP        atm_co_mol, atm_h2_mol, atm_mebr_mol, atm_n2o_mol,            &
!$OMP        c_species, cloud_frac,                                        &
!$OMP        delh2so4_chem, delSO2_wet_H2O2, delSO2_wet_O3, have_nat,      &
!$OMP        ih2so4_hv, ihso3_h2o2, ihso3_o3, iso2_oh, iso3_o3,            &
!$OMP        jpctr, jpdd,  klevel, speci, co2_interactive,                 &
!$OMP        l_asad_use_chem_diags, l_asad_use_drydep, l_co2_interactive,  &
!$OMP        l_asad_use_flux_rxns, l_asad_use_psc_diagnostic,              &
!$OMP        l_asad_use_rxn_rates, l_asad_use_wetdep, l_troposphere,       &
!$OMP        l_ukca_aerchem, l_ukca_chem, l_ukca_het_psc, l_ukca_intdd,    &
!$OMP        l_ukca_nr_aqchem, l_ukca_offline, l_ukca_raq,                 &
!$OMP        l_ukca_raqaero, l_ukca_trop, l_ukca_trophet,                  &
!$OMP        model_levels, n_cf2cl2, n_cfcl3, n_ch4, n_co, n_h2,           &
!$OMP        n_mebr, n_n2o, n_pnts, nlev_with_ddep,                        &
!$OMP        nn_h2o2, nn_h2so4, nn_o1d, nn_o3, nn_o3p, nn_oh, nn_so2,      &
!$OMP        o1d_in_ss, o3p_in_ss, pres, q, qcf, qcl,                      &
!$OMP        row_length, rows, secs_per_step, so4_sa, shno3_3d,            &
!$OMP        temp, tracer, uph2so4inaer, volume, zdryrt, zprt3, zwetrt)

IF (.NOT. ALLOCATED(ystore) .AND. uph2so4inaer == 1)            &
                             ALLOCATE(ystore(model_levels))

!$OMP DO SCHEDULE(STATIC)
DO i=1,rows
  DO j=1,row_length

    ! Copy water vapour and ice field into 1-D arrays
    sph2o(:) = 0.0
    IF (L_ukca_het_psc) THEN
      sph2o(:) = qcf(j,i,:)/c_h2o
    END IF

    zdryrt2 = 0.0
    IF (L_ukca_intdd) THEN
      ! Interactive scheme extracts from levels in boundary layer
      DO l=1,jpdd
        zdryrt2(1:nlev_with_ddep(j,i),l) = zdryrt(j,i,l)
      END DO
    ELSE    ! non-interactive
      zdryrt2(1,:) = zdryrt(j,i,:)
    END IF
   
    !       Put pressure, temperature and tracer mmr into 1-D arrays
    !       for use in ASAD chemical solver
   
    zp(:) = pres(j,i,:)
    zt(:) = temp(j,i,:)
    zq(:) = q(j,i,:)/c_h2o

   IF (ANY(speci(:) == 'CO2       ')) THEN
!  Copy the CO2 concentration into the asad module as VMR
     IF (l_co2_interactive) THEN
       co2_1d(:) = co2_interactive(j,i,:)/c_co2
     ELSE
       co2_1d(:) = rmdi
     END IF

   END IF
   
    zclw(:) = 0.0
    zfcloud(:) = 0.0
    IF (L_UKCA_aerchem) THEN
      zclw(:) = qcl(j,i,:)
      zfcloud(:) = cloud_frac(j,i,:)
    END IF
   

    ! Convert mmr into vmr for tracers
    DO js=1,jpctr
      zftr(:,js) = tracer(j,i,:,js)/c_species(js)
    END DO

    ! Map photolysis rates onto 1-D array.
    IF (l_ukca_offline) THEN
      ! Offline chemistry has no photolysis
      zprt1d(:,:) = 0.0
    ELSE
      ! take column+species info from 4D photolysis array
      zprt1d(:,:) = zprt3(j,i,:,:)
    END IF

    !       Call ASAD routines to do chemistry integration
    !       In lowest levels choose half the dynamical timestep for
    !       chemistry. If dynamical timestep > 20 min, use half and
    !       quarter of dynamical timestep for chemistry.

    IF (.NOT. (L_ukca_trop .OR. L_ukca_aerchem .OR.                  &
               L_ukca_raq .OR. l_ukca_raqaero)) THEN       ! Not B-E

      ! retrieve tropospheric heterogeneous rates from previous time step
      ! for this model level (index k)
      IF (L_ukca_trophet) THEN
        ! N2O5
        l = name2ntpindex(all_ntp, 'het_n2o5  ')
        rc_het(:,1) = all_ntp(l)%data_3d(j,i,:)
        ! HO2+HO2
        l = name2ntpindex(all_ntp, 'het_ho2   ')
        rc_het(:,2) = all_ntp(l)%data_3d(j,i,:)
      ELSE
        rc_het(:,:) = 0.0
      END IF

      ! fill stratospheric flag indicator and SO4 surface area
      stratflag(:) = (.NOT. L_troposphere(j,i,:))
      za(:) = so4_sa(j,i,:)

      ! pass 2D wet-dep field to cdrive
      zwetrt2(:,:) = zwetrt(j,i,:,:)

      IF (uph2so4inaer == 1) THEN
        ! H2SO4 will be updated in MODE, so store old value here
        ystore(:) = y(:,nn_h2so4)
      END IF
    
      ! (re-)initialise DERIV array to 1.0 before each call to ASAD_CDRIVE
      ! to ensure bit-comparability when changing domain decomposition
      deriv(:,:,:) = 1.0

      CALL asad_cdrive(cdot, zftr, zp, zt, zq, co2_1d,              &
                       cloud_frac(j,i,:),qcl(j,i,:),                &
                       j,i,klevel,zdryrt2, zwetrt2, rc_het,         &
                       zprt1d, n_pnts, have_nat(j,i,:), stratflag)

      IF (L_ukca_het_psc) THEN
        ! Save MMR of NAT PSC particles into 3-D array for PSC sedimentation.
        ! Note that sphno3 is NAT in number density of HNO3.
        IF (ANY(sphno3(:) > 0.0)) THEN
          shno3_3d(j,i,:) = (sphno3(:)/tnd(:))*c_hono2
        ELSE
          shno3_3d(j,i,:) = 0.0
        END IF
      END IF

      IF (L_ukca_chem .AND. L_ukca_nr_aqchem) THEN
        ! Calculate chemical fluxes for MODE
        IF (ihso3_h2o2 > 0) delSO2_wet_H2O2(j,i,:) =                &
          delSO2_wet_H2O2(j,i,:) + (rk(:,ihso3_h2o2)*               &
          y(:,nn_so2)*y(:,nn_h2o2))*cdt
        IF (ihso3_o3 > 0) delSO2_wet_O3(j,i,:) =                    &
          delSO2_wet_O3(j,i,:) + (rk(:,ihso3_o3)*                   &
          y(:,nn_so2)*y(:,nn_o3))*cdt
        IF (iso3_o3 > 0) delSO2_wet_O3(j,i,:) =                     &
          delSO2_wet_O3(j,i,:) + (rk(:,iso3_o3)*                    &
          y(:,nn_so2)*y(:,nn_o3))*cdt
        IF (iso2_oh > 0 .AND. ih2so4_hv > 0) THEN  ! net H2SO4 production
          delh2so4_chem(j,i,:) = delh2so4_chem(j,i,:) +             &
           ((rk(:,iso2_oh)*y(:,nn_so2)*y(:,nn_oh)) -                &
            (rk(:,ih2so4_hv)*y(:,nn_h2so4)))*cdt
        ELSE IF (iso2_oh > 0) THEN
          delh2so4_chem(j,i,:) = delh2so4_chem(j,i,:) +             &
            (rk(:,iso2_oh)*y(:,nn_so2)*y(:,nn_oh))*cdt
        END IF
        ! Restore H2SO4 tracer as it will be updated in MODE using 
        ! delh2so4_chem
        IF (uph2so4inaer == 1) y(:,nn_h2so4) = ystore(:)
      END IF

      ! 3D flux diagnostics
      IF (L_asad_use_chem_diags .AND.                               &
           ((L_asad_use_flux_rxns .OR. L_asad_use_rxn_rates) .OR.   &
           (L_asad_use_wetdep .OR. L_asad_use_drydep)))             &
           CALL asad_chemical_diagnostics(row_length,rows,          &
             model_levels,j,i,klevel,secs_per_step,volume,ierr)

      ! PSC diagnostics
      IF (L_asad_use_chem_diags .AND. L_asad_use_psc_diagnostic)    &
        CALL asad_psc_diagnostic(row_length,rows,j,i,klevel,ierr)

      ! Bring results back from vmr to mmr.

      DO l=1,jpctr
        tracer(j,i,:,l) = zftr(:,l) * c_species(l)
      END DO

      ! Set SS species concentrations for output (stratospheric 
      ! configurations)

      ! O1D mmr
      IF (O1D_in_ss) THEN
        l = name2ntpindex(all_ntp, 'O(1D)     ')
        all_ntp(l)%data_3d(j,i,:) = y(:,nn_o1d)/tnd(:)*         &
                                    c_species(nn_o1d)
      END IF

      ! O3P mmr
      IF (O3P_in_ss) THEN
        l = name2ntpindex(all_ntp, 'O(3P)     ')
        all_ntp(l)%data_3d(j,i,:) = y(:,nn_o3p)/tnd(:)*         &
                                    c_species(nn_o3p)
      END IF

      ! First copy the concentations from the zftr array to the 
      ! diag arrays, reshape and convert to moles.

      ! CH4 not masked
      IF (n_ch4 > 0) THEN
         atm_ch4_mol(j,i,:) = zftr(:,n_ch4)*tnd(:)*volume(j,i,:)* &
                              1.0e6/avogadro
      END IF

      ! CO
      IF (n_co > 0) THEN
         atm_co_mol(j,i,:) = zftr(:,n_co)*tnd(:)*volume(j,i,:)*   &
                             1.0e6/avogadro
      END IF
                            
      ! N2O
      IF (n_n2o > 0) THEN
        atm_n2o_mol(j,i,:) = zftr(:,n_n2o)*tnd(:)*volume(j,i,:)*  &
                             1.0e6/avogadro
      END IF

      ! CFC-12
      IF (n_cf2cl2 > 0) THEN
        atm_cf2cl2_mol(j,i,:) = zftr(:,n_cf2cl2)*tnd(:)*           &
                                volume(j,i,:)* 1.0e6/avogadro
      END IF

      ! CFC-11
      IF (n_cfcl3 > 0) THEN
        atm_cfcl3_mol(j,i,:) = zftr(:,n_cfcl3)*tnd(:)*              &
                               volume(j,i,:)* 1.0e6/avogadro
      END IF

      ! CH3Br
      IF (n_mebr > 0) THEN
        atm_mebr_mol(j,i,:) = zftr(:,n_mebr)*tnd(:)*                 &
                              volume(j,i,:)* 1.0e6/avogadro
      END IF

      ! H2
      IF (n_h2 > 0) THEN
        atm_h2_mol(j,i,:) = zftr(:,n_h2)*tnd(:)*                    &
                            volume(j,i,:)* 1.0e6/avogadro
      END IF
    ELSE
      cmessage='Column call is not available for Backward Euler schemes'
      errcode = 5
      CALL ereport('UKCA_CHEMISTRY_CTL_COL',errcode,cmessage)
    END IF
  END DO
END DO ! loop (j,i)
!$OMP END DO

IF (ALLOCATED(ystore)) DEALLOCATE(ystore)

!$OMP END PARALLEL

! Rescale bromine and chlorine tracers to guarantee conservation of total
! chlorine, bromine, and hydrogen over timestep. Only makes sense if at least
! chlorine chemistry is present.

IF (nn_cl > 0) THEN
  CALL ukca_conserve(row_length, rows, model_levels, ntracers,&
       tracer, pres, drain, crain, after_chem)
END IF

IF (L_ukca_strat .OR. L_ukca_stratcfc .OR. L_ukca_strattrop)      &
   THEN
  IF (L_ukca_het_psc) THEN
    ! Do NAT PSC sedimentation

    ! take NAT out of gasphase again
    tracer(:,:,:,n_hono2) = tracer(:,:,:,n_hono2) - shno3_3d

    CALL ukca_sediment(rows, row_length, model_levels, shno3_3d,  &
             qcf, r_theta_levels, mass, secs_per_step,            &
             L_troposphere(:,:,:))

    ! add solid-phase HNO3 back to gasphase HNO3
    tracer(:,:,:,n_hono2) = tracer(:,:,:,n_hono2) + shno3_3d
  END IF

  ! Tracer overwrites required to stop accumulation of tracer mass 
  ! in the uppermost layers.  Exclude water vapour.
  DO i=1,rows
    DO j=1,row_length
      DO k=1,ntracers
        IF (k /= n_h2o) THEN
          tracer(j,i,model_levels  ,k) = tracer(j,i,model_levels-2,k)
          tracer(j,i,model_levels-1,k) = tracer(j,i,model_levels-2,k)
        END IF
      END DO
    END DO
  END DO

  ! Copy NAT MMR into user_diagostics
  IF (L_ukca_het_psc) THEN
      nat_psc(:,:,:)=shno3_3d(:,:,:)
  END IF

ELSE IF (.NOT. l_ukca_offline) THEN   ! tropospheric chemistry
  ! Call routine to overwrite O3, CH4 and NOy species once per day
  ! above level defined by p_above. Only for tropospheric chemistry

  CALL ukca_stratf(i_day_number, row_length,rows, model_levels,  &
                theta_field_size, first_row, global_row_length,  &
                global_rows, jpctr, sinlat, pres,                &
                um_ozone3d, p_above,                             &
                tracer(1:row_length,1:rows,1:model_levels,       &
                       1:jpctr))

END IF     ! L_ukca_strat etc


IF (L_ukca_het_psc .AND. ALLOCATED(shno3_3d)) DEALLOCATE(shno3_3d)

IF (firstcall) firstcall = .FALSE.


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chemistry_ctl_col

END MODULE ukca_chemistry_ctl_col_mod


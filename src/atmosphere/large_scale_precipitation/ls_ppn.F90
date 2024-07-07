! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     SUBROUTINES LS_PPN and LS_PPNC------------------------------------
! Purpose:
!    LS_PPN and LS_PPNC:
!       Calculate large-scale (dynamical) precipitation.
!       LS_PPNC is the gather/scatter routine which then
!       calls LSP_ICE.
! Note: in all cases, level counters (incl subscripts) run from
!       1 (lowest model layer) to tdims%k_end
!       (topmost "wet" model layer)
!
! Programming standard: Unified Model Documentation Paper No 3
!
! Documentation: UM Documentation Paper 26.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation
MODULE ls_ppn_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LS_PPN_MOD'

CONTAINS

SUBROUTINE ls_ppn(                                                            &
  p_layer_boundaries, p_theta_levels, bland,  deltaz,                         &
!-----------------------
! primary fields and
! cloud fractions
!-----------------------
  cf, cfl, cff,                                                               &
  rhcrit,                                                                     &
  lspice_dim1,lspice_dim2,lspice_dim3,                                        &
  rho_r2, q, qcf, qcl, t,                                                     &
  qcf2, qrain, qgraup,                                                        &
!------------------------------------
! Wind field for lateral displacement
! of falling ice by shear
!------------------------------------
  u_on_p, v_on_p,                                                             &
!-----------------------
! aerosol variables
!-----------------------
  sea_salt_film, sea_salt_jet,                                                &
  salt_dim1, salt_dim2, salt_dim3,                                            &
  ukca_cdnc,                                                                  &
  cdnc_dim1, cdnc_dim2, cdnc_dim3,                                            &
  easyaerosol_cdnc,                                                           &
  biogenic,                                                                   &
  snow_depth, land_fract,                                                     &
  so4_ait,                                                                    &
  so4_acc,                                                                    &
  so4_dis,                                                                    &
  bmass_agd,                                                                  &
  bmass_cld,                                                                  &
  ocff_agd,                                                                   &
  ocff_cld,                                                                   &
  nitr_acc,                                                                   &
  nitr_diss,                                                                  &
  aerosol,                                                                    &
  n_arcl_species, n_arcl_compnts, i_arcl_compnts, arcl,                       &
!---------------------------
! Other variables for mphys
!---------------------------
  lsrain,lssnow,lsgraup_mean,                                                 &
  lsrain3d, lssnow3d, lsgraup3d, rainfrac3d,                                  &
  n_drop_tpr, n_drop_3d,                                                      &
  rhc_row_length, rhc_rows,                                                   &
  rhodz_dry, rhodz_moist,                                                     &
! variables for Jules rainfrac code
  ls_rainfrac, land_points, land_index,                                       &
!---------------------------
! COSP variables
!---------------------------
  l_cosp_lsp,                                                                 &
!------------------------------------------------------
! For subgrid orographic rain enhancement
!------------------------------------------------------
  hmteff, zb,                                                                 &
!------------------------------------------------------
! Variables for stochastic physics random parameters2
!------------------------------------------------------
  m_ci,                                                                       &
!-------------
! Error code
!-------------
  error          )

     ! General modules
USE science_fixes_mod,     ONLY: l_mphys_gr_out

USE arcl_mod,              ONLY: npd_arcl_compnts

   ! Stochastic Physics Module
USE stochastic_physics_run_mod, ONLY: x1r_rp, l_rp2, i_rp_scheme,             &
                                      i_rp2b

   ! Microphysics modules
USE mphys_inputs_mod,      ONLY: niters_mp, l_warm_new, l_mcr_qrain,          &
                                 l_mcr_qcf2, l_mcr_qgraup, x1r

USE mphys_bypass_mod,      ONLY: mphys_mod_top, l_last_iter

USE mphys_constants_mod,   ONLY: m_ci_sav, l_calc_mp_const, iter_z,           &
                                 max_it_mlt

USE mphys_diags_mod,       ONLY: l_point_diag, mphys_pts

   ! ENDGame modules
USE level_heights_mod,     ONLY: eta_theta_levels

USE atm_fields_bounds_mod, ONLY: pdims_s, tdims,                              &
                                 pdims

USE jules_hydrology_mod, ONLY: l_var_rainfrac

   ! Dr Hook modules
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

   ! Large scale precipitation modules
USE ls_ppnc_mod,           ONLY: ls_ppnc
USE lsp_taper_ndrop_mod,   ONLY: lsp_taper_ndrop
USE lspcon_mod,            ONLY: lspcon
USE mphys_air_density_mod, ONLY: mphys_air_density
USE lsprec_mod,            ONLY: lsprec_set_reals


USE def_easyaerosol, ONLY: t_easyaerosol_cdnc

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                        &
  rhc_row_length, rhc_rows,                                                   &
  lspice_dim1,lspice_dim2,lspice_dim3,                                        &
! Dimensions for 3D arrays
    salt_dim1,                                                                &
                     ! Array dimensions for sea-salt arrays (equal
    salt_dim2,                                                                &
                     ! either to tdims%i_start:tdims%i_end,
                     ! tdims%j_start:tdims%j_end, and
                     ! 1:tdims%k_end, or
    salt_dim3,                                                                &
                     ! else 1,1,1, depending on L_SEASALT_CCN).
    cdnc_dim1, cdnc_dim2, cdnc_dim3,                                          &
                     ! UKCA cloud drop number concentration dimensions
    n_arcl_species,                                                           &
                      ! Number of requested species within the climatology
    n_arcl_compnts,                                                           &
                      ! Corresponding number of requested components
    i_arcl_compnts(npd_arcl_compnts)
                      ! Array index of each aerosol clim component

REAL ::                                                                       &
  cf( tdims%i_start:tdims%i_end,                                              &
      tdims%j_start:tdims%j_end,                                              &
      1:tdims%k_end ),                                                        &
                                  ! IN Cloud fraction.
  p_theta_levels( tdims%i_start:tdims%i_end,                                  &
                  tdims%j_start:tdims%j_end,                                  &
                  1:tdims%k_end ),                                            &

! Note that the declaration beneath has been written to cope with the
! ENDGame modules, but as ENDGame hasn't yet got a tdims value for 0
! in the k dimension, this may need altering in the future

    p_layer_boundaries( tdims%i_start:tdims%i_end,                            &
                        tdims%j_start:tdims%j_end,                            &
                         0           :tdims%k_end ),                          &

    rhcrit( rhc_row_length, rhc_rows,                                         &
            1:tdims%k_end ),                                                  &
                                                 ! IN Critical humidity
                                                 ! for cloud formation.
   cfl( tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end,                                            &
        1:tdims%k_end ),                                                      &
                                        !IN Cloud liquid fraction.
   cff( tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end,                                            &
        1:tdims%k_end ),                                                      &
                                        !IN Cloud ice fraction.
   rho_r2(pdims_s%i_start:pdims_s%i_end,                                      &
          pdims_s%j_start:pdims_s%j_end,                                      &
          pdims_s%k_start:pdims_s%k_end)
                                ! IN Air density * earth radius**2

LOGICAL :: bland( tdims%i_start : tdims%i_end,                                &
               tdims%j_start : tdims%j_end )
                                   ! IN Land/sea mask

LOGICAL, INTENT(IN) ::  l_cosp_lsp ! Switch for COSP LS diagnostics

REAL, INTENT(INOUT) ::                                                        &
 q( tdims%i_start:tdims%i_end,                                                &
    tdims%j_start:tdims%j_end,                                                &
    1:tdims%k_end ),                                                          &
                                       ! Specific humidity (kg water
 qcf( tdims%i_start:tdims%i_end,                                              &
      tdims%j_start:tdims%j_end,                                              &
      1:tdims%k_end ),                                                        &
                                       ! Cloud ice (kg per kg air).
 qcl( tdims%i_start:tdims%i_end,                                              &
      tdims%j_start:tdims%j_end,                                              &
      1:tdims%k_end ),                                                        &
                                       ! Cloud liquid water (kg per
 qcf2( tdims%i_start:tdims%i_end,                                             &
       tdims%j_start:tdims%j_end,                                             &
       1:tdims%k_end ),                                                       &
                                       ! Ice (kg per kg air)
 qrain( tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end,                                            &
        1:tdims%k_end ),                                                      &
                                       ! Rain (kg per kg air)
 qgraup( tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end,                                           &
         1:tdims%k_end ),                                                     &
                                       ! Graupel (kg per kg air)
 t( tdims%i_start:tdims%i_end,                                                &
    tdims%j_start:tdims%j_end,                                                &
    1:tdims%k_end ),                                                          &
                                       ! Temperature (K).
 aerosol( tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end,                                          &
          1:tdims%k_end )
                                       ! 'Murk' tracer aerosol.

REAL, INTENT(OUT) ::                                                          &
 rhodz_dry( tdims%i_start:tdims%i_end,                                        &
            tdims%j_start:tdims%j_end,                                        &
                        1:tdims%k_end ),                                      &
                                       ! Dry density
 rhodz_moist( tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end,                                      &
                          1:tdims%k_end )
                                       ! Moist density

REAL,INTENT(IN) ::                                                            &
  ! For calculating shear in falling ice cloud fraction calculation.
  u_on_p( pdims%i_start : pdims%i_end,                                        &
          pdims%j_start : pdims%j_end,                                        &
          pdims%k_start : pdims%k_end),                                       &
  v_on_p( pdims%i_start : pdims%i_end,                                        &
          pdims%j_start : pdims%j_end,                                        &
          pdims%k_start : pdims%k_end),                                       &
  ! Aerosol climatology array:
  arcl(   tdims%i_start : tdims%i_end,                                        &
          tdims%j_start : tdims%j_end,                                        &
                      1 : tdims%k_end,                                        &
                      1 : n_arcl_compnts)

REAL,INTENT(IN) ::                                                            &
                                   !Sulphur Cycle tracers (mmr kg/kg)
   so4_ait( tdims%i_start:tdims%i_end,                                        &
            tdims%j_start:tdims%j_end,                                        &
            1:tdims%k_end ),                                                  &
   so4_acc( tdims%i_start:tdims%i_end,                                        &
            tdims%j_start:tdims%j_end,                                        &
            1:tdims%k_end ),                                                  &
   so4_dis( tdims%i_start:tdims%i_end,                                        &
            tdims%j_start:tdims%j_end,                                        &
            1:tdims%k_end ),                                                  &

                                   !Biomass smoke tracers
   bmass_agd( tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end,                                      &
              1:tdims%k_end ),                                                &
   bmass_cld( tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end,                                      &
              1:tdims%k_end ),                                                &

                                   !Fossil-fuel organic carbon tracers
   ocff_agd( tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,                                       &
             1:tdims%k_end ),                                                 &
   ocff_cld( tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,                                       &
             1:tdims%k_end ),                                                 &

                                   !Ammonium nitrate tracers
   nitr_acc( tdims%i_start:tdims%i_end,                                       &
              tdims%j_start:tdims%j_end,                                      &
              1:tdims%k_end ),                                                &
    nitr_diss( tdims%i_start:tdims%i_end,                                     &
               tdims%j_start:tdims%j_end,                                     &
               1:tdims%k_end )

REAL ::                                                                       &
    snow_depth( tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end ),                                  &
                                   ! IN Snow depth (m)
    land_fract( tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end )
                                   ! IN Land fraction

REAL ::                                                                       &
    sea_salt_film(salt_dim1, salt_dim2, salt_dim3),                           &
                                                       ! (m-3)
    sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)  ! (m-3)

REAL ::                                                                       &
    biogenic( tdims%i_start:tdims%i_end,                                      &
              tdims%j_start:tdims%j_end,                                      &
              1:tdims%k_end )
                                                ! (m.m.r.)

REAL :: ukca_cdnc(cdnc_dim1, cdnc_dim2, cdnc_dim3)
!         CDNC from UKCA for 2nd indirect effect (m-3)

TYPE (t_easyaerosol_cdnc), INTENT(IN) :: easyaerosol_cdnc
!         CDNC from EasyAerosol for 2nd indirect effect (m-3)

REAL ::                                                                       &
 lsrain( tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end ),                                         &
                             ! OUT Surface rainfall rate (kg / sq m /
 lssnow( tdims%i_start:tdims%i_end,                                           &
         tdims%j_start:tdims%j_end ),                                         &
                             ! OUT Surface snowfall rate (kg / sq m /
 lsgraup_mean( tdims%i_start:tdims%i_end,                                     &
               tdims%j_start:tdims%j_end ),                                   &
                             ! OUT Surface graupel rate (kg / sq m /
 lsgraup( tdims%i_start:tdims%i_end,                                          &
          tdims%j_start:tdims%j_end )
                             ! Graupel fall rate (kg/m2/s)

REAL, INTENT(OUT) :: n_drop_tpr( tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end ),               &
                     n_drop_3d(  tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )
                     ! Tapered droplet number and droplet number
                     ! from autoconversion scheme

REAL, INTENT(OUT) ::  deltaz( tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

REAL ::                                                                       &
    lsrain3d(lspice_dim1,lspice_dim2,lspice_dim3),                            &
                                                      ! OUT
!                           Rain rate out of each model layer
      lssnow3d(lspice_dim1,lspice_dim2,lspice_dim3),                          &
                                                        ! OUT
!                           Snow rate out of each model layer
     lsgraup3d(lspice_dim1,lspice_dim2,lspice_dim3),                          &
                                                        ! OUT
!                           Graupel rate out of each model layer

      rainfrac3d(lspice_dim1,lspice_dim2,lspice_dim3) ! OUT
!                           rain fraction out of each model layer

! variables for Jules rainfrac code
INTEGER, INTENT(IN) :: land_points
                       ! IN No.of land points being processed, can be 0.
INTEGER, INTENT(IN) ::                                                        &
  land_index( land_points )
REAL, INTENT(OUT) ::                                                          &
  ls_rainfrac(land_points)
                      ! Rain fraction array on land points to be
                      ! passed to atmos_physics2


! Variable for seeder feeder scheme
REAL,    INTENT(IN) :: hmteff(tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end)
REAL,    INTENT(IN) :: zb(tdims%i_start : tdims%i_end,                        &
                          tdims%j_start : tdims%j_end)

! Variables for stochastic physics random parameters
REAL,    INTENT(IN) :: m_ci   ! used to modify ice fall speed.

INTEGER ::    error          ! OUT Return code - 0 if OK,
!                                                - 1 if bad arguments.

!    Workspace usage ---------------------------------------------------

REAL ::  lsrain_mean(tdims%i_start:tdims%i_end,                               &
                     tdims%j_start:tdims%j_end )

REAL ::  lssnow_mean(tdims%i_start:tdims%i_end,                               &
                     tdims%j_start:tdims%j_end )

REAL ::  cf_max(tdims%i_start:tdims%i_end,                                    &
                tdims%j_start:tdims%j_end )

INTEGER ::                                                                    &
 ix(  tdims%i_len  *                                                          &
      tdims%j_len  , 2 ),                                                     &
                               ! Index for compress/expand.
  eta_ratio1               ! ratio of model level height to
                               ! iterative melting height plus 1.

REAL :: vfall( tdims%i_start:tdims%i_end,                                     &
            tdims%j_start:tdims%j_end )
             ! snow fall velocity (m per s).
REAL :: vfall2( tdims%i_start:tdims%i_end,                                    &
             tdims%j_start:tdims%j_end )
             ! fall velocity for qcf2 (m/s)
REAL :: lssnow2( tdims%i_start:tdims%i_end,                                   &
              tdims%j_start:tdims%j_end )
             ! snowfall rate for qcf2
REAL :: droplet_flux( tdims%i_start:tdims%i_end,                              &
                   tdims%j_start:tdims%j_end )
             ! water drop flux / kg m-2 s-1
REAL :: vfall_rain( tdims%i_start:tdims%i_end,                                &
                 tdims%j_start:tdims%j_end )
             ! fall velocity for rain (m/s)
REAL :: vfall_graup( tdims%i_start:tdims%i_end,                               &
                  tdims%j_start:tdims%j_end )
             ! fall vel. for graupel (m/s)
REAL :: cttemp( tdims%i_start:tdims%i_end,                                    &
             tdims%j_start:tdims%j_end )
REAL :: rainfrac( tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end )
REAL :: rainfrac_impr( tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end )
REAL :: frac_ice_above( tdims%i_start:tdims%i_end,                            &
                     tdims%j_start:tdims%j_end )
             ! Cloud ice fraction passed
             ! in layer above
REAL :: layer_thickness( tdims%i_start:tdims%i_end,                           &
                      tdims%j_start:tdims%j_end )

REAL :: iter_eta                      ! eta value at which to
                                       ! start iterative melting

!  Physical constants -------------------------------------------------
REAL :: cfmin
PARAMETER (                                                                   &
 cfmin=1.0e-3                                                                 &
                         ! Used for LS_PPNC  compress.
)
!  Define local variables ----------------------------------------------
INTEGER :: i,k,j,it,                                                          &
                      ! Loop counters: I - horizontal field index;
!                                        K - vertical level index.
             kp1,                                                             &
                        ! Index of level above: k=k+1, apart from
                        ! when k=_dims%end when kp1=k.
             n          ! "nval" for WHEN routine.

REAL :: work
                    ! work variable

LOGICAL :: l_3ddiag ! Flag to determine if we want 3d diagnostics

LOGICAL :: l_use_gridbox( tdims%i_start:tdims%i_end,                          &
                          tdims%j_start:tdims%j_end )
           ! flag for whether grid box is used in microphysics or not

REAL :: one_over_niters_mp  ! 1./niters_mp

!  Variables for Dr Hook:

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LS_PPN'


!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error=0

! Determine if we want 3d diagnosics
l_3ddiag = (lspice_dim1  ==  tdims%i_len  .AND.                               &
   lspice_dim2  ==  tdims%j_len  .AND.                                        &
   lspice_dim3  ==  tdims%k_end  )

!----------------------------------------------------------------------

! Define CX and CONSTP values

! Determine whether m_ci_rp or x1r_rp have been udpated

IF (l_calc_mp_const .OR. m_ci /= m_ci_sav .OR.                                &
    (l_rp2 .AND. i_rp_scheme == i_rp2b .AND. x1r /= x1r_rp)) THEN

  ! If the microphysics constants and number of microphysics iterations
  ! (niters_mp) are not set we can compute them
  ! and save to a module. However, if m_ci (random parameters)
  ! or x1r (random parameters version 2b) changes, some constants
  ! will change with it, so we need to recompute the constants in
  ! this case. Continuation runs (e.g. climate) will also need to
  ! recalculate the constants.
  CALL lspcon( m_ci )

  !Set the values of scheme constants in the lsprec_mod
  CALL lsprec_set_reals()

END IF

! Calculate air density (use a standard routine for both CASIM
! and Wilson and Ballard microphysics )
CALL mphys_air_density( rho_r2, q, qcl, qcf, qcf2, qrain, qgraup,             &
                        rhodz_dry, rhodz_moist, deltaz )


!-------------------------------------------------------------
! Calculation of cloud droplet number. This is now calculated
! here for all models and not in lsp_autoc as was done in
! older versions of the Unified Model
!-------------------------------------------------------------

CALL lsp_taper_ndrop(                                                         &
               ! (Full) Aerosol tracers
                        so4_ait, so4_acc, so4_dis, sea_salt_film,             &
                        biogenic, sea_salt_jet, bmass_agd,                    &
                        bmass_cld, ocff_agd, ocff_cld,                        &
                        nitr_acc, nitr_diss,                                  &
                        n_arcl_species, n_arcl_compnts,                       &
                        i_arcl_compnts, arcl,                                 &
               ! Murk aerosol
                        aerosol,                                              &
               ! CDNC from UKCA
                        ukca_cdnc,                                            &
                        cdnc_dim1, cdnc_dim2, cdnc_dim3,                      &
               ! CDNC from EasyAerosol
                        easyaerosol_cdnc,                           &
               ! Other parameters
                        rhodz_dry,  rhodz_moist,                              &
                        deltaz,                                               &
                        snow_depth, land_fract,                               &
               ! Output parameters
                        n_drop_tpr                                            &
                             )

!-----------------------------------------------------------------------
!  2. Loop round levels from top down (counting bottom level as level 1,
!     as is standard in the Unified model).
!-----------------------------------------------------------------------


!-----------------------------------------------
! Setup for iterative metlting
! Define constants outside of K loop
!-----------------------------------------------

      ! calculate level independent eta value
      ! iter_z is theta_height for level 13 in 38 level set.
      ! hence is now independent of number of levels

iter_eta = iter_z / mphys_mod_top

!-----------------------------------------------------------------------
! Internal structure.
! 2a. Initialise outside of iterative loop.
!-----------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( tdims, lsrain_mean, lssnow_mean, lsgraup_mean, lsrain3d,&
!$OMP&  lssnow3d, lsgraup3d, rainfrac3d, lspice_dim1, lspice_dim2,     &
!$OMP&  lspice_dim3 )  &
!$OMP& PRIVATE( i, j, k )
!$OMP  DO SCHEDULE(STATIC)
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    lsrain_mean(i,j) =0.0
    lssnow_mean(i,j) =0.0
    lsgraup_mean(i,j)=0.0
  END DO ! Loop over points,i
END DO ! Loop over points,j
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = 1, lspice_dim3
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      lsrain3d(i,j,k)=0.0
      lssnow3d(i,j,k)=0.0
      lsgraup3d(i,j,k)=0.0
      rainfrac3d(i,j,k)=0.0
    END DO ! Loop over points,i
  END DO ! Loop over points,j
END DO ! Loop over points,k
!$OMP END DO
!$OMP END PARALLEL
ls_rainfrac=0.0

one_over_niters_mp=1.0/niters_mp

DO it = 1, niters_mp ! Substep outside of column

  IF (it == niters_mp) THEN
    l_last_iter = .TRUE.
  ELSE
    l_last_iter = .FALSE.
  END IF
  !-----------------------------------------------------------------------
  !  Internal structure - moved inside BS iter loop
  !  2.b Initialise rain and snow to zero.
  !   Initialise scavenged amounts of S Cycle tracers to 0 for full field
  !-----------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(NONE)                                         &
!$OMP& SHARED( tdims, lsrain, lssnow, lssnow2, lsgraup, droplet_flux, &
!$OMP& cttemp, rainfrac, frac_ice_above, vfall, vfall2, vfall_rain,   &
!$OMP& vfall_graup, cf_max)                                           &
!$OMP& PRIVATE( i, j, k )
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      lsrain(i,j)=0.0
      lssnow(i,j)=0.0
      lssnow2(i,j)=0.0
      lsgraup(i,j)=0.0
      droplet_flux(i,j)=0.0
      cttemp(i,j)=0.0
      rainfrac(i,j)=0.0
      frac_ice_above(i,j)=0.0
      vfall(i,j)=0.0
      vfall2(i,j)=0.0
      vfall_rain(i,j)=0.0
      vfall_graup(i,j)=0.0
      cf_max(i,j)=-HUGE(1.0)

    END DO ! Loop over points,i
  END DO ! Loop over points,j
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  DO k = tdims%k_end, 1, -1

    ! kp1 is index of model level above, unless we are the model
    ! top in which case it is set to k. Used to calc vertical wind shear.
    IF (k == tdims%k_end) THEN
      kp1=k
    ELSE
      kp1=k+1
    END IF

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE( i, j, work)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        !-----------------------------------------------------------------------
        !  2.5 Form INDEX IX to gather/scatter variables in LS_PPNC
        !-----------------------------------------------------------------------

        ! Initialise n_drop_3d to zero before passing down code tree
        n_drop_3d(i,j,k)=0.0

        !  Set index where cloud fraction > CFMIN or where non-zero pptn
        !  Note: whenimd is functionally equivalent to WHENILE (but autotasks).


        !layer_thickness(i,j) =                                               &
        ! p_layer_boundaries(i,j,k) - p_layer_boundaries(i,j,k-1)

            ! Set up IF statement to determine whether to call the
            ! microphysics code for this grid box (i.e. if there is
            ! already condensate in the grid box or there is
            ! precipitation about to fall into the grid box)
        work = qcf(i,j,k)

            ! Include extra microphysics variables if in use
        IF (l_mcr_qcf2 )  work = work + qcf2(i,j,k) + lssnow2(i,j)
        IF (l_mcr_qrain)  work = work + qrain(i,j,k)
        IF (l_mcr_qgraup) work = work + qgraup(i,j,k)+lsgraup(i,j)
        work = work + droplet_flux(i,j)   ! droplet settling

        IF (cfl(i,j,k) > cfmin .OR.                                           &
           (lsrain(i,j)+lssnow(i,j)) > 0.0 .OR. work > 0.0) THEN
              ! include this grid box.
              ! Strictly speaking the CFL > CFMIN clause is too
              ! restrictive since ice nucleation does not require
              ! liquid water, but the code would be very messy.
          l_use_gridbox(i,j) = .TRUE.
              ! Note that mphys is done on this point if diagnostic is
              ! requested
          IF (l_point_diag) mphys_pts(i,j,k) = .TRUE.
        ELSE
          l_use_gridbox(i,j) = .FALSE.
        END IF

        ! set up improved rain fraction

        rainfrac_impr(i,j) = rainfrac(i,j)

        IF (l_mcr_qrain) THEN
          ! CF_MAX contains the maximum CF values of all the levels
          ! higher than k.
          IF (k /= tdims%k_end) THEN
            cf_max(i,j) = MAX(cf_max(i,j), cf(i,j,k+1))
          END IF

          IF (rainfrac_impr(i,j) == 0.0 .AND. qrain(i,j,k) > 0.0) THEN
            ! if cloud is present, assume this is a good proxy for the rain frac
            rainfrac_impr(i,j) = MAX( cf(i,j,k), cf_max(i,j))
          END IF
          IF (rainfrac_impr(i,j) == 0.0 .AND. qrain(i,j,k) > 0.0) THEN
            ! otherwise set to 0.5
            rainfrac_impr(i,j) = 0.5
            ! set a lower limit for stability
          ELSE IF (rainfrac_impr(i,j) <= 0.01 .AND. qrain(i,j,k) > 0.0) THEN
            rainfrac_impr(i,j) = 0.01
          END IF
        END IF

        IF (l_warm_new) THEN
          ! Always set rain fraction back to improved rain fraction
          rainfrac(i,j) = rainfrac_impr(i,j)
        END IF

      END DO ! Loop over points,i
    END DO ! Loop over points,j
!$OMP END PARALLEL DO

    n=0
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        IF (l_use_gridbox(i,j)) THEN
          n = n + 1
          ix(n,1) = i
          ix(n,2) = j
        END IF

      END DO ! Loop over points,i
    END DO ! Loop over points,j

    IF (n > 0) THEN

      CALL ls_ppnc(k,ix,n,                                                    &
                   lsrain,lssnow,lssnow2,lsgraup,droplet_flux,                &
                   cf(1,1,k),cfl(1,1,k),cff(1,1,k),                           &
                   qcf(1,1,k),qcl(1,1,k),t(1,1,k),                            &
                   qcf2(1,1,k),qrain(1,1,k),qgraup(1,1,k),                    &
                   n_drop_tpr(1,1,k), n_drop_3d(1,1,k),                       &
                   aerosol(1,1,k), land_fract,                                &
                   hmteff, zb,                                                &
                   q(1,1,k), p_theta_levels(1,1,k), layer_thickness,          &
                   deltaz(1,1,k), rhodz_dry(1,1,k), rhodz_moist(1,1,k),       &
                   rhc_row_length, rhc_rows,                                  &
                   bland, rhcrit(1,1,k),                                      &
                   vfall, vfall2, vfall_rain, vfall_graup,                    &
                   frac_ice_above,                                            &
                   cttemp, rainfrac, rainfrac_impr, niters_mp,                &
                   u_on_p(1,1,k),   v_on_p(1,1,k),                            &
                   u_on_p(1,1,kp1), v_on_p(1,1,kp1), k,                       &
                   l_cosp_lsp                                                 &
                   )
    END IF

    ! Copy rainfall and snowfall rates to 3D fields for diagnostic output

    IF ( l_3ddiag ) THEN

      ! Only copy rain and snow to 3D fields if arrays are dimensionalized.

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP& SHARED( tdims, lsrain3d, lsrain, droplet_flux, one_over_niters_mp,&
!$OMP& lssnow3d, lsgraup3d, lssnow, lssnow2, lsgraup, rainfrac3d,        &
!$OMP& rainfrac, k )                                                     &
!$OMP& PRIVATE( i, j )
      DO j = tdims%j_start, tdims%j_end

        DO i = tdims%i_start, tdims%i_end

          lsrain3d(i, j, k)   = lsrain3d(i, j, k) +                           &
             (lsrain(i, j) + droplet_flux(i, j))*one_over_niters_mp

          lssnow3d(i, j, k)   = lssnow3d(i, j, k) +                           &
             (lssnow(i, j) + lssnow2(i, j) + lsgraup(i, j))                   &
             *one_over_niters_mp

          lsgraup3d(i, j, k)  = lsgraup3d(i, j, k) +                          &
                                 ( lsgraup(i, j) * one_over_niters_mp )

          rainfrac3d(i, j, k) = rainfrac3d(i, j, k) +                         &
             rainfrac(i, j)*one_over_niters_mp


        END DO

      END DO
!$OMP END PARALLEL DO

    END IF

  END DO ! Loop over K

! save rain fraction on land points to pass to Jules
! loop over K is backwards, hence rainfrac now contains the level 1
! (surface) value
  IF (l_var_rainfrac) THEN
    DO k = 1, land_points
      j = (land_index(k)-1)/ ( tdims%i_len ) + 1
      i = land_index(k) - (j-1)*( tdims%i_len )
      ls_rainfrac(k) = ls_rainfrac(k) + rainfrac(i,j)*one_over_niters_mp
    END DO
  END IF

  ! If substepping outside of loop over K, then need to accumulate (mean)
  ! precip rates...
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE( i, j )
  DO j = tdims%j_start, tdims%j_end

    DO i = tdims%i_start, tdims%i_end

      lsrain_mean(i,j) = lsrain_mean(i,j)                                     &
         + (lsrain(i, j) + droplet_flux(i, j))*one_over_niters_mp

      ! Add together ice crystals, snow aggregates and graupel
      ! for surface snow rate (kg/m2/s)

      IF (l_mphys_gr_out) THEN

        lssnow_mean(i, j) = lssnow_mean(i, j)                                 &
                        + (lssnow(i, j)  + lssnow2(i, j)                      &
                        + lsgraup(i, j)) * one_over_niters_mp

        lsgraup_mean(i, j) = lsgraup_mean(i, j)                               &
                        + lsgraup(i, j) * one_over_niters_mp

      ELSE

        IF (l_mcr_qcf2) THEN

          lssnow_mean(i, j) = lssnow_mean(i, j)                               &
                            + (lssnow(i, j)  + lssnow2(i, j)                  &
                            + lsgraup(i, j)) * one_over_niters_mp
        ELSE ! l_mcr_qcf2

          lssnow_mean(i, j) = lssnow_mean(i, j)                               &
                            + ( lssnow(i, j) ) * one_over_niters_mp

        END IF ! l_mcr_qcf2

      END IF ! l_mphys_gr_out

    END DO ! tdims&i

  END DO ! tdims%j
!$OMP END PARALLEL DO

END DO ! Outer substepping loop over it

! Recover meaned precip rates
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                      &
!$OMP& SHARED( tdims, lssnow, lssnow_mean, lsrain, lsrain_mean )       &
!$OMP& PRIVATE( i, j )
DO j = tdims%j_start, tdims%j_end

  DO i = tdims%i_start, tdims%i_end

    lssnow(i,j) = lssnow_mean(i,j)
    lsrain(i,j) = lsrain_mean(i,j)

  END DO

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ls_ppn

END MODULE ls_ppn_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Completes lifting of the parcel from layer k to k+1

MODULE convec2_6a_mod

IMPLICIT NONE

!
! Description:
!   Completes lifting of the parcel from layer k to k+1.
!
!   Calls subroutines parcel and environ
!
!   Subroutine parcel carries out the forced
!   detrainment calculation, tests to see if convection is termintating
!   and calculates the precipitation rate from layer k+1.
!
!   Subroutine environ calculates the effect of convection upon the
!   large-scale atmosphere.
!
!   The parcel properties at level k+1 are then fully updated. CAPE and
!   dCAPEbydt are calculated for use in the CAPE closure.
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CONVEC2_6A_MOD'

CONTAINS

SUBROUTINE convec2_6a(k, npnts, np_full, nlev, ntra, ad_on, new_termc,       &
                      start_lev,                                             &
                      pstar, pk, pkp1, delpk,                                &
                      delpkp1, delpkp12, delp_uv_k, delp_uv_kp1,             &
                      exk, exkp1,                                            &
                      thek, thekp1, qek, qekp1,                              &
                      qclek, qclekp1, qcfek, qcfekp1,                        &
                      qsek, qsekp1,                                          &
                      cflek, cflekp1,  cffek,  cffekp1,                      &
                      bcfek,  bcfekp1,                                       &
                      thpk, qpk, qclpk, qcfpk,                               &
                      thpi, qpi, expi,                                       &
                      rbuoyk, rbuoykp1,                                      &
                      watldek, watldekp1, watldpk, watldpkp1,                &
                      Qlkp1, Qfkp1, Frezkp1,                                 &
                      ekp14, ekp34, amdetk, flxk, flx_init,                  &
                      uek, uekp1, vek, vekp1,                                &
                      upk, vpk,                                              &
                      traek, traekp1, trapk,                                 &
                      zkm1, zk, zkp12, zkp1, w2p_km1,                        &
                      l_q_interact, l_mom_gk, l_mom_gk_stable, l_tracer,     &
                      bgmk, bgmkp1, bwk,                                     &
                      bwkp1, blowst, bland,                                  &

                      ! In/out
                      lcbase, lctop,                                         &
                      thpkp1, qpkp1, qclpkp1, qcfpkp1,                       &
                      dthek, dqek, dqclek, dqcfek,                           &
                      tcw, depth, cclwp, lcca,                               &
                      cape, dcpbydt, relh, dptot, max_cfl,                   &
                      eflux_u_ud, eflux_v_ud,                                &
                      duek, dvek,                                            &
                      dtraek, w2p_k,                                         &
                      bterm, blatent, xsbmin,                                &

                      ! Out
                      iccb, icct,                                            &
                      dcflek, dcffek, dbcfek,                                &
                      dthekp1, dqekp1, dqclekp1, dqcfekp1,                   &
                      dcflekp1, dcffekp1, dbcfekp1,                          &
                      prekp1, thrk, qrk, deltak,                             &
                      flxkp12, flxkp1,                                       &
                      cca, ccwkp1,                                           &
                      upkp1, vpkp1,                                          &
                      duekp1, dvekp1,                                        &
                      trapkp1,                                               &
                      dtraekp1, w2p_kp1 )


USE cv_run_mod, ONLY: cnv_wat_load_opt, l_new_dd
USE cv_stash_flg_mod,    ONLY: flg_w_eqn
USE planet_constants_mod, ONLY: cp, c_virtual, r, g
USE water_constants_mod, ONLY: lc, lf
USE cv_derived_constants_mod, ONLY: ls
USE ukca_option_mod,     ONLY: l_ukca, l_ukca_plume_scav 
USE ukca_scavenging_mod, ONLY: ukca_plume_scav

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE parcel_6a_mod, ONLY: parcel_6a
USE calc_w_eqn_mod, ONLY: calc_w_eqn
USE environ_6a_mod
USE pc2_environ_mod

IMPLICIT NONE

! Subroutine arguments

INTEGER,INTENT(IN) :: k             ! present model layer
INTEGER,INTENT(IN) :: npnts         ! Number of points
INTEGER,INTENT(IN) :: np_full       ! Full vector length
INTEGER,INTENT(IN) :: nlev          ! Number of model levels for calculations
INTEGER,INTENT(IN) :: ntra          ! Number of tracer variables
INTEGER,INTENT(IN) :: ad_on         ! Flag for adaptive detrainment
INTEGER,INTENT(IN) :: new_termc     ! Flag for simplified termination of
                                    ! convection

INTEGER,INTENT(IN) :: start_lev(npnts) ! Level at which convection is initiated

REAL,INTENT(IN) :: pstar(npnts)     ! Surface pressure (Pa)
REAL,INTENT(IN) :: pk(npnts)        ! pressure at mid-point of layer k (Pa)
REAL,INTENT(IN) :: pkp1(npnts)      ! pressure at mid-point of layer k+1 (Pa)
REAL,INTENT(IN) :: delpk(npnts)     ! pressure difference across layer k (Pa)
REAL,INTENT(IN) :: delpkp1(npnts)   ! pressure difference across layer k+1 (Pa)
REAL,INTENT(IN) :: delpkp12(npnts)  ! pressure diff. across layer k+1/2 (Pa)
REAL,INTENT(IN) :: delp_uv_k(npnts) ! pressure difference across UV
                                    ! layer k (Pa)
REAL,INTENT(IN) :: delp_uv_kp1(npnts) ! pressure difference across UV
                                    ! layer k+1 (Pa)
REAL,INTENT(IN) :: exk(npnts)       ! Exner ratio at mid-point of layer k
REAL,INTENT(IN) :: exkp1(npnts)     ! Exner ratio at mid-point of layer k+1

REAL,INTENT(IN) :: thek(npnts)      ! Env. pot. temperature in layer k (K)
REAL,INTENT(IN) :: thekp1(npnts)    ! Env. pot. temperature in layer k+1 (K)
REAL,INTENT(IN) :: qek(npnts)       ! Env. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qekp1(npnts)     ! Env. spec. humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qclek(npnts)     ! Env. qcl in layer k (kg/kg)
REAL,INTENT(IN) :: qclekp1(npnts)   ! Env. qcl in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qcfek(npnts)     ! Env. qcf in layer k (kg/kg)
REAL,INTENT(IN) :: qcfekp1(npnts)   ! Env. qcf in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qsek(npnts)      ! Env. saturated specific humidity in
                                    ! in layer k (kg/kg)
REAL,INTENT(IN) :: qsekp1(npnts)    ! Env. saturated specific humidity in
                                    ! in layer k+1 (kg/kg)
REAL,INTENT(IN) :: cflek(npnts)     ! Env. liquid cloud volume fraction
                                    ! in layer k
REAL,INTENT(IN) :: cflekp1(npnts)   ! Env. liquid cloud volume fraction
                                    ! in layer k+1
REAL,INTENT(IN) :: cffek(npnts)     ! Env. frozen cloud volume fraction
                                    ! in layer k
REAL,INTENT(IN) :: cffekp1(npnts)   ! Env. frozen cloud volume fraction
                                    ! in layer k+1
REAL,INTENT(IN) :: bcfek(npnts)     ! Env. total cloud volume fraction
                                    ! in layer k
REAL,INTENT(IN) :: bcfekp1(npnts)   ! Env. total cloud volume fraction
                                    ! in layer k+1
REAL,INTENT(IN) :: thpk(npnts)      ! Par. pot. temperature in layer k (K)
REAL,INTENT(IN) :: qpk(npnts)       ! Par. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qclpk(npnts)     ! Par. qcl in layer k (kg/kg)
REAL,INTENT(IN) :: qcfpk(npnts)     ! Par. qcf in layer k (kg/kg)
REAL,INTENT(IN) :: thpi(npnts)      ! Initial parcel potential temperature (K)
REAL,INTENT(IN) :: qpi(npnts)       ! Initial parcel specific humidity (kg/kg)
REAL,INTENT(IN) :: expi(npnts)      ! Initial parcel Exner pressure
REAL,INTENT(IN) :: rbuoyk(npnts)    ! Par. buoyancy in layer k (K)
REAL,INTENT(IN) :: rbuoykp1(npnts)  ! Par. buoyancy in layer k+1 (K)
REAL,INTENT(IN) :: watldek(npnts)   ! Env. water loading in layer k (kg/kg)
REAL,INTENT(IN) :: watldpk(npnts)   ! Par. water loading in layer k (kg/kg)
REAL,INTENT(IN) :: watldekp1(npnts) ! Env. water loading in layer k+1 (kg/kg)
REAL,INTENT(IN) :: watldpkp1(npnts) ! Par. water loading in layer k+1 (kg/kg)
REAL,INTENT(IN) :: Qlkp1(npnts)     ! Amount of condensation to liquid water
                                    ! in the parcel (kg/kg)
REAL,INTENT(IN) :: Qfkp1(npnts)     ! Amount of deposition to ice water
                                    ! in the parcel (kg/kg)
REAL,INTENT(IN) :: Frezkp1(npnts)   ! Amount of freezing from liquid
                                    ! to ice water in the parcel (kg/kg)
REAL,INTENT(IN) :: ekp14(npnts)     ! Entrainment coefficient at level k+1/4
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: ekp34(npnts)     ! Entrainment coefficient at level k+3/4
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: amdetk(npnts)    ! Mixing detrainment coefficient at level k
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: flxk(npnts)      ! Parcel massflux in layer k (Pa/s)
REAL,INTENT(IN) :: flx_init(npnts)  ! Initial par. massflux at cloud base (Pa/s)
REAL,INTENT(IN) :: uek(npnts)       ! Env. U in layer k (m/s)
REAL,INTENT(IN) :: uekp1(npnts)     ! Env. U in layer k+1 (m/s)
REAL,INTENT(IN) :: vek(npnts)       ! Env. V in layer k (m/s)
REAL,INTENT(IN) :: vekp1(npnts)     ! Env. V in layer k+1 (m/s)
REAL,INTENT(IN) :: upk(npnts)       ! Par. U in layer k (m/s)
REAL,INTENT(IN) :: vpk(npnts)       ! Par. V in layer k (m/s)
REAL,INTENT(IN) :: traek(np_full,ntra)    ! Env. tracer content
                                          ! in layer k (kg/kg)
REAL,INTENT(IN) :: traekp1(np_full,ntra)  ! Env. tracer content
                                          ! in layer k+1 (kg/kg)
REAL,INTENT(IN) :: trapk(np_full,ntra)    ! Par. tracer content
                                          ! in layer k (kg/kg)

! Height above surface of model level ...
REAL,INTENT(IN) :: zkm1(npnts)  ! ...k-1   [m]
REAL,INTENT(IN) :: zk(npnts)    ! ...k     [m]
REAL,INTENT(IN) :: zkp12(npnts) ! ...k+1/2 [m]
REAL,INTENT(IN) :: zkp1(npnts)  ! ...k+1   [m]

REAL,INTENT(IN) :: w2p_km1(npnts)   ! (Parcel vertical velocity)^2
                                    ! on level k-1, [(m/s)^2]

LOGICAL,INTENT(IN) :: l_q_interact  ! True if PC2 is switched on
LOGICAL,INTENT(IN) :: l_mom_gk      ! Switch for inclusion of Gregory-Kershaw
                                    ! CMT
LOGICAL,INTENT(IN) :: l_mom_gk_stable ! Switch for stabilised Gregory-Kershaw
                                      ! CMT
LOGICAL,INTENT(IN) :: l_tracer      ! Switch for tracers

LOGICAL,INTENT(IN) :: bgmk(npnts)   ! mask for parcels which are saturated
                                    ! in layer k
LOGICAL,INTENT(IN) :: bgmkp1(npnts) ! Mask for parcels which are saturated
                                    ! in layer k+1
LOGICAL,INTENT(IN) :: bwk(npnts)    ! mask for parcels which have liquid
                                    ! condensate in layer k
LOGICAL,INTENT(IN) :: bwkp1(npnts)  ! mask for parcels which have liquid
                                    ! condensate in layer k+1
LOGICAL,INTENT(IN) :: blowst(npnts) ! mask for those points at which stability
                                    ! is low enough for convection to occur
LOGICAL,INTENT(IN) :: bland(npnts)  ! Land/sea mask

!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------
INTEGER,INTENT(INOUT) :: lcbase(npnts)! Lowest conv. cloud base level
INTEGER,INTENT(INOUT) :: lctop(npnts) ! Lowest conv. cloud top level

REAL,INTENT(INOUT) :: thpkp1(npnts) ! Par. pot. temperature in layer k+1 (K)
                                    ! IN after entrainment and latent heating
                                    ! OUT after forced detrainment
REAL,INTENT(INOUT) :: qpkp1(npnts)  ! Par. spec. humidity in layer k+1 (kg/kg)
                                    ! IN after entrainment and latent heating
                                    ! OUT after forced detrainment
REAL,INTENT(INOUT) :: qclpkp1(npnts)! Par. qcl in layer k+1 (kg/kg)
                                    ! IN after entrainment and latent heating
                                    ! OUT after forced detrainment
REAL,INTENT(INOUT) :: qcfpkp1(npnts)! Par. qcf in layer k+1 (kg/kg)
                                    ! IN after entrainment and latent heating
                                    ! OUT after forced detrainment

! Convection increments to model fields at level k
! IN before processes at level k. This may be non-zero because of smoothed
! forced detrainment and the initial perturbation
! OUT after processes at level k
REAL,INTENT(INOUT) :: dthek(npnts)  ! Increment to p. temperature
                                    ! in layer k (K/s)
REAL,INTENT(INOUT) :: dqek(npnts)   ! Increment to spec. humidity
                                    ! in layer k (kg/kg/s)
REAL,INTENT(INOUT) :: dqclek(npnts) ! Increment to qcl in layer k (kg/kg/s)
REAL,INTENT(INOUT) :: dqcfek(npnts) ! Increment to qcf in layer k (kg/kg/s)
REAL,INTENT(INOUT) :: dcflek(npnts) ! Increment to liquid cloud volume fraction
                                    ! in layer k
REAL,INTENT(INOUT) :: dcffek(npnts) ! Increment to frozen cloud volume fraction
                                    ! in layer k
REAL,INTENT(INOUT) :: dbcfek(npnts) ! Increment to total cloud volume fraction
                                    ! in layer k
REAL,INTENT(INOUT) :: tcw(npnts)    ! Total condensed water (kg/m**2/s)
                                    ! IN summed to layer k
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: depth(npnts)  ! Depth of convective cloud (m)
                                    ! IN summed to layer k
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: cclwp(npnts)  ! Condensed water path (kg/m**2)
                                    ! IN summed to layer k
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: lcca(npnts)   ! Lowest conv. cloud amount (%)
REAL,INTENT(INOUT) :: cape(npnts)   ! Convective available potential energy
                                    ! IN up to layer k-1
                                    ! OUT up to layer k
REAL,INTENT(INOUT) :: dcpbydt(npnts)! Rate of change of CAPE
                                    ! IN up to layer k-1
                                    ! OUT up to layer k
REAL,INTENT(INOUT) :: relh(npnts)   ! Relative humidity integral (averaged
                                    ! when convection terminates)
REAL,INTENT(INOUT) :: dptot(npnts)  ! Delta P integral
REAL,INTENT(INOUT) :: max_cfl(npnts)! CFL ratio
REAL,INTENT(INOUT) :: eflux_u_ud(npnts) ! Eddy flux of momentum to UD
                                    ! IN at bottom of layer
                                    ! OUT at top of layer
REAL,INTENT(INOUT) :: eflux_v_ud(npnts) ! Eddy flux of momentum to UD
                                    ! IN at bottom of layer
                                    ! OUT at top of layer
! Convection increments to winds and tracers at level k
! IN before processes at level k. This may be non-zero because of smoothed
! forced detrainment and the initial perturbation
! OUT after processes at level k
REAL,INTENT(INOUT) :: duek(npnts)   ! Increment to U in layer k (m/s**2)
REAL,INTENT(INOUT) :: dvek(npnts)   ! Increment to V in layer k (m/s**2)
REAL,INTENT(INOUT) :: dtraek(np_full,ntra)  ! Increment to tracer in layer k
                                            ! (kg/kg/s)

REAL, INTENT(INOUT) :: w2p_k(npnts)     ! (Parcel vertical velocity)^2
                                        ! on level k, [(m/s)^2]

LOGICAL,INTENT(INOUT) :: bterm(npnts)   ! Mask for parcels which terminate
                                        ! in layer k+1
LOGICAL,INTENT(INOUT) :: blatent(npnts) ! Mask for points where latent heat has
                                        ! been released
REAL,INTENT(INOUT) :: xsbmin(npnts)     ! Threshold buoyancy for forced
                                        ! detrainment (K)

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------
INTEGER,INTENT(INOUT) :: iccb(npnts)  ! convective cloud base_level
INTEGER,INTENT(INOUT) :: icct(npnts)  ! convective cloud top level
! Need to be inout as need to be remembered between subsequent calls to convec2

! Convection increments to model fields at level k+1
REAL,INTENT(OUT) :: dthekp1(npnts)  ! Increment to p. temperature
                                      ! in layer k+1 (K/s)
REAL,INTENT(OUT) :: dqekp1(npnts)   ! Increment to spec. humidity
                                      ! in layer k+1 (kg/kg/s)
REAL,INTENT(OUT) :: dqclekp1(npnts) ! Increment to qcl in layer k+1 (kg/kg/s)
REAL,INTENT(OUT) :: dqcfekp1(npnts) ! Increment to qcf in layer k+1 (kg/kg/s)
REAL,INTENT(OUT) :: dcflekp1(npnts) ! Increment to liquid cloud volume fraction
                                      ! in layer k+1
REAL,INTENT(OUT) :: dcffekp1(npnts) ! Increment to frozen cloud volume fraction
                                      ! in layer k+1
REAL,INTENT(OUT) :: dbcfekp1(npnts) ! Increment to total cloud volume fraction
                                      ! in layer k+1
REAL,INTENT(OUT) :: prekp1(npnts)   ! precipitation from parcel as it rises
                                    ! from layer k to k+1 (kg/m**2/s)
REAL,INTENT(OUT) :: thrk(npnts)     ! pot. temperature of forced detrained
                                    ! parcel in layer k (K)
REAL,INTENT(OUT) :: qrk(npnts)      ! Specific humidity of forced detrained
                                    ! parcel in layer k (kg/kg)
REAL,INTENT(OUT) :: deltak(npnts)   ! Parcel forced detrainment rate in
                                    ! layer k multiplied by layer thickness
REAL,INTENT(OUT) :: flxkp12(npnts)  ! parcel massflux in layer k+1/2 (Pa/s)
REAL,INTENT(OUT) :: flxkp1(npnts)   ! parcel massflux in layer k+1 (Pa/s)
REAL,INTENT(OUT) :: cca(npnts)      ! convective cloud amount (%)
REAL,INTENT(OUT) :: ccwkp1(npnts)   ! Total condensate in level k+1 (kg/kg)
REAL,INTENT(OUT) :: upkp1(npnts)    ! Par. U in layer k+1 (m/s)
REAL,INTENT(OUT) :: vpkp1(npnts)    ! Par. V in layer k+1 (m/s)
REAL,INTENT(OUT) :: duekp1(npnts)   ! Increment to U in layer k+1 (m/s**2)
REAL,INTENT(OUT) :: dvekp1(npnts)   ! Increment to V in layer k+1 (m/s**2)
REAL,INTENT(OUT) :: trapkp1(np_full,ntra)  ! Par. tracer content
                                           ! in layer k+1 (kg/kg)
REAL,INTENT(OUT) :: dtraekp1(np_full,ntra) ! Increment to tracer in layer k+1
                                           ! (kg/kg/s)
REAL,INTENT(OUT) :: w2p_kp1(npnts)  ! (Parcel vertical velocity)^2
                                    ! on level k+1, [(m/s)^2]

!-------------------------------------------------------------------------------
! Local variables

INTEGER :: i              ! loop counter
INTEGER :: ktra           ! Loop counter for tracers

REAL :: thvp              ! Virtual tempature of parcel
REAL :: thve              ! Virtual tempature of environment
REAL :: rho               ! Density required in CAPE calculations
REAL :: tmp_dcpbydt       ! Temporary dcpbydt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CONVEC2_6A'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!  Complete lifting parcel to layer k+1.
!
!  SUBROUTINE parcel
!
!  UM Documentation paper 27
!  Sections (5),(6),(7),(8),(9)
! ----------------------------------------------------------------------
!

CALL parcel_6a (k, npnts, np_full, nlev, ntra, ad_on, new_termc,       &
                start_lev,                                             &
                pstar, pk, pkp1, delpkp1, exk, exkp1,                  &
                thek, thekp1, qek, qekp1,                              &
                qclek, qclekp1, qcfek, qcfekp1,                        &
                qsek, qsekp1,                                          &
                thpk, qpk, qclpk, qcfpk,                               &
                thpi, qpi, expi,                                       &
                rbuoyk, rbuoykp1,                                      &
                watldek, watldekp1, watldpk, watldpkp1,                &
                Qlkp1, Qfkp1, Frezkp1,                                 &
                ekp14, ekp34, amdetk, flxk, flx_init,                  &
                uek, uekp1, vek, vekp1,                                &
                upk, vpk,                                              &
                traek, traekp1, trapk,                                 &
                l_q_interact, l_mom_gk, l_tracer,                      &
                bgmk, bgmkp1, bwk, bwkp1, blowst, bland,               &
                ! In/out
                lcbase, lctop,                                         &
                thpkp1, qpkp1, qclpkp1, qcfpkp1,                       &
                tcw, depth, cclwp, lcca,                               &
                bterm, blatent, xsbmin,                                &
                ! Out
                iccb, icct,                                            &
                prekp1, thrk, qrk, deltak, flxkp1, flxkp12,            &
                cca, ccwkp1,                                           &
                upkp1, vpkp1,                                          &
                trapkp1)

! ---------------------------------------------------------------------- 
!  Calculate the change in parcel tracer content due to scavenging 
!  by precipitation for UKCA tracers. 
!  UM Documentation paper 84 
! ---------------------------------------------------------------------- 
 
IF (l_tracer .AND. l_ukca .AND. l_ukca_plume_scav .AND. npnts > 0     &
    .AND. ANY(prekp1 > 1e-30)) THEN 
 
  CALL ukca_plume_scav(k, np_full, ntra, npnts, trapkp1(1:npnts,:),   &
                       qclpkp1, prekp1, flxkp1) 
END IF 


! ----------------------------------------------------------------------
!  Calculate the increments to environmental theta, q, qcl, qcf,
!  winds and tracers.
! ----------------------------------------------------------------------

CALL environ_6a (k, npnts, np_full, ntra,                              &
              delpk, delpkp1, delpkp12, delp_uv_k, delp_uv_kp1,     &
              exk, exkp1,                                           &
              thek, thekp1, qek, qekp1,                             &
              qclek, qclekp1, qcfek, qcfekp1,                       &
              thpk, thpkp1, qpk, qpkp1,                             &
              qclpk, qclpkp1, qcfpk, qcfpkp1,                       &
              thrk, qrk,                                            &
              ekp14, amdetk, deltak, flxk,                          &
              uek, uekp1, vek, vekp1,                               &
              upk, upkp1, vpk, vpkp1,                               &
              traek, traekp1, trapk, trapkp1,                       &
              l_q_interact, l_mom_gk, l_mom_gk_stable, l_tracer,    &
              blowst, bterm,                                        &
              ! In/out
              dthek, dqek, dqclek, dqcfek,                          &
              eflux_u_ud, eflux_v_ud,                               &
              duek, dvek,                                           &
              dtraek,                                               &
              ! Out
              dthekp1, dqekp1, dqclekp1, dqcfekp1,                  &
              duekp1, dvekp1,                                       &
              dtraekp1)

! ----------------------------------------------------------------------
!  Calculate the increments to environmental cloud fractions
! ----------------------------------------------------------------------
CALL pc2_environ (k, npnts,                                         &
              qclek, qclekp1, qcfek, qcfekp1,                       &
              cflek, cflekp1,  cffek,  cffekp1,                     &
              bcfek,  bcfekp1,                                      &
              qclpk, qclpkp1, qcfpk, qcfpkp1,                       &
              dqclek, dqcfek,                                       &
              dqclekp1, dqcfekp1,                                   &
              l_q_interact,                                         &
              bterm,                                                &
              ! Out
              dcflek, dcffek, dbcfek,                               &
              dcflekp1, dcffekp1, dbcfekp1)


! ---------------------------------------------------------------------
!  Calculate contribution to CAPE, the rate of change of CAPE  due to
!  the updraught and Courant number
! ---------------------------------------------------------------------


DO i=1,npnts
  thvp    = thpk(i)*(1.0+c_virtual*qpk(i)-watldpk(i))
  thve    = thek(i)*(1.0+c_virtual*qek(i)-watldek(i))
  rho     = pk(i)/(r*thve*exk(i))

  cape(i) = cape(i)+(thvp-thve)*delpk(i)/(rho*thve)
  relh(i) = relh(i)+(qek(i)/qsek(i))*delpk(i)
  dptot(i)= dptot(i)+delpk(i)

  IF (cnv_wat_load_opt == 0) THEN
    !   Water loading is not included and therefore do not include
    !   the contribution from the condensate increments
    tmp_dcpbydt = ( dthek(i)*(1.0+c_virtual*qek(i))             &
                  + thek(i)*c_virtual*dqek(i) )                 &
                  * (delpk(i)/(rho*thve))
  ELSE
    !   Water loading is on and therefore include the contribution
    !   from the condensate increments. If PC2 is off the condensate
    !   increments will be zero.
    tmp_dcpbydt = ( dthek(i)*(1.0+c_virtual*qek(i)-watldek(i))  &
                  + thek(i)*(c_virtual*dqek(i)                  &
                  - dqclek(i) - dqcfek(i)) )                    &
                  * (delpk(i)/(rho*thve))
  END IF

  IF (tmp_dcpbydt  >   0.0) THEN
    dcpbydt(i) = dcpbydt(i) + tmp_dcpbydt
  END IF

  ! Calculate courant number using mass flux at half level
  max_cfl(i)=MAX(max_cfl(i),                                    &
             flxk(i)/delpk(i)*(1+ekp14(i))*(1.0-deltak(i))*(1.0-amdetk(i)))

END DO


IF (flg_w_eqn .OR. l_new_dd) THEN
  CALL calc_w_eqn ( k, npnts, np_full, nlev, bterm, blowst      &
                  , ekp14, ekp34, zkm1, zk, zkp12, zkp1         &
                  , thek,  thekp1,  thpk,  thpkp1               &
                  , qek,   qekp1,   qpk,   qpkp1                &
                  , qclek, qclekp1, qclpk, qclpkp1              &
                  , qcfek, qcfekp1, qcfpk, qcfpkp1              &
                  , w2p_km1, w2p_k, w2p_kp1 )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE convec2_6a
END MODULE convec2_6a_mod

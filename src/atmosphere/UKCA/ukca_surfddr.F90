! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_surfddr_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_SURFDDR_MOD'

CONTAINS

SUBROUTINE ukca_surfddr( row_length, rows, ntype, npft, sinlat, t0, p0, rh,   &
                         smr, gsf, stcon, t0tile, lai_ft, canwc, so4_vd, rc,  &
                         o3_stom_frac)

USE asad_mod,                ONLY: &
    ndepd,                         &
    nldepd,                        &
    speci

USE ereport_mod,             ONLY: &
    ereport

USE errormessagelength_mod,  ONLY: &
    errormessagelength

USE jules_surface_types_mod, ONLY: &
    brd_leaf,                      &
    brd_leaf_dec,                  &
    brd_leaf_eg_trop,              &
    brd_leaf_eg_temp,              &
    ndl_leaf,                      &
    ndl_leaf_dec,                  &
    ndl_leaf_eg,                   &
    c3_grass,                      &
    c3_crop,                       &
    c3_pasture,                    &
    c4_grass,                      &
    c4_crop,                       &
    c4_pasture,                    &
    shrub,                         &
    shrub_dec,                     &
    shrub_eg,                      &
    urban,                         &
    lake,                          &
    soil,                          &
    ice,                           &
    elev_ice

USE parkind1,                ONLY: &
    jprb,                          &
    jpim

USE science_fixes_mod,       ONLY: &
    l_fix_improve_drydep

USE ukca_constants,          ONLY: &
    rmol,                          &
    rhow

USE ukca_option_mod,         ONLY: &
    jpdd

USE umPrintMgr,              ONLY: &
    PrintStatus,                   &
    PrStatus_Diag,                 &
    umPrint,                       &
    umMessage

USE yomhook,                 ONLY: &
    lhook,                         &
    dr_hook

IMPLICIT NONE

!     Null resistance for deposition (1/r_null ~ 0)

REAL, PARAMETER :: r_null = 1.0e50

INTEGER, INTENT(IN) :: row_length         ! number columns
INTEGER, INTENT(IN) :: rows               ! number of rows
INTEGER, INTENT(IN) :: ntype              ! number of surface types
INTEGER, INTENT(IN) :: npft               ! number of plant functional types
REAL, INTENT(IN) :: sinlat(row_length,rows)
                ! Sine(latitude)
REAL, INTENT(IN) :: t0(row_length,rows)
                ! Surface temperature (K)
REAL, INTENT(IN) :: p0(row_length,rows)
                ! Surface pressure (Pa)
REAL, INTENT(IN) :: rh(row_length,rows)
                ! Relative humidity (fraction)
REAL, INTENT(IN) :: smr(row_length,rows)
                ! Soil moisture content (Fraction by volume)

REAL, INTENT(IN) :: gsf(row_length,rows,ntype)
                ! Global surface fractions
REAL, INTENT(IN) :: t0tile(row_length,rows,ntype)
                ! Surface temperature on tiles (K)

REAL, INTENT(IN) ::  stcon(row_length,rows,npft)
                ! Stomatal conductance (m s-1)
REAL, INTENT(IN) :: lai_ft(row_length,rows,npft)
                ! Leaf area index (m2 leaf m-2)
REAL, INTENT(IN) :: canwc(row_length,rows,npft)
                ! Canopy water content (mm)

REAL, INTENT(IN) :: so4_vd(row_length,rows)
                ! Aerosol deposition velocity (m s-1)
                !(assumed to be the same for SO4 and other aerosols)

!     Surface resistance on tiles (s m-1).

REAL, INTENT(OUT) :: rc(row_length,rows,ntype,jpdd)
REAL, INTENT(OUT) :: o3_stom_frac(row_length,rows)

!     Local variables
INTEGER :: errcode                   ! error code
LOGICAL, SAVE :: first = .TRUE.
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER :: i         ! Loop count over longitudes
INTEGER :: k         ! Loop count over latitudes
INTEGER :: j         ! Loop count over species that deposit
INTEGER :: n         ! Loop count over tiles

REAL :: sm            ! Soil moisture content of gridbox
REAL :: rr            ! General temporary store for resistances.
REAL :: ts            ! Temperature of a particular tile.
REAL :: f             ! Factor to modify CH4 uptake fluxes
!     REAL :: r_cuticle_so2 ! Cuticular resistance for SO2
!     REAL :: r_cuticle_nh3 ! Cuticular resistance for NH3
REAL :: mml = 1.008e5 ! Factor to convert methane flux to dry dep vel.
REAL :: r_wet_o3 = 500.0  ! Wet soil surface resistance for ozone (s m-1)
REAL :: cuticle_o3 = 5000.0 ! Constant for caln of O3 cuticular resistance
REAL :: tundra_s_limit = 0.866 ! Southern limit for tundra (SIN(60))

!     MML - Used to convert methane flux in ug m-2 h-1 to dry dep vel in
!     m s-1; MML=3600*0.016*1.0E9*1.75E-6, where 0.016=RMM methane (kg),
!     1.0E9 converts ug -> kg, 1.75E-6 = assumed CH4 vmr

REAL :: r_cut_o3(row_length,npft) ! Cuticular Resistance for O3
REAL :: r_stom(row_length,rows,npft,5) ! Stomatal resistance
REAL, ALLOCATABLE, SAVE :: rsurf(:,:)        ! Standard surface resistances

!     Scaling of CH4 soil uptake to match present-day TAR value of 30 Tg/year

REAL, PARAMETER :: TAR_scaling = 15.0
REAL, PARAMETER :: glmin       = 1.0e-6  ! Minimum leaf conductance

!     Following arrays used to set up surface resistance array rsurf.
!      Must take the same dimension as ntype
REAL, ALLOCATABLE, SAVE :: zero(:)
REAL, ALLOCATABLE, SAVE :: rooh(:)
REAL, ALLOCATABLE, SAVE :: aerosol(:)
REAL, ALLOCATABLE, SAVE :: tenpointzero(:)

!     CH4 uptake fluxes, in ug m-2 hr-1
REAL, ALLOCATABLE, SAVE :: ch4_up_flux(:)


!     Hydrogen - linear dependence on soil moisture (except savannah)
!      Must take the same dimension as npft
REAL, ALLOCATABLE, SAVE :: h2dd_c(:)
REAL, ALLOCATABLE, SAVE :: h2dd_m(:)
REAL :: h2dd_q = 0.27 ! Quadratic term for H2 loss to savannah

!     Resistances for Tundra if different to standard value in rsurf().

REAL :: r_tundra(5) =                                 &
  (/ 1200.0, 25000.0, 800.0, 3850.0, 1100.0/)
!     NO2      CO      O3      H2     PAN

!     CH4 loss to tundra - Cubic polynomial fit to data.
!     N.B. Loss flux is in units of ug(CH4) m-2 s-1

REAL :: ch4dd_tun(4) = (/ -4.757e-6, 4.0288e-3,       &
                                     -1.13592, 106.636 /)

!     HNO3 dry dep to ice; quadratic dependence

REAL :: hno3dd_ice(3) = (/ -13.57, 6841.9, -857410.6 /)

!     SO2 dry dep to snow/ice; quadratic dependence

REAL :: so2dd_ice(3) = (/ 0.0001, 0.003308, 0.1637 /)

!   Diffusion correction for stomatal conductance Wesley (1989) 
!                                     Atmos. Env. 23, 1293.
REAL :: dif(5) = (/ 1.6, 1.6, 2.6, 1.9, 0.97 /)
!                              NO2   O3  PAN  SO2  NH3

LOGICAL :: todo(row_length,rows,ntype)     ! True if tile fraction > 0.0
LOGICAL :: microb(row_length,rows)   ! True if T > 5 C and RH > 40%
!                  (i.e. microbes in soil are active).

! Number of species for which surface resistance values are not set
INTEGER  :: n_nosurf

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SURFDDR'


! Set up standard resistance array rsurf on first call only
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first) THEN
  ALLOCATE(rsurf(ntype,jpdd))
  ALLOCATE(zero(ntype))
  ALLOCATE(rooh(ntype))
  ALLOCATE(aerosol(ntype))
  ALLOCATE(tenpointzero(ntype))
  ALLOCATE(ch4_up_flux(ntype))
  ALLOCATE(h2dd_c(npft))
  ALLOCATE(h2dd_m(npft))
  rsurf(:,:) = r_null
  zero(:)    = r_null
  tenpointzero(:) = 10.0
  
  ! Check if we have standard JULES tile configuration and fail if not as
  ! Dry deposition in UKCA has not been set up for this.
  ! Use standard tile set up or extend UKCA dry depostion
  WRITE(cmessage,'(A)')                                                       &
       'UKCA does not handle flexible tiles yet: '//                          &
       'Dry deposition needs extending. '         //                          &
       'Please use standard tile configuration'
  SELECT CASE (ntype)
  CASE (9)
    IF ( brd_leaf /= 1  .OR.                                                  &
         ndl_leaf /= 2  .OR.                                                  &
         c3_grass /= 3  .OR.                                                  &
         c4_grass /= 4  .OR.                                                  &
         shrub    /= 5  .OR.                                                  &
         urban    /= 6  .OR.                                                  &
         lake     /= 7  .OR.                                                  &
         soil     /= 8  .OR.                                                  &
         ice      /= 9 ) THEN
        
      ! Tile order changed from standard setup i.e. MOSES-II.
      errcode=1009
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  CASE (13,17,27)
    IF ( brd_leaf_dec     /= 1   .OR.                                         &
         brd_leaf_eg_trop /= 2   .OR.                                         &
         brd_leaf_eg_temp /= 3   .OR.                                         &
         ndl_leaf_dec     /= 4   .OR.                                         &
         ndl_leaf_eg      /= 5   .OR.                                         &
         c3_grass         /= 6 ) THEN 

      ! Tile order does not match that expected
      errcode=1001
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  CASE DEFAULT
    ! ntype must equal 9, 13, 17, 27
    errcode=1002
    CALL ereport(RoutineName,errcode,cmessage)
  END SELECT
  
  SELECT CASE (ntype)
  CASE (13)
    IF ( c4_grass         /= 7  .OR.                                          &
         shrub_dec        /= 8  .OR.                                          &
         shrub_eg         /= 9  .OR.                                          &
         urban            /= 10 .OR.                                          &
         lake             /= 11 .OR.                                          &
         soil             /= 12 .OR.                                          &
         ice              /= 13 ) THEN
        
      ! Tile order does not match that expected
      errcode=1013
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  
  CASE (17, 27)
    IF ( c3_crop          /= 7  .OR.                                          &
         c3_pasture       /= 8  .OR.                                          &
         c4_grass         /= 9  .OR.                                          &
         c4_crop          /= 10 .OR.                                          &
         c4_pasture       /= 11 .OR.                                          &
         shrub_dec        /= 12 .OR.                                          &
         shrub_eg         /= 13 .OR.                                          &
         urban            /= 14 .OR.                                          &
         lake             /= 15 .OR.                                          &
         soil             /= 16 .OR.                                          &
         ice              /= 17 ) THEN
        
      ! Tile order does not match that expected
      errcode=1727
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  END SELECT
  
  IF ( ntype == 27) THEN
    IF ( elev_ice(1)      /= 18 .OR.                                          &
         elev_ice(2)      /= 19 .OR.                                          &
         elev_ice(3)      /= 20 .OR.                                          &
         elev_ice(4)      /= 21 .OR.                                          &
         elev_ice(5)      /= 22 .OR.                                          &
         elev_ice(6)      /= 23 .OR.                                          &
         elev_ice(7)      /= 24 .OR.                                          &
         elev_ice(8)      /= 25 .OR.                                          &
         elev_ice(9)      /= 26 .OR.                                          &
         elev_ice(10)     /= 27 ) THEN
        
      ! Tile order does not match that expected
      errcode=1027
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  END IF
  
  IF (SIZE(rooh) /= ntype) THEN
    errcode = ntype
    cmessage='Changed ntype, initialisation arrays now incorrect'
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
  
  IF (SIZE(h2dd_c) /= npft) THEN
    errcode = npft
    cmessage='Changed npft, initialisation arrays now incorrect'
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
  
  SELECT CASE (ntype)
  CASE (9)
    IF (l_fix_improve_drydep) THEN
      rooh =        (/279.2,238.2,366.3,322.9,362.5,424.9,933.1,585.4,1156.1/)
    ELSE
      rooh =        (/30.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0/)
    END IF
     
    aerosol =     (/r_null,r_null,r_null,r_null,r_null,r_null,                &
                    1000.0,r_null,20000.0/)
     
    ch4_up_flux = (/39.5,50.0,30.0,37.0,27.5,0.0,0.0,27.5,0.0/)
     
    h2dd_c =      (/0.00197,0.00197,0.00177,1.2346,0.0001/)
     
    h2dd_m =      (/-0.00419,-0.00419,-0.00414,-0.472,0.0/)
  CASE (13,17,27)
    rooh(1:6)       = (/   300.3,   270.3,   266.9,   238.0,   238.5,   366.3/)
    aerosol(1:6)    = (/  r_null,  r_null,  r_null,  r_null,  r_null,  r_null/)
    ch4_up_flux(1:6)= (/    39.5,    39.5,    39.5,    50.0,    50.0,    30.0/)
    h2dd_c(1:6)     = (/ 0.00197, 0.00197, 0.00197, 0.00197, 0.00197, 0.00177/)
    h2dd_m(1:6)     = (/-0.00419,-0.00419,-0.00419,-0.00419,-0.00419,-0.00414/) 
  CASE DEFAULT
    errcode=(ntype*100)
    WRITE(cmessage,'(A)')                                                     &
         'UKCA does not handle flexible tiles yet: ' //                       &
         'Dry deposition needs extending. '          //                       &
         'Please use standard tile configuration'
    CALL ereport(RoutineName,errcode,cmessage)
  END SELECT
  
  SELECT CASE (ntype)
  CASE (13)
    rooh(7:13)       = (/   322.9, 332.8, 392.2, 424.9, 933.1, 585.4, 1156.1/)
    aerosol(7:13)    = (/  r_null,r_null,r_null,r_null,1000.0,r_null,20000.0/)
    ch4_up_flux(7:13)= (/    37.0,  27.5,  27.5,   0.0,   0.0,  27.5,    0.0/)
    h2dd_c(7:9)      = (/  1.2346,0.0001,0.0001 /)
    h2dd_m(7:9)      = (/-0.47200,   0.0,   0.0 /)
  CASE (17, 27)
    rooh(7:17)       = (/   366.3,   366.3,   322.9,   322.9,   322.9,        &
                            332.8,   392.2,   424.9,   933.1,   585.4,1156.1 /)
    aerosol(7:17)    = (/  r_null,  r_null,  r_null,  r_null,  r_null,        &
                           r_null,  r_null,  r_null,  1000.0,  r_null,20000.0/)
    ch4_up_flux(7:17)= (/    30.0,    30.0,    37.0,    37.0,    37.0,        &
                             27.5,    27.5,     0.0,     0.0,    27.5,0.0 /)
    h2dd_c(7:13)     = (/ 0.00177, 0.00177,  1.2346,  1.2346,  1.2346,        &
                           0.0001,  0.0001 /)
    h2dd_m(7:13)     = (/-0.00414,-0.00414,-0.47200,-0.47200,-0.47200,        &
                              0.0,     0.0 /)
  END SELECT
  
  IF (ntype == 27) THEN
    rooh(18:27)        = (/  1156.1, 1156.1, 1156.1, 1156.1, 1156.1,          &
                             1156.1, 1156.1, 1156.1, 1156.1, 1156.1 /)
    aerosol(18:27)     = (/ 20000.0,20000.0,20000.0,20000.0,20000.0,          &
                            20000.0,20000.0,20000.0,20000.0,20000.0 /)
    ch4_up_flux(18:27) = (/     0.0,    0.0,    0.0,    0.0,    0.0,          &
                                0.0,    0.0,    0.0,    0.0,    0.0 /)
  END IF
  
  SELECT CASE (ntype)
  CASE (9)
    ! Standard surface resistances (s m-1). Values are for 9 tiles in
    ! order: Broadleaved trees, Needleleaf trees, C3 Grass, C4 Grass,
    ! Shrub, Urban, Water, Bare Soil, Ice.
    !
    ! Nine tile dry deposition surface resistance have been updated April 2018
    ! (behind l_fix_improve_drydep logical) to reflect 13/17/27 tiles
    n_nosurf = 0
    DO n = 1, ndepd
      SELECT CASE (speci(nldepd(n)))
      CASE ('O3        ','O3S       ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 219.3, 233.0, 355.0, 309.3, 358.2,                    &
                        444.4,2000.0, 645.2,2000.0 /)
        ELSE
          rsurf(:,n)=(/ 200.0, 200.0, 200.0, 200.0, 400.0,                    &
                        800.0,2200.0, 800.0,2500.0 /)
        END IF
      CASE ('NO2       ','NO3       ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 364.1, 291.3, 443.8, 386.6, 447.8,                    &
                        555.6,2500.0, 806.5,2500.0 /)
        ELSE
          rsurf(:,n)=(/ 225.0, 225.0, 400.0, 400.0, 600.0,                    &
                       1200.0,2600.0,1200.0,3500.0 /)
        END IF
      CASE ('NO        ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 2184.5, 1747.6, 3662.7, 2319.6, 2686.8,               &
                        3333.3,15000.0, 4838.7,15000.0 /)
        ELSE
          rsurf(:,n)=(/ 1350.0, 1350.0, 2400.0, 2400.0, 3600.0,               &
                       72000.0, r_null,72000.0,21000.0 /)
        END IF
      CASE ('HNO3      ','HONO2     ','B2ndry    ','A2ndry    ','N2O5      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 8.5, 8.4, 13.2, 12.0, 62.3, 12.8, 13.9, 16.0, 19.4 /)
        ELSE
          rsurf(:,n)=tenpointzero
        END IF
      CASE ('HNO4      ','HO2NO2    ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 17.0, 16.8, 26.4, 24.0, 24.9, 25.7, 27.7, 32.1, 38.8 /)
        ELSE
          rsurf(:,n)=tenpointzero
        END IF
      CASE ('ISON      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 582.5, 466.0, 710.1, 618.6, 716.5,                    &
                        888.9,4000.0,1290.3,4000.0 /)
        ELSE
          rsurf(:,n)=tenpointzero
        END IF
      CASE ('HCl       ','HOCl      ','HBr       ','HOBr      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 8.5, 8.4, 13.2, 12.0, 62.3, 12.8, 13.9, 16.0, 19.4 /)
        ELSE
          rsurf(:,n)=zero
        END IF
      CASE ('H2SO4     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)= (/  84.9,  83.8, 131.9, 120.0, 124.4,                   &
                         128.5, 138.6, 160.4, 194.2 /)
        ELSE
          rsurf(:,n)=zero
        END IF
      CASE ('H2O2      ','HOOH      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 84.9,83.8,131.9,120.0,124.4,128.5,138.6,160.4,194.2 /)
        ELSE
          rsurf(:,n)=tenpointzero
        END IF
      CASE ('CH3OOH    ','MeOOH     ','C2H5OOH   ','EtOOH     ',              &
            'n_C3H7OOH ','i_C3H7OOH ','n-PrOOH   ','i-PrOOH   ',              &
            'MeCOCH2OOH','ISOOH     ','MACROOH   ','MeCO3H    ',              &
            'MeCO2H    ','HCOOH     ','PropeOOH  ','MEKOOH    ',              &
            'ALKAOOH   ','AROMOOH   ','BSVOC1    ','BSVOC2    ',              &
            'ASVOC1    ','ASVOC2    ','ISOSVOC1  ','ISOSVOC2  ',              &
            's-BuOOH   ','MVKOOH    ','HACET     ')
        rsurf(:,n)=rooh
      CASE ('PAN       ','PPAN      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 485.4, 388.4, 591.7, 515.5, 597.1,                    &
                        740.7,3333.3,1075.3,3333.3 /)
        ELSE
          rsurf(:,n)=(/ 500.0,  500.0, 500.0,  500.0, 500.0,                  &
                       r_null,12500.0, 500.0,12500.0 /)
        END IF
      CASE ('MPAN      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/  970.9,  776.7, 1183.4, 1030.9, 1194.2,               &
                        1481.5, 6666.7, 2150.5, 6666.7 /)
        ELSE
          rsurf(:,n)=(/ 500.0,  500.0, 500.0,  500.0, 500.0,                  &
                       r_null,12500.0, 500.0,12500.0 /)
        END IF
      CASE ('OnitU     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 582.5,  466.0, 710.1,  618.6, 716.5,                  &
                        888.9, 4000.0,1290.3, 4000.0 /)
        ELSE
          rsurf(:,n)=(/ 500.0,  500.0, 500.0,  500.0, 500.0,                  &
                       r_null,12500.0, 500.0,12500.0 /)
        END IF
      CASE ('NH3       ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 120.0, 130.9, 209.8, 196.1, 191.0,                    &
                        180.7, 148.9, 213.5, 215.1 /)
        ELSE
          rsurf(:,n)=tenpointzero
        END IF
      CASE ('CO        ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 3700.0, 7300.0, 4550.0, 1960.0, 4550.0,               &
                        r_null, r_null, 4550.0, r_null /)
        ELSE
          rsurf(:,n)=(/ 3700.0, 7300.0, 4550.0, 1960.0, 4550.0,               &
                        r_null, r_null, 4550.0, r_null /)
        END IF
        ! Shrub+bare soil set to C3 grass (guess)
      CASE ('CH4       ')
        rsurf(:,n)=zero
      CASE ('HONO      ')
        rsurf(:,n)=zero
      CASE ('H2        ')
        rsurf(:,n)=zero
      CASE ('SO2       ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 120.0, 130.9, 209.8, 196.1, 191.0,                    &
                        180.7, 148.9, 213.5, 215.1 /)
        ELSE
          rsurf(:,n)=(/ 100.0, 100.0, 150.0, 350.0, 400.0,                    &
                        400.0,  10.0, 700.0, r_null /)
        END IF
      CASE ('BSOA      ','ASOA      ','ISOSOA    ')
        rsurf(:,n)=aerosol
      CASE ('ORGNIT    ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 582.5, 466.0, 710.1, 618.6, 716.5,                    &
                        888.9,4000.0,1290.3,4000.0 /)
        ELSE
          rsurf(:,n)=aerosol
        END IF
      CASE ('Sec_Org   ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=aerosol
        ELSE
          rsurf(:,n)=zero
        END IF
      CASE ('HCHO      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 136.0, 143.5, 228.5, 211.6, 210.5,                    &
                        205.1, 182.7, 246.5, 261.8 /)
        ELSE
          rsurf(:,n)=(/ 100.0, 100.0, 150.0, 350.0, 600.0,                    &
                        400.0, 200.0, 700.0, 700.0 /)
        END IF
      CASE ('MeOH      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 187.1, 199.4, 318.3, 295.6, 292.2,                    &
                        282.1, 245.1, 337.2, 352.2 /)
        ELSE
          rsurf(:,n)=zero
        END IF
      CASE ('MeCHO     ','EtCHO     ','MACR      ','NALD      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 5825.2, 4660.3, 7100.6, 6185.6, 7164.8,               &
                        8888.9,40000.0,12903.2,40000.0 /)
        ELSE
          rsurf(:,n)=(/ 1200.0, 1200.0, 1200.0, 1200.0, 1000.0,               &
                        2400.0, r_null, r_null, r_null /)
        END IF
      CASE ('MGLY      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=(/ 12001.2, 13086.3, 20979.0, 19607.8, 19091.9,          &
                        18072.3, 14888.3, 21352.3, 21505.4 /)
        ELSE
          rsurf(:,n)=(/ 1200.0, 1200.0, 1200.0, 1200.0, 1000.0,               &
                        2400.0, r_null, r_null, r_null /)
        END IF
      CASE ('DMSO      ','Monoterp  ')
        rsurf(:,n)=zero
      CASE DEFAULT
        n_nosurf = n_nosurf + 1
        IF (first .AND. PrintStatus == PrStatus_Diag ) THEN
          umMessage = ' Surface resistance values not set for '//             &
                     speci(nldepd(n))
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      END SELECT
    END DO
    ! Fail if surface resistance values unavailable for some species
    IF (first .AND. n_nosurf > 0) THEN
      WRITE(cmessage,'(A,I0,A)')                                              &
        ' Surface resistance values not found for ',n_nosurf,' species.'
      errcode = ABS(n_nosurf)
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  CASE (13, 17, 27)
    ! Standard surface resistances (s m-1). Values are for 13/17/27 tiles in
    ! order: Broadleaf deciduous trees, Broadleaf evergreen tropical trees,
    !        Broadleaf evergreem temperate trees, Needleleaf deciduous trees, 
    !        Needleleaf evergreen trees, C3 Grass,
    !
    n_nosurf = 0 
    DO n = 1, ndepd
      SELECT CASE (speci(nldepd(n)))
      CASE ('O3        ','O3S       ')
        rsurf(1:6,n)=(/ 307.7,285.7,280.4,232.6,233.5,355.0 /)
      CASE ('SO2       ')
        rsurf(1:6,n)=(/ 137.0,111.1,111.9,131.3,130.4,209.8 /)
      CASE ('NO2       ','NO3       ')
        rsurf(1:6,n)=(/ 384.6,357.1,350.5,290.7,291.8,443.8 /)
      CASE ('NO        ')
        rsurf(1:6,n)=(/ 2307.7,2142.9,2102.8,1744.2,1751.0,3662.7 /)
      CASE ('HNO3      ','HONO2     ','B2ndry    ','A2ndry    ','N2O5      ')
        rsurf(1:6,n)=(/ 9.5,8.0,8.0,8.4,8.4,13.2 /)
      CASE ('HCl       ','HOCl      ','HBr       ','HOBr      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(1:6,n)=(/ 9.5,8.0,8.0,8.4,8.4,13.2 /)
        ELSE
          rsurf(:,n)=zero
        END IF
      CASE ('H2SO4     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(1:6,n)=(/ 94.8, 80.0, 80.0, 83.9, 83.7, 131.9 /)
        ELSE
          rsurf(:,n)=zero
        END IF
      CASE ('HNO4      ','HO2NO2    ')
        rsurf(1:6,n)=(/ 19.0,16.0,16.0,16.8,16.7,26.4 /)
      CASE ('HONO      ')
        rsurf(1:6,n)=(/ 47.4,40.0,40.0,42.0,41.8,65.9 /)
      CASE ('H2O2      ','HOOH      ')
        rsurf(1:6,n)=(/ 94.8,80.0,80.0,83.9,83.7,131.9 /)
      CASE ('PAN       ','PPAN      ')
        rsurf(1:6,n)=(/ 512.8,476.2,467.3,387.6,389.1,591.7 /)
      CASE ('MPAN      ')
        rsurf(1:6,n)=(/ 1025.6,952.4,934.6,775.2,778.2,1183.4 /)
      CASE ('OnitU     ','ISON      ','ORGNIT    ')
        rsurf(1:6,n)=(/ 615.4,571.4,560.7,465.1,466.9,710.1 /)
      CASE ('CH3OOH    ','MeOOH     ','C2H5OOH   ','EtOOH     ',              &
            'n_C3H7OOH ','i_C3H7OOH ','n-PrOOH   ','i-PrOOH   ',              &
            'MeCOCH2OOH','ISOOH     ','MACROOH   ','MeCO3H    ',              &
            'MeCO2H    ','HCOOH     ','PropeOOH  ','MEKOOH    ',              &
            'ALKAOOH   ','AROMOOH   ','BSVOC1    ','BSVOC2    ',              &
            'ASVOC1    ','ASVOC2    ','ISOSVOC1  ','ISOSVOC2  ',              &
            's-BuOOH   ','MVKOOH    ','HACET     ')
        rsurf(:,n)=rooh
      CASE ('NH3       ')
        rsurf(1:6,n)=(/ 137.0,111.1,111.9,131.3,130.4,209.8 /)
      CASE ('CO        ')
        rsurf(1:6,n)=(/ 3700.0,3700.0,3700.0,7300.0,7300.0,4550.0 /)  
        ! Shrub+bare soil set to C3 grass (guess)
      CASE ('CH4       ')
        rsurf(:,n)=zero
      CASE ('H2        ')
        rsurf(:,n)=zero
      CASE ('HCHO      ')
        rsurf(1:6,n)=(/ 154.1,126.6,127.2,143.8,143.1,228.5 /)
      CASE ('MeOH      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(1:6,n)=(/ 212.6, 173.9, 174.9, 200.0, 198.0, 318.3 /)
        ELSE
          rsurf(:,n)=zero
        END IF
      CASE ('MeCHO     ','EtCHO     ','MACR      ','NALD      ')
        rsurf(1:6,n)=(/ 6153.8,5714.3,5607.5,4651.2,4669.3,7100.6 /)
      CASE ('MGLY      ')
        rsurf(1:6,n)=(/ 13698.6,11111.1,11194.0,13129.1,13043.5,20979.0 /)
      CASE ('BSOA      ','ASOA      ','ISOSOA    ')
        rsurf(:,n)=aerosol
      CASE ('Sec_Org   ')
        IF (l_fix_improve_drydep) THEN
          rsurf(:,n)=aerosol
        ELSE
          rsurf(:,n)=zero
        END IF
      CASE ('DMSO      ','Monoterp  ')
        rsurf(:,n)=zero
      CASE DEFAULT
        n_nosurf = n_nosurf + 1
        IF (first .AND. PrintStatus == PrStatus_Diag ) THEN
          umMessage = ' Surface resistance values not set for '//             &
                        speci(nldepd(n))
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      END SELECT
    END DO
    ! Fail if surface resistance values unavailable for some species
    IF (first .AND. n_nosurf > 0) THEN
      WRITE(cmessage,'(A,I0,A)')                                              &
      ' Surface resistance values not found for ',n_nosurf,' species.'
      errcode = ABS(n_nosurf)
      CALL ereport(RoutineName,errcode,cmessage)
    END IF
  CASE DEFAULT
    errcode=319
    WRITE(cmessage,'(A)')                                                     &
         'UKCA does not handle flexible tiles yet: '//                        &
         'Dry deposition needs extending. '         //                        &
         'Please use standard tile configuration'
    CALL ereport(RoutineName,errcode,cmessage)
  END SELECT
  
  SELECT CASE (ntype)
  CASE (13)
    ! Standard surface resistances (s m-1). Values are for 13 tiles in
    ! order: C4 Grass, Shrub deciduous, Shrub evergreen,
    !        Urban, Water, Bare Soil, Ice.
    DO n = 1, ndepd
      SELECT CASE (speci(nldepd(n)))
      CASE ('O3        ','O3S       ')
        rsurf(7:13,n)=(/ 309.3,324.3,392.2,444.4,2000.0,645.2,2000.0 /)
      CASE ('SO2       ')
        rsurf(7:13,n)=(/ 196.1,185.8,196.1,180.7,148.9,213.5,215.1 /)
      CASE ('NO2       ','NO3       ')
        rsurf(7:13,n)=(/ 386.6,405.4,490.2,555.6,2500.0,806.5,2500.0 /)
      CASE ('NO        ')
        rsurf(7:13,n)=(/ 2319.6,2432.4,2941.2,3333.3,15000.0,4838.7,15000.0 /)
      CASE ('HNO3      ','HONO2     ','B2ndry    ','A2ndry    ','N2O5      ')
        rsurf(7:13,n)=(/ 12.0,59.1,65.4,12.8,13.9,16.0,19.4 /)
      CASE ('HCl       ','HOCl      ','HBr       ','HOBr      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:13,n)=(/ 12.0,59.1,65.4,12.8,13.9,16.0,19.4 /)
        END IF
      CASE ('H2SO4     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:13,n)=(/ 120.0, 118.1, 130.7, 128.5, 138.6, 160.4, 194.2 /)
        END IF
      CASE ('HNO4      ','HO2NO2    ')
        rsurf(7:13,n)=(/ 24.0,23.6,26.1,25.7,27.7,32.1,38.8 /)
      CASE ('HONO      ')
        rsurf(7:13,n)=(/ 60.0,59.1,65.4,64.2,69.3,80.2,97.1 /)
      CASE ('H2O2      ','HOOH      ')
        rsurf(7:13,n)=(/ 120.0,118.1,130.7,128.5,138.6,160.4,194.2 /)
      CASE ('PAN       ','PPAN      ')
        rsurf(7:13,n)=(/ 515.5,540.5,653.6,740.7,3333.3,1075.3,3333.3 /)
      CASE ('MPAN      ')
        rsurf(7:13,n)=(/ 1030.9,1081.1,1307.2,1481.5,6666.7,2150.5,6666.7 /)
      CASE ('OnitU     ','ISON      ','ORGNIT    ')
        rsurf(7:13,n)=(/ 618.6,648.6,784.3,888.9,4000.0,1290.3,4000.0 /)
      CASE ('NH3       ')
        rsurf(7:13,n)=(/ 196.1,185.8,196.1,180.7,148.9,213.5,215.1 /)
      CASE ('CO        ')
        rsurf(7:13,n)=(/ 1960.0,4550.0,4550.0,r_null,r_null,4550.0,r_null /)  
          ! Shrub+bare soil set to C3 grass (guess)
      CASE ('HCHO      ')
        rsurf(7:13,n)=(/ 211.6,203.1,217.9,205.1,182.7,246.5,261.8 /)
      CASE ('MeOH      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:13,n)=(/ 295.6, 282.7, 301.7, 282.1, 245.1, 337.2, 352.2 /)
        END IF
      CASE ('MeCHO     ','EtCHO     ','MACR      ','NALD      ')
        rsurf(7:13,n)=(/ 6185.6,6486.5,7843.1,8888.9,40000.0,12903.2,40000.0 /)
      CASE ('MGLY      ')
        rsurf(7:13,n)=(/ 19607.8,18575.9,19607.8,18072.3,14888.3,21352.3,     &
                         21505.4 /)
      END SELECT
    END DO
  CASE (17, 27)
    ! Standard surface resistances (s m-1). Values are for 17/27 tiles in
    ! order: C3 Crop, C3 Pasture, C4 Grass, C4 Crop, C4 Pasture,
    !        Shrub deciduous, Shrub evergreen, Urban, Water, Bare Soil, Ice.
    DO n = 1, ndepd
      SELECT CASE (speci(nldepd(n)))
      CASE ('O3        ','O3S       ')
        rsurf(7:17,n)=(/ 355.0,355.0,309.3, 309.3,309.3,                      &
                         324.3,392.2,444.4,2000.0,645.2,2000.0 /)
      CASE ('SO2       ')
        rsurf(7:17,n)=(/ 209.8,209.8,196.1,196.1,196.1,                       &
                         185.8,196.1,180.7,148.9,213.5,215.1 /)
      CASE ('NO2       ','NO3       ')
        rsurf(7:17,n)=(/ 443.8,443.8,386.6, 386.6,386.6,                      &
                         405.4,490.2,555.6,2500.0,806.5,2500.0 /)
      CASE ('NO        ')
        rsurf(7:17,n)=(/ 3662.7,3662.7,2319.6, 2319.6,2319.6,                 &
                         2432.4,2941.2,3333.3,15000.0,4838.7,15000.0 /)
      CASE ('HNO3      ','HONO2     ','B2ndry    ','A2ndry    ','N2O5      ')
        rsurf(7:17,n)=(/ 13.2,13.2,12.0,12.0,12.0,                            &
                         59.1,65.4,12.8,13.9,16.0,19.4 /)
      CASE ('HCl       ','HOCl      ','HBr       ','HOBr      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:17,n)=(/ 13.2,13.2,12.0,12.0,12.0,                          &
                           59.1,65.4,12.8,13.9,16.0,19.4 /)
        END IF
      CASE ('H2SO4     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:17,n)=(/ 131.9, 131.9, 120.0, 120.0, 120.0,                 &
                           118.1, 130.7, 128.5, 138.6, 160.4, 194.2 /)
        END IF
      CASE ('HNO4      ','HO2NO2    ')
        rsurf(7:17,n)=(/ 26.4,26.4,24.0,24.0,24.0,                            &
                         23.6,26.1,25.7,27.7,32.1,38.8 /)
      CASE ('HONO      ')
        rsurf(7:17,n)=(/ 65.9,65.9,60.0,60.0,60.0,                            &
                         59.1,65.4,64.2,69.3,80.2,97.1 /)
      CASE ('H2O2      ','HOOH      ')
        rsurf(7:17,n)=(/ 131.9,131.9,120.0,120.0,120.0,                       &
                         118.1,130.7,128.5,138.6,160.4,194.2 /)
      CASE ('PAN       ','PPAN      ')
        rsurf(7:17,n)=(/ 591.7,591.7,515.5, 515.5, 515.5,                     &
                         540.5,653.6,740.7,3333.3,1075.3,3333.3 /)
      CASE ('MPAN      ')
        rsurf(7:17,n)=(/ 1183.4,1183.4,1030.9,1030.9,1030.9,                  &
                         1081.1,1307.2,1481.5,6666.7,2150.5,6666.7 /)
      CASE ('OnitU     ','ISON      ','ORGNIT    ')
        rsurf(7:17,n)=(/ 710.1,710.1,618.6, 618.6, 618.6,                     &
                         648.6,784.3,888.9,4000.0,1290.3,4000.0 /)
      CASE ('NH3       ')
        rsurf(7:17,n)=(/ 209.8,209.8,196.1,196.1,196.1,                       &
                         185.8,196.1,180.7,148.9,213.5,215.1 /)
      CASE ('CO        ')
        rsurf(7:17,n)=(/ 4550.0,4550.0,1960.0,1960.0,1960.0,                  &
                         4550.0,4550.0,r_null,r_null,4550.0,r_null /)  
          ! Shrub+bare soil set to C3 grass (guess)
      CASE ('HCHO      ')
        rsurf(7:17,n)=(/ 228.5,228.5,211.6,211.6,211.6,                       &
                         203.1,217.9,205.1,182.7,246.5,261.8 /)
      CASE ('MeOH      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(7:17,n)=(/ 318.3, 318.3, 295.6, 295.6, 295.6,                 &
                           282.7, 301.7, 282.1, 245.1, 337.2, 352.2 /)
        END IF
      CASE ('MeCHO     ','EtCHO     ','MACR      ','NALD      ')
        rsurf(7:17,n)=(/ 7100.6,7100.6,6185.6, 6185.6, 6185.6,                &
                         6486.5,7843.1,8888.9,40000.0,12903.2,40000.0 /)
      CASE ('MGLY      ')
        rsurf(7:17,n)=(/ 20979.0,20979.0,19607.8,19607.8,19607.8,             &
                         18575.9,19607.8,18072.3,14888.3,21352.3,21505.4 /)
      END SELECT
    END DO
  END SELECT
  
  IF (ntype ==27) THEN
    ! Standard surface resistances (s m-1). Values are for 27 tiles in
    ! order: Elev_Ice(1-10).
    DO n = 1, ndepd
      SELECT CASE (speci(nldepd(n)))
      CASE ('O3        ','O3S       ')
        rsurf(18:27,n)=(/ 2000.0,2000.0,2000.0,2000.0,2000.0,                 &
                          2000.0,2000.0,2000.0,2000.0,2000.0 /)
      CASE ('SO2       ')
        rsurf(18:27,n)=(/ 215.1,215.1,215.1,215.1,215.1,                      &
                          215.1,215.1,215.1,215.1,215.1 /)
      CASE ('NO2       ','NO3       ')
        rsurf(18:27,n)=(/ 2500.0,2500.0,2500.0,2500.0,2500.0,                 &
                          2500.0,2500.0,2500.0,2500.0,2500.0 /)
      CASE ('NO        ')
        rsurf(18:27,n)=(/ 15000.0,15000.0,15000.0,15000.0,15000.0,            &
                          15000.0,15000.0,15000.0,15000.0,15000.0 /)
      CASE ('HNO3      ','HONO2     ','B2ndry    ','A2ndry    ','N2O5      ')
        rsurf(18:27,n)=(/ 19.4,19.4,19.4,19.4,19.4,                           &
                          19.4,19.4,19.4,19.4,19.4 /)
      CASE ('HCl       ','HOCl      ','HBr       ','HOBr      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(18:27,n)=(/ 19.4,19.4,19.4,19.4,19.4,                         &
                            19.4,19.4,19.4,19.4,19.4 /)
        END IF
      CASE ('H2SO4     ')
        IF (l_fix_improve_drydep) THEN
          rsurf(18:27,n)=(/ 194.2, 194.2, 194.2, 194.2, 194.2,                &
                            194.2, 194.2, 194.2, 194.2, 194.2 /)
        END IF
      CASE ('HNO4      ','HO2NO2    ')
        rsurf(18:27,n)=(/ 38.8,38.8,38.8,38.8,38.8,                           &
                          38.8,38.8,38.8,38.8,38.8 /)
      CASE ('HONO      ')
        rsurf(18:27,n)=(/ 97.1,97.1,97.1,97.1,97.1,                           &
                          97.1,97.1,97.1,97.1,97.1 /)
      CASE ('H2O2      ','HOOH      ')
        rsurf(18:27,n)=(/ 194.2,194.2,194.2,194.2,194.2,                      &
                          194.2,194.2,194.2,194.2,194.2 /)
      CASE ('PAN       ','PPAN      ')
        rsurf(18:27,n)=(/ 3333.3,3333.3,3333.3,3333.3,3333.3,                 &
                          3333.3,3333.3,3333.3,3333.3,3333.3 /)
      CASE ('MPAN      ')
        rsurf(18:27,n)=(/ 6666.7,6666.7,6666.7,6666.7,6666.7,                 &
                          6666.7,6666.7,6666.7,6666.7,6666.7 /)
      CASE ('OnitU     ','ISON      ','ORGNIT    ')
        rsurf(18:27,n)=(/ 4000.0,4000.0,4000.0,4000.0,4000.0,                 &
                          4000.0,4000.0,4000.0,4000.0,4000.0 /)
      CASE ('NH3       ')
        rsurf(18:27,n)=(/ 215.1,215.1,215.1,215.1,215.1,                      &
                          215.1,215.1,215.1,215.1,215.1 /)
      CASE ('CO        ')
        rsurf(18:27,n)=(/ r_null,r_null,r_null,r_null,r_null,                 &
                          r_null,r_null,r_null,r_null,r_null /)  
          ! Shrub+bare soil set to C3 grass (guess)
      CASE ('HCHO      ')
        rsurf(18:27,n)=(/ 261.8,261.8,261.8,261.8,261.8,                      &
                          261.8,261.8,261.8,261.8,261.8 /)
      CASE ('MeOH      ')
        IF (l_fix_improve_drydep) THEN
          rsurf(18:27,n)=(/ 352.2, 352.2, 352.2, 352.2, 352.2,                &
                            352.2, 352.2, 352.2, 352.2, 352.2 /)
        END IF
      CASE ('MeCHO     ','EtCHO     ','MACR      ','NALD      ')
        rsurf(18:27,n)=(/ 40000.0,40000.0,40000.0,40000.0,40000.0,            &
                          40000.0,40000.0,40000.0,40000.0,40000.0 /)
      CASE ('MGLY      ')
        rsurf(18:27,n)=(/ 21505.4,21505.4,21505.4,21505.4,21505.4,            &
                          21505.4,21505.4,21505.4,21505.4,21505.4 /)
      END SELECT
    END DO
  END IF
  first = .FALSE.
END IF !first

o3_stom_frac = 0.0

!     rco3global = 0.0
!     gsfglobal = 0.0


!     Set logical for surface types

DO n = 1, ntype
  DO k = 1, rows
    DO i = 1, row_length
      todo(i,k,n) = (gsf(i,k,n) > 0.0)
    END DO
  END DO
END DO

! Set microb
microb(:,:) = (rh(:,:) > 0.4 .AND. t0(:,:) > 278.0)

!     Set surface resistances to standard values. rsurf is the
!     resistance of the soil, rock, water etc. Set all tiles to
!     standard values. These values will be modified below as
!     necessary. Extra terms for vegetated tiles (stomatal, cuticular)
!     will be added if required. Loop over all parts of array rc to
!     ensure all of it is assigned a value.

DO n = 1, ntype
  DO k = 1, rows
    DO i = 1, row_length
      IF (todo(i,k,n)) THEN
        DO j = 1, ndepd
          rc(i,k,n,j) = rsurf(n,j)
        END DO
      ELSE
        DO j = 1, ndepd
          rc(i,k,n,j) = r_null
        END DO
      END IF
    END DO
  END DO
END DO

!     Calculate stomatal resistances

DO n = 1, npft
  DO k = 1, rows
    DO i = 1, row_length
      IF (todo(i,k,n) .AND. stcon(i,k,n) > glmin) THEN
        r_stom(i,k,n,1) = 1.5 * dif(1) / stcon(i,k,n) ! NO2
        r_stom(i,k,n,2) = dif(2) / stcon(i,k,n)       ! O3
        r_stom(i,k,n,3) = dif(3) / stcon(i,k,n)       ! PAN
        r_stom(i,k,n,4) = dif(4) / stcon(i,k,n)       ! SO2
        r_stom(i,k,n,5) = dif(5) / stcon(i,k,n)       ! NH3
      ELSE
        r_stom(i,k,n,1) = r_null
        r_stom(i,k,n,2) = r_null
        r_stom(i,k,n,3) = r_null
        r_stom(i,k,n,4) = r_null
        r_stom(i,k,n,5) = r_null
      END IF
    END DO
  END DO
END DO

!     Now begin assigning specific surface resistances.

DO j = 1, ndepd

  !       O3: Change land deposition values if surface is wet;
  !       soil moisture value of 0.3 fairly arbitrary.

  IF (speci(nldepd(j)) == 'O3        ') THEN
    DO k = 1, rows
      DO i = 1, row_length
        IF (smr(i,k) > 0.3) THEN
          DO n = 1, npft
            IF (todo(i,k,n)) THEN
              rc(i,k,n,j) = r_wet_o3
            END IF
          END DO
          IF (todo(i,k,soil)) THEN
            rc(i,k,soil,j) = r_wet_o3
          END IF
        END IF
      END DO

      !           Change values for tundra regions

      n = npft
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,n)) rc(i,k,n,j)  = r_tundra(3)
          IF (todo(i,k,soil)) rc(i,k,soil,j)  = r_tundra(3)
        END IF
      END DO

      !           Cuticular resistance for ozone.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n) .AND. lai_ft(i,k,n) > 0.0) THEN
            r_cut_o3(i,n) = cuticle_o3 / lai_ft(i,k,n)
          ELSE
            r_cut_o3(i,n) = r_null
          END IF
        END DO
      END DO

      !           Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,2)) +                                      &
                 (1.0/r_cut_o3(i,n)) +                                        &
                 (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
            o3_stom_frac(i,k) = o3_stom_frac(i,k) +                           &
              gsf(i,k,n) / (rr * r_stom(i,k,n,2))
          END IF
        END DO
      END DO
    END DO

    !       NO2

  ELSE IF (speci(nldepd(j)) == 'NO2       ') THEN

    DO k = 1, rows

      !           Change values for tundra regions

      n = npft
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,n)) rc(i,k,n,j)  = r_tundra(1)
          IF (todo(i,k,soil)) rc(i,k,soil,j)  = r_tundra(1)
        END IF
      END DO

      !           Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,1)) + (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
          END IF
        END DO
      END DO
    END DO

    !       PAN

  ELSE IF (speci(nldepd(j)) == 'PAN       ') THEN

    DO k = 1, rows

      !           Change values for tundra regions

      n = npft
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,n)) rc(i,k,n,j)  = r_tundra(5)
          IF (todo(i,k,soil)) rc(i,k,soil,j)  = r_tundra(5)
        END IF
      END DO

      !           Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,3)) + (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
          END IF
        END DO
      END DO
    END DO

    !         PPAN

  ELSE IF (speci(nldepd(j)) == 'PPAN      ') THEN

    DO k = 1, rows

      !               Change values for tundra regions

      n = npft
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,n)) rc(i,k,n,j)  = r_tundra(5)
          IF (todo(i,k,soil)) rc(i,k,soil,j)  = r_tundra(5)
        END IF
      END DO

      !               Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,3)) + (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
          END IF
        END DO
      END DO
    END DO

    !         MPAN

  ELSE IF (speci(nldepd(j)) == 'MPAN      ') THEN

    DO k = 1, rows

      !             Change values for tundra regions

      n = npft
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,n)) rc(i,k,n,j)  = r_tundra(5)
          IF (todo(i,k,soil)) rc(i,k,soil,j)  = r_tundra(5)
        END IF
      END DO

      !             Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,3)) + (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
          END IF
        END DO
      END DO
    END DO

    !         ONITU

  ELSE IF (speci(nldepd(j)) == 'ONITU     ') THEN

    DO k = 1, rows

      !             Change values for tundra regions

      n = npft
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,n)) rc(i,k,n,j)  = r_tundra(5)
          IF (todo(i,k,soil)) rc(i,k,soil,j)  = r_tundra(5)
        END IF
      END DO

      !             Calculate plant deposition terms.

      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n)) THEN
            rr = (1.0/r_stom(i,k,n,3)) + (1.0/rc(i,k,n,j))
            rc(i,k,n,j) = 1.0 / rr
          END IF
        END DO
      END DO
    END DO

    !       H2 dry dep vel has linear dependence on soil moisture
    !       Limit sm to avoid excessively high deposition velocities

  ELSE IF (speci(nldepd(j)) == 'H2        ') THEN
    DO k = 1, rows
      DO n = 1, npft-2
        DO i = 1, row_length
          IF (todo(i,k,n) .AND. microb(i,k)) THEN
            sm = MAX(smr(i,k),0.1)
            rc(i,k,n,j) = 1.0 / (h2dd_m(n) * sm + h2dd_c(n))
          END IF
        END DO
      END DO

      !           C4 grass: H2 dry dep has quadratic-log dependence on
      !           soil moisture

      n = npft - 1
      DO i = 1, row_length
        IF (todo(i,k,n) .AND. microb(i,k)) THEN
          sm = LOG(MAX(smr(i,k),0.1))
          rr = (h2dd_c(4) + sm*(h2dd_m(4) + sm*h2dd_q)) * 1.0e-4
          IF (rr > 0.00131) rr = 0.00131 ! Conrad/Seiler Max value
          rc(i,k,n,j) = 1.0 / rr
        END IF
      END DO

      !           Shrub: H2 dry dep velocity has no dependence on
      !           soil moisture

      n = npft
      rr = 1.0 / h2dd_c(n)
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit) THEN
          IF (todo(i,k,n) .AND. microb(i,k)) THEN
            rc(i,k,n,j) = rr
            rc(i,k,soil,j) = rr
          END IF
        END IF
      END DO

      !           Bare soil

      DO i = 1, row_length
        IF (sinlat(i,k) <= tundra_s_limit) THEN
          IF (todo(i,k,n) .AND. microb(i,k)) THEN
            sm = MAX(smr(i,k),0.1)
            rc(i,k,n,j) = 1.0 / (h2dd_m(3) * sm + h2dd_c(3))
            rc(i,k,soil,j) = 1.0 / (h2dd_m(3) * sm + h2dd_c(3))
          END IF
        END IF
      END DO
    END DO

    !       CH4

  ELSE IF (speci(nldepd(j)) == 'CH4       ') THEN

    !        CH4: Calculate an uptake flux initially. The uptake
    !        flux depends on soil moisture, based on results of
    !        Reay et al. (2001). Do the PFTs Broadleaf, Needleleaf,
    !        C3 and C4 grasses first.

    DO k = 1, rows
      DO i = 1, row_length
        IF (microb(i,k)) THEN
          sm = smr(i,k)
          IF (sm < 0.16) THEN
            f = sm / 0.16
          ELSE IF (sm > 0.30) THEN
            f = (0.50 - sm) / 0.20
          ELSE
            f = 1.0
          END IF
          f = MAX(f,0.0)
          IF (todo(i,k,1)) rc(i,k,1,j) = ch4_up_flux(1) * f
          IF (todo(i,k,2)) rc(i,k,2,j) = ch4_up_flux(2) * f
          IF (todo(i,k,3)) rc(i,k,3,j) = ch4_up_flux(3) * f
          IF (todo(i,k,4)) rc(i,k,4,j) = ch4_up_flux(4) * f
        END IF
      END DO
    END DO

    !         Now do shrub and bare soil, assumed to be tundra if
    !         latitude > 60N; Ctherwise, calculate an uptake flux
    !         initially. The uptake flux depends on soil moisture,
    !         based on results of Reay et al. (2001).

    n = npft
    DO k = 1, rows
      DO i = 1, row_length
        IF (microb(i,k)) THEN
          IF (sinlat(i,k) > tundra_s_limit) THEN
            ts = t0tile(i,k,n)
            rr = ch4dd_tun(4) + ts * (ch4dd_tun(3) +                          &
              ts * (ch4dd_tun(2) + ts * ch4dd_tun(1)))
            rr = rr * 3600.0 ! Convert from s-1 to h-1
            IF (todo(i,k,n)) rc(i,k,n,j) = MAX(rr,0.0)
            IF (todo(i,k,soil)) rc(i,k,soil,j) = MAX(rr,0.0)
          ELSE
            sm = smr(i,k)
            IF (sm < 0.16) THEN
              f = sm / 0.16
            ELSE IF (sm > 0.30) THEN
              f = (0.50 - sm) / 0.20
            ELSE
              f = 1.0
            END IF
            f = MAX(f,0.0)
            IF (todo(i,k,n)) rc(i,k,n,j) = ch4_up_flux(n) * f
            IF (todo(i,k,soil)) rc(i,k,soil,j) =                              &
              ch4_up_flux(soil) * f
          END IF
        END IF
      END DO
    END DO

    !         Convert CH4 uptake fluxes (ug m-2 h-1) to
    !         resistance (s m-1).

    DO k = 1, rows
      DO n = 1, npft
        DO i = 1, row_length
          IF (todo(i,k,n) .AND. microb(i,k)) THEN
            rr = rc(i,k,n,j)
            IF (rr > 0.0) THEN
              rc(i,k,n,j) = p0(i,k) * mml /                                   &
                (rmol * t0tile(i,k,n) * rr)
            ELSE
              rc(i,k,n,j) = r_null
            END IF
          END IF
        END DO
      END DO
      n = soil
      DO i = 1, row_length
        IF (todo(i,k,n) .AND. microb(i,k)) THEN
          rr = rc(i,k,n,j)
          IF (rr > 0.0) THEN
            rc(i,k,n,j) = p0(i,k) * mml /                                     &
              (rmol * t0tile(i,k,n) * rr * TAR_scaling)
          ELSE
            rc(i,k,n,j) = r_null
          END IF
        END IF
      END DO
    END DO

  ELSE IF (speci(nldepd(j)) == 'CO        ') THEN

    !         Only assign values for CO if microbes are active.

    n = npft
    DO k = 1, rows
      DO i = 1, row_length
        IF (sinlat(i,k) > tundra_s_limit .AND. microb(i,k)) THEN
          IF (todo(i,k,n)) THEN
            rc(i,k,n,j)  = r_tundra(2)       ! CO
          END IF
          IF (todo(i,k,soil)) THEN
            rc(i,k,soil,j)  = r_tundra(2)    ! CO
          END IF
        END IF
      END DO
    END DO

    !       HONO2

    !       Calculate resistances for HONO2 deposition to ice, which
    !       depend on temperature. Ensure resistance for HONO2 does not fall
    !       below 10 s m-1.

  ELSE IF ((speci(nldepd(j)) == 'HNO3      ' .OR.                             &
            speci(nldepd(j)) == 'HONO2     ' .OR.                             &
            speci(nldepd(j)) == 'ISON      ') .OR.                            &
          ((speci(nldepd(j)) == 'HCl       ' .OR.                             &
            speci(nldepd(j)) == 'HOCl      ' .OR.                             &
            speci(nldepd(j)) == 'HBr       ' .OR.                             &
            speci(nldepd(j))== 'HOBr      ') .AND. (l_fix_improve_drydep))) THEN

    n = ntype
    DO k = 1, rows
      DO i = 1, row_length
        IF (todo(i,k,n)) THEN

          !               Limit temperature to a minimum of 252K. Curve used
          !               only used data between 255K and 273K.

          ts = MAX(t0tile(i,k,n), 252.0)
          rr = hno3dd_ice(3) + ts * (hno3dd_ice(2) +                          &
               ts * hno3dd_ice(1))
          rc(i,k,n,j) = MAX(rr,10.0)
        END IF
      END DO
    END DO

    !       ORGNIT is treated as an aerosol.
    !       Assume vd is valid for all land types and aerosol types.

  ELSE IF ((speci(nldepd(j)) == 'ORGNIT    '  .OR.                            &
            speci(nldepd(j)) == 'BSOA      '  .OR.                            &
            speci(nldepd(j)) == 'ASOA      '  .OR.                            &
            speci(nldepd(j)) == 'ISOSOA    ') .OR.                            &
           (speci(nldepd(j)) =='Sec_Org   ' .AND. (l_fix_improve_drydep))) THEN

    DO n = 1, npft+1    !for the 5 functional plant types as
                        !well as for surface type urban (n=6)
      DO k = 1, rows
        DO i = 1, row_length
          IF (so4_vd(i,k) > 0.0 .AND. gsf(i,k,n) > 0.0) THEN
            rr          = 1.0 / so4_vd(i,k)
            rc(i,k,n,j) = rr
          END IF
        END DO
      END DO
    END DO

    n = soil            !for surface type soil (n=8)
    DO k = 1, rows
      DO i = 1, row_length
        IF (so4_vd(i,k) > 0.0 .AND. gsf(i,k,n) > 0.0) THEN
          rr          = 1.0 / so4_vd(i,k)
          rc(i,k,n,j) = rr
        END IF
      END DO
    END DO

  END IF ! End of IF (speci == species name)

END DO  ! End of DO j = 1, ndepd
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_surfddr
END MODULE ukca_surfddr_mod

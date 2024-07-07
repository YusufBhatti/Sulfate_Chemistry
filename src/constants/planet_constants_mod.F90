! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Planet constants

MODULE planet_constants_mod

! Description:
!   Physical constants for a general planet to be read from a namelist.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Constants

USE missing_data_mod, ONLY: imdi, rmdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

! Universal constants
USE conversions_mod, ONLY: pi, rsec_per_day
USE water_constants_mod, ONLY: lc, lf
USE errormessagelength_mod, ONLY: errormessagelength
USE um_types, ONLY: real32

IMPLICIT NONE

!----------------------------------------------------------------------
! Primary planet constants
!----------------------------------------------------------------------

! Use planet orbital parameters
LOGICAL :: l_planet_orbit = .FALSE.

! Set planet rotation from namelist
LOGICAL :: l_set_planet_rotation = .FALSE.

! Use g on model levels for radiative heating
LOGICAL :: l_planet_g = .FALSE.

! Set a grey surface for radiation
LOGICAL :: l_planet_grey_surface = .FALSE.

! Use an intrinsic thermal flux at the lower boundary, set
! with an effective intrinsic temperature planet_t_intrinsic.
LOGICAL :: l_planet_intrinsic_flux = .FALSE.

! Use a constant aerosol mixing ratio
LOGICAL :: l_planet_aerosol = .FALSE.

! Fix the sun at a particular zenith and azimuth angle
LOGICAL :: l_fix_solang = .FALSE.

! Intrinsic planet temperature used to set thermal flux at lower
! boundary if l_planet_intrinsic_flux is .true.
REAL :: planet_t_intrinsic = rmdi

! Effective surface albedo for broadband shortwave radiation
! if l_planet_grey_surface is .true.
REAL :: planet_albedo = rmdi

! Effective surface emissivity for broadband longwave radiation
! if l_planet_grey_surface is .true.
REAL :: planet_emissivity = rmdi

! Constant aerosol mixing ratio
REAL :: planet_aerosol_mmr = rmdi

! Fixed solar zenith angle in radians
REAL :: solar_zenith_angle = rmdi

! Fixed solar azimuth angle measured clockwise from grid north (radians)
REAL :: solar_azimuth_angle = rmdi

! Planet radius in metres
REAL :: planet_radius = rmdi

! Epoch in Julian Days
REAL :: planet_epoch = rmdi

! Eccentricity of the orbit
REAL :: planet_e = rmdi

! Increment to eccentricity per day number from epoch
REAL :: planet_de = rmdi

! Longitude of perihelion in radians
REAL :: planet_lph = rmdi

! Increment to longitude of perihelion per day number from epoch
REAL :: planet_dlph = rmdi

! Obliquity of the orbit in radians
REAL :: planet_oblq = rmdi

! Increment to obliquity of the orbit per day number from epoch
REAL :: planet_doblq = rmdi

! Semi-major axis in AU
REAL :: planet_a = rmdi

! Increment to semi-major axis per day number from epoch
REAL :: planet_da = rmdi

! Mean anomaly at epoch in radians
REAL :: planet_M = rmdi

! Increment to mean anomaly per day number from epoch
REAL :: planet_dM = rmdi

! Planet hour angle at epoch in radians
REAL :: planet_ha = rmdi

! Increment to planet hour angle per day number from epoch
REAL :: planet_dha = rmdi

! Orbital longitude of observer
REAL :: planet_obs_lon = rmdi

! Orbital latitude of observer
REAL :: planet_obs_lat = rmdi

! Stellar radius in metres
REAL :: stellar_radius = rmdi

! Stellar irradiance at 1 astronomical unit (AU) in W/m2
REAL :: sc = rmdi

! Mean acceleration due to gravity at the planet surface
REAL :: g = rmdi

! Gas constant for dry air
REAL :: r = rmdi

! Specific heat of dry air at constant pressure
REAL :: cp = rmdi

! Reference surface pressure
REAL :: pref = rmdi

! Mean scale height for pressure
REAL :: sclht = rmdi

! Near surface environmental lapse rate
REAL :: lapse = rmdi

! Angular speed of planet rotation
REAL :: omega = rmdi

! Selector for formulation of equation of time
INTEGER :: i_eqt = imdi

! Identifiers for equation of time formulation:
INTEGER, PARAMETER :: ip_smart   = 1
!   Equation 29 on page 149 of Smart (1944).
INTEGER, PARAMETER :: ip_mueller = 2 
!   M. Mueller, Acta Physica Polonica A 88 Supplement, S-49 (1995)

!----------------------------------------------------------------------
! Derived planet constants
!----------------------------------------------------------------------

! Angular speed of planet rotation x2
REAL :: two_omega

! Seconds-to-radians converter
REAL :: s2r                 ! planet_dha/rsec_per_day

! Ratio of molecular weights of water and dry air
REAL :: repsilon            ! r/rv
REAL (KIND=real32) :: repsilon_32b

REAL :: p_zero              ! pref
REAL :: recip_p_zero        ! 1.0/pref
REAL :: kappa               ! r/cp
REAL :: recip_kappa         ! 1.0/kappa
REAL :: recip_epsilon       ! 1.0/repsilon
REAL :: c_virtual           ! 1.0/repsilon-1.0
REAL :: one_minus_epsilon   ! 1.0-repsilon
REAL (KIND=real32) :: one_minus_epsilon_32b
REAL :: etar                ! 1.0/(1.0-repsilon)
REAL :: grcp                ! g/cp
REAL :: lcrcp               ! lc/cp
REAL :: lfrcp               ! lf/cp
REAL :: ls                  ! lc+lf
REAL :: lsrcp               ! (lc+lf)/cp
REAL :: cv                  ! cp-r
REAL :: recip_a2            ! 1.0/(planet_radius*planet_radius)
REAL :: g_over_r            ! g/r


!----------------------------------------------------------------------
! Preset planet constants
!----------------------------------------------------------------------

! Selector for preset planet constants
INTEGER :: i_planet = imdi

! Planet identifiers
INTEGER, PARAMETER :: ip_user       = 0
INTEGER, PARAMETER :: ip_mercury    = 1
INTEGER, PARAMETER :: ip_venus      = 2
INTEGER, PARAMETER :: ip_earth      = 3
INTEGER, PARAMETER :: ip_mars       = 4
INTEGER, PARAMETER :: ip_jupiter    = 5
INTEGER, PARAMETER :: ip_saturn     = 6
INTEGER, PARAMETER :: ip_uranus     = 7
INTEGER, PARAMETER :: ip_neptune    = 8
INTEGER, PARAMETER :: ip_pluto      = 9


!----------------------------------------------------------------------
! Earth constants
!----------------------------------------------------------------------

! Earth radius in metres
REAL, PRIVATE, PARAMETER :: earth_radius   = 6371229.0

! Mean acceleration due to gravity at the Earth's surface
REAL, PRIVATE, PARAMETER :: earth_g        = 9.80665

! Gas constant for dry air
REAL, PRIVATE, PARAMETER :: earth_r        = 287.05

! Specific heat of dry air at constant pressure
REAL, PRIVATE, PARAMETER :: earth_cp       = 1005.0

! Reference surface pressure
REAL, PRIVATE, PARAMETER :: earth_pref     = 100000.0

! Mean scale height for pressure
REAL, PRIVATE, PARAMETER :: earth_sclht    = 6.8e+03

! repsilon, ratio of molecular weights of water and dry air
REAL, PRIVATE, PARAMETER :: earth_repsilon = 0.62198

! Near surface environmental lapse rate
REAL, PRIVATE, PARAMETER :: earth_lapse    = 0.0065

! Angular speed of planet rotation (2*pi/siderial day: 23h56m04s)
REAL, PRIVATE, PARAMETER :: earth_omega    = 7.292116e-5

! Increment to Earth's hour angle per day number from epoch
REAL, PRIVATE, PARAMETER :: earth_dha      = 2.0*pi

! Tropopause lapse rate
REAL, PARAMETER :: lapse_trop              = 0.002

! Von Karman's constant
REAL, PARAMETER :: vkman                   = 0.4

! Rv - the gas constant for water vapour
! Note repsilon = r/rv  so to be consistent with above definitions rv=r/repsilon
! Previously convection was using rv=461.1 whereas now rv~461.51 (J/kg/K)
REAL, PARAMETER :: rv                      = earth_r/earth_repsilon

! Solar radius in metres (IAU defined nominal solar radius)
REAL, PRIVATE, PARAMETER :: solar_radius   = 6.957e+08

!----------------------------------------------------------------------
! Planet constants namelist
!----------------------------------------------------------------------

NAMELIST/planet_constants/ i_planet, i_eqt, l_planet_orbit, &
  l_set_planet_rotation, l_planet_g, l_planet_grey_surface, &
  l_planet_intrinsic_flux, l_planet_aerosol, l_fix_solang, &
  planet_t_intrinsic, planet_albedo, planet_emissivity, &
  planet_aerosol_mmr, solar_zenith_angle, solar_azimuth_angle, &
  planet_radius, planet_epoch, &
  planet_e, planet_lph, planet_oblq, planet_a, planet_M, &
  planet_de, planet_dlph, planet_doblq, planet_da, planet_dM, &
  planet_ha, planet_dha, planet_obs_lon, planet_obs_lat, &
  stellar_radius, sc, g, r, cp, pref, omega, sclht, lapse


! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PLANET_CONSTANTS_MOD'

CONTAINS

SUBROUTINE set_planet_constants()

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

INTEGER                              :: icode         ! used for ereport
CHARACTER (LEN=errormessagelength)   :: cmessage      ! used for ereport
CHARACTER (LEN=*), PARAMETER         :: RoutineName = 'set_planet_constants'

SELECT CASE (i_planet)
CASE (ip_user)
  IF (.NOT.l_set_planet_rotation) THEN
    omega = 0.0
  END IF
CASE (ip_earth)
  planet_radius = earth_radius
  g = earth_g
  r = earth_r
  cp = earth_cp
  pref = earth_pref
  sclht = earth_sclht
  lapse = earth_lapse
  repsilon = earth_repsilon
  IF (.NOT.l_set_planet_rotation) THEN
    omega = earth_omega
  END IF
  planet_dha = earth_dha
  stellar_radius = solar_radius
CASE DEFAULT
  WRITE(cmessage,'(a,i2,a)') 'Planet ', i_planet, ' not yet supported'
  icode = 1
  CALL ereport(RoutineName, icode, cmessage)
END SELECT

! Set derived constants:
IF (i_planet /= ip_earth) THEN
  repsilon = r/rv
END IF
two_omega     = 2.0*omega
s2r           = planet_dha/rsec_per_day
p_zero        = pref
recip_p_zero  = 1.0/pref
kappa         = r/cp
recip_kappa   = 1.0/kappa
recip_epsilon = 1.0/repsilon
c_virtual     = 1.0/repsilon-1.0
one_minus_epsilon = 1.0-repsilon
etar          = 1.0/(1.0-repsilon)
grcp          = g/cp
lcrcp         = lc/cp
lfrcp         = lf/cp
ls            = lc+lf
lsrcp         = (lc+lf)/cp
cv            = cp-r
recip_a2      = 1.0/(planet_radius*planet_radius)
g_over_r      = g/r

!Set 32-bit versions as required, eg in qsat_mod
repsilon_32b          = REAL(repsilon,real32)
one_minus_epsilon_32b = REAL(one_minus_epsilon,real32)

END SUBROUTINE set_planet_constants


SUBROUTINE print_nlist_planet_constants()

USE umPrintMgr, ONLY: umPrint, maxLineLen

IMPLICIT NONE
CHARACTER(LEN=maxLineLen) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_PLANET_CONSTANTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist planet_constants',                 &
    src='planet_constants_mod')

WRITE(lineBuffer,'(A,I0)') '  i_planet= ',i_planet
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,I0)') '  i_eqt= ',i_eqt
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_t_intrinsic= ',planet_t_intrinsic
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_albedo= ',planet_albedo
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_emissivity= ',planet_emissivity
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_aerosol_mmr= ',planet_aerosol_mmr
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  solar_zenith_angle= ',solar_zenith_angle
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  solar_azimuth_angle= ',solar_azimuth_angle
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_radius= ',planet_radius
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_epoch= ',planet_epoch
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_e= ',planet_e
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_de= ',planet_de
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_lph= ',planet_lph
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_dlph= ',planet_dlph
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_oblq= ',planet_oblq
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_doblq= ',planet_doblq
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_a= ',planet_a
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_da= ',planet_da
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_M= ',planet_M
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_dM= ',planet_dM
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_ha= ',planet_ha
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_dha= ',planet_dha
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_obs_lon= ',planet_obs_lon
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  planet_obs_lat= ',planet_obs_lat
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  sc= ',sc
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  g= ',g
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  r= ',r
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  cp= ',cp
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  pref= ',pref
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  sclht= ',sclht
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  lapse= ',lapse
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,ES15.7)') '  omega= ',omega
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,L1)') '  l_planet_orbit= ',l_planet_orbit
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,L1)') '  l_set_planet_rotation= ',l_set_planet_rotation
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,L1)') '  l_planet_g= ',l_planet_g
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,L1)') '  l_planet_grey_surface= ',l_planet_grey_surface
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,L1)') ' l_planet_intrinsic_flux = ',l_planet_intrinsic_flux
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,L1)') '  l_planet_aerosol= ',l_planet_aerosol
CALL umPrint(lineBuffer,src='planet_constants_mod')
WRITE(lineBuffer,'(A,L1)') '  l_fix_solang= ',l_fix_solang
CALL umPrint(lineBuffer,src='planet_constants_mod')
CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='planet_constants_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)


END SUBROUTINE print_nlist_planet_constants


#if !defined(LFRIC)
SUBROUTINE read_nml_planet_constants(unit_in)

USE um_parcore, ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_PLANET_CONSTANTS'

! Set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 2
INTEGER, PARAMETER :: n_real = 31
INTEGER, PARAMETER :: n_log = 7

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_planet
  INTEGER :: i_eqt
  REAL :: planet_t_intrinsic
  REAL :: planet_albedo
  REAL :: planet_emissivity
  REAL :: planet_aerosol_mmr
  REAL :: solar_zenith_angle
  REAL :: solar_azimuth_angle
  REAL :: planet_radius
  REAL :: planet_epoch
  REAL :: planet_e
  REAL :: planet_de
  REAL :: planet_lph
  REAL :: planet_dlph
  REAL :: planet_oblq
  REAL :: planet_doblq
  REAL :: planet_a
  REAL :: planet_da
  REAL :: planet_M
  REAL :: planet_dM
  REAL :: planet_ha
  REAL :: planet_dha
  REAL :: planet_obs_lon
  REAL :: planet_obs_lat
  REAL :: stellar_radius
  REAL :: sc
  REAL :: g
  REAL :: r
  REAL :: cp
  REAL :: pref
  REAL :: sclht
  REAL :: lapse
  REAL :: omega               
  LOGICAL :: l_planet_orbit
  LOGICAL :: l_set_planet_rotation
  LOGICAL :: l_planet_g
  LOGICAL :: l_planet_grey_surface
  LOGICAL :: l_planet_intrinsic_flux
  LOGICAL :: l_planet_aerosol
  LOGICAL :: l_fix_solang
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=planet_constants, IOSTAT=ErrorStatus,     &
       IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist planet_constants", iomessage)

  my_nml % i_planet = i_planet
  my_nml % i_eqt    = i_eqt
  ! end of integers
  my_nml % planet_t_intrinsic   = planet_t_intrinsic
  my_nml % planet_albedo        = planet_albedo
  my_nml % planet_emissivity    = planet_emissivity
  my_nml % planet_aerosol_mmr   = planet_aerosol_mmr
  my_nml % solar_zenith_angle   = solar_zenith_angle
  my_nml % solar_azimuth_angle  = solar_azimuth_angle
  my_nml % planet_radius        = planet_radius
  my_nml % planet_epoch         = planet_epoch
  my_nml % planet_e             = planet_e
  my_nml % planet_de            = planet_de
  my_nml % planet_lph           = planet_lph
  my_nml % planet_dlph          = planet_dlph
  my_nml % planet_oblq          = planet_oblq
  my_nml % planet_doblq         = planet_doblq
  my_nml % planet_a             = planet_a
  my_nml % planet_da            = planet_da
  my_nml % planet_M             = planet_M
  my_nml % planet_dM            = planet_dM
  my_nml % planet_ha            = planet_ha
  my_nml % planet_dha           = planet_dha
  my_nml % planet_obs_lon       = planet_obs_lon
  my_nml % planet_obs_lat       = planet_obs_lat
  my_nml % stellar_radius       = stellar_radius
  my_nml % sc                   = sc
  my_nml % g                    = g
  my_nml % r                    = r
  my_nml % cp                   = cp
  my_nml % pref                 = pref
  my_nml % sclht                = sclht
  my_nml % lapse                = lapse
  my_nml % omega                = omega
  ! end of reals
  my_nml % l_planet_orbit           = l_planet_orbit
  my_nml % l_set_planet_rotation    = l_set_planet_rotation
  my_nml % l_planet_g               = l_planet_g
  my_nml % l_planet_grey_surface    = l_planet_grey_surface
  my_nml % l_planet_intrinsic_flux  = l_planet_intrinsic_flux
  my_nml % l_planet_aerosol         = l_planet_aerosol
  my_nml % l_fix_solang             = l_fix_solang

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  i_planet = my_nml % i_planet
  i_eqt    = my_nml % i_eqt
  ! end of integers
  planet_t_intrinsic   = my_nml % planet_t_intrinsic
  planet_albedo        = my_nml % planet_albedo
  planet_emissivity    = my_nml % planet_emissivity
  planet_aerosol_mmr   = my_nml % planet_aerosol_mmr
  solar_zenith_angle   = my_nml % solar_zenith_angle
  solar_azimuth_angle  = my_nml % solar_azimuth_angle
  planet_radius        = my_nml % planet_radius
  planet_epoch         = my_nml % planet_epoch
  planet_e             = my_nml % planet_e
  planet_de            = my_nml % planet_de
  planet_lph           = my_nml % planet_lph
  planet_dlph          = my_nml % planet_dlph
  planet_oblq          = my_nml % planet_oblq
  planet_doblq         = my_nml % planet_doblq
  planet_a             = my_nml % planet_a
  planet_da            = my_nml % planet_da
  planet_M             = my_nml % planet_M
  planet_dM            = my_nml % planet_dM
  planet_ha            = my_nml % planet_ha
  planet_dha           = my_nml % planet_dha
  planet_obs_lon       = my_nml % planet_obs_lon
  planet_obs_lat       = my_nml % planet_obs_lat
  stellar_radius       = my_nml % stellar_radius
  sc                   = my_nml % sc
  g                    = my_nml % g
  r                    = my_nml % r
  cp                   = my_nml % cp
  pref                 = my_nml % pref
  sclht                = my_nml % sclht
  lapse                = my_nml % lapse
  omega                = my_nml % omega
  ! end of reals
  l_planet_orbit           = my_nml % l_planet_orbit
  l_set_planet_rotation    = my_nml % l_set_planet_rotation
  l_planet_g               = my_nml % l_planet_g
  l_planet_grey_surface    = my_nml % l_planet_grey_surface
  l_planet_intrinsic_flux  = my_nml % l_planet_intrinsic_flux
  l_planet_aerosol         = my_nml % l_planet_aerosol
  l_fix_solang             = my_nml % l_fix_solang

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_planet_constants
#endif

END MODULE planet_constants_mod

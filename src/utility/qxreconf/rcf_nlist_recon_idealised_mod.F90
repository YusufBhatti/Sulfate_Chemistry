! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_nlist_recon_idealised_mod

! Description:
!   Defines variables for the recon_idealised namelist,
!   and the namelist definition itself.
!   Also contains routines for reading, printing 
!   and checking the namelist.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

USE rcf_ideal_bubble_constants_mod, ONLY: &
    idl_max_num_bubbles

USE profiles_mod, ONLY: &
    num_data_max

USE missing_data_mod, ONLY: &
    imdi, rmdi

USE filenamelength_mod, ONLY: &
    filenamelength

IMPLICIT NONE

SAVE

!----------------------------------------------------
! Items in the recon_idealised namelist
!----------------------------------------------------

! Horizontal grid definition

REAL :: delta_xi1 = rmdi   ! xi1/lambda grid spacing (m)
REAL :: delta_xi2 = rmdi   ! xi2/phi grid spacing (m)
REAL :: base_xi1  = rmdi   ! xi1/lambda origin ( c.f. first longitude)
REAL :: base_xi2  = rmdi   ! xi2/phi origin (c.f. first latitude)

! Vertical grid definition

INTEGER :: grid_number     = imdi  ! Vertical grid spacing option
INTEGER :: grid_flat       = imdi  ! Vertical flattening above orography option
INTEGER :: big_layers      = imdi  ! Number of big layers in
                                   ! stretched+regular grids
INTEGER :: transit_layers  = imdi  ! Number of transition layers between
                                   ! regular and stretched
INTEGER :: first_constant_r_rho_level_new
                                   ! First constant height level above orography

REAL :: first_theta_height = rmdi  ! Height of lowest theta level (m)
REAL :: thin_theta_height  = rmdi  ! Height of lowest theta level
                                   ! above surface (m)
REAL :: height_domain      = rmdi  ! Height of lid for domain (m)
REAL :: big_factor         = rmdi  ! Big/normal layers scale factor
REAL :: vert_grid_ratio    = rmdi  ! Geometric stretching ratio

! Orography

INTEGER :: surface_type = imdi   ! Choice of orography

REAL :: h_o             = rmdi   ! Maximum orographic height (m)
REAL :: lambda_fraction = rmdi   ! E-W mountain position (normalised 0-1)
REAL :: phi_fraction    = rmdi   ! N-S mountain position (normalised 0-1)
REAL :: half_width_x    = rmdi   ! Mountain E-W half-width
REAL :: half_width_y    = rmdi   ! Mountain N-S half-width

! Initial profiles

INTEGER :: tprofile_number   = imdi   ! Temperature profile
INTEGER :: qprofile_number   = imdi   ! Moisture profile
INTEGER :: pressure_balance  = imdi   ! (i.e. pprofile_number)

REAL :: dtheta_dz1(3)        = rmdi   ! Lapse rate dtheta/dz
REAL :: height_dz1(2)        = rmdi   ! Layer boundaries for dtheta/dz

LOGICAL :: l_init_idealised  = .FALSE.   ! Master idealised initialisation flag.
        ! Causes multiple data fields to be overwritten with idealised profiles.
LOGICAL :: l_constant_dz     = .FALSE.   ! Constant lapse rate T/F
LOGICAL :: l_reset_mixing    = .FALSE.   ! Reset mixing ratios

! Solid rotating body

INTEGER :: AA_jet_m = imdi
INTEGER :: AA_jet_n = imdi

REAL :: AA_jet_A    = rmdi
REAL :: AA_jet_u0   = rmdi
REAL :: grid_NP_lon = rmdi   ! Longitude of north pole when l_rotate_grid=.TRUE.
REAL :: grid_NP_lat = rmdi   ! Latitude of north pole when l_rotate_grid=.TRUE.

! Baroclinic initial profiles

INTEGER :: b_const = imdi         ! A baroclinic wave parameter
INTEGER :: k_const = imdi         ! A baroclinic wave parameter
INTEGER :: initial_profile = imdi ! Choice of initial preset profile

REAL :: t0_p = rmdi          ! T0 pole
REAL :: t0_e = rmdi          ! T0 equator

LOGICAL :: l_rotate_grid    = .FALSE.  ! Use a rotated grid for the solid body
LOGICAL :: l_baro_perturbed = .FALSE.  ! Perturb winds in baroclinic profiles

! Perturbing initial fields

INTEGER :: perturb_type = imdi      ! Perturbation type

REAL :: perturb_magnitude_t = rmdi  ! Maximum magnitude of temp. perturbation
REAL :: perturb_magnitude_q = rmdi  ! Maximum magnitude of moisture perturbation
REAL :: perturb_height(2)   = rmdi  ! Height range for perturbations

LOGICAL :: l_perturb_t = .FALSE.    ! Randomly perturb temperature field
LOGICAL :: l_perturb_q = .FALSE.    ! Randomly perturb moisture field

! Bubble tests

INTEGER :: idl_bubble_option(idl_max_num_bubbles) = imdi  ! Bubble selection

REAL :: idl_bubble_max(idl_max_num_bubbles)     = rmdi
                                             ! Maximum bubble amplitude
REAL :: idl_bubble_xoffset(idl_max_num_bubbles) = rmdi
                                             ! Bubble x-offset (normalised 0-1)
REAL :: idl_bubble_yoffset(idl_max_num_bubbles) = rmdi
                                             ! Bubble y-offset (normalised 0-1)
REAL :: idl_bubble_width(idl_max_num_bubbles)   = rmdi
                                             ! Horizontal bubble width (m)
REAL :: idl_bubble_depth(idl_max_num_bubbles)   = rmdi
                                             ! Bubble depth (m)
REAL :: idl_bubble_height(idl_max_num_bubbles)  = rmdi
                                             ! Height of bubble's centre (m)
LOGICAL :: l_saturate_bubble = .FALSE.       ! Saturate all bubbles if .TRUE.

! Idealised forcing

INTEGER :: num_theta_init_heights = imdi  ! Size of initial pot. temp. profile
INTEGER :: theta_init_field_type  = imdi  ! Type of field for temperature init.
INTEGER :: num_mv_init_heights    = imdi  ! Size of initial vapour profile
INTEGER :: mv_init_field_type     = imdi  ! Type of field for vapour init.
INTEGER :: num_uv_init_heights    = imdi  ! Size of initial horiz. wind profile

REAL :: theta_init_height(num_data_max) = rmdi  ! Heights of temp. data points
REAL :: mv_init_height(num_data_max)    = rmdi  ! Heights of vapour data points
REAL :: uv_init_height(num_data_max)    = rmdi  ! Heights of wind data points
REAL :: theta_init_data(num_data_max)   = rmdi  ! Temp. profile data
REAL :: mv_init_data(num_data_max)      = rmdi  ! Vapour profile data
REAL :: u_init_data(num_data_max)       = rmdi  ! u wind profile data
REAL :: v_init_data(num_data_max)       = rmdi  ! v wind profile data

! Planetary modelling options

REAL :: initial_tolerance = rmdi   ! Initial Newton solver tolerance

CHARACTER(LEN=filenamelength) :: profile_filename = ''
                                   ! NetCDF file containing T-P profile


!----------------------------------------------------
! Namelist recon_idealised
!----------------------------------------------------

NAMELIST /recon_idealised/                                                   &
  ! Horizontal grid
  delta_xi1, delta_xi2, base_xi1, base_xi2,                                  &
  ! Vertical grid
  grid_number, grid_flat, big_layers, transit_layers,                        &
  first_constant_r_rho_level_new, first_theta_height, thin_theta_height,     &
  height_domain, big_factor, vert_grid_ratio,                                &
  ! Orography
  surface_type, h_o, lambda_fraction, phi_fraction, half_width_x,            &
  half_width_y,                                                              &
  ! Initial profiles
  tprofile_number, qprofile_number, pressure_balance,                        &
  dtheta_dz1, height_dz1, l_init_idealised, l_constant_dz, l_reset_mixing,   &
  ! Solid rotating body
  AA_jet_m, AA_jet_n, AA_jet_A, AA_jet_u0,                                   &
  grid_NP_lon, grid_NP_lat,                                                  &
  ! Baroclinic initial profiles
  b_const, k_const, initial_profile, t0_p, t0_e,                             &
  l_rotate_grid, l_baro_perturbed,                                           &
  ! Perturbing initial fields
  perturb_type, perturb_magnitude_t, perturb_magnitude_q,                    &
  perturb_height, l_perturb_t, l_perturb_q,                                  &
  ! Bubble tests
  idl_bubble_option, idl_bubble_max, l_saturate_bubble,                      &
  idl_bubble_xoffset, idl_bubble_yoffset,                                    &
  idl_bubble_width, idl_bubble_depth, idl_bubble_height,                     &
  ! Idealised forcing
  num_theta_init_heights, theta_init_field_type, num_mv_init_heights,        &
  mv_init_field_type, num_uv_init_heights,                                   &
  theta_init_height, mv_init_height, uv_init_height,                         &
  theta_init_data, mv_init_data, u_init_data, v_init_data,                   &
  ! Planetary modelling
  initial_tolerance, profile_filename

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                      & 
  ModuleName='RCF_NLIST_RECON_IDEALISED_MOD'

CONTAINS

SUBROUTINE read_nml_recon_idealised(unit_no)

! Read the recon_idealised namelist.

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

USE um_parcore, ONLY: &
    mype

USE check_iostat_mod, ONLY: &
    check_iostat
USE errormessagelength_mod, ONLY: &
    errormessagelength
USE setup_namelist, ONLY: &
    setup_nml_type

IMPLICIT NONE

! Arguments:
INTEGER, INTENT(IN) :: unit_no   ! Unit number

! Local variables
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode

CHARACTER(LEN=errormessagelength) :: iomessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='READ_NML_RECON_IDEALISED'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 4
INTEGER, PARAMETER :: n_int = 20 + idl_max_num_bubbles
INTEGER, PARAMETER :: n_real = 23 + 1*3 + 2*2 + 6*idl_max_num_bubbles        &
                             + 7*num_data_max
INTEGER, PARAMETER :: n_log = 8
INTEGER, PARAMETER :: n_chars = filenamelength

TYPE my_namelist
  SEQUENCE
  INTEGER :: grid_number
  INTEGER :: grid_flat
  INTEGER :: big_layers
  INTEGER :: transit_layers
  INTEGER :: first_constant_r_rho_level_new
  INTEGER :: surface_type
  INTEGER :: tprofile_number
  INTEGER :: qprofile_number
  INTEGER :: pressure_balance
  INTEGER :: AA_jet_m
  INTEGER :: AA_jet_n
  INTEGER :: b_const
  INTEGER :: k_const
  INTEGER :: initial_profile
  INTEGER :: perturb_type
  INTEGER :: idl_bubble_option(idl_max_num_bubbles)
  INTEGER :: num_theta_init_heights
  INTEGER :: theta_init_field_type
  INTEGER :: num_mv_init_heights
  INTEGER :: mv_init_field_type
  INTEGER :: num_uv_init_heights
  REAL :: delta_xi1
  REAL :: delta_xi2
  REAL :: base_xi1
  REAL :: base_xi2
  REAL :: first_theta_height
  REAL :: thin_theta_height
  REAL :: height_domain
  REAL :: big_factor
  REAL :: vert_grid_ratio
  REAL :: h_o
  REAL :: lambda_fraction
  REAL :: phi_fraction
  REAL :: half_width_x
  REAL :: half_width_y
  REAL :: dtheta_dz1(3)
  REAL :: height_dz1(2)
  REAL :: AA_jet_A
  REAL :: AA_jet_u0
  REAL :: grid_NP_lon
  REAL :: grid_NP_lat
  REAL :: t0_p
  REAL :: t0_e
  REAL :: perturb_magnitude_t
  REAL :: perturb_magnitude_q
  REAL :: perturb_height(2)
  REAL :: idl_bubble_max(idl_max_num_bubbles)
  REAL :: idl_bubble_xoffset(idl_max_num_bubbles)
  REAL :: idl_bubble_yoffset(idl_max_num_bubbles)
  REAL :: idl_bubble_width(idl_max_num_bubbles)
  REAL :: idl_bubble_depth(idl_max_num_bubbles)
  REAL :: idl_bubble_height(idl_max_num_bubbles)
  REAL :: theta_init_height(num_data_max)
  REAL :: mv_init_height(num_data_max)
  REAL :: uv_init_height(num_data_max)
  REAL :: theta_init_data(num_data_max)
  REAL :: mv_init_data(num_data_max)
  REAL :: u_init_data(num_data_max)
  REAL :: v_init_data(num_data_max)
  REAL :: initial_tolerance
  LOGICAL :: l_init_idealised
  LOGICAL :: l_constant_dz
  LOGICAL :: l_reset_mixing
  LOGICAL :: l_rotate_grid
  LOGICAL :: l_baro_perturbed
  LOGICAL :: l_perturb_t
  LOGICAL :: l_perturb_q
  LOGICAL :: l_saturate_bubble
  CHARACTER(LEN=filenamelength) :: profile_filename
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,               &
                    n_real_in=n_real, n_log_in=n_log, n_chars_in=n_chars)

IF (mype == 0) THEN

  READ(UNIT=unit_no, NML=recon_idealised, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(ErrorStatus, "namelist recon_idealised", iomessage)

  ! Integers
  my_nml % grid_number                       = grid_number
  my_nml % grid_flat                         = grid_flat
  my_nml % big_layers                        = big_layers
  my_nml % transit_layers                    = transit_layers
  my_nml % first_constant_r_rho_level_new    = first_constant_r_rho_level_new
  my_nml % surface_type                      = surface_type
  my_nml % tprofile_number                   = tprofile_number
  my_nml % qprofile_number                   = qprofile_number
  my_nml % pressure_balance                  = pressure_balance
  my_nml % AA_jet_m                          = AA_jet_m
  my_nml % AA_jet_n                          = AA_jet_n
  my_nml % b_const                           = b_const
  my_nml % k_const                           = k_const
  my_nml % initial_profile                   = initial_profile
  my_nml % perturb_type                      = perturb_type
  my_nml % idl_bubble_option                 = idl_bubble_option
  my_nml % num_theta_init_heights            = num_theta_init_heights
  my_nml % theta_init_field_type             = theta_init_field_type
  my_nml % num_mv_init_heights               = num_mv_init_heights
  my_nml % mv_init_field_type                = mv_init_field_type
  my_nml % num_uv_init_heights               = num_uv_init_heights
  ! Reals
  my_nml % delta_xi1                         = delta_xi1
  my_nml % delta_xi2                         = delta_xi2
  my_nml % base_xi1                          = base_xi1
  my_nml % base_xi2                          = base_xi2
  my_nml % first_theta_height                = first_theta_height
  my_nml % thin_theta_height                 = thin_theta_height
  my_nml % height_domain                     = height_domain
  my_nml % big_factor                        = big_factor
  my_nml % vert_grid_ratio                   = vert_grid_ratio
  my_nml % h_o                               = h_o
  my_nml % lambda_fraction                   = lambda_fraction
  my_nml % phi_fraction                      = phi_fraction
  my_nml % half_width_x                      = half_width_x
  my_nml % half_width_y                      = half_width_y
  my_nml % dtheta_dz1                        = dtheta_dz1
  my_nml % height_dz1                        = height_dz1
  my_nml % AA_jet_A                          = AA_jet_A
  my_nml % AA_jet_u0                         = AA_jet_u0
  my_nml % grid_NP_lon                       = grid_NP_lon
  my_nml % grid_NP_lat                       = grid_NP_lat
  my_nml % t0_p                              = t0_p
  my_nml % t0_e                              = t0_e
  my_nml % perturb_magnitude_t               = perturb_magnitude_t
  my_nml % perturb_magnitude_q               = perturb_magnitude_q
  my_nml % perturb_height                    = perturb_height
  my_nml % idl_bubble_max                    = idl_bubble_max
  my_nml % idl_bubble_xoffset                = idl_bubble_xoffset
  my_nml % idl_bubble_yoffset                = idl_bubble_yoffset
  my_nml % idl_bubble_width                  = idl_bubble_width
  my_nml % idl_bubble_depth                  = idl_bubble_depth
  my_nml % idl_bubble_height                 = idl_bubble_height
  my_nml % theta_init_height                 = theta_init_height
  my_nml % mv_init_height                    = mv_init_height
  my_nml % uv_init_height                    = uv_init_height
  my_nml % theta_init_data                   = theta_init_data
  my_nml % mv_init_data                      = mv_init_data
  my_nml % u_init_data                       = u_init_data
  my_nml % v_init_data                       = v_init_data
  my_nml % initial_tolerance                 = initial_tolerance
  ! Logicals
  my_nml % l_init_idealised                  = l_init_idealised
  my_nml % l_constant_dz                     = l_constant_dz
  my_nml % l_reset_mixing                    = l_reset_mixing
  my_nml % l_rotate_grid                     = l_rotate_grid
  my_nml % l_baro_perturbed                  = l_baro_perturbed
  my_nml % l_perturb_t                       = l_perturb_t
  my_nml % l_perturb_q                       = l_perturb_q
  my_nml % l_saturate_bubble                 = l_saturate_bubble
  ! Characters
  my_nml % profile_filename                  = profile_filename

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  ! Integers
  grid_number                       = my_nml % grid_number
  grid_flat                         = my_nml % grid_flat
  big_layers                        = my_nml % big_layers
  transit_layers                    = my_nml % transit_layers
  first_constant_r_rho_level_new    = my_nml % first_constant_r_rho_level_new
  surface_type                      = my_nml % surface_type
  tprofile_number                   = my_nml % tprofile_number
  qprofile_number                   = my_nml % qprofile_number
  pressure_balance                  = my_nml % pressure_balance
  AA_jet_m                          = my_nml % AA_jet_m
  AA_jet_n                          = my_nml % AA_jet_n
  b_const                           = my_nml % b_const
  k_const                           = my_nml % k_const
  initial_profile                   = my_nml % initial_profile
  perturb_type                      = my_nml % perturb_type
  idl_bubble_option                 = my_nml % idl_bubble_option
  num_theta_init_heights            = my_nml % num_theta_init_heights
  theta_init_field_type             = my_nml % theta_init_field_type
  num_mv_init_heights               = my_nml % num_mv_init_heights
  mv_init_field_type                = my_nml % mv_init_field_type
  num_uv_init_heights               = my_nml % num_uv_init_heights
  ! Reals
  delta_xi1                         = my_nml % delta_xi1
  delta_xi2                         = my_nml % delta_xi2
  base_xi1                          = my_nml % base_xi1
  base_xi2                          = my_nml % base_xi2
  first_theta_height                = my_nml % first_theta_height
  thin_theta_height                 = my_nml % thin_theta_height
  height_domain                     = my_nml % height_domain
  big_factor                        = my_nml % big_factor
  vert_grid_ratio                   = my_nml % vert_grid_ratio
  h_o                               = my_nml % h_o
  lambda_fraction                   = my_nml % lambda_fraction
  phi_fraction                      = my_nml % phi_fraction
  half_width_x                      = my_nml % half_width_x
  half_width_y                      = my_nml % half_width_y
  dtheta_dz1                        = my_nml % dtheta_dz1
  height_dz1                        = my_nml % height_dz1
  AA_jet_A                          = my_nml % AA_jet_A
  AA_jet_u0                         = my_nml % AA_jet_u0
  grid_NP_lon                       = my_nml % grid_NP_lon
  grid_NP_lat                       = my_nml % grid_NP_lat
  t0_p                              = my_nml % t0_p
  t0_e                              = my_nml % t0_e
  perturb_magnitude_t               = my_nml % perturb_magnitude_t
  perturb_magnitude_q               = my_nml % perturb_magnitude_q
  perturb_height                    = my_nml % perturb_height
  idl_bubble_max                    = my_nml % idl_bubble_max
  idl_bubble_xoffset                = my_nml % idl_bubble_xoffset
  idl_bubble_yoffset                = my_nml % idl_bubble_yoffset
  idl_bubble_width                  = my_nml % idl_bubble_width
  idl_bubble_depth                  = my_nml % idl_bubble_depth
  idl_bubble_height                 = my_nml % idl_bubble_height
  theta_init_height                 = my_nml % theta_init_height
  mv_init_height                    = my_nml % mv_init_height
  uv_init_height                    = my_nml % uv_init_height
  theta_init_data                   = my_nml % theta_init_data
  mv_init_data                      = my_nml % mv_init_data
  u_init_data                       = my_nml % u_init_data
  v_init_data                       = my_nml % v_init_data
  initial_tolerance                 = my_nml % initial_tolerance
  ! Logicals
  l_init_idealised                  = my_nml % l_init_idealised
  l_constant_dz                     = my_nml % l_constant_dz
  l_reset_mixing                    = my_nml % l_reset_mixing
  l_rotate_grid                     = my_nml % l_rotate_grid
  l_baro_perturbed                  = my_nml % l_baro_perturbed
  l_perturb_t                       = my_nml % l_perturb_t
  l_perturb_q                       = my_nml % l_perturb_q
  l_saturate_bubble                 = my_nml % l_saturate_bubble
  ! Characters
  profile_filename                  = my_nml % profile_filename

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_recon_idealised

SUBROUTINE print_nlist_recon_idealised()

! Print the contents of the recon_idealised namelist.
! Some very large arrays are omitted.

USE umPrintMgr, ONLY: umPrint

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer
CHARACTER(LEN=*), PARAMETER  :: RoutineName='PRINT_NLIST_RECON_IDEALISED'
CHARACTER(LEN=15) :: bubble_format = ''

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL umPrint('Contents of namelist recon_idealised', src=RoutineName)

! Horizontal grid
WRITE(lineBuffer,'(A,E12.5)') ' delta_xi1 = ', delta_xi1
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' delta_xi2 = ', delta_xi2
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' base_xi1 = ', base_xi1
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' base_xi2 = ', base_xi2
CALL umPrint(lineBuffer,src=RoutineName)

! Vertical grid
WRITE(lineBuffer,'(A,I0)') ' grid_number = ', grid_number
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' grid_flat = ', grid_flat
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' big_layers = ', big_layers
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' transit_layers = ', transit_layers
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' first_constant_r_rho_level_new = ',             &
                             first_constant_r_rho_level_new
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' first_theta_height = ', first_theta_height
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' thin_theta_height = ', thin_theta_height
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' height_domain = ', height_domain
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' big_factor = ', big_factor
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' vert_grid_ratio = ', vert_grid_ratio
CALL umPrint(lineBuffer,src=RoutineName)

! Orography
WRITE(lineBuffer,'(A,I0)') ' surface_type = ', surface_type
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' h_o = ', h_o
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' lambda_fraction = ', lambda_fraction
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' phi_fraction = ', phi_fraction
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' half_width_x = ', half_width_x
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' half_width_y = ', half_width_y
CALL umPrint(lineBuffer,src=RoutineName)

! Initial profiles
WRITE(lineBuffer,'(A,I0)') ' tprofile_number = ', tprofile_number
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' qprofile_number = ', qprofile_number
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' pressure_balance = ', pressure_balance
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,3(E12.5, 2X))') ' dtheta_dz1 = ', dtheta_dz1(:)
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,2(E12.5, 2X))') ' height_dz1 = ', height_dz1(:)
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,L1)') ' l_init_idealised = ', l_init_idealised
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,L1)') ' l_constant_dz = ', l_constant_dz
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,L1)') ' l_reset_mixing = ', l_reset_mixing
CALL umPrint(lineBuffer,src=RoutineName)

! Solid rotating body
WRITE(lineBuffer,'(A,I0)') ' AA_jet_m = ', AA_jet_m
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' AA_jet_n = ', AA_jet_n
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' AA_jet_A = ', AA_jet_A
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' AA_jet_u0 = ', AA_jet_u0
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' grid_NP_lon = ', grid_NP_lon
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' grid_NP_lat = ', grid_NP_lat
CALL umPrint(lineBuffer,src=RoutineName)

! Baroclinic initial profiles
WRITE(lineBuffer,'(A,I0)') ' b_const = ', b_const
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' k_const = ', k_const
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' initial_profile = ', initial_profile
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' t0_p = ', t0_p
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' t0_e = ', t0_e
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,L1)') ' l_rotate_grid = ', l_rotate_grid
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,L1)') ' l_baro_perturbed = ', l_baro_perturbed
CALL umPrint(lineBuffer,src=RoutineName)

! Perturbing initial fields
WRITE(lineBuffer,'(A,I0)') 'perturb_type = ', perturb_type
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' perturb_magnitude_t = ', perturb_magnitude_t
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,E12.5)') ' perturb_magnitude_q = ', perturb_magnitude_q
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,2(E12.5, 2X))') ' perturb_height = ', perturb_height(:)
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,L1)') ' l_perturb_t = ', l_perturb_t
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,L1)') ' l_perturb_q = ', l_perturb_q
CALL umPrint(lineBuffer,src=RoutineName)

! Bubble tests
! '(A,n(L1,2X))':
WRITE(bubble_format, '(A,I0,A)') '(A,',idl_max_num_bubbles,'(I0,2X))'
WRITE(lineBuffer, bubble_format) ' idl_bubble_option = ', idl_bubble_option
CALL umPrint(lineBuffer,src=RoutineName)
! '(A,n(E12.5,2X))':
WRITE(bubble_format, '(A,I0,A)') '(A,',idl_max_num_bubbles,'(E12.5,2X))'
WRITE(lineBuffer, bubble_format) ' idl_bubble_max = ', idl_bubble_max
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer, bubble_format) ' idl_bubble_xoffset = ', idl_bubble_xoffset
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer, bubble_format) ' idl_bubble_yoffset = ', idl_bubble_yoffset
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer, bubble_format) ' idl_bubble_width = ', idl_bubble_width
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer, bubble_format) ' idl_bubble_depth = ', idl_bubble_depth
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer, bubble_format) ' idl_bubble_height = ', idl_bubble_height
CALL umPrint(lineBuffer,src=RoutineName)
! '(A,n(L1,2X))':
WRITE(bubble_format, '(A,I0,A)') '(A,',idl_max_num_bubbles,'(L1,2X))'
WRITE(lineBuffer,'(A,L1)') ' l_saturate_bubble = ', l_saturate_bubble
CALL umPrint(lineBuffer,src=RoutineName)

! Idealised forcing
WRITE(lineBuffer,'(A,I0)') ' num_theta_init_heights = ', num_theta_init_heights
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' theta_init_field_type = ', theta_init_field_type
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' num_mv_init_heights = ', num_mv_init_heights
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' mv_init_field_type = ', mv_init_field_type
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,I0)') ' num_uv_init_heights = ', num_uv_init_heights
CALL umPrint(lineBuffer,src=RoutineName)
! Heights and data points are large arrays and are not printed.

! Planetary modelling
WRITE(lineBuffer,'(A,E12.5)') ' initial_tolerance = ', initial_tolerance
CALL umPrint(lineBuffer,src=RoutineName)
WRITE(lineBuffer,'(A,A)') ' profile_filename = ', TRIM(profile_filename)
CALL umPrint(lineBuffer,src=RoutineName)

CALL umPrint('- - - - - - end of namelist - - - - - -',    &
    src=RoutineName)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE print_nlist_recon_idealised


SUBROUTINE check_nlist_recon_idealised

! Check the recon_idealised namelist.

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

USE ereport_mod, ONLY: &
    ereport

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE chk_opts_mod, ONLY: &
    chk_var,            &
    def_src

USE dynamics_testing_mod, ONLY: &
    problem_number

USE model_domain_mod, ONLY: &
    model_type,             &
    mt_bi_cyclic_lam,       &
    mt_global,              &
    l_cartesian

IMPLICIT NONE

INTEGER :: i
INTEGER :: ErrorStatus
CHARACTER(LEN=5) :: element = ''

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='CHECK_NLIST_RECON_IDEALISED'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

def_src = ModuleName//':'//RoutineName

IF (grid_number /= 11) THEN
  CALL chk_var(first_constant_r_rho_level_new, &
            'first_constant_r_rho_level_new', '[>=1]')
END IF

CALL chk_var(grid_number, 'grid_number', '[1,3:6,8,10,11,21:23]')

IF (grid_number == 5) THEN
  CALL chk_var(big_factor,     'big_factor', '[0.0:1.0]')
  CALL chk_var(big_layers,     'big_layers', '[>=0]')
  CALL chk_var(transit_layers, 'transit_layers', '[>=0]')
ELSE IF (grid_number == 6) THEN
  CALL chk_var(first_theta_height, 'first_theta_height', '[>=0.0]')
  CALL chk_var(thin_theta_height,  'thin_theta_height',  '[>=0.0]')
ELSE IF (grid_number == 8) THEN
  CALL chk_var(vert_grid_ratio, 'vert_grid_ratio', '[>0.0]')
END IF

IF (grid_number == 1 .OR. grid_number == 3  .OR. grid_number == 8 &
                     .OR. grid_number == 10 .OR. grid_number == 21) THEN
  CALL chk_var(height_domain, 'height_domain', '[>0.0]')
END IF

IF (grid_number <= 23 .AND. grid_number /= 11) THEN
  CALL chk_var(grid_flat, 'grid_flat', '[0,3]')
END IF

IF (model_type == mt_global .AND. l_cartesian) THEN
  errorstatus = 100
  CALL ereport(RoutineName, errorstatus, 'Global model domains '//  &
    'using a Cartesian co-ordinate system are not supported.')
END IF

IF (model_type == mt_bi_cyclic_lam .AND. .NOT. l_cartesian) THEN
  errorstatus = 110
  CALL ereport(RoutineName, errorstatus, 'Bicyclic LAM model domains '//  &
    'require a Cartesian co-ordinate system to be selected.')
END IF


IF (l_cartesian) THEN
  CALL chk_var(delta_xi1, 'delta_xi1', '[>0.0]')
  CALL chk_var(delta_xi2, 'delta_xi2', '[>0.0]')
END IF

IF (l_init_idealised) THEN

  CALL chk_var(surface_type, 'surface_type', '[0,6,8:13]')
  IF (surface_type == 0) THEN

    IF (num_mv_init_heights > 0) THEN
      CALL chk_var(num_mv_init_heights, 'num_mv_init_heights', '[1:100]')
      IF (ANY(mv_init_data(1:num_mv_init_heights) < 0.0)) THEN
        errorstatus = 200
        CALL ereport(RoutineName, errorstatus, 'Negative value in mv_init_data')
      END IF
    END IF
    IF (num_theta_init_heights > 0) THEN
      CALL chk_var(num_theta_init_heights, 'num_theta_init_heights', '[1:100]')
    END IF
    IF (num_uv_init_heights > 0) THEN
      CALL chk_var(num_uv_init_heights, 'num_uv_init_heights', '[1:100]')
    END IF

  ELSE IF (surface_type == 6 .OR. surface_type == 8) THEN
    CALL chk_var(h_o,             'h_o',             '[>=0.0]')
    CALL chk_var(half_width_x,    'half_width_x',    '[>0.0]')
    CALL chk_var(half_width_y,    'half_width_y',    '[>0.0]')
    CALL chk_var(lambda_fraction, 'lambda_fraction', '[0.0:1.0]')
    IF (surface_type == 6) THEN
      CALL chk_var(phi_fraction, 'phi_fraction', '[0.0:1.0]')
    END IF

  ELSE IF (surface_type == 11 .OR. surface_type == 12) THEN
    CALL chk_var(h_o, 'h_o', '[>=0.0]')
  END IF


  ! Tests for specific problem_number cases:
  IF (problem_number == 2) THEN

    CALL chk_var(initial_profile, 'initial_profile','[0:3]')

    IF (l_rotate_grid) THEN
      CALL chk_var(grid_np_lon, 'grid_np_lon', '[0.0:360.0]')
      CALL chk_var(grid_np_lat, 'grid_np_lat', '[-90.0:90.0]')
    END IF

  ELSE IF (problem_number == 3) THEN

    CALL chk_var(pressure_balance, 'pressure_balance', '[0:2]')

    IF (l_perturb_q) THEN
      CALL chk_var(perturb_magnitude_q, 'perturb_magnitude_q', '[>=0.0]')
    END IF
    IF (l_perturb_t) THEN
      CALL chk_var(perturb_magnitude_t, 'perturb_magnitude_t', '[>=0.0]')
    END IF

  ELSE IF (problem_number == 5) THEN

    CALL chk_var(qprofile_number, 'qprofile_number','[0,10,11]')

  END IF

  IF (problem_number == 2 .OR. problem_number == 3) THEN
    DO i = 1, idl_max_num_bubbles
      WRITE(element, '(A,I0,A)') '(',i,')'
      CALL chk_var(idl_bubble_option(i), 'idl_bubble_option'//element, '[0:4]')

      IF (idl_bubble_option(i) > 0) THEN
        CALL chk_var(idl_bubble_depth(i),   &
                    'idl_bubble_depth'//element,   '[>0.0]')
        CALL chk_var(idl_bubble_height(i),  &
                    'idl_bubble_height'//element,  '[>0.0]')
        CALL chk_var(idl_bubble_max(i),     &
                    'idl_bubble_max'//element,     '[>0.0]')
        CALL chk_var(idl_bubble_width(i),   &
                    'idl_bubble_width'//element,   '[>0.0]')
        CALL chk_var(idl_bubble_xoffset(i), &
                    'idl_bubble_xoffset'//element, '[0.0:1.0]')
        CALL chk_var(idl_bubble_yoffset(i), &
                    'idl_bubble_yoffset'//element, '[0.0:1.0]')
      END IF
    END DO
  END IF

  IF (problem_number == 2 .OR. problem_number == 5) THEN

    CALL chk_var(tprofile_number, 'tprofile_number', '[1:4,10:17]')
    IF (tprofile_number == 1) THEN
      IF (ANY(height_dz1 <= 0.0)) THEN
        errorstatus = 300
        CALL ereport(RoutineName, errorstatus, &
          'Negative value of height_dz1 provided')
      END IF

    ELSE IF (tprofile_number == 15) THEN
      IF (LEN_TRIM(profile_filename) == 0) THEN
        errorstatus = 400
        CALL ereport(RoutineName, errorstatus, &
          'tprofile_number = 15 selected but profile_filename is empty')
      END IF
  
    END IF
  END IF

END IF  ! l_init_idealised


def_src = ''

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nlist_recon_idealised
END MODULE rcf_nlist_recon_idealised_mod

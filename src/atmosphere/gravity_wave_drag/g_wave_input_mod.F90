! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE g_wave_input_mod

! Description:
!       Input/namelist control of GWD scheme.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GWD

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.

USE missing_data_mod, ONLY: rmdi,imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! Control logicals
LOGICAL :: l_gwd      = .FALSE.  ! Use orographic drag scheme
LOGICAL :: l_use_ussp = .FALSE.  ! Use spectral GWD scheme, with opaque lid

LOGICAL :: l_gwd_40km = .TRUE.   ! Turn off orographic GWD above 40km

LOGICAL :: l_nonhydro = .FALSE.  ! Turn on nonhydro scheme
LOGICAL :: l_dynbeta  = .FALSE.  ! Turn on dynamically adjusting angle
                                 ! of group vel in nonhydro scheme

LOGICAL :: l_taus_scale = .FALSE.! Froude No. Dependent surface stress

LOGICAL :: l_fix_gwsatn = .TRUE. ! Bug fixes to gwsatn
INTEGER :: sat_scheme = 0        ! Switch to determine whether to use
                                 ! stress based    (sat_scheme=0) or
                                 ! amplitude based (sat_scheme=1)
                                 ! saturation test

! GWD variable inputs

INTEGER :: i_gwd_vn = imdi    ! Switch to determine version of GWD scheme
                              ! 4 => 4A scheme
                              ! 5 => 5A scheme
INTEGER, PARAMETER :: i_gwd_vn_4a = 4
INTEGER, PARAMETER :: i_gwd_vn_5a = 5

REAL :: kay_gwave   =  rmdi   ! Surface stress constant for GW drag
REAL :: gwd_frc     =  rmdi   ! Critical Froude Number for 4A Scheme
REAL :: nsigma      =  rmdi   ! Number of standard deviations above the
                              ! mean orography of top of sub-grid mountains

! Maximum and minimum values for the STPH_RP scheme
! Gravity Wave drag parameters moved to stochastic_physics_run_mod.F90
! Max and min values for the critical froud number
! Max and min values for the gravity wave constant

REAL :: gwd_fsat  = rmdi      ! Critical Froude number for wave breaking
                              ! used in the amplitude based saturation test


!New parameters for 5A version of orographic drag scheme
REAL    :: gsharp  = rmdi     ! Function of mountain sharpness (used to tune
                              ! amplitude of orographic gwd)
REAL    :: fbcd = rmdi        ! Flow blocking drag coefficient

LOGICAL :: l_smooth = .FALSE. ! Turn on smoothing of acceleration over
                              ! a vertical wavelength

LOGICAL :: l_gw_heating(3) = .FALSE. 
                              ! Turn on for heating due to
                              ! (1) flow blocking drag
                              ! (2) mountain-wave dissipation
                              ! (3) non-orographic gravity-wave dissipation

! Non-orographic gravity wave (USSP scheme) variable inputs
!
! Switch to determine version of USSP scheme (presently always '1a')
INTEGER :: i_ussp_vn = imdi
!
! 1 => 1A scheme (Standard global invariant source)
INTEGER, PARAMETER :: i_ussp_vn_1a = 1
!
! Factor enhancement for invariant global wave launch amplitude
REAL    :: ussp_launch_factor = rmdi
!
! Factor enhancement for conversion from convective rain to GW source flux
REAL    :: cgw_scale_factor   = rmdi
!
! Characteristic (spectrum peak) wavelength (m)
REAL    :: wavelstar          = rmdi
!
! Switch for precipitation dependent variable source flux calculation
! F : calculate standard USSP isotropic GW launch flux
! T : calculate variable CGW launch flux
LOGICAL :: l_add_cgw = .FALSE.

NAMELIST/run_gwd/i_gwd_vn                                         &
                ,kay_gwave, gwd_frc, nsigma, gwd_fsat             &
                , gsharp, fbcd, l_smooth, l_gw_heating            &
                ,ussp_launch_factor, wavelstar, l_gwd, l_use_ussp,&
                 i_ussp_vn, l_add_cgw, cgw_scale_factor 

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='G_WAVE_INPUT_MOD'

CONTAINS
SUBROUTINE print_nlist_run_gwd()
USE umPrintMgr, ONLY: umPrint,ummessage,PrNorm
IMPLICIT NONE
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_GWD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_gwd', &
    src='g_wave_input_mod')

WRITE(ummessage,'(A,I0)')' i_gwd_vn = ',i_gwd_vn
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,F14.2)')' kay_gwave = ',kay_gwave
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,F14.2)')' gwd_frc = ',gwd_frc
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,F14.2)')' nsigma = ',nsigma
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,F14.2)')' gwd_fsat = ',gwd_fsat
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,F14.2)')' gsharp = ',gsharp
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,F14.2)')' fbcd = ',fbcd
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,L1)')' l_smooth = ',l_smooth
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,3L1)')' l_gw_heating = ',l_gw_heating
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,F16.4)')' ussp_launch_factor = ',ussp_launch_factor
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,F16.4)')' wavelstar = ',wavelstar
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,L1)')' l_gwd = ',l_gwd
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,L1)')' l_use_ussp = ',l_use_ussp
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,I0)')' i_ussp_vn = ',i_ussp_vn
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,L1)')' l_add_cgw = ',l_add_cgw
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
WRITE(ummessage,'(A,F16.4)')' cgw_scale_factor = ',cgw_scale_factor
CALL umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='g_wave_input_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_gwd

SUBROUTINE read_nml_run_gwd(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_GWD'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 2
INTEGER, PARAMETER :: n_real = 9
INTEGER, PARAMETER :: n_log = 4 + 3

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_gwd_vn
  INTEGER :: i_ussp_vn
  REAL :: kay_gwave
  REAL :: gwd_frc
  REAL :: nsigma
  REAL :: gwd_fsat
  REAL :: gsharp
  REAL :: fbcd
  REAL :: ussp_launch_factor
  REAL :: wavelstar
  REAL :: cgw_scale_factor
  LOGICAL :: l_smooth
  LOGICAL :: l_gw_heating(3)
  LOGICAL :: l_gwd
  LOGICAL :: l_use_ussp
  LOGICAL :: l_add_cgw
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=run_gwd, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_GWD", iomessage)

  my_nml % i_gwd_vn    = i_gwd_vn
  my_nml % i_ussp_vn   = i_ussp_vn
  ! end of integers
  my_nml % kay_gwave   = kay_gwave
  my_nml % gwd_frc     = gwd_frc
  my_nml % nsigma      = nsigma
  my_nml % gwd_fsat    = gwd_fsat
  my_nml % gsharp      = gsharp
  my_nml % fbcd        = fbcd
  my_nml % ussp_launch_factor = ussp_launch_factor
  my_nml % wavelstar   = wavelstar
  my_nml % cgw_scale_factor   = cgw_scale_factor
  ! end of reals
  my_nml % l_smooth      = l_smooth
  my_nml % l_gw_heating  = l_gw_heating
  my_nml % l_gwd         = l_gwd
  my_nml % l_use_ussp    = l_use_ussp
  my_nml % l_add_cgw     = l_add_cgw

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  i_gwd_vn    = my_nml % i_gwd_vn
  i_ussp_vn   = my_nml % i_ussp_vn
  ! end of integers
  kay_gwave   = my_nml % kay_gwave
  gwd_frc     = my_nml % gwd_frc
  nsigma      = my_nml % nsigma
  gwd_fsat    = my_nml % gwd_fsat
  gsharp      = my_nml % gsharp
  fbcd        = my_nml % fbcd
  ussp_launch_factor = my_nml % ussp_launch_factor
  wavelstar   = my_nml % wavelstar
  cgw_scale_factor   = my_nml % cgw_scale_factor
  ! end of reals
  l_smooth      = my_nml % l_smooth
  l_gw_heating  = my_nml % l_gw_heating
  l_gwd         = my_nml % l_gwd
  l_use_ussp    = my_nml % l_use_ussp
  l_add_cgw     = my_nml % l_add_cgw

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_gwd

END MODULE g_wave_input_mod

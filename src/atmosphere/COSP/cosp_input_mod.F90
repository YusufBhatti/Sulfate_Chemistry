! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE cosp_input_mod

! Description:
!       Input/namelist control of COSP simulator.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: COSP

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.

USE missing_data_mod, ONLY: rmdi, imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!- Control variables set in run_cosp namelist
LOGICAL :: l_cosp            = .FALSE.
LOGICAL :: cosp_cloudsat_sim = .FALSE.
LOGICAL :: cosp_lidar_sim    = .FALSE.
LOGICAL :: cosp_isccp_sim    = .FALSE.
LOGICAL :: cosp_misr_sim     = .FALSE.
LOGICAL :: cosp_modis_sim    = .FALSE.
LOGICAL :: cosp_rttov_sim    = .FALSE.
LOGICAL :: cosp_use_vgrid    = .FALSE.
REAL    :: cosp_sr_cloud     = rmdi

!- Defaults for vertical grid if cosp_use_vgrid is set to true
INTEGER :: cosp_nlr        = 40
LOGICAL :: cosp_csat_vgrid = .TRUE.

!- Other control variables
INTEGER,PARAMETER :: cosp_ncolumns_max = 64


!- Inputs related to radar simulations
REAL    :: cosp_radar_freq     = 94.0
INTEGER :: cosp_surface_radar  = 0
INTEGER :: cosp_use_mie_tables = 0
INTEGER :: cosp_use_gas_abs    = 1
INTEGER :: cosp_do_ray         = 0
INTEGER :: cosp_melt_lay       = 0
REAL    :: cosp_k2             = -1
LOGICAL :: cosp_use_reff       = .TRUE.
LOGICAL :: cosp_use_precipitation_fluxes = .FALSE.
!- Inputs related to lidar simulations
INTEGER :: cosp_nprmts_max_hydro = 12
INTEGER :: cosp_naero            = 1
INTEGER :: cosp_nprmts_max_aero  = 1
INTEGER :: cosp_lidar_ice_type   = 0
!- Inputs related to ISCCP simulator
INTEGER :: cosp_overlap                   = 3
INTEGER :: cosp_isccp_topheight           = 1
INTEGER :: cosp_isccp_topheight_direction = 2
REAL    :: cosp_emsfc_lw                  = 0.99
!- Inputs related to RTTOV
INTEGER :: cosp_platform   = 1
INTEGER :: cosp_satellite  = 15
INTEGER :: cosp_instrument = 0
INTEGER,PARAMETER :: COSP_NCHANNELS = 8
INTEGER :: cosp_channels(cosp_nchannels) = (/1,3,5,6,8,10,11,13/)
REAL    :: cosp_surfem(cosp_nchannels) = (/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
REAL    :: cosp_zenang = 50.0
REAL    :: cosp_co2    = 5.241E-04
REAL    :: cosp_ch4    = 9.139E-07
REAL    :: cosp_n2o    = 4.665E-07
REAL    :: cosp_co     = 2.098E-07


NAMELIST/RUN_COSP/ cosp_cloudsat_sim, cosp_lidar_sim, cosp_isccp_sim, &
                   cosp_misr_sim, cosp_modis_sim, cosp_rttov_sim,     &
                   cosp_use_vgrid, l_cosp, cosp_sr_cloud

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='COSP_INPUT_MOD'

CONTAINS

  SUBROUTINE print_nlist_run_cosp()
    USE umPrintMgr, ONLY : umPrint
    IMPLICIT NONE
    CHARACTER(LEN=50000) :: lineBuffer
    REAL(KIND=jprb) :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_COSP'

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    CALL umPrint('Contents of namelist run_cosp', &
        src='cosp_input_mod')

    WRITE(lineBuffer,*)' cosp_cloudsat_sim = ',cosp_cloudsat_sim
    CALL umPrint(lineBuffer,src='cosp_input_mod')
    WRITE(lineBuffer,*)' cosp_lidar_sim = ',cosp_lidar_sim
    CALL umPrint(lineBuffer,src='cosp_input_mod')
    WRITE(lineBuffer,*)' cosp_isccp_sim = ',cosp_isccp_sim
    CALL umPrint(lineBuffer,src='cosp_input_mod')
    WRITE(lineBuffer,*)' cosp_misr_sim = ',cosp_misr_sim
    CALL umPrint(lineBuffer,src='cosp_input_mod')
    WRITE(lineBuffer,*)' cosp_modis_sim = ',cosp_modis_sim
    CALL umPrint(lineBuffer,src='cosp_input_mod')
    WRITE(lineBuffer,*)' cosp_rttov_sim = ',cosp_rttov_sim
    CALL umPrint(lineBuffer,src='cosp_input_mod')
    WRITE(lineBuffer,*)' cosp_use_vgrid = ',cosp_use_vgrid
    CALL umPrint(lineBuffer,src='cosp_input_mod')
    WRITE(lineBuffer,'(A,E10.3)')' cosp_sr_cloud = ',cosp_sr_cloud
    CALL umPrint(lineBuffer,src='cosp_input_mod')

    CALL umPrint('- - - - - - end of namelist - - - - - -', &
        src='cosp_input_mod')

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE print_nlist_run_cosp

  SUBROUTINE read_nml_run_cosp(unit_in)

    USE um_parcore, ONLY: mype

    USE check_iostat_mod, ONLY : check_iostat

    USE setup_namelist, ONLY : setup_nml_type

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: unit_in
    INTEGER :: my_comm
    INTEGER :: mpl_nml_type
    INTEGER :: ErrorStatus
    INTEGER :: icode
    CHARACTER(LEN=errormessagelength) :: iomessage
    REAL(KIND=jprb) :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_COSP'

    ! set number of each type of variable in my_namelist type
    INTEGER, PARAMETER :: no_of_types = 2
    INTEGER, PARAMETER :: n_real = 1
    INTEGER, PARAMETER :: n_log = 8

    TYPE my_namelist
      SEQUENCE
      REAL    ::  cosp_sr_cloud
      LOGICAL ::  cosp_cloudsat_sim
      LOGICAL ::  cosp_lidar_sim
      LOGICAL ::  cosp_isccp_sim
      LOGICAL ::  cosp_misr_sim
      LOGICAL ::  cosp_modis_sim
      LOGICAL ::  cosp_rttov_sim
      LOGICAL ::  cosp_use_vgrid
      LOGICAL ::  l_cosp
    END TYPE my_namelist

    TYPE (my_namelist) :: my_nml

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    CALL gc_get_communicator(my_comm, icode)

    CALL setup_nml_type(no_of_types, mpl_nml_type, n_real_in=n_real,           &
                        n_log_in=n_log)

    IF (mype == 0) THEN

      READ(UNIT=unit_in, NML=RUN_COSP, IOSTAT=ErrorStatus, IOMSG=iomessage)
      CALL check_iostat(errorstatus, "namelist RUN_COSP", iomessage)

      my_nml % cosp_sr_cloud     = cosp_sr_cloud
      my_nml % cosp_cloudsat_sim = cosp_cloudsat_sim
      my_nml % cosp_lidar_sim    = cosp_lidar_sim
      my_nml % cosp_isccp_sim    = cosp_isccp_sim
      my_nml % cosp_misr_sim     = cosp_misr_sim
      my_nml % cosp_modis_sim    = cosp_modis_sim
      my_nml % cosp_rttov_sim    = cosp_rttov_sim
      my_nml % cosp_use_vgrid    = cosp_use_vgrid
      my_nml % l_cosp            = l_cosp

    END IF

    CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

    IF (mype /= 0) THEN

      cosp_sr_cloud     = my_nml % cosp_sr_cloud
      cosp_cloudsat_sim = my_nml % cosp_cloudsat_sim
      cosp_lidar_sim    = my_nml % cosp_lidar_sim
      cosp_isccp_sim    = my_nml % cosp_isccp_sim
      cosp_misr_sim     = my_nml % cosp_misr_sim
      cosp_modis_sim    = my_nml % cosp_modis_sim
      cosp_rttov_sim    = my_nml % cosp_rttov_sim
      cosp_use_vgrid    = my_nml % cosp_use_vgrid
      l_cosp            = my_nml % l_cosp

    END IF

    CALL mpl_type_free(mpl_nml_type,icode)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE read_nml_run_cosp

END MODULE cosp_input_mod

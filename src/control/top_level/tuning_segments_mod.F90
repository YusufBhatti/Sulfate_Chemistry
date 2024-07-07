! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

!  Data module for setting segments sizes used in physics code
! Description:

! Method:
!   segment sizes used by the physics code sections
!   are defined here and assigned default values. These may
!   be overridden by namelist input.

!   Any routine wishing to use these options may do so with the 'USE'
!   statement.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top level Control

! Code Description:
!   Language: FORTRAN 95

MODULE tuning_segments_mod

USE missing_data_mod, ONLY: rmdi, imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength
USE chk_opts_mod, ONLY: chk_var, def_src
USE bl_option_mod, ONLY: i_bl_vn
USE mphys_inputs_mod, ONLY: l_rain, l_casim
USE cv_run_mod, ONLY: i_convection_vn

IMPLICIT NONE

!=======================================================================
! Segments namelist options - ordered as in meta-data sort-key, with namelist
! option and possible values
!=======================================================================

INTEGER :: bl_segment_size = imdi

INTEGER :: gw_seg_size = 32   ! Size of segments for optimising for cache
                              ! use and OpenMP
INTEGER :: ussp_seg_size = imdi ! Segment size for USSP gravity wave drag
INTEGER :: precip_segment_size = imdi

INTEGER :: a_convect_segments = imdi ! No of batches used in convection
                                     ! Original default 1

INTEGER :: a_convect_seg_size = imdi ! Size of convection segments. Can be
                                     ! specified as an alternative to no. of
                                     ! segments. Original default -99

INTEGER :: a_sw_segments = imdi  ! No of batches used in shortwave code
INTEGER :: a_sw_seg_size = -99   ! Size of sw batches, -ve disabled
INTEGER :: a_lw_segments = imdi  ! No of batches used in longwave code
INTEGER :: a_lw_seg_size = -99   ! Size of lw batches, -ve disabled

INTEGER :: ukca_mode_seg_size = imdi  ! No. of columns per chunk

LOGICAL :: l_autotune_segments = .FALSE. ! Automatic segment tuning on/off.
LOGICAL :: l_autotune_verbose  = .FALSE. ! Prints additional information from
                                         ! the tuning algorithm.
INTEGER :: autotune_trial_radius = imdi  ! Maximum change in segment size
                                         ! (in either direction).
INTEGER :: autotune_min_seg_size = imdi  ! The smallest possible segment size,
                                         ! unless smaller than granularity.
INTEGER :: autotune_max_seg_size = imdi  ! The maximum possible segment size.
INTEGER :: autotune_granularity  = imdi  ! Segment sizes will be a multiple of
                                         ! this value.
INTEGER :: autotune_num_leaders = imdi   ! Num. entries in the leader board.
INTEGER :: autotune_cooling_steps  = imdi ! Num. steps over which the
                                          ! ficitious temperature falls to
                                          ! the specified fraction of its
                                          ! initial value.
REAL    :: autotune_cooling_fraction = rmdi ! Fraction of the initial
                                            ! ficitious temperature to be
                                            ! reached after the specified
                                            ! number of steps.
REAL    :: autotune_init_time_coeff = rmdi  ! Proportional increase in initial
                                            ! elapsed time that is to have an
                                            ! acceptance probability of ~50%.

!=======================================================================
!tuning_segments namelist
!=======================================================================

NAMELIST/tuning_segments/ bl_segment_size,gw_seg_size,ussp_seg_size,       &
    precip_segment_size,a_convect_seg_size,a_convect_segments,             &
    a_lw_seg_size,a_lw_segments, a_sw_seg_size,a_sw_segments,              &
    ukca_mode_seg_size, l_autotune_segments, l_autotune_verbose,           &
    autotune_trial_radius, autotune_cooling_fraction,                      &
    autotune_cooling_steps, autotune_min_seg_size, autotune_max_seg_size,  &
    autotune_granularity, autotune_init_time_coeff, autotune_num_leaders

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TUNING_SEGMENTS_MOD'

CONTAINS

SUBROUTINE check_tuning_segments()

! Description:
!   Subroutine to apply logic checks and set control variables based on the
!   options selected in the tuning_segments namelist.

USE ereport_mod,  ONLY: ereport

IMPLICIT NONE

CHARACTER (LEN=errormessagelength) :: cmessage      ! used for ereport
INTEGER                            :: icode         ! used for ereport
CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'CHECK_TUNING_SEGMENTS'
REAL(KIND=jprb)                    :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

IF (i_bl_vn > 1) &
  CALL chk_var(bl_segment_size, 'bl_segment_size', '[1:999]' )
IF (l_rain .AND. .NOT. l_casim) &
  CALL chk_var(precip_segment_size, 'precip_segment_size', '[1:999]' )
CALL chk_var(gw_seg_size, 'gw_seg_size', '[1:999]' )
CALL chk_var(ussp_seg_size, 'ussp_seg_size', '[1:999]' )
IF (i_convection_vn /= 0)  &
  CALL chk_var(a_convect_seg_size, 'a_convect_seg_size', '[<0,1:999]' )
IF (i_convection_vn /= 0 .AND. a_convect_seg_size < 0) &
  CALL chk_var(a_convect_segments, 'a_convect_segments', '[<0,1:999]' )
CALL chk_var(a_lw_seg_size, 'a_lw_seg_size', '[<0,1:999]' )
IF ( a_lw_seg_size < 0) &
  CALL chk_var(a_lw_segments, 'a_lw_segments', '[<0,1:999]' )
CALL chk_var(a_sw_seg_size, 'a_sw_seg_size', '[<0,1:999]' )
IF ( a_sw_seg_size < 0) &
  CALL chk_var(a_sw_segments, 'a_sw_segments', '[<0,1:999]' )
CALL chk_var(ukca_mode_seg_size, 'ukca_mode_seg_size', '[<0,1:216]' )

IF ( l_autotune_segments ) THEN
  CALL chk_var(autotune_trial_radius,     &
                    'autotune_trial_radius'    , '[1:999]')
  CALL chk_var(autotune_granularity,      &
                    'autotune_granularity'     , '[1:999]')
  CALL chk_var(autotune_min_seg_size,     &
                    'autotune_min_seg_size'    , '[1:999]')
  CALL chk_var(autotune_max_seg_size,     &
                    'autotune_max_seg_size'    , '[1:999]')
  CALL chk_var(autotune_num_leaders,      &
                    'autotune_num_leaders'     , '[1:999]')
  CALL chk_var(autotune_cooling_steps,    &
                    'autotune_cooling_steps'   , '[1:999]')
  CALL chk_var(autotune_cooling_fraction, &
                    'autotune_cooling_fraction', '[0.0:1.0]')
  CALL chk_var(autotune_init_time_coeff,  &
                    'autotune_init_time_coeff' , '[1.0:2.0]')
END IF

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_tuning_segments

SUBROUTINE print_nlist_tuning_segments()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_TUNING_SEGMENTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist tuning_segments',                    &
    src='tuning_segments_mod')

!Segment size / number of segments
WRITE(lineBuffer,'(A,I0)') 'bl_segment_size = ',bl_segment_size
CALL umPrint(lineBuffer,src='tuning_segments_mod')
WRITE(lineBuffer,'(A,I0)') 'gw_seg_size = ',gw_seg_size
CALL umPrint(lineBuffer,src='tuning_segments_mod')
WRITE(lineBuffer,'(A,I0)') 'ussp_seg_size = ',ussp_seg_size
CALL umPrint(lineBuffer,src='tuning_segments_mod')
WRITE(lineBuffer,'(A,I0)') 'precip_segment_size = ',precip_segment_size
CALL umPrint(lineBuffer,src='tuning_segments_mod')
WRITE(lineBuffer,'(A,I0)') 'a_convect_seg_size = ',a_convect_seg_size
CALL umPrint(lineBuffer,src='tuning_segments_mod')
WRITE(lineBuffer,'(A,I0)') 'a_convect_segments = ',a_convect_segments
CALL umPrint(lineBuffer,src='tuning_segments_mod')
WRITE(lineBuffer,'(A,I0)') 'a_lw_seg_size = ',a_lw_seg_size
CALL umPrint(lineBuffer,src='tuning_segments_mod')
WRITE(lineBuffer,'(A,I0)') 'a_lw_segments = ',a_lw_segments
CALL umPrint(lineBuffer,src='tuning_segments_mod')
WRITE(lineBuffer,'(A,I0)') 'a_sw_seg_size = ',a_sw_seg_size
CALL umPrint(lineBuffer,src='tuning_segments_mod')
WRITE(lineBuffer,'(A,I0)') 'a_sw_segments = ',a_sw_segments
CALL umPrint(lineBuffer,src='tuning_segments_mod')
WRITE(lineBuffer,'(A,I0)') 'ukca_mode_seg_size = ',ukca_mode_seg_size
CALL umPrint(lineBuffer,src='tuning_segments_mod')

!Automatic segment size tuning
WRITE(lineBuffer,'(A,L1)') 'l_autotune_segments = ',l_autotune_segments
CALL umPrint(lineBuffer,src='tuning_segments_mod')

IF (l_autotune_segments) THEN
  WRITE(lineBuffer,'(A,L1)') 'l_autotune_verbose = ',l_autotune_verbose
  CALL umPrint(lineBuffer,src='tuning_segments_mod')
  WRITE(lineBuffer,'(A,I0)') 'autotune_trial_radius = ',autotune_trial_radius
  CALL umPrint(lineBuffer,src='tuning_segments_mod')
  WRITE(lineBuffer,'(A,I0)') 'autotune_min_seg_size = ',autotune_min_seg_size
  CALL umPrint(lineBuffer,src='tuning_segments_mod')
  WRITE(lineBuffer,'(A,I0)') 'autotune_max_seg_size = ',autotune_max_seg_size
  CALL umPrint(lineBuffer,src='tuning_segments_mod')
  WRITE(lineBuffer,'(A,I0)') 'autotune_granularity = ',autotune_granularity
  CALL umPrint(lineBuffer,src='tuning_segments_mod')
  WRITE(lineBuffer,'(A,I0)') 'autotune_num_leaders = ',autotune_num_leaders
  CALL umPrint(lineBuffer,src='tuning_segments_mod')
  WRITE(lineBuffer,'(A,I0)') 'autotune_cooling_steps = ',autotune_cooling_steps
  CALL umPrint(lineBuffer,src='tuning_segments_mod')
  WRITE(lineBuffer,'(A,ES14.5)') 'autotune_cooling_fraction = ', &
        autotune_cooling_fraction
  CALL umPrint(lineBuffer,src='tuning_segments_mod')
  WRITE(lineBuffer,'(A,ES14.5)') 'autotune_init_time_coeff = ', &
        autotune_init_time_coeff
  CALL umPrint(lineBuffer,src='tuning_segments_mod')
END IF

CALL umPrint('- - - - - - end of namelist - - - - - -',                 &
    src='tuning_segments_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_tuning_segments

#if !defined(LFRIC)
SUBROUTINE read_nml_tuning_segments(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: errorstatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_TUNING_SEGMENTS'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int  = 17
INTEGER, PARAMETER :: n_real = 2
INTEGER, PARAMETER :: n_log  = 2

TYPE my_namelist
  SEQUENCE
  INTEGER :: bl_segment_size
  INTEGER :: gw_seg_size
  INTEGER :: ussp_seg_size
  INTEGER :: precip_segment_size
  INTEGER :: a_convect_seg_size
  INTEGER :: a_convect_segments
  INTEGER :: a_lw_seg_size
  INTEGER :: a_lw_segments 
  INTEGER :: a_sw_seg_size
  INTEGER :: a_sw_segments 
  INTEGER :: ukca_mode_seg_size
  INTEGER :: autotune_trial_radius
  INTEGER :: autotune_min_seg_size
  INTEGER :: autotune_max_seg_size
  INTEGER :: autotune_granularity
  INTEGER :: autotune_num_leaders
  INTEGER :: autotune_cooling_steps
  REAL    :: autotune_cooling_fraction
  REAL    :: autotune_init_time_coeff
  LOGICAL :: l_autotune_segments
  LOGICAL :: l_autotune_verbose
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type,   &
                    n_int_in=n_int, n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=tuning_segments, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist TUNING_SEGMENTS", iomessage)

  my_nml % bl_segment_size            = bl_segment_size
  my_nml % gw_seg_size                = gw_seg_size
  my_nml % ussp_seg_size              = ussp_seg_size
  my_nml % precip_segment_size        = precip_segment_size
  my_nml % a_convect_seg_size         = a_convect_seg_size
  my_nml % a_convect_segments         = a_convect_segments
  my_nml % a_lw_seg_size              = a_lw_seg_size
  my_nml % a_lw_segments              = a_lw_segments
  my_nml % a_sw_seg_size              = a_sw_seg_size
  my_nml % a_sw_segments              = a_sw_segments
  my_nml % ukca_mode_seg_size         = ukca_mode_seg_size
  my_nml % autotune_trial_radius      = autotune_trial_radius
  my_nml % autotune_min_seg_size      = autotune_min_seg_size
  my_nml % autotune_max_seg_size      = autotune_max_seg_size
  my_nml % autotune_granularity       = autotune_granularity
  my_nml % autotune_num_leaders       = autotune_num_leaders
  my_nml % autotune_cooling_steps     = autotune_cooling_steps
  ! end of integers

  my_nml % autotune_cooling_fraction  = autotune_cooling_fraction
  my_nml % autotune_init_time_coeff   = autotune_init_time_coeff
  ! end of reals

  my_nml % l_autotune_segments        = l_autotune_segments
  my_nml % l_autotune_verbose         = l_autotune_verbose
  ! end of logicals

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  bl_segment_size            = my_nml % bl_segment_size
  gw_seg_size                = my_nml % gw_seg_size
  ussp_seg_size              = my_nml % ussp_seg_size
  precip_segment_size        = my_nml % precip_segment_size
  a_convect_seg_size         = my_nml % a_convect_seg_size
  a_convect_segments         = my_nml % a_convect_segments
  a_lw_seg_size              = my_nml % a_lw_seg_size
  a_lw_segments              = my_nml % a_lw_segments
  a_sw_seg_size              = my_nml % a_sw_seg_size
  a_sw_segments              = my_nml % a_sw_segments
  ukca_mode_seg_size         = my_nml % ukca_mode_seg_size
  autotune_trial_radius      = my_nml % autotune_trial_radius
  autotune_min_seg_size      = my_nml % autotune_min_seg_size
  autotune_max_seg_size      = my_nml % autotune_max_seg_size
  autotune_granularity       = my_nml % autotune_granularity
  autotune_num_leaders       = my_nml % autotune_num_leaders
  autotune_cooling_steps     = my_nml % autotune_cooling_steps
  ! end of integers

  autotune_cooling_fraction = my_nml % autotune_cooling_fraction
  autotune_init_time_coeff  = my_nml % autotune_init_time_coeff
  ! end of reals

  l_autotune_segments        = my_nml % l_autotune_segments
  l_autotune_verbose         = my_nml % l_autotune_verbose
  ! end of logicals

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_tuning_segments
#endif

END MODULE tuning_segments_mod

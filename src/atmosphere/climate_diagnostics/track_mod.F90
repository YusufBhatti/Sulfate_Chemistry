! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Global data module for switches and options for FV-TRACK.
!
!     Code Owner: Please refer to the UM file CodeOwners.txt
!     This file belongs in section: Climate Diagnostics

MODULE track_mod

USE missing_data_mod,        ONLY: rmdi, imdi
USE umPrintMgr,              ONLY: umprint, ummessage
USE errormessagelength_mod,  ONLY: errormessagelength
USE yomhook,                 ONLY: lhook, dr_hook
USE parkind1,                ONLY: jprb, jpim

IMPLICIT NONE

! Params
LOGICAL, PARAMETER :: tc_tracking_ON  = .TRUE.
LOGICAL, PARAMETER :: tc_tracking_OFF = .FALSE.

! Module variables
LOGICAL ::                                                               &
    l_hoskins = .FALSE.
        ! Flag to activate hoskins filter
INTEGER ::                                                               &
    nbot_850 = imdi                                                      &
        ! Truncate all wavenumbers below this one for 850hPa field
,   ntop_850 = imdi                                                      &
        ! Truncate all wavenumbers above this value
,   ntop_tc = imdi                                                       &
        ! Truncate all wavenumbers above this value
,   nlevs_avg = imdi
        ! Number of vorticity levels (pressure) to average
        ! for the TC tracking  at ntop_tc.

REAL :: sm = rmdi ! Hoskin filter value at truncation


! Define Namelist RUN_track
NAMELIST/Run_Track/                                                      &
l_hoskins, nbot_850, ntop_850, ntop_tc, nlevs_avg, sm

!-------------------------------------------------------------------------
!   winds for new TC method
REAL, ALLOCATABLE ::                                                     &
       u_tc_new(:,:,:)                                                   &
        ! U component to store the nlevs_avg levels for the TC-new method
,      v_tc_new(:,:,:)                                                   &
        ! V component to store the nlevs_avg levels for the TC-new method
,      Ymn(:,:,:)                                                        &
        ! Matrix to store the Associate Legendre Polynomials for each lat.
,      Fm(:,:,:)                                                         &
        ! Matrix for SH projections for 850hPa filtering
,      Fm_TC(:,:,:)                                                      &
        ! Matrix for SH projections for TC filtering
,      cos_ki(:,:)                                                       &
        ! Cosine functions for freq and longitudes
,      sin_ki(:,:)                                                       &
        ! Sine functions for freq and longitudes
,      cos_ki_T(:,:)                                                     &
        ! Transpose of Cosine functions
,      sin_ki_T(:,:)
        !  Transpose of Sine functions

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TRACK_MOD'

CONTAINS

SUBROUTINE print_nlist_run_track()
! Print out the run_track namelist on the out file.

IMPLICIT NONE

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_TRACK'
!End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist Run_Track', src='track_mod')

WRITE(umMessage,'(A,L3)')'l_hoskins = ',l_hoskins
CALL umPrint(umMessage,src='track_mod')
WRITE(umMessage,'(A,I3)')'nbot_850 = ',nbot_850
CALL umPrint(umMessage,src='track_mod')
WRITE(umMessage,'(A,I3)')'ntop_850 = ',ntop_850
CALL umPrint(umMessage,src='track_mod')
WRITE(umMessage,'(A,I3)')'ntop_tc = ',ntop_tc
CALL umPrint(umMessage,src='track_mod')
WRITE(umMessage,'(A,I3)')'nlevs_avg = ',nlevs_avg
CALL umPrint(umMessage,src='track_mod')
WRITE(umMessage,'(A,ES12.4)') 'sm = ',sm
CALL umPrint(umMessage,src='track_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', src='track_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_track


SUBROUTINE read_nml_run_track(unit_in)
! Read in namelist values from run_track input array
USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE
! Argument variables
INTEGER, INTENT(IN) :: unit_in

! MPL and error handling variables
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
! DrHook handle
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_TRACK'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 4
INTEGER, PARAMETER :: n_real = 1
INTEGER, PARAMETER :: n_log = 1

! Define my_namelist
TYPE my_namelist
  SEQUENCE
  INTEGER :: nbot_850
  INTEGER :: ntop_850
  INTEGER :: ntop_tc
  INTEGER :: nlevs_avg
  REAL :: sm
  LOGICAL :: l_hoskins
END TYPE my_namelist

TYPE (my_namelist) :: my_nml
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                        n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_track, IOSTAT=ErrorStatus,            &
       IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_track", iomessage)

  my_nml % nbot_850     = nbot_850
  my_nml % ntop_850     = ntop_850
  my_nml % ntop_tc      = ntop_tc
  my_nml % nlevs_avg    = nlevs_avg
  my_nml % sm           = sm
  my_nml % l_hoskins    = l_hoskins
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN
  nbot_850     = my_nml % nbot_850
  ntop_850     = my_nml % ntop_850
  ntop_tc      = my_nml % ntop_tc
  nlevs_avg    = my_nml % nlevs_avg
  sm           = my_nml % sm
  l_hoskins    = my_nml % l_hoskins
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_track

SUBROUTINE destroy_track_arr()
! Deallocates saved arrays employed by FV-TRACK

! DrHook modules
USE parkind1,      ONLY: jpim, jprb
USE yomhook,       ONLY: lhook, dr_hook

IMPLICIT NONE
!Local variables (DrHook only)
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DESTROY_TRACK_ARR'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Deallocate Fm matrices
IF (ALLOCATED(Fm))    DEALLOCATE ( Fm )
IF (ALLOCATED(Fm_TC)) DEALLOCATE ( Fm_TC )

! Deallocate trigonometric matrices for fft_track
IF (ALLOCATED(sin_ki)) DEALLOCATE ( sin_ki )
IF (ALLOCATED(cos_ki)) DEALLOCATE ( cos_ki )

IF (ALLOCATED(sin_ki_T)) DEALLOCATE ( sin_ki_T )
IF (ALLOCATED(cos_ki_T)) DEALLOCATE ( cos_ki_T )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE destroy_track_arr

END MODULE track_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for LBC reading

! Description:
!   Module containing data concerning LBC reading
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: LBC input

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE lbc_read_data_mod

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE lbc_mod,  ONLY: nrim_timesa
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Items from the retired cbound.h include file

INTEGER, PARAMETER :: max_rimwidth = 30    ! Maximum rim width for LBCs

INTEGER :: current_lbc_step         ! Timestep at which LBCs were last updated
INTEGER :: albc_num                 ! Number of atmos boundary file currently
                                    !  in use
INTEGER :: albc2_starttime_steps    ! VT of first block of data in 2nd atmos
                                    !  boundary file, in steps from start of run

INTEGER :: albc_swapstep            ! Step on which to swap to second atmos
                                    ! boundary_file

!-----------------------------------------------------
! Items set in namelist
!-----------------------------------------------------

! Rim weights
REAL :: rimweightsa(max_rimwidth)

! Switch for interpolated winds in lbcs
! True for advecting winds interpolated in boundary zone
LOGICAL :: l_int_uvw_lbc = .FALSE.

! Limited area model uses lateral boundary tendencies
LOGICAL :: l_lateral_boundary = .FALSE.

LOGICAL :: L_old_lbc_file = .FALSE.  ! .true. if using a 13 field
                                      ! LBC file - ensures legacy LBC
                                      ! files will still work with advected
                                      ! wind changes

NAMELIST /lbc_options/ l_int_uvw_lbc, l_lateral_boundary,  &
                       rimweightsa, l_old_lbc_file, nrim_timesa

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LBC_READ_DATA_MOD'

CONTAINS

SUBROUTINE print_nlist_lbc_options()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_LBC_OPTIONS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist lbc_options', &
     src='lbc_read_data_mod')

WRITE(lineBuffer,*) ' l_int_uvw_lbc = ',l_int_uvw_lbc
CALL umPrint(lineBuffer,src='lbc_read_data_mod')
WRITE(lineBuffer,*) ' l_lateral_boundary = ',l_lateral_boundary
CALL umPrint(lineBuffer,src='lbc_read_data_mod')
WRITE(lineBuffer,*) ' rimweightsa = ',rimweightsa
CALL umPrint(lineBuffer,src='lbc_read_data_mod')
WRITE(lineBuffer,*) ' l_old_lbc_file = ',l_old_lbc_file
CALL umPrint(lineBuffer,src='lbc_read_data_mod')
WRITE(lineBuffer,*) ' nrim_timesa = ',nrim_timesa
CALL umPrint(lineBuffer,src='lbc_read_data_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
     src='lbc_read_data_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_lbc_options

SUBROUTINE read_nml_lbc_options(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_LBC_OPTIONS'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_real = max_rimwidth
INTEGER, PARAMETER :: n_log = 3
INTEGER, PARAMETER :: n_int = 1

CHARACTER(LEN=errormessagelength) :: iomessage

TYPE my_namelist
  SEQUENCE
  REAL :: rimweightsa(max_rimwidth)
  LOGICAL :: l_int_uvw_lbc
  LOGICAL :: l_lateral_boundary
  LOGICAL :: l_old_lbc_file
  INTEGER :: nrim_timesa
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,             &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=lbc_options, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist lbc_options", iomessage)

  my_nml % rimweightsa        = rimweightsa
  my_nml % l_int_uvw_lbc      = l_int_uvw_lbc
  my_nml % l_lateral_boundary = l_lateral_boundary
  my_nml % l_old_lbc_file     = l_old_lbc_file
  my_nml % nrim_timesa        = nrim_timesa

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  rimweightsa        = my_nml % rimweightsa
  l_int_uvw_lbc      = my_nml % l_int_uvw_lbc
  l_lateral_boundary = my_nml % l_lateral_boundary
  l_old_lbc_file     = my_nml % l_old_lbc_file
  nrim_timesa        = my_nml % nrim_timesa

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_lbc_options

END MODULE lbc_read_data_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

MODULE nlcfiles_namelist_mod

! Description:
!   Module which defines the nlcfiles namelist and associated read
!   routine
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc

USE filenamelength_mod, ONLY: filenamelength
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=filenamelength), SAVE :: &
    alabcin1 = "", &  ! LBC input file #1
    alabcin2 = "", &  ! LBC input file #2
    astart = "",   &  ! Atmos start dump / Recon output dump
    atmanl = "",   &  ! Atmos analysis dump
    iau_inc = "",  &  ! IAU increment file
    obs01 = "",    &  ! Observation file directory #1
    obs02 = "",    &  ! Observation file directory #2
    obs03 = "",    &  ! Observation file directory #3
    obs04 = "",    &  ! Observation file directory #4
    obs05 = "",    &  ! Observation file directory #5
    rp2_seed = "", &  ! Random seed for RP2 (stochastic physics)
    streqlog = ""     ! STASH Requests log file

NAMELIST / nlcfiles /  &
    alabcin1, &
    alabcin2, &
    astart,   &
    atmanl,   &
    iau_inc,  &
    obs01,    &
    obs02,    &
    obs03,    &
    obs04,    &
    obs05,    &
    rp2_seed, & 
    streqlog

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NLCFILES_NAMELIST_MOD'

CONTAINS

SUBROUTINE read_nml_nlcfiles(unithist)

USE um_parcore, ONLY: mype
USE setup_namelist, ONLY: setup_nml_type
USE check_iostat_mod, ONLY: check_iostat
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: unithist
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: errorstatus
INTEGER :: icode

CHARACTER(LEN=errormessagelength) :: iomessage

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLCFILES'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_chars =  filenamelength * 12

TYPE my_namelist
  SEQUENCE
  CHARACTER(LEN=filenamelength) :: alabcin1
  CHARACTER(LEN=filenamelength) :: alabcin2
  CHARACTER(LEN=filenamelength) :: astart
  CHARACTER(LEN=filenamelength) :: atmanl
  CHARACTER(LEN=filenamelength) :: iau_inc
  CHARACTER(LEN=filenamelength) :: obs01
  CHARACTER(LEN=filenamelength) :: obs02
  CHARACTER(LEN=filenamelength) :: obs03
  CHARACTER(LEN=filenamelength) :: obs04
  CHARACTER(LEN=filenamelength) :: obs05
  CHARACTER(LEN=filenamelength) :: rp2_seed
  CHARACTER(LEN=filenamelength) :: streqlog
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_chars_in=n_chars)

IF (mype == 0) THEN

  READ(UNIT=unithist, NML=nlcfiles, IOSTAT=errorstatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist NLCFILES", iomessage)
  
  my_nml % alabcin1 = alabcin1
  my_nml % alabcin2 = alabcin2
  my_nml % astart   = astart
  my_nml % atmanl   = atmanl
  my_nml % iau_inc  = iau_inc
  my_nml % obs01    = obs01
  my_nml % obs02    = obs02
  my_nml % obs03    = obs03
  my_nml % obs04    = obs04
  my_nml % obs05    = obs05
  my_nml % rp2_seed = rp2_seed
  my_nml % streqlog = streqlog

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  alabcin1 = my_nml % alabcin1
  alabcin2 = my_nml % alabcin2
  astart   = my_nml % astart
  atmanl   = my_nml % atmanl
  iau_inc  = my_nml % iau_inc
  obs01    = my_nml % obs01
  obs02    = my_nml % obs02
  obs03    = my_nml % obs03
  obs04    = my_nml % obs04
  obs05    = my_nml % obs05
  rp2_seed = my_nml % rp2_seed
  streqlog = my_nml % streqlog

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlcfiles


SUBROUTINE print_nlist_nlcfiles()

USE umprintMgr, ONLY: umprint

IMPLICIT NONE

CHARACTER(LEN=50000) :: linebuffer
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_NLCFILES'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umprint('Contents of namelist nlcfiles', src=modulename)

WRITE(linebuffer,"(A,A)")' alabcin1 = ', TRIM(alabcin1)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' alabcin2 = ', TRIM(alabcin2)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' astart   = ', TRIM(astart)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' atmanl   = ', TRIM(atmanl)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' iau_inc  = ', TRIM(iau_inc)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' obs01    = ', TRIM(obs01)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' obs02    = ', TRIM(obs02)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' obs03    = ', TRIM(obs03)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' obs04    = ', TRIM(obs04)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' obs05    = ', TRIM(obs05)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' rp2_seed = ', TRIM(rp2_seed)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' streqlog = ', TRIM(streqlog)
CALL umprint(linebuffer,src=modulename)

CALL umprint('- - - - - - end of namelist - - - - - -', src=modulename)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_nlcfiles

END MODULE nlcfiles_namelist_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Define headers namelist

MODULE rcf_headers_mod

! Description:
!    Define headers namelist
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE missing_data_mod, ONLY: imdi,rmdi

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE

PUBLIC ::   read_nml_headers,      &
            rcf_override_headers,  &
            print_nlist_headers,   &
            check_nml_headers,     &
            inthd,                 &
            relhd,                 &
            fixhd

! options for how to override the date and time
INTEGER :: i_override_date_time = imdi
INTEGER, PARAMETER :: date_from_dump = 0
INTEGER, PARAMETER :: year_only = 1
INTEGER, PARAMETER :: full_date_and_time = 2

! new date and time to use
INTEGER :: new_date_time(6) = imdi

NAMELIST /Headers/ i_override_date_time, new_date_time

! Items which may be set via code branch
INTEGER :: FixHd(256) = imdi
INTEGER :: IntHd(100) = imdi
REAL    :: RelHd(100) = rmdi

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_HEADERS_MOD'

CONTAINS

SUBROUTINE rcf_override_headers()

IMPLICIT NONE

SELECT CASE(i_override_date_time)

CASE(year_only)
  ! year
  FixHD(21) = new_date_time(1)
  FixHD(28) = new_date_time(1)

CASE(full_date_and_time)
  ! year
  FixHD(21) = new_date_time(1)
  FixHD(28) = new_date_time(1)
  ! month
  FixHD(22) = new_date_time(2)
  FixHD(29) = new_date_time(2)
  ! day
  FixHD(23) = new_date_time(3)
  FixHD(30) = new_date_time(3)
  ! hour
  FixHD(24) = new_date_time(4)
  FixHD(31) = new_date_time(4)
  ! minute
  FixHD(25) = new_date_time(5)
  FixHD(32) = new_date_time(5)
  ! second
  FixHD(26) = new_date_time(6)
  FixHD(33) = new_date_time(6)

END SELECT

END SUBROUTINE rcf_override_headers

SUBROUTINE print_nlist_headers()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist headers', &
        src='rcf_headers_mod')

WRITE(lineBuffer,*)' i_override_date_time = ',i_override_date_time
CALL umPrint(lineBuffer,src='rcf_headers_mod')
WRITE(lineBuffer,*)' new_date_time = ',new_date_time
CALL umPrint(lineBuffer,src='rcf_headers_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
        src='rcf_headers_mod')

END SUBROUTINE print_nlist_headers

SUBROUTINE read_nml_headers(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat
USE errormessagelength_mod, ONLY: errormessagelength

USE setup_namelist, ONLY: setup_nml_type


IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode

CHARACTER(LEN=errormessagelength) :: iomessage
    
! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_int = 1 + 6

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_override_date_time
  INTEGER :: new_date_time(6)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

CHARACTER (LEN=*), PARAMETER  :: RoutineName='READ_NML_HEADERS'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=headers, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist headers", iomessage)

  my_nml % i_override_date_time = i_override_date_time
  my_nml % new_date_time        = new_date_time

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  i_override_date_time  = my_nml % i_override_date_time
  new_date_time         = my_nml % new_date_time

END IF

CALL mpl_type_free(mpl_nml_type,icode)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_headers

SUBROUTINE check_nml_headers()
USE chk_opts_mod, ONLY: chk_var, def_src

IMPLICIT NONE
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'CHECK_NML_HEADERS'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

! Check option for overriding date and time values in dump
CALL chk_var( i_override_date_time, 'i_override_date_time',             &
              [date_from_dump, year_only, full_date_and_time],          &
              cmessage='Unrecognised option for overriding time.')

! Which parts of the new time need to be valid varies depending on what
! parts of the timestamps are being overridden. If time taken from dump,
! the values in the namelist aren't used and don't need to be checked.
SELECT CASE(i_override_date_time)

CASE(year_only)
  CALL chk_var( new_date_time(1), 'new_date_time(1)', '[>0]',           &
       cmessage='Please enter a valid year')

CASE(full_date_and_time)
  CALL chk_var( new_date_time(1), 'new_date_time(1)', '[>0]',           &
       cmessage='Please enter a valid year')
  CALL chk_var( new_date_time(2), 'new_date_time(2)', '[1:12]',         &
       cmessage='Please enter a valid month')
  CALL chk_var( new_date_time(3), 'new_date_time(3)', '[1:31]',         &
       cmessage='Please enter a valid day')
  CALL chk_var( new_date_time(4), 'new_date_time(4)', '[0:23]',         &
       cmessage='Please enter a valid hour')
  CALL chk_var( new_date_time(5), 'new_date_time(5)', '[0:59]',         &
       cmessage='Please enter a valid minute')
  CALL chk_var( new_date_time(6), 'new_date_time(6)', '[0:59]',         &
       cmessage='Please enter a valid second')
END SELECT

def_src = ''

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nml_headers

END MODULE rcf_headers_mod

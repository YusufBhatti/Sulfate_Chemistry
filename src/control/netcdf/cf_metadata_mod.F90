! *****************************COPYRIGHT*******************************
! 
! Copyright 2017-2018 University of Reading
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice, 
! this list of conditions and the following disclaimer in the documentation 
! and/or other materials provided with the distribution.
! 
! 3. Neither the name of the copyright holder nor the names of its contributors 
! may be used to endorse or promote products derived from this software without 
! specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: NetCDF output
!
MODULE cf_metadata_mod

USE mpl,                    ONLY: mpl_character
USE yomhook,                ONLY: lhook,dr_hook
USE parkind1,               ONLY: jprb,jpim
USE um_parcore,             ONLY: mype
USE umprintmgr,             ONLY: umprint,ummessage,newline,printstatus,prdiag
USE ereport_mod,            ONLY: ereport
USE get_env_var_mod,        ONLY: get_env_var
USE version_mod,            ONLY: nsectp,nitemp
USE filenamelength_mod,     ONLY: filenamelength
USE errormessagelength_mod, ONLY: errormessagelength
USE file_manager,           ONLY: assign_file_unit,release_file_unit
USE stash_model_mod,        ONLY: h_a_polelat,h_a_polelong
USE submodel_mod,           ONLY: atmos_im
USE um_version_mod,         ONLY: um_version_int

IMPLICIT NONE

INTEGER, PARAMETER :: maxfieldlen = 100 ! Maximum length of field entry

CHARACTER(LEN=*), PARAMETER  :: conventions = 'CF-1.6'
                                ! Conventions attribute string
CHARACTER(LEN=maxfieldlen)   :: standard_name(0:nsectp,nitemp)
                                ! standard_name attribute string
CHARACTER(LEN=maxfieldlen)   :: units(0:nsectp,nitemp)
                                ! units attribute string
CHARACTER(LEN=maxfieldlen)   :: cf_extra_info(0:nsectp,nitemp)
                                ! Extra CF information string

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CF_METADATA_MOD'

CONTAINS

! Read stash2cf file and extract cf metadata

SUBROUTINE get_std_name_units()

IMPLICIT NONE

!
!  Local variables
!
INTEGER, PARAMETER :: nseparator = 5    ! Number of separators per line
INTEGER, PARAMETER :: maxlinelen = 256  ! Maximum line length in STASH2CF file
INTEGER, PARAMETER :: maxsectionlen = 3 ! Maximum length of stash code entry
INTEGER, PARAMETER :: maxitemlen = 3    ! Maximum length of stash code entry

INTEGER :: ErrorStatus           ! Error return code
INTEGER :: i                     ! Loop index
INTEGER :: my_comm               ! MPI communicator
INTEGER :: seppos(nseparator)    ! Array containing position of field separators
INTEGER :: pos1                  ! Variable to temporarily hold line position
INTEGER :: iunit                 ! Assigned Fortran unit number
INTEGER :: section               ! STASH section code
INTEGER :: item                  ! STASH item code

LOGICAL :: lrotgrid ! TRUE if rotated grid

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'GET_STD_NAME_UNITS'
CHARACTER(LEN=errormessagelength) :: message   ! Used for ereport
CHARACTER(LEN=errormessagelength) :: iomessage ! IO message from file open
CHARACTER(LEN=filenamelength)  :: stash2cf     ! STASH to CF meta-data file name
CHARACTER(LEN=maxlinelen)      :: line         ! Line read from meta-data file
CHARACTER(LEN=maxsectionlen)   :: csection     ! STASH section number string
CHARACTER(LEN=maxitemlen)      :: citem        ! STASH item number string
CHARACTER(LEN=maxfieldlen)     :: pp_extra_info! Info on valid grid rotation

! Dr-Hook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Default values of CF strings
standard_name = 'undefined'
units = 'undefined'
cf_extra_info = ''

! Read file on PE 0 only
IF (mype == 0) THEN

  lrotgrid = .NOT. (h_a_polelat == 90.0 .AND. h_a_polelong == 0.0)

  ! Get name for stash2cf meta-data file
  CALL get_env_var("STASH2CF",stash2cf)

  WRITE(umMessage,'(A,A)') 'Reading STASH2CF file ',TRIM(stash2cf)
  CALL umPrint(umMessage,src=RoutineName)

  ! Open stash2cf meta-data file
  CALL assign_file_unit(stash2cf,iunit,handler="fortran")
  OPEN(UNIT=iunit,FILE=stash2cf,ACTION='READ',STATUS='OLD', &
       IOSTAT=ErrorStatus,IOMSG=iomessage)
  IF (ErrorStatus >  0) THEN
    WRITE(message,'(A,A,A,A)') &
          'Failed in OPEN of stash2cf meta-data file',newline, &
          'IoMsg: ',TRIM(iomessage)
    CALL ereport(RoutineName,ErrorStatus,message)
  ELSE IF (ErrorStatus <  0) THEN
    WRITE(message,'(A,A,A,A)') &
          'Warning message on OPEN of stash2cf meta-data file',newline, &
          'IoMsg: ',TRIM(iomessage)
    CALL ereport(RoutineName,ErrorStatus,message)
  END IF

  ! Loop over all entries in STASH2CF file
  DO
    READ(UNIT=iunit,FMT='(A)',IOSTAT=ErrorStatus)line
    IF (ErrorStatus /= 0) EXIT

    IF (line(1:1) == '#') CYCLE ! Comment line

    ! Save position of field separators
    pos1 = 0
    DO i=1,nseparator
      seppos(i) = pos1 + SCAN(line(pos1+1:),'!')
      pos1 = seppos(i)
    END DO

    ! Extract fields necesssary to determine if
    ! this entry is needed in this model configuaration
    pp_extra_info = line(seppos(5)+1:)

    ! Either this entry applies to both rotated and unrotated grids or
    ! only select entry for correct grid rotation
    IF (pp_extra_info == '' .OR. &
        ((lrotgrid .AND. pp_extra_info == 'rotated_latitude_longitude') .OR. &
        (.NOT. lrotgrid .AND. pp_extra_info == 'true_latitude_longitude'))) THEN

      ! Extract STASH section and item numbers
      csection = line(1:seppos(1)-1)
      READ(csection,'(i3)') section
      citem = line(seppos(1)+1:seppos(2)-1)
      READ(citem,'(i3)') item

      ! Extract standard_name, units and cf_extra values
      IF (seppos(3)-seppos(2) > 1) THEN
        standard_name(section,item) = TRIM(line(seppos(2)+1:seppos(3)-1))
      END IF
      IF (seppos(4)-seppos(3) > 1) THEN
        units(section,item) = TRIM(line(seppos(3)+1:seppos(4)-1))
      END IF
      IF (seppos(5)-seppos(4) > 1) THEN
        cf_extra_info(section,item) = TRIM(line(seppos(4)+1:seppos(5)-1))
      END IF
    END IF
  END DO

  ! Close stash2cf meta-data file
  CLOSE(UNIT=iunit)
  CALL release_file_unit(iunit, HANDLER="fortran")

END IF

CALL gc_get_communicator(my_comm,ErrorStatus)

! Broadcast CF meta-data information to all other PEs
CALL mpl_bcast(standard_name,maxfieldlen*(nsectp+1)*nitemp,mpl_character, &
               0,my_comm,ErrorStatus)
CALL mpl_bcast(units,maxfieldlen*(nsectp+1)*nitemp,mpl_character, &
               0,my_comm,ErrorStatus)
CALL mpl_bcast(cf_extra_info,maxfieldlen*(nsectp+1)*nitemp,mpl_character, &
               0,my_comm,ErrorStatus)

IF (printstatus >= prdiag) THEN
  DO section=0,nsectp
    DO item=1,nitemp
      IF (standard_name(section,item) /= 'undefined') THEN
        WRITE(umMessage,'(I4,I4,A)')section,item, &
                                    ' '//TRIM(standard_name(section,item))// &
                                    ' '//TRIM(units(section,item))// &
                                    ' '//TRIM(cf_extra_info(section,item))
        CALL umPrint(umMessage,src=RoutineName)
      END IF
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_std_name_units

END MODULE cf_metadata_mod


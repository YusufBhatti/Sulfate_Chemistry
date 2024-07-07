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
! Description
! When a stream is due to be reinitialised, the old file is closed,
! the post-processing system is informed if requested, and a new file
! is initialised.
!
!  Method : 1. Files to be processed on each call are controlled by
!              the file unit reinitialisation switch set by ncfile_init 
!              and at regular intervals thereafter.
!           2. Completed files are closed, and if post-processing is
!              selected a zero-length indicator .arch file is created
!           3. Where streams are to be reinitialised, a new file is
!              created with an appropriate name
!
! Code Description:
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

MODULE ncfile_reinit_mod

USE yomhook,                 ONLY: lhook,dr_hook
USE parkind1,                ONLY: jprb,jpim
USE um_parcore,              ONLY: mype
USE umnetcdf_mod,            ONLY:                                     &
    nc_file_open,nc_file_close,nc_is_file_open,nc_mode_create
USE filenamelength_mod,      ONLY: filenamelength
USE io_configuration_mod,    ONLY: l_postp
USE submodel_mod,            ONLY: atmos_im,internal_id_max
USE iau_mod,                 ONLY: l_iau
USE umPrintMgr,              ONLY: umPrint,umMessage
USE filename_generation_mod, ONLY: get_filename
USE nlstcall_mod,            ONLY: model_analysis_mins,l_fastrun   
USE model_time_mod,          ONLY: secs_per_stepim
USE file_manager,            ONLY:                                     &
    init_file_loop,um_file_type,assign_file_unit,release_file_unit
USE init_nc_mod,             ONLY: init_nc
USE init_stash_nc_mod,       ONLY: init_stash_nc

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NCFILE_REINIT_MOD'

CONTAINS

SUBROUTINE ncfile_reinit

IMPLICIT NONE

!
! Local variables
!
INTEGER :: name_offset ! Character offset in file name generation
INTEGER :: arch_unit   ! Fortran unit for archiving flag file

CHARACTER(LEN=filenamelength) :: oldncfile ! Previous netCDF file on unit
CHARACTER(LEN=filenamelength) :: archfile  ! Archiving file

CHARACTER(LEN=*), PARAMETER :: RoutineName = "NCFILE_REINIT"

TYPE(um_file_type), POINTER :: um_file ! Pointer to file objects

! Dr-Hook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0) THEN

  ! If the run is using an analysis calculate the offset to use in any
  ! relative filenames
  IF (l_fastrun .OR. l_iau) THEN
    name_offset = - (model_analysis_mins * 60)/secs_per_stepim(atmos_im)
  ELSE
    name_offset = 0
  END IF

  NULLIFY(um_file)
  um_file => init_file_loop(handler="netcdf")
  DO WHILE (ASSOCIATED(um_file))

    IF (um_file % meta % initialise) THEN

      !  Save the name of the existing file before regenerating it
      oldncfile = um_file % filename

      ! Generate the archiving flag file
      IF (l_postp) THEN
        IF (INDEX(oldncfile, "Reserved unit") == 0) THEN
          archfile = TRIM(oldncfile(INDEX(oldncfile,'/',.TRUE.)+1:))//'.arch'
          CALL assign_file_unit(archfile, arch_unit, HANDLER="fortran")
          OPEN(UNIT=arch_unit, FILE=archfile)
          CLOSE(UNIT=arch_unit)
          CALL release_file_unit(arch_unit, HANDLER="fortran")
        END IF
      END IF

      !  If a previous file existed it will be open
      !  so make sure we close it before opening the new file
      IF (nc_is_file_open(um_file)) THEN
        CALL nc_file_close(um_file)
      ELSE
        WRITE(umMessage,'(A,I4,A)')                                           &
          TRIM(RoutineName)//': Warning tried to close unopen file on unit ', &
          um_file % UNIT,' name='//TRIM(oldncfile)
        CALL umPrint(umMessage,src=RoutineName)
      END IF

      ! Construct filename.
      CALL get_filename(um_file % meta % filename_base,             &
                        um_file % filename,                         &
                        relative_offset_steps = name_offset,        &
                        reinit_steps = um_file % meta % init_steps)

      ! Open the netCDF file and write global netCDF attributes
      WRITE(umMessage,'(A,I3)')TRIM(RoutineName)//': Opening new file '//  &
          TRIM(um_file % filename)//' on unit ', um_file % UNIT
      CALL umPrint(umMessage,src=RoutineName)

      CALL nc_file_open(um_file,nc_mode_create)

      IF (um_file % meta % is_output_file) THEN
        WRITE(umMessage,'(A,I3)')                                     &
            TRIM(RoutineName)//': Initialising new file on unit ',    &
            um_file % UNIT
        CALL umPrint(umMessage,src=RoutineName)
        CALL init_nc(um_file)
      END IF

    END IF ! File needs to be initialised

    ! Increment the pointer to the next file
    um_file => um_file % next
  END DO

  ! Initialise all variables/dimensions in NetCDF files using STASH
  CALL init_stash_nc()

END IF ! mype == 0


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_reinit

END MODULE ncfile_reinit_mod

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
! Description:
!   Initialises the setup of fieldsfiles at the beginning of an NRUN or CRUN
!
! Method:
!           1. Opens all active input and output netCDF files on initial
!              call to routine.
!           2. Files to be processed on each call are controlled by
!              the file unit reinitialisation switch set on step 0 and 
!              at regular intervals thereafter.
!
! Code Description:
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

MODULE ncfile_init_mod

USE yomhook,                 ONLY: lhook,dr_hook
USE parkind1,                ONLY: jprb,jpim
USE um_parcore,              ONLY: mype,nproc
USE umnetcdf_mod,            ONLY: nc_file_open,nc_mode_create,nc_set_data_type
USE submodel_mod,            ONLY: atmos_im,internal_id_max
USE iau_mod,                 ONLY: l_iau
USE umPrintMgr,              ONLY: umPrint,umMessage
USE file_manager,            ONLY: init_file_loop,um_file_type
USE filename_generation_mod, ONLY: get_filename
USE nlstcall_mod,            ONLY: model_analysis_mins,l_fastrun   
USE model_time_mod,          ONLY: secs_per_stepim
USE init_nc_crun_mod,        ONLY: init_nc_crun
USE init_nc_mod,             ONLY: init_nc
USE cf_metadata_mod,         ONLY: get_std_name_units
USE init_stash_nc_mod,       ONLY: init_stash_nc

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NCFILE_INIT_MOD'

CONTAINS

SUBROUTINE ncfile_init

IMPLICIT NONE

!
! Local variables
!
INTEGER :: name_offset ! Character offset in file name generation
INTEGER :: ini_nc_file ! = 0 if there are netCDF files to initialise
INTEGER :: info        ! Return value from gc_ibcast

CHARACTER(LEN=*), PARAMETER :: RoutineName = "NCFILE_INIT"

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
  IF (ASSOCIATED(um_file)) THEN
     ini_nc_file = 0
  ELSE
     ini_nc_file = 1
  END IF
  DO WHILE (ASSOCIATED(um_file))

    IF (um_file % meta % partially_written) THEN

      ! A CRUN will start with active files that exist on disk. This
      ! call reads the netCDF file metadata to initialise variables needed
      ! to write further data
      CALL init_nc_crun(um_file)

    ELSE

      ! Otherwise check whether or not the file needs to be initialised at the 
      ! start of the run
      IF (um_file % meta % initialise) THEN

        IF (um_file % meta % init_steps /= 0) THEN
          ! For re-initialising output streams a suitable filename must be 
          ! generated (if not the fixed filename is already set)
          CALL get_filename(um_file % meta % filename_base,             &
                            um_file % filename,                         &
                            relative_offset_steps = name_offset,        &
                            reinit_steps = um_file % meta % init_steps)
        END IF

        ! Open the netCDF file and write global netCDF attributes
        WRITE(umMessage,'(A,I3)')TRIM(RoutineName)//': Opening new file '// &
              TRIM(um_file % filename)//' on unit ', um_file % UNIT
        CALL umPrint(umMessage,src=RoutineName)

        CALL nc_file_open(um_file,nc_mode_create)

        CALL init_nc(um_file)

      END IF ! File needs to be initialised
    END IF ! File is partially written

    ! Increment the pointer to the next file
    um_file => um_file % next
  END DO

END IF ! mype == 0

! Broadcast ini_nc_file to all PEs
CALL gc_ibcast(201,1,0,nproc,info,ini_nc_file)

IF (ini_nc_file == 0) THEN

  ! Get CF standard names and associated units
  CALL get_std_name_units()

  ! Select correct precision for dimensions and variables
  CALL nc_set_data_type()

  IF (mype == 0) THEN
    ! Initialise all variables/dimensions in NetCDF files using STASH
    CALL init_stash_nc()
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_init

END MODULE ncfile_init_mod

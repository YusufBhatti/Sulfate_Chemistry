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
!  Purpose: Gets the first and last timestep values which can be written
!           to re-initialised file um_file.


MODULE reinit_file_times_mod

USE yomhook,        ONLY: lhook,dr_hook
USE parkind1,       ONLY: jprb,jpim
USE umprintmgr,     ONLY: umprint,ummessage,printstatus,prdiag
USE file_manager,   ONLY: um_file_type
USE submodel_mod,   ONLY: atmos_im
USE nlstcall_mod,   ONLY: lcal360
USE model_time_mod, ONLY: i_year,i_month,basis_time_days,basis_time_secs,stepim
USE nlstgen_mod,    ONLY: steps_per_periodim,secs_per_periodim

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='REINIT_FILE_TIMES_MOD'

CONTAINS

SUBROUTINE reinit_file_times(um_file,filetime1,filetime2)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(OUT) :: filetime1 ! First timestep in current output file
INTEGER, INTENT(OUT) :: filetime2 ! Last timestep in current output file

!
!  Local variables
!
INTEGER :: init_steps         ! Frequency of file initialisation in timesteps
INTEGER :: init_start_step    ! Start timestep for file initialisation
INTEGER :: init_end_step      ! End timestep for file initialisation
INTEGER :: year               ! Year of filetime2 as a date
INTEGER :: month              ! Month of filetime2 as a date
INTEGER :: day                ! Day of filetime2 as a date
INTEGER :: hour               ! Hour of filetime2 as a date
INTEGER :: minute             ! Minute of filetime2 as a date
INTEGER :: second1            ! Second of filetime2 as a date
INTEGER :: elapsed_days       ! Elapsed days since basis time of filetime2
INTEGER :: elapsed_secs       ! Elapsed seconds in part of day since basis time
                              ! of filetime2

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'REINIT_FILE_TIMES'

! Dr-Hook
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Frequency, start and end of netCDF file initialisation in timesteps
init_steps = um_file % meta % init_steps
init_start_step = um_file % meta % init_start_step
init_end_step = um_file % meta % init_end_step

IF (printstatus >= prdiag) THEN
  CALL umPrint(TRIM(RoutineName)//' variables:',src=RoutineName)
  WRITE(umMessage,'(A,3(I6))') &
            'init_start_step,init_end_step,init_steps = ', &
             init_start_step,init_end_step,init_steps
  CALL umPrint(umMessage,src=RoutineName)
END IF

! First time value for current output file is current timestep
filetime1 = stepim(atmos_im)

IF (init_steps < 0) THEN ! Update frequency is in units of real months
  ! Calculate date of last valid time in output file
  ! Date will be 1st of the month, -init_steps months from current date
  ! init_step should have values -1,-3 or -12 see subroutine settsctl
  month = MOD(i_month + (-init_steps) - 1, 12) + 1
  year = i_year + (i_month + (-init_steps) - 1)/12
  day = 1
  hour = 0
  minute = 0
  second1 = 0

  ! Calculate elapsed days/seconds of last valid time in output file 
  ! relative to basis time
  ! DEPENDS ON: time2sec
  CALL time2sec(year,month,day,hour,minute,second1, &
                basis_time_days,basis_time_secs, &
                elapsed_days,elapsed_secs,lcal360)

  ! Convert elapsed time into timestep value for 
  ! last valid time in output file (filetime2)
  ! DEPENDS ON: tim2step
  CALL tim2step(elapsed_days,elapsed_secs, &
                steps_per_periodim(atmos_im),secs_per_periodim(atmos_im), &
                filetime2)
ELSE
  ! Initial value of last time value for current output file
  filetime2 = init_start_step + init_steps
  ! Find last time value for current output file
  DO WHILE (filetime2 < filetime1)
    filetime2 = filetime2 + init_steps
  END DO
END IF
! Reset filetime2 if value is greater than end time of reinitilised files
IF (init_end_step > 0) THEN
  filetime2 = MIN(init_end_step,filetime2)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE reinit_file_times

END MODULE reinit_file_times_mod

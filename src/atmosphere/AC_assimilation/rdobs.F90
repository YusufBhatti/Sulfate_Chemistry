! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Purpose : Read from ACOBS Files,reformat and place OBS header
!              details in COMOBS. The bulk of the required OBS data
!              is put into dynamic work array OBS for transmission via
!              argument list to GETOBS. OBS is written out to a cache
!              file for subsequent reading at later timesteps.
!              Thus reread of ACOBS files only required intermittently
!
!
!    Programming standard: Unified Model Documentation Paper No. 3
!
!    Project Task : P3
!
!    External documentation:
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE rdobs_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RDOBS_MOD'

CONTAINS

SUBROUTINE rdobs (nfiles,timestep_no,timestep,obs,obs_flag,       &
                  tndv,tnobs,                                     &
                  lambda_p,phi_p,p_rows,row_length,               &
                  icode,cmessage)
! ----------------------------------------------------------------------
!    INTENT IN:
!       NFILES   : NO OF AC OBSERVATION FILES TO BE READ
!                : SEE ACP NAMELIST PARAMETER NO_OBS_FILES
!       TIMESTEP_NO : TIMESTEP NO
!       TIMESTEP : TIMESTEP IN SECONDS
!    INTENT OUT:
!       TNDV     : ACTUAL SIZE OF OBS ARRAY
!       TNOBS    : ACTUAL NO OF OBS
!       OBS      : OBS array
!       OBS_FLAG : do not use flags
!       ICODE/CMESSAGE: for error processing
! ----------------------------------------------------------------------

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE comobs_mod, ONLY: nobtypmx, obs_info
USE rdobs2_mod, ONLY: rdobs2
USE ac_control_mod
USE num_obs_mod, ONLY: ac_num_obs_max, ac_tot_obs_size_max
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------
!     ARGUMENTS
INTEGER :: tndv,tnobs
INTEGER :: nfiles
INTEGER :: timestep_no
REAL :: timestep
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER :: p_rows,row_length
REAL :: lambda_p(1-halo_i:row_length+halo_i)
REAL :: phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)

!-----------------------------------------------------------------------
!     LOCAL VARIABLES
INTEGER :: j
REAL :: timerel
!-----------------------------------------------------------------------
!     Dynamic allocation
INTEGER :: obs_flag(ac_num_obs_max)
REAL :: obs(ac_tot_obs_size_max)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RDOBS'

!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ltimer_ac) CALL timer ('RDOBS   ',3)

!-----------------------------------------------------------------------
!            SECTION 1: COPY INPUT FILES TO WORK AREA, ETC.
!-----------------------------------------------------------------------

!     Find relative time for this timestep

timerel = NINT ((timestep_no-1)*timestep) / 60.0
IF (ABS(obs_info % timenext-timerel) <  0.016) timerel = NINT (timerel)

!     Determine whether AC Observation files need
!     to be read in on this timestep.

IF (timerel >=  obs_info % timenext) THEN

  !       Read in the AC Observation Files
  CALL rdobs2(nfiles,timestep_no,obs,obs_flag,tndv,               &
              tnobs,timerel,                                      &
              lambda_p,phi_p,p_rows,row_length,                   &
              icode,cmessage)
  IF (icode >  0) GO TO 999

END IF

999   CONTINUE
IF (ltimer_ac) CALL timer ('RDOBS   ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rdobs
END MODULE rdobs_mod

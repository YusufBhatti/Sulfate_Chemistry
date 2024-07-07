! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Increment the age-of-air tracer in the ukcatracers 
! structure
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: Fortran
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
MODULE ukca_age_air_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_AGE_AIR_MOD'

CONTAINS 

SUBROUTINE ukca_age_air(timestep, z_top_of_model, all_tracers)

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,               ONLY: jprb, jpim
USE ukca_d1_defs,           ONLY: Nukca_D1items, UkcaD1Codes, UKCA_sect
USE ukca_nmspec_mod,        ONLY: nm_spec
USE ukca_cspecies,          ONLY: n_age
USE ukca_option_mod,        ONLY: i_ageair_reset_method, reset_age_by_level,  &
                               max_ageair_reset_level, max_ageair_reset_height
USE nlsizes_namelist_mod, ONLY: model_levels
USE level_heights_mod,    ONLY: eta_theta_levels
USE UmPrintMgr

USE yomhook,                ONLY: lhook, dr_hook

IMPLICIT NONE

REAL, INTENT(IN)    :: timestep
REAL, INTENT(IN)    :: z_top_of_model ! Model top height
REAL, INTENT(INOUT) :: all_tracers(:,:,:,:)

! Model level upto which to reset the tracer values
INTEGER, SAVE :: max_zero_age
LOGICAL, SAVE :: first=.TRUE.
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message

! loop counters
INTEGER :: i, k
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_AGE_AIR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set n_age and maximum reset level on first call 
IF (first) THEN
  IF ( n_age == 0 )  THEN
    cmessage = ' Error setting index for Age-of-Air tracer.'
    errcode = 1
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF
  max_zero_age = -1
  IF ( i_ageair_reset_method == reset_age_by_level ) THEN
    ! Check input namelist value
    IF ( max_ageair_reset_level > model_levels ) THEN
      WRITE(cmessage,'(A,2I8)')'Age-of-air tracer reset level > model_levels ',&
            max_ageair_reset_level, model_levels
      errcode = 2
      CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
    END IF
    max_zero_age = max_ageair_reset_level   ! Up to user-specified level
  ELSE          
   ! Reset based on user-specified height
    DO k = 1, model_levels-1
      IF ( (eta_theta_levels(k) * z_top_of_model)   &
            <= max_ageair_reset_height ) max_zero_age = k
    END DO
  END IF   ! reset method
  IF ( max_zero_age < 0 )  THEN
    cmessage = ' Error setting maximum reset level for Age-of-Air tracer.'
    errcode = 3
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  ELSE
    WRITE(cmessage,'(2(A,I5))') 'UKCA AGE-OF-AIR: Reset method= ', &
       i_ageair_reset_method,'. Tracer will be reset upto level ', &
       max_zero_age
    CALL umPrint(cmessage,src=ModuleName//':'//RoutineName)
  END IF
  first =.FALSE.
END IF   ! l_first

! Increment the 'age-of-air' tracer by the length of the current
! timestep. This makes this tracer equal to the time since
! the air was last in the lowest model levels
all_tracers(:,:,1:model_levels,n_age)=                          &
    all_tracers(:,:,1:model_levels,n_age) + timestep

! set age in lower levels to zero on all timesteps
DO k = 1, max_zero_age
    all_tracers(:,:,k,n_age) = 0.0
END DO

! enforce upper boundary condition
all_tracers(:,:,model_levels,n_age) =                           &
    all_tracers(:,:,model_levels-1,n_age)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_age_air

END MODULE ukca_age_air_mod

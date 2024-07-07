! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    MODULE GETOB2_MOD
!
!    SUBROUTINE GETOB2
!
!    Purpose : Extract appropriate observations from the OBS array  
!              (passed from RDOBS via argument list)
!              GETOBS called from AC gets list of obs in time window
!              GETOB2 called from AC2 gets lat, long, time and
!                     model obs type for an obs type.
!              GETOB3 called from VERTANL gets data specific to a type
!                     (eg data values and assoc error ratios)
!
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: AC Assimilation
MODULE getob2_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GETOB2_MOD'

CONTAINS

SUBROUTINE getob2 (kact,obs,inobs,obs_lat,obs_long,obs_time,      &
                   model_obs_type,obs_no,lenobt,                  &
                   icode,cmessage)
!
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx, obs_info
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------
INTEGER :: kact                     ! IN  Obs type index
INTEGER :: inobs                    ! IN  No of obs for this type
INTEGER :: lenobt                   ! IN  No of chosen obs
INTEGER :: obs_no (lenobt)          ! IN  pointers to chosen obs
INTEGER :: model_obs_type (lenobt)  ! OUT model obs types

REAL :: obs(inobs,*)             ! IN  obs from RDOBS
REAL :: obs_lat (lenobt)         ! OUT latitudes
REAL :: obs_long(lenobt)         ! OUT longitudes
REAL :: obs_time(lenobt)         ! OUT times

INTEGER :: icode                    ! OUT error code
CHARACTER(LEN=errormessagelength) :: cmessage           ! OUT error message
!-----------------------------------------------------------------------
!     Local variables
INTEGER :: job       !  Loop counter over obs
INTEGER :: ip_lat    !  Pointer to obs latitude
INTEGER :: ip_long   !  Pointer to obs longitude
INTEGER :: ip_time   !  Pointer to obs times
INTEGER :: ip_type   !  Pointer to model obs type numbers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GETOB2'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ltimer_ac) CALL timer('GETOB2  ',3)

!     Set up pointers to data in OBS array
ip_lat  = 1   !  Latitude
ip_long = 2   !  Longitude
ip_time = 3   !  Time
ip_type = 4   !  Model observation type (**)

!     Get latitudes, longitudes, time and model obs type numbers
DO job=1,lenobt
  obs_lat(job)        = obs(obs_no(job) - obs_info % obs_no_st(kact), &
                            ip_lat)
  obs_long(job)       = obs(obs_no(job) - obs_info % obs_no_st(kact), &
                            ip_long)
  obs_time(job)       = obs(obs_no(job) - obs_info % obs_no_st(kact), &
                            ip_time)
  model_obs_type(job) = NINT(obs(obs_no(job) - obs_info % obs_no_st(kact),&
                             ip_type))
END DO


IF (ltimer_ac) CALL timer('GETOB2  ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE getob2
END MODULE getob2_mod

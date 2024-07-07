! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE DEF_TYPE ----------------------------------------------
!
!    Purpose : Initialise Observation Type Dependent Arrays
!
!    Programming Standard : UM Doc Paper No 3 ; Version 4 ; 5/2/92
!
!    Project Task : P3
!
!
!
!    ARGUMENTS:---------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE def_type_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DEF_TYPE_MOD'

CONTAINS

SUBROUTINE def_type (icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx
USE free_tracers_inputs_mod, ONLY: a_max_trvars
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

INTEGER :: icode           !  Return Code
CHARACTER(LEN=errormessagelength) :: cmessage  !  Error Message


!    WORKSPACE USAGE:---------------------------------------------------
!     Local variables
INTEGER :: i,j,jobt  !  Array index, loop counters.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEF_TYPE'


!     Arrays initialised for each observation type.
!     MASTER_AC_TYPES   - AC Observation Types known in AC Scheme
!     DEF_TIMEB         - Insertion period before observation time
!     DEF_TIMEA         - Insertion period after  observation time
!     DEF_TGETOBB       - Time Window before obs time to fetch obs
!     DEF_TGETOBA       - Time Window after  obs time to fetch obs
!     DEF_RADINF        - Maximum Normalised Influence Radius.
!     DEF_OBTHIN        - Observation thinning factor
!     DEF_CSCALE_START  - Correlation Scale at start of insertion period
!     DEF_CSCALE_OBTIME -                      observation time
!     DEF_CSCALE_END    -                      end of insertion period.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
i=0

IF (model_type == mt_global) THEN

  i = i + 1              !  Defaults for Type 101 pstar
  master_ac_types(i)   = 101
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 201 sonde temp
  master_ac_types(i)   = 201
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 202 ship surf temp
  master_ac_types(i)   = 202
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 203 airep temp
  master_ac_types(i)   = 203
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 204 synop surf temp
  master_ac_types(i)   = 204
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 200.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0

  i = i + 1              !  Defaults for Type 205 LASS temp
  master_ac_types(i)   = 205
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 3
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 206 SATEM temp
  master_ac_types(i)   = 206
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 207 SAT120 temp
  master_ac_types(i)   = 207
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 3
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 208 constrained LASS
  master_ac_types(i)   = 208
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 3
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 209 BOGUS 1000-500
  master_ac_types(i)   = 209
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 211 UARS Temp
  master_ac_types(i)   = 211
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 301 Sonde winds
  master_ac_types(i)   = 301
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 302 Ship surf wind
  master_ac_types(i)   = 302
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 303 Airep wind
  master_ac_types(i)   = 303
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 304 Synop surf wind
  master_ac_types(i)   = 304
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 240.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0

  i = i + 1              !  Defaults for Type 305 Scatwind
  master_ac_types(i)   = 305
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 5
  def_cscale_start(i)  = 240.0
  def_cscale_obtime(i) = 120.0
  def_cscale_end(i)    = 120.0

  i = i + 1              !  Defaults for Type 306 Drifter winds
  master_ac_types(i)   = 306
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 311 UARS wind
  master_ac_types(i)   = 311
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 200.0
  def_cscale_end(i)    = 200.0

  i = i + 1              !  Defaults for Type 401 Sonde RH
  master_ac_types(i)   = 401
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0

  i = i + 1              !  Defaults for Type 402 Ship surf RH
  master_ac_types(i)   = 402
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0

  i = i + 1              !  Defaults for Type 403 Bogus RH
  master_ac_types(i)   = 403
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0

  i = i + 1              !  Defaults for Type 404 Synop surf RH
  master_ac_types(i)   = 404
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 200.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0

  i = i + 1              !  Defaults for Type 405 LASS RH
  master_ac_types(i)   = 405
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 3
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0

  i = i + 1              !  Defaults for Type 406 MOPS RH
  master_ac_types(i)   = 406
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 0.01
  def_obthin(i)        = 1
  def_cscale_start(i)  = 999.0
  def_cscale_obtime(i) = 999.0
  def_cscale_end(i)    = 999.0

  i = i + 1              !  Defaults for Type 407 Cloud histograms
  master_ac_types(i)   = 407
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 100.0
  def_cscale_obtime(i) = 100.0
  def_cscale_end(i)    = 100.0

  i = i + 1              !  Defaults for type 506 MOPS precip
  master_ac_types(i)   = 506
  def_timeb(i)         = 180.0
  def_timea(i)         = 180.0
  def_tgetobb(i)       = 180.0
  def_tgetoba(i)       = 180.0
  def_radinf(i)        = 0.01
  def_obthin(i)        = 1
  def_cscale_start(i)  = 999.0
  def_cscale_obtime(i) = 999.0
  def_cscale_end(i)    = 999.0

  !  Defaults for Type 6nn (tracers)
  !  Allow for A_MAX_TRVARS (currently 29) possible tracers
  DO j=1,a_max_trvars
    i = i + 1
    master_ac_types(i) = 600+j
    def_timeb(i)       = 240.0
    def_timea(i)       = 60.0
    def_tgetobb(i)     = 240.0
    def_tgetoba(i)     = 60.0
    def_radinf(i)      = 3.5
    def_obthin(i)      = 1
    def_cscale_start(i) = 360.0
    def_cscale_obtime(i) = 200.0
    def_cscale_end(i)  = 200.0
  END DO

  i = i + 1              !  Defaults for Type 901 LOG Visibility
  master_ac_types(i)   = 901
  def_timeb(i)         = 240.0
  def_timea(i)         = 60.0
  def_tgetobb(i)       = 240.0
  def_tgetoba(i)       = 60.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 360.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0

  !  Defaults for Limited Area.
ELSE

  i = i + 1              !  Defaults for Type 101 pstar
  master_ac_types(i)   = 101
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 225.0
    def_cscale_obtime(i) = 150.0
    def_cscale_end(i)    = 165.0
    def_timeb(i)         = 120.0
    def_timea(i)         = 24.0
    def_tgetobb(i)       = 120.0
    def_tgetoba(i)       = 24.0
    def_radinf(i)        = 1.75
  END IF

  i = i + 1              !  Defaults for Type 201 sonde temp
  master_ac_types(i)   = 201
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 120.0
    def_cscale_obtime(i) = 100.0
    def_cscale_end(i)    = 100.0
  END IF

  i = i + 1              !  Defaults for Type 202 ship surf temp
  master_ac_types(i)   = 202
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 225.0
    def_cscale_obtime(i) = 150.0
    def_cscale_end(i)    = 165.0
  END IF

  i = i + 1              !  Defaults for Type 203 airep temp
  master_ac_types(i)   = 203
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 120.0
    def_cscale_obtime(i) = 100.0
    def_cscale_end(i)    = 100.0
  END IF

  i = i + 1              !  Defaults for Type 204 synop surf temp
  master_ac_types(i)   = 204
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 150.0
  def_cscale_obtime(i) = 120.0
  def_cscale_end(i)    = 120.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 105.0
    def_cscale_obtime(i) = 85.0
    def_cscale_end(i)    = 85.0
    def_timeb(i)         = 120.0
    def_timea(i)         = 24.0
    def_tgetobb(i)       = 120.0
    def_tgetoba(i)       = 24.0
    def_radinf(i)        = 1.75
  END IF

  i = i + 1              !  Defaults for Type 205 LASS temp
  master_ac_types(i)   = 205
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 120.0
    def_cscale_obtime(i) = 100.0
    def_cscale_end(i)    = 100.0
  END IF

  i = i + 1              !  Defaults for Type 206 SATEM temp
  master_ac_types(i)   = 206
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 120.0
    def_cscale_obtime(i) = 100.0
    def_cscale_end(i)    = 100.0
  END IF

  i = i + 1              !  Defaults for Type 207 SAT120 temp
  master_ac_types(i)   = 207
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 120.0
    def_cscale_obtime(i) = 100.0
    def_cscale_end(i)    = 100.0
  END IF

  i = i + 1              !  Defaults for Type 208 constrained LASS
  master_ac_types(i)   = 208
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0

  i = i + 1              !  Defaults for Type 209 BOGUS 1000-500
  master_ac_types(i)   = 209
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 120.0
    def_cscale_obtime(i) = 100.0
    def_cscale_end(i)    = 100.0
  END IF

  i = i + 1              !  Defaults for Type 211 UARS Temp
  master_ac_types(i)   = 211
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0

  i = i + 1              !  Defaults for Type 301 Sonde winds
  master_ac_types(i)   = 301
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 120.0
    def_cscale_obtime(i) = 100.0
    def_cscale_end(i)    = 100.0
  END IF

  i = i + 1              !  Defaults for Type 302 Ship surf wind
  master_ac_types(i)   = 302
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  =  95.0
    def_cscale_obtime(i) =  75.0
    def_cscale_end(i)    =  75.0
  END IF

  i = i + 1              !  Defaults for Type 303 Airep wind
  master_ac_types(i)   = 303
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 120.0
    def_cscale_obtime(i) = 100.0
    def_cscale_end(i)    = 100.0
  END IF

  i = i + 1              !  Defaults for Type 304 Synop surf wind
  master_ac_types(i)   = 304
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 240.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0
  IF (lac_mes) THEN
    def_cscale_start(i)  =  95.0
    def_cscale_obtime(i) =  75.0
    def_cscale_end(i)    =  75.0
    def_timeb(i)         = 120.0
    def_timea(i)         = 24.0
    def_tgetobb(i)       = 120.0
    def_tgetoba(i)       = 24.0
    def_radinf(i)        = 1.75
  END IF

  i = i + 1              !  Defaults for Type 305 Scatwind
  master_ac_types(i)   = 305
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 240.0
  def_cscale_obtime(i) = 120.0
  def_cscale_end(i)    = 120.0

  i = i + 1              !  Defaults for Type 306 Ship surf wind
  master_ac_types(i)   = 306
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0
  IF (lac_mes) THEN
    def_cscale_start(i)  =  95.0
    def_cscale_obtime(i) =  75.0
    def_cscale_end(i)    =  75.0
  END IF

  i = i + 1              !  Defaults for Type 311 UARS wind
  master_ac_types(i)   = 311
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 3.5
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 180.0
  def_cscale_end(i)    = 180.0

  i = i + 1              !  Defaults for Type 401 Sonde RH
  master_ac_types(i)   = 401
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 115.0
    def_cscale_obtime(i) = 85.0
    def_cscale_end(i)    = 85.0
  END IF

  i = i + 1              !  Defaults for Type 402 Ship surf RH
  master_ac_types(i)   = 402
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 240.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 225.0
    def_cscale_obtime(i) = 150.0
    def_cscale_end(i)    = 165.0
  END IF

  i = i + 1              !  Defaults for Type 403 Bogus RH
  master_ac_types(i)   = 403
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 115.0
    def_cscale_obtime(i) =  85.0
    def_cscale_end(i)    =  85.0
  END IF

  i = i + 1              !  Defaults for Type 404 Synop surf RH
  master_ac_types(i)   = 404
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 150.0
  def_cscale_obtime(i) = 120.0
  def_cscale_end(i)    = 120.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 100.0
    def_cscale_obtime(i) =  75.0
    def_cscale_end(i)    =  75.0
    def_timeb(i)         = 120.0
    def_timea(i)         = 24.0
    def_tgetobb(i)       = 120.0
    def_tgetoba(i)       = 24.0
    def_radinf(i)        = 1.75
  END IF

  i = i + 1              !  Defaults for Type 405 LASS RH
  master_ac_types(i)   = 405
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0

  i = i + 1              !  Defaults for Type 406 MOPS RH
  master_ac_types(i)   = 406
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 0.01
  def_obthin(i)        = 1
  def_cscale_start(i)  = 999.0
  def_cscale_obtime(i) = 999.0
  def_cscale_end(i)    = 999.0

  i = i + 1              !  Defaults for Type 407 Cloud histograms
  master_ac_types(i)   = 407
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 50.0
  def_cscale_obtime(i) = 50.0
  def_cscale_end(i)    = 50.0

  i = i + 1              !  Defaults for type 506 MOPS precip
  master_ac_types(i)   = 506
  def_timeb(i)         = 180.0
  def_timea(i)         = 180.0
  def_tgetobb(i)       = 180.0
  def_tgetoba(i)       = 180.0
  def_radinf(i)        = 0.01
  def_obthin(i)        = 1
  def_cscale_start(i)  = 999.0
  def_cscale_obtime(i) = 999.0
  def_cscale_end(i)    = 999.0

  i = i + 1              !  Defaults for Type 901 LOG Visibility
  master_ac_types(i)   = 901
  def_timeb(i)         = 150.0
  def_timea(i)         = 30.0
  def_tgetobb(i)       = 150.0
  def_tgetoba(i)       = 30.0
  def_radinf(i)        = 2.25
  def_obthin(i)        = 1
  def_cscale_start(i)  = 300.0
  def_cscale_obtime(i) = 150.0
  def_cscale_end(i)    = 150.0
  IF (lac_mes) THEN
    def_cscale_start(i)  = 100.0
    def_cscale_obtime(i) =  75.0
    def_cscale_end(i)    =  75.0
    def_timeb(i)         = 120.0
    def_timea(i)         = 24.0
    def_tgetobb(i)       = 120.0
    def_tgetoba(i)       = 24.0
    def_radinf(i)        = 1.75
  END IF

  !  Defaults for Type 6nn (tracers)
  !  Allow for A_MAX_TRVARS (currently 29) possible tracers
  DO j=1,a_max_trvars
    i = i + 1
    master_ac_types(i) = 600+j
    def_timeb(i)       = 150.0
    def_timea(i)       = 30.0
    def_tgetobb(i)     = 150.0
    def_tgetoba(i)     = 30.0
    def_radinf(i)      = 3.5
    def_obthin(i)      = 1
    def_cscale_start(i) = 300.0
    def_cscale_obtime(i) = 180.0
    def_cscale_end(i)  = 180.0
  END DO

END IF  !  if GLOBAL

!     Modify defaults for UARS assimilation.
IF (lac_uars) THEN
  DO jobt=1,i
    def_cscale_start(jobt)  = 600.0
    def_cscale_obtime(jobt) = 400.0
    def_cscale_end(jobt)    = 400.0
    IF (master_ac_types(jobt) == 207) THEN
      def_obthin(jobt)      = 2
    ELSE
      def_obthin(jobt)      = 1
    END IF
  END DO
END IF

!     Initialise rest of arrays.
IF (i <  nobtypmx) THEN
  DO jobt = i+1,nobtypmx
    master_ac_types(jobt)   = 0
    def_timeb(jobt)         = 0.0
    def_timea(jobt)         = 0.0
    def_tgetobb(jobt)       = 0.0
    def_tgetoba(jobt)       = 0.0
    def_radinf(jobt)        = 0.0
    def_obthin(jobt)        = 0
    def_cscale_start(jobt)  = 0.0
    def_cscale_obtime(jobt) = 0.0
    def_cscale_end(jobt)    = 0.0
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE def_type
END MODULE def_type_mod

!
! MODULE SCM_UTILS--------------------------------------------------------------
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! ******************************COPYRIGHT***************************************


MODULE scm_utils

!-------------------------------------------------------------------------------
! Description:
!   SCM UTILITY MODULE; For use with SCM in routines with "USE scm_utils"
!
!   Contains subroutines:
!   scm_message: To send messages to stdout or unit number with timestep
!                information. Calls with no arguments will output the current
!                model timestep.
!
!   scm_trap_nan:
!   alloc_common_scm:
!   dealloc_common_scm:
!-------------------------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
USE um_types
USE missing_data_mod
USE jules_surface_types_mod, ONLY: ntype, npft, nnvg, soil
USE s_maxdim
USE scm_cntl_mod, ONLY: lsname
USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
IMPLICIT NONE

! Missing data indicators

  ! Precision of integers passed in/out out netcdf routines
INTEGER, PARAMETER, PRIVATE :: incdf = integer32
INTEGER           , PRIVATE :: nan_list_count = 0 ! Updated by scm trap nan
INTEGER, PARAMETER, PRIVATE :: max_nans = 50      ! Kill the run if the number
                                                  ! of diagnostics containing
                                                  ! NaNs reaches this value.
CHARACTER(LEN=lsname) , PRIVATE :: nan_list(max_nans) ! Updated by scm trap nan

INTEGER :: time_info(7)

! Integer id for time and date variables in netcdf files
INTEGER(incdf) :: scm_nc_time_id
INTEGER(incdf) :: scm_nc_date_id
INTEGER(incdf) :: scm_nc_local_id
INTEGER(incdf) :: scm_nc_hour_id
INTEGER        :: scm_timestep_count ! Number of current timestep
INTEGER        :: daycount
REAL           :: scm_timestep       ! SCM run timestep

INTEGER ::    &
  rw          &! Rows
, rw_lng      &! Row length
, nbl_lv      &! No. of boundary layer levels
, o3_lv       &! No. of ozone levels
, nlnd        &! No. of land points
, ntile       &! No. of tiles per land point
, ntyp        &! No. of plant profiles
, st_lv       &! No. of soil temperature levels
, sm_lv       &! No. of soil moisture levels
, ntr_lv      &! No. of tracer levels
, ntr_var      ! No. of tracers

INTEGER ::    &
  nml_nmod_lv &! No. of model levels specified in obs forcing
, nobs        &! No. of observational forcing profiles in namelist
, obs_t0       ! Observational profile at which to begin run

REAL, ALLOCATABLE :: &
  z_th(:)     &! Height of theta levels
, z_rh(:)     &! Height of rho   levels
, nml_z_th(:) &
, nml_z_rh(:)

LOGICAL ::              &
  cv_run_mess  = .FALSE. &! Display runtime convective messages
, nc_obsfor    = .FALSE. &!
, old_nml      = .FALSE. &! Old namelist format (pre vn7.7)
, old_rlx      = .FALSE. &! Old forcing relaxation code (pre vn7.7)
, old_vertadv  = .FALSE. &
              ! Old interactive vertical advection code (pre vn7.7)
, interpolate  = .FALSE. &! Auto-interpolation (unsupported)
, l_emcorr_opt = .FALSE.  ! Energy correction logical

! Parameters for forcing relaxation case options
INTEGER, PARAMETER ::   &
  rlx_none      = 0     &! No relaxation
, rlx_init      = 1     &! Relax to initial profile across tau
, rlx_bgrd      = 2     &! Relax to background profile across tau
, rlx_inst_init = 3     &! Reset to initial profile
, rlx_inst_bgrd = 4      ! Reset to background profile


! Note: the number of levels in this super-array follows its use
! for saving fields in subroutine restart_dump (called from scm_main),
! and loading fields in subroutine dumpinit (called from run_init).
INTEGER, PARAMETER :: nprimvars = 9*mx_mod_lv  & ! u,v,w,t,theta,q,qcl,qcf,cf
                                + mx_st_lv     & ! t_deep_soil
                                + 3*mx_mod_lv  & ! p, rho, cca
                                + mx_sm_lv     & ! smcl
                                + 4*mx_mod_lv  & ! cca_dp,cca_md,cca_sh,bl_w_var
                                + 8              ! smc,canopy_gb,snodep,tstar,
                                                 ! zh,z0msea,rccb,rcct


CHARACTER(LEN=8) :: time_string ! Actual time at end of timestep

! Array for holding dumps/primary variables/convection sub-step data
REAL, ALLOCATABLE:: resdump(:,:,:)


! Dr Hook Parameters
!=================================
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SCM_UTILS'

CONTAINS

SUBROUTINE scm_message(string, UNIT)

IMPLICIT NONE

CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: string
INTEGER,      OPTIONAL, INTENT(IN) :: UNIT

CHARACTER(LEN=10) :: ts_str
INTEGER       :: unitno

! Dr Hook
!==============================
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SCM_MESSAGE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(UNIT)) THEN
  unitno = UNIT
ELSE
  unitno = 6
END IF

WRITE(ts_str,'(A4,I4,A2)') '[TS ', scm_timestep_count,'] '

IF (PRESENT(string)) THEN
  WRITE(unitno,*) ts_str, TRIM(ADJUSTL(string))
ELSE
  WRITE(unitno,*) ts_str
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE scm_message

!-----------------------------------------------------------------------------

SUBROUTINE scm_trap_nan(string, routine)

USE umPrintMgr, ONLY: umPrint, umMessage, newline
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: string  ! Diagnostic with NaN
CHARACTER(LEN=*), INTENT(IN) :: routine ! Calling routine
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: icode
INTEGER :: i

! Dr Hook
!==============================
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SCM_TRAP_NAN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (nan_list_count <= max_nans) THEN

  ! Check whether we've already detected NaNs in this diagnostic;
  ! Do nothing if so
  DO i=1, nan_list_count
    IF (nan_list(i) == string) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, &
                              zhook_handle)
      RETURN
    END IF
  END DO

  ! Raise a warning if we've found NaNs in a diagnostic that we haven't
  ! found them in before
  icode = -777
  WRITE(cmessage,"(A,I4,A)") '[TS ', scm_timestep_count,'] : ' //  &
                             'NaN detected: '//string//' from '//routine
  CALL ereport(routinename, icode, cmessage)

  ! Store the latest NaN-diag's name and increment the counter.
  nan_list(nan_list_count+1) = string
  nan_list_count = nan_list_count + 1

  IF (nan_list_count == max_nans) THEN
    ! Raise a fatal error if detected NaNa in max_nans diagnostics.
    icode = 777
    WRITE(cmessage,"(A,I4,A,I3,A)") '[TS ', scm_timestep_count,'] : ',  &
                     max_nans, '+ Diagnostics contain NaNs.  Stopping run.'
    CALL ereport(routinename, icode, cmessage)
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE scm_trap_nan

!=============================================================================

SUBROUTINE alloc_common_scm

USE nlsizes_namelist_mod, ONLY: model_levels

IMPLICIT NONE

! Dr Hook
!==============================
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOC_COMMON_SCM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE                    &
  ( z_th(model_levels)      &
  , z_rh(model_levels+1)    &
  , nml_z_th(nml_nmod_lv)   &
  , nml_z_rh(nml_nmod_lv+1) )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE alloc_common_scm

!-------------------------------------------------------------------------------

SUBROUTINE dealloc_common_scm

IMPLICIT NONE

! Dr Hook
!=============================================
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEALLOC_COMMON_SCM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DEALLOCATE   &
  ( z_th     &
  , z_rh     &
  , nml_z_th &
  , nml_z_rh )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE dealloc_common_scm


!-------------------------------------------------------------------------------
END MODULE scm_utils

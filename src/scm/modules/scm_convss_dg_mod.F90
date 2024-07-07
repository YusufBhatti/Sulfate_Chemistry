! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module containing a type definition for SCM diagnostics from convection,
!  and subroutines to initialise, output and deallocate the diagnostic arrays.

MODULE scm_convss_dg_mod


IMPLICIT NONE

!
! Description:
!   Defines a derived type structure containing all the SCM diagnostics that
!   need to be passed up from inside the convection scheme (e.g. diagnostics
!   based on variables which are only available inside deep_conv_mod).
!   SCMoutput can't be called from inside the convection scheme, because
!   it is only called on some timesteps at some points, using compressed
!   arrays.  SCMoutput assumes the data passed to it are of dimension
!   (row_length, rows).  Therefore, the type defined here is used to declare 
!   a structure, containing numerous diagnostics, with the full dimensions
!   in ni_conv_ctl, which is passed down into compressed dimensions in
!   glueconv, deep_conv_mod, etc.
!   This module also contains subroutines to allocate and initialise the
!   diagnostics to zero, and output them via scmoutput.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName =                    &
                                        'SCM_CONVSS_DG_MOD'



!------------------------------------------------------------------------
! Type definition for SCM convection diagnostics
!------------------------------------------------------------------------
! NOTE: If adding more diagnostics to this list, they also need to be
! added to the equivalent type definition in the stub routine
! src/scm/stub/scm_convss_dg_mod_stub.F90.
! Otherwise, where the new components are referenced in the UM code,
! the full UM will fail at compile time due to the stated fields not
! actually being components of this type.
TYPE scm_convss_dg_type

  ! Note that making all the fields below allocatable is unnecessary,
  ! given that the first routine that uses this type already knows that
  ! all the arrays should have the same size (model_levels).
  ! This structure is itself declared as an allocatable array; no
  ! memory is used until that is allocated, so there is no advantage
  ! to using allocatable fields in here.
  ! This would be an ideal time to use the fortran 2003 parameterised
  ! derived type feature (i.e. it allows all the arrays to be dimensioned
  ! with model_levels from the outset, by making model_levels an argument
  ! to this type).  We then wouldn't need the allocate / deallocate
  ! subroutines below.  However, UMDP3 states that f2003 parameterised
  ! derived types are not allowed, so we'll just have to do it this way
  ! for now.

  ! Environment profile seen by the convective main ascent
  REAL, ALLOCATABLE :: env_theta(:)
  REAL, ALLOCATABLE :: env_q(:)
  REAL, ALLOCATABLE :: env_qcl(:)
  REAL, ALLOCATABLE :: env_qcf(:)
  REAL, ALLOCATABLE :: env_thetav(:)

  ! Parcel properties from the main ascent
  REAL, ALLOCATABLE :: par_theta(:)
  REAL, ALLOCATABLE :: par_q(:)
  REAL, ALLOCATABLE :: par_qcl(:)
  REAL, ALLOCATABLE :: par_qcf(:)
  REAL, ALLOCATABLE :: par_thetav(:)

  ! Parcel buoyancy excess
  REAL, ALLOCATABLE :: par_thetav_excess(:)

  ! Sub-step mass-flux profile before and after closure
  REAL, ALLOCATABLE :: up_flx_guess(:)
  REAL, ALLOCATABLE :: up_flx(:)

  ! Fractional entrainment and mixing detrainment rates in layer
  REAL, ALLOCATABLE :: ekp14(:)
  REAL, ALLOCATABLE :: ekp34(:)
  REAL, ALLOCATABLE :: amdetk(:)

  ! Diagnostics for adaptive detrainment
  REAL, ALLOCATABLE :: deltak(:)
  REAL, ALLOCATABLE :: rbuoy_star(:)
  REAL, ALLOCATABLE :: xsbmin(:)
  REAL, ALLOCATABLE :: thrk(:)
  REAL, ALLOCATABLE :: qrk(:)
  REAL, ALLOCATABLE :: thvrk_excess(:)

  ! Status of deep, shallow and mid-level convection
  ! (0 = not diagnosed, 1 = failed ascent, 2 = zero closure, 3 = real convection
  INTEGER :: status_deep
  INTEGER :: status_shallow
  INTEGER :: status_mid

  ! Surface precip from each convection type
  REAL :: precip_deep
  REAL :: precip_shallow
  REAL :: precip_mid

END TYPE scm_convss_dg_type


CONTAINS



!------------------------------------------------------------------------
! Subroutine to allocate the SCM diagnostic fields within the above structure
!------------------------------------------------------------------------
SUBROUTINE scm_convss_dg_allocate( scm_convss_dg )

USE nlsizes_namelist_mod,     ONLY: model_levels

USE yomhook,                  ONLY: lhook, dr_hook
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

TYPE(scm_convss_dg_type), INTENT(INOUT) :: scm_convss_dg

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SCM_CONVSS_DG_ALLOCATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Environment profile seen by the convective main ascent
ALLOCATE( scm_convss_dg % env_theta(model_levels) )
ALLOCATE( scm_convss_dg % env_q(model_levels) )
ALLOCATE( scm_convss_dg % env_qcl(model_levels) )
ALLOCATE( scm_convss_dg % env_qcf(model_levels) )
ALLOCATE( scm_convss_dg % env_thetav(model_levels) )

! Parcel properties from the main ascent
ALLOCATE( scm_convss_dg % par_theta(model_levels) )
ALLOCATE( scm_convss_dg % par_q(model_levels) )
ALLOCATE( scm_convss_dg % par_qcl(model_levels) )
ALLOCATE( scm_convss_dg % par_qcf(model_levels) )
ALLOCATE( scm_convss_dg % par_thetav(model_levels) )

! Parcel buoyancy excess
ALLOCATE( scm_convss_dg % par_thetav_excess(model_levels) )

! Sub-step mass-flux profile before and after closure
ALLOCATE( scm_convss_dg % up_flx_guess(model_levels) )
ALLOCATE( scm_convss_dg % up_flx(model_levels) )

! Fractional entrainment and mixing detrainment rates in layer
ALLOCATE( scm_convss_dg % ekp14(model_levels) )
ALLOCATE( scm_convss_dg % ekp34(model_levels) )
ALLOCATE( scm_convss_dg % amdetk(model_levels) )

! Diagnostics for adaptive detrainment
ALLOCATE( scm_convss_dg % deltak(model_levels) )
ALLOCATE( scm_convss_dg % rbuoy_star(model_levels) )
ALLOCATE( scm_convss_dg % xsbmin(model_levels) )
ALLOCATE( scm_convss_dg % thrk(model_levels) )
ALLOCATE( scm_convss_dg % qrk(model_levels) )
ALLOCATE( scm_convss_dg % thvrk_excess(model_levels) )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE scm_convss_dg_allocate



!------------------------------------------------------------------------
! Subroutine to deallocate the above SCM diagnostics
!------------------------------------------------------------------------
SUBROUTINE scm_convss_dg_deallocate( scm_convss_dg )

USE yomhook,                  ONLY: lhook, dr_hook
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

TYPE(scm_convss_dg_type), INTENT(INOUT) :: scm_convss_dg

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SCM_CONVSS_DG_DEALLOCATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Environment profile seen by the convective main ascent
DEALLOCATE( scm_convss_dg % env_theta )
DEALLOCATE( scm_convss_dg % env_q )
DEALLOCATE( scm_convss_dg % env_qcl )
DEALLOCATE( scm_convss_dg % env_qcf )
DEALLOCATE( scm_convss_dg % env_thetav )

! Parcel properties from the main ascent
DEALLOCATE( scm_convss_dg % par_theta )
DEALLOCATE( scm_convss_dg % par_q )
DEALLOCATE( scm_convss_dg % par_qcl )
DEALLOCATE( scm_convss_dg % par_qcf )
DEALLOCATE( scm_convss_dg % par_thetav )

! Parcel buoyancy excess
DEALLOCATE( scm_convss_dg % par_thetav_excess )

! Sub-step mass-flux profile before and after closure
DEALLOCATE( scm_convss_dg % up_flx_guess )
DEALLOCATE( scm_convss_dg % up_flx )

! Fractional entrainment and mixing detrainment rates in layer
DEALLOCATE( scm_convss_dg % ekp14 )
DEALLOCATE( scm_convss_dg % ekp34 )
DEALLOCATE( scm_convss_dg % amdetk )

! Diagnostics for adaptive detrainment
DEALLOCATE( scm_convss_dg % deltak )
DEALLOCATE( scm_convss_dg % rbuoy_star )
DEALLOCATE( scm_convss_dg % xsbmin )
DEALLOCATE( scm_convss_dg % thrk )
DEALLOCATE( scm_convss_dg % qrk )
DEALLOCATE( scm_convss_dg % thvrk_excess )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE scm_convss_dg_deallocate



!------------------------------------------------------------------------
! Subroutine to initialise the above SCM convection diagnostics to zero
!------------------------------------------------------------------------
SUBROUTINE scm_convss_dg_initzero( scm_convss_dg )

USE nlsizes_namelist_mod,     ONLY: row_length, rows, model_levels

USE yomhook,                  ONLY: lhook, dr_hook
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

TYPE(scm_convss_dg_type), INTENT(INOUT) :: scm_convss_dg(row_length,rows)

INTEGER :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SCM_CONVSS_DG_INITZERO'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise all the diagnostics to zero
DO j = 1, rows
  DO i = 1, row_length

    ! 3D diagnostics
    DO k = 1, model_levels

      ! Environment profile seen by the convective main ascent
      scm_convss_dg(i,j) % env_theta(k) = 0.0
      scm_convss_dg(i,j) % env_q(k) = 0.0
      scm_convss_dg(i,j) % env_qcl(k) = 0.0
      scm_convss_dg(i,j) % env_qcf(k) = 0.0
      scm_convss_dg(i,j) % env_thetav(k) = 0.0

      ! Parcel properties from the main ascent
      scm_convss_dg(i,j) % par_theta(k) = 0.0
      scm_convss_dg(i,j) % par_q(k) = 0.0
      scm_convss_dg(i,j) % par_qcl(k) = 0.0
      scm_convss_dg(i,j) % par_qcf(k) = 0.0
      scm_convss_dg(i,j) % par_thetav(k) = 0.0

      ! Parcel buoyancy excess
      scm_convss_dg(i,j) % par_thetav_excess(k) = 0.0

      ! Sub-step mass-flux profile before and after closure
      scm_convss_dg(i,j) % up_flx_guess(k) = 0.0
      scm_convss_dg(i,j) % up_flx(k) = 0.0

      ! Fractional entrainment and mixing detrainment rates in layer
      scm_convss_dg(i,j) % ekp14(k) = 0.0
      scm_convss_dg(i,j) % ekp34(k) = 0.0
      scm_convss_dg(i,j) % amdetk(k) = 0.0

      ! Diagnostics for adaptive detrainment
      scm_convss_dg(i,j) % deltak(k) = 0.0
      scm_convss_dg(i,j) % rbuoy_star(k) = 0.0
      scm_convss_dg(i,j) % xsbmin(k) = 0.0
      scm_convss_dg(i,j) % thrk(k) = 0.0
      scm_convss_dg(i,j) % qrk(k) = 0.0
      scm_convss_dg(i,j) % thvrk_excess(k) = 0.0

    END DO

    ! 2D diagnostics

    ! Status of deep, shallow and mid-level convection
    ! (0 = not diagnosed, 1 = failed ascent, 2 = zero closure, 3 = real conv
    scm_convss_dg(i,j) % status_deep = 0
    scm_convss_dg(i,j) % status_shallow = 0
    scm_convss_dg(i,j) % status_mid = 0

    ! Surface precip from each convection type
    scm_convss_dg(i,j) % precip_deep = 0.0
    scm_convss_dg(i,j) % precip_shallow = 0.0
    scm_convss_dg(i,j) % precip_mid = 0.0

  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE scm_convss_dg_initzero



!------------------------------------------------------------------------
! Subroutine to output the above diagnostics using scmoutput.
! This routine is called from ni_conv_ctl
!------------------------------------------------------------------------
SUBROUTINE scm_convss_dg_output( scm_convss_dg, call_number )

USE nlsizes_namelist_mod,     ONLY: row_length, rows, model_levels
USE s_scmop_mod,              ONLY: default_streams, t_avg, d_all,  &
                                    t_inst, d_point

USE scmoutput_mod,            ONLY: scmoutput

USE yomhook,                  ONLY: lhook, dr_hook
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE


! Structure declared in ni_conv_ctl, containing all the SCM diagnostics
! on the i,j grid (not compressed)
TYPE(scm_convss_dg_type), INTENT(INOUT) :: scm_convss_dg(row_length,rows)

! Current convective substep number
INTEGER, INTENT(IN) :: call_number

! Array with correct dimensions, for copying diagnostics into for output
REAL :: tmpscm3d(row_length,rows,model_levels)
REAL :: tmpscm2d(row_length,rows)

! call_number stored in a character string
CHARACTER (LEN=8) :: substepchar
CHARACTER (LEN=12) :: substepchar_long

INTEGER :: i,j,k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SCM_CONVSS_DG_OUTPUT'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Copy the current substep number into a character string, to be appended
! to all the diagnostic names.  This will make the diagnostics from
! subsequent substeps appear to be distinct in scmoutput, preventing
! the error with duplicated diagnostics in the same timestep.
WRITE( substepchar, "(a,i1,a)" ) "convss", call_number, "_"

! Long version for adding to the diagnostic description
WRITE( substepchar_long, "(a,i1,a)" ) ' (substep ', call_number, ')'


! Copy diagnostics into tmpscm3d, and output...



! Environment profile seen by the convective main ascent

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % env_theta(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'env_theta',                     &
     'Convection environment potential temperature' // substepchar_long,  &
     'K', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % env_q(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'env_q',                         &
     'Convection environment specific humidity' // substepchar_long,      &
     'kg kg-1', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % env_qcl(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'env_qcl',                       &
     'Convection environment cloud liquid water' // substepchar_long,     &
     'kg kg-1', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % env_qcf(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'env_qcf',                       &
     'Convection environment cloud ice' // substepchar_long,              &
     'kg kg-1', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % env_thetav(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'env_thetav',                    &
     'Convection environment virt. pot. temp.' // substepchar_long,       &
     'K', t_avg, d_all, default_streams, '', RoutineName )



  ! Parcel properties from the main ascent

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % par_theta(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'par_theta',                     &
     'Convection parcel potential temperature' // substepchar_long,       &
     'K', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % par_q(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'par_q',                         &
     'Convection parcel specific humidity' // substepchar_long,           &
     'kg kg-1', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % par_qcl(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'par_qcl',                       &
     'Convection parcel cloud liquid water' // substepchar_long,          &
     'kg kg-1', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % par_qcf(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'par_qcf',                       &
     'Convection parcel cloud ice' // substepchar_long,                   &
     'kg kg-1', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % par_thetav(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'par_thetav',                    &
     'Convection parcel virt. pot. temp.' // substepchar_long,            &
     'K', t_avg, d_all, default_streams, '', RoutineName )



! Parcel buoyancy excess

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % par_thetav_excess(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'par_thetav_excess',             &
     'Convection parcel buoyancy w.r.t. environment' // substepchar_long, &
     'K', t_avg, d_all, default_streams, '', RoutineName )



! Sub-step mass-flux profile before and after closure

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % up_flx_guess(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'up_flx_guess',                  &
     'Convection updraft mass-flux before closure' // substepchar_long,   &
     'Pa s-1', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % up_flx(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'up_flx',                        &
     'Convection updraft mass-flux after closure' // substepchar_long,    &
     'Pa s-1', t_avg, d_all, default_streams, '', RoutineName )



! Fractional entrainment and mixing detrainment rates in layer

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % ekp14(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'ekp14',                         &
     'Fractional entrainment between k and k+1/2' // substepchar_long,    &
     'dimensionless', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % ekp34(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'ekp34',                         &
     'Fractional entrainment between k+1/2 and k+1' // substepchar_long,  &
     'dimensionless', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % amdetk(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'amdetk',                        &
     'Fractional mixing detrainment between k and k+1' //substepchar_long,&
     'dimensionless', t_avg, d_all, default_streams, '', RoutineName )



! Diagnostics for adaptive detrainment

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % deltak(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'deltak',                        &
     'Fractional forced detrainment between k and k+1' //substepchar_long,&
     'dimensionless', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % rbuoy_star(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'rbuoy_star',                    &
     'Parcel buoyancy before forced detrainment' // substepchar_long,     &
     'K', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % xsbmin(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'xsbmin',                        &
     'Target min buoyancy after forced detrainment' // substepchar_long,  &
     'K', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % thrk(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'thrk',                          &
     'Potential temperature of forced detrained air' // substepchar_long, &
     'K', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % qrk(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'qrk',                           &
     'Specific humidity of forced detrained air' // substepchar_long,     &
     'kg kg-1', t_avg, d_all, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm3d(i,j,:) = scm_convss_dg(i,j) % thvrk_excess(:)
  END DO
END DO
CALL scmoutput( tmpscm3d, substepchar // 'thvrk_excess',                  &
     'Buoyancy excess of forced detrained air' // substepchar_long,       &
     'K', t_avg, d_all, default_streams, '', RoutineName )



! Status of deep, shallow and mid-level convection

DO j = 1, rows
  DO i = 1, row_length
    tmpscm2d(i,j) = REAL( scm_convss_dg(i,j) % status_deep )
  END DO
END DO
CALL scmoutput( tmpscm2d, substepchar // 'status_deep',                   &
     'Deep convection status: 0=undiagnosed, 1=failed-ascent, ' //        &
     '2=zero-closure, 3=active'//substepchar_long,                        &
     'dimensionless', t_inst, d_point, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm2d(i,j) = REAL( scm_convss_dg(i,j) % status_shallow )
  END DO
END DO
CALL scmoutput( tmpscm2d, substepchar // 'status_shallow',                &
     'Shallow convection status: 0=undiagnosed, 1=failed-ascent, ' //     &
     '2=zero-closure, 3=active'//substepchar_long,                        &
     'dimensionless', t_inst, d_point, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm2d(i,j) = REAL( scm_convss_dg(i,j) % status_mid )
  END DO
END DO
CALL scmoutput( tmpscm2d, substepchar // 'status_mid',                    &
     'Mid-level convection status: 0=undiagnosed, 1=failed-ascent, ' //   &
     '2=zero-closure, 3=active'//substepchar_long,                        &
     'dimensionless', t_inst, d_point, default_streams, '', RoutineName )



! Surface precip from each convection type

DO j = 1, rows
  DO i = 1, row_length
    tmpscm2d(i,j) = scm_convss_dg(i,j) % precip_deep
  END DO
END DO
CALL scmoutput( tmpscm2d, substepchar // 'precip_deep',                   &
     'Deep convection precipitation rate' // substepchar_long,            &
     'kg m-2 s-1', t_inst, d_point, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm2d(i,j) = scm_convss_dg(i,j) % precip_shallow
  END DO
END DO
CALL scmoutput( tmpscm2d, substepchar // 'precip_shallow',                &
     'Shallow convection precipitation rate' // substepchar_long,         &
     'kg m-2 s-1', t_inst, d_point, default_streams, '', RoutineName )

DO j = 1, rows
  DO i = 1, row_length
    tmpscm2d(i,j) = scm_convss_dg(i,j) % precip_mid
  END DO
END DO
CALL scmoutput( tmpscm2d, substepchar // 'precip_mid',                    &
     'Mid-level convection precipitation rate' // substepchar_long,       &
     'kg m-2 s-1', t_inst, d_point, default_streams, '', RoutineName )



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE scm_convss_dg_output



END MODULE scm_convss_dg_mod

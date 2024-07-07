! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module containing a type definition for SCM diagnostics from convection,
!  and subroutines to initialise and output the diagnostics
!  (stub versions which do nothing)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model


MODULE scm_convss_dg_mod

IMPLICIT NONE

!------------------------------------------------------------------------
! Type definition for SCM convection diagnostics
!------------------------------------------------------------------------
! NOTE: even though these arrays will never be allocated or used in the 
! full UM, the components need to be declared in order to prevent an error
! during compile
TYPE scm_convss_dg_type

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


SUBROUTINE scm_convss_dg_allocate( scm_convss_dg )
IMPLICIT NONE
TYPE(scm_convss_dg_type) :: scm_convss_dg
RETURN
END SUBROUTINE scm_convss_dg_allocate


SUBROUTINE scm_convss_dg_deallocate( scm_convss_dg )
IMPLICIT NONE
TYPE(scm_convss_dg_type) :: scm_convss_dg
RETURN
END SUBROUTINE scm_convss_dg_deallocate


SUBROUTINE scm_convss_dg_initzero( scm_convss_dg )
USE nlsizes_namelist_mod, ONLY: row_length, rows
IMPLICIT NONE
TYPE(scm_convss_dg_type) :: scm_convss_dg(row_length,rows)
RETURN
END SUBROUTINE scm_convss_dg_initzero


SUBROUTINE scm_convss_dg_output( scm_convss_dg, call_number )
USE nlsizes_namelist_mod, ONLY: row_length, rows
IMPLICIT NONE
TYPE(scm_convss_dg_type) :: scm_convss_dg(row_length,rows)
INTEGER :: call_number
RETURN
END SUBROUTINE scm_convss_dg_output


END MODULE scm_convss_dg_mod

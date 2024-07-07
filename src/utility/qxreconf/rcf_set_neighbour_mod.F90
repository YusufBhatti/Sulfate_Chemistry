! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Parallel RCF: sets array of tids of adjacent processes.

MODULE Rcf_Set_Neighbour_Mod

!  Subroutine Rcf_Set_Neighbour - What are my neighbouring processors?
!
! Description:
! This routine finds the tids of the North, South, East and West
! neighbouring processes. It takes account of the boundary
! condition in each dimension (X and Y) which can be either:
! cyclic : wrap around
! static : no wrap around
!
! Method:
! The tid of each neighbouring process is calculated (taking into
! account the relevant boundary conditions) and placed in the
! neighbour array.
!
! Derived from UM4.5 code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SET_NEIGHBOUR_MOD'

CONTAINS

SUBROUTINE Rcf_Set_Neighbour(decomp_type)

USE UM_ParVars
USE Decomp_DB
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine Arguments:
INTEGER :: decomp_type  ! Decomposition type to update neighbours for

CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SET_NEIGHBOUR'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ------------------------------------------------------------------

! Set Northern neighbour

IF ( decompDB (decomp_type) % g_gridpos (2,mype) > 0) THEN
  !       This processor is not at the top of the LPG
  decompDB (decomp_type) % neighbour (PNorth) = &
      mype - decompDB (decomp_type) % gridsize (1)
ELSE IF (decompDB (decomp_type) % bound (2) == bc_cyclic) THEN
  !       This processor at the top of LPG, and has cyclic BCs.
  decompDB (decomp_type) % neighbour (PNorth) = &
      mype + decompDB (decomp_type) % nproc  -  &
      decompDB (decomp_type) % gridsize (1)
ELSE
  !       This processor at top of LPG and has static BCs
  decompDB (decomp_type) % neighbour (PNorth) =NoDomain
END IF

! Set Southern neighbour

IF ( decompDB (decomp_type) % g_gridpos (2,mype) < &
    (decompDB (decomp_type) % gridsize (2)-1) ) THEN
  !       This processor is not at the bottom of the LPG
  decompDB (decomp_type) % neighbour (PSouth) =    &
      mype + decompDB (decomp_type) % gridsize (1)
ELSE IF (decompDB (decomp_type) % bound (2) == bc_cyclic) THEN
  !       This processor at the bottom of LPG, and has cyclic BCs.
  decompDB (decomp_type) % neighbour (PSouth) =    &
      mype - decompDB (decomp_type) % nproc  +     &
      decompDB (decomp_type) % gridsize (1)
ELSE
  !       This processor at top of LPG and has static BCs
  decompDB (decomp_type) % neighbour (PSouth) =NoDomain
END IF

! Set Western neighbour

IF ( decompDB (decomp_type) % g_gridpos (1,mype) > 0) THEN
  !       This processor is not at the left of the LPG
  decompDB (decomp_type) % neighbour (PWest) =     &
      mype - 1
ELSE IF (decompDB (decomp_type) % bound (1) == bc_cyclic) THEN
  !       This processor at the left of the LPG, and has cyclic BCs.
  decompDB (decomp_type) % neighbour (PWest) =     &
      mype + decompDB (decomp_type) % gridsize (1) - 1
ELSE
  !       This processor at top of LPG and has static BCs
  decompDB (decomp_type) % neighbour (PWest) =NoDomain
END IF

! Set Eastern neighbour
IF ( decompDB (decomp_type) % g_gridpos (1,mype) < &
    (decompDB (decomp_type) % gridsize (1)-1) ) THEN
  !       This processor is not at the right of the LPG
  decompDB (decomp_type) % neighbour (PEast) =     &
      mype + 1
ELSE IF (decompDB (decomp_type) % bound (1) == bc_cyclic) THEN
  !       This processor at the left of the LPG, and has cyclic BCs.
  decompDB (decomp_type) % neighbour (PEast) =     &
      mype - decompDB (decomp_type) % gridsize (1) + 1
ELSE
  !       This processor at top of LPG and has static BCs
  decompDB (decomp_type) % neighbour (PEast) =NoDomain
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Set_Neighbour

END MODULE Rcf_Set_Neighbour_Mod



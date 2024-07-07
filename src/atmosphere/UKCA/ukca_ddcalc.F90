! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Combine surface resistance rc with aerodynamic resistance ra and
!  quasi-laminar resistance rb to get overall dry deposition velocity.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_DDEPCTL.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_ddcalc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_DDCALC_MOD'

CONTAINS

SUBROUTINE ukca_ddcalc(row_length, rows, bl_levels, timestep,                 &
                       dzl, zbl, gsf, ra, rb, rc, nlev_with_ddep, zdryrt)


USE asad_mod,                ONLY: &
    ndepd

USE ereport_mod,             ONLY: &
    ereport

USE errormessagelength_mod,  ONLY: &
    errormessagelength

USE jules_surface_types_mod, ONLY: &
    ntype,                         &
    npft

USE parkind1,                ONLY: &
    jprb,                          &
    jpim

USE ukca_option_mod,         ONLY: &
    jpdd, l_ukca_ddep_lev1

USE yomhook,                 ONLY: &
    lhook,                         &
    dr_hook


IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length, rows
INTEGER, INTENT(IN) :: bl_levels

REAL, INTENT(IN) :: timestep
REAL, INTENT(IN) :: zbl(row_length,rows)             ! boundary layer depth
REAL, INTENT(IN) :: dzl(row_length, rows, bl_levels) ! thickness of BL levels
REAL, INTENT(IN) :: gsf(row_length,rows,ntype)       ! surface heat flux
REAL, INTENT(IN) :: ra(row_length,rows,ntype)        ! aerodynamic resistance
REAL, INTENT(IN) :: rb(row_length,rows,jpdd)         ! quasi-laminar resistance
REAL, INTENT(IN) :: rc(row_length,rows,ntype,jpdd)   ! surface resistance

! no of levels over which dry deposition acts
INTEGER, INTENT(OUT) :: nlev_with_ddep(row_length,rows)
REAL, INTENT(OUT)    :: zdryrt(row_length,rows,jpdd) ! dry deposition rate

INTEGER :: i, j, k, l, n
INTEGER :: errcode    ! Error code for ereport
CHARACTER (LEN=errormessagelength) :: cmessage  ! error message
REAL :: dd
REAL :: layer_depth(row_length,rows)
REAL :: vd(row_length,rows,ntype,jpdd)    ! deposition velocity

REAL :: r_nodep = 1.0e40

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DDCALC'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!     Set all arrays to zero
DO j = 1, jpdd
  DO k = 1, rows
    DO i = 1, row_length
      zdryrt(i,k,j) = 0.0
    END DO
  END DO
END DO

DO j = 1, jpdd
  DO n = 1, ntype
    DO k = 1, rows
      DO i = 1, row_length
        vd(i,k,n,j) = 0.0
      END DO
    END DO
  END DO
END DO


! Two options:
!
! (A) dry deposition confined to lowest model layer:
IF (l_ukca_ddep_lev1) THEN 
    nlev_with_ddep(:,:) = 1
    layer_depth(:,:) = dzl(:,:,1)
ELSE

! (B) If dry deposition applied everywhere within the boundary
!     layer then look for the highest model level completely
!     contained in it.
    DO k = 1, rows
      DO i = 1, row_length
        layer_depth(i,k) = dzl(i,k,1)
        nlev_with_ddep(i,k) = 1
      END DO
    END DO

    DO l = 2, bl_levels
      DO k = 1, rows
        DO i = 1, row_length
          dd = layer_depth(i,k) + dzl(i,k,l)
          IF (dd < zbl(i,k)) THEN
            layer_depth(i,k) = dd
            nlev_with_ddep(i,k) = l
          END IF
        END DO
      END DO
    END DO

END IF 

!     Calculate overall dry deposition velocity [vd = 1/(ra + rb + rc)]
!     Do vegetated tiles first. Quasi-laminar resistance pre-multiplied by
!     ln[z0m/z0] = 2.0 for vegetated areas, or 1.0 for smooth surfaces
!     See Ganzeveld & Lelieveld, JGR 1995 Vol.100 No. D10 pp.20999-21012.

DO k = 1, rows
  DO j = 1, ndepd
    DO n = 1, npft
      DO i = 1, row_length
        IF (rc(i,k,n,j) < r_nodep .AND. gsf(i,k,n) > 0.0) THEN
          vd(i,k,n,j) = 1.0 /                                    &
            (ra(i,k,n) + 2.0*rb(i,k,j) + rc(i,k,n,j))
        END IF
      END DO
    END DO

    !         Now do calculation for non-vegetated tiles

    DO n = npft+1, ntype
      DO i = 1, row_length
        IF (rc(i,k,n,j) < r_nodep .AND. gsf(i,k,n) > 0.0) THEN
          vd(i,k,n,j) = 1.0 /                                    &
            (ra(i,k,n) + rb(i,k,j) + rc(i,k,n,j))
        END IF
      END DO
    END DO
  END DO
END DO

!     VD() now contains dry deposition velocities for each tile
!     in each grid sq. Calculate overall first-order loss rate
!     over time "timestep" for each tile and sum over all tiles
!     to obtain overall first-order loss rate zdryrt().

DO n = 1, ntype
  DO j = 1, ndepd
    DO k = 1, rows
      DO i = 1, row_length
        IF (vd(i,k,n,j) > 0.0) THEN
          zdryrt(i,k,j) = zdryrt(i,k,j) + gsf(i,k,n) *           &
            (1.0-EXP(-vd(i,k,n,j) * timestep / layer_depth(i,k)))
        END IF
      END DO
    END DO
  END DO
END DO

!     ZDRYRT() contains loss rate over time "timestep".
!     Divide by timestep to get rate in s-1.

DO j = 1, ndepd
  DO k = 1, rows
    DO i = 1, row_length
      zdryrt(i,k,j) = -LOG(1.0 - zdryrt(i,k,j)) / timestep
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_ddcalc
END MODULE ukca_ddcalc_mod

! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  Calculate lots of different partitions for convective analysis

MODULE crmstyle_cal_weights_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_CAL_WEIGHTS_MOD'

CONTAINS
! ------------------------------------------------------------------------------
! Description:
!  Calculate weights for vertical interpolation to a set of fixed heights
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CRMstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------
SUBROUTINE crmstyle_cal_weights(ncols,nrows,nlevs,r_in,orog,r_out,             &
                       k_m_level,weights)

USE planet_constants_mod, ONLY: planet_radius
USE word_sizes_mod, ONLY: iwp, wp

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: &
  ncols                & ! columns
 ,nrows                & ! rows
 ,nlevs                  ! levels

REAL, INTENT(IN) ::          &
  r_in(ncols,nrows,nlevs)      ! Input heights (full precision must be used)

REAL(wp), INTENT(IN) ::      &
  orog(ncols,nrows)          & ! orography height
 ,r_out(nlevs)                 ! Required fixed output heights


INTEGER(iwp), INTENT(OUT) ::    &
  k_m_level(ncols,nrows,nlevs)    ! level below
REAL(wp), INTENT(OUT) ::        &
  weights(ncols,nrows,nlevs)      ! Weights for interpolation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

INTEGER :: i,j,k,kk
! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CRMSTYLE_CAL_WEIGHTS'

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
!-------------------------------------------------------------------------
! Find level below


!$OMP PARALLEL DO PRIVATE(i,j,k,kk) DEFAULT(NONE)                      &
!$OMP& SHARED(nlevs, nrows, ncols, k_m_level, r_out, orog, r_in,       &
!$OMP&  weights, planet_radius)

DO k=1,nlevs               ! output levels
  DO j=1,nrows
    DO i=1,ncols
      k_m_level(i,j,k) = -99     ! indicate below model lowest level
      weights(i,j,k)   = 0.0     ! initialise to zero
    END DO
  END DO

  DO j=1,nrows
    DO i=1,ncols
      IF (r_out(k) >  (orog(i,j)+planet_radius) ) THEN

        IF (r_out(k) > (orog(i,j)+planet_radius)                    &
                            .AND. r_out(k) <= r_in(i,j,1) ) THEN
          k_m_level(i,j,k) = 0
        END IF

        DO kk=1,nlevs-1          ! input levels
          IF (r_out(k) > r_in(i,j,kk) .AND. r_out(k) <= r_in(i,j,kk+1) ) THEN
            k_m_level(i,j,k) = kk
            weights(i,j,k) = (r_out(k)      -r_in(i,j,kk))/            &
                             (r_in(i,j,kk+1)-r_in(i,j,kk))
          END IF
        END DO

        ! New level above top  - Only likely to occur for 1 level
        IF (r_out(k) > r_in(i,j,nlevs)) THEN
          k_m_level(i,j,k) = nlevs
        END IF

      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_cal_weights
END MODULE crmstyle_cal_weights_mod

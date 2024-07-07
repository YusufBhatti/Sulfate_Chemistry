! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  put input field on fixed height from hybrid heights.

MODULE put_on_fixed_heights_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PUT_ON_FIXED_HEIGHTS_MOD'

CONTAINS

SUBROUTINE put_on_fixed_heights(ncols,nrows,nlevs,k_m_level,f_in,weights,&
                               l_inter,f_out)

USE word_sizes_mod, ONLY: iwp,wp

USE missing_data_mod, ONLY: rmdi

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
! ------------------------------------------------------------------------------
! Description:
!   Read in pp file of model output
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utilty - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------

INTEGER, INTENT(IN) :: &
  ncols                & ! columns
 ,nrows                & ! rows
 ,nlevs                  ! levels

INTEGER(iwp), INTENT(IN) ::     &
  k_m_level(ncols,nrows,nlevs)    ! level below

REAL(wp), INTENT(IN) ::         &
  f_in(ncols,nrows,nlevs)       & ! Input field
 ,weights(ncols,nrows,nlevs)      ! Weights for interpolation

LOGICAL, INTENT(IN) ::          &
  l_inter                         ! .true. input to be interpolated
                                  ! .false. copy original to output field
                                  ! used for specail case of flat grid, theta
                                  ! level fields.

REAL(wp), INTENT(OUT) ::        &
  f_out(ncols,nrows,nlevs)        ! Output field

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
INTEGER :: i,j,k,kk

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PUT_ON_FIXED_HEIGHTS'

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

IF (l_inter) THEN
  ! Find level below

!$OMP PARALLEL DO PRIVATE(i,j,k, kk) DEFAULT(NONE)                 &
!$OMP& SHARED(nlevs, nrows, ncols, k_m_level, f_in, f_out, weights)

  DO k=1,nlevs               ! output levels

    DO j=1,nrows
      DO i=1,ncols
        IF (k_m_level(i,j,k) == 0) THEN    ! above surface

          f_out(i,j,k) = f_in(i,j,1)       ! same as level 1 value

        ELSE IF ( k_m_level(i,j,k) > 0 .AND. k_m_level(i,j,k) < nlevs) THEN

          kk=k_m_level(i,j,k)
          f_out(i,j,k) = f_in(i,j,kk)*(1.0-weights(i,j,k)) +            &
                                  f_in(i,j,kk+1)*weights(i,j,k)

        ELSE IF ( k_m_level(i,j,k) == nlevs) THEN  ! can       get a value

          kk=k_m_level(i,j,k)
          f_out(i,j,k) = f_in(i,j,kk)   ! set to top level value

        ELSE           ! below surface

          f_out(i,j,k) = rmdi

        END IF
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE ! copy field to output field

  DO k=1,nlevs               ! output levels
    DO j=1,nrows
      DO i=1,ncols
        f_out(i,j,k) = f_in(i,j,k)
      END DO
    END DO
  END DO

END IF
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------

RETURN
END SUBROUTINE put_on_fixed_heights

END MODULE put_on_fixed_heights_mod

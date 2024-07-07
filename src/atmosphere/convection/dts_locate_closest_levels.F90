! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Locates the element of 'array' that is closest to 'val'
!
MODULE dts_locate_closest_levels_mod

IMPLICIT NONE

CHARACTER(LEN=*),                                                           &
      PARAMETER, PRIVATE :: ModuleName = 'DTS_LOCATE_CLOSEST_LEVELS_MOD'
CONTAINS

SUBROUTINE dts_locate_closest_levels(nx,ny,nv,nl,l_increase,array,vals      &
                                     ,klev,klevm,klevp)
! Modules none used


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr
IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! Locates the element of 'array' that is closest to 'val'
! Options are: (nx=1, nv=nl) or (nv=1, nx=nl) or (nx=nv=nl)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
! ------------------------------------------------------------------------------

! Subroutine arguments

INTEGER, INTENT(IN) :: &
  nx,ny,nv,nl            ! Dimensions

LOGICAL, INTENT(IN) :: &
  l_increase             ! does array increase or decrease with height
                         !(true = increase, false = decrease, eg temperature)

REAL, INTENT(IN) ::    &
  array(nx,ny)         & ! Array of values
  ,vals(nv)              ! Values to find in array

INTEGER, INTENT(OUT) :: &
  klev(nl)              & ! closest level k
 ,klevm(nl)             & ! k-1
 ,klevp(nl)               ! k+1

! Local variables

INTEGER ::        &
  i,k                ! loop counters

REAL  ::          &
  EPSILON(nl)     &  ! difference array - val
 ,besteps(nl)        ! smallest Abs |epsilon|

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DTS_LOCATE_CLOSEST_LEVELS'

!------------------------------------------------------------------------------
! Initialise besteps and klev

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (nv == nx .AND. nx /= nl) THEN  
  WRITE(umMessage,*) 'ERROR in locate_clsest'
  CALL umPrint(umMessage,src='dts_locate_closest_levels')
END IF


! attempt to speed up code
klev(:) = 1

IF (nx == 1) THEN

  besteps(:) = ABS(array(1,1)-vals(:))
  DO k=2,ny
    DO i=1,nl
      EPSILON(i) = (array(1,k)-vals(i))
      IF (ABS(EPSILON(i)) < besteps(i)) THEN
        besteps(i) = ABS(EPSILON(i))
        klev(i) = k
        ! depending on whether the field increases with height or decreases, if
                         ! the best fit value is negative or positive,  it
                         ! tells us whether the closest level is above or below it
        IF ((EPSILON(i) < 0.0 .AND. l_increase) .OR.                 &
           (EPSILON(i) > 0.0 .AND. .NOT. l_increase) ) THEN
          klevm(i) = k
        ELSE
          klevm(i) = k-1
        END IF
        klevp(i) = klevm(i)+1
      END IF
    END DO
  END DO

ELSE IF (nv == 1) THEN

  besteps(:) = ABS(array(:,1)-vals(1))
  DO k=2,ny
    DO i=1,nl

      EPSILON(i) = (array(i,k)-vals(1))

      IF (ABS(EPSILON(i)) < besteps(i)) THEN
        besteps(i) = ABS(EPSILON(i))
        klev(i) = k
        ! depending on whether the field increases with height or decreases, if
                         ! the best fit value is negative or positive,  it
                         ! tells us whether the closest level is above or below it
        IF ((EPSILON(i) < 0.0 .AND. l_increase) .OR.                 &
           (EPSILON(i) > 0.0 .AND. .NOT. l_increase) ) THEN
          klevm(i) = k
        ELSE
          klevm(i) = k-1
        END IF
        klevp(i) = klevm(i)+1
      END IF
    END DO
  END DO

ELSE IF (nv == nx .AND. nx == nl) THEN

  besteps(:) = ABS(array(:,1)-vals(:))

  DO k=2,ny
    DO i=1,nl
      EPSILON(i) = (array(i,k)-vals(i))
      IF (ABS(EPSILON(i)) < besteps(i)) THEN
        besteps(i) = ABS(EPSILON(i))
        klev(i) = k
        ! depending on whether the field increases with height or decreases, if
                         ! the best fit value is negative or positive,  it
                         ! tells us whether the closest level is above or below it
        IF ((EPSILON(i) < 0.0 .AND. l_increase) .OR.                 &
           (EPSILON(i) > 0.0 .AND. .NOT. l_increase) ) THEN
          klevm(i) = k
        ELSE
          klevm(i) = k-1
        END IF
        klevp(i) = klevm(i)+1
      END IF
    END DO
  END DO
END IF    ! tests input dimensions
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN


END SUBROUTINE dts_locate_closest_levels
END MODULE dts_locate_closest_levels_mod

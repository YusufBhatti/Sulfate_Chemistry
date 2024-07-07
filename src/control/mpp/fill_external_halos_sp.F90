! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Routine for filling in External Halos in single precision

SUBROUTINE fill_external_halos_sp(                                 &
  field,row_length,rows,levels,                                   &
  halo_x, halo_y)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParParams, ONLY: nodomain, peast, pwest, pnorth, psouth
USE um_types,              ONLY: real32

IMPLICIT NONE



! Purpose:
!  Fills external halos (those around the edge of model domain) with
!  sensible (copy of interior points) numbers.
!  This is useful for LAM models when SWAPBOUNDS doesn't fill these
!  halos as they are normally filled with LBC data. However, not all
!  fields have LBC data applied to them, and such fields need to have
!  sensible numbers put in these external halo regions.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Arguments:

INTEGER, INTENT(IN) :: row_length ! IN: number of points on a row
                                  !     (not including halos)
INTEGER, INTENT(IN) :: rows       ! IN: number of rows in a theta field
                                  !     (not including halos)
INTEGER, INTENT(IN) :: levels     ! IN: number of model levels
INTEGER, INTENT(IN) :: halo_x     ! IN: size of halo in "i" direction
INTEGER, INTENT(IN) :: halo_y     ! IN: size of halo in "j" direction

REAL(KIND=real32), INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,  &
                             1-halo_y:rows+halo_y,                &
                                                 levels)
                                  ! IN/OUT : Field to have its halos updated

! Local variables

INTEGER  ::  i,j,k            ! loop indicies
INTEGER  ::  j_start,j_stop   ! loop indicies

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILL_EXTERNAL_HALOS_SP'

!----------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Loop limits changed to stop use of undefined data
IF (neighbour(pwest) == nodomain .OR. neighbour(peast) == nodomain) THEN
  IF (neighbour(psouth)  ==  nodomain) THEN
    j_start=1
  ELSE
    j_start=1-halo_y
  END IF
  IF (neighbour(pnorth)  ==  nodomain) THEN
    j_stop=rows
  ELSE
    j_stop=rows+halo_y
  END IF
END IF

! parameters used pwest, peast, psouth, pnorth, nodomain
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                          &
!$OMP& SHARED(neighbour,levels,rows,row_length,halo_y,halo_x,         &
!$OMP& field,j_start,j_stop)

! Western halo region
IF (neighbour(pwest)  ==  nodomain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = j_start, j_stop
      DO i = 1-halo_x, 0
        field(i,j,k)=field(1,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PWest)  ==  NoDomain)

! Eastern halo region
IF (neighbour(peast)  ==  nodomain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = j_start, j_stop
      DO i = row_length+1, row_length+halo_x
        field(i,j,k)=field(row_length,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PEast)  ==  NoDomain)

! Northern halo region
IF (neighbour(pnorth)  ==  nodomain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = rows+1, rows+halo_y
      DO i = 1-halo_x, row_length+halo_x
        field(i,j,k)=field(i,rows,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PNorth)  ==  NoDomain)

! Southern halo region
IF (neighbour(psouth)  ==  nodomain) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = 1-halo_y, 0
      DO i = 1-halo_x, row_length+halo_x
        field(i,j,k)=field(i,1,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

END IF ! IF (neighbour(PSouth)  ==  NoDomain)

! End of OpenMP parallel region
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fill_external_halos_sp

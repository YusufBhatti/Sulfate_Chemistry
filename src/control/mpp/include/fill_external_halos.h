! -- begin fill_external_halos.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

INTEGER, INTENT(IN) :: row_length ! IN: number of points on a row
                                  !     (not including halos)
INTEGER, INTENT(IN) :: rows       ! IN: number of rows in a theta field
                                  !     (not including halos)
INTEGER, INTENT(IN) :: levels     ! IN: number of model levels
INTEGER, INTENT(IN) :: halo_x     ! IN: size of halo in "i" direction
INTEGER, INTENT(IN) :: halo_y     ! IN: size of halo in "j" direction

#if   defined(FILL_EXTERNAL_HALOS_INTEGER)
INTEGER,               INTENT(INOUT) ::                                    &
#elif defined(FILL_EXTERNAL_HALOS_LOGICAL)
LOGICAL,               INTENT(INOUT) ::                                    &
#else
REAL(KIND=field_kind), INTENT(INOUT) ::                                    &
#endif
                                        field(1-halo_x:row_length+halo_x,  &
                                              1-halo_y:rows+halo_y,        &
                                              levels)
                                  ! IN/OUT : Field to have its halos updated


! Local variables

INTEGER  ::  i,j,k            ! loop indicies
INTEGER  ::  j_start,j_stop   ! loop indicies

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

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
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                          &
!$OMP SHARED(neighbour,levels,rows,row_length,halo_y,halo_x,         &
!$OMP field,j_start,j_stop)

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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
! -- end fill_external_halos.h --

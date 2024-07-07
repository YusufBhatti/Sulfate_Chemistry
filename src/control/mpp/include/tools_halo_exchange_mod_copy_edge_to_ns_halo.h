! -- begin tools_halo_exchange_mod_copy_edge_to_ns_halo.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP


INTEGER, INTENT(IN) :: row_length       ! number of points on a row
                                        !     (not including halos)
INTEGER, INTENT(IN) :: rows             ! number of rows in a theta field
                                        !     (not including halos)
INTEGER, INTENT(IN) :: levels           ! number of model levels
INTEGER, INTENT(IN) :: halo_x           ! size of halo in "i" direction
INTEGER, INTENT(IN) :: halo_y           ! size of halo in "j" direction

LOGICAL, INTENT(IN) :: edge_south       ! Copy south edge into halo
LOGICAL, INTENT(IN) :: edge_north       ! Copy north edge into halo



#if   defined(TOOLS_HALO_EXCHANGE_INTEGER)
INTEGER, INTENT(INOUT) ::                                                 &
  field(1-halo_x:row_length+halo_x,                                       &
        1-halo_y:rows+halo_y,                                             &
        levels)
#elif defined(TOOLS_HALO_EXCHANGE_LOGICAL)
LOGICAL, INTENT(INOUT) ::                                                 &
  field(1-halo_x:row_length+halo_x,                                       &
        1-halo_y:rows+halo_y,                                             &
        levels)
#else
REAL(KIND=field_kind), INTENT(INOUT) ::                                   &
  field(1-halo_x:row_length+halo_x,                                       &
        1-halo_y:rows+halo_y,                                             &
        levels)

#endif

! Local variables

INTEGER :: i,j,k

! DrHook

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)



IF (edge_south) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, halo_x, rows, row_length, field)
  DO k = 1, levels
    DO j = 1, halo_y
      DO i = 1-halo_x, row_length+halo_x
        field(i,j-halo_y,k)=field(i,1,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO   
  
END IF

IF (edge_north) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, halo_x, rows, row_length, field)
  DO k = 1, levels
    DO j = 1, halo_y
      DO i = 1-halo_x, row_length+halo_x
        field(i,rows+j,k)=field(i,rows,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO   

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
! -- end  tools_halo_exchange_mod_copy_edge_to_ns_halo.h --

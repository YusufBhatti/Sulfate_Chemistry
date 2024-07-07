! -- begin tools_halo_exchange_mod_copy_cr_halo_to_buffer.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP


INTEGER, INTENT(IN) :: row_length ! number of points on a row
                                  !     (not including halos)
INTEGER, INTENT(IN) :: rows       ! number of rows in a theta field
                                  !     (not including halos)
INTEGER, INTENT(IN) :: levels     ! number of model levels
INTEGER, INTENT(IN) :: halo_x     ! size of halo in "i" direction
INTEGER, INTENT(IN) :: halo_y     ! size of halo in "j" direction

INTEGER, INTENT(IN) :: disp_x     ! displacement in "i" direction
INTEGER, INTENT(IN) :: disp_y     ! displacement in "j" direction

#if   defined(TOOLS_HALO_EXCHANGE_INTEGER)
INTEGER, INTENT(IN) ::                                                    &
  field(1-halo_x:row_length+halo_x,                                       &
        1-halo_y:rows+halo_y,                                             &
        levels)
#elif defined(TOOLS_HALO_EXCHANGE_LOGICAL)
LOGICAL, INTENT(IN) ::                                                    &
  field(1-halo_x:row_length+halo_x,                                       &
        1-halo_y:rows+halo_y,                                             &
        levels)
#else
REAL(KIND=field_kind), INTENT(IN) ::                                      &
  field(1-halo_x:row_length+halo_x,                                       &
        1-halo_y:rows+halo_y,                                             &
        levels)
#endif

#if   defined(TOOLS_HALO_EXCHANGE_INTEGER)
INTEGER, INTENT(INOUT) :: buffer(halo_x,halo_y,levels) 
#elif defined(TOOLS_HALO_EXCHANGE_LOGICAL)
LOGICAL, INTENT(INOUT) :: buffer(halo_x,halo_y,levels) 
#else
REAL(KIND=field_kind), INTENT(INOUT) :: buffer(halo_x,halo_y,levels) 
#endif
  

! Local variables

INTEGER :: i,j,k

! DrHook

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, halo_x, field, buffer, disp_x, disp_y)
    DO k=1,levels
      DO j=1,halo_y
        DO i=1,halo_x
          buffer(i,j,k) = field(i+disp_x,j+disp_y,k)
        END DO ! I
      END DO ! J
    END DO ! K
!$OMP END PARALLEL DO 


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
! -- end  tools_halo_exchange_mod_copy_cr_halo_to_buffer.h --

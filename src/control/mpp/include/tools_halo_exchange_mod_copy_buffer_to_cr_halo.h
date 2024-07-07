! -- begin tools_halo_exchange_mod_copy_buffer_to_cr_halo.h --
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

INTEGER, INTENT(IN) :: disp_x           ! displacement in "i" direction
INTEGER, INTENT(IN) :: disp_y           ! displacement in "j" direction

LOGICAL, INTENT(IN) :: overpole_south   ! Invert the halos if at the poles
LOGICAL, INTENT(IN) :: overpole_north   ! Invert the halos if at the poles
LOGICAL, INTENT(IN) :: change_sign      ! Change the sign

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

#if   defined(TOOLS_HALO_EXCHANGE_INTEGER)
INTEGER, INTENT(IN) :: buffer( halo_x, halo_y, levels) 
#elif defined(TOOLS_HALO_EXCHANGE_LOGICAL)
LOGICAL, INTENT(IN) :: buffer( halo_x, halo_y, levels) 
#else
REAL(KIND=field_kind), INTENT(IN) :: buffer( halo_x, halo_y, levels) 
#endif


! Local variables

INTEGER :: i,j,k

! DrHook

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (overpole_south) THEN

  IF (change_sign) THEN
#if !defined(TOOLS_HALO_EXCHANGE_LOGICAL)
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, halo_x, field, buffer, disp_x)
    DO k=1,levels
      DO j=1,halo_y
        DO i=1,halo_x
          
          field(i+disp_x,1-j,k) = -buffer(i,j,k)
          
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO 
#endif
  ELSE ! don't change sign
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, halo_x, field, buffer, disp_x)
    DO k=1,levels
      DO j=1,halo_y    
        DO i=1,halo_x
          
          field(i+disp_x,1-j,k) = buffer(i,j,k)
          
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO 
  END IF !  IF (change_sign)
        
END IF

IF (overpole_north) THEN

  IF (change_sign) THEN
#if !defined(TOOLS_HALO_EXCHANGE_LOGICAL)
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, halo_x, field, buffer, disp_x, rows)
    DO k=1,levels
      DO j=1,halo_y
        DO i=1,halo_x
          
          field(i+disp_x,j+rows,k)= -buffer(i,halo_y-j+1,k) 
          
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO 
#endif
  ELSE ! don't change sign
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, halo_x, field, buffer, disp_x, rows)
    DO k=1,levels
      DO j=1,halo_y
        DO i=1,halo_x
          
          field(i+disp_x,j+rows,k)= buffer(i,halo_y-j+1,k)
          
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO 
  END IF ! IF (change_sign)

END IF

! If it is not a overpole copy

IF( .NOT. (overpole_south .OR. overpole_north) ) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, halo_x, field, buffer, disp_x, disp_y)
  DO k=1,levels
    DO j=1,halo_y
      DO i=1,halo_x                        
        field(i+disp_x,j+disp_y,k) = buffer(i,j,k) 
      END DO ! I
    END DO ! J
  END DO ! K
!$OMP END PARALLEL DO 

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
! -- end  tools_halo_exchange_mod_copy_buffer_to_cr_halo.h --

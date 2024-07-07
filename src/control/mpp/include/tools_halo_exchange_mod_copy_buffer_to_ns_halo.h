! -- begin tools_halo_exchange_mod_copy_buffer_to_ns_halo.h --
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

INTEGER, INTENT(IN) :: disp_y           ! displacement in "j" direction

INTEGER, INTENT(IN) :: halo_type        ! Full or normal halo

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

! If full halo, halo_type is equal 1, otherwise is 0
#if   defined(TOOLS_HALO_EXCHANGE_INTEGER)
INTEGER, INTENT(IN) ::                                                    &
  buffer(1-halo_x*halo_type:row_length+halo_x*halo_type,                  & 
         halo_y,levels) 
#elif defined(TOOLS_HALO_EXCHANGE_LOGICAL)
LOGICAL, INTENT(IN) ::                                                    &
  buffer(1-halo_x*halo_type:row_length+halo_x*halo_type,                  & 
         halo_y,levels) 
#else
REAL(KIND=field_kind), INTENT(IN) ::                                      &
  buffer(1-halo_x*halo_type:row_length+halo_x*halo_type,                  & 
         halo_y,levels) 
#endif

! Local variables

INTEGER :: i,j,k
INTEGER :: is,ie
INTEGER :: error_code

! DrHook

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


SELECT CASE(halo_type)
CASE(normal_halo)
  is     = 1
  ie     = row_length
CASE(full_halo)
  is     = 1 - halo_x
  ie     = row_length  + halo_x
CASE DEFAULT
  error_code = 10
  CALL ereport(RoutineName, error_code, 'Unknown halo type')
END SELECT


IF (overpole_south) THEN

  IF (change_sign) THEN
#if !defined(TOOLS_HALO_EXCHANGE_LOGICAL)
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, is, ie, field, buffer)
    DO k=1,levels
      DO j=1,halo_y
        DO i=is,ie
          
          field(i,1-j,k) = -buffer(i,j,k)
          
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO 
#endif
  ELSE ! don't change sign
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, is, ie, field, buffer)
    DO k=1,levels
      DO j=1,halo_y    
        DO i=is,ie
          
          field(i,1-j,k) = buffer(i,j,k)
          
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
!$OMP SHARED(levels, halo_y, is, ie, field, rows, buffer)
    DO k=1,levels
      DO j=1,halo_y
        DO i=is,ie
          
          field(i,j+rows,k)= -buffer(i,halo_y-j+1,k) 
          
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO 
#endif
  ELSE ! don't change sign
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, is, ie, field, rows, buffer)
    DO k=1,levels
      DO j=1,halo_y
        DO i=is,ie
          
          field(i,j+rows,k)= buffer(i,halo_y-j+1,k)
          
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
!$OMP SHARED(levels, halo_y, is, ie, field, disp_y, buffer)
  DO k=1,levels
    DO j=1,halo_y
      DO i=is,ie                        

        field(i,j+disp_y,k) = buffer(i,j,k) 

      END DO ! I
    END DO ! J
  END DO ! K
!$OMP END PARALLEL DO 
  
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
! -- end  tools_halo_exchange_mod_copy_buffer_to_ns_halo.h --

! -- begin tools_halo_exchange_mod_copy_ns_halo_single.h --
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

INTEGER, INTENT(IN) :: north_off        ! offsets to use when copying data
INTEGER, INTENT(IN) :: south_off        ! to send around poles

LOGICAL, INTENT(IN) :: overpoles        ! North/South halos need to be
                                        ! communicated over the poles
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


! Local variables

INTEGER :: i,j,k
INTEGER :: is,ie
INTEGER :: half_full_row_length
INTEGER :: half_row_length

! DrHook

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

half_full_row_length = (row_length + 2*halo_x) / 2
half_row_length      = (row_length) / 2


IF (overpole_south) THEN

  IF (change_sign) THEN
#if !defined(TOOLS_HALO_EXCHANGE_LOGICAL)
!$OMP PARALLEL DO SCHEDULE(STATIC)                           &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k)                           &
!$OMP SHARED(levels, halo_y, halo_x, half_full_row_length,   & 
!$OMP        half_row_length, field, rows, south_off)
    DO k=1,levels
      DO j=1,halo_y
        DO i=1,half_full_row_length
          field(half_row_length+i,1-j,k)=                    &
            -field(i,j+south_off,k)
          field(i-halo_x,1-j,k)=                             &
            -field(half_row_length+i-halo_x,j+south_off,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO 
#endif
  ELSE ! don't change sign
!$OMP PARALLEL DO SCHEDULE(STATIC)                           &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k)                           &
!$OMP SHARED(levels, halo_y, halo_x, half_full_row_length,   & 
!$OMP        half_row_length, field, rows, south_off)
    DO k=1,levels
      DO j=1,halo_y
        DO i=1,half_full_row_length
          field(half_row_length+i,1-j,k)=                    &
            field(i,j+south_off,k)
          field(i-halo_x,1-j,k)=                             &
            field(half_row_length+i-halo_x,j+south_off,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO 
  END IF !  IF (change_sign)
        
END IF

IF (overpole_north) THEN

  IF (change_sign) THEN
#if !defined(TOOLS_HALO_EXCHANGE_LOGICAL)
!$OMP PARALLEL DO SCHEDULE(STATIC)                           &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k)                           &
!$OMP SHARED(levels, halo_y, halo_x, half_full_row_length,   & 
!$OMP        half_row_length, field, rows, north_off)
    DO k=1,levels
      DO j=1,halo_y
        DO i=1,half_full_row_length
          
          field(half_row_length+i,rows+j,k)=                 &
            -field(i,rows-j+1-north_off,k)
          field(i-halo_x,rows+j,k)=                          &
            -field(half_row_length+i-halo_x,                 &
            rows-j+1-north_off,k)
          
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO 
#endif
  ELSE ! don't change sign
!$OMP PARALLEL DO SCHEDULE(STATIC)                           &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k)                           &
!$OMP SHARED(levels, halo_y, halo_x, half_full_row_length,   & 
!$OMP        half_row_length, field, rows, north_off)
    DO k=1,levels
      DO j=1,halo_y
        DO i=1,half_full_row_length
          
          field(half_row_length+i,rows+j,k)=                 &
            field(i,rows-j+1-north_off,k)
          field(i-halo_x,rows+j,k)=                          &
            field(half_row_length+i-halo_x,                  &
            rows-j+1-north_off,k)
          
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO 
  END IF ! IF (change_sign)

END IF

! The case of a N/S cyclic LAM on a single N/S PE.

IF (.NOT. overpoles) THEN 

!$OMP PARALLEL                                               &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k)                           &
!$OMP SHARED(levels, halo_y, halo_x, rows, row_length, field)

  ! Fill southern halo
!$OMP DO SCHEDULE(STATIC)
  DO k=1,levels
    DO j=1,halo_y
      DO i=1-halo_x,row_length+halo_x
        
        field(i,j-halo_y,k)=field(i,rows-halo_y+j,k)
      END DO
    END DO 
  END DO
!$OMP END DO NOWAIT

  ! Fill northern halo  
!$OMP DO SCHEDULE(STATIC)
  DO k=1,levels
    DO j=1,halo_y
      DO i=1-halo_x,row_length+halo_x

        field(i,rows+j,k)=field(i,j,k)
      END DO 
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
! -- end  tools_halo_exchange_mod_copy_ns_halo_single.h --

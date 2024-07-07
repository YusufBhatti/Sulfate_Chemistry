! -- begin halo_exchange_mod_swap_bounds_hub_2D.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP


INTEGER, INTENT(IN) :: row_length ! Number of points on a row 
!                                 ! (not including halos)
INTEGER, INTENT(IN) :: rows       ! Number of rows in a theta field
!                                 ! (not including halos)
INTEGER, INTENT(IN) :: levels     ! Number of model levels
INTEGER, INTENT(IN) :: halo_x     ! Size of halo in "i" direction
INTEGER, INTENT(IN) :: halo_y     ! Size of halo in "j" direction
INTEGER, INTENT(IN) :: field_type ! Defines the grid interpolation type
!                                 !     of the input FIELD (u,v or w)
LOGICAL, INTENT(IN) :: l_vector   ! TRUE:  Horizontal vector
!                                 ! FALSE: Scalar

#if   defined(HALO_EXCHANGE_INTEGER)
INTEGER, INTENT(INOUT) ::                                                  &
#elif defined(HALO_EXCHANGE_LOGICAL)
LOGICAL, INTENT(INOUT) ::                                                  &
#else
REAL(KIND=field_kind), INTENT(INOUT) ::                                    &
#endif
  field(1-halo_x:row_length+halo_x,                                        &
        1-halo_y:rows+halo_y)                                              

! The next set of logicals control whether swaps are performed for a
! given direction. If none are present, all are done. If any are
! present, only the ones given TRUE values will be done.
LOGICAL, INTENT(IN), OPTIONAL :: do_east_arg
LOGICAL, INTENT(IN), OPTIONAL :: do_west_arg
LOGICAL, INTENT(IN), OPTIONAL :: do_south_arg
LOGICAL, INTENT(IN), OPTIONAL :: do_north_arg

! The next logical controls if corners are swapped or not. 
LOGICAL, INTENT(IN), OPTIONAL :: do_corners_arg

! Local variables

! Control which swap bounds routine is used to perform the halo
! exchange.
INTEGER :: swap_bounds_method

INTEGER :: error_code

REAL(KIND=jprb)               :: zhook_handle

LOGICAL :: do_east    ! perform east swap
LOGICAL :: do_west    ! perform west swap
LOGICAL :: do_south   ! perform south swap
LOGICAL :: do_north   ! perform north swap
LOGICAL :: do_corners ! do corner swaps

!------------------------------------------------------------------
! If there are no halos, do not perform the swap bounds

IF(halo_x == 0 .AND. halo_y == 0) THEN
  RETURN
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Setup the logicals
do_east    = .TRUE.
do_west    = .TRUE.
do_south   = .TRUE.
do_north   = .TRUE.
do_corners = .TRUE.

IF (PRESENT(do_corners_arg)) do_corners=do_corners_arg

IF (PRESENT(do_east_arg)  .OR. PRESENT(do_west_arg) .OR.  &
    PRESENT(do_south_arg) .OR. PRESENT(do_north_arg) ) THEN
  do_east  = .FALSE.
  do_west  = .FALSE.
  do_south = .FALSE.
  do_north = .FALSE.

  IF (PRESENT(do_east_arg))  do_east  = do_east_arg
  IF (PRESENT(do_west_arg))  do_west  = do_west_arg
  IF (PRESENT(do_south_arg)) do_south = do_south_arg
  IF (PRESENT(do_north_arg)) do_north = do_north_arg
END IF

IF (halo_x == 0) THEN
  do_east    = .FALSE.
  do_west    = .FALSE.
  do_corners = .FALSE.
END IF

IF(halo_y == 0) THEN
  do_south   = .FALSE.
  do_north   = .FALSE.
  do_corners = .FALSE.
END IF

!------------------------------------------------------------------
! CAll the 3D swap bounds hub routine

#if   defined(HALO_EXCHANGE_INTEGER)

CALL swap_bounds_hub_int( field, row_length, rows, levels,                 & 
    halo_x, halo_y, field_type, l_vector,                                  & 
    do_east, do_west, do_south, do_north, do_corners  )

#elif defined(HALO_EXCHANGE_LOGICAL)

CALL swap_bounds_hub_log( field, row_length, rows, levels,                 & 
    halo_x, halo_y, field_type, l_vector,                                  & 
    do_east, do_west, do_south, do_north, do_corners  )

#elif defined(HALO_EXCHANGE_SINGLE)

CALL swap_bounds_hub_sp( field, row_length, rows, levels,                  & 
    halo_x, halo_y, field_type, l_vector,                                  & 
    do_east, do_west, do_south, do_north, do_corners  )

#else

CALL swap_bounds_hub_dp( field, row_length, rows, levels,                  & 
    halo_x, halo_y, field_type, l_vector,                                  & 
    do_east, do_west, do_south, do_north, do_corners  )

#endif


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

! -- end  halo_exchange_mod_swap_bounds_hub_2D.h --

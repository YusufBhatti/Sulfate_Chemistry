! -- begin non_blocking_halo_exchange_mod_begin_swap_bounds_hub.h --
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

REAL(KIND=field_kind), INTENT(INOUT) ::                                    &
  field(1-halo_x:row_length+halo_x,                                        &
        1-halo_y:rows+halo_y,                                              &
        levels)                   

! The next set of logicals control whether swaps are performed for a
! given direction. If none are present, all are done. If any are
! present, only the ones given TRUE values will be done.
LOGICAL, INTENT(IN), OPTIONAL :: do_east_arg
LOGICAL, INTENT(IN), OPTIONAL :: do_west_arg
LOGICAL, INTENT(IN), OPTIONAL :: do_south_arg
LOGICAL, INTENT(IN), OPTIONAL :: do_north_arg

! The next logical controls if corners are swapped or not. 
LOGICAL, INTENT(IN), OPTIONAL :: do_corners_arg

! The optional argument with the field node 
TYPE(nb_field_node_type), INTENT(IN), OPTIONAL, POINTER  :: field_node_arg


! Local variables

LOGICAL :: do_east    ! perform east swap
LOGICAL :: do_west    ! perform west swap
LOGICAL :: do_south   ! perform south swap
LOGICAL :: do_north   ! perform north swap
LOGICAL :: do_corners ! do corner swaps

! Old ddt variables.
INTEGER, SAVE  :: old_rows       = imdi
INTEGER, SAVE  :: old_row_length = imdi
INTEGER, SAVE  :: old_halo_y     = imdi       
INTEGER, SAVE  :: old_halo_x     = imdi       
INTEGER, SAVE  :: old_mpl_type   = imdi       
INTEGER, SAVE  :: old_byte_width = imdi       

! Used by the ddt methods
TYPE(ddt_type),SAVE :: ddt = ddt_type(imdi,imdi,imdi,imdi,imdi, & 
                                      imdi,imdi,imdi,imdi,imdi)

! The field node with all the information
TYPE(nb_field_node_type), POINTER       :: field_node

INTEGER(KIND=keytype) :: key

INTEGER :: error_code

REAL(KIND=jprb)               :: zhook_handle

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
! Do blocking swap-bounds? 
! Call the default swap bounds routine and exit

IF (nb_swap_bounds_method == nb_swap_bounds_off) THEN
  
  CALL swap_bounds( field, row_length, rows, levels,    & 
                    halo_x, halo_y,                     &
                    field_type, l_vector,               &
                    do_east, do_west,                   & 
                    do_south, do_north,                 & 
                    do_corners                          & 
                    )

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)  
  RETURN
END IF

!------------------------------------------------------------------
! Retrieve node

IF ( PRESENT(field_node_arg)) THEN

  field_node => field_node_arg

ELSE

  ! Retrieve the node of the given key, i.e the addess location.
  key = get_key( field )
  field_node => get_nb_field_node(nb_field_list, key )

END IF  

! Check that the status of the communication is ok,
! i.e. there is not any non-blocking swap bounds happening.
IF(field_node % status /= nb_ok) THEN
  error_code = 99
  CALL ereport(RoutineName,error_code,                                     &
               ' wrong staus of non-blocking communication')  
END IF


!------------------------------------------------------------------
! Save the field information

field_node % rows       =  rows
field_node % row_length =  row_length
field_node % levels     =  levels
field_node % halo_x     =  halo_x
field_node % halo_y     =  halo_y
field_node % field_type =  field_type
field_node % l_vector   =  l_vector
field_node % do_east    =  do_east
field_node % do_west    =  do_west
field_node % do_south   =  do_south
field_node % do_north   =  do_north
field_node % do_corners =  do_corners

!------------------------------------------------------------------
! If there are no halos, do not perform the swap bounds

IF(halo_x == 0 .AND. halo_y == 0) THEN

  ! The node status need to be started so that the
  ! "end_swap_bounds_hub" works.
  field_node % status = nb_started

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)  
  RETURN  
END IF


!------------------------------------------------------------------
! Run the selected method

SELECT CASE(nb_swap_bounds_method)

CASE(ddt_aao)

  !****************************!
  !          DDT AAO           !
  !****************************!
  
  !------------------------------------------------------------------
  ! Initialise the derived data types(ddts) if the ones already saved
  ! cannot be reused since the dimensions are different. The number of
  ! levels do not impact the ddts.

  IF ((rows      /=  old_rows         ) .OR.                     &
    (row_length  /=  old_row_length   ) .OR.                     &
    (halo_y      /=  old_halo_y       ) .OR.                     &
    (halo_x      /=  old_halo_x       ) .OR.                     &
    (mpl_type    /=  old_mpl_type     ) .OR.                     &
    (byte_width  /=  old_byte_width    )) THEN

    CALL free_ddt(ddt)

    CALL create_ddt_ew( row_length, rows, levels,     & 
      halo_x,halo_y, mpl_type, byte_width,            &
      normal_halo, ddt)                               
    CALL create_ddt_ns( row_length, rows, levels,     & 
      halo_x,halo_y, mpl_type, byte_width,            &
      normal_halo, ddt)                               
    CALL create_ddt_cr( row_length, rows, levels,     & 
      halo_x,halo_y, mpl_type, byte_width,            &
      ddt)                                       

    old_rows       = rows
    old_row_length = row_length
    old_halo_y     = halo_y 
    old_halo_x     = halo_x
    old_mpl_type   = mpl_type
    old_byte_width = byte_width

  END IF

  CALL begin_swap_bounds_ddt_aao(                                          &
    field,                                                                 &
    row_length, rows,     levels,                                          &
    halo_x    , halo_y,                                                    &
    field_type, l_vector,                                                  &
    do_east   , do_west,                                                   &
    do_south  , do_north,                                                  &
    do_corners, ddt,                                                       &
    field_node % nreq_recv,  field_node % ireq_recv(1:8),                  &
    field_node % tag_index * nb_tag_base)    
  
CASE(ddt_nsew_wa)

  !****************************!
  !        DDT NSEW WA         !
  !****************************!
  
  !------------------------------------------------------------------
  ! Initialise the derived data types(ddts) if the ones already saved
  ! cannot be reused since the dimensions are different. The number of
  ! levels do not impact the ddts.

  IF ((rows      /=  old_rows         ) .OR.                     &
    (row_length  /=  old_row_length   ) .OR.                     &
    (halo_y      /=  old_halo_y       ) .OR.                     &
    (halo_x      /=  old_halo_x       ) .OR.                     &
    (mpl_type    /=  old_mpl_type     ) .OR.                     &
    (byte_width  /=  old_byte_width    )) THEN

    CALL free_ddt(ddt)

    CALL create_ddt_ew( row_length, rows, levels,     & 
      halo_x,halo_y, mpl_type, byte_width,            &
      full_halo, ddt)   ! < --- full halo
    CALL create_ddt_ns( row_length, rows, levels,     & 
      halo_x,halo_y, mpl_type, byte_width,            &
      normal_halo, ddt) ! < --- normal halo                               

    old_rows       = rows
    old_row_length = row_length
    old_halo_y     = halo_y 
    old_halo_x     = halo_x
    old_mpl_type   = mpl_type
    old_byte_width = byte_width

  END IF


  CALL begin_swap_bounds_ddt_nsew_wa(                                      &
    field,                                                                 &
    row_length, rows, levels,                                              &
    halo_x    , halo_y,                                                    &
    field_type, l_vector,                                                  &
    do_east   , do_west,                                                   &
    do_south  , do_north,                                                  &
    do_corners, ddt,                                                       &
    field_node % ns_halos_recv,                                            &
    field_node % nreq_recv,  field_node % ireq_recv(1:4),                  &
    field_node % tag_index * nb_tag_base)    

CASE DEFAULT
  error_code = 10
  CALL ereport(RoutineName, error_code,                                    & 
    'Unknown or not supported swap bounds method')

END SELECT

!------------------------------------------------------------------
! Set the status
field_node % status = nb_started

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

! -- end non_blocking_halo_exchange_mod_begin_swap_bounds_hub.h --


! -- begin non_blocking_halo_exchange_mod_end_swap_bounds_hub.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Assumed shape field to have its halos updated
REAL(KIND=field_kind), INTENT(INOUT) ::                                    &
  field(:,:,:)                   

! The optional argument with the field node 
TYPE(nb_field_node_type), INTENT(IN), OPTIONAL, POINTER  :: field_node_arg

! Local variables

INTEGER :: row_length ! Number of points on a row 
!                                 ! (not including halos)
INTEGER :: rows       ! Number of rows in a theta field
!                     ! (not including halos)
INTEGER :: levels     ! Number of model levels
INTEGER :: halo_x     ! Size of halo in "i" direction
INTEGER :: halo_y     ! Size of halo in "j" direction
INTEGER :: field_type ! Defines the grid interpolation type
!                     !     of the input FIELD (u,v or w)
LOGICAL :: l_vector   ! TRUE:  Horizontal vector
!                     ! FALSE: Scalar
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
! Do blocking swap-bounds? 
! This routine does nothing.
IF (nb_swap_bounds_method == nb_swap_bounds_off) THEN  
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

! Check that the swap bounds was started
IF(field_node % status /= nb_started) THEN
  error_code = 99
  CALL ereport(RoutineName,error_code,                                     &
               ' wrong staus of non-blocking communication')  
END IF

!------------------------------------------------------------------
! Retrieve the field information

rows       = field_node % rows       
row_length = field_node % row_length 
levels     = field_node % levels     
halo_x     = field_node % halo_x     
halo_y     = field_node % halo_y     
field_type = field_node % field_type 
l_vector   = field_node % l_vector   
do_east    = field_node % do_east    
do_west    = field_node % do_west    
do_south   = field_node % do_south   
do_north   = field_node % do_north   
do_corners = field_node % do_corners 

!------------------------------------------------------------------
! If there are no halos, do not perform the swap bounds

IF(halo_x == 0 .AND. halo_y == 0) THEN

  ! Remove the node from the list of active node.
  field_node % status  = nb_ok
  CALL remove_nb_field_node( nb_field_list, field_node)

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

  CALL end_swap_bounds_ddt_aao(                                            &
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


  CALL end_swap_bounds_ddt_nsew_wa(                                        &
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

IF (sb_model_type == mt_lam) THEN
  CALL fill_external_halos(field, row_length, rows, levels, halo_x, halo_y)
END IF

!------------------------------------------------------------------
! Remove the node from the list of active node.
field_node % status  = nb_ok
CALL remove_nb_field_node( nb_field_list, field_node)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

! -- end non_blocking_halo_exchange_mod_end_swap_bounds_hub.h --


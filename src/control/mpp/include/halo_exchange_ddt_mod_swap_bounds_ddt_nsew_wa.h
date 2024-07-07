! -- begin halo_exchange_ddt_mod_swap_bounds_ddt_nsew_wa.h --
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

#if   defined(HALO_EXCHANGE_DDT_INTEGER)
INTEGER, INTENT(INOUT) ::                                                  &
#elif defined(HALO_EXCHANGE_DDT_LOGICAL)
LOGICAL, INTENT(INOUT) ::                                                  &
#else
REAL(KIND=field_kind), INTENT(INOUT) ::                                    &
#endif
  field(1-halo_x:row_length+halo_x,                                        &
        1-halo_y:rows+halo_y,                                              &
        levels)                   
  
! The next set of logicals control whether swaps are performed for a
! given direction. 
LOGICAL, INTENT(IN) :: do_east_arg
LOGICAL, INTENT(IN) :: do_west_arg
LOGICAL, INTENT(IN) :: do_south_arg
LOGICAL, INTENT(IN) :: do_north_arg
LOGICAL, INTENT(IN) :: do_corners_arg

! Local variables

TYPE(ddt_type),SAVE :: ddt = ddt_type(imdi,imdi,imdi,imdi,imdi, & 
                                      imdi,imdi,imdi,imdi,imdi)

! Old ddt variables.
INTEGER, SAVE  :: old_rows       = imdi
INTEGER, SAVE  :: old_row_length = imdi
INTEGER, SAVE  :: old_halo_y     = imdi       
INTEGER, SAVE  :: old_halo_x     = imdi       
INTEGER, SAVE  :: old_mpl_type   = imdi       
INTEGER, SAVE  :: old_byte_width = imdi       


INTEGER  :: ireq_recv(8)
INTEGER  :: nreq_recv
INTEGER  :: ns_halos_recv             ! Number of N/S halos received
INTEGER  :: tag_base

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Integer and logical byte size is fixed at compile time.
#if   defined(HALO_EXCHANGE_DDT_INTEGER)
byte_width = umFortranIntegerSize()
#elif defined(HALO_EXCHANGE_DDT_LOGICAL)
byte_width = umFortranLogicalSize()
#endif

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

! Use always the same tag
tag_base = 0

CALL begin_swap_bounds_ddt_nsew_wa(               &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg, ddt, ns_halos_recv,             &  
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

CALL end_swap_bounds_ddt_nsew_wa(                 &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg, ddt, ns_halos_recv,             &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

! -- end halo_exchange_ddt_mod_swap_bounds_ddt_nsew_wa.h --

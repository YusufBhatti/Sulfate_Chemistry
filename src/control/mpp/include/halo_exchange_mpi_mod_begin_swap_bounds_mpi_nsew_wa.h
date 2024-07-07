! -- begin halo_exchange_mpi_mod_begin_swap_bounds_mpi_nsew_wa.h --
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

#if   defined(HALO_EXCHANGE_MPI_INTEGER)
INTEGER, INTENT(INOUT) ::                                                  &
#elif defined(HALO_EXCHANGE_MPI_LOGICAL)
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

INTEGER, INTENT(INOUT) :: ns_halos_recv           ! Number of N/S halos received
INTEGER, INTENT(INOUT) :: nreq_recv               ! number of requests
INTEGER, INTENT(INOUT) :: ireq_recv(4)            ! MPL request

INTEGER, INTENT(IN) :: tag_base ! Used to identify unique tags

! The buffers used to receive the data
#if   defined(HALO_EXCHANGE_MPI_INTEGER)
INTEGER, INTENT(INOUT) ::                                                  &
#elif defined(HALO_EXCHANGE_MPI_LOGICAL)
LOGICAL, INTENT(INOUT) ::                                                  &
#else
REAL(KIND=field_kind), INTENT(INOUT) ::                                    &
#endif
  recv_buffer_e(halo_x, 1-halo_y:rows+halo_y, levels),                     &
  recv_buffer_w(halo_x, 1-halo_y:rows+halo_y, levels),                     &
  recv_buffer_n(row_length, halo_y, levels),                               &
  recv_buffer_s(row_length, halo_y, levels)

! Local variables

! The buffers used to send the data
#if   defined(HALO_EXCHANGE_MPI_INTEGER)
INTEGER ::                                                                 &
#elif defined(HALO_EXCHANGE_MPI_LOGICAL)
LOGICAL ::                                                                 &
#else
REAL(KIND=field_kind) ::                                                   &
#endif
  send_buffer_e(halo_x, 1-halo_y:rows+halo_y, levels),                     &
  send_buffer_w(halo_x, 1-halo_y:rows+halo_y, levels),                     &
  send_buffer_n(row_length, halo_y, levels),                               &
  send_buffer_s(row_length, halo_y, levels)


INTEGER  :: ierr                    ! MPI error code

LOGICAL  :: send_south, send_north, send_east, send_west
LOGICAL  :: recv_south, recv_north, recv_east, recv_west

INTEGER  ::                &
  full_row_length,         & ! length of row including halos
  full_rows,               & ! number of rows including halo rows
  ew_halo_size,            & ! size of ew_halo exchanged
  ns_halo_size,            & ! size of ns_halo exchanged
  north_off,               & ! Offsets to use when copying data
  south_off                  ! to send around poles


! Displacement when halo data is copied into buffer data
INTEGER ::                           &
  from_halo_disp_x_east,             &
  from_halo_disp_x_west,             &
  from_halo_disp_y_north,            &
  from_halo_disp_y_south          

! Tags of the north/south messagges
INTEGER :: tag_send_pnorth
INTEGER :: tag_send_psouth
INTEGER :: tag_recv_pnorth
INTEGER :: tag_recv_psouth

! If .TRUE. the North/South halos need to be communicated over the poles
LOGICAL  :: overpoles

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Initialize output variables and compute sizes

ns_halos_recv = 0
nreq_recv = 0
ireq_recv(:) = mpl_request_null


full_row_length = row_length + 2*halo_x
full_rows       = rows + 2*halo_y
  
ew_halo_size    = halo_x * full_rows * levels
ns_halo_size    = row_length * halo_y * levels

!------------------------------------------------------------------
! Over the poles set up

IF(bound(2)  /=  bc_cyclic .AND. sb_model_type == mt_global) THEN
  overpoles = .TRUE.
ELSE
  overpoles = .FALSE.
END IF

recv_south = do_south_arg
send_north = do_south_arg

recv_north = do_north_arg
send_south = do_north_arg

recv_east  = do_east_arg
send_west  = do_east_arg

recv_west  = do_west_arg
send_east  = do_west_arg

north_off=0
south_off=0
  
tag_send_pnorth = tag_to_n
tag_recv_pnorth = tag_from_n

tag_send_psouth = tag_to_s
tag_recv_psouth = tag_from_s


IF (overpoles) THEN

  ! Set up the offsets. When copying data that is to be passed over
  ! the pole, on wind (u or v) grid, then copy data one row away
  ! from the pole
  
  IF ( field_type  ==  fld_type_v) THEN
    IF (at_extremity(pnorth)) north_off=1
    IF (at_extremity(psouth)) south_off=1
  END IF
    
  ! If the processor is at the poles, change the tags

  IF (at_extremity(pnorth)) THEN
    
    tag_send_pnorth = tag_to_n_p
    tag_recv_pnorth = tag_from_n_p
        
    recv_north = do_north_arg
    send_north = do_north_arg

  END IF

  IF (at_extremity(psouth)) THEN

    tag_send_psouth = tag_to_s_p
    tag_recv_psouth = tag_from_s_p
    
    recv_south = do_south_arg
    send_south = do_south_arg
    
  END IF

END IF

!------------------------------------------------------------------
! Set displacements

from_halo_disp_x_east  = (row_length-halo_x)
from_halo_disp_x_west  = (0)
from_halo_disp_y_north = (rows-halo_y) - north_off           
from_halo_disp_y_south = (0)           + south_off

!---------------------------------------
! Receive the data

IF (he_neighbour(pnorth)  /=  nodomain .AND.                      &
    he_neighbour(pnorth)  /=  mype .AND.                          &
    recv_north) THEN

  nreq_recv=nreq_recv+1
  ! Receive from north
  CALL mpl_irecv(recv_buffer_n, ns_halo_size,                     &
    mpl_type, he_neighbour(pnorth), tag_base + tag_recv_pnorth,   & 
    halos_comm, ireq_recv(nreq_recv), ierr)      

ELSE

  ! Since we are not going to receive this halo
  ! Increase the number of N/S halos received
  ns_halos_recv = ns_halos_recv + 1

END IF

IF (he_neighbour(psouth)  /=  nodomain .AND.                      &
    he_neighbour(psouth)  /=  mype .AND.                          &
    recv_south) THEN

  nreq_recv=nreq_recv+1
  ! Receive from south
  CALL mpl_irecv(recv_buffer_s, ns_halo_size,                     &
    mpl_type, he_neighbour(psouth), tag_base + tag_recv_psouth,   &
    halos_comm, ireq_recv(nreq_recv), ierr)

ELSE

  ! Since we are not going to receive this halo
  ! Increase the number of N/S halos received
  ns_halos_recv = ns_halos_recv + 1

END IF

IF (he_neighbour(peast)  /=  nodomain .AND. recv_east) THEN
  ! Receive from East
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(recv_buffer_e, ew_halo_size,                     &
    mpl_type, he_neighbour(peast), tag_base + tag_from_e,         &
    halos_comm, ireq_recv(nreq_recv), ierr)

END IF

IF (he_neighbour(pwest)  /=  nodomain .AND. recv_west) THEN
  ! Receive from West
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(recv_buffer_w, ew_halo_size,                     &
    mpl_type, he_neighbour(pwest), tag_base + tag_from_w,         &
    halos_comm, ireq_recv(nreq_recv), ierr)

END IF

!------------------------------
! North-South communications

IF (he_neighbour(psouth)  /=  nodomain .AND.                      &
    he_neighbour(psouth)  /=  mype .AND.                          &
    send_south) THEN

  CALL copy_ns_halo_to_buffer(field,row_length,rows,levels,       &
    halo_x,halo_y, send_buffer_s,                                 &
    from_halo_disp_y_south, normal_halo)

  ! Send south
  CALL mpl_send(send_buffer_s, ns_halo_size,                      &
    mpl_type, he_neighbour(psouth), tag_base + tag_send_psouth,   &
    halos_comm, ierr)

END IF

IF (he_neighbour(pnorth)  /=  nodomain .AND.                      &
    he_neighbour(pnorth)  /=  mype .AND.                          &
    send_north) THEN

  CALL copy_ns_halo_to_buffer(field,row_length,rows,levels,       &
    halo_x,halo_y, send_buffer_n,                                 &
    from_halo_disp_y_north, normal_halo)

  ! Send north
  CALL mpl_send(send_buffer_n, ns_halo_size,                      &
    mpl_type, he_neighbour(pnorth), tag_base + tag_send_pnorth,   &
    halos_comm, ierr)

END IF

! If neither North or South halos need to be exchanged or the corners
! do not need to be exchanged, then send immediately the eventual East
! and West halos.
IF (ns_halos_recv==2 .OR. (.NOT. do_corners_arg)) THEN
  
  !------------------------------------------------------------------
  ! East-West communications
  
  IF (he_neighbour(pwest)  /=  nodomain .AND. send_west) THEN

    CALL copy_ew_halo_to_buffer(field,row_length,rows,levels,     &
      halo_x,halo_y, send_buffer_w,                               &
      from_halo_disp_x_west, full_halo)

    ! Send West
    CALL mpl_send(send_buffer_w, ew_halo_size,                    &
      mpl_type,he_neighbour(pwest), tag_base + tag_to_w,          &
      halos_comm, ierr)

  END IF

  IF (he_neighbour(peast)  /=  nodomain .AND. send_east) THEN

    CALL copy_ew_halo_to_buffer(field,row_length,rows,levels,     &
      halo_x,halo_y, send_buffer_e,                               &
      from_halo_disp_x_east, full_halo)

    ! Send East
    CALL mpl_send(send_buffer_e, ew_halo_size,                    &
      mpl_type, he_neighbour(peast), tag_base + tag_to_e,         &
      halos_comm, ierr)

  END IF

  ! We want to do this only once
  ns_halos_recv = -1
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

! -- end halo_exchange_mpi_mod_begin_swap_bounds_mpi_nsew_wa.h --


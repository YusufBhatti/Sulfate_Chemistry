! -- begin halo_exchange_ddt_mod_begin_swap_bounds_ddt_aao.h --
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

! The derived that types. This is expected to be already initialised.
TYPE(ddt_type), INTENT(IN) :: ddt 

INTEGER, INTENT(INOUT) :: nreq_recv               ! number of requests
INTEGER, INTENT(INOUT) :: ireq_recv(8)            ! MPL request

INTEGER, INTENT(IN) :: tag_base ! Used to identify unique tags

! Local variables

INTEGER  :: ierr                    ! MPI error code

LOGICAL  :: send_south, send_north, send_east, send_west
LOGICAL  :: recv_south, recv_north, recv_east, recv_west
LOGICAL  :: do_corners

INTEGER  ::                &
  north_off,               & ! Offsets to use when copying data
  south_off                  ! to send around poles

! Tags of the north/south messagges
INTEGER :: tag_send_pnorth
INTEGER :: tag_send_psouth
INTEGER :: tag_recv_pnorth
INTEGER :: tag_recv_psouth
INTEGER :: tag_send_pne
INTEGER :: tag_send_pse
INTEGER :: tag_send_psw
INTEGER :: tag_send_pnw
INTEGER :: tag_recv_pne
INTEGER :: tag_recv_pse
INTEGER :: tag_recv_psw
INTEGER :: tag_recv_pnw

! If .TRUE. the North/South halos need to be communicated over the poles
LOGICAL  :: overpoles

! Displacement for North/South halo data
INTEGER :: disp_y_north
INTEGER :: disp_y_south          

! The derived data type used for the over poles communication. If
! halo_y is greater than one then the type is the reverse halo type,
! otherwise is the normal.
INTEGER  :: ns_type_south
INTEGER  :: ns_type_north  
INTEGER  :: cr_type_south
INTEGER  :: cr_type_north


REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Initialize variables and compute sizes

nreq_recv = 0
ireq_recv(:) = mpl_request_null
  
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

tag_send_pne    = tag_to_ne
tag_send_pse    = tag_to_se
tag_send_psw    = tag_to_sw
tag_send_pnw    = tag_to_nw

tag_recv_pne    = tag_from_ne
tag_recv_pse    = tag_from_se
tag_recv_psw    = tag_from_sw
tag_recv_pnw    = tag_from_nw

disp_y_north = rows + 1
disp_y_south = 1-halo_y

ns_type_south    = ddt%ns_halo_ext_type
ns_type_north    = ddt%ns_halo_ext_type
cr_type_south    = ddt%cr_halo_ext_type
cr_type_north    = ddt%cr_halo_ext_type

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
    
    tag_send_pne    = tag_to_ne_p
    tag_recv_pne    = tag_from_ne_p

    tag_send_pnw    = tag_to_nw_p
    tag_recv_pnw    = tag_from_nw_p
    
    recv_north = do_north_arg
    send_north = do_north_arg

    disp_y_north    = rows + halo_y

    IF (halo_y > 1) THEN
      ns_type_north    = ddt%rev_ns_halo_ext_type
      cr_type_north    = ddt%rev_cr_halo_ext_type
    END IF

  END IF

  IF (at_extremity(psouth)) THEN

    tag_send_psouth = tag_to_s_p
    tag_recv_psouth = tag_from_s_p

    tag_send_pse    = tag_to_se_p
    tag_recv_pse    = tag_from_se_p

    tag_send_psw    = tag_to_sw_p
    tag_recv_psw    = tag_from_sw_p
    
    recv_south = do_south_arg
    send_south = do_south_arg

    disp_y_south    = 0            

    IF (halo_y > 1) THEN
      ns_type_south    = ddt%rev_ns_halo_ext_type
      cr_type_south    = ddt%rev_cr_halo_ext_type
    END IF
    
  END IF

END IF

! If one direction is not swapped, then do not swap the corners.
IF ( (.NOT. do_east_arg ) .OR. (.NOT. do_west_arg) .OR. & 
     (.NOT. do_south_arg) .OR. (.NOT. do_north_arg)) THEN
  do_corners = .FALSE.
ELSE
  do_corners = do_corners_arg
END IF

!---------------------------------------
! Receive the data

IF (he_neighbour(pwest)  /=  nodomain .AND. recv_west) THEN

  ! Receive from West
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(field(1-halo_x,1,1), levels,                                  &
    ddt%ew_halo_ext_type, he_neighbour(pwest), tag_base + tag_from_w,          &
    halos_comm, ireq_recv(nreq_recv), ierr)

END IF

IF (he_neighbour(peast)  /=  nodomain .AND. recv_east) THEN

  ! Receive from East
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(field(row_length+1,1,1), levels,                              &
    ddt%ew_halo_ext_type, he_neighbour(peast), tag_base + tag_from_e,          &
    halos_comm, ireq_recv(nreq_recv), ierr)

END IF


IF (he_neighbour(pnorth)  /=  nodomain .AND.                    &
    he_neighbour(pnorth)  /=  mype .AND.                        &
    recv_north) THEN

  ! Receive from north
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(field(1,disp_y_north,1), levels,                              &
    ns_type_north, he_neighbour(pnorth), tag_base + tag_recv_pnorth,           &
    halos_comm, ireq_recv(nreq_recv), ierr)      

END IF

IF (he_neighbour(psouth)  /=  nodomain .AND.                    &
    he_neighbour(psouth)  /=  mype .AND.                        &
    recv_south) THEN

  ! Receive from south
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(field(1,disp_y_south,1), levels,                              &
    ns_type_south, he_neighbour(psouth), tag_base + tag_recv_psouth,           &
    halos_comm, ireq_recv(nreq_recv), ierr)

END IF

! Corners

IF (he_neighbour(pnortheast)  /=  nodomain .AND.                &
    he_neighbour(pnorth)  /=  mype         .AND.                &
    do_corners) THEN

  ! Receive from north east
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(field(row_length+1,disp_y_north,1), levels,                   &
    cr_type_north, he_neighbour(pnortheast), tag_base + tag_recv_pne,          &
    halos_comm, ireq_recv(nreq_recv), ierr)      

END IF

IF (he_neighbour(psoutheast)  /=  nodomain .AND.                &
    he_neighbour(psouth)  /=  mype         .AND.                &
    do_corners) THEN

  ! Receive from south east
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(field(row_length+1,disp_y_south,1), levels,                   &
    cr_type_south, he_neighbour(psoutheast), tag_base + tag_recv_pse,          &
    halos_comm, ireq_recv(nreq_recv), ierr)      

END IF

IF (he_neighbour(pnorthwest)  /=  nodomain .AND.                &
    he_neighbour(pnorth)  /=  mype         .AND.                &
    do_corners) THEN

  ! Receive from north west
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(field(1-halo_x,disp_y_north,1), levels,                       &
    cr_type_north, he_neighbour(pnorthwest), tag_base + tag_recv_pnw,          &
    halos_comm, ireq_recv(nreq_recv), ierr)      

END IF

IF (he_neighbour(psouthwest)  /=  nodomain .AND.                &
    he_neighbour(psouth)  /=  mype         .AND.                &
    do_corners) THEN

  ! Receive from south west
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(field(1-halo_x,disp_y_south,1), levels,                       &
    cr_type_south, he_neighbour(psouthwest), tag_base + tag_recv_psw,          &
    halos_comm, ireq_recv(nreq_recv), ierr)      

END IF

!---------------------------------------
! Send the data

IF (he_neighbour(psouth)  /=  nodomain .AND.                      &
    he_neighbour(psouth)  /=  mype .AND.                          &
    send_south) THEN

  ! Send south
  CALL mpl_send(field(1,1+south_off,1), levels,                                &
    ddt%ns_halo_ext_type, he_neighbour(psouth), tag_base + tag_send_psouth,    &
    halos_comm, ierr)

END IF

IF (he_neighbour(pnorth)  /=  nodomain .AND.                      &
    he_neighbour(pnorth)  /=  mype .AND.                          &
    send_north) THEN

  ! Send north
  CALL mpl_send(field(1,rows-halo_y+1-north_off,1), levels,                    &
    ddt%ns_halo_ext_type, he_neighbour(pnorth), tag_base + tag_send_pnorth,    &
    halos_comm, ierr)

END IF
  
IF (he_neighbour(pwest)  /=  nodomain .AND. send_west) THEN

  ! Send West
  CALL mpl_send(field(1,1,1), levels,                                          &
    ddt%ew_halo_ext_type, he_neighbour(pwest), tag_base + tag_to_w,            &
    halos_comm, ierr)

END IF

IF (he_neighbour(peast)  /=  nodomain .AND. send_east) THEN
  
  ! Send East
  CALL mpl_send(field(row_length-halo_x+1,1,1), levels,                        &
    ddt%ew_halo_ext_type, he_neighbour(peast), tag_base + tag_to_e,            &
    halos_comm, ierr)

END IF

!------------------------------
! Corner communications

IF (he_neighbour(pnortheast)  /=  nodomain .AND.                &
    he_neighbour(pnorth)  /=  mype         .AND.                &
    do_corners) THEN

  ! Send north east
  CALL mpl_send(field(row_length+1-halo_x,rows-halo_y+1-north_off,1), levels,  &
    ddt%cr_halo_ext_type, he_neighbour(pnortheast), tag_base + tag_send_pne,   &
    halos_comm, ierr)

END IF

IF (he_neighbour(psoutheast)  /=  nodomain .AND.                &
    he_neighbour(psouth)  /=  mype         .AND.                &
    do_corners) THEN

  ! Send south east
  CALL mpl_send(field(row_length+1-halo_x,1+south_off,1), levels,              &
    ddt%cr_halo_ext_type, he_neighbour(psoutheast), tag_base + tag_send_pse,   &
    halos_comm, ierr)

END IF

IF (he_neighbour(pnorthwest)  /=  nodomain .AND.                &
    he_neighbour(pnorth)  /=  mype         .AND.                &
    do_corners) THEN

  ! Send north west
  CALL mpl_send(field(1,rows-halo_y+1-north_off,1), levels,                    &
    ddt%cr_halo_ext_type, he_neighbour(pnorthwest), tag_base + tag_send_pnw,   &
    halos_comm, ierr)

END IF

IF (he_neighbour(psouthwest)  /=  nodomain .AND.                &
    he_neighbour(psouth)  /=  mype         .AND.                &
    do_corners) THEN

  ! Send south west
  CALL mpl_send(field(1,1+south_off,1), levels,                                &
    ddt%cr_halo_ext_type, he_neighbour(psouthwest), tag_base + tag_send_psw,   &
    halos_comm, ierr)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

! -- end halo_exchange_ddt_mod_begin_swap_bounds_ddt_aao.h --

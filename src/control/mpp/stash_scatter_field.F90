! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Scatters STASHed data from one processor to many processors

MODULE stash_scatter_field_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'STASH_SCATTER_FIELD_MOD'

CONTAINS
! Subroutine interface:
SUBROUTINE stash_scatter_field (                                  &
  local_field , global_field ,                                    &
  local_size, global_size, levels,                                &
  global_north , global_east_in , global_south , global_west,     &
  gridtype_code, halo_type,                                       &
  scatter_pe,                                                     &
  icode, cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE gcom_mod
USE UM_ParVars
USE UM_ParCore,  ONLY: mype, nproc_max
USE Field_Types, ONLY: fld_type_u, fld_type_v, fld_type_unknown
USE mpp_conf_mod, ONLY: include_halos_ew, include_halos_ns
USE cppxref_mod, ONLY:                                            &
    ppx_atm_tzonal
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE missing_data_mod, ONLY: imdi
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE



! Description:
! Takes a decomposed, STASH processed field and gathers
! it to a single processor, ready for I/O,

! Method:
! See in-line documentation

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Subroutine arguments:

INTEGER, INTENT(IN) :: local_size      ! IN: size of level of LOCAL_FIELD
INTEGER, INTENT(IN) :: global_size     ! IN: size of level of GLOBAL_FIELD
INTEGER, INTENT(IN) :: levels          ! IN: number of levels
INTEGER, INTENT(IN) :: global_north    ! IN: specification of subdomain boundaries
INTEGER, INTENT(IN) :: global_east_in  ! IN: ""
INTEGER, INTENT(IN) :: global_south    ! IN: ""
INTEGER, INTENT(IN) :: global_west     ! IN: ""
INTEGER, INTENT(IN) :: gridtype_code   ! IN: indicates the type of grid output
INTEGER, INTENT(IN) :: halo_type       ! IN: type of halo on this field
INTEGER, INTENT(IN) :: scatter_pe      ! IN: the PE to scatter global field from

INTEGER, INTENT(OUT) :: icode          ! OUT: return code, 0=OK

REAL, INTENT(OUT) :: local_field(local_size,levels)
                                       ! OUT : local scattered data
REAL, INTENT(IN)  :: global_field(global_size,levels)
                                       ! IN : (PE SCATTER_PE only) - full field

CHARACTER(LEN=errormessagelength) :: cmessage  
                                       ! OUT: Error message if ICODE  /=  0

! Local variables

INTEGER  ::  global_east        ! copy of GLOBAL_EAST_IN with wrap around s.t.
                                !         GLOBAL_EAST > GLOBAL_ROW_LEN
INTEGER  ::  global_x           ! size of global data EW
INTEGER  ::  global_y           ! size of global data NS
INTEGER  ::  fld_type           ! indicates if field is on P or U grid
INTEGER  ::  level              ! loop index for loop over levels
INTEGER  ::  i                  ! loop counter
INTEGER  ::  proc_topleft_x,proc_topleft_y    ! processors at corners of
INTEGER  ::  proc_botright_x,proc_botright_y  ! the subarea
INTEGER  ::  dummy1,dummy2      ! ignored return arguments
INTEGER  ::  procx,procy        ! loop indexes for loops over processors
INTEGER  ::  eff_procx          ! real x co-ord of processor column procx
INTEGER  ::  procid             ! processor id of (procx,procy)
INTEGER  ::  local_xstart,local_xend          ! boundaries of subdomain for
INTEGER  ::  local_ystart,local_yend          ! processor procid
INTEGER  ::  local_start_row    ! first row to receive on procid
INTEGER  ::  local_start_col    ! first column to receive on procid
INTEGER  ::  sendsize_x         ! number of points on each row to send to procid
INTEGER  ::  nrows_to_send      ! number of rows to send to procid
INTEGER  ::  local_row_length   ! size of receiving array EW
INTEGER  ::  global_start_row   ! first row to send on SCATTER_PE
INTEGER  ::  global_start_col   ! first col. to send on SCATTER_PE
INTEGER  ::  global_row_length  ! size of sending array EW
INTEGER  ::  flag,info          ! GCOM arguments

LOGICAL  ::  l_vec              ! Indicates if a field is a vector quantity
                                ! (SWAPBOUNDS argument)
! Copies of arguments / variables used to decide if we can use the
! send/receive maps used in the last call

! Save all the variables that may be used in the next call

INTEGER, SAVE ::                                                  &
  old_local_size , old_global_size                                &
, old_global_north , old_global_east_in                           &
, old_global_south , old_global_west                              &
, old_gridtype_code , old_scatter_pe                              &
, old_halo_type,old_current_decomp_type

INTEGER, SAVE ::                                                  &
! variables defining send and receive maps to be passed to
! GC_RALL_TO_ALL, defining the data transposition
  receive_map(7,2)                                                &
, n_sends,n_recvs  ! number of sends and receives
INTEGER, SAVE, ALLOCATABLE ::   send_map(:,:)


LOGICAL :: wrap          ! if the subdomain wraps over 0 degree meridion
LOGICAL :: wrap_special  ! if there is a wrap around, which starts and
                         ! ends on the same processor
LOGICAL :: zonal_data    ! if this is a zonal data grid
LOGICAL :: fullfield     ! if this is a full field - NOT a subarea


! Set all the old_* variables to a number indicating they've
! not been used yet

DATA                                                              &
  old_local_size , old_global_size                                &
, old_global_north , old_global_east_in                           &
, old_global_south , old_global_west                              &
, old_gridtype_code , old_scatter_pe                              &
, old_halo_type,old_current_decomp_type                           &
  / -1,-1,-1,-1,-1,-1,-1,-1,-1,-1 /

! Functions

INTEGER :: get_fld_type

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STASH_SCATTER_FIELD'
! ------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(send_map)) &
  ALLOCATE(send_map(7,2*nproc_max))

! DEPENDS ON: get_fld_type
fld_type=get_fld_type(gridtype_code)

l_vec=((fld_type  ==  fld_type_u) .OR.                            &
       (fld_type  ==  fld_type_v))

! See if there is wrap around over meridion, and if so make
! sure that GLOBAL_EAST is > glsize(1)

global_east=global_east_in
IF (fld_type  ==  fld_type_unknown) THEN
  WRITE(umMessage,*)                                                      &
  'STASH_SCATTER_FIELD cannot process field with ppx gridtype ',  &
  gridtype_code
  CALL umPrint(umMessage,src='stash_scatter_field')
  icode=1
  cmessage='STASH_SCATTER_FIELD : Incompatible GRIDTYPE code'
  GO TO 9999
END IF

IF (global_east  >   glsize(1,fld_type)) THEN
  wrap=.TRUE.
ELSE IF (global_east  <   global_west) THEN
  wrap=.TRUE.
  global_east=global_east_in+glsize(1,fld_type)
ELSE
  wrap=.FALSE.
END IF


IF (gridtype_code  ==  ppx_atm_tzonal) THEN
  zonal_data=.TRUE.
  global_x=1
ELSE

  ! This is a normal field

  zonal_data=.FALSE.
  global_x=glsize(1,fld_type)

END IF

global_y=glsize(2,fld_type)

! Set up logical indicating if this is a full field, or just
! a subdomain

IF (zonal_data) THEN

  fullfield= ( ( global_north  ==  global_y) .AND.                &
             ( global_south  ==  1))

ELSE

  fullfield = (( global_west  ==  1) .AND.                        &
               ( global_east  ==  global_x) .AND.                 &
               ( global_north  ==  global_y) .AND.                &
               ( global_south  ==  1))

END IF

! Dealing with fields not in model grid

IF ((global_x == 0) .OR. (global_x == imdi)) THEN
  DO level=1,levels
    DO i=1,global_size
      local_field(i,level)=global_field(i,level)
    END DO
  END DO
ELSE

  ! If this is a fullfield, we can simply use the standard
  ! SCATTER_FIELD routine

  IF (fullfield) THEN

    IF (zonal_data) THEN

      ! DEPENDS ON: scatter_zonal_field
      CALL scatter_zonal_field( local_field,global_field,           &
                                lasize(2,fld_type,halo_type),       &
                                global_y,                           &
                                levels,gridtype_code,fld_type,      &
                                halo_type,                          &
                                scatter_pe)

    ELSE

      DO level=1,levels

        ! DEPENDS ON: scatter_field
        CALL scatter_field( local_field(1,level) ,                  &
                            global_field(1,level),                  &
                            lasize(1,fld_type,halo_type),           &
                            lasize(2,fld_type,halo_type),           &
                            global_x,global_y,                      &
                            fld_type,halo_type,                     &
                            scatter_pe,gc_all_proc_group )

      END DO
    END IF
  ELSE

    ! for subdomains, life is not so easy - we must explicitly
    ! calculate our own send and receive maps, and use GCG_RALLTOALLE
    ! to shift the data around.

    ! If the same arguments are used as were used in the last call
    ! to this routine, we can just use the previously calculated
    ! send and receive maps, otherwise we need to calculate new maps

    IF (.NOT. (                                                     &
      (local_size  ==  old_local_size) .AND.                        &
      (global_size  ==  old_global_size) .AND.                      &
      (global_north  ==  old_global_north) .AND.                    &
      (global_east_in  ==  old_global_east_in) .AND.                &
      (global_south  ==  old_global_south) .AND.                    &
      (global_west  ==  old_global_west) .AND.                      &
      (gridtype_code  ==  old_gridtype_code) .AND.                  &
      (halo_type  ==  old_halo_type) .AND.                          &
      (scatter_pe  ==  old_scatter_pe) .AND.                        &
      (current_decomp_type  ==  old_current_decomp_type ))) THEN

      old_local_size=local_size
      old_global_size=global_size
      old_global_north=global_north
      old_global_east_in=global_east_in
      old_global_south=global_south
      old_global_west=global_west
      old_gridtype_code=gridtype_code
      old_halo_type=halo_type
      old_scatter_pe=scatter_pe
      old_current_decomp_type=current_decomp_type

      ! Find out what the boundaries of the subdomain area

      ! DEPENDS ON: global_to_local_rc
      CALL global_to_local_rc(gridtype_code,halo_type,              &
                              global_west,global_north,             &
                              proc_topleft_x,proc_topleft_y,        &
                              dummy1,dummy2)
      ! DEPENDS ON: global_to_local_rc
      CALL global_to_local_rc(gridtype_code,halo_type,              &
                              global_east,global_south,             &
                              proc_botright_x,proc_botright_y,      &
                              dummy1,dummy2)

      ! Ensure that the processor x co-ords are such that the botright_x is
      ! always greater than (or equal to) top_left_x.
      IF (wrap) proc_botright_x=gridsize(1)+proc_botright_x

      ! wrap_special is set to true if there is a wrap around which starts
      ! and ends on the same processor. This case requires extra work as
      ! the processor in question
      IF (wrap .AND. (proc_topleft_x+gridsize(1)  ==                &
                      proc_botright_x)) THEN
        wrap_special=.TRUE.
      ELSE
        wrap_special=.FALSE.
      END IF

      n_sends=0
      n_recvs=0

      DO procy=proc_botright_y,proc_topleft_y
        DO procx=proc_topleft_x,proc_botright_x

          eff_procx=MOD(procx,gridsize(1))
          procid=eff_procx+procy*gridsize(1)

          ! DEPENDS ON: global_to_local_subdomain
          CALL global_to_local_subdomain(                           &
            include_halos_ew, include_halos_ns,                     &
            gridtype_code,halo_type,procid,                         &
            global_south,global_east,                               &
            global_north,global_west,                               &
            local_ystart,local_xend,                                &
            local_yend  ,local_xstart)

          ! Calculate the shape of the arrays, and where to start sending/
          ! receiving data, and how many rows to send

          local_start_row=1
          nrows_to_send=local_yend-local_ystart+1

          global_start_row=g_datastart(2,procid)+local_ystart -     &
                           halosize(2,halo_type) -                  &
                           global_south
          global_row_length=global_east-global_west+1

          ! Calculate the following variables:
          ! local_row_length : the X dimension size of the local array
          ! local_send_offx  : the offset into each row to start sending from
          ! sendsize_x       : the number of points on each row to send
          ! The calculation of these numbers is different for processors
          ! at the start and end of a wrap_special case

          IF (wrap_special .AND. procx  ==  proc_topleft_x) THEN
            local_row_length=g_blsize(1,fld_type,procid) +          &
                             local_xend - local_xstart + 1
            local_start_col=1
            sendsize_x=g_lasize(1,fld_type,halo_type,procid) -      &
                       local_xstart
            global_start_col=1

          ELSE IF (wrap_special .AND. procx  ==  proc_botright_x)    &
          THEN
            local_row_length=g_blsize(1,fld_type,procid) +          &
                             local_xend - local_xstart + 1
            local_start_col=local_row_length - local_xend +         &
                            halosize(1,halo_type) + 1
            sendsize_x=local_xend - halosize(1,halo_type)
            global_start_col=global_row_length-sendsize_x+1

          ELSE
            local_row_length=local_xend-local_xstart+1
            local_start_col=1
            sendsize_x=local_xend-local_xstart+1
            global_start_col=local_xstart -                         &
                             (halosize(1,halo_type) + 1 ) +         &
                             g_datastart(1,procid)-global_west+1
          END IF

          IF (global_start_col  <   0) THEN
            ! Wrapped around field, but this processor is not start or end
            ! processor
            global_start_col=global_start_col+glsize(1,fld_type)
          END IF

          ! Now we can set up the send and receive map entries for the data on
          ! this processor

          IF (mype  ==  procid) THEN  ! I need to receive some data

            n_recvs=n_recvs+1

            receive_map(r_source_pe,n_recvs) = scatter_pe
            receive_map(r_base_address_in_recv_array,n_recvs) =     &
              (local_start_row-1)*local_row_length +                &
              local_start_col
            receive_map(r_number_of_elements_in_item,n_recvs) =     &
              nrows_to_send
            receive_map(r_stride_in_recv_array,n_recvs) =           &
              local_row_length
            receive_map(r_element_length,n_recvs) = sendsize_x
            receive_map(r_base_address_in_send_array,n_recvs) =     &
              (global_start_row-1)*global_row_length +              &
              global_start_col
            receive_map(r_stride_in_send_array,n_recvs) =           &
              global_row_length

          END IF ! if I'm receiving data

          IF (mype  ==  scatter_pe) THEN ! I need to send data

            n_sends=n_sends+1

            send_map(s_destination_pe,n_sends) = procid
            send_map(s_base_address_in_send_array,n_sends) =        &
              (global_start_row-1)*global_row_length +              &
              global_start_col
            send_map(s_number_of_elements_in_item,n_sends) =        &
              nrows_to_send
            send_map(s_stride_in_send_array,n_sends) =              &
              global_row_length
            send_map(s_element_length,n_sends) = sendsize_x
            send_map(s_base_address_in_recv_array,n_sends) =        &
              (local_start_row-1)*local_row_length +                &
              local_start_col
            send_map(s_stride_in_recv_array,n_sends) =              &
              local_row_length

          END IF ! if I'm sending data

        END DO ! procx : loop along processor row

      END DO ! procy : loop down processor column

    END IF ! if I need to recalculate my send/receive maps

    ! Send / receive the data using GCG_RALLTOALLE

    DO level=1,levels

      flag=0  ! This is currently ignored at GCG v1.1
      info=gc_none

      CALL gcg_ralltoalle(                                          &
        global_field(1,level)  ,                                    &
        send_map    , n_sends  ,global_size  ,                      &
        local_field(1,level) ,                                      &
        receive_map , n_recvs , local_size ,                        &
        gc_all_proc_group , flag, info)

    END DO

  END IF ! if this is a full or extracted field

END IF

9999 CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stash_scatter_field

END MODULE stash_scatter_field_mod

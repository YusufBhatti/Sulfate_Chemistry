! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Gathers STASHed data from many processors to one processor

MODULE stash_gather_field_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'STASH_GATHER_FIELD_MOD'
  
CONTAINS

SUBROUTINE stash_gather_field (                                   &
   local_field , global_field ,                                    &
   local_size, global_size, levels,                                &
   global_north , global_east_in , global_south , global_west,     &
   gridtype_code ,halo_type,                                       &
   gather_pe,                                                      &
   data_extracted,                                                 &
   packing, im_ident, packing_type,                          &
   num_out,                                                        &
   comp_accrcy, loc_rmdi,                                          &
   icode, cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE gcom_mod
USE UM_ParVars
USE UM_ParCore,   ONLY: mype, nproc, nproc_max
USE UM_ParParams, ONLY: halo_type_no_halo
USE Field_Types,  ONLY: fld_type_p, fld_type_u
USE mpp_conf_mod, ONLY: include_halos_ew, include_halos_ns
USE cppxref_mod, ONLY:                                            &
    ppx_atm_tzonal, ppx_atm_uzonal
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE missing_data_mod, ONLY: imdi 
USE errormessagelength_mod, ONLY: errormessagelength
USE mpp_conf_mod, ONLY: gather_sync

IMPLICIT NONE


! Description:
! Takes a decomposed, STASH processed field and gathers
! it to a single processor, ready for I/O,

! Method:
! See in-line documentation

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Subroutine arguments:


INTEGER, INTENT(IN) ::                                            &
  local_size                                                      &
                  ! IN: size of level of LOCAL_FIELD
, global_size                                                     &
                  ! IN: size of level of GLOBAL_FIELD
, levels                                                          &
                  ! IN: number of levels
, global_north                                                    &
                  ! IN: specification of subdomain boundaries
, global_east_in                                                  &
                  ! IN: ""
, global_south                                                    &
                  ! IN: ""
, global_west                                                     &
                  ! IN: ""
, gridtype_code                                                   &
                  ! IN: indicates the type of grid output
, halo_type                                                       &
                  ! IN: type of halo on this field
, gather_pe       ! IN: the PE to gather the global field to

INTEGER, INTENT(OUT) ::                                           &
  icode           ! OUT: return code, 0=OK

! Optional Arguments to handle the WGDOS packing if necessary

LOGICAL, INTENT(IN), OPTIONAL ::                                  &
  packing
                  ! IN: Set .true. if packing of the input
                  !     field is to be packed!

INTEGER, INTENT(IN), OPTIONAL ::                                  &
  im_ident        ! IN: Internal model identifier

INTEGER, INTENT(INOUT), OPTIONAL ::                               &
  packing_type    ! IN/OUT: This flag is zero on input,
                  !         then stash packing is selected,
                  !         and the routine returns the
                  !         packing flag.
                  !
                  !         If the variable is set to 1 on input
                  !         then 32-bit packing for dumpfiles
                  !         is selected

INTEGER, INTENT(OUT), OPTIONAL ::                                 &
  num_out         ! OUT: Number of 32-bit IBM words in the Packed
                  !      field for WDGOS packing

INTEGER, INTENT(IN), OPTIONAL ::                                  &
  comp_accrcy     ! IN: Packing Accuracy in Power of 2

REAL, INTENT(IN), OPTIONAL ::                                     &
  loc_rmdi        ! IN: Missing data indicator

! Remaining Non-Optional Arguments

LOGICAL, INTENT(IN) ::                                            &
  data_extracted  ! IN: TRUE if the data in LOCAL_FIELD has
                  !     already been extracted, or FALSE if
                  !     the extraction must be done here.

REAL, INTENT(IN) ::                                               &
  local_field(local_size,levels)
                  ! IN : local data

REAL, INTENT(OUT) ::                                              &
  global_field(global_size,levels)
                  ! OUT : (PE GATHER_PE only) - full gathered
                  !       field

CHARACTER(LEN=errormessagelength), INTENT(INOUT) ::               &
  cmessage        ! INOUT: Error message if ICODE  /=  0

! Local variables

INTEGER    ::                                                     &
  global_east                                                     &
                ! copy of GLOBAL_EAST_IN with wrap around s.t.
!                     ! GLOBAL_EAST > GLOBAL_ROW_LEN
, global_x                                                        &
                ! size of global data EW
, global_y                                                        &
                ! size of global data NS
, fld_type                                                        &
                ! indicates if field is on P or U grid
, level                                                           &
                ! loop index for loop over levels
, i                                                               &
                ! loop counter
, proc_topleft_x,proc_topleft_y                                   &
                                  ! processors at corners of
, proc_botright_x,proc_botright_y                                 &
                                  ! the subarea
, dummy1,dummy2                                                   &
                ! ignored return arguments
, procx,procy                                                     &
                ! loop indexes for loops over processors
, eff_procx                                                       &
                ! real x co-ord of processor column procx
, procid                                                          &
                ! processor id of (procx,procy)
, local_xstart,local_xend                                         &
                           ! boundaries of subdomain for
, local_ystart,local_yend                                         &
                           ! processor procid
, local_start_row                                                 &
                    ! first row to send from procid
, local_start_col                                                 &
                    ! first column to send from procid
, sendsize_x                                                      &
                    ! number of points on each row to send
, nrows_to_send                                                   &
                    ! number of rows to send from procid
, local_row_length                                                &
                    ! size of sending array EW
, global_start_row                                                &
                    ! first row to receive at on GATHER_PE
, global_start_col                                                &
                    ! first col. to recieve on GATHER_PE
, global_row_length                                               &
                    ! size of receiving array EW
, flag,info         ! GCOM arguments

! Copies of arguments / variables used to decide if we can use the
! send/receive maps used in the last call

! Save all the variables that may be used in the next call

INTEGER, SAVE ::                                                  &
  old_local_size , old_global_size                                &
, old_global_north , old_global_east_in                           &
, old_global_south , old_global_west                              &
, old_gridtype_code , old_gather_pe                               &
, old_current_decomp_type, old_halo_type

! variables defining send and receive maps to be passed to
! GC_RALL_TO_ALL, defining the data transposition
INTEGER, SAVE :: send_map(7,2), n_sends, n_recvs
INTEGER, SAVE, ALLOCATABLE :: receive_map(:,:)


LOGICAL  ::                                                       &
  wrap                                                            &
               ! if the subdomain wraps over 0 degree meridion
, wrap_special                                                    &
               ! if there is a wrap around, which starts and
!                      ends on the same processor
, zonal_data                                                      &
               ! if this is a zonal data grid
, fullfield                                                       &
               ! if this is a full field - NOT a subarea
, l_packing    ! if packing of data is required



! Set all the old_* variables to a number indicating they've
! not been used yet

DATA                                                              &
  old_local_size , old_global_size                                &
, old_global_north , old_global_east_in                           &
, old_global_south , old_global_west                              &
, old_gridtype_code , old_gather_pe                               &
, old_current_decomp_type, old_halo_type                          &
  / -1,-1,-1,-1,-1,-1,-1,-1,-1,-1 /

! Functions

INTEGER :: get_fld_type

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STASH_GATHER_FIELD'
! ------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(receive_map)) &
  ALLOCATE(receive_map(7,2*nproc_max))

icode=0
IF (PRESENT(packing)) THEN
  l_packing = packing
ELSE
  l_packing = .FALSE.
END IF

! DEPENDS ON: get_fld_type
fld_type=get_fld_type(gridtype_code)
! See if there is wrap around over meridion, and if so make
! sure that GLOBAL_EAST is > glsize(1)

global_east=global_east_in
IF (global_east  >   glsize(1,fld_type)) THEN
  wrap=.TRUE.
ELSE IF (global_east  <   global_west) THEN
  wrap=.TRUE.
  global_east=global_east_in+glsize(1,fld_type)
ELSE
  wrap=.FALSE.
END IF

IF ((gridtype_code  ==  ppx_atm_tzonal) .OR.   & ! Atmos T zonal
   ( gridtype_code  ==  ppx_atm_uzonal) )  THEN ! Atmos U zonal
  ! This is a zonal field
  zonal_data=.TRUE.
  global_x=1

  IF (gridtype_code  ==  ppx_atm_tzonal) THEN ! Atmos T zonal
    fld_type=fld_type_p
  ELSE
    fld_type=fld_type_u
  END IF
ELSE

  ! This is a normal field

  zonal_data=.FALSE.
  global_x=glsize(1,fld_type)
END IF

global_y=glsize(2,fld_type)

! Set up logical indicating if this is a full field, or just
! a subdomain

IF (zonal_data) THEN

  fullfield= ( ( global_south  ==  1) .AND.                       &
             ( global_north  ==  global_y))

ELSE

  fullfield = (( global_west  ==  1) .AND.                        &
               ( global_east  ==  global_x) .AND.                 &
               ( global_south  ==  1) .AND.                       &
               ( global_north  ==  global_y))

END IF

! Dealing with fields not in model grid

IF ((global_x == 0) .OR. (global_x == imdi)) THEN
  WRITE(umMessage,*)'local_size=',local_size
  CALL umPrint(umMessage,src='stash_gather_field')
  WRITE(umMessage,*)'global_size=',global_size
  CALL umPrint(umMessage,src='stash_gather_field')
  DO level=1,levels
    DO i=1,global_size
      global_field(i,level)=local_field(i,level)
    END DO
  END DO
ELSE

  ! If this a fullfield, we can simply use the standard
  ! GATHER_FIELD routine

  IF (fullfield) THEN

    IF (zonal_data) THEN

      ! DEPENDS ON: gather_zonal_field
      CALL gather_zonal_field( local_field,global_field,            &
                            lasize(2,fld_type,halo_type),global_y,  &
                              levels,gridtype_code,fld_type,        &
                              halo_type,gather_pe)

    ELSE

      DO level=1,levels

        IF (l_packing) THEN
          ! DEPENDS ON: gather_pack_field
          CALL gather_pack_field(                                   &
                          local_field(1,level),                     &
                          global_field(1,level),                    &
                          lasize(1,fld_type,halo_type),             &
                          lasize(2,fld_type,halo_type),             &
                          global_x,global_y,                        &
                          fld_type,halo_type,                       &
                          gather_pe,gc_all_proc_group,              &
                          packing, im_ident, packing_type,    &
                          num_out,                                  &
                          comp_accrcy, loc_rmdi)
        ELSE
          ! DEPENDS ON: gather_field
          CALL gather_field( local_field(1,level),                  &
                             global_field(1,level),                 &
                             lasize(1,fld_type,halo_type),          &
                             lasize(2,fld_type,halo_type),          &
                             global_x, global_y,                    &
                             fld_type, halo_type,                   &
                             gather_pe, gc_all_proc_group)
        END IF

        IF (gather_sync) THEN
          ! A synchronisation can help with message congestion
          CALL gc_gsync(nproc,icode)
        END IF

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
      (gather_pe  ==  old_gather_pe) .AND.                          &
      (halo_type  ==  old_halo_type) .AND.                          &
      (current_decomp_type  ==  old_current_decomp_type ))) THEN

      old_local_size=local_size
      old_global_size=global_size
      old_global_north=global_north
      old_global_east_in=global_east_in
      old_global_south=global_south
      old_global_west=global_west
      old_gridtype_code=gridtype_code
      old_gather_pe=gather_pe
      old_current_decomp_type=current_decomp_type
      old_halo_type=halo_type

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

          IF (data_extracted) THEN
            local_start_row=1
          ELSE
            local_start_row=local_ystart
          END IF
          nrows_to_send=local_yend-local_ystart+1

          global_start_row=g_datastart(2,procid)+local_ystart-       &
                            halosize(2,halo_type) - global_south
          global_row_length=global_east-global_west+1

          ! Calculate the following variables:
          ! local_row_length : the X dimension size of the local array
          ! local_send_offx  : the offset into each row to start sending from
          ! sendsize_x       : the number of points on each row to send
          ! The calculation of these numbers is different for processors
          ! at the start and end of a wrap_special case
          ! Note that when DATA_EXTRACTED is true, then local_field has no
          ! halos.

          IF (wrap_special .AND. procx  ==  proc_topleft_x) THEN
            IF (data_extracted) THEN
              local_row_length = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
                                 + local_xend - local_xstart + 1
              sendsize_x       = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
                                 - local_xstart + 1
              local_start_col  = 1
            ELSE
              local_row_length = g_lasize(1,fld_type,halo_type,procid)
              sendsize_x       = g_lasize(1,fld_type,halo_type,procid)          &
                                 - local_xstart
              local_start_col  = local_xstart

            END IF
            global_start_col=1

          ELSE IF (wrap_special .AND. procx  ==  proc_botright_x)    &
          THEN
            IF (data_extracted) THEN
              local_row_length = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
                                 + local_xend - local_xstart + 1
              local_start_col  = local_row_length - local_xend + 1
              sendsize_x       = local_xend
            ELSE
              local_row_length = g_lasize(1,fld_type,halo_type,procid)
              local_start_col  = offx + 1
              sendsize_x       = local_xend - offx
            END IF
            global_start_col=global_row_length-sendsize_x+1

          ELSE
            IF (data_extracted) THEN
              local_row_length=local_xend-local_xstart+1
              local_start_col=1
            ELSE
              local_row_length=g_lasize(1,fld_type,halo_type,procid)
              local_start_col=local_xstart
            END IF
            sendsize_x=local_xend-local_xstart+1
            global_start_col=local_xstart-halosize(1,halo_type)+    &
                             g_datastart(1,procid)-global_west
          END IF

          IF (global_start_col  <   0) THEN
            ! Wrapped around field, but this processor is not start or end
            ! processor
            global_start_col=global_start_col+glsize(1,fld_type)
          END IF

          ! Now we can set up the send and receive map entries for the data on
          ! this processor

          IF (mype  ==  procid) THEN  ! I need to send some data


            n_sends=n_sends+1

            send_map(s_destination_pe,n_sends) = gather_pe
            send_map(s_base_address_in_send_array,n_sends) =        &
              (local_start_row-1)*local_row_length +                &
              local_start_col
            send_map(s_number_of_elements_in_item,n_sends) =        &
              nrows_to_send
            send_map(s_stride_in_send_array,n_sends) =              &
              local_row_length
            send_map(s_element_length,n_sends) = sendsize_x
            send_map(s_base_address_in_recv_array,n_sends) =        &
              (global_start_row-1)*global_row_length +              &
              global_start_col
            send_map(s_stride_in_recv_array,n_sends) =              &
              global_row_length

          END IF ! if I'm sending data

          IF (mype  ==  gather_pe) THEN  ! I need to receive data

            n_recvs=n_recvs+1

            receive_map(r_source_pe,n_recvs) = procid
            receive_map(r_base_address_in_recv_array,n_recvs) =     &
              (global_start_row-1)*global_row_length +              &
              global_start_col
            receive_map(r_number_of_elements_in_item,n_recvs) =     &
              nrows_to_send
            receive_map(r_stride_in_recv_array,n_recvs) =           &
              global_row_length
            receive_map(r_element_length,n_recvs) = sendsize_x
            receive_map(r_base_address_in_send_array,n_recvs) =     &
              (local_start_row-1)*local_row_length +                &
              local_start_col
            receive_map(r_stride_in_send_array,n_recvs) =           &
              local_row_length

          END IF ! if I'm receiving data

        END DO ! procx : loop along processor row

      END DO ! procy : loop down processor column

    END IF ! if I need to recalculate my send/receive maps

    ! Send / receive the data using GCG_RALLTOALLE


    DO level=1,levels

      flag=0  ! This is currently ignored at GCG v1.1
      info=gc_none

      CALL gcg_ralltoalle(                                          &
        local_field(1,level)  ,                                     &
        send_map    , n_sends  ,local_size  ,                       &
        global_field(1,level) ,                                     &
        receive_map , n_recvs , global_size ,                       &
        gc_all_proc_group , flag, info)

    END DO

  END IF ! if this is a full or extracted field

END IF


9999 CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stash_gather_field

END MODULE stash_gather_field_mod

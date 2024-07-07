! -- begin halo_exchange_mod_swap_bounds_RB_hub.h --
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

LOGICAL, INTENT(IN) :: red        ! If .TRUE. the red elements need to
                                  ! be swapped

REAL(KIND=field_kind), INTENT(INOUT) ::                                    &
  field(0:row_length+1,0:rows+1,levels)  


! Local variables

LOGICAL  :: l_vector
LOGICAL  :: do_east, do_west, do_south, do_north, do_corners

REAL(KIND=jprb)               :: zhook_handle


! Check if this routine can perform the swap bounds, 
! if not call the base routine.

IF (nproc_x==1 .OR. nproc_y==1 .OR. &
  (sb_model_type /= mt_global .AND. sb_model_type /= mt_lam) .OR. &
  (sb_model_type == mt_global .AND. MOD(glsize(1,fld_type_p),2)/=0)) THEN

  l_vector   = .FALSE.
  do_east    = .TRUE.  
  do_west    = .TRUE. 
  do_south   = .TRUE. 
  do_north   = .TRUE. 
  do_corners = .FALSE.

  CALL swap_bounds(                     &
    field, row_length, rows, levels,    & 
    1, 1,                               & 
    fld_type_p, l_vector,               & 
    do_east, do_west,                   & 
    do_south, do_north,                 & 
    do_corners                          & 
    )
  RETURN
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Only one choice.
CALL swap_bounds_mpi_rb(field, row_length, rows, levels, red)


IF (sb_model_type == mt_lam) THEN
  CALL fill_external_halos(field, row_length, rows, levels, 1, 1)
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

! -- end  halo_exchange_mod_swap_bounds_RB_hub.h --

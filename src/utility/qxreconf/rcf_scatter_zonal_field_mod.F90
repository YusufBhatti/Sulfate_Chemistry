! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! + Scatters zonal field from one processor to many processors

MODULE Rcf_Scatter_Zonal_Field_Mod

IMPLICIT NONE
!  Subroutine Rcf_Scatter_Zonal - scatters a zonal field
!
! Description:
! Takes a zonal field on a single processor, and decomposes it over
! many processors.
!
! Method:
!   Use gcg_ralltoalle to scatter data to processes by utilising
!   the decomposition information (sizes and data starts) previously
!   computed and stored in UM_parvars, etc.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SCATTER_ZONAL_FIELD_MOD'

CONTAINS

SUBROUTINE Rcf_Scatter_Zonal_Field (  local_field , global_field , &
                                      local_size  , global_size  , &
                                      levels, grid_type , &
                                      scatter_pe)

USE UM_ParVars     !  Lots used.
USE UM_ParCore, ONLY: &
    mype,             &
    nproc,            &
    nproc_max
USE UM_ParParams, ONLY: &
    halo_type_single
USE Field_Types, ONLY: &
    fld_type_p,        &
    fld_type_u
USE gcom_mod
USE cppxref_mod, ONLY: &
    ppx_atm_tzonal
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)   :: Local_Size    ! size of local_field level
INTEGER, INTENT(IN)   :: Global_Size   ! size of global_filed level
INTEGER, INTENT(IN)   :: Levels        ! number of levels in field
INTEGER, INTENT(IN)   :: Grid_Type     ! grid type of field
INTEGER, INTENT(IN)   :: Scatter_PE    ! PE on which Global_field lives

                                       ! Field to scatter
REAL, INTENT(IN)      :: Global_Field( Global_Size, Levels )
                                       ! Local part of field
REAL, INTENT(OUT)     :: Local_Field( Local_Size, Levels )

! Local variables
INTEGER               :: fld_type            ! P or U
INTEGER               :: info                ! return from GCOM
INTEGER, ALLOCATABLE, SAVE  :: send_map(:,:)       ! send map
INTEGER               :: receive_map(7,1)    ! receive map
INTEGER               :: flag                ! dummy arg for GCOM
INTEGER               :: n_mess_to_send      ! number of messages
INTEGER               :: k                   ! loop counter
INTEGER               :: iproc               ! loop counter

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_SCATTER_ZONAL_FIELD'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!====================================================================

IF (.NOT. ALLOCATED(send_map)) &
    ALLOCATE(send_map(7,nproc_max))

IF (grid_type  ==  ppx_atm_tzonal) THEN
  fld_type=fld_type_p
ELSE
  fld_type=fld_type_u
END IF

!--------------------------------------------------------------------

n_mess_to_send=0

IF (mype  ==  scatter_pe) THEN
  DO iproc=0,nproc-1
    send_map(s_destination_pe,iproc+1) = iproc
    send_map(s_base_address_in_send_array,iproc+1) = &
             g_datastart(2,iproc)
    send_map(s_number_of_elements_in_item,iproc+1) = 1
    send_map(s_stride_in_send_array,iproc+1) = 0
    IF (fld_type  ==  fld_type_p) THEN
      send_map(s_element_length,iproc+1) = g_blsizep(2,iproc)
    ELSE
      send_map(s_element_length,iproc+1) = g_blsizeu(2,iproc)
    END IF
    send_map(s_base_address_in_recv_array,iproc+1) = Offy+1
    send_map(s_stride_in_recv_array,iproc+1) = 0
  END DO
  n_mess_to_send=nproc
END IF

receive_map(r_source_pe,1) = scatter_pe
receive_map(r_base_address_in_recv_array,1) = Offy+1
receive_map(r_number_of_elements_in_item,1) = 1
receive_map(r_stride_in_recv_array,1) = 0
IF (fld_type  ==  fld_type_p) THEN
  receive_map(r_element_length,1) = blsizep(2)
ELSE
  receive_map(r_element_length,1) = blsizeu(2)
END IF
receive_map(r_base_address_in_send_array,1) = datastart(2)
receive_map(r_stride_in_send_array,1) = 0

DO k=1,levels

  info=gc_none
  flag=gc_none

  IF (fld_type  ==  fld_type_p) THEN
    CALL gcg_ralltoalle( global_field(1,k),    send_map,        &
                         n_mess_to_send,       glsizep(2),       &
                         local_field(1,k),     receive_map,     &
                         1,                                     &
                         lasize(2,fld_type_p,halo_type_single), &
                         gc_all_proc_group,    flag,    info)
  ELSE
    CALL gcg_ralltoalle( global_field(1,k),     send_map,       &
                         n_mess_to_send,        glsizep(2)-1,    &
                         local_field(1,k),      receive_map,    &
                         1,                                     &
                         lasize(2,fld_type_p,halo_type_single), &
                         gc_all_proc_group,     flag,    info)
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Scatter_Zonal_Field
END MODULE Rcf_Scatter_Zonal_Field_Mod



! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! + Gathers a zonal field from many processors to one processor

MODULE Rcf_Gather_Zonal_Field_Mod

!  Subroutine Rcf_GAther_Zonal_Field -  Gathers a field onto 1 pe
!
! Description:
! Takes a zonal field decomposed on many processors and gathers
! it onto a single specified PE.
!
! Method:
!  Calculates send and receive maps to use with GCG_RALLETOALLE
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GATHER_ZONAL_FIELD_MOD'

CONTAINS

SUBROUTINE Rcf_Gather_Zonal_Field ( local_field, global_field,  &
                                    local_size,  global_size,   &
                                    levels,      grid_type,     &
                                    gather_pe)

USE UM_ParVars    ! Most of this used
USE gcom_mod
USE cppxref_mod, ONLY: ppx_atm_tzonal, ppx_atm_ozone

USE UM_ParCore, ONLY:  &
    mype,      nproc,  &
    nproc_max

USE UM_ParParams, ONLY: &
    halo_type_single

USE Field_Types, ONLY: &
    fld_type_p,  fld_type_u

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)   :: Local_Size
INTEGER, INTENT(IN)   :: Global_Size
INTEGER, INTENT(IN)   :: Levels
INTEGER, INTENT(IN)   :: Grid_Type
INTEGER, INTENT(IN)   :: Gather_PE
REAL, INTENT(IN)      :: Local_Field ( Local_Size, Levels )
REAL, INTENT(OUT)     :: Global_Field( Global_Size, Levels )

! Local variables
INTEGER               :: k                 ! looper
INTEGER               :: flag              ! Gcom argument
INTEGER               :: fld_type          ! P or U field
INTEGER               :: info              ! Gcom return code
INTEGER               :: iproc             ! loop counter
INTEGER               :: send_map(7,1)
INTEGER, ALLOCATABLE,SAVE  :: receive_map(:,:)
INTEGER               :: n_mess_to_send
INTEGER               :: n_mess_to_receive
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_GATHER_ZONAL_FIELD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!====================================================================

IF (.NOT. ALLOCATED(receive_map)) &
    ALLOCATE(receive_map(7,nproc_max))

send_map(:,:)    = 0
receive_map(:,:) = 0

! Note no V zonal grids....
IF ((grid_type  ==  ppx_atm_tzonal) .OR. &
    (grid_type  ==  ppx_atm_ozone)) THEN
  fld_type=fld_type_p
ELSE
  fld_type=fld_type_u
END IF

!--------------------------------------------------------------------

n_mess_to_receive=0

IF (mype  ==  gather_pe) THEN
  DO iproc=0,nproc-1
    IF (g_gridpos(1,iproc)  ==  0) THEN
      !           Only one processor per LPG row needs to send the data
      !           as it will be the same for each processor along the
      !           row.
      receive_map(r_source_pe,n_mess_to_receive+1) = iproc
      receive_map(r_base_address_in_recv_array, n_mess_to_receive+1) = &
                  g_datastart(2,iproc)
      receive_map(r_number_of_elements_in_item, n_mess_to_receive+1) = 1
      receive_map(r_stride_in_recv_array, n_mess_to_receive+1) = 0
      IF (fld_type  ==  fld_type_p) THEN
        receive_map(r_element_length,n_mess_to_receive+1) = &
                    g_blsizep(2,iproc)
      ELSE
        receive_map(r_element_length,n_mess_to_receive+1) = &
                    g_blsizeu(2,iproc)
      END IF
      receive_map(r_base_address_in_send_array, &
                  n_mess_to_receive+1) = Offy+1
      receive_map(r_stride_in_send_array, n_mess_to_receive+1) = 0
      n_mess_to_receive=n_mess_to_receive+1
    END IF
  END DO
END IF

n_mess_to_send=0
IF (atwest) THEN ! only processors at the left of the LPG will
                 ! send anything
  send_map(s_destination_pe,1) = gather_pe
  send_map(s_base_address_in_send_array,1) = Offy+1
  send_map(s_number_of_elements_in_item,1) = 1
  send_map(s_stride_in_send_array,1) = 0
  IF (fld_type  ==  fld_type_p) THEN
    send_map(s_element_length,1) = blsizep(2)
  ELSE
    send_map(s_element_length,1) = blsizeu(2)
  END IF
  send_map(s_base_address_in_recv_array,1) = datastart(2)
  send_map(s_stride_in_recv_array,1) = 0

  n_mess_to_send=1
END IF


DO k=1,levels

  info=gc_none
  flag=gc_none

  IF (fld_type  ==  fld_type_p) THEN
    CALL gcg_ralltoalle(                                 &
        local_field(1,k),send_map,n_mess_to_send,        &
        lasize(2,fld_type_p,halo_type_single),          &
        global_field(1,k),receive_map,n_mess_to_receive, &
        glsizep(2),gc_all_proc_group,flag,info)
  ELSE
    CALL gcg_ralltoalle(                                 &
        local_field(1,k),send_map,n_mess_to_send,        &
        lasize(2,fld_type_p,halo_type_single),          &
        global_field(1,k),receive_map,n_mess_to_receive, &
        glsizep(2)-1,gc_all_proc_group,flag,info)
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Gather_Zonal_Field
END MODULE Rcf_Gather_Zonal_Field_Mod



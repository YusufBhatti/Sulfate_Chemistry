! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Gathers a field from many processors to one processor

MODULE Rcf_Gather_Field_Mod

!  Subroutine Rcf_Gather_field_real  - gather real fields
!  Subroutine Rcf_Gather_field_log   - gathers logical fields
!
! Description:
!  Takes a model field that has been decomposed over a group of
!  processors, and gathers the data together so that one processor
!  contains the entire global field.
!
! Method:
!  A send and receive map is constructed which instructs the GCOM
!  permute operation to do a gather from all processors in the
!  group to the GATHER_PE
!--------------------------------------------------------------------
!   ************************************************************
!      This subroutine has an overloaded interface so as to
!      cope cleanly with a number of different data types.
!   ***********************************************************
!--------------------------------------------------------------------
!
! Derived from UM4.5 code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb

INTERFACE Rcf_Gather_Field
MODULE PROCEDURE Rcf_Gather_Field_Log, Rcf_Gather_Field_Real
END INTERFACE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GATHER_FIELD_MOD'

CONTAINS

!--------------------------------------------------------------------
!   Thus subroutine copes with the REAL case
!--------------------------------------------------------------------
SUBROUTINE Rcf_Gather_Field_Real(local_field,    global_field, &
                                 local_row_len,  local_rows,   &
                                 global_row_len, global_rows,  &
                                 gather_pe,      proc_group )

USE Ereport_Mod, ONLY: &
    Ereport

USE UM_ParVars   ! Use a lot of this

USE UM_ParCore, ONLY:  &
    mype,       nproc, &
    nproc_max

USE Field_Types, ONLY:      &
    fld_type_p, fld_type_r, &
    fld_type_u, fld_type_v

USE gcom_mod

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN)    :: Local_Row_Len  ! row length in local field
INTEGER, INTENT(IN)    :: Global_Row_Len ! row length in global field
INTEGER, INTENT(IN)    :: Local_Rows     ! number of rows locally
INTEGER, INTENT(IN)    :: Global_Rows    ! number of rows globally
INTEGER, INTENT(IN)    :: Gather_PE      ! PE on which to gather data
INTEGER, INTENT(IN)    :: Proc_Group     ! group of PEs involved

REAL, INTENT(IN)       :: Local_Field( Local_Row_Len * Local_Rows )
REAL, INTENT(OUT)      :: Global_Field( Global_Row_Len * Global_Rows)

! Local variables
INTEGER                :: Info        ! return code from Gcom
INTEGER                :: Fld_Type    ! P, U or V
INTEGER                :: Flag        ! dummy for Gcom
INTEGER                :: iproc       ! processor number
INTEGER                :: ErrorStatus
CHARACTER (LEN=errormessagelength)     :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_GATHER_FIELD_REAL'

INTEGER, SAVE          :: n_mess_to_rec
INTEGER, SAVE          :: send_map(7,1)
INTEGER, ALLOCATABLE, SAVE :: receive_map(:,:)
INTEGER, SAVE          :: Old_Global_Row_Len       = -1234  ! old values
INTEGER, SAVE          :: Old_Global_Rows          = -1234  ! from
INTEGER, SAVE          :: Old_Proc_Group           = -1234  ! previous
INTEGER, SAVE          :: Old_Gather_PE            = -1234  ! calls to
INTEGER, SAVE          :: Old_Decomp               = -1234  ! routine
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------
! 0.0 Can we use the same send/receive map that we calculated
!     last time round?
!-----------------------------------------------------------------

IF (.NOT. ALLOCATED(receive_map)) &
    ALLOCATE(receive_map(7,nproc_max))

IF ((global_row_len  /=  old_GLOBAL_ROW_LEN) .OR. &
    (global_rows     /=  old_GLOBAL_ROWS   ) .OR. &
    (proc_group      /=  old_PROC_GROUP    ) .OR. &
    (gather_pe       /=  old_GATHER_PE     ) .OR. &
    (current_decomp_type  /=  old_DECOMP   )) THEN
  !       Different arguments from the last call so we need
  !       to calculate a new send/receive map

  ! 1.0 Find the type of field (P or U) being done

  IF (global_rows     ==  glsizep(2) .AND. &
      global_row_len  ==  glsizep(1) ) THEN
    ! This should cover P field and U (global) field for C grid
    fld_type=fld_type_p
  ELSE IF (global_rows     ==  glsizev(2) .AND. &
          global_row_len  ==  glsizev(1) ) THEN
    fld_type=fld_type_v
  ELSE IF (global_rows     ==  glsizeu(2) .AND. &
          global_row_len  ==  glsizeu(1) ) THEN
    fld_type=fld_type_u
  ELSE IF (global_rows     ==  glsizer(2) .AND. &
          global_row_len  ==  glsizer(1) ) THEN
    fld_type=fld_type_r
  ELSE
    Cmessage = 'Unrecognised field type'
    ErrorStatus = 10
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF


  ! 2.0 Set up send map

  send_map(s_destination_pe,1) = gather_pe !  processor to send to

  send_map(s_base_address_in_send_array,1) = Offy*local_row_len+1+Offx
  !       first data to send

  IF (fld_type  ==  fld_type_p) THEN
    send_map(s_number_of_elements_in_item,1)= blsizep(2)
    send_map(s_stride_in_send_array,1) = blsizep(1) + 2 * Offx
    send_map(s_element_length,1) = blsizep(1)
    send_map(s_base_address_in_recv_array,1) = &
           datastart(1)+(datastart(2)-1)*global_row_len
  ELSE IF (fld_type  ==  fld_type_u) THEN
    send_map(s_number_of_elements_in_item,1)= blsizeu(2)
    send_map(s_stride_in_send_array,1) = blsizeu(1) + 2 * Offx
    send_map(s_element_length,1) = blsizeu(1)
    send_map(s_base_address_in_recv_array,1) = &
           datastart(1)+(datastart(2)-1)*global_row_len
  ELSE IF (fld_type  ==  fld_type_v) THEN
    send_map(s_number_of_elements_in_item,1)= blsizev(2)
    send_map(s_stride_in_send_array,1) = blsizev(1) + 2 * Offx
    send_map(s_element_length,1) = blsizev(1)
    send_map(s_base_address_in_recv_array,1) = &
           datastart(1)+(datastart(2)-1)*global_row_len
  ELSE IF (fld_type  ==  fld_type_r) THEN
    send_map(s_number_of_elements_in_item,1)= blsizer(2)
    send_map(s_stride_in_send_array,1) = blsizer(1) + 2 * Offx
    send_map(s_element_length,1) = blsizer(1)
    send_map(s_base_address_in_recv_array,1) = &
           datastartr(1)+(datastartr(2)-1)*global_row_len
  END IF

  send_map(s_stride_in_recv_array,1) = global_row_len
  !       stride between rows in global data


  ! 3.0 Set up the receive map (for PE GATHER_PE only)

  ! Assume here that this group consists of all processors
  ! We'll get some new GCG functionality soon to improve this

  n_mess_to_rec=0

  IF (mype  ==  gather_pe) THEN
    DO iproc=0,nproc-1
      receive_map(r_source_pe,iproc+1) = iproc

      IF (fld_type  ==  fld_type_p) THEN
        receive_map(r_base_address_in_recv_array,iproc+1) = &
                g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizep(1)
        receive_map(r_number_of_elements_in_item,iproc+1) = &
                g_blsizep(2,iproc)
        receive_map(r_element_length,iproc+1) = g_blsizep(1,iproc)
        receive_map(r_stride_in_send_array,iproc+1) = g_blsizep(1,iproc)
        receive_map(r_base_address_in_send_array,iproc+1) = &
                Offy*g_blsizep(1,iproc)+Offx+1

      ELSE IF (fld_type  ==  fld_type_u) THEN
        receive_map(r_base_address_in_recv_array,iproc+1) = &
                g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizeu(1)
        receive_map(r_number_of_elements_in_item,iproc+1) = &
                g_blsizeu(2,iproc)
        receive_map(r_element_length,iproc+1) = g_blsizeu(1,iproc)
        receive_map(r_stride_in_send_array,iproc+1) = &
                g_blsizeu(1,iproc)
        receive_map(r_base_address_in_send_array,iproc+1) = &
                Offy*g_blsizeu(1,iproc)+Offx+1

      ELSE IF (fld_type  ==  fld_type_v) THEN
        receive_map(r_base_address_in_recv_array,iproc+1) = &
                g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizev(1)
        receive_map(r_number_of_elements_in_item,iproc+1) = &
                g_blsizev(2,iproc)
        receive_map(r_element_length,iproc+1) = g_blsizev(1,iproc)
        receive_map(r_stride_in_send_array,iproc+1) = &
                g_blsizev(1,iproc)
        receive_map(r_base_address_in_send_array,iproc+1) = &
                Offy*g_blsizev(1,iproc)+Offx+1

      ELSE IF (fld_type  ==  fld_type_r) THEN
        receive_map(r_base_address_in_recv_array,iproc+1) = &
                  g_datastartr(1,iproc) + (g_datastartr(2,iproc) -1 ) &
                  * glsizer(1)
        receive_map(r_number_of_elements_in_item,iproc+1) = &
                g_blsizer(2,iproc)
        receive_map(r_element_length,iproc+1) = g_blsizer(1,iproc)
        receive_map(r_stride_in_send_array,iproc+1) = g_blsizer(1,iproc)
        receive_map(r_base_address_in_send_array,iproc+1) = &
                Offy*g_blsizer(1,iproc)+Offx+1
      END IF

      ! Should the following go in the above list?
      receive_map(r_stride_in_recv_array,iproc+1) = global_row_len


    END DO
    n_mess_to_rec=nproc
  END IF

  old_GLOBAL_ROW_LEN=global_row_len
  old_GLOBAL_ROWS=global_rows
  old_PROC_GROUP=proc_group
  old_GATHER_PE=gather_pe
  old_DECOMP=current_decomp_type

END IF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

flag=0  ! This is currently ignored at GCG v1.1
info=gc_none

CALL gcg_ralltoalle(local_field,send_map,1, &
                    local_row_len*local_rows, &
                    global_field,receive_map,n_mess_to_rec, &
                    global_row_len*global_rows, &
                    proc_group,flag,info)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Gather_Field_Real

!--------------------------------------------------------------------
!   Thus subroutine copes with the LOGICAL case
!--------------------------------------------------------------------
SUBROUTINE Rcf_Gather_Field_Log(local_field,    global_field, &
                                local_row_len,  local_rows,   &
                                global_row_len, global_rows,  &
                                gather_pe,      proc_group )

USE Ereport_Mod, ONLY: &
    Ereport

USE UM_ParVars   ! Use a lot of this

USE UM_ParCore, ONLY:  &
    mype,       nproc, &
    nproc_max

USE Field_Types, ONLY:      &
    fld_type_p, fld_type_r, &
    fld_type_u, fld_type_v

USE gcom_mod

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN)    :: Local_Row_Len  ! row length in local field
INTEGER, INTENT(IN)    :: Global_Row_Len ! row length in global field
INTEGER, INTENT(IN)    :: Local_Rows     ! number of rows locally
INTEGER, INTENT(IN)    :: Global_Rows    ! number of rows globally
INTEGER, INTENT(IN)    :: Gather_PE      ! PE on which to gather data
INTEGER, INTENT(IN)    :: Proc_Group     ! group of PEs involved

LOGICAL, INTENT(IN)    :: Local_Field( Local_Row_Len * Local_Rows )
LOGICAL, INTENT(OUT)   :: Global_Field( Global_Row_Len * Global_Rows)

! Local variables
INTEGER                :: Info        ! return code from Gcom
INTEGER                :: Fld_Type    ! P, U or V
INTEGER                :: Flag        ! dummy for Gcom
INTEGER                :: iproc       ! processor number
INTEGER                :: ErrorStatus
CHARACTER (LEN=errormessagelength)     :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_GATHER_FIELD_LOG'

INTEGER, SAVE          :: n_mess_to_rec
INTEGER, SAVE          :: send_map(7,1)
INTEGER, ALLOCATABLE, SAVE :: receive_map(:,:)
INTEGER, SAVE          :: Old_Global_Row_Len       = -1234  ! old values
INTEGER, SAVE          :: Old_Global_Rows          = -1234  ! from
INTEGER, SAVE          :: Old_Proc_Group           = -1234  ! previous
INTEGER, SAVE          :: Old_Gather_PE            = -1234  ! calls to
INTEGER, SAVE          :: Old_Decomp               = -1234  ! routine
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------
! 0.0 Can we use the same send/receive map that we calculated
!     last time round?
!-----------------------------------------------------------------

IF (.NOT. ALLOCATED(receive_map)) &
    ALLOCATE(receive_map(7,nproc_max))

IF ((global_row_len  /=  old_GLOBAL_ROW_LEN) .OR. &
    (global_rows     /=  old_GLOBAL_ROWS   ) .OR. &
    (proc_group      /=  old_PROC_GROUP    ) .OR. &
    (gather_pe       /=  old_GATHER_PE     ) .OR. &
    (current_decomp_type  /=  old_DECOMP   )) THEN
  !       Different arguments from the last call so we need
  !       to calculate a new send/receive map

  ! 1.0 Find the type of field (P or U) being done

  IF (global_rows     ==  glsizep(2) .AND. &
      global_row_len  ==  glsizep(1) ) THEN
    ! This should cover P field and U (global) field for C grid
    fld_type=fld_type_p
  ELSE IF (global_rows     ==  glsizev(2) .AND. &
          global_row_len  ==  glsizev(1) ) THEN
    fld_type=fld_type_v
  ELSE IF (global_rows     ==  glsizeu(2) .AND. &
          global_row_len  ==  glsizeu(1) ) THEN
    fld_type=fld_type_u
  ELSE IF (global_rows     ==  glsizer(2) .AND. &
          global_row_len  ==  glsizer(1) ) THEN
    fld_type=fld_type_r
  ELSE
    Cmessage = 'Unrecognised field type'
    ErrorStatus = 10
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF


  ! 2.0 Set up send map

  send_map(s_destination_pe,1) = gather_pe !  processor to send to

  send_map(s_base_address_in_send_array,1) = Offy*local_row_len+1+Offx
  !       first data to send

  IF (fld_type  ==  fld_type_p) THEN
    send_map(s_number_of_elements_in_item,1)= blsizep(2)
    send_map(s_stride_in_send_array,1) = blsizep(1) + 2 * Offx
    send_map(s_element_length,1) = blsizep(1)
    send_map(s_base_address_in_recv_array,1) = &
           datastart(1)+(datastart(2)-1)*global_row_len
  ELSE IF (fld_type  ==  fld_type_u) THEN
    send_map(s_number_of_elements_in_item,1)= blsizeu(2)
    send_map(s_stride_in_send_array,1) = blsizeu(1) + 2 * Offx
    send_map(s_element_length,1) = blsizeu(1)
    send_map(s_base_address_in_recv_array,1) = &
           datastart(1)+(datastart(2)-1)*global_row_len
  ELSE IF (fld_type  ==  fld_type_v) THEN
    send_map(s_number_of_elements_in_item,1)= blsizev(2)
    send_map(s_stride_in_send_array,1) = blsizev(1) + 2 * Offx
    send_map(s_element_length,1) = blsizev(1)
    send_map(s_base_address_in_recv_array,1) = &
           datastart(1)+(datastart(2)-1)*global_row_len
  ELSE IF (fld_type  ==  fld_type_r) THEN
    send_map(s_number_of_elements_in_item,1)= blsizer(2)
    send_map(s_stride_in_send_array,1) = blsizer(1) + 2 * Offx
    send_map(s_element_length,1) = blsizer(1)
    send_map(s_base_address_in_recv_array,1) = &
           datastartr(1)+(datastartr(2)-1)*global_row_len
  END IF


  send_map(s_stride_in_recv_array,1) = global_row_len
  !       stride between rows in global data


  ! 3.0 Set up the receive map (for PE GATHER_PE only)

  ! Assume here that this group consists of all processors
  ! We'll get some new GCG functionality soon to improve this

  n_mess_to_rec=0

  IF (mype  ==  gather_pe) THEN
    DO iproc=0,nproc-1
      receive_map(r_source_pe,iproc+1) = iproc

      IF (fld_type  ==  fld_type_p) THEN
        receive_map(r_base_address_in_recv_array,iproc+1) = &
                g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizep(1)
        receive_map(r_number_of_elements_in_item,iproc+1) = &
                g_blsizep(2,iproc)
        receive_map(r_element_length,iproc+1) = g_blsizep(1,iproc)
        receive_map(r_stride_in_send_array,iproc+1) = g_blsizep(1,iproc)
        receive_map(r_base_address_in_send_array,iproc+1) = &
                Offy*g_blsizep(1,iproc)+Offx+1

      ELSE IF (fld_type  ==  fld_type_u) THEN
        receive_map(r_base_address_in_recv_array,iproc+1) = &
                g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizeu(1)
        receive_map(r_number_of_elements_in_item,iproc+1) = &
                g_blsizeu(2,iproc)
        receive_map(r_element_length,iproc+1) = g_blsizeu(1,iproc)
        receive_map(r_stride_in_send_array,iproc+1) = &
                g_blsizeu(1,iproc)
        receive_map(r_base_address_in_send_array,iproc+1) = &
                Offy*g_blsizeu(1,iproc)+Offx+1

      ELSE IF (fld_type  ==  fld_type_v) THEN
        receive_map(r_base_address_in_recv_array,iproc+1) = &
                g_datastart(1,iproc)+(g_datastart(2,iproc)-1)*glsizev(1)
        receive_map(r_number_of_elements_in_item,iproc+1) = &
                g_blsizev(2,iproc)
        receive_map(r_element_length,iproc+1) = g_blsizev(1,iproc)
        receive_map(r_stride_in_send_array,iproc+1) = &
                g_blsizev(1,iproc)
        receive_map(r_base_address_in_send_array,iproc+1) = &
                Offy*g_blsizev(1,iproc)+Offx+1

      ELSE IF (fld_type  ==  fld_type_r) THEN
        receive_map(r_base_address_in_recv_array,iproc+1) = &
                  g_datastartr(1,iproc) + (g_datastartr(2,iproc) -1 ) &
                  * glsizer(1)
        receive_map(r_number_of_elements_in_item,iproc+1) = &
                g_blsizer(2,iproc)
        receive_map(r_element_length,iproc+1) = g_blsizer(1,iproc)
        receive_map(r_stride_in_send_array,iproc+1) = g_blsizer(1,iproc)
        receive_map(r_base_address_in_send_array,iproc+1) = &
                Offy*g_blsizer(1,iproc)+Offx+1
      END IF

      ! Should the following go in the above list?
      receive_map(r_stride_in_recv_array,iproc+1) = global_row_len


    END DO
    n_mess_to_rec=nproc
  END IF

  old_GLOBAL_ROW_LEN=global_row_len
  old_GLOBAL_ROWS=global_rows
  old_PROC_GROUP=proc_group
  old_GATHER_PE=gather_pe
  old_DECOMP=current_decomp_type

END IF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

flag=0  ! This is currently ignored at GCG v1.1
info=gc_none

CALL gcg_ralltoalle(local_field,send_map,1, &
                    local_row_len*local_rows, &
                    global_field,receive_map,n_mess_to_rec, &
                    global_row_len*global_rows, &
                    proc_group,flag,info)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Gather_Field_Log

END MODULE Rcf_Gather_Field_Mod



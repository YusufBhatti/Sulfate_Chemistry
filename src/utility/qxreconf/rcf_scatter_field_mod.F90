! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Scatters a field from one processor to many processors

MODULE Rcf_Scatter_Field_Mod

!  Subroutine Rcf_Scatter_Field_Real - for real data
!  Subroutine Rcf_Scatter_field_Log  - for logical data
!
! Description:
!  Takes a model field which is stored entirely on one processor
!  and distributes it over a group of processors.
!
! Method:
!  A send and receive map is constructed which instructs the GCOM
!  permute operation to do a scatter to all processors in the
!  group from the SCATTER_PE
!
!  Derived from UM4.5 code.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

!--------------------------------------------------------------------
!   ************************************************************
!      This subroutine has an overloaded interface so as to
!      cope cleanly with a number of different data types.
!   ***********************************************************
!--------------------------------------------------------------------
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb

INTERFACE Rcf_Scatter_Field
MODULE PROCEDURE Rcf_Scatter_Field_Log, Rcf_Scatter_Field_Real
END INTERFACE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SCATTER_FIELD_MOD'

CONTAINS

!-----------------------------------------------------------------
! This deals with real data
!-----------------------------------------------------------------

SUBROUTINE Rcf_Scatter_Field_Real(local_field,    global_field, &
                                  local_row_len,  local_rows,   &
                                  global_row_len, global_rows,  &
                                  scatter_pe,     proc_group    )

USE UM_ParVars

USE UM_ParCore, ONLY: &
    mype,             &
    nproc,            &
    nproc_max

USE Field_Types, ONLY:      &
    fld_type_p, fld_type_r, &
    fld_type_u, fld_type_v

USE gcom_mod

USE Ereport_Mod, ONLY: &
    Ereport

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE


! Subroutine Arguments:
INTEGER, INTENT(IN)   :: Local_Row_Len   ! len of rows in local field
INTEGER, INTENT(IN)   :: Local_Rows      ! num of rows in local field
INTEGER, INTENT(IN)   :: Global_Row_Len  ! len of rows in global field
INTEGER, INTENT(IN)   :: Global_Rows     ! num of rows in global field
INTEGER, INTENT(IN)   :: Scatter_PE      ! processor to scatter from
INTEGER, INTENT(IN)   :: Proc_Group      ! group ID of involved PEs

REAL, INTENT(IN)      :: Global_Field( global_row_len * global_rows )
REAL, INTENT(OUT)     :: Local_Field ( local_row_len * local_rows )

! Local variables
INTEGER               :: Info            ! return code from GCOM
INTEGER               :: fld_type
INTEGER               :: iproc
INTEGER               :: flag
INTEGER               :: ErrorStatus
CHARACTER (LEN=errormessagelength)    :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SCATTER_FIELD_REAL'

INTEGER, SAVE         :: n_mess_to_send
INTEGER, ALLOCATABLE, SAVE :: send_map(:,:)
INTEGER, SAVE         :: receive_map(7,1)
INTEGER, SAVE         :: Old_Global_Row_Len   = -1234   ! old values
INTEGER, SAVE         :: Old_Global_Rows      = -1234   ! from previous
INTEGER, SAVE         :: Old_Proc_Group       = -1234   ! calls to this
INTEGER, SAVE         :: Old_Scatter_PE       = -1234   ! routine
INTEGER, SAVE         :: Old_Decomp           = -1234   ! ...
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(send_map)) &
    ALLOCATE(send_map(7,nproc_max))

!-------------------------------------------------------

! 0.0 Can we use the same send/receive map that we calculated
!     last time round?

IF ((global_row_len  /=  old_GLOBAL_ROW_LEN) .OR. &
    (global_rows     /=  old_GLOBAL_ROWS   ) .OR. &
    (proc_group      /=  old_PROC_GROUP    ) .OR. &
    (scatter_pe      /=  old_SCATTER_PE    ) .OR. &
    (current_decomp_type  /=  old_DECOMP  )) THEN
  !       Different arguments from the last call so we need
  !       to calculate a new send/receive map

  ! 1.0 Find the type of field (P or U) being done

  IF (global_rows         ==  glsizep(2) .AND. &
      global_row_len      ==  glsizep(1) ) THEN
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
    WRITE(umMessage,'(A,I4,I4)')'Rows/Row length provided: ',global_rows,global_row_len
    CALL umPrint(umMessage,src='rcf_scatter_field_mod')
    WRITE(umMessage,'(A,I4,I4)')'P - fields:', glsizep(2), glsizep(1)
    CALL umPrint(umMessage,src='rcf_scatter_field_mod')
    WRITE(umMessage,'(A,I4,I4)')'V - fields:', glsizev(2), glsizev(1)
    CALL umPrint(umMessage,src='rcf_scatter_field_mod')
    WRITE(umMessage,'(A,I4,I4)')'U - fields:', glsizeu(2), glsizeu(1)
    CALL umPrint(umMessage,src='rcf_scatter_field_mod')
    WRITE(umMessage,'(A,I4,I4)')'R - fields:', glsizer(2), glsizer(1)
    CALL umPrint(umMessage,src='rcf_scatter_field_mod')
    Cmessage = 'Unable to determine field type from supplied rows/row_length'
    ErrorStatus = 10
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  ! 2.0 Set up the send map (for PE SCATTER_PE only)

  ! Assume here that this group consists of all processors
  ! We'll get some new GCG functionality soon to improve this

  n_mess_to_send=0

  IF (mype  ==  scatter_pe) THEN
    DO iproc=0,nproc-1
      send_map(s_destination_pe,iproc+1) = iproc

      IF (fld_type  ==  fld_type_p) THEN
        send_map(s_base_address_in_send_array,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*global_row_len
        send_map(s_number_of_elements_in_item,iproc+1) = &
                 g_blsizep(2,iproc)
        send_map(s_element_length,iproc+1) = g_blsizep(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizep(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = &
                 g_blsizep(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizep(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = g_blsizep(1,iproc)

      ELSE IF (fld_type  ==  fld_type_u) THEN
        send_map(s_base_address_in_send_array,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*global_row_len
        send_map(s_number_of_elements_in_item,iproc+1) = &
                 g_blsizeu(2,iproc)
        send_map(s_element_length,iproc+1) = g_blsizeu(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizeu(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = &
                 g_blsizeu(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizeu(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = g_blsizeu(1,iproc)

      ELSE IF (fld_type  ==  fld_type_v) THEN
        send_map(s_base_address_in_send_array,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*global_row_len
        send_map(s_number_of_elements_in_item,iproc+1) = &
                 g_blsizev(2,iproc)
        send_map(s_element_length,iproc+1) = g_blsizev(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizev(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = g_blsizev(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizev(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = g_blsizev(1,iproc)

      ELSE IF (fld_type  ==  fld_type_r) THEN
        send_map(s_base_address_in_send_array,iproc+1) = &
                 g_datastartr(1,iproc)+                  &
                (g_datastartr(2,iproc)-1)*global_row_len
        send_map(s_number_of_elements_in_item,iproc+1) = &
                 g_blsizer(2,iproc)
        send_map(s_element_length,iproc+1) = g_blsizer(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizer(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = &
                 g_blsizer(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizer(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = g_blsizer(1,iproc)

      END IF

      send_map(s_stride_in_send_array,iproc+1) = global_row_len

    END DO
    n_mess_to_send=nproc
  END IF

  ! 3.0 Set up the receive map

  receive_map(r_source_pe,1) = scatter_pe
  receive_map(r_base_address_in_recv_array,1) = &
              Offy*local_row_len+1+Offx

  IF (fld_type  ==  fld_type_p) THEN
    receive_map(r_number_of_elements_in_item,1) =  blsizep(2)
    receive_map(r_stride_in_recv_array,1) = blsizep(1) + 2 * Offx
    receive_map(r_element_length,1) = blsizep(1)
    receive_map(r_base_address_in_send_array,1) = &
              datastart(1)+(datastart(2)-1)*global_row_len
  ELSE IF (fld_type  ==  fld_type_u) THEN
    receive_map(r_number_of_elements_in_item,1) = blsizeu(2)
    receive_map(r_stride_in_recv_array,1) = blsizeu(1) + 2 * Offx
    receive_map(r_element_length,1) = blsizeu(1)
    receive_map(r_base_address_in_send_array,1) = &
              datastart(1)+(datastart(2)-1)*global_row_len
  ELSE IF (fld_type  ==  fld_type_v) THEN
    receive_map(r_number_of_elements_in_item,1) = blsizev(2)
    receive_map(r_stride_in_recv_array,1) = blsizev(1) + 2 * Offx
    receive_map(r_element_length,1) = blsizev(1)
    receive_map(r_base_address_in_send_array,1) = &
              datastart(1)+(datastart(2)-1)*global_row_len
  ELSE IF (fld_type  ==  fld_type_r) THEN
    receive_map(r_number_of_elements_in_item,1) = blsizer(2)
    receive_map(r_stride_in_recv_array,1) = blsizer(1) + 2 * Offx
    receive_map(r_element_length,1) = blsizer(1)
    receive_map(r_base_address_in_send_array,1) = &
              datastartr(1)+(datastartr(2)-1)*global_row_len
  END IF

  receive_map(r_stride_in_send_array,1) = global_row_len


  old_GLOBAL_ROW_LEN=global_row_len
  old_GLOBAL_ROWS=global_rows
  old_PROC_GROUP=proc_group
  old_SCATTER_PE=scatter_pe
  old_DECOMP=current_decomp_type

END IF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

flag=gc_none  ! This is currently ignored at GCG v1.1
info=gc_none

CALL gcg_ralltoalle(global_field, send_map,n_mess_to_send, &
                    global_row_len*global_rows,            &
                    local_field,  receive_map,  1,         &
                    local_row_len*local_rows,              &
                    proc_group, flag, info)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Scatter_Field_Real


!------------------------------------------------------------------
! This deals with logical fields
!------------------------------------------------------------------
SUBROUTINE Rcf_Scatter_Field_Log(local_field,    global_field, &
                                 local_row_len,  local_rows,   &
                                 global_row_len, global_rows,  &
                                 scatter_pe,     proc_group    )

USE UM_ParVars   !

USE UM_ParCore, ONLY: &
    mype,             &
    nproc,            &
    nproc_max

USE Field_Types, ONLY:      &
    fld_type_p, fld_type_r, &
    fld_type_u, fld_type_v

USE gcom_mod

USE Ereport_Mod, ONLY: &
    Ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Subroutine Arguments:
INTEGER, INTENT(IN)   :: Local_Row_Len   ! len of rows in local field
INTEGER, INTENT(IN)   :: Local_Rows      ! num of rows in local field
INTEGER, INTENT(IN)   :: Global_Row_Len  ! len of rows in global field
INTEGER, INTENT(IN)   :: Global_Rows     ! num of rows in global field
INTEGER, INTENT(IN)   :: Scatter_PE      ! processor to scatter from
INTEGER, INTENT(IN)   :: Proc_Group      ! group ID of involved PEs

LOGICAL, INTENT(IN)      :: Global_Field( global_row_len * global_rows )
LOGICAL, INTENT(OUT)     :: Local_Field ( local_row_len * local_rows )

! Local variables
INTEGER               :: Info            ! return code from GCOM
INTEGER               :: fld_type
INTEGER               :: iproc
INTEGER               :: flag
INTEGER               :: ErrorStatus
CHARACTER (LEN=errormessagelength)    :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SCATTER_FIELD_LOG'

INTEGER, SAVE         :: n_mess_to_send
INTEGER, ALLOCATABLE, SAVE :: send_map(:,:)
INTEGER, SAVE         :: receive_map(7,1)
INTEGER, SAVE         :: Old_Global_Row_Len   = -1234   ! old values
INTEGER, SAVE         :: Old_Global_Rows      = -1234   ! from previous
INTEGER, SAVE         :: Old_Proc_Group       = -1234   ! calls to this
INTEGER, SAVE         :: Old_Scatter_PE       = -1234   ! routine
INTEGER, SAVE         :: Old_Decomp           = -1234   ! ...
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(send_map)) &
    ALLOCATE(send_map(7,nproc_max))

!-------------------------------------------------------

! 0.0 Can we use the same send/receive map that we calculated
!     last time round?

IF ((global_row_len  /=  old_GLOBAL_ROW_LEN) .OR. &
    (global_rows     /=  old_GLOBAL_ROWS   ) .OR. &
    (proc_group      /=  old_PROC_GROUP    ) .OR. &
    (scatter_pe      /=  old_SCATTER_PE    ) .OR. &
    (current_decomp_type  /=  old_DECOMP  )) THEN
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

  ! 2.0 Set up the send map (for PE SCATTER_PE only)

  ! Assume here that this group consists of all processors
  ! We'll get some new GCG functionality soon to improve this

  n_mess_to_send=0

  IF (mype  ==  scatter_pe) THEN
    DO iproc=0,nproc-1
      send_map(s_destination_pe,iproc+1) = iproc

      IF (fld_type  ==  fld_type_p) THEN
        send_map(s_base_address_in_send_array,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*global_row_len
        send_map(s_number_of_elements_in_item,iproc+1) = &
                 g_blsizep(2,iproc)
        send_map(s_element_length,iproc+1) = g_blsizep(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizep(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = &
                 g_blsizep(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizep(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = g_blsizep(1,iproc)

      ELSE IF (fld_type  ==  fld_type_u) THEN
        send_map(s_base_address_in_send_array,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*global_row_len
        send_map(s_number_of_elements_in_item,iproc+1) = &
                 g_blsizeu(2,iproc)
        send_map(s_element_length,iproc+1) = g_blsizeu(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizeu(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = &
                 g_blsizeu(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizeu(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = g_blsizeu(1,iproc)

      ELSE IF (fld_type  ==  fld_type_v) THEN
        send_map(s_base_address_in_send_array,iproc+1) = &
                 g_datastart(1,iproc)+                   &
                (g_datastart(2,iproc)-1)*global_row_len
        send_map(s_number_of_elements_in_item,iproc+1) = &
                 g_blsizev(2,iproc)
        send_map(s_element_length,iproc+1) = g_blsizev(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizev(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = g_blsizev(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizev(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = g_blsizev(1,iproc)

      ELSE IF (fld_type  ==  fld_type_r) THEN
        send_map(s_base_address_in_send_array,iproc+1) = &
                 g_datastartr(1,iproc)+                  &
                (g_datastartr(2,iproc)-1)*global_row_len
        send_map(s_number_of_elements_in_item,iproc+1) = &
                 g_blsizer(2,iproc)
        send_map(s_element_length,iproc+1) = g_blsizer(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizer(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = g_blsizer(1,iproc)
        send_map(s_base_address_in_recv_array,iproc+1) = &
                 Offy*g_blsizer(1,iproc)+Offx+1
        send_map(s_stride_in_recv_array,iproc+1) = g_blsizer(1,iproc)

      END IF

      send_map(s_stride_in_send_array,iproc+1) = global_row_len

    END DO
    n_mess_to_send=nproc
  END IF

  ! 3.0 Set up the receive map

  receive_map(r_source_pe,1) = scatter_pe
  receive_map(r_base_address_in_recv_array,1) = &
              Offy*local_row_len+1+Offx

  IF (fld_type  ==  fld_type_p) THEN
    receive_map(r_number_of_elements_in_item,1) =  blsizep(2)
    receive_map(r_stride_in_recv_array,1) = blsizep(1) + 2 * Offx
    receive_map(r_element_length,1) = blsizep(1)
  ELSE IF (fld_type  ==  fld_type_u) THEN
    receive_map(r_number_of_elements_in_item,1) = blsizeu(2)
    receive_map(r_stride_in_recv_array,1) = blsizeu(1) + 2 * Offx
    receive_map(r_element_length,1) = blsizeu(1)
  ELSE
    receive_map(r_number_of_elements_in_item,1) = blsizev(2)
    receive_map(r_stride_in_recv_array,1) = blsizev(1) + 2 * Offx
    receive_map(r_element_length,1) = blsizev(1)
  END IF

  receive_map(r_base_address_in_send_array,1) = &
              datastart(1)+(datastart(2)-1)*global_row_len
  receive_map(r_stride_in_send_array,1) = global_row_len


  old_GLOBAL_ROW_LEN=global_row_len
  old_GLOBAL_ROWS=global_rows
  old_PROC_GROUP=proc_group
  old_SCATTER_PE=scatter_pe
  old_DECOMP=current_decomp_type

END IF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

flag=gc_none  ! This is currently ignored at GCG v1.1
info=gc_none

CALL gcg_ralltoalle(global_field, send_map,n_mess_to_send, &
                    global_row_len*global_rows,            &
                    local_field,  receive_map,  1,         &
                    local_row_len*local_rows,              &
                    proc_group, flag, info)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Scatter_Field_Log

END MODULE Rcf_Scatter_Field_Mod



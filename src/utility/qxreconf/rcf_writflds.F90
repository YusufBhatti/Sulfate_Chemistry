! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Description:
!  This routine writes a number of fields out to dump.
!
! Method:
!  Setpos is used to set the position in the dump to write too.
!  The addressing is calculated and a number of variables
!  (grid code and lbc_levels) required by lower level routines are set.
!  Data is written by rcf_write_multi


!-------------------------------------------------------------------
! Note that this routine does not have a module as it is called
! with different type arguents for D1
!-------------------------------------------------------------------

SUBROUTINE Rcf_WritFlds( nftout,   number_of_fields,                           &
    POSITION, lookup,len1_lookup,                                              &
    d1,   len_buf,  fixhd,                                                     &
    errorstatus,  cmessage, multi_pe)

USE UM_ParVars, ONLY:                                                       &
    blsizeu,                                                                 &
    blsizep,             blsizev

USE umPrintMgr, ONLY:                                                       &
    umPrint,                        &
    umMessage,                      &
    PrintStatus,                                                             &
    PrStatus_Diag

USE Rcf_Exppx_Mod, ONLY:                                                    &
    Rcf_Exppx

USE Rcf_Ppx_Info_Mod, ONLY:                                                 &
    STM_record_type


USE Rcf_Write_Multi_Mod, ONLY:                                              &
    Rcf_Write_Multi

USE Ereport_Mod, ONLY:                                                      &
    Ereport

USE Rcf_Level_Code_Mod, ONLY:                                               &
    Rcf_Level_Code

USE Rcf_Grid_Type_Mod, ONLY:                                                &
    Output_Grid

USE Rcf_Global_To_Local_Mod, ONLY:                                          &
    Rcf_Get_Fld_Type

USE io
USE lookup_addresses

USE cppxref_mod, ONLY:                                                      &
     ppx_atm_ozone,                                                          &
     ppx_atm_tzonal

USE nlstcall_mod, ONLY: ltimer

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)  :: nftout            ! Unit number for I/O
INTEGER, INTENT(IN)  :: number_of_fields ! No of fields to be
                                        ! written
INTEGER, INTENT(IN)  :: len_buf           ! Length of I/O buffer
INTEGER, INTENT(IN)  :: POSITION          ! Field number from
                                        ! which to begin I/O
INTEGER, INTENT(IN)  :: fixhd(*)          ! fixed length header
INTEGER, INTENT(IN)  :: Len1_Lookup       ! 1st dim of lookup
INTEGER, INTENT(IN)  :: Lookup(Len1_Lookup,*) ! lookup table
REAL,    INTENT(IN)  :: d1(*)             ! Data to write


INTEGER, INTENT(OUT)           :: ErrorStatus
CHARACTER(LEN=errormessagelength), INTENT(OUT) :: Cmessage

LOGICAL, INTENT(IN)            :: Multi_PE

! Local variables
INTEGER                :: k            ! index
INTEGER                :: D1_Off       ! Offset in D1 array
INTEGER                :: len_io
INTEGER                :: Field_Start  ! word address to begin I/O
INTEGER                :: Data_Write_Size ! data to write
INTEGER                :: Fld_Type
INTEGER                :: local_len    ! local size of field
INTEGER                :: grid_type
INTEGER                :: LBC_levels
INTEGER                :: field_model
INTEGER                :: field_sect
INTEGER                :: field_item
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_WRITFLDS'

REAL                   :: a_io
TYPE (stm_record_type) :: Stash_Record

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (ltimer) CALL timer( routinename, 3 )

errorstatus = 0
cmessage    = ' '

!  2. Buffer out NUMBER_OF_FIELDS fields of data:
D1_Off = 0
DO k = POSITION, POSITION + number_of_fields-1

  ! Location on disk from which to begin I/O
  Field_Start=lookup(lbegin,k)

  Data_Write_Size = Lookup( lbnrec, k)

  ! Position file pointer

  CALL setpos( nftout, Field_Start, ErrorStatus)


  ! Get some information about this field
  field_item  = MOD(lookup(42,k),1000)
  field_sect  = (lookup(42,k)-field_item)/1000
  field_model = lookup(45,k)

  stash_record = rcf_exppx( field_model, field_sect, field_item )
  grid_type = stash_record % grid_type

  ! For atmosphere zonal ozone fields - set to zonal grid type
  IF ( grid_type == ppx_atm_ozone .AND. lookup(lbnpt,k) == 1) THEN
    Stash_Record % grid_type = ppx_atm_tzonal
  END IF

  CALL Rcf_Write_Multi( nftout, d1( D1_Off + 1), Data_Write_Size,            &
      len_io, local_len, a_io,                                               &
      lookup(1,k), fixhd(12),                                                &
      Stash_Record, Multi_PE )

  ! Reset the STASHmaster grid-type to what it really should be
  Stash_Record % grid_type = grid_type

  ! Check for I/O errors
  IF (a_io /= -1.0 .OR. len_io /= Data_Write_Size ) THEN
    WRITE(umMessage,'('' *ERROR* Writing field no'',I5)')k
    CALL umPrint(umMessage,src='rcf_writflds')
    IF (fixhd(5) < 6 .OR. fixhd(5) > 10) THEN ! Not AC/Cx/Cov/ObSt
      ! DEPENDS ON: rcf_pr_look
      CALL rcf_pr_look(lookup,lookup,len1_lookup,k)
    END IF

    ! DEPENDS ON: ioerror
    CALL ioerror('rcf_writfields: buffer out of real data', &
        a_io, len_io,Data_Write_Size)
    ErrorStatus = NINT(a_io)+1
    cmessage    = 'Rcf_WRITFLDS:I/O error'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  ! Data summary used to be here - removed for time being

  ! Increment offset by size of level written out
  D1_Off = D1_Off + Local_Len

END DO

IF (ltimer) CALL timer( routinename, 4 )

IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_WritFlds



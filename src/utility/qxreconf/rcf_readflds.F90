! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! This routine reads in a number of fields from the dump.

! Description:
!   reads a number of fields from a dump into the D1 array.
!
! Method:
!   UMDP F3 is the relevant documentation.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


! No module for this routine as it is called with real/logical/integer
! arguments and this could cause problems. Also called from old
! F77 code and this would also cause problems here.

SUBROUTINE Rcf_ReadFlds(nftin, number_of_fields,        &
                        POSITION, lookup, len1_lookup,  &
                        d1, len_buf, fixhd,             &
                        icode, cmessage)

USE Field_Types, ONLY:      &
    fld_type_p,             &
    fld_type_u,             &
    fld_type_v

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Diag

USE Rcf_Exppx_Mod, ONLY: &
    Rcf_Exppx

USE Rcf_Ppx_Info_Mod, ONLY: &
    STM_record_type

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Read_Multi_Mod, ONLY: &
    Rcf_Read_Multi

USE Rcf_Level_Code_Mod, ONLY: &
    Rcf_Level_Code

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid

USE Rcf_Global_To_Local_Mod, ONLY: &
    Rcf_Get_Fld_Type

USE rcf_headaddress_mod, ONLY: &
    fh_lookupsize2

USE io
USE lookup_addresses
USE cppxref_mod, ONLY:  &
    ppx_atm_compressed, &
    ppx_atm_tall,       &
    ppx_atm_cuall,      &
    ppx_atm_cvall,      &
    ppx_atm_ozone,      &
    ppx_atm_tzonal

USE missing_data_mod, ONLY: imdi

USE nlstcall_mod, ONLY: ltimer

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

 
USE timer_mod, ONLY: timer
 
USE packing_codes_mod, ONLY: PC_No_CompressType, PC_WGDOS_Packing,             &
    PC_Cray32_Packing

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)  :: nftin                 ! Unit no. for I/O
INTEGER, INTENT(IN)  :: number_of_fields      ! No of fields to read
INTEGER, INTENT(IN)  :: len_buf               ! Length of I/O buffer
INTEGER, INTENT(IN)  :: POSITION              ! Field number from which
                                              ! to begin I/O
INTEGER, INTENT(IN)  :: fixhd(*)              ! Fixed lenght header
INTEGER, INTENT(IN)  :: len1_lookup           ! 1st dimension of lookup
INTEGER, INTENT(IN)  :: lookup(len1_lookup,*) ! PP Lookup table
REAL,    INTENT(OUT) :: d1(*)                 ! data space to be filled

INTEGER, INTENT(OUT)            :: icode      ! error code
CHARACTER (LEN=errormessagelength), INTENT(OUT) :: Cmessage   ! Error Message

! Local variables:---------------------------------------------
INTEGER, PARAMETER    :: unset = -1     ! flag
INTEGER   :: D1_Off         ! Offset in D1 array
INTEGER   :: k              ! index
INTEGER   :: len_io         ! Length of I/O returned by LENGTH
INTEGER   :: Data_Size      ! No of words of data on disk (inc WFIO pad)
INTEGER   :: Field_Start    ! word address to begin I/O
INTEGER   :: local_len      ! length of local part of field read in
INTEGER   :: Data_Read_Size ! No of words to read from disk
INTEGER   :: field_model    !  Model code for field
INTEGER   :: field_sect     !  Section code for field
INTEGER   :: field_item     !  Item code for field
INTEGER   :: grid_type      !  Grid code from stash
INTEGER   :: ErrorStatus
INTEGER   :: Fld_Type
INTEGER   :: expand_size    ! expanded data size

REAL      :: a_io

CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'RCF_READFLDS'

TYPE (STM_record_type)        :: Stash_Record
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! -------------------------------------------------------------

IF (ltimer) CALL timer( routinename, 3 )

icode = 0
Cmessage = ' '

!  Buffer in NUMBER_OF_FIELDS fields of real data:
D1_Off = 0
DO  k=POSITION,position+number_of_fields-1

  ! Location on disk from which to begin I/O
  Field_Start=lookup(lbegin,k)

  ! data_size contains the number of words of data used to store the
  ! field on disk
  IF (MOD(lookup(lbpack,k),10) == PC_Cray32_Packing) THEN
    Data_Size = (lookup(lblrec,k)+1)/2    ! 32 bit packed field
  ELSE
    Data_Size = lookup(lblrec,k)
  END IF

  ! data_read_size contains the number of words to data that need to
  ! be read in for a field. Each field has extra words of dummy data
  ! added at the end to ensure each field starts on a disk sector
  ! boundary.
  data_read_size = lookup(lbnrec,k)

  ! The last field on a dump does not have these extra words
  ! added.  So check against number of lookups and if next lookup is -99 assume
  ! we are at the end.  Otherwise assume we are the last field.
  IF (k < fixhd(fh_lookupsize2)) THEN
    IF (lookup(lbyr,k+1) == -99) THEN
      data_read_size = data_size
    END IF
  ELSE
    data_read_size = data_size
  END IF


  ! We will only deal with well-formed files. Thus an error is required
  ! if an old dump is encountered
  IF ( lookup(lblrec,k) /= 0 .AND.                             &
      (lookup(lbnrec,k) == 0 .OR. lookup(lbnrec,k) == imdi) ) THEN
    ErrorStatus = 10
    Cmessage = 'Invalid dump addressing: Possibly data file is '//&
        'ill-formed'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  ! Position file pointer

  CALL Setpos(nftin, Field_Start, icode)

  ! Get some information about this field
  field_item  = MOD(lookup(42,k),1000)
  field_sect  = (lookup(42,k)-field_item)/1000
  field_model = lookup(45,k)

  Stash_Record = Rcf_Exppx( field_model, field_sect, field_item )
  grid_type = Stash_Record % grid_type

  fld_type = Rcf_Get_Fld_Type( Stash_Record % grid_type)

  ! Special case for d1_grid_type for ancillary fields
  IF ( FixHd(5) == 4) THEN     ! Ancil data
    IF ( (MOD( Lookup(lbpack,k)/10, 10 ) == PC_No_CompressType )  .AND. &
          grid_type == ppx_atm_compressed ) THEN

      ! Compressed in stashmaster but uncompressed in header
      SELECT CASE( fld_type )
      CASE ( fld_type_p )
        Stash_Record % grid_type = ppx_atm_tall

      CASE ( fld_type_u )
        Stash_Record % grid_type = ppx_atm_cuall

      CASE ( fld_type_v )
        Stash_Record % grid_type = ppx_atm_cvall

      CASE DEFAULT
        ErrorStatus = 20
        Cmessage = 'Unknown grid type for Ancil data field.'
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )

      END SELECT
    END IF
  END IF

  ! For atmosphere zonal ozone fields - set to zonal grid type
  IF (grid_type  ==  ppx_atm_ozone .AND.  lookup(lbnpt,k) == 1) THEN
    Stash_Record % grid_type = ppx_atm_tzonal
  END IF

  ! Set the size of the expanded data buffer
  IF ( MOD (Lookup(lbpack, k), 10) == PC_WGDOS_Packing) THEN
    expand_size = Lookup(lbrow, k) * Lookup(lbnpt, k)
  ELSE
    expand_size = 2 * Data_Read_Size
  END IF

  CALL Rcf_Read_Multi(nftin, d1(D1_Off+1), Data_Read_Size,    &
                      expand_size, len_io, local_len, a_io,   &
                      lookup(1,k), fixhd(12), Stash_Record )

  ! Change the grid-type in the STASHmaster back to what it should
  ! be...
  Stash_Record % grid_type = grid_type

  ! Check for I/O errors
  IF (a_io  /=  -1.0 .OR. len_io  /=  Data_Read_Size) THEN
    WRITE(umMessage,'('' *ERROR* Reading field no'',I5)')k
    CALL umPrint(umMessage,src='rcf_readflds')
    IF (fixhd(5) <  6 .OR. fixhd(5) >  10) THEN ! Not AC/Cx/Cov/ObSt
      ! DEPENDS ON: rcf_pr_look
      CALL rcf_PR_Look( lookup,lookup,len1_lookup,k)
    END IF

    ! DEPENDS ON: ioerror
    CALL ioerror('buffer in of real data',a_io,len_io, Data_Read_Size)

    icode=NINT(a_io)+1
    Cmessage = 'Rcf_READFLDS:I/O error'
    CALL Ereport( RoutineName, Icode, Cmessage )
  END IF

  ! Data summary used to be here - removed for the time being

  ! increment size of level read in
  D1_Off = D1_Off + local_len

END DO

IF (ltimer) CALL timer( routinename, 4 )

IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Readflds

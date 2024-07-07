! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Parallel RCF interface to BUFFOUT

MODULE Rcf_Write_Multi_Mod

!  Subroutine Rcf_Write_Multi - parallel interface to buffout
!
! Description:
!  This routine provides an interface to BUFFOUT for the parallel
!  Reconfiguration. It is used where each process must write out a
!  local section of a global field.
!
! Method:
!  Each processor sends its local part of the global field to PE 0
!  which assembles all the parts, and then writes them to disk.
!  Fields which are compressed to land points are expanded before
!  sending to PE 0, PE 0 then compresses the global field before
!  writing it to disk.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_WRITE_MULTI_MOD'

CONTAINS

SUBROUTINE Rcf_Write_Multi(  nft, d1, isize, len_io, local_len,              &
    iostatus, lookup, fixhd12,                                               &
    stash_record, Multi_PE )

USE UM_ParVars, ONLY:                                                     &
    gc_all_proc_group,                                                    &
    glsize,                                                               &
    g_blsize

USE UM_ParCore, ONLY:                                                     &
    mype

USE UM_ParParams, ONLY:                                                   &
    halo_type_single

USE cppxref_mod, ONLY:  &
    ppx_atm_compressed, &
    ppx_atm_tzonal,     &
    ppx_atm_uzonal

USE Rcf_General_Gather_Field_Mod, ONLY:                                   &
    Rcf_General_Gather_Field

USE Rcf_Ppx_Info_Mod, ONLY:                                               &
    STM_Record_Type

USE Rcf_Global_To_Local_Mod, ONLY: &
    Rcf_Get_Fld_Type

USE io, ONLY: buffout

USE lookup_addresses

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

USE packing_codes_mod, ONLY: PC_Cray32_Packing

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN)    :: Nft         ! Fortran unit number
INTEGER, INTENT(IN)    :: Isize       ! no. of words to write out
INTEGER, INTENT(IN)    :: Fixhd12     ! 12th element of fixed hdr.
INTEGER, INTENT(IN)    :: Lookup(64)  ! Lookup table
INTEGER, INTENT(OUT)   :: Len_IO      ! no. or words written out
INTEGER, INTENT(OUT)   :: Local_Len   ! size of local field written out

REAL, INTENT(IN)       :: d1(*)       ! Array to write out
REAL, INTENT(OUT)      :: iostatus    ! return code

TYPE( STM_Record_Type ), INTENT(IN) :: stash_record

LOGICAL, INTENT(IN)    :: Multi_PE    ! flag to confirm 'gather' req'd

! Local variables
INTEGER                :: i                 ! looper
INTEGER                :: info              ! Gcom code
INTEGER                :: local_x           ! local x dimension
INTEGER                :: local_y           ! local y dimension
INTEGER                :: global_x          ! global x dimension
INTEGER                :: global_y          ! global y dimension
INTEGER                :: fld_type          ! p/u/v field
LOGICAL                :: l_compressed      ! Have we compressed data already?
REAL                   :: CompBuf(lookup(lbnrec)) ! compression buffer
REAL                   :: buf( isize * 2 )  ! buffer for data
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_WRITE_MULTI'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! ------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

iostatus     = -1.0
len_io       = isize
local_len    = 0
l_compressed = .FALSE.

! If data spread across multiple PE's then gather it.
IF (Multi_PE) THEN
  IF (MOD(lookup(lbpack),10) == PC_Cray32_Packing  .AND.                       &
      (stash_record % grid_type /= ppx_atm_compressed  .AND.      &
       stash_record % grid_type /= ppx_atm_tzonal      .AND.      &
       stash_record % grid_type /= ppx_atm_uzonal ) )  THEN

    ! Field needs to be 32-bit packed and we can use a 32-bit communication
    fld_type = Rcf_Get_Fld_Type(stash_record % grid_type)

    global_x = glsize(1, fld_type)
    global_y = glsize(2, fld_type)
    local_x  = g_blsize(1, fld_type, mype)
    local_y  = g_blsize(2, fld_type, mype)

    local_len = local_x * local_y

    ! Do the packing (from d1 into buf), then the gather (from buf into
    ! compbuf)
    CALL pack21( local_len, d1, buf)

    ! DEPENDS ON: gather_field_mpl32
    CALL gather_field_mpl32( buf,      compbuf,               &
                             local_x,  local_y,               &
                             global_x, global_y,              &
                             fld_type, halo_type_single,      &
                             0,        gc_all_proc_group)

    l_compressed = .TRUE.


  ELSE         ! Not 32-bit packed

    ! Gather the field from the local D1 array to buf
    CALL Rcf_General_Gather_Field( d1,        buf,                           &
                                   local_len, lookup(lblrec),                &
                                   stash_record,   0 )
  END IF

ELSE
  ! data is already held as a full copy on each PE.
  buf(: lookup(lblrec)) = d1(: lookup(lblrec))
END IF

IF (mype == 0) THEN
  !       Does this field need to be compressed?
  IF (MOD((lookup(lbpack)),10) == PC_Cray32_Packing .AND. .NOT. l_compressed)  &
  THEN
    IF (lookup(data_type) == 1) THEN
      ! DEPENDS ON: pack21
      CALL pack21( lookup(lblrec), buf, compbuf )
    END IF
  ELSE IF (.NOT. l_compressed) THEN ! no compression required - just do a copy
    DO i=1, lookup( lblrec )
      compbuf(i)=buf(i)
    END DO
  END IF

  ! Now write out the global field

  IF (MOD((lookup(lbpack)),10) == PC_Cray32_Packing) THEN
    IF (lookup(data_type) == 1) THEN
      ! Data is packed using CRAY 32 bit method - note that we need to write
      ! out 2*ISIZE 32 bit words using BUFFOUT32_F77

      ! DEPENDS ON: buffout32_f77
      CALL buffout32_f77( nft,compbuf,2*isize,len_io,iostatus )

      ! And then halve LEN_IO to satisfy tests against ISIZE
      len_io = len_io/2
    END IF
  ELSE
    ! For non-packed data
    CALL buffout( nft, compbuf, isize, len_io, iostatus )
  END IF

END IF ! am I PE 0 ?

! If the field was compressed for writing on disk, we need to compress
! and expand the field in memory. This ensures the same field exists in
! memory that would exist if this dump was read back in.

IF (MOD((lookup(lbpack)),10) == PC_Cray32_Packing) THEN
  IF (lookup(data_type) == 1) THEN
    ! DEPENDS ON: pack21
    CALL pack21(   local_len, d1, compbuf )
    ! DEPENDS ON: expand21
    CALL expand21( local_len, compbuf, d1 )
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Write_Multi
END MODULE Rcf_Write_Multi_Mod

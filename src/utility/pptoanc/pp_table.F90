! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Works out the lookup tables for the dump/ancillary    *
!           file header from the pp fields                       *
!
! Subroutine interface:
SUBROUTINE pp_table(pp_int,pp_real,nfields,lookup,rlookup,        &
 fieldsize,n,levn,m,runtot,number_of_codes,field_code,            &
 stash_code,add_wrap_pts,compress,pack32,wave,len1_levdepc,       &
 len2_levdepc,lev_dep_consts,len_realc,real_const,icode)

USE c_model_id_mod, ONLY: model_id
USE lookup_addresses

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE missing_data_mod, ONLY: &
    rmdi 

USE um_version_mod, ONLY: um_version_int

USE Packing_Codes_Mod, ONLY:                                                   &
    PC_LandMask_Compression,                                                   &
    PC_FieldIDX_CompressType,                                                  &
    PC_SeaMask_Compression,                                                    &
    PC_BitMask_CompressType,                                                   &
    PC_No_Packing,                                                             &
    PC_Cray32_Packing

IMPLICIT NONE

!
! Description:
!            Works out the lookup tables for the dump/ancillary
!            file header from the pp fields
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
! 1.0 Global variables 

! Subroutine arguments
!   Scalar arguments with intent(in):

INTEGER :: nfields   ! dimension for lookup tables

INTEGER :: fieldsize  ! size of field to be stored in anc file
INTEGER :: n          ! field number
INTEGER :: levn       ! level number
INTEGER :: m          ! field type
INTEGER :: number_of_codes
INTEGER :: len1_levdepc
INTEGER :: len2_levdepc
INTEGER :: len_realc

LOGICAL :: add_wrap_pts
LOGICAL :: compress
LOGICAL :: pack32
LOGICAL :: wave  ! T for wave dump creation


!   Array  arguments with intent(in):

INTEGER :: pp_int(45)
REAL :: pp_real(19)
INTEGER :: field_code(number_of_codes)  !stash and field codes
INTEGER :: stash_code(number_of_codes)  !input by user


!   Scalar arguments with intent(in/out):
INTEGER :: runtot  ! start address for this field on input and
                ! for next field on output
INTEGER :: icode

!   Array  arguments with intent(in/out):

INTEGER :: lookup(45,nfields)
REAL :: rlookup(46:64,nfields)

REAL :: lev_dep_consts(1+len1_levdepc*len2_levdepc)
REAL :: real_const(len_realc)


! Local Scalars
INTEGER :: i

!- End of header

DO i=1,45
  lookup(i,n) = 0
END DO

DO i=46,64
  rlookup(i,n) = 0.0
END DO

lookup(lbyr,n)   = pp_int(1)   ! lbyr
lookup(lbmon,n)  = pp_int(2)   ! lbmon
lookup(lbdat,n)  = pp_int(3)   ! lbdat
lookup(lbhr,n)   = pp_int(4)   ! lbhr
lookup(lbmin,n)  = pp_int(5)   ! lbmin
lookup(lbsec,n)  = pp_int(6)   ! lbsec
lookup(lbyrd,n)  = pp_int(7)   ! lbyrd
lookup(lbmond,n) = pp_int(8)   ! lbmond
lookup(lbdatd,n) = pp_int(9)   ! lbdatd
lookup(lbhrd,n)  = pp_int(10)  ! lbhrd
lookup(lbmind,n) = pp_int(11)  ! lbmind
lookup(lbsecd,n) = pp_int(12)  ! lbsecd

lookup(lbtim,n)  = pp_int(13) ! lbtim
lookup(lbft,n)   = pp_int(14) ! lbft
lookup(lbcode,n) = pp_int(16) ! lbcode
lookup(lbhem,n)  = pp_int(17) ! lbhem

!  1.0 Obtain rows and columns depending on compress and
!  add_wrap_pts.

IF (add_wrap_pts) THEN
  IF (compress) THEN
    lookup(lbrow,n) = 0  ! no rows if data is compressed
    lookup(lbnpt,n) = 0  ! no columns if data compressed
  ELSE
    lookup(lbrow,n) = pp_int(18) ! lbrow
    lookup(lbnpt,n) = pp_int(19)+2 ! lbnpt
  END IF

ELSE

  IF (compress) THEN
    lookup(lbrow,n) = 0  ! no rows if data is compressed
    lookup(lbnpt,n) = 0  ! no columns if data compressed
  ELSE
    lookup(lbrow,n) = pp_int(18) ! lbrow
    lookup(lbnpt,n) = pp_int(19) ! lbnpt
  END IF

END IF
!
IF (compress) THEN

  IF (.NOT. wave) THEN

    !C      compression using CFI
    lookup(lbpack,n) = 100*PC_LandMask_Compression + 10*PC_FieldIDX_CompressType

  ELSE

    !C      for wave dump compression using ls mask (compression to sea points)
    lookup(lbpack,n) = 100*PC_SeaMask_Compression + 10*PC_BitMask_CompressType

  END IF

ELSE
  lookup(lbpack,n) = PC_No_Packing  ! no compression
END IF

IF (pack32) THEN
  lookup(lbpack,n) = lookup(lbpack,n) + PC_Cray32_Packing  ! lbpack
END IF

lookup(lblrec,n) = fieldsize  ! lblrec
lookup(lbext,n)  = pp_int(20) ! lbext
lookup(lbrel,n)  = 3          ! lbrel
lookup(lbfc,n)   = pp_int(23) ! lbfc
lookup(lbproc,n) = pp_int(25) ! lbproc
lookup(lbvc,n)   = pp_int(26) ! lbvc
lookup(lbegin,n) = runtot     ! lbegin

lookup(lblev,n)  = levn       ! lblev (level number)

!MH for spectral wave energy the required value is already in pp-header
IF (pp_int(23) == 351) THEN
  lookup(lblev,n) = pp_int(33) ! wave model freq number
  lookup(44,n)    = pp_int(44) ! wave model dir number
END IF

lookup(lbproj,n) = pp_int(lbproj)
lookup(lbtyp,n)  = pp_int(lbtyp)
lookup(lblev,n)  = pp_int(lblev)
lookup(lbplev,n) = pp_int(lbplev)

lookup(lbsrce,n) = um_version_int*10000 + model_id
lookup(data_type,n) = pp_int(data_type)  !   data type

IF (lookup(data_type,n) <  1 .OR. lookup(data_type,n) >  3) THEN
  WRITE(umMessage,*) '********** WARNING ****************** '
  CALL umPrint(umMessage,src='pp_table')
  WRITE(umMessage,*) ' Data Type= ',lookup(data_type,n),' for field ',n, &
              ' is not recognised.'
  CALL umPrint(umMessage,src='pp_table')
  WRITE(umMessage,*) ' Either correct PP Header or set it through',      &
              ' the HEADER_DATA namelist.'
  CALL umPrint(umMessage,src='pp_table')
  WRITE(umMessage,*) '********** WARNING ****************** '
  CALL umPrint(umMessage,src='pp_table')
END IF

lookup(naddr,n) = runtot     ! start address in data
lookup(item_code,n) = pp_int(item_code)

lookup(model_code,n) = pp_int(model_code) ! sub model identifier

WRITE(umMessage,*) 'Field No ',n,' PP Field Code = ',lookup(lbfc,n),     &
                          ' Stash Code  = ',lookup(item_code,n)
CALL umPrint(umMessage,src='pp_table')

runtot=runtot+fieldsize

rlookup(blev,n)   = pp_real(52-45) ! blev / hybrid lev 'B' value
rlookup(brlev,n)  = pp_real(53-45) ! brlev
rlookup(bhlev,n)  = pp_real(54-45) ! bhlev / hybrid lev 'A' value
rlookup(bhrlev,n) = pp_real(55-45) ! bhrlev
rlookup(bplat,n)  = pp_real(56-45) ! bplat
rlookup(bplon,n)  = pp_real(57-45) ! bplon
rlookup(bgor,n)   = pp_real(58-45) ! bgor
rlookup(bzy,n)    = pp_real(59-45) ! bzy
rlookup(bdy,n)    = pp_real(60-45) ! bdy
rlookup(bzx,n)    = pp_real(61-45) ! bzx
rlookup(bdx,n)    = pp_real(62-45) ! bdx

rlookup(bmdi,n) = rmdi           ! bmdi
rlookup(bmks,n) = pp_real(64-45) ! bmks

! for spectral wave energy the required value is set from real_const
IF (pp_int(23) == 351) THEN
  rlookup(blev,n) = lev_dep_consts(pp_int(33)) ! wave model freq
  rlookup(bhlev,n)= (pp_int(44)-1)*real_const(13) ! direction
END IF

9999 CONTINUE
RETURN
END SUBROUTINE pp_table

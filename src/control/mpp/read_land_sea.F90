! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Parallel UM : Reads in the local section of Land-Sea Mask.

! Subroutine Interface:
SUBROUTINE read_land_sea(nft,icode,lookup,loc_len1_lookup,       &
                         loc_len2_lookup,fixhd,loc_len_fixhd)


USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE io
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParCore,        ONLY: mype, nproc
USE UM_ParParams,      ONLY: halo_type_no_halo
USE Field_Types,       ONLY: fld_type_p
USE atm_land_sea_mask, ONLY: atmos_landmask, atmos_landmask_local,             &
                             atmos_number_of_landpts,                          &
                             atmos_number_of_landpts_proc
USE lookup_addresses

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY: &
    land_field

USE missing_data_mod, ONLY: imdi
USE errormessagelength_mod, ONLY: errormessagelength
USE Packing_Codes_Mod, ONLY: PC_Cray32_Packing
IMPLICIT NONE

! Description:
!  This routine reads the land-sea mask (LSM) from the dump.
!  It is required for unpacking and packing fields which are 
!  stored compressed to land points.

! Method:
!  The position of the LSM within the dump is found from examining
!  the LOOKUP headers, it is then read in, and the relevant part
!  of the field sent to each processor. The local number of land
!  points is counted, and the LAND_FIELD variable is reset to this
!  new value.
!  Note : Halos can contain land points - but only those halos
!         which are updated by SWAPBNDS.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Subroutine Arguments:

INTEGER, INTENT(IN) :: nft              ! IN : FORTRAN unit number
INTEGER, INTENT(IN) :: loc_len1_lookup  ! IN : Dimension of the LOOKUP array
INTEGER, INTENT(IN) :: loc_len2_lookup  ! IN : Dimension of the LOOKUP array
INTEGER, INTENT(IN) :: loc_len_fixhd    ! IN : Dimension of the FIXHD array

INTEGER, INTENT(IN) :: lookup(loc_len1_lookup,loc_len2_lookup)
                                            ! IN : LOOKUP array from dump header
INTEGER, INTENT(IN) :: fixhd(loc_len_fixhd) ! IN : FIXHD array from dump header

REAL, INTENT(OUT) ::   icode               ! OUT : Return code

! Local variables

INTEGER :: i,j,k,word_address,ipts,iproc,info,len_io
INTEGER :: landpts_local,local_off,global_off,local_landmask_size

LOGICAL::read_mask
LOGICAL, ALLOCATABLE :: atmos_landmask_tmp(:) ! Space to read LSM from disk


INTEGER                      ::  errorstatus = 0
CHARACTER (LEN=errormessagelength)           ::  cmessage = ' '
INTEGER                      ::  i_icode
REAL :: io_ret(3)

! Local Parameters
CHARACTER (LEN=*), PARAMETER ::  RoutineName = 'READ_LAND_SEA'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! --------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
icode=-1.0
ipts=0
len_io=0

! Allocate arrays for global and processor land-sea mask
ALLOCATE ( atmos_landmask(glsize(1,fld_type_p)*glsize(2,fld_type_p)) )
ALLOCATE ( atmos_landmask_local(lasize(1,fld_type_p,halo_type_no_halo)*        &
                                lasize(2,fld_type_p,halo_type_no_halo)) )
ALLOCATE ( atmos_number_of_landpts_proc(0:nproc-1) )

! Find location of LSM in the dump

! disable read-broadcast for buffin by pe 0
CALL set_unit_bcast_flag(nft)
IF (mype  ==  0) THEN

  read_mask=.FALSE.

  DO i=1,loc_len2_lookup
    IF (lookup(item_code,i)  ==  30) THEN
      read_mask=.TRUE.
      EXIT
    END IF
  END DO

  IF (.NOT. read_mask) THEN
    WRITE(umMessage,'(a,i4)')                                               &
    'Error in READ_LAND_SEA_MASK: Missing Field of Type ',item_code
    CALL umPrint(umMessage,src='read_land_sea')
    errorstatus = 1
    cmessage = 'Missing Field or Wrong Dumpfile'

    CALL ereport(routinename,errorstatus,cmessage)
  END IF

  k=i
  word_address=1
  ! Old Format dumpfiles
  IF ((lookup(lbnrec,k) == 0) .OR.                                 &
    ! Prog lookups in dump before vn3.2:
        ((lookup(lbnrec,k) == imdi) .AND. (fixhd(12) <= 301))) THEN
    ! Dump and ancillary files
    word_address=1
    IF (i  >   1) THEN
      DO k=2,i
        IF (MOD((lookup(lbpack,k-1)),10) == PC_Cray32_Packing) THEN
          ipts=(lookup(lblrec,k-1)+1)/2
        ELSE
          ipts=(lookup(lblrec,k-1))
        END IF
        word_address=word_address+ipts
      END DO
    END IF
    word_address=fixhd(160)+word_address-2
    ipts=lookup(lblrec, i)
  ELSE
    ! PP type files and new format Dumpfiles (vn4.4 onwards)
    word_address=lookup(lbegin,i)
    ! Use the stored round-up value
    ipts=lookup(lbnrec,i)
  END IF

  CALL setpos(nft,word_address,i_icode)
  IF (i_icode  /=  0) THEN
    WRITE(umMessage,*) 'READ_LAND_SEA: Error Return from SETPOS',  &
               ' Status is ',i_icode
    CALL umPrint(umMessage,src='read_land_sea')

    errorstatus = 3
    cmessage = 'Error from SETPOS - see output for Status'

    CALL ereport(routinename,errorstatus,cmessage)

  END IF

  ! Read the LSM in to PE 0
  ! -----------------------
  ! Allocate enough space to read 'well formed' field from disc.
  ALLOCATE ( atmos_landmask_tmp(ipts) )

  CALL buffin(nft,atmos_landmask_tmp,ipts,                     &
                     len_io,icode)

  ! copy just the field data.
  atmos_landmask(:) =  &
                 atmos_landmask_tmp(1:glsize(1,fld_type_p)*glsize(2,fld_type_p))
  DEALLOCATE ( atmos_landmask_tmp )

END IF   ! (mype == 0)
CALL clear_unit_bcast_flag(nft)! Restore broadcast flag

! Broadcast the I/O Status to the other PE's

io_ret(1)=len_io
io_ret(2)=icode
io_ret(3)=ipts

CALL gc_rbcast(99, 3, 0, nproc, info, io_ret)

len_io=NINT(io_ret(1))
icode=io_ret(2)
ipts=io_ret(3)

! Check the I/O Status on all PE'e

IF ((icode /= -1.0) .OR. (len_io /= ipts)) THEN
  WRITE(umMessage,*)'ERROR READING DUMP ON UNIT ',nft
  CALL umPrint(umMessage,src='read_land_sea')
  ! DEPENDS ON: ioerror
  CALL ioerror('BUFFER IN FROM READ_LAND_SEA_MASK',               &
   icode,len_io,ipts)
  errorstatus = 1
  WRITE(cmessage,'(A,I5)') 'Error reading dump on using ',nft
  CALL ereport(routinename,errorstatus,cmessage)
END IF

! Broadcast the global LSM to all processors

CALL gc_ibcast(100,glsize(1,fld_type_p)*glsize(2,fld_type_p),     &
               0,nproc,info,atmos_landmask)

! Copy my local part of the full LSM into atmos_landmask_local

atmos_landmask_local(:)=.FALSE.

DO j=1,blsize(2,fld_type_p)

  local_off=(j-1)*lasize(1,fld_type_p,halo_type_no_halo)

  global_off=(j-1+datastart(2)-1)*glsize(1,fld_type_p)  &
              +datastart(1)-1

  DO i=1,blsize(1,fld_type_p)

    atmos_landmask_local(local_off+i)=atmos_landmask(global_off+i)

  END DO ! i
END DO ! j

! Count the number of global land points

atmos_number_of_landpts=0
DO i=1,glsize(1,fld_type_p)*glsize(2,fld_type_p)
  IF (atmos_landmask(i))                                          &
      atmos_number_of_landpts=atmos_number_of_landpts+1
END DO

! Do a swap to get land points in halo areas

landpts_local=0
DO i=1,lasize(1,fld_type_p,halo_type_no_halo)*                    &
       lasize(2,fld_type_p,halo_type_no_halo)
  IF (atmos_landmask_local(i))                                    &
    landpts_local=landpts_local+1
END DO

IF (landpts_local  /=  land_field) THEN
  WRITE(umMessage,*) 'PE ',mype,' : LAND_FIELD is being reset from ',     &
             land_field,' to ',landpts_local
  CALL umPrint(umMessage,src='read_land_sea')
  land_field=landpts_local
END IF

atmos_number_of_landpts_proc(mype) = landpts_local

DO i=0,nproc-1
  CALL gc_ibcast(101,1,i,nproc,info,atmos_number_of_landpts_proc(i))
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_land_sea


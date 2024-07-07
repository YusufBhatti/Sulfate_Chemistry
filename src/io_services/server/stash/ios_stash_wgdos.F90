! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

!
! An F90 wrapper to WGDOS packing routines
!

MODULE ios_stash_wgdos

USE IOS_print_mgr, ONLY:                                                    &
    IOS_print,                                                               &
    IOS_Verbosity,                                                           &
    IOS_message
IMPLICIT NONE


INTEGER, PARAMETER :: packedheaderwords32 = 3



CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IOS_STASH_WGDOS'

CONTAINS

SUBROUTINE ios_stash_pack_wgdos(                                             &
    field,                                                                   &
    compressedField,                                                         &
    num,                                                                     &
    accuracy,                                                                &
    StashItem,                                                               &
    missingDataIndicator)

USE ios
USE IOS_Common
USE UM_Types
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength
USE wgdos_packing_mod, ONLY: wgdos_compress_field

IMPLICIT NONE

INTEGER, INTENT(IN)  :: accuracy ! accuracy is accuracy as a power of two
INTEGER, INTENT(IN)  :: StashItem! 1000*Section+Item STASH code 
REAL,    INTENT(IN)  :: field(:,:)          ! is the field to compress
REAL,    INTENT(IN)  :: MissingDataIndicator! parameter to the algorithm
INTEGER, INTENT(OUT) :: num
REAL(KIND=real32),                                                         &
    INTENT(OUT)      :: compressedField(:)  ! is the output buffer
INTEGER(KIND=integer64) :: compressedField_int(SIZE(compressedField)/2)

INTEGER              :: ix  ! ix,iy are the field dimensions
INTEGER              :: iy
INTEGER              :: icode               ! is the error code
INTEGER              :: i
CHARACTER(LEN=errormessagelength)       :: message     ! is the error message


! params for dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOS_STASH_PACK_WGDOS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode=0

! Note that input buffer is 64 bit, but output buffer is 32 bit, leading to
! occasional odd factors of two in specifying buffer sizes.

ix=SIZE(field,1)
iy=SIZE(field,2)

IF ( SIZE(compressedField) < 2*(ix + packedheaderwords32)*iy +             &
    packedheaderwords32 ) THEN  
!$OMP CRITICAL(internal_write)
  WRITE(message,'(A,I10,A,I10,A)')                                         &
      'size of compressedField array (',SIZE(compressedField),             &
      ') smaller than potential need  (',2*(ix + packedheaderwords32 )*iy,')'
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport("IO SERVER WGDOS MODULE",99,message)
END IF

CALL wgdos_compress_field(field, SIZE(field), compressedfield_int,          &
    SIZE(compressedfield_int), ix, num, accuracy, MissingDataIndicator, &
    StashItem, icode)
compressedField = TRANSFER(compressedField_int, compressedField)

! Compression factor - note factor of 2 as num is in 32 bit words
IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
  IF (ix*iy>0) THEN !Avoid div0
!$OMP CRITICAL(internal_write)
    WRITE(ios_message,'(A,F5.2,A,I4)')                                     &
        'Info: Stash Server: WGDOS Compression is ',                       &
        (100.0*num)/(2*ix*iy),'%, for accuracy of ',accuracy
!$OMP END CRITICAL(internal_write)
    CALL ios_print(ios_message,src='ios_stash_wgdos')
  ELSE
    WRITE(ios_message,'(A)')                                               &
        'Info: Stash Server: WGDOS Compression: No Data'
    CALL ios_print(ios_message,src='ios_stash_wgdos')
  END IF
END IF

IF (num >= SIZE(compressedField)) THEN
  WRITE(ios_message,'(A)')'A buffer overrun may have occurred'
  CALL IOS_print(IOS_message,src='ios_stash_wgdos')
  WRITE(ios_message,'(A)')' '
  CALL IOS_print(IOS_message,src='ios_stash_wgdos')
  WRITE(ios_message,'(A)') 'Arguments after call to wgdos_compress_field:'
  CALL IOS_print(IOS_message,src='ios_stash_wgdos')
!$OMP CRITICAL(internal_write)
  WRITE(ios_message,'(A,I6,A,I12)') 'StashItem: ',StashItem,'data size: ',ix*iy
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_stash_wgdos')
!$OMP CRITICAL(internal_write)
  WRITE(ios_message,'(A,I12)')   'compressed field size: ',                &
      SIZE(compressedField)
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_stash_wgdos')
!$OMP CRITICAL(internal_write)
  WRITE(ios_message,'(A,I12)')   'x:',ix
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_stash_wgdos')
!$OMP CRITICAL(internal_write)
  WRITE(ios_message,'(A,I12)')   'y:',iy
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_stash_wgdos')
!$OMP CRITICAL(internal_write)
  WRITE(ios_message,'(A,I12)')   'num (32 bit words after compression:',num
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_stash_wgdos')
!$OMP CRITICAL(internal_write)
  WRITE(ios_message,'(A,I12)')   'accuracy:',accuracy
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_stash_wgdos')
!$OMP CRITICAL(internal_write)
  WRITE(ios_message,'(A,F32.16)')'MDI:',MissingDataIndicator
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_stash_wgdos')
!$OMP CRITICAL(internal_write)
  WRITE(ios_message,'(A,I12)')   'icode:',icode
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_stash_wgdos')

  CALL IOS_Ereport("IO SERVER WGDOS MODULE",99,                            &
      'Potential buffer overrun due to -ve compression')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE  ios_stash_pack_wgdos

END MODULE ios_stash_wgdos

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
! Purpose: To repack data from the input array FIELD and return

SUBROUTINE re_pack(pack_type,full_size,field,num_cray_words,      &
                   ilabel,rlabel,pp_fixhd,icode,cmessage)
USE wgdos_packing_mod, ONLY: wgdos_compress_field
USE lookup_addresses
USE errormessagelength_mod, ONLY: errormessagelength
USE packing_codes_mod, ONLY: PC_No_Packing, PC_WGDOS_Packing,                  &
    PC_Cray32_Packing, PC_RunLength_Packing

IMPLICIT NONE
INTEGER ::                                                        &
     pack_type                                                    &
                          !IN  The type of packing used
    ,full_size                                                    &
                          !IN  The full unpacked size of a field
    ,ilabel(45)                                                   &
                          !OUT holds integer part of LOOKUP
    ,icode                                                        &
                          !OUT Non zero for any error
    ,pp_fixhd(*)          !IN  PPfile fixed length header
REAL ::                                                           &
     field(full_size)                                             &
                          !INOUT On Input contains data.On output
    ,rlabel(19)           !    holds real part of LOOKUP
CHARACTER(LEN=errormessagelength) :: cmessage    
                          !OUT Will contain any error mesages.
!
!     LOCAL  VARIABLES
REAL ::                                                           &
     work_array(full_size)                                        &
                            ! WORK array used for packing
    ,amdi                   ! Missing data indicator.
INTEGER ::                                                        &
     int_work_array(full_size)                                    &
                            ! INTEGER WORK array used for packing
    ,ixx                                                          &
                            ! X dimension for packing
    ,iyy                                                          &
                            ! Y dimension for packing
    ,isc                                                          &
                            ! Accuracy required for WGDOS packing
    ,idum                                                         &
                            ! Dummy variable
    ,num_cray_words                                               &
                            ! IN no of values in an input field
    ,num_unpack_values                                            &
                            ! Number of numbers originally packed
    ,i                                                            &
                            ! Loop counter
    ,pack_code                                                    &
                            ! Packing actually used
    ,grib_packing           ! OUT - profile for packing
!
!
amdi=rlabel(18)
! Packing actually used (might be different then packing wanted)
pack_code=pack_type

IF (pack_type == PC_WGDOS_Packing) THEN     ! WGDOS packing
  ixx=ilabel(lbnpt)
  isc=NINT(rlabel(6))
  CALL wgdos_compress_field(field,full_size,int_work_array,full_size    &
      ,ixx,num_cray_words,isc,amdi,ilabel(item_code),icode)
  work_array = TRANSFER(int_work_array, work_array)
ELSE IF (pack_type == PC_Cray32_Packing) THEN !  32 Bit CRAY packing

ELSE IF (pack_type == PC_RunLength_Packing) THEN ! Run length encoding
  ixx=ilabel(lbrow)
  iyy=ilabel(lbnpt)
  ! DEPENDS ON: runlen_encode
  CALL runlen_encode(field,ixx*iyy,work_array,ixx*iyy,            &
                     num_cray_words,amdi,icode,cmessage)
  ! Size of run length encoded data is greater than unpacked
  ! field therefore leave field unpacked.
  IF (num_cray_words  >=  ixx*iyy) THEN
    pack_code = PC_No_Packing
    DO i=1,ixx*iyy
      work_array(i) = field(i)
    END DO
    num_cray_words = ixx*iyy
  END IF

ELSE
  icode=6
  cmessage=' UNPACK - packing type not yet supported'
END IF
DO i=1,NUM_cray_words
  field(i)=work_array(i)
END DO
ilabel(data_type)=1  ! The data type must now be real
ilabel(lbpack)=ilabel(lbpack)+pack_code ! data now packed
RETURN
END SUBROUTINE re_pack

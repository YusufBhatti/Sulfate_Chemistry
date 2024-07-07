! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Packing Tools

MODULE Packing_Codes_Mod

! This module contains parameters for the different packing codes used to make
! up the lbpack lookup header value for fields.

! See also: UMDP F03

IMPLICIT NONE

PRIVATE

INTEGER, PUBLIC, PARAMETER ::                                                  &

  ! N1 Packing Codes [Packing]
  PC_No_Packing = 0,                                                           &
  PC_WGDOS_Packing = 1,                                                        &
  PC_Cray32_Packing = 2,                                                       &
  PC_GRIB_Packing = 3,                                                         &
  PC_RunLength_Packing = 4,                                                    &

  ! N2 Packing Codes [Data Compression Type]
  PC_No_CompressType = 0,                                                      &
  PC_FieldIDX_CompressType = 1,                                                &
  PC_BitMask_CompressType = 2,                                                 &

  ! N3 Packing Codes [Compression]
  PC_No_Compression = 0,                                                       &
  PC_LandMask_Compression = 1,                                                 &
  PC_SeaMask_Compression = 2,                                                  &

  ! N4 Packing Codes [Number Format]
  PC_Native_Format = 0,                                                        &
  PC_IBM_Format = 1,                                                           &
  PC_Cray_Format = 2,                                                          &
  PC_IEEE_Format = 3,                                                          &
  PC_GRIB_Format = 4,                                                          &
  PC_VAX_Format = 5

END MODULE
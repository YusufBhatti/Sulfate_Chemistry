! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE NEWPACK:--------------------------------------
!
!    Purpose: Packing codes stored in LOOKUP(21,K) & LOOKUP(39,K)
!             are changed from pre vn2.8 values to
!             specification required at release 2.8
!
!
!    Documentation: UM Documentation Paper F3
!
!      -----------------------------------------------------------------
!
!      Code Owner: Please refer to the UM file CodeOwners.txt
!      This file belongs in section: Dump I/O

SUBROUTINE newpack                                                &
(lookup,len1_lookup,len2_lookup)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE lookup_addresses, ONLY: lbpack

USE Packing_Codes_Mod, ONLY: PC_No_Packing, PC_No_CompressType,                &
    PC_No_Compression, PC_Cray32_Packing, PC_FieldIDX_CompressType,            &
    PC_BitMask_CompressType, PC_LandMask_Compression

IMPLICIT NONE

INTEGER ::                                                        &
 len1_lookup                                                      &
,len2_lookup                                                      &
,lookup(len1_lookup,len2_lookup)

INTEGER ::                                                        &
 n1                                                               &
,n2                                                               &
,n3                                                               &
,k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NEWPACK'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO k=1,len2_lookup
  n1=PC_No_Packing
  n2=PC_No_CompressType
  n3=PC_No_Compression
  IF (lookup(lbpack,k) == -2) n1=PC_Cray32_Packing
  ! Ocean field packed using index array
  IF (lookup(lbpack,k) >  9 .AND. lookup(lbpack,k) <  100) THEN
    n2=PC_FieldIDX_CompressType
    n3=lookup(lbpack,k)-10
  END IF
  ! Ocean field compressed using bit mask
  IF (lookup(lbpack,k) > 99) THEN
    n2=PC_BitMask_CompressType
    n3=lookup(lbpack,k)-100
  END IF
  ! Real field stored at land pts
  IF (lookup(39,k) == 4) THEN
    lookup(39,k)=1
    n2=PC_BitMask_CompressType
    n3=PC_LandMask_Compression
  END IF
  ! Integer field stored at land pts
  IF (lookup(39,k) == 5) THEN
    lookup(39,k)=2
    n2=PC_BitMask_CompressType
    n3=PC_LandMask_Compression
  END IF

  lookup(lbpack,k)=100*n3+10*n2+n1
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE newpack

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Data module defining um_header_type

MODULE Rcf_UMhead_Mod

! Description:
!  Contains parameters used in UM dumps,
!  a structure to hold a UM header
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

!   Global Constants
INTEGER, PARAMETER :: LenFixHd = 256   ! Length of Fixed Length Header.

!   Global Type Definitions:
TYPE UM_header_type
  SEQUENCE

  INTEGER ::      LenIntC         ! Length of Integer Constants array.
  INTEGER ::      LenRealC        ! Length of Real Constants array.
  INTEGER ::      Len1LevDepC     ! 1st dimension \  of Level Dependent
  INTEGER ::      Len2LevDepC     ! 2nd dimension /  Constants array.
  INTEGER ::      Len1RowDepC     ! 1st dimension \  of Row Dependent
  INTEGER ::      Len2RowDepC     ! 2nd dimension /  Constants array.
  INTEGER ::      Len1ColDepC     ! 1st dimension \  of Column Dependent
  INTEGER ::      Len2ColDepC     ! 2nd dimension /  Constants array.
  INTEGER ::      Len1FldsOfC     ! 1st dimension \  of Fields of
  INTEGER ::      Len2FldsOfC     ! 2nd dimension /  Constants array.
  INTEGER ::      LenExtraC       ! Length of Extra Constants array.
  INTEGER ::      LenHistFile     ! Length of History File.
  INTEGER ::      LenCompFldI1    ! Length of Compressed Field Index 1.
  INTEGER ::      LenCompFldI2    ! Length of Compressed Field Index 2.
  INTEGER ::      LenCompFldI3    ! Length of Compressed Field Index 3.
  INTEGER ::      Len1Lookup      ! 1st dimension of Lookup table.
  INTEGER ::      Len2Lookup      ! 2nd dimension of Lookup table.
  INTEGER ::      LenData         ! Length of Data array.
  INTEGER ::      StartData       ! Position of start of Data array.
  INTEGER ::      NumFlds         ! Number of data fields.
  INTEGER ::      MaxFldSize      ! Maximum size of field.
  INTEGER ::      UnitNum         ! Unit number associated with UM dump.

  INTEGER, ALLOCATABLE :: FixHd    (:)  ! Fixed length header.
  INTEGER, ALLOCATABLE :: IntC     (:)  ! Integer Constants array.
  INTEGER, ALLOCATABLE :: CompFldI1(:)  ! Compressed Field Index array 1.
  INTEGER, ALLOCATABLE :: CompFldI2(:)  ! Compressed Field Index array 2.
  INTEGER, ALLOCATABLE :: CompFldI3(:)  ! Compressed Field Index array 3.
  INTEGER, ALLOCATABLE :: Lookup(:,:)   ! Lookup table.

  REAL, ALLOCATABLE    :: RealC  (:)    ! Real Constants array.
  REAL, ALLOCATABLE    :: LevDepC(:,:)  ! Level Dependent Constants array.

  !Use 2-dimentional Row and Column dependent consts.
  REAL, ALLOCATABLE    :: RowDepC(:,:)  ! Row Dependent Const. array.
  REAL, ALLOCATABLE    :: ColDepC(:,:)  ! Column Dependent Const. array.

  REAL, ALLOCATABLE    :: FldsOfC (:)   ! Field Dependent Constants array.
  REAL, ALLOCATABLE    :: ExtraC  (:)   ! Extra Constants array.
  REAL, ALLOCATABLE    :: HistFile(:)   ! History File.

END TYPE UM_header_type

END MODULE Rcf_UMhead_Mod

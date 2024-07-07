! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! A module to contain information about submodels

MODULE Submodel_Mod

! Description:
!   Data module to contain information about submodels. Largely
!   unnecessary but used for consistency with UM.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

!  1. Maximum internal model/submodel array sizes for this version.
! Max no. of internal models
INTEGER, PARAMETER  ::  N_Internal_Model_Max = 2

! Max no. of subm. dump parts
INTEGER, PARAMETER  ::  N_Submodel_Partition_Max = 1

! Max value of int. model id
INTEGER, PARAMETER  ::  Internal_Id_Max = N_Internal_Model_Max

! Max value of subm. dump id
INTEGER, PARAMETER  ::  Submodel_Id_Max = N_Submodel_Partition_Max


!  2. Atmos Internal Model identifiers
INTEGER, PARAMETER  ::  atmos_im = 1


! Fieldcalc (retired)
INTEGER, PARAMETER  ::  fieldcalc_im = 10

! Atmos Submodel partition identifiers
INTEGER, PARAMETER  ::  atmos_sm = 1

!  3. Lists of internal models and their submodel dump partitions -
INTEGER             :: N_Internal_Model  = 1    ! No. of internal models
                     ! For all forecast jobs this is 1 (atmos).
                     ! Not able to set as parameter due to use in utilities
                     ! and SCM.
INTEGER, SAVE  :: Internal_Model_List(N_Internal_Model_Max)

! Submodel identifier for each internal model in list
INTEGER, SAVE       :: Submodel_For_IM(N_Internal_Model_Max)

! Submodel number for each submodel id
INTEGER, SAVE       :: Submodel_For_SM(N_Internal_Model_Max)

!  4. Lists calculated in model from user interface supplied arrays -
!     - experiment specific.
! No of internal models in  each submodel partition indexed by sm identifier
INTEGER, SAVE  :: N_Internal_For_SM(Submodel_ID_Max)

! List of submodel partition identifiers
INTEGER, SAVE  :: Submodel_Partition_List(N_Submodel_Partition_Max)

! Submodel partition identifier indexed by internal model identifier
INTEGER, SAVE  :: Submodel_Partition_Index(Internal_ID_Max)

END MODULE Submodel_Mod

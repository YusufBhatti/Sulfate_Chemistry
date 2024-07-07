! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Output dump addressing arrays

MODULE Rcf_NRecon_Mod

! Description:
!   Arrays used for output dump addressing.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


USE version_mod, ONLY: &
    Nitemp,                   &
    NSectP

USE Submodel_Mod, ONLY:    &
    N_Internal_Model_Max

IMPLICIT NONE

TYPE :: Type_Recondat
  INTEGER                    :: sec_item
  INTEGER                    :: rlevs
  INTEGER                    :: LEN
  INTEGER                    :: raddress
  INTEGER                    :: rplevs
END TYPE Type_Recondat

TYPE :: Type_Recondat_Node
  TYPE ( Type_Recondat ), POINTER    :: recondat_info => NULL()
  TYPE ( Type_Recondat_Node ), POINTER :: next => NULL()
END TYPE Type_Recondat_Node

TYPE ( Type_Recondat_Node ), POINTER :: recondat_node => NULL()

TYPE ( Type_Recondat_Node ), SAVE, TARGET &
  :: RecondatList( N_Internal_Model_Max, 0 : NSectP )
INTEGER                     :: PrimDataLen ( N_Internal_Model_Max )
INTEGER                     :: DumpProgLevs( N_Internal_Model_Max )

END MODULE Rcf_NRecon_Mod

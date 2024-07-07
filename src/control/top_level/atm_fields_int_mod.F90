! "Integer" and "logical" arrays in atm_fields_mod
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE atm_fields_int_mod

! The arrays in this module are passed from atm_step into lower level routines
! (atmos_physics1/2, etc.) as arguments, and are declared as integer/logical 
! in those routines. However, they are also pointers to D1, and so have to be
! declared real in this module. At UM 10.2 and earlier they were declared real
! in atm_fields_mod but as integer/logical in typ_atm_fields.h, which was 
! #included in atm_step (while atm_fields_mod was USEd above atm_step level). 
! These arrays thus undergo a type change at the interface of atm_step. 
! As typ_atm_fields has been removed, it is therefore necessary to 
! pass these arrays as explicit arguments into atm_step - in view of the type 
! change they cannot be USEd from the module, as this causes compile errors. 
! In order to avoid the need for very long "ONLY" clauses in 
! "USE atm_fields_mod" it was decided to split atm_fields_mod into two separate
! modules, one for the real arrays and this one for the non-real arrays.

IMPLICIT NONE

! Arrays which are integer below atm_step level

! Number of model levels in the  turbulently mixed layer
REAL, POINTER :: ntml(:,:)

! Top level for turb mixing in any decoupled Sc layer
REAL, POINTER :: ntdsc(:,:)

! Bottom level for turb mixing in any decoupled Sc layer
REAL, POINTER :: nbdsc(:,:)

REAL, POINTER :: ccb(:,:)            ! Convective cloud base
REAL, POINTER :: cct(:,:)            ! Convective cloud top

! Arrays which are logical below atm_step level

! Land sea mask
REAL, POINTER :: land(:,:)
        
! Boundary layer convection flag
REAL, POINTER :: cumulus(:,:)      

END MODULE atm_fields_int_mod


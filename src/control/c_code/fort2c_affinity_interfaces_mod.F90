! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Provides Fortran2003 C-interoperable interfaces to affinity-related
! C-routines.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C Code
!

#if defined(__linux__)
MODULE fort2c_affinity_interfaces_mod

! DEPENDS ON: c_affinity.o

USE io_dependencies
IMPLICIT NONE
PRIVATE 

!-------------------------------------------------------------------------------
! C-interoperable interfaces
!-------------------------------------------------------------------------------

INTERFACE
  FUNCTION fort2c_affinity_max_available_cpus()     &
           bind(c,name="c_max_available_cpus")
    USE, INTRINSIC    :: iso_c_binding, ONLY: c_int
    INTEGER(c_int)    :: fort2c_affinity_max_available_cpus
  END FUNCTION fort2c_affinity_max_available_cpus

  FUNCTION fort2c_affinity_num_available_cpus()     &
           bind(c,name="c_num_available_cpus")
    USE, INTRINSIC    :: iso_c_binding, ONLY: c_int
    INTEGER(c_int)    :: fort2c_affinity_num_available_cpus
  END FUNCTION fort2c_affinity_num_available_cpus

  FUNCTION fort2c_affinity_running_on_core()        &
           bind(c,name="c_running_on_core")
    USE, INTRINSIC    :: iso_c_binding, ONLY: c_int
    INTEGER(c_int)    :: fort2c_affinity_running_on_core
  END FUNCTION fort2c_affinity_running_on_core
END INTERFACE 

!-------------------------------------------------------------------------------
! Public interfaces
!-------------------------------------------------------------------------------

PUBLIC :: fort2c_affinity_max_available_cpus
PUBLIC :: fort2c_affinity_num_available_cpus
PUBLIC :: fort2c_affinity_running_on_core

END MODULE fort2c_affinity_interfaces_mod
#endif


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Module to contain photolysis rate arrays
!   Contains subroutines: ukca_strat_photol_init
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE ukca_dissoc


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
IMPLICIT NONE

REAL, ALLOCATABLE, SAVE :: ajhno3 (:)
REAL, ALLOCATABLE, SAVE :: ajpna  (:)
REAL, ALLOCATABLE, SAVE :: ajh2o2 (:)
REAL, ALLOCATABLE, SAVE :: aj2a   (:)
REAL, ALLOCATABLE, SAVE :: aj2b   (:)
REAL, ALLOCATABLE, SAVE :: aj3    (:)
REAL, ALLOCATABLE, SAVE :: aj3a   (:)
REAL, ALLOCATABLE, SAVE :: ajcnita(:)
REAL, ALLOCATABLE, SAVE :: ajcnitb(:)
REAL, ALLOCATABLE, SAVE :: ajbrno3(:)
REAL, ALLOCATABLE, SAVE :: ajbrcl (:)
REAL, ALLOCATABLE, SAVE :: ajoclo (:)
REAL, ALLOCATABLE, SAVE :: ajcl2o2(:)
REAL, ALLOCATABLE, SAVE :: ajhocl (:)
REAL, ALLOCATABLE, SAVE :: ajno   (:)
REAL, ALLOCATABLE, SAVE :: ajno2  (:)
REAL, ALLOCATABLE, SAVE :: ajn2o5 (:)
REAL, ALLOCATABLE, SAVE :: ajno31 (:)
REAL, ALLOCATABLE, SAVE :: ajno32 (:)
REAL, ALLOCATABLE, SAVE :: ajbro  (:)
REAL, ALLOCATABLE, SAVE :: ajhcl  (:)
REAL, ALLOCATABLE, SAVE :: ajn2o  (:)
REAL, ALLOCATABLE, SAVE :: ajhobr (:)
REAL, ALLOCATABLE, SAVE :: ajf11  (:)
REAL, ALLOCATABLE, SAVE :: ajf12  (:)
REAL, ALLOCATABLE, SAVE :: ajh2o  (:)
REAL, ALLOCATABLE, SAVE :: ajccl4 (:)
REAL, ALLOCATABLE, SAVE :: ajf113 (:)
REAL, ALLOCATABLE, SAVE :: ajf22  (:)
REAL, ALLOCATABLE, SAVE :: ajch3cl(:)
REAL, ALLOCATABLE, SAVE :: ajc2oa (:)
REAL, ALLOCATABLE, SAVE :: ajc2ob (:)
REAL, ALLOCATABLE, SAVE :: ajmhp  (:)
REAL, ALLOCATABLE, SAVE :: ajch3br(:)
REAL, ALLOCATABLE, SAVE :: ajmcfm (:)
REAL, ALLOCATABLE, SAVE :: ajch4  (:)
REAL, ALLOCATABLE, SAVE :: ajf12b1(:)
REAL, ALLOCATABLE, SAVE :: ajf13b1(:)
REAL, ALLOCATABLE, SAVE :: ajcof2 (:)
REAL, ALLOCATABLE, SAVE :: ajcofcl(:)
REAL, ALLOCATABLE, SAVE :: ajco2  (:)
REAL, ALLOCATABLE, SAVE :: ajcos  (:)
REAL, ALLOCATABLE, SAVE :: ajhono (:)
REAL, ALLOCATABLE, SAVE :: ajmena (:)
REAL, ALLOCATABLE, SAVE :: ajchbr3(:)
REAL, ALLOCATABLE, SAVE :: ajdbrm (:)
REAL, ALLOCATABLE, SAVE :: ajcs2  (:)
REAL, ALLOCATABLE, SAVE :: ajh2so4(:)
REAL, ALLOCATABLE, SAVE :: ajso3  (:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_DISSOC'

CONTAINS

!  UKCA stratospheric photolysis module
!  Test version
! ######################################################################
!
! Subroutine Interface:
!
!---------------------------------------------------------------------------
! Subroutine UKCA_STRAT_PHOTOL_INIT
!------------------------------------------------------------------------
!
! This routine computes stratospheric photolysis rates and merges the
! rates, where necessary, with the tropospheric rates. This is done for
! one level at a time. The stratospheric photolysis routines are taken
! from SLIMCAT.

SUBROUTINE ukca_strat_photol_init
USE UM_ParVars
USE nlsizes_namelist_mod, ONLY: &
    theta_field_size

IMPLICIT NONE


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_STRAT_PHOTOL_INIT'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(ajhno3 )) ALLOCATE(ajhno3 (theta_field_size))
IF (.NOT. ALLOCATED(ajpna  )) ALLOCATE(ajpna  (theta_field_size))
IF (.NOT. ALLOCATED(ajh2o2 )) ALLOCATE(ajh2o2 (theta_field_size))
IF (.NOT. ALLOCATED(aj2a   )) ALLOCATE(aj2a   (theta_field_size))
IF (.NOT. ALLOCATED(aj2b   )) ALLOCATE(aj2b   (theta_field_size))
IF (.NOT. ALLOCATED(aj3    )) ALLOCATE(aj3    (theta_field_size))
IF (.NOT. ALLOCATED(aj3a   )) ALLOCATE(aj3a   (theta_field_size))
IF (.NOT. ALLOCATED(ajcnita)) ALLOCATE(ajcnita(theta_field_size))
IF (.NOT. ALLOCATED(ajcnitb)) ALLOCATE(ajcnitb(theta_field_size))
IF (.NOT. ALLOCATED(ajbrno3)) ALLOCATE(ajbrno3(theta_field_size))
IF (.NOT. ALLOCATED(ajbrcl )) ALLOCATE(ajbrcl (theta_field_size))
IF (.NOT. ALLOCATED(ajoclo )) ALLOCATE(ajoclo (theta_field_size))
IF (.NOT. ALLOCATED(ajcl2o2)) ALLOCATE(ajcl2o2(theta_field_size))
IF (.NOT. ALLOCATED(ajhocl )) ALLOCATE(ajhocl (theta_field_size))
IF (.NOT. ALLOCATED(ajno   )) ALLOCATE(ajno   (theta_field_size))
IF (.NOT. ALLOCATED(ajno2  )) ALLOCATE(ajno2  (theta_field_size))
IF (.NOT. ALLOCATED(ajn2o5 )) ALLOCATE(ajn2o5 (theta_field_size))
IF (.NOT. ALLOCATED(ajno31 )) ALLOCATE(ajno31 (theta_field_size))
IF (.NOT. ALLOCATED(ajno32 )) ALLOCATE(ajno32 (theta_field_size))
IF (.NOT. ALLOCATED(ajbro  )) ALLOCATE(ajbro  (theta_field_size))
IF (.NOT. ALLOCATED(ajhcl  )) ALLOCATE(ajhcl  (theta_field_size))
IF (.NOT. ALLOCATED(ajn2o  )) ALLOCATE(ajn2o  (theta_field_size))
IF (.NOT. ALLOCATED(ajhobr )) ALLOCATE(ajhobr (theta_field_size))
IF (.NOT. ALLOCATED(ajf11  )) ALLOCATE(ajf11  (theta_field_size))
IF (.NOT. ALLOCATED(ajf12  )) ALLOCATE(ajf12  (theta_field_size))
IF (.NOT. ALLOCATED(ajh2o  )) ALLOCATE(ajh2o  (theta_field_size))
IF (.NOT. ALLOCATED(ajccl4 )) ALLOCATE(ajccl4 (theta_field_size))
IF (.NOT. ALLOCATED(ajf113 )) ALLOCATE(ajf113 (theta_field_size))
IF (.NOT. ALLOCATED(ajf22  )) ALLOCATE(ajf22  (theta_field_size))
IF (.NOT. ALLOCATED(ajch3cl)) ALLOCATE(ajch3cl(theta_field_size))
IF (.NOT. ALLOCATED(ajc2oa )) ALLOCATE(ajc2oa (theta_field_size))
IF (.NOT. ALLOCATED(ajc2ob )) ALLOCATE(ajc2ob (theta_field_size))
IF (.NOT. ALLOCATED(ajmhp  )) ALLOCATE(ajmhp  (theta_field_size))
IF (.NOT. ALLOCATED(ajch3br)) ALLOCATE(ajch3br(theta_field_size))
IF (.NOT. ALLOCATED(ajmcfm )) ALLOCATE(ajmcfm (theta_field_size))
IF (.NOT. ALLOCATED(ajch4  )) ALLOCATE(ajch4  (theta_field_size))
IF (.NOT. ALLOCATED(ajf12b1)) ALLOCATE(ajf12b1(theta_field_size))
IF (.NOT. ALLOCATED(ajf13b1)) ALLOCATE(ajf13b1(theta_field_size))
IF (.NOT. ALLOCATED(ajcof2 )) ALLOCATE(ajcof2 (theta_field_size))
IF (.NOT. ALLOCATED(ajcofcl)) ALLOCATE(ajcofcl(theta_field_size))
IF (.NOT. ALLOCATED(ajco2  )) ALLOCATE(ajco2  (theta_field_size))
IF (.NOT. ALLOCATED(ajcos  )) ALLOCATE(ajcos  (theta_field_size))
IF (.NOT. ALLOCATED(ajhono )) ALLOCATE(ajhono (theta_field_size))
IF (.NOT. ALLOCATED(ajmena )) ALLOCATE(ajmena (theta_field_size))
IF (.NOT. ALLOCATED(ajchbr3)) ALLOCATE(ajchbr3(theta_field_size))
IF (.NOT. ALLOCATED(ajdbrm )) ALLOCATE(ajdbrm (theta_field_size))
IF (.NOT. ALLOCATED(ajcs2  )) ALLOCATE(ajcs2  (theta_field_size))
IF (.NOT. ALLOCATED(ajh2so4)) ALLOCATE(ajh2so4(theta_field_size))
IF (.NOT. ALLOCATED(ajso3  )) ALLOCATE(ajso3  (theta_field_size))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_strat_photol_init

END MODULE ukca_dissoc

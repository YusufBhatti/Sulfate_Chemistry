! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt 
! This file belongs in section: Top Level


MODULE atm_boundary_headers_mod

USE lbc_mod,              ONLY: rim_lookupsa, bound_lookupsa
USE nlsizes_namelist_mod, ONLY: len_fixhd, a_len_inthd, len1_lookup,  &
                                len1_lbc_comp_lookup, a_len_realhd
USE parkind1,             ONLY: jprb, jpim
USE yomhook,              ONLY: lhook, dr_hook
USE ereport_mod,          ONLY: ereport
USE umPrintMgr,           ONLY: umMessage, newline


IMPLICIT NONE

SAVE


! Headers from atmosphere boundary data sets
! Second index of header arrays = 1 Lateral boundary data
! = 2 Lower boundary data
INTEGER, ALLOCATABLE :: fixhd_bounda  (:,:)   ! Fixed header
INTEGER, ALLOCATABLE :: inthd_bounda  (:,:)   ! Integer header
REAL,    ALLOCATABLE :: realhd_bounda (:,:)   ! Real header


! Lookups for the first set of LBCs in the LBC file
INTEGER, ALLOCATABLE :: lookup_bounda (:,:)

! Varying items from the LOOKUP table for the entire LBC file
INTEGER, ALLOCATABLE :: lookup_comp_bounda(:,:)

! Control data calculated from namelist
INTEGER :: rim_stepsa      ! Set by IN_BOUND
INTEGER :: nbound_lookup(2)


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ATM_BOUNDARY_HEADERS_MOD'

CONTAINS

SUBROUTINE allocate_atm_boundary_headers()

IMPLICIT NONE

INTEGER                       :: icode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ALLOCATE_ATM_BOUNDARY_HEADERS'
CHARACTER(LEN=*), PARAMETER :: err_str = 'Failure in allocating array '

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode=0


IF (ALLOCATED (lookup_comp_bounda )) DEALLOCATE(lookup_comp_bounda)
IF (ALLOCATED (lookup_bounda )) DEALLOCATE(lookup_bounda)
IF (ALLOCATED (realhd_bounda )) DEALLOCATE(realhd_bounda)
IF (ALLOCATED (inthd_bounda )) DEALLOCATE(inthd_bounda)
IF (ALLOCATED (fixhd_bounda )) DEALLOCATE(fixhd_bounda)


ALLOCATE( fixhd_bounda(len_fixhd,1)                                            &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//":"//routinename, icode,                             &
            err_str//'[fixhd_bounda]'//                               newline//&
            TRIM(umMessage) )

ALLOCATE( inthd_bounda(a_len_inthd,1)                                          &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//":"//routinename, icode,                             &
            err_str//'[inthd_bounda]'//                               newline//&
            TRIM(umMessage) )

ALLOCATE( realhd_bounda(a_len_realhd,1)                                        &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//":"//routinename, icode,                             &
            err_str//'[realhd_bounda]'//                              newline//&
            TRIM(umMessage) )

ALLOCATE( lookup_bounda(len1_lookup,rim_lookupsa)                              &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//":"//routinename, icode,                             &
            err_str//'[lookup_bounda]'//                              newline//&
            TRIM(umMessage) )

ALLOCATE( lookup_comp_bounda(len1_lbc_comp_lookup,bound_lookupsa)              &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//":"//routinename, icode,                             &
            err_str//'[lookup_comp_bounda]'//                         newline//&
            TRIM(umMessage) )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_atm_boundary_headers

END MODULE atm_boundary_headers_mod

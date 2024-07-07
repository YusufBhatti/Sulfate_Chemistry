! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt 
! This file belongs in section: Ancillaries

MODULE ancil_headers_mod

USE ancil_mod,             ONLY: num_ancil_files
USE nlsizes_namelist_mod,  ONLY: len_fixhd, a_len_inthd, len1_lookup,  &
                                 a_len_realhd
USE ancilcta_namelist_mod, ONLY: nancil_lookupsa

USE parkind1,              ONLY: jprb, jpim
USE yomhook,               ONLY: lhook, dr_hook
USE ereport_mod,           ONLY: ereport
USE umPrintMgr,            ONLY: umMessage, newline

IMPLICIT NONE

SAVE

! Headers from ancillary datasets

INTEGER, ALLOCATABLE :: fixhd_ancila(:,:)
INTEGER, ALLOCATABLE :: inthd_ancila(:,:)
REAL,    ALLOCATABLE :: realhd_ancila(:,:)
INTEGER, ALLOCATABLE :: lookup_ancila(:,:)
INTEGER, ALLOCATABLE :: lookup_start_ancila(:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ANCIL_HEADERS_MOD'

CONTAINS

SUBROUTINE allocate_ancil_headers()

IMPLICIT NONE

INTEGER                       :: icode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ALLOCATE_ANCIL_HEADERS'
CHARACTER(LEN=*), PARAMETER :: err_str = 'Failure in allocating array '

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode=0

IF (ALLOCATED (lookup_start_ancila)) DEALLOCATE(lookup_start_ancila)
IF (ALLOCATED (lookup_ancila))       DEALLOCATE(lookup_ancila)
IF (ALLOCATED (realhd_ancila))       DEALLOCATE(realhd_ancila)
IF (ALLOCATED (inthd_ancila))        DEALLOCATE(inthd_ancila)
IF (ALLOCATED (fixhd_ancila))        DEALLOCATE(fixhd_ancila)

ALLOCATE( fixhd_ancila(len_fixhd,num_ancil_files)                      &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                        &
CALL ereport( modulename//'::'//routinename, icode,                    &
            err_str//'[fixhd_ancila]'//                       newline//&
            TRIM(umMessage) )

ALLOCATE( inthd_ancila(a_len_inthd,num_ancil_files)                    &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                        &
CALL ereport( modulename//'::'//routinename, icode,                    &
            err_str//'[inthd_ancila]'//                       newline//&
            TRIM(umMessage) )

ALLOCATE(realhd_ancila(a_len_realhd,num_ancil_files)                   &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                        &
CALL ereport( modulename//'::'//routinename, icode,                    &
            err_str//'[realhd_ancila]'//                      newline//&
            TRIM(umMessage) )

ALLOCATE( lookup_ancila(len1_lookup,nancil_lookupsa)                   &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                        &
CALL ereport( modulename//'::'//routinename, icode,                    &
            err_str//'[lookup_ancila]'//                      newline//&
            TRIM(umMessage) )

ALLOCATE( lookup_start_ancila(num_ancil_files, 2)                      &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                        &
CALL ereport( modulename//'::'//routinename, icode,                    &
            err_str//'[lookup_start_ancila]'//                newline//&
            TRIM(umMessage) )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_ancil_headers


END MODULE ancil_headers_mod

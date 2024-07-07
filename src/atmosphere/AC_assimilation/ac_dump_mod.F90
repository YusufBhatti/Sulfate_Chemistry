! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation

MODULE ac_dump_mod

USE nlsizes_namelist_mod, ONLY: len_fixhd
USE parkind1,     ONLY: jprb, jpim
USE yomhook,      ONLY: lhook, dr_hook
USE ereport_mod,  ONLY: ereport
USE umPrintMgr,   ONLY: umMessage, newline

IMPLICIT NONE

INTEGER :: len_inthd
INTEGER :: len_realhd
INTEGER :: len1_levdepc
INTEGER :: len2_levdepc
INTEGER :: len1_rowdepc
INTEGER :: len2_rowdepc
INTEGER :: len1_coldepc 
INTEGER :: len2_coldepc
INTEGER :: len1_flddepc
INTEGER :: len2_flddepc
INTEGER :: len_extcnst
INTEGER :: len_dumphist
INTEGER :: len_cfi1 
INTEGER :: len_cfi2 
INTEGER :: len_cfi3
INTEGER :: len1_lookup_obs 
INTEGER :: len2_lookup_obs

INTEGER :: fixhd (len_fixhd)

INTEGER, ALLOCATABLE :: inthd(:)
REAL, ALLOCATABLE :: realhd(:)
REAL, ALLOCATABLE :: levdepc(:,:)
REAL, ALLOCATABLE :: rowdepc(:,:)
REAL, ALLOCATABLE :: coldepc(:,:)
REAL, ALLOCATABLE :: flddepc(:,:)
REAL, ALLOCATABLE :: extcnst(:)
REAL, ALLOCATABLE :: dumphist(:)
INTEGER, ALLOCATABLE :: cfi1(:)
INTEGER, ALLOCATABLE :: cfi2(:)
INTEGER, ALLOCATABLE :: cfi3(:)
INTEGER, ALLOCATABLE :: lookup(:,:)


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='AC_DUMP_MOD'

CONTAINS

SUBROUTINE allocate_ac_dump_headers()

IMPLICIT NONE

INTEGER                       :: icode
CHARACTER(LEN=*),   PARAMETER :: err_str = 'Failure in allocating array'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*),   PARAMETER :: routinename = 'ALLOCATE_AC_DUMP_HEADERS'

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_in,zhook_handle)

CALL deallocate_ac_dump_headers()

ALLOCATE( inthd(len_inthd)                                            &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[inthd]'//                             newline//&
            TRIM(umMessage) )

ALLOCATE( realhd(len_realhd)                                          &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[realhd]'//                            newline//&
            TRIM(umMessage) )

ALLOCATE( levdepc(len1_levdepc,len2_levdepc)                          &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[levdepc]'//                           newline//&
            TRIM(umMessage) )

ALLOCATE( rowdepc(len1_rowdepc,len2_rowdepc)                          &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[rowdepc]'//                           newline//&
            TRIM(umMessage) )

ALLOCATE( coldepc(len1_coldepc,len2_coldepc)                          &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[coldepc]'//                           newline//&
            TRIM(umMessage) )

ALLOCATE( flddepc(len1_flddepc,len2_flddepc)                          &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[flddepc]'//                           newline//&
            TRIM(umMessage) )

ALLOCATE( extcnst(len_extcnst)                                        &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[extcnst]'//                           newline//&
            TRIM(umMessage) )

ALLOCATE( dumphist(len_dumphist)                                      &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[dumphist]'//                          newline//&
            TRIM(umMessage) )

ALLOCATE( cfi1(len_cfi1)                                              &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[cfi1]'//                              newline//&
            TRIM(umMessage) )

ALLOCATE( cfi2(len_cfi2)                                              &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[cfi2]'//                              newline//&
            TRIM(umMessage) )

ALLOCATE( cfi3(len_cfi3)                                              &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[cfi3]'//                              newline//&
            TRIM(umMessage) )

ALLOCATE( lookup(len1_lookup_obs,len2_lookup_obs)                     &
        , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                       &
CALL ereport( modulename//":"//routinename, icode,                    &
            err_str//'[lookup]'//                            newline//&
            TRIM(umMessage) )

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_out,zhook_handle)

END SUBROUTINE allocate_ac_dump_headers

SUBROUTINE deallocate_ac_dump_headers()

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*),   PARAMETER :: routinename = 'DEALLOCATE_AC_DUMP_HEADERS'

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_in,zhook_handle)

IF (ALLOCATED(inthd)) DEALLOCATE(inthd)
IF (ALLOCATED(realhd)) DEALLOCATE(realhd)
IF (ALLOCATED(levdepc)) DEALLOCATE(levdepc)
IF (ALLOCATED(rowdepc)) DEALLOCATE(rowdepc)
IF (ALLOCATED(coldepc)) DEALLOCATE(coldepc)
IF (ALLOCATED(flddepc)) DEALLOCATE(flddepc)
IF (ALLOCATED(extcnst)) DEALLOCATE(extcnst)
IF (ALLOCATED(dumphist)) DEALLOCATE(dumphist)
IF (ALLOCATED(cfi1)) DEALLOCATE(cfi1)
IF (ALLOCATED(cfi2)) DEALLOCATE(cfi2) 
IF (ALLOCATED(cfi3)) DEALLOCATE(cfi3)
IF (ALLOCATED(lookup)) DEALLOCATE(lookup)

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_out,zhook_handle)

END SUBROUTINE deallocate_ac_dump_headers

END MODULE ac_dump_mod

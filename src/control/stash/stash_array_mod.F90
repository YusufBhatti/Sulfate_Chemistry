! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: STASH related variables/arrays for describing output requests
!          and space management.
!
! Programming Standard : Unified Model Documentation paper number 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH
!
!
MODULE stash_array_mod

USE submodel_mod, ONLY: n_internal_model
USE version_mod,  ONLY: nelemp
USE parkind1,     ONLY: jprb, jpim
USE yomhook,      ONLY: lhook, dr_hook
USE ereport_mod,  ONLY: ereport
USE umPrintMgr,   ONLY: umMessage, newline

IMPLICIT NONE

SAVE

INTEGER, PARAMETER :: len_stlist   = nelemp ! No of items in stlist

! No of items per timeseries recd
INTEGER, PARAMETER :: time_series_rec_len = 9

INTEGER :: nsects                ! Max no of diagnostic sections
INTEGER :: n_req_items           ! Max item number in any section
INTEGER :: nitems                ! No of distinct items requested
INTEGER :: n_ppxrecs             ! No of PP_XREF records this run
INTEGER :: totitems              ! Total no of processing requests
INTEGER :: nsttims               ! Max no of STASHtimes in a table
INTEGER :: nsttabl               ! No of STASHtimes tables
INTEGER :: num_stash_levels      ! Max no of levels in a levelslist
INTEGER :: num_level_lists       ! No of levels lists
INTEGER :: num_stash_pseudo      ! Max no of pseudo-levs in a list
INTEGER :: num_pseudo_lists      ! No of pseudo-level lists
INTEGER :: nstash_series_block   ! No of blocks of timeseries recds
INTEGER :: nstash_series_records ! Total no of timeseries records


! This file is needed to get ppxref_codelen to dimension PP_XREF
! sizes in STASH used for defining local array dimensions at a
! lower level.
INTEGER :: max_stash_levs  ! Max no of output levels for any diag
INTEGER :: pp_len2_lookup  ! Max no of lookups needed in stwork


! stashflag (.TRUE. for processing this timestep). sf(0,is) .FALSE.
! if no flags on for section is.
LOGICAL, ALLOCATABLE :: sf (:,:)

! Whether a calculation is needed for sf above
LOGICAL, ALLOCATABLE :: sf_calc (:,:)

! STASH list index
INTEGER, ALLOCATABLE :: stindex (:,:,:,:)

! List of STASH output requests
INTEGER, ALLOCATABLE :: stlist (:,:)

! Address of item from generating plug compatible routine (often workspace)
INTEGER, ALLOCATABLE :: si (:,:,:)

! STASH times tables
INTEGER, ALLOCATABLE :: sttabl (:,:)

! Length of STASH workspace required in each section
INTEGER, ALLOCATABLE :: stash_maxlen        (:,:)
INTEGER, ALLOCATABLE :: ppindex             (:,:)
INTEGER, ALLOCATABLE :: stash_levels        (:,:)
INTEGER, ALLOCATABLE :: stash_pseudo_levels (:,:)
INTEGER, ALLOCATABLE :: stash_series        (:,:)
INTEGER, ALLOCATABLE :: stash_series_index  (:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'STASH_ARRAY_MOD'

PRIVATE :: n_internal_model

CONTAINS

SUBROUTINE allocate_stash_arrays()

! Purpose: Allocate STASH arrays based on current values of
!          the dimensioning variables!

IMPLICIT NONE

INTEGER                       :: icode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ALLOCATE_STASH_ARRAYS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode = 0

IF (ALLOCATED( sf                  )) DEALLOCATE(sf)
IF (ALLOCATED( sf_calc             )) DEALLOCATE(sf_calc)
IF (ALLOCATED( stindex             )) DEALLOCATE(stindex)
IF (ALLOCATED( stlist              )) DEALLOCATE(stlist)
IF (ALLOCATED( si                  )) DEALLOCATE(si)
IF (ALLOCATED( sttabl              )) DEALLOCATE(sttabl)
IF (ALLOCATED( stash_maxlen        )) DEALLOCATE(stash_maxlen)
IF (ALLOCATED( ppindex             )) DEALLOCATE(ppindex)
IF (ALLOCATED( stash_levels        )) DEALLOCATE(stash_levels)
IF (ALLOCATED( stash_pseudo_levels )) DEALLOCATE(stash_pseudo_levels)
IF (ALLOCATED( stash_series        )) DEALLOCATE(stash_series)
IF (ALLOCATED( stash_series_index  )) DEALLOCATE(stash_series_index)


ALLOCATE( sf (0:nitems, 0:nsects)                                              &
        , stat=icode, errmsg=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [sf]'//                    newline//&
              TRIM(umMessage) )


ALLOCATE( sf_calc (0:nitems, 0:nsects)                                         &
        , stat=icode, errmsg=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [sf_calc]'//               newline//&
              TRIM(umMessage) )


ALLOCATE( stindex (2, nitems, 0:nsects, n_internal_model)                      &
        , stat=icode, errmsg=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [sf_calc]'//               newline//&
              TRIM(umMessage) )


ALLOCATE( stlist (len_stlist, totitems)                                        &
        , stat=icode, errmsg=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stlist]'//                newline//&
              TRIM(umMessage) )


ALLOCATE( si (nitems, 0:nsects, n_internal_model)                              &
        , stat=icode, errmsg=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [si]'//                    newline//&
              TRIM(umMessage) )


ALLOCATE( sttabl (nsttims, nsttabl)                                            &
        , stat=icode, errmsg=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [sttabl]'//                newline//&
              TRIM(umMessage) )


ALLOCATE( stash_maxlen (0:nsects, n_internal_model)                            &
         , stat=icode, errmsg=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stash_maxlen]'//          newline//&
              TRIM(umMessage) )


ALLOCATE( ppindex (nitems, n_internal_model)                                   &
        , stat=icode, errmsg=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [ppindex]'//               newline//&
              TRIM(umMessage) )


ALLOCATE( stash_levels (num_stash_levels+1, num_level_lists)                   &
        , stat=icode, errmsg=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stash_levels]'//          newline//&
              TRIM(umMessage) )


ALLOCATE( stash_pseudo_levels (num_stash_pseudo+1, num_pseudo_lists)           &
        , stat=icode, errmsg=umMessage)
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stash_pseudo_levels]'//   newline//&
              TRIM(umMessage) )


ALLOCATE( stash_series (time_series_rec_len, nstash_series_records)            &
        , stat=icode, errmsg=umMessage)
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stash_series]'//          newline//&
              TRIM(umMessage) )


ALLOCATE( stash_series_index (2, nstash_series_block)                          &
        , stat=icode, errmsg=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [stash_series_index]'//    newline//&
              TRIM(umMessage) )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_stash_arrays

END MODULE stash_array_mod

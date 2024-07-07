! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Control routine for processing of basis library STASH file
!
! Subroutine Interface:

SUBROUTINE stash_proc(ErrorStatus, cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE ppxlook_mod

USE version_mod, ONLY:                                            &
    nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
    nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
    nlevp, npslevp, npslistp

USE submodel_mod, ONLY: n_internal_model_max, n_submodel_partition_max
USE stextend_mod, ONLY: llistty, indx_s, in_s, ppind_s,     &
                        npos_ts, nrecs_ts, list_s, itim_s,  &
                        rlevlst_s, levlst_s, lenplst
USE cstash_mod, ONLY: pslist_d, npslists,  &
                      nseries

USE h_vers_mod, ONLY: h_vers

USE missing_data_mod, ONLY: imdi, rmdi

USE errormessagelength_mod, ONLY: errormessagelength

USE stash_model_mod, ONLY:                                                     &
    nserrec_s, nserblk_s, len_prim, len_dump, global_len_prim, global_len_dump,&
    len_secd, len_extra, len_work, nheadsub, nhead, len_primim, len_dumpim,    &
    global_len_primim, global_len_dumpim, len_secdim, nrecs_s, ntimes_s,       &
    nlevl_s, item_max_all, nmaxlev_s, npslists_s, nmaxpsl_s

IMPLICIT NONE
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!

! Subroutine arguments

!   Array arguments with intent(out):
CHARACTER(LEN=errormessagelength) :: cmessage    ! Error return message

!   Error status:
INTEGER ::     ErrorStatus ! +ve = fatal error

! Local scalars
INTEGER :: i,j,l,ii
INTEGER :: nlevels
INTEGER :: nrecs
INTEGER :: ntimes

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STASH_PROC'

!- End of Header ---------------------------------------------------

!Initialisation
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
npslists    =0   ! Counter for no. of pseudo levels lists
nseries     =0   ! Time series block counter
nserrec_s   =0   ! Total no. of time series records
nserblk_s   =0   ! Total no. of time series blocks
ErrorStatus =0

DO i=1,nprofdp
  nrecs_ts(i)=0
  npos_ts(i)=0
END DO

! Initialisation of data length arrays
DO i=1,n_submodel_partition_max
  len_prim(i)        = 0 ! Length of primary data
  len_dump(i)        = 0 ! Length of dump extension (diagnostic)
  global_len_prim(i) = 0 ! Length of global primary data
  global_len_dump(i) = 0 ! Length of global dump extension
  len_secd(i)        = 0 ! Length of secondary atmos
  len_extra(i)       = 0 ! Length of space that is addressed in model
  len_work(i)        = 0 ! Length of work
  NHeadSub(i)        = 0 ! No. of pp headers for each submodel
END DO
DO i=1,n_internal_model_max
  nhead(i)             = 0 ! No. of pp headers for each internal model
  len_primim(i)        = 0 ! Primary data length for each internal model
  len_dumpim(i)        = 0 ! Diagnostic do.
  global_len_primim(i) = 0 ! Global primary data length for each
                           ! internal model
  global_len_dumpim(i) = 0 ! Dlobal diagnostic do.
  len_secdim(i)        = 0 ! Secondary do.
END DO

DO i=1,nprofdp
  llistty(i)=' '
END DO


DO j=1,nlevp_s
  DO i=1,nlevlstsp
    rlevlst_s(j,i)=rmdi
    levlst_s(j,i)=imdi
  END DO
END DO

DO l=1,n_internal_model_max
  DO j=0,nsectp
    DO i=1,nitemp
      in_s(1,l,j,i)=0
      in_s(2,l,j,i)=0
      indx_s(1,l,j,i)=0
      indx_s(2,l,j,i)=0
    END DO
  END DO
END DO

DO j=1,nelemp+1
  DO i=1,nrecdp
    list_s(j,i)=0
  END DO
END DO

DO j=1,  ntimep
  DO i=1,2*nproftp+2
    itim_s(j,i)=-1
  END DO
END DO

DO j=1,n_internal_model_max
  DO i=1,nitemp
    ppind_s(j,i)=0
  END DO
END DO

DO i=1,nprofdp
  nrecs_ts(i)=0
  npos_ts(i)=0
END DO

DO i=1,npslistp
  lenplst(i)=0
END DO

DO i=1,n_internal_model_max
  DO j=0,nsectp
    h_vers(i,j)=0
  END DO
END DO

! Read stash basis file from job library
! DEPENDS ON: rdbasis
CALL rdbasis(cmessage,ErrorStatus)
! Adjust stash time series records (if any)
IF (nseries >  0) THEN
  ! DEPENDS ON: timser
  CALL timser(cmessage,ErrorStatus)
END IF
IF (ErrorStatus  /=  0) GO TO 9999


CALL read_atmos_stashmaster()

! Define submodel and section/version configuration
! DEPENDS ON: setmodl
CALL setmodl(ErrorStatus,cmessage)
IF (ErrorStatus  /=  0) GO TO 9999

nrecs=0
ntimes=0
nlevels=0
DO i=1,npslistp
  lenplst(i)=0
END DO


! Construct preliminary STASH list
! DEPENDS ON: prelim
CALL prelim(nrecs,                                              &
            ntimes,nlevels,ErrorStatus,cmessage)
IF (ErrorStatus >  0) GO TO 9999

! REORDER STASH LIST & SET UP INDEX

! DEPENDS ON: order
CALL order(nrecs)
! DEPENDS ON: sindx
CALL sindx(nrecs)

! DELETE DUPLIC ENTRIES, CONCATENATE OVERLAP LEVELS ETC
! DELETE DUPLICATE STASH_TIMES, REORDER

! DEPENDS ON: duplic
CALL duplic(nrecs,ntimes,nlevels)

! REORDER STASH LIST & SET UP INDEX

! DEPENDS ON: order
CALL order(nrecs)
! DEPENDS ON: sindx
CALL sindx(nrecs)

! CHANGE POINTER SYSTEM, ADD ADDRESSES, LENGTHS AND INPUT LEVELS
! DEPENDS ON: pointr
CALL pointr(nrecs)

! OUTPUT LENGTH
! DEPENDS ON: outptl
CALL outptl(                                                    &
            nrecs,ErrorStatus,cmessage)
IF (ErrorStatus  /=  0) GO TO 9999

! INPUT LENGTH AND INPUT LEVELS, SET STLIST(NELEMP+1,I) TO MODEL_ST
! ALSO INPUT PSEUDO LEVELS
! DEPENDS ON: inputl
CALL inputl(nrecs,                                              &
            nlevels,ErrorStatus,cmessage)
IF (ErrorStatus  /=  0) GO TO 9999


!     ADDRESSING
! DEPENDS ON: addres
CALL addres(                                                      &
            nrecs,ErrorStatus,cmessage)
IF (ErrorStatus  /=  0) GO TO 9999

! SET RETURN VALUES FOR OTHER FILES & write out STASH list

nrecs_s=nrecs
ntimes_s=ntimes
nlevl_s=nlevels
item_max_all=nitemp
!       ITEM_MAX_REQ IS DONE IN WSTLIST
nmaxlev_s=1
DO i =1,nlevels
  nmaxlev_s=MAX(nmaxlev_s,levlst_s(1,i))
END DO

npslists_s=npslists
nmaxpsl_s=1
DO i =1,npslists
  nmaxpsl_s=MAX(nmaxpsl_s,lenplst(i))
END DO

! Write output file (for checking purposes).

! DEPENDS ON: wstlst
CALL wstlst(nrecs,ntimes,nlevels)

9999  CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stash_proc

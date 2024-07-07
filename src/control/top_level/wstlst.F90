! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:
SUBROUTINE wstlst(nrecs,ntimes,nlevels)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr
USE UM_ParParams
USE UM_ParCore, ONLY: mype

USE ppxlook_mod, ONLY: ppxref_sections
USE version_mod, ONLY: nitemp, nelemp, nsectp, ntimep


USE stextend_mod, ONLY: llistty, indx_s, in_s, ppind_s,     &
                        npos_ts, nrecs_ts, list_s, itim_s,  &
                        rlevlst_s, levlst_s, lenplst
USE cstash_mod, ONLY: pslist_d, npslists, tlim_ts, blim_ts,  &
                      wlim_ts, nlim_ts, ndprof, elim_ts, slim_ts
USE submodel_mod, ONLY: n_internal_model_max, atmos_sm, atmos_im
USE stash_array_mod, ONLY:                                                     &
    n_req_items, num_pseudo_lists, totitems, num_stash_pseudo, nitems,         &
    nstash_series_records, num_level_lists, nstash_series_block,               &
    num_stash_levels, nsttims, nsttabl, nsects, n_ppxrecs

USE file_manager, ONLY: get_file_unit_by_id, init_file_loop, um_file_type

USE nlsizes_namelist_mod, ONLY: &
    a_len2_lookup, a_len_d1, a_len_data, global_a_len_data, len_tot

USE stash_model_mod, ONLY:                                                     &
    nserblk_s, nserrec_s, item_max_req, nhead, len_primim, len_dumpim,         &
    global_len_prim, global_len_primim, global_len_dumpim, len_secdim,         &
    len_extra, nmaxlev_s, nlevl_s, nmaxpsl_s, npslists_s

IMPLICIT NONE
! Description: Print STASH control arrays [but also initialises some
!              STASH control variables]
!
! Method:  Simple interception of control arrays for printing.
!          [ Note that modularity would be improved if the function of
!            simply printing out variables was separated from the
!            initialisation of a no. of derived STASH variable
!            interspersed through the routine.]
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
! Global variables:


INTEGER :: number_times
INTEGER :: i,j,k
INTEGER :: ft_unit
INTEGER :: lend1a
INTEGER :: nlevels
INTEGER :: nrecs
INTEGER :: ntimes
INTEGER :: BlkId           !Time series block identifier
INTEGER :: BlkSt           !Start position of ts block data
INTEGER :: idp             !Domain profile loop counter
INTEGER :: ipos
INTEGER :: Nrecs_prev      !No of recs in previous time ser block

TYPE(um_file_type), POINTER :: pp_file
INTEGER :: stash_unit    ! Stash requests log file id
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='WSTLST'

NAMELIST/stsizes/                                                 &
        a_len2_lookup,a_len_data,a_len_d1,                        &
        len_tot,                                                  &
        nsects,n_req_items,nitems,totitems,                       &
        nsttabl,num_stash_levels,num_level_lists,                 &
        num_stash_pseudo,num_pseudo_lists,                        &
        nsttims,nstash_series_block,                              &
        nstash_series_records,n_ppxrecs
!
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
lend1a=0

IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN

  stash_unit = get_file_unit_by_id("stash_req_log", handler="fortran")

  WRITE(stash_unit,'(I5,A)') nrecs, ' STASH LIST RECORDS'
  DO j=1,nrecs
    WRITE(stash_unit,'(A)') ' '
    WRITE(stash_unit,'(8I10)')(list_s(i,j),i=1,nelemp)
  END DO

  WRITE(stash_unit,'(I4,A)') ntimes,' STASH TIMES'
  DO i=1,ntimes
    number_times=0
    DO j=1,ntimep
      IF (itim_s(j,i) >  0) THEN
        number_times=number_times+1
      END IF
    END DO
    WRITE(stash_unit,'(100I5)') (itim_s(j,i),j=1,number_times)
  END DO

  WRITE(stash_unit,'(I4,A)') nlevels, ' STASH LEVELS LIST(S)'
  DO i=1,nlevels
    WRITE(stash_unit,'(I5,A1)') levlst_s(1,i),llistty(i)
    IF (llistty(i) == 'I') THEN
      WRITE(stash_unit,'(16I4)') (levlst_s(j,i),j=2,levlst_s(1,i)+1)
    ELSE
      WRITE(stash_unit,'(6F12.3)') (rlevlst_s(j,i),j=2,levlst_s(1,i)+1)
    END IF
  END DO


  WRITE(stash_unit,'(I4,A)') npslists, ' STASH PSEUDO DIMENSION LIST(S)'
  DO i=1,npslists
    WRITE(stash_unit,'(I5)') lenplst(i)
    WRITE(stash_unit,'(16I4)') (pslist_d(j,i),j=1,lenplst(i))
  END DO
  !
  WRITE(stash_unit,'(I4,A)') nserblk_s, ' STASH TIME SERIES BLOCKS'
  !
  BlkSt =1
  DO idp=1,ndprof
    IF (npos_ts(idp) >  0) THEN
      BlkId = npos_ts (idp)
      IF (BlkId >  1) THEN
        BlkSt=BlkSt+Nrecs_prev
      END IF
      WRITE(stash_unit,'(A,I4,A,I4,A)') 'SERIES NUMBER', npos_ts(idp),     &
                              ' WITH ', nrecs_ts(npos_ts(idp)), ' RECORDS'
      WRITE(stash_unit,'(A)')                                              &
          '   NORTH   SOUTH    EAST    WEST  BOTTOM     TOP'
      DO ipos=BlkSt,BlkSt+nrecs_ts(npos_ts(idp))-1
        WRITE(stash_unit,'(6I8)') nlim_ts(ipos),slim_ts(ipos),             &
            elim_ts(ipos), wlim_ts(ipos),blim_ts(ipos),tlim_ts(ipos)
      END DO
      WRITE(stash_unit,'(A)') ' '
      Nrecs_prev=nrecs_ts(npos_ts(idp)) ! For next TS block
    END IF
  END DO

  WRITE(stash_unit,'(A,A)')                                                &
  ' MODL SECT ITEM   IN_S(1)   IN_S(2) INDX_S(1) INDX_S(2)',               &
  '   PPIND_S'
  WRITE(stash_unit,'(A)')                                                  &
  '                  St addr   St  len StListPos StListNum '

END IF ! on PrintStatus and mype=0


n_ppxrecs=0
item_max_req=1

DO k=1,n_internal_model_max
  DO j=0,ppxref_sections
    DO i=1,nitemp
      IF (in_s(1,k,j,i) /= 0) THEN
        n_ppxrecs=n_ppxrecs+1
        item_max_req=MAX(j,item_max_req)

        IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
          IF (j == 0) THEN
            WRITE(stash_unit,'(3I5,5I10)')                               &
              k,j,i,  in_s(1,k,j,i),  in_s(2,k,j,i),                     &
                    indx_s(1,k,j,i),indx_s(2,k,j,i),ppind_s(k,i)
          ELSE
            WRITE(stash_unit,'(3I5,5I10)')                               &
              k,j,i,  in_s(1,k,j,i),  in_s(2,k,j,i),                     &
                    indx_s(1,k,j,i),indx_s(2,k,j,i)
          END IF
        END IF ! on PrintStatus and mype=0

      END IF
    END DO ! over i
  END DO ! over j
END DO ! over k

IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  i=-1
  WRITE(stash_unit,'(3I5,5I10)') i,i,i,i,i,i
END IF ! on PrintStatus and mype=0

! Variables in COMMON STSIZES - for UMINDEX routine.
!          N_PPXRECS was obtained in the loop above.

a_len2_lookup = nhead(atmos_im)
a_len_data    = MAX( 1, len_primim(atmos_im) + len_dumpim(atmos_im) )
lend1a        = len_primim(atmos_im) + len_dumpim(atmos_im) +  &
                                      len_secdim(atmos_im) + len_extra(atmos_sm)
len_tot       = lend1a
a_len_d1      = lend1a

global_a_len_data = MAX( 1, global_len_primim(atmos_im) +    &
                                                   global_len_dumpim(atmos_im) )

nsects       =nsectp
n_req_items  =item_max_req
nitems       =nitemp
totitems     =MAX(1,nrecs)
nsttabl      =MAX(1,ntimes)

num_stash_levels=MAX(1,nmaxlev_s)
num_level_lists =MAX(1,nlevl_s)
num_stash_pseudo=MAX(1,nmaxpsl_s)
num_pseudo_lists=MAX(1,npslists_s)
nsttims         =ntimep

nstash_series_block =MAX(1,nserblk_s)
nstash_series_records=MAX(1,nserrec_s)

!Assign number of reserved fields, and designate output files
NULLIFY(pp_file)
pp_file => init_file_loop(handler="portio")
DO WHILE (ASSOCIATED(pp_file))
  IF (pp_file % pp_meta % reserved_headers < &
      MAX(4096, pp_file % pp_meta % reserved_headers_calc)) THEN
    ! Previously a loop over nrecs below would go through the output
    ! stream units again and set them all back to 4096.  Since the
    ! comparison above would previously take place *after* this point
    ! we instead have to ignore the results of the calculation here
    ! and stick with the user set value or 4096
    pp_file % pp_meta % reserved_headers = &
        MAX(4096, pp_file % pp_meta % reserved_headers)

    ! The line below would replicate the (intended?) behaviour and
    ! have the results of the calculation take precedence over the
    ! user supplied value
    !pp_file % pp_meta % reserved_headers = &
    !    MAX(4096, pp_file % pp_meta % reserved_headers_calc)
  END IF

  ! Increment loop to the next file
  pp_file => pp_file % next

END DO

IF (PrintStatus >= PrStatus_Oper .AND. mype==0) THEN
  WRITE(stash_unit,stsizes)
END IF ! on PrintStatus and mype=0

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE wstlst

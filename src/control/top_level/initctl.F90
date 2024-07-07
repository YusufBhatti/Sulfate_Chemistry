! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE INITCTL-------------------------------------------------
!
!    PROGRAMMING STANDARD: UNIFIED MODEL DP NO. 3, VERSION 3
!
!    SYSTEM TASK: C4
!
!    SYSTEM COMPONENTS: C30, C40
!
!    PURPOSE: Initialises STASH control arrays from STASH control file.
!
!    EXTERNAL DOCUMENTATION: UMDP C4
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

SUBROUTINE initctl(                                               &
                 icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE filenamelength_mod, ONLY:                                     &
    filenamelength
USE version_mod, ONLY: nprofdp, npslevp
USE stextend_mod, ONLY: llistty, indx_s, in_s, ppind_s, levlst_s, &
           rlevlst_s, npos_ts, nrecs_ts, list_s, itim_s, lenplst
USE umPrintMgr
USE UM_ParParams
USE lookup_addresses
USE submodel_mod, ONLY:                                                        &
    internal_model_list, n_internal_model_max, n_internal_model

USE stash_array_mod, ONLY:                                                     &
    len_stlist, num_pseudo_lists, stash_pseudo_levels, pp_len2_lookup,         &
    totitems, stash_series, num_stash_pseudo, stash_maxlen, ppindex, nitems,   &
    stindex, nstash_series_records, stlist, num_level_lists,                   &
    nstash_series_block, num_stash_levels, nsttims, stash_series_index,        &
    max_stash_levs, stash_levels, nsttabl, si, sttabl, nsects,                 &
    time_series_rec_len, n_ppxrecs
USE stparam_mod, ONLY: st_input_code, st_output_code, st_pseudo_in,&
                    st_input_bottom, st_input_top, st_special_code,&
                    st_output_type, st_fieldsfile

USE ppxlook_mod, ONLY: exppxc, exppxi
USE cppxref_mod, ONLY:                                            &
    ppxref_codelen, ppx_model_number,                             &
    ppx_data_type, ppx_grid_type,                                 &
    ppx_field_code, ppx_cf_levelcode, ppx_cf_fieldcode,           &
    ppx_lv_code, ppx_pack_acc
USE version_mod, ONLY:                                            &
    nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
    nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
    nlevp, npslevp, npslistp
USE cstash_mod, ONLY: elim_ts, wlim_ts, slim_ts, nlim_ts, blim_ts,&
                      tlim_ts, ndprof, nseries, pslist_d, ig_ts,  &
                      i51_ts, i1_ts

USE file_manager, ONLY: get_file_by_unit, um_file_type

USE nlsizes_namelist_mod, ONLY: &
    len_tot, n_obj_d1_max

USE errormessagelength_mod, ONLY: errormessagelength

USE diagdesc_mod, ONLY: diagdesc

IMPLICIT NONE

!  Arguments


INTEGER ::                                                        &
   icode                  ! OUT: Error return code
!
CHARACTER(LEN=errormessagelength) ::                              &
   cmessage               ! OUT: Error return message

!  Local arrays

!  STASH input lengths
INTEGER :: si_len(nitems,0:nsects,n_internal_model_max)
INTEGER :: ppxref_dat(ppxref_codelen)

! Local variables

CHARACTER(LEN=36) :: NAME
REAL ::                                                           &
      real_levels(num_stash_levels,num_level_lists)

INTEGER ::                                                        &
       num_lists,                                                &
       num_levels,                                               &
       num_pseudo_levels,                                        &
       n_tables,                                                 &
       i,                                                        &
       ipk,                                                      &
       kk,                                                       &
       l,                                                        &
       IS,                                                       &
       ie,                                                       &
       ii,                                                       &
       sm,                                                       &
       iobj,                                                     &
       isec,                                                     &
       itm,                                                      &
       Im_ident,                                                 &
       Sm_ident,                                                 &
       ix,                                                       &
       istep,                                                    &
       il,                                                       &
       im,                                                       &
       jj,                                                       &
       input_length

INTEGER :: ir1               ! loop start
INTEGER :: ir,j,k            ! loop count
CHARACTER(LEN=1) :: var_type
CHARACTER(LEN=filenamelength) :: filename
LOGICAL :: int_mod_included  ! Flag to indicate whether a particular
                          !   internal model is included
INTEGER :: pp_len2_look
TYPE(um_file_type), POINTER :: pp_file

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INITCTL'

!-----------------------------------------------------------------------

!  1. Assign STASHlist and associated lists to appropriate UM arrays

!     Initialise STLIST to zero

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO   ii = 1,totitems
  DO ie = 1,len_stlist
    stlist(ie,ii)=0
  END DO
END DO

!     Assign STASH list to STLIST

DO   i = 1,totitems
  DO j = 1,len_stlist
    stlist(j,i) = list_s(j,i)
  END DO
END DO

!     Assign STASH times tables to STTABL

DO   i = 1,nsttabl
  DO j = 1,nsttims
    sttabl(j,i) = itim_s(j,i)
  END DO
END DO

!     Assign STASH levels lists to STASH_LEVELS

DO   ii=1,num_level_lists      ! Initialise STASH_LEVELS to -99
  DO jj=1,num_stash_levels+1
    stash_levels(jj,ii)=-99
  END DO
END DO

DO   i = 1,num_level_lists
  num_levels        = levlst_s(1,i)
  stash_levels(1,i) = num_levels
  DO j = 1,num_levels
    IF (llistty(i) == 'R') THEN
      real_levels (j  ,i) = rlevlst_s  (j+1,i)
      IF (real_levels (j  ,i) >= 0.0) THEN
        stash_levels(j+1,i) =(real_levels(j  ,i)+0.0001)*1000.0
      ELSE
        stash_levels(j+1,i) =(real_levels(j  ,i)-0.0001)*1000.0
      END IF
    ELSE IF (llistty(i) == 'I') THEN
      stash_levels(j+1,i) = levlst_s   (j+1,i)
    END IF
  END DO
END DO

!     Store STASH pseudo levels lists in STASH_PSEUDO_LEVELS

DO   ii=1,num_pseudo_lists       ! Initialise STASH_PSEUDO_LEVELS
  DO jj=1,num_stash_pseudo+1     !   to -99
    stash_pseudo_levels(jj,ii)=-99
  END DO
END DO

DO   i = 1,num_pseudo_lists
  num_levels = lenplst(i)
  stash_pseudo_levels(1,i) = num_levels
  DO j = 1,num_levels
    stash_pseudo_levels(j+1,i) = pslist_d(j,i)
  END DO
END DO

!Transfer time series data to STASH_SERIES array
IF (nseries >  0) THEN
  !There are timeseries domains
  !  Loop over STASHC domain profiles
  DO i=1,ndprof
    IF (npos_ts(i)  >   0) THEN
      !  This domain profile has a block of time series domains
      !    J=time series block identifier (pointer):
      j=npos_ts(i)
      !    STASH_SERIES_INDEX(1,J)=sequence no. of first record for
      !                            ts block J in STASH_SERIES
      IF (j == 1) THEN
        stash_series_index(1,j)=1
      ELSE
        stash_series_index(1,j)= stash_series_index(1,j-1)        &
                               +stash_series_index(2,j-1)
      END IF
      !    STASH_SERIES_INDEX(2,J)=no. of records in ts block J
      stash_series_index(2,j)=nrecs_ts(j)
      ir1  =stash_series_index(1,j)
      DO ir=ir1,ir1+nrecs_ts(j)-1
        stash_series(1,ir)=ig_ts
        stash_series(2,ir)=i1_ts
        stash_series(3,ir)=i51_ts
        stash_series(4,ir)=nlim_ts(ir)
        stash_series(5,ir)=slim_ts(ir)
        stash_series(6,ir)=wlim_ts(ir)
        stash_series(7,ir)=elim_ts(ir)
        stash_series(8,ir)=blim_ts(ir)
        stash_series(9,ir)=tlim_ts(ir)
      END DO
    END IF
  END DO
END IF

!     Initialise STINDEX and SI

DO ie=1,nitems
  DO IS=0,nsects
    DO im=1,n_internal_model
      DO ii=1,2
        stindex(ii,ie,IS,im)=0
      END DO
      si    (ie,IS,im)=1
      si_len(ie,IS,im)=0
      IF (IS == 0) THEN
        ppindex(ie,im)=0
      END IF
    END DO
  END DO
END DO

! 2. Read STASHindex and compute STASHWORK array lengths.
!    The Lth. row in STINDEX, SI, SI_LEN, PPINDEX, corresponds to the
!        Lth. internal model in INTERNAL_MODEL_LIST.
!    Output a formatted description of the selected diagnostics.

ii  =0       ! Counter for checking no. of diags. printed
l   =0       ! Counter for rows in STINDEX, SI, etc.
DO k=1,n_internal_model_max
  int_mod_included=.FALSE.
  ! Find out whether int. model  'K' is included. If it is:
  !  Set logical flag, increment row number
  DO kk=1,n_internal_model_max
    IF (internal_model_list(kk) == k) THEN
      int_mod_included=.TRUE.
      l = l + 1
    END IF
  END DO
  IF (int_mod_included) THEN
    DO j=0,nsects                        ! NSECTS=NSECTP (WSTLST)
      DO i=1,nitems                        ! NITEMS=NITEMP (WSTLST)
        IF (in_s(1,k,j,i) >= 1) THEN        ! Entry in STASH list
          ii = ii + 1
          stindex(1,i,j,l) = indx_s (1,k,j,i)  ! STASH index
          stindex(2,i,j,l) = indx_s (2,k,j,i)
          si     (  i,j,l) = in_s   (1,k,j,i)  ! STASH lengths and
          si_len (  i,j,l) = in_s   (2,k,j,i)  !   addresses in D1
          IF (j == 0 .AND. ppind_s(k,i) /= 0) THEN
            ppindex(i,  l) = ppind_s(  k,  i)  ! Index for pp header
          END IF                               !   array
          ! Extract ppxref information to be passed into
          !   diagnostic description routine
          NAME = exppxc(k,j,i,                                          &
                       icode,cmessage)
          ppxref_dat(ppx_model_number) = k
          ppxref_dat(ppx_field_code  ) = exppxi(k,j,i,ppx_field_code  , &
                                               icode,cmessage)
          ppxref_dat(ppx_data_type   ) = exppxi(k,j,i,ppx_data_type   , &
                                               icode,cmessage)
          ppxref_dat(ppx_grid_type   ) = exppxi(k,j,i,ppx_grid_type   , &
                                               icode,cmessage)
          ppxref_dat(ppx_lv_code     ) = exppxi(k,j,i,ppx_lv_code     , &
                                               icode,cmessage)
          ppxref_dat(ppx_cf_levelcode) = exppxi(k,j,i,ppx_cf_levelcode, &
                                               icode,cmessage)
          ppxref_dat(ppx_cf_fieldcode) = exppxi(k,j,i,ppx_cf_fieldcode, &
                                               icode,cmessage)
          DO ipk = 0,9
            ppxref_dat(ppx_pack_acc+ipk) = exppxi(k,j,i,ppx_pack_acc+ipk, &
                                                 icode,cmessage)
          END DO

          IF (PrintStatus >= PrStatus_Normal) THEN
            !  Write a formatted description of the diagnostic to output file.
            DO ix=stindex(1,i,j,l),                                       &
                                                             ! Loop over
                 stindex(1,i,j,l)+stindex(2,i,j,l)-1        !   entries
              CALL diagdesc(ix,NAME,stlist(1,ix),ppxref_dat(1),             &
             stash_levels,num_stash_levels,num_level_lists,                &
             stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,        &
             sttabl,nsttims,nsttabl,                                       &
             stash_series,time_series_rec_len,nstash_series_records,       &
             stash_series_index,nstash_series_block)
            END DO
          END IF  ! PrintStatus test

        ELSE IF (in_s(1,k,j,i) == -1) THEN
          ii = ii + 1
        END IF                       !  INDX_S(1,K,J,I) >= 1
        IF (icode >  0) GO TO 9999
      END DO                       !  Items
    END DO                       !  Sections

  END IF                       !  INT_MOD_INCLUDED
END DO                       !  Models

IF (ii /= n_ppxrecs) THEN
  WRITE(umMessage,*) ' Error in INITCTL: N_PPXRECS not correct  ',ii
  CALL umPrint(umMessage,src='initctl')
  cmessage='INITCTL  : N_PPXRECS not correct               '
  icode=1
  GO TO 9999
END IF

!  2.1 Find the max length in STASH_WORK and store in STASH_MAXLEN

DO im=1,n_internal_model
  DO IS=1,nsects  !  Note not section Zero as the data is in D1
    stash_maxlen(IS,im)=1
    DO ie=1,nitems  ! Again only data not in D1
      IF (stindex(1,ie,IS,im) /= 0) THEN
        !         Item is in STASHlist ...
        IF (stlist(st_input_code,stindex(1,ie,IS,im)) == 1) THEN
          !           ...  and input from STASHwork
          !             ... input length not from ST_LIST as this is post STOCGT
          input_length=si_len(ie,IS,im)
          stash_maxlen(IS,im)=stash_maxlen(IS,im)+input_length
        END IF
      END IF
    END DO
  END DO
END DO

!
!
! ----------------------------------------------------------------------
!   3.   Set derived control variables for use in STASH/STWORK
!
!        Set pp_len2_lookUP to maximum pp_len2_look value for any PP
!        unit referenced in the STASHlist (minimum value possible is 8).
!        Set MAX_STASH_LEVS to the maximum possible no of output levels
!        for any diagnostic, allowing for possible pseudo-levels.
!
pp_len2_lookUP=8
max_stash_levs=1
DO ii=1,totitems
  IF (stlist(st_output_code,ii) <  0 .AND. &
      stlist(st_output_type,ii) == st_fieldsfile) THEN
    ! output is to PP file
    NULLIFY(pp_file)
    pp_file => get_file_by_unit(-stlist(st_output_code,ii), &
                                handler="portio")
    pp_len2_look = pp_file % pp_meta % reserved_headers
    IF (pp_len2_look > pp_len2_lookUP) THEN
      pp_len2_lookUP = pp_len2_look
    END IF
  END IF

  !       Input levels list/range is always longer than output
  IF (stlist(st_input_bottom,ii) == st_special_code) THEN
    !          On special level
    num_levels=1
  ELSE IF (stlist(st_input_bottom,ii) <  0) THEN
    !          Using levels list, element 1 holds length.
    num_levels=stash_levels(1,-stlist(st_input_bottom,ii))
  ELSE
    !          Range
    num_levels=                                                  &
      stlist(st_input_top,ii)-stlist(st_input_bottom,ii)+1
  END IF
  IF (stlist(st_pseudo_in,ii) /= 0) THEN
    !          On pseudo levels
    num_pseudo_levels=stash_pseudo_levels(1,                     &
                     stlist(st_pseudo_in,ii))
  ELSE
    !          Not on pseudo levels
    num_pseudo_levels=1
  END IF
  max_stash_levs=MAX(max_stash_levs,num_levels*num_pseudo_levels)
END DO
! Round pp_len2_lookUP up to a multiple of 8
pp_len2_lookUP=((pp_len2_lookUP+7)/8)*8


! DEPENDS ON: fill_d1_array
CALL fill_d1_array(                                               &
                 icode,cmessage)
IF (icode /= 0) GO TO 9999

!----------------------------------------------------------------------
9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE initctl

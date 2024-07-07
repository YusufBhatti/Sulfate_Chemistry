! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine STASH --------------------------------------------------
!
! Purpose: Control routine for diagnostic processing step-by-step.
!          Called after each code section to process diagnostic fields
!          from workspace STASH_WORK to their final destination in D1
!          or PP file.  This routine loops over raw input fields and
!          calls a service routine STWORK to do the actual processing.
!
! Programming standard : UM Doc Paper no 3 vn8.3
!
! Project task : C4
!
! External documentation : UMDP no C4
!
! Interface and arguments --------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: stash

SUBROUTINE stash(sm_ident,im_ident,IS,stash_work,                 &
                 icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    printstatus,            &
    PrStatus_Diag
USE io_configuration_mod, ONLY: io_field_padding
USE UM_ParVars
USE decomp_params, ONLY: decomp_standard_atmos
USE Control_Max_Sizes
USE Decomp_DB
USE lookup_addresses
USE dump_headers_mod, ONLY: a_fixhd, a_inthd, a_realhd, a_levdepc, a_lookup
USE submodel_mod, ONLY: atmos_sm
USE d1_array_mod, ONLY: d1
USE stash_array_mod, ONLY:                                                     &
    len_stlist, num_pseudo_lists, stash_pseudo_levels, pp_len2_lookup,         &
    totitems, stash_series, num_stash_pseudo, stash_maxlen, nitems, stindex,   &
    nstash_series_records, stlist, num_level_lists, nstash_series_block,       &
    num_stash_levels, nsttims, stash_series_index, max_stash_levs,             &
    stash_levels, nsttabl, si, sttabl, nsects, time_series_rec_len, sf
USE stparam_mod, ONLY: st_proc_no_code, st_output_length, &
                       st_dump_level_output_length
USE cderived_mod, ONLY: elf

USE nlsizes_namelist_mod, ONLY:                                     &
    land_field, len1_lookup, len_fixhd, len_tot, model_levels,      &
    n_obj_d1_max, n_rows, river_row_length, river_rows, row_length, &
    rows, theta_field_size, v_field_size, a_len_inthd, a_len_realhd,&
    a_len1_levdepc, a_len2_levdepc, a_len2_lookup

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


INTEGER ::                                                        &
    sm_ident                                                      &
                   !IN      Submodel identifier
   ,im_ident                                                      &
                   !IN      Internal model identifier
   ,is                                                            &
                   !IN      Section number
   ,icode          !OUT     Return code

REAL ::                                                           &
    stash_work(*)  !IN     Area holding the data if not in D1
CHARACTER(LEN=errormessagelength) ::                              &
    cmessage       !OUT     ANY ERROR MESSAGE PASSED BACK

!---------------------------------------------------------------------
!
! Local variables.
!
LOGICAL ::                                                        &
   lcyclic                       ! TRUE if submodel is cyclic

INTEGER ::                                                        &
   ie,                                                            &
                                 ! Index over items in section
   ilend,                                                         &
                                 ! End point in STASHlist
   ilstart,                                                       &
                                 ! Start point in STASHlist
   il,                                                            &
                                 ! STASHlist index
   im,                                                            &
                                 ! Item number in section
   ippx,                                                          &
                                 ! Index to record in PP_XREF
   lenout                                                         &
                                 ! Maximum output length
   ,im_index                                                      &
                           ! Internal model index number
   ,num_rows1                                                     &
                           ! Number of rows in field type 1
   ,num_rows2                                                     &
                           ! Number of rows in field type 1
   ,row_len                                                       &
                           ! Row length
   ,field_len1                                                    &
                           ! Length of field type 1
   ,field_len2                                                    &
                           ! Length of field type 2
   ,num_levels                                                    &
                           ! Number of levels
   ,orig_decomp                                                   &
                           ! Decomposition on entry
,   global_LENOUT          ! Size of output field on disk

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STASH'

!--------------------------------------------------------------------
!0. Initialise.
!   Set LCYCLIC to indicate EW boundary condition.
!
! Find current decomposition, so we can return to this after STASH
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
orig_decomp=current_decomp_type

IF (sm_ident == atmos_sm) lcyclic = .NOT. elf

im_index = 1

icode = 0

!--------------------------------------------------------------------
!1. Loop over items within this section and call STWORK with
!   appropriate argument list depending on whether atmosphere/ocean
!
DO ie=1,nitems    ! max no of items in this section
  IF (stindex(2,ie,IS,im_index) >  0) THEN ! Entries for this
                                          ! SECTION/ITEM
    ilstart=stindex(1,ie,IS,im_index)
    ilend=ilstart+stindex(2,ie,IS,im_index)-1
    im=stlist(1,ilstart)             ! item number

    IF (sf(im,IS)) THEN       ! item/section reqd for this t/s
      IF (stlist(st_proc_no_code,stindex(1,ie,IS,im_index))       &
                                                 /= 0) THEN
        ! required by STASH so continue.
        ! It should not be possible to have any multiple entries for a given
        ! ITEM SECTION number and one of those not be required by STASH

        IF (printstatus >= prstatus_diag) THEN
          WRITE(umMessage,'(1X,A,I4,A,I4,A)') 'STASH: Item',im,          &
          ' Section',IS,' required by stash'
          CALL umPrint(umMessage,src='stash')
        END IF

        ! NB: max poss output length for the item/sect can be LONGER than
        !     the input length in the case of timeseries.
        lenout=0
        global_LENOUT=0
        DO il=ilstart,ilend
          lenout=MAX(lenout,stlist(st_output_length,il))
          global_LENOUT=                                        &
            MAX(global_LENOUT,                                  &
                stlist(st_dump_level_output_length,il))
        END DO
        ! Add on an extra io_field_padding on the end, ensuring array is
        ! big enough for data+extra space to round up to the next
        ! io_field_padding
        global_LENOUT=global_LENOUT+io_field_padding
        ! Make sure global_lenout at least as large as lenout
        ! so that time series data are not corrupted
        global_LENOUT = MAX(global_LENOUT,lenout)

        !    1 : Use PPXREF file to control packing
        !    2 : Do not pack prognostics, as 1 for diagnostics
        !    3 : Do not pack prognostics or diagnostics

        ! Make superarrays to pass into STWORK
        IF (im_ident  ==  atmos_sm) THEN
          ! Change to atmosphere decomposition
          IF (current_decomp_type /= decomp_standard_atmos) THEN
            CALL change_decomposition(decomp_standard_atmos,icode)
          END IF
          IF (icode  /=  0) THEN
            cmessage='STASH : Unsupported MPP submodel : atmos'
            GO TO 9999
          END IF
          num_rows1  = rows
          num_rows2  = n_rows
          row_len    = row_length
          field_len1 = theta_field_size
          field_len2 = v_field_size
          num_levels = model_levels

          ! DEPENDS ON: stwork
          CALL stwork(                                            &
           d1,len_tot,stash_work,stash_maxlen(IS,im_index),lenout,&
           global_LENOUT,                                         &
           IS,im,ilstart,ilend,                                   &
           stlist,len_stlist,totitems,si,nsects,nitems,           &
           stash_levels,num_stash_levels,num_level_lists,         &
           stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists, &
           max_stash_levs,sttabl,nsttims,nsttabl,                 &
           stash_series,nstash_series_records,time_series_rec_len,&
           stash_series_index,nstash_series_block,                &
           a_fixhd, a_inthd,                                      &
           a_realhd, len_fixhd, a_len_inthd, a_len_realhd,        &
           a_levdepc, a_len1_levdepc, a_len2_levdepc,             &
           a_lookup, a_lookup,                                    &
                           ! 2nd copy used as REAL in PP_HEAD
           len1_lookup, a_len2_lookup, pp_len2_lookUP,            &
           lcyclic,num_rows1,num_rows2,                           &
           row_len,field_len1,field_len2,num_levels,              &
           river_rows, river_row_length,                          &
           elf,                                                   &
           sm_ident,im_ident,                                     &
           icode,cmessage)
        END IF

      END IF

    END IF

  END IF
  ! Handle warning conditions on return from STWORK
  IF (icode <  0) THEN
    WRITE(umMessage,*) &
        'STASH    : Warning processing diagnostic section ',&
        IS,', item ',im,', code ',icode
    CALL umPrint(umMessage,src='stash')
    CALL umPrint(cmessage,src='stash')
    icode=0
  END IF
END DO

9999 CONTINUE
IF (current_decomp_type  /=  orig_decomp) THEN
  CALL change_decomposition(orig_decomp,                          &
                                      icode)
  IF (icode  /=  0) THEN
    cmessage='STASH : Unsupported MPP submodel'
  END IF
END IF
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stash

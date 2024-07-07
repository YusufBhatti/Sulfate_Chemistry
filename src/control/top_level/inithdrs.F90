! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE INITHDRS -----------------------------------------------
!
!    PURPOSE:   Initialises dump LOOKUP headers reserved for diagnostic
!               fields with size and other basic information to allow
!               dump IO routines to work correctly before STASH has
!               updated the addressed fields.
!
!    EXTERNAL DOCUMENTATION: UMDP NO. C4
!
!     -------------------------------------------------------------

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

SUBROUTINE inithdrs(                                              &
                  icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE lookup_addresses
USE cppxref_mod, ONLY: ppx_dump_packing, ppx_data_type
USE nlstgen_mod, ONLY: dump_packim
USE submodel_mod, ONLY: atmos_sm, submodel_partition_index
USE dump_headers_mod, ONLY: a_fixhd, a_lookup
USE stash_array_mod, ONLY: stash_pseudo_levels, totitems, stlist, stash_levels
USE stparam_mod, ONLY: st_output_code, st_dump, st_model_code,        &
    st_sect_no_code, st_item_code, st_proc_no_code, st_lookup_ptr,    &
    st_gridpoint_code, block_size, st_output_bottom,                  &
    st_time_series_code, st_time_series_mean, st_append_traj_code,    &
    st_output_top, st_pseudo_out, st_dump_output_length,              &
    st_dump_output_addr, vert_mean_base, global_mean_base
USE ppxlook_mod, ONLY: exppxi
USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lookup,                            &
    len_dumphist, len_fixhd, mpp_len1_lookup

USE missing_data_mod, ONLY: imdi, rmdi

USE errormessagelength_mod, ONLY: errormessagelength

USE Packing_Codes_Mod, ONLY: PC_Cray_Format, PC_No_Packing


IMPLICIT NONE

!   Arguments
!
INTEGER, INTENT(OUT) :: icode                 ! Error return code

CHARACTER(LEN=errormessagelength), INTENT(OUT)  ::  cmessage  
                                              !  Error return message

! Local variables

INTEGER :: i_eqv_rmdi
INTEGER :: disk_address                    ! Current rounded disk address
INTEGER :: number_of_data_words_on_disk    ! Number of data words on disk
INTEGER :: number_of_data_words_in_memory  ! Number of Data Words in memory
INTEGER :: i
INTEGER :: ii
INTEGER :: IS
INTEGER :: im
INTEGER :: iproc
INTEGER :: ilookup
INTEGER :: iheaders
INTEGER :: ilength
INTEGER :: imean
INTEGER :: j
INTEGER :: im_ident ! internal model identifier
INTEGER :: sm_ident ! submodel partition (dump) identifier

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INITHDRS'

! ----------------------------------------------------------------------
!   1. Set dump LOOKUP headers with basic information needed by
!      READDUMP and WRITDUMP, by scanning STASHlist items for
!      diagnostics destined for dump addresses.  NB: timeseries
!      fields cannot be 32-bit packed in dumps as the extra data
!      will contain integers.
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
sm_ident   = imdi     ! make sure this is always initialised

DO ii=1,totitems
  IF (stlist(st_output_code,ii) == st_dump) THEN
    !       output is to addressed D1

    im_ident  =stlist(st_model_code,ii)          ! internal model
    sm_ident  =submodel_partition_index(im_ident)! submodel

    IS     =stlist(st_sect_no_code,ii)
    im     =stlist(st_item_code,ii)
    iproc  =stlist(st_proc_no_code,ii)
    ilookup=stlist(st_lookup_ptr,ii)

    IF (sm_ident /= atmos_sm) THEN
      icode=111
      cmessage='INITHDRS: Only Atmos diagnostic requests possible'
      GO TO 9999
    END IF

    !         Calculate the total number of headers required by this
    !         stashlist record.
    imean=(stlist(st_gridpoint_code,ii)/block_size)*block_size
    IF (stlist(st_output_bottom,ii) == 100) THEN
      !           single level
      iheaders=1
    ELSE IF (stlist(st_proc_no_code,ii)  ==                       &
              st_time_series_code .OR.                             &
      stlist(st_proc_no_code,ii) == st_time_series_mean) THEN
      !           time series
      iheaders=1
    ELSE IF (stlist(st_proc_no_code,ii)  ==                       &
             st_append_traj_code) THEN
      !           append trajectories
      iheaders=1
    ELSE IF (imean == vert_mean_base) THEN
      !           vertical mean
      iheaders=1
    ELSE IF (imean == global_mean_base) THEN
      !           total 3-D mean
      iheaders=1
    ELSE IF (stlist(st_output_bottom,ii) <  0) THEN
      !           level list, not vertical mean.
      iheaders=stash_levels(1, -stlist(st_output_bottom,ii) )
    ELSE
      !           level range, not vertical mean.
      iheaders=stlist(st_output_top,ii)-                          &
               stlist(st_output_bottom,ii)+1
    END IF
    IF (stlist(st_pseudo_out,ii) >  0) THEN !Output pseudo levs
      iheaders=iheaders*                                          &
      stash_pseudo_levels(1,stlist(st_pseudo_out,ii))
    END IF

    ilength=stlist(st_dump_output_length,ii) / iheaders
    !         Loop down the headers.
    DO i=0,iheaders-1
      IF (sm_ident == atmos_sm) THEN
        DO j=1,len1_lookup
          a_lookup(j,ilookup+i)=imdi
        END DO
        DO j=46,len1_lookup
          a_lookup(j,ilookup+i)=TRANSFER(rmdi,i_eqv_rmdi)
        END DO
        a_lookup(lbnrec   ,ilookup+i)=0
        a_lookup(item_code,ilookup+i)=IS*1000+im
        a_lookup(model_code,ilookup+i)=im_ident
        a_lookup(data_type,ilookup+i)=                            &
                          exppxi(im_ident,IS,im,ppx_data_type,    &
                                             icode,cmessage)
        a_lookup(lblrec,ilookup+i)=ilength
        a_lookup(naddr ,ilookup+i)=                               &
          stlist(st_dump_output_addr  ,ii)+( ilength * i )
        IF (iproc == st_time_series_code .OR.                     &
            iproc == st_time_series_mean .OR.                     &
        iproc == st_append_traj_code) THEN
          a_lookup(lbpack,ilookup+i)=1000 * PC_Cray_Format
        ELSE
          a_lookup(lbpack,ilookup+i)=1000 * PC_Cray_Format +                   &
                                      exppxi(im_ident,IS,im,ppx_dump_packing,  &
                                             icode,cmessage)
          IF (dump_packim(sm_ident) == 3 ) THEN
            ! Do not pack data ; Override PPXREF packing indicator
            a_lookup(lbpack,ilookup+i) =                          &
           (a_lookup(lbpack,ilookup+i)/10)*10 + PC_No_Packing
          END IF
        END IF
      END IF
    END DO ! I, Loop over headers for this STASHlist entry

  END IF
END DO  ! II LOOP OVER TOTITEMS

!--reset the disk addresses and lengths for well-formed I/O
! DEPENDS ON: set_dumpfile_address
CALL set_dumpfile_address(a_fixhd, len_fixhd,                   &
                          a_lookup, len1_lookup,                &
                          a_len2_lookup,                        &
                          number_of_data_words_in_memory,       &
                          number_of_data_words_on_disk,         &
                          disk_address)
9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE inithdrs

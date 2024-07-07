! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the STASH addresses for D1
! Subroutine Interface:
SUBROUTINE addres(nrecs,                                          &
                  ErrorStatus,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE decomp_params, ONLY: decomp_standard_atmos
USE Decomp_DB
USE ppxlook_mod, ONLY: ppxref_items, ppxref_sections,             &
                       stashmaster_record_exists, exppxi
USE cppxref_mod, ONLY:                                            &
    ppx_version_mask, ppx_grid_type, ppx_halo_type, ppx_lev_flag, &
    ppx_lb_code, ppx_lt_code, ppx_lv_code,                        &
    ppx_pf_code, ppx_pl_code, ppx_pt_code,                        &
    ppx_ptr_code, ppx_opt_code, ppx_space_code

USE submodel_mod, ONLY:                                                        &
    internal_model_list, n_internal_model_max, n_internal_model,               &
    submodel_for_im, submodel_partition_index, n_submodel_partition_max,       &
    atmos_im

USE um_stashcode_mod, ONLY: stashcode_prog_sec, stashcode_lbc_input_sec, &
    stashcode_lbc_output_sec, stashcode_tracer_sec, stashcode_ukca_sec,  &
    stashcode_tracer_lbc_sec, stashcode_ukca_lbc_sec,                    &
    stashcode_glomap_clim_sec
USE stparam_mod, ONLY: st_model_code, st_output_code, st_output_addr,&
    st_dump_output_addr, st_output_length, st_dump_output_length,    &
    st_output_bottom, st_series_ptr, st_gridpoint_code,st_output_top,&
    st_pseudo_out, st_lookup_ptr, st_freq_code, st_start_time_code,  &
    st_end_time_code, st_sect_no_code, st_item_code, st_output_type
USE stextend_mod, ONLY: indx_s, in_s, d1_paddr, n_obj_d1,  &
           max_d1_len, list_s, itim_s, levlst_s, lenplst,  &
               d1_type, d1_im, d1_extra_info, diag, seco
USE cstash_mod, ONLY: iopn, ipseudo, iflag, ptr_prog, halo_type,  &
                      ipfirst, iplast, igp, ispace, vmsk, itop,  &
                      ibot, ilev, stsh_hours
USE totimp_mod, ONLY: totimp
USE glomap_clim_option_mod,     ONLY: l_glomap_mode_clim
USE ukca_nmspec_mod, ONLY: ukca_set_nmspec
USE ukca_scavenging_mod, ONLY: ukca_set_conv_indices
USE file_manager, ONLY: get_file_by_unit, um_file_type

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlstcall_mod, ONLY: run_target_end
USE nlsizes_namelist_mod, ONLY: &
    a_prog_len, a_prog_lookup, n_obj_d1_max

USE errormessagelength_mod, ONLY: errormessagelength

USE version_mod, ONLY: ntimep
USE stash_model_mod, ONLY:                                                     &
    len_extra, nhead, len_primim, len_prim, len_dump, global_len_prim,         &
    len_dumpim, global_len_dump, global_len_dumpim, nheadsub, len_secd,        &
    len_secdim, len_work

IMPLICIT NONE

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    Language: Fortran 95.
!    Written to UMDP3 programming standards.
!
! Global variables:

! Subroutine arguments:
!   Scalar arguments with intent(in):
INTEGER,INTENT(IN) :: nrecs

!   Scalar arguments with intent(out):
CHARACTER(LEN=errormessagelength),INTENT(OUT) :: cmessage

! ErrorStatus:
INTEGER,INTENT(OUT) :: ErrorStatus

! Local scalars:
INTEGER :: Im_ident  !Internal model identifier (absolute)
INTEGER :: Im_index  !Internal model index (expt. dependent)
INTEGER :: Sm_ident  !Submodel identifier (absolute)
INTEGER :: isec
INTEGER :: isec_loop
INTEGER :: iitm
INTEGER :: rlevs
INTEGER :: raddress
INTEGER :: PIrow
INTEGER :: i,j
INTEGER :: ifirst
INTEGER :: ifreq
INTEGER :: ihours
INTEGER :: ilast
INTEGER :: irec
INTEGER :: ih,il,ip,it
INTEGER :: len_work_s(n_submodel_partition_max)
INTEGER :: icode ! return from CHANGE_DECOMPOSITION

TYPE(um_file_type), POINTER :: pp_file

! Local arrays:
!    Submodel definitions array: stores list of Im_index's
!     for each submodel partition
INTEGER :: SM_def(n_submodel_partition_max,n_internal_model_max)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ADDRES'

!- End of Header ----------------------------------------------------


! 1.  Set STASHIN addresses and input lengths for primary fields

!   The address loop for primary fields is performed for each
!   internal model in turn. Hence, each internal model's primary
!   data occupies a contiguous block in D1. The order of these blocks
!   is the same as the order of the internal models given in the
!   array INTERNAL_MODEL_LIST.
!   User-defined prognostics are included in this primary addressing
!   routine, since they are incorporated into the ppxref lookup
!   arrays PPXI, PPXC in routine GETPPX.

!   Initialisation
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
n_obj_d1_max=0
DO i = 1,n_submodel_partition_max
  n_obj_d1(i)=0
END DO

DO i = 1,n_submodel_partition_max
  DO j = 1,n_internal_model_max
    SM_def(i,j) = 0
  END DO
END DO

!   Obtain submodel definitions and store in SMdef array
DO Im_index = 1,n_internal_model
  !   Submodel ident.
  Sm_ident =   submodel_for_im(Im_index)
  !   Internal model index
  SM_def(Sm_ident,Im_index) = Im_index
END DO

!   Primary address loop

!     Loop over submodel partitions
DO Sm_ident = 1,n_submodel_partition_max

  !       Initialise len_extra
  len_extra(Sm_ident)=0

  !     Initialise address for reconfiguration
  raddress = 1

  !      Loop over internal models for each SM partition
  DO Im_index = 1,n_internal_model

    !       Test whether current SM contains this IM
    IF (SM_def(Sm_ident,Im_index) >  0) THEN

      !        Obtain internal model identifier
      Im_ident   = internal_model_list(Im_index)

      ! Set the correct decomposition in PARVARS

      icode=0

      IF (Im_ident  ==  atmos_im) THEN
        IF (current_decomp_type  /=  decomp_standard_atmos) &
        CALL change_decomposition(decomp_standard_atmos,icode)

      ELSE  ! unsupported decomposition type
        WRITE(umMessage,'(A,A)') 'ADDRES1 : Error - Only atmosphere ',   &
             'submodel is currently supported for UM code.'
        CALL umPrint(umMessage,src='addres')
        ErrorStatus=-1
        cmessage='Unsupported submodel for UM code'
        GO TO 9999
      END IF

      IF (icode  /=  0) THEN
        WRITE(umMessage,'(A,A)') 'ADDRES1 : Error - Could not set ',     &
          'decomposition for selected submodel.'
        CALL umPrint(umMessage,src='addres')
        ErrorStatus=-2
        cmessage='Unsupported decomposition selected for UM code'
        GO TO 9999
      END IF

      !        Initialise primary data lengths
      IF (Im_ident == atmos_im) a_prog_len=0
      IF (Im_ident == atmos_im) a_prog_lookup=0
      PIrow  = 0
      DO ISEC_loop = 1,8   ! Currently there are eight sections
                           ! that contain "primary" type fields
        IF (ISEC_loop == 1) isec=stashcode_prog_sec        ! section 0 primary
        IF (ISEC_loop == 2) isec=stashcode_lbc_input_sec   ! LBC input
        IF (ISEC_loop == 3) isec=stashcode_lbc_output_sec  ! LBC output
        IF (ISEC_loop == 4) isec=stashcode_tracer_sec      ! Free tracers 
        IF (ISEC_loop == 5) isec=stashcode_ukca_sec        ! UKCA tracers
        IF (ISEC_loop == 6) isec=stashcode_tracer_lbc_sec  ! Free tracer lbcs
        IF (ISEC_loop == 7) isec=stashcode_ukca_lbc_sec    ! UKCA tracer lbcs
        IF (ISEC_loop == 8) isec=stashcode_glomap_clim_sec ! GLOMAP NWP
                                                           !   climatology
        !       Loop over primary section items
        DO iitm   = 1,ppxref_items
          !   Check whether there is a primary field corresponding
          !         to this item number
          IF (stashmaster_record_exists(isec,iitm)) THEN
            vmsk    = exppxi(Im_ident,isec,iitm,ppx_version_mask, &
                        ErrorStatus,cmessage)
            ispace  = exppxi(Im_ident,isec,iitm,ppx_space_code,   &
                        ErrorStatus,cmessage)
            igp     = exppxi(Im_ident,isec,iitm,ppx_grid_type,    &
                        ErrorStatus,cmessage)
            ilev    = exppxi(Im_ident,isec,iitm,ppx_lv_code,      &
                        ErrorStatus,cmessage)
            ibot    = exppxi(Im_ident,isec,iitm,ppx_lb_code,      &
                        ErrorStatus,cmessage)
            itop    = exppxi(Im_ident,isec,iitm,ppx_lt_code,      &
                        ErrorStatus,cmessage)
            DO i=1,6
              iopn(i)=exppxi(Im_ident,isec,iitm,ppx_opt_code+i-1, &
                        ErrorStatus,cmessage)
            END DO
            iflag   = exppxi(Im_ident,isec,iitm,ppx_lev_flag,     &
                        ErrorStatus,cmessage)
            ipseudo = exppxi(Im_ident,isec,iitm,ppx_pt_code,      &
                        ErrorStatus,cmessage)
            ipfirst = exppxi(Im_ident,isec,iitm,ppx_pf_code,      &
                        ErrorStatus,cmessage)
            iplast  = exppxi(Im_ident,isec,iitm,ppx_pl_code,      &
                        ErrorStatus,cmessage)
            halo_type = exppxi(Im_ident,isec,iitm,ppx_halo_type,  &
                        ErrorStatus,cmessage)
            IF ((ispace == 2) .OR. (ispace == 3) .OR. (ispace == 9)    &
           .OR. (ispace == 4) .OR. (ispace == 5) .OR. (ispace == 10)   &
           .OR. (ispace == 8)) THEN ! Primary variable
              ! DEPENDS ON: primary
              CALL primary(isec,iitm,Im_index,Im_ident,Sm_ident,  &
                   rlevs,raddress,PIrow,ErrorStatus,cmessage)
            END IF
          END IF  !  stashmaster_record_exists
        END DO    !  Loop over items
      END DO      !  ISEC_loop : Loop over sections
    END IF        !  test whether SM contains IM
  END DO          !  Loop over Im_index
END DO            !  Loop over SM partitions

! LOOKUP array lengths
a_prog_lookup = nhead(atmos_im)
! Primary data lengths
a_prog_len = len_primim(atmos_im)
CALL umPrint('',src='addres')
CALL umPrint(' ***********************************',src='addres')
WRITE(umMessage,'(A,I0)') ' ADDRES : A_PROG_LOOKUP = ',a_prog_lookup
CALL umPrint(umMessage,src='addres')
WRITE(umMessage,'(A,I0)') ' ADDRES : A_PROG_LEN    = ',a_prog_len
CALL umPrint(umMessage,src='addres')
WRITE(umMessage,'(A,2I6,2I9)')                                                &
  ' ADDRES : NHEAD, len_primim = ', nhead, len_primim
CALL umPrint(umMessage,src='addres')

! 2. Loop through stash list to set output addresses and
!                 header positions for diagnostics
DO irec=1,nrecs

  ! Read internal model number from stash list. Stash list has already
  ! been ordered by internal model, section, item. Thus, all the atmos
  ! diagnostic addressing will be done first, followed by the slab
  ! addressing in the case of a slab model.
  Im_ident = list_s(st_model_code,irec)
  ! Obtain submodel partition id.
  Sm_ident = submodel_partition_index(Im_ident)

  ! Set output address relative to D1
  IF (list_s(st_output_code,irec) == 1) THEN

    ! Diagnostic output to dump rather than direct output pp file
    !   Add the output length for this diag to len_dump; total length of
    !   dump so far = len_prim + len_dump; hence obtain the start address for
    !   the output from the next diagnostic to be stored in dump.

    list_s(st_output_addr,irec)                                   &
             = len_prim(Sm_ident)+len_dump(Sm_ident)+1
    ! Information for preliminary D1 addressing array
    n_obj_d1(Sm_ident)     =n_obj_d1(Sm_ident)+1
    IF (n_obj_d1(Sm_ident) <= max_d1_len) THEN
      d1_paddr(d1_type,n_obj_d1(Sm_ident),Sm_ident)=diag
      d1_paddr(d1_im,n_obj_d1(Sm_ident),Sm_ident)=Im_ident
      d1_paddr(d1_extra_info,n_obj_d1(Sm_ident),Sm_ident)=irec
    END IF
    list_s(st_dump_output_addr,irec)=                             &
             global_len_prim(Sm_ident)+global_len_dump(Sm_ident)+1
    len_dump(Sm_ident)                                               &
             = len_dump(Sm_ident)+list_s(st_output_length,irec)
    len_dumpim(Im_ident)                                             &
             = len_dumpim(Im_ident)+list_s(st_output_length,irec)
    global_len_dump(Sm_ident)=                                       &
      global_len_dump(Sm_ident)+list_s(st_dump_output_length,irec)
    global_len_dumpim(Sm_ident)=                                     &
      global_len_dumpim(Im_ident)+list_s(st_dump_output_length,irec)

    IF (list_s(st_output_bottom,irec) == 100) THEN
      ! Special levels
      rlevs=1
    ELSE IF (list_s(st_series_ptr,irec) /= 0) THEN
      ! Time series domain
      rlevs=1
    ELSE IF (list_s(st_gridpoint_code,irec) >= 10                  &
       .AND. list_s(st_gridpoint_code,irec) <  20) THEN
      ! Vertical ave.
      rlevs=1
    ELSE IF (list_s(st_output_bottom,irec) <  0) THEN
      ! Levels list
      rlevs=levlst_s(1,-list_s(st_output_bottom,irec))
    ELSE
      ! Range of model levels
      rlevs=list_s(st_output_top   ,irec)                         &
           -list_s(st_output_bottom,irec)+1
    END IF

    IF (list_s(st_pseudo_out,irec) >  0) THEN
      ! Pseudo levels
      rlevs=rlevs*lenplst(list_s(st_pseudo_out,irec))
    END IF

    ! Set position of pp lookup header in the dump
    list_s(st_lookup_ptr,irec)=NHeadSub(Sm_ident)+1

    ! Increment NHEAD (there is one pp header for each level at
    !  which a diagnostic is output
    nhead   (Im_ident)=nhead   (Im_ident)+rlevs
    NHeadSub(Sm_ident)=NHeadSub(Sm_ident)+rlevs

  ELSE IF (list_s(st_output_code,irec) == 2) THEN

    ! Secondary data in D1.
    ! Compute and store secondary data lengths. Start address for
    ! secondary data is determined below, after total dump
    ! diagnostic length has been found.

    list_s(st_output_addr,irec)=len_secd(Sm_ident)+1
    len_secd(Sm_ident)                                               &
   =len_secd(Sm_ident)+list_s(st_output_length,irec)
    len_secdim(Im_ident)                                             &
   =len_secdim(Im_ident)+list_s(st_output_length,irec)
    ! Set pointer for pp header
    list_s(st_lookup_ptr,irec)=-1

  ELSE IF (list_s(st_output_code,irec) <  0 .AND.                    &
           list_s(st_output_type,irec) == 1) THEN

    ! Diagnostic output to PP file

    ! Compute no. of pp headers for this diagnostic
    !   = output levels * pseudo output levels * output times

    !   No. of levels
    IF (list_s(st_output_bottom,irec) == 100) THEN
      ! Special levels
      il=1
    ELSE IF (list_s(st_series_ptr,irec) /= 0) THEN
      ! Time series dom
      il=1
    ELSE IF (list_s(st_gridpoint_code,irec) >= 10                  &
      .AND. list_s(st_gridpoint_code,irec) <  20) THEN
      ! Vertical average
      il=1
    ELSE IF (list_s(st_output_bottom,irec) <  0) THEN
      ! Levels list
      il=levlst_s(1,-list_s(st_output_bottom,irec))
    ELSE
      ! Range of mod levs
      il=list_s(st_output_top,irec)                               &
       -list_s(st_output_bottom,irec)+1
    END IF

    !   No. of pseudo levels
    IF (list_s(st_pseudo_out,irec) >  0) THEN
      ip=lenplst(list_s(st_pseudo_out,irec))
    ELSE
      ip=1
    END IF

    !   No. of output times
    IF (list_s(st_freq_code,irec) >  0) THEN
      ifirst=list_s(st_start_time_code,irec)
      ifreq =list_s(st_freq_code      ,irec)
      IF (list_s(st_end_time_code,irec) == -1) THEN
        ! Output to continues to end of run
        ihours=1+8760*run_target_end(1)                           &
                + 744*run_target_end(2)                           &
                +  24*run_target_end(3)                           &
                +     run_target_end(4)
        ilast=totimp(ihours,stsh_hours,Im_ident)
        IF (ilast  ==  -999) THEN
          errorStatus = 1
          cmessage = 'TOTIMP:UNEXPECTED TIME UNIT '//             &
              'or IRREGULAR DUMPS FOR DUMP FREQUENCY'
          GO TO 9999
        END IF
      ELSE
        ! Last output time before end of run
        ilast=list_s(st_end_time_code,irec)
      END IF

      it= 1 + (ilast-ifirst)/ifreq
      IF (it <  0) THEN
        it=0
        CALL umPrint(' Output time error detected in routine ADDRESS:', &
            src='addres')
        CALL umPrint(' Output time starts after specified end of run', &
            src='addres')
        WRITE(umMessage,'(A,4I6)')                                       &
      ' STASH record no.,MODEL,SECTION,ITEM as follows: ',               &
                        irec, list_s(st_model_code,irec),                &
                              list_s(st_sect_no_code,irec),              &
                              list_s(st_item_code,irec)
        CALL umPrint(umMessage,src='addres')
        WRITE(umMessage,'(A,I6)') 'OUTPUT CODE: ',                       &
                              list_s(st_output_code,irec)
        CALL umPrint(umMessage,src='addres')
      END IF
    ELSE
      ! Times table in STASH_times array
      it=1
      DO i=1,ntimep
        IF (itim_s(i,-list_s(st_freq_code,irec)) == -1) THEN
          it=i-1
          GO TO 260
        END IF
      END DO
      260        CONTINUE
    END IF
    ! No. of output "headers" - (levels)*(pseudo-levels)*(output times)
    ih=il*ip*it
    ! Assign output unit no. (nn) to (st_output_addr)
    list_s(st_output_addr,irec)=-list_s(st_output_code,irec)
    ! Accumulate no. of output headers to corresponding file
    NULLIFY(pp_file)
    pp_file => get_file_by_unit(list_s(st_output_addr, irec), &
                                handler="portio")
    pp_file % pp_meta % reserved_headers_calc = &
        pp_file % pp_meta % reserved_headers_calc + ih

  ELSE IF (list_s(st_output_code,irec) <  0 .AND.                    &
           list_s(st_output_type,irec) == 2) THEN

    ! Diagnostic output to NetCDF file

    ! Assign output unit no. (nn) to (st_output_addr)
    list_s(st_output_addr,irec)=-list_s(st_output_code,irec)

  ELSE IF (list_s(st_output_code,irec) == 0) THEN
    ! Inactive record, not output
    list_s(st_output_addr,irec)=-list_s(st_output_code,irec)
  ELSE
    CALL umPrint( 'ERROR detected in routine ADDRESS ',src='addres')
    CALL umPrint( 'ILLEGAL OUTPUT CODE FOR STASH RECORD ',src='addres')
    WRITE(umMessage,'(A,4I6)')                                           &
  ' STASH record no.,MODEL,SECTION,ITEM as follows: ',                   &
                   irec, list_s(st_model_code,irec),                     &
                         list_s(st_sect_no_code,irec),                   &
                         list_s(st_item_code,irec)
    CALL umPrint(umMessage,src='addres')
  END IF

END DO      ! End of loop over records for D1 addressing


!     Correct the addressing of SPACE=9 items from being relative
!     to start of len_extra space to being relative to start of dump

!     Loop over submodel partitions
DO  Sm_ident = 1,n_submodel_partition_max

  !       Loop over internal models for each SM partition
  DO Im_index = 1,n_internal_model

    !         Test whether current SM contains this IM
    IF (SM_def(Sm_ident,Im_index) >  0) THEN

      !           Obtain internal model identifier
      Im_ident   = internal_model_list(Im_index)

      DO ISEC_loop=1,8
        IF (ISEC_loop == 1) isec=stashcode_prog_sec        ! section 0 primary
        IF (ISEC_loop == 2) isec=stashcode_lbc_input_sec   ! LBC input
        IF (ISEC_loop == 3) isec=stashcode_lbc_output_sec  ! LBC output
        IF (ISEC_loop == 4) isec=stashcode_tracer_sec      ! Free tracers 
        IF (ISEC_loop == 5) isec=stashcode_ukca_sec        ! UKCA tracers
        IF (ISEC_loop == 6) isec=stashcode_tracer_lbc_sec  ! Free tracer lbcs
        IF (ISEC_loop == 7) isec=stashcode_ukca_lbc_sec    ! UKCA tracer lbcs
        IF (ISEC_loop == 8) isec=stashcode_glomap_clim_sec ! GLOMAP NWP
                                                           !   climatology
        DO iitm   = 1,ppxref_items
          !             Check whether there is a primary field corresponding
          IF (stashmaster_record_exists(isec,iitm) ) THEN
            ispace  = exppxi(im_ident,isec,iitm,ppx_space_code,     &
                             errorstatus,cmessage)
            IF (in_s(1,im_ident,isec,iitm) /= 0                     &
                .AND. ispace == 9) THEN  ! item is active
              in_s(1,im_ident,isec,iitm)=in_s(1,im_ident,isec,iitm)+&
                                         len_prim(sm_ident)+len_dump(sm_ident)
            END IF
          END IF
        END DO
      END DO ! ISEC_loop
    END IF
  END DO
END DO


! Set secondary data addresses relative to start of D1
DO irec=1,nrecs
  Im_ident = list_s(st_model_code,irec)
  Sm_ident = submodel_partition_index(Im_ident)

  IF (list_s(st_output_code,irec) == 2) THEN
    list_s(st_output_addr,irec)  =list_s(st_output_addr,irec)     &
  + len_prim(Sm_ident)+len_dump(Sm_ident)+len_extra(Sm_ident)
    ! Information for preliminary D1 addressing array
    n_obj_d1(Sm_ident)     =n_obj_d1(Sm_ident)+1
    IF (n_obj_d1(Sm_ident) <= max_d1_len) THEN
      d1_paddr(d1_type,n_obj_d1(Sm_ident),Sm_ident)=seco
      d1_paddr(d1_im,n_obj_d1(Sm_ident),Sm_ident)=Im_ident
      d1_paddr(d1_extra_info,n_obj_d1(Sm_ident),Sm_ident)=irec
    END IF
  END IF
END DO

! 3.  Set input addresses and work lengths for non-primary
!            fields (i.e., ISPACE=0 or 7)
DO Im_ident=1,n_internal_model_max
  Sm_ident=  submodel_partition_index(Im_ident)
  DO isec    =0,ppxref_sections
    ! Re-initialise sectional work lengths
    DO i=1,n_submodel_partition_max
      len_work_s(i)=0
    END DO
    DO iitm  =1,ppxref_items
      IF (indx_s(2,Im_ident,isec,iitm) >  0) THEN
        ! Item in STASH list
        ! Obtain space code & section zero point-back code
        !   from ppxref lookup array
        ispace  = exppxi(Im_ident,isec,iitm,ppx_space_code   ,        &
                                              ErrorStatus,cmessage)
        ptr_prog= exppxi(Im_ident,isec,iitm,ppx_ptr_code     ,        &
                                              ErrorStatus,cmessage)
        ! Compute length of work space required
        IF (ispace == 0) THEN
          ! STASH_WORK address & length
          in_s(1,Im_ident,isec,iitm)=len_work_s(Sm_ident)+1
          len_work_s(Sm_ident)=len_work_s(Sm_ident)               &
                            +in_s(2,Im_ident,isec,iitm)
        ELSE IF (ispace == 7) THEN
          ! Point-back to primary space in section 0
          in_s(1,Im_ident,isec,iitm    )                          &
         =in_s(1,Im_ident,0   ,ptr_prog)
          in_s(2,Im_ident,isec,iitm    )                          &
         =in_s(2,Im_ident,0   ,ptr_prog)
        END IF
      END IF
    END DO   ! Items

    ! Find max sectional work length for each submodel partition
    DO i=1,n_submodel_partition_max
      len_work(i)=MAX(len_work(i),len_work_s(i))
    END DO

  END DO     ! Sections
  IF (Sm_ident /= 0) THEN
    !       Save the maximum value for dimensioning full D1 address array
    n_obj_d1_max=MAX(n_obj_d1_max,n_obj_d1(Sm_ident))
    WRITE(umMessage,'(I12,A,I6)') n_obj_d1(Sm_ident),                    &
                          ' D1 items in submodel ',Sm_ident
    CALL umPrint(umMessage,src='addres')
  END IF
END DO     ! Models
IF (n_obj_d1_max >  max_d1_len) THEN
  CALL umPrint('ADDRES1: No of items in D1 exceeds maximum allowed:')
  WRITE(umMessage,'(A,I12,A,I12)') 'Number allowed ',max_d1_len,         &
                           ' Number requested ',n_obj_d1_max
  CALL umPrint(umMessage,src='addres')
  CALL umPrint('Modify the module STEXTEND_MOD to increase',src='addres')
  CALL umPrint('MAX_D1_LEN parameter as required',src='addres')
  CALL umPrint('Such a change can be safely made',src='addres')
  cmessage='ADDRES1: No of D1 items exceeds max: See output'
  ErrorStatus=1
END IF

CALL umPrint( '********************************************'// &
    '***********************************',src='addres')
CALL umPrint('',src='addres')

! 4. Initialise UKCA nm_spec array and mode indexing arrays
CALL ukca_set_nmspec()
CALL ukca_set_conv_indices()

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE addres

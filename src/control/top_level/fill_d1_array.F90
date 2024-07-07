! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE FILL_D1_ARRAY------------------------------------------
!
!    PURPOSE: Fill D1 addressing array with useful information.
!
!             Code Owner: Please refer to the UM file CodeOwners.txt
!             This file belongs in section: Top Level

SUBROUTINE fill_d1_array(                                         &
                  icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr
USE UM_ParParams
USE cppxref_mod, ONLY: ppx_grid_type, ppx_halo_type

USE submodel_mod, ONLY: submodel_for_sm
USE d1_array_mod, ONLY: d1_list_len, d1_object_type, d1_imodl, d1_section, &
                        d1_item, d1_address, d1_length, d1_grid_type,      &
                        d1_no_levels, d1_stlist_no, d1_lookup_ptr,         &
                        d1_north_code, d1_south_code, d1_east_code,        &
                        d1_west_code, d1_gridpoint_code, d1_proc_no_code,  &
                        d1_halo_type, prognostic, diagnostic, secondary,   &
                        other, d1_addr, no_obj_d1
USE stash_array_mod, ONLY: stlist
USE ppxlook_mod, ONLY: exppxi
USE stparam_mod, ONLY: st_output_addr, st_sect_no_code, st_item_code,&
    st_output_length, st_position_in_d1, st_north_code,st_south_code,&
    st_east_code, st_west_code, s_grid, s_proc, st_output_bottom,    &
    st_series_ptr, st_gridpoint_code, st_output_top, st_pseudo_out
USE stextend_mod, ONLY: in_s, d1_paddr, n_obj_d1,             &
            d1_type, d1_im, d1_extra_info, d1_levs, d1_sect,  &
            prog, diag, seco, extra_d1, levlst_s, lenplst

USE nlsizes_namelist_mod, ONLY: &
    len_tot, n_obj_d1_max

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER ::                                                        &
  ii,                                                             &
           ! Addresses preliminary array
  sm,                                                             &
           ! Addresses final array=1 for 1st submod =2 for 2nd
           ! submodel etc
  TYPE,                                                           &
           ! Code for prognostic, diagnostic, secondary or other
  iobj,                                                           &
           ! Addresses final array
  isec,                                                           &
           ! Section number
  itm,                                                            &
           ! Item number
  levs,                                                           &
           ! No of levels
  inf,                                                            &
        ! Diagnostic STASHlist number or prognosic item number
  Im_ident,                                                       &
  Sm_ident,                                                       &
  lookup_ptr,                                                     &
              ! Pointer to lookup table
  ext_addr,                                                       &
            ! Temporary pointer
  icode                   ! OUT: Error return code
!
CHARACTER(LEN=errormessagelength) ::                              &
    cmessage               ! OUT: Error return message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILL_D1_ARRAY'

! Initialise array
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
Sm_ident=1
DO ii=1,n_obj_d1_max
  DO inf=1,d1_list_len
    d1_addr(inf,ii,Sm_ident)=-1
    no_obj_d1(Sm_ident)=0
  END DO
END DO

IF (PrintStatus >= PrStatus_Oper) THEN
  ! Set up addressing of D1
  CALL umPrint('',src='fill_d1_array')
  CALL umPrint( '********************************************'// &
                   '***********************************', &
                   src='fill_d1_array')
  CALL umPrint('Addressing of D1 array',src='fill_d1_array')
  CALL umPrint('Key to Type:',src='fill_d1_array')
  CALL umPrint('Type=0: Prognostic',src='fill_d1_array')
  CALL umPrint('Type=1: Diagnostics in dump',src='fill_d1_array')
  CALL umPrint('Type=2: Secondary diagnostics',src='fill_d1_array')
  CALL umPrint('Type=3: Others (eg P_EXNER)',src='fill_d1_array')
END IF  ! PrintStatus test
sm=0
! Remove superfluous DO Sm_ident=1,N_SUBMODEL_PARTITION_MAX
Sm_ident=1   ! N_SUBMODEL_PARTITION_MAX is PARAMETER in submodel_mod
iobj=0
sm=submodel_for_sm(Sm_ident)
IF (sm /= 0) THEN
  IF (no_obj_d1(sm) == 0) THEN
    no_obj_d1(sm)=n_obj_d1(Sm_ident)
    WRITE(umMessage,'(A,I4)')'Submodel id ',Sm_ident
    CALL umPrint(umMessage,src='fill_d1_array')
    WRITE(umMessage,'(A,I4)')'Submodel Number ',sm
    CALL umPrint(umMessage,src='fill_d1_array')
    WRITE(umMessage,'(A,I6)') &
        'No of objects in this submodel: ',no_obj_d1(sm)
    CALL umPrint(umMessage,src='fill_d1_array')
    ! Address if submodel not empty and not already addressed
    DO ii=1,no_obj_d1(sm)
      !           Preliminary array held in D1_PADDR - full array in D1_ADDR
      !           Index II in D1_PADDR goes into index IOBJ of D1_ADDR
      !           First add prognostics followed by diagnostics...
      Im_ident=d1_paddr(d1_im,ii,Sm_ident)
      inf=d1_paddr(d1_extra_info,ii,Sm_ident)
      isec=d1_paddr(d1_sect,ii,Sm_ident)
      TYPE=d1_paddr(d1_type,ii,Sm_ident)
      IF (TYPE == prog) THEN
        iobj=iobj+1
        d1_addr(d1_stlist_no,iobj,sm)=inf
        d1_addr(d1_section,iobj,sm)=isec
        d1_addr(d1_no_levels,iobj,sm)=d1_paddr(d1_levs,ii,Sm_ident)
        d1_addr(d1_object_type,iobj,sm)=prognostic
        d1_addr(d1_imodl,iobj,sm)  = Im_ident
        d1_addr(d1_address,iobj,sm)= in_s(1,Im_ident,isec,inf)
      ELSE IF (TYPE == diag) THEN
        iobj=iobj+1
        d1_addr(d1_stlist_no,iobj,sm)=inf
        d1_addr(d1_object_type,iobj,sm)=diagnostic
        d1_addr(d1_imodl,iobj,sm)  = Im_ident
        d1_addr(d1_address,iobj,sm)= stlist(st_output_addr,inf)
      END IF
    END DO

    !         Extra data between primary and secondary diagnostics
    DO ii=1,no_obj_d1(sm)
      isec=d1_paddr(d1_sect,ii,Sm_ident)
      TYPE=d1_paddr(d1_type,ii,Sm_ident)
      IF (TYPE == extra_d1) THEN
        Im_ident=d1_paddr(d1_im,ii,Sm_ident)
        inf=d1_paddr(d1_extra_info,ii,Sm_ident)
        iobj=iobj+1
        d1_addr(d1_stlist_no,iobj,sm)=inf
        d1_addr(d1_section,iobj,sm)=isec
        d1_addr(d1_no_levels,iobj,sm)=d1_paddr(d1_levs,ii,Sm_ident)
        d1_addr(d1_object_type,iobj,sm)=other
        d1_addr(d1_imodl,iobj,sm)  = Im_ident
        !               NOT OCEAN: Address was calculated in ADDRES
        d1_addr(d1_address,iobj,sm)=in_s(1,Im_ident,isec,inf)
      END IF
    END DO
    !         Finally add secondary diagnostics
    DO ii=1,no_obj_d1(sm)
      !           Preliminary array held in D1_PADDR - full array in D1_ADDR
      Im_ident=d1_paddr(d1_im,ii,Sm_ident)
      inf=d1_paddr(d1_extra_info,ii,Sm_ident)
      TYPE=d1_paddr(d1_type,ii,Sm_ident)
      IF (TYPE == seco) THEN
        iobj=iobj+1
        d1_addr(d1_stlist_no,iobj,sm)=inf
        d1_addr(d1_object_type,iobj,sm)=secondary
        d1_addr(d1_imodl,iobj,sm)  = Im_ident
        d1_addr(d1_address,iobj,sm)= stlist(st_output_addr,inf)
      END IF
    END DO

    lookup_ptr=0
    DO ii=1,no_obj_d1(sm)
      TYPE= d1_addr(d1_object_type,ii,sm)
      isec= d1_addr(d1_section,ii,sm)
      inf = d1_addr(d1_stlist_no,ii,sm)
      Im_ident = d1_addr(d1_imodl,ii,sm)
      IF ((TYPE == prognostic) .OR. (TYPE == other)) THEN
        ! Prognostics don't have STASHlist numbers
        d1_addr(d1_stlist_no,ii,sm)= -1
        d1_addr(d1_item,ii,sm)   = inf
        d1_addr(d1_length,ii,sm) = in_s(2,Im_ident,isec,inf)
        isec = d1_addr(d1_section,ii,sm)
        itm  = inf
        !-------------------------------------------------------------------
        ! Prognostic items:
        ! Additional items can be added to the array here. Its code (eg
        ! d1_item, d1_levels) should be added to d1_array_mod and
        ! set as a parameter. The D1_LIST_LEN parameter should be changed
        ! as required
        !-------------------------------------------------------------------
      ELSE
        d1_addr(d1_section,ii,sm)= stlist(st_sect_no_code,inf)
        d1_addr(d1_item,ii,sm)   = stlist(st_item_code,inf)
        d1_addr(d1_length,ii,sm) = stlist(st_output_length,inf)
        isec=d1_addr(d1_section,ii,sm)
        itm=d1_addr(d1_item,ii,sm)
        ! STASH list pointer to D1 address information
        stlist(st_position_in_d1,inf) = ii
        !-------------------------------------------------------------------
        ! Diagnostic items
        ! Add items as per prognostics
        !-------------------------------------------------------------------
        d1_addr(d1_north_code,ii,sm)    =stlist(st_north_code,inf)
        d1_addr(d1_south_code,ii,sm)    =stlist(st_south_code,inf)
        d1_addr(d1_east_code,ii,sm)     =stlist(st_east_code,inf)
        d1_addr(d1_west_code,ii,sm)     =stlist(st_west_code,inf)
        d1_addr(d1_gridpoint_code,ii,sm)=stlist(s_grid,inf)
        d1_addr(d1_proc_no_code,ii,sm)  =stlist(s_proc,inf)
        ! 1. Number of levels
        IF (stlist(st_output_bottom,inf) == 100) THEN
          ! Special levels
          levs=1
        ELSE IF (stlist(st_series_ptr,inf) /= 0) THEN
          ! Time series domain
          levs=1
        ELSE IF (stlist(st_gridpoint_code,inf) >= 10               &
            .AND. stlist(st_gridpoint_code,inf) <  20) THEN
          ! Vertical ave.
          levs=1
        ELSE IF (stlist(st_output_bottom,inf) <  0) THEN
          ! Levels list
          levs=levlst_s(1,-stlist(st_output_bottom,inf))
        ELSE
          ! Range of model levels
          levs=stlist(st_output_top   ,inf)                       &
            -stlist(st_output_bottom,inf)+1
        END IF

        IF (stlist(st_pseudo_out,inf) >  0) THEN
          ! Pseudo levels
          levs=levs*lenplst(stlist(st_pseudo_out,inf))
        END IF
        d1_addr(d1_no_levels,ii,sm) = levs
      END IF
      !-------------------------------------------------------------------
      ! Items whose settings are common to progs and diags (eg from PPXREF)
      ! Add items as per prognostics
      ! ISEC and ITM set above
      !-------------------------------------------------------------------
      d1_addr(d1_grid_type,ii,sm) =                               &
        exppxi(Im_ident,isec,itm,ppx_grid_type,                   &
              icode, cmessage)
      d1_addr(d1_halo_type,ii,sm) =                               &
        exppxi(Im_ident,isec,itm,ppx_halo_type,                   &
              icode, cmessage)
      lookup_ptr=lookup_ptr+d1_addr(d1_no_levels,ii,sm)
      d1_addr(d1_lookup_ptr,ii,sm)=lookup_ptr
    END DO
    IF (PrintStatus >= PrStatus_Normal) THEN
      WRITE(umMessage,'(A1,5A6,A12,A11,A8,2A10)' ) '*', '#', 'Type', 'Modl',   &
        'Sect', 'Item', 'Address', 'Length', 'Levels', 'Gridtype', 'Halotype'
      CALL umPrint(umMessage,src='fill_d1_array')

      DO ii=1,no_obj_d1(sm)
        WRITE(umMessage,'(A1,5I6,I12,I11,I8,2I10)')                            &
          '*', ii, d1_addr(d1_object_type,ii,sm), d1_addr(d1_imodl,iobj,sm),   &
          d1_addr(d1_section,ii,sm), d1_addr(d1_item,ii,sm),                   &
          d1_addr(d1_address,ii,sm), d1_addr(d1_length,ii,sm),                 &
          d1_addr(d1_no_levels,ii,sm), d1_addr(d1_grid_type,ii,sm),            &
          d1_addr(d1_halo_type,ii,sm)
        CALL umPrint(umMessage,src='fill_d1_array')
      END DO
    END IF  ! PrintStatus test
  END IF ! IF (NO_OBJ_D1(SM) == 0) THEN
END IF


CALL umPrint( '********************************************'// &
    '***********************************',src='fill_d1_array')
CALL umPrint('',src='fill_d1_array')
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fill_d1_array

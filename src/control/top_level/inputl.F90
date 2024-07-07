! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Interface:

SUBROUTINE inputl(nrecs,                                          &
                   nlevels,ErrorStatus,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParParams,  ONLY: local_data
USE decomp_params, ONLY: decomp_standard_atmos
USE mpp_conf_mod, ONLY: include_halos_ew, include_halos_ns
USE Decomp_DB

USE ppxlook_mod, ONLY: ppxref_sections, ppxref_items, exppxi
USE cppxref_mod, ONLY:                                            &
    ppx_halo_type, ppx_grid_type, ppx_lev_flag,                   &
    ppx_pf_code, ppx_pl_code, ppx_pt_code, ppx_lv_code
USE version_mod, ONLY: nelemp, nlevlstsp, npslistp
USE submodel_mod, ONLY: n_internal_model_max, atmos_sm,           &
                        submodel_partition_index
USE stparam_mod, ONLY: st_output_bottom, st_input_bottom, st_pseudo_in,&
                       st_output_top, st_pseudo_out, st_input_top, &
                       st_input_code, st_input_length, st_output_length
USE stextend_mod, ONLY: llistty, indx_s, in_s, list_s,  &
                        rlevlst_s, levlst_s, lenplst
USE cstash_mod,   ONLY: pslist_d, npslists, ipfirst, iplast,  &
                        halo_type, ipseudo, igp, ilev, iflag

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

USE levsrt_mod, ONLY: levsrt

IMPLICIT NONE

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 90
!    Written to UMDP3 programming standards version 8.3.
!
! Global variables:

! Subroutine arguments:
!   Scalar arguments with intent(in):
INTEGER :: nrecs

!   Scalar arguments with intent(out):
INTEGER :: nlevels   ! total no. of sets of stash levels

!   Scalar arguments with intent(out):

! ErrorStatus:
INTEGER :: ErrorStatus

! Local scalars:
LOGICAL :: model_lev
CHARACTER(LEN=errormessagelength) :: cmessage
LOGICAL :: ladd
LOGICAL :: ldupll
INTEGER :: i,il,ilin
INTEGER :: istart,iend
INTEGER :: modl
INTEGER :: isec
INTEGER :: iitm
INTEGER :: ip_in
INTEGER :: ix1,ix2
INTEGER :: iy1,iy2
INTEGER :: iz_in
INTEGER :: len_in
INTEGER :: len_primin
INTEGER :: ndupll
INTEGER :: nlevin
INTEGER :: leno
INTEGER :: ipf,ipl
! local versions of the global subdomain limits
INTEGER :: local_IX1,local_IX2,local_IY1,local_IY2
INTEGER ::                                                        &
        orig_decomp                                               &
                         ! MPP decomposition before start
       ,decomp_type                                               &
                         ! decomposition type
       ,sm_ident         ! submodel identifier

! Function and subroutine calls:
LOGICAL :: disct_lev

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INPUTL'

!- End of Header ----------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

orig_decomp = current_decomp_type

DO modl=1,n_internal_model_max
  !
  !    Ensure that domain decomposition is consistent with submodel
  !
  sm_ident = submodel_partition_index(modl)
  IF (sm_ident == atmos_sm) THEN
    decomp_type = decomp_standard_atmos
  ELSE                            ! No decomposition defined:
    decomp_type = orig_decomp    !  return to original
  END IF

  CALL change_decomposition(decomp_type,ErrorStatus)

  IF (ErrorStatus >  0) THEN
    cmessage='INPUTL: ERROR in changing MPP decomposition'
    WRITE(umMessage,*) cmessage
    CALL umPrint(umMessage,src=RoutineName)
    GO TO 9999
  END IF
  DO isec=0,ppxref_sections
    DO iitm=1,ppxref_items
      IF (indx_s(2,modl,isec,iitm) >= 1) THEN
        ! At least one stash rec
        igp     = exppxi(modl,isec,iitm,ppx_grid_type     ,             &
                                               ErrorStatus,cmessage)
        ilev    = exppxi(modl,isec,iitm,ppx_lv_code       ,             &
                                               ErrorStatus,cmessage)
        iflag   = exppxi(modl,isec,iitm,ppx_lev_flag      ,             &
                                               ErrorStatus,cmessage)
        ipseudo = exppxi(modl,isec,iitm,ppx_pt_code       ,             &
                                               ErrorStatus,cmessage)
        ipfirst = exppxi(modl,isec,iitm,ppx_pf_code       ,             &
                                               ErrorStatus,cmessage)
        iplast  = exppxi(modl,isec,iitm,ppx_pl_code       ,             &
                                               ErrorStatus,cmessage)
        halo_type  = exppxi(modl,isec,iitm,ppx_halo_type,               &
                                               ErrorStatus,cmessage)

        istart=       indx_s(1,modl,isec,iitm)   ! Pos of 1st rec
        iend  =istart+indx_s(2,modl,isec,iitm)-1 ! Pos of last rec
        ! Diagnostics with input on levels list (IFLAG=1),
        !  rather than on all possible levels
        IF ((iflag  == 1   ) .AND.                                      &
                                                  !Input on lev list
           (istart == iend) .AND.                                      &
                                                  !Only 1 stash rec
           (list_s(st_output_bottom,istart) <  0))                    &
                                                  !Output on lev list
            THEN
          ! Only one stash record for this m,s,i - output levels list is
          !  the same as the input levels list
          list_s(st_input_bottom ,istart)=                            &
          list_s(st_output_bottom,istart)
        ELSE IF (iflag == 1 .AND. ilev /= 5) THEN
          ! Input on levels list & more than one stash request -
          !  construct input levels list
          nlevels=nlevels+1
          IF (nlevels >  nlevlstsp) THEN
            WRITE(umMessage,*) 'ERROR IN ROUTINE INPUTL:'
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,*) 'TOO MANY STASH LEVELS LISTS REQUESTED ',      &
                       'ARRAYS WILL BE OVERWRITTEN'
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,*) 'REDUCE NUMBER OF LEVELS LISTS'
            CALL umPrint(umMessage,src=RoutineName)
            ErrorStatus=1
            GO TO 9999
          END IF
          ! Construct input levels list: this is the combined list of all
          !  the output levels for all the stash requests for this m,s,i
          nlevin=1
          ! Set levels list type - real or integer
          ! DEPENDS ON: disct_lev
          model_lev=disct_lev(ilev,ErrorStatus,cmessage)
          IF (.NOT. model_lev) THEN
            ! Non-model levels - real
            llistty(nlevels)='R'
            ! Model levels - integer
          ELSE
            llistty(nlevels)='I'
          END IF
          ! Loop over stash recs for this m,s,i
          DO i=istart,iend
            ! Pointer for input level list
            list_s(st_input_bottom,i)=-nlevels
            IF (list_s(st_output_bottom,i) <  0) THEN
              ! There is an output levels list:
              !  For each of the levels in the output levels lists for the stash
              !   record I, find out whether this level is already present in the
              !   input levels list NLEVELS constructed so far.
              !   If it is, set LADD=F. Otherwise, LADD=T.
              ! Loop over output levels and check each one
              DO il=2,levlst_s(1,-list_s(st_output_bottom,i))+1
                ladd=.TRUE.
                IF (nlevin >  1) THEN
                  DO ilin=2,nlevin
                    IF (list_s(st_output_top,i) /= 1) THEN
                      ! Non-model levels: real
                      IF (rlevlst_s(il,-list_s(st_output_bottom,i))    &
                       ==                                             &
                         rlevlst_s(ilin,nlevels)) ladd=.FALSE.
                    ELSE
                      ! Model levels: integer
                      IF ( levlst_s(il,-list_s(st_output_bottom,i))    &
                       ==                                             &
                          levlst_s(ilin,nlevels)) ladd=.FALSE.
                    END IF
                  END DO
                END IF

                ! If LADD=T, add level 'IL' from stash record 'I' output levels list
                !  to input levels list NLEVELS
                IF (ladd) THEN
                  nlevin=nlevin+1
                  IF (list_s(st_output_top,i) /= 1) THEN
                    rlevlst_s(nlevin,nlevels)=                        &
                    rlevlst_s(il,-list_s(st_output_bottom,i))
                  ELSE
                    levlst_s(nlevin,nlevels)=                         &
                    levlst_s(il,-list_s(st_output_bottom,i))
                  END IF
                END IF
              END DO     ! Loop over levels

            ELSE
              ! Contiguous range of model levels for output, rather than list
              !  Compare output levels range for stash record I with input levs
              !   range NLEVELS. Any of the output levels not already present
              !   in the input range is added to the input list.
              DO il=list_s(st_output_bottom,i),                       &
                    list_s(st_output_top   ,i)
                ladd=.TRUE.
                DO ilin=2,nlevin
                  IF (il == levlst_s(ilin,nlevels)) ladd=.FALSE.
                END DO
                IF (ladd) THEN
                  nlevin=nlevin+1
                  levlst_s(nlevin,nlevels)=il
                END IF
              END DO
            END IF   !  Levels list/range
          END DO     !  Loop over stash recs

          ! Record no. of levels in input list just constructed
          levlst_s(1,nlevels)=nlevin-1

          IF (nlevin-1 == 0) THEN
            WRITE(umMessage,*) 'ORDINARY LEVEL'
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,*) 'ISEC=',isec
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,*) 'IITM=',iitm
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,*) 'NLEVELS=',nlevels
            CALL umPrint(umMessage,src=RoutineName)
            DO i=istart,iend
              WRITE(umMessage,*) 'I=',i
              CALL umPrint(umMessage,src=RoutineName)
              WRITE(umMessage,*) 'LIST_S(st_output_bottom)=',                 &
                          list_s(st_output_bottom,i)
              CALL umPrint(umMessage,src=RoutineName)
              IF (list_s(st_output_bottom,i) <  0) THEN
                DO il=2,levlst_s(1,-list_s(st_output_bottom,i))+1
                  WRITE(umMessage,*) 'IL=',il
                  CALL umPrint(umMessage,src=RoutineName)
                  WRITE(umMessage,*)                                          &
                  'LEVLST',levlst_s(il,-list_s(st_output_bottom,i))
                  CALL umPrint(umMessage,src=RoutineName)
                END DO
              ELSE
                WRITE(umMessage,*)                                            &
                'LIST_S(st_output_top=',list_s(st_output_top,i)
                CALL umPrint(umMessage,src=RoutineName)
              END IF
            END DO
          END IF
          ! Sort levels list
          CALL levsrt(llistty(  nlevels), levlst_s(1,nlevels),        &
                     levlst_s(2,nlevels),rlevlst_s(2,nlevels))

          ! Determine whether this levels list is a duplicate of another list
          ! DEPENDS ON: duplevl
          CALL duplevl(nlevels,ldupll,ndupll)
          IF (ldupll) THEN
            ! Duplicate list at NDUPLL - reset pointer and reduce NLEVELS by 1
            nlevels=nlevels-1
            DO i=istart,iend
              list_s(st_input_bottom,i)=-ndupll
            END DO
          END IF
        END IF   !Levels lists

        ! Pseudo levels lists
        IF ((iflag  == 1   ) .AND.                                      &
          ((istart == iend) .OR. (ipseudo == 0)) ) THEN
          ! Either no pseudo levels or only one request:
          ! Input pseudo levels list equals output list
          list_s(st_pseudo_in,istart)=list_s(st_pseudo_out,istart)
        ELSE IF (iflag == 1) THEN
          ! Input pseudo levels list with more than one request
          npslists=npslists+1

          IF (npslists >  npslistp) THEN
            WRITE(umMessage,'(A,A)') RoutineName,                       &
            ': ERROR at npslists check for iflag=1'
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,'(A,I0)')                                   &
            'STASH pseudo levels lists request exceeds limit of ', npslistp
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,'(A,A)') 'Reduce pseudo levels lists ',     & 
            'or increase npslistp in version_mod'
            CALL umPrint(umMessage,src=RoutineName)
            ErrorStatus=1
            GO TO 9999
          END IF
          ! Construct input pseudo list: combined list of all output
          !  pseudo levels for all stash requests for this m,s,i
          nlevin=0
          DO i=istart,iend
            list_s(st_pseudo_in,i)=npslists
            DO il=1,lenplst(list_s(st_pseudo_out,i))
              ladd=.TRUE.
              IF (nlevin >  0) THEN
                DO ilin=1,nlevin
                  IF ( pslist_d(il,list_s(st_pseudo_out,i)) ==         &
                      pslist_d(ilin,npslists)) ladd=.FALSE.
                END DO
              END IF
              IF (ladd) THEN
                nlevin=nlevin+1
                pslist_d(nlevin,npslists)=                            &
                pslist_d(il,list_s(st_pseudo_out,i))
              END IF
            END DO
          END DO
          lenplst(npslists)=nlevin

          IF (nlevin == 0) THEN
            WRITE(umMessage,*) 'PSEUDO LEVEL'
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,*) 'ISEC=',isec
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,*) 'IITM=',iitm
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,*) 'NPSLISTS=',npslists
            CALL umPrint(umMessage,src=RoutineName)
            DO i=istart,iend
              WRITE(umMessage,*) 'I=',i
              CALL umPrint(umMessage,src=RoutineName)
              WRITE(umMessage,*) 'LENPLST=',lenplst(list_s(st_pseudo_out,i))
              CALL umPrint(umMessage,src=RoutineName)
              DO il=1,lenplst(list_s(st_pseudo_out,i))
                WRITE(umMessage,*) 'IL=',il
                CALL umPrint(umMessage,src=RoutineName)
                WRITE(umMessage,*)                                            &
                'PSLIST_D',pslist_d(il,list_s(st_pseudo_out,i))
                CALL umPrint(umMessage,src=RoutineName)
              END DO
            END DO
          END IF
          ! Sort input pseudo levels list.  The REAL argument is really just a dummy
          ! since we are processing integer levels here.
          CALL levsrt('I',lenplst(npslists),pslist_d(1,npslists),     &
                      REAL(pslist_d(1:lenplst(npslists),npslists)))
          ! Find out if duplicate
          ! DEPENDS ON: duppsll
          CALL duppsll(ldupll,ndupll)
          IF (ldupll) THEN
            ! Duplicate pseudo list at NDUPLL
            npslists=npslists-1
            DO i=istart,iend
              list_s(st_pseudo_in,i)=ndupll
            END DO
          END IF
        ELSE IF (iflag == 0 .AND. ipseudo /= 0) THEN
          ! Input pseudo levels list contains all possible pseudo levels for
          !  this diagnostic
          npslists=npslists+1

          IF (npslists >  npslistp) THEN
            WRITE(umMessage,'(A,A)') RoutineName, &
            ': ERROR at npslists check for iflag=0'
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,'(A,I0)')                                   &
            'STASH pseudo levels lists request exceeds limit of ', npslistp
            CALL umPrint(umMessage,src=RoutineName)
            WRITE(umMessage,'(A,A)') 'Reduce pseudo levels lists ',     & 
            'or increase npslistp in version_mod'
            CALL umPrint(umMessage,src=RoutineName)
            ErrorStatus=1
            GO TO 9999
          END IF

          DO i=istart,iend
            list_s(st_pseudo_in,i)=npslists
            ! Decode first & last pseudo level codes from stash master
            ! DEPENDS ON: pslevcod
            CALL pslevcod(ipfirst,ipf,'F',ErrorStatus,cmessage)
            ! DEPENDS ON: pslevcod
            CALL pslevcod(iplast ,ipl,'L',ErrorStatus,cmessage)
            ! Construct list
            DO nlevin = ipf,ipl
              pslist_d(nlevin,npslists)=nlevin
            END DO
          END DO
          lenplst(npslists)=ipl-ipf+1
        END IF   ! Pseudo levels

        ! Calculate horizontal factor for input length
        ! DEPENDS ON: lltorc
        CALL lltorc(igp,90,-90,0,360,iy1,iy2,ix1,ix2)

        ! Convert from global to local subdomain limits
        ! DEPENDS ON: global_to_local_subdomain
        CALL global_to_local_subdomain( include_halos_ew,include_halos_ns,&
                                        igp,halo_type,mype,             &
                                        iy1,ix2,iy2,ix1,                &
                                        local_IY1,local_IX2,            &
                                        local_IY2,local_IX1)
        ix1=local_IX1
        ix2=local_IX2
        iy1=local_IY1
        iy2=local_IY2
        ! All sub-model grids: atmos/ocean/wave now ordered S->N
        len_in=(ix2-ix1+1)*(iy2-iy1+1)

        ! Calculate vertical levels factor for input length
        IF (ilev /= 5) THEN
          ! More than one level
          IF (list_s(st_input_bottom,istart) <  0) THEN
            ! Level list
            iz_in=levlst_s(1,-list_s(st_input_bottom,istart))
          ELSE
            ! Range of model levs
            iz_in=list_s(st_input_top   ,istart)-                     &
                  list_s(st_input_bottom,istart)+1
          END IF
        ELSE
          ! Single level input
          iz_in=1
        END IF

        ! Calculate pseudo levels factor for input length
        IF (ipseudo /= 0) THEN
          ip_in=lenplst(list_s(st_pseudo_in,istart))
        ELSE
          ip_in=1
        END IF

        ! Calculate input length for this diag. and store in LIST_S
        ! Input_code <  0 means that a diag already processed into D1 is being
        !   reprocessed, so input length of child diag equals output length of
        !   parent.
        ! Otherwise, the input len is given by the product of the appropriate
        !   x-,y-,z-, and p-dimensions.
        DO i=istart,iend
          IF (list_s(st_input_code  ,i) >= 0) THEN
            list_s(st_input_length,i)=len_in*iz_in*ip_in
          ELSE
            list_s(st_input_length ,i)=                              &
            list_s(st_output_length,-list_s(st_input_code,i))
          END IF
          ! Store model no. in last element of LIST_S - for ADDRES
          list_s(nelemp+1,i)=modl
        END DO

        ! Recalculate input length for non-primary (length unchanged for
        ! most cases) and store in IN_S array.
        IF (isec /= 0) THEN
          IF ((igp /= 31) .AND. (igp /= 32)) THEN
            ! DEPENDS ON: addrln
            CALL addrln(igp,halo_type,len_primin,local_data)
            in_s(2,modl,isec,iitm)=len_primin*iz_in*ip_in
          END IF
        END IF

      END IF ! At least one stash record for m,s,i

    END DO   ! Items
  END DO   ! Sections
END DO   ! Models
!
CALL change_decomposition(orig_decomp,ErrorStatus)

IF (ErrorStatus >  0) THEN
  cmessage='INPUTL: ERROR in original MPP decomposition'
  WRITE(umMessage,*) cmessage
  CALL umPrint(umMessage,src=RoutineName)
  GO TO 9999
END IF

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE inputl

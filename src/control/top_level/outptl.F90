! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calc stash list output lens; reset boundary spec for full area output.

! Subroutine Interface:

SUBROUTINE outptl(                                                &
                  nrecs,ErrorStatus,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE decomp_params, ONLY: decomp_standard_atmos, decomp_unset
USE mpp_conf_mod, ONLY: include_halos_ew, include_halos_ns
USE Decomp_DB

USE cppxref_mod, ONLY:                                            &
    ppx_grid_type, ppx_halo_type, ppx_pt_code, ppx_lv_code
USE version_mod, ONLY:                                            &
    nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
    nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
    nlevp, npslevp, npslistp
USE cstash_mod, ONLY: ipseudo, halo_type, ilev, igp
USE submodel_mod, ONLY: atmos_im
USE stextend_mod, ONLY: nrecs_ts, list_s, levlst_s, lenplst
USE stparam_mod, ONLY: st_model_code, st_sect_no_code, st_item_code,&
    st_proc_no_code, st_output_length, st_input_code,st_period_code,&
    st_freq_code, st_series_ptr, st_south_code, st_east_code,       &
    st_north_code, st_west_code, st_gridpoint_code,st_output_bottom,&
    st_output_top, st_pseudo_out, st_dump_output_length,            &
    st_dump_level_output_length

USE umPrintMgr, ONLY:      &
    umPrint,               &
    umMessage
USE ppxlook_mod, ONLY: exppxi

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 90
!    Written to UMDP3 programming standards version 8.3.
!
! Subroutine arguments
!   Array arguments with intent(in):
INTEGER :: nrecs

!   Array arguments with intent(out):
CHARACTER(LEN=errormessagelength) :: cmessage

! ErrorStatus:
INTEGER :: ErrorStatus

! Local variables
INTEGER :: output_length
INTEGER :: ie
INTEGER :: IN
INTEGER :: ip_dim
INTEGER :: irec
INTEGER :: IS
INTEGER :: modl
INTEGER :: isec
INTEGER :: item
INTEGER :: it_dim
INTEGER :: iw
INTEGER :: ix_dim
INTEGER :: iy_dim
INTEGER :: iz_dim
INTEGER ::                                                        &
! local versions of the global subdomain boundaries
        local_north,local_east,local_south,local_west                   &
      , local_IN,local_IE,local_IS,local_IW                             &
! global versions of the X and Y horizontal dimensions, and
! total output size
      , global_IX_DIM,global_IY_DIM,global_output_length                &
! variables indicating the decomposition type at various stages
      , orig_decomp,decomp_type

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OUTPTL'

!- End of Header --------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
orig_decomp=current_decomp_type

! Loop over STASH records
DO irec=1,nrecs

  ! Obtain model, section, item for this record
  modl = list_s(st_model_code  ,irec)
  isec = list_s(st_sect_no_code,irec)
  item = list_s(st_item_code   ,irec)

  ! Set the correct decomposition type for this model
  IF (modl  ==  atmos_im) THEN
    decomp_type=decomp_standard_atmos
  ELSE
    !       Shouldn't get to this
    decomp_type=decomp_unset
    CALL umPrint( 'OUTPTL : Error',src='outptl')
    WRITE(umMessage,'(A,I4,A)') 'Unsupported Model ',modl,' for MPP code'
    CALL umPrint(umMessage,src='outptl')
    cmessage='Unsupported Model for MPP code'
    ErrorStatus=-1
    GO TO 9999
  END IF

  IF (current_decomp_type  /=  decomp_type) THEN
    CALL change_decomposition(decomp_type,ErrorStatus)
    IF (ErrorStatus  /=  0) THEN
      CALL umPrint( 'OUTPUTL : Error',src='outptl')
      WRITE(umMessage,'(A,A,I4)') 'Call to CHANGE_DECOMPOSITION failed ',  &
                 'with decomposition type ',decomp_type
      CALL umPrint(umMessage,src='outptl')
      cmessage='Unsupported decomposition for MPP code'
      GO TO 9999
    END IF
  END IF


  ! Extract level code, grid type code from ppx lookup array
  ilev    = exppxi(modl,isec,item,ppx_lv_code,                    &
                                ErrorStatus,cmessage)
  igp     = exppxi(modl,isec,item,ppx_grid_type,                  &
                                ErrorStatus,cmessage)
  ipseudo = exppxi(modl ,isec ,item,ppx_pt_code      ,            &
                                ErrorStatus,cmessage)
  halo_type = exppxi(modl,isec,item,ppx_halo_type,                &
                                ErrorStatus,cmessage)
  IF (list_s(st_proc_no_code,irec) == 0) THEN
    ! Dummy record - output length zero
    list_s(st_output_length,irec)=0
  ELSE IF (list_s(st_input_code,irec) <  0 .AND.                    &
         list_s(st_proc_no_code,irec) /= 8) THEN
    ! Child record - get output length from parent
    list_s(st_output_length,irec)=                                &
    list_s(st_output_length,-list_s(st_input_code,irec))

  ELSE
    ! Neither dummy nor child - calculate output length
    !   T dimension (equals 1 except for the time series case)
    IF ((list_s(st_proc_no_code,irec) == 1) .OR.                    &
       (list_s(st_proc_no_code,irec) == 2) .OR.                    &
       (list_s(st_proc_no_code,irec) == 3) .OR.                    &
       (list_s(st_proc_no_code,irec) == 5) .OR.                    &
       (list_s(st_proc_no_code,irec) == 6)) THEN     
      it_dim=1
    ELSE IF (list_s(st_proc_no_code,irec) == 4 .OR.                &
            list_s(st_proc_no_code,irec) == 8) THEN
      ! Time series case
      it_dim=                                                     &
      list_s(st_period_code,irec)/list_s(st_freq_code,irec)

    ELSE
      WRITE(umMessage,'(A,A,I4,A,I6)')'OUTPTL: ',                        &
      'ERROR UNEXPECTED PROCESSING CODE',                                &
      list_s(st_proc_no_code,irec),'   FOR RECORD ',irec
      CALL umPrint(umMessage,src='outptl')
      CALL umPrint(umMessage,src='outptl')
    END IF

    IF (list_s(st_series_ptr,irec) == 0) THEN
      ! Set up local versions of the boundaries of the subdomain

      ! DEPENDS ON: global_to_local_subdomain
      CALL global_to_local_subdomain(include_halos_ew,include_halos_ns,&
                               igp,halo_type,mype,                &
                               list_s(st_south_code,irec),        &
                               list_s(st_east_code,irec),         &
                               list_s(st_north_code,irec),        &
                               list_s(st_west_code,irec),         &
                               local_south,local_east,            &
                               local_north,local_west)

      ! Not a time series profile
      !   X dimension
      IF (list_s(st_gridpoint_code,irec) <  0) THEN
        WRITE(umMessage,'(A,A,I4,A,I6)')'OUTPTL: ',                      &
        'ERROR UNEXPECTED GRIDPOINT CODE',                               &
        list_s(st_gridpoint_code,irec),'   FOR RECORD ',irec
        CALL umPrint(umMessage,src='outptl')
        CALL umPrint(umMessage,src='outptl')
      ELSE IF (list_s(st_gridpoint_code,irec) <  20) THEN
        ix_dim=local_east-local_west+1
        global_IX_DIM=list_s(st_east_code,irec)-                  &
                      list_s(st_west_code,irec)+1
      ELSE IF (list_s(st_gridpoint_code,irec) <  30) THEN
        ix_dim=1
        global_IX_DIM=1
      ELSE IF (list_s(st_gridpoint_code,irec) <  40) THEN
        ix_dim=local_east-local_west+1
        global_IX_DIM=list_s(st_east_code,irec)-                  &
                      list_s(st_west_code,irec)+1
      ELSE IF (list_s(st_gridpoint_code,irec) <= 43) THEN
        ix_dim=1
        global_IX_DIM=1
      ELSE
        WRITE(umMessage,'(A,A,I4,A,I6)')'OUTPTL: ',                      &
        'ERROR UNEXPECTED GRIDPOINT CODE',                               &
        list_s(st_gridpoint_code,irec),'   FOR RECORD ',irec
        CALL umPrint(umMessage,src='outptl')
        CALL umPrint(umMessage,src='outptl')
      END IF  ! X dim

      IF (ix_dim <  1) THEN
        ! Area cut by global model
        ! DEPENDS ON: lltorc
        CALL lltorc(igp,90,-90,0,360,IN,IS,iw,ie)
        ! DEPENDS ON: global_to_local_subdomain
        CALL global_to_local_subdomain(                           &
          include_halos_ew,include_halos_ns, igp ,halo_type, mype,&
          IN,ie,IS,iw,                                            &
          local_IN,local_IE,local_IS,local_IW)

        ix_dim=ix_dim+local_IE-2*halosize(1,halo_type)
        ! Subtract two halos, because we don't want wrap around to include
        ! the halo at the end, and the beginning of field

      END IF
      IF (global_IX_DIM <  1) THEN
        ! DEPENDS ON: lltorc
        CALL lltorc(igp,90,-90,0,360,IN,IS,iw,ie)
        global_IX_DIM=global_IX_DIM+ie
      END IF

      !   Y dimension
      IF (list_s(st_gridpoint_code,irec) <  0) THEN
        WRITE(umMessage,'(A,A,I4,A,I6)')'OUTPTL: ',                      &
        'ERROR UNEXPECTED GRIDPOINT CODE',                               &
        list_s(st_gridpoint_code,irec),'   FOR RECORD ',irec
        CALL umPrint(umMessage,src='outptl')
        CALL umPrint(umMessage,src='outptl')
      ELSE IF (list_s(st_gridpoint_code,irec) <  30) THEN
        ! Atmos grid - first lat is southern most
        iy_dim=local_north-local_south+1
        global_IY_DIM=list_s(st_north_code,irec)-                 &
                      list_s(st_south_code,irec)+1
      ELSE IF (list_s(st_gridpoint_code,irec) <= 40) THEN
        iy_dim=1
        global_IY_DIM=1
      ELSE IF (list_s(st_gridpoint_code,irec) <= 43) THEN
        iy_dim=1
        global_IY_DIM=1
      ELSE
        WRITE(umMessage,'(A,A,I4,A,I6)')'OUTPTL: ',                      &
        'ERROR UNEXPECTED GRIDPOINT CODE',                               &
        list_s(st_gridpoint_code,irec),'   FOR RECORD ',irec
        CALL umPrint(umMessage,src='outptl')
        CALL umPrint(umMessage,src='outptl')
      END IF  ! Y dim

      !   Z dimension
      IF (list_s(st_gridpoint_code,irec) <  0) THEN
        WRITE(umMessage,'(A,A,I4,A,I6)')'OUTPTL: ',                      &
        'ERROR UNEXPECTED GRIDPOINT CODE',                               &
        list_s(st_gridpoint_code,irec),'   FOR RECORD ',irec
        CALL umPrint(umMessage,src='outptl')
        CALL umPrint(umMessage,src='outptl')
      ELSE IF (list_s(st_gridpoint_code,irec) <  10) THEN
        IF (ilev == 5) THEN
          iz_dim=1
        ELSE IF (list_s(st_output_bottom,irec) <  0) THEN
          iz_dim=levlst_s(1,-list_s(st_output_bottom,irec))
        ELSE
          iz_dim=list_s(st_output_top,irec)-                      &
                 list_s(st_output_bottom,irec)+1
        END IF
      ELSE IF (list_s(st_gridpoint_code,irec) <  20) THEN
        iz_dim=1
      ELSE IF (list_s(st_gridpoint_code,irec) <= 43) THEN
        IF (ilev == 5) THEN
          iz_dim=1
        ELSE IF (list_s(st_output_bottom,irec) <  0) THEN
          iz_dim=levlst_s(1,-list_s(st_output_bottom,irec))
        ELSE
          iz_dim=list_s(st_output_top,irec)-                      &
                 list_s(st_output_bottom,irec)+1
        END IF
      ELSE
        WRITE(umMessage,'(A,A,I4,A,I6)')'OUTPTL: ',                      &
        'ERROR UNEXPECTED GRIDPOINT CODE',                               &
        list_s(st_gridpoint_code,irec),'   FOR RECORD ',irec
        CALL umPrint(umMessage,src='outptl')
        CALL umPrint(umMessage,src='outptl')
      END IF  ! Z dim

      !   P dimension - pseudo levels
      IF (ipseudo >  0) THEN
        ip_dim=lenplst(list_s(st_pseudo_out,irec))
      ELSE
        ip_dim=1
      END IF

      ! Output length - total number of points
      output_length = it_dim*ix_dim*iy_dim*iz_dim*ip_dim

      global_output_length =                                      &
        it_dim*global_IX_DIM*global_IY_DIM*iz_dim*ip_dim
      list_s(st_output_length,irec) = output_length
      list_s(st_dump_output_length,irec) =                        &
        global_output_length
      list_s(st_dump_level_output_length,irec) =                  &
        global_IX_DIM*global_IY_DIM  ! size of horizontal field
    ELSE    ! Time series profile
      list_s(st_output_length,irec)=                              &
       nrecs_ts(list_s(st_series_ptr,irec))*it_dim+               &
      (nrecs_ts(list_s(st_series_ptr,irec))+1)*6

      list_s(st_dump_output_length,irec)=                         &
        list_s(st_output_length,irec)
    END IF
  END IF    ! Neither dummy nor child
END DO      ! Loop over STASH records
IF ((orig_decomp  /=  current_decomp_type) .AND.                  &
    (orig_decomp  /=  decomp_unset)) THEN
  CALL change_decomposition(orig_decomp,ErrorStatus)
END IF

9999  CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE outptl

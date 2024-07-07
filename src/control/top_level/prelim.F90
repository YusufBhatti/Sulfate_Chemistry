! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Construct preliminary STASH list of user requests

! Subroutine Interface:

SUBROUTINE prelim(nrecs,                                                      &
  ntimes,nlevels,errorstatus,cmessage)

USE rad_input_mod, ONLY: a_lw_radstep_diag, a_sw_radstep_diag,                &
                         a_lw_radstep_prog, a_sw_radstep_prog
USE asad_mod,      ONLY: interval

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE UM_ParParams
USE cv_param_mod,       ONLY: a_conv_step
USE river_inputs_mod, ONLY: river_step
USE jules_vegetation_mod, ONLY: phenol_period, triffid_period
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
USE submodel_mod, ONLY: atmos_im

USE stextend_mod, ONLY: llistty, npos_ts, list_s, itim_s,                     &
                        rlevlst_s, levlst_s, lenplst
USE stparam_mod, ONLY: st_model_code, st_sect_no_code, st_item_code,          &
    st_lookup_ptr, st_output_addr, st_input_code, st_gridpoint_code,          &
    st_weight_code, st_south_code, st_north_code, st_west_code,               &
    st_east_code, st_input_bottom, st_input_top, st_output_bottom,            &
    st_output_top, st_pseudo_in, st_pseudo_out, st_series_ptr,                &
    st_output_code, st_macrotag, st_proc_no_code, st_offset_code,             &
    st_period_code, st_freq_code, st_start_time_code,                         &
    st_end_time_code, st_output_ts0_code, st_diag_address,                    &
    st_domain_code, st_time_unit2_code, st_time_unit3_code,                   &
    st_output_type, st_domain_global, st_domain_n_hemisphere,                 &
    st_domain_s_hemisphere, st_domain_30_to_90_N, st_domain_30_to_90_S,       &
    st_domain_0_to_30_N, st_domain_0_to_30_S, st_domain_30_S_to_30_N,         &
    st_domain_whole_degrees, st_domain_gridpoints, st_levels_model_rho,       &
    st_levels_model_theta, st_levels_deep_soil, st_levels_single

USE cppxref_mod, ONLY:                                                        &
    ppx_version_mask, ppx_grid_type, ppx_lev_flag,                            &
    ppx_timavail_code, ppx_space_code, ppx_opt_code,                          &
    ppx_lbvc_code, ppx_ptr_code, ppx_item_number,                             &
    ppx_lb_code, ppx_lt_code, ppx_lv_code,                                    &
    ppx_pf_code, ppx_pl_code, ppx_pt_code

USE version_mod, ONLY: nelemp, nrecdp, nlevlstsp, nproftp

USE cstash_mod, ONLY: iest_d, iwst_d, isth_d, inth_d, ts_d,                   &
                      pllen_d, plpos_d, plt_d, imsk_d, ptr_prog,              &
                      iplast, imn_d, iopa_d, iwt_d, iend_t,                   &
                      istr_t, isdt_t, iedt_t, isam_t, unt2_t,                 &
                      ioff_t, iunt_u, locn_u, pslist_d, ityp_t,               &
                      ifre_t, unt1_t, intv_t, levlst_d, levt_d,               &
                      iopl_d, rlevlst_d, isec_b, modl_b, ndiag,               &
                      itim_t, modl_t, iopt_t, iser_t, levb_d,                 &
                      ndprof, unt3_t, item_b, itop, ibot, ilev,               &
                      iflag, ipfirst, ipseudo, iopn, itim_b,                  &
                      iuse_b, idom_b, vmsk, igp, itima, ispace,               &
                      stsh_hours, lts0_t, lnetcdf_u

USE nlstcall_mod, ONLY: model_basis_time, lmean, lcal360

USE ppxlook_mod, ONLY: exppxi
USE totimp_mod,  ONLY: totimp

USE umPrintMgr, ONLY: umPrint, umMessage, newline
USE conversions_mod, ONLY: isec_per_day

USE errormessagelength_mod, ONLY: errormessagelength

USE levsrt_mod, ONLY: levsrt

IMPLICIT NONE

!  Description:
!  Constructs a preliminary STASH list of user requests. Uses interim
!  pointer system, by means of the "extra entry" NELEMP+1 in the LIST_S
!  array. At this stage, the input levels encompass all possible levels.
!  Called by STPROC.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

!  Code description:
!    FORTRAN 90
!    Written to UMDP3 programming standards

! Subroutine arguments

!   Scalar arguments with intent(out):

INTEGER :: nrecs
INTEGER :: ntimes
INTEGER :: nlevels ! Total no. of sets of levs for diags (inpt+outp)
CHARACTER(LEN=errormessagelength) :: cmessage

! ErrorStatus:
INTEGER :: errorstatus

! Local scalars:
LOGICAL ::   model_lev
LOGICAL ::   lmask
LOGICAL ::   loffset
INTEGER ::   i
INTEGER ::   ibot1
INTEGER ::   idiag
INTEGER ::   idomlev
INTEGER ::   idom_l
LOGICAL ::   ldum
INTEGER ::   ifirst
INTEGER ::   ifirst1
INTEGER ::   ilast
INTEGER ::   ilast1
INTEGER ::   im
INTEGER ::   imd
INTEGER ::   iplof
INTEGER ::   modl_l
INTEGER ::   isec_l
INTEGER ::   item_l
INTEGER ::   itim_l
INTEGER ::   itim
INTEGER ::   itop1
INTEGER ::   iuse_l
INTEGER ::   ix1
INTEGER ::   ix2
INTEGER ::   iy1
INTEGER ::   iy2
INTEGER ::   jlev
INTEGER ::   lev_offset
INTEGER ::   lbvc
INTEGER ::   imax          ! to find max of times-table
INTEGER ::   itimlst       ! column of times-table
INTEGER ::   item_chk
INTEGER ::   river_step_ts     ! river routing period in timesteps
INTEGER ::   phenol_period_ts  ! phenol period in timesteps
INTEGER ::   triffid_period_ts ! triffid period in timesteps
INTEGER ::   start_days,start_secs
INTEGER ::   end_days,end_secs
INTEGER ::   ref_days,ref_secs
INTEGER ::   period

INTEGER :: offset_size ! offset for diagnostics not output every timestep

CHARACTER (LEN=errormessagelength)          :: cmessage2
CHARACTER (LEN=*), PARAMETER :: RoutineName='PRELIM'

CHARACTER (LEN=30) :: FieldStr

! Function and subroutine calls:
LOGICAL :: disct_lev

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of Header ------------------------------------------------------

! river_step is in seconds - calculate in timesteps

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
river_step_ts = river_step*(REAL(steps_per_periodim(atmos_im))                &
                  /REAL(secs_per_periodim(atmos_im)))

! triffid_period and phenol_period are in days - calculate in timesteps
triffid_period_ts = triffid_period * (REAL(steps_per_periodim(atmos_im))      &
     / REAL(secs_per_periodim(atmos_im))) * isec_per_day

phenol_period_ts = phenol_period * (REAL(steps_per_periodim(atmos_im))        &
     / REAL(secs_per_periodim(atmos_im))) * isec_per_day

! 0.1  Store output-times tables in array ITIM_S

IF (ntimes == 0) THEN
  DO i = 1, nproftp
    IF (iopt_t(i) == 2 .AND. modl_t(i) >  0) THEN

      ! Profile has output times list
      !  MODL_T(I) labels internal model for times list
      DO itim = 1, itim_t(i)

        itim_s(itim,i) = totimp(iser_t(itim,i), unt3_t(i), modl_t(i))
        IF (itim_s(itim,i)  ==  -999) THEN
          errorstatus = 100
          WRITE (cmessage,'(a,a,i3)')                                         &
            'PRELIM:TOTIMP:Error in time period conversion',                  &
            ' output times table no.=', i
          WRITE(umMessage,*) cmessage
          CALL umPrint(umMessage,src='prelim')
          GO TO 9999
        END IF
      END DO
      itim_s(itim_t(i)+1,i) = -1
    ELSE
      itim_s(1,i) = -1
    END IF
  END DO
  ntimes = nproftp
END IF

! 0.2  Store output levels lists in array LEVLST_S

lev_offset = nlevels ! Initialised to 0 before entering this routine

! Loop over domain profiles in STASH basis file
DO i = 1, ndprof
  IF (levb_d(i) == -1) THEN

    ! There is a levels list for this dom prof
    IF (iopl_d(i) == st_levels_model_rho   .OR.                                &
        iopl_d(i) == st_levels_model_theta .OR.                                &
        iopl_d(i) == st_levels_deep_soil) THEN

      ! Levs list contains model levs - list type is integer
      llistty(i+lev_offset) = 'I'
    ELSE

      ! Not model levs - list type real
      llistty(i+lev_offset) = 'R'
    END IF

    ! LEVT_D(I) = no. of levs in list 'I'
    levlst_s(1, i+lev_offset) = levt_d(i)

    ! Levels list 'I' was read into (R)LEVLST_D(J,I), J=1,LEVT_D(I),
    !  by RDBASIS.
    !  Transfer this levels list to (R)LEVLST_S(J,I+LEV_OFFSET),
    !  J=2,LEVT_D(I)+1.

    DO jlev = 1, levt_d(i)
      IF (iopl_d(i) == st_levels_model_rho   .OR.                              &
          iopl_d(i) == st_levels_model_theta .OR.                              &
          iopl_d(i) == st_levels_deep_soil) THEN
  
        !         Model levels
        levlst_s(jlev+1, i+lev_offset) = levlst_d(jlev, i)
      ELSE IF (iopl_d(i) /= st_levels_single) THEN

        !         Real levels
        rlevlst_s(jlev+1, i+lev_offset) = rlevlst_d(jlev, i)
      END IF
    END DO

    iplof = i + lev_offset

    !   Sort this levels list into correct order (if not already in order)
    CALL levsrt( llistty(   iplof),  levlst_s(1, iplof),                      &
                levlst_s(2, iplof), rlevlst_s(2, iplof))
  ELSE

    ! No levels list, i.e., the output from this diag. is on a
    !    contiguous range of model levels
    levlst_s(1, i+lev_offset) = 0
    rlevlst_s(1, i+lev_offset) = 0
  END IF
END DO  !  Domain profiles

nlevels = ndprof+lev_offset  ! NDPROF = no. of sets of input levels

IF (nlevels > nlevlstsp) THEN
  WRITE(umMessage,'(A)')                                                      &
    'PRELIM: Too many levels lists, arrays overwritten.'
  CALL umPrint(umMessage,src='prelim')
  cmessage =                                                                  &
    'PRELIM: Too many levels lists, arrays overwritten.'
  GO TO 9999
END IF

! Section 1. MAIN LOOP - loop over diag requests in STASH basis file

IF (ndiag > 0) THEN

  DO idiag = 1, ndiag

    modl_l = modl_b(idiag)
    isec_l = isec_b(idiag)
    item_l = item_b(idiag)
    idom_l = idom_b(idiag)
    iuse_l = iuse_b(idiag)
    itim_l = itim_b(idiag)

    item_chk = 0

    item_chk = exppxi(modl_l, isec_l, item_l, ppx_item_number,                &
                      errorstatus, cmessage)

    WRITE(FieldStr,'(2(A,I0))')                                               &
      'Field - Section:', isec_l, ', Item:', item_l

    IF (item_chk /= item_l) THEN
      WRITE(cmessage2,'(A)')                                         newline//&
        TRIM(ADJUSTL(FieldStr)) // ' discarded.'//                   newline//&
        'No stashmaster record.'

      errorstatus = -10

      CALL ereport(routinename, errorstatus, cmessage2)

      ! Make this item null
      item_b(idiag) = 0
      GO TO 999
    END IF
    IF (itim_l /= 0) THEN       ! If the diag is not a null request

      ! Section 1.0  Extract data required for STASH processing from PPXI

      IF (nrecs == nrecdp) THEN
        WRITE(cmessage2,'(A)')                                       newline//&
          TRIM(ADJUSTL(FieldStr)) // ' request denied.'//            newline//&
          'Too many stash list entries.'

        errorstatus = -20

        CALL ereport(routinename, errorstatus, cmessage2)
        GO TO 999
      END IF

      vmsk    = exppxi(modl_l ,isec_l ,item_l,ppx_version_mask ,              &
                       errorstatus, cmessage)
      ispace  = exppxi(modl_l ,isec_l ,item_l,ppx_space_code   ,              &
                       errorstatus, cmessage)
      itima   = exppxi(modl_l ,isec_l ,item_l,ppx_timavail_code,              &
                       errorstatus, cmessage)
      igp     = exppxi(modl_l ,isec_l ,item_l,ppx_grid_type    ,              &
                       errorstatus, cmessage)
      ilev    = exppxi(modl_l ,isec_l ,item_l,ppx_lv_code      ,              &
                       errorstatus, cmessage)
      ibot    = exppxi(modl_l ,isec_l ,item_l,ppx_lb_code      ,              &
                       errorstatus, cmessage)
      itop    = exppxi(modl_l ,isec_l ,item_l,ppx_lt_code      ,              &
                       errorstatus, cmessage)
      iflag   = exppxi(modl_l ,isec_l ,item_l,ppx_lev_flag     ,              &
                       errorstatus, cmessage)
      DO i = 1, 6

        iopn(i) = exppxi(modl_l ,isec_l ,item_l, ppx_opt_code+i-1 ,           &
                         errorstatus, cmessage)
      END DO
      ipseudo = exppxi(modl_l ,isec_l ,item_l,ppx_pt_code         ,           &
                       errorstatus, cmessage)
      ipfirst = exppxi(modl_l ,isec_l ,item_l,ppx_pf_code         ,           &
                       errorstatus, cmessage)
      iplast  = exppxi(modl_l ,isec_l ,item_l,ppx_pl_code         ,           &
                       errorstatus, cmessage)
      ptr_prog= exppxi(modl_l ,isec_l ,item_l,ppx_ptr_code        ,           &
                       errorstatus, cmessage)
      lbvc    = exppxi(modl_l ,isec_l ,item_l,ppx_lbvc_code       ,           &
                       errorstatus, cmessage)

      ! Check availability of diagnostic
      ! DEPENDS ON: tstmsk
      CALL tstmsk(modl_l, isec_l, lmask, ldum, errorstatus, cmessage)
      IF (.NOT. lmask) THEN
        WRITE(cmessage2,'(A)')                                       newline//&
          TRIM(ADJUSTL(FieldStr)) // ' request denied.' //           newline//&
          'Unavailable to this model version.'

        errorstatus = -30

        CALL ereport(routinename, errorstatus, cmessage2)
        GO TO 999
      END IF

      nrecs = nrecs+1

      list_s(st_diag_address, nrecs) = idiag
      list_s(st_model_code  , nrecs) = modl_l
      list_s(st_sect_no_code, nrecs) = isec_l
      list_s(st_item_code   , nrecs) = item_l

      ! Prelim pointer for 'child' records
      list_s(nelemp+1       , nrecs) = nrecs
      list_s(st_lookup_ptr  , nrecs) = -1

      ! Set input code for STASH requests:
      !  =0 Use primary or secondary field:       D1(SI(item,section,model))
      !  =1 Use field in diagnostic space: STASHwork(SI(item,section,model))
      !  =-j Use diagnostic at D1(LIST_S(st_output_addr,j))
      IF ( (ispace == 2) .OR. (ispace == 4)                                   &
        .OR. (ispace == 7) .OR. (ispace == 8) .OR. (ispace == 9)) THEN
        list_s(st_input_code, nrecs) = 0
      ELSE
        list_s(st_input_code, nrecs) = 1
      END IF


      ! 1.1   Expand the domain profile ---------------------------

      !   Averaging and Weighting
      im = imsk_d(idom_l)
      IF ( (igp ==  2) .OR. (igp ==  3)    .OR.                               &
           (igp == 12) .OR. (igp == 13)) THEN

        ! Diags only available over land/sea
        IF ( (imsk_d(idom_l)  ==  1)         .AND.                            &
             (igp == 3 .OR. igp == 13) ) THEN

          ! Diag requested over land+sea, only available over sea
          im = 3
        ELSE IF ( (imsk_d(idom_l)  ==  1)     .AND.                           &
                  (igp == 2 .OR. igp == 12) ) THEN

          ! Diag requested over land+sea, only available over land
          im = 2
        ELSE IF ( (imsk_d(idom_l)  ==  2)    .AND.                            &
                (igp == 3 .OR. igp == 13) ) THEN

          ! Diag requested over land, only available over sea
          WRITE(umMessage,'(A)')                                              &
            'PRELIM: Changed to sea diagnostic. ' // TRIM(ADJUSTL(FieldStr))
          CALL umPrint(umMessage,src='prelim')
          im = 3
        ELSE IF ( (imsk_d(idom_l)  ==  3)    .AND.                            &
                  (igp == 2 .OR. igp == 12) ) THEN

          ! Diag requested over sea, only available over land
          WRITE(umMessage,'(A)')                                              &
            'PRELIM: Changed to land diagnostic. ' // TRIM(ADJUSTL(FieldStr))
          CALL umPrint(umMessage,src='prelim')
          im = 2
        END IF
      END IF

      IF ( (imsk_d(idom_l)  ==  4) ) THEN

        ! Diag requested for processing of ALL non-mdi values
        WRITE(umMessage,'(A)')                                                &
          'PRELIM: All non-mdi diagnostic. ' // TRIM(ADJUSTL(FieldStr))
        CALL umPrint(umMessage,src='prelim')
        im = 4
      END IF

      list_s(st_gridpoint_code, nrecs) = im+10*imn_d(idom_l)
      list_s(st_weight_code   , nrecs) =       iwt_d(idom_l)
      list_s(st_domain_code   , nrecs) =      iopa_d(idom_l)

      !   Horizontal area
      !    - convert lat/long spec to row/column numbers if appropriate;
      !    - convert lat/long spec to equatorial lat/long if appropriate.
      IF (iopa_d(idom_l) == st_domain_global) THEN

        ! DEPENDS ON: lltorc
        CALL lltorc(igp, 90, -90, 0, 360,                                     &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),     &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == st_domain_n_hemisphere ) THEN

        ! DEPENDS ON: lltorc
        CALL lltorc(igp, 90, 0, 0, 360,                                       &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),     &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == st_domain_s_hemisphere ) THEN

        ! DEPENDS ON: lltorc
        CALL lltorc(igp, 0, -90, 0, 360,                                      &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),     &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == st_domain_30_to_90_N ) THEN

        ! DEPENDS ON: lltorc
        CALL lltorc(igp, 90, 30, 0, 360,                                      &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),     &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == st_domain_30_to_90_S ) THEN

        ! DEPENDS ON: lltorc
        CALL lltorc(igp, -30, -90, 0, 360,                                    &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),     &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == st_domain_0_to_30_N ) THEN

        ! DEPENDS ON: lltorc
        CALL lltorc(igp, 30, 00, 0, 360,                                      &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),     &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == st_domain_0_to_30_S ) THEN

        ! DEPENDS ON: lltorc
        CALL lltorc(igp, 00, -30, 0, 360,                                     &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),     &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == st_domain_30_S_to_30_N ) THEN

        ! DEPENDS ON: lltorc
        CALL lltorc(igp, 30, -30, 0, 360,                                     &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),     &
              list_s(st_west_code, nrecs), list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == st_domain_whole_degrees ) THEN

        ! DEPENDS ON: lltorc
        CALL lltorc(igp, inth_d(idom_l), isth_d(idom_l),                      &
              iwst_d(idom_l), iest_d(idom_l),                                 &
              list_s(st_south_code, nrecs), list_s(st_north_code, nrecs),     &
              list_s(st_west_code, nrecs),  list_s(st_east_code, nrecs))
      ELSE IF (iopa_d(idom_l) == st_domain_gridpoints ) THEN

        ! DEPENDS ON: lltorc
        CALL lltorc(igp, 90, -90, 0, 360, iy1, iy2, ix1, ix2)
        list_s(st_north_code, nrecs) = MIN(inth_d(idom_l), iy2)
        list_s(st_south_code, nrecs) = MIN(isth_d(idom_l), iy2)
        list_s(st_west_code , nrecs) = MIN(iwst_d(idom_l), ix2)
        list_s(st_east_code , nrecs) = MIN(iest_d(idom_l), ix2)
      ELSE
        WRITE(cmessage2,'(A,I0)')                                    newline//&
          TRIM(ADJUSTL(FieldStr)) //                                 newline//&
          'Invalid domain area option = ', iopa_d(idom_l)

        errorstatus = -35

        CALL ereport(routinename, errorstatus, cmessage2)
        nrecs = nrecs-1
        GO TO 999
      END IF

      ! Input level setting
      ! DEPENDS ON: disct_lev
      model_lev = disct_lev(ilev, errorstatus, cmessage)
      IF (model_lev) THEN

        ! Model levels
        ! Set bottom level
        ! DEPENDS ON: levcod
        CALL levcod(ibot, ibot1, errorstatus, cmessage)

        ! Set top level
        ! DEPENDS ON: levcod
        CALL levcod(itop, itop1, errorstatus, cmessage)

        ! Contig. range of model levels
        IF (iflag == 0) THEN
          list_s(st_input_bottom, nrecs) = ibot1
          list_s(st_input_top   , nrecs) = itop1

          ! Non-contig. levels list
        ELSE IF (iflag == 1) THEN
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  1
        END IF
      ELSE

        ! Non-model levels
        IF (ilev == 3) THEN

          !  Pressure levels
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  2
        ELSE IF (ilev == 4) THEN

          !  Height levels
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  3
        ELSE IF (ilev == 5) THEN

          !  Special levels
          list_s(st_input_bottom, nrecs) = 100
          list_s(st_input_top   , nrecs) = lbvc
        ELSE IF (ilev == 7) THEN

          !  Theta levels
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  4
        ELSE IF (ilev == 8) THEN

          !  PV levels
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  5
        ELSE IF (ilev == 9) THEN

          !  Cloud threshold levels
          list_s(st_input_bottom, nrecs) = -1
          list_s(st_input_top   , nrecs) =  6
        END IF
      END IF

      ! Output level specification
      ! DEPENDS ON: disct_lev
      model_lev = disct_lev(ilev, errorstatus, cmessage)
      IF (model_lev) THEN

        ! Model levels
        IF (levb_d(idom_l) >= 0) THEN

          ! Contiguous range of model levels
          list_s(st_output_bottom, nrecs) = MAX(levb_d(idom_l), ibot1)
          list_s(st_output_top   , nrecs) = MIN(levt_d(idom_l), itop1)
          IF ( (levb_d(idom_l) <  ibot1)   .OR.                               &
               (levt_d(idom_l) >  itop1) ) THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) //                             newline//&
              'Field has level range out of bounds; Corrected'
            errorstatus = -40

            CALL ereport(routinename, errorstatus, cmessage2)
          END IF
          IF ( (  ts_d(idom_l) )      .AND.                                   &
               ( (levb_d(idom_l) <  ibot1)    .OR.                            &
                 (levt_d(idom_l) >  itop1) )  ) THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Time series domain has inconsistent levels.'

            errorstatus = -50

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs-1
            GO TO 999
          END IF
          IF ( (levt_d(idom_l) <  ibot1)  .OR.                                &
               (levb_d(idom_l) >  itop1)) THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Field has inconsistent top and bottom levels.'

            errorstatus = -60

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs-1
            GO TO 999
          END IF
        ELSE

          ! Non-contig. list of model levels
          list_s(st_output_bottom, nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   , nrecs) = 1
        END IF
      ELSE

        ! Non-model levels
        IF (ilev == 5) THEN

          ! Special level
          list_s(st_output_bottom, nrecs) = 100
          list_s(st_output_top   , nrecs) = lbvc
        ELSE IF (ilev == 3) THEN

          ! Pressure levels
          list_s(st_output_bottom,nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   ,nrecs) = 2
        ELSE IF (ilev == 4) THEN

          ! Height levels
          list_s(st_output_bottom, nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   , nrecs) = 3
        ELSE IF (ilev == 7 ) THEN

          ! Theta levels
          list_s(st_output_bottom, nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   , nrecs) = 4
        ELSE IF (ilev == 8 ) THEN

          ! PV levels
          list_s(st_output_bottom, nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   , nrecs) = 5
        ELSE IF (ilev == 9 ) THEN

          ! Cloud threshold levels
          list_s(st_output_bottom, nrecs) = -(idom_l + lev_offset)
          list_s(st_output_top   , nrecs) = 6
        ELSE
          WRITE(cmessage2,'(A,I0)')                                  newline//&
            TRIM(ADJUSTL(FieldStr)) // ' ignored.'//                 newline//&
            'Domain level option = ', iopl_d(idom_l)

          errorstatus=-70

          CALL ereport(routinename, errorstatus, cmessage2)
          nrecs = nrecs-1
          GO TO 999
        END IF
      END IF

      ! Output pseudo-levels level setting
      IF (ipseudo /= plt_d(idom_l)) THEN
        WRITE(cmessage2,'(A)')                                       newline//&
          TRIM(ADJUSTL(FieldStr)) // ' ignored.'//                   newline//&
          'Invalid pseudo-level type.'

        errorstatus = -80

        CALL ereport(routinename, errorstatus, cmessage2)
        nrecs = nrecs-1
        GO TO 999
      END IF
      list_s(st_pseudo_in,nrecs) = 0  !(This is set in INPUTL)
      IF (ipseudo >  0) THEN

        ! Pseudo levels list for this diagnostic
        list_s(st_pseudo_out, nrecs) = plpos_d(idom_l)
        lenplst(plpos_d(idom_l))     = pllen_d(idom_l)
        ifirst = pslist_d(1, plpos_d(idom_l))
        ilast  = pslist_d(pllen_d(idom_l), plpos_d(idom_l))

        ! Check pseudo level limits
        ! DEPENDS ON: pslims
        CALL pslims(ipfirst, iplast, ifirst1, ilast1)
        IF (ifirst <  ifirst1) THEN
          WRITE(cmessage2,'(A)')                                     newline//&
            TRIM(ADJUSTL(FieldStr)) // ' ignored.'//                 newline//&
            'Field has first pseudo-level set too low.'

          errorstatus = -90

          CALL ereport(routinename, errorstatus, cmessage2)
          nrecs = nrecs-1
          GO TO 999
        END IF
        IF (ilast >  ilast1) THEN
          WRITE(cmessage2,'(A)')                                     newline//&
            TRIM(ADJUSTL(FieldStr)) // ' ignored.'//                 newline//&
            'Field has last pseudo-level set too high.'

          errorstatus = -95

          CALL ereport(routinename, errorstatus, cmessage2)
          nrecs = nrecs-1
          GO TO 999
        END IF
      ELSE
        list_s(st_pseudo_out,nrecs) = 0
      END IF

      ! Time-series domain profiles
      IF (ts_d(idom_l) ) THEN

        ! Pointer for location of time series
        list_s(st_series_ptr, nrecs) = npos_ts(idom_l)
      ELSE
        list_s(st_series_ptr, nrecs) = 0
      END IF

      ! 1.2   Expand the useage profile --------------------------

      IF (locn_u(iuse_l) == 3) THEN                ! PP file

        list_s(st_output_code, nrecs) = -iunt_u(iuse_l)
        list_s(st_macrotag,    nrecs) = 0

        IF (lnetcdf_u(iuse_l)) THEN
          list_s(st_output_type, nrecs) = 2
        ELSE
          list_s(st_output_type, nrecs) = 1
        END IF

      ELSE IF (locn_u(iuse_l) == 1) THEN ! Dump store: set user tag

        list_s(st_output_code, nrecs) = 1
        list_s(st_macrotag,    nrecs) = iunt_u(iuse_l)
        list_s(st_output_type, nrecs) = 0

      ELSE IF (locn_u(iuse_l) == 6) THEN ! Secondary dump store:
                                         !             set user tag
        list_s(st_output_code, nrecs) = 2
        list_s(st_macrotag,    nrecs) = iunt_u(iuse_l)
        list_s(st_output_type, nrecs) = 0

      ELSE IF (locn_u(iuse_l) == 2) THEN ! Climate mean: tag set
                                         !   1000*(time mean tag)

        list_s(st_output_code, nrecs) = 1
        list_s(st_macrotag,    nrecs) = iunt_u(iuse_l) * 1000
        list_s(st_output_type, nrecs) = 0

      ELSE

        WRITE(cmessage2,'(A,I0)')                                   newline//&
          TRIM(ADJUSTL(FieldStr)) // ' ignored.'//                  newline//&
          'Invalid usage option = ', locn_u(iuse_l)

        errorstatus = -120

        CALL ereport(routinename, errorstatus, cmessage2)
        nrecs = nrecs-1
        GO TO 999

      END IF

      ! 1.3   Expand the time profile ------------------------------

      ! Initialise as single time field

      !   Set time processing record

      list_s(st_proc_no_code, nrecs) = ityp_t(itim_l)

      ! Set units for sampling period (frequency)
      list_s(st_time_unit2_code, nrecs) = unt2_t(itim_l)

      ! Set units for output times of diagnostic
      list_s(st_time_unit3_code, nrecs) = unt3_t(itim_l)

      ! Initialise offset to 0
      list_s(st_offset_code,nrecs)=0

      ! Time series currently not supported for netCDF output
      IF (list_s(st_output_type, nrecs) == 2 .AND. &
           (list_s(st_series_ptr, nrecs) > 0 .OR. &
            list_s(st_proc_no_code, nrecs) == 4 .OR. &
            list_s(st_proc_no_code, nrecs) == 8)) THEN

        WRITE(cmessage2,'(A)')                                      newline//&
          TRIM(ADJUSTL(FieldStr)) // ' ignored.'//                  newline//&
          'Time series currently not supported for netCDF output.'

        errorstatus = -130

        CALL ereport(routinename, errorstatus, cmessage2)
        nrecs = nrecs-1
        CYCLE

      END IF

      !   Set period record
      IF (ityp_t(itim_l) == 1) THEN        ! No period
        list_s(st_period_code, nrecs) = 0
      ELSE IF ( (intv_t(itim_l) == -1)  .AND.                                 &
                (ityp_t(itim_l) == 2) ) THEN
        list_s(st_period_code, nrecs) = -1
      ELSE
        list_s(st_period_code, nrecs) =                                       &
                  totimp(intv_t(itim_l), unt1_t(itim_l), modl_l)
        IF (list_s(st_period_code,nrecs)  ==  -999) THEN
          errorstatus = 101
          WRITE (cmessage,'(A)')                                              &
            TRIM(ADJUSTL(FieldStr)) //                               newline//&
            'PRELIM:TOTIMP:Error in time period conversion'
          WRITE(umMessage,*) cmessage
          CALL umPrint(umMessage,src='prelim')
          GO TO 9999
        END IF
      END IF

      IF (iopt_t(itim_l) == 1 .OR. iopt_t(itim_l) == 3) THEN

        ! Regular output times
        list_s(st_freq_code, nrecs)=                                          &
                  totimp(ifre_t(itim_l), unt3_t(itim_l), modl_l)
        IF (list_s(st_freq_code, nrecs)  ==  -999) THEN
          errorstatus = 102
          WRITE (cmessage,'(A)') TRIM(ADJUSTL(FieldStr)) //          newline//&
            'PRELIM:TOTIMP:Error in time period conversion'
          WRITE(umMessage,*) cmessage
          CALL umPrint(umMessage,src='prelim')
          GO TO 9999
        END IF

        ! Time profile is specified with a start and end date, between which,
        ! the diagnostic is written out.
        IF (iopt_t(itim_l) == 3) THEN

          ! Calculate the difference in (integer number of) hours between the
          ! start date in the time profile and the model basis time.
          ! DEPENDS ON: time2sec
          CALL time2sec(isdt_t(1, itim_l), isdt_t(2, itim_l),                 &
                        isdt_t(3, itim_l), isdt_t(4, itim_l),                 &
                        isdt_t(5, itim_l), isdt_t(6, itim_l),                 &
                        0, 0, start_days, start_secs, lcal360)

          ! DEPENDS ON: time2sec
          CALL time2sec(model_basis_time(1), model_basis_time(2),             &
                        model_basis_time(3), model_basis_time(4),             &
                        model_basis_time(5), model_basis_time(6),             &
                        0, 0, ref_days, ref_secs, lcal360)

          period = (start_days-ref_days) * 24 +                               &
                       INT((start_secs-ref_secs) / 3600.0)

          ! Convert the start date (hours) into timesteps since basis time.
          list_s(st_start_time_code, nrecs) = totimp(period,stsh_hours, modl_l)

        ELSE

          list_s(st_start_time_code, nrecs) =                                 &
                      totimp(istr_t(itim_l), unt3_t(itim_l), modl_l)
        END IF
        IF (list_s(st_start_time_code, nrecs)  ==  -999) THEN
          errorstatus = 103
          WRITE (cmessage,'(A)') TRIM(ADJUSTL(FieldStr)) //          newline//&
            'PRELIM:TOTIMP:Error in time period conversion'
          WRITE(umMessage,*) cmessage
          CALL umPrint(umMessage,src='prelim')
          GO TO 9999
        END IF
        IF (iend_t(itim_l) == -1) THEN
          list_s(st_end_time_code, nrecs) = -1
        ELSE
          IF (iopt_t(itim_l) == 3) THEN

            ! Calculate the difference in (integer number of) hours between the
            ! end date in the time profile and the model basis time.
            ! DEPENDS ON: time2sec
            CALL time2sec(iedt_t(1, itim_l), iedt_t(2, itim_l),               &
                          iedt_t(3, itim_l), iedt_t(4, itim_l),               &
                          iedt_t(5, itim_l), iedt_t(6, itim_l),               &
                          0, 0, end_days, end_secs, lcal360)

            ! DEPENDS ON: time2sec
            CALL time2sec(model_basis_time(1), model_basis_time(2),           &
                          model_basis_time(3), model_basis_time(4),           &
                          model_basis_time(5), model_basis_time(6),           &
                          0, 0, ref_days, ref_secs, lcal360)

            ! End date uses ceiling to calculate the next largest integer secs
            ! so that the end date specified in the namelist is inclusive.
            period = (end_days-ref_days) * 24 +                               &
                       CEILING((end_secs-ref_secs) / 3600.0)

            ! Convert the end date (hours) into timesteps since basis time.
            list_s(st_end_time_code, nrecs)= totimp(period,stsh_hours, modl_l)
          ELSE

            list_s(st_end_time_code, nrecs)=                                  &
                          totimp(iend_t(itim_l), unt3_t(itim_l), modl_l)
          END IF
        END IF
        IF (list_s(st_end_time_code, nrecs)  ==  -999) THEN
          errorstatus = 104
          WRITE (cmessage,'(A)') TRIM(ADJUSTL(FieldStr)) //          newline//&
            'PRELIM:TOTIMP:Error in time period conversion'
          WRITE(umMessage,*) cmessage
          CALL umPrint(umMessage,src='prelim')
          GO TO 9999
        END IF

        ! Set end time to -1 if output requested to end of run

        ! Correct start time for radiation, periodic convection, leaf
        ! phenology and vegetation competition and river routing
        ! most of these diagnostics are output at the beginning of
        ! timestep so:
        offset_size = 1
        IF ( (itima == 2) .AND. (a_lw_radstep_diag /= 1) ) THEN
          imd = MOD( list_s(st_start_time_code, nrecs), a_lw_radstep_diag )
          list_s(st_start_time_code, nrecs) =                                 &
            list_s(st_start_time_code, nrecs) + 1 - imd
          loffset = .TRUE.
        ELSE IF ( (itima == 3) .AND. (a_sw_radstep_diag /= 1) ) THEN
          imd = MOD( list_s(st_start_time_code, nrecs), a_sw_radstep_diag )
          list_s(st_start_time_code, nrecs) =                                 &
            list_s(st_start_time_code, nrecs) + 1 - imd
          loffset = .TRUE.
        ELSE IF ( (itima == 13) .AND. (a_conv_step /= 1) ) THEN
          imd = MOD( list_s(st_start_time_code, nrecs), a_conv_step )
          list_s(st_start_time_code, nrecs) =                                 &
            list_s(st_start_time_code, nrecs) + 1 - imd
          loffset = .TRUE.
        ELSE IF ( (itima == 14) .AND. (phenol_period_ts /= 1) ) THEN
          offset_size = phenol_period_ts
          imd = MOD(list_s(st_start_time_code, nrecs), phenol_period_ts)
          list_s(st_start_time_code, nrecs) =                                 &
            list_s(st_start_time_code, nrecs) - imd + offset_size
          loffset = .TRUE.
        ELSE IF ( (itima == 15) .AND. (triffid_period_ts /= 1) ) THEN
          offset_size = triffid_period_ts
          imd = MOD(list_s(st_start_time_code, nrecs), triffid_period_ts)
          list_s(st_start_time_code, nrecs)=                                  &
            list_s(st_start_time_code, nrecs) - imd + offset_size
          loffset = .TRUE.
        ELSE IF ( (itima == 16) .AND. (river_step_ts /= 1)) THEN
          ! river routing diagnostics are output at end of river routing
          ! timestep
          offset_size = river_step_ts
          imd=MOD(list_s(st_start_time_code,nrecs),river_step_ts)
          list_s(st_start_time_code,nrecs)=                                   &
            list_s(st_start_time_code,nrecs)-imd+offset_size
          loffset=.TRUE.
        ELSE IF ((itima == 17) .AND. (interval /= 1)) THEN
          offset_size = interval
          imd = MOD( list_s(st_start_time_code, nrecs), interval)
          list_s(st_start_time_code, nrecs) =                                 &
             list_s(st_start_time_code, nrecs)-imd+offset_size
          loffset = .TRUE.
        ELSE IF ((itima == 18) .AND. (a_lw_radstep_prog /= 1)) THEN
          imd = MOD( list_s(st_start_time_code, nrecs), a_lw_radstep_prog)
          list_s(st_start_time_code, nrecs) =                                 &
             list_s(st_start_time_code, nrecs) + 1 - imd
          loffset = .TRUE.
        ELSE IF ((itima == 19) .AND. (a_sw_radstep_prog /= 1)) THEN
          imd = MOD( list_s(st_start_time_code, nrecs), a_sw_radstep_prog )
          list_s(st_start_time_code, nrecs) =                                 &
             list_s(st_start_time_code, nrecs) + 1 - imd
          loffset = .TRUE.
        ELSE
          loffset = .FALSE.
        END IF
      ELSE IF (iopt_t(itim_l) == 2) THEN

        ! List of specified output times
        list_s(st_freq_code, nrecs) = -itim_l
      ELSE
        WRITE(cmessage2,'(A)')                                       newline//&
          TRIM(ADJUSTL(FieldStr)) // ' ignored.'//                   newline//&
          'Invalid output times code.'

        errorstatus = -150

        CALL ereport(routinename, errorstatus, cmessage2)
        nrecs = nrecs-1
        GO TO 999
      END IF

      ! Set flag for PP output at TimeStep 0
      IF (lts0_t(itim_l)) THEN
        list_s(st_output_ts0_code, nrecs) = 0
      ELSE
        list_s(st_output_ts0_code, nrecs) = 1
      END IF

      IF ( (list_s(st_proc_no_code, nrecs) >  1)   .AND.                      &
           (list_s(st_proc_no_code, nrecs) <= 6) ) THEN

        ! Other than single time field
        IF (nrecs >= nrecdp) THEN
          WRITE(cmessage2,'(A)')                                     newline//&
            TRIM(ADJUSTL(FieldStr)) // ' ignored.'//                 newline//&
            'Too many s_list requests.'

          errorstatus = -160

          CALL ereport(routinename, errorstatus, cmessage2)
          nrecs = nrecs-1
          GO TO 999
        END IF

        DO i = 1, nelemp + 1          ! Copy stash list forward
          list_s(i, nrecs+1) = list_s(i, nrecs)
        END DO

        IF (loffset) THEN
          ! Rad, conv, or river routing timesteps,
          ! offset already added
          list_s(st_start_time_code, nrecs+1) =                               &
             list_s(st_start_time_code, nrecs+1) - offset_size
          IF (list_s(st_period_code,nrecs) /= -1) THEN

            ! Offsets are added to start time
            list_s(st_offset_code, nrecs) =                                   &
                          totimp(ioff_t(itim_l), unt2_t(itim_l), modl_l)
            IF (list_s(st_offset_code, nrecs)  ==  -999) THEN
              errorstatus = 1
              WRITE (cmessage,'(A)') TRIM(ADJUSTL(FieldStr)) //      newline//&
                'PRELIM:TOTIMP:Error in time period conversion'
              WRITE(umMessage,*) cmessage
              CALL umPrint(umMessage,src='prelim')
              GO TO 9999
            END IF
            list_s(st_start_time_code, nrecs) =                               &
              list_s(st_start_time_code, nrecs) -                             &
              list_s(st_period_code, nrecs)     +                             &
              list_s(st_offset_code, nrecs)
          ELSE
            list_s(st_start_time_code,nrecs) = 1
          END IF

        ELSE

          IF (list_s(st_period_code, nrecs) /= -1) THEN

            ! Offsets are added to start time
            list_s(st_offset_code, nrecs) =                                   &
                          totimp(ioff_t(itim_l), unt2_t(itim_l), modl_l)
            IF (list_s(st_offset_code, nrecs)  ==  -999) THEN
              errorstatus = 1
              WRITE (cmessage,'(A)') TRIM(ADJUSTL(FieldStr)) //      newline//&
                'PRELIM:TOTIMP:Error in time period conversion'
              WRITE(umMessage,*) cmessage
              CALL umPrint(umMessage,src='prelim')
              GO TO 9999
            END IF
            list_s(st_start_time_code, nrecs) =                               &
              list_s(st_start_time_code, nrecs) -                             &
              list_s(st_period_code, nrecs) + 1 +                             &
              list_s(st_offset_code, nrecs)
          ELSE
            list_s(st_start_time_code, nrecs) = 1
          END IF

        END IF

        IF (list_s(st_start_time_code, nrecs) <  1) THEN
          WRITE(cmessage2,'(A)')                                     newline//&
            TRIM(ADJUSTL(FieldStr)) //                               newline//&
            'Start time before period, setting to 1.'
          errorstatus = -170

          CALL ereport(routinename, errorstatus, cmessage2)
          list_s(st_start_time_code, nrecs) = 1
        END IF

        ! Check if offset corresponds to a valid timestep.
        ! If not, reject the diagnostic.
        SELECT CASE (itima)
        CASE (2)   ! Long-Wave Radiation
          IF (MOD( list_s(st_offset_code, nrecs), a_lw_radstep_diag ) /= 0)   &
          THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with lw_rad step.'

            errorstatus = -180

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF

        CASE (3)   ! Short-Wave Radiation
          IF (MOD( list_s(st_offset_code, nrecs), a_sw_radstep_diag ) /= 0)   &
          THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with sw_rad step.'

            errorstatus = -190

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        CASE (13)   ! Convection
          IF (MOD( list_s(st_offset_code, nrecs), a_conv_step ) /= 0) THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with convect step.'

            errorstatus = -200

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF

        CASE (14)   ! Leaf Phenology
          IF (MOD( list_s(st_offset_code, nrecs), phenol_period_ts ) /= 0) THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with phenol period.'

            errorstatus = -210

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF

        CASE (15)   ! Triffid
          IF (MOD( list_s(st_offset_code, nrecs), triffid_period_ts ) /= 0) THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with triffid period timesteps.'

            errorstatus = -220

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF

        CASE (16)   ! Trip river routing model
          IF (MOD(list_s(st_offset_code,nrecs),river_step_ts)/= 0) THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with river routing period.'
            errorstatus = -230

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs=nrecs-1
            GO TO 999
          END IF

        CASE (17)   ! UKCA
          IF (MOD(list_s(st_offset_code,nrecs),interval)/= 0) THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with chemistry period.'
            errorstatus = -230

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs=nrecs-1
            GO TO 999
          END IF

        CASE (18)   ! Long-Wave Radiation Field
          IF (MOD( list_s(st_offset_code, nrecs), a_lw_radstep_prog ) /= 0)   &
          THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with lw_rad step.'

            errorstatus = -180

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF

        CASE (19)   ! Short-Wave Radiation Field
          IF (MOD( list_s(st_offset_code, nrecs), a_sw_radstep_prog ) /= 0)   &
          THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with sw_rad step.'

            errorstatus = -190

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        END SELECT

        list_s(st_proc_no_code ,nrecs+1) = 1

        list_s(st_input_bottom ,nrecs+1) = list_s(st_output_bottom, nrecs)

        list_s(st_input_top    ,nrecs+1) = list_s(st_output_top   , nrecs)

        list_s(st_input_code   ,nrecs+1) = -nrecs
        list_s(st_output_code  ,nrecs)   = 1
        list_s(st_series_ptr   ,nrecs+1) = 0
        list_s(nelemp+1        ,nrecs+1) = nrecs + 1

        ! Frequency
        list_s(st_freq_code,nrecs)=                                           &
                  totimp(isam_t(itim_l), unt2_t(itim_l), modl_l)
        IF (list_s(st_freq_code, nrecs)  ==  -999) THEN
          errorstatus = 105
          WRITE (cmessage,'(A)') TRIM(ADJUSTL(FieldStr)) //          newline//&
            'PRELIM:TOTIMP:Error in time period conversion'

          WRITE(umMessage,*) cmessage
          CALL umPrint(umMessage,src='prelim')
          GO TO 9999
        END IF

        !   Correct frequency for radiation, periodic convection, leaf
        !   phenology and vegetation competition

        IF (itima == 2) THEN
          IF (list_s(st_freq_code, nrecs) == 1) THEN
            list_s(st_freq_code, nrecs) = a_lw_radstep_diag
          ELSE IF (MOD(list_s(st_freq_code, nrecs), a_lw_radstep_diag) /= 0)  &
          THEN
            WRITE(cmessage2,'(A,I0)')                                newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Incorrect sampling for lw_radstep. Frequency = ',              &
              list_s(st_freq_code, nrecs)

            errorstatus = -225

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        ELSE IF (itima == 3) THEN
          IF (list_s(st_freq_code, nrecs) == 1) THEN
            list_s(st_freq_code, nrecs) = a_sw_radstep_diag
          ELSE IF (MOD(list_s(st_freq_code, nrecs), a_sw_radstep_diag) /= 0)  &
          THEN
            WRITE(cmessage2,'(A,I0)')                                newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Incorrect sampling for sw_radstep. Frequency = ',              &
              list_s(st_freq_code, nrecs)

            errorstatus = -230

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        ELSE IF (itima == 13) THEN
          IF (list_s(st_freq_code, nrecs) == 1) THEN
            list_s(st_freq_code, nrecs) = a_conv_step
          ELSE IF (MOD(list_s(st_freq_code, nrecs), a_conv_step) /= 0) THEN
            WRITE(cmessage2,'(A,I0)')                                newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'             //  newline//&
              'Incorrect sampling for conv_step. Frequency = ',               &
              list_s(st_freq_code, nrecs)

            errorstatus = -240

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        ELSE IF (itima == 14) THEN
          IF (list_s(st_freq_code, nrecs) == 1) THEN
            list_s(st_freq_code, nrecs) = phenol_period_ts
          ELSE IF (MOD(list_s(st_freq_code, nrecs), phenol_period_ts) /= 0) THEN
            WRITE(cmessage2,'(A,I0)')                                newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Incorrect sampling for phenol_period (timesteps).'//  newline//&
              'Frequency = ', list_s(st_freq_code,nrecs)

            errorstatus = -250

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        ELSE IF (itima == 15) THEN
          IF (list_s(st_freq_code, nrecs) == 1) THEN
            list_s(st_freq_code, nrecs) = triffid_period_ts
          ELSE IF                                                             &
            (MOD(list_s(st_freq_code, nrecs), triffid_period_ts) /= 0) THEN
            WRITE(cmessage2,'(A,I0)')                                newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Incorrect sampling for triffid_period (timesteps).'// newline//&
              'Frequency = ', list_s(st_freq_code,nrecs)

            errorstatus = -260

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        ELSE IF (itima == 16) THEN
          IF (list_s(st_freq_code,nrecs) == 1) THEN
            list_s(st_freq_code,nrecs)=river_step_ts
          ELSE IF                                                             &
            (MOD(list_s(st_freq_code,nrecs),river_step_ts) /= 0) THEN
            WRITE(cmessage2,'(A,I0)')                                newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Incorrect sampling for river_step_ts.'//              newline//&
              'Frequency = ', list_s(st_freq_code,nrecs)

            errorstatus=-270

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs=nrecs-1
            GO TO 999
          END IF
        ELSE IF (itima == 17) THEN
          IF (list_s(st_freq_code,nrecs) == 1) THEN
            list_s(st_freq_code,nrecs)=interval
          ELSE IF                                                             &
             (MOD(list_s(st_freq_code,nrecs),interval) /= 0) THEN
            WRITE(cmessage2,'(A,I0)')                                newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Incorrect sampling for chemistry timestep. '//        newline//&
              'Frequency = ', list_s(st_freq_code,nrecs)

            errorstatus=-270

            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs=nrecs-1
            GO TO 999
          END IF
        END IF

        ! Correct frequency of the radiation fields, and UKCA diagnostic
        SELECT CASE (itima)
        CASE (17) ! UKCA diagnostics
          IF (MOD( list_s(st_offset_code, nrecs), interval ) /= 0)       &
             THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with chemistry step.'
            errorstatus = -180
            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        CASE (18) ! Long-Wave radiation field
          IF (MOD( list_s(st_offset_code, nrecs), a_lw_radstep_prog ) /= 0)   &
             THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with lw_rad step.'
            errorstatus = -180
            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        CASE (19) ! Short-Wave radiation field
          IF (MOD( list_s(st_offset_code, nrecs), a_sw_radstep_prog ) /= 0)   &
             THEN
            WRITE(cmessage2,'(A)')                                   newline//&
              TRIM(ADJUSTL(FieldStr)) // ' ignored.'//               newline//&
              'Offset does not agree with sw_rad step.'
            errorstatus = -180
            CALL ereport(routinename, errorstatus, cmessage2)
            nrecs = nrecs - 1
            GO TO 999
          END IF
        END SELECT

        ! For the NRECS item an end_time_code needs to be set if we
        ! are dealing with a times table rather than  regular diagn.
        ! This should be the maximum timestep in the time list. The list
        ! should be ready sorted (and thus maximum is last member) but
        ! will run through and find maximum to be on the safe side.

        imax = 0
        itimlst = -list_s(st_freq_code, nrecs+1)

        IF (itimlst  >   0) THEN      ! List *not* regular
          DO i = 1, itim_t(itimlst)
            IF (imax  <   itim_s(i, itimlst)) THEN
              imax = itim_s(i, itimlst)
            END IF
          END DO

          list_s(st_end_time_code, nrecs) = imax

        END IF

        !   Period

        IF ((intv_t(itim_l) == -1) .AND. (ityp_t(itim_l) == 2)) THEN
          list_s(st_period_code, nrecs) = -1
        ELSE
          list_s(st_period_code, nrecs)=                                      &
                      totimp(intv_t(itim_l), unt1_t(itim_l), modl_l)
          IF (list_s(st_period_code, nrecs)  ==  -999) THEN
            errorstatus = 106
            WRITE (cmessage,'(A)') TRIM(ADJUSTL(FieldStr)) //        newline//&
              'PRELIM:TOTIMP:Error in time period conversion'
            WRITE(umMessage,*) cmessage
            CALL umPrint(umMessage,src='prelim')
            GO TO 9999
          END IF
        END IF

        !   Add the record - unless the output destination is the dump,
        !                      and output at the accumulating period
        IF (    locn_u(iuse_l) >  2                                           &
            .OR.                                                              &
             ( (list_s(st_freq_code  , nrecs+1 )    /=                        &
                list_s(st_period_code, nrecs   ) )                            &
                                                   .AND.                      &
               (list_s(st_start_time_code, nrecs+1) /=                        &
                list_s(st_end_time_code  , nrecs+1) )   )  ) THEN

          ! No tag for parent
          list_s(st_macrotag, nrecs) = 0
          nrecs = nrecs + 1
        END IF

      ELSE IF (list_s(st_proc_no_code, nrecs) == 8) THEN

        ! Option of "daily" mean timeseries
        IF (nrecs >= nrecdp) THEN
          WRITE(cmessage2,'(A)')                                     newline//&
            TRIM(ADJUSTL(FieldStr)) // ' ignored.' //                newline//&
            'Too many s_list requests.'

          errorstatus = -270

          CALL ereport(routinename, errorstatus, cmessage2)
          nrecs = nrecs - 1
          GO TO 999
        END IF

        ! Special case where 2 extra records required
        !  Record 1 - time mean only no spatial processing
        !  Record 2 - timeseries formed extracting from record 1
        !  Record 3 - extract timeseries from dump ie record 2

        DO i = 1, nelemp + 1      ! Copy stash list forward
          list_s(i, nrecs+1) = list_s(i, nrecs)
          list_s(i, nrecs+2) = list_s(i, nrecs)
        END DO

        IF (loffset) THEN
          ! Rad, conv, or river routing timesteps,
          ! offset allready added
          list_s(st_start_time_code, nrecs+2) =                               &
             list_s(st_start_time_code, nrecs+2) - offset_size
          IF (list_s(st_period_code, nrecs) /= -1) THEN
            list_s(st_start_time_code, nrecs) =                               &
               list_s(st_start_time_code, nrecs) -                            &
               list_s(st_period_code, nrecs)
          ELSE
            list_s(st_start_time_code, nrecs) = 1
          END IF

        ELSE

          IF (list_s(st_period_code, nrecs) /= -1) THEN
            list_s(st_start_time_code, nrecs) =                               &
              list_s(st_start_time_code, nrecs) -                             &
              list_s(st_period_code, nrecs) + 1
          ELSE
            list_s(st_start_time_code, nrecs) = 1
          END IF

        END IF

        IF (list_s(st_start_time_code, nrecs) <  1) THEN
          WRITE(cmessage2,'(A)')                                     newline//&
            TRIM(ADJUSTL(FieldStr)) //                               newline//&
            'Start time before period, setting to 1'

          errorstatus = -280

          CALL ereport(routinename, errorstatus, cmessage2)
          list_s(st_start_time_code, nrecs) = 1
        END IF

        list_s(st_proc_no_code ,nrecs)   = 3  ! time mean
        list_s(st_proc_no_code ,nrecs+1) = 8  ! timseries special case
        list_s(st_proc_no_code ,nrecs+2) = 1  !  extract

        ! Reset first record to no area weight or spatial processing
        ! ie first record just controls time meaning

        list_s(st_gridpoint_code, nrecs) = 1
        list_s(st_weight_code, nrecs) = 0

        list_s(st_input_bottom ,nrecs+1) =                                    &
          list_s(st_output_bottom, nrecs)
        list_s(st_input_bottom ,nrecs+2) =                                    &
          list_s(st_output_bottom, nrecs+1)

        list_s(st_input_top    ,nrecs+1) =                                    &
          list_s(st_output_top   , nrecs)
        list_s(st_input_top    ,nrecs+2) =                                    &
          list_s(st_output_top   , nrecs+1)

        list_s(st_input_code   ,nrecs+1) = -nrecs
        list_s(st_input_code   ,nrecs+2) = -nrecs - 1
        list_s(st_output_code  ,nrecs  ) = 1
        list_s(st_output_code  ,nrecs+1) = 1  ! dump
        list_s(st_series_ptr   ,nrecs+2) = 0
        list_s(st_series_ptr   ,nrecs)   = 0
        list_s(nelemp+1        ,nrecs+1) = nrecs + 1
        list_s(nelemp+1        ,nrecs+2) = nrecs + 2

        !  definition 8 implies frequency of time mean over every timestep
        list_s(st_freq_code, nrecs) = 1
        ! Frequency
        list_s(st_freq_code, nrecs+1) =                                       &
                    totimp(isam_t(itim_l), unt2_t(itim_l), modl_l)
        IF (list_s(st_freq_code, nrecs)  ==  -999) THEN
          errorstatus = 107
          WRITE (cmessage,'(A)') TRIM(ADJUSTL(FieldStr)) //          newline//&
            'PRELIM:TOTIMP:Error in time period conversion'
          GO TO 9999
        END IF


        !   Correct frequency for radiation, periodic convection, leaf
        !   phenology and vegetation competition

        IF (itima == 2) THEN
          list_s(st_freq_code, nrecs) = a_lw_radstep_diag
        ELSE IF (itima == 3) THEN
          list_s(st_freq_code, nrecs) = a_sw_radstep_diag
        ELSE IF (itima == 13) THEN
          list_s(st_freq_code, nrecs) = a_conv_step
        ELSE IF (itima == 14) THEN
          list_s(st_freq_code, nrecs) = phenol_period_ts
        ELSE IF (itima == 15) THEN
          list_s(st_freq_code, nrecs) = triffid_period_ts
        ELSE IF (itima == 16) THEN
          list_s(st_freq_code,nrecs)=river_step_ts
        ELSE IF (itima == 17) THEN
          list_s(st_freq_code,nrecs)=interval
        ELSE IF (itima == 18) THEN
          list_s(st_freq_code, nrecs) = a_lw_radstep_prog
        ELSE IF (itima == 19) THEN
          list_s(st_freq_code, nrecs) = a_sw_radstep_prog
        END IF

        !   Period
        ! time mean over sampling period
        list_s(st_period_code, nrecs) =                                       &
                  totimp(isam_t(itim_l), unt2_t(itim_l), modl_l)
        IF (list_s(st_period_code, nrecs)  ==  -999) THEN
          errorstatus = 108
          WRITE (cmessage,'(A)') TRIM(ADJUSTL(FieldStr)) //          newline//&
            'PRELIM:TOTIMP:Error in time period conversion'
          GO TO 9999
        END IF

        ! period for timeseries recycle period
        list_s(st_period_code, nrecs+1) =                                     &
                  totimp(intv_t(itim_l), unt1_t(itim_l), modl_l)
        IF (list_s(st_period_code, nrecs+1)  ==  -999) THEN
          errorstatus = 109
          WRITE (cmessage,'(A)') TRIM(ADJUSTL(FieldStr)) //          newline//&
            'PRELIM:TOTIMP:Error in time period conversion'
          GO TO 9999
        END IF


        ! st_start_time for 2 record should be period for first record
        ! unless offset from start of run. Note value independent of logical
        !  OFFSET

        IF (list_s(st_period_code, nrecs) /= -1) THEN
          IF (loffset) THEN
            list_s(st_start_time_code, nrecs+1) =                             &
              list_s(st_start_time_code, nrecs+1) -                           &
              list_s(st_period_code, nrecs+1)     +                           &
              list_s(st_freq_code, nrecs+1) - 1
          ELSE
            list_s(st_start_time_code, nrecs+1) =                             &
              list_s(st_start_time_code, nrecs+1) -                           &
              list_s(st_period_code, nrecs+1)     +                           &
              list_s(st_freq_code, nrecs+1)
          END IF
        ELSE
          list_s(st_start_time_code, nrecs+1) = 1
        END IF


        !   Add both record
        list_s(st_macrotag, nrecs)=0
        nrecs = nrecs + 2

      END IF       ! Other than single time field

    END IF         ! Diag request not null - ITIM_L /= 0
    999 CONTINUE
  END DO         ! Loop over diagnostic requests

END IF         ! NDIAG >  0

! DEPENDS ON: pslcom
CALL pslcom(nrecs)    ! Compress out unused pseudo levels lists

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE prelim

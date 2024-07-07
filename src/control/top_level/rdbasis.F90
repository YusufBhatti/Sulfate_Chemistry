! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Read the basis file information
! Subroutine Interface:

SUBROUTINE rdbasis                                                      &
   (cmessage,ErrorStatus)

USE ereport_mod, ONLY: ereport
USE check_iostat_mod, ONLY: check_iostat
USE yomhook, ONLY: lhook, dr_hook
USE setup_namelist, ONLY: setup_nml_type
USE parkind1, ONLY: jprb, jpim
USE filenamelength_mod, ONLY: filenamelength
USE stextend_mod, ONLY: npos_ts, nrecs_ts
USE lbc_mod, ONLY: rimwidtha, nrim_timesa, rim_lookupsa, bound_lookupsa
USE missing_data_mod, ONLY:  rmdi, imdi
USE version_mod, ONLY: nrecdp, nproftp, nprofdp, nprofup, nprofpp,    &
                       ndiagpm, ntimep, NTimSerP, nlevp, npslevp,     &
                       npslistp, tsdp
USE cstash_mod, ONLY: pslist_d, npslists, imn_d, imsk_d, iwst_d,      &
                      iwt_d, plpos_d, pllen_d, ts_d, plt_d, levt_d,   &
                      levb_d, iopa_d, iest_d, isth_d, inth_d, ilev_d, &
                      tlimr_ts, blimr_ts, tlim_ts, nlim_ts, wlim_ts,  &
                      elim_ts, slim_ts, usepro, rlevlst_d, levlst_d,  &
                      locn_u, blim_ts, nseries, iunt_u, iuse_b,       &
                      idom_b, itim_b, timpro, unt1_t, intv_t, ityp_t, &
                      ndprof, ntprof, ndiag, nuprof, item_b, isec_b,  &
                      modl_b, isam_t, itim_t, unt3_t, ioff_t, modl_t, &
                      iopl_d, dompro, iser_t, istr_t, iopt_t, unt2_t, &
                      iend_t, ifre_t, iedt_t, isdt_t, lts0_t,         &
                      pckpro, npprof, lnetcdf_u
USE ukca_option_mod,   ONLY: tc_lbc_ukca
USE ukca_tracer_stash, ONLY: a_max_ukcavars
USE umPrintMgr, ONLY: umPrint, umMessage, newline
USE free_tracers_inputs_mod, ONLY: a_max_trvars, i_free_tracer_lbc
USE um_parcore, ONLY: mype
USE mpl, ONLY: mpl_integer, mpl_real, mpl_address_kind, mpl_character
USE rimtypes, ONLY: nrima_max
USE nlsizes_namelist_mod, ONLY: &
    intf_lookupsa, tr_lbc_ukca, tr_lbc_vars, tr_ukca, tr_vars
USE file_manager, ONLY: assign_file_unit, release_file_unit, &
    get_file_unit_by_id
USE get_env_var_mod, ONLY: get_env_var
USE errormessagelength_mod, ONLY: errormessagelength
USE stash_model_mod, ONLY: nserrec_s, nserblk_s

USE model_domain_mod, ONLY: model_type, mt_global, mt_lam, mt_cyclic_lam       &
                          , mt_bi_cyclic_lam

USE profilename_length_mod, ONLY: profilename_length, fileid_length,    &
                                  packagename_length

USE stparam_mod, ONLY:                                                         &
  st_domain_whole_degrees, st_domain_gridpoints, st_levels_single

USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: IOSTAT_END

IMPLICIT NONE

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Code description:
! FORTRAN 90
! Written to UM programming standards UMDP3  vn 8.3

! Subroutine arguments

! Scalar arguments with intent(out):
CHARACTER(LEN=errormessagelength) :: cmessage      ! Error return message

! Local variables:
LOGICAL :: model_lev  !TRUE for model levels
LOGICAL :: lswap ! Indicates unfinished sort in the bubble sort block.
LOGICAL :: lswap_diag  ! Current element requires a swap in the bubble sort.
INTEGER :: iu       ! Unit no. of stash basis file
INTEGER :: i,j,k,l
INTEGER :: idum, imod
INTEGER :: IOSTAT
INTEGER :: ntsrecs    !Counter for time series records
INTEGER :: list_req_t1(3) ! Temporary list for request's STASHMaster data.
INTEGER :: list_req_t2(3) ! Temporary list for request's STASHMaster data.
INTEGER :: list_n_pro_t1(3) ! Temporary list for request's profile indices.
INTEGER :: list_n_pro_t2(3) ! Temporary list for request's profile indices.

! Temporary list for request's profiles.
CHARACTER(LEN=profilename_length) :: list_pro_t1(3) 
CHARACTER(LEN=profilename_length) :: list_pro_t2(3) 

CHARACTER (LEN=filenamelength) :: filename = "dummy file"

! Error status:
INTEGER :: ErrorStatus ! Error return code
INTEGER :: Readstatus  ! read  return code
CHARACTER(LEN=errormessagelength) :: iomessage

! Namelist UMSTASH_STREQ: STASH requests
INTEGER     :: isec, item       ! stash section and item
CHARACTER(LEN=profilename_length) :: dom_name         ! domain profile name
CHARACTER(LEN=profilename_length) :: use_name         ! usage profile name
CHARACTER(LEN=profilename_length) :: tim_name         ! time profile name
CHARACTER(LEN=packagename_length) :: package          ! name of package

NAMELIST/umstash_streq/ isec,item,dom_name,tim_name,use_name,package

! Namelist UMSTASH_TIME: Time profiles
INTEGER      :: unt1,unt2,unt3
INTEGER     :: ityp,isam,intv,iopt,itim
INTEGER     :: istr,iend,ifre,ioff,itimes,iser(ntimep)
INTEGER     :: isdt(6), iedt(6)
LOGICAL     :: lts0                 ! Disable/Enable pp output at TimeStep 0

NAMELIST/umstash_time/ityp,isam,intv,unt1  ,unt2,unt3,iopt            &
                     ,istr,iend,ifre,ioff,itimes,iser,tim_name        &
                     ,isdt,iedt,lts0

! Namelist UMSTASH_DOMAIN: Domain profiles
INTEGER     ::  iopl                !Level type code
INTEGER     ::  ilevs               !Flag for range/selected model levels
INTEGER     ::  levb,levt           !Bottom/top levels for range
INTEGER     ::  plt                 !Pseudo level type code
INTEGER     ::  iopa                !Horizontal domain type code
INTEGER     ::  inth,isth,iest,iwst !Horiz domain limits (IOPA=9,10)
INTEGER     ::  imsk                !Grid type code
INTEGER     ::  imn                 !Spatial meaning code
INTEGER     ::  iwt                 !Weighting code
INTEGER     ::  levlst (nlevp)      !Levels lists array: integer
REAL        ::  rlevlst(nlevp)      !real
INTEGER     ::  pslist (npslevp)    !pseudo
LOGICAL     ::  ts                        !Flag for time series domain
INTEGER     ::  tsnum                     !No. of time ser doms in prof
INTEGER     ::  tblim (tsdp),ttlim (tsdp) !TS dom limits (top/bottom)
REAL        ::  tblimr(tsdp),ttlimr(tsdp) !ditto for real levels
INTEGER     ::  tnlim (tsdp),tslim (tsdp) !TS dom limits (N/S)
INTEGER     ::  twlim (tsdp),telim (tsdp) !TS dom limits (E/W)
! Single Point Multiple Level (SPML) Time Series
LOGICAL     ::  l_spml_ts = .FALSE.
INTEGER     ::  spml_ns
INTEGER     ::  spml_ew
INTEGER     ::  spml_bot
INTEGER     ::  spml_top

NAMELIST/umstash_domain/iopl ,ilevs ,levb  ,levt  ,plt    ,iopa ,imsk ,       &
                              imn  ,iwt   ,ts    ,levlst,rlevlst,dom_name ,   &
                              inth ,isth  ,iest  ,iwst  ,pslist ,             &
                              tsnum,tblim ,ttlim ,tnlim ,tslim  ,telim,twlim, &
                              tblimr,ttlimr, l_spml_ts, spml_ns, spml_ew,     &
                              spml_bot, spml_top

! Namelist UMSTASH_USE: Usage profiles
INTEGER :: locn
INTEGER :: macrotag
CHARACTER(LEN=fileid_length) :: file_id
NAMELIST/umstash_use/locn,file_id,use_name,macrotag

! Namelist exclude_package: a package of stash requests to ignore 
NAMELIST/exclude_package/package

! variables for mpi communication of namelists
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

INTEGER :: information(10+nrima_max)

! set number of each type of variable in my_namelist type

INTEGER, PARAMETER :: no_of_types = 4
INTEGER, PARAMETER :: n_int = 7 + 6*ndiagpm + (25+ntimep)*nproftp   &
 + 6*ntimserp + (17+nlevp)*nprofdp  + npslevp*npslistp + 2*nprofup 
INTEGER, PARAMETER :: n_real = 2*ntimserp + nlevp*nprofdp
INTEGER, PARAMETER :: n_log = nprofdp + nproftp + nprofup
INTEGER, PARAMETER :: n_chars = profilename_length*nproftp +  &
                                profilename_length*nprofdp  +  &
                                profilename_length*nprofup +  &
                                packagename_length*nprofpp

! Function and subroutine calls:
LOGICAL :: disct_lev
LOGICAL :: diag_exclude

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RDBASIS'

!- End of Header ------------------------------------------------------

! setup namelist type - variables of the same type in the namelist
! should be together
! each variable in the namelist should be declared in the TYPE statement

TYPE my_namelist
  SEQUENCE
  INTEGER :: ndiag
  INTEGER :: ntprof
  INTEGER :: nseries
  INTEGER :: ndprof
  INTEGER :: nuprof
  INTEGER :: npprof
  INTEGER :: modl_b(ndiagpm)
  INTEGER :: isec_b(ndiagpm)
  INTEGER :: item_b(ndiagpm)
  INTEGER :: itim_b(ndiagpm)
  INTEGER :: idom_b(ndiagpm)
  INTEGER :: iuse_b(ndiagpm)
  INTEGER :: unt1_t(nproftp)
  INTEGER :: unt2_t(nproftp)
  INTEGER :: unt3_t(nproftp)
  INTEGER :: ityp_t(nproftp)
  INTEGER :: intv_t(nproftp)
  INTEGER :: isam_t(nproftp)
  INTEGER :: iopt_t(nproftp)
  INTEGER :: istr_t(nproftp)
  INTEGER :: iend_t(nproftp)
  INTEGER :: isdt_t(6,nproftp)
  INTEGER :: iedt_t(6,nproftp)
  INTEGER :: ifre_t(nproftp)
  INTEGER :: ioff_t(nproftp)
  INTEGER :: itim_t(nproftp)
  INTEGER :: iser_t(ntimep,nproftp)
  INTEGER :: modl_t(nproftp)
  INTEGER :: iopl_d(nprofdp)
  INTEGER :: levb_d(nprofdp)
  INTEGER :: levt_d(nprofdp)
  INTEGER :: iopa_d(nprofdp)
  INTEGER :: inth_d(nprofdp)
  INTEGER :: isth_d(nprofdp)
  INTEGER :: iest_d(nprofdp)
  INTEGER :: iwst_d(nprofdp)
  INTEGER :: imsk_d(nprofdp)
  INTEGER :: imn_d(nprofdp)
  INTEGER :: iwt_d(nprofdp)
  INTEGER :: blim_ts(ntimserp)
  INTEGER :: tlim_ts(ntimserp)
  INTEGER :: nlim_ts(ntimserp)
  INTEGER :: slim_ts(ntimserp)
  INTEGER :: elim_ts(ntimserp)
  INTEGER :: wlim_ts(ntimserp)
  INTEGER :: ilev_d(nprofdp)
  INTEGER :: levlst_d(nlevp,nprofdp )
  INTEGER :: plt_d(nprofdp)
  INTEGER :: pllen_d(nprofdp)
  INTEGER :: plpos_d(nprofdp)
  INTEGER :: pslist_d(npslevp,npslistp)
  INTEGER :: npslists
  INTEGER :: locn_u(nprofup)
  INTEGER :: iunt_u(nprofup)
  INTEGER :: npos_ts(nprofdp)
  INTEGER :: nrecs_ts(nprofdp)
  REAL :: blimr_ts(ntimserp)
  REAL :: tlimr_ts(ntimserp)
  REAL :: rlevlst_d(nlevp,nprofdp)
  LOGICAL :: lts0_t(nproftp)
  LOGICAL :: ts_d(nprofdp)
  LOGICAL :: lnetcdf_u(nprofup)
  CHARACTER (LEN=profilename_length) :: timpro(nproftp)
  CHARACTER (LEN=profilename_length) :: dompro(nprofdp)
  CHARACTER (LEN=profilename_length) :: usepro(nprofup)
  CHARACTER (LEN=packagename_length) :: pckpro(nprofpp)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int, &
                 n_real_in=n_real, n_log_in=n_log, n_chars_in=n_chars)

! Initialisation

information(:)=0

! All the reading in of information from namelists is done on task 0.
! Then the derived arrays and any other variables are broadcast to the
! other task. If introducing a new variable make sure it is broadcast.
! N.B. we are broadcasting the derived arrays not the namelists.

IF (mype==0) THEN

  ndiag   =0
  ntprof  =0
  ndprof  =0
  nuprof  =0
  npprof  =0
  ntsrecs =0

  DO i = 1,ndiagpm
    modl_b(i)=0
    isec_b(i)=0
    item_b(i)=0
    itim_b(i)=0
    idom_b(i)=0
    iuse_b(i)=0
  END DO

  itimes   =0
  DO i=1,nproftp
    timpro(i)=''
    ityp_t(i)= imdi
    intv_t(i)=imdi
    unt1_t(i)=imdi
    isam_t(i)=imdi
    unt2_t(i)=imdi
    iopt_t(i)=imdi
    istr_t(i)=0
    iend_t(i)=0
    lts0_t(i)=.FALSE.
    DO j=1,6
      isdt_t(j, i)=0
      iedt_t(j, i)=0
    END DO
    ifre_t(i)=imdi
    ioff_t(i)=imdi
    unt3_t(i)=imdi
    itim_t(i)=imdi
    modl_t(i)=1          ! hard wired to set model to atmosphere
    DO j = 1,ntimep
      iser_t(j,i)=imdi
    END DO
  END DO

  DO i=1,nprofdp
    dompro  (i)=''
    iopl_d  (i)=0
    levb_d  (i)=0
    levt_d  (i)=0
    plt_d   (i)=0
    iopa_d  (i)=0
    inth_d  (i)=0
    isth_d  (i)=0
    iest_d  (i)=0
    iwst_d  (i)=0
    imsk_d  (i)=0
    imn_d   (i)=0
    iwt_d   (i)=0
    ts_d    (i)=.FALSE.
    pllen_d (i)=0
    plpos_d (i)=0
    ilev_d  (i)=0
  END DO
  DO i = 1,nlevp
    DO j = 1,nprofdp
      levlst_d (i,j)=0
      rlevlst_d(i,j)=0.0
    END DO
  END DO
  DO i = 1,npslevp
    DO j = 1,npslistp
      pslist_d(i,j)=0
    END DO
  END DO

  DO i = 1,nprofup
    usepro(i)=''
    locn_u(i)=imdi
    iunt_u(i)=imdi
    lnetcdf_u(i)=.FALSE.
  END DO

  DO i = 1,nprofpp
    pckpro(i)=''
  END DO

  ! Get name for stash control file
  CALL get_env_var("STASHC",filename)

  ! Open stash control file
  CALL assign_file_unit(filename, iu, handler="fortran")
  OPEN(iu,FILE=filename, ACTION='READ', IOSTAT=ErrorStatus, IOMSG=iomessage)
  IF (ErrorStatus >  0) THEN
    CALL umPrint('RDBASIS : Failed in OPEN of Stash Control File:'// newline//&
                 'IoMsg: '//TRIM(iomessage), src='rdbasis')
    information(1)=ErrorStatus
    GO TO 9999
  ELSE IF (ErrorStatus <  0) THEN
    CALL umPrint('RDBASIS : Warning message on OPEN of Stash Control File', &
       src='rdbasis')
    WRITE(umMessage,'(a,i5)')'IOSTAT= ',ErrorStatus
    CALL umPrint(umMessage,src='rdbasis')
    WRITE(umMessage,'(a,i5)')'IOMSG= ',TRIM(iomessage)
    CALL umPrint(umMessage,src='rdbasis')
  END IF

  ! ---------------------------------------------

  ReadStatus = 0

  ! read in the domain profiles
  DO

    ! Initialise
    dom_name =''
    iopl =imdi
    levb =imdi
    levt =imdi
    ilevs=imdi
    DO j = 1,nlevp
      levlst (j)= imdi
      rlevlst(j)= rmdi
    END DO
    DO j = 1,npslevp
      pslist(j)=imdi
    END DO
    DO j = 1,tsdp
      tblim (j)=imdi
      ttlim (j)=imdi
      tblimr(j)=rmdi
      ttlimr(j)=rmdi
      tnlim (j)=imdi
      tslim (j)=imdi
      telim (j)=imdi
      twlim (j)=imdi
    END DO

    READ (UNIT=iu, NML=umstash_domain, IOSTAT=ReadStatus, IOMSG=iomessage)
    ! only call check_iostat if ReadStatus is > 0
    ! we do not want EOF to report a warning.

    IF (Readstatus == IOSTAT_END) THEN
      EXIT
    ELSE IF (Readstatus == 0) THEN

      ndprof = ndprof + 1

      IF (ndprof >  nprofdp) THEN
        CALL umPrint( 'ERROR IN STASHC:',src='rdbasis')
        WRITE(umMessage,'(a,i5)')'NUMBER OF DOMAIN PROFILES EXCEEDS LIMIT OF ',&
              nprofdp
        CALL umPrint(umMessage,src='rdbasis')
        CALL umPrint( 'ARRAYS WILL BE OVERWRITTEN',src='rdbasis')
        ErrorStatus=1
        information(1)=ErrorStatus
        GO TO 9999
      END IF

      ! Check for errors in levels lists
      ! DEPENDS ON: disct_lev
      model_lev=disct_lev(iopl,ErrorStatus,cmessage)

      IF (model_lev) THEN
        ! Model levels
        IF (ilevs == 2) THEN
          IF ( levlst(1) == imdi ) THEN
            WRITE(umMessage,'(a,i5,a)')                               &
               'ERROR,RDBASIS: LEVELS LIST IN DOMAIN PROFILE '        &
               ,ndprof,' HAS NO ENTRIES'
            CALL umPrint(umMessage,src='rdbasis')
            cmessage='ERROR,RDBASIS: LEVELS LIST HAS NO ENTRIES'
            ErrorStatus=1
            information(1)=ErrorStatus
            GO TO 9999
          END IF
        END IF
      ELSE IF (iopl /= st_levels_single) THEN
        IF (rlevlst(1) == rmdi) THEN
          WRITE(umMessage,'(a,i5,a)')                                 &
             'ERROR,RDBASIS: LEVELS LIST IN DOMAIN PROFILE '          &
             ,ndprof,' HAS NO ENTRIES'
          CALL umPrint(umMessage,src='rdbasis')
          cmessage='ERROR,RDBASIS: LEVELS LIST HAS NO ENTRIES'
          ErrorStatus=1
          information(1)=ErrorStatus
          GO TO 9999
        END IF
      END IF

      ! Profile name
      dompro(ndprof)=dom_name
      ! Store level type code in IOPL_D
      iopl_d(ndprof)=iopl
      ! DEPENDS ON: disct_lev
      model_lev=disct_lev(iopl,ErrorStatus,cmessage)

      IF (model_lev) THEN
        ! Integer levels
        ilev_d(ndprof)=ilevs
        IF (ilevs == 1) THEN
          !  Range of model levels
          levb_d(ndprof)=levb
          levt_d(ndprof)=levt
        END IF
        IF (ilevs == 2) THEN
          ! List of selected model levels
          levb_d(ndprof)=-1
          DO j=1,nlevp
            levlst_d(j,ndprof) = levlst(j)
            IF (levlst(j) >  0) THEN
              ! Store no. of levels in LEVT_D(ndprof)
              levt_d(ndprof)=levt_d(ndprof)+1
            END IF
          END DO
        END IF
      ELSE IF (iopl /= st_levels_single) THEN
        ! Real levels
        levb_d(ndprof)=-1
        DO j=1,nlevp
          rlevlst_d(j,ndprof) = rlevlst(j)
          IF (rlevlst(j) /=  rmdi) THEN
            ! Store no. of levels in LEVT_D(NDPROF)
            levt_d(ndprof)=levt_d(ndprof)+1
          END IF
        END DO
      END IF
      ! Store pseudo level type code in PLT_D
      plt_d (ndprof)=plt
      IF (plt >  0) THEN
        ! Domain profile 'NDPROF' has pseudo levels list
        ! Count total no. of pseudo levels lists
        npslists = npslists + 1
        IF (npslists >  npslistp) THEN
          WRITE(umMessage,'(a,a,i3,a,i3)')                          &
             'MESSAGE FROM ROUTINE RDBASIS: ',                    &
             'no. of pseudo levels lists ',npslists,              &
             ' exceeds allowed limit of ',npslistp
          CALL umPrint(umMessage,src='rdbasis')
          cmessage='ERROR,RDBASIS: NUMBER OF PSEUDO LEVELS LISTS EXCEEDS LIMIT'
          ErrorStatus=10
          CALL ereport("rdbasis",errorstatus,cmessage)
        END IF
        ! Store list in column 'NPSLISTS' of PSLIST_D
        DO j=1,npslevp
          pslist_d(j,npslists) = pslist (j)
          ! PPLEN_D(NDPROF) stores length of ps.lev.list for domain 'ndprof'
          IF (pslist(j) >  0) THEN
            pllen_d (ndprof)        = pllen_d(ndprof) + 1
          END IF
        END DO
        ! PLPOS(NDPROF) stores the column no. in PSLIST_D for dom. prof. 'NDPROF'
        plpos_d(ndprof) = npslists
      END IF
      ! Store horizontal domain type in IOPA_D
      iopa_d(ndprof)=iopa
      IF (iopa == st_domain_whole_degrees .OR.                                 &
               iopa == st_domain_gridpoints) THEN
        ! Specified area
        inth_d(ndprof)=inth
        isth_d(ndprof)=isth
        iest_d(ndprof)=iest
        iwst_d(ndprof)=iwst
      END IF
      imsk_d(ndprof)=imsk ! Gridpoint option
      imn_d (ndprof)=imn  ! Meaning option
      iwt_d (ndprof)=iwt  ! Weighting option
      ts_d  (ndprof)=ts   ! Time series domain
      IF (ts_d(ndprof)) THEN
        ! This domain profile has a block of time series domains
        ! Store time series data for current profile in _TS arrays
        nseries    = nseries+1        ! Time series block number:
        npos_ts(ndprof) = nseries          !      used as a pointer
        IF (l_spml_ts) THEN
           ! SPML: there is a time series for each level
           tsnum = spml_top - spml_bot+1
        END IF
        IF (tsnum < 1) THEN
          WRITE(cmessage,'(A,I0,A)')                                 &
              'Number of time series requested (tsnum = ',           &
              tsnum,                                                 &
              ') must be greater or equal to unity when ts == .true.'
          errorstatus=9
          CALL ereport("rdbasis",errorstatus,cmessage)
        END IF
        nrecs_ts(nseries) = tsnum     ! No. of records in ts block
        nserrec_s  = nserrec_s+tsnum  ! Cumulative total ts records

        IF (nserrec_s <= ntimserp) THEN
          DO j = 1,tsnum
            IF (j <= tsdp) THEN
              ntsrecs = ntsrecs+1
! DEPENDS ON: disct_lev
              model_lev=disct_lev(iopl,ErrorStatus,cmessage)
              
              IF (model_lev .AND. l_spml_ts .AND. j>=spml_bot         &
                  .AND. j<=spml_top) THEN
                blim_ts (ntsrecs)= j
                tlim_ts (ntsrecs)= j
              END IF
                
              IF (model_lev .AND. .NOT. l_spml_ts) THEN
                blim_ts (ntsrecs)= tblim (j)
                tlim_ts (ntsrecs)= ttlim (j)
              ELSE IF (iopl /= st_levels_single .AND. .NOT. l_spml_ts) THEN
                blimr_ts(ntsrecs)= tblimr(j)
                tlimr_ts(ntsrecs)= ttlimr(j)
              END IF

              IF (l_spml_ts) THEN
                nlim_ts(ntsrecs) = spml_ns
                slim_ts(ntsrecs) = spml_ns
                elim_ts(ntsrecs) = spml_ew
                wlim_ts(ntsrecs) = spml_ew
              ELSE
                nlim_ts(ntsrecs) = tnlim(j)
                slim_ts(ntsrecs) = tslim(j)
                elim_ts(ntsrecs) = telim(j)
                wlim_ts(ntsrecs) = twlim(j)
              END IF
            ELSE
              WRITE(umMessage,'(a,a,i5,a,i5,a,a)')                    &
                 'MESSAGE FROM ROUTINE RDBASIS: ',                    &
                 'no. of time series in domain profile ',ndprof,      &
                 ' exceeds allowed limit of ',tsdp,' some will be',   &
                 ' ignored'
              CALL umPrint(umMessage,src='rdbasis')
            END IF
          END DO
        ELSE
          WRITE(umMessage,'(a,a,i5,a)')                               &
             'TIMSER: total no. of time series requested exceeds ',   &
             'allowed limit of ',ntimserp,'; some will be ignored.'
          CALL umPrint(umMessage,src='rdbasis')
        END IF
      ELSE
        npos_ts  (ndprof)=-1
      END IF
    ELSE
      CALL check_iostat(ReadStatus, "namelist umstash_domain", iomessage )
    END IF
  END DO

  nserblk_s = nseries

  ! Rewind stash control file ahead of next namelist searches.
  REWIND(iu)

  ReadStatus = 0

  ! Time profile namelists
  DO 
    tim_name=''
    isam=imdi
    intv=imdi
    iopt=imdi
    istr=imdi
    iend=imdi
    ifre=imdi
    ioff=imdi
    unt1=imdi
    unt2=imdi
    unt3=imdi
    DO j = 1,ntimep
      iser(j)=imdi
    END DO
    DO j = 1,6
      isdt(j) = imdi
      iedt(j) = imdi
    END DO
    lts0=.FALSE.

    READ (UNIT=iu, NML=umstash_time, IOSTAT=ReadStatus, IOMSG=iomessage)
    ! only call check_iostat if ReadStatus is > 0
    ! we do not want EOF to report a warning.

    IF (Readstatus == IOSTAT_END) THEN
      EXIT

    ELSE IF (readstatus == 0) THEN

      ntprof = ntprof + 1

      IF (ntprof >  nproftp) THEN
        CALL umPrint( 'ERROR IN STASHC:',src='rdbasis')
        WRITE(umMessage,'(a,i5)') 'NUMBER OF TIME PROFILES EXCEEDS LIMIT OF ',  &
               nproftp
        CALL umPrint(umMessage,src='rdbasis')
        CALL umPrint( 'ARRAYS WILL BE OVERWRITTEN',src='rdbasis')
        ErrorStatus=1
        information(1)=ErrorStatus
        GO TO 9999
      END IF

      timpro(ntprof)=tim_name
      ityp_t(ntprof)=ityp
      IF (ityp /= 1) THEN
        ! Diagnostic output is time-processed
        isam_t(ntprof)=isam  !Sampling frequency
        unt2_t(ntprof)=unt2
        intv_t(ntprof)=intv  !Processing interval
        unt1_t(ntprof)=unt1
      END IF
      ! Diag. output time option
      iopt_t(ntprof)=iopt
      unt3_t(ntprof)=unt3
      IF (iopt == 1) THEN
        ! Regular output time interval
        istr_t(ntprof)=istr
        iend_t(ntprof)=iend
        ifre_t(ntprof)=ifre
        ioff_t(ntprof)=ioff
      ELSE IF (iopt == 2) THEN
        ! Specified list of output times
        ! Length of times table
        itim_t(ntprof)=itimes
        ! Times table
        DO j = 1,itimes
          iser_t(j,ntprof)=iser(j)
        END DO
      ELSE IF (iopt == 3) THEN
        DO j = 1,6
          isdt_t(j,ntprof)=isdt(j)
          iedt_t(j,ntprof)=iedt(j)
        END DO
        ifre_t(ntprof)=ifre
      END IF
      IF (iopt == 1 .AND. istr_t(ntprof) == 0) THEN
        lts0_t(ntprof)=lts0
      END IF
    ELSE
       CALL check_iostat(ReadStatus, "namelist umstash_time", iomessage )
    END IF
  END DO


  ! Rewind stash control file
  REWIND(iu)

  ReadStatus = 0

  !Usage profile namelists
  DO 
    use_name=''
    locn=imdi
    macrotag=imdi
    file_id=""

    READ (UNIT=iu, NML=umstash_use, IOSTAT=ReadStatus, IOMSG=iomessage)
    ! only call check_iostat if ReadStatus is > 0
    ! we do not want EOF to report a warning.

    IF (Readstatus == IOSTAT_END) THEN
      EXIT

    ELSE IF ( readstatus == 0) THEN

      nuprof = nuprof + 1

      IF (nuprof >  nprofup) THEN
        CALL umPrint( 'ERROR IN STASHC:',src='rdbasis')
        WRITE(umMessage,'(a,i5)') 'NUMBER OF USAGE PROFILES EXCEEDS LIMIT OF ',&
              nprofup
        CALL umPrint(umMessage,src='rdbasis')
        CALL umPrint( 'ARRAYS WILL BE OVERWRITTEN',src='rdbasis')
        ErrorStatus=1
        information(1)=ErrorStatus
        GO TO 9999
      END IF

      usepro(nuprof)=use_name
      locn_u(nuprof)=locn

      IF (locn == 3) THEN
        ! If this request is destined for an output file, the file id should
        ! have been provided - convert this to the correct unit number
        iunt_u(nuprof) = get_file_unit_by_id(file_id, handler="portio", &
                                             ignore_missing=.TRUE.)
        IF (iunt_u(nuprof) == imdi) THEN
          ! Output file not handled by "portio", next try "netcdf" file handler
          iunt_u(nuprof) = get_file_unit_by_id(file_id, handler="netcdf", &
                                               ignore_missing=.TRUE.)
          IF (iunt_u(nuprof) == imdi) THEN
            ! Error file id not found
            WRITE(cmessage,'(A)') &
               "Unable to find file id ("//TRIM(file_id)// &
               ") in portio or netcdf file lists"
            errorstatus=14
            CALL ereport("rdbasis",errorstatus,cmessage)
          ELSE
            ! file id associated with netCDF format output file
            lnetcdf_u(nuprof) = .TRUE.
          END IF
        END IF
      ELSE
        ! Otherwise iunt_u will be storing the value of the macro tag
        iunt_u(nuprof) = macrotag
      END IF
    ELSE
      CALL check_iostat(ReadStatus, "namelist umstash_use", iomessage )
    END IF
  END DO

  ! Rewind stash control file
  REWIND(iu)

  
  ReadStatus = 0

  !package profile exclusion namelist
  DO 

    package=''

    READ (UNIT=iu, NML=exclude_package, IOSTAT=ReadStatus, IOMSG=iomessage)
    ! only call check_iostat if ReadStatus is > 0
    ! we do not want EOF to report a warning.

    IF (Readstatus == IOSTAT_END) THEN
      EXIT

    ELSE IF ( readstatus == 0) THEN

      npprof = npprof + 1

      IF (npprof >  nprofpp) THEN
        CALL umPrint( 'ERROR IN STASHC:',src='rdbasis')
        WRITE(umMessage,'(a,i5)')                         &
              'NUMBER OF EXCLUDED PACKAGES EXCEEDS LIMIT OF ', nprofpp
        CALL umPrint(umMessage,src='rdbasis')
        CALL umPrint( 'ARRAYS WILL BE OVERWRITTEN',src='rdbasis')
        ErrorStatus=1
        information(1)=ErrorStatus
        GO TO 9999
      END IF

      pckpro(npprof)=package
    ELSE
      CALL check_iostat(ReadStatus, "namelist exclude_package", iomessage )
    END IF
  END DO

  ! Rewind stash control file
  REWIND(iu)

  !--------------------------------------------------

  ReadStatus = 0

  imod=1   ! hard-wired to atmosphere as only sub model

  ! read in the stash requests
  DO 

    isec=0
    item=0
    package=''

    READ (UNIT=iu, NML=umstash_streq, IOSTAT=ReadStatus, IOMSG=iomessage)
    ! only call check_iostat if ReadStatus is > 0
    ! we do not want EOF to report a warning.

    IF (Readstatus == IOSTAT_END) THEN
      EXIT

    ELSE IF (readstatus ==0 ) THEN

      ndiag = ndiag + 1

      IF (ndiag  >  nrecdp ) THEN
        WRITE(cmessage,'(a,a,i5,a)') 'NUMBER OF DIAGNOSTIC REQUESTS ',  &
                                   'EXCEEDS LIMIT OF ',               &
                                   nrecdp ,' SOME HAVE BEEN IGNORED'
        errorstatus=-15
        CALL ereport("rdbasis",errorstatus,cmessage)
        EXIT
      END IF

      modl_b(ndiag)=imod
      isec_b(ndiag)=isec
      item_b(ndiag)=item
      

      diag_exclude = .FALSE.
      DO j=1,npprof
        IF (package == pckpro(j)) THEN
          diag_exclude = .TRUE.
          EXIT
        END IF
      END DO

      IF (diag_exclude) THEN
        ndiag = ndiag - 1
        WRITE(cmessage,'(a,a,a,i5,i5)') 'Package ', TRIM(package),    &
                 ' is set to be excluded, excluding stash request', isec, item
        errorstatus=-1
        CALL ereport("rdbasis",errorstatus,cmessage)
        CYCLE
      END IF


      DO j=1,ntprof
        IF (tim_name == timpro(j)) THEN
          itim_b(ndiag)= j
          EXIT
        END IF
      END DO

      ! check time profile has been set and is one of the time profiles
      ! already read in from the STASHC file.

      IF (itim_b(ndiag)== 0) THEN
        WRITE(cmessage,'(a,a)') 'UNRECOGNISED time profile request ',   &
                             tim_name
        errorstatus=16
        CALL ereport("rdbasis",errorstatus,cmessage)
      END IF

      DO j=1,ndprof
        IF (dom_name == dompro(j)) THEN
          idom_b(ndiag)= j
          EXIT
        END IF
      END DO

      ! check domain profile has been set and is one of the domain profiles
      ! already read in from the STASHC file.

      IF (idom_b(ndiag)== 0) THEN
        WRITE(cmessage,'(a,a)') 'UNRECOGNISED domain profile request ', &
                             dom_name
        errorstatus=17
        CALL ereport("rdbasis",errorstatus,cmessage)
      END IF


      DO j=1,nuprof
        IF (use_name == usepro(j)) THEN
          iuse_b(ndiag)= j
          EXIT
        END IF
      END DO

      ! check usage profile has been set and is one of the usage profiles
      ! already read in from the STASHC file.

      IF (iuse_b(ndiag)== 0) THEN
        WRITE(cmessage,'(a,a)') 'UNRECOGNISED usage profile request ',  &
                             use_name
        errorstatus=18
        CALL ereport("rdbasis",errorstatus,cmessage)
      END IF

    ELSE
      CALL check_iostat(ReadStatus, "namelist umstash_streq", iomessage )
    END IF

  END DO

  ! Rewind stash control file
  REWIND(iu)


  ! Now sort requests by modl_b, isec_b, item_b, and profile names.
  ! Use a bubble sort.

  lswap = .TRUE.

  DO WHILE (lswap)

    ! Loop over USTASH_STREQ records and reorder based on their data.
    ! Continue until no more swaps occur (lswap is .FALSE.).
    lswap = .FALSE.

    DO i = ndiag-1, 1, -1
      j = i + 1

      ! Construct temporary lists used for sorting and temporary storage.
      ! list_req_t1, list_req_t2 hold STASHMaster request data
      ! list_pro_t1, list_pro_t2 hold profile name data
      ! list_n_pro_t1, list_n_pro_t2 hold the profile index data

      list_req_t1 = (/ modl_b(i), isec_b(i), item_b(i) /)
      list_req_t2 = (/ modl_b(j), isec_b(j), item_b(j) /)
      list_n_pro_t1 = (/ idom_b(i), itim_b(i), iuse_b(i) /)
      list_n_pro_t2 = (/ idom_b(j), itim_b(j), iuse_b(j) /)
      list_pro_t1 = (/ dompro(idom_b(i)), &
                       timpro(itim_b(i)), &
                       usepro(iuse_b(i)) /)
      list_pro_t2 = (/ dompro(idom_b(j)), &
                       timpro(itim_b(j)), &
                       usepro(iuse_b(j)) /)

      lswap_diag = .FALSE.
      DO k = 1, 6
        IF (k < 4) THEN
          IF (list_req_t1(k) > list_req_t2(k)) THEN
            lswap_diag = .TRUE.
          ELSE IF (list_req_t1(k) < list_req_t2(k)) THEN
            EXIT
          END IF
        ELSE IF (k > 3) THEN
          IF (list_pro_t1(k-3) > list_pro_t2(k-3)) THEN
            lswap_diag = .TRUE.
          ELSE IF (list_pro_t1(k-3) < list_pro_t2(k-3)) THEN
            EXIT
          END IF
        END IF
      END DO ! End loop over sort properties

      IF (lswap_diag) THEN
        lswap=.TRUE.
        modl_b(i) = list_req_t2(1)
        modl_b(j) = list_req_t1(1)
        isec_b(i) = list_req_t2(2)
        isec_b(j) = list_req_t1(2)
        item_b(i) = list_req_t2(3)
        item_b(j) = list_req_t1(3)
        idom_b(i) = list_n_pro_t2(1)
        idom_b(j) = list_n_pro_t1(1)
        itim_b(i) = list_n_pro_t2(2)
        itim_b(j) = list_n_pro_t1(2)
        iuse_b(i) = list_n_pro_t2(3)
        iuse_b(j) = list_n_pro_t1(3)
      END IF
    END DO  ! End loop over modl_b, etc elements

  END DO ! End while loop

  ! Set the number of tracers being used by counting the non-zero elements
  ! of the arrays that list them
  tr_lbc_vars=0
  tr_ukca=0
  tr_lbc_ukca=0
  DO i=1,a_max_trvars
    IF (i_free_tracer_lbc(i) /= 0) tr_lbc_vars = tr_lbc_vars + 1
  END DO
  DO i=1,a_max_ukcavars
    IF (tc_lbc_ukca(i) /= 0) tr_lbc_ukca = tr_lbc_ukca + 1
  END DO

  ! Set rim parameters to 1 unless this is a limited-area run
  IF (model_type == mt_global) THEN
    rimwidtha=1
    nrim_timesa=1
  ELSE IF (model_type == mt_lam .OR. model_type == mt_cyclic_lam .OR. &
      model_type == mt_bi_cyclic_lam) THEN
    nrim_timesa=MAX(nrim_timesa,1)
  END IF

  ! Update RIM_LOOKUPSA with the number of tracer LBCs
  rim_lookupsa = rim_lookupsa + tr_lbc_vars + tr_lbc_ukca
  ! Calculate total number of headers in the LBC file
  bound_lookupsa=rim_lookupsa+(nrim_timesa-1)*(rim_lookupsa-1)
  ! Number of atmosphere model interface lookups
  intf_lookupsa = intf_lookupsa + tr_vars + tr_ukca

  CLOSE(UNIT=iu,STATUS='KEEP',IOSTAT=iostat)
  CALL release_file_unit(iu, handler="fortran")

  information(2)=tr_lbc_vars
  information(3)=tr_ukca
  information(4)=tr_lbc_ukca
  information(5)=nserblk_s
  information(6)=nserrec_s
  information(7)=nrim_timesa
  information(8)=rim_lookupsa
  information(9)=bound_lookupsa
  information(10)=intf_lookupsa
  information(11:nrima_max+10) =rimwidtha(1:nrima_max)

  !!! set up mpi type array
  my_nml % ndiag    = ndiag
  my_nml % ntprof   = ntprof
  my_nml % nseries  = nseries
  my_nml % ndprof   = ndprof
  my_nml % nuprof   = nuprof
  my_nml % npprof   = npprof
  my_nml % modl_b   = modl_b
  my_nml % isec_b   = isec_b
  my_nml % item_b   = item_b
  my_nml % itim_b   = itim_b
  my_nml % idom_b   = idom_b
  my_nml % iuse_b   = iuse_b
  my_nml % unt1_t   = unt1_t
  my_nml % unt2_t   = unt2_t
  my_nml % unt3_t   = unt3_t
  my_nml % ityp_t   = ityp_t
  my_nml % intv_t   = intv_t
  my_nml % isam_t   = isam_t
  my_nml % iopt_t   = iopt_t
  my_nml % istr_t   = istr_t
  my_nml % iend_t   = iend_t
  my_nml % isdt_t   = isdt_t
  my_nml % iedt_t   = iedt_t
  my_nml % ifre_t   = ifre_t
  my_nml % ioff_t   = ioff_t
  my_nml % itim_t   = itim_t
  my_nml % iser_t   = iser_t
  my_nml % modl_t   = modl_t
  my_nml % iopl_d   = iopl_d
  my_nml % levb_d   = levb_d
  my_nml % levt_d   = levt_d
  my_nml % iopa_d   = iopa_d
  my_nml % inth_d   = inth_d
  my_nml % isth_d   = isth_d
  my_nml % iest_d   = iest_d
  my_nml % iwst_d   = iwst_d
  my_nml % imsk_d   = imsk_d
  my_nml % imn_d    = imn_d
  my_nml % iwt_d    = iwt_d
  my_nml % blim_ts  = blim_ts
  my_nml % tlim_ts  = tlim_ts
  my_nml % nlim_ts  = nlim_ts
  my_nml % slim_ts  = slim_ts
  my_nml % elim_ts  = elim_ts
  my_nml % wlim_ts  = wlim_ts
  my_nml % ilev_d   = ilev_d
  my_nml % levlst_d = levlst_d
  my_nml % plt_d    = plt_d
  my_nml % pllen_d  = pllen_d
  my_nml % plpos_d  = plpos_d
  my_nml % pslist_d = pslist_d
  my_nml % npslists = npslists
  my_nml % locn_u   = locn_u
  my_nml % iunt_u   = iunt_u
  my_nml % npos_ts  = npos_ts
  my_nml % nrecs_ts = nrecs_ts
  ! end of integers
  my_nml % blimr_ts  = blimr_ts
  my_nml % tlimr_ts  = tlimr_ts
  my_nml % rlevlst_d = rlevlst_d
  my_nml % lts0_t    = lts0_t
  my_nml % ts_d      = ts_d
  my_nml % lnetcdf_u = lnetcdf_u
  my_nml % timpro    = timpro
  my_nml % dompro    = dompro
  my_nml % usepro    = usepro
  my_nml % pckpro    = pckpro

  9999 CONTINUE

END IF ! mype==0

!! End of mype == 0 section, now broadcast all information to other processors

CALL mpl_bcast(information,10+nrima_max,mpl_integer,0,my_comm,icode)

IF (mype/=0) THEN
  errorstatus            = information(1)
  tr_lbc_vars            = information(2)
  tr_ukca                = information(3)
  tr_lbc_ukca            = information(4)
  nserblk_s              = information(5)
  nserrec_s              = information(6)
  nrim_timesa            = information(7)
  rim_lookupsa           = information(8)
  bound_lookupsa         = information(9)
  intf_lookupsa          = information(10)
  rimwidtha(1:nrima_max) = information(11:nrima_max+10)
END IF

IF (errorstatus <=0) THEN

  CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

  !!!  if(mype/=0)then unpopulate mpi type
  IF (mype /= 0) THEN
    ndiag     = my_nml % ndiag
    ntprof    = my_nml % ntprof
    nseries   = my_nml % nseries
    ndprof    = my_nml % ndprof
    nuprof    = my_nml % nuprof
    npprof    = my_nml % npprof
    modl_b    = my_nml % modl_b
    isec_b    = my_nml % isec_b
    item_b    = my_nml % item_b
    itim_b    = my_nml % itim_b
    idom_b    = my_nml % idom_b
    iuse_b    = my_nml % iuse_b
    unt1_t    = my_nml % unt1_t
    unt2_t    = my_nml % unt2_t
    unt3_t    = my_nml % unt3_t
    ityp_t    = my_nml % ityp_t
    intv_t    = my_nml % intv_t
    isam_t    = my_nml % isam_t
    iopt_t    = my_nml % iopt_t
    istr_t    = my_nml % istr_t
    iend_t    = my_nml % iend_t
    isdt_t    = my_nml % isdt_t
    iedt_t    = my_nml % iedt_t
    ifre_t    = my_nml % ifre_t
    ioff_t    = my_nml % ioff_t
    itim_t    = my_nml % itim_t
    iser_t    = my_nml % iser_t
    modl_t    = my_nml % modl_t
    iopl_d    = my_nml % iopl_d
    levb_d    = my_nml % levb_d
    levt_d    = my_nml % levt_d
    iopa_d    = my_nml % iopa_d
    inth_d    = my_nml % inth_d
    isth_d    = my_nml % isth_d
    iest_d    = my_nml % iest_d
    iwst_d    = my_nml % iwst_d
    imsk_d    = my_nml % imsk_d
    imn_d     = my_nml % imn_d
    iwt_d     = my_nml % iwt_d
    blim_ts   = my_nml % blim_ts
    tlim_ts   = my_nml % tlim_ts
    nlim_ts   = my_nml % nlim_ts
    slim_ts   = my_nml % slim_ts
    elim_ts   = my_nml % elim_ts
    wlim_ts   = my_nml % wlim_ts
    ilev_d    = my_nml % ilev_d
    levlst_d  = my_nml % levlst_d
    plt_d     = my_nml % plt_d
    pllen_d   = my_nml % pllen_d
    plpos_d   = my_nml % plpos_d
    pslist_d  = my_nml % pslist_d
    npslists  = my_nml % npslists
    locn_u    = my_nml % locn_u
    iunt_u    = my_nml % iunt_u
    npos_ts   = my_nml % npos_ts
    nrecs_ts  = my_nml % nrecs_ts
    blimr_ts  = my_nml % blimr_ts
    tlimr_ts  = my_nml % tlimr_ts
    rlevlst_d = my_nml % rlevlst_d
    lts0_t    = my_nml % lts0_t
    ts_d      = my_nml % ts_d
    lnetcdf_u = my_nml % lnetcdf_u
    timpro    = my_nml % timpro
    dompro    = my_nml % dompro
    usepro    = my_nml % usepro
    pckpro    = my_nml % pckpro
  END IF

END IF   !! errorstatus<=0

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE rdbasis

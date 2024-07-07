! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!    Purpose : Read from ACOBS Files,reformat and place OBS header
!              details in COMOBS. The bulk of the required OBS data
!              is put into dynamic work array OBS for transmission via
!              argument list to GETOBS. OBS is written out to a cache
!              file for subsequent reading at later timesteps.
!              Thus reread of ACOBS files only required intermittently
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE rdobs2_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RDOBS2_MOD'

CONTAINS

SUBROUTINE rdobs2(nfiles,timestep_no,obs,obs_flag,tndv,           &
                  tnobs,timerel,                                  &
                  lambda_p,phi_p,p_rows,row_length,               &
                  icode,cmessage)
! ----------------------------------------------------------------------
!    INTENT IN:
!       NFILES   : NO OF AC OBSERVATION FILES TO BE READ
!       TIMESTEP_NO : TIMESTEP NUMBER
!       TIMEREL     : relative time for this timestep
!    INTENT OUT:
!       TNDV     : ACTUAL SIZE OF OBS ARRAY
!       TNOBS    : ACTUAL NO OF OBS
!       OBS      : OBS array
!       OBS_FLAG : do not use flags
!       ICODE/CMESSAGE: for error processing
!
!  ONLY APPLICABLE TO LAMS NOW
!  MOPS DATA ONLY - NO WIND OBS
! ----------------------------------------------------------------------
!
!
USE conversions_mod, ONLY: pi_over_180
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io
USE model_file, ONLY: model_file_open, model_file_close
USE UM_ParVars
USE UM_ParCore,   ONLY: mype, nproc
USE UM_ParParams, ONLY: peast, pwest, pnorth, psouth
USE comobs_mod, ONLY: nobtypmx, ndatavmx, obs_info, used_files
USE days_mod, ONLY: days
USE rdobs3_mod, ONLY: rdobs3
USE setdac_mod, ONLY: setdac
USE latlon_eq_rotation_mod, ONLY: rotate_latlon_to_eq
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE umPrintMgr, ONLY:      &
    str,                    &
    umPrint,                &
    umMessage
USE ac_control_mod
USE nlsizes_namelist_mod, ONLY: model_levels, len_fixhd
USE num_obs_mod, ONLY: ac_num_obs_max, ac_tot_obs_size_max
USE ac_dump_mod, ONLY: len_inthd, len_realhd, len1_levdepc, len2_levdepc, &
                       len1_rowdepc, len2_rowdepc, len1_coldepc,          &
                       len2_coldepc, len1_flddepc, len2_flddepc,          &
                       len_extcnst, len_dumphist, len_cfi1, len_cfi2,     &
                       len_cfi3, len1_lookup_obs, len2_lookup_obs, fixhd, &
                       allocate_ac_dump_headers,                          &
                       deallocate_ac_dump_headers
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: l_regular, model_type, mt_global

IMPLICIT NONE

INTEGER :: common_length

!-----------------------------------------------------------------------
!     ARGUMENTS
!-----------------------------------------------------------------------
INTEGER :: tndv,tnobs,max_ndv
INTEGER :: nfiles
INTEGER :: timestep_no
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER :: p_rows,row_length
REAL :: lambda_p(1-halo_i:row_length+halo_i)
REAL :: phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)

!-----------------------------------------------------------------------
!     LOCAL VARIABLES
!-----------------------------------------------------------------------
INTEGER :: ktype
INTEGER :: io_stat
INTEGER :: jf,jobt,jobt2,jdv,j,job,jtyp,jlev,jact,jmot
INTEGER :: itot0,itot1,itotal
INTEGER :: irday,ifday
INTEGER :: ipt,ipc,ifile,i,iworksp,iobt,ijf
INTEGER :: ipt_this_file
INTEGER :: max_ndatav
INTEGER :: indvmax
REAL ::                                                           &
 zzlatmx,zzlatmn,zzlongmx,zzlongmn,                               &
 timerel,timeadj,twstart,twend,tgetobb,tgetoba
INTEGER :: iproc,istat,imsg
REAL :: w_limit,e_limit,n_limit,s_limit,delta_phi,delta_lambda
!-----------------------------------------------------------------------
REAL :: datalevs (model_levels+1,nobtypmx)
INTEGER :: iref     (nobtypmx,ndatavmx,nfiles)
INTEGER :: inobs    (nobtypmx,nfiles)
INTEGER :: iobstyp  (nobtypmx,nfiles)
INTEGER :: indatav  (nobtypmx,nfiles)
INTEGER :: inoblev  (nobtypmx,nfiles)
INTEGER :: ioblvtyp (nobtypmx,nfiles)
INTEGER :: indvhdr  (nfiles)
INTEGER :: inobtyp  (nfiles)
INTEGER :: imaxnlev (nfiles)
REAL :: plevels     (model_levels+1,nobtypmx,nfiles)
INTEGER :: kobstyp  (nobtypmx)
LOGICAL :: lempty   (nfiles)
!-----------------------------------------------------------------------
INTEGER :: obs_file_yy,obs_file_mm,obs_file_dd
INTEGER :: obs_file_hh,obs_file_min,obs_file_sec
INTEGER :: ip_lat,ip_long,ip_time,ip_type,ip_mot,ip_u,ip_v
INTEGER :: ip_mot1,ip_mot2,ip_type1,ip_type2
INTEGER :: idummy
!-----------------------------------------------------------------------

INTEGER :: len_data
!-----------------------------------------------------------------------
!     DYNAMIC ALLOCATION
INTEGER :: obs_flag(ac_num_obs_max)
INTEGER :: iwork(ac_num_obs_max)
REAL :: obs(ac_tot_obs_size_max)
REAL :: work(ac_tot_obs_size_max+2048)
REAL :: u_wrk(ac_num_obs_max),v_wrk(ac_num_obs_max)
REAL :: wrklat(ac_num_obs_max),wrklon(ac_num_obs_max)
REAL :: coeff1(ac_num_obs_max),coeff2(ac_num_obs_max)
INTEGER :: ispare(7)
INTEGER :: envvar

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RDOBS2'


! Calculate limits of domain in degrees
IF (l_regular) THEN

  n_limit=lat_n
  s_limit=lat_s-0.5*dlat
  w_limit=long_w_model
  e_limit=long_e_model+0.5*dlong

ELSE

  delta_phi=(phi_p(1,2)-phi_p(1,1))/pi_over_180
  delta_lambda=(lambda_p(row_length)-lambda_p(row_length-1)) &
           /pi_over_180

  n_limit=phi_p(1,p_rows)/pi_over_180
  s_limit=phi_p(1,1)/pi_over_180-0.5*delta_phi
  w_limit=lambda_p(1)/pi_over_180
  e_limit=lambda_p(row_length)/pi_over_180+0.5*delta_lambda

END IF

!-----------------------------------------------------------------------
!            SECTION 1: COPY INPUT FILES TO WORK AREA, ETC.
!-----------------------------------------------------------------------

!       Read in the AC Observation files and merge the observations.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (diag_rdobs >= 1 .AND. mype == 0) THEN
  CALL umPrint(' ',src='rdobs2')
  WRITE(umMessage,*) 'READ IN AC OBS FILES - TIMESTEP : ',timestep_no
  CALL umPrint(umMessage,src='rdobs2')
END IF

!       Set up time for next read of observation files
obs_info % timenext = timerel + obs_info % timeint

!       Initilaise to zero to allow for empty files.
DO jf=1,nfiles
  lempty  (jf) = .FALSE.
  indvhdr (jf) = 0
  inobtyp (jf) = 0
  imaxnlev(jf) = 0
  DO jobt=1,nobtypmx
    iobstyp  (jobt,jf) = 0
    indatav  (jobt,jf) = 0
    inobs    (jobt,jf) = 0
    inoblev  (jobt,jf) = 0
    ioblvtyp (jobt,jf) = 0
    DO jlev=1,model_levels+1
      plevels(jlev,jobt,jf) = 0.0
    END DO
  END DO
END DO

!       Initialise IREF
DO jf = 1,nfiles
  DO jdv = 1,ndatavmx
    DO jobt = 1,nobtypmx
      iref(jobt,jdv,jf) = 0
    END DO
  END DO
END DO

!       IPT IS THE AMOUNT OF WORK SPACE USED IN THE ARRAY WORK SO FAR
ipt=0

!       JF LOOP: ONE CYCLE FOR EACH INPUT FILE.
DO jf=1,obs_info % num_used_files

  CALL assign_file_unit(used_files(jf), ifile, handler="portio")

  !       Get unit number of observation file

  envvar = 1
  CALL model_file_open(ifile,used_files(jf),&
                       obs_info % filename_len(jf),0,envvar,icode)


  !       READ FILE CONTENTS TO BUFFER.
  !       IF FILE IS EMPTY - PROCEED TO READ NEXT FILE

  obs_info % nobtyp = 0
  tndv   = 0
  tnobs  = 0
  obs_info % ndvhdr = 0
  DO jobt=1,nobtypmx
    obs_info % obstyp  (jobt) = 0
    obs_info % ndatav  (jobt) = 0
    obs_info % noblev  (jobt) = 0
    obs_info % oblevtyp(jobt) = 0
    obs_info % nobs    (jobt) = 0
  END DO

  !-----------------------------------------------------------------------
  IF (obs_format == 2 .OR. obs_format == 3) THEN

    !         Go to start of obs file
    ! read in fixed length header

    ! DEPENDS ON: read_flh
    CALL read_flh (ifile,fixhd,len_fixhd,icode,cmessage)
    IF (icode >  0) GO TO 9999

    !         Get dimensions of data set components in Fixed length header
    ! DEPENDS ON: get_dim
    CALL get_dim (fixhd,                                                    &
                  len_fixhd, len_inthd, len_realhd,                         &
                  len1_levdepc, len2_levdepc, len1_rowdepc, len2_rowdepc,   &
                  len1_coldepc, len2_coldepc, len1_flddepc,len2_flddepc,    &
                  len_extcnst,  len_dumphist,                               &
                  len_cfi1, len_cfi2, len_cfi3,                             &
                  len1_lookup_obs, len2_lookup_obs,                         &
                  len_data)

    CALL allocate_ac_dump_headers()

    obs_file_yy  = fixhd(21)   ! )
    obs_file_mm  = fixhd(22)   ! ) Date/time
    obs_file_dd  = fixhd(23)   ! ) of
    obs_file_hh  = fixhd(24)   ! ) Observation
    obs_file_min = fixhd(25)   ! ) File
    obs_file_sec = fixhd(26)   ! )

    max_ndv      = fixhd(162)  !  Max total no of data values

    !         If MAX_NDV = 0, set to 1. (Obs file with no observations)
    !         Prevents dynamic allocation with 0 in READACOBS
    IF (max_ndv == 0) max_ndv = 1

    !
    !  check dimension of WORK
    !
    common_length = 2048+ipt+ obs_info % per_file_tndvmax(jf)
    ! 2048 allows space for headers when no or very few obs; was 3000.
    IF (common_length  >   ac_tot_obs_size_max+2048) THEN
      WRITE(umMessage,*)'RDOBS2 - common_length,IPT,PER_FILE_TNDVMAX,' &
          ,'AC_TOT_OBS_SIZE_MAX+2048 ',common_length,ipt,              &
          obs_info % per_file_tndvmax,  ac_tot_obs_size_max+2048
      CALL umPrint(umMessage,src='rdobs2')
      !
      ! important failure messages to PEn,OUT and OUT2 output streams
      WRITE(umMessage,*)' RDOBS2 : Insufficient space in WORK array'
      CALL umPrint(umMessage,src='rdobs2')
      WRITE(umMessage,*) ' Recode or rerun with fewer obs'
      CALL umPrint(umMessage,src='rdobs2')
      WRITE(umMessage,*)'AC_TOT_OBS_SIZE_MAX',ac_tot_obs_size_max,     &
                        ' should be >= ', common_length-2048
      CALL umPrint(umMessage,src='rdobs2')
      icode = 1
      cmessage =' RDOBS2 : Insufficient space in WORK array'
      GO TO 9999
    ELSE
      IF (diag_rdobs >= 1) THEN
        ! diagnostic message directed to every PEn output stream
        ipc=(common_length*100)/(ac_tot_obs_size_max+2048)
        WRITE(umMessage,*) &
            'RDOBS2:% of WORK used so far for reading obs=',ipc
        CALL umPrint(umMessage,src='rdobs2')
      END IF
    END IF
    CALL rdobs3(ifile,nobtypmx,obs_info % nobtyp,obs_info % obstyp,  &
                obs_info % nobs, obs_info % ndatav,                  &
                obs_info % noblev,obs_info % oblevtyp,datalevs,      &
                work(ipt+1),                                         &
                len_data,max_ndv,tnobs,obs_info % missd,             &
                icode,cmessage                                       &
                         ,ipt                                        &
                             )

    IF (icode >  0) GO TO 9999

    obs_info % ndvhdr = 5                 ! No of data values in header
    obs_info % maxnlev1 = len1_levdepc-2  ! Max no of levels + 1

  ELSE
    icode=1
    cmessage='RDOBS2: ILLEGAL OBS_FORMAT'
    GO TO 9999
  END IF

  !       Convert any old obs type numbers to new obs type numbers.
  !       For version 2.6 onwards.
  IF (obs_info % nobtyp >  0) THEN
    DO jobt=1,obs_info % nobtyp
      IF (obs_info % obstyp(jobt) == 501) THEN
        obs_info % obstyp(jobt) = 302
        CALL umPrint('Type 501 in Obs file changed to 302', &
            src='rdobs2',pe=0)
      END IF
      IF (obs_info % obstyp(jobt) == 502) THEN
        obs_info % obstyp(jobt) = 305
        CALL umPrint('Type 502 in Obs file changed to 305', &
            src='rdobs2',pe=0)
      END IF
    END DO
  END IF

  !       PRINT CONTENTS OF FILE HEADER.
  IF (diag_rdobs >= 1 .AND. mype == 0) THEN
    CALL umPrint(' ',src='rdobs2')
    CALL umPrint(' AC OBS FILE - UNIT NO :'// &
        TRIM(str(ifile)),src='rdobs2')
    CALL umPrint(' NO OF OBS TYPES       :'// &
        TRIM(str(obs_info % nobtyp)),src='rdobs2')
    CALL umPrint(' TOTAL NO OF OBS       :'// &
        TRIM(str(tnobs)),src='rdobs2')
    CALL umPrint(' ',src='rdobs2')
    WRITE(umMessage,'(A,I3.2,A,I2.2,A,I4)')                              &
        ' DATE :',obs_file_dd,'/',obs_file_mm,'/',obs_file_yy
    CALL umPrint(umMessage,src='rdobs2')
    WRITE(umMessage,'(A,I3.2,I2.2,A)')                                   &
        ' TIME :',obs_file_hh,obs_file_min,'Z'
    CALL umPrint(umMessage,src='rdobs2')
    CALL umPrint(' ',src='rdobs2')

    IF (obs_info % nobtyp >  0) THEN

      WRITE(umMessage,'(A,T30,18I7/(T30,18I7))')                         &
          ' AC OBS TYPES      :',(obs_info % obstyp(i),          &
          i=1,obs_info % nobtyp)
      CALL umPrint(umMessage,src='rdobs2')
      WRITE(umMessage,'(A,T30,18I7/(T30,18I7))')                         &
          ' NO OF LEVELS      :',(obs_info % noblev(i),         &
          i=1,obs_info % nobtyp)
      CALL umPrint(umMessage,src='rdobs2')
      WRITE(umMessage,'(A,T30,18I7/(T30,18I7))')                         &
          ' NO OF DATA VALUES :',(obs_info % ndatav(i),         &
          i=1,obs_info % nobtyp)
      CALL umPrint(umMessage,src='rdobs2')
      WRITE(umMessage,'(A,T30,18I7/(T30,18I7))')                         &
          ' OBS LEVEL TYPE    :',(obs_info % oblevtyp(i),       &
          i=1,obs_info % nobtyp)
      CALL umPrint(umMessage,src='rdobs2')
      WRITE(umMessage,'(A,T30,18I7/(T30,18I7))')                         &
          ' NO OF OBS         :',                                &
          (obs_info % nobs(i),i=1,obs_info % nobtyp)
      CALL umPrint(umMessage,src='rdobs2')

    END IF

  END IF

  IF (timestep_no == 1 .AND. tnobs == 0) THEN
    IF (obs_info % obstyp(1) == 406) THEN
      WRITE (umMessage,'(A,I3.2,I2.2,A,I3.2,A,I2.2,A,I4,A)')             &
          ' AC Observation File (MOPS) - ',                         &
          obs_file_hh,obs_file_min,'Z',                             &
          obs_file_dd,'/',obs_file_mm,'/',obs_file_yy,              &
          ' - No MOPS ACOBS data'
      CALL umPrint(umMessage,src='rdobs2')
    ELSE
      WRITE (umMessage,'(A,I3.2,I2.2,A,I3.2,A,I2.2,A,I4,A)')             &
          ' AC Observation File - ',                                &
          obs_file_hh,obs_file_min,'Z',                             &
          obs_file_dd,'/',obs_file_mm,'/',obs_file_yy,              &
          ' - No standard ACOBS data'
      CALL umPrint(umMessage,src='rdobs2')
    END IF
  END IF

  !       Check that no of observation types in file (NOBTYP) does not
  !       exceed maximum allowed. (NOBTYPMX)
  IF (obs_info % nobtyp >  nobtypmx) THEN
    icode = 1
    cmessage = ' RDOBS2 : TOO MANY OBSERVATION TYPES IN OBS FILE'
    IF (mype == 0) THEN
      CALL umPrint(' RDOBS2 : TOO MANY OBSERVATION TYPES IN OBS FILE', &
          src='rdobs2')
      CALL umPrint(' NO OF OBS TYPES = '//TRIM(str(obs_info % nobtyp)), &
          src='rdobs2')
      CALL umPrint(' MAXIMUM ALLOWED = '//TRIM(str(nobtypmx)), &
          src='rdobs2')
    END IF
    GO TO 9999
  END IF

  !       Check that no of data values for each obs type (NDATAV) does
  !       not exceed maximum allowed. (NDATAVMX)
  IF (obs_info % nobtyp >  0) THEN
    DO jobt=1,obs_info % nobtyp
      IF (obs_info % ndatav(jobt) >  ndatavmx) THEN
        icode = 1
        cmessage = ' RDOBS2 : TOO MANY DATA VALUES IN OBS FILE'
        IF (mype == 0) THEN
          WRITE(umMessage,*) ' RDOBS2 : Too many Data values in Obs file'
          CALL umPrint(umMessage,src='rdobs2')
          WRITE(umMessage,*) ' Observation Types =  ',obs_info % obstyp(jobt)
          CALL umPrint(umMessage,src='rdobs2')
          WRITE(umMessage,*) ' No of Data Values =  ',obs_info % ndatav(jobt)
          CALL umPrint(umMessage,src='rdobs2')
          WRITE(umMessage,*) ' Maximum allowed   =  ',ndatavmx
          CALL umPrint(umMessage,src='rdobs2')
        END IF
        GO TO 9999
      END IF
    END DO
  END IF

  !       Store information for this obs file
  inobtyp(jf) = obs_info % nobtyp
  indvhdr(jf) = obs_info % ndvhdr
  imaxnlev(jf) = obs_info % maxnlev1
  DO jobt=1,obs_info % nobtyp
    iobstyp(jobt,jf)  = obs_info % obstyp(jobt)
    indatav(jobt,jf)  = obs_info % ndatav(jobt)
    inoblev(jobt,jf)  = obs_info % noblev(jobt)
    ioblvtyp(jobt,jf) = obs_info % oblevtyp(jobt)
    DO jlev=1, obs_info % maxnlev1
      plevels(jlev,jobt,jf) = datalevs(jlev,jobt)
    END DO
  END DO

  !       EVALUATE POINTERS TO EACH SUB-VECTOR.
  !       IPT=0 BEFORE READING FIRST FILE

  ipt_this_file = ipt

  DO jobt=1,obs_info % nobtyp
    DO jdv=1,obs_info % ndatav(jobt)
      iref(jobt,jdv,jf)=ipt
      ipt = ipt+ obs_info % nobs(jobt)
    END DO
  END DO

  !       CONVERT REFERENCE DATE & FILE DATE TO ELAPSED DAYS

  CALL days (obs_info % obs_ref_dd, obs_info % obs_ref_mm, &
             obs_info % obs_ref_yy, irday)
  CALL days (obs_file_dd,obs_file_mm,obs_file_yy,ifday)

  !       FIND RELATIVE TIME ADJUSTMENT.

  timeadj = 1440*(ifday-irday) +                                  &
   60*(obs_file_hh-obs_info % obs_ref_hh) + obs_file_min - obs_info % obs_ref_min

  IF (obs_info % nobtyp >  0) THEN

    !         LOOP OVER OBSERVATION TYPES

    DO jobt=1,obs_info % nobtyp

      tgetobb = 0.0
      tgetoba = 0.0
      DO jobt2=1,nobtypmx
        IF (obs_info % obstyp(jobt) == master_ac_types(jobt2)) THEN
          tgetobb = def_tgetobb(jobt2)
          tgetoba = def_tgetoba(jobt2)
        END IF
      END DO
      IF (tgetobb == 0.0 .OR. tgetoba == 0.0) THEN
        icode = 1
        cmessage = 'RDOBS2 : TGETOBB and/or TGETOBA = 0.0 ?'
        GO TO 9999
      END IF

      !         Set up time window to get observations required for assim.
      twstart  = timerel - tgetoba
      twend    = timerel + obs_info % timeint + tgetobb

      zzlatmx = n_limit
      zzlatmn = s_limit
      IF (at_extremity(PNorth))zzlatmx= obs_info % obs_lat_n
      IF (at_extremity(PSouth))zzlatmn= obs_info % obs_lat_s
      zzlongmx= e_limit
      zzlongmn= w_limit
      IF (at_extremity(PEast))zzlongmx= obs_info % obs_long_e
      IF (at_extremity(PWest))zzlongmn= obs_info % obs_long_w
      IF (zzlongmn >= 360.0)zzlongmn=zzlongmn-360.0
      IF (zzlongmx >= 360.0)zzlongmx=zzlongmx-360.0

      ip_lat  = 1
      ip_long = 2
      ip_time = 3

      IF (obs_info % nobs(jobt) >  0) THEN
        !
        !  Ensure boundaries are numerically the same.

        ! process PEs not along eastern boundary
        IF (.NOT. at_extremity(peast)) THEN

          DO iproc=1,nproc-1

            imsg=jobt*100+iproc  ! message tag

            IF (mype == iproc) THEN
              IF (.NOT. at_extremity(pwest)) THEN
                ! this PE receives if not along western boundary
                CALL gc_rrecv(imsg,1,iproc-1,istat,zzlongmn,zzlongmx)
              END IF
            ELSE IF (mype == iproc-1) THEN
              ! PEs not along eastern boundary send to adjacent PE
              CALL gc_rsend(imsg,1,iproc,istat,zzlongmn,zzlongmx)
            END IF

          END DO

        ELSE IF (.NOT. at_extremity(pwest)) THEN
          ! PEs along eastern boundary can receive only

          DO iproc=1,nproc-1
            IF (mype == iproc) THEN        ! this PE receives
              imsg=jobt*100+iproc
              CALL gc_rrecv(imsg,1,iproc-1,istat,zzlongmn,zzlongmx)
            END IF
          END DO

        END IF

        ! Now along southern/northern boundaries

        ! process PEs not along southern boundary
        IF (.NOT. at_extremity(psouth)) THEN

          DO iproc=nproc-1,0,-1

            imsg=jobt*1000+iproc  ! message tag

            IF (mype == iproc-nproc_x) THEN
              IF (.NOT. at_extremity(pnorth)) THEN
                ! this PE receives if not along northern boundary
                CALL gc_rrecv(imsg,1,iproc,istat,zzlatmx,zzlatmn)
              END IF
            ELSE IF (mype == iproc) THEN
              ! PEs not along southern boundary send to adjacent PE
              CALL gc_rsend(imsg,1,iproc-nproc_x,istat,zzlatmx,zzlatmn)
            END IF

          END DO

        ELSE
          ! PEs along southern boundary can receive only

          DO iproc=nproc-1,0,-1
            IF (mype == iproc-nproc_x) THEN        ! this PE receives
              imsg=jobt*1000+iproc
              CALL gc_rrecv(imsg,1,iproc,istat,zzlatmx,zzlatmn)
            END IF
          END DO

        END IF


        IF (model_type /= mt_global) THEN
          !       Rotate real Lat/Lon of obs to ELF co-ordinates.
          !       Write transformed Lat,Lon back to same area of work array.
          !       Output Longitudes from rotate_latlon_to_eq are in range 
          !       0-->360 degrees.
          CALL rotate_latlon_to_eq                                        &
          (work(iref(jobt,ip_lat,jf)+1),work(iref(jobt,ip_long,jf)+1),    &
           work(iref(jobt,ip_lat,jf)+1),work(iref(jobt,ip_long,jf)+1),    &
           elfplat,elfplon,obs_info % nobs(jobt) )
        END IF

        DO job=1,obs_info % nobs(jobt)

          !       Reset Observation time so times are relative to start of assm.

          work(iref(jobt,ip_time,jf)+job) =                               &
          work(iref(jobt,ip_time,jf)+job)+timeadj

          !       Test if observation in area and time window.
          !   -0.5 on timewindow helps make results same on different machines

          IF ( zzlongmn  <   zzlongmx ) THEN

            IF ( work(iref(jobt,ip_lat, jf)+job)  <   zzlatmx .AND.       &
                 work(iref(jobt,ip_lat, jf)+job)  >=  zzlatmn .AND.       &
                 work(iref(jobt,ip_long,jf)+job)  <   zzlongmx  .AND.     &
                 work(iref(jobt,ip_long,jf)+job)  >=  zzlongmn  .AND.     &
                 work(iref(jobt,ip_time,jf)+job)  <   twend-0.5 .AND.     &
                 work(iref(jobt,ip_time,jf)+job)  >   twstart+0.5  ) THEN

              !           Count no of obs in area and time window and
              !           record those observations.

              inobs(jobt,jf)        = inobs(jobt,jf)+1
              iwork(inobs(jobt,jf)) = job

            END IF

          ELSE

            IF ( work(iref(jobt,ip_lat, jf)+job)  <   zzlatmx  .AND.      &
                 work(iref(jobt,ip_lat, jf)+job)  >=  zzlatmn  .AND.      &
               ( work(iref(jobt,ip_long,jf)+job)  <   zzlongmx .OR.       &
                 work(iref(jobt,ip_long,jf)+job)  >=  zzlongmn ) .AND.    &
                 work(iref(jobt,ip_time,jf)+job)  <   twend-0.5 .AND.     &
                 work(iref(jobt,ip_time,jf)+job)  >   twstart+0.5  ) THEN

              !           Count no of obs in area and time window and
              !           record those observations.

              inobs(jobt,jf)        = inobs(jobt,jf)+1
              iwork(inobs(jobt,jf)) = job

            END IF
          END IF

        END DO

        !       Compress out observations not required.
        DO jdv=1,obs_info % ndatav(jobt)
          DO job=1,inobs(jobt,jf)
            work(iref(jobt,jdv,jf)+job)=work(iref(jobt,jdv,jf)+iwork(job))
          END DO
        END DO

      END IF

    END DO ! jobt

  END IF

  !       Print no of observations in time window for file JF.
  DO jobt=1,obs_info % nobtyp
    CountA(jobt)=inobs(jobt,jf)
  END DO

  IF (mype == 0) THEN
    DO jobt=1,obs_info % nobtyp
      CountC(jobt)=0
    END DO
  END IF

  DO iproc=0,nproc-1
    IF (mype == 0) THEN
      CALL gc_irecv(iproc,obs_info % nobtyp,iproc,istat,countb,counta)
      DO job=1,obs_info % nobtyp
        CountC(job)=CountC(job)+CountB(job)
      END DO
    ELSE IF (mype == iproc) THEN
      CALL gc_isend(iproc,obs_info % nobtyp,0,istat,countb,counta)
    END IF
    !          IF (DIAG_RDOBS >  1.AND.mype == 0) THEN
    WRITE(umMessage,'(A,I5,T30,18I7/(T30,18I7))')                      &
          ' NO OF OBS IN T.W: pe=',iproc,                              &
          (CountB(jobt),jobt=1,obs_info % nobtyp)
    CALL umPrint(umMessage,src='rdobs2')
    !          ENDIF
  END DO
  !        If(mype == 0) THEN
  CALL umPrint(' ',src='rdobs2')
  WRITE(umMessage,'(A,T30,18I7/(T30,18I7))')                           &
      ' NO OF OBS IN T.W: total:',                                     &
      (CountC(jobt),jobt=1,obs_info % nobtyp)
  CALL umPrint(umMessage,src='rdobs2')
  !        ENDIF

  !       Compress out unused work space
  !       ==============================

  !       Following commented line may be brought back in future

  ipt = ipt_this_file
  DO jobt=1,obs_info % nobtyp
    DO jdv=1,obs_info % ndatav(jobt)
      IF (inobs(jobt,jf) >  0) THEN
        DO job=1,inobs(jobt,jf)
          work(ipt+job) = work(iref(jobt,jdv,jf)+job)
        END DO
      END IF
      iref(jobt,jdv,jf)=ipt
      ipt = ipt+inobs(jobt,jf)
    END DO
  END DO
  IF (diag_rdobs >= 1) THEN
    ipc = (ipt*100)/ac_tot_obs_size_max
    ! diagnostic message directed to every PEn output stream
    WRITE(umMessage,*)'RDOBS2:% of OBS array required=',ipc
    CALL umPrint(umMessage,src='rdobs2')
  END IF

  CALL model_file_close(ifile,used_files(jf),obs_info % filename_len(jf),&
                  envvar,0,icode)
  CALL release_file_unit(ifile, handler="portio")

END DO ! JF=1,NUM_USED_FILES

CALL deallocate_ac_dump_headers()

!-----------------------------------------------------------------------

!         DETERMINE LIST OF AC OBS TYPES FOR MERGED FILES
!         NOBTYP IS NOW THE NO OF AC OBS TYPES IN THE MERGED LIST

!         INITIAL LIST IS LIST FOR FIRST FILE
obs_info % nobtyp=inobtyp(1)
IF (obs_info % nobtyp >  0) THEN
  DO jobt=1,obs_info % nobtyp
    kobstyp(jobt) = iobstyp(jobt,1)
  END DO
END IF

IF (nfiles >  1) THEN

  !         GET FULL LIST FROM OTHER FILES
  DO jf=2,nfiles
    IF (.NOT. lempty(jf)) THEN
      DO jtyp=1,inobtyp(jf)
        DO jobt=1,obs_info % nobtyp
          IF (iobstyp(jtyp,jf) == kobstyp(jobt)) GO TO 2020
        END DO
        obs_info % nobtyp = obs_info % nobtyp+1
        IF (obs_info % nobtyp <= nobtypmx) THEN
          kobstyp(obs_info % nobtyp) = iobstyp(jtyp,jf)
        ELSE
          icode = 1
          cmessage =                                                  &
          ' RDOBS2 : MAX NO OF AC OBS TYPES REACHED FOR MERGED FILES'
          CALL umPrint( cmessage, src='rdobs2',pe=0)
          GO TO 9999
        END IF

        2020     CONTINUE
      END DO ! jtyp
    END IF
  END DO

END IF

!       NOBTYP is now set to the no of obs types to be assimilated
!       which is the same as NACT.
obs_info % nobtyp = nact

DO jobt=1,obs_info % nobtyp
  obs_info % obstyp(jobt)   = -1
  obs_info % nobs  (jobt)   =  0
  obs_info % ndatav(jobt)   = -1
  obs_info % noblev(jobt)   = -1
  obs_info % oblevtyp(jobt) = -1
  DO jlev=1, obs_info % maxnlev1
    datalevs(jlev,jobt) = obs_info % missd
  END DO
END DO

IF (obs_info % nobtyp >  0) THEN

  !         OBSTYP is now the list of obs types to be assimilated
  !         which is the same as LACT.
  DO jact=1,nact
    obs_info % obstyp(jact) = lact(jact)
  END DO

  DO jobt=1,obs_info % nobtyp
    DO jf=1,nfiles
      IF (.NOT. lempty(jf)) THEN
        DO jtyp=1,inobtyp(jf)
          IF (iobstyp(jtyp,jf) == obs_info % obstyp(jobt)) THEN

            obs_info % ndatav(jobt)   = indatav(jtyp,jf)
            obs_info % noblev(jobt)   = inoblev(jtyp,jf)
            obs_info % oblevtyp(jobt) = ioblvtyp(jtyp,jf)

            DO jlev=1, obs_info % maxnlev1
              datalevs(jlev,jobt) = plevels(jlev,jtyp,jf)
            END DO
            GO TO 2110

          END IF
        END DO
      END IF
    END DO
    2110       CONTINUE
  END DO

END IF

!       Check that data is consistent for obs types being assimilated.

IF (nfiles >  1 .AND. obs_info % nobtyp >  0) THEN

  DO jf=1,nfiles

    !         ChecK NDVHDR

    IF (.NOT. lempty(jf) .AND. indvhdr(jf) /= obs_info % ndvhdr) THEN
      icode = 1
      cmessage =' RDOBS2 : MIS-MATCH IN AC OBS FILES - NDVHDR'
      CALL umPrint( ' ',src='rdobs2',pe=0)
      WRITE(umMessage,*)' RDOBS2 : DIFFERENT VALUES OF NDVHDR ?'
      CALL umPrint(umMessage,src='rdobs2',pe=0)
      WRITE(umMessage,'(A,5I5)')' NDVHDR =',(indvhdr(ijf),ijf=1,nfiles)
      CALL umPrint(umMessage,src='rdobs2',pe=0)
      GO TO 9999
    END IF

    !         Check MAXNLEV1

    IF (.NOT. lempty(jf) .AND. imaxnlev(jf) /= obs_info % maxnlev1) THEN
      icode = 1
      cmessage =' RDOBS2 : MIS-MATCH IN AC OBS FILES - MAXNLEV1'
      CALL umPrint(' ',src='rdobs2',pe=0)
      WRITE(umMessage,*)' RDOBS2 : DIFFERENT VALUES OF MAXNLEV1'
      CALL umPrint(umMessage,src='rdobs2',pe=0)
      WRITE(umMessage,'(A,5I5)') &
          ' MAXNLEV1 =',(imaxnlev(ijf),ijf=1,nfiles)
      CALL umPrint(umMessage,src='rdobs2',pe=0)
      GO TO 9999
    END IF

  END DO

  DO jobt=1,obs_info % nobtyp
    DO jf=1,nfiles
      IF (inobtyp(jf) >  0) THEN
        DO jtyp=1,inobtyp(jf)

          IF (iobstyp(jtyp,jf) == obs_info % obstyp(jobt)) THEN

            !               Check no of data values (NDATAV)

            IF (indatav(jtyp,jf) /= obs_info % ndatav(jobt)) THEN
              icode = 1
              cmessage =' RDOBS2 : MIS-MATCH IN OBS FILES - NDATAV'
              WRITE(umMessage,*) ' RDOBS2 : Different No of Data Values.'
              CALL umPrint(umMessage,src='rdobs2',pe=0)
              WRITE(umMessage,*) '        : See Obs Type ', &
                  obs_info % obstyp(jobt)
              CALL umPrint(umMessage,src='rdobs2',pe=0)
              GO TO 9999
            END IF

            !               Check no of observation levels (NOBLEV)

            IF (inoblev(jtyp,jf) /= obs_info % noblev(jobt)) THEN
              icode = 1
              cmessage =' RDOBS2 : MIS-MATCH IN OBS FILES - NOBLEV'
              WRITE(umMessage,*) ' RDOBS2 : Different No of Obs levels.'
              CALL umPrint(umMessage,src='rdobs2',pe=0)
              WRITE(umMessage,*) '        : See Obs Type ', &
                  obs_info % obstyp(jobt)
              CALL umPrint(umMessage,src='rdobs2',pe=0)
              GO TO 9999
            END IF

            !               Check Observation level type (OBLEVTYP)

            IF (ioblvtyp(jtyp,jf) /= obs_info % oblevtyp(jobt)) THEN
              icode = 1
              cmessage =                                            &
              ' RDOBS2 : MIS-MATCH IN OBS FILES - OBLEVTYP'
              WRITE(umMessage,*)                                    &
                  ' RDOBS2 : Different Observation Level Type'
              CALL umPrint(umMessage,src='rdobs2',pe=0)
              WRITE(umMessage,*) '        : See Obs Type ',         &
                  obs_info % obstyp(jobt)
              CALL umPrint(umMessage,src='rdobs2',pe=0)
              GO TO 9999
            END IF

            !               Check pressure levels of observations (DATALEVS)

            DO jlev=1, obs_info % maxnlev1
              IF (plevels(jlev,jtyp,jf) /= datalevs(jlev,jobt)) THEN
                icode = 1
                cmessage =                                           &
                ' RDOBS2 : MIS-MATCH IN OBS FILES - DATALEVS'
                WRITE(umMessage,*) ' RDOBS2 : Different levels in Obs files'
                CALL umPrint(umMessage,src='rdobs2',pe=0)
                WRITE(umMessage,*) '        : See Obs Type ', &
                    obs_info % obstyp(jobt)
                CALL umPrint(umMessage,src='rdobs2',pe=0)
                WRITE(umMessage,'(A,I5,A,2F10.1)') ' LEVEL',jlev,           &
                    ' ; Pressures =',plevels(jlev,jtyp,jf),              &
                    datalevs(jlev,jobt)
                CALL umPrint(umMessage,src='rdobs2',pe=0)
                GO TO 9999
              END IF
            END DO

          END IF

        END DO  !  End of loop over JTYP
      END IF
    END DO  !  End of loop over JF
  END DO  !  End of loop over JOBT

END IF

!-----------------------------------------------------------------------
!            SECTION 2: MERGE INPUT DATA & SET UP HEADER.
!-----------------------------------------------------------------------

!       Merge observation data and put into output buffer OBS
tndv=0
DO jobt=1,obs_info % nobtyp
  DO jdv=1,obs_info % ndatav(jobt)
    DO jf=1,nfiles
      DO jtyp=1,inobtyp(jf)
        IF (iobstyp(jtyp,jf) == obs_info % obstyp(jobt) .AND.   &
            inobs(jtyp,jf)   >  0) THEN

          IF (tndv+inobs(jtyp,jf) <= ac_tot_obs_size_max) THEN

            DO job=1,inobs(jtyp,jf)
              obs(tndv+job) = work(iref(jtyp,jdv,jf)+job)
            END DO
            tndv = tndv+inobs(jtyp,jf)

          ELSE

            icode = 1
            cmessage =                                          &
            ' RDOBS2 : Insufficient space in OBS array'
            ! important failure messages to PEn,OUT and OUT2 output streams
            CALL umPrint(cmessage,src='rdobs2')
            WRITE(umMessage,*) ' Recode or rerun with fewer obs'
            CALL umPrint(umMessage,src='rdobs2')
            GO TO 9999

          END IF
        END IF
      END DO   !  Loop over JTYP
    END DO   !  Loop over JF
  END DO   !  Loop over JDV
END DO   !  Loop over JOBT
WRITE(umMessage,*)'temp check OBS ',tndv,ac_tot_obs_size_max,          &
                                   (tndv*ac_tot_obs_size_max)*100.0,'%'
CALL umPrint(umMessage,src='rdobs2')

!       Count total no of obs for each obs type
DO jobt=1,obs_info % nobtyp
  obs_info % nobs(jobt)=0
  DO jf=1,nfiles
    DO jtyp=1,inobtyp(jf)
      IF (iobstyp(jtyp,jf) == obs_info % obstyp(jobt)) THEN
        obs_info % nobs(jobt) = obs_info % nobs(jobt)+inobs(jtyp,jf)
      END IF
    END DO
  END DO
END DO

!       Get total no of obs (TNOBS)
tnobs = 0
DO jobt=1,obs_info % nobtyp
  tnobs = tnobs + obs_info % nobs(jobt)
END DO

!       Set up pointers to start of data for each obs type
obs_info % mdispobt(1)=0
DO jobt=2,obs_info % nobtyp
  obs_info % mdispobt(jobt) = obs_info % mdispobt(jobt-1) +       &
                              obs_info % nobs(jobt-1) *           &
                              obs_info % ndatav(jobt-1)
END DO

!       Offset to first obs for each type
obs_info % obs_no_st(1)=0
DO jobt=2,obs_info % nobtyp
  obs_info % obs_no_st(jobt) = obs_info % obs_no_st(jobt-1) +   &
                               obs_info % nobs(jobt-1)
END DO


!       Check no of observations against maximum allowed (AC_NUM_OBS_MAX)
IF (tnobs >  ac_num_obs_max) THEN
  icode = 1
  cmessage = ' RDOBS2 : Insufficient space in OBS_FLAG array'
  ! check in every pe?
  WRITE(umMessage,*)' RDOBS2 : Insufficient space in OBS_FLAG array'
  CALL umPrint(umMessage,src='rdobs2')
  WRITE(umMessage,*)' Recode or rerun with fewer obs'
  CALL umPrint(umMessage,src='rdobs2')
  WRITE(umMessage,'(A,I8)')' NO OF OBS   = ',tnobs
  CALL umPrint(umMessage,src='rdobs2')
  WRITE(umMessage,'(A,I8)')' MAX ALLOWED = ',ac_num_obs_max
  CALL umPrint(umMessage,src='rdobs2')
  GO TO 9999
ELSE
  IF (diag_rdobs >= 1) THEN
    ipc = (tnobs*100)/ac_num_obs_max
    ! diagnostic message directed to every PEn output stream
    WRITE(umMessage,*)'RDOBS2:% of OBS_FLAG array required=',ipc
    CALL umPrint(umMessage,src='rdobs2')
  END IF
END IF

!       Print summary of output file.
IF (diag_rdobs >  0) THEN
  DO jobt=1,obs_info % nobtyp
    CountA(jobt)=obs_info % nobs(jobt)
  END DO

  IF (mype == 0) THEN
    CALL umPrint(' ',src='rdobs2')
    WRITE(umMessage,'(A,I3.2,A,I2.2,A,I4)')                              &
    ' REF DATE :',obs_info % obs_ref_dd,'/',                      &
                  obs_info % obs_ref_mm,'/',                      &
                  obs_info % obs_ref_yy
    CALL umPrint(umMessage,src='rdobs2')
    WRITE(umMessage,'(A,I3.2,I2.2,A)')                                   &
    ' REF TIME :',obs_info % obs_ref_hh, obs_info % obs_ref_min,'Z'
    CALL umPrint(umMessage,src='rdobs2')
    WRITE(umMessage,'(A,F8.1,A)') ' REL TIME      :',timerel,'M'
    CALL umPrint(umMessage,src='rdobs2')
    WRITE(umMessage,'(A,F8.1,A)') ' REL T.W START :',twstart,'M'
    CALL umPrint(umMessage,src='rdobs2')
    WRITE(umMessage,'(A,F8.1,A)') ' REL T.W END   :',twend,'M'
    CALL umPrint(umMessage,src='rdobs2')
    CALL umPrint(' ',src='rdobs2')
    WRITE(umMessage,'(A,T30,18I6/(T30,18I6))')                           &
        ' AC OBS TYPES :',(obs_info % obstyp(jobt),             &
        jobt=1,obs_info % nobtyp)
    CALL umPrint(umMessage,src='rdobs2')
    CALL umPrint(' ',src='rdobs2')
    DO jobt=1,obs_info % nobtyp
      CountC(jobt)=0
    END DO
  END IF

  DO iproc=0,nproc-1
    IF (mype == 0) THEN
      CALL gc_irecv(iproc,obs_info % nobtyp,iproc,istat,countb,counta)
      DO job=1,5
        CountC(job)=CountC(job)+CountB(job)
      END DO
    ELSE IF (mype == iproc) THEN
      CALL gc_isend(iproc,obs_info % nobtyp,0,istat,countb,counta)
    END IF

    IF (diag_rdobs >  1) THEN
      WRITE(umMessage,'(A,I5,T30,18I6/(T30,18I6))')                        &
          ' NO OF OBS    : pe=',iproc,(CountB(jobt),jobt=1,         &
          obs_info % nobtyp)
      CALL umPrint(umMessage,src='rdobs2')
    END IF
  END DO
  WRITE(umMessage,'(A,T30,18I6/(T30,18I6))')                  &
      ' NO OF OBS    : total:',(CountC(jobt),jobt=1,           &
      obs_info % nobtyp)
  CALL umPrint(umMessage,src='rdobs2',pe=0)

  !        Get total no of obs (all pe's)
  IF (obs_info % nobtyp >= 2) THEN
    DO jobt=2,obs_info % nobtyp
      CountC(1) = CountC(1) + CountC(jobt)
    END DO
  END IF  !NOBTYP >= 2
  IF (mype == 0) THEN
    CALL umPrint(' ',src='rdobs2')
    WRITE(umMessage,'(A,I8,A)') ' TOTAL NO OF OBS IN T.W :',CountC(1)
    CALL umPrint(umMessage,src='rdobs2')
    !        Any observations to assimilate ?
    IF (CountC(1) == 0) THEN
      CALL umPrint(' ',src='rdobs2')
      CALL umPrint('Timestep '//TRIM(str(timestep_no)),src='rdobs2')
      CALL umPrint('There are no observations to be assimilated.', &
          src='rdobs2')
    END IF !TNOBS == 0
  END IF

END IF

IF (obs_info % nobtyp >  0) THEN

  DO jobt=1,obs_info % nobtyp

    IF (obs_info % nobs(jobt) >  0) THEN

      !             Convert Obs Latitudes (degrees) to Co-latitudes (radians)
      ip_lat = obs_info % mdispobt(jobt)
      DO job=1,obs_info % nobs(jobt)
        obs(ip_lat+job) = (90.0-obs(ip_lat+job))*pi_over_180
      END DO

      !             Convert Obs Longitudes to Radians in range 0 to 2*PI.
      ip_long = obs_info % mdispobt(jobt) + obs_info % nobs(jobt)
      DO job=1,obs_info % nobs(jobt)
        IF (obs(ip_long+job)  <   0.0)                          &
        obs(ip_long+job) = obs(ip_long+job)+360.0
        obs(ip_long+job) = obs(ip_long+job)*pi_over_180
      END DO

    END IF

  END DO

  !         SET UP OBLEVELS AND OBLAYERB (LEVELS IN PASCALS)
  !         ------------------------------------------------
  DO jobt=1,obs_info % nobtyp

    IF (obs_info % oblevtyp(jobt) == 3 .OR.                       &
        obs_info % oblevtyp(jobt) == 4) THEN

      IF (obs_info % oblevtyp(jobt) == 3) THEN

        DO jlev=1,obs_info % noblev(jobt)
          obs_info % oblevels(jlev,jobt) = datalevs(jlev,jobt)
        END DO

        DO jlev=2,obs_info % noblev(jobt)
          obs_info % oblayerb(jlev,jobt) =                          &
          SQRT ( obs_info % oblevels(jlev-1,jobt) *                 &
                 obs_info % oblevels(jlev,jobt) )
        END DO

        obs_info % oblayerb(1,jobt) =  obs_info % oblevels(1,jobt) *&
         obs_info % oblevels(1,jobt) / obs_info % oblayerb(2,jobt)

        obs_info % oblayerb(obs_info % noblev(jobt)+1,jobt) =       &
               obs_info % oblevels(obs_info % noblev(jobt),jobt) *  &
               obs_info % oblevels(obs_info % noblev(jobt),jobt) /  &
               obs_info % oblayerb(obs_info % noblev(jobt),jobt)

      END IF

      IF (obs_info % oblevtyp(jobt) == 4) THEN

        DO jlev=1,obs_info % noblev(jobt)+1
          obs_info % oblayerb(jlev,jobt) = datalevs(jlev,jobt)
        END DO

        DO jlev=1,obs_info % noblev(jobt)
          obs_info % oblevels(jlev,jobt) =                          &
          SQRT ( obs_info % oblayerb(jlev,jobt) *                   &
                 obs_info % oblayerb(jlev+1,jobt) )
        END DO

      END IF

      IF (mype == 0) THEN
        CALL umPrint(' ',src='rdobs2')
        WRITE(umMessage,*) ' Observation Type ',obs_info % obstyp(jobt)
        CALL umPrint(umMessage,src='rdobs2')
        CALL umPrint(' ',src='rdobs2')
        CALL umPrint('  Levels (mb) =',src='rdobs2')
        WRITE(umMessage, '(1X,5F8.1)')                   &
            (obs_info % oblevels(jlev,jobt)*0.01,                  &
            jlev=1,obs_info % noblev(jobt))
        CALL umPrint(umMessage,src='rdobs2')
        CALL umPrint('  Layer boundaries (mb) =',src='rdobs2')
        WRITE(umMessage,'(1X,5F8.1)')        &
            (obs_info % oblayerb(jlev,jobt)*0.01,                  &
            jlev=1,obs_info % noblev(jobt)+1)
        CALL umPrint(umMessage,src='rdobs2')

      END IF
    END IF

  END DO

END IF
!     ====================================================
!     SET UP ARRAY NERLEV1 WHICH POINTS THE FIRST DATA VALUE
!     CORRESPONDING TO THE FIRST LEVEL OF ERROR RATIO FOR THE OBS TYPE

IF (obs_info % nobtyp >  0) THEN
  DO jobt=1,obs_info % nobtyp
    obs_info % nerlev1(jobt) = obs_info % ndatav(jobt) -          &
                               obs_info % noblev(jobt)+1
  END DO
END IF
!     ====================================================
IF (nact >  0) THEN

  WRITE(umMessage,'(A,15I5)')'0TYPES TO BE PROCESSED : LACT    =', &
      (lact(j),j=1,nact)
  CALL umPrint(umMessage,src='rdobs2',pe=0)
  WRITE(umMessage,'(A,15I5)')'          NO OF LEVELS : NOBLEV  =', &
      (obs_info % noblev(j),j=1,nact)
  CALL umPrint(umMessage,src='rdobs2',pe=0)
  WRITE(umMessage,'(A,15I5)')'  NO OF DATA VARIABLES : NDATAV  =', &
      (obs_info % ndatav(j),j=1,nact)
  CALL umPrint(umMessage,src='rdobs2',pe=0)
  WRITE(umMessage,'(A,15I5)')'       FIRST ERROR LEV : NERLEV1 =', &
      (obs_info % nerlev1(j),j=1,nact)
  CALL umPrint(umMessage,src='rdobs2',pe=0)

END IF
!     ====================================================
!     CALL SETDAC TO SET UP FOR ANY DIAGNOSTICS REQUIRED
CALL setdac (obs,tndv)
!     ====================================================

!     RE-USE IWORK IN THIS SECTION

IF (tnobs >  0) THEN

  DO job=1,tnobs
    obs_flag(job)=0
  END DO

  !       Convert Analysis Types into INTEGER.

  ip_type = 4
  DO jobt=1,obs_info % nobtyp
    IF (obs_info % nobs(jobt) >  0) THEN
      ip_mot = obs_info % mdispobt(jobt) +                        &
               (ip_type-1) * obs_info % nobs(jobt)
      DO job=1,obs_info % nobs(jobt)
        iwork(obs_info % obs_no_st(jobt)+job) = NINT( obs(ip_mot+job) )
      END DO
    END IF
  END DO

  itot0=0
  itot1=0
  DO jmot=1,nanaltyp
    IF (iomitobs(jmot) >  0) THEN
      DO job=1,tnobs
        IF (iwork(job) == iomitobs(jmot)) obs_flag(job)=1
      END DO
      IF (diag_rdobs >= 1) THEN
        itot1=0
        DO job=1,tnobs
          IF (obs_flag(job) == 1) itot1=itot1+1
        END DO
        itotal=itot1-itot0
        ! in every pe?
        WRITE(umMessage,'(A,I5,A,I6)')' ANALYSIS TYPE ',iomitobs(jmot),  &
            ' : NO OF OBSERVATIONS FLAGGED (DO NOT USE)  -',itotal
        CALL umPrint(umMessage,src='rdobs2')
        itot0=itot1
      END IF
    END IF
  END DO


  !       OPTIONAL PRINT OUT.
  IF (diag_rdobs >= 1 .AND. itot1 /= 0) THEN
    ! in every pe?
    WRITE(umMessage, '(A,I6)')                                           &
        ' TOTAL NO OF OBSERVATIONS FLAGGED (DO NOT USE)  -',itot1
    CALL umPrint(umMessage,src='rdobs2')

  END IF

END IF
!-----------------------------------------------------------------------


9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rdobs2
END MODULE rdobs2_mod

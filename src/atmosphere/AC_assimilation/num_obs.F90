! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE NUM_OBS  -----------------------------------------------
!
!    Purpose : Calculate no of observations and no of observation
!              values in AC observation files. These values are
!              used to dimension any observation data arrays in the
!              assimilation code.
!
!    Programming standard: Unified Model Documentation Paper No. 3
!
!    Project Task : P3
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE num_obs_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

INTEGER :: ac_num_obs_max         ! Total no of obs in obs files
INTEGER :: ac_tot_obs_size_max ! Total no of data values in obs files

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NUM_OBS_MOD'

CONTAINS

SUBROUTINE num_obs (nfiles,                                       &
                    p_rows,row_length,ak,bk,realhd1,realhd2,      &
                    realhd3,realhd4,realhd5,realhd6,              &
                    icode,cmessage)
! ----------------------------------------------------------------------

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io, ONLY: setpos
USE io_constants, ONLY: ioOpenReadOnly, ioNameProvided, ioAllLocal, ioNoDelete
USE model_file, ONLY: model_file_open, model_file_close
USE UM_ParVars
USE UM_ParCore,   ONLY: mype, nproc
USE comobs_mod
USE num_obs2_mod, ONLY: num_obs2
USE umPrintMgr, ONLY:      &
    umPrint,                &
    str,                    &
    umMessage
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE nlcfiles_namelist_mod, ONLY: &
    obs01_file => obs01,  &
    obs02_file => obs02,  &
    obs03_file => obs03,  &
    obs04_file => obs04,  &
    obs05_file => obs05
USE filenamelength_mod, ONLY: filenamelength
USE ac_control_mod
USE ac_dump_mod, ONLY: len_inthd, len_realhd, len1_levdepc, len2_levdepc, &
                       len1_rowdepc, len2_rowdepc, len1_coldepc,          &
                       len2_coldepc, len1_flddepc, len2_flddepc,          &
                       len_extcnst, len_dumphist, len_cfi1, len_cfi2,     &
                       len_cfi3, len1_lookup_obs, len2_lookup_obs, fixhd, &
                       allocate_ac_dump_headers,                          &
                       deallocate_ac_dump_headers
USE nlsizes_namelist_mod, ONLY: model_levels, len_fixhd

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------
!     ARGUMENTS
!-----------------------------------------------------------------------
INTEGER :: nfiles        ! IN  No of obs files to be read
INTEGER :: p_rows        ! IN  No of model rows
INTEGER :: row_length    ! IN  No of points on row
REAL :: ak(model_levels) ! IN  Vertical grid
REAL :: bk(model_levels) !
REAL :: realhd1,realhd2  ! IN  Horizontal grid
REAL :: realhd3,realhd4
REAL :: realhd5,realhd6
INTEGER :: icode          ! OUT Return code
CHARACTER(LEN=errormessagelength) :: cmessage ! OUT Error message
!-----------------------------------------------------------------------
!     LOCAL VARIABLES
!-----------------------------------------------------------------------
INTEGER :: idummy
INTEGER :: io_stat
REAL :: datalevs(model_levels+1,nobtypmx)
INTEGER :: jf,jobt,jlev   !  Loop counters
INTEGER :: unit_no        !  Unit no for obs file
INTEGER :: ispare(7)
INTEGER :: istat          ! GCOM status
!-----------------------------------------------------------------------
INTEGER :: tndv, indv  (max_num_acob_files)
INTEGER :: tnobs,inobs (max_num_acob_files)
!-----------------------------------------------------------------------
INTEGER :: len_data
!-----------------------------------------------------------------------
INTEGER :: obs_file_yy,obs_file_mm,obs_file_dd
INTEGER :: obs_file_hh,obs_file_min,obs_file_sec
INTEGER :: dirnum, filenum
CHARACTER(LEN=filenamelength) :: obs_file_name
CHARACTER(LEN=filenamelength) :: obs_dir_name
LOGICAL :: lfile_exists
INTEGER :: len_obs_file_name
INTEGER :: ifound
INTEGER :: ibcast_buf(4)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUM_OBS'


!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ltimer_ac) CALL timer ('NUMOBS  ',3)
!-----------------------------------------------------------------------

! All obs files opened on the same unit number.
CALL assign_file_unit("observation_files", unit_no, handler="portio")

IF (mype == 0) THEN
  CALL umPrint('READ IN AC OBS FILES - Headers only',src='num_obs')

  IF (obs_format  ==  3) THEN
    nfiles = nfiles * num_ob_file_types
  END IF

  obs_info % num_used_files = 0

  DO jf=1,nfiles   !  Loop over observation files

    IF (obs_format  ==  3) THEN
      dirnum=INT((jf-1)/(num_ob_file_types))+1
      filenum=MOD(jf-1,num_ob_file_types)+1

      SELECT CASE(dirnum)
      CASE (1)
        obs_dir_name = obs01_file
      CASE (2)
        obs_dir_name = obs02_file
      CASE (3)
        obs_dir_name = obs03_file
      CASE (4)
        obs_dir_name = obs04_file
      CASE (5)
        obs_dir_name = obs05_file
      CASE DEFAULT
        cmessage = &
          'NUM_OBS: directory number not supported'
        GO TO 9999
      END SELECT

      IF (LEN_TRIM(obs_dir_name) + LEN_TRIM(ob_file_type(filenum)) + 6 &
            > filenamelength) THEN
        icode = 10
        cmessage = &
          'NUM_OBS: Full filename of obs_dir_name too long to store'
        GO TO 9999
      END IF
      obs_file_name = TRIM(obs_dir_name) // '/' // &
                      TRIM(ob_file_type(filenum)) // '.acobs'
      len_obs_file_name = LEN_TRIM(obs_file_name)
    ELSE

      SELECT CASE(jf)
      CASE (1)
        obs_file_name = obs01_file
      CASE (2)
        obs_file_name = obs02_file
      CASE (3)
        obs_file_name = obs03_file
      CASE (4)
        obs_file_name = obs04_file
      CASE (5)
        obs_file_name = obs05_file
      CASE DEFAULT
        cmessage = &
          'NUM_OBS: file number not supported'
        GO TO 9999
      END SELECT
      len_obs_file_name = LEN_TRIM(obs_file_name)
    END IF

    INQUIRE(FILE=obs_file_name,IOSTAT=icode,EXIST=lfile_exists)
    IF (icode  /=  0) GO TO 9999

    IF (lfile_exists) THEN
      WRITE(umMessage,*)'Obs File name: ',TRIM(obs_file_name), &
                         len_obs_file_name
      CALL umPrint(umMessage,src='num_obs')

      obs_info % num_used_files = obs_info % num_used_files + 1
      used_files(obs_info % num_used_files) = obs_file_name
      obs_info % filename_len(obs_info % num_used_files) =  &
                 len_obs_file_name

      ! Specify the locality of the file such as to avoid read broadcast
      ! replication because this code is in a mype==0 clause
      CALL model_file_open(unit_no, obs_file_name, &
          read_write=ioOpenReadOnly, error=icode, iolocality=ioAllLocal)

      inobs(obs_info % num_used_files) = 0
      indv (obs_info % num_used_files) = 0

      CALL umPrint(' AC OBS FILE - UNIT NO :'//TRIM(str(unit_no)), &
          src='num_obs')

      obs_info % nobtyp = 0
      !       Read File contents
      !       If FILE is empty - proceed to read next file.

      !         Go to start of obs file

      CALL setpos (unit_no,0,icode)

      !         Read in fixed length header (FLH)
      ! DEPENDS ON: read_flh
      CALL read_flh (unit_no,fixhd,len_fixhd,                       &
          icode,cmessage)
      IF (icode >  0) GO TO 9999

      !         Get dimensions of all data set components from FLH
      ! DEPENDS ON: get_dim
      CALL get_dim (fixhd,                                         &
          len_fixhd, len_inthd, len_realhd,                        &
          len1_levdepc, len2_levdepc, len1_rowdepc, len2_rowdepc,  &
          len1_coldepc, len2_coldepc, len1_flddepc,len2_flddepc,   &
          len_extcnst,  len_dumphist,                              &
          len_cfi1, len_cfi2, len_cfi3,                            &
          len1_lookup_obs, len2_lookup_obs,                        &
          len_data)

      CALL allocate_ac_dump_headers()

      !         Get date/time of observation file
      obs_file_yy     = fixhd(21)
      obs_file_mm     = fixhd(22)
      obs_file_dd     = fixhd(23)
      obs_file_hh     = fixhd(24)
      obs_file_min    = fixhd(25)
      obs_file_sec    = fixhd(26)


      CALL num_obs2 (unit_no,p_rows,row_length,   &
                   ak,bk,realhd1,realhd2,realhd3,realhd4,         &
                   realhd5,realhd6,                               &
                   len_data,obs_info % nobtyp,tnobs,tndv,         &
                   icode,cmessage)

      IF (icode >  0) GO TO 9999

      !       Record values for this file
      inobs(obs_info % num_used_files) = tnobs
      indv (obs_info % num_used_files) = tndv

      WRITE(umMessage,*) ' '
      CALL umPrint(umMessage,src='num_obs')
      WRITE(umMessage,'(A,I3.2,I2.2,A,I4.2,A,I2.2,A,I4)')        &
          ' TIME and DATE         :',                  &
          obs_file_hh,obs_file_min,'Z',                &
          obs_file_dd,'/',obs_file_mm,'/',obs_file_yy
      CALL umPrint(umMessage,src='num_obs')
      WRITE(umMessage,'(A,I8)') ' No of Obs Types       :', &
          obs_info % nobtyp
      CALL umPrint(umMessage,src='num_obs')
      WRITE(umMessage,'(A,I8)') ' No of Observations    :',tnobs
      CALL umPrint(umMessage,src='num_obs')
      WRITE(umMessage,'(A,I8)') ' No of Data Values     :',tndv
      CALL umPrint(umMessage,src='num_obs')

      CALL model_file_close(unit_no, obs_file_name,             &
             delete=ioNoDelete, error=icode)
    END IF
    nfiles=obs_info % num_used_files

  END DO

  WRITE(umMessage,'('' Obs file no       '',5I8)')                       &
      (jf,jf=1,nfiles)
  CALL umPrint(umMessage,src='num_obs')
  WRITE(umMessage,'(  '' No of obs         '',5I8)')                     &
      (inobs(jf),jf=1,nfiles)
  CALL umPrint(umMessage,src='num_obs')
  WRITE(umMessage,'(  '' No of data values '',5I8)')                     &
      (indv(jf),jf=1,nfiles)
  CALL umPrint(umMessage,src='num_obs')

  !     Add up INOBS and INDV to get AC_NUM_OBS_MAX and AC_TOT_OBS_SIZE_MAX

  ac_num_obs_max = 0
  ac_tot_obs_size_max = 0

  DO jf=1,nfiles
    IF (inobs(jf) >= 0) THEN
      ac_num_obs_max = ac_num_obs_max + inobs(jf)
    END IF
    IF (indv(jf) >= 0) THEN
      ac_tot_obs_size_max = ac_tot_obs_size_max + indv (jf)
    END IF
    obs_info % per_file_tndvmax(jf)=indv (jf)
  END DO

  IF (ac_num_obs_max == 0) THEN
    ac_num_obs_max = 1   !  Reset to 1 for dynamic allocation
    WRITE(umMessage,*) ' AC_NUM_OBS_MAX = ',ac_num_obs_max
    CALL umPrint(umMessage,src='num_obs')
    WRITE(umMessage,*) ' Reset to 1 prevent allocation problems.'
    CALL umPrint(umMessage,src='num_obs')
  END IF

  IF (ac_tot_obs_size_max == 0) THEN
    ac_tot_obs_size_max = 1   !  Reset to 1 for dynamic allocation
    WRITE(umMessage,*) ' AC_TOT_OBS_SIZE_MAX = ',ac_tot_obs_size_max
    CALL umPrint(umMessage,src='num_obs')
    WRITE(umMessage,*) ' Reset to 1 prevent allocation problems.'
    CALL umPrint(umMessage,src='num_obs')
  END IF

  CALL umPrint(' Values calculated by NUM_OBS ',src='num_obs')
  CALL umPrint(' Total no of obs (AC_NUM_OBS_MAX)         = '//      &
      TRIM(str(ac_num_obs_max)),src='num_obs')
  CALL umPrint(' Total no of data values (AC_TOT_OBS_SIZE_MAX) = '// &
      TRIM(str(ac_tot_obs_size_max)),src='num_obs')

  9999    CONTINUE
  ibcast_buf(1)=nfiles
  ibcast_buf(2)=ac_num_obs_max
  ibcast_buf(3)=ac_tot_obs_size_max
  ibcast_buf(4)=icode

END IF ! if(mype == 0)

! Release the file unit again
CALL release_file_unit(unit_no, handler="portio")

CALL deallocate_ac_dump_headers()

! Broadcast the NUM_OBS output arguments to all pes
CALL gc_ibcast(1,4,0,nproc,istat,ibcast_buf)

nfiles              = ibcast_buf(1)
ac_num_obs_max      = ibcast_buf(2)
ac_tot_obs_size_max = ibcast_buf(3)
icode               = ibcast_buf(4)

CALL gc_ibcast(2, obs_info_int_len, 0, nproc, istat,              &
               obs_info % maxnlev1)

CALL gc_ibcast(3, obs_info_real_len, 0, nproc, istat,             &
               obs_info % missd)

IF (ltimer_ac) CALL timer ('NUMOBS  ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE num_obs
!-----------------------------------------------------------------------
END MODULE num_obs_mod

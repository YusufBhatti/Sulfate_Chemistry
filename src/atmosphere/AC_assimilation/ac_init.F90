! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE AC_INIT ------------------------------------------------
!
!    Purpose : Initialise and set up for assimilation.
!            : via calls to ACP_NAMEL,ACDIAG_NAMEL,NUM_OBS,SETTPS
!
!
!    Programming standard: Unified Model Documentation Paper No. 3
!
!    Project Task : P3
!
!
!    Arguments:---------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE ac_init_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='AC_INIT_MOD'

CONTAINS

SUBROUTINE ac_init (bl_levels, tr_levels,                         &
                    p_rows, u_rows, row_length,                   &
                    timestep,                                     &
                    basis_time_yy, basis_time_mm, basis_time_dd,  &
                    basis_time_hh, basis_time_min,                &
                    realhd1, realhd2, realhd3,                    &
                    realhd4, realhd5, realhd6,                    &
                    ak, bk, lambda_p,phi_p, no_lambda_p,          &
                    icode, cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE conversions_mod, ONLY: pi_over_180

USE filenamelength_mod, ONLY:                                          &
    filenamelength

USE UM_ParVars
USE UM_ParCore,      ONLY: mype, nproc
USE UM_ParParams,    ONLY: halo_type_no_halo
USE Field_Types,     ONLY: fld_type_p
USE atmos_max_sizes, ONLY: row_length_max
USE Control_Max_Sizes
USE comobs_mod, ONLY: nobtypmx, obs_info
USE acdiag_namel_mod, ONLY: acdiag_namel
USE acp_namel_mod, ONLY: acp_namel
USE num_obs_mod, ONLY: num_obs
USE settps_mod, ONLY: settps
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE ac_control_mod

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: l_regular, model_type, mt_global
USE nlsizes_namelist_mod, ONLY: model_levels

IMPLICIT NONE


INTEGER :: bl_levels
INTEGER :: tr_levels
INTEGER :: p_rows
INTEGER :: u_rows
INTEGER :: row_length
REAL :: timestep
INTEGER :: basis_time_yy
INTEGER :: basis_time_mm
INTEGER :: basis_time_dd
INTEGER :: basis_time_hh
INTEGER :: basis_time_min
REAL :: realhd1,realhd2,realhd3,realhd4,realhd5,realhd6
REAL :: ak(model_levels),bk(model_levels)
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: icode

INTEGER :: no_lambda_p
REAL :: lambda_p(no_lambda_p)
REAL :: phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)

!-INTENT=IN--------------------------------------------------------
!     BL_LEVELS     : TOTAL NUMBER OF boundary layer LEVELS
!     TR_LEVELS     : TOTAL NUMBER OF tracer LEVELS
!     ROW_LENGTH    : NUMBER OF POINTS ON ROW
!     P_ROWS        : NUMBER OF ROWS (FOR PSTAR)
!     U_ROWS        : NUMBER OF ROWS (FOR wind)
!     TIMESTEP      : TIMESTEP IN SECONDS
!     BASIS_TIME_## : DEFINES DATA TIME
!     REALHD1-6     : DEFINES HORIZONTAL GRID
!     AK,BK         : DEFINES VERTICAL GRID
!-INTENT=OUT-----------------------------------------------------
!     ICODE         : NON ZERO FOR FAILURE
!     CMESSAGE      : REASON FOR FAILURE
! ---------------------------------------------------------------------

INTEGER :: jfile
INTEGER :: iccode
CHARACTER(LEN=filenamelength) :: filename
REAL, ALLOCATABLE :: phi_p_global(:,:)
INTEGER :: istat,i,j
REAL :: phi_p_local(row_length,p_rows)
REAL :: xlats,xlonge ! used in variable grid runs
REAL :: dlat_n,dlat_s,dlon_w,dlon_e
INTEGER :: row_length_global,p_rows_global

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='AC_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

row_length_global=glsize(1,fld_type_p)
p_rows_global=glsize(2,fld_type_p)


!-----------------------------------------------------------------------

CALL umPrint(' IN AC_INIT',src='ac_init,pe=0')

!  2. Set up COMMG Variables and check dimensions
! For regular grids read values from header otherwise from coordinate arrays
IF (l_regular) THEN

  dlong  = realhd1
  dlat   = realhd2
  xlatn  = realhd3 + (p_rows_global-1)*dlat
  xlongw = realhd4

ELSE

  ! Gather phi and lambda coordinates for entire domain
  DO i=1,p_rows
    DO j=1,row_length
      phi_p_local(j,i)=phi_p(j,i)
    END DO
  END DO

  ALLOCATE (phi_p_global(row_length_global,p_rows_global))

  ! DEPENDS ON: gather_field
  CALL gather_field(phi_p_local,   phi_p_global,                  &
                    row_length,    p_rows,                        &
                    row_length_global,                            &
                    p_rows_global,                                &
                    fld_type_p,    halo_type_no_halo,             &
                    0,             gc_all_proc_group )

  ! In a variable grid set boundaries and DLAT and DLONG to be
  ! outer values
  ! Calculate on PE0 then broadcast
  ! DLAT and DLON are the outermost grid lengths
  IF (mype == 0) THEN
    dlong = (lambda_p(2)-lambda_p(1))/pi_over_180
    dlat = (phi_p_global(1,2)-phi_p_global(1,1))/pi_over_180
    xlatn = phi_p_global(1,p_rows_global)/pi_over_180
    xlongw = lambda_p(1)/pi_over_180
    xlats = phi_p_global(1,1)/pi_over_180
    xlonge = lambda_p(row_length_global)/pi_over_180
  END IF

  ! Broadcast to all PEs
  CALL gc_rbcast(1000,1,0,nproc,istat,dlong)
  CALL gc_rbcast(1001,1,0,nproc,istat,dlat)
  CALL gc_rbcast(1002,1,0,nproc,istat,xlatn)
  CALL gc_rbcast(1003,1,0,nproc,istat,xlats)
  CALL gc_rbcast(1004,1,0,nproc,istat,xlongw)
  CALL gc_rbcast(1005,1,0,nproc,istat,xlonge)

END IF

IF (.NOT. l_regular) THEN
  dlat_n=(phi_p(1,p_rows)-phi_p(1,p_rows-1))/pi_over_180
  dlat_s=(phi_p(1,2)-phi_p(1,1))/pi_over_180
  dlon_w=(lambda_p(2)-lambda_p(1))/pi_over_180
  dlon_e=(lambda_p(row_length)-lambda_p(row_length-1))/pi_over_180
END IF

IF (model_levels >   model_levels_max) THEN
  icode=1
  cmessage = ' ACINIT: MODEL_LEVELS_MAX too small'
  WRITE(umMessage,*)  cmessage,                                &
      ' MODEL_LEVELS=',model_levels,' MODEL_LEVELS_MAX=',model_levels_max
  CALL umPrint(umMessage,src='ac_init,pe=0')
  GO TO 999
END IF

IF (p_rows >   rows_max) THEN
  icode=1
  cmessage = ' ACINIT: ROWS_MAX too small'
  WRITE(umMessage,*)cmessage,                                &
      ' P_ROWS=',p_rows,' ROWS_MAX=',rows_max
  CALL umPrint(umMessage,src='ac_init,pe=0')
  GO TO 999
END IF

IF (row_length >   row_length_max) THEN
  icode=1
  cmessage = ' ACINIT: ROW_LENGTH_MAX too small'
  WRITE(umMessage,*)cmessage,                                &
      ' ROW_LENGTH=',row_length,' ROW_LENGTH_MAX=',row_length_max
  CALL umPrint(umMessage,src='ac_init,pe=0')
  GO TO 999
END IF

! Make sure western boundary longitude in range 0 - 360 degrees
IF (xlongw <  0.0) xlongw = xlongw + 360.0

! Real lat/lon of pseudo N. pole in degrees

IF (model_type /= mt_global) THEN
  elfplat = realhd5
  elfplon = realhd6
END IF  ! .NOT. GLOBAL

CALL umPrint('',src='ac_init')
WRITE(umMessage,'(A,(T25,10F10.6))') ' DLAT,DLONG,XLATN,XLONGW',       &
                                  dlat,dlong,xlatn,xlongw
CALL umPrint(umMessage,src='ac_init')

!  3. ACP namelist. Set defaults, read in and process.
CALL acp_namel ( bl_levels, tr_levels,         &
                 p_rows, u_rows, row_length,   &
                 timestep, icode, cmessage)

IF (icode >  0) GO TO 999

!  4. ADIAG namelist. Set defaults, read in and process.
CALL acdiag_namel (icode,cmessage)
IF (icode >  0) GO TO 999

!  6. Read in AC Obs Files and compute number of observations

CALL num_obs (no_obs_files,                                       &
              p_rows,row_length,ak,bk,realhd1,realhd2,            &
              realhd3,realhd4,realhd5,realhd6,                    &
              icode,cmessage)
IF (icode >  0) GO TO 999

!  7. Set up list of obs types and groups to be processed this run
CALL settps (icode,cmessage)
IF (icode >  0) GO TO 999

!  8. Set up Analysis Grid in COMAG
! Initialise COMAG variables which don't change during a run.

! Model Grid Spacings
dlatmg  = dlat *pi_over_180
dlongmg = dlong*pi_over_180

! First row/point on p*/theta grid
row1mgth = ( 90.0 -xlatn )*pi_over_180
pt1mgth  = xlongw*pi_over_180

! Convert Lat at which pts/row starts decreasing to co-lat/radians.
aglatdec = (90.0-aglatdec)*pi_over_180

! Width of area according to model grid dimensions.
agrowlen = dlongmg*(row_length_global-1)


! The COMAG variables DEF_AGRES_ROWS and DEF_AGRES_PTS are
! initialised in DEF_GROUP. Use &ACP namelist arrays AGRES_ROWS
! and AGRES_PTS to change initialised values.

! The remaining COMAG variables are set in SETAG.

!  9. Set up COMOBS Variables
! Time Interval (mins) between reading AC Obs files
obs_info % timeint = 180.0

! Reference Time/Date which is start of assimilation
obs_info % obs_ref_yy  = basis_time_yy
obs_info % obs_ref_mm  = basis_time_mm
obs_info % obs_ref_dd  = basis_time_dd
obs_info % obs_ref_hh  = basis_time_hh
obs_info % obs_ref_min = basis_time_min

! Set up Latitudes and Longitudes of grid boundaries
! for observations to be used in assimilation.

! For Limited area assimilations, reject observations within
! one grid length of boundary.
IF (l_regular) THEN
  obs_info % obs_lat_n  = xlatn + 0.01 - dlat
  obs_info % obs_lat_s  = xlatn - 0.01 - (p_rows_global-2)*dlat
ELSE
  obs_info % obs_lat_n = xlatn + 0.01 - dlat_n
  obs_info % obs_lat_s = xlats - 0.01 + dlat_s
END IF
IF (l_regular) THEN
  obs_info % obs_long_w = xlongw - 0.01 + dlong
  obs_info % obs_long_e = xlongw + 0.01 + (row_length_global-2)*dlong
ELSE
  obs_info % obs_long_w = xlongw - 0.01 + dlon_w
  obs_info % obs_long_e = xlonge + 0.01 - dlon_e
END IF

! After rotation of the Lat/Long of Obs to ELF co-ords,
! longitude values will be in range 0-360 degrees, so
! make sure that boundary values are consistent. Note that
! if this leaves ZLONMN > ZLONMX, the test for obs in the area
! will assume that the Limited Area grid straddles the Meridian
! between ZLONMN and ZLONMX.
IF (obs_info % obs_long_w  <   0.0) THEN
  obs_info % obs_long_w = obs_info % obs_long_w + 360.0

ELSE IF (obs_info % obs_long_w  >   360.0) THEN
  obs_info % obs_long_w = obs_info % obs_long_w - 360.0

END IF

IF (obs_info % obs_long_e  <   0.0) THEN
  obs_info % obs_long_e = obs_info % obs_long_e + 360.0

ELSE IF (obs_info % obs_long_e  >   360.0) THEN
  obs_info % obs_long_e = obs_info % obs_long_e - 360.0

END IF

! Initialise time for next read in of observation files.
! Forces read on first timestep.
obs_info % timenext = -1440.0

999  CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ac_init

END MODULE ac_init_mod

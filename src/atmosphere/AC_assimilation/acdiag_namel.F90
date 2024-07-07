! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    MODULE ACDIAG_NAMEL_MOD  -------------------------------
!
!    Purpose : Read in ACDIAG Namelist and process.
!
!   ACDIAG_NAMEL:
!   Set defaults for ACDIAG namelist variables.
!                    Read in namelist and process.
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Logical system components covered:
!
!    Project Task : P3
!
!    External documentation:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE acdiag_namel_mod

USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

NAMELIST /acdiag/ ldiagac, lldac, lrms, ltemp, lverif

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACDIAG_NAMEL_MOD'

CONTAINS

SUBROUTINE acdiag_namel (icode,cmessage)

USE conversions_mod, ONLY: pi_over_180

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE comobs_mod, ONLY: nobtypmx
USE latlon_eq_rotation_mod, ONLY: rotate_latlon_to_eq

USE umPrintMgr, ONLY:      &
    umPrint,               &
    umMessage
USE check_iostat_mod
USE ereport_mod, ONLY: ereport

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE


! Global parameters:

! Exported variables (INTENT=OUT)
INTEGER ::    icode              ! Non zero for failure

CHARACTER(LEN=errormessagelength) :: cmessage           ! Reason for failure

REAL :: lat(4),lon(4)
REAL :: DAGLAT_temp(1)
REAL :: DAGLON_temp(1)
REAL :: xtemp
REAL :: zlat,zlong
INTEGER :: j

INTEGER :: errorstatus
CHARACTER (LEN=errormessagelength) :: iomessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'ACDIAG_NAMEL'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Hard-wired settings

lnormf   = .TRUE.

dlatn    =  55.0
dlats    =  50.0
dlongw   = -10.0
dlonge   = - 5.0

lldag0   = .FALSE.
lldag(1) = .FALSE.
lldag(2) = .FALSE.
lldag(3) = .FALSE.
lldag(4) = .FALSE.
lldag(5) = .FALSE.

DO j = 1, ndacop
  mdaco(j) = 0

END DO
modaco = 1

DO j = 1, nobtypmx
  mdiagtps(j) = 0

END DO

llband  = .FALSE.
mdiag   = 0
ndacprt = 10

daglat  = 57.0
daglon  = 340.0

ndgprt  = 6

!  3. Process the namelist
! Process ACDIAG namelist parameters and set up variables in
! COMACDG for diagnostic control in the assimilation.

! Check MODACO
IF (modaco <  1 .OR. modaco >  3) THEN
  lldac(1) = .FALSE.
  WRITE(umMessage,*)' Invalid value for MODACO - ',modaco,    &
      ' ; LLDAC(1) set to F'
  CALL umPrint(umMessage,src='acdiag_namel',pe=0)
END IF

! If 'obs-model' statistics required for verification purposes
! then reset various parameters in ACP and ACDIAG namelist.
IF (lverif) THEN
  lgeo    = .FALSE.   ! Do not calculate geostrophic increments
  lhydr   = .FALSE.   ! Do not calculate hydrostatic increments
  lhydrol = .FALSE.   ! Do not calculate hydrology   increments
  l_lhn = .FALSE.     ! Do not calculate LHN increments
  lnormf  = .FALSE.   ! Set normalisation factors to 1.0
  lwbal_sf= .FALSE.   ! Do not calculate P* & theta from windincs
  lwbal_ua= .FALSE.   !   "        "        "         "     "
  lrms    = .TRUE.    ! Print RMS values.
  ltemp   = .TRUE.    ! Use Temp Increments

  DO j = 1, nobtypmx
    def_tgetobb(j) = 180.0   ! Time Window of Obs - Before
    def_tgetoba(j) = 181.0   ! Time Window of Obs - After
    def_obthin(j) = 1        ! Dont reduce data vols by thinning

  END DO
END IF

IF (model_type /= mt_global) THEN
  lat(1) = dlatn
  lon(1) = dlongw
  lat(2) = lat(1)
  lon(2) = dlonge
  lat(3) = dlats
  lon(3) = lon(2)
  lat(4) = lat(3)
  lon(4) = lon(1)

  ! Transform corners of diagnostic area to elf co-ordinates
  CALL rotate_latlon_to_eq(lat,lon,lat,lon,elfplat,elfplon,4)

  ! Adjust area to be rectangular in elf lat/lon
  dlats  = 0.5* ( lat(3) + lat(4) )
  dlongw = 0.5* ( lon(1) + lon(4) )
  dlatn  = 0.5* ( lat(1) + lat(2) )
  dlonge = 0.5* ( lon(2) + lon(3) )

  ! Make sure dlatn >  dlats in new elf co-ords
  IF (dlats  >   dlatn) THEN
    xtemp = dlats
    dlats = dlatn
    dlatn = xtemp

  END IF

  ! The search for obs permits dlongw >  dlonge by assuming
  ! that an area crossing the meridian is intended, so no test
  ! as for latitudes.

  ! Transform diagnostic point to elf co-ordinates
  ! Move DAGLAT/LON into temp array as rotate_latlon_to_eq expects arrays
  DAGLAT_temp(1) = daglat
  DAGLON_temp(1) = daglon
  CALL rotate_latlon_to_eq(DAGLAT_temp,DAGLON_temp,DAGLAT_temp,DAGLON_temp, &
                           elfplat,elfplon,1)
  daglat = DAGLAT_temp(1)
  daglon = DAGLON_temp(1)
END IF  ! .NOT. GLOBAL

! Convert Latitudes into Co-Latitudes (0-180) (NP-SP)
dlats = 90.0 - dlats
dlatn = 90.0 - dlatn

! Convert Co-Latitudes to Radians
dlats = dlats * pi_over_180
dlatn = dlatn * pi_over_180

! Check Longitudes within Range (0-360 DEGREES)
IF (dlongw <  0.0) dlongw = dlongw+360.0
IF (dlonge <  0.0) dlonge = dlonge+360.0

! Convert Longitudes to Radians
dlongw = dlongw * pi_over_180
dlonge = dlonge * pi_over_180

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE acdiag_namel


SUBROUTINE print_nlist_acdiag()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist acdiag', &
    src='acdiag_namel')

WRITE(lineBuffer,*)' LDIAGAC = ',ldiagac
CALL umPrint(lineBuffer,src='acdiag_namel')
WRITE(lineBuffer,*)' LLDAC = ',lldac
CALL umPrint(lineBuffer,src='acdiag_namel')
WRITE(lineBuffer,*)' LRMS = ',lrms
CALL umPrint(lineBuffer,src='acdiag_namel')
WRITE(lineBuffer,*)' LTEMP = ',ltemp
CALL umPrint(lineBuffer,src='acdiag_namel')
WRITE(lineBuffer,*)' LVERIF = ',lverif
CALL umPrint(lineBuffer,src='acdiag_namel')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='acdiag_namel')

END SUBROUTINE print_nlist_acdiag

! ------------------------------------------------------------------

SUBROUTINE read_nml_acdiag(nml_unit)

USE um_parcore, ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE setup_namelist, ONLY: setup_nml_type
USE parkind1, ONLY: jprb,jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Argument
INTEGER :: nml_unit

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_ACDIAG'

INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode

CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_log = nldacp + 4

TYPE my_namelist
  SEQUENCE
  LOGICAL :: ldiagac
  LOGICAL :: lldac(nldacp)
  LOGICAL :: lrms
  LOGICAL :: ltemp
  LOGICAL :: lverif 
END TYPE my_namelist
        
TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_log_in=n_log)

!     Namelist defaults for ACDIAG

ldiagac   =.FALSE.
lldac(:) = .FALSE.
lrms     = .FALSE.
ltemp    = .FALSE.
lverif   = .FALSE.

IF (mype == 0) THEN

  READ (UNIT=nml_unit, NML=acdiag, IOSTAT=errorstatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist ACDIAG", iomessage)

  my_nml % ldiagac = ldiagac
  my_nml % lldac(:)= lldac(:)
  my_nml % lrms    = lrms
  my_nml % ltemp   = ltemp
  my_nml % lverif  = lverif

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  ldiagac = my_nml % ldiagac
  lldac(:)= my_nml % lldac(:)  
  lrms    = my_nml % lrms   
  ltemp   = my_nml % ltemp  
  lverif  = my_nml % lverif

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_acdiag

END MODULE acdiag_namel_mod

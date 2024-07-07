! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module containing all routines relating to 2D photolysis
!  used in UKCA sub-model.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_phot2d

USE asad_mod,            ONLY: jpspj
USE ukca_chem_defs_mod,  ONLY: ratj_t, ratj_defs
USE ukca_flupj_mod,      ONLY: ukca_flupj
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

USE conversions_mod,     ONLY: pi_over_180
USE umPrintMgr,          ONLY: PrintStatus, PrStatus_Diag, newline, &
                               umPrint, umMessage
USE ereport_mod,         ONLY: ereport
USE UM_ParVars
USE UM_ParCore,          ONLY: mype, nproc
USE file_manager,        ONLY: assign_file_unit, release_file_unit
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
SAVE
PRIVATE


INTEGER              :: myjppj  ! set = jppj below

INTEGER,PARAMETER,PUBLIC        :: nolat=19   ! no of 2D lats
INTEGER,PARAMETER,PUBLIC        :: nolev=17   ! no of 2D levs
INTEGER,PARAMETER,PUBLIC        :: nlphot=51  ! no of 2D photol levs
INTEGER,PARAMETER,PUBLIC        :: ntphot=3   ! no of times of day

!     Number of time intervals in data (74 corresponds to
!     a 5 day interval, with one extra set).

INTEGER,PARAMETER               :: n_data_times=74
REAL,   PARAMETER               :: data_interval=5.0 ! interval in days

LOGICAL, PARAMETER              :: L_optimised=.TRUE.  ! T for optimised read MONSooN

!     Photolysis rates from 2-D model which have been interpolated onto 3-D

REAL,ALLOCATABLE,PUBLIC          :: pjin(:,:,:,:)

! 2D model level pressures
REAL,PUBLIC :: pr2d(nolev)

PUBLIC ukca_photin, ukca_curve, ukca_inpr2d

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_PHOT2D'

CONTAINS

!-----------------------------------------------------------------------

SUBROUTINE ukca_photin(i_day_number,row_lengthda,tot_p_rows,     &
                p_levelsda,p_fieldda,first_row, glon, glat,      &
                sinlat,pl,jppj)

! Purpose: Subroutine to read in 2D photolysis rates, and interpolate on
!          3-d latitudes, heights. reconstruct daily curves (every 5 day
!          plus code to account for hour angle (longitude)
!
!          3 values are stored for each level and each latitude
!          symmetrical distribution plus zero values at dawn and
!          dusk gives 7 points altogether.
!          Based on PHOTIN.F from Cambridge TOMCAT model.

! dataset: 51(levels)x19(latitudes)x3(points)x74(every 5days)
!
!          Called from UKCA_CHEMISTRY_CTL.
!
! ---------------------------------------------------------------------
!
USE ereport_mod, ONLY: ereport
USE UM_ParVars
IMPLICIT NONE

INTEGER, INTENT(IN) :: jppj
INTEGER, INTENT(IN) :: i_day_number
INTEGER, INTENT(IN) :: row_lengthda
INTEGER, INTENT(IN) :: tot_p_rows
INTEGER, INTENT(IN) :: p_levelsda
INTEGER, INTENT(IN) :: p_fieldda
INTEGER, INTENT(IN) :: first_row
INTEGER, INTENT(IN) :: glon
INTEGER, INTENT(IN) :: glat

REAL, INTENT(IN) :: sinlat(:)
REAL, INTENT(IN) :: pl(:,:,:)

! Local variables

INTEGER :: i,j
INTEGER :: idofy                         ! Day number
INTEGER :: ipos
LOGICAL,SAVE  :: first_call=.TRUE.

REAL :: fpos
REAL :: pr2dj(nlphot)                    ! 2D photolysis level p
REAL :: zmean_pl(tot_p_rows,p_levelsda)  ! zonal mean of 3D press
REAL,ALLOCATABLE    :: pjin2d(:,:,:,:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PHOTIN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (first_call) THEN  
  myjppj = jppj
  CALL phot2d_allocate_memory (tot_P_ROWS,p_levelsda,ntphot)
  first_call = .FALSE.
END IF

ALLOCATE(pjin2d(nolat,nlphot,ntphot,myjppj))

! Find nearest day in photolysis dataset - day 1 is 31. December

idofy = i_day_number
fpos  = idofy/data_interval + 1.0
ipos  = NINT(fpos)
IF ((fpos-ipos*1.0) < 0.0) ipos = ipos-1

! Read in photolysis data

IF (L_optimised) THEN
  CALL read2d_opt(ipos,fpos,pjin2d)
ELSE
  CALL read2d_orig(ipos,fpos,pjin2d)
END IF

! Set up 2-D pressure arrays

CALL ukca_inpr2d(pr2d,pr2dj)

! Interpolate 2D photolysis rates onto 3-D levels and
! latitudes. Longitude comes later.

CALL ukca_calc_zmean_press(row_lengthda, tot_p_rows,           &
                         p_levelsda, glon, pl, zmean_pl)

CALL ukca_interpj(pjin2d,pr2dj,zmean_pl,row_lengthda,          &
             tot_p_rows, p_levelsda,                           &
             p_fieldda,sinlat,pi_over_180 )

DEALLOCATE(pjin2d)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_photin

!---------------------------------------------------------------------

SUBROUTINE ukca_curve(                                         &
               pjinda,tloc,dayl,p_field,p_levels,              &
               tot_p_rows,row_length,wks)
!
! Purpose: Subroutine to interpolate tropospheric photolysis rates
!          in time. Based on curve.F from Cambridge TOMCAT model.
!
!          Called from UKCA_CHEMISTRY_CTL.
!
! ---------------------------------------------------------------------
!
USE ereport_mod, ONLY: ereport
USE UM_ParVars
IMPLICIT NONE

INTEGER, INTENT(IN) :: p_field                   ! no of points
INTEGER, INTENT(IN) :: p_levels                  ! no of vert
INTEGER, INTENT(IN) :: tot_p_rows                ! no of rows
INTEGER, INTENT(IN) :: row_length                ! no of cols

REAL, INTENT(IN) :: dayl(:)           ! day length
REAL, INTENT(IN) :: tloc(:)           ! local time
REAL, INTENT(IN) :: pjinda(:,:,:)     ! 2D photolys

REAL, INTENT(OUT) :: wks(:,:)         ! interpolated

! Local variables

INTEGER :: i                                       ! loop variab
INTEGER :: j                                       ! loop variables
INTEGER :: k                                       ! loop variables
INTEGER :: jr                                      ! loop variables
INTEGER :: errcode                                 ! Variable passed to ereport

REAL, PARAMETER :: tfrac1 = 0.04691008             ! determines
REAL, PARAMETER :: tfrac2 = 0.23076534             ! 2D photolys

REAL :: dawn                ! time of dawn
REAL :: dusk                ! time of dusk
REAL :: timel               ! local time
REAL :: slope               ! slope used in linear interpolation
REAL :: const               ! intercept used in linear interpola
REAL :: fgmt(7)             ! times at which photol rates are va


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CURVE'


! Initialise wks

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
wks = 0.0

! Calculate rates using a simple linear interpolation.

DO j = 1,tot_p_rows
  DO i = 1,row_length
    k = i+(j-1)*row_length

    ! Non-Polar night

    IF (dayl(k) > 0.0) THEN
      dawn = 12.00 - (dayl(k)/2.0)
      dusk = 12.00 + (dayl(k)/2.0)

      fgmt(1) = dawn
      fgmt(2) = dawn + tfrac1*dayl(k)
      fgmt(3) = dawn + tfrac2*dayl(k)
      fgmt(4) = 12.00
      fgmt(5) = dawn + (1.0-tfrac2)*dayl(k)
      fgmt(6) = dawn + (1.0-tfrac1)*dayl(k)
      fgmt(7) = dusk

      timel = tloc(k)

      IF (timel > 24.0) timel = timel-24.0

      timel = MIN(timel,24.0)
      timel = MAX(timel, 0.0)

      ! Local Night-time

      IF (timel < dawn .OR. timel > dusk) THEN

        wks(k,1:MYjppj) = 0.0

        ! For the time between dawn and PJIN(1) or PJIN(5) and dusk

      ELSE IF ((timel >= dawn   .AND. timel < fgmt(2))         &
           .OR. (timel > fgmt(6) .AND. timel <= dusk)) THEN

        IF (timel > fgmt(6)) timel = 24.00 - timel

        ! trap for -ve (timel-fgmt(1))
        IF ((fgmt(1) - timel) < 1.0e-6) timel = fgmt(1)

        DO jr = 1, myjppj
          slope = pjinda(j,1,jr)/(fgmt(2) - fgmt(1))
          wks(k,jr) = slope*(timel - fgmt(1))

          IF (wks(k,jr) < 0.0) THEN
            WRITE(umMessage,*) 'negative wks in ukca_curve 1',wks(k,jr)
            CALL umPrint(umMessage,src='ukca_phot2d')
            WRITE(umMessage,*) i,j,k,jr,slope,fgmt(1),fgmt(2),timel
            CALL umPrint(umMessage,src='ukca_phot2d')
            WRITE(umMessage,*) pjinda(j,1,jr),fgmt(1)-timel
            CALL umPrint(umMessage,src='ukca_phot2d')
            WRITE(umMessage,*) fgmt,dawn,dusk
            CALL umPrint(umMessage,src='ukca_phot2d')
            errcode = jr
            CALL ereport('UKCA_CURVE',errcode,' Negative photolysis')
          END IF
        END DO

        ! For the time between PJIN(1) and PJIN(2) or PJIN(4) and PJIN(5)

      ELSE IF ((timel >= fgmt(2) .AND. timel < fgmt(3))        &
           .OR. (timel >  fgmt(5) .AND. timel <= fgmt(6))) THEN

        IF (timel > fgmt(5)) timel = 24.00 - timel

        DO jr = 1, myjppj
          slope = (pjinda(j,2,jr)-pjinda(j,1,jr))/             &
               (fgmt(3) - fgmt(2))
          const = pjinda(j,1,jr)- slope* fgmt(2)
          wks(k,jr)= slope*timel + const
        END DO

        ! For the time between PJIN(2), PJIN(3) and PJIN(4)

      ELSE IF (timel >= fgmt(3) .AND. timel <= fgmt(5)) THEN

        IF (timel > fgmt(4)) timel = 24.00 - timel

        DO jr = 1, myjppj
          slope = (pjinda(j,3,jr)-pjinda(j,2,jr))/             &
              (fgmt(4) - fgmt(3))
          const = pjinda(j,2,jr)- slope* fgmt(3)
          wks(k,jr)= slope*timel + const
        END DO

      END IF    ! end of IF (timel < dawn .OR. timel > dusk)

      ! End of the condition on the non polar night

    ELSE

      wks(k,1:myjppj) = 0.0

    END IF     ! end of IF (dayl(k) > 0.0)

  END DO       ! end of looping over row length
END DO         ! end of looping over rows

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_curve

! ----------------------------------------------------------------------

SUBROUTINE read2d_opt(ipos,fpos,pjin2d)

! ----------------------------------------------------------------------
! Purpose: Subroutine to read 2D photolysis rates.
!          Based on read2d.F from Cambridge TOMCAT model.
!
!          Called from UKCA_PHOTIN.
!
! ----------------------------------------------------------------------

USE ukca_option_mod, ONLY: jppj, phot2d_dir
USE ereport_mod, ONLY: ereport
USE UM_ParVars
IMPLICIT NONE



INTEGER, INTENT(INOUT) :: ipos     ! integer position
REAL,    INTENT(IN)    :: fpos     ! real position
REAL,    INTENT(OUT)   :: pjin2d(:,:,:,:)   ! Photol rate

! Local variables
INTEGER :: ifinx(myjppj)
INTEGER :: j                             ! Loop variable
INTEGER :: jr                            ! Loop variable
INTEGER :: k                             ! Loop variable
INTEGER :: kk                            ! Loop variable
INTEGER :: IN                            ! Index
INTEGER :: ierror                        ! Error flag
INTEGER :: info                          ! Tag for communication
INTEGER :: errcode                       ! Variable passed to ereport
INTEGER :: ukca2pho_unit                 ! File unit

REAL :: delpos
! Photol rates at all times
REAL, ALLOCATABLE, SAVE :: pjin2da(:,:,:,:,:)
REAL :: pr(myjppj,3)
REAL :: pfrac(1,3)

LOGICAL, SAVE :: L_first=.TRUE.        ! Logical for firstcall

CHARACTER(LEN=256):: file2             ! Dir for photol files
CHARACTER(LEN=256):: fnme(myjppj)      ! inc length of string
CHARACTER(LEN=10) :: csp(myjppj,jpspj)
CHARACTER(LEN=10) :: cmnt(myjppj)      ! Comment line
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message
CHARACTER(LEN=errormessagelength) :: iomessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ2D_OPT'


! 1.  Determine filenames containing specified photolysis rates


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
pjin2d = 0.0
! Reset ipos if nearest day=365 or zero
delpos = fpos - ipos*1.0
IF (ipos == 74) ipos=1

! 1.1  Read filenames from photolysis ratefile on first call

IF (L_first) THEN
  ALLOCATE(pjin2da(myjppj,n_data_times,nlphot,ntphot,nolat))

  !       use module to get cmnt

  IF (SIZE(ratj_defs) /= jppj) THEN
    cmessage='size of ratj_defs is not equal to jppj'
    errcode=1
    CALL ereport('UKCA_PHOT2D.UKCA_READ2D',errcode,cmessage)
  END IF
  DO k=1,jppj
    cmnt(k)=TRIM(ADJUSTL(ratj_defs(k)%fname))
  END DO

  !       1.2  Add '.bin' extension

  file2 = TRIM(phot2d_dir)//'/'
  DO jr = 1, MYjppj
    IN = INDEX(cmnt(jr),' ') - 1
    IF ( IN < 0 ) IN = 10
    fnme(jr) = TRIM(ADJUSTL(file2))//                             &
         TRIM(ADJUSTL(cmnt(jr)(1:IN)))//'.bin'
    IF (mype == 0 .AND. PrintStatus >= PrStatus_Diag) THEN
      WRITE(umMessage,'(3A)') 'fnme =', fnme(jr), cmnt(jr)
      CALL umPrint(umMessage,src='ukca_phot2d')
    END IF
  END DO

  ! 2.  Read specified photolysis rates

  IF (mype == 0) THEN
    DO jr = 1, MYjppj

      CALL assign_file_unit(fnme(jr),ukca2pho_unit,handler="fortran")
      OPEN(ukca2pho_unit, FILE=fnme(jr), FORM='UNFORMATTED',        &
           ACTION='READ', IOSTAT=ierror, IOMSG=iomessage)
      IF (ierror /= 0) THEN
        cmessage =                                         newline//&
        'Error opening file:'//                            newline//&
        TRIM(fnme(jr))//                                   newline//&
        'IoMsg: '//TRIM(iomessage)
        errcode = jr 
        CALL ereport('UKCA_PHOT2D.READ2D_OPT',errcode,cmessage)
      END IF

      READ(ukca2pho_unit) pjin2da(jr,:,:,:,:)
      CLOSE(ukca2pho_unit)
      CALL release_file_unit(ukca2pho_unit,handler="fortran")

    END DO   ! jr
  END IF  ! mype
  l_first=.FALSE.
END IF  ! l_first

! Interpolate in time using the saved values
IF (mype == 0) THEN
  DO jr = 1, MYjppj
    DO k = 1,nolat
      DO kk = 1,ntphot
        DO j = 1,nlphot
          pjin2d(k,j,kk,jr) =                                  &
                (pjin2da(jr,ipos+1,j,kk,k)-                    &
                 pjin2da(jr,ipos,j,kk,k))*delpos +             &
                 pjin2da(jr,ipos,j,kk,k)
        END DO
      END DO
    END DO
  END DO  ! jr
END IF   ! mype

CALL gc_rbcast (1,nlphot*ntphot*nolat*MYjppj,0,nproc             &
                 ,info,pjin2d)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read2d_opt

! ---------------------------------------------------------------------

SUBROUTINE read2d_orig(ipos,fpos,pjin2d)

! Purpose: Subroutine to read 2D photolysis rates.
!          Based on read2d.F from Cambridge TOMCAT model.
!
!          Called from UKCA_PHOTIN.
!
! ---------------------------------------------------------------------

USE ukca_option_mod, ONLY: jppj, phot2d_dir
USE ereport_mod, ONLY: ereport
USE UM_ParVars
IMPLICIT NONE



REAL, INTENT(IN) :: fpos

INTEGER, INTENT(INOUT) :: ipos

REAL, INTENT(OUT) :: pjin2d(nolat,nlphot,ntphot,myjppj) ! Photol rate

! Local variables

INTEGER :: ifinx(myjppj)
INTEGER :: ij                            ! Loop variable
INTEGER :: j                             ! Loop variable
INTEGER :: jr                            ! Loop variable
INTEGER :: k                             ! Loop variable
INTEGER :: kk                            ! Loop variable
INTEGER :: IN                            ! Index
INTEGER :: ierror                        ! Error flag
INTEGER :: info                          ! Tag for communication
INTEGER :: errcode                       ! Variable passed to ereport
INTEGER :: ukca2pho_unit                 ! File unit

REAL :: delpos
REAL :: pjin2d1(nlphot,ntphot,nolat)! Photol rates straddling time
REAL :: pjin2d2(nlphot,ntphot,nolat)! Photol rates straddling time
REAL :: pr(myjppj,3)
REAL :: pfrac(1,3)


CHARACTER(LEN=80) :: file2             ! Dir for photol files
CHARACTER(LEN=60) :: fnme(myjppj)      ! inc length of string
CHARACTER(LEN=10) :: csp(myjppj,jpspj)
CHARACTER(LEN=10) :: cmnt(myjppj)      ! Comment line
CHARACTER(LEN=errormessagelength) :: cmessage          ! Error message
CHARACTER(LEN=errormessagelength) :: iomessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ2D_ORIG'


! 1.  Determine filenames from module info

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (SIZE(ratj_defs) /= jppj) THEN
  cmessage='size of ratj_defs is not equal to jppj'
  errcode=1
  CALL ereport('UKCA_PHOT2D.UKCA_READ2D',errcode,cmessage)
END IF

DO k=1,jppj
  cmnt(k)=TRIM(ADJUSTL(ratj_defs(k)%fname))
END DO


! 1.2  Add '.d' extension

file2 = TRIM(phot2d_dir)//'/'

DO jr = 1, myjppj
  IN = INDEX(cmnt(jr),' ') - 1
  IF ( IN < 0 ) IN = 10
  fnme(jr) = TRIM(ADJUSTL(file2))//TRIM(ADJUSTL(cmnt(jr)(1:IN)))  &
       //'.dat'
  IF (mype == 0 .AND. PrintStatus >= PrStatus_Diag) THEN
    WRITE(umMessage,'(3A)') 'fnme =', fnme(jr), cmnt(jr)
    CALL umPrint(umMessage,src='ukca_phot2d')
  END IF
END DO

! Reset ipos if nearest day=365 or zero

delpos = fpos - ipos*1.0
IF (ipos == 74) ipos=1

! 2.  Read specified photolysis rates

pjin2d = 0.0
IF (mype == 0) THEN
  DO jr = 1, myjppj

    CALL assign_file_unit(fnme(jr), ukca2pho_unit, handler="fortran")
    OPEN(ukca2pho_unit, FILE=fnme(jr), FORM='FORMATTED',                &
         ACTION='READ', IOSTAT=ierror, IOMSG=iomessage)
    IF (ierror /= 0) THEN
      cmessage =                                               newline//&
        'Error reading file:'//                                newline//&
        TRIM(fnme(jr))//                                       newline//&
        'IoMsg: '//TRIM(iomessage)
      errcode = jr
      CALL ereport('UKCA_PHOT2D.READ2D_ORIG',errcode,cmessage)
    END IF

    DO ij = 1,ipos
      READ(ukca2pho_unit,*) (((pjin2d1(j,kk,k),j=1,51),kk=1,3)          &
                     ,k=nolat,1,-1)
    END DO

    READ(ukca2pho_unit,*) (((pjin2d2(j,kk,k),j=1,51),kk=1,3)            &
                  ,k=nolat,1,-1)

    DO k = 1,nolat
      DO kk = 1,ntphot
        DO j = 1,nlphot
          pjin2d(k,j,kk,jr) =                                    &
                   (pjin2d2(j,kk,k)-pjin2d1(j,kk,k))*            &
                    delpos + pjin2d1(j,kk,k)
        END DO
      END DO
    END DO

    CLOSE(ukca2pho_unit)
    CALL release_file_unit(ukca2pho_unit, handler="fortran")

  END DO
END IF      ! End of IF mype statement

CALL gc_rbcast (1,nlphot*ntphot*nolat*myjppj,0,nproc             &
                 ,info,pjin2d)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read2d_orig

!------------------------------------------------------------

SUBROUTINE ukca_inpr2d(pr2d,pr2dj)
!
! Purpose: Subroutine to calculate the pressure levels for the
!          2D photolysis rates. Original version taken from the
!          Cambridge TOMCAT model.
!
!          Called from UKCA_PHOTIN.
!
! ---------------------------------------------------------------------
!
USE ereport_mod, ONLY: ereport
USE UM_ParVars
IMPLICIT NONE

!       Local variables

INTEGER, PARAMETER :: maxlev = 30

INTEGER :: j                        ! Loop variable
INTEGER :: ij                       ! Loop variable

REAL, PARAMETER :: fac  = 1.284025417  ! Used to ensure that pre
REAL, PARAMETER :: psur = 1.0e5        ! Surface pressure in Pas
REAL, PARAMETER :: ares = 3.0          ! Factor related to verti
REAL, PARAMETER :: eps  = 1.0e-10      ! Factor used to calculat

REAL :: ee                       ! Temporary store
REAL :: fj                       ! Factor related to vertical re
REAL :: za                       ! Factor related to vertical re
REAL :: pr2d(nolev)              ! Pressures of 2D model levels
REAL :: pr2dj(nlphot)            ! Pressures of 2D photolysis le
REAL :: pp(nlphot+1)
REAL :: pres(maxlev)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INPR2D'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO j = 1, nolev
  ee = EXP((j-1)/2.0)
  pr2d(j) = psur/(ee * fac)
END DO

!       2D pressure levels - normal 2-D order - up to level 30
!       for photolysis

DO j = 1,maxlev-1
  ee = EXP((j-1) / 2.0)
  pres(j) = psur / (ee * fac)
END DO
pres(maxlev) = pres(maxlev-1)

fj    = 2.0/ares
pp(1) = (1.0-fj)*LOG(psur)+fj*LOG(pres(1))
pp(1) = EXP(pp(1))

DO ij = 2,nlphot+1
  za     = ij/(2.0*ares)
  fj     = 2.0*za+0.5+eps
  j      = INT(fj)
  fj     = fj-j-eps
  j      = j+1
  pp(ij) = (1.0-fj)*LOG(pres(j-1))+ fj*LOG(pres(j))
  pp(ij) = EXP(pp(ij))
END DO

pr2dj(1)=(psur+pp(1))*0.5

DO ij = 2,nlphot
  pr2dj(ij) = (pp(ij)+pp(ij-1))*0.5
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_inpr2d

!----------------------------------------------------------------

SUBROUTINE ukca_interpj(pjin2d,pr2dj,zm_pl,lon,lat,lev,         &
                        p_field,sinlat,degrad)
!
! Purpose: Subroutine to interpolate photolysis rates
!
!          Called from UKCA_PHOTIN.
!
! ---------------------------------------------------------------------
USE ereport_mod, ONLY: ereport
USE UM_ParVars
IMPLICIT NONE

INTEGER :: lon                               ! No of longitudes
INTEGER :: lat                               ! No of latitudes
INTEGER :: lev                               ! No of levels
INTEGER :: p_field                           ! No of spatial poi

REAL, INTENT(IN)       :: pjin2d(:,:,:,:)  ! 2D photo
REAL, INTENT(IN)       :: zm_pl(:,:)   ! zmean 3D press
REAL, INTENT(IN)       :: pr2dj(:)   ! 2D photo
REAL, INTENT(IN)       :: sinlat(:)  ! Sine (3D
REAL, INTENT(IN)                     :: degrad  ! To conve

!       Local variables

INTEGER :: i                     ! Loop variable
INTEGER :: j                     ! Loop variable
INTEGER :: jr                    ! Loop variable
INTEGER :: k                     ! Loop variable
INTEGER :: kk                    ! Loop variable
INTEGER :: l                     ! Loop variable

REAL, PARAMETER :: Npole = 90.0

REAL :: p2d(nlphot)
REAL :: lat2d(nolat)          ! 2D model latitudes
REAL :: wks1(lat,nlphot)      ! Working array
REAL :: wks2(nolat)           ! Working array
REAL :: wks3(nlphot)          ! Working array
REAL :: wks4(lat,lev)         ! Working array
REAL :: lati
REAL :: press
REAL :: delphi

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INTERPJ'


!       Set up 2D latitudes. lat=1 pole nord
!       LAT2D() is the latitude in the centre of the 2D box in radians

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
delphi = Npole/nolat
DO i = 1,nolat
  lat2d(i)=SIN((90.0 - (2*i-1)*delphi)*degrad)
END DO

DO jr = 1,MYjppj              ! Loop over photolysis reactions

  !         Interpolate linearly in Sin(lat) (KK is the point through the

  DO kk = 1,ntphot          ! Loop over times of day

    DO j = 1,nlphot
      DO i = 1,nolat
        wks2(i) = pjin2d(i,j,kk,jr)
      END DO

      DO k = 1,lat
        lati = sinlat((k-1)*lon+1)
        lati = MAX(lati, lat2d(nolat))
        lati = MIN(lati, lat2d(    1))
        wks1(k,j) = ukca_flupj(lati,lat2d,wks2,nolat)
        wks1(k,j) = MAX(wks1(k,j), 0.0)
      END DO
    END DO

    !           Interpolate linearly in log(P)

    DO k = 1,lat
      DO j = 1,nlphot
        wks3(j) = wks1(k,j)
        p2d(j)  = LOG(pr2dj(j))
      END DO

      DO l = 1,lev
        press = LOG(zm_pl(k,l))
        press = MAX(press, p2d(nlphot))
        press = MIN(press, p2d(1))
        wks4(k,l) = ukca_flupj(press,p2d,wks3,nlphot)
        wks4(k,l) = MAX(wks4(k,l), 0.0)
      END DO
    END DO

    DO l = 1,lev
      DO k = 1,lat
        pjin(k,l,kk,jr) = wks4(k,l)
      END DO
    END DO


  END DO     ! End of loop over times of day
END DO       ! End of loop over photolysis reactions

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_interpj

!-----------------------------------------------------------------

SUBROUTINE phot2d_allocate_memory(lat,lev,ntphot)
!
! Purpose: Subroutine to allocate array to hold
!          2D photolysis rates
!
!          Called from UKCA_PHOTIN.
!
! ---------------------------------------------------------------------
USE ereport_mod, ONLY: ereport
USE UM_ParVars
IMPLICIT NONE
INTEGER,INTENT(IN)         :: lat,lev,ntphot

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHOT2D_ALLOCATE_MEMORY'


!       lat is the number of latitude on current PE, NOT nolat!
!       lev is numer of model level, NOT nolev

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
ALLOCATE(pjin(lat,lev,ntphot,MYjppj))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE phot2d_allocate_memory

!----------------------------------------------------------------------

SUBROUTINE ukca_calc_zmean_press(lon, lat, lev,              &
                                 glon,pl, zmean_pl)
!
! Purpose: Subroutine to calculate zonal mean pressure profile
!
!          Called from UKCA_PHOTIN.
!
! ---------------------------------------------------------------------
USE global_2d_sums_mod, ONLY: global_2d_sums
USE ereport_mod, ONLY: ereport
USE UM_ParVars

IMPLICIT NONE


INTEGER, INTENT(IN) :: lon        ! No of longitudes
INTEGER, INTENT(IN) :: lat        ! No of latitudes
INTEGER, INTENT(IN) :: lev        ! No of levels
INTEGER, INTENT(IN) :: glon       ! No of global lons

REAL, INTENT(IN)    :: pl(lon,lat,lev)   ! Pressure

REAL, INTENT(OUT)   :: zmean_pl (lat,lev)

!     Local variables

INTEGER :: j,k                 ! Loop variables
INTEGER :: istat               ! Output status from GCOM routine

REAL    :: fac                 ! Factor used to calculate zmean
REAL    :: sumpl_row(lat,lev)  ! Global sum along rows

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALC_ZMEAN_PRESS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
zmean_pl = 0.0

CALL global_2d_sums(pl, lon, 1, 0, 0, lev*lat,                    &
                    sumpl_row, gc_proc_row_group)

fac      = 1.0/REAL(glon)
zmean_pl = sumpl_row*fac

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calc_zmean_press

!----------------------------------------------------------------------

END MODULE ukca_phot2d

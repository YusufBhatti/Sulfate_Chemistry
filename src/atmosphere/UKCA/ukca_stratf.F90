! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Purpose: Subroutine to overwrite values at top of model
!          using interpolated 5-day fields from the 2-d model.
!          Based on STRATF.F from Cambridge TOMCAT model and
!          modified by Olaf Morgenstern to allow for flexible
!          positioning of the NOy species.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_CHEMISTRY_CTL.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_stratf_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_STRATF_MOD'

CONTAINS

SUBROUTINE ukca_stratf(i_day_number, row_length, rows,            &
                       model_levels, theta_field_size,            &
                       first_row, global_row_length,              &
                       global_rows, ntracer, sinlat,              &
                       pres, um_ozone3d, p_above, tracer)

USE asad_mod,          ONLY: advt
USE ukca_tropopause,   ONLY: tropopause_level
USE ukca_cspecies,     ONLY: n_o3, n_o3s, n_hono2, n_ch4, n_no,   &
                             n_no2, n_no3, n_n2o5, n_ho2no2
USE ukca_constants,    ONLY: c_o3, c_n, c_ch4, c_hno3, c_no,      &
                             c_no3, c_no2, c_n2o5, c_hono2,       &
                             c_ho2no2
USE ukca_option_mod,   ONLY: L_ukca_use_2dtop, strat2d_dir
USE parkind1,          ONLY: jprb, jpim
USE yomhook,           ONLY: lhook, dr_hook
USE UM_ParVars
USE UM_ParCore,        ONLY: mype, nproc
USE ereport_mod,       ONLY: ereport
USE umPrintMgr
USE param2d_mod, ONLY: nolat,nolev,nlphot,ntphot,jpphio,jpphin

USE errormessagelength_mod, ONLY: errormessagelength

USE ukca_2d_bc_read_interp_mod, ONLY: ukca_2d_bc_read_interp
USE ukca_calc_noy_zmeans_mod, ONLY: ukca_calc_noy_zmeans
USE ukca_interp_mod, ONLY: ukca_interp
IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length         ! No of longitudes
INTEGER, INTENT(IN) :: rows               ! No of latitudes
INTEGER, INTENT(IN) :: model_levels       ! No of levels
INTEGER, INTENT(IN) :: theta_field_size   ! No of spatial points
INTEGER, INTENT(IN) :: first_row          ! First global latitude
INTEGER, INTENT(IN) :: global_row_length  ! No of global longitudes
INTEGER, INTENT(IN) :: global_rows        ! No of global latitudes
INTEGER, INTENT(IN) :: ntracer            ! No of chemical tracers

REAL, INTENT(IN) :: sinlat(theta_field_size)           ! Sine(latitude)
REAL, INTENT(IN) :: pres(row_length,rows,model_levels) ! Model pressures
REAL, INTENT(IN) :: um_ozone3d(row_length,rows,model_levels) ! UM O3
REAL, INTENT(IN) :: p_above               ! above which tbc are applied

! Tracer concentrations in mass mixing ratio
REAL, INTENT(INOUT) :: tracer(row_length,rows,model_levels,ntracer)

!     Local variables

INTEGER :: i_day_number   ! Day number
INTEGER :: info           ! Tag used in communication
INTEGER :: i,ij,j,k,l     ! Loop variables
INTEGER :: m              ! Loop variable

LOGICAL :: mask(row_length,rows,model_levels) ! mask to identify stratosphere

REAL, PARAMETER :: o3_hno3_ratio = 1.0/1000.0 ! kg[N]/kg[O3] from
                                              ! Murphy and Fahey 1994

REAL :: noy(nolat,nolev)                      ! 2D field interpolated in time
REAL :: o3(nolat,nolev)                       ! 2D field interpolated in time
REAL :: ch4(nolat,nolev)                      ! 2D field interpolated in time
REAL :: o33d(row_length,rows,model_levels)    ! 2D field interpolated onto 3D
REAL :: noy3d(row_length,rows,model_levels)   ! 2D field interpolated onto 3D
REAL :: ch43d(row_length,rows,model_levels)   ! 2D field interpolated onto 3D
REAL :: hno33d(row_length,rows,model_levels)  ! 3D field from fixed o3:hno3 ratio
REAL :: noyzm(rows,model_levels)              ! NOy zonal mean
REAL :: hno3zm(rows,model_levels)             ! HNO3 zonal mean
REAL :: hno4zm(rows,model_levels)             ! HNO4 zonal mean
REAL :: nozm(rows,model_levels)               ! NOz zonal mean
REAL :: no2zm(rows,model_levels)              ! NO2 zonal mean
REAL :: no3zm(rows,model_levels)              ! NO3 zonal mean
REAL :: n2o5zm(rows,model_levels)             ! N2O5 zonal mean

! Parameter to use UM ancil ozone instead of 2-D O3
LOGICAL, PARAMETER :: L_use_UMO3 = .TRUE.

! Parameter to use fixed o3 to hno3 ratio rather than 2d NOy
LOGICAL, PARAMETER :: L_use_O3HNO3ratio = .TRUE.

! Parameter to overwrite stratosphere (fixed no of levels above tropopause)
LOGICAL, PARAMETER :: L_all_strat   = .TRUE.
INTEGER, PARAMETER :: no_above_trop1 = 3        ! Suitable for L38/L60
INTEGER, PARAMETER :: no_above_trop2 = 10       ! Suitable for L63/L85
INTEGER :: no_above_trop

! Parameter to overwrite CH4 with 2D boundary conditions
LOGICAL, PARAMETER :: L_overwrite_CH4 = .FALSE.

INTEGER           :: errcode                    ! Error code: ereport
CHARACTER(LEN=errormessagelength) :: cmessage                   ! Error message
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_STRATF'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set number of levels above tropopause depending on vertical resolution
IF (model_levels == 38 .OR. model_levels == 60) THEN
  no_above_trop = no_above_trop1
ELSE IF (model_levels == 63 .OR. model_levels == 70 .OR.          &
         model_levels == 85) THEN
  no_above_trop = no_above_trop2
ELSE
  errcode = 1
  cmessage = 'Levels above tropopause not set at this resolution'
  CALL ereport('UKCA_STRATF',errcode,cmessage)
END IF

IF (printstatus == Prstatus_Diag) THEN
  WRITE(umMessage,*) 'UKCA_STRATF: Logicals in use:'
  CALL umPrint(umMessage,src='ukca_stratf')
  WRITE(umMessage,*) 'L_overwrite_CH4: ',L_overwrite_CH4
  CALL umPrint(umMessage,src='ukca_stratf')
  WRITE(umMessage,*) 'L_all_strat: ',L_all_strat
  CALL umPrint(umMessage,src='ukca_stratf')
  WRITE(umMessage,*) 'L_use_O3HNO3ratio: ',L_use_O3HNO3ratio
  CALL umPrint(umMessage,src='ukca_stratf')
  WRITE(umMessage,*) 'L_use_UMO3: ',L_use_UMO3
  CALL umPrint(umMessage,src='ukca_stratf')
  WRITE(umMessage,*) 'L_ukca_use_2dtop: ',L_ukca_use_2dtop
  CALL umPrint(umMessage,src='ukca_stratf')
  WRITE(umMessage,*) 'no_above_trop: ',no_above_trop
  CALL umPrint(umMessage,src='ukca_stratf')
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='ukca_stratf')
END IF

IF (L_ukca_use_2dtop) THEN
  ! Using data from Cambridge 2D model
  ! Get 2D data for current time then  interpolate in space

  CALL ukca_2d_bc_read_interp(mype,nproc,i_day_number,strat2d_dir, &
                              noy,o3,ch4)
END IF    ! End of IF(L_ukca_use_2dtop) statement

! Interpolate onto 3-D lats and levs and convert to mass mixing ratio (O3,CH4)
!  or calculate as VMR.
IF (L_use_UMO3) THEN
  o33d(:,:,:) = um_ozone3d(:,:,:)
ELSE
  CALL ukca_interp(o3, pres, row_length, rows, model_levels,      &
                   theta_field_size, sinlat, o33d)
  o33d(:,:,:) = o33d(:,:,:)*c_o3
END IF

IF (L_use_O3HNO3ratio) THEN
  hno33d(:,:,:) = o33d(:,:,:)*o3_hno3_ratio*c_hno3/c_n
ELSE
  CALL ukca_interp(noy, pres, row_length, rows, model_levels,     &
                   theta_field_size, sinlat, noy3d)
END IF

IF (L_overwrite_CH4 ) THEN
  CALL ukca_interp(ch4, pres, row_length, rows, model_levels,     &
                   theta_field_size, sinlat, ch43d)
  ch43d(:,:,:) = ch43d(:,:,:)*c_ch4
END IF

IF (.NOT. L_use_O3HNO3ratio) THEN
  ! Calculate zonal means for NOy species
  CALL ukca_calc_noy_zmeans(row_length, rows, model_levels,      &
                     global_row_length, global_rows, first_row,  &
                     tracer, noyzm, nozm, no2zm, no3zm,          &
                     n2o5zm, hno3zm, hno4zm)
END IF

! Overwrite values in tracer array

!       New N fields = zonal 3D field * (2d NOy/3d NOy) - at present only
!       scale to total NOy - the partitioning between NOy species is
!       not really correct because the species are treated differently
!       in the two models. This can be improved. --- Comments from TOMCAT

IF (L_use_O3HNO3ratio .AND. (.NOT. L_all_strat)) THEN

  ! Overwrite o3, hno3, and ch4 vmr at all gridboxes above p_above
  !  O3   - either UM ancillary or Cambridge 2d
  !  HNO3 - using fixed o3:hno3 ratio
  !  CH4  - from Cambridge 2D model

  mask(:,:,:) = .FALSE.
  WHERE (pres <= p_above) mask = .TRUE.

  WHERE (mask(:,:,:))
    tracer(:,:,:,n_o3)    = o33d(:,:,:)
    tracer(:,:,:,n_hono2) = hno33d(:,:,:)
  END WHERE
  IF (n_o3s > 0) THEN
    WHERE (mask(:,:,:))
      tracer(:,:,:,n_o3s) = o33d(:,:,:)
    END WHERE
  END IF

  IF (L_overwrite_ch4) THEN
    WHERE (mask(:,:,:))
      tracer(:,:,:,n_ch4) = ch43d(:,:,:)
    END WHERE
  END IF

ELSE IF (L_use_O3HNO3ratio .AND. L_all_strat) THEN

  ! Overwrite o3, hno3, and ch4 at all gridboxes a fixed
  !  number of model levels above the tropopause
  !  O3   - either UM ancillary or Cambridge 2d
  !  HNO3 - using fixed o3:hno3 ratio
  !  CH4  - from Cambridge 2D model

  mask(:,:,:) = .FALSE.
  DO l = model_levels,1,-1
    mask(:,:,l) = tropopause_level(:,:)+no_above_trop <= l
  END DO

  WHERE (mask(:,:,:))
    tracer(:,:,:,n_o3)    = o33d(:,:,:)
    tracer(:,:,:,n_hono2) = hno33d(:,:,:)
  END WHERE
  IF (n_o3s > 0) THEN
    WHERE (mask(:,:,:))
      tracer(:,:,:,n_o3s) = o33d(:,:,:)
    END WHERE
  END IF

  IF (L_overwrite_ch4) THEN
    WHERE (mask(:,:,:))
      tracer(:,:,:,n_ch4) = ch43d(:,:,:)
    END WHERE
  END IF
ELSE

  ! Overwrite O3, NOy and CH4 at all gridboxes above p_above
  !  O3  - use either UM ancillary or Cambridge 2D
  !  NOy - use Cambridge 2D
  !  CH4 - use Cambridge 2D (if L_overwrite_ch4)

  mask(:,:,:) = .FALSE.
  WHERE (pres <= p_above) mask = .TRUE.

  DO k = 1,ntracer
    SELECT CASE (advt(k))
    CASE ('O3        ')
      WHERE (mask) tracer(:,:,:,n_o3) = o33d(:,:,:)
    CASE ('O3S       ')
      WHERE (mask) tracer(:,:,:,n_o3s) = o33d(:,:,:)
    CASE ('CH4       ')
      IF (L_overwrite_ch4) THEN
        WHERE (mask) tracer(:,:,:,n_ch4) = ch43d(:,:,:)
      END IF
    CASE ('NO        ')
      WHERE (mask) tracer(:,:,:,n_no) = noy3d(:,:,:)*            &
          SPREAD(nozm/noyzm,dim=1,ncopies=row_length)*c_no
    CASE ('NO3       ')
      WHERE (mask) tracer(:,:,:,n_no3) = noy3d(:,:,:)*           &
          SPREAD(no3zm/noyzm,dim=1,ncopies=row_length)*c_no3
    CASE ('NO2       ')
      WHERE (mask) tracer(:,:,:,n_no2) = noy3d(:,:,:)*           &
          SPREAD(no2zm/noyzm,dim=1,ncopies=row_length)*c_no2
    CASE ('N2O5      ')
      WHERE (mask) tracer(:,:,:,n_n2o5) = noy3d(:,:,:)*          &
          SPREAD(n2o5zm/noyzm,dim=1,ncopies=row_length)*c_n2o5
    CASE ('HO2NO2    ')
      WHERE (mask) tracer(:,:,:,n_ho2no2) = noy3d(:,:,:)*        &
          SPREAD(hno4zm/noyzm,dim=1,ncopies=row_length)*c_ho2no2
    CASE ('HONO2     ')
      WHERE (mask) tracer(:,:,:,n_hono2) = noy3d(:,:,:)*         &
          SPREAD(hno3zm/noyzm,dim=1,ncopies=row_length)*c_hono2
    END SELECT
  END DO

END IF    ! End of IF(L_use_O3HNO3ratio and L_all_strat) statement

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_stratf
END MODULE ukca_stratf_mod

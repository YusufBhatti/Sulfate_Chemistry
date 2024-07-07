! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Classify columns as convective or stratiform based on surface precip

MODULE crmstyle_class_col_mod

IMPLICIT NONE
! ---------------------------------------------------------------------------
! Description:
! 
! This subroutine determines the column classification for convective
! and stratiform regions. (Code taken from Met Office LEM model.)
! Also additional partitions of trailing  anvil regions, shallow clouds
! and ice anvils.
!
! The convective and stratiform partition that is coded by default is
! based on the technique of Steiner et al (1995,J.Appl.Meteor.,1978-2007)
! This is based on 2km radar data and therefore the precipitation data should
! be mapped onto a 2km grid.
! The criteria are (P=surface preciptation rate):
!    - Intensity: Any grid point with P >= 10 mm/h is classified as CONVECTIVE
!    - Peakedness: Any grid point with P that exceeds the background rate
!        P_b by at least a factor of 2 is CONVECTIVE. P_b is the average of
!        the rates in the rainy (P > 0.01mm/h) grid points within an 11km
!        radius.
!    - Surrounding Area: For each grid point identified as a convective centre
!        by at least one of the above criteria, all surrounding grid points
!        within an intensity dependent "convective radius" about the grid
!        point are defined as CONVECTIVE. This is based on the background
!        precipitation rate, P_b.
!               P_b         Convective Radius
!             < 1.2 mm/h          1 km
!           1.2 -> 2.5 mm/h       2 km
!           2.5 -> 5.0 mm/h       3 km
!           5.0 -> 10.0 mm/h      4 km
!              >= 10.0 mm/h       5 km
!
!    - Stratiform Column: All rainy points not convective are defined as
!        STRATIFORM.
! Note the code is setup treat input resolutions between 1-2km as if they are
! 2km data. This will allow the 1.5km data to be processed.
!
! The other column classifications come from EUCREM (?) and are based on
! the following. For all points which are neither CONVECTIVE or STRATIFORM
!    If IWP > 1.0 kg/m2 the point is in a TRAILING REGION
!    Elseif IWP > 0.01 kg/m2 the point is in an ICE ANVIL
!    Elseif LWP > 0.02 kg/m2 the point is in a SHALLOW CLOUD
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 95.
!   This code is written to UMDP 3 standards
! ---------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_CLASS_COL_MOD'

CONTAINS

SUBROUTINE crmstyle_class_col( )

USE crmstyle_cntl_mod, ONLY:                                             &
   in_cols,in_rows, new_res, mlevs, l_qgraup

USE crmstyle_grid_info_mod, ONLY:                                       &
   nprocs, local_row_len, local_rows, local_new_x, local_new_y

USE hires_data_mod, ONLY:                                                &
  iclass_col, precip, density, qcl, qcf, qrain, qgraup 

USE crmstyle_sample_arrays_mod, ONLY:                                    &
  fract_conv, fract_strat, prec_conv, prec_strat

USE crmwork_arrays_mod, ONLY:                                            &
  prec_full, mask, h_theta_sea, dz_theta_sea

USE crmstyle_pp_data_mod, ONLY:                                          &
  bdy, bdx

USE word_sizes_mod, ONLY: iwp,wp    ! Allows use of 4 byte words 

! UM parallel info and grids
USE UM_ParVars, ONLY: gc_all_proc_group
USE UM_ParParams, ONLY: halo_type_no_halo, halo_type_single, fld_type_p 
USE UM_ParCore, ONLY: mype

USE umPrintMgr                              ! for writing output

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments  - None
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! local variables
!----------------------------------------------------------------------

INTEGER ::             &
  i,j,k,ii,jj,ij       & ! loop counters
 , iii, jjj            & ! loop counters
 , nx, ny, ix, iy      & ! more loop counters
 , ncon, n             & !
 , ierror              & ! Error flag if grid spacing > radgrid
 , nxpts, nypts        & ! No. of model grid points in radar grid box
 , nradxpts, nradypts  & ! No. of radar grid boxes in domain.
 , nxrange             & ! guess of no. of radar grid points
 , nyrange             & ! guess of no. of radar grid points
                         ! within background radius
 , icode               & ! return code
 , icount                ! No. of rainy radar gridpoints in background rad.

INTEGER, PARAMETER ::  &
  nconradp = 5           ! number of thresholds  and distances

INTEGER(iwp), ALLOCATABLE :: &
  iclasstemp(:,:)              ! Temporary array for provisional 
                               ! classifications

LOGICAL ::             &
  l_notfit             & ! local logical true if grid-spacing not
                         ! integer multiple of RADGRID
 ,l_toosmall             ! local logical true if domain size LT backrad

REAL ::                                      &
  iclass_full(in_cols,in_rows)               & ! Column classification full 
                                               ! grid 
 ,field_out_local(local_row_len,local_rows)  & ! work array
 ,iwpath(local_row_len,local_rows)           & ! ice water path
 ,lwpath(local_row_len,local_rows)             ! liquid water path

REAL(wp), ALLOCATABLE  :: &
  pr_b(:,:)               & ! Array Containing The Background Precip. Rate 
, pr_rad(:,:)               ! Array Containing precip. rate on radar grid 

REAL ::                   &
  dist                    & ! Distance from a point (changes in code)
 ,radius                  & ! A tempory store of the convective radius
 ,rnpts                   & ! real number of points in xy plane
 ,ryrange                 & ! real no. of radar points within radius
 ,qicetemp                & ! sum of q ice components
 ,qltemp                  & ! sum of q liquid components
 ,nvalues                 & ! number of values
 ,rnvalues                  ! 1/ nvalues

REAL, PARAMETER ::        &
  radgrid_deg = 0.018     & ! 2km approximate resolution in degrees. 
 ,backrad_deg = 0.099     & ! Background Radius in degrees  (was originally
                            ! 11000. m)
 ,deg_1km   = 0.009         ! 1km resolution in degrees

REAL, PARAMETER ::        &
  con_intense  = 10.0     & ! ppt rate for Convective Point (Intensity)
 ,con_peak     =  2.0     & ! ppt rate > than background (Peakednesss)
 ,rainy_val    =  0.01    & ! ppt rate for rainy gridpoint (mm/hr)
 ,trail_iwp    =  1.0     & ! IWP value for trailing region (kg/m2)
 ,anvil_iwp    =  0.01    & ! IWP value for ice anvil (kg/m2)
 ,shallow_lwp  =  0.02      ! LWP value for shallow clouds (kg/m2)

! conrad_val contains the convective radius values & conrad_thres the 
! precipitation thresholds for each radius value
! Note the original radii were in m but these have been changed to degrees.
!  conrad_val(nconradp)   = (/1000., 2000., 3000., 4000., 5000./) 

REAL                                                                 &
! The following variables cannot be parameters when compiling with gfortran
! due to a compiler bug affecting parameter arrays in OpenMP regions.
#if !defined(GNU_FORTRAN)
, PARAMETER                                                          &
#endif
 ::  conrad_deg(nconradp)   = (/0.009, 0.018, 0.027, 0.036, 0.045/)  & ! degrees
    ,conrad_thres(nconradp) = (/1.2,   2.5,   5.0,   10.0,  1.0e5/)    ! mm/hr  


CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_CLASS_COL"
CHARACTER(LEN=errormessagelength) :: cmessage

! Required by Dr hook 
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

nypts = 0
nxpts = 0
l_notfit = .FALSE.

! Work out whether input data needs to be regridded to 2km or whether
! resolution close enough to 2km so use as is.
! At present assuming dx and dy same? 
! Error flags if grid spacing greater than RADGRID or
! if grid-spacing not integer multiple of RADGRID (Used MOD)
! or the domain is less than backrad

IF (bdy == radgrid_deg) THEN     ! 2 km input data
  nypts = 1
ELSE IF (bdy > deg_1km .AND. bdy < radgrid_deg) THEN
  nypts = 1
ELSE IF (bdy <= deg_1km) THEN
  ! Hope exact divisor of 2km specified in degrees
  nypts = NINT(radgrid_deg/bdy)
  l_notfit = MOD(radgrid_deg, bdy) /= 0.0 
ELSE
  l_notfit = .TRUE.      ! Data not suitable as input grid too large
END IF

IF (bdx == radgrid_deg) THEN     ! 2 km input data
  nxpts = 1
ELSE IF (bdx > deg_1km .AND. bdx < radgrid_deg) THEN
  nxpts = 1
ELSE IF (bdx <= deg_1km) THEN
  ! Hope exact divisor of 2km specified in degrees
  nxpts = NINT(radgrid_deg/bdx)
  l_notfit = MOD(radgrid_deg, bdx) /= 0.0 
ELSE
  l_notfit = .TRUE.      ! Data not suitable as input grid too large
END IF

l_toosmall = ( bdx*in_cols < backrad_deg)  .OR. (bdy*in_rows < backrad_deg)

! Exit routine if there is an error 
IF ((nypts == 0 .OR. nxpts == 0 ) .OR. l_notfit .OR. l_toosmall) THEN
  ierror = 1 
ELSE 
  ierror = 0
END IF
IF (ierror /= 0) THEN 
  WRITE(umMessage,'(A)')         &
        '***WARNING:COLUMN CLASSIFCIATION NOT VALID FOR THIS RUN***'
  CALL umPrint(umMessage)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! number of points for a 2km grid

IF (nypts /= 1) THEN
  nradypts = NINT((bdy*in_rows)/radgrid_deg)
ELSE
  nradypts = in_rows
END IF
IF (nypts /= 1) THEN
  nradxpts = NINT((bdx*in_cols)/radgrid_deg)
ELSE
  nradxpts = in_cols
END IF

! Only have the full precipitation array on mype 0

IF (mype == 0) THEN

  ! Initialise arrays

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                                &
!$OMP& SHARED(in_rows, in_cols , iclass_full, prec_full)
  DO j=1,in_rows
    DO i=1,in_cols
      iclass_full(i,j) = 0.0 
      ! Convert surface precipitation from mm/s to mm/hour for use here
      prec_full(i,j) = prec_full(i,j)*3600.0
    END DO
  END DO
!$OMP END PARALLEL DO

    ! Allocate work arrays for 2km grid 

  ALLOCATE(iclasstemp(nradxpts,nradypts))
  ALLOCATE(pr_b(nradxpts,nradypts))
  ALLOCATE(pr_rad(nradxpts,nradypts))

  ! regrid to 2km
  ! Loop over all points required in PR_RAD

  rnpts = 1.0/(REAL(nxpts*nypts))

!$OMP PARALLEL DO PRIVATE(ny,nx,iy,ix, i,j) DEFAULT(NONE)                  &
!$OMP& SHARED(nradypts, nradxpts, nypts, nxpts, pr_rad, pr_b, prec_full,   &
!$OMP& rnpts)
  DO ny = 1,nradypts
    DO nx = 1,nradxpts
      pr_rad(nx,ny) = 0.0
      pr_b(nx,ny)   = 0.0
      ! For each point loop over each set of points for RADGRID square
      DO iy = 1,nypts
        DO ix = 1,nxpts
          i = (nx-1)*nxpts + ix
          j = (ny-1)*nypts + iy
          pr_rad(nx,ny) = pr_rad(nx,ny) + prec_full(i,j)
        END DO ! ix
      END DO   ! iy
      ! Form mean precip for new grid
      pr_rad(nx,ny) = pr_rad(nx,ny)*rnpts
    END DO     ! nx
  END DO       ! ny
!$OMP END PARALLEL DO

    ! Precipitation data now averaged onto the required resolution
    ! First of all work out the precipitation background field.
    ! Select the range over which to search for points within
    ! the  background radius

  nyrange = INT(backrad_deg/radgrid_deg) + 1
  nxrange = nyrange

!$OMP PARALLEL DO PRIVATE(ny,nx,iy,ix, i,j, dist, icount) DEFAULT(NONE)    &
!$OMP& SHARED(nradypts, nradxpts, nypts, nxpts, nyrange, nxrange, pr_rad,  &
!$OMP& pr_b, prec_full, rnpts, iclasstemp)
  DO ny = 1,nradypts
    DO nx = 1,nradxpts
      icount = 0
      DO iy = ny - nyrange, ny + nyrange
        DO ix = nx - nxrange,nx + nxrange

          ! Check that point is within BACKRAD of the centre

          dist = SQRT( ((nx - ix)*radgrid_deg)**2 +       &
                          ((ny - iy)*radgrid_deg)**2 ) 
          IF (dist < backrad_deg .AND. dist > 0.0 ) THEN

            ! Do wrap around on IX and IY

            i = ix
            IF (i < 1)        i = i + nradxpts
            IF (i > nradxpts) i = i - nradxpts
            j = iy
            IF (j < 1)        j = j + nradypts
            IF (j > nradypts) j = j - nradypts

            ! Calc. the background PPN (points must be "RAINY")
 
            IF (pr_rad(i,j) > rainy_val) THEN
              icount = icount + 1
              pr_b(nx,ny) = pr_b(nx,ny) + pr_rad(i,j)
            END IF
          END IF ! dist < backrad .AND. dist > 0.0
        END DO   ! ix
      END DO     ! iy
      IF (icount > 0) THEN
        pr_b(nx,ny) = pr_b(nx,ny)/REAL(icount)
      END IF

      ! PR_B now contains the background precipitation
      ! rate & we can determine which points are initially
      ! defined as convective. Put this information in ICLASSTEMP

      IF (pr_rad(nx,ny) > con_intense) THEN ! intensity
        iclasstemp(nx,ny) = 1
      ELSE IF (pr_rad(nx,ny) > con_peak*pr_b(nx,ny)) THEN ! peakednesss
        iclasstemp(nx,ny) = 1
      ELSE
        iclasstemp(nx,ny) = 0
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO

    ! Next stage is to assign iclasstemp = 1 within convective radius
    ! thresholds of the convective grid points defined above
    ! here a code of -1 will mean convective from the convective radius
    ! when we transfer back to the normal grid this will be set back to 1

!$OMP PARALLEL DO PRIVATE(ny,nx,ncon,iy,ix, i,j, radius, ryrange, nyrange,  &
!$OMP& nxrange, dist) DEFAULT(NONE) SHARED(                                 &
#if defined(GNU_FORTRAN)
! Make these two arrays SHARED if they don't have the PARAMETER attribute:
!$OMP& conrad_thres, conrad_deg,                                            &
#endif
!$OMP& nradypts, nradxpts, pr_b, iclasstemp)
  DO ny = 1,nradypts
    DO nx = 1,nradxpts
      IF (iclasstemp(nx,ny) == 1) THEN

        radius = 0.0
        DO ncon=1,nconradp
          IF (pr_b(nx,ny) < conrad_thres(ncon) .AND. radius == 0.0) THEN
            radius = conrad_deg(ncon)
          END IF
        END DO

        ! Now lets work out the grid point range as before

        ryrange = radius/radgrid_deg
        nyrange = INT(ryrange) + 1
        nxrange = nyrange

        DO iy = ny - nyrange, ny + nyrange
          DO ix = nx - nxrange, nx + nxrange

            ! Check that point is within radius of the centre

            dist = SQRT( ((nx - ix)*radgrid_deg)**2 +    &
                            ((ny - iy)*radgrid_deg)**2 )
            IF (dist < radius .AND. dist > 0.0 ) THEN

              ! Do wrap around on IX and IY. Only correct if bicyclic
              ! Otherwise hope discarded area of limited area large enough
              ! for this approach to make no difference to results

              i = ix
              IF (i < 1)        i = i + nradxpts
              IF (i > nradxpts) i = i - nradxpts
              j = iy
              IF (j < 1)        j = j + nradypts
              IF (j > nradypts) j = j - nradypts

              ! Set ICLASSTEMP = -1 if not = 1 already

              IF (iclasstemp(i,j) /= 1) THEN
                iclasstemp(i,j) = -1
              END IF
            END IF  ! dist < radius
          END DO    ! IX
        END DO      ! IY

      END IF        ! iclasstemp == 1
    END DO          ! NY
  END DO            ! NX

  ! If everything's gone fine then iclasstemp should have all
  ! the CONVECTIVE points identified. now do STRATIFORM points 

!$OMP PARALLEL DO PRIVATE(ny,nx,iy,ix,i,j) DEFAULT(NONE)                   &
!$OMP& SHARED(nradypts, nradxpts, nypts, nxpts, iclasstemp, pr_rad,        &
!$OMP& iclass_full)
  DO ny=1,nradypts
    DO nx=1,nradxpts
      IF (iclasstemp(nx,ny) == 0 .AND. pr_rad(nx,ny) > rainy_val) THEN
        iclasstemp(nx,ny) = 2
      END IF


      ! Now we must transfer all this data back onto the model grid
      !  i.e. into ICLASS - follow roughly the same procedure for
      ! putting it onto the RADGRID grid in the first place 
      IF (iclasstemp(nx,ny) /= 0) THEN
        DO iy=1,nypts
          DO ix=1,nxpts
            i = (nx-1)*nxpts + ix
            j = (ny-1)*nypts + iy
            iclass_full(i,j) = REAL(ABS(iclasstemp(nx,ny)))
          END DO
        END DO
      END IF       ! ICLASSTEMP /= 0 

    END DO         ! NY
  END DO           ! NX
!$OMP END PARALLEL DO

  DEALLOCATE(pr_b)
  DEALLOCATE(pr_rad)
  DEALLOCATE(iclasstemp)

END IF    ! mype == 0

!---------------------------------------------------------------------------
! Force synchronisation before trying to scatter back fields
!---------------------------------------------------------------------------
icode = 0
CALL  gc_gsync(nprocs,icode)

! Now scatter back column classification to Nodes (only valid for real arrays)
! DEPENDS ON: scatter_field
CALL scatter_field(field_out_local, iclass_full,                      &
                  local_row_len,local_rows,                           &
                  in_cols,in_rows,                                    &
                  fld_type_p,halo_type_no_halo,                       &
                  0,gc_all_proc_group)


!---------------------------------------------------------------------------
! Convert back to an integer from a real
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                                &
!$OMP& SHARED(local_rows, local_row_len , iclass_col, iwpath, lwpath,       &
!$OMP& field_out_local)    
DO j = 1,local_rows
  DO i = 1,local_row_len
    iclass_col(i,j) = INT(field_out_local(i,j)) 
    iwpath(i,j) = 0.0
    lwpath(i,j) = 0.0
  END DO
END DO
!$OMP END PARALLEL DO

! Calculate liquid and ice water paths on hires grid
! Correct for specific humidities and wet density as variables

DO k=1,mlevs
  DO j = 1,local_rows
    DO i = 1,local_row_len
      IF (mask(i,j,k)) THEN  ! above surface
        qicetemp=  qcf(i,j,k)
        qltemp=  qcl(i,j,k)+qrain(i,j,k)
        IF (l_qgraup) THEN
          qicetemp= qicetemp +qgraup(i,j,k)
        END IF 
        iwpath(i,j) = iwpath(i,j) + qicetemp*density(i,j,k)*dz_theta_sea(k)
        lwpath(i,j) = lwpath(i,j) + qltemp*density(i,j,k)*dz_theta_sea(k)
      END IF 
    END DO
  END DO
END DO

! ICLASS now contains all the convective(1) and stratiform(2) points
! perform the rest of the column classification
 
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                                &
!$OMP& SHARED(local_rows, local_row_len , iclass_col, iwpath, lwpath)
DO j=1,local_rows
  DO i=1,local_row_len
    IF (iclass_col(i,j) == 0) THEN 
      ! 3 - trailing region
      IF (iwpath(i,j) > trail_iwp) THEN
        iclass_col(i,j) = 3
        ! 4 - ice anvil
      ELSE IF (iwpath(i,j) > anvil_iwp) THEN
        iclass_col(i,j) = 4 
        ! 5 - shallow cloud - non precipitating at surface?
      ELSE IF (lwpath(i,j) > shallow_lwp) THEN
        iclass_col(i,j) = 5
      END IF
    END IF
  END DO
END DO
!$OMP END PARALLEL DO

nvalues = REAL(new_res(1)*new_res(1)) 
rnvalues = 1.0/nvalues

! Classification means for single level fields
!$OMP PARALLEL DO PRIVATE(i,j,ii,jj,iii,jjj) DEFAULT(NONE)                  &
!$OMP& SHARED(local_new_x, local_new_y, new_res, rnvalues, iclass_col,      &
!$OMP& precip, fract_conv,fract_strat, prec_conv, prec_strat)    
DO jj=1,local_new_y
  DO ii= 1,local_new_x
    fract_conv(ii,jj)  = 0.0
    fract_strat(ii,jj) = 0.0
    prec_conv(ii,jj)  = 0.0
    prec_strat(ii,jj) = 0.0

    DO j=1,new_res(1)
      jjj = (jj-1)*new_res(1)+j
      DO i=1,new_res(1)
        iii = (ii-1)*new_res(1)+i
        IF (iclass_col(iii,jjj) == 1) THEN         ! Convective
          fract_conv(ii,jj)  = fract_conv(ii,jj) + 1.0
          prec_conv(ii,jj)   = prec_conv(ii,jj) + precip(iii,jjj)
        ELSE IF (iclass_col(iii,jjj) == 2) THEN    ! Stratiform
          fract_strat(ii,jj)  = fract_strat(ii,jj) + 1.0
          prec_strat(ii,jj)   = prec_strat(ii,jj) + precip(iii,jjj)
        END IF
      END DO
    END DO

    ! Calculate prec means first - gives means calculated by dividing by all
    ! points so see proportion of total mean which is conv or stratiform
    prec_conv(ii,jj)   = prec_conv(ii,jj) *rnvalues
    prec_strat(ii,jj)  = prec_strat(ii,jj)*rnvalues
    fract_conv(ii,jj)  = fract_conv(ii,jj) *rnvalues
    fract_strat(ii,jj) = fract_strat(ii,jj)*rnvalues

  END DO
END DO
!$OMP END PARALLEL DO

!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!--------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_class_col

END MODULE crmstyle_class_col_mod

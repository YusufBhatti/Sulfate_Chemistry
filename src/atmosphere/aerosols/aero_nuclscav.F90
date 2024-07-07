! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Perform nucleation scavenging of aged aerosol to cloud aerosol. 
! (aerosol either bmas or ocff)
MODULE aero_nuclscav_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='AERO_NUCLSCAV_MOD'

CONTAINS

SUBROUTINE aero_nuclscav(                                                      &
  ! Arguments IN
  rows, row_length, off_x, off_y, halo_i, halo_j,                              &
  cloudf, qcl, qcf,                                                            &
  aerosol_agd, aerosol_cld,                                                    &
  ! Arguments OUT
  delta_aeronuclscav                                                           &
  )

! Purpose:
!  To perform nucleation scavenging of aged aerosol to form
!   the third mode of aerosol, aerosol in cloud water.
!
!   Called by Aero_ctl
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
!
! Code Description:
!  Language: Fortran 90.
!  This code is written to UMDP3 v8 programming standards
!
! Documentation: UMDP20

USE timestep_mod, ONLY: timestep

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE nlsizes_namelist_mod, ONLY: model_levels

USE c_aero_chm_mod, ONLY: cloudtau, evaptau, nuctau, thold
IMPLICIT NONE

!  includes parameters for rate of aerosol nucleation/evaporation
!  includes gas constant, R

! Arguments with intent IN:

INTEGER :: rows                 ! no. of rows
INTEGER :: row_length           ! no. of pts along a row
INTEGER :: off_x                ! size of small halo in i
INTEGER :: off_y                ! size of small halo in j
INTEGER :: halo_i               ! EW halo size
INTEGER :: halo_j               ! NS halo size

REAL :: cloudf(row_length,rows,model_levels)        ! Decimal cloud fraction
! Cloud liquid water
REAL :: qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,model_levels)
! Cloud frozen water
REAL :: qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,model_levels)

! Arguments with intent IN:
REAL ::aerosol_agd(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)
                                !mmr of aged aerosol
! mmr of aerosol-in-cloud
REAL ::aerosol_cld(1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)

! Arguments with intent OUT:

! cloud aerosol increment due to nucleation scavenging
REAL :: delta_aeronuclscav(row_length,rows,model_levels)

! Local variables:

INTEGER :: i,j,k  ! loop counters

REAL :: delta_nuc(row_length,rows,model_levels)   ! Increment to cloud aerosol
REAL :: delta_evap(row_length,rows,model_levels)  ! Increment to aged aerosol
REAL :: qctotal(row_length,rows,model_levels)     ! Total cloud water
REAL :: clear_frac                 ! Clear fraction of grid box (1.0 - cloudf)

REAL :: evaptime     ! timescale for cloud droplets to evaporate
REAL :: nuctime      ! timescale for particles to enter a cloud and nucleate.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='AERO_NUCLSCAV'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!-----------------------------------------------------------------------
! 1. Initialise increments to zero
!-----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,evaptime,clear_frac,nuctime)
!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      delta_nuc(i,j,k)=0.0
      delta_evap(i,j,k)=0.0
    END DO
  END DO
END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 2. Release aerosol from evaporating cloud droplets in partly cloudy grid
!    boxes. Also release any aerosol-in-cloud in cloud-free grid boxes.
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      qctotal(i,j,k)=qcl(i,j,k) + qcf(i,j,k)
      IF (qctotal(i,j,k)  <   thold) THEN
        delta_evap(i,j,k) = aerosol_cld(i,j,k)
        !                    evaporate all the cloud aerosol in this grid box
      ELSE IF (cloudf(i,j,k)  <   0.95) THEN
        evaptime=evaptau + 0.5*cloudtau
        delta_evap(i,j,k) = (1.0 - EXP(-timestep/evaptime)) *                  &
          aerosol_cld(i,j,k)
      ELSE
        delta_evap(i,j,k) = 0.0
      END IF
    END DO
  END DO
END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 3. In-cloud scavenging of aged aerosol particles by cloud droplets.
!    It is assumed that aged aerosol particles can nucleate cloud
!    droplets.
!-----------------------------------------------------------------------


!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length

      clear_frac = 1.0 - cloudf(i,j,k)

      IF ((qctotal(i,j,k)  >=  thold) .AND.                                    &
        (cloudf(i,j,k)  >  0.0)) THEN

        nuctime=nuctau+((clear_frac*cloudtau)/(2.0*cloudf(i,j,k)))
        delta_nuc(i,j,k)=(1.0-EXP(-timestep/nuctime))*aerosol_agd(i,j,k)

      END IF ! Test on qctotal and thold
    END DO
  END DO
END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 4. Calculate total increment for output
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      delta_aeronuclscav(i,j,k) =                                             &
        delta_nuc(i,j,k) - delta_evap(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE aero_nuclscav
END MODULE aero_nuclscav_mod

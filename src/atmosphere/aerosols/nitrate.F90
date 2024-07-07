! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE nitrate_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NITRATE_MOD'

CONTAINS

SUBROUTINE nitrate(                                                &
  ! Arguments IN
  halo_i, halo_j, off_x, off_y,                                    &
  row_length, rows,                                                &
  tstep, cloudf, p, t, q, qcl, qcf,                                &
  ! Arguments IN/OUT
  hno3, nh3, nitr_acc, nitr_diss,                                  &
  ! Arguments OUT
  delta_n_chem, delta_n_evap, delta_n_nuc  )

USE nlsizes_namelist_mod, ONLY: model_levels

USE ukca_constants, ONLY: c_n, c_hno3

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE c_aero_chm_mod, ONLY: cloudtau, evaptau, nuctau, thold
USE c_rmol, ONLY: rmol
USE science_fixes_mod, ONLY: l_fix_nh4no3_equilibrium
USE chemistry_constants_mod, ONLY: avogadro

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    l_new_qsat_chem_aero !Currently defaults to FALSE

IMPLICIT NONE

!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Aerosols

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NITRATE'

! Arguments with intent IN:

INTEGER :: row_length          !no. of pts along a row
INTEGER :: rows                !no. of rows
INTEGER :: halo_i              !EW halo size
INTEGER :: halo_j              !NS halo size
INTEGER :: off_x
INTEGER :: off_y

! temperature on theta levels
REAL :: t(row_length,rows,model_levels)
! pressure on theta levels
REAL :: p(1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)
!chemistry tstep: LE physics tstep
REAL :: tstep
!cloud fraction (0-1)
REAL :: cloudf(row_length,rows,model_levels)
! specific humidity
REAL :: q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j, model_levels)
! cloud liquid water
REAL :: qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,model_levels)
! cloud frozen water
REAL :: qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,model_levels)
REAL :: drh                     !deliquescence relative humidity (as a fraction)

! Arguments with intent IN/OUT:
! HNO3 MMR from UKCA (kg[HNO3]/kg[air])
REAL :: hno3(row_length,rows,model_levels)
! NH3 MMR from CLASSIC scheme (kg[N]/kg[air])
REAL :: nh3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)
! MMR (kg[N]/kg[air])
REAL :: nitr_acc(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)
! MMR (kg[N]/kg[air])
REAL :: nitr_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)

! Arguments with intent OUT (diagnostics):

! Flux thro' reaction HNO3 + NH3 --> NH4NO3
REAL :: delta_n_chem(row_length, rows, model_levels)
! nitr_diss released by evapn of cloud droplets to nitr_diss per ts
REAL    ::delta_n_evap(row_length, rows, model_levels)
! nitr_acc transferred by nucleation to nitr_diss per timestep.
REAL    ::delta_n_nuc(row_length, rows, model_levels)

! Local variables:

INTEGER :: i,j,k               ! loop variables
REAL    :: lnp1, lnp2, lnp3, p1, p2, p3
REAL    :: kp1, kp
REAL    :: total_nh3  ! nitr_acc + NH3
REAL    :: total_hno3 ! nitr_acc + HNO3
REAL    :: old_nh4no3
REAL    :: kp_rt2
REAL    :: evaptime    ! timescale for cloud droplets to evaporate
REAL    :: nuctime     ! timescale for particles to enter a cloud and nucleate.
!saturation water vapour mixing ratio
REAL    :: q_saturation(row_length, rows, model_levels)
! Conversion factor MMR --> concentration for NH3
REAL     :: convert_nh3(row_length, rows, model_levels)
! Conversion factor MMR --> concentration for NH4NO3
REAL     :: convert_nh4no3(row_length, rows, model_levels)
! Conversion factor MMR --> concentration for HNO3
REAL     :: convert_hno3(row_length, rows, model_levels)
REAL     :: air_density(row_length, rows, model_levels)
! relative humidity as a fraction
REAL     :: rel_humid_frac(row_length, rows, model_levels)
! pressure (no halo)
REAL     :: p_array(row_length, rows, model_levels)
! specific humidity (no halo)
REAL     :: q_array(row_length, rows, model_levels)
! QCL, no halo
REAL     :: qcl_array(row_length, rows, model_levels)
! QCF, no halo
REAL     :: qcf_array(row_length, rows, model_levels)
! NH3 (no halo)
REAL     :: nh3_array(row_length, rows, model_levels)
! accumulation-mode nitrate (no halo)
REAL     :: nitr_acc_array(row_length, rows, model_levels)
! dissolved nitrate (no halo)
REAL     :: nitr_diss_array(row_length, rows, model_levels)
!clear air fraction
REAL     :: clearf(row_length,rows,model_levels)
REAL     :: qctotal(row_length,rows,model_levels)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

p_array(:,:,:) = p(1:row_length, 1:rows, :)
q_array(:,:,:) =  q(1:row_length, 1:rows, :)
qcl_array(:,:,:) =  qcl(1:row_length, 1:rows, :)
qcf_array(:,:,:) =  qcf(1:row_length, 1:rows, :)
nh3_array(:,:,:) = nh3(1:row_length, 1:rows, :)
nitr_acc_array(:,:,:) = nitr_acc(1:row_length, 1:rows, :)
nitr_diss_array(:,:,:) = nitr_diss(1:row_length, 1:rows, :)

!Initialise relative humidity to zero everywhere:
rel_humid_frac(:,:,:) = 0.0

! Initialise arrays to zero everywhere:
q_saturation(:,:,:) = 0.0
convert_nh3(:,:,:) = 0.0
convert_nh4no3(:,:,:) = 0.0
convert_hno3(:,:,:) = 0.0
delta_n_chem(:,:,:) = 0.0
delta_n_evap(:,:,:) = 0.0
delta_n_nuc(:,:,:) = 0.0

IF (l_new_qsat_chem_aero) THEN
  DO k=1,model_levels
    ! Calculate saturation mixing ratio:
    CALL qsat_new(q_saturation(:,:,k),t(:,:,k),p_array(:,:,k),row_length,rows)
    ! Calculate relative humidity as a fraction on model_levels:
    rel_humid_frac(:,:,k) = q_array(:,:,k) / q_saturation(:,:,k)
  END DO
ELSE
  ! DEPENDS ON: qsat
  DO k=1,model_levels
    ! Calculate saturation mixing ratio:
    CALL qsat(q_saturation(:,:,k),t(:,:,k),                                   &
      p_array(:,:,k), SIZE(t(:,:,k)))
    ! Calculate relative humidity as a fraction on model_levels:
    rel_humid_frac(:,:,k) = q_array(:,:,k)                                    &
      / q_saturation(:,:,k)
  END DO
END IF
!
WHERE (rel_humid_frac < 0.0)
  rel_humid_frac=0.0           ! remove negatives
END WHERE
WHERE (rel_humid_frac >= 1.0)
  rel_humid_frac=0.999         ! remove values >= 1.0
END WHERE

! Calculate air density for use in MMR --> concentration coversion factors
air_density(:,:,:) = p_array(:,:,:) * avogadro                                 &
  * 1.0e-06 / (rmol * t(:,:,:))

! ---------------------------------------------------------------------
! Calculate MMR --> concentration coversion factors
convert_nh3(:,:,:) = air_density(:,:,:) / c_n

! Prior to UM 10.3, there has always been a factor 2.0 in the
! denominator of the following calculation. However that is wrong.
! Note that nitr_acc_array, which will be converted from MMR to
! molec cm-3 by the application of the conversion factor convert_nh4no3,
! is the MMR of N in NO3 (i.e. only 1 N atom) instead of the MMR 
! of N in NH4NO3 (i.e. 2 N atoms). Therefore the factor 2.0 is not
! needed.
! This incorrect behaviour is currently maintained by default, but can
! be fixed by setting the temporary logical l_fix_nh4no3_equilibrium
! to TRUE.
IF (l_fix_nh4no3_equilibrium) THEN
  convert_nh4no3(:,:,:) = air_density(:,:,:) / c_n
ELSE
  convert_nh4no3(:,:,:) = air_density(:,:,:) / (2.0 * c_n)
END IF

convert_hno3(:,:,:) = air_density(:,:,:) / (c_hno3)

!----------------------------------------------------------------------
! 1.   Calculation of production of NH4NO3 in reaction of NH3 with HNO3
!----------------------------------------------------------------------

DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length

      ! Deliquescence RH:
      drh = EXP((618.3 / t(i,j,k)) - 2.551)
      ! Calculate dissociation constant:
      IF (rel_humid_frac(i,j,k) < drh) THEN
        kp = (t(i,j,k)**(-6.025))                                              &
          * EXP(118.87 - (24084.0 / t(i,j,k)))
      ELSE

        ! p1 = EXP(-135.94 + (8763.0 / T(i,j,k))) * (T(i,j,k)**19.12):
        lnp1 = -135.94 + (8763.0 / t(i,j,k))                                   &
          + (19.12 * LOG(t(i,j,k)))
        p1 = EXP(lnp1)

        ! p2 = EXP(-122.65 + (9969.0 / T(i,j,k))) * (T(i,j,k)**16.22):
        lnp2 = -122.65 + (9969.0 / t(i,j,k))                                   &
          + (16.22 * LOG(t(i,j,k)))
        p2 = EXP(lnp2)

        ! p3 = EXP(-182.61 + (13875.0 / T(i,j,k))) * (T(i,j,k)**24.46):
        lnp3 = -182.61 + (13875.0 / t(i,j,k))                                  &
          + (24.46 * LOG(t(i,j,k)))
        p3 = EXP(lnp3)

        kp1 = (t(i,j,k)**(-6.025))                                             &
          * EXP(118.87 - (24084.0 / t(i,j,k)))

        kp = (p1 - p2 * (1.0 - rel_humid_frac(i,j,k))                          &
          + p3 * ((1.0 - rel_humid_frac(i,j,k))**2))                         &
          * ((1.0 - rel_humid_frac(i,j,k))**1.75) * kp1
      END IF

      ! Convert from nbar^2 to (molec.cm^-3)^2 units
      kp_rt2 = kp * 1.0e-8 / ((rmol * 1.0e6 * t(i,j,k)                         &
        / avogadro)**2)

      total_nh3 = (nitr_acc_array(i,j,k)                                       &
        * convert_nh4no3(i,j,k))                                               &
        + (nh3_array(i,j,k) * convert_nh3(i,j,k))
      total_hno3 = (nitr_acc_array(i,j,k)                                      &
        * convert_nh4no3(i,j,k))                                               &
        + (hno3(i,j,k) * convert_hno3(i,j,k))

      old_nh4no3 = nitr_acc_array(i,j,k) ! MMR

      IF (total_nh3 * total_hno3 > kp_rt2) THEN
        nitr_acc_array(i,j,k) =  0.5 * (total_nh3 + total_hno3                 &
          - SQRT(((total_nh3+total_hno3)**2)                                   &
          - 4.0*((total_nh3*total_hno3)-kp_rt2)))
        nh3_array(i,j,k) = (total_nh3 - nitr_acc_array(i,j,k))                 &
          / convert_nh3(i,j,k)
        hno3(i,j,k) =  (total_hno3 - nitr_acc_array(i,j,k))                    &
          / convert_hno3(i,j,k)
        nitr_acc_array(i,j,k) = nitr_acc_array(i,j,k)                          &
          / convert_nh4no3(i,j,k)
      ELSE
        nitr_acc_array(i,j,k) =  0.0
        nh3_array(i,j,k) = total_nh3 / convert_nh3(i,j,k)
        hno3(i,j,k) =  total_hno3 / convert_hno3(i,j,k)
      END IF
      !
      nitr_acc_array(i,j,k) = MAX(nitr_acc_array(i,j,k), 0.0)
      nh3_array(i,j,k) = MAX(nh3_array(i,j,k), 0.0)
      hno3(i,j,k) = MAX(hno3(i,j,k), 0.0)

      delta_n_chem(i,j,k) = nitr_acc_array(i,j,k) - old_nh4no3

    END DO        !End i loop
  END DO           !End j loop
END DO              !End k loop

!---------------------------------------------------------------------
! 2. Release of aerosol from evaporating cloud droplets:
!    if no condensed water (liquid + ice) in grid box, release
!    dissolved sulphate as accumulation mode aerosol.
!     If cloud fraction less than 0.95, release some in clear  air.
!--------------------------------------------------------------------

DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      qctotal(i,j,k)= qcl_array(i,j,k) + qcf_array(i,j,k)
      IF (qctotal(i,j,k)  <   thold) THEN
        delta_n_evap(i,j,k) = nitr_diss_array(i,j,k)
      ELSE IF (cloudf(i,j,k) <  0.95) THEN
        evaptime = evaptau + 0.5 * cloudtau
        delta_n_evap(i,j,k)=                                                   &
          (1.0-EXP(-tstep / evaptime)) * nitr_diss_array(i,j,k)
      ELSE
        delta_n_evap(i,j,k)=0.0
      END IF
    END DO
  END DO
END DO

!---------------------------------------------------------------------
! 3. Nucleation of accumulation mode aerosol forming dissolved nitrate
!---------------------------------------------------------------------

!    THIS CODE ASSUMES THAT THE PARAMETER NUCTAU, WHICH IS THE
!    TIMESCALE FOR NUCLEATION ONCE A PARTICLE ENTERS A CLOUD, IS
!    VERY SHORT COMPARED WITH CLOUDTAU.

!cdir collapse
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      clearf(i,j,k) = 1.0 - cloudf(i,j,k) ! Calculate clear air fraction
      IF ((qctotal(i,j,k) >= thold) .AND.                                      &
        (cloudf(i,j,k) >  0.0)) THEN  
        nuctime = nuctau + ((clearf(i,j,k) * cloudtau) /                       &
          (2.0 * cloudf(i,j,k)))
        delta_n_nuc(i,j,k)=                                                    &
          (1.0-EXP(-tstep / nuctime)) * nitr_acc_array(i,j,k)
      END IF

    END DO
  END DO
END DO

!----------------------------------------------------------------
! 4. UPDATE NH3, HNO3, accumulation nitrate and dissolved nitrate
!----------------------------------------------------------------

nitr_acc_array(:,:,:) = nitr_acc_array(:,:,:)                                  &
  + delta_n_evap(:,:,:)                                                        &
  - delta_n_nuc(:,:,:)

nitr_diss_array(:,:,:) = nitr_diss_array(:,:,:)                                &
  + delta_n_nuc(:,:,:)                                                         &
  - delta_n_evap(:,:,:)

nh3(1:row_length, 1:rows, :) = nh3_array(:,:,:)
nitr_acc(1:row_length, 1:rows, :) = nitr_acc_array(:,:,:)
nitr_diss(1:row_length, 1:rows, :) = nitr_diss_array(:,:,:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE nitrate
END MODULE nitrate_mod

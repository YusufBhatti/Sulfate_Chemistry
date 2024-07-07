! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sets up variables before the call to the physics routines.
!
! Subroutine Interface:

SUBROUTINE pre_physics                                                        &
   ! (In)
   ( row_length, rows, nfor, ichgf, qcl, qcf                                  &
   , ch_ug, ch_vg, ilscnt, f_coriolis, lcal360, daycount, stepcount           &
   , a_sw_radstep_diag, a_sw_radstep_prog, l_triffid, npft                    &
   ! (InOut)
   , u, v, ug_scm, vg_scm, npft_trif                                          &
   ! (Out)
   , co2_mmr )

USE s_main_force,             ONLY: ug,                vg,                    &
                                    timestep,          obs,                   &
                                    geoforce,          co2start,              &
                                    co2rate,           co2end,                &
                                    ntrad1,            l_geo_centred

USE rad_input_mod,            ONLY: l_radiation
USE set_rad_steps_mod,        ONLY: set_l_rad_step

USE conversions_mod,          ONLY: isec_per_day

USE nlsizes_namelist_mod,     ONLY: model_levels

USE parkind1,                 ONLY: jpim, jprb
USE yomhook,                  ONLY: lhook, dr_hook

IMPLICIT NONE
!
! Description: This routine sets several variables before the call to
!              the two physics routines.
! Method: It takes the relevant parts from the two physics routines
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: FORTRAN 90
!
!  ARGUMENTS WITH INTENT In

! Parameters

! Model dimensions
INTEGER ::               &
  row_length             &
, rows                   &
, npft

INTEGER ::               &
  nfor                   &! Number of obs. forcing
, ichgf                  &
, ilscnt

! Arrays
REAL ::                                      &
  qcl(row_length,rows,model_levels)          &
, qcf(row_length,rows,model_levels)          &
, ch_ug(row_length,rows,nfor-1,model_levels) &
, ch_vg(row_length,rows,nfor-1,model_levels) &
, f_coriolis(row_length,rows)


! time information
INTEGER ::               &
  daycount               &
, stepcount              &
, a_sw_radstep_diag      &
, a_sw_radstep_prog


! Logicals

LOGICAL ::               &
  lcal360                &
, l_triffid

! ARGUMENTS WITH INTENT In/Out

REAL ::                                   &
  u(row_length,rows,model_levels)         &
, v(row_length,rows,model_levels)         &
, ug_scm(row_length,rows,model_levels)    &
, vg_scm(row_length,rows,model_levels)

INTEGER ::               &
  npft_trif

! ARGUMENTS WITH INTENT Out

REAL ::                  &
  co2_mmr

! local variables.

! loop counters
INTEGER ::               &
  i, j, k
REAL ::  utmp      ! Temporary u-wind

! Dr Hook
!=============================================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRE_PHYSICS'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! External routines:

! --------------------------------------------------------------------
IF (l_triffid) THEN
  npft_trif = npft
ELSE
  npft_trif = 1
END IF

!
!---------------------------------------------------------------------
!     If Geostrophic forcing is chosen
!---------------------------------------------------------------------
!
IF (geoforce) THEN

  IF (ilscnt > 0) THEN
    ug_scm (:,:,:) = ug_scm(:,:,:) + timestep * ch_ug(:,:,ilscnt,:)
    vg_scm (:,:,:) = vg_scm(:,:,:) + timestep * ch_vg(:,:,ilscnt,:)
  END IF

  IF (l_geo_centred) THEN
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          utmp = u(i,j,k)

          u(i,j,k) = ug_scm(i,j,k)                                     &
              + (  (u(i,j,k) - ug_scm(i,j,k))                          &
                 * (1.0 - (0.5*f_coriolis(i,j)*timestep)**2)           &
                 + (v(i,j,k) - vg_scm(i,j,k))                          &
                 * (f_coriolis(i,j)*timestep) )                        &
              / (1.0 + (0.5*f_coriolis(i,j)*timestep)**2)

          v(i,j,k) = vg_scm(i,j,k)                                     &
              + (  (v(i,j,k) - vg_scm(i,j,k))                          &
                 * (1.0 - (0.5*f_coriolis(i,j)*timestep)**2)           &
                 - (utmp-ug_scm(i,j,k))                                &
                 * (f_coriolis(i,j)*timestep) )                        &
              / (1.0 + (0.5*f_coriolis(i,j)*timestep)**2)
        END DO
      END DO
    END DO
  ELSE
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length

          !         Store current u-wind to ensure temporally consistent
          !         updating of v-wind.
          utmp = u(i,j,k)
          u(i,j,k) = u(i,j,k)                                          &
                    - f_coriolis(i,j) * timestep                       &
                                      * (vg_scm(i,j,k)-v(i,j,k))
          v(i,j,k) = v(i,j,k)                                          &
                    + f_coriolis(i,j) * timestep                       &
                                      * (ug_scm(i,j,k)-utmp)
        END DO
      END DO
    END DO
  END IF            ! l_geo_centred
END IF            ! geoforce


!---------------------------------------------------------------------
!     Calculate the CO2 mass mixing ratio using the rate of change
!     (per year)
!---------------------------------------------------------------------

IF (lcal360) THEN
  co2_mmr = co2start + co2rate                                               &
    * ((daycount-1)*isec_per_day + (stepcount-1)*timestep) / 360*isec_per_day
ELSE
  co2_mmr = co2start + co2rate                                               &
    * ((daycount-1)*isec_per_day + (stepcount-1)*timestep) / 365*isec_per_day
END IF

IF (co2_mmr > co2end) THEN
  co2_mmr = co2end
END IF

!---------------------------------------------------------------------
!     Is this a radiation timestep?
!---------------------------------------------------------------------

! If using radiation, call subroutine to set whether this is a radiation
! timestep
IF (l_radiation)  CALL set_l_rad_step(stepcount)


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE pre_physics


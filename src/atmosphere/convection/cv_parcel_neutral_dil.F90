! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate level of neutral buoyancy for a dilute and an undilute ascent.
!
MODULE cv_parcel_neutral_dil_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CV_PARCEL_NEUTRAL_DIL_MOD'
CONTAINS

SUBROUTINE cv_parcel_neutral_dil(nunstable,                                   &
               nlcl_c,k_plume, l_dilute,                                      &
               z_lcl_c, thv_pert, z_full_c, z_half_c, exner_theta_levels_c,   &
               buoyancy, buoyancy_dil, t_dens_env, dqsatdz,                   &
               zh_c,                                                          &
               k_max,k_max_dil,k_neutral,k_neutral_dil,                       &
               max_buoy,max_buoy_dil,ql_ad_c,delthvu_c,                       &
               cape_c,cin_c)

USE planet_constants_mod, ONLY: g

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE nlsizes_namelist_mod,  ONLY: model_levels

IMPLICIT NONE
!
! Description:
!   Calculate the level of neutral buoyancy for a dilute and an undilute parcel
!   ascent.
!   Also calcuate the CAPE and CIN of the ascent
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER, INTENT(IN) :: &
  nunstable              ! Number of parcel ascents

INTEGER, INTENT(IN) :: &
  nlcl_c(nunstable)    & ! Level number of LCL
 ,k_plume(nunstable)     ! start level for surface-driven plume

LOGICAL, INTENT(IN) :: &
  l_dilute               ! true if dilute ascent as well as undilutre ascent

REAL, INTENT(IN)    ::   &
  z_lcl_c(nunstable)     & ! LCL height rounded to nearest model level
 ,thv_pert(nunstable)      ! threshold thv of parcel  (K)

REAL, INTENT(IN)    ::                         &
  z_full_c(nunstable, model_levels)            & ! Height theta lev (m)
 ,z_half_c(nunstable, model_levels)            & ! Height uv lev    (m)
 ,exner_theta_levels_c(nunstable, model_levels)& ! Exner on theta lev
 ,buoyancy(nunstable, model_levels)            & ! undilute parcel buoyancy (K)
 ,buoyancy_dil(nunstable, model_levels)        & ! dilute parcel buoyancy (K)
 ,t_dens_env(nunstable, model_levels)          & ! Density potential temperature
                                                 ! of environment. (K)
 ,dqsatdz(nunstable, model_levels)               ! dqsat/dz along adiabat


REAL, INTENT(INOUT)    :: &
  zh_c(nunstable)           ! BL depth compressed

INTEGER, INTENT(OUT) ::     &
  k_max(nunstable)          & ! level of max parcel buoyancy
 ,k_max_dil(nunstable)      & ! level of max parcel buoyancy dilute
 ,k_neutral(nunstable)      & ! level of neutral parcel buoyancy for undilute
 ,k_neutral_dil(nunstable)    ! level of neutral parcel buoyancy for dilute

REAL, INTENT(OUT)    ::     &
  max_buoy(nunstable)       & ! Maximum buoyancy undilute ascent
 ,max_buoy_dil(nunstable)     ! Maximum buoyancy dilute ascent

REAL, INTENT(OUT)   ::      &
  cape_c(nunstable)         & ! CAPE from undilute parcel ascent (m2/s2)
 ,cin_c(nunstable)          & ! CIN from undilute parcel ascent (m2/s2)
 ,delthvu_c(nunstable)      & ! Integral of undilute parcel buoyancy
                              ! over convective cloud layer (for convection)
 ,ql_ad_c(nunstable)          ! adiabatic liquid water content at inversion
                              ! or cloud top (kg/kg)

!-------------------------------------------------------------------------
! Local variables

INTEGER :: ii, k            !Loop counters


REAL    :: &
  inc      & ! CAPE in layer
 ,dz       & ! layer thickness
 ,factor     ! multiplies thv_pert

REAL    ::                &
  dtv_min(nunstable)      & ! min Tv of parcel in cld layer 1 virtual
                            ! temperature (K).
 ,dtv_min_dil(nunstable)    ! min Tv of parcel in cld layer dilute ascent but
                            ! undilute value
! UNUSED AT PRESENT - left in case required in future
!REAL    ::                &
! ,cape2_c(nunstable)      & ! undilute CAPE from dilute parcel ascent (m2/s2)
! ,cin2_c(nunstable)         ! undilute CIN from dilute parcel ascent (m2/s2)

LOGICAL ::                &
  topprof(nunstable)      & ! Flag set when top of ascent is reached.
 ,topprof_dil(nunstable)    ! Flag set when top of ascent is reached.



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CV_PARCEL_NEUTRAL_DIL'

!-------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_dilute) THEN

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(ii, k, factor, dz, inc)      &  
!$OMP& SHARED(nunstable,k_max,k_neutral,max_buoy,ql_ad_c,delthvu_c,   &
!$OMP&        CAPE_c,CIN_c,dtv_min,topprof,k_max_dil,k_neutral_dil,   &
!$OMP&        max_buoy_dil,dtv_min_dil,topprof_dil,model_levels,      &
!$OMP&        nlcl_c,k_plume,z_lcl_c,thv_pert,  z_full_c,z_half_c,    &
!$OMP&        exner_theta_levels_c,buoyancy, buoyancy_dil,t_dens_env, &
!$OMP&        dqsatdz,zh_c, g) SCHEDULE(STATIC)
! The main index of the loop is "ii" so that it can be
! parallelised using OpenMP. 
  DO ii=1,nunstable

  !-----------------------------------------------------------------------
  ! Initialise output arrays

    topprof(ii)     = .FALSE.

    k_max(ii)          = 1
    k_neutral(ii)      = 1

    dtv_min(ii)   = 0.0
    max_buoy(ii)  = 0.0
    delthvu_c(ii) = 0.0
    ql_ad_c(ii)   = 0.0
    CAPE_c(ii)    = 0.0
    CIN_c(ii)     = 0.0

  ! initialise more arrays

    topprof_dil(ii) = .FALSE.

    k_max_dil(ii)      = 1
    k_neutral_dil(ii)  = 1

    dtv_min_dil(ii)   = 0.0
    max_buoy_dil(ii)  = 0.0
    !    CAPE2_c(ii)    = 0.0
    !    CIN2_c(ii)     = 0.0

    DO  k = 2,model_levels

      !-----------------------------------------------------------------------
      ! Only perform tests if parcel ascent If unstable
      !-----------------------------------------------------------------------
      ! No flag for above_lcl required. Reduce thv_pert by a factor dependent
      ! on height relative to LCL.

      IF (k-1 >  nlcl_c(ii)+1                                            &
                        .AND. z_full_c(ii,k-1) >  1.1*z_lcl_c(ii)) THEN

        ! decrease thv_pert by exp(-(z-zlcl)/1000.)

        factor = EXP( (z_lcl_c(ii)-z_full_c(ii,k))*1.0e-3)

        ! set to zero if z-zlcl >1000?

        IF ((z_full_c(ii,k)-z_lcl_c(ii)) >   1000.0) THEN
          factor =0.0
        END IF

      ELSE
        factor= 1.0
      END IF

      !-----------------------------------------------------------------------
      ! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
      !-----------------------------------------------------------------------
      ! Not reached LNB continue testing

      IF ( .NOT. topprof(ii) .AND. k >  k_plume(ii) ) THEN

        IF (buoyancy(ii,k) >  max_buoy(ii)) THEN
          max_buoy(ii) = buoyancy(ii,k)
          k_max(ii)    = k
        END IF

        ! Is parcel still buoyant ?

        IF ( (buoyancy(ii,k)  <=  - thv_pert(ii))                        &
          !                      or reached top of model
                   .OR. (k  >   model_levels-1)  ) THEN

          k_neutral(ii) = k-1
          topprof(ii) = .TRUE.

        END IF
      END IF
      !-----------------------------------------------------------------------
      ! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
      ! Dilute plume
      !-----------------------------------------------------------------------
      ! Not reached LNB continue testing

      IF ( .NOT. topprof_dil(ii) .AND. k >  k_plume(ii) ) THEN
        IF (buoyancy_dil(ii,k) >  max_buoy(ii)) THEN
          max_buoy_dil(ii) = buoyancy_dil(ii,k)
          k_max_dil(ii)    = k
        END IF

        ! Is parcel still buoyant ?

        IF ( (buoyancy_dil(ii,k)  <=  - thv_pert(ii)*factor)       &
          !                      or reached top of model
                   .OR. (k  >   model_levels-1)  ) THEN

          k_neutral_dil(ii) = k-1
          topprof_dil(ii) = .TRUE.
          zh_c(ii) = z_half_c(ii,k)

          IF ( delthvu_c(ii)  >   0.0) THEN
            ! Note delthvu_c is an undilute value calculated to top of dilute
            ! plume.
            ! compensate for any negative buoyancy of parcel in cloud layer
            delthvu_c(ii) = delthvu_c(ii) - dtv_min_dil(ii) *      &
                                      ( z_half_c(ii,k) - z_lcl_C(ii) )
          END IF
        END IF
      END IF

      !-----------------------------------------------------------------------
      ! While doing parcel ascent
      ! (a) find minimum buoyancy
      ! (b) integrate CAPE over the ascent
      !-----------------------------------------------------------------------

      IF (k > nlcl_c(ii) .AND. k < model_levels ) THEN

        dz = z_half_c(ii,k+1) - z_half_c(ii,k)
        inc = g *  buoyancy(ii,k) * dz/t_dens_env(ii,k)


        !----------------------------------------------------------
        ! If not reached level of neutral buoyancy (top of ascent)
        !----------------------------------------------------------

        IF (.NOT. topprof(ii)) THEN

          ! Note only calculating CIN and CAPE from ascents reaching
          ! level of neutral buoyancy. This may not always correspond
          ! to the diagnosed top for the convection scheme.

          IF (inc <  0.0) THEN
            CIN_c(ii)  = CIN_c(ii) + inc
          ELSE      ! CAPE holds only postive part
            CAPE_c(ii) = CAPE_c(ii) + inc
          END IF
          ! adiabatic liquid water content at cloud top
          !                            = -dqsat/dz * zcld

          dtv_min(ii) = MIN( dtv_min(ii),                         &
                           buoyancy(ii,k)/exner_theta_levels_c(ii,k)  )

        END IF    ! test on topprof

        ! Not reached top of ascent - dilute ascent
        IF (.NOT. topprof_dil(ii)) THEN
          ! require undilute buoyancy here as calculating delthvu as undilute
          dtv_min_dil(ii) = MIN( dtv_min_dil(ii),                        &
                            buoyancy(ii,k)/exner_theta_levels_c(ii,k) )

          dz = z_half_c(ii,k+1) - z_half_c(ii,k)
          ! undilute value
          delthvu_c(ii) = delthvu_c(ii) + buoyancy(ii,k)* dz             &
                                            /exner_theta_levels_c(ii,k)

          ! calculation of CIN and CAPE from profiles - undilute values
          ! NOT USED at present so commenting out - may want in future?
          !          inc =g * buoyancy(ii,k) * dz /t_dens_env(ii,k)

          !          IF (inc <  0.0) THEN
          !            CIN2_c(ii)  = CIN2_c(ii) + inc
          !          ELSE      ! CAPE holds only positive part
          !            CAPE2_c(ii) = CAPE2_c(ii) + inc
          !          END IF

                    ! ql_ad = -dqsat/dz * zcld

          ql_ad_c(ii) = -1.0* dqsatdz(ii,k)                                &
                          *(z_half_c(ii,k+1) - z_half_c(ii,nlcl_c(ii)+1))
        END IF    ! test on topprof

      END IF

    END DO     ! level loop

  END DO   ! ii loop
!$OMP END PARALLEL DO

  !     write(6,*) ' k neutral ',k_neutral(i),k_neutral_dil(i)
  !     write(6,*) ' k max buoy ',k_max(i),k_max_dil(i)
  !     write(6,*) 'max buoy ',max_buoy(i),max_buoy_dil(i)
ELSE     ! undilute ascent only

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(ii, k, dz, inc)             &  
!$OMP& SHARED(nunstable,k_max,k_neutral,max_buoy,ql_ad_c,delthvu_c,  &
!$OMP&        CAPE_c,CIN_c,dtv_min,topprof,model_levels,             &
!$OMP&        nlcl_c,k_plume,z_lcl_c,thv_pert,  z_full_c,z_half_c,   &
!$OMP&        exner_theta_levels_c,buoyancy, t_dens_env,             &
!$OMP&        dqsatdz,zh_c, g) SCHEDULE(STATIC)
! The main index of the loop is "ii" so that it can be
! parallelised using OpenMP. 
  DO ii=1,nunstable

  !-----------------------------------------------------------------------
  ! Initialise output arrays

    topprof(ii)     = .FALSE.

    k_max(ii)          = 1
    k_neutral(ii)      = 1

    dtv_min(ii)   = 0.0
    max_buoy(ii)  = 0.0
    delthvu_c(ii) = 0.0
    ql_ad_c(ii)   = 0.0
    CAPE_c(ii)    = 0.0
    CIN_c(ii)     = 0.0

    DO  k = 2,model_levels

      !-----------------------------------------------------------------------
      ! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
      !-----------------------------------------------------------------------
      ! Not reached LNB continue testing

      IF ( .NOT. topprof(ii) .AND. k >  k_plume(ii) ) THEN

        IF (buoyancy(ii,k) >  max_buoy(ii)) THEN
          max_buoy(ii) = buoyancy(ii,k)
          k_max(ii)    = k
        END IF

        ! Is parcel still buoyant ?

        IF ( (buoyancy(ii,k)  <=  - thv_pert(ii))                        &
          !                      or reached top of model
                   .OR. (k  >   model_levels-1)  ) THEN

          k_neutral(ii) = k-1
          topprof(ii) = .TRUE.
          zh_c(ii) = z_half_c(ii,k)

          IF ( delthvu_c(ii)  >   0.0) THEN
            ! compensate for any negative buoyancy of parcel in cloud layer
            delthvu_c(ii) = delthvu_c(ii) - dtv_min(ii) *                 &
                                      ( z_half_c(ii,k) - z_lcl_c(ii) )
          END IF
        END IF
      END IF

      !-----------------------------------------------------------------------
      ! While doing parcel ascent
      ! (a) find minimum buoyancy
      ! (b) integrate CAPE over the ascent
      !-----------------------------------------------------------------------

      IF (k > nlcl_c(ii) .AND. k < model_levels ) THEN

        dz = z_half_c(ii,k+1) - z_half_c(ii,k)
        inc = g *  buoyancy(ii,k) * dz/t_dens_env(ii,k)


        !----------------------------------------------------------
        ! If not reached level of neutral buoyancy (top of ascent)
        !----------------------------------------------------------

        IF (.NOT. topprof(ii)) THEN

          ! Note only calculating CIN and CAPE from ascents reaching
          ! level of neutral buoyancy. This may not always correspond
          ! to the diagnosed top for the convection scheme.

          IF (inc <  0.0) THEN
            CIN_c(ii)  = CIN_c(ii) + inc
          ELSE      ! CAPE holds only postive part
            CAPE_c(ii) = CAPE_c(ii) + inc
          END IF
          ! adiabatic liquid water content at cloud top
          !                            = -dqsat/dz * zcld

          ql_ad_c(ii) = -1.0* dqsatdz(ii,k)                        &
                       *(z_half_c(ii,k+1) - z_half_c(ii,nlcl_c(ii)+1))

          dtv_min(ii) = MIN( dtv_min(ii),                         &
                           buoyancy(ii,k)/exner_theta_levels_c(ii,k)  )

          delthvu_c(ii) = delthvu_c(ii) + buoyancy(ii,k)* dz      &
                                            /exner_theta_levels_c(ii,k)

        END IF    ! test on topprof

      END IF

    END DO     ! level loop
  END DO   ! ii loop
!&OMP END PARALLEL DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE cv_parcel_neutral_dil
END MODULE cv_parcel_neutral_dil_mod

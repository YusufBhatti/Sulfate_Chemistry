! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Cloud Aerosol Interacting Microphysics (CASIM)
! Module for allocating and storing prognostics required for the CASIM
! Microphysics

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation (CASIM)

MODULE casim_prognostics

USE casim_switches, ONLY:                                               &
       l_mp_CloudNumber, l_mp_RainNumber,                               &
       l_mp_Rain3mom, l_mp_IceNumber,                                   &
       l_mp_SnowNumber, l_mp_Snow3mom,                                  &
       l_mp_GraupNumber, l_mp_Graup3mom,                                &
       ! Switches for the activated aerosols
       l_mp_ActiveSolLiquid, l_mp_ActiveSolRain,                        &
       l_mp_ActiveInsolIce, l_mp_ActiveSolIce,                          &
       l_mp_ActiveInsolLiquid, l_mp_ActiveSolNumber,                    &
       l_mp_ActiveInSolNumber

USE atm_fields_bounds_mod, ONLY: tdims

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
SAVE

! Description
! Declare the new prognostics and increments
! for the CASIM microphysics scheme
!
! The new prognostics are declared as 3D pointer arrays which
! points to the D1 array (normally this is done in
! control/top_level/set_atm_fields, the pointer is a 1D array
! which is modified to 3D array in the argument list to atmos_phys1.
!
! This module also contains the allocatable arrays for the
! increments which do not exist in D1 but are required during 
! the timestep
!
! Set 3D pointer arrays for the new cloud and ice progs in D1
REAL, POINTER :: CloudNumber(:,:,:)         ! Cloud drop number concentration
REAL, POINTER :: RainNumber(:,:,:)          ! Rain drop number concentration
REAL, POINTER :: Rain3mom(:,:,:)            ! 3rd moment of rain PSD
REAL, POINTER :: IceNumber(:,:,:)           ! Ice particle number concentration
REAL, POINTER :: SnowNumber(:,:,:)          ! Snow particle number concentration
REAL, POINTER :: Snow3mom(:,:,:)            ! 3rd moment of snow PSD
REAL, POINTER :: GraupNumber(:,:,:)         ! Graupel number concentration
REAL, POINTER :: Graup3mom(:,:,:)           ! 3rd moment of graupel PSD

! Set 3D pointer arrays for the activated aerosol and IN prognostic
REAL, POINTER :: ActiveSolLiquid(:,:,:)     ! act soluble aerosol in liquid
REAL, POINTER :: ActiveSolRain(:,:,:)       ! act soluble aerosol in rain
REAL, POINTER :: ActiveInsolIce(:,:,:)      ! act insoluble aerosol in ice
REAL, POINTER :: ActiveSolIce(:,:,:)        ! act soluble aerosol in ice
REAL, POINTER :: ActiveInsolLiquid(:,:,:)   ! act insoluble aerosol in liquid
REAL, POINTER :: ActiveSolNumber(:,:,:)     ! act soluble aerosol number
REAL, POINTER :: ActiveInSolNumber(:,:,:)   ! act insoluble aerosol number

! Set allocatable arrays for the mphys increments
REAL, ALLOCATABLE :: cloudnumber_inc(:,:,:)
REAL, ALLOCATABLE :: RainNumber_inc(:,:,:)
REAL, ALLOCATABLE :: Rain3mom_inc(:,:,:)
REAL, ALLOCATABLE :: IceNumber_inc(:,:,:)
REAL, ALLOCATABLE :: SnowNumber_inc(:,:,:)
REAL, ALLOCATABLE :: Snow3mom_inc(:,:,:)
REAL, ALLOCATABLE :: GraupNumber_inc(:,:,:)
REAL, ALLOCATABLE :: Graup3mom_inc(:,:,:)
! Set allocatable arrys for activated aerosol and IN prognostic mphys increment
REAL, ALLOCATABLE :: dActiveSolLiquid(:,:,:)
REAL, ALLOCATABLE :: dActiveSolRain(:,:,:)
REAL, ALLOCATABLE :: dActiveInsolIce(:,:,:)
REAL, ALLOCATABLE :: dActiveSolIce(:,:,:)
REAL, ALLOCATABLE :: dActiveInsolLiquid(:,:,:)
REAL, ALLOCATABLE :: dActiveSolNumber(:,:,:)
REAL, ALLOCATABLE :: dActiveInSolNumber(:,:,:)

INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

!==============================================================================

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CASIM_PROGNOSTICS'

CONTAINS

SUBROUTINE casim_alloc_increments

IMPLICIT NONE

REAL(KIND=jprb) :: zhook_handle ! For Dr Hook

CHARACTER(LEN=*), PARAMETER :: RoutineName='CASIM_ALLOC_INCREMENTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_mp_cloudnumber) THEN

  ALLOCATE( cloudnumber_inc ( tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end ) )

  cloudnumber_inc(:,:,:) = 0.0

ELSE 

  ALLOCATE( cloudnumber_inc(1,1,1) )

END IF ! l_mp_cloudnumber

IF (l_mp_rainnumber) THEN

  ALLOCATE( rainnumber_inc (  tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end ) )

  rainnumber_inc(:,:,:) = 0.0

ELSE 

  ALLOCATE( rainnumber_inc(1,1,1) )

END IF 

IF (l_mp_rain3mom) THEN

  ALLOCATE( rain3mom_inc (    tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end ) )

  rain3mom_inc(:,:,:) = 0.0

ELSE 

  ALLOCATE( rain3mom_inc(1,1,1) )

END IF 

IF (l_mp_icenumber) THEN

  ALLOCATE( icenumber_inc (  tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end ) )

  icenumber_inc(:,:,:) = 0.0

ELSE

  ALLOCATE( icenumber_inc(1,1,1) )

END IF 

IF (l_mp_snownumber) THEN

  ALLOCATE( snownumber_inc (  tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end ) )

  snownumber_inc(:,:,:) = 0.0

ELSE 

  ALLOCATE( snownumber_inc(1,1,1) )

END IF 

IF (l_mp_snow3mom) THEN

  ALLOCATE( snow3mom_inc (    tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end ) )

  snow3mom_inc(:,:,:) = 0.0

ELSE 

  ALLOCATE( snow3mom_inc(1,1,1) )

END IF 

IF (l_mp_graupnumber) THEN

 ALLOCATE( graupnumber_inc (  tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end ) )

  graupnumber_inc(:,:,:) = 0.0

ELSE 

  ALLOCATE( graupnumber_inc(1,1,1) )

END IF 

IF (l_mp_graup3mom) THEN

  ALLOCATE( graup3mom_inc (   tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end ) )

  graup3mom_inc(:,:,:) = 0.0

ELSE 

  ALLOCATE( graup3mom_inc(1,1,1) )

END IF 

IF (l_mp_activesolliquid) THEN

  ALLOCATE( dactivesolliquid ( tdims%i_start : tdims%i_end,                   &
                               tdims%j_start : tdims%j_end,                   &
                                           1 : tdims%k_end ) )

  dactivesolliquid(:,:,:) = 0.0

ELSE 

  ALLOCATE( dactivesolliquid(1,1,1) ) 

END IF 

IF (l_mp_activesolrain) THEN

  ALLOCATE( dactivesolrain (     tdims%i_start : tdims%i_end,                &
                                 tdims%j_start : tdims%j_end,                &
                                             1 : tdims%k_end ) )

  dactivesolrain(:,:,:) = 0.0


ELSE 

  ALLOCATE( dactivesolrain(1,1,1) ) 

END IF 

IF (l_mp_activeinsolice) THEN

  ALLOCATE( dactiveinsolice (    tdims%i_start : tdims%i_end,                &
                                 tdims%j_start : tdims%j_end,                &
                                             1 : tdims%k_end ) )

  dactiveinsolice(:,:,:) = 0.0

ELSE 

  ALLOCATE( dactiveinsolice(1,1,1) ) 

END IF 

IF (l_mp_activesolice) THEN

  ALLOCATE( dactivesolice (      tdims%i_start : tdims%i_end,                &
                                 tdims%j_start : tdims%j_end,                &
                                             1 : tdims%k_end ) )

  dactivesolice(:,:,:) = 0.0

ELSE 

  ALLOCATE( dactivesolice(1,1,1) ) 

END IF 

IF (l_mp_activeinsolliquid) THEN

  ALLOCATE( dactiveinsolliquid ( tdims%i_start : tdims%i_end,                &
                                 tdims%j_start : tdims%j_end,                &
                                             1 : tdims%k_end ) )

  dactiveinsolliquid(:,:,:) = 0.0

ELSE 

  ALLOCATE( dactiveinsolliquid(1,1,1) )

END IF 

IF (l_mp_activesolnumber) THEN 

  ALLOCATE( dactivesolnumber (   tdims%i_start : tdims%i_end,                &
                                 tdims%j_start : tdims%j_end,                &
                                             1 : tdims%k_end ) )

  dactivesolnumber(:,:,:) = 0.0

ELSE 

  ALLOCATE( dactivesolnumber(1,1,1) )

END IF 

IF (l_mp_activeinsolnumber) THEN

  ALLOCATE( dactiveinsolnumber ( tdims%i_start : tdims%i_end,                &
                                 tdims%j_start : tdims%j_end,                &
                                             1 : tdims%k_end ) )

  dactiveinsolnumber(:,:,:) = 0.0

ELSE

  ALLOCATE( dactiveinsolnumber(1,1,1) )

END IF 

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE casim_alloc_increments

!==============================================================================

SUBROUTINE casim_prognostics_update

IMPLICIT NONE

REAL(KIND=jprb) :: zhook_handle ! for Dr Hook

CHARACTER(LEN=*), PARAMETER :: RoutineName='CASIM_PROGNOSTICS_UPDATE'

INTEGER :: i, j, k ! loop counters

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_mp_cloudnumber) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        cloudnumber(i,j,k) = cloudnumber(i,j,k) + cloudnumber_inc(i,j,k)
      END DO
    END DO
  END DO

  cloudnumber(:,:,0) = cloudnumber(:,:,1)

END IF ! l_mp_cloudnumber

IF (l_mp_rainnumber) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        rainnumber(i,j,k) = rainnumber(i,j,k) + rainnumber_inc(i,j,k)
      END DO
    END DO
  END DO

  rainnumber(:,:,0) = rainnumber(:,:,1)

END IF ! l_mp_rainnumber

IF (l_mp_rain3mom) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        rain3mom(i,j,k) = rain3mom(i,j,k) + rain3mom_inc(i,j,k)
      END DO
    END DO
  END DO

  rain3mom(:,:,0) = rain3mom(:,:,1)

END IF ! l_mp_rain3mom

IF (l_mp_icenumber) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        icenumber(i,j,k) = icenumber(i,j,k) + icenumber_inc(i,j,k)
      END DO
    END DO
  END DO

  icenumber(:,:,0) = icenumber(:,:,1)

END IF ! l_mp_icenumber

IF (l_mp_snownumber) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        snownumber(i,j,k) = snownumber(i,j,k) + snownumber_inc(i,j,k)
      END DO
    END DO
  END DO

  snownumber(:,:,0) = snownumber(:,:,1)

END IF ! l_mp_snownumber

IF (l_mp_snow3mom) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        snow3mom(i,j,k) = snow3mom(i,j,k) + snow3mom_inc(i,j,k)
      END DO
    END DO
  END DO

  snow3mom(:,:,0) = snow3mom(:,:,1)

END IF ! l_mp_snow3mom

IF (l_mp_graupnumber) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        graupnumber(i,j,k) = graupnumber(i,j,k) + graupnumber_inc(i,j,k)
      END DO
    END DO
  END DO

  graupnumber(:,:,0) = graupnumber(:,:,1)

END IF ! l_mp_graupnumber

IF (l_mp_graup3mom) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        graup3mom(i,j,k) = graup3mom(i,j,k) + graup3mom_inc(i,j,k)
      END DO
    END DO
  END DO

  graup3mom(:,:,0) = graup3mom(:,:,1)

END IF ! l_mp_graup3mom


IF (l_mp_activesolliquid) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        activesolliquid(i,j,k) = activesolliquid(i,j,k)                       &
                               + dactivesolliquid(i,j,k)
      END DO
    END DO
  END DO

  activesolliquid(:,:,0) = activesolliquid(:,:,1)

END IF ! l_mp_activesolliquid


IF (l_mp_activesolrain) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        activesolrain(i,j,k) = activesolrain(i,j,k)                           &
                             + dactivesolrain(i,j,k)
      END DO
    END DO
  END DO

  activesolrain(:,:,0) = activesolrain(:,:,1)

END IF ! l_mp_activesolrain

IF (l_mp_activeinsolice) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        activeinsolice(i,j,k) = activeinsolice(i,j,k)                         &
                              + dactiveinsolice(i,j,k)
      END DO
    END DO
  END DO

  activeinsolice(:,:,0) = activeinsolice(:,:,1)

END IF ! l_mp_activeinsolice

IF (l_mp_activesolice) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        activesolice(i,j,k) = activesolice(i,j,k)                             &
                            + dactivesolice(i,j,k)
      END DO
    END DO
  END DO

  activesolice(:,:,0) = activesolice(:,:,1)

END IF ! l_mp_activesolice

IF (l_mp_activeinsolliquid) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        activeinsolliquid(i,j,k) = activeinsolliquid(i,j,k)                   &
                                 + dactiveinsolliquid(i,j,k)
      END DO
    END DO
  END DO

  activeinsolliquid(:,:,0) = activeinsolliquid(:,:,1)

END IF ! l_mp_activeinsolliquid

IF (l_mp_activesolnumber) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        activesolnumber(i,j,k) = activesolnumber(i,j,k)                       &
                               + dactivesolnumber(i,j,k)
      END DO
    END DO
  END DO

  activesolnumber(:,:,0) = activesolnumber(:,:,1)

END IF ! l_mp_activesolnumber

IF (l_mp_activeinsolnumber) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        activeinsolnumber(i,j,k) = activeinsolnumber(i,j,k)                   &
                                 + dactiveinsolnumber(i,j,k)
      END DO
    END DO
  END DO

  activeinsolnumber(:,:,0) = activeinsolnumber(:,:,1)

END IF ! l_mp_activeinsolnumber

!------------------------------------------------------------------------
! Finally, deallocate the increments in the reverse order to which they
! were allocated
!------------------------------------------------------------------------
DEALLOCATE( dactiveinsolnumber )
DEALLOCATE( dactivesolnumber   )
DEALLOCATE( dactiveinsolliquid )
DEALLOCATE( dactivesolice      )
DEALLOCATE( dactiveinsolice    )
DEALLOCATE( dactivesolrain     )
DEALLOCATE( dactivesolliquid   )
DEALLOCATE( graup3mom_inc      )
DEALLOCATE( graupnumber_inc    )
DEALLOCATE( snow3mom_inc       )
DEALLOCATE( snownumber_inc     )
DEALLOCATE( icenumber_inc      )
DEALLOCATE( rain3mom_inc       )
DEALLOCATE( rainnumber_inc     )
DEALLOCATE( cloudnumber_inc    )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE casim_prognostics_update

!==============================================================================
END MODULE casim_prognostics

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Create a NULL cosp_gridbox variable to prevent runaway memory usage when
! the UM is run with COSP inactive
!
! Description:
!   A module to allocate and deallocate the cosp gridbox variable during
!   model runs in which COSP is not activated. This is required to ensure
!   efficicent use of memory during model runs
!
! Method:
!   Subroutine allocate_null_gbx allocates all dimensions of all arrays to
!   have a length of 1
!   Subroutine free_null_gbx deallocates these arrays
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 95
!   This code is written to UMDP3 standards

MODULE null_cosp_gridbox_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NULL_COSP_GRIDBOX_MOD'

CONTAINS

SUBROUTINE allocate_null_gbx()

! Description:
!   Allocate all arrays in the cosp gridbox to have lengths 1 in all
!   dimensions

USE cosp_variable_mod, ONLY: cosp_gbx
USE parkind1, ONLY: jpim, jprb
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_NULL_GBX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! allocate the arrays
ALLOCATE(cosp_gbx%zlev(1,1), cosp_gbx%zlev_half(1,1), &
   cosp_gbx%dlev(1,1), cosp_gbx%p(1,1), cosp_gbx%ph(1,1), &
   cosp_gbx%t(1,1), cosp_gbx%q(1,1), cosp_gbx%sh(1,1), &
   cosp_gbx%dtau_s(1,1), cosp_gbx%dtau_c(1,1), &
   cosp_gbx%dem_s(1,1), cosp_gbx%dem_c(1,1), &
   cosp_gbx%tca(1,1), cosp_gbx%cca(1,1), &
   cosp_gbx%rain_ls(1,1), cosp_gbx%rain_cv(1,1), &
   cosp_gbx%grpl_ls(1,1), cosp_gbx%snow_ls(1,1),  &
   cosp_gbx%snow_cv(1,1),cosp_gbx%mr_ozone(1,1))

ALLOCATE(cosp_gbx%toffset(1), cosp_gbx%longitude(1),cosp_gbx%latitude(1), &
   cosp_gbx%psfc(1), cosp_gbx%land(1), cosp_gbx%sunlit(1),cosp_gbx%skt(1), &
   cosp_gbx%u_wind(1),cosp_gbx%v_wind(1))

ALLOCATE(cosp_gbx%mr_hydro(1,1,1), &
   cosp_gbx%dist_prmts_hydro(1,1), &
   cosp_gbx%Reff(1,1,1), &
   cosp_gbx%Np(1,1,1))

ALLOCATE(cosp_gbx%conc_aero(1,1,1), cosp_gbx%dist_type_aero(1), &
   cosp_gbx%dist_prmts_aero(1,1,1,1))

ALLOCATE(cosp_gbx%ichan(1),cosp_gbx%surfem(1))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE allocate_null_gbx



SUBROUTINE free_null_gbx()

USE cosp_variable_mod, ONLY: cosp_gbx
USE parkind1, ONLY: jpim, jprb
USE yomhook, ONLY: lhook, dr_hook

! Description:
!   Frees the allocated items in cosp_gridbox

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FREE_NULL_GBX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DEALLOCATE(cosp_gbx%zlev, cosp_gbx%zlev_half, &
   cosp_gbx%dlev, cosp_gbx%p, cosp_gbx%ph, &
   cosp_gbx%t, cosp_gbx%q, cosp_gbx%sh, &
   cosp_gbx%dtau_s, cosp_gbx%dtau_c, &
   cosp_gbx%dem_s, cosp_gbx%dem_c, &
   cosp_gbx%tca, cosp_gbx%cca, &
   cosp_gbx%rain_ls, cosp_gbx%rain_cv, &
   cosp_gbx%grpl_ls, cosp_gbx%snow_ls,  &
   cosp_gbx%snow_cv,cosp_gbx%mr_ozone)

DEALLOCATE(cosp_gbx%toffset, cosp_gbx%longitude,cosp_gbx%latitude,&
   cosp_gbx%psfc, cosp_gbx%land, cosp_gbx%sunlit,cosp_gbx%skt, &
   cosp_gbx%u_wind,cosp_gbx%v_wind)

DEALLOCATE(cosp_gbx%mr_hydro, &
   cosp_gbx%dist_prmts_hydro, &
   cosp_gbx%Reff, &
   cosp_gbx%Np)

DEALLOCATE(cosp_gbx%conc_aero, cosp_gbx%dist_type_aero, &
   cosp_gbx%dist_prmts_aero)

DEALLOCATE(cosp_gbx%ichan,cosp_gbx%surfem)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE free_null_gbx


END MODULE null_cosp_gridbox_mod

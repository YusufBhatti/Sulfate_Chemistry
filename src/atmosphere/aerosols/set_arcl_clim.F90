! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine for the aerosol climatology for NWP
!
! Purpose:
!  Copy individual fields of the aerosol climatology for NWP into
!  a single array. Aerosol species that are not requested are not
!  copied. The array index corresponding to each requested component
!  are stored for latter access to the corresponding field.
!
! Method:
!   Straightforward.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
!
! Description of code:
!   FORTRAN 90
!- ---------------------------------------------------------------------
MODULE set_arcl_clim_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SET_ARCL_CLIM_MOD'

CONTAINS

SUBROUTINE set_arcl_clim(                                                      &
                                ! Array dimensions
  n_arcl_compnts,                                                              &
                                ! Internal model switches
  l_use_arcl,                                                                  &
                                ! Climatologies from ancillary files:
                                !    biomass-burning
  arclbiom_fr, arclbiom_ag, arclbiom_ic,                                       &
                                !    black-carbon
  arclblck_fr, arclblck_ag,                                                    &
                                !    sea-salt
  arclsslt_fi, arclsslt_jt,                                                    &
                                !    sulphate
  arclsulp_ac, arclsulp_ak, arclsulp_di,                                       &
                                !    mineral dust
  arcldust_b1, arcldust_b2, arcldust_b3,                                       &
  arcldust_b4, arcldust_b5, arcldust_b6,                                       &
                                !    fossil-fuel organic carbon
  arclocff_fr, arclocff_ag, arclocff_ic,                                       &
                                !    delta aerosol
  arcldlta_dl,                                                                 &
                                ! Internal climatology array
  arcl,                                                                        &
                                ! Component array indices
  i_arcl_compnts                                                               &
  )

USE arcl_mod,        ONLY: npd_arcl_compnts, npd_arcl_species, ip_arcl_sulp,   &
                           ip_arcl_dust, ip_arcl_sslt, ip_arcl_blck,           &
                           ip_arcl_biom, ip_arcl_ocff, ip_arcl_dlta,           &
                           ip_arcl_sulp_ac, ip_arcl_sulp_ak, ip_arcl_sulp_di,  &
                           ip_arcl_dust_b1, ip_arcl_dust_b2, ip_arcl_dust_b3,  &
                           ip_arcl_dust_b4, ip_arcl_dust_b5, ip_arcl_dust_b6,  &
                           ip_arcl_sslt_fi, ip_arcl_sslt_jt, ip_arcl_blck_fr,  &
                           ip_arcl_blck_ag, ip_arcl_biom_fr, ip_arcl_biom_ag,  &
                           ip_arcl_biom_ic, ip_arcl_ocff_fr, ip_arcl_ocff_ag,  &
                           ip_arcl_ocff_ic, ip_arcl_dlta_dl

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: tdims

USE science_fixes_mod, ONLY: l_fix_arcl_eg_levs

IMPLICIT NONE

!
! Arguments with intent(in)
!

!
! Array dimensions
!
INTEGER, INTENT(IN) ::                                                         &
  n_arcl_compnts

!
! Internal model switches
!
LOGICAL, INTENT(IN) :: l_use_arcl(npd_arcl_species)

!
! Climatologies from ancillary files
!
!    Three components of biomass-burning (fresh, aged, in-cloud):
!
REAL, INTENT(IN) :: arclbiom_fr(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arclbiom_ag(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arclbiom_ic(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end)
!
!    Two components of black-carbon (fresh, aged)
!
REAL, INTENT(IN) :: arclblck_fr(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arclblck_ag(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end)
!
!    Two components of sea-salt (film and jet mode)
!
REAL, INTENT(IN) :: arclsslt_fi(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arclsslt_jt(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end)
!
!    Three components of sulphate (accumulation mode, aitken mode, dissolved)
!
REAL, INTENT(IN) :: arclsulp_ac(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arclsulp_ak(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arclsulp_di(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end)
!
!    Six components of mineral dust (six size bins)
!
REAL, INTENT(IN) :: arcldust_b1(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arcldust_b2(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arcldust_b3(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arcldust_b4(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arcldust_b5(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arcldust_b6(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end)
!
!    Three components of fossil-fuel organic carbon (fresh, aged, in-cloud)
!
REAL, INTENT(IN) :: arclocff_fr(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arclocff_ag(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end),                    &
                    arclocff_ic(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end)
!
!    One component of delta aerosol (for use as an additional aerosol)
!
REAL, INTENT(IN) :: arcldlta_dl(tdims%i_start:tdims%i_end,                     &
                                tdims%j_start:tdims%j_end,                     &
                                tdims%k_start:tdims%k_end)

!
! Arguments with intent(out)
!

!
! Internal climatology array
!
REAL, INTENT(OUT) :: arcl(tdims%i_start:tdims%i_end,                           &
                          tdims%j_start:tdims%j_end,                           &
                                      1:tdims%k_end,                           &
                                      1:n_arcl_compnts)
!
! Component array indices
!
INTEGER, INTENT(OUT) :: i_arcl_compnts(npd_arcl_compnts)

!
! Local variables
!

INTEGER :: i_cmp

INTEGER :: i, j, k, level_offset

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_ARCL_CLIM'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! for endgame, there is a level indexing bug that looks at level 0, instead
! of level 1, which can be reproduced by setting level_offset to 1:
IF (l_fix_arcl_eg_levs) THEN
  level_offset=0
ELSE
  level_offset=1
END IF


i_cmp = 1 ! since this routine has been called, there is at least
! one component to process.

!
! For each requested species, copy the corresponding components
! into the gathering array.
! We also keep track of the index where component has been put.
! An index of -1 denotes components that are not used.
!
IF (l_use_arcl(ip_arcl_sulp)) THEN

  i_arcl_compnts(ip_arcl_sulp_ac) = i_cmp
  i_arcl_compnts(ip_arcl_sulp_ak) = i_cmp + 1
  i_arcl_compnts(ip_arcl_sulp_di) = i_cmp + 2

#if defined (INTEL_FORTRAN) && (INTEL_FORTRAN < 15000000)
! Usage of OpenMP directives is suppressed due to a bug in old Intel compilers
! See UM ticket #2370 for details
#else
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, arcl, i_cmp, arclsulp_ac, level_offset,           &
!$OMP         arclsulp_ak, arclsulp_di )                               &
!$OMP PRIVATE( i, j, k )
#endif
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        arcl(i, j, k, i_cmp  ) = arclsulp_ac(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+1) = arclsulp_ak(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+2) = arclsulp_di(i, j, k-level_offset)

      END DO ! i
    END DO ! j
  END DO ! k
#if defined (INTEL_FORTRAN) && (INTEL_FORTRAN < 15000000)
! Usage of OpenMP directives is suppressed due to a bug in old Intel compilers
#else
!$OMP END PARALLEL DO
#endif

  i_cmp = i_cmp + 3

ELSE

  i_arcl_compnts(ip_arcl_sulp_ac) = -1
  i_arcl_compnts(ip_arcl_sulp_ak) = -1
  i_arcl_compnts(ip_arcl_sulp_di) = -1

END IF

IF (l_use_arcl(ip_arcl_dust)) THEN

  i_arcl_compnts(ip_arcl_dust_b1) = i_cmp
  i_arcl_compnts(ip_arcl_dust_b2) = i_cmp + 1
  i_arcl_compnts(ip_arcl_dust_b3) = i_cmp + 2
  i_arcl_compnts(ip_arcl_dust_b4) = i_cmp + 3
  i_arcl_compnts(ip_arcl_dust_b5) = i_cmp + 4
  i_arcl_compnts(ip_arcl_dust_b6) = i_cmp + 5

#if defined (INTEL_FORTRAN) && (INTEL_FORTRAN < 15000000)
! Usage of OpenMP directives is suppressed due to a bug in old Intel compilers
! See UM ticket #2370 for details
#else
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, arcl, i_cmp, arcldust_b1, level_offset,           &
!$OMP         arcldust_b2, arcldust_b3, arcldust_b4, arcldust_b5,      &
!$OMP         arcldust_b6 )                                            &
!$OMP PRIVATE( i, j, k )
#endif
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        arcl(i, j, k, i_cmp  ) = arcldust_b1(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+1) = arcldust_b2(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+2) = arcldust_b3(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+3) = arcldust_b4(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+4) = arcldust_b5(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+5) = arcldust_b6(i, j, k-level_offset)

      END DO ! i
    END DO ! j
  END DO ! k
#if defined (INTEL_FORTRAN) && (INTEL_FORTRAN < 15000000)
! Usage of OpenMP directives is suppressed due to a bug in old Intel compilers
#else
!$OMP END PARALLEL DO
#endif

  i_cmp = i_cmp + 6

ELSE

  i_arcl_compnts(ip_arcl_dust_b1) = -1
  i_arcl_compnts(ip_arcl_dust_b2) = -1
  i_arcl_compnts(ip_arcl_dust_b3) = -1
  i_arcl_compnts(ip_arcl_dust_b4) = -1
  i_arcl_compnts(ip_arcl_dust_b5) = -1
  i_arcl_compnts(ip_arcl_dust_b6) = -1

END IF

IF (l_use_arcl(ip_arcl_sslt)) THEN

  i_arcl_compnts(ip_arcl_sslt_fi) = i_cmp
  i_arcl_compnts(ip_arcl_sslt_jt) = i_cmp + 1

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, arcl, i_cmp, arclsslt_fi, level_offset,           &
!$OMP         arclsslt_jt )                                            &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        arcl(i, j, k, i_cmp  ) = arclsslt_fi(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+1) = arclsslt_jt(i, j, k-level_offset)

      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END PARALLEL DO

  i_cmp = i_cmp + 2

ELSE

  i_arcl_compnts(ip_arcl_sslt_fi) = -1
  i_arcl_compnts(ip_arcl_sslt_jt) = -1

END IF

IF (l_use_arcl(ip_arcl_blck)) THEN

  i_arcl_compnts(ip_arcl_blck_fr) = i_cmp
  i_arcl_compnts(ip_arcl_blck_ag) = i_cmp + 1

#if defined (INTEL_FORTRAN) && (INTEL_FORTRAN < 15000000)
! Usage of OpenMP directives is suppressed due to a bug in old Intel compilers
! See UM ticket #2370 for details
#else
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, arcl, i_cmp, arclblck_fr, level_offset,           &
!$OMP         arclblck_ag )                                            &
!$OMP PRIVATE( i, j, k )
#endif
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        arcl(i, j, k, i_cmp  ) = arclblck_fr(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+1) = arclblck_ag(i, j, k-level_offset)

      END DO ! i
    END DO ! j
  END DO ! k
#if defined (INTEL_FORTRAN) && (INTEL_FORTRAN < 15000000)
! Usage of OpenMP directives is suppressed due to a bug in old Intel compilers
#else
!$OMP END PARALLEL DO
#endif

  i_cmp = i_cmp + 2

ELSE

  i_arcl_compnts(ip_arcl_blck_fr) = -1
  i_arcl_compnts(ip_arcl_blck_ag) = -1

END IF

IF (l_use_arcl(ip_arcl_biom)) THEN

  i_arcl_compnts(ip_arcl_biom_fr) = i_cmp
  i_arcl_compnts(ip_arcl_biom_ag) = i_cmp + 1
  i_arcl_compnts(ip_arcl_biom_ic) = i_cmp + 2

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, arcl, i_cmp, arclbiom_fr, level_offset,           &
!$OMP         arclbiom_ag, arclbiom_ic )                               &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        arcl(i, j, k, i_cmp  ) = arclbiom_fr(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+1) = arclbiom_ag(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+2) = arclbiom_ic(i, j, k-level_offset)

      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END PARALLEL DO

  i_cmp = i_cmp + 3

ELSE

  i_arcl_compnts(ip_arcl_biom_fr) = -1
  i_arcl_compnts(ip_arcl_biom_ag) = -1
  i_arcl_compnts(ip_arcl_biom_ic) = -1

END IF

IF (l_use_arcl(ip_arcl_ocff)) THEN

  i_arcl_compnts(ip_arcl_ocff_fr) = i_cmp
  i_arcl_compnts(ip_arcl_ocff_ag) = i_cmp + 1
  i_arcl_compnts(ip_arcl_ocff_ic) = i_cmp + 2

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, arcl, i_cmp, arclocff_fr, level_offset,           &
!$OMP         arclocff_ag, arclocff_ic )                               &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        arcl(i, j, k, i_cmp  ) = arclocff_fr(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+1) = arclocff_ag(i, j, k-level_offset)
        arcl(i, j, k, i_cmp+2) = arclocff_ic(i, j, k-level_offset)

      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END PARALLEL DO

  i_cmp = i_cmp + 3

ELSE

  i_arcl_compnts(ip_arcl_ocff_fr) = -1
  i_arcl_compnts(ip_arcl_ocff_ag) = -1
  i_arcl_compnts(ip_arcl_ocff_ic) = -1

END IF

IF (l_use_arcl(ip_arcl_dlta)) THEN

  i_arcl_compnts(ip_arcl_dlta_dl) = i_cmp

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        arcl(i, j, k, i_cmp  ) = arcldlta_dl(i, j, k-level_offset)

      END DO ! i
    END DO ! j
  END DO ! k

  i_cmp = i_cmp + 1

ELSE

  i_arcl_compnts(ip_arcl_dlta_dl) = -1

END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE set_arcl_clim
END MODULE set_arcl_clim_mod

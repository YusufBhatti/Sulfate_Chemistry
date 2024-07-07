! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
MODULE eg_ungroup_tracers_mod
IMPLICIT NONE

! Description:
!      This routine ungroup/unpack a superarray into individual
!      tracers array -- this is the reverse operation of the
!      the routine "eg_group_tracers".
!
! Method: ENDGame formulation version 3.02
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Tracer Advection
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_UNGROUP_TRACERS_MOD'

CONTAINS
SUBROUTINE eg_ungroup_tracers(                                  &
                    super_array_size,                           &
                    super_array,                                &
                    L_CO2_interactive, co2,                     &
                    L_Murk_advect, murk,                        &
                    L_Soot, soot_new, soot_agd, soot_cld,       &
                    l_sulpc_so2, so2, SO4_aitken, so4_accu,     &
                                 so4_diss,                      &
                    L_sulpc_nh3, nh3,                           &
                    L_sulpc_dms, dms,                           &
                    L_dust, dust_div1, dust_div2, dust_div3,    &
                            dust_div4, dust_div5, dust_div6,    &
                    L_biomass, bmass_new, bmass_agd, bmass_cld, &
                    L_ocff, ocff_new, ocff_agd, ocff_cld,       &
                    L_nitrate, nitr_acc, nitr_diss,             &
                    l_use_cariolle, ozone_tracer,               &
                    tracers, tr_vars,                           &
                    tr_ukca, tracer_ukca,                       &
                    L_twobin_dust                               )

USE atm_fields_bounds_mod
USE ukca_option_mod,  ONLY: l_conserve_ukca_with_tr
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_UNGROUP_TRACERS'

INTEGER, INTENT(IN) :: super_array_size
INTEGER, INTENT(IN) :: tr_ukca, tr_vars

REAL, INTENT(OUT)   ::  co2      (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  murk     (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  soot_new (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  soot_agd (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  soot_cld (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  so2      (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  so4_aitken                              &
                                 (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  so4_accu (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  so4_diss (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  nh3      (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  dms      (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  dust_div1(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  dust_div2(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  dust_div3(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  dust_div4(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  dust_div5(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  dust_div6(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  bmass_new(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  bmass_agd(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  bmass_cld(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  ocff_new (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  ocff_agd (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  ocff_cld (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  nitr_acc (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  nitr_diss(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  ozone_tracer                            &
                                 (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)   ::  tracer_ukca                             &
                                 (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end,&
                                  tr_ukca )
REAL, INTENT(OUT)   ::  tracers                                 &
                                 (tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end,&
                                  tr_vars )

REAL, INTENT(IN)  ::  super_array                               &
                       (tdims_s%i_start:tdims_s%i_end,          &
                        tdims_s%j_start:tdims_s%j_end,          &
                        tdims_s%k_start:tdims_s%k_end,          &
                        super_array_size)

LOGICAL, INTENT(IN) ::                                           &
           L_CO2_interactive,L_Murk_advect,L_soot,l_sulpc_so2,   &
           L_sulpc_nh3,L_sulpc_dms,L_biomass,L_dust, L_ocff,     &
           l_use_cariolle, L_nitrate, L_twobin_dust

INTEGER ::  i,j,k,l,tr_count,array_size_count

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

array_size_count=0

! ----------------------------------------------------------------------
! Section 1  carbon cycle.
! ----------------------------------------------------------------------
! Note that throughout this code the OpenMP makes heavy use of implicit
! barriers through END DO and END SINGLE

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,tr_count)

IF (L_CO2_interactive) THEN
!$OMP SINGLE
  array_size_count = array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        co2(i,j,k) = super_array(i,j,k,array_size_count)
      END DO
    END DO
  END DO
!$OMP END DO
END IF   ! L_CO2_interactive

! ----------------------------------------------------------------------
! Section 2  Soot cycle.
! ----------------------------------------------------------------------
IF (L_soot) THEN
!$OMP SINGLE
  array_size_count = array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        soot_new(i,j,k) = super_array(i,j,k,array_size_count)
        soot_agd(i,j,k) = super_array(i,j,k,array_size_count+1)
        soot_cld(i,j,k) = super_array(i,j,k,array_size_count+2)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP SINGLE
  array_size_count=array_size_count + 2
!$OMP END SINGLE
END IF    ! L_soot

! ----------------------------------------------------------------------
! Section 3  Biomass aerosol.
! ----------------------------------------------------------------------
IF (l_Biomass) THEN
!$OMP SINGLE
  array_size_count = array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        bmass_new(i,j,k) = super_array(i,j,k,array_size_count)
        bmass_agd(i,j,k) = super_array(i,j,k,array_size_count+1)
        bmass_cld(i,j,k) = super_array(i,j,k,array_size_count+2)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP SINGLE
  array_size_count=array_size_count + 2
!$OMP END SINGLE
END IF  ! l_biomass

! ----------------------------------------------------------------------
! Section 4  sulphur cycle.
! ----------------------------------------------------------------------
IF (l_sulpc_so2) THEN
!$OMP SINGLE
  array_size_count = array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        so4_aitken(i,j,k) = super_array(i,j,k,array_size_count)
        so4_accu(i,j,k) = super_array(i,j,k,array_size_count+1)
        so4_diss(i,j,k) = super_array(i,j,k,array_size_count+2)
        so2(i,j,k) = super_array(i,j,k,array_size_count+3)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP SINGLE
  array_size_count=array_size_count + 3
!$OMP END SINGLE

  IF (L_sulpc_nh3) THEN
!$OMP SINGLE
    array_size_count = array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims_s%k_start, tdims_s%k_end
      DO j = tdims_s%j_start, tdims_s%j_end
        DO i = tdims_s%i_start, tdims_s%i_end
          nh3(i,j,k) = super_array(i,j,k,array_size_count)
        END DO
      END DO
    END DO
!$OMP END DO
  END IF  ! L_sulpc_nh3

  IF (L_sulpc_dms) THEN
!$OMP SINGLE
    array_size_count = array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims_s%k_start, tdims_s%k_end
      DO j = tdims_s%j_start, tdims_s%j_end
        DO i = tdims_s%i_start, tdims_s%i_end
          dms(i,j,k) = super_array(i,j,k,array_size_count)
        END DO
      END DO
    END DO
!$OMP END DO
  END IF  ! L_sulpc_dms
END IF  ! L_sulpc_SO2

! ----------------------------------------------------------------------
! Section 5  Mineral dust.
! ----------------------------------------------------------------------
IF (L_dust) THEN
!$OMP SINGLE
  array_size_count=array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        dust_div1(i,j,k) = super_array(i,j,k,array_size_count  )
        dust_div2(i,j,k) = super_array(i,j,k,array_size_count+1)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP SINGLE
  array_size_count = array_size_count + 1
!$OMP END SINGLE
  IF (.NOT. L_twobin_dust) THEN
!$OMP SINGLE
    array_size_count = array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims_s%k_start, tdims_s%k_end
      DO j = tdims_s%j_start, tdims_s%j_end
        DO i = tdims_s%i_start, tdims_s%i_end
          dust_div3(i,j,k) = super_array(i,j,k,array_size_count  )
          dust_div4(i,j,k) = super_array(i,j,k,array_size_count+1)
          dust_div5(i,j,k) = super_array(i,j,k,array_size_count+2)
          dust_div6(i,j,k) = super_array(i,j,k,array_size_count+3)
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP SINGLE
    array_size_count=array_size_count + 3
!$OMP END SINGLE
  END IF
END IF    ! L_dust

! ----------------------------------------------------------------------
! Section 6  Fossil-fuel organic carbon aerosol
! ----------------------------------------------------------------------
IF (L_ocff) THEN
!$OMP SINGLE
  array_size_count=array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        ocff_new(i,j,k) = super_array(i,j,k,array_size_count)
        ocff_agd(i,j,k) = super_array(i,j,k,array_size_count+1)
        ocff_cld(i,j,k) = super_array(i,j,k,array_size_count+2)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP SINGLE
  array_size_count=array_size_count + 2
!$OMP END SINGLE
END IF    ! L_ocff

! ----------------------------------------------------------------------
! Section 7  Cariolle ozone tracer.
! ----------------------------------------------------------------------
IF (l_use_cariolle) THEN
!$OMP SINGLE
  array_size_count=array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        ozone_tracer(i,j,k) = super_array(i,j,k,array_size_count)
      END DO
    END DO
  END DO
!$OMP END DO
END IF    ! L_USE_CARIOLLE

! ----------------------------------------------------------------------
! Section 8  Ammonium nitrate aerosol.
! ----------------------------------------------------------------------
IF (L_nitrate) THEN
!$OMP SINGLE
  array_size_count=array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        nitr_acc(i,j,k) = super_array(i,j,k,array_size_count)
        nitr_diss(i,j,k) = super_array(i,j,k,array_size_count+1)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP SINGLE
  array_size_count=array_size_count + 1
!$OMP END SINGLE
END IF    ! L_nitrate

! add further stracers if needed
! ----------------------------------------------------------------------
! Section 9 Free tracers 
! ----------------------------------------------------------------------

IF ( tr_vars > 0 ) THEN 
  DO tr_count=1,tr_vars 
!$OMP SINGLE
    array_size_count=array_size_count + 1 
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims_s%k_start,tdims_s%k_end 
      DO j = tdims_s%j_start,tdims_s%j_end 
        DO i = tdims_s%i_start,tdims_s%i_end 
          tracers(i,j,k,tr_count) =                              &
                             super_array(i,j,k,array_size_count) 
        END DO 
      END DO 
    END DO 
!$OMP END DO
  END DO 
END IF 

! ----------------------------------------------------------------------
! Section 10  UKCA tracers
! ----------------------------------------------------------------------
IF ( tr_ukca > 0 .AND. l_conserve_ukca_with_tr ) THEN
  DO tr_count=1,tr_ukca
!$OMP SINGLE
    array_size_count=array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims_s%k_start,tdims_s%k_end
      DO j = tdims_s%j_start,tdims_s%j_end
        DO i = tdims_s%i_start,tdims_s%i_end
          tracer_ukca(i,j,k,tr_count) =                            &
                 super_array(i,j,k,array_size_count)
        END DO
      END DO
    END DO
!$OMP END DO
  END DO

END IF  ! tr_ukca > 0

! ----------------------------------------------------------------------
! Section 11  Murk cycle
! ----------------------------------------------------------------------
IF (L_Murk_advect) THEN
!$OMP SINGLE
  array_size_count = array_size_count + 1
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        murk(i,j,k) = super_array(i,j,k,array_size_count)
      END DO
    END DO
  END DO
!$OMP END DO
END IF  ! L_Murk_advect

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE eg_ungroup_tracers

END MODULE eg_ungroup_tracers_mod


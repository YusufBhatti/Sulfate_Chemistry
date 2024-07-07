! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE tr_set_phys_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TR_SET_PHYS_MOD'

CONTAINS

SUBROUTINE tr_set_phys_4A(                                      &
                    super_array_size, super_tracer_phys,        &
                    l_co2_interactive, co2,                     &
                    l_murk_advect, murk,                        &
                    l_soot, soot_new, soot_agd, soot_cld,       &
                    l_sulpc_so2, so2, so4_aitken, so4_accu,     &
                                 so4_diss,                      &
                    l_sulpc_nh3, nh3,                           &
                    l_sulpc_dms, dms,                           &
                    l_dust, dust_div1, dust_div2, dust_div3,    &
                            dust_div4, dust_div5, dust_div6,    &
                    l_biomass, bmass_new, bmass_agd, bmass_cld, &
                    l_ocff, ocff_new, ocff_agd, ocff_cld,       &
                    l_nitrate, nitr_acc, nitr_diss,             &
                    l_use_cariolle, ozone_tracer,               &
                    tracer, tracer_ukca,                        &
                    row_length, rows,                           &
                    model_levels, tr_levels, tr_vars, tr_ukca,  &
                    l_init, supertrdims )

! Purpose: Interface routine to initialise tracer fields
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Tracer Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE dust_parameters_mod, ONLY: l_twobin_dust
USE murk_inputs_mod,     ONLY: l_murk_lbc
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE field_types, ONLY: fld_type_p
USE atm_fields_bounds_mod, ONLY: array_dims, tdims_s
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE model_domain_mod, ONLY: model_type, mt_lam

IMPLICIT NONE


TYPE (array_dims)  supertrdims

! Arguments with INTENT IN. ie: Input variables.

! Arguments with INTENT IN. ie: Input variables.

INTEGER ::                                                      &
        ! model dimensions
  row_length                                                    &
                   ! number of points on a row
, rows                                                          &
                   ! number of rows in a theta field
, model_levels                                                  &
                   ! number of model levels
, tr_levels                                                     &
                   ! number of tracer levels
, tr_vars                                                       &
                   ! number of tracers
, tr_ukca                                                       &
                   ! number of ukca tracers
, super_array_size


REAL, INTENT(IN) :: co2                                         &
                  (tdims_s%i_start:tdims_s%i_end,               &
                   tdims_s%j_start:tdims_s%j_end,               &
                   tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: murk                                       &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: soot_new                                   &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: soot_agd                                   &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: soot_cld                                   &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so2                                        &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so4_aitken                                 &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so4_accu                                   &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: so4_diss                                   &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN)  :: nh3                                        &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dms                                        &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN)  :: dust_div1                                  &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div2                                  &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div3                                  &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div4                                  &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div5                                  &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: dust_div6                                  &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: bmass_new                                  &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: bmass_agd                                  &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: bmass_cld                                  &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN)  :: ocff_new                                   &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: ocff_agd                                   &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: ocff_cld                                   &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: nitr_acc                                   &
                 (tdims_s%i_start:tdims_s%i_end,                &
                  tdims_s%j_start:tdims_s%j_end,                &
                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: nitr_diss                                  &
                 (tdims_s%i_start:tdims_s%i_end,                &
                  tdims_s%j_start:tdims_s%j_end,                &
                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN) ::                                             &
  tracer      (tdims_s%i_start:tdims_s%i_end,                   &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end, 1:tr_vars ),      &
  tracer_ukca(tdims_s%i_start:tdims_s%i_end,                    &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end,1:tr_ukca)

! Add cariolle specific parameters for ozone tracer
REAL, INTENT(IN) ::                                             &
  ozone_tracer(tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,             &
                     tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT) :: super_tracer_phys                          &
                    (supertrdims%i_start:supertrdims%i_end,     &
                     supertrdims%j_start:supertrdims%j_end,     &
                     supertrdims%k_start:supertrdims%k_end,     &
                     super_array_size)

LOGICAL ::                                                      &
  l_co2_interactive,                                            &
  l_murk_advect,                                                &
  l_soot,                                                       &
  l_sulpc_so2,                                                  &
  l_sulpc_nh3,                                                  &
  l_sulpc_dms,                                                  &
  l_biomass,                                                    &
  l_dust,                                                       &
  l_ocff,                                                       &
  l_use_cariolle,                                               &
  l_nitrate

LOGICAL :: L_init   ! flag for setting halo values

!      local variables
INTEGER ::                                                      &
  i,j,k,counter        !loop variables

INTEGER :: array_size_count

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TR_SET_PHYS_4A'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
array_size_count=0

! ----------------------------------------------------------------------
! Section 1.1  carbon cycle.
! ----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)

IF (L_CO2_interactive) THEN

!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = 1, rows
      DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) = co2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
END IF   ! L_CO2_interactive

! ----------------------------------------------------------------------
! Section 1.2  Soot cycle.
! ----------------------------------------------------------------------
IF (l_soot) THEN
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = 1, rows
      DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count)   =  soot_new(i,j,k)
        super_tracer_phys(i,j,k,array_size_count+1) =  soot_agd(i,j,k)
        super_tracer_phys(i,j,k,array_size_count+2) =  soot_cld(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP MASTER
  array_size_count=array_size_count +2
!$OMP END MASTER

!$OMP BARRIER
END IF    ! L_soot


! ----------------------------------------------------------------------
! Section 1.3  Biomass aerosol.
! ----------------------------------------------------------------------
IF (l_biomass) THEN
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = 1, rows
      DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count)   =  bmass_new(i,j,k)
        super_tracer_phys(i,j,k,array_size_count+1) =  bmass_agd(i,j,k)
        super_tracer_phys(i,j,k,array_size_count+2) =  bmass_cld(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP MASTER
  array_size_count=array_size_count +2
!$OMP END MASTER

!$OMP BARRIER
END IF  ! l_biomass

! ----------------------------------------------------------------------
! Section 1.4  sulphur cycle.
! ----------------------------------------------------------------------
IF (l_sulpc_so2) THEN
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = 1, rows
      DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count)   =  so4_aitken(i,j,k)
        super_tracer_phys(i,j,k,array_size_count+1) =  so4_accu(i,j,k)
        super_tracer_phys(i,j,k,array_size_count+2) =  so4_diss(i,j,k)
        super_tracer_phys(i,j,k,array_size_count+3) =  so2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP MASTER
  array_size_count=array_size_count +3
!$OMP END MASTER

!$OMP BARRIER

  IF (l_sulpc_nh3) THEN
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
    DO k = tdims_s%k_start, tdims_s%k_end
      DO j = 1, rows
        DO i = 1, row_length
          super_tracer_phys(i,j,k,array_size_count) =  nh3(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO
  END IF  ! L_sulpc_nh3

  IF (l_sulpc_dms) THEN
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
    DO k = tdims_s%k_start, tdims_s%k_end
      DO j = 1, rows
        DO i = 1, row_length
          super_tracer_phys(i,j,k,array_size_count) =  dms(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO
  END IF  ! L_sulpc_dms
END IF  ! L_sulpc_SO2

! ----------------------------------------------------------------------
! Section 1.5  Mineral dust.
! ----------------------------------------------------------------------
IF (l_dust) THEN
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = 1, rows
      DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  dust_div1(i,j,k)
        super_tracer_phys(i,j,k,array_size_count+1) =  dust_div2(i,j,k)
        IF (.NOT. l_twobin_dust) THEN
          super_tracer_phys(i,j,k,array_size_count+2) =  dust_div3(i,j,k)
          super_tracer_phys(i,j,k,array_size_count+3) =  dust_div4(i,j,k)
          super_tracer_phys(i,j,k,array_size_count+4) =  dust_div5(i,j,k)
          super_tracer_phys(i,j,k,array_size_count+5) =  dust_div6(i,j,k)
        END IF
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP MASTER
  IF (l_twobin_dust) THEN
    array_size_count=array_size_count + 1
  ELSE
    array_size_count=array_size_count + 5
  END IF
!$OMP END MASTER

!$OMP BARRIER

END IF    ! L_dust

! ----------------------------------------------------------------------
! Section 1.6  Fossil-fuel organic carbon aerosol
! ----------------------------------------------------------------------
IF (l_ocff) THEN
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = 1, rows
      DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  ocff_new(i,j,k)
        super_tracer_phys(i,j,k,array_size_count+1) =  ocff_agd(i,j,k)
        super_tracer_phys(i,j,k,array_size_count+2) =  ocff_cld(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP MASTER
  array_size_count=array_size_count +2
!$OMP END MASTER

!$OMP BARRIER
END IF    ! L_ocff


! ----------------------------------------------------------------------
! Section 1.7  Cariolle ozone tracer.
! ----------------------------------------------------------------------
IF (l_use_cariolle) THEN
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = 1, rows
      DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  ozone_tracer(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
END IF    ! L_USE_CARIOLLE

! ----------------------------------------------------------------------
! Section 1.8  Ammonium nitrate aerosol.
! ----------------------------------------------------------------------
IF (l_nitrate) THEN
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = 1, rows
      DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) = nitr_acc(i,j,k)
        super_tracer_phys(i,j,k,array_size_count+1) = nitr_diss(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER
END IF    ! L_nitrate

! ----------------------------------------------------------------------
!  ANY NEW NAMED TRACER SPECIES SHOULD BE ADDED HERE
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Section 1.33  Free tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------

IF (model_levels==tr_levels .AND. tr_vars>0) THEN
  DO counter=1,tr_vars
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          super_tracer_phys(i,j,k,array_size_count) =          &
                                            tracer(i,j,k,counter)
        END DO
      END DO
    END DO
!$OMP END DO 
  END DO

END IF  ! tr_vars>0


! ----------------------------------------------------------------------
! Section 1.34  UKCA tracers
! ----------------------------------------------------------------------
IF ( tr_ukca > 0) THEN

  DO counter=1,tr_ukca

!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
    DO k = tdims_s%k_start,tdims_s%k_end
      DO j = 1, rows
        DO i = 1, row_length
          super_tracer_phys(i,j,k,array_size_count) =          &
               tracer_ukca(i,j,k,counter)
        END DO
      END DO
    END DO
!$OMP END DO 

  END DO

END IF     ! tr_ukca>0

! ----------------------------------------------------------------------
! Section 1.99  Murk cycle.  This must be the last Full level field in
!                            the super_array
! ----------------------------------------------------------------------

IF (l_murk_advect) THEN
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = 1, rows
      DO i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  murk(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
END IF  ! L_Murk_advect

!$OMP DO SCHEDULE(STATIC)
  DO k = 1,super_array_size
    DO j = supertrdims%j_start,supertrdims%j_end
      DO i = supertrdims%i_start,supertrdims%i_end
        super_tracer_phys(i,j,0,k) =  super_tracer_phys(i,j,1,k) 
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

IF (l_init) THEN
  ! ----------------------------------------------------------------------
  ! Call SWAPBOUNDS and set external halos for all arrays if required
  ! ----------------------------------------------------------------------

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds(                                             &
       super_tracer_phys,                                       &
        row_length, rows,                                       &
   (supertrdims%k_end-supertrdims%k_start+1)*array_size_count,  &
    supertrdims%halo_i,supertrdims%halo_j, fld_type_p,swap_field_is_scalar)

  IF (model_type == mt_lam) THEN
    ! DEPENDS ON: set_external_halos
    CALL set_external_halos(                                    &
        super_tracer_phys,                                      &
        row_length, rows,                                       &
        (supertrdims%k_end-supertrdims%k_start+1)*array_size_count,  &
        supertrdims%halo_i, supertrdims%halo_j, 0.0)
  END IF

END IF  !L_INIT
IF (l_murk_advect .AND. .NOT.l_murk_lbc) THEN
  IF (model_type == mt_lam) THEN
    ! DEPENDS ON: fill_external_halos
    CALL fill_external_halos(                                &
     super_tracer_phys(supertrdims%i_start,                  &
                       supertrdims%j_start,                  &
                       supertrdims%k_start,                  &
                       array_size_count)                     &
                 , row_length, rows,                         &
                   supertrdims%k_end-supertrdims%k_start+1,  &
                   supertrdims%halo_i, supertrdims%halo_j)
  END IF
END IF  ! L_Murk_advect

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE TR_Set_Phys_4A
END MODULE tr_set_phys_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SL_Tracer1
!
MODULE SL_Tracer1_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SL_TRACER1_MOD'

CONTAINS

SUBROUTINE SL_Tracer1_4A(                                         &
                      super_array,                                &
                      super_array_size,                           &
                      eta_theta_levels,                           &
                      row_length, rows, n_rows, model_levels,     &
                      g_i_pe, g_j_pe,                             &
                      high_order_scheme, monotone_scheme,         &
                      l_high, l_mono,                             &
                      depart_xi1, depart_xi2, depart_xi3,         &
                      co2, L_CO2_interactive,                     &
                      Murk, L_murk_advect,                        &
                      dust_div1,dust_div2,dust_div3,              &
                      dust_div4,dust_div5,dust_div6,              &
                      L_dust,                                     &
                      soot_new, soot_agd, soot_cld, L_soot,       &
                      bmass_new, bmass_agd, bmass_cld,            &
                      l_biomass,                                  &
                      ocff_new, ocff_agd, ocff_cld, l_ocff,       &
                      so2, so4_aitken, so4_accu,                  &
                      so4_diss, nh3, dms,                         &
                      L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms,      &
                      nitr_acc, nitr_diss, L_nitrate,             &
                      tracers, tr_vars,                           &
                      tracer_ukca, tr_ukca,                       &
                      L_use_cariolle, ozone_tracer,               &
                      rho_n, rho_np1, Error_Code)

! Purpose:
!          Performs semi-Lagrangian advection of tracers
!          (based on sl_thermo)
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Tracer Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE dust_parameters_mod, ONLY: l_twobin_dust
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE um_parcore, ONLY: mype, nproc
USE um_parvars, ONLY: halo_i, halo_j, offx, offy, datastart,     &
                      at_extremity,                              &
                      gc_proc_row_group, gc_proc_col_group,      &
                      nproc_x, nproc_y
USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows
USE atm_fields_bounds_mod
USE murk_inputs_mod,    ONLY: l_murk_lbc
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE Field_Types
USE horiz_grid_mod,     ONLY: xi1_p, xi2_p, xi1_p, xi2_p
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

USE model_domain_mod, ONLY: model_type, mt_lam
USE eg_zlf_conservation_mod, ONLY: eg_zlf_conservation
USE eg_zlf_mod, ONLY: zlf_cfl_top_level_tracers
USE dynamics_input_mod, ONLY: l_sl_bc_correction,                &
                              zlf_conservation_tracers_option

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER ::                                                        &
  row_length                                                      &
               ! number of points on a row
, rows                                                            &
               ! number of rows.
, n_rows                                                          &
               ! number of v-rows.
, model_levels                                                    &
               ! Number of model levels.
, g_i_pe(1-halo_i:global_row_length+halo_i)                       &
                       ! processor on my processor-row
                       ! holding a given value in i direction
, g_j_pe(1-halo_j:global_rows      +halo_j)
                       ! processor on my processor-column
                       ! holding a given value in j direction

INTEGER ::                                                        &
  tr_vars                                                         &
, tr_ukca              ! No of UKCA tracers

INTEGER ::                                                        &
  high_order_scheme                                               &
                           ! a code saying which high order
                           ! scheme to use for interpolation
, monotone_scheme
                        ! a code saying which monotone
                        ! scheme to use for interpolation

INTEGER, INTENT(IN) :: super_array_size

LOGICAL ::                                                        &
  l_high                                                          &
                 ! True, if high order interpolation required
, l_mono
                 ! True, if interpolation required to be monotone

LOGICAL, INTENT(IN) ::                                            &
  L_CO2_interactive                                               &
, L_murk_advect                                                   &
, l_dust                                                          &
, L_Soot, L_biomass, L_ocff, L_nitrate                            &
, L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms                           &
, l_use_cariolle

REAL, INTENT(INOUT)  ::                                         &
  super_array(       tdims_l%i_start:tdims_l%i_end,             &
                     tdims_l%j_start:tdims_l%j_end,             &
                     tdims_l%k_start:tdims_l%k_end,             &
                     super_array_size)

REAL, INTENT(OUT) :: co2                                      &
                  (tdims_s%i_start:tdims_s%i_end,               &
                   tdims_s%j_start:tdims_s%j_end,               &
                   tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: murk                                    &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: soot_new                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: soot_agd                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: soot_cld                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: so2                                     &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: so4_aitken                              &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: so4_accu                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: so4_diss                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: nh3                                     &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: dms                                     &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: dust_div1                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: dust_div2                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: dust_div3                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: dust_div4                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: dust_div5                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: dust_div6                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: bmass_new                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: bmass_agd                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: bmass_cld                               &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)  :: ocff_new                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: ocff_agd                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: ocff_cld                                &
                   (tdims_s%i_start:tdims_s%i_end,              &
                    tdims_s%j_start:tdims_s%j_end,              &
                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  :: nitr_acc                                &
                 (tdims_s%i_start:tdims_s%i_end,                &
                  tdims_s%j_start:tdims_s%j_end,                &
                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(OUT)  :: nitr_diss                               &
                 (tdims_s%i_start:tdims_s%i_end,                &
                  tdims_s%j_start:tdims_s%j_end,                &
                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)  ::                                         &
 tracers(tdims_s%i_start:tdims_s%i_end,                         &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end,tr_vars)           &
, tracer_ukca(tdims_s%i_start:tdims_s%i_end,                    &
              tdims_s%j_start:tdims_s%j_end,                    &
              tdims_s%k_start:tdims_s%k_end,tr_ukca)            &
! Add cariolle specific parameters for ozone tracer
            , ozone_tracer(tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end)

REAL ::                                                         &
  depart_xi1    (wdims%i_start:wdims%i_end,                     &
                 wdims%j_start:wdims%j_end,                     &
                 wdims%k_start:wdims%k_end)                     &
                                                ! Lambda
                                                ! co-ordinate of
                                                ! departure point.
, depart_xi2 (wdims%i_start:wdims%i_end,                        &
              wdims%j_start:wdims%j_end,                        &
              wdims%k_start:wdims%k_end)                        &
                                              ! Phi Co-ordinate of
                                              ! co-ordinate of
                                              ! departure point.
, depart_xi3(wdims%i_start:wdims%i_end,                         &
             wdims%j_start:wdims%j_end,                         &
             wdims%k_start:wdims%k_end)     ! Vertical
                                           ! co-ordinate of
                                           ! departure point.

REAL ::                                                           &
  eta_theta_levels(0:model_levels)

! Arguments with INTENT OUT. ie: Output variables.

INTEGER ::                                                        &
  Error_Code     ! Non-zero on exit if error detected.


REAL, INTENT(IN) ::   rho_n(pdims_s%i_start:pdims_s%i_end,        &
                            pdims_s%j_start:pdims_s%j_end,        &
                            pdims_s%k_start:pdims_s%k_end)
REAL, INTENT(IN) :: rho_np1(pdims_s%i_start:pdims_s%i_end,        &
                            pdims_s%j_start:pdims_s%j_end,        &
                            pdims_s%k_start:pdims_s%k_end)

! Local Variables.

! scalars

INTEGER ::                                                        &
  i, j, k                                                         &
              ! Loop indices
, temp                                                            &
, counter                                                         &
, tr_start

INTEGER :: k_int_linear ! Linear interpolation is used at departure
                        ! points in this layer and below.
                        ! (Optional argument for subroutine
                        !  eg_interpolation_eta_pmf.)
INTEGER :: number_of_inputs_zlf   ! number of tracers to be used with ZLF

! arrays

INTEGER ::  array_size_count

REAL ::                                                         &
  data_out_super(    tdims%i_start:tdims%i_end,                 &
                     tdims%j_start:tdims%j_end,                 &
                     tdims%k_start:tdims%k_end,                 &
                     super_array_size)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


CHARACTER(LEN=*), PARAMETER :: RoutineName='SL_TRACER1_4A'


! ----------------------------------------------------------------------
!  Section 0.    Initialise array_size_count
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
array_size_count = super_array_size


! ----------------------------------------------------------------------
! Call SWAPBOUNDS and set external halos for all arrays
! ----------------------------------------------------------------------

! DEPENDS ON: swap_bounds
CALL Swap_Bounds(                                                 &
       super_array,                                               &
       row_length, rows,                                          &
       tdims%k_len*array_size_count,                              &
       halo_i, halo_j, fld_type_p,  swap_field_is_scalar )

! For tracers without LBCs we need to do something:
! murk uses fill_external_halos (so already done by swap_bounds)
! to copy interior values to the exterior.
! For other tracers we set the halos to zeros, set_external_halo
! (this is why murk is the last field in the super_array).

IF (model_type == mt_lam) THEN

  IF (l_murk_advect .AND. .NOT.l_murk_lbc) THEN

    ! DEPENDS ON: set_external_halos
    CALL set_external_halos(                                    &
        super_array,                                            &
        row_length, rows,                                       &
        tdims_l%k_len*(array_size_count-1),                     &
        halo_i, halo_j, 0.0)
  ELSE

    ! DEPENDS ON: set_external_halos
    CALL set_external_halos(                                    &
        super_array,                                            &
        row_length, rows,                                       &
        tdims_l%k_len*array_size_count,                         &
        halo_i, halo_j, 0.0)

  END IF  ! L_Murk_advect

END IF

!------------------------------------------------------------
! having defined array now do interpolation
!------------------------------------------------------------


IF (super_array_size >= 1) THEN

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)
!$OMP DO SCHEDULE(STATIC)
  DO k = 1,super_array_size        
    DO j = tdims_l%j_start,tdims_l%j_end
      DO i = tdims_l%i_start, tdims_l%i_end
        super_array(i,j,0,k) =  super_array(i,j,1,k) 
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  IF ( zlf_conservation_tracers_option > 0 ) THEN
     !
     ! Here we can add or remove any tracer to the ZLF scheme 
     ! from 0 to all tracers (default here is treat all with zlf)
     !
     ! Note also non-ZLF tracers should be grouped first in the 
     !      super_array followed by the ZLF tracers
     !
     
    number_of_inputs_zlf   = super_array_size
     
    CALL eg_zlf_conservation(super_array, data_out_super,             &
          super_array_size, number_of_inputs_zlf,                     &
          row_length, rows, n_rows, model_levels, halo_i, halo_j,     &
          0, 0, datastart, g_i_pe, high_order_scheme,                 &
          monotone_scheme,  l_high, l_mono,                           &
          zlf_conservation_tracers_option, zlf_cfl_top_level_tracers, &
          rho_n, rho_np1, error_code                                  )
 
  ELSE
 
      ! Set layers over which linear interpolation is used
    IF (l_sl_bc_correction) THEN
      k_int_linear=2
    ELSE
      k_int_linear=1
    END IF

    CALL eg_interpolation_eta_pmf(                           &
            eta_theta_levels,fld_type_w,                     &
            super_array_size,                                &
            row_length, rows,model_levels+1,                 &
            rows,                                            &
            row_length, rows,model_levels+1,                 &
            high_order_scheme, monotone_scheme,              &
            l_high, l_mono,                                  &
            depart_xi3,depart_xi1,depart_xi2,                &
            mype, nproc, nproc_x, nproc_y,                   &
            halo_i, halo_j,                                  &
            global_row_length, datastart, at_extremity,      &
            g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
            0, 0, error_code,                                &
            super_array,data_out_super,                      &
            k_int_linear_in=k_int_linear)
  END IF  
   
END IF


! ----------------------------------------------------------------------
! Section 3.2.1  carbon cycle.
! ----------------------------------------------------------------------

array_size_count=0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, counter)

IF (l_CO2_interactive) THEN

!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        co2(i,j,k) = data_out_super(i,j,k,array_size_count)
      END DO
    END DO
  END DO
!$OMP END DO
END IF  ! l_CO2_interactive

! ----------------------------------------------------------------------
! Section 3.2.2.1  Soot cycle.
! ----------------------------------------------------------------------

IF (l_Soot) THEN

!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        soot_new(i,j,k) = data_out_super(i,j,k,array_size_count)
        soot_agd(i,j,k) = data_out_super(i,j,k,array_size_count+1)
        soot_cld(i,j,k) = data_out_super(i,j,k,array_size_count+2)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP MASTER
  array_size_count=array_size_count +2
!$OMP END MASTER
!$OMP BARRIER
END IF  ! l_soot

! ----------------------------------------------------------------------
! Section 3.2.2.2  Biomass aerosol.
! ----------------------------------------------------------------------
IF (l_Biomass) THEN

!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        bmass_new(i,j,k) = data_out_super(i,j,k,array_size_count)
        bmass_agd(i,j,k) = data_out_super(i,j,k,array_size_count+1)
        bmass_cld(i,j,k) = data_out_super(i,j,k,array_size_count+2)
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
! Section 3.2.3  sulphur cycle.
! ----------------------------------------------------------------------
IF (l_Sulpc_so2) THEN

!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        so4_aitken(i,j,k) = data_out_super(i,j,k,array_size_count)
        so4_accu(i,j,k) = data_out_super(i,j,k,array_size_count+1)
        so4_diss(i,j,k) = data_out_super(i,j,k,array_size_count+2)
        so2(i,j,k) = data_out_super(i,j,k,array_size_count+3)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP MASTER
  array_size_count=array_size_count +3
!$OMP END MASTER
!$OMP BARRIER

  IF (L_sulpc_nh3) THEN

!$OMP MASTER
    array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          nh3(i,j,k) = data_out_super(i,j,k,array_size_count)
        END DO
      END DO
    END DO
!$OMP END DO
  END IF
  IF (L_sulpc_dms) THEN

!$OMP MASTER
    array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          dms(i,j,k) = data_out_super(i,j,k,array_size_count)
        END DO
      END DO
    END DO
!$OMP END DO
  END IF

END IF  ! l_sulpc_so2

! ----------------------------------------------------------------------
! Section 3.2.4  mineral dust.
! ----------------------------------------------------------------------
IF (l_dust) THEN
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        dust_div1(i,j,k) = data_out_super(i,j,k,array_size_count)
        dust_div2(i,j,k) = data_out_super(i,j,k,array_size_count+1)
        IF (.NOT. l_twobin_dust) THEN
          dust_div3(i,j,k) = data_out_super(i,j,k,array_size_count+2)
          dust_div4(i,j,k) = data_out_super(i,j,k,array_size_count+3)
          dust_div5(i,j,k) = data_out_super(i,j,k,array_size_count+4)
          dust_div6(i,j,k) = data_out_super(i,j,k,array_size_count+5)
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
END IF ! L_DUST

! ----------------------------------------------------------------------
! New addition  Fossil-fuel organic carbon aerosol.
! ----------------------------------------------------------------------
IF (L_ocff) THEN

!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        ocff_new(i,j,k) = data_out_super(i,j,k,array_size_count)
        ocff_agd(i,j,k) = data_out_super(i,j,k,array_size_count+1)
        ocff_cld(i,j,k) = data_out_super(i,j,k,array_size_count+2)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP MASTER
  array_size_count=array_size_count +2
!$OMP END MASTER
!$OMP BARRIER

END IF  ! L_ocff

! ----------------------------------------------------------------------
! Section 3.2.4.1  Cariolle ozone tracer.
! ----------------------------------------------------------------------
IF (l_use_cariolle) THEN

!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        ozone_tracer(i,j,k) = data_out_super(i,j,k,array_size_count)
      END DO
    END DO
  END DO
!$OMP END DO
END IF ! L_USE_CARIOLLE

! ----------------------------------------------------------------------
! New addition  Ammonium nitrate aerosol
! ----------------------------------------------------------------------
IF (L_nitrate) THEN

!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        nitr_acc(i,j,k) = data_out_super(i,j,k,array_size_count)
        nitr_diss(i,j,k)= data_out_super(i,j,k,array_size_count+1)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

END IF  ! L_nitrate

! ----------------------------------------------------------------------
! ANY NEW NAMED AEROSOL SPECIES SHOULD BE ADDED HERE ^^^
! and in the same location/order in tr_set_phys
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 3.2.5  free tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
IF ( tr_vars > 0 ) THEN

  DO counter=1,tr_vars
!$OMP MASTER
    array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          tracers(i,j,k,counter) = data_out_super(i,j,k,array_size_count)
        END DO
      END DO
    END DO
!$OMP END DO 
  END DO
END IF  ! tr_vars > 0

! ----------------------------------------------------------------------
! Section 3.2.6  UKCA tracers  (tr_levels = model_levels)
! ----------------------------------------------------------------------
IF ( tr_ukca > 0 ) THEN

  DO counter = 1, tr_ukca
!$OMP MASTER
    array_size_count=array_size_count +1
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          tracer_ukca(i,j,k,counter) =                             &
                           data_out_super(i,j,k,array_size_count)
        END DO
      END DO
    END DO  ! Loop over levels
!$OMP END DO 
  END DO !  counter = 1, tr_ukca
END IF  ! tr_ukca > 0

! ----------------------------------------------------------------------
! Section 3.2.99  Murk cycle. This must be the last Full level field in
!                            the super_array
! ----------------------------------------------------------------------

IF (L_Murk_advect) THEN

!$OMP MASTER
  array_size_count=array_size_count +1
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        murk(i,j,k) = data_out_super(i,j,k,array_size_count)
      END DO
    END DO
  END DO
!$OMP END DO
END IF  ! L_Murk_advect

!$OMP END PARALLEL

! END of routine.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE SL_Tracer1_4A

END MODULE SL_Tracer1_mod

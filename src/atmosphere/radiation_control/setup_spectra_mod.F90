! *****************************COPYRIGHT*******************************
!(C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE setup_spectra_mod

IMPLICIT NONE

!
!  Subroutine SETUP_SPECTRA
!
! Purpose:
!   Takes the data from spectral files read in on processor 0 and   
!   broadcasts it to all the other processors
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
!*
!* CAUTION - This routine needs to mirror any change made in the 
!*           def_spectrum module: 
!* see the Reading_namelists wiki page for guidance.
!* 
!- ---------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SETUP_SPECTRA_MOD'
CONTAINS

SUBROUTINE setup_spectra(con, spec, n_call)

USE def_control,    ONLY: StrCtrl
USE def_spectrum,   ONLY: StrSpecData, n_dim, allocate_spectrum
USE um_parcore,     ONLY: mype
USE setup_namelist, ONLY: setup_nml_type
USE mpl,            ONLY: mpl_packed, mpl_integer, mpl_real, mpl_character,    &
                          mpl_logical     
USE filenamelength_mod, ONLY: filenamelength

IMPLICIT NONE

INTEGER, INTENT(IN)  :: n_call
TYPE(StrSpecData)    :: spec(n_call)
TYPE(StrCtrl)        :: con(n_call)

!  local variables
INTEGER              :: i
INTEGER, PARAMETER   :: no_of_types = 1
INTEGER              :: my_comm
INTEGER              :: icode
INTEGER              :: mpl_spectraldims_type
INTEGER              :: POSITION
INTEGER              :: count_int
INTEGER              :: count_real
INTEGER              :: count_log
INTEGER              :: k1
INTEGER              :: k2
INTEGER              :: k3
INTEGER              :: k
CHARACTER, ALLOCATABLE :: buffer(:)

CALL gc_get_communicator(my_comm, icode)

! the line to setup the spectral file names is now inside an mype==0 
! block in initphys so is done in this subroutine for all other pe's
CALL mpl_bcast(con%spectral_file, filenamelength, mpl_character,  &
             0, my_comm, icode)

! broadcast spectral file dimensions to rest of pe's
CALL setup_nml_type(no_of_types, mpl_spectraldims_type, n_int_in=n_dim)
CALL mpl_bcast(spec(:)%dim, n_call, mpl_spectraldims_type, 0, my_comm, icode)

! allocate arrays on rest of pe's
IF (mype /= 0) THEN
  DO i = 1, n_call
    CALL allocate_spectrum(spec(i))
  END DO
END IF !!mype/=0

CALL mpl_type_free(mpl_spectraldims_type,icode)

! now arrays have been allocated, can calculate size of buffer required
! need to calculate number of each type of element

count_log = 0
count_int = 0
count_real = 0

DO i = 1, n_call
  count_log  = count_log  + spec(i)%Dim%nd_alloc_log
  count_int  = count_int  + spec(i)%Dim%nd_alloc_int
  count_real = count_real + spec(i)%Dim%nd_alloc_real
END DO
 
CALL mpl_pack_size(count_log, mpl_logical, my_comm, k1, icode)
CALL mpl_pack_size(count_int, mpl_integer, my_comm, k2, icode)
CALL mpl_pack_size(count_real, mpl_real, my_comm, k3, icode)
k=k1+k2+k3
k= INT( (11*k) /10 )    !! add 10% to buffer for safety

ALLOCATE(buffer(k))

! pack spectra data on pe 0 into buffer

POSITION=0

IF ( mype == 0 ) THEN
  DO i = 1, n_call
    CALL mpl_pack(spec(i)%Basic% l_present, SIZE(spec(i)%Basic%l_present),   &
                 mpl_logical, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Basic%n_band, 1, mpl_integer, buffer, k, POSITION, &
                 my_comm, icode)
    CALL mpl_pack(spec(i)%Basic%wavelength_long , spec(i)%Dim%nd_band,       &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Basic%wavelength_short , spec(i)%Dim%nd_band,      &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Basic%n_band_exclude, spec(i)%Dim%nd_band,         &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Basic%index_exclude,                               &
                  SIZE(spec(i)%Basic%index_exclude),                         &
                  mpl_integer, buffer, k, POSITION, my_comm, icode)

    CALL mpl_pack(spec(i)%Solar%solar_flux_band, spec(i)%Dim%nd_band,        &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Solar%solar_flux_band_ses,                         &
                  SIZE(spec(i)%Solar%solar_flux_band_ses),                   &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Solar%weight_blue, spec(i)%Dim%nd_band,            &
                 mpl_real, buffer, k, POSITION, my_comm, icode)

    CALL mpl_pack(spec(i)%Rayleigh%i_rayleigh_scheme, 1,                     &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Rayleigh%rayleigh_coeff, spec(i)%Dim%nd_band,      &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Rayleigh%n_gas_rayleigh, 1,                        &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Rayleigh%index_rayleigh, spec(i)%Dim%nd_species,   &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Rayleigh%rayleigh_coeff_gas,                       &
                 SIZE(spec(i)%Rayleigh%rayleigh_coeff_gas),                  &
                 mpl_real, buffer, k, POSITION, my_comm, icode)

    CALL mpl_pack(spec(i)%Gas%n_absorb, 1,                                   &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%n_absorb_sb, 1,                                &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%n_gas_frac, 1,                                 &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%n_band_absorb, spec(i)%Dim%nd_band,            &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%index_absorb, SIZE(spec(i)%Gas%index_absorb),  &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%index_sb, SIZE(spec(i)%Gas%index_sb),          &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%type_absorb, spec(i)%Dim%nd_species,           &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%n_mix_gas, spec(i)%Dim%nd_band,                &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%index_mix_gas,                                 &
                 SIZE(spec(i)%Gas%index_mix_gas),                            &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%num_mix, spec(i)%Dim%nd_band,                  &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%mix_gas_band, spec(i)%Dim%nd_band,             &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%num_ref_p, SIZE(spec(i)%Gas%num_ref_p),        &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%num_ref_t, SIZE(spec(i)%Gas%num_ref_t),        &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%i_band_k, SIZE(spec(i)%Gas%i_band_k),          &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%i_band_k_ses, spec(i)%Dim%nd_band,             &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%i_scale_k, SIZE(spec(i)%Gas%i_scale_k),        &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%i_scale_fnc, SIZE(spec(i)%Gas%i_scale_fnc),    &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%i_scat, SIZE(spec(i)%Gas%i_scat),              &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%l_self_broadening,                             &
             SIZE(spec(i)%Gas%l_self_broadening),                            &
                 mpl_logical, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%k, SIZE(spec(i)%Gas%k),                        &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%w, SIZE(spec(i)%Gas%w),                        &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%scale, SIZE(spec(i)%Gas%scale),                &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%p_ref, SIZE(spec(i)%Gas%p_ref),                &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%t_ref, SIZE(spec(i)%Gas%t_ref),                &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%p_lookup, spec(i)%Dim%nd_pre,                  &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%t_lookup, SIZE(spec(i)%Gas%t_lookup),          &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%gf_lookup, SIZE(spec(i)%Gas%gf_lookup),        &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%k_lookup, SIZE(spec(i)%Gas%k_lookup),          &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%k_lookup_sb, SIZE(spec(i)%Gas%k_lookup_sb),    &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%w_ses, SIZE(spec(i)%Gas%w_ses),                &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%k_mix_gas, SIZE(spec(i)%Gas%k_mix_gas),        &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%f_mix, spec(i)%Dim%nd_band,                    &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%l_doppler, spec(i)%Dim%nd_species,             &
                 mpl_logical, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Gas%doppler_cor, spec(i)%Dim%nd_species,           &
                 mpl_real, buffer, k, POSITION, my_comm, icode)

    CALL mpl_pack(spec(i)%Planck%n_deg_fit, 1, mpl_integer, buffer, k,       &
                 POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Planck% thermal_coeff,                             &
                 SIZE(spec(i)%Planck% thermal_coeff)  ,                      &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Planck% theta_planck_tbl,                          &
                 SIZE(spec(i)%Planck% theta_planck_tbl),                     &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Planck%t_ref_planck, 1, mpl_real, buffer, k,       &
                 POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Planck%l_planck_tbl, 1, mpl_logical, buffer, k,    &
                 POSITION, my_comm, icode)

    CALL mpl_pack(spec(i)%Cont%n_band_continuum, spec(i)%Dim%nd_band,        &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Cont%index_continuum ,                             &
                 SIZE(spec(i)%Cont%index_continuum),                         &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Cont%index_water, 1,                               &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Cont%i_scale_fnc_cont,                             &
                 SIZE(spec(i)%Cont%i_scale_fnc_cont),                        &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Cont%k_cont, SIZE(spec(i)%Cont%k_cont),            &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Cont%scale_cont, SIZE(spec(i)%Cont%scale_cont),    &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Cont%p_ref_cont, SIZE(spec(i)%Cont%p_ref_cont),    &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Cont%t_ref_cont, SIZE(spec(i)%Cont%t_ref_cont),    &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Cont%k_cont_ses, SIZE(spec(i)%Cont%k_cont_ses),    &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Cont%k_h2oc, SIZE(spec(i)%Cont%k_h2oc),            &
                 mpl_real, buffer, k, POSITION, my_comm, icode)

    CALL mpl_pack(spec(i)%ContGen%n_cont, 1,                                 &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%n_band_cont,                               &
             SIZE(spec(i)%ContGen%n_band_cont),                              &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%index_cont,                                &
             SIZE(spec(i)%ContGen%index_cont),                               &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%index_cont_gas_1,                          &
             SIZE(spec(i)%ContGen%index_cont_gas_1),                         &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%index_cont_gas_2,                          &
             SIZE(spec(i)%ContGen%index_cont_gas_2),                         &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%i_band_k_cont,                             &
             SIZE(spec(i)%ContGen%i_band_k_cont),                            &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%i_cont_overlap_band,                       &
             SIZE(spec(i)%ContGen%i_cont_overlap_band),                      &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%i_scat_cont,                               &
             SIZE(spec(i)%ContGen%i_scat_cont),                              &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%l_cont_major,                              &
             SIZE(spec(i)%ContGen%l_cont_major),                             &
                 mpl_logical, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%k_cont,                                    &
             SIZE(spec(i)%ContGen%k_cont),                                   &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%w_cont,                                    &
             SIZE(spec(i)%ContGen%w_cont),                                   &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%n_t_lookup_cont, 1,                        &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%t_lookup_cont,                             &
             SIZE(spec(i)%ContGen%t_lookup_cont),                            &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%ContGen%k_lookup_cont,                             &
             SIZE(spec(i)%ContGen%k_lookup_cont),                            &
                 mpl_real, buffer, k, POSITION, my_comm, icode)


    CALL mpl_pack(spec(i)%Drop%l_drop_type, spec(i)%Dim%nd_drop_type,        &
                 mpl_logical, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Drop%i_drop_parm, spec(i)%Dim%nd_drop_type,        &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Drop%n_phf, spec(i)%Dim%nd_drop_type,              &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Drop%parm_list, SIZE(spec(i)%Drop% parm_list),     &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Drop%parm_min_dim, spec(i)%Dim%nd_drop_type,       &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Drop%parm_max_dim, spec(i)%Dim%nd_drop_type,       &
                 mpl_real, buffer, k, POSITION, my_comm, icode)

    CALL mpl_pack(spec(i)%Aerosol%l_aero_spec,                               &
                  spec(i)%Dim%nd_aerosol_species,                            &
                 mpl_logical, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%n_aerosol, 1,                              &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%n_aerosol_mr, 1,                           &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%type_aerosol,                              &
                  spec(i)%Dim%nd_aerosol_species,                            &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%i_aerosol_parm,                            &
                  spec(i)%Dim%nd_aerosol_species,                            &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%n_aerosol_phf_term,                        &
                  spec(i)%Dim%nd_aerosol_species,                            &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%nhumidity,                                 &
                  spec(i)%Dim%nd_aerosol_species,                            &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%abs, SIZE(spec(i)%Aerosol%abs),            &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%scat, SIZE(spec(i)%Aerosol%scat),          &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%phf_fnc, SIZE(spec(i)%Aerosol%phf_fnc),    &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%humidities,                                &
                 SIZE(spec(i)%Aerosol%humidities),                           &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%n_aod_wavel, 1,                            &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%i_aod_type,                                &
                  spec(i)%Dim%nd_aerosol_species,                            &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%aod_wavel, spec(i)%Dim%nd_aod_wavel,       &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%aod_abs, SIZE(spec(i)%Aerosol%aod_abs),    &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Aerosol%aod_scat, SIZE(spec(i)%Aerosol%aod_scat),  &
                 mpl_real, buffer, k, POSITION, my_comm, icode)

    CALL mpl_pack(spec(i)%Ice%l_ice_type, spec(i)%Dim%nd_ice_type,           &
                 mpl_logical, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Ice%i_ice_parm, spec(i)%Dim%nd_ice_type,           &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Ice%n_phf, spec(i)%Dim%nd_ice_type,                &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Ice%parm_list, SIZE(spec(i)%Ice% parm_list),       &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Ice%parm_min_dim, spec(i)%Dim%nd_ice_type,         &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Ice%parm_max_dim, spec(i)%Dim%nd_ice_type,         &
                 mpl_real, buffer, k, POSITION, my_comm, icode)

    CALL mpl_pack(spec(i)%Var%n_sub_band, 1,                                 &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Var%n_times, 1,                                    &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Var%n_repeat_times, 1,                             &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Var%n_rayleigh_coeff, 1,                           &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Var%index_sub_band,                                &
                 SIZE(spec(i)%Var%index_sub_band),                           &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Var%wavelength_sub_band,                           &
                 SIZE(spec(i)%Var%wavelength_sub_band),                      &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Var%time, SIZE(spec(i)%Var%time),                  &
                 mpl_integer, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Var%total_solar_flux,                              &
                 SIZE(spec(i)%Var%total_solar_flux),                         &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Var%solar_flux_sub_band,                           &
                 SIZE(spec(i)%Var%solar_flux_sub_band),                      &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
    CALL mpl_pack(spec(i)%Var%rayleigh_coeff,                                &
                 SIZE(spec(i)%Var%rayleigh_coeff),                           &
                 mpl_real, buffer, k, POSITION, my_comm, icode)
  END DO
END IF  ! mype==0

! broadcast buffer
CALL mpl_bcast(buffer, k, mpl_packed, 0, my_comm, icode)

! unpack data on pe/=0
IF ( mype /= 0 ) THEN
  POSITION=0
  DO i = 1, n_call
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Basic%l_present,            &
          SIZE(spec(i)%Basic% l_present), mpl_logical, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Basic%n_band,               &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Basic%wavelength_long,      &
          spec(i)%Dim%nd_band, mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Basic%wavelength_short,     &
          spec(i)%Dim%nd_band, mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Basic%n_band_exclude,       &
          spec(i)%Dim%nd_band, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Basic%index_exclude,        &
          SIZE(spec(i)%Basic%index_exclude), mpl_integer, my_comm, icode)

    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Solar% solar_flux_band,     &
          spec(i)%Dim%nd_band, mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Solar% solar_flux_band_ses, &
          SIZE(spec(i)%Solar%solar_flux_band_ses), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Solar% weight_blue,         &
          spec(i)%Dim%nd_band, mpl_real, my_comm, icode)

    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Rayleigh%i_rayleigh_scheme, &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Rayleigh%rayleigh_coeff,    &
          spec(i)%Dim%nd_band, mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Rayleigh%n_gas_rayleigh,    &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Rayleigh%index_rayleigh,    &
          spec(i)%Dim%nd_species, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Rayleigh%rayleigh_coeff_gas,&
          SIZE(spec(i)%Rayleigh%rayleigh_coeff_gas), mpl_real, my_comm, icode)

    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%n_absorb,               &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%n_absorb_sb,            &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%n_gas_frac,             &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%n_band_absorb,          &
          spec(i)%Dim%nd_band, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%index_absorb,           &
          SIZE(spec(i)%Gas%index_absorb), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%index_sb,               &
          SIZE(spec(i)%Gas%index_sb), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%type_absorb,            &
          spec(i)%Dim%nd_species, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%n_mix_gas,              &
          spec(i)%Dim%nd_band, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%index_mix_gas,          &
          SIZE(spec(i)%Gas%index_mix_gas), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%num_mix,                &
          spec(i)%Dim%nd_band, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%mix_gas_band,           &
          spec(i)%Dim%nd_band, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%num_ref_p,              &
          SIZE(spec(i)%Gas%num_ref_p), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%num_ref_t,              &
          SIZE(spec(i)%Gas%num_ref_t), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%i_band_k,               &
          SIZE(spec(i)%Gas%i_band_k), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%i_band_k_ses,           &
          spec(i)%Dim%nd_band, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%i_scale_k,              &
          SIZE(spec(i)%Gas%i_scale_k), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%i_scale_fnc,            &
          SIZE(spec(i)%Gas%i_scale_fnc), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%i_scat,                 &
          SIZE(spec(i)%Gas%i_scat), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%l_self_broadening,      &
          SIZE(spec(i)%Gas%l_self_broadening), mpl_logical, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%k,                      &
          SIZE(spec(i)%Gas%k), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%w,                      &
          SIZE(spec(i)%Gas%w), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%scale,                  &
          SIZE(spec(i)%Gas%scale), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%p_ref,                  &
          SIZE(spec(i)%Gas%p_ref), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%t_ref,                  &
          SIZE(spec(i)%Gas%t_ref), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%p_lookup,               &
          spec(i)%Dim%nd_pre, mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%t_lookup,               &
          SIZE(spec(i)%Gas%t_lookup), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%gf_lookup,              &
          SIZE(spec(i)%Gas%gf_lookup), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%k_lookup,               &
          SIZE(spec(i)%Gas%k_lookup), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%k_lookup_sb,            &
          SIZE(spec(i)%Gas%k_lookup_sb), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%w_ses,                  &
          SIZE(spec(i)%Gas%w_ses), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%k_mix_gas,              &
          SIZE(spec(i)%Gas%k_mix_gas), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%f_mix,                  &
          spec(i)%Dim%nd_band, mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%l_doppler,              &
          spec(i)%Dim%nd_species, mpl_logical, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Gas%doppler_cor,            &
          spec(i)%Dim%nd_species, mpl_real, my_comm, icode)

    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Planck%n_deg_fit,           &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Planck%thermal_coeff,       &
          SIZE(spec(i)%Planck% thermal_coeff), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Planck%theta_planck_tbl,    &
          SIZE(spec(i)%Planck% theta_planck_tbl), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Planck%t_ref_planck,        &
          1, mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Planck%l_planck_tbl,        &
          1, mpl_logical, my_comm, icode)

    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Cont%n_band_continuum,      &
          spec(i)%Dim%nd_band, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Cont%index_continuum ,      &
          SIZE(spec(i)%Cont%index_continuum), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Cont%index_water,           &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Cont%i_scale_fnc_cont,      &
          SIZE(spec(i)%Cont%i_scale_fnc_cont), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Cont%k_cont,                &
          SIZE(spec(i)%Cont%k_cont), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Cont%scale_cont,            &
          SIZE(spec(i)%Cont%scale_cont), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Cont%p_ref_cont,            &
          SIZE(spec(i)%Cont%p_ref_cont), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Cont%t_ref_cont,            &
          SIZE(spec(i)%Cont%t_ref_cont), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Cont%k_cont_ses,            &
          SIZE(spec(i)%Cont%k_cont_ses), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Cont%k_h2oc,                &
          SIZE(spec(i)%Cont%k_h2oc), mpl_real, my_comm, icode)

    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%n_cont,             &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%n_band_cont,        &
          SIZE(spec(i)%ContGen%n_band_cont), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%index_cont,         &
          SIZE(spec(i)%ContGen%index_cont), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%index_cont_gas_1,   &
          SIZE(spec(i)%ContGen%index_cont_gas_1), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%index_cont_gas_2,   &
          SIZE(spec(i)%ContGen%index_cont_gas_2), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%i_band_k_cont,      &
          SIZE(spec(i)%ContGen%i_band_k_cont), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%i_cont_overlap_band,&
          SIZE(spec(i)%ContGen%i_cont_overlap_band), mpl_integer, my_comm,icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%i_scat_cont,        &
          SIZE(spec(i)%ContGen%i_scat_cont), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%l_cont_major,       &
          SIZE(spec(i)%ContGen%l_cont_major), mpl_logical, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%k_cont,             &
          SIZE(spec(i)%ContGen%k_cont), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%w_cont,             &
          SIZE(spec(i)%ContGen%w_cont), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%n_t_lookup_cont,    &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%t_lookup_cont,      &
          SIZE(spec(i)%ContGen%t_lookup_cont), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%ContGen%k_lookup_cont,      &
          SIZE(spec(i)%ContGen%k_lookup_cont), mpl_real, my_comm, icode)

    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Drop%l_drop_type,           &
          spec(i)%Dim%nd_drop_type, mpl_logical, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Drop%i_drop_parm,           &
          spec(i)%Dim%nd_drop_type, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Drop%n_phf,                 &
          spec(i)%Dim%nd_drop_type, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Drop%parm_list,             &
          SIZE(spec(i)%Drop%parm_list), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Drop%parm_min_dim,          &
          spec(i)%Dim%nd_drop_type, mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Drop%parm_max_dim,          &
          spec(i)%Dim%nd_drop_type, mpl_real, my_comm, icode)

    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%l_aero_spec,        &
          spec(i)%Dim%nd_aerosol_species, mpl_logical, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%n_aerosol,          &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%n_aerosol_mr,       &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%type_aerosol,       &
          spec(i)%Dim%nd_aerosol_species, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%i_aerosol_parm ,    &
          spec(i)%Dim%nd_aerosol_species, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%n_aerosol_phf_term, &
          spec(i)%Dim%nd_aerosol_species, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%nhumidity,          &
          spec(i)%Dim%nd_aerosol_species, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%abs,                &
          SIZE(spec(i)%Aerosol%abs), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%scat,               &
          SIZE(spec(i)%Aerosol%scat), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%phf_fnc,            &
          SIZE(spec(i)%Aerosol%phf_fnc), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%humidities,         &
          SIZE(spec(i)%Aerosol%humidities), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%n_aod_wavel,        &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%i_aod_type,         &
          spec(i)%Dim%nd_aerosol_species, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%aod_wavel,          &
          spec(i)%Dim%nd_aod_wavel, mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%aod_abs,            &
          SIZE(spec(i)%Aerosol%aod_abs), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Aerosol%aod_scat,           &
          SIZE(spec(i)%Aerosol%aod_scat), mpl_real, my_comm, icode)

    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Ice%l_ice_type,             &
          spec(i)%Dim%nd_ice_type, mpl_logical, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Ice%i_ice_parm,             &
          spec(i)%Dim%nd_ice_type, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Ice%n_phf,                  &
          spec(i)%Dim%nd_ice_type, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Ice%parm_list,              &
          SIZE(spec(i)%Ice% parm_list), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Ice%parm_min_dim,           &
          spec(i)%Dim%nd_ice_type, mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Ice%parm_max_dim,           &
          spec(i)%Dim%nd_ice_type, mpl_real, my_comm, icode)

    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Var%n_sub_band,             &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Var%n_times,                &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Var%n_repeat_times,         &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Var%n_rayleigh_coeff,       &
          1, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Var%index_sub_band,         &
          SIZE(spec(i)%Var%index_sub_band), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Var%wavelength_sub_band,    &
          SIZE(spec(i)%Var%wavelength_sub_band), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Var%time,                   &
          SIZE(spec(i)%Var%time), mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Var%total_solar_flux,       &
          SIZE(spec(i)%Var%total_solar_flux), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Var%solar_flux_sub_band,    &
          SIZE(spec(i)%Var%solar_flux_sub_band), mpl_real, my_comm, icode)
    CALL mpl_unpack(buffer, k, POSITION, spec(i)%Var%rayleigh_coeff,         &
          SIZE(spec(i)%Var%rayleigh_coeff), mpl_real, my_comm, icode)
  END DO
END IF !!mype/=0

DEALLOCATE(buffer)

END SUBROUTINE setup_spectra

END MODULE setup_spectra_mod

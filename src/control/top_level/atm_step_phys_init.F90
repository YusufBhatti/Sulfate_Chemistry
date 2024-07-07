! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Contains various chunks of code from atm_step - the purpose of each
! section is indicated at the head of the section
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Subroutine Interface:
SUBROUTINE atm_step_phys_init( &
r_v, r_u, theta_star, q_star, qcl_star, qcf_star, cf_star,  &
cfl_star,   cff_star, exner_lbc_real_tend, w_lbc_real_tend, &
errorstatus, flag )

USE jules_surface_types_mod
USE atm_step_local
USE atm_fields_bounds_mod, ONLY : o3dims2,                              &
                                  tdims, tdims_s, udims, udims_s,       &
                                  vdims, vdims_s, wdims_s
USE turb_diff_mod, ONLY: L_subfilter_horiz, L_subfilter_vert
USE turb_diff_ctl_mod, ONLY:                                            &
    visc_m, visc_h, rneutml_sq, shear, delta_smag, max_diff

USE mpp_conf_mod,  ONLY: swap_field_is_scalar

USE jules_snow_mod, ONLY: nsmax
USE timestep_mod

USE rad_input_mod, ONLY:                                                &
    l_use_seasalt_indirect,l_use_seasalt_direct,                        &
    i_ozone_int, l_use_cariolle, l_use_ozoneinrad, lexpand_ozone
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE UM_ParVars
USE UM_ParCore,   ONLY: nproc
USE UM_ParParams, ONLY: halo_type_extended
USE Field_Types,  ONLY: fld_type_p
USE Control_Max_Sizes
USE rimtypes
USE lbc_mod
USE o3intp_mod, ONLY: io3_3dspec, io3_2dspec, io3_2dmasscon,            &
                      io3_trop_map, io3_trop_map_masscon
USE problem_mod, ONLY: standard, dynamical_core,                        &
                       idealised_problem, idealised_planet

USE dynamics_testing_mod, ONLY:  L_Backwards, problem_number
USE dyncore_ctl_mod, ONLY: friction_level
USE idealise_run_mod, ONLY: base_frictional_timescale,                  &
                            suhe_sigma_cutoff, suhe_fric

USE carbon_options_mod,  ONLY: l_co2_interactive
USE run_aerosol_mod, ONLY: l_dms, l_dms_em, l_dms_em_inter,  &
     l_use_seasalt_sulpc

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain, &
                            l_use_seasalt_autoconv
USE stash_array_mod, ONLY: sf_calc
USE jules_vegetation_mod, ONLY: l_triffid

USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first

USE jules_sea_seaice_mod, ONLY: nice, nice_use

USE nlsizes_namelist_mod, ONLY:                                            &
    a_len_inthd, a_len_realhd, bl_levels, global_row_length, land_field,   &
    len1_lbc_comp_lookup, len1_lookup, len_fixhd, model_levels, n_cca_lev, &
    n_rows, ntiles, river_row_length, river_rows,                          &
    row_length, rows, sm_levels, theta_off_size, tr_levels, tr_ukca, tr_vars

USE atm_fields_mod, ONLY: exner_lbc, exner_lbc_tend, ozone_tracer, w_lbc, &
    w_lbc_tend, p, pstar, u, v, soil_carb1, rsp_s_acc1, soil_carb, rsp_s_acc

USE atm_boundary_headers_mod, ONLY: rim_stepsa

IMPLICIT NONE


! Subroutine arguments

REAL, TARGET :: r_u(udims_s%i_start:udims_s%i_end,                     &
                    udims_s%j_start:udims_s%j_end,                     &
                    udims_s%k_start:udims_s%k_end)
REAL, TARGET :: r_v(vdims_s%i_start:vdims_s%i_end,                     &
                    vdims_s%j_start:vdims_s%j_end,                     &
                    vdims_s%k_start:vdims_s%k_end)

REAL :: theta_star(tdims%i_start:tdims%i_end,                          &
                   tdims%j_start:tdims%j_end,                          &
                   tdims%k_start:tdims%k_end)
REAL :: q_star  (tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)
REAL :: qcl_star(tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)
REAL :: qcf_star(tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)

REAL :: cf_star (tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)
REAL :: cfl_star(tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)
REAL :: cff_star(tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)

REAL :: exner_lbc_real_tend(lenrima(fld_type_p,halo_type_extended,     &
                                    rima_type_norm),model_levels+1)
REAL :: w_lbc_real_tend (lenrima(fld_type_p,halo_type_extended,     &
                                 rima_type_norm),wdims_s%k_start:wdims_s%k_end)

INTEGER :: errorstatus

! Local variables

REAL :: lbc_fac

CHARACTER(LEN=8) :: flag

REAL :: temp_ozone(tdims_s%i_start:tdims_s%i_end,                      &
                   tdims_s%j_start:tdims_s%j_end,                      &
                   tdims_s%k_start:tdims_s%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATM_STEP_PHYS_INIT'

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (flag == 'uvw_btnd') THEN

  ! Obtain tendencies in the boundary zone for u, v, w and use to make.
  ! lbcs for u_adv, v_adv, w_adv

  l_do_halos=.TRUE.
  l_do_boundaries=.TRUE.

  IF (rim_stepsa  ==  0) THEN
    increment_factor=0.0
  ELSE
    increment_factor=1.0/ (rim_stepsa-MOD(timestep_number-1,rim_stepsa))
  END IF

  lbc_size=lenrima(fld_type_p,halo_type_extended, rima_type_norm)

  lbc_fac = 0.5

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,k) SCHEDULE(STATIC)    &
!$OMP SHARED(wdims_s,lbc_size,w_lbc_real_tend,w_lbc,lbc_fac,     &
!$OMP        increment_factor,w_lbc_tend)
  DO k=wdims_s%k_start,wdims_s%k_end
    DO i=1,lbc_size
      w_lbc_real_tend(i,k) = w_lbc(i,k) +  lbc_fac * increment_factor *&
                                ( w_lbc_tend(i,k) - w_lbc(i,k) )
    END DO
  END DO
!$OMP END PARALLEL DO

  ! ----------------------------------------------------------------------

ELSE IF (flag == 'zeroincs') THEN

  ! initialise arrays that hold physics increments to zero, only needs
  ! doing at non-halo points
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                     &
!$OMP SHARED(tdims,theta_star,q_star,qcl_star,qcf_star,cf_star, &
!$OMP cfl_star,cff_star,l_mcr_qcf2,qcf2_star,l_mcr_qrain,qrain_star, &
!$OMP l_mcr_qgraup,qgraup_star,udims,r_u,vdims,r_v)

!$OMP DO SCHEDULE(STATIC)
  DO k=tdims%k_start, tdims%k_end
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        theta_star(i,j,k) = 0.0
        q_star(i,j,k) = 0.0
        qcl_star(i,j,k) = 0.0
        qcf_star(i,j,k) = 0.0
        cf_star(i,j,k) = 0.0
        cfl_star(i,j,k) = 0.0
        cff_star(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  ! Initialise additional microphysics variables if in use
  IF (l_mcr_qcf2) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k=tdims%k_start, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          qcf2_star(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_mcr_qrain) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k=tdims%k_start, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          qrain_star(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_mcr_qgraup) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k=tdims%k_start, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          qgraup_star(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

!$OMP DO SCHEDULE(STATIC)
  DO k=udims%k_start, udims%k_end
    DO j=udims%j_start, udims%j_end
      DO i=udims%i_start, udims%i_end
        r_u(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k=vdims%k_start, vdims%k_end
    DO j=vdims%j_start, vdims%j_end
      DO i=vdims%i_start, vdims%i_end
        r_v(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  ! -----------------------------------------------------------


ELSE IF (flag == 'friction') THEN

  ! Run with a positive timestep if integrating backwards.
  IF (l_backwards) timestep = pos_timestep
  !  When no physics applied then may need a simple friction
  !   if non-inviscid (idealised_problem) then DO NOT apply friction
  !   if we are modeling an idealised planet we also do not apply 
  !   friction.

  IF (problem_number  /=  idealised_problem .AND.                   &
       problem_number /= idealised_planet) THEN
    !  standard simple friction being used

    IF (problem_number  ==  dynamical_core) THEN

      ! DEPENDS ON: idl_friction_suarez_held
      CALL idl_friction_suarez_held(row_length, rows, n_rows,       &
                          model_levels, timestep,                   &
                          offx, offy,                               &
                          global_row_length,nproc, nproc_y,         &
                          gc_proc_row_group,at_extremity,           &
                          friction_level,                           &
                          base_frictional_timescale,                &
                          suhe_sigma_cutoff, suhe_fric,             &
                          p, pstar, u, v, r_u, r_v )

    ELSE     ! problem_number  /=  dynamical_core

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SCHEDULE(STATIC)   &
!$OMP SHARED(bl_levels,j_begin,j_end,udims,r_u,timestep,friction_level, &
!$OMP u,vdims,r_v,v)
      DO k = 1, bl_levels
        DO j = j_begin, j_end
          DO i = udims%i_start, udims%i_end
            r_u(i,j,k) = r_u(i,j,k) -timestep *friction_level(k) *u(i,j,k)
          END DO
        END DO

        DO j = vdims%j_start, vdims%j_end
          DO i = vdims%i_start, vdims%i_end
            r_v(i,j,k) = r_v(i,j,k) -timestep *friction_level(k) *v(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

    END IF    ! problem_number  ==  dynamical_core

  END IF    ! problem_number  /=  idealised_problem

  ! Go back to negative timestep if integrating backwards.
  IF (l_backwards) timestep = neg_timestep


  ! ----------------------------------------------------------------------


ELSE IF (flag == 'microphy') THEN

  ! Microphysics and ozone code from atm_step

    ! Set LAND_PTS_TRIF and NPFT_TRIF according to TRIFFID on/off
    ! and set pointers to single or multi pool soil carbon, and
    ! soil carbon dimensions to be used in subroutines.
  IF (l_triffid) THEN
    land_pts_trif = land_field
    npft_trif = npft
    cs => soil_carb1
    rsa => rsp_s_acc1
    dim_cs1 = 4
    dim_cs2 = land_field
  ELSE
    land_pts_trif = 1
    npft_trif = 1
    cs => soil_carb
    rsa => rsp_s_acc
    dim_cs1 = 1
    dim_cs2 = 1
  END IF
  !
  !  set up CO2 field to be passed down
  !
  IF (l_co2_interactive) THEN
    co2_dim_len = row_length
    co2_dim_row = rows
    co2_dim_lev = model_levels
  ELSE
    co2_dim_len = 1
    co2_dim_row = 1
    co2_dim_lev = 1
  END IF

  IF (l_use_seasalt_autoconv .OR. l_use_seasalt_sulpc .OR.          &
      l_use_seasalt_indirect .OR. l_use_seasalt_direct) THEN
    salt_dim1=row_length
    salt_dim2=rows
    salt_dim3=model_levels
  ELSE
    salt_dim1=1
    salt_dim2=1
    salt_dim3=1
  END IF

  IF (l_use_seasalt_sulpc .OR. l_dms_em_inter) THEN
    aero_dim1=row_length
    aero_dim2=rows
  ELSE
    aero_dim1=1
    aero_dim2=1
  END IF

  IF (l_use_seasalt_sulpc) THEN
    aero_dim3=model_levels
  ELSE
    aero_dim3=1
  END IF

  ! Allocate additional microphysics variables to full size
  ! IF in use, otherwise allocate minimum amount of space

  IF (l_mcr_qcf2) THEN  ! Allocate second cloud ice
    ALLOCATE ( qcf2_star(tdims%i_start:tdims%i_end,              &
                         tdims%j_start:tdims%j_end,              &
                         tdims%k_start:tdims%k_end) )
  ELSE
    ALLOCATE ( qcf2_star(1,1,1) )
  END IF

  IF (l_mcr_qrain) THEN  ! Allocate rain
    ALLOCATE ( qrain_star(tdims%i_start:tdims%i_end,             &
                          tdims%j_start:tdims%j_end,             &
                          tdims%k_start:tdims%k_end) )
  ELSE
    ALLOCATE ( qrain_star(1,1,1) )
  END IF

  IF (l_mcr_qgraup) THEN  ! Allocate graupel
    ALLOCATE ( qgraup_star(tdims%i_start:tdims%i_end,            &
                           tdims%j_start:tdims%j_end,            &
                           tdims%k_start:tdims%k_end) )
  ELSE
    ALLOCATE ( qgraup_star(1,1,1) )
  END IF

  ! The _star fields are used to store the increments to theta, q, qcl,
  !  and qcf


  ! ------------------------------------------------------------


ELSE IF (flag == 'ozoninit') THEN

  ! Code for ozone initialisation in atm_step

    ! IF the ozone_tracer is initialised to 0 then it should be reset to
    ! something realistic so set it to climatology expanded to 3D

  IF (l_use_cariolle) THEN
    IF (ozone_tracer(1,1,1) == 0.0) THEN
      WRITE(umMessage,*)'O3 tracer must not be set to 0 reset to clim'
      CALL umPrint(umMessage,src='atm_step_phys_init')

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)        &
!$OMP SHARED(tdims,ozone_tracer,ozone3d)
      DO k=tdims%k_start, tdims%k_end
        DO j=tdims%j_start, tdims%j_end
          DO i=tdims%i_start, tdims%i_end
            ozone_tracer(i,j,k)=ozone3d(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

      ! Halos updated
      ! DEPENDS ON: swap_bounds
      CALL swap_bounds(ozone_tracer,                             &
                       row_length, rows,                         &
                       tdims%k_len,                              &
                       offx, offy, fld_type_p, swap_field_is_scalar)
    ELSE
      WRITE(umMessage,*) 'At least first row of Ozone tracer is not zero'
      CALL umPrint(umMessage,src='atm_step_phys_init')
    END IF ! End of check for zero ozone
  END IF

  !   After a successful call ozone3D will contain the mixing ratios of
  !   ozone on the model's grid. Uses ozone tracer from the last timestep
  !   in the radiation scheme. Radiation is done first so all variables used
  !   are from the last timestep.
  !

  IF (l_use_ozoneinrad) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)        &
!$OMP SHARED(tdims,ozone_tracer,temp_ozone,ozone3d,rows,row_length)
    DO k=tdims%k_start, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          temp_ozone(i,j,k)=ozone_tracer(i,j,k)

          !  Use the 3D ozone in radiation but check for zero or -ve ozone as some
          !  tropospheric analysed ozone may be -ve; do not use -ve values.

          IF (temp_ozone(i,j,k) > 0.0) THEN
            ozone3d(i,j,k) = temp_ozone(i,j,k)
          ELSE
            IF (j > 1 .AND. j < rows .AND.     &
                i > 1 .AND. i < row_length) THEN       
              IF (ozone_tracer(i,j-1,k) > 0.0) THEN   
                ozone3d(i,j,k) = ozone_tracer(i,j-1,k)
              ELSE
                IF (ozone_tracer(i,j+1,k) > 0.0) THEN
                  ozone3d(i,j,k) = ozone_tracer(i,j+1,k)
                ELSE
                  IF (ozone_tracer(i-1,j,k) > 0.0) THEN    
                    ozone3d(i,j,k) = ozone_tracer(i-1,j,k)
                  ELSE
                    IF (ozone_tracer(i+1,j,k) > 0.0) THEN
                      ozone3d(i,j,k) = ozone_tracer(i+1,j,k)
                    END IF   ! Check on i+1>0
                  END IF    ! Check on i-1>0
                END IF    ! Check on j+1>0
              END IF    ! Check on j-1>0
            END IF    ! Check on j (1:rows) and i (1:row_length)
          END IF    ! Check on temp_ozone > 0
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF


  ! ----------------------------------------------------------------------


ELSE IF (flag =='MRtoSpHm' ) THEN

  ! Allocate mix_star variables to hold mixing ratio increments
  ! from atmos_physics1

  ALLOCATE ( mix_v_star (tdims%i_start:tdims%i_end,            &
                         tdims%j_start:tdims%j_end,            &
                         tdims%k_start:tdims%k_end) )
  ALLOCATE ( mix_cl_star(tdims%i_start:tdims%i_end,            &
                         tdims%j_start:tdims%j_end,            &
                         tdims%k_start:tdims%k_end) )
  ALLOCATE ( mix_cf_star(tdims%i_start:tdims%i_end,            &
                         tdims%j_start:tdims%j_end,            &
                         tdims%k_start:tdims%k_end) )

  IF (l_mcr_qcf2) THEN
    ALLOCATE ( mix_cf2_star(tdims%i_start:tdims%i_end,         &
                            tdims%j_start:tdims%j_end,         &
                            tdims%k_start:tdims%k_end) )
  ELSE
    ALLOCATE ( mix_cf2_star(1,1,1) )
  END IF

  IF (l_mcr_qrain) THEN
    ALLOCATE ( mix_rain_star(tdims%i_start:tdims%i_end,        &
                             tdims%j_start:tdims%j_end,        &
                             tdims%k_start:tdims%k_end) )
  ELSE
    ALLOCATE ( mix_rain_star(1,1,1) )
  END IF

  IF (l_mcr_qgraup) THEN
    ALLOCATE ( mix_graup_star(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end) )
  ELSE
    ALLOCATE ( mix_graup_star(1,1,1) )
  END IF

  ! q_star currently holds mixing ratio increments from physics1
  ! but we wish it to hold specific humidity increments.
  ! Start by copying q_star into mix_v_star etc.

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)          &
!$OMP SHARED(tdims,mix_v_star,mix_cl_star,mix_cf_star,q_star,qcl_star, &
!$OMP qcf_star,mix_cf2_star,qcf2_star,l_mcr_qcf2,l_mcr_qrain,          &
!$OMP mix_rain_star,qrain_star,l_mcr_qgraup,mix_graup_star,qgraup_star)

!$OMP DO SCHEDULE(STATIC)
  DO k=tdims%k_start, tdims%k_end
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        mix_v_star(i,j,k)  = q_star(i,j,k)     ! Vapour
        mix_cl_star(i,j,k) = qcl_star(i,j,k)   ! Liquid
        mix_cf_star(i,j,k) = qcf_star(i,j,k)   ! Ice
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  IF (l_mcr_qcf2) THEN  ! Ice2
!$OMP DO SCHEDULE(STATIC)
    DO k=tdims%k_start, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          mix_cf2_star(i,j,k) = qcf2_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_mcr_qrain) THEN  ! Rain
!$OMP DO SCHEDULE(STATIC)
    DO k=tdims%k_start, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          mix_rain_star(i,j,k) = qrain_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_mcr_qgraup) THEN  ! Graupel
!$OMP DO SCHEDULE(STATIC)
    DO k=tdims%k_start, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          mix_graup_star = qgraup_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

!$OMP END PARALLEL

  ! ----------------------------------------------------------------------


ELSE IF (flag == 'tropinit') THEN

  ! Set logicals for tropopause diagnostics

  ! i_ozone_int  : cntlatm
  ! IO3_TROP_MAP : o3intp.h

  l_o3_trop_level  = ( sf_calc(280,2) ) .AND.                          &
                     ( ( i_ozone_int  ==  io3_trop_map) .OR.      &
                       ( i_ozone_int  ==  io3_trop_map_masscon) )
  l_o3_trop_height = ( sf_calc(281,2) ) .AND.                          &
                     ( ( i_ozone_int  ==  io3_trop_map) .OR.      &
                       ( i_ozone_int  ==  io3_trop_map_masscon) )
  l_t_trop_level   = ( sf_calc(282,2) ) .AND.                          &
                     ( ( i_ozone_int  ==  io3_trop_map) .OR.      &
                       ( i_ozone_int  ==  io3_trop_map_masscon) )
  l_t_trop_height  = ( sf_calc(283,2) ) .AND.                          &
                     ( ( i_ozone_int  ==  io3_trop_map) .OR.      &
                       ( i_ozone_int  ==  io3_trop_map_masscon) )

  ! Allocate space for tropopause diagnostics
  IF ( l_o3_trop_level ) THEN
    ALLOCATE ( o3_trop_level(row_length,rows) )
  ELSE
    ALLOCATE ( o3_trop_level(1,1) )
  END IF
  IF ( l_o3_trop_height ) THEN
    ALLOCATE ( o3_trop_height(row_length,rows) )
  ELSE
    ALLOCATE ( o3_trop_height(1,1) )
  END IF
  IF ( l_t_trop_level ) THEN
    ALLOCATE ( t_trop_level(row_length,rows) )
  ELSE
    ALLOCATE ( t_trop_level(1,1) )
  END IF
  IF ( l_t_trop_height ) THEN
    ALLOCATE ( t_trop_height(row_length,rows) )
  ELSE
    ALLOCATE ( t_trop_height(1,1) )
  END IF

  ! Check whether ozone option is applicable to the specified
  ! ozone ancillary file. If not output error message.
  IF (printstatus  >=   prstatus_diag) THEN
    WRITE(umMessage,*) 'Atm Step: Lexpand_ozone', lexpand_ozone
    CALL umPrint(umMessage,src='atm_step_phys_init')
  END IF

  IF ((lexpand_ozone) .AND. (i_ozone_int == 1)) THEN
    errorstatus=123

    CALL ereport("ATM_STEP", errorstatus,                        &
       "A 2D ozone ancillary has been specified with a" //       &
       "3D ozone option.")
  END IF

  IF (.NOT. lexpand_ozone) THEN
    IF (i_ozone_int == 2) THEN
      errorstatus=123

      CALL ereport("ATM_STEP", errorstatus,                     &
         "A 3D ozone ancillary has been specified with a" //    &
         "2D ozone option.")
    END IF
    IF (i_ozone_int == 5) THEN
      errorstatus=123

      CALL ereport("ATM_STEP", errorstatus,                     &
        "A 3D ozone ancillary has been specified with a" //    &
        "2D ozone option.")
    END IF
  END IF

  ! Convert the ozone field supplied from the dump to a full 3-D
  ! array.



  IF (lexpand_ozone) THEN
    nd_o3=rows*o3dims2%k_len
  ELSE
    nd_o3=row_length*rows*o3dims2%k_len
  END IF

  ! Allocate workspace, required only in atmos_physics1
  ALLOCATE (ozone3d (row_length, rows, o3dims2%k_start:o3dims2%k_end) )


  ! ----------------------------------------------------------------------


ELSE IF (flag == 'bc_exner') THEN

  ! Calculate EXNER_LBCs at next time level

  ! Change halo size for exner_lbc
  lbc_size=lenrima(fld_type_p,halo_type_extended,rima_type_norm)

  IF (rim_stepsa  ==  0) THEN

    ! No LBC updating
    DO k=1,model_levels+1
      DO i=1,lbc_size
        exner_lbc_real_tend(i,k)=exner_lbc(i,k)
      END DO !i
    END DO !k

  ELSE IF (MOD(timestep_number,rim_stepsa)  ==  0) THEN

    ! End of current LBC period
    DO k=1,model_levels+1
      DO i=1,lbc_size
        exner_lbc_real_tend(i,k)= exner_lbc_tend(i,k)
      END DO !i
    END DO !k

  ELSE ! Just a normal timestep during a LBC period

    increment_factor=1.0/ (rim_stepsa-MOD(timestep_number-1,rim_stepsa))

    DO k=1,model_levels+1
      DO i=1,lbc_size
        exner_lbc_real_tend(i,k) = exner_lbc(i,k) + increment_factor *  &
                                  (exner_lbc_tend(i,k) - exner_lbc(i,k))
      END DO !i
    END DO !k

  END IF ! End of current LBC period?

  !  n_rims_to_do set at begining of subroutine. Only wts=1 rims
  !  updated here since derived p fields are only being used for
  !  diagnostics. They get properly weighted at start of next timestep
  l_do_halos=.FALSE.
  l_do_boundaries=.TRUE.

  ! ----------------------------------------------------------------------

ELSE IF (flag == 'turb_cof') THEN

  ! Code to prepare for Smagorinsky turbulence scheme

  IF (L_subfilter_horiz .OR. L_subfilter_vert) THEN

    ALLOCATE (visc_h(1-halo_i:row_length+halo_i,                       &
                     1-halo_j:rows+halo_j, model_levels) )
    ALLOCATE (visc_m(1-halo_i:row_length+halo_i,                       &
                     1-halo_j:rows+halo_j, model_levels) )
    ALLOCATE (rneutml_sq(row_length, rows, model_levels-1))
    ALLOCATE (shear     (row_length, rows, model_levels-1))
    ALLOCATE (max_diff  (row_length, rows) )
    ALLOCATE (delta_smag(row_length, rows) )

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)  SHARED(halo_i,row_length, &
!$OMP halo_j,rows,model_levels,visc_m,visc_h,delta_smag)

!$OMP DO SCHEDULE(STATIC)
    DO k=1,model_levels
      DO j=1-halo_j,rows+halo_j
        DO i=1-halo_i,row_length+halo_i
          visc_m(i,j,k)  = 0.0
          visc_h(i,j,k)  = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO j=1,rows
      DO i=1,row_length
        delta_smag(i,j)= 1.0
      END DO
    END DO
!$OMP END DO

!$OMP END PARALLEL

  ELSE

    ALLOCATE ( visc_h(1,1,1) )
    ALLOCATE ( visc_m(1,1,1) )
    ALLOCATE ( rneutml_sq(1,1,1) )
    ALLOCATE ( shear(1,1,1) )
    ALLOCATE ( max_diff(1,1) )
    ! A full array is required in convection and boundary layer so lets
    ! just allocate to the required sized (even if not used).
    ALLOCATE (delta_smag(row_length, rows) )

  END IF !  L_subfilter_horiz .OR. L_subfilter_vert

END IF ! flag

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE atm_step_phys_init

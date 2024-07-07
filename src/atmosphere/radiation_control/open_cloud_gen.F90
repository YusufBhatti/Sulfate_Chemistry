! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
MODULE open_cloud_gen_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'OPEN_CLOUD_GEN_MOD'
CONTAINS

SUBROUTINE open_cloud_gen(                                              &
!                       Parallel variables
  global_row_length, global_rows                                        &
, me, n_proc, at_extremity                                              &
!                       Model Dimensions
, row_length, rows                                                      &
, n_profile, p_temp, offx, offy                                         &
!                       Properties of clouds
, w_cloud1, dp_corr_strat, cct, cloud_levels                            &
, cca, n_cca_levels, ccw, qcl_n, qcf_n, xx_cos_theta_latitude           &
!                       Model switches
, l_rad_step_diag, l_rad_step_prog                                      &
!                       in time stepping information.
, val_year, val_day_number, val_hour, val_minute                        &
, val_second                                                            &
!                       error information
, ierr)

!  Subroutine to generate sub-columns for McICA

! Method:
!       Allocates arrays to store sub-columns, rearranges cloud fields
!       so that they are in suitable order for generator and calls
!       generator which fills these arrays. Copies cloudy sub-columns
!       if more are required

USE rad_pcf
USE rand_no_mcica
USE mcica_mod
USE sw_control_struct
USE rad_input_mod, ONLY: rad_mcica_sampling, rad_mcica_sigma,           &
                         two_d_fsd_factor, i_fsd
USE cld_generator_mod
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE UM_ParVars, ONLY: g_datastart
USE cloud_inputs_mod, ONLY: l_add_cca_to_mcica
USE level_heights_mod, ONLY: r_theta_levels
USE fsd_parameters_mod, ONLY: f_arr, f_cons, ip_fsd_constant,           &
                              fsd_eff_lam, fsd_eff_phi
USE nlsizes_namelist_mod, ONLY: model_levels

USE model_domain_mod, ONLY: model_type, mt_global, mt_single_column

IMPLICIT NONE



!                       Dummy arguments.
INTEGER, INTENT(INOUT) ::                                             &
  ierr
!           Error flag
!
!                       Parallel Variables
INTEGER, INTENT(IN) ::                                                &
  global_row_length                                                   &
!           Number of global points on a row
  , global_rows                                                         &
!           Number of global rows
  , n_proc                                                              &
!           Total number of processors
  , me
!           My processor number
!
LOGICAL, INTENT(IN) ::                                                &
  at_extremity(4)
!           Indicates if this processor is at north, south
!           east or west of the processor grid
!
!
!                       Dimensions of arrays
INTEGER, INTENT(IN) ::                                                &
  row_length                                                          &
!           Number of points on a row
  , rows                                                                &
!           Number of rows
  , n_profile                                                           &
!           Size allocated for atmospheric profiles
  , n_cca_levels                                                        &
!           number of convective cloud amount levels
  , offx                                                                &
!           Size of small halo in i
  , offy
!           Size of small halo in i

!
!                       Information about model time
INTEGER, INTENT(IN) ::                                                &
  val_year                                                            &
!           time information for current timestep
  , val_day_number                                                      &
!           time information for current timestep
  , val_hour                                                            &
!           time information for current timestep
  , val_minute                                                          &
!           time information for current timestep
  , val_second
!           time information for current timestep

REAL, INTENT(IN) ::                                                   &
  p_temp(n_profile, 0:model_levels)                                   &
!           Pressure
  , w_cloud1(n_profile, model_levels)                                   &
!           Amount of cloud
  , cca(n_profile, n_cca_levels)                                        &
!           Convective cloud amount
  , ccw(n_profile, model_levels)                                        &
!           Convective in-cloud water content (kg/kg)
  , qcl_n(n_profile, model_levels), qcf_n(n_profile, model_levels)      &

  , xx_cos_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)
!           Finite volume cosine of latitude.
!
!                       Properties of clouds
REAL, INTENT(IN) ::                                                   &
  dp_corr_strat
!           Decorrelation pressure scale for large scale cloud!

INTEGER, INTENT(IN) ::                                                &
  cct(n_profile)
!           Level of top of convective cloud

INTEGER, INTENT(OUT) ::                                               &
  cloud_levels
!           Number of global cloud levels

!                       Model Switches
LOGICAL ::                                                            &
  l_rad_step_diag                                                     &
!           true if fast radiation timestep    (3C)
  , l_rad_step_prog
!           true if slow radiation timestep    (3C)
!
!
!
!                       Local variables.
INTEGER ::                                                            &
  i                                                                   &
!           Loop variable
  , j                                                                   &
!           Loop variable
  , k                                                                   &
!           Loop variable
  , l                                                                   &
!           Loop variable
  , ll                                                                  &
!           Loop variable
  , info                                                                &
!           Loop variable
  , list(n_profile)                                                     &
!           Global location of each point on the processor.
  , random_dummy_init                                                   &
!           Seed for generating seeds for random numbers used in the
!           generator
  , random_dummy(n_profile)                                             &
!           Seed for generating seeds for random numbers used in the
!           generator
  , first_row                                                           &
!           index of first row
  , last_row                                                            &
!           index of last row
  , cloud_top
!           Top global cloudy layer

REAL ::                                                               &
  p(n_profile, model_levels)                                          &
!           Pressure
  , eps                                                                 &
!           small number to prevent division by zero.
  , cloud_scale                                                         &
!           cloud fraction times gridbox size
  , thickness_part(n_profile, model_levels)                             &
!           part of FSD param related to layer thickness
  , dp_corr_cloud(n_profile,model_levels)                               &
!           Cloud fraction decorrelation length
  , dp_corr_cond(n_profile,model_levels)                                &
!           Cloud condensate decorrelation length
  , sigma_qcw(n_profile,model_levels)                                   &
!           Normalized cloud condensate std. dev
  , w_cloud(n_profile, 0:model_levels)                                  &
!           Amount of cloud
  , c_cloud(n_profile, model_levels)                                    &
!           Amount of convective cloud
  , c_ratio(n_profile, model_levels)                                    &
!           Ratio of convective cloud condensate to mean condensate
  , ls_ratio(n_profile, model_levels)                                   &
!           Ratio of large-scale cloud condensate to mean condensate
  , mix_ratio(n_profile, model_levels)                                  &
!           Mean in-cloud condensate mixing ratio
  , x_in_km
!           grid-box size in km

REAL,PARAMETER ::                                                     &
  one_third = 1.0/3.0
!           Constant for use in FSD parametrization

LOGICAL ::                                                            &
  l_layer_clear
!           Flag for layer free of clouds
INTEGER ::                                                            &
  k_dum
!           dummy storage for k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OPEN_CLOUD_GEN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ((l_rad_step_prog) .OR. (l_rad_step_diag)) THEN

  IF ((sw_control(1)%i_cloud == ip_cloud_mcica) .OR.                  &
      (i_fsd /= ip_fsd_constant)) THEN

    eps=EPSILON(eps)

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(j, i) SCHEDULE(STATIC)       &
!$OMP& SHARED(model_levels,n_profile,l_add_cca_to_mcica,n_cca_levels, &
!$OMP&    cca,ccw,w_cloud,w_cloud1,c_cloud,c_ratio,ls_ratio,mix_ratio,&
!$OMP&    qcl_n,qcf_n,eps)
    DO j=1,model_levels
      DO i=1,n_profile
        IF (l_add_cca_to_mcica .AND. j <= n_cca_levels .AND.          &
            cca(i,j) > 0.0 .AND. ccw(i,j) > 0.0) THEN
          w_cloud(i,model_levels+1-j)=cca(i,j)+(1.0-cca(i,j))*w_cloud1(i,j)
        ELSE
          c_cloud(i,model_levels+1-j) = 0.0
          c_ratio(i,model_levels+1-j) = 0.0
          ls_ratio(i,model_levels+1-j) = 1.0
          w_cloud(i,model_levels+1-j)=w_cloud1(i,j)
        END IF
        IF (w_cloud(i,model_levels+1-j) > cut) THEN
          IF (l_add_cca_to_mcica .AND. j <= n_cca_levels .AND.        &
              cca(i,j) > 0.0 .AND. ccw(i,j) > 0.0) THEN
            c_cloud(i,model_levels+1-j)=cca(i,j)
            mix_ratio(i,model_levels+1-j) = ( qcl_n(i,j)+qcf_n(i,j) &
              +cca(i,j)*ccw(i,j) ) / w_cloud(i,model_levels+1-j)
          c_ratio(i,model_levels+1-j) = ccw(i,j) / mix_ratio(i,model_levels+1-j)
            ls_ratio(i,model_levels+1-j) = (qcl_n(i,j)+qcf_n(i,j)) &
              / (w_cloud1(i,j)*mix_ratio(i,model_levels+1-j) + eps)
          END IF
        ELSE
          w_cloud(i,model_levels+1-j) = 0.0
        END IF
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  IF (i_fsd == ip_fsd_constant) THEN
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                      &
!$OMP& SHARED( model_levels, n_profile, sigma_qcw, rad_mcica_sigma )   &
!$OMP& PRIVATE( i, j )
    DO j=1, model_levels
      DO i=1, n_profile
        sigma_qcw(i,j)=rad_mcica_sigma
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE
    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          l=(j-1)*row_length+i
          ll=model_levels+1-k
          x_in_km = 0.001*SQRT ( r_theta_levels(i,j,k) * fsd_eff_lam    &
                               * r_theta_levels(i,j,k) * fsd_eff_phi    &
                               * xx_cos_theta_latitude(i,j)  )
          IF (w_cloud(l,ll) == 1.0) THEN
            sigma_qcw(l,ll) = two_d_fsd_factor * f_arr(1,i,j,k)       &
            * (x_in_km ** one_third)                                  &
            * ((((f_cons(1) * x_in_km) ** f_cons(2)) + 1)             &
            ** (f_cons(3)))
          ELSE
            cloud_scale=w_cloud(l,ll)*x_in_km
            sigma_qcw(l,ll) = two_d_fsd_factor *                      &
              (f_arr(2,i,j,k)-(f_arr(3,i,j,k)*w_cloud(l,ll)))         &
              * (cloud_scale ** one_third)                            &
              * ((((f_cons(1) * cloud_scale) ** f_cons(2)) + 1)       &
              ** (f_cons(3)))
          END IF
        END DO
      END DO
    END DO

    ALLOCATE(cloud_inhom_param_full(n_profile,1:model_levels))

    DO j=1, model_levels
      DO i=1, n_profile
        ! Calculate scaling parameter, =EXP(MEAN(LOG(water content)))/MEAN(WC)
        ! from standard deviation of water content, assuming a log-normal
        ! distribution to simplify the maths.
        cloud_inhom_param_full(i,j)=(1.0/SQRT((sigma_qcw(i,j)**2)+1))
      END DO
    END DO
  END IF


  IF (sw_control(1)%i_cloud == ip_cloud_mcica) THEN

    cloud_levels = 1
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                      &
!$OMP& SHARED( model_levels, n_profile, w_cloud )                      &
!$OMP& PRIVATE( i, j ) REDUCTION(MAX:cloud_levels)
    DO j=1,model_levels
      DO i=1,n_profile
        IF (w_cloud(i,model_levels+1-j)  >   cut) THEN
          cloud_levels=MAX(cloud_levels, j)
        END IF
      END DO
    END DO
!$OMP END PARALLEL DO

    IF (model_type /= mt_single_column) THEN
      ! To obtain reproducible results independent of the
      ! decomposition of the domain used on an MPP machine a global
      ! value for the topmost cloudy layer is used.
      CALL gc_imax(1, n_proc, info, cloud_levels)
    END IF

    cloud_top=model_levels+1-cloud_levels

    IF (l_rad_step_prog) THEN

      ! Set the SW and LW values of subcol_k (the number of sub-columns
      ! each k-term "sees") and subcol_reorder (a reordering of the
      ! sub-columns so that each sub-column is equivalently important in
      ! the SW and LW).
      !
      SELECT CASE (rad_mcica_sampling)

      CASE (ip_mcica_full_sampling)
        subcol_need=tot_subcol_gen
        sw_subcol_k=tot_subcol_gen
        lw_subcol_k=tot_subcol_gen

        ALLOCATE(lw_subcol_reorder(subcol_need))
        DO i=1,subcol_need
          lw_subcol_reorder(i)=i
        END DO

      CASE (ip_mcica_single_sampling)
        subcol_need=subcol_need_single
        sw_subcol_k=1
        lw_subcol_k=1

        ALLOCATE(lw_subcol_reorder(subcol_need))
        DO i=1,subcol_need
          lw_subcol_reorder(i) =                                      &
            MOD(lw_subcol_reorder_single(i),tot_subcol_gen)
          IF (lw_subcol_reorder(i) == 0) THEN
            lw_subcol_reorder(i) = tot_subcol_gen
          END IF
        END DO

      CASE (ip_mcica_optimal_sampling)
        subcol_need=subcol_need_optimal
        !          sw_subcol_k and lw_subcol_k have been read from data file

        ALLOCATE(lw_subcol_reorder(subcol_need))
        DO i=1,subcol_need
          lw_subcol_reorder(i) =                                      &
            MOD(lw_subcol_reorder_optimal(i),tot_subcol_gen)
          IF (lw_subcol_reorder(i) == 0) THEN
            lw_subcol_reorder(i) = tot_subcol_gen
          END IF
        END DO

      END SELECT

      ALLOCATE(sw_subcol_reorder(subcol_need))
      DO i=1,subcol_need
        sw_subcol_reorder(i) = MOD(i, tot_subcol_gen)
        IF (sw_subcol_reorder(i) == 0) THEN
          sw_subcol_reorder(i) = tot_subcol_gen
        END IF
      END DO

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(j, i)                         &
!$OMP& SHARED(model_levels,n_profile,p,p_temp,dp_corr_cloud,        &
!$OMP&    dp_corr_strat,dp_corr_cond)

!$OMP DO SCHEDULE(STATIC)  
      DO j=1,model_levels
        DO i=1,n_profile
          p(i,model_levels+1-j)=p_temp(i,j)
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)  
      DO j=1, model_levels
        DO i=1, n_profile
          dp_corr_cloud(i,j)=dp_corr_strat
          dp_corr_cond(i,j)=dp_corr_cloud(i,j)*0.5
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      ALLOCATE(clw_sub_full(n_profile,cloud_top:model_levels,tot_subcol_gen))
      !  ALLOCATE(cic_sub_full(n_profile,cloud_top:model_levels,tot_subcol_gen))
      ALLOCATE(frac_cloudy_full(n_profile))
      ALLOCATE(ncldy(n_profile))

      k = 0
      first_row=0
      last_row=rows-1
 
      ! included additional lines with k_dum usage have replaced k=k+1
      ! to remove dependency in the loop and to enable threading
      k_dum = k
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                      &
!$OMP& SHARED( first_row, last_row, row_length, k_dum, list, me,       &
!$OMP&         global_row_length )                                     &
!$OMP& PRIVATE( i, j, k )
      DO j=first_row, last_row
        DO i=0, row_length-1
          k = k_dum + row_length*(j-first_row) + i + 1
          list(k) = (g_datastart(2,me)+j-1)*global_row_length         &
                     + g_datastart(1,me)+i
        END DO
      END DO
!$OMP END PARALLEL DO
      k = k_dum + row_length*(last_row-first_row+1)


      random_dummy_init=(ABS(val_year-2000)*366*24*60*60)             &
        +(val_day_number*24*60*60)+(val_hour*60*60)                   &
        +(val_minute*60)+val_second

      DO i=1,n_profile
        random_dummy(i)=random_dummy_init+list(i)
      ! Zero out fields
        ncldy(i)  = 0
      END DO

      ! The initial values of random dummy are successive integers, which
      ! result in random numbers that are close to each other.This first
      ! call to mcica_rand_no is to ensure that the seed for each
      ! profile is itself random.

      CALL mcica_rand_no(random_dummy,n_profile,tot_subcol_gen)

      ALLOCATE(rand_seed_x(n_profile,tot_subcol_gen))

      DO i=1,tot_subcol_gen
        CALL mcica_rand_no(random_dummy,n_profile,1)
        DO j=1, n_profile
          rand_seed_x(j,i)=random_dummy(j)
        END DO
      END DO

      ALLOCATE(rand_seed_y(n_profile,tot_subcol_gen))

      DO i=1,tot_subcol_gen
        CALL mcica_rand_no(random_dummy,n_profile,1)
        DO j=1, n_profile
          rand_seed_y(j,i)=random_dummy(j)
        END DO
      END DO

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(k, j, i) SCHEDULE(STATIC)   &
!$OMP& SHARED(tot_subcol_gen,cloud_top,model_levels,n_profile,       &
!$OMP&    clw_sub_full)
      DO k = 1, tot_subcol_gen
        DO j = cloud_top, model_levels
          DO i = 1, n_profile
            !              cic_sub_full(i,j,k) = 0.0
            clw_sub_full(i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

      ! Set the overlap used in the cloud generator
      ioverlap=sw_control(1)%i_overlap

      CALL cld_generator(model_levels, cloud_top, n_profile          &
      , tot_subcol_gen, dp_corr_cloud, dp_corr_cond, sigma_qcw       &
      , w_cloud, c_cloud, c_ratio, ls_ratio, p, 1, n_profile)

      IF (ALLOCATED(rand_seed_y)) DEALLOCATE(rand_seed_y)
      IF (ALLOCATED(rand_seed_x)) DEALLOCATE(rand_seed_x)

      IF (rad_mcica_sampling /= ip_mcica_full_sampling) THEN
        DO i=1,n_profile
          frac_cloudy_full(i)=REAL(ncldy(i))/REAL(tot_subcol_gen)
        END DO
      ELSE IF (rad_mcica_sampling == ip_mcica_full_sampling) THEN
        ! In this case we treat the clear sub-columns as cloudy sub-columns
        ! so the clear-sky fraction is implicit in the summing of the
        ! "cloudy" sub-columns
        DO i=1,n_profile
          frac_cloudy_full(i)=1.0
        END DO
      END IF

      IF (rad_mcica_sampling /= ip_mcica_full_sampling) THEN
        ! For the case where are less cloudy subcolumns than required,
        ! copy cloudy values to ensure enough cloudy subcolumns
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE (i, j, k, l)               &
!$OMP& SHARED(n_profile,ncldy,subcol_need,tot_subcol_gen,cloud_top, &
!$OMP&    model_levels,clw_sub_full)
        DO i=1, n_profile
          IF (ncldy(i) < subcol_need .AND. ncldy(i) > 0) THEN
            DO j=ncldy(i)+1,MIN(subcol_need,tot_subcol_gen)
              DO k=cloud_top,model_levels
                l=j-ncldy(i)
                clw_sub_full(i,k,j)=clw_sub_full(i,k,l)
                !                  cic_sub_full(i,k,j)=cic_sub_full(i,k,l)
              END DO
            END DO
          END IF
        END DO
!$OMP END PARALLEL DO
      END IF
    END IF! l_rad_step_prog

  ELSE! If not using McICA in either SW or LW

    l_layer_clear=.TRUE.
    DO j=model_levels,1,-1
      DO i=1,n_profile
        l_layer_clear = l_layer_clear .AND.                             &
          (w_cloud1(i, j) <= 0.0e+00) .AND. (cct(i) < j-1)
      END DO
      IF (.NOT. l_layer_clear) THEN
        cloud_levels = j
        EXIT
      END IF
    END DO

    IF (model_type /= mt_single_column) THEN
      !     To obtain reproducible results independent of the
      !     decomposition of the domain used on an MPP machine a global
      !     value for the topmost cloudy layer is used.
      CALL gc_imax(1, n_proc, info, cloud_levels)
    END IF

  END IF! not mcica

END IF! l_rad_step_prog .OR. l_rad_step_diag

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE open_cloud_gen
END MODULE open_cloud_gen_mod

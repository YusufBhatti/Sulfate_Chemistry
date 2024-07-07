! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! PC2 homogenous forcing for convective increments and optional erosion

MODULE pc2_from_conv_ctl_mod

IMPLICIT NONE
!---------------------------------------------------------------------------
! Description: Called from ni_conv_ctl to do either PC2 homogenous forcing
! using convection increments or  PC2 homogenous forcing plus erosion

! method:

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large scale cloud

! code description:
!   language: fortran 90
!   This code is written to umdp3 programming standards version 10.4
!---------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PC2_FROM_CONV_CTL_MOD'

CONTAINS

SUBROUTINE  pc2_from_conv_ctl(                                            &
                 row_length, rows, n_conv_levels,                         &
                 ! Model switches
                 l_calc_dxek, l_q_interact, l_mixing_ratio,               &
                 ! Coordinate info
                 exner_theta_levels, p_layer_centres,                     &  
                 ! in/out Convection increments
                 theta_inc, q_inc, qcl_inc, qcf_inc,                      &
                 cf_liquid_inc, cf_frozen_inc, bulk_cf_inc,               &
                 ! in/out values with all increments added
                 theta_star, q_star, qcl_star, qcf_star,                  &
                 cf_liquid_star, cf_frozen_star, bulk_cf_star             &
                     )

USE atm_fields_bounds_mod, ONLY:  tdims, tdims_s

! Model time stepping information
USE timestep_mod, ONLY: timestep

USE cloud_inputs_mod, ONLY: i_pc2_erosion_method, l_micro_eros,   &
                            dbsdtbs_turb_0
USE pc2_constants_mod,ONLY: pc2eros_hybrid_allfaces, dbsdtbs_turb_1

USE cv_run_mod,  ONLY: l_pc2_diag_sh

! Convective diagnostic output arrays
USE cv_diagnostic_array_mod, ONLY:                                      &
        theta_incr_diag_conv                                            &
       ,T_incr_diag_conv ,q_incr_diag_conv ,qcl_incr_diag_conv          &
       ,qcf_incr_diag_conv ,cf_liquid_incr_diag_conv                    &
       ,cf_frozen_incr_diag_conv ,bulk_cf_incr_diag_conv                &
       ,qcl_incr_inhom_diag ,qcf_incr_inhom_diag                        &
       ,bulk_cf_incr_inhom_diag ,cf_liquid_incr_inhom_diag              &
       ,cf_frozen_incr_inhom_diag
USE cv_stash_flg_mod, ONLY:                                             &
        l_theta_incr_conv                                               &
       ,l_T_incr_conv, l_q_incr_conv, l_qcl_incr_conv, l_qcf_incr_conv  &
       ,l_cfl_incr_conv, l_bcf_incr_conv

! Add convection tendencies to physics_tendencies module SPT scheme
USE physics_tendencies_mod,  ONLY:                                      &
    l_retain_conv_tendencies, dt_conv, dq_conv


USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

INTEGER, INTENT(IN)  ::      &
  row_length                 & ! Row length
, rows                       & ! Number of rows
, n_conv_levels                ! number of convection levels

! Model switches
LOGICAL, INTENT(IN)  ::  &
  l_calc_dxek            & ! Switch for calculation of condensate increment
, l_q_interact           & ! Switch allows overwriting of parcel variables
                           ! when calculating condensate increments.
, l_mixing_ratio           ! Use mixing ratio formulation

REAL, INTENT(IN) ::                                  &
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,  & ! Exner on theta level
                     tdims_s%j_start:tdims_s%j_end,  & !
                     tdims_s%k_start:tdims_s%k_end)  &
, p_layer_centres(tdims%i_end,                       & ! pressure at layer
                  tdims%j_end,                       & ! centres (Pa)
                  0:tdims%k_end)

! '_inc' - convection + PC2 increments to a field
REAL,INTENT(INOUT) ::                                                    &
  theta_inc(tdims%i_end,tdims%j_end,tdims%k_end)                         &
, q_inc(tdims%i_end,tdims%j_end,tdims%k_end)                             &
, qcl_inc(tdims%i_end,tdims%j_end,tdims%k_end)                           &
, qcf_inc(tdims%i_end,tdims%j_end,tdims%k_end)                           &
, cf_liquid_inc(tdims%i_end,tdims%j_end,tdims%k_end)                     &
, cf_frozen_inc(tdims%i_end,tdims%j_end,tdims%k_end)                     &
, bulk_cf_inc(tdims%i_end,tdims%j_end,tdims%k_end)

! arguments with intent in/out. ie: input variables changed on output.
! '_star' = '_n' + all increments up to now in time step
! _star variables held on same levels as prognostics in the vertical for
! ENDGame. 

REAL,INTENT(INOUT) ::                                                    &
  theta_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
                         1:tdims%k_end)                                  &
, q_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
                     1:tdims%k_end)                                      &
, qcl_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
                       1:tdims%k_end)                                    &
, qcf_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
                       1:tdims%k_end)                                    &
, cf_liquid_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                             1:tdims%k_end)                              &
, cf_frozen_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                             1:tdims%k_end)                              &
, bulk_cf_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
               1:tdims%k_end)


!----------------------------------------------------------------------
! Local variables
!----------------------------------------------------------------------
INTEGER ::       &
  i, j, k        & ! Loop counters
 ,kp1, km1         ! Level numbers plus one, minus one

LOGICAL ::       &
  l_full_zero      ! True if a dummy zero full field is required

LOGICAL ::            &
 l_pc2_prod_qcl_mp

! Dummy array to pass to pc2_turbulence_ctl.
REAL ::                                                                 &
 dummy_zeros(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             1:tdims%k_end)

! Allocatable arrays for PC2 scheme increment calculation - used to save
! memory when increment calculations are not requested.
REAL,ALLOCATABLE:: &
  t_earliest (:,:,:)          & !  temperature at start of convection
, t_inc_latest (:,:,:)        & !  temperature increment on input to PC2 homog
, q_earliest (:,:,:)          & !  humidity    at start of convection
, qcl_earliest (:,:,:)        & !  qCL         at start of convection
, cfl_earliest (:,:,:)        & !  cf_liquid   at start of convection
, cff_earliest (:,:,:)        & !  cf_frozen   at start of convection
, bcf_earliest (:,:,:)        & !  bulk cloud  at start of convection
, theta_inc_pc2(:,:,:)        & !  pot temperature increment due to PC2 homog
, q_inc_pc2 (:,:,:)           & !  humidity        increment due to PC2 homog
, qcl_inc_pc2(:,:,:)          & !  qCL             increment due to PC2 homog
, cfl_inc_pc2(:,:,:)          & !  cf_liquid       increment due to PC2 homog
, bcf_inc_pc2(:,:,:)          & !  bulk cloud      increment due to PC2 homog
, bcf_above(:,:,:)            & !  Bulk cloud fraction in layer above
, bcf_below(:,:,:)              !  Bulk cloud fraction in layer below

REAL,ALLOCATABLE:: full_zero(:,:,:)        !  a dummy array for a zero field


CHARACTER(LEN=*), PARAMETER :: RoutineName='PC2_FROM_CONV_CTL'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------

l_full_zero = l_calc_dxek

! L_full_zero_if1:
IF (l_full_zero) THEN

  !  Set up a field of zero values if dummy field is required

  ALLOCATE ( full_zero(row_length,rows,1) )

  DO j=1, rows
    DO i=1, row_length
      full_zero(i,j,1)  = 0.0
    END DO  ! I
  END DO  ! J

END IF  ! L_full_zero_if1

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(k, j, i, t_inc_latest,            &
!$OMP& theta_inc_pc2, q_inc_pc2, qcl_inc_pc2, cfl_inc_pc2, bcf_inc_pc2, &
!$OMP& t_earliest, q_earliest, qcl_earliest, bcf_earliest,              &
!$OMP& cfl_earliest, cff_earliest, bcf_above, bcf_below, km1, kp1,      &
!$OMP& l_pc2_prod_qcl_mp)                                               &
!$OMP& SHARED(n_conv_levels,rows,row_length,i_pc2_erosion_method,       &
!$OMP&     theta_inc,exner_theta_levels,theta_star,q_star,q_inc,        &
!$OMP&     qcl_star,qcl_inc,bulk_cf_star,bulk_cf_inc,cf_liquid_star,    &
!$OMP&     cf_liquid_inc,cf_frozen_star,cf_frozen_inc,l_pc2_diag_sh,    &
!$OMP&     l_micro_eros,p_layer_centres,timestep,full_zero,             &
!$OMP&     cf_liquid_incr_inhom_diag,l_mixing_ratio,dbsdtbs_turb_0,     &
!$OMP&     l_q_interact,                                                &
!$OMP&     theta_incr_diag_conv,t_incr_diag_conv,q_incr_diag_conv,      &
!$OMP&     qcl_incr_diag_conv,cf_liquid_incr_diag_conv,                 &
!$OMP&     bulk_cf_incr_diag_conv,l_retain_conv_tendencies,tdims,       &
!$OMP&     dt_conv,dq_conv,dummy_zeros,                                 &
!$OMP&     l_theta_incr_conv, l_T_incr_conv,                            &
!$OMP&     l_q_incr_conv, l_qcl_incr_conv,                              &
!$OMP&     l_cfl_incr_conv, l_bcf_incr_conv )

ALLOCATE ( t_inc_latest(row_length,rows,1) )

ALLOCATE ( theta_inc_pc2(row_length,rows,1) )
ALLOCATE ( q_inc_pc2(row_length,rows,1) )
ALLOCATE ( qcl_inc_pc2(row_length,rows,1) )
ALLOCATE ( cfl_inc_pc2(row_length,rows,1) )
ALLOCATE ( bcf_inc_pc2(row_length,rows,1) )

ALLOCATE ( t_earliest(row_length,rows,1) )
ALLOCATE ( q_earliest(row_length,rows,1) )
ALLOCATE ( qcl_earliest(row_length,rows,1) )
ALLOCATE ( cfl_earliest(row_length,rows,1) )
ALLOCATE ( cff_earliest(row_length,rows,1) )
ALLOCATE ( bcf_earliest(row_length,rows,1) )

IF (i_pc2_erosion_method == pc2eros_hybrid_allfaces) THEN
  ALLOCATE ( bcf_above(row_length,rows,1) )
  ALLOCATE ( bcf_below(row_length,rows,1) )
ELSE
  ALLOCATE ( bcf_above(1,1,1) )
  ALLOCATE ( bcf_below(1,1,1) )
END IF

! Because of the switches used for the following pc2 calls, a
! dummy array input is needed in place of a variable that is not
! used.
!$OMP DO SCHEDULE(STATIC)
DO k=1, tdims%k_end
  DO j=tdims%j_start, tdims%j_end
    DO i=tdims%i_start, tdims%i_end
      dummy_zeros(i,j,k)=0.0
    END DO
  END DO
END DO
!$OMP END DO

! Convert potential temperature increment to temperature increment

!$OMP  DO SCHEDULE(DYNAMIC)
   
DO k=1, n_conv_levels

  DO j=1, rows
    DO i=1, row_length
      t_inc_latest(i,j,1) = theta_inc(i,j,k) * exner_theta_levels(i,j,k)

      t_earliest(i,j,1)   = theta_star(i,j,k)*exner_theta_levels(i,j,k)  &
                                                    - t_inc_latest(i,j,1)
      q_earliest(i,j,1)   = q_star(i,j,k)       - q_inc(i,j,k)
      qcl_earliest(i,j,1) = qcl_star(i,j,k)     - qcl_inc(i,j,k)
      bcf_earliest(i,j,1) = bulk_cf_star(i,j,k) - bulk_cf_inc(i,j,k)
      cfl_earliest(i,j,1) = cf_liquid_star(i,j,k) - cf_liquid_inc(i,j,k)
      cff_earliest(i,j,1) = cf_frozen_star(i,j,k) - cf_frozen_inc(i,j,k)
    END DO ! i
  END DO ! j

  IF (i_pc2_erosion_method == pc2eros_hybrid_allfaces) THEN

    IF (k==1) THEN
      km1=k
    ELSE
      km1=k-1
    END IF

    IF (k==n_conv_levels) THEN
      kp1=k
    ELSE
      kp1=k+1
    END IF

    DO j=1, rows
      DO i=1, row_length
        bcf_above(i,j,1)    = bulk_cf_star(i,j,kp1) - bulk_cf_inc(i,j,kp1)
        bcf_below(i,j,1)    = bulk_cf_star(i,j,km1) - bulk_cf_inc(i,j,km1)
      END DO ! i
    END DO ! j

  END IF

  IF (l_pc2_diag_sh .OR. l_micro_eros) THEN

    ! The following logical is set to indicate that pc2_hom_conv is not 
    ! being called by microphys_ctl. 
    l_pc2_prod_qcl_mp = .FALSE.

    ! Call pc2_hom_conv without erosion term i.e. dbsdtbs_turb_0=0.0
    ! and dbsdtbs_turb_1=0.0.

    ! DEPENDS ON: pc2_hom_conv
    CALL pc2_hom_conv(                                                    &
        ! Input variables
               p_layer_centres(1,1,k), 1                                  &
         ,     timestep                                                   &
        ! INput variables
         ,     t_earliest(1,1,1), q_earliest(1,1,1), qcl_earliest(1,1,1)  &
         ,     bcf_earliest(1,1,1),cfl_earliest(1,1,1),cff_earliest(1,1,1)&
         ,     t_inc_latest,q_inc(1,1,k),full_zero(1,1,1),full_zero(1,1,1)&
         ,     cf_liquid_inc(1,1,k)                                       &
         ,     bcf_above(1,1,1),bcf_below(1,1,1),dummy_zeros(1,1,1)       &
        ! OUTput variables
         ,     theta_inc_PC2,q_inc_PC2,qcl_inc_PC2,bcf_inc_PC2,cfl_inc_PC2&
        ! INput variables (other quantities)
         ,     0.0, 0.0                                                   &
        ! Model switches
         ,     l_mixing_ratio, l_pc2_prod_qcl_mp)
  ELSE
    ! call pc2_hom_conv with erosion term

    ! set to indicate that pc2_hom_conv is not being called by microphys_ctl
    !
    l_pc2_prod_qcl_mp = .FALSE.

    ! DEPENDS ON: pc2_hom_conv
    CALL pc2_hom_conv(                                                   &
        ! Input variables
              p_layer_centres(1,1,k), 1                                  &
        ,     timestep                                                   &
        ! INput variables
        ,     t_earliest(1,1,1), q_earliest(1,1,1), qcl_earliest(1,1,1)  &
        ,     bcf_earliest(1,1,1),cfl_earliest(1,1,1),cff_earliest(1,1,1)&
        ,     t_inc_latest,q_inc(1,1,k),full_zero(1,1,1),full_zero(1,1,1)&
        ,     cf_liquid_inc(1,1,k)                                       &
        ,     bcf_above(1,1,1),bcf_below(1,1,1),dummy_zeros(1,1,1)       &
        ! OUTput variables
        ,     theta_inc_PC2,q_inc_PC2,qcl_inc_PC2,bcf_inc_PC2,cfl_inc_PC2&
        ! INput variables (other quantities)
        ,     dbsdtbs_turb_0, dbsdtbs_turb_1                             &
        ! Model switches
        ,     l_mixing_ratio, l_pc2_prod_qcl_mp)
  END IF
  !
  ! Calculate potential temperature increment (on convect levels only)
  ! from temperature increment output by PC2_Homog.

  DO j=1,rows
    DO i=1,row_length
      theta_inc_pc2(i,j,1) = theta_inc_pc2(i,j,1)/exner_theta_levels(i,j,k)
    END DO ! i
  END DO ! j

  ! L_q_interact_if1:
  IF (l_q_interact) THEN  

    DO j = 1, rows
      DO i = 1, row_length

        ! Update increments to theta, moisture and cloud fields with additional
        ! increments from the response to environment changes (homogenous),

        theta_inc(i,j,k) = theta_inc(i,j,k) + theta_inc_pc2(i,j,1)
        q_inc(i,j,k)   =   q_inc(i,j,k) +   q_inc_pc2(i,j,1)
        qcl_inc(i,j,k) = qcl_inc(i,j,k) + qcl_inc_pc2(i,j,1)
        ! Not updated   qcf_inc(i,j,k) = qcf_inc(i,j,k)
        cf_liquid_inc(i,j,k) = cf_liquid_inc(i,j,k) + cfl_inc_pc2(i,j,1)
        ! Not updated   cf_frozen_inc(i,j,k) = cf_frozen_inc(i,j,k)
        bulk_cf_inc(i,j,k)   = bulk_cf_inc(i,j,k) +  bcf_inc_pc2(i,j,1)

        ! ... and update working version of theta, moisture and cloud fields.

        theta_star(i,j,k) = theta_star(i,j,k) + theta_inc_pc2(i,j,1)
        q_star(i,j,k)     =     q_star(i,j,k) +     q_inc_pc2(i,j,1)
        qcl_star(i,j,k)   =   qcl_star(i,j,k) +   qcl_inc_pc2(i,j,1)
        ! Not updated   qcf_star(i,j,k)   =   qcf_star(i,j,k)
        cf_liquid_star(i,j,k) = cf_liquid_star(i,j,k) + cfl_inc_pc2(i,j,1)
        ! Not updated   cf_frozen_n(i,j,k) = cf_frozen_n(i,j,k)
        bulk_cf_star(i,j,k)   = bulk_cf_star(i,j,k) +           &
                                  bcf_inc_pc2(i,j,1)

      END DO  ! i loop
    END DO  ! j

    ! ... and update the diagnostics.
    IF ( l_theta_incr_conv ) THEN
      DO j = 1, rows
        DO i = 1, row_length
          theta_incr_diag_conv(i,j,k) = theta_incr_diag_conv(i,j,k)     &
                + theta_inc_pc2(i,j,1)
        END DO
      END DO
    END IF
    IF ( l_T_incr_conv ) THEN
      DO j = 1, rows
        DO i = 1, row_length
          t_incr_diag_conv(i,j,k) = t_incr_diag_conv(i,j,k)     &
                + exner_theta_levels(i,j,k) * theta_inc_pc2(i,j,1)
        END DO
      END DO
    END IF
    IF ( l_q_incr_conv ) THEN
      DO j = 1, rows
        DO i = 1, row_length
          q_incr_diag_conv(i,j,k) = q_incr_diag_conv(i,j,k)     &
                                    + q_inc_pc2(i,j,1)
        END DO
      END DO
    END IF
    IF ( l_qcl_incr_conv ) THEN
      DO j = 1, rows
        DO i = 1, row_length
          qcl_incr_diag_conv(i,j,k)= qcl_incr_diag_conv(i,j,k)  &
                                     + qcl_inc_pc2(i,j,1)
        END DO
      END DO
    END IF
    ! Not updated   qcf_incr_diag_conv(i,j,k)= qcf_incr_diag_conv(i,j,k)
    IF ( l_cfl_incr_conv ) THEN
      DO j = 1, rows
        DO i = 1, row_length
          cf_liquid_incr_diag_conv(i,j,k)                        &
              = cf_liquid_incr_diag_conv(i,j,k) + cfl_inc_pc2(i,j,1)
        END DO
      END DO
    END IF
    ! Not updated   cf_frozen_incr_diag_conv(i,j,k)
    !    &        = cf_frozen_incr_diag_conv(i,j,k)
    IF ( l_bcf_incr_conv ) THEN
      DO j = 1, rows
        DO i = 1, row_length
          bulk_cf_incr_diag_conv(i,j,k)                          &
              = bulk_cf_incr_diag_conv(i,j,k)   + bcf_inc_pc2(i,j,1)
        END DO
      END DO
    END IF



    ! Add PC2 increments from convection to physics_tendencies_mod 
    ! arrays
    IF (l_retain_conv_tendencies) THEN
      ! Add SPT increments for T
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dt_conv(i,j,k) = dt_conv(i,j,k) +                 &
                       exner_theta_levels(i,j,k) * theta_inc_pc2(i,j,1)
        END DO
      END DO

      ! Add SPT increments for q
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dq_conv(i,j,k) = dq_conv(i,j,k) + q_inc_pc2(i,j,1)
        END DO
      END DO

    END IF
    ! End adding increments to physics_tendencies_mod 
    ! from PC2 inc from conv.

  END IF  ! L_q_interact_if1
END DO  ! k
!$OMP  END DO

DEALLOCATE ( bcf_below )
DEALLOCATE ( bcf_above )

DEALLOCATE ( t_earliest )
DEALLOCATE ( q_earliest )
DEALLOCATE ( qcl_earliest )
DEALLOCATE ( cfl_earliest )
DEALLOCATE ( cff_earliest )
DEALLOCATE ( bcf_earliest )
DEALLOCATE ( t_inc_latest )

DEALLOCATE ( theta_inc_pc2 )
DEALLOCATE ( q_inc_pc2 )
DEALLOCATE ( qcl_inc_pc2 )
DEALLOCATE ( cfl_inc_pc2 )
DEALLOCATE ( bcf_inc_pc2 )

!$OMP  END PARALLEL

IF (l_full_zero) THEN
  DEALLOCATE (full_zero )
END IF
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!----------------------------------------------------------------------
RETURN
END SUBROUTINE pc2_from_conv_ctl

END MODULE pc2_from_conv_ctl_mod

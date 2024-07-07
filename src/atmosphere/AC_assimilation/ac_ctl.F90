! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine AC_CTL   -----------------------------------------------
!
!   programming standard : unified model documentation paper No 3
!
!   Logical components covered : P1
!
!   Project task : P0
!
!   Arguments
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC assimilation

! AC scheme variables are documented in file acp_namel.F90

MODULE ac_ctl_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='AC_CTL_MOD'

CONTAINS

SUBROUTINE ac_ctl(int18, p_field, theta_star_no_halo,             &
 q_star, qcl_star, qcf_star,                                      &
 obs_flag,obs,obs_flag_len,obs_len,                               &
 p, p_theta_levels, exner_theta_levels,                           &
 fv_cos_theta_latitude,                                           &
 cf_area, cf_bulk, cf_liquid, cf_frozen,                          &
 pstar, ntml, cumulus,                                            &
 STASHwork,lambda_p,phi_p,                                        &
                  l_mr_physics, icode, cmessage)

USE atm_fields_bounds_mod

USE planet_constants_mod, ONLY: cp

USE ac_diagnostics_mod

USE cloud_inputs_mod, ONLY: rhcrit, i_cld_area,  i_cld_vn
USE pc2_constants_mod, ONLY: acf_brooks, i_cld_pc2
USE water_constants_mod, ONLY: lc

USE missing_data_mod, ONLY: rmdi, imdi

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParParams, ONLY: peast, pwest, pnorth, psouth
USE Control_Max_Sizes
USE comobs_mod, ONLY: nobtypmx
USE ac_mod, ONLY: ac

USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, stash_pseudo_levels, num_stash_pseudo, stindex, stlist,        &
    num_stash_levels, stash_levels, si, sf
USE nlstcall_mod, ONLY: ltimer

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, bl_levels, len1_lookup,                 &
    len_dumphist, len_fixhd, model_levels,                             &
    mpp_len1_lookup, row_length, rows, theta_field_size

USE ac_control_mod

USE errormessagelength_mod, ONLY: errormessagelength
      
USE model_time_mod, ONLY: &
    secs_per_stepim, stepim

USE model_domain_mod, ONLY: model_type, mt_lam, mt_cyclic_lam

USE ls_acf_brooks_mod, ONLY: ls_acf_brooks
USE ls_arcld_mod, ONLY: ls_arcld
USE pc2_assim_mod, ONLY: pc2_assim
IMPLICIT NONE

INTEGER ::    int18        ! Dummy variable for STASH_MAXLEN(18)
INTEGER ::    p_field
INTEGER ::    icode        ! Return code : 0   = Normal Exit
!                          !             : > 0 = Error

REAL :: theta_star_no_halo(row_length,rows,model_levels) ! latest theta
REAL :: q_star(p_field,model_levels)           ! specific humidity,
REAL :: qcl_star(p_field,model_levels)         ! cloud liquid,
REAL :: qcf_star(p_field,model_levels)         ! cloud ice content
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message if return code >0
INTEGER :: obs_flag_len,obs_len
INTEGER :: obs_flag(obs_flag_len)
REAL    :: obs(obs_len)
LOGICAL, INTENT(IN)::                                             &
 l_mr_physics            ! Use mixing ratio (if code available)

REAL, INTENT (INOUT) ::                                           &
  p(pdims_s%i_start:pdims_s%i_end,                                &
    pdims_s%j_start:pdims_s%j_end,                                &
    pdims_s%k_start:pdims_s%k_end+1)                              &

, p_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
                 tdims_s%j_start:tdims_s%j_end,                   &
                 tdims_s%k_start:tdims_s%k_end)                   &

, exner_theta_levels(tdims_s%i_start:tdims_s%i_end,               &
                     tdims_s%j_start:tdims_s%j_end,               &
                     tdims_s%k_start:tdims_s%k_end)               &

, fv_cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,            &
                        tdims_s%j_start:tdims_s%j_end)            &

, cf_area(tdims%i_start:tdims%i_end,                              &
          tdims%j_start:tdims%j_end,                              &
                        tdims%k_end)                              &

, cf_bulk(tdims_l%i_start:tdims_l%i_end,                          &
          tdims_l%j_start:tdims_l%j_end,                          &
          tdims_l%k_start:tdims_l%k_end)                          &

, cf_liquid(tdims_l%i_start:tdims_l%i_end,                        &
            tdims_l%j_start:tdims_l%j_end,                        &
            tdims_l%k_start:tdims_l%k_end)                        &

, cf_frozen(tdims_l%i_start:tdims_l%i_end,                        &
            tdims_l%j_start:tdims_l%j_end,                        &
            tdims_l%k_start:tdims_l%k_end)                        &

, pstar(pdims%i_start:pdims%i_end,                                &
        pdims%j_start:pdims%j_end)

INTEGER, INTENT(INOUT) ::                                        &
 ntml(pdims%i_start:pdims%i_end,                                 &
      pdims%j_start:pdims%j_end)

LOGICAL, INTENT(INOUT) ::                                        &
 cumulus(tdims%i_start:tdims%i_end,                              &
         tdims%j_start:tdims%j_end)

REAL :: lambda_p(tdims_l%i_start:tdims_l%i_end)
REAL :: phi_p   (tdims_l%i_start:tdims_l%i_end,                      &
              tdims_l%j_start:tdims_l%j_end)

REAL :: lcrcp

! PC2 options are only important if i_cld_vn = i_cld_pc2
! Seek advice from the PC2 team before altering these parameters
! from .true., you will need to have put in place large amounts
! of extra code first.
LOGICAL,PARAMETER:: L_pc2_cond=.TRUE.  ! Do condensation and
! liquid cloud fraction changes.
LOGICAL,PARAMETER:: L_pc2_cfl =.TRUE.  ! Do liquid cloud
! fraction changes. This requires that the condensation as a
! result of assimilation is calculated directly.
! Note: One must not try to run with condensation on but liquid cloud
! fraction off with this code.
LOGICAL,PARAMETER:: L_pc2_cff =.TRUE.  ! Do ice cloud fraction
! changes. This requires that the ice increment from assimilation is
! calculated directly.

!  Dynamically allocated area for stash processing
REAL :: stashwork(*)


INTEGER ::                                                        &
       stashmacro_tag,                                            &
                                 ! STASHmacro tag number
       mdi,                                                       &
                                 ! Missing data indicator
       k, error                  ! do loop variable/ error

INTEGER :: j                        ! DO Loop Variable.
INTEGER :: i, iind                  ! temporary scalars
INTEGER :: im_index                 ! Internal model index
INTEGER :: rhc_row_length, rhc_rows   ! rhcrit dimensions
INTEGER :: int_p                      ! field pointers
INTEGER :: j_start, j_end             ! start and end indices
INTEGER :: i_start, i_end             ! for processors
INTEGER :: ji                         ! 2D array index for halo
                                   ! i/j variables
INTEGER :: levels_per_level, large_levels ! needed by ls_arcld
REAL :: rhcpt (1,1,model_levels)        ! rhcrit array
REAL :: p_layer_centres(row_length, rows, 0:model_levels)
!          pressure at layer centres. Same as p_theta_levels
!          except bottom level = pstar, and at top = 0
REAL :: p_layer_boundaries(row_length, rows, 0:model_levels)
!                pressure at layer boundaries. Same as p except at
!                bottom level = pstar, and at top = 0.
REAL :: exner_layer_centres(row_length, rows, model_levels)
!          exner at layer centres. Same as exner_theta_levels
REAL :: work(p_field,model_levels)    ! array for large-scale
                                   ! latent heating
                                   ! or ls_cld dummy output
REAL :: work2(p_field,model_levels)   ! convective heating
                                   ! or ls_cld dummy output
REAL :: theta_star(p_field,model_levels) ! non-halo theta
REAL :: bulk_cloud_nohalo(row_length,rows,model_levels)

REAL :: dummy_field(p_field)

! These variables are for PC2 calculations
REAL, ALLOCATABLE::                             &
           q_work (:,:,:)                                                &
                               ! Vapour work array (kg kg-1)
,          qcl_work(:,:,:)                                               &
                               ! Liquid work array (kg kg-1)
,          qcf_work(:,:,:)                                               &
                               ! Ice work array    (kg kg-1)
,          t_work(:,:,:)                                                 &
                               ! Temperature/theta work array
                               ! please see inline comments (K)
,          bulk_cloud_fraction(:,:,:)                                    &
                               ! Bulk cloud fraction (no units)
,          cloud_fraction_liquid(:,:,:)                                  &
                               ! Liquid cloud fraction  "
,          cloud_fraction_frozen(:,:,:)                                  &
                               ! Ice cloud fraction     "
,          delta_q(:,:,:)                                                &
                               ! Change in q to force condensation
,          delta_qcl(:,:,:)                                              &
                               ! Change in qcl to force condensation
,          delta_qcf(:,:,:)                                              &
                               ! Change in qcf
,          delta_t(:,:,:)                                                &
                               ! Change in t to force condensation
,          pc2_work(:,:,:)     ! A work array for PC2.

! Declare allocatable arrays for passing cloud fractions
! to LS_ACF_Brooks
REAL, ALLOCATABLE::                            &
 cf_bulk_nohalo(:,:,:), cf_liquid_nohalo(:,:,:), cf_frozen_nohalo(:,:,:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='AC_CTL'
!
! Start Routine
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

lcrcp=lc/cp

DO j=1,p_field
  dummy_field(j) =0.0
END DO
!
!  1.0 Get address for each field from its STASH section/item code
!      and STASHmacro tag  (searching only on STASHmacro tag)
mdi            = imdi
stashmacro_tag = 30

! Initialise STASHWORK for section 18.
DO j = 1, int18
  stashwork(j) = rmdi

END DO

! Create local versions of exner and p on layer centres/boundaries.
! These local arrays have NO halo.
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    p_layer_centres(i,j,0) = pstar(i,j)
    p_layer_boundaries(i,j,0) = pstar(i,j)
  END DO
END DO

DO k = 1, pdims%k_end - 1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      p_layer_boundaries(i,j,k) = p(i,j,k+1)
      p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
      exner_layer_centres(i,j,k)=exner_theta_levels(i,j,k)
    END DO
  END DO
END DO
k=model_levels
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    p_layer_boundaries(i,j,k) = 0.0
    p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
    exner_layer_centres(i,j,k)=exner_theta_levels(i,j,k)
  END DO
END DO

!  1.5 large scale rainfall rate LSRR stored in ac_diagnostics module

IF (lsrr(1)  ==  rmdi) THEN
  icode    = 4203
  cmessage = "AC_CTL: large scale rainfall rate not available"
  WRITE(umMessage,*)'AC_CTL 4203 ',icode,cmessage
  CALL umPrint(umMessage,src='ac_ctl')
END IF

IF (icode  >   0) GO TO 9999

!  1.6 large scale snowfall rate LSSR stored in ac_diagnostics module

IF (lssr(1)  ==  rmdi) THEN
  icode    = 4204
  cmessage = "AC_CTL: large scale snowfall rate not available"
  WRITE(umMessage,*)'AC_CTL 4204 ',icode,cmessage
  CALL umPrint(umMessage,src='ac_ctl')
END IF

IF (icode  >   0) GO TO 9999

IF (use_conv_in_mops) THEN

  !  1.7 convective rainfall rate CVRR stored in ac_diagnostics module

  IF (cvrr(1)  ==  rmdi) THEN
    icode    = 5205
    cmessage = "AC_CTL: convective rainfall rate not available"

  END IF

  IF (icode  >   0) GO TO 9999

  !  1.8 convective snowfall rate CVSR stored in ac_diagnostics module

  IF (cvsr(1)  ==  rmdi) THEN
    icode    = 5206
    cmessage = "AC_CTL: convective snowfall rate not available"

  END IF

  IF (icode  >   0) GO TO 9999

  !  1.10 convective cloud cover on each model level CONVCC stored in module

  IF (convcc(1,1)  ==  rmdi) THEN
    icode    = 5212
    cmessage = "AC_CTL: convective cloud amount not available"

  END IF

  IF (icode  >   0) GO TO 9999

  !  1.11 bulk cloud fraction (liquid+ice) after large scale cloud CF_LSC

  IF (cf_lsc(1,1)  ==  rmdi) THEN
    icode    = 9201
    cmessage = "AC_CTL: cf after large scale cloud not available"
  END IF

  IF (icode  >   0) GO TO 9999

END IF  !(IF USE_CONV_IN_MOPS)

IF ( l_lhn ) THEN
  !  seek convective heating rate and
  !  diagnostics for calculating large-scale latent heating rate

  IF (use_conv_in_mops) THEN
    !  1.13 temperature increment across convection TINC_CVN stored in module

    IF (tinc_cvn(1,1)  ==  rmdi) THEN
      icode    = 5181
      cmessage = "AC_CTL: temp incrs across conv'n not available"
    END IF

    IF (icode  >   0) GO TO 9999

  END IF

  !  1.14 temperature increment across large scale precipitation TINC_PPN

  IF (tinc_ppn(1,1)  ==  rmdi) THEN
    icode    = 4181
    cmessage = "AC_CTL: temp incrs across ls_ppn not available"
  END IF

  IF (icode  >   0) GO TO 9999

  !  1.15 cloud liquid water after large scale cloud QCL_LSC stored in module

  IF (qcl_lsc(1,1)  ==  rmdi) THEN
    icode    = 9206
    cmessage = "AC_CTL: qcl after large scale cloud not available"
  END IF

  IF (icode  >   0) GO TO 9999

  !  1.16 cloud liquid water after advection QCL_ADV stored in ac_diagnostics module

  IF (qcl_adv(1,1)  ==  rmdi) THEN
    icode    = 12254
    cmessage = "AC_CTL: qcl after advection not available"
  END IF

  IF (icode  >   0) GO TO 9999


  ! 1.25  Calculate 'large-scale' latent heating contributions
  !       ----------------------------------------------------
  DO k=1,model_levels
    DO j=1,p_field
      work(j,k) = tinc_ppn(j,k) +                                   &
                lcrcp*( qcl_lsc(j,k)- qcl_adv(j,k)  )
    END DO
  END DO
  !  large scale latent heating currently dT/dt in K/timestep
  !  as is convective heating - convert both to dtheta/dt in K/s
  IF (use_conv_in_mops) THEN
    DO k = 1, model_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          int_p = (j-1) * row_length + i
          work(int_p,k)  =  work(int_p,k) /                          &
             (exner_layer_centres(i,j,k) * secs_per_stepim(atmos_im))
          work2(int_p,k) =  tinc_cvn(int_p,k) /                      &
             (exner_layer_centres(i,j,k) * secs_per_stepim(atmos_im))
        END DO
      END DO
    END DO

  ELSE     ! IF USE_CONV_IN_MOPS is false

    DO k = 1, model_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          int_p = (j-1) * row_length + i
          work(int_p,k)  =  work(int_p,k) /                          &
             (exner_layer_centres(i,j,k) * secs_per_stepim(atmos_im))
          work2(int_p,k) = 0.0
        END DO
      END DO
    END DO
  END IF

ELSE     !  if LHN not selected
  !                 initialise dummy heating rate array to pass to AC
  DO k=1,model_levels
    DO j=1,p_field
      work(j,k) =  0.0
      work2(j,k) =  0.0
    END DO
  END DO

END IF   !  L_LHN

! ----------------------------------------------------------------------
!  2. --- Section 18 Data Assimilation ------

IF (ltimer) THEN
  CALL timer('AC      ', 3)

END IF

im_index=1


!  copy non halo values from theta_star_no_halo to
!  2-d array theta_star
DO k =             1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      int_p = (j-1) * row_length + i
      theta_star(int_p,k) = theta_star_no_halo(i,j,k)
    END DO
  END DO
END DO

IF (i_cld_vn == i_cld_pc2 .AND.                                   &
   (L_pc2_cond .OR. L_pc2_cfl .OR. L_pc2_cff)   ) THEN
  ! Input values needed for the PC2 scheme
  ALLOCATE ( q_work  (tdims%i_start:tdims%i_end,                  &
                      tdims%j_start:tdims%j_end,                  &
                      model_levels) )
  ALLOCATE ( qcl_work(tdims%i_start:tdims%i_end,                  &
                      tdims%j_start:tdims%j_end,                  &
                      model_levels) )
  ALLOCATE ( qcf_work(tdims%i_start:tdims%i_end,                  &
                      tdims%j_start:tdims%j_end,                  &
                      model_levels) )
  ALLOCATE ( t_work  (tdims%i_start:tdims%i_end,                  &
                      tdims%j_start:tdims%j_end,                  &
                      model_levels) )
  ALLOCATE ( pc2_work(tdims%i_start:tdims%i_end,                  &
                      tdims%j_start:tdims%j_end,                  &
                      model_levels) )
  DO k = 1, model_levels
    DO j = tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        int_p = (j-1) * row_length + i
        ! Copy input values into the work arrays
        q_work  (i,j,k) = q_star  (int_p,k)
        qcl_work(i,j,k) = qcl_star(int_p,k)
        qcf_work(i,j,k) = qcf_star(int_p,k)
        t_work  (i,j,k) = theta_star(int_p,k)
        ! t_work is now theta before AC assimilation
        pc2_work(i,j,k) = exner_layer_centres(i,j,k)
        ! pc2_work is now exner before AC assimilation
      END DO
    END DO
  END DO
END IF  ! i_cld_pc2 .and. (L_pc2_cond.or.L_pc2_cfl.or.L_pc2_cff)

IF (use_conv_in_mops) THEN

  CALL ac (                                                         &
    row_length, rows, bl_levels,                                    &
    theta_field_size,                                               &
    stepim(atmos_im), secs_per_stepim(atmos_im),                    &
    exner_layer_centres, pstar,                                     &
    p_layer_centres(1,1,1),                                         &
    theta_star, q_star, qcl_star, qcf_star,                         &
    convcc, lsrr, lssr, cvrr, cvsr,                                 &
    cf_lsc,work2,work, rhcrit,                                      &
    obs_flag,obs,                                                   &
    stindex(1,1,18,im_index),                                       &
    stlist, len_stlist, si(1,18,im_index), sf(1,18),                &
    stashwork, stash_levels,                                        &
    num_stash_levels, stash_pseudo_levels, num_stash_pseudo,        &
    lambda_p,phi_p,                                                 &
    icode, cmessage)
ELSE
  CALL ac (                                                         &
    row_length, rows, bl_levels,                                    &
    theta_field_size,                                               &
    stepim(atmos_im), secs_per_stepim(atmos_im),                    &
    exner_layer_centres, pstar,                                     &
    p_layer_centres(1,1,1),                                         &
    theta_star, q_star, qcl_star, qcf_star,                         &
    work2, lsrr, lssr, dummy_field, dummy_field,                    &
    cf_lsc,work2,work, rhcrit,                                      &
    obs_flag,obs,                                                   &
    stindex(1,1,18,im_index),                                       &
    stlist, len_stlist, si(1,18,im_index), sf(1,18),                &
    stashwork, stash_levels,                                        &
    num_stash_levels, stash_pseudo_levels, num_stash_pseudo,        &
    lambda_p,phi_p,                                                 &
    icode, cmessage)
END IF

! For PC2 we need to calculate the forcing values of qT and TL
IF (i_cld_vn == i_cld_pc2) THEN

  IF (L_pc2_cond .OR. L_pc2_cfl .OR. L_pc2_cff) THEN

    ALLOCATE ( delta_q  (tdims%i_start:tdims%i_end,                  &
                         tdims%j_start:tdims%j_end,                  &
                         model_levels) )
    ALLOCATE ( delta_qcl(tdims%i_start:tdims%i_end,                  &
                         tdims%j_start:tdims%j_end,                  &
                         model_levels) )
    ALLOCATE ( delta_qcf(tdims%i_start:tdims%i_end,                  &
                         tdims%j_start:tdims%j_end,                  &
                         model_levels) )
    ALLOCATE ( delta_t  (tdims%i_start:tdims%i_end,                  &
                         tdims%j_start:tdims%j_end,                  &
                         model_levels) )
    ALLOCATE ( bulk_cloud_fraction(tdims%i_start:tdims%i_end,        &
                         tdims%j_start:tdims%j_end,                  &
                         model_levels)  )
    ALLOCATE ( cloud_fraction_liquid(tdims%i_start:tdims%i_end,      &
                                     tdims%j_start:tdims%j_end,      &
                                     model_levels))
    ALLOCATE ( cloud_fraction_frozen(tdims%i_start:tdims%i_end,      &
                                     tdims%j_start:tdims%j_end,      &
                                     model_levels))
    DO k = 1, model_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          int_p = (j-1) * row_length + i
          ! Calculate the increments across the AC scheme
          delta_q  (i,j,k) =  q_star  (int_p,k) - q_work  (i,j,k)
          delta_qcl(i,j,k) =  qcl_star(int_p,k) - qcl_work(i,j,k)
          delta_qcf(i,j,k) =  qcf_star(int_p,k) - qcf_work(i,j,k)

          ! t_work is currently theta before AC assimilation (K)
          delta_t  (i,j,k) =  exner_layer_centres(i,j,k)            &
                           * (theta_star(int_p,k) - t_work(i,j,k) )
          ! deltat is now change in temperature associated with
          ! AC assimilation minus the adiabatic contribution
          ! associated with the change in pressure. This last
          ! part is dealt with in pc2_pressure_ctl.
          ! Pc2_work currently contains exner before AC assim
          t_work  (i,j,k) =  pc2_work(i,j,k) * t_work(i,j,k)
          ! t_work now contains temperature before AC assimilation

          pc2_work(i,j,k) =  0.0
          ! pc2_work now contains an array of zeros

          ! Now copy cloud fraction information from the prognostics
          ! arrays
          bulk_cloud_fraction  (i,j,k) = cf_bulk(i,j,k)
          cloud_fraction_liquid(i,j,k) = cf_liquid(i,j,k)
          cloud_fraction_frozen(i,j,k) = cf_frozen(i,j,k)

        END DO
      END DO
    END DO

    ! Now force condensation and cloud fraction updates.
    ! PC2 is not set up for prognostic area_cloud_fraction.

    CALL pc2_assim( secs_per_stepim(atmos_im), l_pc2_cfl, l_pc2_cff &
  ,                 l_mr_physics                                  &
  ,                 t_work(1,1,1), bulk_cloud_fraction(1,1,1)       &
  ,                 cloud_fraction_liquid(1,1,1)                    &
  ,                 cloud_fraction_frozen(1,1,1)                    &
  ,                 q_work(1,1,1), qcl_work(1,1,1), qcf_work(1,1,1) &
  ,              p_layer_centres(1:row_length,1:rows,1:model_levels) &
  ,                 delta_t(1,1,1), delta_q(1,1,1), delta_qcl(1,1,1)&
  ,                 delta_qcf(1,1,1),pc2_work(1,1,1)                &
                    )

    ! q_work, qcl_work and t_work have now been updated by the
    ! forcing and condensation terms together.
    ! cloud_fraction_liquid (and _frozen and _bulk) has been
    ! updated by the condensation. qcf_work and p_work are not
    ! updated.

    ! Now copy q_work, qcl_work and t_work variables back
    ! into the inout variables if we haven't a confident
    ! direct estimate of condensation from another source.
    IF (L_pc2_cond) THEN
      DO k =             1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            int_p = (j-1) * row_length + i
            q_star  (int_p,k) = q_work  (i,j,k)
            qcl_star(int_p,k) = qcl_work(i,j,k)
            ! Remember the inout variable is theta, not temp.
            theta_star(int_p,k) = t_work(i,j,k)                     &
                                / exner_layer_centres(i,j,k)
          END DO
        END DO
      END DO
    END IF  ! L_pc2_cond

    ! Update the D1 cloud fractions
    DO k =             1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          cf_bulk(i,j,k) = bulk_cloud_fraction(i,j,k)
          cf_liquid(i,j,k) = cloud_fraction_liquid(i,j,k)
          cf_frozen(i,j,k) = cloud_fraction_frozen(i,j,k)
        END DO
      END DO
    END DO

    IF (i_cld_area == acf_brooks) THEN
      ALLOCATE ( cf_bulk_nohalo  (tdims%i_start:tdims%i_end,      &
                                  tdims%j_start:tdims%j_end,      &
                                              1:tdims%k_end) )
      ALLOCATE ( cf_liquid_nohalo(tdims%i_start:tdims%i_end,      &
                                  tdims%j_start:tdims%j_end,      &
                                              1:tdims%k_end) )
      ALLOCATE ( cf_frozen_nohalo(tdims%i_start:tdims%i_end,      &
                                  tdims%j_start:tdims%j_end,      &
                                              1:tdims%k_end) )

      i_start = 1
      i_end = row_length
      j_start = 1
      j_end = rows
      IF (model_type == mt_lam) THEN
        IF (at_extremity(PSouth)) j_start = 2
        IF (at_extremity(PNorth)) j_end = rows-1
        IF (at_extremity(PWest)) i_start = 2
        IF (at_extremity(PEast)) i_end = row_length-1
      END IF
      IF (model_type == mt_cyclic_lam) THEN
        IF (at_extremity(PSouth)) j_start = 2
        IF (at_extremity(PNorth)) j_end = rows-1
      END IF

      ! Place bulk, liquid and frozen cloud fractions in halo-free arrays
      DO k = 1, tdims%k_end
        DO j = j_start, j_end
          DO i = i_start, i_end
            ji = i+halo_i + (j+halo_j-1) * (row_length+2*halo_i)
            cf_bulk_nohalo(i,j,k)  = cf_bulk(i,j,k)
            cf_liquid_nohalo(i,j,k)= cf_liquid(i,j,k)
            cf_frozen_nohalo(i,j,k)= cf_frozen(i,j,k)
          END DO
        END DO
      END DO

      CALL LS_ACF_Brooks (                                      &
           FV_cos_theta_latitude                                &
          ,cf_bulk_nohalo, cf_liquid_nohalo                     &
          ,cf_frozen_nohalo, cumulus                            &
          ,cf_area )

      DEALLOCATE ( cf_bulk_nohalo )
      DEALLOCATE ( cf_liquid_nohalo )
      DEALLOCATE ( cf_frozen_nohalo )

    END IF ! i_cld_area

    ! Now deallocate the arrays
    DEALLOCATE(q_work  )
    DEALLOCATE(qcl_work)
    DEALLOCATE(qcf_work)
    DEALLOCATE(t_work  )
    DEALLOCATE(pc2_work)
    DEALLOCATE(delta_q  )
    DEALLOCATE(delta_qcl)
    DEALLOCATE(delta_qcf)
    DEALLOCATE(delta_t  )
    DEALLOCATE(bulk_cloud_fraction  )
    DEALLOCATE(cloud_fraction_liquid)
    DEALLOCATE(cloud_fraction_frozen)

  END IF !  L_pc2_cond.or.L_pc2_cfl.or.L_pc2_cff

ELSE  ! i_cld_pc2

  !  2.1 call cloud scheme to 'rebalance' thermodynamic fields
  !  calculate Tl and qt
  DO k=1,model_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        int_p = (j-1) * row_length + i
        theta_star(int_p,k) = theta_star(int_p,k)*               &
            exner_layer_centres(i,j,k) - lcrcp*qcl_star(int_p,k)
        q_star(int_p,k) = q_star(int_p,k) + qcl_star(int_p,k)
      END DO
    END DO
  END DO
  ! set up some arguments for ls_cld

  rhc_row_length = 1
  rhc_rows       = 1

  DO k = 1, model_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        bulk_cloud_nohalo(i,j,k) = cf_bulk(i,j,k)
      END DO
    END DO
  END DO
  DO k = 1, model_levels
    rhcpt (1,1,k) = rhcrit(k)
  END DO
  ! Determine number of sublevels for vertical gradient area cloud
  ! Want an odd number of sublevels per level: 3 is hardwired in do loops
  levels_per_level = 3
  large_levels = ((model_levels - 2)*levels_per_level) + 2

  CALL ls_arcld( p_layer_centres,rhcpt(:,:,1:model_levels),       &
                 p_layer_boundaries,                              &
                 rhc_row_length, rhc_rows, bl_levels,             &
                 levels_per_level, large_levels,                  &
                 FV_cos_theta_latitude,                           &
                 ntml, cumulus, l_mr_physics,                   &
                 qcf_star(:,1:),theta_star(:,1:), q_star(:,1:),   &
                 qcl_star(:,1:),                                  &
                 cf_area(:,:,1:), bulk_cloud_nohalo(:,:,1:),      &
                 work2(:,1:), work(:,1:),                         &
                 error)

  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        cf_bulk(i,j,k) = bulk_cloud_nohalo(i,j,k)
      END DO
    END DO
  END DO

  ! "1.24.4"  convert t back to theta
  DO k=1,model_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        int_p = (j-1) * row_length + i
        theta_star(int_p,k) = theta_star(int_p,k)/               &
                        exner_layer_centres(i,j,k)
      END DO
    END DO
  END DO

END IF  ! i_cld_pc2

!  copy back theta_star into theta_star_no_halo
DO k =             1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      int_p = (j-1) * row_length + i
      theta_star_no_halo(i,j,k) = theta_star(int_p,k)
    END DO
  END DO
END DO

IF (ltimer) THEN
  CALL timer('AC      ', 4)

END IF

9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ac_ctl
END MODULE ac_ctl_mod

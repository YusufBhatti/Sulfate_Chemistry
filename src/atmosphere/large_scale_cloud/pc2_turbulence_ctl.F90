! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the effect of changing the width of the PDF used in PC2.
!
SUBROUTINE pc2_turbulence_ctl (                                         &
! Primary fields passed in/out - these will be unchanged if calling from
! AP1 and incremented if calling from microphys_ctl
 t, q, qcl, cf, cfl, cff,                                               &
 p_theta_levels,                                                        &
 dqcl_mp,                                                               &
! diagnostic info
 STASHwork4,                                                            &
! SCM diagnostics switches (dummy in full UM)
 nSCMDpkgs, L_SCMDiags,                                                 &
! increments passed in/out - these are unchanged if calling
! from microphys_ctl but updated if calling from AP1
 T_inc, q_inc, qcl_inc, cf_inc, cfl_inc, l_pc2_prod_qcl_mp)
!
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE cloud_inputs_mod, ONLY: dbsdtbs_turb_0, l_micro_eros,               &
                            i_pc2_erosion_method
USE pc2_constants_mod, ONLY: dbsdtbs_turb_1, pc2eros_hybrid_allfaces
USE atm_fields_bounds_mod, ONLY: tdims, pdims
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE timestep_mod, ONLY: timestep
USE um_parvars,   ONLY: at_extremity
USE ereport_mod,  ONLY: ereport
USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                              &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE s_scmop_mod,  ONLY: default_streams,                                &
                        t_avg, d_wet, d_all, scmdiag_pc2
USE scmoutput_mod,ONLY: scmoutput

USE errormessagelength_mod, ONLY: errormessagelength
USE model_domain_mod, ONLY: model_type, mt_single_column

IMPLICIT NONE
!
! Description: Condensation and cloud fraction changes in the PC2
!              framework as a result of changing the width of the
!              moisture PDF without changing its shape.
!
! Method:      Uses the equations outlined in the PC2 cloud scheme
!              documentation UMDP029A
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Cloud
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! Primary fields passed in
REAL ::                                                               &
 t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tdims%k_end),  &
 q(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tdims%k_end),  &
 qcl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
       tdims%k_end),                                                  &
 cf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
      tdims%k_end),                                                   &
 cfl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
       tdims%k_end),                                                  &
 cff(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
       tdims%k_end),                                                  &
 p_theta_levels(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                pdims%k_end)
!
!
! Diagnostics info
REAL ::                                                               &
 STASHwork4(*)     ! STASH workspace
!
! Additional variables for SCM diagnostics
INTEGER ::                                                            &
 nSCMDpkgs              ! No of SCM diagnostics packages
!
LOGICAL ::                                                            &
 L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages
!
REAL ::                                                               &
 T_inc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
       tdims%k_end),                                                  &
 q_inc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
       tdims%k_end),                                                  &
 qcl_inc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
         tdims%k_end),                                                &
 cf_inc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
        tdims%k_end),                                                 &
 cfl_inc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
         tdims%k_end),                                                &
 dqcl_mp(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
         tdims%k_end)

!
! Local variables
!
CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'PC2_TURBULENCE_CTL'
CHARACTER(LEN=errormessagelength) :: cmessage
!
INTEGER ::                                                            &
 i,j,k,                             & ! Loop counters
 icode,item,im_index
!
INTEGER, PARAMETER :: sect = 4        ! STASH section for microphys
!
REAL ::                                                               &
 T_work(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
        tdims%k_end),                                                 &
 q_work(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
        tdims%k_end),                                                 &
 qcl_work(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
          tdims%k_end),                                               &
 cf_work(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
         tdims%k_end),                                                &
 cfl_work(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
          tdims%k_end),                                               &
 zeros(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tdims%k_end)
!
REAL,ALLOCATABLE::                                                    &
 cf_above(:,:,:), cf_below(:,:,:)
!
LOGICAL ::                                                            &
                       !, INTENT(IN)
 l_pc2_prod_qcl_mp        ! Called to apply dqcl from mp calc

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!
!- End of header
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
!
! Define zeros array
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      zeros(i,j,k) = 0.0
    END DO
  END DO
END DO
!
IF (i_pc2_erosion_method == pc2eros_hybrid_allfaces) THEN
  ALLOCATE(cf_above(tdims%i_start:tdims%i_end,                        &
                    tdims%j_start:tdims%j_end,tdims%k_end))
  ALLOCATE(cf_below(tdims%i_start:tdims%i_end,                        &
                    tdims%j_start:tdims%j_end,tdims%k_end))
  k = 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      cf_above(i,j,k) = cf(i,j,k+1)
      cf_below(i,j,k) = cf(i,j,k)
    END DO
  END DO
  DO k = 2, tdims%k_end-1
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        cf_above(i,j,k) = cf(i,j,k+1)
        cf_below(i,j,k) = cf(i,j,k-1)
      END DO
    END DO
  END DO
  k = tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      cf_above(i,j,k) = cf(i,j,k)
      cf_below(i,j,k) = cf(i,j,k-1)
    END DO
  END DO
ELSE
  ALLOCATE(cf_above(1,1,1))
  ALLOCATE(cf_below(1,1,1))
END IF
!
! Call homogenous forcing routine - this does not add the increments
! and outputs them separately
! DEPENDS ON: pc2_hom_conv
CALL pc2_hom_conv(p_theta_levels,tdims%k_end,timestep,                &
                  t, q, qcl, cf, cfl, cff,                            &
                  zeros, zeros, zeros, zeros, zeros,                  &
                  cf_above, cf_below,                                 &
                  dqcl_mp,                                            &
                  t_work, q_work, qcl_work, cf_work, cfl_work,        &
                  dbsdtbs_turb_0, dbsdtbs_turb_1, l_mr_physics,       &
                  l_pc2_prod_qcl_mp)
!
DEALLOCATE(cf_below)
DEALLOCATE(cf_above)
!
IF (l_micro_eros) THEN
  IF ( l_pc2_prod_qcl_mp ) THEN   
    ! if calling from microphys ctl to apply result of mp calc,
    ! update increment fields
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q_inc(i,j,k)   = q_work(i,j,k)   + q_inc(i,j,k)
          qcl_inc(i,j,k) = qcl_work(i,j,k) + qcl_inc(i,j,k)
          cf_inc(i,j,k)  = cf_work(i,j,k)  + cf_inc(i,j,k)
          cfl_inc(i,j,k) = cfl_work(i,j,k) + cfl_inc(i,j,k)
          T_inc(i,j,k)   = T_work(i,j,k)   + T_inc(i,j,k)
          q(i,j,k)       = q(i,j,k)        + q_work(i,j,k)
          qcl(i,j,k)     = qcl(i,j,k)      + qcl_work(i,j,k)
          cf(i,j,k)      = cf(i,j,k)       + cf_work(i,j,k)
          cfl(i,j,k)     = cfl(i,j,k)      + cfl_work(i,j,k)
          t(i,j,k)       = t(i,j,k)        + t_work(i,j,k)
        END DO
      END DO
    END DO 
  ELSE
    ! if calling from microphys ctl to do erosion, update the main fields
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q(i,j,k)   = q(i,j,k)   + q_work(i,j,k)
          qcl(i,j,k) = qcl(i,j,k) + qcl_work(i,j,k)
          cf(i,j,k)  = cf(i,j,k)  + cf_work(i,j,k)
          cfl(i,j,k) = cfl(i,j,k) + cfl_work(i,j,k)
          t(i,j,k)   = t(i,j,k)   + t_work(i,j,k)
        END DO
      END DO
    END DO 
  END IF
ELSE
  ! if calling from AP1,  update the increment fields
  ! or if called from microphys ctl to apply result of mp calc
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q_inc(i,j,k)   = q_work(i,j,k)   + q_inc(i,j,k)
        qcl_inc(i,j,k) = qcl_work(i,j,k) + qcl_inc(i,j,k)
        cf_inc(i,j,k)  = cf_work(i,j,k)  + cf_inc(i,j,k)
        cfl_inc(i,j,k) = cfl_work(i,j,k) + cfl_inc(i,j,k)
        T_inc(i,j,k)   = T_work(i,j,k)   + T_inc(i,j,k)
        q(i,j,k)       = q(i,j,k)        + q_work(i,j,k)
        qcl(i,j,k)     = qcl(i,j,k)      + qcl_work(i,j,k)
        cf(i,j,k)      = cf(i,j,k)       + cf_work(i,j,k)
        cfl(i,j,k)     = cfl(i,j,k)      + cfl_work(i,j,k)
        t(i,j,k)       = t(i,j,k)        + t_work(i,j,k)
      END DO
    END DO
  END DO
END IF
!
! ----------------------------------------------------------------------
! Output Diagnostics
! ----------------------------------------------------------------------
!
SELECT CASE (model_type)
CASE DEFAULT

  icode = 0 ! Initialise error status
  im_index = 1
  ! ----------------------------------------------------------------------
  ! DIAG.04281 Copy hom_conv qcl incr to stashwork
  ! ----------------------------------------------------------------------
  item = 281
  IF (sf(item,sect)) THEN
    !
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork4(si(item,sect,im_index)),t_work,    &
    tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
    stash_levels,num_stash_levels+1,                                &
    atmos_im,sect,item,                                             &
    icode,cmessage)
    !
    IF (icode /=  0) THEN
      CALL ereport(routinename,icode,cmessage)
      icode = 0
    END IF
  END IF
  !
  ! ----------------------------------------------------------------------
  ! DIAG.04282 Copy hom_conv qcl incr to stashwork
  ! ----------------------------------------------------------------------
  item = 282
  IF (sf(item,sect)) THEN

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork4(si(item,sect,im_index)),q_work,         &
         tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,      &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
         stash_levels,num_stash_levels+1,                                &
         atmos_im,sect,item,                                             &
         icode,cmessage)

    IF (icode /=  0) THEN
      CALL ereport(routinename,icode,cmessage)
      icode = 0
    END IF
  END IF
  !
  ! ----------------------------------------------------------------------
  ! DIAG.04283 Copy hom_conv qcl incr to stashwork
  ! ----------------------------------------------------------------------
  item = 283
  IF (sf(item,sect)) THEN

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork4(si(item,sect,im_index)),qcl_work,  &
    tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
    stash_levels,num_stash_levels+1,                                &
    atmos_im,sect,item,                                             &
    icode,cmessage)

    IF (icode /=  0) THEN
      CALL ereport(routinename,icode,cmessage)
      icode = 0
    END IF
  END IF

  ! ----------------------------------------------------------------------
  ! DIAG.04292 Copy hom_conv qcl incr to stashwork
  ! ----------------------------------------------------------------------
  item = 292
  IF (sf(item,sect)) THEN

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork4(si(item,sect,im_index)),cf_work,   &
    tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
    stash_levels,num_stash_levels+1,                                &
    atmos_im,sect,item,                                             &
    icode,cmessage)

    IF (icode /=  0) THEN
      CALL ereport(routinename,icode,cmessage)
      icode = 0
    END IF
  END IF

  ! ----------------------------------------------------------------------
  ! DIAG.04293 Copy hom_conv cfl incr to stashwork
  ! ----------------------------------------------------------------------
  item = 293
  IF (sf(item,sect)) THEN

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork4(si(item,sect,im_index)),cfl_work,  &
         tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity, &
         stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
         stash_levels,num_stash_levels+1,                           &
         atmos_im,sect,item,                                        &
         icode,cmessage)

    IF (icode /=  0) THEN
      CALL ereport(routinename,icode,cmessage)
      icode = 0
    END IF
  END IF


  !-----------------------------------------------------------------------
  ! SCM PC2 Diagnostics Package
  !-----------------------------------------------------------------------
CASE (mt_single_column)
  IF (l_scmdiags(scmdiag_pc2)) THEN

    !       Stash 4 281
    CALL scmoutput(T_work,'dt_pc2turb',                             &
         'Temperature increment PC2 turbulence','K',                &
         t_avg,d_all,default_streams,'',routinename)

    !       Stash 4 282
    CALL scmoutput(q_work,'dq_pc2turb',                             &
         'Specific humidity increment PC2 turbulence','kg/kg',      &
         t_avg,d_wet,default_streams,'',routinename)

    !       Stash 4,283
    CALL scmoutput(qcl_work,'dqcl_pc2turb',                         &
         'QCL increment PC2 turbulence','kg/kg',                    &
         t_avg,d_wet,default_streams,'',routinename)

    !       Stash 4, 292
    CALL scmoutput(cf_work,'dbcf_pc2turb',                          &
         'Bulk cloud fraction increment PC2 turbulence','Fraction', &
         t_avg,d_wet,default_streams,'',routinename)

    !       Stash 4,293
    CALL scmoutput(cfl_work,'dcfl_pc2turb',                           &
         'Liquid cloud fraction increment PC2 turbulence','Fraction', &
         t_avg,d_wet,default_streams,'',routinename)

  END IF ! scmdiag_pc2

END SELECT ! model_type

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pc2_turbulence_ctl

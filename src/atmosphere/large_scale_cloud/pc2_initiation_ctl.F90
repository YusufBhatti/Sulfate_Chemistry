! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initiate cloud and liquid water within the PC2 Cloud Scheme.

MODULE pc2_initiation_ctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'PC2_INITIATION_CTL_MOD'
CONTAINS

SUBROUTINE pc2_initiation_ctl (                                         &
! Dimensions of Rh crit array
  rhc_row_length, rhc_rows, zlcl_mixed,                                 &

! Model switches
  l_mixing_ratio,                                                       &

! SCM diagnostics switches
  nSCMDpkgs,L_SCMDiags,                                                 &

! Primary fields passed IN/OUT
  t,q,qcl,qcf,cf,cfl,cff,rhts,tlts,qtts,ptts,cf_area,                   &

! Primary fields passed IN
  p,pstar,p_theta_levels,ccb,cumulus,rhcrit,                            &

! Output increments
  t_work,q_work,qcl_work,qcf_work,cf_work,cfl_work,cff_work             &
 )

USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: pdims, tdims, pdims_s,                 &
                                 ScmRowLen, ScmRow
USE cloud_inputs_mod,      ONLY: i_cld_area
USE pc2_constants_mod,     ONLY: acf_cusack
USE level_heights_mod,     ONLY: r_theta_levels

USE s_scmop_mod,           ONLY: default_streams,                     &
                                 t_inst, d_wet, d_all, scmdiag_pc2
USE scmoutput_mod,         ONLY: scmoutput

USE model_domain_mod, ONLY: model_type, mt_single_column

USE pc2_arcld_mod, ONLY: pc2_arcld
USE pc2_checks_mod, ONLY: pc2_checks
USE pc2_checks2_mod, ONLY: pc2_checks2
USE pc2_hom_arcld_mod, ONLY: pc2_hom_arcld
USE pc2_initiate_mod, ONLY: pc2_initiate
IMPLICIT NONE

! Description:
!   Initiate a small amount of cloud fraction and liquid
!   water content for the PC2 Cloud Scheme. Check that moisture
!   variables are consistent with each other.
!
! Method:
!   See the PC2 documentation.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Cloud
!
! Code Description:
!   Language: fortran 77 + common extensions
!   This code is written to UMDP3 v6 programming standards.

! Declarations:
! Arguments with intent in. ie: input variables.

! Model dimensions
INTEGER ::         &
  rhc_row_length,  & ! Dimensions of RHcrit variable
  rhc_rows           ! Dimensions of RHcrit variable
!
LOGICAL ::                            &
  cumulus(tdims%i_start:tdims%i_end,  &  
          tdims%j_start:tdims%j_end), & ! Convection is occurring.
  l_mixing_ratio                        ! Use mixing ratio formulation

! height of lcl in a well-mixed BL (types 3 or 4), 0 otherwise
REAL :: zlcl_mixed(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! Primary fields passed in/out
REAL ::                                                               &
  t(                 tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  q(                 tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  qcl(               tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  qcf(               tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  cf(                tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  cfl(               tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  cff(               tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  cf_area(           tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end)

! Primary fields passed in
REAL ::                                                               &
  p(                 pdims_s%i_start:pdims_s%i_end,                   &
                     pdims_s%j_start:pdims_s%j_end,                   &
                     pdims_s%k_start:pdims_s%k_end+1),                &
  pstar(             pdims%i_start  :pdims%i_end,                     &
                     pdims%j_start  :pdims%j_end),                    &
  p_theta_levels(    pdims%i_start:pdims%i_end,                       &
                     pdims%j_start:pdims%j_end,                       &
                     pdims%k_start:pdims%k_end),                      &
  rhcrit(            rhc_row_length,                                  &
                     rhc_rows,                                        &
                                 1:tdims%k_end),                      &
  rhts(              tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  tlts(              tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
!       TL at start of timestep
    qtts(              tdims%i_start:tdims%i_end,                       &
                       tdims%j_start:tdims%j_end,                       &
                                   1:tdims%k_end),                      &
!       qT at start of timestep
    ptts(              pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end,                       &
                       pdims%k_start:pdims%k_end)
!       Pressure at theta levels at start of timestep

INTEGER ::                                                            &
  ccb(               tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end)

INTEGER ::                                                            &
  nSCMDpkgs             ! No of SCM diagnostics packages

LOGICAL ::                                                            &
  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

! Local variables


! Output increment diagnostics
REAL ::                                                               &
  t_work(            tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  q_work(            tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  qcl_work(          tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  qcf_work(          tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  cf_work(           tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  cfl_work(          tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  cff_work(          tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                                 1:tdims%k_end),                      &
  p_layer_boundaries(pdims%i_start:pdims%i_end,                       &
                     pdims%j_start:pdims%j_end,                       &
                                 0:pdims%k_end),                      &
!       Pressure at layer boundaries. Same as p except at
!       bottom level = pstar, and at top = 0.
    p_layer_centres(   pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end,                       &
                                   0:pdims%k_end)
!       Pressure at layer centres. Same as p_theta_levels
!       except bottom level = pstar, and at top = 0.

INTEGER ::                                                            &
  i,j,k,iScm,jScm,                                                    &
!       Loop counters
    large_levels,                                                       &
!       Total no. of sub-levels being processed by cloud scheme.
!       Currently ((model_levels - 2)*levels_per_level) + 2
    levels_per_level
!       No. of sub-levels being processed by area cloud scheme.
!       Want an odd number of sublevels per level.
!       NB: levels_per_level = 3 is currently hardwired in the do loops

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'PC2_INITIATION_CTL'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL   (KIND=jprb)            :: zhook_handle

! External Functions:

!- End of header

REAL :: t_initinc    (ScmRowLen, ScmRow, tdims%k_end)
REAL :: q_initinc    (ScmRowLen, ScmRow, tdims%k_end)
REAL :: qcl_initinc  (ScmRowLen, ScmRow, tdims%k_end)
REAL :: qcf_initinc  (ScmRowLen, ScmRow, tdims%k_end)
REAL :: cf_initinc   (ScmRowLen, ScmRow, tdims%k_end)
REAL :: cfl_initinc  (ScmRowLen, ScmRow, tdims%k_end)
REAL :: cff_initinc  (ScmRowLen, ScmRow, tdims%k_end)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      ! Work fields are set to starting fields so diagnostics
      ! can be calculated

      q_work(i,j,k)   = q(i,j,k)
      qcl_work(i,j,k) = qcl(i,j,k)
      qcf_work(i,j,k) = qcf(i,j,k)
      cf_work(i,j,k)  = cf(i,j,k)
      cfl_work(i,j,k) = cfl(i,j,k)
      cff_work(i,j,k) = cff(i,j,k)

    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      t_work(i,j,k)   = t(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! Call checking routine

CALL pc2_checks(p_theta_levels,                                       &
    t, cf, cfl, cff, q, qcl, qcf,                                     &
    l_mixing_ratio)


!-----------------------------------------------------------------------
!       SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------
IF ( l_scmdiags(scmdiag_pc2) .AND.                                    &
     model_type == mt_single_column ) THEN

  ! Initialise initinc arrays to calculate net increment from initiation
  DO k=1, tdims%k_end
    DO j=tdims%j_start, tdims%j_end
      jScm = j - tdims%j_start + 1
      DO i=tdims%i_start, tdims%i_end
        iScm = i - tdims%i_start + 1
        q_initinc   (iScm,jScm,k) = q   (i,j,k)
        qcl_initinc (iScm,jScm,k) = qcl (i,j,k)
        qcf_initinc (iScm,jScm,k) = qcf (i,j,k)
        cf_initinc  (iScm,jScm,k) = cf  (i,j,k)
        cfl_initinc (iScm,jScm,k) = cfl (i,j,k)
        cff_initinc (iScm,jScm,k) = cff (i,j,k)
      END DO
    END DO
  END DO

  DO k=1, tdims%k_end
    DO j=tdims%j_start, tdims%j_end
      jScm = j - tdims%j_start + 1
      DO i=tdims%i_start, tdims%i_end
        iScm = i - tdims%i_start + 1
        t_initinc(iScm,jScm,k) = t(i,j,k)
      END DO
    END DO
  END DO

END IF ! scmdiag_pc2 / model_type


! Call initiation routine

! using area cloud parameterisation
IF (i_cld_area == acf_cusack) THEN

  ! set p at layer boundaries.
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      p_layer_boundaries(i,j,0) = pstar(i,j)
      p_layer_centres(i,j,0) = pstar(i,j)
    END DO
  END DO
  DO k = pdims%k_start, pdims%k_end - 1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        p_layer_boundaries(i,j,k) = p(i,j,k+1)
        p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
      END DO
    END DO
  END DO
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      p_layer_boundaries(i,j,pdims%k_end) = 0.0
      p_layer_centres(i,j,pdims%k_end) =                           &
                       p_theta_levels(i,j,pdims%k_end)
    END DO
  END DO

  ! Determine number of sublevels for vertical gradient area cloud
  ! Want an odd number of sublevels per level: 3 is hardwired in do loops
  levels_per_level = 3
  large_levels = ((pdims%k_end - 2)*levels_per_level) + 2

  CALL pc2_arcld(p_layer_centres,p_layer_boundaries,                &
    ccb,cumulus,rhcrit,                                             &
    rhc_row_length,rhc_rows,zlcl_mixed,                             &
    large_levels,levels_per_level,cf_area,                          &
    t,cf,cfl,cff,q,qcl,qcf,rhts,tlts,qtts,ptts,l_mixing_ratio)

ELSE !i_cld_area

  CALL pc2_initiate(p_theta_levels,ccb,cumulus,rhcrit,                &
    tdims%k_end, rhc_row_length,rhc_rows,zlcl_mixed,r_theta_levels,   &
    t,cf,cfl,cff,q,qcl,rhts,l_mixing_ratio)

END IF !i_cld_area



!-----------------------------------------------------------------------
!       SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------
! Update initinc arrays to hold net increment from initiation

IF ( l_scmdiags(scmdiag_pc2) .AND.                                    &
     model_type == mt_single_column ) THEN

  DO k=1, tdims%k_end
    DO j=tdims%j_start, tdims%j_end
      jScm = j - tdims%j_start + 1
      DO i=tdims%i_start, tdims%i_end
        iScm = i - tdims%i_start + 1

        q_initinc   (iScm,jScm,k) = q   (i,j,k) - q_initinc   (iScm,jScm,k)
        qcl_initinc (iScm,jScm,k) = qcl (i,j,k) - qcl_initinc (iScm,jScm,k)
        qcf_initinc (iScm,jScm,k) = qcf (i,j,k) - qcf_initinc (iScm,jScm,k)
        cf_initinc  (iScm,jScm,k) = cf  (i,j,k) - cf_initinc  (iScm,jScm,k)
        cfl_initinc (iScm,jScm,k) = cfl (i,j,k) - cfl_initinc (iScm,jScm,k)
        cff_initinc (iScm,jScm,k) = cff (i,j,k) - cff_initinc (iScm,jScm,k)
      END DO
    END DO
  END DO

  DO k=1, tdims%k_end
    DO j=tdims%j_start, tdims%j_end
      jScm = j - tdims%j_start + 1
      DO i=tdims%i_start, tdims%i_end
        iScm = i - tdims%i_start + 1
        t_initinc(iScm,jScm,k)   = t(i,j,k)   - t_initinc(iScm,jScm,k)
      END DO
    END DO
  END DO



  !         Stash 16,161
  CALL scmoutput(t_initinc,'dt_pc2init',                              &
       'Temperature increment PC2 init','K',                          &
       t_inst,d_all,default_streams,'',routinename)

  !         Stash 16,162
  CALL scmoutput(q_initinc,'dq_pc2init',                              &
       'Specific humidity increment PC2 init','kg/kg',                &
       t_inst,d_wet,default_streams,'',routinename)

  !         Stash 16,163
  CALL scmoutput(qcl_initinc,'dqcl_pc2init',                          &
       'QCL increment PC2 init','kg/kg',                              &
       t_inst,d_wet,default_streams,'',routinename)

  !         Stash 16,164
  CALL scmoutput(qcf_initinc,'dqcf_pc2init',                          &
       'QCF increment PC2 init','kg/kg',                              &
       t_inst,d_wet,default_streams,'',routinename)

  !         Stash 16,172
  CALL scmoutput(cf_initinc,'dbcf_pc2init',                           &
       'Bulk cloud fraction increment PC2 init','Fraction',           &
       t_inst,d_wet,default_streams,'',routinename)

  !         Stash 16,173
  CALL scmoutput(cfl_initinc,'dcfl_pc2init',                          &
       'Liquid cloud fraction increment PC2 init','Fraction',         &
       t_inst,d_wet,default_streams,'',routinename)

  !         Stash 16,174
  CALL scmoutput(cff_initinc,'dcff_pc2init',                          &
       'Frozen cloud fraction increment PC2 init','Fraction',         &
       t_inst,d_wet,default_streams,'',routinename)

END IF ! scmdiag_pc2 / model_type


! Call second checking routine

CALL pc2_checks2(p_theta_levels,rhcrit,                               &
    rhc_row_length,rhc_rows,                                          &
    t, cf, cfl, cff, q, qcl, l_mixing_ratio)

! Call first checking routine again

CALL pc2_checks(p_theta_levels,                                       &
    t, cf, cfl, cff, q, qcl, qcf,                                     &
    l_mixing_ratio)

! using area cloud parameterisation
IF (i_cld_area == acf_cusack) THEN
  CALL pc2_hom_arcld(p_layer_centres,p_layer_boundaries,            &
     large_levels,levels_per_level,                                 &
     cf_area,t,cf,cfl,cff,q,qcl,qcf,                                &
     l_mixing_ratio)
END IF    ! i_cld_area

! Update work array to hold net increment from the above routines

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      q_work(i,j,k)   = q(i,j,k)   - q_work(i,j,k)
      qcl_work(i,j,k) = qcl(i,j,k) - qcl_work(i,j,k)
      qcf_work(i,j,k) = qcf(i,j,k) - qcf_work(i,j,k)
      cf_work(i,j,k)  = cf(i,j,k)  - cf_work(i,j,k)
      cfl_work(i,j,k) = cfl(i,j,k) - cfl_work(i,j,k)
      cff_work(i,j,k) = cff(i,j,k) - cff_work(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      t_work(i,j,k)   = t(i,j,k)   - t_work(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

!-----------------------------------------------------------------------
!       SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------
IF ( l_scmdiags(scmdiag_pc2) .AND.                                    &
     model_type == mt_single_column ) THEN

  !         Stash 16,161
  CALL scmoutput(t_work,'dt_pc2initchk',                              &
       'Temperature increment PC2 init+chks','K',                     &
       t_inst,d_all,default_streams,'',routinename)

  !         Stash 16,162
  CALL scmoutput(q_work,'dq_pc2initchk',                              &
       'Specific humidity increment PC2 init+chks','kg/kg',           &
       t_inst,d_wet,default_streams,'',routinename)

  !         Stash 16,163
  CALL scmoutput(qcl_work,'dqcl_pc2initchk',                          &
       'QCL increment PC2 init+chks','kg/kg',                         &
       t_inst,d_wet,default_streams,'',routinename)

  !         Stash 16,164
  CALL scmoutput(qcf_work,'dqcf_pc2initchk',                          &
       'QCF increment PC2 init+chks','kg/kg',                         &
       t_inst,d_wet,default_streams,'',routinename)

  !         Stash 16,172
  CALL scmoutput(cf_work,'dbcf_pc2initchk',                           &
       'Bulk cloud fraction increment PC2 init+chks','Fraction',      &
       t_inst,d_wet,default_streams,'',routinename)

  !         Stash 16,173
  CALL scmoutput(cfl_work,'dcfl_pc2initchk',                          &
       'Liquid cloud fraction increment PC2 init+chks','Fraction',    &
       t_inst,d_wet,default_streams,'',routinename)

  !         Stash 16,174
  CALL scmoutput(cff_work,'dcff_pc2initchk',                          &
       'Frozen cloud fraction increment PC2 init+chks','Fraction',    &
       t_inst,d_wet,default_streams,'',routinename)

END IF ! scmdiag_pc2 / model_type


! End of routine initiation_ctl

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pc2_initiation_ctl
END MODULE pc2_initiation_ctl_mod

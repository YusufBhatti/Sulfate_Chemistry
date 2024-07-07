! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To set the turbulent mixing coefficients KM and KH
!           (Note: should be used after any vertical interpolation
!                  but before any horizontal interpolation.)

!  Programming standard: UMDP 3

!  Documentation: UMDP No.24

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
MODULE kmkh_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'KMKH_MOD'
CONTAINS

SUBROUTINE kmkh (                                                       &
! IN data
 bl_levels,BL_diag,nSCMDpkgs,L_SCMDiags,                                &
 ntml,cumulus,ntdsc,sml_disc_inv,dsc_disc_inv,                          &
 weight_1dbl,weight_1dbl_rho,                                           &
! INOUT data
 rhokm,rhokh,rhokmz,rhokhz,rhokm_top,rhokh_top                          &
 )

USE atm_fields_bounds_mod, ONLY: pdims,pdims_s,tdims,ScmRowLen, ScmRow, &
     tdims
USE bl_option_mod, ONLY:                                                &
    Keep_Ri_FA, off, on, kprof_cu, except_disc_inv
USE bl_diags_mod, ONLY: strnewbldiag
USE cv_run_mod, ONLY: l_param_conv
USE model_domain_mod, ONLY: model_type, mt_single_column
USE s_scmop_mod,      ONLY: default_streams,                            &
                            t_avg, d_bl, scmdiag_bl
USE scmoutput_mod,    ONLY: scmoutput

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! IN arguments
INTEGER, INTENT(IN) ::                                                  &
 bl_levels

LOGICAL, INTENT(IN) ::                                                  &
 cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! IN flag for Cu in the bl

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER,INTENT(IN) ::                                                   &
 nSCMDpkgs              ! No of SCM diagnostics packages

LOGICAL,INTENT(IN) ::                                                   &
 L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages

! Declaration of new BL diagnostics.
TYPE (strnewbldiag), INTENT(INOUT) :: BL_diag

INTEGER, INTENT(IN) ::                                                  &
 ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                              ! IN Number of model levels in the
!                                   turbulently mixed layer.
   ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                ! IN Top level for turb mixing in
!                                   cloud layer
   sml_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                ! IN Flags for whether discontinuous
   dsc_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                ! IN inversions are diagnosed

REAL, INTENT(IN) ::                                                     &
 weight_1dbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             bl_levels),                                                &
 weight_1dbl_rho(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,   &
             bl_levels)
                            ! IN Weighting applied to 1D BL scheme,
                            !    to blend with Smagorinsky scheme,
                            !    on theta and rho levels
! INOUT arguments
REAL, INTENT(INOUT) ::                                                  &
 rhokmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        2:bl_levels),                                                   &
                            ! INOUT Non-local turbulent mixing
!                                  coefficient for momentum.
   rhokhz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2:bl_levels),                                                 &
                              ! INOUT Non-local turbulent mixing
!                                  coefficient for heat and moisture
   rhokm_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             2:bl_levels),                                              &
                              ! INOUT Non-local top-down turbulent
                              !    mixing coefficient for momentum.
   rhokh_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             2:bl_levels)
                              ! INOUT Non-local top-down turbulent
                              !    mixing coefficient for heat
                              !    and moisture.
REAL,INTENT(INOUT) ::                                                   &
 rhokm(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,     &
       bl_levels),                                                      &
                            ! INOUT Layer K-1 - to - layer K
!                                     turbulent mixing coefficient
!                                     for momentum.
   rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                              ! INOUT Layer K-1 - to - layer K
!                                     turbulent mixing coefficient
!                                     for heat and moisture.
!----------------------------------------------------------------------
!  Define local storage.

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'KMKH'
! Scm diags
REAL :: rhokh_diag(ScmRowLen,ScmRow,bl_levels)
REAL :: rhokm_diag(ScmRowLen,ScmRow,bl_levels)
                           ! diffusivity of heat and momentum kg/(ms)

INTEGER ::                                                              &
 i,j,iScm,jScm,                                                         &
                     ! Loop counter (horizontal field index).
 k             ! Loop counter (vertical level index).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (model_type == mt_single_column) THEN
  DO k=1, bl_levels
    DO j=pdims%j_start, pdims%j_end
      jScm = j - pdims%j_start + 1
      DO i=pdims%i_start, pdims%i_end
        iScm = i - pdims%i_start + 1
        rhokh_diag(iScm,jScm,k) = rhokh(i,j,k)
        rhokm_diag(iScm,jScm,k) = rhokm(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
END IF ! model_type

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,j,i)
IF (Keep_Ri_FA == on) THEN
  !-----------------------------------------------------------------------
  ! Set local K's to zero at the LCL in cumulus and at the
  ! top of a turbulent layer with a well-defined inversion
  !-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF ( ( cumulus(i,j) .OR. sml_disc_inv(i,j) == 1) .AND.          &
             (k == ntml(i,j)+1 .OR. k == ntml(i,j)+2) ) THEN
          rhokh(i,j,k) = 0.0
          rhokm(i,j,k) = 0.0
        END IF

        IF ( dsc_disc_inv(i,j)  ==  1 .AND.                             &
             (k == ntdsc(i,j)+1 .OR. k == ntdsc(i,j)+2) ) THEN
          rhokh(i,j,k) = 0.0
          rhokm(i,j,k) = 0.0
        END IF

      END DO ! P_POINTS,i
    END DO ! P_POINTS,j
  END DO ! BL_LEVELS
!$OMP END DO NOWAIT

ELSE IF (Keep_Ri_FA == except_disc_inv) THEN
  !-----------------------------------------------------------------------
          ! Reduce local K's only at the top of a turbulent
          ! layer with a well-defined inversion
  !-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF ( sml_disc_inv(i,j)  ==  1 .AND.                             &
             (k == ntml(i,j)+1 .OR. k == ntml(i,j)+2) ) THEN
          rhokh(i,j,k) = (1.0-weight_1dbl(i,j,k))*rhokh(i,j,k)
          rhokm(i,j,k) = (1.0-weight_1dbl(i,j,k))*rhokm(i,j,k)
        END IF

        IF ( dsc_disc_inv(i,j)  ==  1 .AND.                             &
             (k == ntdsc(i,j)+1 .OR. k == ntdsc(i,j)+2) ) THEN
          rhokh(i,j,k) = (1.0-weight_1dbl(i,j,k))*rhokh(i,j,k)
          rhokm(i,j,k) = (1.0-weight_1dbl(i,j,k))*rhokm(i,j,k)
        END IF

      END DO ! P_POINTS,i
    END DO ! P_POINTS,j
  END DO ! BL_LEVELS
!$OMP END DO NOWAIT

ELSE
  !-----------------------------------------------------------------------
  ! Set local K's to zero from the LCL in cumulus and from the
  ! top of a turbulent layer with a well-defined inversion
  !-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (cumulus(i,j) .AND. ( (l_param_conv .AND. k >  ntml(i,j))    &
           .OR. (.NOT. l_param_conv .AND. k >= ntml(i,j)) )) THEN
          rhokh(i,j,k)=0.0
          rhokm(i,j,k)=0.0
        END IF

        IF ( dsc_disc_inv(i,j)  ==  1 .AND. k  >   ntdsc(i,j) ) THEN
          rhokh(i,j,k) = 0.0
          rhokm(i,j,k) = 0.0
        END IF

        IF ( sml_disc_inv(i,j)  ==  1 .AND. k  >   ntml(i,j) ) THEN
            !   This also means no local mixing within any DSC layer
          rhokh(i,j,k) = 0.0
          rhokm(i,j,k) = 0.0
        END IF

      END DO ! P_POINTS,i
    END DO ! P_POINTS,j
  END DO ! BL_LEVELS
!$OMP END DO NOWAIT

END IF ! test on Keep_Ri_FA

IF (kprof_cu == off) THEN
  !-----------------------------------------------------------------------
  ! Set non-local K's to zero at the LCL in cumulus layers,
  ! including level NTML if not l_param_conv convection scheme
  !-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        IF (cumulus(i,j) .AND. ( (l_param_conv .AND. k == ntml(i,j)+1)  &
             .OR. (.NOT. l_param_conv .AND.                             &
                         k >= ntml(i,j) .AND. k <  ntml(i,j)+2) )) THEN
          rhokhz(i,j,k)=0.0
          rhokmz(i,j,k)=0.0
          rhokh_top(i,j,k)=0.0
          rhokm_top(i,j,k)=0.0
        END IF
      END DO ! P_POINTS,i
    END DO ! P_POINTS,j
  END DO ! BL_LEVELS
!$OMP END DO NOWAIT

END IF  ! test on kprof_cu

! Need a barrier to ensure all previous possible loops have completed
!$OMP BARRIER
!-----------------------------------------------------------------------
! Save diffusion coefficients for diagnostics
!-----------------------------------------------------------------------
IF (BL_diag%l_rhokmloc) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%rhokmloc(i,j,k)=rhokm(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (BL_diag%l_rhokhloc) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%rhokhloc(i,j,k)=rhokh(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (BL_diag%l_rhokmsurf) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%rhokmsurf(i,j,k)=weight_1dbl(i,j,k)*rhokmz(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (BL_diag%l_rhokhsurf) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%rhokhsurf(i,j,k)=weight_1dbl_rho(i,j,k)*rhokhz(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (BL_diag%l_rhokmsc) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%rhokmsc(i,j,k)=weight_1dbl(i,j,k)*rhokm_top(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (BL_diag%l_rhokhsc) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%rhokhsc(i,j,k)=weight_1dbl_rho(i,j,k)*rhokh_top(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
!-----------------------------------------------------------------------
! Set KM and KH to be the maximum of the local and non-local
! values andstore RHO_KM on P-grid for output.
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      rhokh(i,j,k) = MAX( rhokh(i,j,k) ,                                &
             weight_1dbl_rho(i,j,k)*(rhokhz(i,j,k)+rhokh_top(i,j,k)) )
      rhokm(i,j,k) = MAX( rhokm(i,j,k) ,                                &
             weight_1dbl(i,j,k)*(rhokmz(i,j,k)+rhokm_top(i,j,k)) )

    END DO ! P_POINTS,i
  END DO ! P_POINTS,j
END DO ! BL_LEVELS
!$OMP END DO

!$OMP END PARALLEL

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
IF ( l_scmdiags(scmdiag_bl) .AND.                                       &
     (model_type == mt_single_column) ) THEN

  !     Output the diffusivities.
  CALL scmoutput(rhokm,'momdif',                                        &
       'Diffusivity of momentum','kg/(ms)',                             &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(rhokh,'htdiff',                                        &
       'Diffusivity of heat','kg/(ms)',                                 &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(rhokm_diag,'KM_local',                                 &
       'Diffusivity of momentum','kg/(ms)',                             &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(rhokh_diag,'KH_local',                                 &
       'Diffusivity of heat','kg/(ms)',                                 &
       t_avg,d_bl,default_streams,'',routinename)

  rhokh_diag(:,:,1) = 0.0
  rhokm_diag(:,:,1) = 0.0

  DO k=2, bl_levels
    DO j=pdims%j_start, pdims%j_end
      jScm = j - pdims%j_start + 1
      DO i=pdims%i_start, pdims%i_end
        iScm = i - pdims%i_start + 1

        rhokh_diag(iScm,jScm,k) = weight_1dbl_rho(i,j,k)*rhokhz(i,j,k)
        rhokm_diag(iScm,jScm,k) = weight_1dbl(i,j,k)*rhokmz(i,j,k)
      END DO ! k
    END DO ! i
  END DO ! j

  CALL scmoutput(rhokm_diag,'KM_surf',                                  &
       'Diffusivity of momentum','kg/(ms)',                             &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(rhokh_diag,'KH_surf',                                  &
       'Diffusivity of heat','kg/(ms)',                                 &
       t_avg,d_bl,default_streams,'',routinename)

  DO k=2, bl_levels
    DO j=pdims%j_start, pdims%j_end
      jScm = j - pdims%j_start + 1
      DO i=pdims%i_start, pdims%i_end
        iScm = i - pdims%i_start + 1
        rhokh_diag(iScm,jScm,k) = weight_1dbl_rho(i,j,k)*rhokh_top(i,j,k)
        rhokm_diag(iScm,jScm,k) = weight_1dbl(i,j,k)*rhokm_top(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(rhokm_diag,'KM_top',                                   &
       'Diffusivity of momentum','kg/(ms)',                             &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(rhokh_diag,'KH_top',                                   &
       'Diffusivity of heat','kg/(ms)',                                 &
       t_avg,d_bl,default_streams,'',routinename)

END IF ! scmdiag_bl / model_type


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE kmkh
END MODULE kmkh_mod

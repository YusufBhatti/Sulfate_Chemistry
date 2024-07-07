! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Purpose: Calculate explicit fluxes of TL and QT

!  Programming standard: UMDP 3

!  Documentation: UM Documentation Paper No 24.

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
MODULE ex_flux_tq_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'EX_FLUX_TQ_MOD'
CONTAINS

SUBROUTINE ex_flux_tq (                                                 &
! IN levels etc
  bl_levels,nSCMDpkgs,L_SCMDiags,                                       &
! IN fields
  tl, qw, rdz, rhokh, rhokhz, grad_t_adj, grad_q_adj, rhof2, rhofsc,    &
  ft_nt, fq_nt, ft_nt_dscb, fq_nt_dscb, tothf_zh, tothf_zhsc, totqf_zh, &
  totqf_zhsc,  weight_1dbl, ntml, ntdsc, nbdsc,                         &
! INOUT fields
  ftl, fqw                                                              &
  )

USE atm_fields_bounds_mod, ONLY: pdims, tdims
USE bl_option_mod, ONLY: flux_grad, LockWhelan2006
USE planet_constants_mod, ONLY: cp, grcp
USE s_scmop_mod,   ONLY: default_streams,                               &
    t_avg, d_bl, scmdiag_bl
USE scmoutput_mod, ONLY: scmoutput
USE model_domain_mod, ONLY: model_type, mt_single_column

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

INTEGER, INTENT(IN) ::                                                  &
 bl_levels,                                                             &
                            ! IN No. of atmospheric levels for which
!                                boundary layer fluxes are calculated.
   ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                              ! IN Number of model layers in turbulent
!                                  mixing layer.
   ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! IN Top level for turb mixing in any
!                                  decoupled Sc layer
   nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! IN lowest flux level in DSC layer

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::                                                  &
  nSCMDpkgs             ! No of SCM diagnostics packages

LOGICAL, INTENT(IN) ::                                                  &
  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

REAL, INTENT(IN) ::                                                     &
  tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),   &
                            ! IN Liquid/frozen water temperture (K)
  qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),   &
                            ! IN Total water content (kg/kg)
  rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        bl_levels),                                                     &
                            ! IN Exchange coeffs for scalars
  rhokhz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         2:bl_levels),                                                  &
                            ! IN Non-local turbulent mixing
                            !    coefficient for heat and moisture
  weight_1dbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              bl_levels),                                               &
                            ! IN Weighting applied to 1D BL scheme,
                            !    to blend with Smagorinsky scheme
  rdz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),   &
                            ! IN RDZ(,1) is the reciprocal
                            !     height of level 1, i.e. of the
                            !     middle of layer 1.  For K > 1,
                            !     RDZ(,K) is the reciprocal of the
                            !     vertical distance from level
                            !     K-1 to level K.
  grad_t_adj(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                            ! IN Temperature gradient adjustmenent
                            !    for non-local mixing in unstable
                            !    turbulent boundary layer.
  grad_q_adj(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! IN Humidity gradient adjustment
!                                  for non-local mixing in unstable
!                                  turbulent boundary layer.

REAL, INTENT(IN) ::                                                     &
  ft_nt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels+1),                                                   &
                            ! IN Non-turbulent heat and moisture flux
  fq_nt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels+1)        !    (on rho levels, surface flux(K=1)=0)
REAL, INTENT(IN) ::                                                     &
  rhof2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        2:bl_levels),                                                   &
                            ! IN f2 and fsc term shape profiles
  rhofsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         2:bl_levels)
REAL, INTENT(IN) ::                                                     &
  tothf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                            ! IN Total heat fluxes at inversions
  tothf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &

  totqf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                            ! IN Total moisture fluxes at inversions
  totqf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &

  ft_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                            ! IN Non-turbulent heat and moisture flux
  fq_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            !      at the base of the DSC layer.

! ARGUMENTS WITH INTENT INOUT.
REAL, INTENT(INOUT) ::                                                  &
  ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, bl_levels),  &
                           ! INOUT FTL(,K) contains net turb
!                                   sensible heat flux into layer K
!                                   from below; so FTL(,1) is the
!                                   surface sensible heat, H. (W/m2)
    fqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, bl_levels)
                             ! INOUT Moisture flux between layers
!                                   (kg per square metre per sec).
!                                   FQW(,1) is total water flux
!                                   from surface, 'E'.

! LOCAL VARIABLES.

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'EX_FLUX_TQ'

INTEGER ::                                                              &
  i, j, k

REAL ::                                                                 &
 grad_ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels),                                                   &
                                          ! K*dth/dz
 grad_fqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels),                                                   &
                                          ! K*dq/dz
 non_grad_ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              bl_levels),                                               &
                                              ! Non-gradient flux
 non_grad_fqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              bl_levels),                                               &
                                              ! Non-gradient flux
 f2_ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
                                          ! Heat flux: f2 term
 f2_fqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
                                          ! Moisture flux: f2 term
 fsc_ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         bl_levels),                                                    &
                                          ! Heat flux: fsc term
 fsc_fqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         bl_levels)  ! Moisture flux: fsc term

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
!$OMP DO SCHEDULE(STATIC)
DO k = 1, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      grad_ftl(i,j,k)=0.0
      grad_fqw(i,j,k)=0.0
      non_grad_ftl(i,j,k)=0.0
      non_grad_fqw(i,j,k)=0.0
      f2_ftl(i,j,k) =0.0
      f2_fqw(i,j,k) =0.0
      fsc_ftl(i,j,k)=0.0
      fsc_fqw(i,j,k)=0.0
    END DO
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  !-----------------------------------------------------------------------
  ! 1. "Explicit" fluxes of TL and QW, on P-grid.
  !-----------------------------------------------------------------------
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      grad_ftl(i,j,k)= - rhokh(i,j,k) *                                 &
        ( ( ( tl(i,j,k) - tl(i,j,k-1) ) * rdz(i,j,k) ) + grcp )
      grad_fqw(i,j,k)= - rhokh(i,j,k) *                                 &
            ( qw(i,j,k) - qw(i,j,k-1) ) * rdz(i,j,k)
        !-------------------------------------------------------------
        ! Entrainment fluxes were specified directly in FTL,FQW in
        ! KMKHZ so now add on gradient-dependent fluxes (note that
        ! RHOKH should be zero at these levels).
        !-------------------------------------------------------------
      ftl(i,j,k) = weight_1dbl(i,j,k)*ftl(i,j,k) + grad_ftl(i,j,k)
      fqw(i,j,k) = weight_1dbl(i,j,k)*fqw(i,j,k) + grad_fqw(i,j,k)
        !-------------------------------------------------------------
        !  Add surface-drive gradient adjustment terms to fluxes
        !  within the surface-based mixed layer.
        !-------------------------------------------------------------
      IF (k  <=  ntml(i,j) ) THEN
        non_grad_ftl(i,j,k) = weight_1dbl(i,j,k) *                      &
                                   rhokhz(i,j,k) * grad_t_adj(i,j)
        non_grad_fqw(i,j,k) = weight_1dbl(i,j,k) *                      &
                                   rhokhz(i,j,k) * grad_q_adj(i,j)
        ftl(i,j,k) = ftl(i,j,k) + non_grad_ftl(i,j,k)
        fqw(i,j,k) = fqw(i,j,k) + non_grad_fqw(i,j,k)
      END IF
    END DO
  END DO
END DO
!$OMP END DO
!-----------------------------------------------------------------------
! 2. Lock and Whelan revised non-gradient formulation
!-----------------------------------------------------------------------
IF (flux_grad  ==  LockWhelan2006) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        IF ( k  <=  ntml(i,j) ) THEN
          f2_ftl(i,j,k)  = rhof2(i,j,k)  * tothf_zh(i,j)
          fsc_ftl(i,j,k) = rhofsc(i,j,k) * tothf_zh(i,j)
          ftl(i,j,k)   = ftl(i,j,k) + weight_1dbl(i,j,k) *              &
                     ( f2_ftl(i,j,k) + fsc_ftl(i,j,k) - ft_nt(i,j,k) )

          f2_fqw(i,j,k)  = rhof2(i,j,k)  * totqf_zh(i,j)
          fsc_fqw(i,j,k) = rhofsc(i,j,k) * totqf_zh(i,j)
          fqw(i,j,k)   = fqw(i,j,k) + weight_1dbl(i,j,k) *              &
                     ( f2_fqw(i,j,k) + fsc_fqw(i,j,k) - fq_nt(i,j,k) )
        END IF

        IF ( k  <=  ntdsc(i,j) .AND. k >= nbdsc(i,j) ) THEN

          f2_ftl(i,j,k)  = rhof2(i,j,k)  *                              &
                            ( tothf_zhsc(i,j)-ft_nt_dscb(i,j) )
          fsc_ftl(i,j,k) = rhofsc(i,j,k) *                              &
                            ( tothf_zhsc(i,j)-ft_nt_dscb(i,j) )
          ftl(i,j,k)   = ftl(i,j,k) + weight_1dbl(i,j,k) *              &
                     ( f2_ftl(i,j,k) + fsc_ftl(i,j,k)                   &
                         - ( ft_nt(i,j,k)-ft_nt_dscb(i,j) ) )

          f2_fqw(i,j,k)  = rhof2(i,j,k)  *                              &
                            ( totqf_zhsc(i,j)-fq_nt_dscb(i,j) )
          fsc_fqw(i,j,k) = rhofsc(i,j,k) *                              &
                            ( totqf_zhsc(i,j)-fq_nt_dscb(i,j) )
          fqw(i,j,k)   = fqw(i,j,k) + weight_1dbl(i,j,k) *              &
                     ( f2_fqw(i,j,k) + fsc_fqw(i,j,k)                   &
                         - ( fq_nt(i,j,k)-fq_nt_dscb(i,j) ) )
        END IF
      END DO
    END DO
  END DO
!$OMP END DO
END IF   ! FLUX_GRAD
!$OMP END PARALLEL

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
IF ( l_scmdiags(scmdiag_bl) .AND.                                       &
     model_type == mt_single_column ) THEN

  DO k=1, bl_levels
    DO j=pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        grad_ftl(i,j,k)= cp * grad_ftl(i,j,k)
        non_grad_ftl(i,j,k) = cp * non_grad_ftl(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(grad_ftl,'Grad_ftl',                                   &
       'Down gradient flux of TL','W/m2',                               &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(non_grad_ftl,'Surf_ng_ftl',                            &
       'Surface driven non-gradient flux of TL','W/m2',                 &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(grad_fqw,'Grad_fqw',                                   &
       'Down-gradient flux of QW','kg/m2/s',                            &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(non_grad_fqw,'Surf_ng_fqw',                            &
       'Surface driven non-gradient flux of QW','kg/m2/s',              &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(f2_ftl,'f2_ftl',                                       &
       'FTL: f2 term','W/m2',                                           &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(f2_fqw,'f2_fqw',                                       &
       'FQW: f2 term','kg/m2/s',                                        &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(fsc_ftl,'fSc_ftl',                                     &
       'FTL: fSc term','W/m2',                                          &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(fsc_fqw,'fSc_fqw',                                     &
       'FQW: fSc term','kg/m2/s',                                       &
       t_avg,d_bl,default_streams,'',routinename)

END IF ! scmdiag_bl / model_type


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ex_flux_tq
END MODULE ex_flux_tq_mod

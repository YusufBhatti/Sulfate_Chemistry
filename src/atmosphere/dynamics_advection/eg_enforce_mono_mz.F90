! *****************************COPYRIGHT******************************* 
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine enforce_mono_mz
!
MODULE enforce_mono_mz_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ENFORCE_MONO_MZ_MOD'
!
CONTAINS
SUBROUTINE enforce_mono_PMF(f_h,f_l,ise,ife,jse,jfe,kse,kfe,    &
                                    isr,ifr,jsr,jfr,ksr,kfr,nt, &
                                    l_mono_stringent            )
!======================================================================
! This routine to enforce the PMF monotonicity scheme
! The method needs a high-order and low-order solutions
! The method detect non-monotonic behaviour by looking at product of
! slopes around the target point then if it detect ocillatory (zigzag)
! behaviour replaces the solution with the low-order one.
!
! For details see references:
! M. Zerroukat, J. Comput. Physics, vol.229, pp 9011-9019 (2010)
! M. Zerroukat and T. Allen, J. Comput. Physics,302, pp 285-299 (2015)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!======================================================================

    USE mpp_conf_mod,  ONLY: swap_field_is_scalar

    USE Field_Types, ONLY: fld_type_p
    USE parkind1, ONLY: jpim, jprb       !DrHook
    USE yomhook,  ONLY: lhook, dr_hook   !DrHook      

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ise,ife,jse,jfe,kse,kfe,isr,ifr,   &
                           jsr,jfr,ksr,kfr,nt
    REAL, INTENT(INOUT) :: f_h(ise:ife,jse:jfe,kse:kfe,nt)
    REAL, INTENT(IN)    :: f_l(ise:ife,jse:jfe,kse:kfe,nt)
    LOGICAL, INTENT(IN) :: l_mono_stringent
! local
    INTEGER, PARAMETER :: nmon=2
    REAL :: f_he(isr-nmon:ifr+nmon,jsr-nmon:jfr+nmon,ksr:kfr,nt)
       
    INTEGER :: i,j,k,kk
    REAL    :: t1x, t2x, t3x, t1y, t2y, t3y, t1z, t2z, t3z
    REAL, PARAMETER  :: eps = 1.0e-15
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle
    CHARACTER(LEN=*), PARAMETER :: RoutineName='ENFORCE_MONO_PMF'
   
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,kk)                           &
!$OMP          SHARED(nt, ksr, kfr, jsr, jfr, isr, ifr, f_he, f_h)
DO kk = 1, nt
!$OMP DO SCHEDULE(STATIC)
  DO k = ksr, kfr
    DO j = jsr,jfr
      DO i = isr,ifr
        f_he(i,j,k,kk) = f_h(i,j,k,kk)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END DO
!$OMP END PARALLEL
   
CALL Swap_Bounds(f_he,ifr-isr+1,jfr-jsr+1,nt*(kfr-ksr+1),  &
                  nmon, nmon, fld_type_p, swap_field_is_scalar   )
!
! l_mono_stringent=T: Means the answer at any point HAS to lie inside
!                     the range of the surounding points.                    
!
! l_mono_stringent=F: Allows the solution to be oustside the range of
!                     the neighbouring points IF the slopes suggests the
!                     target point is a resolved turning point 
!        This is a less diffusive option compared to l_mono_stringent=T
!                                    
IF ( l_mono_stringent ) THEN
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,kk,t1x,t1y,t1z)                &
!$OMP          SHARED(nt, ksr, kfr, jsr, jfr, isr, ifr, f_he, f_l, f_h)
  DO kk = 1, nt
!$OMP DO SCHEDULE(STATIC)
    DO k = ksr, kfr
      DO j = jsr,jfr
        DO i = isr,ifr
             
          t1x = (f_he(i  ,j,k,kk)-f_he(i-1,j,k,kk)) *  &
                (f_he(i+1,j,k,kk)-f_he(i  ,j,k,kk))
                 
          t1y = (f_he(i,j  ,k,kk)-f_he(i,j-1,k,kk)) *  &
                (f_he(i,j+1,k,kk)-f_he(i,j  ,k,kk))
                       
          t1z = (f_he(i,j,k,kk)-f_he(i,j,MAX(1,k-1) ,kk)) *  &
                (f_he(i,j,MIN(kfr,k+1),kk)-f_he(i,j,k ,kk))
                                                 
          IF ( (t1x < eps) .OR. (t1y < eps) .OR. (t1z < eps) ) THEN
                       
            f_h(i,j,k,kk) = f_l(i,j,k,kk)
                       
          END IF
                 
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL
ELSE
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,kk,t1x,t1y,t1z,t2x,t2y,t2z,  &
!$OMP                                t3x,t3y,t3z)                       &
!$OMP          SHARED(nt, ksr, kfr, jsr, jfr, isr, ifr, f_he, f_l, f_h)
  DO kk = 1, nt
!$OMP DO SCHEDULE(STATIC)
    DO k = ksr, kfr
      DO j = jsr,jfr
        DO i = isr,ifr
             
          t1x = (f_he(i  ,j,k,kk)-f_he(i-1,j,k,kk)) *  &
                (f_he(i+1,j,k,kk)-f_he(i  ,j,k,kk))
                 
          t2x = (f_he(i  ,j,k,kk)-f_he(i-1,j,k,kk)) *  &
                (f_he(i-1,j,k,kk)-f_he(i-2,j,k,kk))
               
          t3x = (f_he(i+1,j,k,kk)-f_he(i  ,j,k,kk)) *  &
                (f_he(i+2,j,k,kk)-f_he(i+1,j,k,kk))
                       
          t1y = (f_he(i,j  ,k,kk)-f_he(i,j-1,k,kk)) *  &
                (f_he(i,j+1,k,kk)-f_he(i,j  ,k,kk))
                       
          t2y = (f_he(i,j  ,k,kk)-f_he(i,j-1,k,kk)) *  &
                (f_he(i,j-1,k,kk)-f_he(i,j-2,k,kk))
               
          t3y = (f_he(i,j+1,k,kk)-f_he(i,j  ,k,kk)) *  &
                (f_he(i,j+2,k,kk)-f_he(i,j+1,k,kk))
                       
          t1z = (f_he(i,j,k,kk)-f_he(i,j,MAX(1,k-1) ,kk)) *   &
                (f_he(i,j,MIN(kfr,k+1),kk)-f_he(i,j,k,kk))
                       
          t2z = (f_he(i,j,k,kk)-f_he(i,j,MAX(ksr,k-1),kk)) *  &
                (f_he(i,j,MAX(ksr,k-1),kk)-f_he(i,j,MAX(ksr,k-2),kk))
               
          t3z = (f_he(i,j,MIN(kfr,k+1),kk)-f_he(i,j,k,kk)) *  &
                (f_he(i,j,MIN(kfr,k+2),kk)-f_he(i,j,MIN(kfr,k+1),kk))
                                 
          IF ( ( t1x < eps .AND. (t2x < eps .OR. t3x < eps ) )  .OR.  &
               ( t1y < eps .AND. (t2y < eps .OR. t3y < eps ) )  .OR.  &
               ( t1z < eps .AND. (t2z < eps .OR. t3z < eps ) ) ) THEN
                                                                           
             f_h(i,j,k,kk) = f_l(i,j,k,kk)
                       
          END IF
                 
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE enforce_mono_PMF
!====================================================================
SUBROUTINE enforce_conserv_OCF(f_h,f_l,f_n,rho_n,rho_np1,      &
                               ise,ife,jse,jfe,kse,kfe,        &
                               isr,ifr,jsr,jfr,ksr,kfr,nt      )
!====================================================================
! This routine to enforce the OCF conservative scheme
! The method needs a high-order and low-order solutions
! The method find an optimal solution between the high and low that
! satisfies the conservation constraint:
!
!    SUM{f_h * rho_np1 * volumes} = SUM{f_n * rho_n * volumes}
! 
! See details in:
! M. Zerroukat and T. Allen, J. Comput. Physics, 302, pp 285-299 (2015)
!
!====================================================================
                               
    USE Field_Types, ONLY: fld_type_p
    USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims 
    USE parkind1, ONLY: jpim, jprb       !DrHook
    USE yomhook,  ONLY: lhook, dr_hook   !DrHook
   
    USE eg_helmholtz_mod, only : ec_vol
    USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels
    USE global_2d_sums_mod,    ONLY: global_2d_sums
   
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ise,ife,jse,jfe,kse,kfe,isr,ifr,jsr,  &
                           jfr,ksr,kfr,nt
    REAL, INTENT(INOUT) :: f_h(ise:ife,jse:jfe,kse:kfe,nt)
    REAL, INTENT(IN)    :: f_l(ise:ife,jse:jfe,kse:kfe,nt),      &
                           f_n(ise:ife,jse:jfe,kse:kfe,nt)
    REAL, INTENT(IN)  ::   rho_n(pdims_s%i_start:pdims_s%i_end,  &
                                 pdims_s%j_start:pdims_s%j_end,  &
                                 pdims_s%k_start:pdims_s%k_end)
    REAL, INTENT(IN)  :: rho_np1(pdims_s%i_start:pdims_s%i_end,  &
                                 pdims_s%j_start:pdims_s%j_end,  &
                                 pdims_s%k_start:pdims_s%k_end)
! local
         
    INTEGER :: i,j,k,kk
    REAL    :: t1x, t2x, t3x, t1y, t2y, t3y, t1z, t2z, t3z
    REAL, PARAMETER  :: eps = 1.0e-15
      
    REAL :: alfa_za(pdims%k_start:pdims%k_end)
    REAL :: beta_za(pdims%k_start:pdims%k_end)
    REAL :: du(isr:ifr,jsr:jfr,ksr:kfr,nt)
    REAL :: local_sum(isr:ifr,jsr:jfr,nt,4)
    REAL :: w1(nt), w2(nt)
    REAL :: s1, s2, s3, deno, temp
    REAL :: global_sum(nt,4)
    REAL ::   psi_n(tdims%i_start:tdims%i_end, &
                    tdims%j_start:tdims%j_end, &
                    tdims%k_start:tdims%k_end  ) 
    REAL :: psi_np1(tdims%i_start:tdims%i_end, &
                    tdims%j_start:tdims%j_end, &
                    tdims%k_start:tdims%k_end  )
                    
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='ENFORCE_CONSERV_OCF'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
DO k = pdims%k_start+1, pdims%k_end
      alfa_za(k) = (   eta_rho_levels(k) - eta_theta_levels(k-1) ) /   &
                   ( eta_theta_levels(k) - eta_theta_levels(k-1) )
      beta_za(k) = 1.0 - alfa_za(k) 
END DO
alfa_za(pdims%k_start) = 1.0
beta_za(pdims%k_start) = 0.0 

!======================================================================
! Compute the 3D array psi(i,j,k) = rho*volume*average
!   so mass_tracer = SUM[psi*tracer] = SUM[av(tracer)*rho*vol]
!======================================================================

psi_n(:,:,tdims%k_start) = 0.0
psi_np1(:,:,tdims%k_start) = 0.0
  
DO k = tdims%k_start + 1, tdims%k_end - 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      psi_n(i,j,k)  =  rho_n(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )      &
                    +  rho_n(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)
      psi_np1(i,j,k)  =  rho_np1(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )  &
                      +  rho_np1(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)
    END DO
  END DO
END DO

k = tdims%k_end
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
      psi_n(i,j,k)   =     rho_n(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
    psi_np1(i,j,k)   =   rho_np1(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
  END DO
END DO

!====================================================
! compute some global sums needed for mass conservation
!====================================================

DO kk = 1, nt
  DO i = isr, ifr
    DO j = jsr, jfr
      DO k = ksr+1, kfr
        du(i,j,k,kk) = psi_np1(i,j,k) * (f_h(i,j,k,kk)-f_l(i,j,k,kk))
      END DO
    END DO
  END DO
END DO
    
global_sum = 0.0
local_sum(:,:,:,:) = 0.0 

DO kk = 1, nt
  DO i = isr, ifr
    DO j = jsr, jfr
      DO k = ksr+1, kfr
        local_sum(i,j,kk,1) = local_sum(i,j,kk,1)  +  &
                             psi_n(i,j,k)*f_n(i,j,k,kk)
      END DO
    END DO
  END DO
END DO 
       
DO kk = 1, nt
  DO i = isr, ifr
    DO j = jsr, jfr
      DO k = ksr+1, kfr               
        local_sum(i,j,kk,2) = local_sum(i,j,kk,2) +  &
                          psi_np1(i,j,k)*f_l(i,j,k,kk)
        IF (du(i,j,k,kk) < 0.0 ) THEN
          local_sum(i,j,kk,3) = local_sum(i,j,kk,3) - du(i,j,k,kk)
        ELSE
          local_sum(i,j,kk,4) = local_sum(i,j,kk,4) + du(i,j,k,kk)
        END IF
      END DO
    END DO
  END DO
END DO
    
CALL global_2d_sums(local_sum,ifr-isr+1,jfr-jsr+1,0,0,4*nt,global_sum)
    
DO kk = 1, nt
  s1 = global_sum(kk,4)
  s2 = global_sum(kk,3)
  s3 = global_sum(kk,1) - global_sum(kk,2)
  deno = MAX(s1*s1 + s2*s2, eps)
  w1(kk) = ( s2*s2 + s1*s3 + s2*s1 )/ deno
  w2(kk) = ( s1*s1 - s2*s3 + s2*s1 )/ deno
END DO
DO kk = 1, nt
  IF ( ABS(w1(kk)) < eps .AND. ABS(w2(kk)) < eps ) THEN 
    !the special case of w1=w2=0
    s1 = global_sum(kk,1)/MAX(global_sum(kk,2),eps)
    f_h(:,:,:,kk) = f_h(:,:,:,kk) * s1
  ELSE
    DO k = ksr+1, kfr
      DO j = jsr, jfr
        DO i = isr, ifr
          IF (du(i,j,k,kk) < 0.0 ) THEN
            f_h(i,j,k,kk) = w2(kk)*f_h(i,j,k,kk) +  &
                           (1.0-w2(kk))*f_l(i,j,k,kk)
          ELSE
            f_h(i,j,k,kk) = w1(kk)*f_h(i,j,k,kk) +  &
                           (1.0-w1(kk))*f_l(i,j,k,kk)
          END IF
        END DO
      END DO
    END DO
  END IF     ! ABS(w1) < eps
END DO       ! Loop kk
f_h(:,:,ksr,:) = f_h(:,:,ksr+1,:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE enforce_conserv_OCF

!======================================================================
SUBROUTINE enforce_conserv_OCF2(f_h,f_l,f_n,f_s,lb_flux,        &
                                psi_n, psi_np1,                 &
                                ise,ife,jse,jfe,kse,kfe,        &
                                isr,ifr,jsr,jfr,ksr,kfr,nt      )
!=====================================================================
! This is the same routine as above with exception there are extra
! sources. In this case the mass conservation constraint is:
!
!    SUM{f_h * rho_np1 * volumes} = SUM{f_n * rho_n   * volumes} +
!                                   SUM{f_s * rho_np1 * volumes} +
!                                   lb_flux
!
! lb_flux = refers to net mass through the boundary of the 
! domain (e.g., LAM)
! f_s     = sources added at the end of time step (e.g., physics2) 
!
! See details in:
! M. Zerroukat and T. Allen, J. Comput. Physics, 302, pp 285-299 (2015)
!
!=====================================================================
    USE Field_Types, ONLY: fld_type_p
    USE parkind1, ONLY: jpim, jprb       !DrHook
    USE yomhook,  ONLY: lhook, dr_hook   !DrHook
    USE global_2d_sums_mod,    ONLY: global_2d_sums
   
    IMPLICIT NONE
    INTEGER, INTENT(IN) ::  ise,ife,jse,jfe,kse,kfe,isr,ifr,jsr, &
                            jfr,ksr,kfr,nt
    REAL, INTENT(INOUT) :: f_h(ise:ife,jse:jfe,kse:kfe,nt)
    REAL, INTENT(IN)    :: f_l(ise:ife,jse:jfe,kse:kfe,nt), &
                           f_n(ise:ife,jse:jfe,kse:kfe,nt), &
                           f_s(ise:ife,jse:jfe,kse:kfe,nt)
    REAL, INTENT(IN)    :: psi_n(isr:ifr,jsr:jfr,ksr:kfr),  &
                           psi_np1(isr:ifr,jsr:jfr,ksr:kfr)
    REAL, INTENT(IN)    :: lb_flux(nt)
   
! local
         
    INTEGER :: i,j,k,kk
    REAL    :: t1x, t2x, t3x, t1y, t2y, t3y, t1z, t2z, t3z
    REAL, PARAMETER  :: eps = 1.0e-15
 
    REAL :: du(isr:ifr,jsr:jfr,ksr:kfr,nt)
    REAL :: local_sum(isr:ifr,jsr:jfr,nt,4)
    REAL :: w1(nt), w2(nt)
    REAL :: s1, s2, s3, deno, temp
    REAL :: global_sum(nt,4)
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle
    CHARACTER(LEN=*), PARAMETER :: RoutineName='ENFORCE_CONSERV_OCF2'
   
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
!====================================================
! compute some global sums needed for mass conservation
!====================================================
   
DO kk = 1, nt
  DO i = isr, ifr
    DO j = jsr, jfr
      DO k = ksr+1, kfr
        du(i,j,k,kk) = psi_np1(i,j,k) * (f_h(i,j,k,kk)-f_l(i,j,k,kk))
      END DO
    END DO
  END DO
END DO
    
global_sum = 0.0
local_sum(:,:,:,:) = 0.0 

DO kk = 1, nt
  DO i = isr, ifr
    DO j = jsr, jfr
      DO k = ksr+1, kfr
        local_sum(i,j,kk,1) = local_sum(i,j,kk,1) +        &
                               psi_n(i,j,k)*f_n(i,j,k,kk)  &
                            +  psi_np1(i,j,k)*f_s(i,j,k,kk)
      END DO
    END DO
  END DO
END DO 
       
DO kk = 1, nt
  DO i = isr, ifr
    DO j = jsr, jfr
       DO k = ksr+1, kfr               
          local_sum(i,j,kk,2) = local_sum(i,j,kk,2) +  &
                            psi_np1(i,j,k)*f_l(i,j,k,kk)
          IF (du(i,j,k,kk) < 0.0 ) THEN
            local_sum(i,j,kk,3) = local_sum(i,j,kk,3) - du(i,j,k,kk)
          ELSE
            local_sum(i,j,kk,4) = local_sum(i,j,kk,4) + du(i,j,k,kk)
          END IF
       END DO
    END DO
  END DO
END DO
    
CALL global_2d_sums(local_sum,ifr-isr+1,jfr-jsr+1,0,0,4*nt,global_sum)
DO kk = 1, nt
  global_sum(kk,1) = global_sum(kk,1) + lb_flux(kk)
END DO
    
DO kk = 1, nt
  s1 = global_sum(kk,4)
  s2 = global_sum(kk,3)
  s3 = global_sum(kk,1) - global_sum(kk,2)
  deno = MAX(s1*s1 + s2*s2, eps)
  w1(kk) = ( s2*s2 + s1*s3 + s2*s1 )/ deno
  w2(kk) = ( s1*s1 - s2*s3 + s2*s1 )/ deno
END DO
DO kk = 1, nt
  IF ( ABS(w1(kk)) < eps .AND. ABS(w2(kk)) < eps ) THEN 
    ! the special case of w1=w2=0
    s1 = global_sum(kk,1)/MAX(global_sum(kk,2),eps)
    f_h(:,:,:,kk) = f_h(:,:,:,kk) * s1
  ELSE
    DO k = ksr+1, kfr
      DO j = jsr, jfr
        DO i = isr, ifr
           IF (du(i,j,k,kk) < 0.0 ) THEN
             f_h(i,j,k,kk) = w2(kk)*f_h(i,j,k,kk) +  &
                            (1.0-w2(kk))*f_l(i,j,k,kk)
           ELSE
             f_h(i,j,k,kk) = w1(kk)*f_h(i,j,k,kk) +  &
                            (1.0-w1(kk))*f_l(i,j,k,kk)
           END IF
        END DO
      END DO
    END DO
  END IF     ! ABS(w1) < eps
END DO       ! Loop kk
f_h(:,:,ksr,:) = f_h(:,:,ksr+1,:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)  
END SUBROUTINE enforce_conserv_OCF2

END MODULE enforce_mono_mz_mod

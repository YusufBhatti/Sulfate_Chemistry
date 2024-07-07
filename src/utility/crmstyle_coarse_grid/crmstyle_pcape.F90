! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  Calculate PCAPE

MODULE crmstyle_pcape_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_PCAPE_MOD'

CONTAINS

! ------------------------------------------------------------------------------
! Description:
!  Calculate PCAPE using cloud base and cloud top based on partition
! Using a definition where water loading ignored i.e. Tv = T(1+cvq)
! Could be called for bcu, bcw, ucu
! Needs total tendencies
!
!  PCAPE      = sum (Tv_up - Tv_mean)dp/Tv_mean (J/m3)
!  dPCAPE/dt (J/m3/s)   - See P Bechtold et al 2014
!  cloud base   - lowest cloud 
!  cloud top    - lowest cloud top
!  dPCAPE_BL/dt -  from surface to cloud base sum dTv/dt * dp   
!
!  dilCAPE  =1/g * sum(Tv_up - Tv_mean)dz/Tv_mean 
!  dCAPE/dt = As dPCAPE/dt but not mass weighted. 
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------


SUBROUTINE crmstyle_pcape( )


USE crmstyle_cntl_mod, ONLY:                                             &
   mlevs, iprint, l_bcu, l_bcw, l_ppd, timestep

USE crmstyle_grid_info_mod, ONLY:                                        &
   local_new_x, local_new_y

USE crmstyle_sample_arrays_mod, ONLY:                                     &
  all_th, all_thv, all_rho, all_pstar, all_orog, zlcl,                    &
  all_q, all_qcl, all_qcf, all_qrain, all_qgraup,                         &
  all_ptheta, all_a,                                                      &
  all_dt1, all_dt2, all_dt4, all_dt9, all_dt12,all_dt30,                  &
  all_dq4, all_dq9, all_dq12, all_dq30,                                   &
  all_dqcl4, all_dqcl9, all_dqcl12, all_dqcl30,                           &
  all_dqcf4, all_dqcf3, all_dqcf12, all_dqcf30,                           &
  all_dqrain30, all_dqgr30, all_drho,                                     &
  all_exner, all_t, all_tv,                                               &
  bcu_h_w, bcu_h_u, bcu_h_v, bcu_h_th, bcu_h_thv, bcu_h_rho,              &
  bcu_h_rh, bcu_h_a, bcu_h_dpx, bcu_h_dpy,                                &
  bcu_h_q,   bcu_h_qcl, bcu_h_qcf, bcu_h_qrain, bcu_h_qgraup,             &
  bcu_h_dt1, bcu_h_dt2, bcu_h_dt4, bcu_h_dt9, bcu_h_dt12, bcu_h_dt30,     &
  bcu_h_dq4, bcu_h_dq9, bcu_h_dq12, bcu_h_dq30,                           &
  bcu_h_dqcl4, bcu_h_dqcl9, bcu_h_dqcl12, bcu_h_dqcl30,                   &
  bcu_h_dqcf4, bcu_h_dqcf3, bcu_h_dqcf12, bcu_h_dqcf30,                   &
  bcu_h_dqrain30, bcu_h_dqgr30, bcu_h_drho,                               &
  bcu_pcape, bcu_dpcapedt, bcu_dpcapedt_bl,                               &
  bcu_dilcape, bcu_dcapedt,                                               &
  bcw_h_w, bcw_h_th, bcw_h_thv, bcw_h_rho,                                &
  bcw_h_rh, bcw_h_a, bcw_h_dpx, bcw_h_dpy,                                &
  bcw_h_q,   bcw_h_qcl, bcw_h_qcf, bcw_h_qrain, bcw_h_qgraup,             &
  bcw_h_dt1, bcw_h_dt2, bcw_h_dt4, bcw_h_dt9, bcw_h_dt12, bcw_h_dt30,     &
  bcw_h_dq4, bcw_h_dq9, bcw_h_dq12, bcw_h_dq30,                           &
  bcw_h_dqcl4, bcw_h_dqcl9, bcw_h_dqcl12, bcw_h_dqcl30,                   &
  bcw_h_dqcf4, bcw_h_dqcf3, bcw_h_dqcf12, bcw_h_dqcf30,                   &
  bcw_h_dqrain30, bcw_h_dqgr30, bcw_h_drho,                               &
  bcw_pcape, bcw_dpcapedt, bcw_dpcapedt_bl,                               &
  bcw_dilcape, bcw_dcapedt,                                               &
  ppd_h_w, ppd_h_th, ppd_h_thv, ppd_h_rho,                                &
  ppd_h_rh, ppd_h_a, ppd_h_dpx, ppd_h_dpy,                                &
  ppd_h_q,   ppd_h_qcl, ppd_h_qcf, ppd_h_qrain, ppd_h_qgraup,             &
  ppd_h_dt1, ppd_h_dt2, ppd_h_dt4, ppd_h_dt9, ppd_h_dt12, ppd_h_dt30,     &
  ppd_h_dq4, ppd_h_dq9, ppd_h_dq12, ppd_h_dq30,                           &
  ppd_h_dqcl4, ppd_h_dqcl9, ppd_h_dqcl12, ppd_h_dqcl30,                   &
  ppd_h_dqcf4, ppd_h_dqcf3, ppd_h_dqcf12, ppd_h_dqcf30,                   &
  ppd_h_dqrain30, ppd_h_dqgr30, ppd_h_drho,                               &
  ppd_dpcapedt_bl

USE crmwork_arrays_mod, ONLY: h_theta_sea, dz_theta_sea

USE word_sizes_mod, ONLY: iwp,wp   ! Allows use of 4 byte words to reduce
                                    ! memory

USE planet_constants_mod, ONLY: r, kappa, c_virtual, g


USE missing_data_mod, ONLY: rmdi

USE umPrintMgr                              ! for writing output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE


!---------------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------------

INTEGER :: i,j,k                      ! loop counters


INTEGER(iwp) ::                     &
  k_base(local_new_x,local_new_y)   & ! Cloud base 
 ,k_top(local_new_x,local_new_y)    & ! Cloud top
 ,k_lcl(local_new_x,local_new_y)      ! level less than or equal to zlcl

REAL(wp) ::                            &
  p_rho(local_new_x,local_new_y,mlevs) & ! p half way between theta levels (Pa)
 ,dp(local_new_x,local_new_y,mlevs)      ! estimate of theta layer thickness(Pa)

REAL(wp) ::                                &
  all_dtvdt(local_new_x,local_new_y,mlevs)   ! d(thetav)/dt   (K/s)

REAL(wp) ::    &
  rtimestep    & ! 1/timestep    (/s)
 ,dtvupdt      & ! d(thetav)/dt in buoyant updraughts   (K/s)
 ,t_up         & ! Temperature in updraught (K)
 ,tv_up        & ! Virtual Temperature in updraught (K)
 ,dtvdddt      & ! dTvDD/dt 
 ,t_down       & ! Temperature in downdraught (K)
 ,tv_down        ! Virtual Temperature in downdraught (K)


CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_PCAPE"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

rtimestep = 1./timestep

!-------------------------------------------------------------------------------
! dp required scale with area of grid above surface
!-------------------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(k, i, j) DEFAULT(NONE)                  &
!$OMP& SHARED(mlevs, local_new_y, local_new_x, all_a, all_ptheta, p_rho)
DO k=1,mlevs-1

  DO j=1,local_new_y
    DO i=1,local_new_x
      IF (all_a(i,j,k) > 0.0) THEN
        p_rho(i,j,k) = 0.5*(all_ptheta(i,j,k+1)+all_ptheta(i,j,k))
      END IF
    END DO
  END DO

END DO
!$OMP END PARALLEL DO

! Want dp positive values
k=1
  DO j=1,local_new_y
    DO i=1,local_new_x
      IF (all_a(i,j,k) > 0.0) THEN
        dp(i,j,k) = (all_pstar(i,j) - p_rho(i,j,k))*all_a(i,j,k)
      END IF
    END DO
  END DO

!$OMP PARALLEL DO PRIVATE(k, i, j) DEFAULT(NONE)                  &
!$OMP& SHARED(mlevs, local_new_y, local_new_x, all_a, dp, p_rho, all_pstar)
DO k=2,mlevs-1
  DO j=1,local_new_y
    DO i=1,local_new_x
      IF (all_a(i,j,k) > 0.0) THEN
        IF (all_a(i,j,k-1) > 0.0) THEN
          dp(i,j,k) = (p_rho(i,j,k-1) - p_rho(i,j,k))*all_a(i,j,k)
        ELSE
          dp(i,j,k) = (all_pstar(i,j) - p_rho(i,j,k))*all_a(i,j,k)
        END IF
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

k=mlevs
  DO j=1,local_new_y
    DO i=1,local_new_x
      dp(i,j,k) = (p_rho(i,j,k-1) - all_ptheta(i,j,k))*all_a(i,j,k)
    END DO
  END DO

!-------------------------------------------------------------------------------
! exner, Tv  and dTv/dt
!  p = rhoRT
! dp/dt  = RT drho/dt + R rho dT/dt
!        = R*all_t *all_drho + R all_rho *all_dt30 
!-------------------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(k, i, j) DEFAULT(NONE)                           &
!$OMP& SHARED(mlevs, local_new_y, local_new_x, all_a, all_dtvdt, all_dt30, & 
!$OMP& c_virtual, all_q, all_dq30, all_t, rtimestep)
DO k=1,mlevs
  DO j=1,local_new_y
    DO i=1,local_new_x
      IF (all_a(i,j,k) > 0.0) THEN
        ! Require value per second rather than per timestep
        all_dtvdt(i,j,k) = (all_dt30(i,j,k)*(1.0+ c_virtual*all_q(i,j,k)) + &
                            all_dq30(i,j,k)*c_virtual*all_t(i,j,k))*rtimestep
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! find k lcl
DO j=1,local_new_y
  DO i=1,local_new_x
    k_lcl(i,j)  = 1 
  END DO
END DO
DO k=2,mlevs-1
  DO j=1,local_new_y
    DO i=1,local_new_x
      IF (all_a(i,j,k) > 0.0) THEN
        IF (h_theta_sea(k) <= zlcl(i,j) .AND.                            & 
               h_theta_sea(k+1) > zlcl(i,j)) THEN
          k_lcl(i,j)  = k 
        END IF 
      END IF
    END DO
  END DO
END DO

!-------------------------------------------------------------------------------
! Rest depends on partition 
!-------------------------------------------------------------------------------
IF (l_bcu) THEN
  DO j=1,local_new_y
    DO i=1,local_new_x
      k_base(i,j) = 1
      k_top(i,j)  = 0 
      bcu_pcape(i,j)  = 0.0
      bcu_dilcape(i,j)  = 0.0
      bcu_dcapedt(i,j)  = 0.0
      bcu_dpcapedt(i,j)  = 0.0
      bcu_dpcapedt_bl(i,j)  = 0.0
      IF (l_ppd) THEN
        ppd_dpcapedt_bl(i,j)  = 0.0   
      END IF
    END DO
  END DO


  ! find cloud base and top assuming defined by partition a>0.0 and 
  ! continuously so
  DO k=2,mlevs-1
    DO j=1,local_new_y
      DO i=1,local_new_x
        IF (all_a(i,j,k) > 0.0) THEN
          ! Find first level with a non zero area
          IF (bcu_h_a(i,j,k) > 0.0 .AND. k_base(i,j) == 1 ) THEN
            k_base(i,j) = k
          END IF
          IF (bcu_h_a(i,j,k) > 0.0 .AND. bcu_h_a(i,j,k+1) < 0.0001         &
                 .AND. k_top(i,j) == 0) THEN
            k_top(i,j)  = k 
          END IF 
        END IF
      END DO
    END DO
  END DO

  ! Check that base is not more than one level below the LCL if so alter it to
  ! the level just below.

  DO j=1,local_new_y
    DO i=1,local_new_x
      k=k_base(i,j)             ! Cloud base level
      IF (all_a(i,j,k) > 0.0) THEN
        IF (h_theta_sea(k+1) < zlcl(i,j) ) THEN
          k_base(i,j) = k_lcl(i,j)     
        END IF
      END IF
    END DO
  END DO

  ! Restrict printout 
  IF (iprint == 2) THEN 
    WRITE(umMessage,'(A)') ' bcu_h_a (1,1)       '
    CALL umPrint(umMessage,src=RoutineName)
    DO k = 1,mlevs
      WRITE(umMessage,'(i4,10E10.3)') k,bcu_h_a(1,1,k),bcw_h_a(1,1,k),   &
       all_exner(1,1,k),all_dt30(1,1,k),dp(1,1,k),all_tv(1,1,k),all_dtvdt(1,1,k)
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF
 
! Force cloud base to be level 2 so some BL below
  DO j=1,local_new_y
    DO i=1,local_new_x
      IF (k_base(i,j) == 1) THEN
        k_base(i,j)  = 2
      END IF
    END DO
  END DO

  DO k=1,mlevs
    DO j=1,local_new_y
      DO i=1,local_new_x
        IF (all_a(i,j,k) > 0.0) THEN
          t_up = (bcu_h_th(i,j,k)+all_th(i,j,k))*all_exner(i,j,k)
          tv_up = t_up*(1.0+c_virtual*bcu_h_q(i,j,k))  
          IF ( k >= k_base(i,j) .AND. k <= k_top(i,j)) THEN
            ! Terms 1 and 2 in eq 3  Bechtold et al 2104 
            bcu_pcape(i,j) = bcu_pcape(i,j) +                                &
                 (tv_up - all_tv(i,j,k))*dp(i,j,k)/all_tv(i,j,k)
            bcu_dilcape(i,j) = bcu_dilcape(i,j) + g*(tv_up - all_tv(i,j,k))* &
                                 dz_theta_sea(k)*all_a(i,j,k)/all_tv(i,j,k)

            dtvupdt = ( bcu_h_dt30(i,j,k)*(1.0+ c_virtual*bcu_h_q(i,j,k)) +  &
                            bcu_h_dq30(i,j,k)*c_virtual*t_up)*rtimestep
            bcu_dpcapedt(i,j) = bcu_dpcapedt(i,j) -                          &
                          (all_dtvdt(i,j,k) - dtvupdt)*dp(i,j,k)/all_tv(i,j,k)
            bcu_dcapedt(i,j) = bcu_dcapedt(i,j) -                            &
                                   g*(all_dtvdt(i,j,k) - dtvupdt)*           &
                                  dz_theta_sea(k)*all_a(i,j,k)/all_tv(i,j,k)
          END IF
          ! Cannot estimate term 3 as don't know how the cloud base height
          ! is changing so assuming term 3 is zero.

          IF (k < k_base(i,j)) THEN
            ! Definition   dPCAPE/dt = -(1/T*) sum below pbase (dTv/dt) dp
            ! dp is negative (but we have calculated |dp| a positive value)
            ! Peter Bechtold takes T* as 1K. Peter's paper implies
            ! this should not include the convective part. Here we have 
            ! calculated the total including the convective downdraught
            ! contribution to the BL.
            bcu_dpcapedt_bl(i,j) = bcu_dpcapedt_bl(i,j) +            &
                                     all_dtvdt(i,j,k)*dp(i,j,k)

            IF (l_ppd) THEN
              t_down = (ppd_h_th(i,j,k)+all_th(i,j,k))*all_exner(i,j,k)
!              tv_down = t_down*(1.0+c_virtual*ppd_h_q(i,j,k))  
              dtvdddt = ( ppd_h_dt30(i,j,k)*(1.0+ c_virtual*ppd_h_q(i,j,k)) + &
                           ppd_h_dq30(i,j,k)*c_virtual*t_down)*rtimestep

              ! Integral over just the DD part so adding an area factor
              ppd_dpcapedt_bl(i,j) = ppd_dpcapedt_bl(i,j) +            &
                                        dtvdddt*dp(i,j,k)*ppd_h_a(i,j,k)

            END IF
          END IF
        END IF   ! all_a > 0.0
      END DO
    END DO
  END DO

  ! Restrict printout for checking only 
  IF (iprint == 2) THEN
    WRITE(umMessage,'(A)') ' BCU PCAPE         '
    CALL umPrint(umMessage,src=RoutineName)
    DO j = 1,local_new_y
      WRITE(umMessage,'(20F10.0)') (bcu_pcape(i,j),i=1,local_new_x)
      CALL umPrint(umMessage,src=RoutineName)
    END DO

    WRITE(umMessage,'(A)') ' k_base        '
    CALL umPrint(umMessage,src=RoutineName)
    DO j = 1,local_new_y
      WRITE(umMessage,'(20i10)') (k_base(i,j),i=1,local_new_x)
      CALL umPrint(umMessage,src=RoutineName)
    END DO

    WRITE(umMessage,'(A)') ' k_top        '
    CALL umPrint(umMessage,src=RoutineName)
    DO j = 1,local_new_y
      WRITE(umMessage,'(20i10)') (k_top(i,j),i=1,local_new_x)
      CALL umPrint(umMessage,src=RoutineName)
    END DO

    WRITE(umMessage,'(A)') ' BCU dPCAPE/dt         '
    CALL umPrint(umMessage,src=RoutineName)
    DO j = 1,local_new_y
      WRITE(umMessage,'(20F10.3)') (bcu_dpcapedt(i,j),i=1,local_new_x)
      CALL umPrint(umMessage,src=RoutineName)
    END DO

    WRITE(umMessage,'(A)') ' BCU dPCAPE/dt BL         '
    CALL umPrint(umMessage,src=RoutineName)
    DO j = 1,local_new_y
      WRITE(umMessage,'(20F10.3)') (bcu_dpcapedt_bl(i,j),i=1,local_new_x)
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF
END IF
!-------------------------------------------------------------------------------
!  BCw 
!-------------------------------------------------------------------------------
IF (l_bcw) THEN
  DO j=1,local_new_y
    DO i=1,local_new_x
      k_base(i,j) = 1 
      k_top(i,j)  = 0 
      bcw_pcape(i,j)  = 0.0
      bcw_dilcape(i,j)  = 0.0
      bcw_dpcapedt(i,j)  = 0.0
      bcw_dcapedt(i,j)  = 0.0
      bcw_dpcapedt_bl(i,j)  = 0.0
    END DO
  END DO
  ! find cloud base and top assuming defined by partition a>0.0 and 
  ! continuously so
  DO k=2,mlevs-1
    DO j=1,local_new_y
      DO i=1,local_new_x
        IF (all_a(i,j,k) > 0.0) THEN
          IF (bcw_h_a(i,j,k) > 0.0 .AND. k_base(i,j) == 1 ) THEN
            k_base(i,j) = k
          END IF
          IF (bcw_h_a(i,j,k) > 0.0 .AND. bcw_h_a(i,j,k+1) < 0.0001         &
                 .AND. k_top(i,j) == 0) THEN
            k_top(i,j)  = k 
          END IF 
        END IF
      END DO
    END DO
  END DO
  ! Check that base is not more than one level below the LCL if so alter it to
  ! the level just below.

  DO j=1,local_new_y
    DO i=1,local_new_x
      k=k_base(i,j)             ! Cloud base level
      IF (all_a(i,j,k) > 0.0) THEN
        IF (h_theta_sea(k) < zlcl(i,j) .AND.               &
                h_theta_sea(k+1) < zlcl(i,j)) THEN  
          k_base(i,j) = k_lcl(i,j)     
        END IF
      END IF
    END DO
  END DO

! Force cloud base to be level 2 so some BL below
  DO j=1,local_new_y
    DO i=1,local_new_x
      IF (k_base(i,j) == 1) THEN
        k_base(i,j)  = 2
      END IF
    END DO
  END DO

  DO k=1,mlevs
    DO j=1,local_new_y
      DO i=1,local_new_x
        IF (all_a(i,j,k) > 0.0) THEN
          t_up = (bcw_h_th(i,j,k)+all_th(i,j,k))*all_exner(i,j,k)
          tv_up = t_up*(1.0+c_virtual*bcw_h_q(i,j,k))  
          IF ( k >= k_base(i,j) .AND. k <= k_top(i,j)) THEN
            ! Terms 1 and 2 in eq 3  Bechtold et al 2104 
            bcw_pcape(i,j) = bcw_pcape(i,j) +                                &
              (tv_up - all_tv(i,j,k))*dp(i,j,k)/all_tv(i,j,k)
            bcw_dilcape(i,j) = bcw_dilcape(i,j) + g*(tv_up - all_tv(i,j,k))  &
                                    *dz_theta_sea(k)*all_a(i,j,k)/all_tv(i,j,k)

            dtvupdt = ( bcw_h_dt30(i,j,k)*(1.0+ c_virtual*bcw_h_q(i,j,k)) + &
                            bcw_h_dq30(i,j,k)*c_virtual*t_up)*rtimestep
            bcw_dpcapedt(i,j) = bcw_dpcapedt(i,j) -                         &
                     (all_dtvdt(i,j,k) - dtvupdt)*dp(i,j,k)/all_tv(i,j,k)
            bcw_dcapedt(i,j) = bcw_dcapedt(i,j) -                           &
                                    g*(all_dtvdt(i,j,k) - dtvupdt)*         &
                                dz_theta_sea(k)*all_a(i,j,k)/all_tv(i,j,k)
          END IF
          IF (k < k_base(i,j)) THEN
          ! includes convective downdraughts as total
          ! Should we be removing a downdraught part?
            bcw_dpcapedt_bl(i,j) = bcw_dpcapedt_bl(i,j) +                 &
                                  all_dtvdt(i,j,k)*dp(i,j,k)
          END IF
        END IF
      END DO
    END DO
  END DO

  IF (iprint == 2) THEN
    WRITE(umMessage,'(A)') ' k_base        '
    CALL umPrint(umMessage,src=RoutineName)
    DO j = 1,local_new_y
      WRITE(umMessage,'(20i10)') (k_base(i,j),i=1,local_new_x)
      CALL umPrint(umMessage,src=RoutineName)
    END DO

    WRITE(umMessage,'(A)') ' k_top        '
    CALL umPrint(umMessage,src=RoutineName)
    DO j = 1,local_new_y
      WRITE(umMessage,'(20i10)') (k_top(i,j),i=1,local_new_x)
      CALL umPrint(umMessage,src=RoutineName)
    END DO

    WRITE(umMessage,'(A)') ' BCW PCAPE         '
    CALL umPrint(umMessage,src=RoutineName)
    DO j = 1,local_new_y
      WRITE(umMessage,'(20F10.0)') (bcw_pcape(i,j),i=1,local_new_x)
      CALL umPrint(umMessage,src=RoutineName)
    END DO

    WRITE(umMessage,'(A)') ' BCW dPCAPE/dt         '
    CALL umPrint(umMessage,src=RoutineName)
    DO j = 1,local_new_y
      WRITE(umMessage,'(20F10.0)') (bcw_dpcapedt(i,j),i=1,local_new_x)
      CALL umPrint(umMessage,src=RoutineName)
    END DO

    WRITE(umMessage,'(A)') ' BCW dPCAPE/dt BL         '
    CALL umPrint(umMessage,src=RoutineName)

    DO j = 1,local_new_y
      WRITE(umMessage,'(20F10.3)') (bcw_dpcapedt_bl(i,j),i=1,local_new_x)
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF

END IF

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_pcape

END MODULE crmstyle_pcape_mod

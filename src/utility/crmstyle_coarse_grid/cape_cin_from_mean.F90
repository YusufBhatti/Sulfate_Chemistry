! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculate CAPE and CIN from an undilute ascent

MODULE cape_cin_from_mean_mod

IMPLICIT NONE
!---------------------------------------------------------------------------
! Description: Calculate CAPE and CIN using mean profiles calculated
! for CRMstyle_coarse_grid 

! method: Does an undilute parcel ascent from the near surface or higher
!         to the level of neutral buoyancy. Method being used is a copy
!         of the code in dts_cape (the deep turbulence scheme). 

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid

! code description:
!   language: fortran 90
!   This code is written to umdp3 programming standards version 10.4
!---------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CAPE_CIN_FROM_MEAN_MOD'

CONTAINS

SUBROUTINE  cape_cin_from_mean( )

USE crmstyle_cntl_mod, ONLY:                                         &
   mlevs, iprint, l_cape

USE crmstyle_grid_info_mod, ONLY:                                    &
   local_new_x, local_new_y

USE crmstyle_sample_arrays_mod, ONLY:                                &
      all_pstar, all_orog, all_zh, all_a, cape, cin, zneutral, zlcl, &
      all_th, all_t, all_exner, all_tv, all_q, all_rho, all_ptheta,  &
      all_qcl, all_qcf, zfree, all_tstar

USE crmwork_arrays_mod, ONLY:                                         &
      h_theta_sea, h_rho_sea

USE planet_constants_mod, ONLY: r, cp, kappa, c_virtual, rv, pref,   &
                    repsilon, recip_kappa

USE cv_diag_param_mod, ONLY:                                             &
    a_bolton, b_bolton, c_bolton, d_bolton

USE water_constants_mod, ONLY: lc, lf, tm

USE missing_data_mod, ONLY: rmdi

USE word_sizes_mod, ONLY: iwp,wp       ! Allows use of 4 byte words to reduce
                                        ! memory

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE qsat_mod, ONLY: qsat, qsat_mix

USE umPrintMgr                        ! Required for writing output

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments

!----------------------------------------------------------------------
! Local variables

INTEGER :: i,j,k    & ! loop counters
 ,n_itime           & ! number of points with itime=1
 ,kmin              & ! minimum klev across all points
 ,kmax              & ! maximum ktop across all points
 ,npnts               ! total number of points

INTEGER(iwp) ::                    &
  iitime(local_new_x,local_new_y)   & ! Counter
 ,k_plume(local_new_x,local_new_y) & ! start level for LCL calculation
 ,kstart(local_new_x,local_new_y)  & ! ascent start level
 ,ktop(local_new_x,local_new_y)      ! ascent end level

LOGICAL ::         &
  l_mr_physics     & ! .TRUE. if prognostics are mixing ratios.
 ,l_itime            ! no itime values still zero


REAL(wp) ::                              &
  p_rho(local_new_x,local_new_y,mlevs)   & ! pressure on theta layer boundaries
 ,exner_rho(local_new_x,local_new_y,mlevs) ! Exner on rho levels

! qsat calls expect double precision input variables

REAL ::                             &
  delp(local_new_x,local_new_y)     & ! pressure diff. between one level and
 ,delp_cen(local_new_x,local_new_y) & !  and next (Pa)
 ,pnb(local_new_x,local_new_y)      & ! pressure at level of neutral buoyancy
 ,tparcel(local_new_x,local_new_y)  & ! temperature of parcel  (K)
 ,thparcel(local_new_x,local_new_y) & ! potential temperature of parcel  (K)
 ,p(local_new_x,local_new_y)        & ! pressure of parcel  (Pa)
 ,ptemp(local_new_x,local_new_y)    & ! copy of pressure   (Pa)
 ,ttemp(local_new_x,local_new_y)    & ! copy of temperature   (K)
 ,qvbot(local_new_x,local_new_y)    & ! initial parcel water vapour (kg/kg)
 ,qlbot(local_new_x,local_new_y)    & ! initial parcel liquid water  (kg/kg)
 ,qlval(local_new_x,local_new_y)    & ! liquid water content inside parcel
 ,q0(local_new_x,local_new_y)       & ! local qsat variables
 ,q1(local_new_x,local_new_y)       & ! local qsat variables
 ,q2(local_new_x,local_new_y)       & ! local qsat variables
 ,qnext(local_new_x,local_new_y)    & ! local qsat variables
 ,qvparcel(local_new_x,local_new_y) & ! water vapour in parcel
 ,qcl_plume(local_new_x,local_new_y) & ! Estimate of max allowed plume water
 ,qse(local_new_x,local_new_y)      & ! Qsat for environment
 ,p_lcl(local_new_x,local_new_y)    & ! pressure of LCL (Pa)
 ,T_lcl(local_new_x,local_new_y)    & ! temperature at LCL (K)
 ,exner_lcl(local_new_x,local_new_y) & ! Exner at LCL 
 ,exner_surf(local_new_x,local_new_y)& ! Surface Exner
 ,zstart(local_new_x,local_new_y)     ! Height at start of ascent

REAL ::     &
  rbb2      &  ! 1/bb2
, rdeltat   &  ! 1/deltat
, dp_o_p    &  ! dp /p
, gdz_o_thv    ! dp/(thv rho)

REAL ::          &
  thvparcel      & ! virtual potential temp of parcel 
 ,vap_press      & ! vapour pressure
 ,z_surf         & ! 
 ,zlcl_surf      & ! Height of LCL above surface (m)
 ,thetav         & ! Environmental thetav value
 ,factor         & ! 
 ,capecom        & ! contribution to cape at a given level
 ,vlt            & ! latent heat of deposition or condensation
 ,t0,t1,t2       & ! local temp variables
 ,alpha          & !
 ,beta           & !
 ,bb1            & !
 ,bb2            & ! local variables
 ,drs0           & !
 ,dt             & !
 ,dq             & !
 ,deltat         & ! local variables
 ,qlcrit         & !
 ,qlmin          & ! maximum liquid water content of parcel -- beyond
 ,smallcape      & ! minimum value of cape before negative buoyancy counts 
                   ! as level of neutral buoyancy rather than cin
 ,smallpressure    ! minimum value of pressure -- set to 1 Pa

CHARACTER(LEN=*), PARAMETER :: RoutineName='CAPE_CIN_FROM_MEAN'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------

npnts=local_new_y*local_new_x

l_mr_physics = .FALSE.   ! assumes prognostic input specific quantities


!----------------------------------------------------------------------
! Lifting condensation level 
!----------------------------------------------------------------------

! What level to try to start ascent?
DO j=1,local_new_y
  DO i=1,local_new_x
    k_plume(i,j) = 1
    z_surf = all_orog(i,j) + 0.1*all_zh(i,j)
    DO WHILE (h_theta_sea(k_plume(i,j)) < z_surf)      
      k_plume(i,j) = k_plume(i,j) + 1
    END DO
  END DO
END DO

IF (l_mr_physics) THEN   ! expression for mixing ratio

!$OMP PARALLEL DO PRIVATE(i, j, vap_press) DEFAULT(NONE)                     &
!$OMP& SHARED(local_new_y, local_new_x, all_q, all_ptheta,                   &
!$OMP&  k_plume, repsilon, recip_kappa, T_lcl, p_lcl, all_pstar, all_t)
  DO j=1,local_new_y
    DO i=1,local_new_x

      vap_press = 0.01*all_q(i,j,k_plume(i,j)) * all_ptheta(i,j,k_plume(i,j)) &
                                      / (repsilon+all_q(i,j,k_plume(i,j)) )
      IF (vap_press  >   0.0) THEN
        T_lcl(i,j) = a_bolton + b_bolton/                                     &
                              (c_bolton*LOG(all_t(i,j,k_plume(i,j)))          &
                                         - LOG(vap_press) - d_bolton )

        p_lcl(i,j) = all_ptheta(i,j,k_plume(i,j)) *                           &
                     ( T_lcl(i,j) / all_t(i,j,k_plume(i,j)) )**recip_kappa

      ELSE
        p_lcl(i,j) = all_pstar(i,j)
      END IF

    END DO
  END DO
!$OMP END PARALLEL DO

ELSE          ! expression for specific humidity
!$OMP PARALLEL DO PRIVATE(i, j, vap_press) DEFAULT(NONE)                      &
!$OMP& SHARED(local_new_y, local_new_x, all_q, all_ptheta,                    &
!$OMP&  k_plume, repsilon, recip_kappa, T_lcl, p_lcl, all_pstar, all_t)
  DO j=1,local_new_y
    DO i=1,local_new_x

      vap_press = all_q(i,j,k_plume(i,j)) * all_ptheta(i,j,k_plume(i,j)) &
                                      / (repsilon*100.0 )
      IF (vap_press  >   0.0) THEN
        T_lcl(i,j) = a_bolton + b_bolton/                                     &
                              (c_bolton*LOG(all_t(i,j,k_plume(i,j)))          &
                                         - LOG(vap_press) - d_bolton )

        p_lcl(i,j) = all_ptheta(i,j,k_plume(i,j)) *                           &
                     ( T_lcl(i,j) / all_t(i,j,k_plume(i,j)) )**recip_kappa

      ELSE
        p_lcl(i,j) = all_pstar(i,j)
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! test on l_mr_physics

! Exner and pressure on rho levels
DO j=1,local_new_y
  DO i=1,local_new_x
    exner_lcl(i,j)  = (p_lcl(i,j)/pref)**kappa
    exner_surf(i,j) = (all_pstar(i,j)/pref)**kappa
    p_rho(i,j,1) = all_pstar(i,j)       ! surface value as used for layer cal
  END DO
END DO

! Work out exner at rho from exner at theta
!$OMP PARALLEL DO PRIVATE(k, i, j, factor) DEFAULT(NONE)                   &
!$OMP& SHARED(mlevs, local_new_y, local_new_x, all_a,  exner_rho,          &
!$OMP&  p_rho, h_theta_sea, h_rho_sea, all_exner, kappa, pref, all_pstar)
DO k=2,mlevs-1
  DO j=1,local_new_y
    DO i=1,local_new_x
      IF (all_a(i,j,k) > 0.0 .AND.  all_a(i,j,k-1) > 0.0) THEN  
        factor = (h_theta_sea(k)-h_rho_sea(k))/h_theta_sea(k)
        exner_rho(i,j,k) = all_exner(i,j,k-1)*factor                    &
                             +(1.-factor)*all_exner(i,j,k)
        p_rho(i,j,k) = (exner_rho(i,j,k)**(1.0/kappa))*pref
      ELSE
        p_rho(i,j,k) = all_pstar(i,j)  ! surface  
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO


!-----------------------------------------------------------------------
! Accurate calculation of height of LCL using exner_lcl rather than p_lcl
! as UM interpolates exner linearly in height but not pressure.
!-----------------------------------------------------------------------
DO j=1,local_new_y
  DO i=1,local_new_x

    k = 1

    IF ( exner_lcl(i,j) >= exner_surf(i,j)) THEN
      zlcl(i,j) = 0.0           ! at or below surface
    ELSE IF (exner_lcl(i,j) < exner_surf(i,j)                                 &
                   .AND. exner_lcl(i,j) > exner_rho(i,j,k)) THEN
      factor= (exner_rho(i,j,k) - exner_lcl(i,j))/                            &
                      (exner_rho(i,j,k) - exner_surf(i,j))
      zlcl(i,j) = (1.0-factor)*h_rho_sea(k)
    END IF

    DO k=2,mlevs
      IF (exner_lcl(i,j) >= exner_rho(i,j,k)                                  &
                        .AND. exner_lcl(i,j) < exner_rho(i,j,k-1) ) THEN
        factor= (exner_rho(i,j,k) - exner_lcl(i,j))/                          &
                        (exner_rho(i,j,k) - exner_rho(i,j,k-1))
        zlcl(i,j) = (1.0-factor)*h_rho_sea(k)+factor*h_rho_sea(k-1)
      END IF
    END DO         ! level loop

  END DO
END DO

IF (iprint >1) THEN
  WRITE(umMessage,'(A)') '  orog        '
  CALL umPrint(umMessage,src=RoutineName)
  DO j = 1,local_new_y
    WRITE(umMessage,'(20F6.0)') (all_orog(i,j),i=1,local_new_x)
    CALL umPrint(umMessage,src=RoutineName)
  END DO

  WRITE(umMessage,'(A)') '  zlcl        '
  CALL umPrint(umMessage,src=RoutineName)
  DO j = 1,local_new_y
    WRITE(umMessage,'(20F6.0)') (zlcl(i,j),i=1,local_new_x)
    CALL umPrint(umMessage,src=RoutineName)
  END DO

  WRITE(umMessage,'(A)') '  tlcl        '
  CALL umPrint(umMessage,src=RoutineName)
  DO j = 1,local_new_y
    WRITE(umMessage,'(20F6.0)') (t_lcl(i,j),i=1,local_new_x)
    CALL umPrint(umMessage,src=RoutineName)
  END DO

  WRITE(umMessage,'(A)') '  plcl        '
  CALL umPrint(umMessage,src=RoutineName)
  DO j = 1,local_new_y
    WRITE(umMessage,'(20F8.0)') (p_lcl(i,j),i=1,local_new_x)
    CALL umPrint(umMessage,src=RoutineName)
  END DO
END IF
!-----------------------------------------------------------------------
! Section to do a parcel ascent from near surface to calculate CAPE & CIN
!-----------------------------------------------------------------------

IF (l_cape) THEN

  n_itime = 0
  l_itime = .TRUE.
  qlcrit = 1.0e-3 ! kg/kg
  qlmin  = 2.0e-4 ! kg/kg
  deltat   = 0.2 ! K iteration for temperature interval
  rdeltat  = 1.0/0.2 ! K iteration for temperature interval
  smallcape= 1.0 ! J/kg a small value for cape (somewhat arbitrary)
  smallpressure = 1.0 ! Pa

  ! Initiliase arrays
  DO j=1,local_new_y
    DO i=1,local_new_x
      pnb(i,j) = 0.0
      cin(i,j) = 0.0
      cape(i,j) = 0.0
      qvparcel(i,j) = 0.0
      ktop(i,j) = 0
      iitime(i,j) = 0 ! initialise counter
      qlval(i,j) = 0.0
      zneutral(i,j) = 0.0

      ! LCL level relative to surface
      zlcl_surf = zlcl(i,j) - all_orog(i,j)
      IF (zlcl_surf < 0.0) THEN
         zlcl_surf = 5.0    ! lowest possible layer depth
      END IF  
      ! Deep turbulence ascent starts from zlcl/2.
      zstart(i,j) = zlcl_surf*0.5 + all_orog(i,j)
      kstart(i,j) = 1   ! Initialise start level for ascent
      zfree(i,j) = 0.0
    END DO
  END DO

  !-----------------------------------------------------------------------
  ! Starting level for ascent
  !-----------------------------------------------------------------------

  DO k=2, mlevs-1
    DO j=1,local_new_y
      DO i=1,local_new_x
        IF (zstart(i,j) >= h_theta_sea(k) .AND.          &
              zstart(i,j) < h_theta_sea(k+1)) THEN
          kstart(i,j) = k
        END IF
      END DO
    END DO
  END DO

  kmin = mlevs-1      ! lowest start level
  kmax = mlevs-1    ! maximum level for ascent
  DO j=1,local_new_y
    DO i=1,local_new_x
      IF (kstart(i,j) < kmin ) THEN
        kmin = kstart(i,j)
      END IF
    END DO
  END DO

  IF (iprint >1) THEN
    WRITE(umMessage,'(A)') '  kstart        '
    CALL umPrint(umMessage,src=RoutineName)

    DO j = 1,local_new_y
      WRITE(umMessage,'(20i10)') (kstart(i,j),i=1,local_new_x)
      CALL umPrint(umMessage,src=RoutineName)
    END DO

    WRITE(umMessage,'(A,2I6)') ' kmin, kmax  ',kmin,kmax
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  !----------------------------------------------------------------------
  ! Parcel ascent from   zlcl_surf/2  or bottom level ?
  !----------------------------------------------------------------------

  ! Use method from deep turbulence scheme for doing ascent
  DO j=1,local_new_y
    DO i=1,local_new_x
      thparcel(i,j) = all_th(i,j,kstart(i,j))
      tparcel(i,j)  = thparcel(i,j) *all_exner(i,j,kstart(i,j))
      ! use large-scale value
      qlbot(i,j)    = all_qcl(i,j,kstart(i,j))+ all_qcf(i,j,kstart(i,j))
      p(i,j) = all_ptheta(i,j,kstart(i,j))
    END DO
  END DO

  ! Find the saturation value of this initial parcel
  IF (l_mr_physics) THEN
    CALL qsat_mix(q0,tparcel,p,local_new_x,local_new_y)
  ELSE
    CALL qsat(q0,tparcel,p,local_new_x,local_new_y)
  END IF

  DO j=1,local_new_y
    DO i=1,local_new_x
      IF (h_theta_sea(k_plume(i,j)) <= zlcl(i,j)) THEN
        ! if the parcel starts beneath the lcl, then give it the water vapour
        ! value for that level multiplying by an extra factor of 1.2
        qvbot(i,j) = all_q(i,j,k_plume(i,j))*1.2
      ELSE
        ! if the parcel starts above the lcl, then say that it's already
        ! saturated and give it the plume liquid value
        qvbot(i,j) = q0(i,j)
        qlbot(i,j) = 0.5*q0(i,j)   ! approx qcl_plume
        IF (qlbot(i,j) > qlcrit) THEN
          qlbot(i,j) = qlcrit   
        END IF
      END IF
    END DO
  END DO

  !==================================================================
  ! Start parcel ascent -- increasing height by vertical levels
  !==================================================================
  ! loop dependent on k-1 values so cannot use openMP on K
  DO k = kmin,kmax
    DO j=1,local_new_y
      DO i=1,local_new_x
        IF (all_a(i,j,k) > 0.0) THEN  ! safe should be values 
          p(i,j) = all_ptheta(i,j,k)
          ttemp(i,j) = all_t(i,j,k)
          delp_cen(i,j) = all_ptheta(i,j,k)- all_ptheta(i,j,k+1)
          delp(i,j)     =( p_rho(i,j,k) - p_rho(i,j,k+1)) *all_a(i,j,k)
        ELSE  ! All points in mean below model surface
          p(i,j) = all_pstar(i,j)   ! setting to suface pressure for qsat
          ttemp(i,j) = all_tstar(i,j) ! setting to suface temp for qsat
          delp(i,j) = 0.0
          delp_cen(i,j) = 0.0
        END IF  
      END DO
    END DO

    ! calculate how much water vapour would be needed at this T and p in
    ! order to be saturated
    ! for later moist adiabatic ascent, will need gradient as a function of T
    IF (l_mr_physics) THEN
      CALL qsat_mix(qse,ttemp,ptemp,local_new_x,local_new_y)
      CALL qsat_mix(q0,tparcel,p,local_new_x,local_new_y)
    ELSE
      CALL qsat(qse,ttemp,ptemp,local_new_x,local_new_y)
      CALL qsat(q0,tparcel,p,local_new_x,local_new_y)
    END IF


    DO j=1,local_new_y
      DO i=1,local_new_x
        qvparcel(i,j) = q0(i,j) ! first suggest saturated

        ! Plume cloud water  
        qcl_plume(i,j) = 0.5*qse(i,j)
        IF (qcl_plume(i,j) > qlcrit) THEN
          qcl_plume(i,j) = qlcrit
        ELSE IF (qcl_plume(i,j) < qlmin) THEN
          qcl_plume(i,j) = qlmin
        END IF

        IF (k >= kstart(i,j) .AND. pnb(i,j) < smallpressure) THEN

          ! liquid water in parcel is total water minus saturation value
          qlval(i,j) = qlbot(i,j)+qvbot(i,j)-qvparcel(i,j)

          IF (qlval(i,j) < 0.0) THEN
            qlval(i,j) = 0.0 ! ensure >= 0
          ELSE IF (qlval(i,j) > qcl_plume(i,j)) THEN
            qlval(i,j) = qcl_plume(i,j)
          END IF
          ! If the parcel was unsaturated, then its value stays that of the 
          ! initial parcel
          IF (qlval(i,j) < 1.0e-10) THEN
            qvparcel(i,j)  = qvbot(i,j)
          END IF

          ! evaluation of parcel virt.pot.temp keeping all water loading
          thvparcel = thparcel(i,j)*(1.0+c_virtual*qvparcel(i,j)-qlval(i,j))

          thetav = all_th(i,j,k)*(1.0+c_virtual*all_q(i,j,k))
          ! evaluation of parcel virt.pot.temp losing all liquid water
          !thvparcel = thparcel(i,j)*(1.0+c_virtual*qvparcel(i,j))

          ! contribution to cape from this level:

          gdz_o_thv =  delp(i,j)/(thetav*all_rho(i,j,k))
          capecom = gdz_o_thv*(thvparcel-thetav)

          IF (capecom > 0.0) THEN
            cape(i,j)=cape(i,j)+capecom 
          END IF

          ! if negative and cape is still small, then count as cin
          IF (capecom < 0.0 .AND. cape(i,j) <= smallcape) THEN
            cin(i,j)=cin(i,j)+capecom
          END IF

          ! Level of free convection 
          IF (zfree(i,j) == 0.0 .AND. h_theta_sea(k) >= zlcl(i,j) .AND. &
             capecom > 0.0 ) THEN
            zfree(i,j) = h_theta_sea(k)
          END IF

          ! define level of neutral buoyancy as the first level at which
          ! negatively buoyant once cape is greater than a certain value

          IF (capecom < 0.0 .AND. cape(i,j) > smallcape) THEN
            pnb(i,j)=p(i,j)
            zneutral(i,j) = h_theta_sea(k)
            ktop(i,j) = k
          END IF

        END IF !k >= kstart(i,j) .and. pnb(i,j) == 0.0
      END DO 
    END DO 

    ! Determine qsat for the level below starting point
    ! So only call if still points yet to stat ascent.
    IF (l_itime) THEN
      DO j=1,local_new_y
        DO i=1,local_new_x
          IF (k == 1) THEN
            ptemp(i,j) = all_pstar(i,j)
          ELSE
            ptemp(i,j) = all_ptheta(i,j,k-1)
          END IF
        END DO 
      END DO 
    END IF

    DO j=1,local_new_y
      DO i=1,local_new_x
        ptemp(i,j) = all_ptheta(i,j,k+1)
      END DO 
    END DO 

    ! break out of points loop to dO qsat calculations
    IF (l_mr_physics) THEN
      CALL qsat_mix(qnext,tparcel,ptemp,local_new_x,local_new_y)
      CALL qsat_mix(q1,tparcel-deltat,ptemp,local_new_x,local_new_y)
    ELSE
      CALL qsat(qnext,tparcel,ptemp,local_new_x,local_new_y)
      CALL qsat(q1,tparcel-deltat,ptemp,local_new_x,local_new_y)
    END IF

    DO j=1,local_new_y
      DO i=1,local_new_x
        IF (k >= kstart(i,j) .AND. pnb(i,j) == 0.0) THEN

          ! determine a local gradient for the qsaturation
          ! (assume a local linear relationship between T and q for a given p)

          IF (qlval(i,j) > 0.0) THEN ! if it's saturated then do moist ascent
            t0 = tparcel(i,j)
            !        t1 = t0 - deltat   ! not required as use (t0-t1)=deltat

            ! choose between L_vap and L_fus tm = 273.15 K
            IF (tparcel(i,j) >= tm) THEN
              vlt = lc
            ELSE
              vlt = lc+lf
            END IF
            ! determine gradient of q vs T line
            alpha = (qnext(i,j)-q1(i,j))*rdeltat
            beta = q0(i,j)-alpha*t0

            dp_o_p = delp_cen(i,j)/p(i,j)

            ! Calculation of bb1 including 2nd order term
            bb1 =  -(r*tparcel(i,j)*dp_o_p)*(1.0-(kappa-1.0)*dp_o_p*0.5)

            !        bb2 = 1.0 + q0(i,j)
            rbb2 = 1.0/(1.0 + q0(i,j))
            drs0 = qnext(i,j)-q0(i,j)

            !        dt = (bb1-vlt*drs0/bb2)/(cp+alpha*vlt/bb2)
            dt = (bb1-vlt*drs0*rbb2)/(cp+alpha*vlt*rbb2)

            dq = drs0 + alpha*dt

            tparcel(i,j) = tparcel(i,j) + dt
            iitime(i,j) = iitime(i,j) + 1
            n_itime = n_itime + 1

            ! this is the potential temperature at the next height up
            thparcel(i,j) = tparcel(i,j)/all_exner(i,j,k+1)

          ELSE

            !deduce T given new pressure
            tparcel(i,j) = thparcel(i,j)*all_exner(i,j,k+1)

          END IF

        END IF ! (k >= kstart(i,j) .and. pnb(i,j) == 0.0)
      END DO
    END DO

    IF (n_itime == npnts) THEN  ! All ascents reached pnb
      l_itime = .FALSE.
    END IF

  END DO ! k

  ! Need to set missing data for ascent where unable to calculate
  ! CAPE and CIN 

  DO j=1,local_new_y
    DO i=1,local_new_x
      IF (zneutral(i,j) == h_theta_sea(kmax) ) THEN
        ! Gone to top of model - a problem
        zneutral(i,j) = rmdi
        zfree(i,j)    = rmdi
        cin(i,j)      = rmdi
        cape(i,j)     = rmdi
      ELSE IF ( zneutral(i,j) < zlcl(i,j)) THEN 
        ! No CAPE and undilute parcel top  
        zneutral(i,j) = rmdi
        zfree(i,j)    = rmdi
        cape(i,j)     = rmdi
        cin(i,j)      = rmdi
      ELSE IF ( zfree(i,j) > 6000.) THEN 
        ! Do I want some check on whether zfree is too high?
        ! Say 6000m or more
        zneutral(i,j) = rmdi
        zfree(i,j)    = rmdi
        cape(i,j)     = rmdi
        cin(i,j)      = rmdi
      END IF 
    END DO 
  END DO 

END IF   ! l_cape
!==============================================================
! End of parcel ascent
!==============================================================
! Restrict printout for checking only
IF (iprint == 2) THEN
  WRITE(umMessage,'(A)') ' CAPE         '
  CALL umPrint(umMessage,src=RoutineName)
  DO j = 1,local_new_y
    WRITE(umMessage,'(20F10.0)') (cape(i,j),i=1,local_new_x)
    CALL umPrint(umMessage,src=RoutineName)
  END DO

  WRITE(umMessage,'(A)') '  CIN        '
  CALL umPrint(umMessage,src=RoutineName)
  DO j = 1,local_new_y
    WRITE(umMessage,'(20F10.3)') (cin(i,j),i=1,local_new_x)
    CALL umPrint(umMessage,src=RoutineName)
  END DO

  WRITE(umMessage,'(A)') '  zfree        '
  CALL umPrint(umMessage,src=RoutineName)
  DO j = 1,local_new_y
    WRITE(umMessage,'(20F7.0)') (zfree(i,j),i=1,local_new_x)
    CALL umPrint(umMessage,src=RoutineName)
  END DO

  WRITE(umMessage,'(A)') '  zneutral        '
  CALL umPrint(umMessage,src=RoutineName)
  DO j = 1,local_new_y
    WRITE(umMessage,'(20F7.0)') (zneutral(i,j),i=1,local_new_x)
    CALL umPrint(umMessage,src=RoutineName)
  END DO
END IF

!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!----------------------------------------------------------------------
RETURN
END SUBROUTINE cape_cin_from_mean

END MODULE cape_cin_from_mean_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate the scaled divergence of the Eliassen-Palm flux from
! u, v, w and temperature on TEM pressure levels.
! Residual mean meridional circulation, Meridional and Vertical
! components of Eliassen-Palm, Meridional heat flux and momentum
! fluxes are returned also.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Climate Diagnostics
MODULE calc_div_ep_flux_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_DIV_EP_FLUX_MOD'

CONTAINS


SUBROUTINE Calc_Div_EP_Flux(                                            &
    rows, n_rows, row_length, global_row_length, sin_v_latitude,        &
    qvstarbar_p, qwstarbar_p, qepy_p, qepz_p, qdivep_p,                 &
    qzvptp_p, qzupvp_p, tem_press, n_tem_levs,                          &
    u_tem_p, v_tem_p, w_tem_p, temp_tem_p,                              &
    vstarbar, wstarbar, Fy, Fz, divF, zvptp, zupvp)

USE planet_constants_mod, ONLY: planet_radius, omega, &
                                kappa, p_zero, sclht

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE atm_fields_bounds_mod, ONLY:            &
                 udims, udims_s, vdims, vdims_s

USE UM_ParVars
USE UM_ParParams,      ONLY: pnorth, psouth
USE Field_Types,       ONLY: fld_type_v
USE zonal_average_mod, ONLY: zonal_average
USE mpp_conf_mod,      ONLY: swap_field_is_vector, swap_field_is_scalar

USE ereport_mod,             ONLY: ereport
USE errormessagelength_mod,  ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: rows, n_rows, row_length, global_row_length

! No of levels for input u_tem_p, v_tem_p, w_tem_p and temp_tem_p
!          and for ouputs
INTEGER, INTENT(IN) :: n_tem_levs

! Surface density at reference latitude
REAL, PARAMETER  :: rho_surface = 1.212

REAL :: sin_v_latitude(vdims%i_start:vdims%i_end,                      &
                    vdims%j_start:vdims%j_end)

REAL :: zlogc(n_tem_levs)  ! "log-pressure" coordinate

REAL, INTENT(IN) ::                                                 &
     u_tem_p(udims%i_start:udims%i_end,                                 &
         vdims%j_start:vdims%j_end,n_tem_levs)                        &
                                       ! u at selected pressures
    ,v_tem_p(udims%i_start:udims%i_end,                                 &
         vdims%j_start:vdims%j_end,n_tem_levs)                        &
                                       ! v at selected pressures
    ,w_tem_p(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,       &
         n_tem_levs)
                                       ! w at selected pressures

REAL, INTENT(IN) :: tem_press(n_tem_levs)

LOGICAL, INTENT(IN) ::                                              &
     qvstarbar_p                                                    &
             ! Flag for residual mean meridional circulation
    ,qwstarbar_p                                                    &
             ! Flag for residual mean meridional circulation
    ,qepy_p                                                         &
             ! Flag for Eliassen Palm flux (phi component)
    ,qepz_p                                                         &
             ! Flag for Eliassen Palm flux (vertical component)
    ,qdivep_p                                                       &
             ! Flag for divergence of Eliassen Palm flux
    ,qzvptp_p                                                       &
             ! Flag for meridional heat flux
    ,qzupvp_p                                                       
             ! Flag for meridional momentum flux

INTEGER :: start_j, end_j

REAL ::  uv(udims%i_start:udims%i_end,                              &
            vdims%j_start:vdims%j_end, n_tem_levs)                    &
       , theta_p(udims%i_start:udims%i_end,                         &
                 vdims%j_start:vdims%j_end, n_tem_levs)               &
       , vT(udims%i_start:udims%i_end,                              &
            vdims%j_start:vdims%j_end, n_tem_levs)                    &
       , vtheta(udims%i_start:udims%i_end,                          &
                vdims%j_start:vdims%j_end, n_tem_levs)                &
       , temp_tem_p(udims%i_start:udims%i_end,                     &
                     vdims%j_start:vdims%j_end,  n_tem_levs)          &
       , uw(udims%i_start:udims%i_end,                              &
            vdims%j_start:vdims%j_end, n_tem_levs)                    &
       , zonal_u_halo(vdims_s%j_start:vdims_s%j_end, n_tem_levs)      &
       , zonal_u(vdims%j_start:vdims%j_end, n_tem_levs)               &
       , zonal_v(vdims%j_start:vdims%j_end, n_tem_levs)               &
       , zonal_uv(vdims%j_start:vdims%j_end, n_tem_levs)              &
       , zonal_theta(vdims%j_start:vdims%j_end, n_tem_levs)           &
       , zonal_t(vdims%j_start:vdims%j_end, n_tem_levs)               &
       , zonal_vT(vdims%j_start:vdims%j_end, n_tem_levs)              &
       , zonal_vtheta(vdims%j_start:vdims%j_end, n_tem_levs)          &
       , zonal_w(vdims%j_start:vdims%j_end, n_tem_levs)               &
       , zonal_uw(vdims%j_start:vdims%j_end, n_tem_levs)              &
       , zvpthp(vdims_s%j_start:vdims_s%j_end, n_tem_levs)            &
       , zupwp(vdims%j_start:vdims%j_end, n_tem_levs)                 &
       , dthdz(vdims_s%j_start:vdims_s%j_end, n_tem_levs)             &
       , phi(vdims_s%j_start:vdims_s%j_end)                         &
       , rho0(n_tem_levs)                                             &
       , dudz(vdims%j_start:vdims%j_end, n_tem_levs)                  &
       , ducosdphi(vdims%j_start:vdims%j_end, n_tem_levs)             &
       , divFz(vdims%j_start:vdims%j_end, n_tem_levs)                 &
       , divFy(vdims%j_start:vdims%j_end, n_tem_levs)                 &
       , Fy_halo(vdims_s%j_start:vdims_s%j_end, n_tem_levs)

REAL, INTENT(OUT) ::                                                &
     vstarbar(vdims%j_start:vdims%j_end, n_tem_levs)                  &
    ,wstarbar(vdims%j_start:vdims%j_end, n_tem_levs)                  &
    ,divF(vdims%j_start:vdims%j_end, n_tem_levs)                      &
    ,Fy(vdims%j_start:vdims%j_end, n_tem_levs)                        &
    ,Fz(vdims%j_start:vdims%j_end, n_tem_levs)                        &
    ,zvptp(vdims%j_start:vdims%j_end, n_tem_levs)                     &
    ,zupvp(vdims%j_start:vdims%j_end, n_tem_levs)

INTEGER :: x,y,z
INTEGER :: top

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
REAL :: work_zupvp(vdims%j_start:vdims%j_end, n_tem_levs)

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_DIV_EP_FLUX'

! ----------------------------------------------------------------------
! STASH items 310 311: residual circulation on pressure surfaces/B grid
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (qvstarbar_p .OR. qwstarbar_p .OR. qepy_p .OR. qepz_p &
     .OR. qdivep_p .OR. qzvptp_p .OR. qzupvp_p) THEN
  !
  ! Calculate the "log-pressure" coordinate z=-Hln(p/ps) and basic density
  ! from surface density at reference latitude.
  !
  DO z=1,n_tem_levs
    zlogc(z) = - sclht * LOG(tem_press(z)/(p_zero*0.01))
    rho0(z)      = rho_surface * tem_press(z) / (p_zero*0.01)
  END DO

  DO z=1,n_tem_levs
    DO y= vdims%j_start,vdims%j_end      ! latitude loop
      DO x=udims%i_start,udims%i_end     ! longitude loop
        theta_p(x,y,z) = temp_tem_p(x,y,z) *                     &
                            ((p_zero*0.01) / tem_press(z))**kappa
      END DO
    END DO
  END DO

  DO z=1,n_tem_levs
    DO y= vdims%j_start,vdims%j_end      ! latitude loop
      DO x= udims%i_start,udims%i_end    ! longitude loop
        uv(x,y,z)     = u_tem_p(x,y,z) * v_tem_p(x,y,z)
        vtheta(x,y,z) = v_tem_p(x,y,z) * theta_p(x,y,z)
        vT(x,y,z)     = v_tem_p(x,y,z) * temp_tem_p(x,y,z)
        uw(x,y,z)     = u_tem_p(x,y,z) * w_tem_p(x,y,z)
      END DO
    END DO
  END DO

  CALL zonal_average(w_tem_p, zonal_w,udims%i_start,udims%i_end,      &
                                  vdims%j_start,vdims%j_end,      &
                                  global_row_length,              &
                                  n_tem_levs, gc_proc_row_group)

  CALL zonal_average(u_tem_p, zonal_u,udims%i_start,udims%i_end,      &
                                  vdims%j_start,vdims%j_end,      &
                                  global_row_length, &
                                  n_tem_levs, gc_proc_row_group)

  CALL zonal_average(v_tem_p, zonal_v,udims%i_start,udims%i_end,      &
                                  vdims%j_start,vdims%j_end,      &
                                  global_row_length, &
                                  n_tem_levs, gc_proc_row_group)

  CALL zonal_average(uv, zonal_uv,udims%i_start,udims%i_end,      &
                                  vdims%j_start,vdims%j_end,      &
                                  global_row_length, &
                                  n_tem_levs, gc_proc_row_group)

  CALL zonal_average(theta_p, zonal_theta,                        &
                                  udims%i_start,udims%i_end,      &
                                  vdims%j_start,vdims%j_end,      &
                                  global_row_length, &
                                  n_tem_levs, gc_proc_row_group)

  CALL zonal_average(temp_tem_p, zonal_t,                        &
                                  udims%i_start,udims%i_end,      &
                                  vdims%j_start,vdims%j_end,      &
                                  global_row_length, &
                                  n_tem_levs, gc_proc_row_group)

  CALL zonal_average(vtheta, zonal_vtheta,                        &
                                  udims%i_start,udims%i_end,      &
                                  vdims%j_start,vdims%j_end,      &
                                  global_row_length, &
                                  n_tem_levs, gc_proc_row_group)

  CALL zonal_average(uw, zonal_uw,udims%i_start,udims%i_end,      &
                                  vdims%j_start,vdims%j_end,      &
                                  global_row_length, &
                                  n_tem_levs, gc_proc_row_group)

  CALL zonal_average(vT, zonal_vT,udims%i_start,udims%i_end,      &
                                  vdims%j_start,vdims%j_end,      &
                                  global_row_length, &
                                  n_tem_levs, gc_proc_row_group)

  ! Now work out the zonal average of the prime quantities needed
  ! for EP flux. eg {u'v'}={uv}-{u}{v}, where {}=zonal average.
  ! Calculate the meridional heat (v'T') and momentum (u'v') fluxes
  !
  DO z=1,n_tem_levs
    DO y=vdims%j_start,vdims%j_end
      work_zupvp(y,z) =zonal_uv(y,z) - zonal_u(y,z) * zonal_v(y,z)
      zvpthp(y,z)=zonal_vtheta(y,z) - zonal_v(y,z) * zonal_theta(y,z)
      zupwp(y,z) =zonal_uw(y,z) - zonal_u(y,z) * zonal_w(y,z)
    END DO
  END DO

  ! Copy values to meridional momentum flux (STASH item 316) when requested
  IF (qzupvp_p) THEN
    DO z=1,n_tem_levs
      DO y=vdims%j_start,vdims%j_end
        zupvp(y,z) = work_zupvp(y,z)
      END DO
    END DO
  END IF

  ! Calculate meridional heat flux (STASH item 315) when requested
  IF (qzvptp_p) THEN
    DO z=1,n_tem_levs
      DO y=vdims%j_start,vdims%j_end
        zvptp(y,z)=zonal_vT(y,z) - zonal_v(y,z) * zonal_T(y,z)
      END DO
    END DO
  END IF

  !
  !  Calculate dtheta/dz
  !
  ! For lowest level
  ! 
  DO y=vdims%j_start,vdims%j_end
    dthdz(y,1) = (zonal_theta(y,2) -                           &
                  zonal_theta(y,1)) /  (zlogc(2) - zlogc(1))
  END DO

  ! For highest level
  !
  DO y=vdims%j_start,vdims%j_end
    dthdz(y,n_tem_levs) = (zonal_theta(y,n_tem_levs) -             &
                           zonal_theta(y,n_tem_levs-1)) /          &
                          (zlogc(n_tem_levs) - zlogc(n_tem_levs-1))
  END DO

  DO z=2,n_tem_levs-1
    DO y=vdims%j_start,vdims%j_end
      dthdz(y,z) = (zonal_theta(y,z+1) - zonal_theta(y,z-1)) / &
                         (zlogc(z+1)-zlogc(z-1))
    END DO
  END DO
END IF

IF (qvstarbar_p) THEN
  !
  !  Calculate vstarbar. STASH item 310
  ! ----------------------------------------------------------------------
  !  v* = v_bar - (d/dz [rho0] / rho0 * v'th'/dthdz) - d/dz [ v'th'/dthdz ]
  !
  top = n_tem_levs
  ! For bottom level
  DO y=vdims%j_start,vdims%j_end
    vstarbar(y,1) = zonal_v(y,1) +                         &
                 ((zvpthp(y,1) / dthdz(y,1)) / sclht) -    &
                 (                                         &
                   ((zvpthp(y,2) / dthdz(y,2)) -           &
                    (zvpthp(y,1) / dthdz(y,1)) )/          &
                   (zlogc(2) - zlogc(1))                   &
                 )
  END DO

  ! For top level
  DO y=vdims%j_start,vdims%j_end
    vstarbar(y,top) = zonal_v(y,top) +                     &
              ((zvpthp(y,top) / dthdz(y,top)) / sclht) -   &
              (                                            &
                  ( (zvpthp(y,top) / dthdz(y,top) )-       &
                   ( zvpthp(y,top-1) / dthdz(y,top-1))) /  &
                   (zlogc(top) -zlogc(top-1))              &
              )
  END DO

  DO z=2,top-1
    DO y=vdims%j_start,vdims%j_end
      vstarbar(y,z) = zonal_v(y,z) +                         &
                      ((zvpthp(y,z) / dthdz(y,z)) / sclht) - &
                      (                                      &
                      (  (zvpthp(y,z+1) / dthdz(y,z+1) )-    &
                       (  zvpthp(y,z-1) / dthdz(y,z-1)))/    &
                        (zlogc(z+1) - zlogc(z-1))            &
                      )
    END DO
  END DO
END IF ! on STASHflag

IF (qwstarbar_p .OR. qepy_p .OR. qepz_p .OR. qdivep_p ) THEN
  !
  ! Obtain latitude needed in calculation of wstarbar and Eliassen-Palm flux.
  !
  DO y=vdims%j_start,vdims%j_end
    phi(y) = ASIN(sin_v_latitude(1,y))
  END DO
END IF

IF (qwstarbar_p) THEN
  !
  !  Calculate wstarbar. STASH item 311
  ! ----------------------------------------------------------------------
  !  w* = w_bar - 1/(a cos{phi}) * d/dphi [ cos{phi} * v'th'/dthdz ]
  !
  ! DEPENDS ON: swap_bounds
  CALL Swap_bounds(zvpthp,1,n_rows,n_tem_levs,                    &
                   0,offy,fld_type_v,swap_field_is_vector)
  ! DEPENDS ON: swap_bounds
  CALL Swap_bounds(dthdz,1,n_rows,n_tem_levs,                     &
                   0,offy,fld_type_v, swap_field_is_scalar)
  ! DEPENDS ON: swap_bounds
  CALL Swap_bounds(phi,1,n_rows,1,                              &
                   0,offy,fld_type_v, swap_field_is_scalar)

  DO z=n_tem_levs,1,-1
    start_j = vdims%j_start
    end_j = vdims%j_end

    IF (at_extremity(PSouth)) THEN
      start_j = vdims%j_start+1
      wstarbar(vdims%j_start,z) = zonal_w(vdims%j_start,z) +    &
             (1.0 / (planet_radius * COS(phi(vdims%j_start)))) *&
               ( -SIN(phi(vdims%j_start))                       &
      * zvpthp(vdims%j_start,z) / dthdz(vdims%j_start,z) +      &
                 COS(phi(vdims%j_start)) *                      &
                 (                                              &
     (zvpthp(vdims%j_start+1,z) / dthdz(vdims%j_start+1,z)) -   &
     (zvpthp(vdims%j_start,z) / dthdz(vdims%j_start,z))         &
                 )                                              &
         / (phi(vdims%j_start+1) - phi(vdims%j_start))          &
               )
    END IF

    IF (at_extremity(PNorth)) THEN
      end_j = vdims%j_end-1
      wstarbar(vdims%j_end,z) = zonal_w(vdims%j_end,z) +        &
               (1.0 / (planet_radius * COS(phi(vdims%j_end)))) * &
               (                                                &
       -SIN(phi(vdims%j_end)) * zvpthp(vdims%j_end,z)           &
                                  / dthdz(vdims%j_end,z) +      &
                 COS(phi(vdims%j_end)) *                        &
                 (                                              &
         (zvpthp(vdims%j_end,z) / dthdz(vdims%j_end,z)) -       &
         (zvpthp(vdims%j_end-1,z) / dthdz(vdims%j_end-1,z))     &
                 )                                              &
        / (phi(vdims%j_end) - phi(vdims%j_end-1))               &
               )
    END IF

    !CDIR NOUNROLL
    DO y=start_j,end_j
      wstarbar(y,z) = zonal_w(y,z) +                            &
               (1.0 / (planet_radius * COS(phi(y)))) *          &
               (                                                &
                 -SIN(phi(y)) * zvpthp(y,z) / dthdz(y,z) +      &
                 COS(phi(y)) *                                  &
                 (                                              &
                   (zvpthp(y+1,z) / dthdz(y+1,z)) -             &
                   (zvpthp(y-1,z) / dthdz(y-1,z))               &
                 )                                              &
                  / (phi(y+1) - phi(y-1))                       &
               )
    END DO
  END DO
END IF

IF (qepy_p .OR. qdivep_p) THEN
  !
  !  Calculate dubar/dz if required
  !
  ! bottom level
  DO y=vdims%j_start,vdims%j_end
    dudz(y,1) = (zonal_u(y,2) - zonal_u(y,1)) /                 &
                       (zlogc(2) - zlogc(1))
  END DO

  ! top level
  DO y=vdims%j_start,vdims%j_end
    dudz(y,n_tem_levs) = (zonal_u(y,n_tem_levs) -                   &
                        zonal_u(y,n_tem_levs-1)) /                &
                        (zlogc(n_tem_levs) - zlogc(n_tem_levs-1))
  END DO

  DO z=2,n_tem_levs-1
    DO  y=vdims%j_start,vdims%j_end
      dudz(y,z) = (zonal_u(y,z+1) - zonal_u(y,z-1)) /           &
                          (zlogc(z+1) - zlogc(z-1))
    END DO
  END DO
END IF

IF (qepz_p .OR. qdivep_p) THEN
  !
  !  Calculate d(ubar*cos(phi))/dphi * (1./radius*cos(phi))
  !
  DO z=1,n_tem_levs
    DO y=vdims%j_start,vdims%j_end
      zonal_u_halo(y,z) = zonal_u(y,z)
    END DO
  END DO

  ! DEPENDS ON: swap_bounds
  CALL Swap_bounds(zonal_u_halo,1,n_rows,n_tem_levs,            &
                   0,offy,fld_type_v,swap_field_is_vector)

  DO z=1,n_tem_levs
    start_j = vdims%j_start
    end_j = vdims%j_end
    IF (at_extremity(PSouth)) THEN
      start_j = vdims%j_start+1
      ducosdphi(vdims%j_start,z) = (                               &
(zonal_u_halo(vdims%j_start,z) - zonal_u_halo(vdims%j_start+1,z)) /  &
(phi(vdims%j_start) - phi(vdims%j_start+1)) -                       &
 zonal_u_halo(vdims%j_start,z) * TAN(phi(vdims%j_start))           &
                       ) / planet_radius
    END IF

    IF (at_extremity(PNorth)) THEN
      end_j = vdims%j_end-1
      ducosdphi(vdims%j_end,z) = (                                 &
  (zonal_u_halo(vdims%j_end-1,z) - zonal_u_halo(vdims%j_end,z)) /  &
   (phi(vdims%j_end-1) - phi(vdims%j_end)) -                       &
    zonal_u_halo(vdims%j_end,z) * TAN(phi(vdims%j_end))            &
                            ) / planet_radius
    END IF
    !CDIR NOUNROLL
    DO y=start_j,end_j
      ducosdphi(y,z) = (                                              &
                        (zonal_u_halo(y-1,z) - zonal_u_halo(y+1,z)) / &
                         (phi(y-1) - phi(y+1)) -                      &
                         zonal_u_halo(y,z) * TAN(phi(y))              &
                       ) / planet_radius
    END DO
  END DO
END IF
! ----------------------------------------------------------------------
! STASH items 312 313: Eliassen-Palm flux (phi & z components)
!                      on pressure surfaces/B grid
! ----------------------------------------------------------------------
!
! Calculate Eliassen-Palm flux
! ----------------------------------------------------------------------
!  F_phi = rho0 * a cos{phi} * ( dubardz * v'th'/dthdz - u'v' )
!  F_z   = rho0 * a cos{phi} * ( (2 Omega sin{phi} - 1/(a cos{phi}) *
!               d/dphi [ubar cos{phi}] ) * v'th'/dthdz - u'w' )
!
IF (qepy_p .OR. qdivep_p) THEN
  DO z=1,n_tem_levs
    DO y=vdims%j_start,vdims%j_end
      Fy(y,z) = rho0(z) * planet_radius * COS(phi(y)) *      &
                 (                                           &
                   dudz(y,z) * zvpthp(y,z) / dthdz(y,z)      &
                   - work_zupvp(y,z)                         &
                 )
    END DO
  END DO
END IF

IF (qepz_p .OR. qdivep_p) THEN
  DO z=1,n_tem_levs
    DO y=vdims%j_start,vdims%j_end
      Fz(y,z) = rho0(z) * planet_radius * COS(phi(y)) *      &
                 (                                           &
                   (2.0 * omega * SIN(phi(y)) -               &
                    ducosdphi(y,z)) *                        &
                    (zvpthp(y,z) / dthdz(y,z)) - zupwp(y,z)  &
                 )
    END DO
  END DO
END IF

IF (qdivep_p) THEN
  !
  ! Calculate Div F. STASH item 314
  ! ----------------------------------------------------------------------
  ! 1/(a cos{phi}) * d/dphi [ cos{phi} F_phi ] + d/dz [F_z]
  !
  DO z=1,n_tem_levs
    DO y=vdims%j_start,vdims%j_end
      Fy_halo(y,z) = Fy(y,z)
    END DO
  END DO

  ! DEPENDS ON: swap_bounds
  CALL Swap_bounds(Fy_halo,1,n_rows,n_tem_levs,                   &
                   0,offy,fld_type_v,swap_field_is_vector)

  ! bottom level
  DO y=vdims%j_start,vdims%j_end
    divFz(y,1) = (Fz(y,2) - Fz(y,1)) /                          &
                      (zlogc(2) - zlogc(1))
  END DO

  ! top level
  DO y=vdims%j_start,vdims%j_end
    divFz(y,n_tem_levs) = (Fz(y,n_tem_levs) -                       &
                                Fz(y,n_tem_levs-1)) /             &
                          (zlogc(n_tem_levs) - zlogc(n_tem_levs-1))
  END DO

  DO z=2,n_tem_levs-1
    !CDIR NOUNROLL
    DO y=vdims%j_start,vdims%j_end
      divFz(y,z) = (Fz(y,z+1) - Fz(y,z-1)) /                  &
                    (zlogc(z+1) - zlogc(z-1))
    END DO
  END DO

  DO z=1,n_tem_levs
    start_j = vdims%j_start
    end_j = vdims%j_end

    IF (at_extremity(PSouth)) THEN
      start_j = vdims%j_start+1
      divFy(vdims%j_start,z) = (                                &
  (Fy_halo(vdims%j_start,z) - Fy_halo(vdims%j_start+1,z)) /     &
   (phi(vdims%j_start) - phi(vdims%j_start+1)) -                &
          Fy_halo(vdims%j_start,z) * TAN(phi(vdims%j_start))    &
                   ) / planet_radius
    END IF

    IF (at_extremity(PNorth)) THEN
      end_j = vdims%j_end-1
      divFy(vdims%j_end,z) = (                                  &
         (Fy_halo(vdims%j_end-1,z) - Fy_halo(vdims%j_end,z)) /  &
            (phi(vdims%j_end-1) - phi(vdims%j_end)) -           &
             Fy_halo(vdims%j_end,z) * TAN(phi(vdims%j_end))     &
                   ) / planet_radius
    END IF
    !CDIR NOUNROLL
    DO y=start_j,end_j
      divFy(y,z) = (                                      &
                    (Fy_halo(y-1,z) - Fy_halo(y+1,z)) /   &
                     (phi(y-1) - phi(y+1)) -              &
                     Fy_halo(y,z) * TAN(phi(y))           &
                   ) / planet_radius
    END DO
    !
    ! Calculate scaled Div F
    !----------------------------------------------------------------------
    ! Div.F / (ro0 * a cos{phi}
    ! Add code to cover for issues when cos(phi(y))=0 as at ENDGame poles.
    ! set poles to zero to avoid divde by tiny values...

    start_j = vdims%j_start
    end_j   = vdims%j_end

    IF (at_extremity(PSouth)) THEN
      start_j = vdims%j_start+1
      divF(vdims%j_start,z) = 0.0
    END IF

    IF (at_extremity(PNorth)) THEN
      end_j = vdims%j_end-1
      divF(vdims%j_end,z) = 0.0
    END IF
 
    !CDIR NOUNROLL
    DO y=start_j,end_j

      divF(y,z) = (divFy(y,z) + divFz(y,z)) /                &
                   (rho0(z) * planet_radius * COS(phi(y)))
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE Calc_Div_EP_Flux

END MODULE calc_div_ep_flux_mod

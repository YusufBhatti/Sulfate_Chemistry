! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE eg_mass_conserv_mod
IMPLICIT NONE

! Description:
! ----------------------------------------------------------------------
!  This routine assumes rho_species = rho * q
!
!  therefore, it should work for both:
!
!            (a) rho = rho_dry & q = mixing_ratios, and
!            (b) rho = rho_wet & q = specific humidity
!
!
! The quantity we want to conserve is:
!
! (1) SUM{ volume(i,j,k)*rho_np1(i,j,k)*qbar_np1(i,j,k) } = C (known value)
!
! where volume(i,j,k) is the volume of the box surrounding a rho-point
! or rho(i,j,k) is the centre of volume(i,j,k), and,
!
! (2) qbar(i,j,k) = alfa_za(k) * q(i,j,k) + beta_za(k) * q(i,j,k-1)
!
! Using (1) and (2) we can write the SUM in (1) as:
!
! (3) SUM( psi_np1(i,j,k) * q_np1(i,j,k) )  = C
!
! and
!
! (4) C = SUM(   psi_n(i,j,k) *    qs_n(i,j,k) )  +
!         SUM( psi_np1(i,j,k) *    qs_s(i,j,k) )
!         + pseudo_lbflux
!
!     Note that the source_term qs_s are at (n+1) and refers to Physics2
!     the sources due to physics1 are contained within qs_n
!
!  Note also that this routine intrinsically impose the assumption of
!  constant fields at the bottom layer: i.e.,
!
!  (5)   f(:,:,tdims_s%k_start,:) = f(:,:,tdims_s%k_start+1,:),
!
! where f = qs_n, qs_np1 or qs_s, irrespective of whether they are set
! or not [ i.e., the routine doesn't involve f(:,:,tdims_s%k_start,:) ]
!
! ----------------------------------------------------------------------
!
! Method: ENDGame formulation version 3.02
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_MASS_CONSERV_MOD'

CONTAINS

SUBROUTINE eg_mass_conservation_fix(rho_n, rho_np1, qs_n, qs_np1, qs_s,     &
                                    qsmin, number_qs, L_conserv_smooth_lap, &
                                    pseudo_lbflux  )


USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims, tdims_s
USE eg_helmholtz_mod,      ONLY: ec_vol
USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels
USE horiz_grid_mod,        ONLY: xi1_p, xi2_p
USE global_2d_sums_mod,    ONLY: global_2d_sums
USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE Field_Types
USE eg_total_mass_region_mod, ONLY: eg_total_mass_region
USE mpp_conf_mod,          ONLY: swap_field_is_scalar


IMPLICIT NONE
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_MASS_CONSERVATION_FIX'

INTEGER, INTENT(IN) :: number_qs

REAL,    INTENT(IN) ::   qs_n(tdims_s%i_start:tdims_s%i_end,          &
                              tdims_s%j_start:tdims_s%j_end,          &
                              tdims_s%k_start:tdims_s%k_end,          &
                              number_qs)

REAL,    INTENT(IN) ::   qs_s(tdims_s%i_start:tdims_s%i_end,          &
                              tdims_s%j_start:tdims_s%j_end,          &
                              tdims_s%k_start:tdims_s%k_end,          &
                              number_qs)

REAL,    INTENT(INOUT) :: qs_np1(tdims_s%i_start:tdims_s%i_end,       &
                                 tdims_s%j_start:tdims_s%j_end,       &
                                 tdims_s%k_start:tdims_s%k_end,       &
                                 number_qs)

REAL,    INTENT(IN)    ::  rho_n(pdims_s%i_start:pdims_s%i_end,       &
                                 pdims_s%j_start:pdims_s%j_end,       &
                                 pdims_s%k_start:pdims_s%k_end)

REAL,    INTENT(IN)    :: rho_np1(pdims_s%i_start:pdims_s%i_end,      &
                                  pdims_s%j_start:pdims_s%j_end,      &
                                  pdims_s%k_start:pdims_s%k_end)

REAL,    INTENT(IN) :: qsmin (number_qs)
LOGICAL, INTENT(IN) :: L_conserv_smooth_lap
REAL,    INTENT(IN), OPTIONAL :: pseudo_lbflux(number_qs)


! locals

REAL :: alfa_za(pdims%k_start:pdims%k_end)
REAL :: beta_za(pdims%k_start:pdims%k_end)

REAL :: psi_n  (tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                tdims%k_start:tdims%k_end)
REAL :: psi_np1(tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                tdims%k_start:tdims%k_end)

REAL :: local_sums(tdims%i_start:tdims%i_end,                         &
                   tdims%j_start:tdims%j_end,                         &
                   number_qs,3)

REAL :: global_sums(number_qs,3)

REAL :: local_sums2(tdims%i_start:tdims%i_end,                        &
                    tdims%j_start:tdims%j_end,                        &
                    number_qs+1)

REAL :: global_sums2(number_qs+1)

REAL    :: dfdxp,dfdxm,dfdyp,dfdym,dfdzp,dfdzm,delsq,gama

REAL    :: dx2(tdims%i_start  :tdims%i_end)
REAL    :: dx3(tdims%i_start-1:tdims%i_end)
REAL    :: dy2(tdims%j_start  :tdims%j_end)
REAL    :: dy3(tdims%j_start-1:tdims%j_end)
REAL    :: dz2(tdims%k_start  :tdims%k_end)
REAL    :: dz3(tdims%k_start  :tdims%k_end)

REAL    :: mass_deficit, temp1, temp2
REAL    :: one_over_total_lap

REAL    :: small_tol = 1.0e-30

REAL    :: lap (tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                tdims%k_start:tdims%k_end,number_qs)

REAL    ::  qs_lap(tdims%i_start-2:tdims%i_end+2,                     &
                   tdims%j_start-2:tdims%j_end+2,                     &
                   tdims%k_start:tdims%k_end,                         &
                   number_qs)

INTEGER :: i,j,k,kk
INTEGER :: filter_data = 1 ! default
INTEGER :: is, ie, js, je
LOGICAL :: l_exclude_rim


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ------------------------------------------------------------------
! Section 0 : set loop indicies, rho_n, rho_np1 set to be zero for LAM case
! ------------------------------------------------------------------
IF (PRESENT(pseudo_lbflux)) THEN
  l_exclude_rim = .TRUE.
ELSE
  l_exclude_rim = .FALSE.
END IF
CALL eg_total_mass_region(is, ie, js, je, l_exclude_rim)

! ------------------------------------------------------------------
! Section 1 : make qs(n+1) >= qsmin
! ------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,kk)                          &
!$OMP SHARED(tdims,psi_n,rho_n,ec_vol,alfa_za,beta_za,psi_np1,rho_np1,  &
!$OMP        tdims_s, qs_np1, qsmin, pdims, eta_rho_levels,             &
!$OMP        eta_theta_levels, number_qs, dx2, dx3, dy2, dy3, dz2, dz3, &
!$OMP        xi1_p, xi2_p, lap)
DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        qs_np1(i,j,k,kk) = MAX(qs_np1(i,j,k,kk),qsmin(kk))
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END DO

!------------------------------------------------------------------
! section 2: Compute the vertical averaging weigths for Eq.(2) above
!------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO k = pdims%k_start+1, pdims%k_end
  alfa_za(k) = (   eta_rho_levels(k) - eta_theta_levels(k-1) ) /    &
               ( eta_theta_levels(k) - eta_theta_levels(k-1) )
  beta_za(k) = 1.0 - alfa_za(k)
END DO
!$OMP END DO
!$OMP SINGLE
alfa_za(pdims%k_start) = 1.0
beta_za(pdims%k_start) = 0.0
!$OMP END SINGLE

!-------------------------------------------------------------
! section 3: Compute the 3D array psi(i,j,k) for Eqs (3) & (4)
!-------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO k = tdims%k_start + 1, tdims%k_end - 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      psi_n(i,j,k)  =  rho_n(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )    &
                    +  rho_n(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)

      psi_np1(i,j,k) = rho_np1(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )  &
                     + rho_np1(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

k = tdims%k_end
!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    psi_n(i,j,k) =   rho_n(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
    psi_np1(i,j,k) = rho_np1(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
  END DO
END DO
!$OMP END DO NOWAIT

! --------------------------------------------------------------------
!  section 4:  Calculate the laplacian (div) of qs. This div is computed
!              like Cartesian (and this consistent with interpolation).
!              These calculations can be simplified for constant mesh
!              (here a general(variable) mesh is assumed)
! --------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO i = tdims%i_start, tdims%i_end
  dx2(i) = 0.5 * ( xi1_p(i+1) - xi1_p(i-1) )
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO i = tdims%i_start - 1 , tdims%i_end
  dx3(i) = xi1_p(i+1) - xi1_p(i)
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start, tdims%j_end
  dy2(j) = 0.5 * ( xi2_p(j+1) - xi2_p(j-1) )
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start - 1, tdims%j_end
  dy3(j) = xi2_p(j+1) - xi2_p(j)
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = tdims%k_start+1 , tdims%k_end - 1
  dz2(k) = 0.5 * ( eta_theta_levels(k+1) - eta_theta_levels(k-1) )
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = tdims%k_start , tdims%k_end - 1
  dz3(k) = eta_theta_levels(k+1) - eta_theta_levels(k)
END DO
!$OMP END DO NOWAIT

DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        lap(i,j,k,kk) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END DO
!$OMP END PARALLEL

IF ( filter_data == 1 .AND. L_conserv_smooth_lap ) THEN  

  ! Compute the laplacian of filtered qs_np1 (using 1-2-1 filter in
  ! all directions). This is equivalent to apply the Laplacian to
  ! a coarse grid (similar to Multigrid).
  ! These filtered data are used only in computing "lap" as
  ! a measure of smoothness -- this is to avoid giving too much
  ! weight to sub-grid noise

!$OMP PARALLEL DEFAULT(NONE)                                    &
!$OMP& PRIVATE(i,j,k,kk) SHARED(tdims,qs_lap,qs_np1,number_qs)
  DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qs_lap(i,j,k,kk) = qs_np1(i,j,k,kk)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL

  i =  tdims%i_len
  j =  tdims%j_len
  k =  tdims%k_len*number_qs

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds(qs_lap,i,j,k,2,2,fld_type_p, swap_field_is_scalar)


!$OMP PARALLEL DEFAULT(NONE)                                            &
!$OMP& PRIVATE(i,j,k,kk,dfdxp,dfdxm,dfdyp,dfdym,dfdzp,dfdzm,delsq)      &
!$OMP& SHARED(tdims,qs_lap,dx3,dx2,dy3,dy2,dz3,dz2,lap,qs_np1,psi_np1,  &
!$OMP& number_qs)
  DO kk = 1, number_qs

    k = tdims%k_start + 1
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        dfdxp  = (qs_lap(i+2,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
        dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-2,j,k,kk))/dx3(i-1)

        dfdyp  = (qs_lap(i,j+2,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
        dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-2,k,kk))/dy3(j-1)

        dfdzp  = (qs_lap(i,j,k+1,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
        dfdzm  = 0.0

        delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                 (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                 (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k)

        lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)

      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start + 2, tdims%k_end - 2
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          dfdxp  = (qs_lap(i+2,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-2,j,k,kk))/dx3(i-1)

          dfdyp  = (qs_lap(i,j+2,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-2,k,kk))/dy3(j-1)

          dfdzp  = (qs_lap(i,j,k+2,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
          dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-2,kk))/dz3(k-1)

          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                   (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k)

          lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT


    k = tdims%k_end - 1
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        dfdxp  = (qs_lap(i+2,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
        dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-2,j,k,kk))/dx3(i-1)

        dfdyp  = (qs_lap(i,j+2,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
        dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-2,k,kk))/dy3(j-1)

        dfdzp  = (qs_lap(i,j,k+1,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
        dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-1,kk))/dz3(k-1)

        delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                 (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                 (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k)

        lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)

      END DO
    END DO
!$OMP END DO NOWAIT

    k = tdims%k_end
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
        dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))/dx3(i-1)

        dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
        dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))/dy3(j-1)

        delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) + &
                 (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j)

        lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)

      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL


ELSE

  IF ( filter_data == 1  ) THEN  

    ! filter qs_np1 (apply 1-2-1 filter in all directions)

    DO kk = 1, number_qs
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                      &
!$OMP& PRIVATE(i,j,k) SHARED(tdims,qs_lap,qs_np1,kk)
      DO k = tdims%k_start+1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            qs_lap(i,j,k,kk) = ( qs_np1(i-1,j-1,k,kk)                   &
                                     + qs_np1(i+1,j-1,k,kk)             &
                                     + qs_np1(i-1,j+1,k,kk)             &
                                     + qs_np1(i+1,j+1,k,kk)             &
                                     + 2.0*( qs_np1(i,j-1,k,kk)         &
                                           + qs_np1(i,j+1,k,kk)         &
                                           + qs_np1(i-1,j,k,kk)         &
                                           + qs_np1(i+1,j,k,kk) )       &
                                     + 4.0*qs_np1(i,j,k,kk) )/16.0
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    END DO


    i = tdims%i_len
    j = tdims%j_len
    k = tdims%k_len*number_qs

    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(qs_lap,i,j,k,2,2,fld_type_p, swap_field_is_scalar)

  ELSE !  use unfiltred raw data

!$OMP PARALLEL DEFAULT(NONE)                                             &
!$OMP& PRIVATE(i,j,k,kk) SHARED(tdims,qs_lap,qs_np1,number_qs)
    DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
      DO k = tdims%k_start+1, tdims%k_end
        DO j = tdims%j_start-1, tdims%j_end+1
          DO i = tdims%i_start-1, tdims%i_end+1
            qs_lap(i,j,k,kk) = qs_np1(i,j,k,kk)
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
!$OMP END PARALLEL

  END IF

!$OMP PARALLEL DEFAULT(NONE)                                            &
!$OMP& PRIVATE(i,j,k,kk,dfdxp,dfdxm,dfdyp,dfdym,dfdzp,dfdzm,            &
!$OMP& delsq) SHARED(tdims,qs_lap,dx3,dx2,dy3,dy2,dz3,                  &
!$OMP& dz2,lap,qs_np1,psi_np1,number_qs)

  DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start-2, tdims%j_end+2
      DO i = tdims%i_start-2, tdims%i_end+2
        qs_lap(i,j,tdims%k_start,kk)=qs_lap(i,j,tdims%k_start+1,kk)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO

! Need a barrier as dependencies from above loop set for the next
!$OMP BARRIER

  DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start + 1, tdims%k_end - 1
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))/dx3(i-1)

          dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))/dy3(j-1)

          dfdzp  = (qs_lap(i,j,k+1,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
          dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-1,kk))/dz3(k-1)

          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                   (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k)

          lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    k = tdims%k_end
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
        dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))/dx3(i-1)

        dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
        dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))/dy3(j-1)

        delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) + &
                 (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j)

        lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)

      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL

END IF

! ----------------------------------------------------------------------
!  section 5:  Calculate the global sums in (3) & (4) and total lap
! ----------------------------------------------------------------------

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,k,kk) DEFAULT(NONE)    &
!$OMP& SHARED(local_sums,tdims,psi_n,psi_np1,qs_n,qs_s,qs_np1,lap,    &
!$OMP& number_qs,is,ie,js,je)
DO kk = 1, number_qs

  local_sums(:,:,kk,1) = 0.0
  local_sums(:,:,kk,2) = 0.0
  local_sums(:,:,kk,3) = 0.0

  DO k = tdims%k_start+1, tdims%k_end
    DO j = js, je
      DO i = is, ie

        local_sums(i,j,kk,1) =  local_sums(i,j,kk,1)                &
                                +    psi_n(i,j,k) * qs_n(i,j,k,kk)  &
                                +  psi_np1(i,j,k) * qs_s(i,j,k,kk)

        local_sums(i,j,kk,2) =  local_sums(i,j,kk,2)                &
                                +  psi_np1(i,j,k) * qs_np1(i,j,k,kk)

        local_sums(i,j,kk,3) =  local_sums(i,j,kk,3)                &
                                  +  lap(i,j,k,kk)
      END DO
    END DO
  END DO

END DO
!$OMP END PARALLEL DO


CALL global_2d_sums(local_sums, tdims%i_len,          &
                                       tdims%j_len,   &
                                       0, 0, 3*number_qs, global_sums )

IF (PRESENT(pseudo_lbflux)) THEN
  global_sums(:,1) = global_sums(:,1) + pseudo_lbflux(:)
END IF

! ----------------------------------------------------------------------
!  section 6:  Modify qs to achieve mass conservation
! ----------------------------------------------------------------------

DO kk = 1, number_qs
  IF ( ABS(global_sums(kk,3)) > small_tol ) THEN
    mass_deficit = global_sums(kk,1) -  global_sums(kk,2)
    one_over_total_lap = 1.0/global_sums(kk,3)
!$OMP PARALLEL DO  SCHEDULE(STATIC) DEFAULT(NONE)                     &
!$OMP& PRIVATE(i,j,k, gama) SHARED(tdims,lap,mass_deficit,            &
!$OMP& one_over_total_lap,qs_np1,psi_np1,kk,is,ie,js,je)
    DO k = tdims%k_start + 1, tdims%k_end
      DO j = js, je
        DO i = is, ie
          gama = lap(i,j,k,kk) * mass_deficit * one_over_total_lap
          qs_np1(i,j,k,kk) = qs_np1(i,j,k,kk) + gama/psi_np1(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END IF
END DO


!=======================================================
!  section 7: Apply a second correction to force q >= qmin
!             while maintaining mass conservation
!=======================================================

! clip the field points that are below qsmin

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,kk)                          &
!$OMP SHARED(tdims_s, qs_np1, qsmin, number_qs, local_sums2, js,je,is,  &
!$OMP        ie, tdims, psi_np1)
DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        qs_np1(i,j,k,kk) = MAX(qs_np1(i,j,k,kk),qsmin(kk))
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END DO

! compute new masses at (n+1) after cliping
DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
  DO j = js, je
    DO i = is, ie
      local_sums2(i,j,kk) = 0.0
    END DO
  END DO
!$OMP END DO

  DO k = tdims%k_start+1, tdims%k_end
!$OMP DO SCHEDULE(STATIC)
    DO j = js, je
      DO i = is, ie
        local_sums2(i,j,kk) = local_sums2(i,j,kk)                    &
                            + psi_np1(i,j,k) * qs_np1(i,j,k,kk)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
END DO

! this sum doesn't invlove qs (it's the sum above with q=1)
!$OMP DO SCHEDULE(STATIC)
DO j = js, je
  DO i = is, ie
    local_sums2(i,j,number_qs+1) = 0.0
  END DO
END DO
!$OMP END DO
DO k = tdims%k_start+1, tdims%k_end
!$OMP DO SCHEDULE(STATIC)
  DO j = js, je
    DO i = is, ie
      local_sums2(i,j,number_qs+1) = local_sums2(i,j,number_qs+1)    &
                                   + psi_np1(i,j,k)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO

!$OMP END PARALLEL

CALL global_2d_sums(local_sums2, tdims%i_len,       &
                    tdims%j_len,                    &
                    0, 0, number_qs+1, global_sums2                 )

DO kk = 1, number_qs

  mass_deficit = global_sums2(kk)  - qsmin(kk)*global_sums2(number_qs+1)

  IF ( ABS(mass_deficit) > small_tol ) THEN
    temp1 = ( global_sums(kk,1)-qsmin(kk)*global_sums2(number_qs+1) )/&
            mass_deficit
    temp2 = qsmin(kk)*(1.0 - temp1)

!$OMP PARALLEL DO DEFAULT(NONE)          &
!$OMP& PRIVATE(k, j, i)                  &
!$OMP& SHARED(tdims, js, je, is, ie, kk, &
!$OMP&        qs_np1, temp1, temp2)
    DO k = tdims%k_start, tdims%k_end
      DO j = js, je
        DO i = is, ie
          qs_np1(i,j,k,kk) = temp1 * qs_np1(i,j,k,kk) + temp2
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END IF
END DO

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,kk)                          &
!$OMP          SHARED(number_qs, tdims_s, qs_np1, tdims)
DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims_s%j_start, tdims_s%j_end
    DO i = tdims_s%i_start, tdims_s%i_end
      qs_np1(i,j,tdims%k_start,kk) = qs_np1(i,j,tdims%k_start+1,kk)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE eg_mass_conservation_fix

!====================================================================
!
! Version to be used if l_fix_conserv = .FALSE.
!
SUBROUTINE eg_mass_conservation(rho_n, rho_np1, qs_n, qs_np1, qs_s, &
                             qsmin, qswitches, number_qs,           &
                             L_conserv_smooth_lap, pseudo_lbflux    )


USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims, tdims_s
USE eg_helmholtz_mod,      ONLY: ec_vol
USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels
USE horiz_grid_mod,        ONLY: xi1_p, xi2_p
USE global_2d_sums_mod,    ONLY: global_2d_sums
USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE Field_Types
USE eg_total_mass_region_mod, ONLY: eg_total_mass_region
USE mpp_conf_mod,          ONLY: swap_field_is_scalar

IMPLICIT NONE
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_MASS_CONSERVATION'

INTEGER, INTENT(IN) :: number_qs

REAL,    INTENT(IN) ::   qs_n(tdims_s%i_start:tdims_s%i_end,          &
                              tdims_s%j_start:tdims_s%j_end,          &
                              tdims_s%k_start:tdims_s%k_end,          &
                              number_qs)

REAL,    INTENT(IN) ::   qs_s(tdims_s%i_start:tdims_s%i_end,          &
                              tdims_s%j_start:tdims_s%j_end,          &
                              tdims_s%k_start:tdims_s%k_end,          &
                              number_qs)

REAL,    INTENT(INOUT) :: qs_np1(tdims_s%i_start:tdims_s%i_end,       &
                                 tdims_s%j_start:tdims_s%j_end,       &
                                 tdims_s%k_start:tdims_s%k_end,       &
                                 number_qs)

REAL,    INTENT(IN)    ::  rho_n(pdims_s%i_start:pdims_s%i_end,       &
                                 pdims_s%j_start:pdims_s%j_end,       &
                                 pdims_s%k_start:pdims_s%k_end)

REAL,    INTENT(IN)    :: rho_np1(pdims_s%i_start:pdims_s%i_end,      &
                                  pdims_s%j_start:pdims_s%j_end,      &
                                  pdims_s%k_start:pdims_s%k_end)

REAL,    INTENT(IN) :: qsmin (number_qs)
INTEGER, INTENT(IN) :: qswitches (number_qs)
LOGICAL, INTENT(IN) :: L_conserv_smooth_lap
REAL,    INTENT(IN), OPTIONAL :: pseudo_lbflux(number_qs)

! locals

REAL :: alfa_za(pdims%k_start:pdims%k_end)
REAL :: beta_za(pdims%k_start:pdims%k_end)

REAL :: psi_n  (tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                tdims%k_start:tdims%k_end)
REAL :: psi_np1(tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                tdims%k_start:tdims%k_end)

REAL :: local_sums(tdims%i_start:tdims%i_end,                         &
                   tdims%j_start:tdims%j_end,                         &
                   number_qs,3)

REAL :: global_sums(number_qs,3)

REAL :: local_sums2(tdims%i_start:tdims%i_end,                        &
                    tdims%j_start:tdims%j_end,                        &
                    number_qs+1)

REAL :: global_sums2(number_qs+1)

REAL    :: dfdxp,dfdxm,dfdyp,dfdym,dfdzp,dfdzm,delsq,gama

REAL    :: dx2(tdims%i_start  :tdims%i_end)
REAL    :: dx3(tdims%i_start-1:tdims%i_end)
REAL    :: dy2(tdims%j_start  :tdims%j_end)
REAL    :: dy3(tdims%j_start-1:tdims%j_end)
REAL    :: dz2(tdims%k_start  :tdims%k_end)
REAL    :: dz3(tdims%k_start  :tdims%k_end)

REAL    :: mass_deficit, temp1, temp2
REAL    :: one_over_total_lap

REAL    :: small_tol = 1.0e-30

REAL    :: lap (tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                tdims%k_start:tdims%k_end,number_qs)

REAL    ::  qs_lap(tdims%i_start-2:tdims%i_end+2,                     &
                   tdims%j_start-2:tdims%j_end+2,                     &
                   tdims%k_start:tdims%k_end,                         &
                   number_qs)

INTEGER :: i,j,k,kk
INTEGER :: filter_data = 1 ! default
INTEGER :: IS, ie, js, je
LOGICAL :: l_exclude_rim

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ------------------------------------------------------------------
! Section 0 : set loop indicies, rho_n, rho_np1 set to be zero for LAM case
! ------------------------------------------------------------------
IF (PRESENT(pseudo_lbflux)) THEN
  l_exclude_rim = .TRUE.
ELSE
  l_exclude_rim = .FALSE.
END IF
CALL eg_total_mass_region(IS, ie, js, je, L_exclude_rim)

!------------------------------------------------------------------
! section 1: Compute the vertical averaging weigths for Eq.(2) above
!------------------------------------------------------------------

DO k = pdims%k_start, pdims%k_end

  alfa_za(k) = (   eta_rho_levels(k) - eta_theta_levels(k-1) ) /    &
               ( eta_theta_levels(k) - eta_theta_levels(k-1) )
  beta_za(k) = 1.0 - alfa_za(k)

END DO

! ------------------------------------------------------------------
! Section 2 : make qs(n+1) > qsmin if necessary
! ------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,kk)                      &
!$OMP SHARED(tdims,psi_n,rho_n,ec_vol,alfa_za,beta_za,psi_np1,      &
!$OMP     rho_np1,qs_np1,qsmin,number_qs)


DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    qs_np1(:,:,k,kk) = MAX(qs_np1(:,:,k,kk),qsmin(kk))
  END DO
!$OMP END DO NOWAIT
END DO

!-------------------------------------------------------------
! section 3: Compute the 3D array psi(i,j,k) for Eqs (3) & (4)
!-------------------------------------------------------------

k = tdims%k_start

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    psi_n(i,j,k) =   rho_n(i,j,k+1) * ec_vol(i,j,k+1) * beta_za(k+1)
    psi_np1(i,j,k) = rho_np1(i,j,k+1) * ec_vol(i,j,k+1) * beta_za(k+1)
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = tdims%k_start + 1, tdims%k_end - 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      psi_n(i,j,k)  =  rho_n(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )    &
                    +  rho_n(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)

      psi_np1(i,j,k) = rho_np1(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )  &
                     + rho_np1(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

k = tdims%k_end

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    psi_n(i,j,k)   =   rho_n(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
    psi_np1(i,j,k) = rho_np1(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

! --------------------------------------------------------------------
!  section 4:  Calculate the laplacian (div) of qs. This div is computed
!              like Cartesian (and this consistent with interpolation).
!              These calculations can be simplified for constant mesh
!              (here a general(variable) mesh is assumed)
! --------------------------------------------------------------------

DO i = tdims%i_start, tdims%i_end
  dx2(i) = 0.5 * ( xi1_p(i+1) - xi1_p(i-1) )
END DO

DO i = tdims%i_start - 1 , tdims%i_end
  dx3(i) = xi1_p(i+1) - xi1_p(i)
END DO

DO j = tdims%j_start, tdims%j_end
  dy2(j) = 0.5 * ( xi2_p(j+1) - xi2_p(j-1) )
END DO

DO j = tdims%j_start - 1, tdims%j_end
  dy3(j) = xi2_p(j+1) - xi2_p(j)
END DO

DO k = tdims%k_start+1 , tdims%k_end - 1
  dz2(k) = 0.5 * ( eta_theta_levels(k+1) - eta_theta_levels(k-1) )
END DO

DO k = tdims%k_start , tdims%k_end - 1
  dz3(k) = eta_theta_levels(k+1) - eta_theta_levels(k)
END DO

! Reset array
lap = 0.0

IF ( filter_data == 1 .AND. L_conserv_smooth_lap ) THEN  

  ! Compute the laplacian of filtered qs_np1 (using 1-2-1 filter in
  ! all directions). This is equivalent to apply the Laplacian to
  ! a coarse grid (similar to Multigrid).
  ! These filtered data are used only in computing "lap" as
  ! a measure of smoothness -- this is to avoid giving too much
  ! weight to sub-grid noise

  DO kk = 1, number_qs
!$OMP PARALLEL DO  SCHEDULE(STATIC) DEFAULT(NONE)               &
!$OMP PRIVATE(i,j,k) SHARED(tdims,qs_lap,qs_np1,kk)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qs_lap(i,j,k,kk) = qs_np1(i,j,k,kk)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO

  i =  tdims%i_len
  j =  tdims%j_len
  k =  tdims%k_len*number_qs

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds(qs_lap,i,j,k,2,2,fld_type_p, swap_field_is_scalar)


  DO kk = 1, number_qs

    IF ( qswitches(kk) == 1 ) THEN

      k = tdims%k_start
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          dfdxp  = (qs_lap(i+2,j,k,kk) - qs_lap(i  ,j,k,kk))/ dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-2,j,k,kk))/ dx3(i-1)

          dfdyp  = (qs_lap(i,j+2,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-2,k,kk))/dy3(j-1)

          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +         &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j)

          lap(i,j,k,kk) = ABS(delsq) * psi_np1(i,j,k)

        END DO
      END DO

      k = tdims%k_start + 1
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          dfdxp  = (qs_lap(i+2,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-2,j,k,kk))/dx3(i-1)

          dfdyp  = (qs_lap(i,j+2,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-2,k,kk))/dy3(j-1)

          dfdzp  = (qs_lap(i,j,k+1,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
          dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-1,kk))/dz3(k-1)

          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                   (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k)

          lap(i,j,k,kk) = ABS(delsq) * psi_np1(i,j,k)

        END DO
      END DO


!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP PRIVATE(i,j,k,dfdxp,dfdxm,dfdyp,dfdym,dfdzp,dfdzm,delsq)         &
!$OMP SHARED(tdims,qs_lap,dx3,dx2,dy3,dy2,dz3,dz2,lap,qs_np1,psi_np1,  &
!$OMP kk)

      DO k = tdims%k_start + 2, tdims%k_end - 2
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            dfdxp  = (qs_lap(i+2,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
            dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-2,j,k,kk))/dx3(i-1)

            dfdyp  = (qs_lap(i,j+2,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
            dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-2,k,kk))/dy3(j-1)

            dfdzp  = (qs_lap(i,j,k+2,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
            dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-2,kk))/dz3(k-1)

            delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                     (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                     (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k)

            lap(i,j,k,kk) = ABS(delsq) * psi_np1(i,j,k)

          END DO
        END DO
      END DO

!$OMP END PARALLEL DO


      k = tdims%k_end - 1
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          dfdxp  = (qs_lap(i+2,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-2,j,k,kk))/dx3(i-1)

          dfdyp  = (qs_lap(i,j+2,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-2,k,kk))/dy3(j-1)

          dfdzp  = (qs_lap(i,j,k+1,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
          dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-1,kk))/dz3(k-1)

          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                   (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k)

          lap(i,j,k,kk) = ABS(delsq) * psi_np1(i,j,k)

        END DO
      END DO

      k = tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))/dx3(i-1)

          dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))/dy3(j-1)

          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) + &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j)

          lap(i,j,k,kk) = ABS(delsq) * psi_np1(i,j,k)

        END DO
      END DO

    END IF
  END DO


ELSE

!$OMP PARALLEL DEFAULT(NONE)                                            &
!$OMP PRIVATE(i,j,k, dfdxp,dfdxm,dfdyp,dfdym,dfdzp,dfdzm,               &
!$OMP delsq,kk)                                                         &
!$OMP SHARED(tdims,qs_lap,number_qs,filter_data,dx3,dx2,dy3,dy2,dz3,    &
!$OMP    dz2,lap,qs_np1,psi_np1,qswitches)

  IF ( filter_data == 1  ) THEN  

    ! filter qs_np1 (apply 1-2-1 filter in all directions)

    DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
      DO k = tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            qs_lap(i,j,k,kk) = ( qs_np1(i-1,j-1,k,kk)                   &
                                     + qs_np1(i+1,j-1,k,kk)             &
                                     + qs_np1(i-1,j+1,k,kk)             &
                                     + qs_np1(i+1,j+1,k,kk)             &
                                     + 2.0*( qs_np1(i,j-1,k,kk)         &
                                           + qs_np1(i,j+1,k,kk)         &
                                           + qs_np1(i-1,j,k,kk)         &
                                           + qs_np1(i+1,j,k,kk) )       &
                                     + 4.0*qs_np1(i,j,k,kk) )/16.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO

! Wait for qs_lap
!$OMP BARRIER
!$OMP MASTER

    i = tdims%i_len
    j = tdims%j_len
    k = tdims%k_len*number_qs

    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(qs_lap,i,j,k,2,2,fld_type_p, swap_field_is_scalar)
!$OMP END MASTER

  ELSE !  use unfiltred raw data

    DO kk = 1, number_qs
!$OMP DO  SCHEDULE(STATIC)
      DO k = tdims%k_start, tdims%k_end
        DO j = tdims%j_start-1, tdims%j_end+1
          DO i = tdims%i_start-1, tdims%i_end+1
            qs_lap(i,j,k,kk) = qs_np1(i,j,k,kk)
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO

  END IF

! Wait for qs_lap
!$OMP BARRIER

  DO kk = 1, number_qs

    IF ( qswitches(kk) == 1 ) THEN

      k = tdims%k_start
!$OMP DO  SCHEDULE(STATIC)
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          
          dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))/dx3(i-1)
          
          dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))/dy3(j-1)
          
          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +           &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j)

          lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)

        END DO
      END DO
!$OMP END DO NOWAIT
    
!$OMP DO SCHEDULE(STATIC)
      DO k = tdims%k_start + 1, tdims%k_end - 1
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
            dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))/dx3(i-1)

            dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
            dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))/dy3(j-1)

            dfdzp  = (qs_lap(i,j,k+1,kk) - qs_lap(i,j,k  ,kk))/dz3(k  )
            dfdzm  = (qs_lap(i,j,k  ,kk) - qs_lap(i,j,k-1,kk))/dz3(k-1)

            delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) +  &
                     (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j) +  &
                     (dz3(k)*dz3(k-1)) * (dfdzp-dfdzm)/dz2(k)

            lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)

          END DO
        END DO
      END DO
!$OMP END DO NOWAIT

      k = tdims%k_end
!$OMP DO  SCHEDULE(STATIC)
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          dfdxp  = (qs_lap(i+1,j,k,kk) - qs_lap(i  ,j,k,kk))/dx3(i  )
          dfdxm  = (qs_lap(i  ,j,k,kk) - qs_lap(i-1,j,k,kk))/dx3(i-1)

          dfdyp  = (qs_lap(i,j+1,k,kk) - qs_lap(i,j  ,k,kk))/dy3(j  )
          dfdym  = (qs_lap(i,j  ,k,kk) - qs_lap(i,j-1,k,kk))/dy3(j-1)

          delsq  = (dx3(i)*dx3(i-1)) * (dfdxp-dfdxm)/dx2(i) + &
                   (dy3(j)*dy3(j-1)) * (dfdyp-dfdym)/dy2(j)

          lap(i,j,k,kk) = ABS(delsq)*ABS(qs_np1(i,j,k,kk))*psi_np1(i,j,k)

        END DO
      END DO
!$OMP END DO NOWAIT

    END IF
  END DO ! kk

!$OMP END PARALLEL

END IF

! ----------------------------------------------------------------------
!  section 5:  Calculate the global sums in (3) & (4) and total lap
! ----------------------------------------------------------------------

!$OMP PARALLEL PRIVATE(i,j,k,kk) DEFAULT(NONE)                        &
!$OMP SHARED(local_sums,tdims,psi_n,psi_np1,qs_n,qs_s,qs_np1,lap,     &
!$OMP qswitches,number_qs,is,ie,js,je)

DO kk = 1, number_qs

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      local_sums(i,j,kk,1) = 0.0
      local_sums(i,j,kk,2) = 0.0
      local_sums(i,j,kk,3) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT
END DO

!$OMP BARRIER

DO kk = 1, number_qs
  IF ( qswitches(kk) == 1 ) THEN
    DO k = tdims%k_start, tdims%k_end
!$OMP DO SCHEDULE(STATIC) 
      DO j = js, je
        DO i = IS, ie

          local_sums(i,j,kk,1) =  local_sums(i,j,kk,1)                &
                                  +    psi_n(i,j,k) * qs_n(i,j,k,kk)  &
                                  +  psi_np1(i,j,k) * qs_s(i,j,k,kk)

          local_sums(i,j,kk,2) =  local_sums(i,j,kk,2)                &
                                  +  psi_np1(i,j,k) * qs_np1(i,j,k,kk)

          local_sums(i,j,kk,3) =  local_sums(i,j,kk,3)                &
                                    +  lap(i,j,k,kk)
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
  END IF
END DO

!$OMP END PARALLEL


CALL global_2d_sums(local_sums, tdims%i_len,          &
                                       tdims%j_len,   &
                                       0, 0, 3*number_qs, global_sums )

IF (PRESENT(pseudo_lbflux)) THEN
  global_sums(:,1) = global_sums(:,1) + pseudo_lbflux(:)
END IF

! ----------------------------------------------------------------------
!  section 6:  Modify qs to achieve mass conservation
! ----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,gama,kk,mass_deficit,      &
!$OMP one_over_total_lap,temp1,temp2)                                 &
!$OMP SHARED(tdims,lap,number_qs,global_sums,small_tol,qsmin,         &
!$OMP    qs_np1,psi_np1,is,ie,js,je,local_sums2,global_sums2,         &
!$OMP    qswitches)

DO kk = 1, number_qs
  IF ( qswitches(kk) == 1 .AND. ABS(global_sums(kk,3)) > small_tol ) THEN
    mass_deficit = global_sums(kk,1) -  global_sums(kk,2)
    one_over_total_lap = 1.0/global_sums(kk,3)
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = js, je
        DO i = IS, ie
          gama = lap(i,j,k,kk) * mass_deficit * one_over_total_lap
          qs_np1(i,j,k,kk) = qs_np1(i,j,k,kk) + gama/psi_np1(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
END DO

!=======================================================
!  section 7: Apply a second correction to force q >= qmin
!      while maintaining mass conservation
!=======================================================

! clip the field points that are below qsmin

DO kk = 1, number_qs
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    qs_np1(:,:,k,kk) = MAX(qs_np1(:,:,k,kk),qsmin(kk))
  END DO
!$OMP END DO NOWAIT
END DO

! compute new masses at (n+1) after cliping

DO kk = 1, number_qs+1
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      local_sums2(i,j,kk) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT
END DO

! Wait for local_sums2 and qs_np1
!$OMP BARRIER

DO kk = 1, number_qs
  DO k = tdims%k_start, tdims%k_end
!$OMP DO SCHEDULE(STATIC)
    DO j = js, je
      DO i = IS, ie

        local_sums2(i,j,kk) = local_sums2(i,j,kk)                    &
                            + psi_np1(i,j,k) * qs_np1(i,j,k,kk)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
END DO

! this sum doesn't invlove qs (it's the sum above with q=1)

DO k = tdims%k_start, tdims%k_end
!$OMP DO SCHEDULE(STATIC)
  DO j = js, je
    DO i = IS, ie
      local_sums2(i,j,number_qs+1) = local_sums2(i,j,number_qs+1)    &
                                   + psi_np1(i,j,k)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO

!$OMP BARRIER
!$OMP MASTER
CALL global_2d_sums(local_sums2, tdims%i_len,         &
                    tdims%j_len,                      &
                    0, 0, number_qs+1, global_sums2 )
!$OMP END MASTER
!$OMP BARRIER

DO kk = 1, number_qs

  mass_deficit = global_sums2(kk)  - qsmin(kk)*global_sums2(number_qs+1)

  IF ( qswitches(kk) == 1 .AND. ABS(mass_deficit) > small_tol ) THEN
    temp1 = ( global_sums(kk,1)-qsmin(kk)*global_sums2(number_qs+1) )/&
            mass_deficit
    temp2 = qsmin(kk)*(1.0 - temp1)

!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = js, je
        DO i = IS, ie
          qs_np1(i,j,k,kk) = temp1 * qs_np1(i,j,k,kk) + temp2
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

  END IF
END DO

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE eg_mass_conservation

END MODULE eg_mass_conserv_mod


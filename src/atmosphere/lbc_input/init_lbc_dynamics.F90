! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE init_lbc_dynamics(u,v,w, etadot,                                 &
                             u_lbc, v_lbc, w_lbc,                           &
                             u_lbc_tnd, v_lbc_tnd, w_lbc_tnd,               &
                             inc_fact, rim_stepsa, timestep_number,         &
                             row_length, rows, n_rows, model_levels,        &
                             offx, offy, halo_i, halo_j,                    &
                             lenrim, rimwidth, rimweights,                  &
                             lbc_size, lbc_start,                           &
                             at_extremity,                                  &
                             L_do_boundaries, L_do_halos                    &
                            )

USE um_parparams, ONLY: Nhalo_max, halo_type_extended, Pnorth, Psouth,      &
                        Peast, Pwest
USE field_types, ONLY: nfld_max, fld_type_u, fld_type_v, fld_type_w,        &
                       fld_type_p
USE atm_fields_bounds_mod, ONLY: wdims, udims_s, vdims_s, wdims_s
USE metric_terms_mod, ONLY: deta_xi3_theta, h1_p_eta, h2_p_eta, h3_p_eta,   &
                            dxi1_xi3,dxi2_xi3
USE horiz_grid_mod, ONLY: intw_rho2w, intw_u2p, intw_v2p
USE lam_lbc_weights_mod, ONLY: ns_weights, ew_weights
USE mpp_conf_mod,   ONLY: swap_field_is_scalar

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!
! Description: Update tendancies for LBC's and copy into field
!              (used by ENDGame to simplify atm_step).
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: LBC Input
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

INTEGER, INTENT(IN)    :: row_length, rows, n_rows, model_levels
INTEGER, INTENT(IN)    :: offx, offy, halo_i, halo_j
INTEGER, INTENT(IN)    :: rim_stepsa, timestep_number
INTEGER, INTENT(IN)    :: rimwidth, lenrim(Nfld_max,NHalo_max),             &
                                lbc_size(4,Nfld_max,NHalo_max),             &
                               lbc_start(4,Nfld_max,NHalo_max)

REAL,    INTENT(IN)    :: rimweights(rimwidth)

REAL,    INTENT(INOUT) :: inc_fact
LOGICAL, INTENT(INOUT) :: at_extremity(4)
LOGICAL, INTENT(INOUT) :: L_do_boundaries, L_do_halos

REAL,    INTENT(INOUT) :: u(udims_s%i_start:udims_s%i_end,                  &
                            udims_s%j_start:udims_s%j_end,                  &
                            udims_s%k_start:udims_s%k_end)
REAL,    INTENT(INOUT) :: v(vdims_s%i_start:vdims_s%i_end,                  &
                            vdims_s%j_start:vdims_s%j_end,                  &
                            vdims_s%k_start:vdims_s%k_end)
REAL,    INTENT(INOUT) :: w(wdims_s%i_start:wdims_s%i_end,                  &
                            wdims_s%j_start:wdims_s%j_end,                  &
                            wdims_s%k_start:wdims_s%k_end),                 &
                     etadot(wdims_s%i_start:wdims_s%i_end,                  &
                            wdims_s%j_start:wdims_s%j_end,                  &
                            wdims_s%k_start:wdims_s%k_end)

REAL,    INTENT (IN)   :: u_lbc(lenrim(fld_type_u,halo_type_extended),      &
                                udims_s%k_start:udims_s%k_end),             &
                      u_lbc_tnd(lenrim(fld_type_u,halo_type_extended),      &
                                udims_s%k_start:udims_s%k_end)
REAL,    INTENT (IN)   :: v_lbc(lenrim(fld_type_v,halo_type_extended),      &
                                vdims_s%k_start:vdims_s%k_end),             &
                      v_lbc_tnd(lenrim(fld_type_v,halo_type_extended),      &
                                vdims_s%k_start:vdims_s%k_end)
REAL,    INTENT (IN)   :: w_lbc(lenrim(fld_type_p,halo_type_extended),      &
                                wdims_s%k_start:wdims_s%k_end),             &
                      w_lbc_tnd(lenrim(fld_type_p,halo_type_extended),      &
                                wdims_s%k_start:wdims_s%k_end)

INTEGER                :: i, j, k, lbc_len
INTEGER                :: ii, jj


REAL                   :: u_tend(lenrim(fld_type_u,halo_type_extended),     &
                                 udims_s%k_start:udims_s%k_end),            &
                          v_tend(lenrim(fld_type_v,halo_type_extended),     &
                                 vdims_s%k_start:vdims_s%k_end),            &
                          w_tend(lenrim(fld_type_p,halo_type_extended),     &
                                 wdims_s%k_start:wdims_s%k_end)

REAL                   :: etadot_rim

REAL                   :: u_at_w, v_at_w

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_LBC_DYNAMICS'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

l_do_halos      = .TRUE.
l_do_boundaries = .TRUE.

IF (rim_stepsa == 0) THEN
  inc_fact = 0.0
ELSE
  inc_fact = 1.0/ (rim_stepsa-MOD(timestep_number-1,rim_stepsa))
END IF

! update u, v, w

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,k,lbc_len)                     &
!$OMP SHARED(lenrim, u_tend,u_lbc,inc_fact,u_lbc_tnd,                 &
!$OMP                v_tend,v_lbc,v_lbc_tnd,  model_levels,           &
!$OMP                w_tend,w_lbc,         w_lbc_tnd)

lbc_len = lenrim(fld_type_u,halo_type_extended)

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO i = 1, lbc_len
    u_tend(i,k) = u_lbc(i,k) + inc_fact*(u_lbc_tnd(i,k) - u_lbc(i,k))
  END DO
END DO
!$OMP END DO NOWAIT

lbc_len = lenrim(fld_type_v,halo_type_extended)

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO i = 1, lbc_len
    v_tend(i,k) = v_lbc(i,k) + inc_fact*(v_lbc_tnd(i,k) - v_lbc(i,k))
  END DO
END DO
!$OMP END DO NOWAIT

lbc_len = lenrim(fld_type_p,halo_type_extended)

!$OMP DO SCHEDULE(STATIC)
DO k = 0, model_levels
  DO i = 1, lbc_len
    w_tend(i,k) = w_lbc(i,k) + inc_fact*(w_lbc_tnd(i,k) - w_lbc(i,k))
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL


IF ( rimwidth > 0 ) THEN
! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(row_length, rows, offx, offy,             &
       model_levels, fld_type_u, u,                                     &
       lenrim(fld_type_u,halo_type_extended),                           &
       lbc_size(1,fld_type_u,halo_type_extended),                       &
       lbc_start(1,fld_type_u,halo_type_extended),                      &
       halo_i, halo_j,                                                  &
       u_tend, rimwidth, rimwidth,  rimweights, at_extremity,           &
       L_do_boundaries, L_do_halos)

! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(row_length, n_rows, offx, offy,           &
       model_levels, fld_type_v, v,                                     &
       lenrim(fld_type_v,halo_type_extended),                           &
       lbc_size(1,fld_type_v,halo_type_extended),                       &
       lbc_start(1,fld_type_v,halo_type_extended),                      &
       halo_i, halo_j,                                                  &
       v_tend, rimwidth, rimwidth,  rimweights, at_extremity,           &
       L_do_boundaries, L_do_halos)

! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(row_length, rows, offx, offy,             &
       model_levels+1, fld_type_p, w,                                   &
       lenrim(fld_type_p,halo_type_extended),                           &
       lbc_size(1,fld_type_p,halo_type_extended),                       &
       lbc_start(1,fld_type_p,halo_type_extended),                      &
       halo_i, halo_j,                                                  &
       w_tend, rimwidth, rimwidth,  rimweights, at_extremity,           &
       L_do_boundaries, L_do_halos)
END IF


! Updated u, v & w => need to change etadot

IF ( at_extremity(PSouth) ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,etadot_rim,u_at_w,v_at_w) &
!$OMP SHARED(wdims,rimwidth,etadot,ns_weights,u,v,w,intw_rho2w,         &
!$OMP intw_u2p,intw_v2p,h3_p_eta,dxi1_xi3,dxi2_xi3,h1_p_eta,h2_p_eta,   &
!$OMP deta_xi3_theta) SCHEDULE(STATIC) 
  DO k = wdims%k_start+1, wdims%k_end-1
    DO j = wdims%j_start, wdims%j_start + rimwidth-1
      DO i = wdims%i_start, wdims%i_end
        u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +         &
                                   intw_u2p(i,2)*u(i,j,k+1) ) +         &
                 intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k)   +         &
                                   intw_u2p(i,2)*u(i,j,k) )
        v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +         &
                                   intw_v2p(j,2)*v(i,j,k+1) ) +         &
                 intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +           &
                                   intw_v2p(j,2)*v(i,j,k) )
        etadot_rim = ( w(i,j,k)/h3_p_eta(i,j,k) -                       &
                       u_at_w*dxi1_xi3(i,j,k)/h1_p_eta(i,j,k) -         &
                       v_at_w*dxi2_xi3(i,j,k)/h2_p_eta(i,j,k)           &
                       )/deta_xi3_theta(i,j,k)
        etadot(i,j,k) = (1.0 - ns_weights(i,j))*etadot(i,j,k) +         &
                               ns_weights(i,j)*etadot_rim
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

IF ( at_extremity(PWest) ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,etadot_rim,u_at_w,v_at_w) &
!$OMP SHARED(wdims,rimwidth,etadot,ew_weights,u,v,w,intw_rho2w,         &
!$OMP intw_u2p,intw_v2p,h3_p_eta,h1_p_eta,h2_p_eta,dxi1_xi3,dxi2_xi3,   &
!$OMP deta_xi3_theta) SCHEDULE(STATIC) 
  DO k = wdims%k_start+1, wdims%k_end-1
    DO j = wdims%j_start, wdims%j_end
      DO i = wdims%i_start, wdims%i_start + rimwidth-1
        u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +         &
                                   intw_u2p(i,2)*u(i,j,k+1) ) +         &
                 intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k)   +         &
                                   intw_u2p(i,2)*u(i,j,k) )
        v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +         &
                                   intw_v2p(j,2)*v(i,j,k+1) ) +         &
                 intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +           &
                                   intw_v2p(j,2)*v(i,j,k) )
        etadot_rim = ( w(i,j,k)/h3_p_eta(i,j,k) -                       &
                       u_at_w*dxi1_xi3(i,j,k)/h1_p_eta(i,j,k) -         &
                       v_at_w*dxi2_xi3(i,j,k)/h2_p_eta(i,j,k)           &
                       )/deta_xi3_theta(i,j,k)
        etadot(i,j,k) =  (1.0 - ew_weights(i,j))*etadot(i,j,k) +        &
                                ew_weights(i,j)*etadot_rim
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

IF ( at_extremity(PNorth) ) THEN
!$OMP PARALLEL DO PRIVATE(ii,jj,i,j,k,etadot_rim,u_at_w,v_at_w)         &
!$OMP SHARED(wdims,rimwidth,etadot,ns_weights,u,v,w,intw_rho2w,         &
!$OMP intw_u2p,intw_v2p,h3_p_eta,h1_p_eta,h2_p_eta,dxi1_xi3,dxi2_xi3,   &
!$OMP deta_xi3_theta) SCHEDULE(STATIC) DEFAULT(NONE)
  DO k = wdims%k_start+1, wdims%k_end-1
    jj = rimwidth + 1
    DO j = wdims%j_end-rimwidth+1, wdims%j_end
      jj = jj - 1
      DO i = wdims%i_start, wdims%i_end
        u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +         &
                                   intw_u2p(i,2)*u(i,j,k+1) ) +         &
                 intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k)   +         &
                                   intw_u2p(i,2)*u(i,j,k) )
        v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +         &
                                   intw_v2p(j,2)*v(i,j,k+1) ) +         &
                 intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +           &
                                   intw_v2p(j,2)*v(i,j,k) )
        etadot_rim = ( w(i,j,k)/h3_p_eta(i,j,k) -                       &
                       u_at_w*dxi1_xi3(i,j,k)/h1_p_eta(i,j,k) -         &
                       v_at_w*dxi2_xi3(i,j,k)/h2_p_eta(i,j,k)           &
                       )/deta_xi3_theta(i,j,k)
        etadot(i,j,k) =  (1.0 - ns_weights(i,jj))*etadot(i,j,k) +       &
                                ns_weights(i,jj)*etadot_rim
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

IF ( at_extremity(PEast) ) THEN
!$OMP PARALLEL DO PRIVATE(ii,jj,i,j,k,etadot_rim,u_at_w,v_at_w)         &
!$OMP SHARED(wdims,rimwidth,etadot,ew_weights,u,v,w,intw_rho2w,         &
!$OMP intw_u2p,intw_v2p,h3_p_eta,h1_p_eta,h2_p_eta,dxi1_xi3,dxi2_xi3,   &
!$OMP deta_xi3_theta) SCHEDULE(STATIC) DEFAULT(NONE) 
  DO k = wdims%k_start+1, wdims%k_end-1
    DO j = wdims%j_start, wdims%j_end
      ii = rimwidth + 1
      DO i = wdims%i_end-rimwidth+1, wdims%i_end
        ii = ii - 1
        u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +         &
                                   intw_u2p(i,2)*u(i,j,k+1) ) +         &
                 intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k)   +         &
                                   intw_u2p(i,2)*u(i,j,k) )
        v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +         &
                                   intw_v2p(j,2)*v(i,j,k+1) ) +         &
                 intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +           &
                                   intw_v2p(j,2)*v(i,j,k) )
        etadot_rim = ( w(i,j,k)/h3_p_eta(i,j,k) -                       &
                       u_at_w*dxi1_xi3(i,j,k)/h1_p_eta(i,j,k) -         &
                       v_at_w*dxi2_xi3(i,j,k)/h2_p_eta(i,j,k)           &
                       )/deta_xi3_theta(i,j,k)
        etadot(i,j,k) =  (1.0 - ew_weights(ii,j))*etadot(i,j,k) +       &
                                ew_weights(ii,j)*etadot_rim
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! DEPENDS ON: swap_bounds
CALL swap_bounds(etadot, row_length, rows, model_levels+1,              &
                 offx, offy, fld_type_p,swap_field_is_scalar)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_lbc_dynamics

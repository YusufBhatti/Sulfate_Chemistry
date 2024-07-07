! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SETCONA_CTL (ENDGAME version)
!
! Purpose: Interface routine to SETCONA to allow SETCONA to be called
!          from INITAL. Also performs some intialisation of the
!          atmosphere STASH array that was previously carried out
!          in INITDUMP
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

SUBROUTINE setcona_ctl_4a(                                        &
  icode,cmessage)

USE jules_sea_seaice_mod, ONLY: l_ctile
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: tdims 
USE atm_d1_indices_mod, ONLY: jetatheta, jetarho, jzseak_theta, jzseak_rho, &
                              jlambda_input_p, jlambda_input_u,             &
                              jphi_input_p, jphi_input_v
USE atm_fields_mod, ONLY: rho, exner_rho_levels, exner_theta_levels,        &
                          exner_lbc, p_theta_levels, orography, orog_lbc,   &
                          orog_grad_x, orog_grad_y,orog_unfilt,             &
                          vol_smc_sat, land, frac_land, p, pstar
USE dump_headers_mod, ONLY: ih_height_gen, ih_1_c_rho_level, rh_deltaEW,    &
                            rh_deltaNS, rh_baselat, rh_baselong, rh_rotlat, &
                            rh_rotlong, rh_z_top_theta, fh_RowDepCStart,    &
                            a_fixhd, a_inthd, a_realhd, a_levdepc,          &
                            a_rowdepc, a_coldepc, sea_mask, fh_CoordSystem
USE mask_compression, ONLY: expand_from_mask
USE Control_Max_Sizes
USE UM_ParVars
USE atm_land_sea_mask, ONLY: atmos_landmask_local
USE lbc_mod

USE lbc_read_data_mod, ONLY: rimweightsa
USE nlsizes_namelist_mod, ONLY:                                             &
    a_len2_coldepc, a_len2_rowdepc, bl_levels, cloud_levels,                &
    global_row_length, global_rows, land_field, model_levels,               &
    n_rows, ozone_levels, row_length, rows,st_levels, theta_field_size,     &
    tr_levels, tr_ukca, tr_vars
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


INTEGER ::                                                        &
  icode      ! Return code

CHARACTER(LEN=errormessagelength) ::                              &
  cmessage   ! Error message

! Local variables

INTEGER ::                                                        &
  i,j,ij                                                          &
             ! loop counters
, dummy      ! Dummy integer

REAL ::                                                           &
  fland(theta_field_size)
                           ! Land fraction un-compressed

LOGICAL ::                                                        &
  lsea                     ! True if any of gridbox is sea
LOGICAL :: log_val         ! used for transfers
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETCONA_CTL_4A'

!----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! 1.0 Call SETCONA

! DEPENDS ON: setcona_4A
CALL setcona_4a( a_levdepc(jetatheta), a_levdepc(jetarho),        &
                 a_levdepc(Jzseak_theta), a_levdepc(Jzseak_rho),  &
              vol_smc_sat, land, orography,                       &
              orog_grad_x,orog_grad_y, orog_unfilt,               &
              rho, exner_rho_levels,                              &
              orog_lbc,exner_lbc,                                 &
              lenrima,lbc_sizea,lbc_starta,                       &
              rimwidtha, rimweightsa,                             &
              global_row_length, global_rows,                     &
              rows, n_rows, row_length,                           &
              land_field, bl_levels,                              &
              a_inthd(ih_1_c_rho_level), cloud_levels,            &
              a_realhd(rh_z_top_theta),                           &
              ozone_levels,bl_levels,tr_levels, tr_vars, tr_ukca, &
              a_inthd(ih_height_gen),                             &
              a_realhd(rh_deltaEW), a_realhd(rh_deltaNS),         &
              a_realhd(rh_baselat), a_realhd(rh_baselong),        &
              a_realhd(rh_rotlat), a_realhd(rh_rotlong),          &
              a_coldepc(jlambda_input_p),                         &
              a_coldepc(jlambda_input_u),                         &
              a_rowdepc(jphi_input_p),                            &
              a_rowdepc(jphi_input_v),                            &
              a_fixhd(fh_CoordSystem),                            &
              a_fixhd(fh_RowDepCStart),                           &
              exner_theta_levels,                                 &
              p_theta_levels,                                     &
              p, pstar,                                           &
              a_len2_coldepc, a_len2_rowdepc,                     &
              icode, cmessage)

IF (l_ctile) THEN
  CALL expand_from_mask(fland,frac_land,                     &
      atmos_landmask_local, theta_field_size, dummy)
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      ij= (i-1) + (j-1)*tdims%i_len
      lsea = (fland(ij+1) /= 1.0)
      sea_mask(i,j) = lsea
    END DO
  END DO
ELSE
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      lsea = (.NOT. TRANSFER(land(i,j),log_val))
      sea_mask(i,j) = lsea
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE setcona_ctl_4a

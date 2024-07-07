! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to compute Leonard term vertical fluxes of meridional wind

MODULE leonard_term_vert_v_mod

IMPLICIT NONE

!
! Description:
!   Calculates Leonard term vertical fluxes and increments for
!   the meridional wind v.  The increments are passed out.
!   Also the increments and fluxes maybe output as diagnostics by
!   calls to STASH routines.
!
! Method:
!   The flux of v is computed at points above and below the v-points.
!   The Leonard term flux is proportional to the correlation of the
!   horizontal gradient of v with the horizontal gradient of w
!   (see Moeng et al 2010).
!   The increment is then computed by differencing the flux in the
!   vertical (with appropriate scaling by rho * r^2).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Diffusion and filtering
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName =                    &
                                        'LEONARD_TERM_VERT_V_MOD'


CONTAINS


SUBROUTINE leonard_term_vert_v( v, w, rho_wet_rsq_tq,                   &
                                dtrdz_v, v_inc_leonard, stashwork3 )

USE atm_fields_bounds_mod, ONLY: tdims_s, vdims, vdims_s, wdims_s
USE nlsizes_namelist_mod,  ONLY: bl_levels
USE level_heights_mod,     ONLY: r_at_v, r_at_v_w
USE timestep_mod,          ONLY: timestep
USE turb_diff_mod,         ONLY: leonard_kl
USE p_to_v_mod,            ONLY: p_to_v

USE ereport_mod,           ONLY: ereport
USE errormessagelength_mod,ONLY: errormessagelength

USE stash_array_mod,       ONLY: len_stlist, stindex, stlist,           &
                                 num_stash_levels, stash_levels, si, sf
USE submodel_mod,          ONLY: atmos_im
USE um_parvars,            ONLY: at_extremity
USE um_stashcode_mod,      ONLY: stashcode_bl_sec

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Meridional wind field for which the vertical flux is to be computed
REAL, INTENT(IN) :: v (     vdims_s%i_start:vdims_s%i_end,              &
                            vdims_s%j_start:vdims_s%j_end,              &
                            vdims_s%k_start:vdims_s%k_end )

! Vertical velocity used to compute the flux of v
REAL, INTENT(IN) :: w     ( wdims_s%i_start:wdims_s%i_end,              &
                            wdims_s%j_start:wdims_s%j_end,              &
                            bl_levels )

! Wet (total) density * r^2 on theta-levels, with halos
REAL, INTENT(IN) :: rho_wet_rsq_tq ( tdims_s%i_start:tdims_s%i_end,     &
                                     tdims_s%j_start:tdims_s%j_end,     &
                                     bl_levels )

! Array of timestep / ( dz * rho * r^2 ), precalculated in the BL scheme
REAL, INTENT(IN) :: dtrdz_v        ( vdims%i_start:vdims%i_end,         &
                                     vdims%j_start:vdims%j_end,         &
                                     bl_levels )

! Output increment to u (gets added to latest field array at higher level)
REAL, INTENT(OUT) :: v_inc_leonard ( vdims%i_start:vdims%i_end,         &
                                     vdims%j_start:vdims%j_end,         &
                                     bl_levels )

! STASH workspace
REAL, INTENT(INOUT) :: stashwork3(*)


! Array to store vertical flux of u at theta-levels above and below u
REAL :: flux            ( vdims%i_start:vdims%i_end,                    &
                          vdims%j_start:vdims%j_end,                    &
                          0:bl_levels )

! density * r^2 horizontally interpolated to theta-levels above and below v
REAL :: rho_wet_rsq_vth ( vdims%i_start:vdims%i_end,                    &
                          vdims%j_start:vdims%j_end,                    &
                          bl_levels )

! Leonard term parameter at theta-levels above and below v
REAL :: kl_v            ( vdims%i_start:vdims%i_end,                    &
                          vdims%j_start:vdims%j_end,                    &
                          bl_levels )


! Loop counters
INTEGER :: i, j, k

! Stuff for STASH calls
INTEGER :: icode, item
INTEGER, PARAMETER :: im_index = 1
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LEONARD_TERM_VERT_V'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! ----------------------------------------------------------------------
! Horizontal interpolation of density*r^2 to above and below v-points.
! ----------------------------------------------------------------------

CALL p_to_v ( rho_wet_rsq_tq,                                           &
              tdims_s%i_start,tdims_s%i_end,                            &
              tdims_s%j_start,tdims_s%j_end,                            &
              vdims%i_start,vdims%i_end,                                &
              vdims%j_start,vdims%j_end,                                &
              1,bl_levels,rho_wet_rsq_vth )

! ----------------------------------------------------------------------
! Calculate kl at theta-levels above and below v points,
! accounting for stability limit
! ----------------------------------------------------------------------

DO k = 1, bl_levels
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end

      ! Leonard term parameter is the min of the input leonard_kl
      ! and the max stable value   6 dz / ( dt dw )

      ! For dw we use the maximum horizontal finite difference that
      ! contributes to the flux at each point.

      kl_v(i,j,k) = MIN( leonard_kl,                                    &
          6.0 * ( r_at_v(i,j,k+1) - r_at_v(i,j,k) )                     &
              / ( timestep * MAX(                                       &
                    ! Difference either side (point lies exactly between these):
                    ABS( w(i,j+1,k)   - w(i,j,k)     ),                 &
                    ! Neighbouring differences to the west:
                    ABS( w(i,j,k)     - w(i-1,j,k)   ),                 &
                    ABS( w(i,j+1,k)   - w(i-1,j+1,k) ),                 &
                    ! Neighbouring differences to the east:
                    ABS( w(i+1,j,k)   - w(i,j,k)   ),                   &
                    ABS( w(i+1,j+1,k) - w(i,j+1,k) )                    &
                                )                                       &
                )      )

    END DO
  END DO
END DO


! ----------------------------------------------------------------------
! Calculate flux on theta levels above and below v points
! ----------------------------------------------------------------------

DO k = 1, bl_levels-1 
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      ! Leonard term sub-grid vertical flux is proportional to the
      ! product of the grid-scale horizontal differences in v and w:
      flux(i,j,k) = (kl_v(i,j,k) / 12.0)                                &
      ! 8 terms contribute to each of x and y direction, so scale by 1/8
                  * (1.0 / 8.0) *(                                      &
      !
      ! Terms from gradient in y-direction...
          ! ( point is half-way between 2 w-points in y-direction,
          !   so only one dw_y contributes ):
          2.0 * ( (v(i,j+1,k)  -v(i,j-1,k)  )                           &
                + (v(i,j+1,k+1)-v(i,j-1,k+1)) ) * (w(i,j+1,k)-w(i,j,k)) &
      !
      ! Terms from gradient in x-direction...
          ! from difference to the west:
        + ( (v(i,j,k)  -v(i-1,j,k)  )                                   &
          + (v(i,j,k+1)-v(i-1,j,k+1)) ) * ( (w(i,j,k)  -w(i-1,j,k)  )   &
                                          + (w(i,j+1,k)-w(i-1,j+1,k)) ) &
          ! from difference to the east:
        + ( (v(i+1,j,k)-v(i,j,k)    )                                   &
          + (v(i+1,j,k+1)-v(i,j,k+1)) ) * ( (w(i+1,j,k)-w(i,j,k)    )   &
                                          + (w(i+1,j+1,k)-w(i,j+1,k)) ) &
      !
                                 )                                      &
                       * rho_wet_rsq_vth(i,j,k)
      ! Flux now includes density * r^2 factor
    END DO
  END DO
END DO

! Set flux to zero at top and bottom
DO j = vdims%j_start, vdims%j_end
  DO i = vdims%i_start, vdims%i_end
    flux(i,j,0) = 0.0
    flux(i,j,bl_levels) = 0.0
  END DO
END DO


! ----------------------------------------------------------------------
! Difference the flux in the vertical to get increment on rho levels
! ----------------------------------------------------------------------

DO k = 1, bl_levels
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      v_inc_leonard(i,j,k) = -( flux(i,j,k) - flux(i,j,k-1) )           &
                           * dtrdz_v(i,j,k)
      ! Note: both flux and dtrdz_v contain the rho * r^2 factor
    END DO
  END DO
END DO


! ----------------------------------------------------------------------
! Output increment and flux diagnostics to STASH
! ----------------------------------------------------------------------

icode=0

! Increment diagnostic
item = 196
IF (icode <= 0 .AND. sf(item,stashcode_bl_sec)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork3(si(item,stashcode_bl_sec,im_index)),     &
                    v_inc_leonard,                                      &
                    vdims%i_len, vdims%j_len, bl_levels,                &
                    0, 0, 0, 0, at_extremity,                           &
                    stlist(1,stindex(1,item,stashcode_bl_sec,im_index)),&
                    len_stlist, stash_levels, num_stash_levels+1,       &
                    atmos_im, stashcode_bl_sec, item,                   &
                    icode, cmessage )
  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item)"//TRIM(cmessage)
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF

! Flux diagnostic
item = 554
IF (icode <= 0 .AND. sf(item,stashcode_bl_sec)) THEN

  ! Remove factor of r^2 from the flux for diagnostic output
  DO k = 1, bl_levels-1
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        flux(i,j,k) = flux(i,j,k) / ( r_at_v_w(i,j,k)*r_at_v_w(i,j,k) )
      END DO
    END DO
  END DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork3(si(item,stashcode_bl_sec,im_index)),     &
                    flux(:,:,1:bl_levels),                              &
                    vdims%i_len, vdims%j_len, bl_levels,                &
                    0, 0, 0, 0, at_extremity,                           &
                    stlist(1,stindex(1,item,stashcode_bl_sec,im_index)),&
                    len_stlist, stash_levels, num_stash_levels+1,       &
                    atmos_im, stashcode_bl_sec, item,                   &
                    icode, cmessage )
  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item)"//TRIM(cmessage)
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF

! Limited Kl diagnostic
item = 551
IF (icode <= 0 .AND. sf(item,stashcode_bl_sec)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork3(si(item,stashcode_bl_sec,im_index)),     &
                    kl_v,                                               &
                    vdims%i_len, vdims%j_len, bl_levels,                &
                    0, 0, 0, 0, at_extremity,                           &
                    stlist(1,stindex(1,item,stashcode_bl_sec,im_index)),&
                    len_stlist, stash_levels, num_stash_levels+1,       &
                    atmos_im, stashcode_bl_sec, item,                   &
                    icode, cmessage )
  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item)"//TRIM(cmessage)
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE Leonard_term_vert_v

END MODULE leonard_term_vert_v_mod

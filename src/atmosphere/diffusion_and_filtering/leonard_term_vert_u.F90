! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to compute Leonard term vertical fluxes of zonal wind

MODULE leonard_term_vert_u_mod

IMPLICIT NONE

!
! Description:
!   Calculates Leonard term vertical fluxes and increments for
!   the zonal wind u.  The increments are passed out.
!   Also the increments and fluxes maybe output as diagnostics by
!   calls to STASH routines.
!
! Method:
!   The flux of u is computed at points above and below the u-points.
!   The Leonard term flux is proportional to the correlation of the
!   horizontal gradient of u with the horizontal gradient of w
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
                                        'LEONARD_TERM_VERT_U_MOD'


CONTAINS


SUBROUTINE leonard_term_vert_u( u, w, rho_wet_rsq_tq,                   &
                                dtrdz_u, u_inc_leonard, stashwork3 )

USE atm_fields_bounds_mod, ONLY: tdims_s, udims, udims_s, wdims_s
USE nlsizes_namelist_mod,  ONLY: bl_levels
USE level_heights_mod,     ONLY: r_at_u, r_at_u_w 
USE timestep_mod,          ONLY: timestep
USE turb_diff_mod,         ONLY: leonard_kl
USE p_to_u_mod,            ONLY: p_to_u

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

! Zonal wind field for which the vertical flux is to be computed
REAL, INTENT(IN) :: u (     udims_s%i_start:udims_s%i_end,              &
                            udims_s%j_start:udims_s%j_end,              &
                            udims_s%k_start:udims_s%k_end )

! Vertical velocity used to compute the flux of u
REAL, INTENT(IN) :: w     ( wdims_s%i_start:wdims_s%i_end,              &
                            wdims_s%j_start:wdims_s%j_end,              &
                            bl_levels )

! Wet (total) density * r^2 on theta-levels, with halos
REAL, INTENT(IN) :: rho_wet_rsq_tq ( tdims_s%i_start:tdims_s%i_end,     &
                                     tdims_s%j_start:tdims_s%j_end,     &
                                     bl_levels )

! Array of timestep / ( dz * rho * r^2 ), precalculated in the BL scheme
REAL, INTENT(IN) :: dtrdz_u        ( udims%i_start:udims%i_end,         &
                                     udims%j_start:udims%j_end,         &
                                     bl_levels )

! Output increment to u (gets added to latest field array at higher level)
REAL, INTENT(OUT) :: u_inc_leonard ( udims%i_start:udims%i_end,         &
                                     udims%j_start:udims%j_end,         &
                                     bl_levels )

! STASH workspace
REAL, INTENT(INOUT) :: stashwork3(*)


! Array to store vertical flux of u at theta-levels above and below u
REAL :: flux            ( udims%i_start:udims%i_end,                    &
                          udims%j_start:udims%j_end,                    &
                          0:bl_levels )

! density * r^2 horizontally interpolated to theta-levels above and below u
REAL :: rho_wet_rsq_uth ( udims%i_start:udims%i_end,                    &
                          udims%j_start:udims%j_end,                    &
                          bl_levels )

! Leonard term parameter at theta-levels above and below u
REAL :: kl_u            ( udims%i_start:udims%i_end,                    &
                          udims%j_start:udims%j_end,                    &
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

CHARACTER(LEN=*), PARAMETER :: RoutineName='LEONARD_TERM_VERT_U'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! ----------------------------------------------------------------------
! Horizontal interpolation of density*r^2 to above and below u-points.
! ----------------------------------------------------------------------

CALL p_to_u ( rho_wet_rsq_tq,                                           &
              tdims_s%i_start,tdims_s%i_end,                            &
              tdims_s%j_start,tdims_s%j_end,                            &
              udims%i_start,udims%i_end,                                &
              udims%j_start,udims%j_end,                                &
              1,bl_levels,rho_wet_rsq_uth )

! ----------------------------------------------------------------------
! Calculate kl at theta-levels above and below u points,
! accounting for stability limit
! ----------------------------------------------------------------------

DO k = 1, bl_levels
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end

      ! Leonard term parameter is the min of the input leonard_kl
      ! and the max stable value   6 dz / ( dt dw )

      ! For dw we use the maximum horizontal finite difference that
      ! contributes to the flux at each point.

      kl_u(i,j,k) = MIN( leonard_kl,                                    &
          6.0 * ( r_at_u(i,j,k+1) - r_at_u(i,j,k) )                     &
              / ( timestep * MAX(                                       &
                    ! Difference either side (point lies exactly between these):
                    ABS( w(i+1,j,k)   - w(i,j,k)     ),                 &
                    ! Neighbouring differences to the south:
                    ABS( w(i,j,k)     - w(i,j-1,k)   ),                 &
                    ABS( w(i+1,j,k)   - w(i+1,j-1,k) ),                 &
                    ! Neighbouring differences to the north:
                    ABS( w(i,j+1,k)   - w(i,j,k)   ),                   &
                    ABS( w(i+1,j+1,k) - w(i+1,j,k) )                    &
                                )                                       &
                )      )

    END DO
  END DO
END DO


! ----------------------------------------------------------------------
! Calculate flux on theta levels above and below u points
! ----------------------------------------------------------------------

DO k = 1, bl_levels-1 
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      ! Leonard term sub-grid vertical flux is proportional to the
      ! product of the grid-scale horizontal differences in u and w:
      flux(i,j,k) = (kl_u(i,j,k) / 12.0)                                &
      ! 8 terms contribute to each of x and y direction, so scale by 1/8
                  * (1.0 / 8.0) *(                                      &
      !
      ! Terms from gradient in x-direction...
          ! ( point is half-way between 2 w-points in x-direction,
          !   so only one dw_x contributes ):
          2.0 * ( (u(i+1,j,k)  -u(i-1,j,k)  )                           &
                + (u(i+1,j,k+1)-u(i-1,j,k+1)) ) * (w(i+1,j,k)-w(i,j,k)) &
      !
      ! Terms from gradient in y-direction...
          ! from difference to the south:
        + ( (u(i,j,k)  -u(i,j-1,k)  )                                   &
          + (u(i,j,k+1)-u(i,j-1,k+1)) ) * ( (w(i,j,k)  -w(i,j-1,k)  )   &
                                          + (w(i+1,j,k)-w(i+1,j-1,k)) ) &
          ! from difference to the north:
        + ( (u(i,j+1,k)-u(i,j,k)    )                                   &
          + (u(i,j+1,k+1)-u(i,j,k+1)) ) * ( (w(i,j+1,k)-w(i,j,k)    )   &
                                          + (w(i+1,j+1,k)-w(i+1,j,k)) ) &
      !
                                 )                                      &
                       * rho_wet_rsq_uth(i,j,k)
      ! Flux now includes density * r^2 factor
    END DO
  END DO
END DO

! Set flux to zero at top and bottom
DO j = udims%j_start, udims%j_end
  DO i = udims%i_start, udims%i_end
    flux(i,j,0) = 0.0
    flux(i,j,bl_levels) = 0.0
  END DO
END DO


! ----------------------------------------------------------------------
! Difference the flux in the vertical to get increment on rho levels
! ----------------------------------------------------------------------

DO k = 1, bl_levels
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      u_inc_leonard(i,j,k) = -( flux(i,j,k) - flux(i,j,k-1) )           &
                           * dtrdz_u(i,j,k)
      ! Note: both flux and dtrdz_u contain the rho * r^2 factor
    END DO
  END DO
END DO


! ----------------------------------------------------------------------
! Output increment and flux diagnostics to STASH
! ----------------------------------------------------------------------

icode=0

! Increment diagnostic
item = 195
IF (icode <= 0 .AND. sf(item,stashcode_bl_sec)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork3(si(item,stashcode_bl_sec,im_index)),     &
                    u_inc_leonard,                                      &
                    udims%i_len, udims%j_len, bl_levels,                &
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
item = 553
IF (icode <= 0 .AND. sf(item,stashcode_bl_sec)) THEN

  ! Remove factor of r^2 from the flux for diagnostic output
  DO k = 1, bl_levels-1
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        flux(i,j,k) = flux(i,j,k) / ( r_at_u_w(i,j,k)*r_at_u_w(i,j,k) )
      END DO
    END DO
  END DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork3(si(item,stashcode_bl_sec,im_index)),     &
                    flux(:,:,1:bl_levels),                              &
                    udims%i_len, udims%j_len, bl_levels,                &
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
item = 550
IF (icode <= 0 .AND. sf(item,stashcode_bl_sec)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork3(si(item,stashcode_bl_sec,im_index)),     &
                    kl_u,                                               &
                    udims%i_len, udims%j_len, bl_levels,                &
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

END SUBROUTINE Leonard_term_vert_u

END MODULE leonard_term_vert_u_mod

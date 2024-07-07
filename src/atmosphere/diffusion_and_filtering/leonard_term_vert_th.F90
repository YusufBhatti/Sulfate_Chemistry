! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to compute Leonard term vertical fluxes of fields at theta-points

MODULE leonard_term_vert_th_mod

IMPLICIT NONE

!
! Description:
!   Calculates Leonard term vertical fluxes and increments for any
!   field defined at theta-points.  The increments are passed out.
!   Also the increments and fluxes maybe output as diagnostics by
!   calls to STASH routines.
!
! Method:
!   The flux of the input field is computed at rho-points (above and
!   below the theta points where the field is defined).  The Leonard
!   term flux is proportional to the correlation of the horizontal
!   gradient of the field with the horizontal gradient of w
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
                                        'LEONARD_TERM_VERT_TH_MOD'


CONTAINS


SUBROUTINE leonard_term_vert_th( field, w, kl,                          &
                                 rho_wet_rsq, dtrdz_charney_grid,       &
                                 field_inc_leonard, exner_theta_levels, &
                                 exner_rho_levels,                      &
                                 item_inc, item_flx, stashwork3 )

USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims, tdims_s, wdims_s
USE level_heights_mod,     ONLY: r_rho_levels
USE nlsizes_namelist_mod,  ONLY: bl_levels
USE planet_constants_mod,  ONLY: cp

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

! Field for which the vertical flux is to be computed
REAL, INTENT(IN) :: field ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            bl_levels )

! Vertical velocity used to compute the flux of field
REAL, INTENT(IN) :: w     ( wdims_s%i_start:wdims_s%i_end,              &
                            wdims_s%j_start:wdims_s%j_end,              &
                            bl_levels )

! Leonard term parameter used to scale the fluxes
! (limiting for numerical stability already aplied at a higher level)
REAL, INTENT(IN) :: kl    ( tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            1:bl_levels )

! Wet (total) density * r^2 on rho-levels
REAL, INTENT(IN) :: rho_wet_rsq        ( tdims_s%i_start:tdims_s%i_end, &
                                         tdims_s%j_start:tdims_s%j_end, &
                                         bl_levels )
! exner pressure on theta and rho levels
REAL, INTENT(IN) :: exner_theta_levels ( tdims_s%i_start:tdims_s%i_end, &
                                         tdims_s%j_start:tdims_s%j_end, &
                                         tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(IN) :: exner_rho_levels   ( pdims_s%i_start:pdims_s%i_end, &
                                         pdims_s%j_start:pdims_s%j_end, &
                                         pdims_s%k_start:pdims_s%k_end + 1)

! Array of timestep / ( dz * rho * r^2 ), precalculated in the BL scheme
REAL, INTENT(IN) :: dtrdz_charney_grid ( pdims%i_start:pdims%i_end,     &
                                         pdims%j_start:pdims%j_end,     &
                                         bl_levels )

! Output increment to field (gets added to latest field array at higher level)
! Declare the increment array with the TARGET attribute, 
! so that we can point a pointer at it for diagnostics
REAL, TARGET, INTENT(OUT) :: field_inc_leonard ( tdims%i_start:tdims%i_end,  &
                                                 tdims%j_start:tdims%j_end,  &
                                                 bl_levels )

! Declare a pointer for the increment diagnostic actually output.  
! For qt and w increments, this just points to the main increment array,
! but for theta, we want to scale by exner to get a temperature increment.
REAL, POINTER :: field_inc_diag(:,:,:) => NULL()

! STASH item numbers for the increment and flux diagnostics for current field
INTEGER, INTENT(IN) :: item_inc
INTEGER, INTENT(IN) :: item_flx

! STASH workspace
REAL, INTENT(INOUT) :: stashwork3(*)


! Array to store vertical flux of field on rho-levels
REAL :: flux ( tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               bl_levels+1 )


! Loop counters
INTEGER :: i, j, k

! Stuff for STASH calls
INTEGER :: icode
INTEGER, PARAMETER :: im_index = 1
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LEONARD_TERM_VERT_TH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! ----------------------------------------------------------------------
! Calculate vertical flux on rho levels
! ----------------------------------------------------------------------

DO k = 2, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ! Leonard term sub-grid vertical flux is proportional to the
      ! product of the grid-scale horizontal differences in field and w:
      flux(i,j,k) = (Kl(i,j,k) / 12.0)                                  &
      ! 4 terms contribute to each of x and y direction, so scale by 1/4
                  * (1.0 / 4.0) *(                                      &
      !
      ! Terms from gradient in x-direction...
          ! from theta-level below:
          (field(i,j,k-1)-field(i-1,j,k-1)) * (w(i,j,k-1)-w(i-1,j,k-1)) &
        + (field(i+1,j,k-1)-field(i,j,k-1)) * (w(i+1,j,k-1)-w(i,j,k-1)) &
          ! from theta-level above:
        + (field(i,j,k)  -field(i-1,j,k)  ) * (w(i,j,k)  -w(i-1,j,k)  ) &
        + (field(i+1,j,k)-field(i,j,k)    ) * (w(i+1,j,k)-w(i,j,k)    ) &
      !
      ! Terms from gradient in y-direction...
          ! from theta-level below:
        + (field(i,j,k-1)-field(i,j-1,k-1)) * (w(i,j,k-1)-w(i,j-1,k-1)) &
        + (field(i,j+1,k-1)-field(i,j,k-1)) * (w(i,j+1,k-1)-w(i,j,k-1)) &
          ! from theta-level above:
        + (field(i,j,k)  -field(i,j-1,k)  ) * (w(i,j,k)  -w(i,j-1,k)  ) &
        + (field(i,j+1,k)-field(i,j,k)    ) * (w(i,j+1,k)-w(i,j,k)    ) &
      !
                                 )                                      &
                       * rho_wet_rsq(i,j,k)
      ! Flux now includes density * r^2 factor
    END DO
  END DO
END DO

! Set flux to zero at top and bottom
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    flux(i,j,1) = 0.0
    flux(i,j,bl_levels+1) = 0.0
  END DO
END DO


! ----------------------------------------------------------------------
! Difference the flux in the vertical to get increment on theta-levels
! ----------------------------------------------------------------------

DO k = 1, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      field_inc_leonard(i,j,k) = -( flux(i,j,k+1) - flux(i,j,k) )       &
                               * dtrdz_charney_grid(i,j,k)
      ! Note: both flux and dtrdz_charney_grid contain the rho * r^2 factor
    END DO
  END DO
END DO

! ----------------------------------------------------------------------
! Output increment and flux diagnostics to STASH
! Multiply theta increment by exner_theta_levels to get an increment in
! temperature rather than theta
! Multiply theta flux by cp * exner_rho_levels to get heat flux in W/m2
! ----------------------------------------------------------------------

icode=0

! If increment diagnostic requested...
IF (icode <= 0 .AND. sf(item_inc,stashcode_bl_sec)) THEN

  ! If this is the theta increment, we need to convert to temperature increment
  IF (item_inc == 197) THEN

    ! Allocate separate space for the diag for converting to T increment
    ALLOCATE( field_inc_diag( tdims%i_start:tdims%i_end,               &
                              tdims%j_start:tdims%j_end,               &
                              bl_levels ) )
    ! Convert to T increment
    DO k = 1, bl_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          field_inc_diag(i,j,k) = field_inc_leonard(i,j,k) *           &
                                  exner_theta_levels(i,j,k)
        END DO
      END DO
    END DO

  ! If this is not theta, we can just output the increment field as-is
  ELSE

    ! Diag field is just a pointer to the main increment field
    field_inc_diag => field_inc_leonard

  END IF

  ! Output the increment diagnostic to STASH
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork3(si(item_inc,stashcode_bl_sec,im_index)),      &
                    field_inc_diag,                                          &
                    tdims%i_end, tdims%j_end, bl_levels,                     &
                    0, 0, 0, 0, at_extremity,                                &
                    stlist(1,stindex(1,item_inc,stashcode_bl_sec,im_index)), &
                    len_stlist, stash_levels, num_stash_levels+1,            &
                    atmos_im, stashcode_bl_sec, item_inc,                    &
                    icode, cmessage )
  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item)"//TRIM(cmessage)
    CALL ereport(RoutineName,icode,cmessage)
  END IF

  ! Deallocate / nullify the temporary diagnostic array
  IF (item_inc == 197)  DEALLOCATE(field_inc_diag)
  field_inc_diag => NULL()

END IF

! If flux diagnostic requested...
IF (icode <= 0 .AND. sf(item_flx,stashcode_bl_sec)) THEN

  ! Remove factor of r^2 from the flux for diagnostic output
  DO k = 2, bl_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        flux(i,j,k) = flux(i,j,k) / ( r_rho_levels(i,j,k)*r_rho_levels(i,j,k) )
      END DO
    END DO
  END DO

  ! If this is the flux of theta_l, convert to a heat flux in Wm-2
  IF ( item_flx == 556 ) THEN
    DO k = 2, bl_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          flux(i,j,k) = flux(i,j,k) * cp * exner_rho_levels(i,j,k)
        END DO
      END DO
    END DO
  END IF

  ! Output the flux diagnostic to STASH
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork3(si(item_flx,stashcode_bl_sec,im_index)),      &
                    flux(:,:,1:bl_levels),                                   &
                    tdims%i_end, tdims%j_end, bl_levels,                     &
                    0, 0, 0, 0, at_extremity,                                &
                    stlist(1,stindex(1,item_flx,stashcode_bl_sec,im_index)), &
                    len_stlist, stash_levels, num_stash_levels+1,            &
                    atmos_im, stashcode_bl_sec, item_flx,                    &
                    icode, cmessage )
  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item)"//TRIM(cmessage)
    CALL ereport(RoutineName,icode,cmessage)
  END IF

END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE leonard_term_vert_th

END MODULE leonard_term_vert_th_mod

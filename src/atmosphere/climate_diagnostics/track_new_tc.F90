! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
! This routine average the vorticity of 'nlevs_avg' levels,
! requested in winds_for_track.F90 and send them to
! the Filter_2D to do the filtering to ntop_tc.
!
!     Code Owner: Please refer to the UM file CodeOwners.txt
!     This file belongs in section: Climate Diagnostics
MODULE track_new_TC_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TRACK_NEW_TC_MOD'

CONTAINS

SUBROUTINE track_new_TC(nbot, ntop, delta_lambda, delta_phi, tc_activate,  &
                        vort, vort_T)

! Use array dimensions
USE atm_fields_bounds_mod,   ONLY: udims, vdims, udims_s,                  &
                                   vdims_s, array_dims
! Load logicals and arrays for 'new-TC' method
USE track_mod,               ONLY: nlevs_avg, u_tc_new, v_tc_new
! planet radius to compute vorticity
USE planet_constants_mod,    ONLY: planet_radius
! Model grid trigonometry
USE trignometric_mod,        ONLY: cos_v_latitude
! Use swap bounds
USE halo_exchange,           ONLY: swap_bounds
USE mpp_conf_mod,            ONLY: swap_field_is_scalar
! get Global dimension
USE nlsizes_namelist_mod,    ONLY: global_row_length
! Use FV-TRACK 2D filter
USE filter_2D_mod,           ONLY: filter_2D

! UM settings for swap and handling at poles
USE UM_ParVars,              ONLY: at_extremity
USE UM_ParParams,            ONLY: psouth, pnorth
USE Field_Types,             ONLY: fld_type_v
! UM error and printout variables
USE umPrintMgr,              ONLY: umprint,ummessage
USE ereport_mod,             ONLY: ereport
USE errormessagelength_mod,  ONLY: errormessagelength

! DrHook variables
USE yomhook,                 ONLY: lhook, dr_hook
USE parkind1,                ONLY: jprb, jpim

IMPLICIT NONE
!  Argument variables
INTEGER, INTENT(IN) :: nbot         ! Lower truncation wavenumber
INTEGER, INTENT(IN) :: ntop         ! Upper truncation wavenumber
REAL,    INTENT(IN) :: delta_lambda ! grid longitude spacing in radians
REAL,    INTENT(IN) :: delta_phi    ! grid latitude  spacing in radians
LOGICAL, INTENT(IN) :: tc_activate  ! Activates TC-tracking (uses different
                                    ! settings for filtering)

REAL, INTENT(INOUT) ::                                                     &
    vort(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end)              &
   ! Vorticity field
,   vort_T(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end)
   ! Filtered Vorticity

! Local variables
INTEGER :: i,j,k       ! indexing types
INTEGER, SAVE :: nlim  ! total number of longitude points div. by 2

REAL ::                                                                    &
    u_s(udims_s%i_start:udims_s%i_end,                                     &
        vdims_s%j_start:vdims_s%j_end, nlevs_avg)                          &
   ! U with haloes to compute vorticity
 ,  v_s(udims_s%i_start:udims_s%i_end,                                     &
        vdims_s%j_start:vdims_s%j_end, nlevs_avg)                          &
   ! V with haloes to compute vorticity

 ,  vort_all_levs(udims%i_start:udims%i_end,                               &
                  vdims%j_start:vdims%j_end,                               &
                  nlevs_avg)
   ! Vorticity in all levels requested


LOGICAL, SAVE :: declare = .TRUE. ! to save arrays

!Error message handlers
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER  :: RoutineName='TRACK_NEW_TC'

! DrHook handlers
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up vor dimensions
IF (declare) THEN

  ! Initialize variables from UM data for the sph.harm calculations
  ! For the UM, not being an spectral model, NLIM:
  IF (MOD(global_row_length,2) == 0) THEN
    nlim=global_row_length/2
  ELSE
    nlim=(global_row_length+1)/2
  END IF

  ! Stop model  if ntop >= nlim
  IF (ntop >= nlim ) THEN
    WRITE(umMessage,'(A)')'  **********************************  '
    CALL umPrint(umMessage,src='track_new_TC')
    WRITE(umMessage,'(A)')' NTOP is higher than Truncation N'
    CALL umPrint(umMessage,src='track_new_TC')
    WRITE(umMessage,'(A)')'(Please choose a lower truncation number' //&
                        'in sect 30 panel)'
    CALL umPrint(umMessage,src='track_new_TC')
    icode   = 1
    WRITE (cmessage,'(A)')'Ntop >= Nlim, choose a lower truncation'
    CALL ereport(routinename, icode, cmessage)
  END IF

  ! switch off the definition of arrays
  declare=.FALSE.
END IF


! Fill u and v haloes with u_850 and v_850
! for u (loop over levels too)

DO k=1,nlevs_avg
  DO j = vdims%j_start, vdims%j_end
    DO i = udims%i_start, udims%i_end
      u_s(i,j,k)=u_tc_new(i,j,k)
    END DO
  END DO

  ! for v
  DO j = vdims%j_start, vdims%j_end
    DO i = udims%i_start, udims%i_end
      v_s(i,j,k)=v_tc_new(i,j,k)
    END DO
  END DO
END DO

IF (ALLOCATED(u_tc_new)) DEALLOCATE(u_tc_new)
IF (ALLOCATED(v_tc_new)) DEALLOCATE(v_tc_new)

! swap bounds for u_s and v_s
CALL swap_bounds(u_s,                                                    &
                   udims_s%i_len - 2*udims_s%halo_i,                     &
                   vdims_s%j_len - 2*vdims_s%halo_j,                     &
                   nlevs_avg, udims_s%halo_i, vdims_s%halo_j,            &
                   fld_type_v,swap_field_is_scalar)

CALL swap_bounds(v_s,                                                    &
                   udims_s%i_len - 2*udims_s%halo_i,                     &
                   vdims_s%j_len - 2*vdims_s%halo_j,                     &
                   nlevs_avg, udims_s%halo_i, vdims_s%halo_j,            &
                   fld_type_v,swap_field_is_scalar)

! Compute vorticity

! On a B grid, each U,V is located in each corner so in order to compute
! vorticity we take the left and right value -> Thus we divide by 2 delta_X
!
!     v(i-1) ....... v(i) ...... v(i+j)

! loop over levels
DO k=1,nlevs_avg

  DO j = vdims%j_start, vdims%j_end
    DO i = udims%i_start, udims%i_end
      vort_all_levs(i,j,k) = (v_s(i+1,j,k) - v_s(i-1,j,k))/             &
                             (2.0*delta_lambda) -                       &
                             (u_s(i,j+1,k) * cos_v_latitude(i,j+1) -    &
                              u_s(i,j-1,k) * cos_v_latitude(i,j-1))     &
                             /(2.0*delta_phi)
      ! Compute the latitude dependency (using rho as
      ! the Planet radius (approximation).
      vort_all_levs(i,j,k) = vort_all_levs(i,j,k)/                      &
                             (planet_radius*cos_v_latitude(i,j))
    END DO
  END DO

  ! Set vorticity as 0 at the poles.
  IF (at_extremity(psouth)) THEN
    DO i = udims%i_start,udims%i_end
      vort_all_levs(i,vdims%j_start,k) = 0.0
    END DO
  END IF

  IF (at_extremity(pnorth)) THEN
    DO i = udims%i_start,udims%i_end
      vort_all_levs(i,vdims%j_end,k) = 0.0
    END DO
  END IF

END DO

! Do averaging bewteen levels

!Initialize the vorticity array
DO j = vdims%j_start, vdims%j_end
  DO i = udims%i_start, udims%i_end
    vort(i,j) = 0.0
  END DO
END DO

! Do the Summatory over all levels
DO j = vdims%j_start, vdims%j_end
  DO i = udims%i_start, udims%i_end
    DO k = 1,nlevs_avg
      vort(i,j) = vort(i,j) + vort_all_levs(i,j,k)
    END DO
  END DO
END DO

! Divide by the number of levels.
DO j = vdims%j_start, vdims%j_end
  DO i = udims%i_start, udims%i_end
    vort(i,j) = vort(i,j) / REAL(nlevs_avg)
  END DO
END DO

! Do filtering using spherical harmonic truncation
CALL filter_2D(nbot,ntop,nlim,delta_phi,vort,tc_activate,vort_T)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE track_new_TC

END MODULE track_new_TC_mod

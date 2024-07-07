! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE pws_rel_vortic_diag_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_REL_VORTIC_DIAG_MOD'

CONTAINS
!
!  Routine to calculate relative vorticity diagnostic

SUBROUTINE pws_rel_vortic_diag(UField,          &  ! in
                               VField,          &  ! in
                               number_plev,     &  ! in
                               Rel_vortic,      &  ! out
                               icode,cmessage )    ! inout

! Description: Routine to calculate relative vorticity
!
! Method: u,v-winds on B grid and pressure level
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP 003

USE atm_fields_bounds_mod
USE um_stashcode_mod, ONLY: stashcode_pws_sec,                    &
    stashcode_pws_rel_vorticity

USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod, ONLY: swap_field_is_vector
USE missing_data_mod, ONLY: rmdi
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE planet_constants_mod, ONLY: planet_radius
USE trignometric_mod, ONLY: cos_v_latitude
USE um_parvars
USE um_parparams,     ONLY: pnorth, psouth
USE Field_Types,      ONLY: fld_type_u, fld_type_v
USE errormessagelength_mod, ONLY: errormessagelength
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: number_plev ! no.of pressure levels
! U Component of wind
REAL, INTENT(IN)    :: UField(udims%i_start:udims%i_end,               &
                              vdims%j_start:vdims%j_end,number_plev)  
! V Component of wind
REAL, INTENT(IN)    :: VField(udims%i_start:udims%i_end,               &
                              vdims%j_start:vdims%j_end,number_plev)
! rel vorticity on p levels
REAL, INTENT(OUT)   :: Rel_vortic(udims%i_start:udims%i_end,               &
                       vdims%j_start:vdims%j_end,number_plev)
INTEGER :: icode ! Return code : 0 Normal exit, >0 Error exit
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message if ICODE > 0

! Local Constants:
REAL :: HdelX  ! 1/2*W-E grid spacing
REAL :: HdelY  ! 1/2*S-N grid spacing
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PWS_REL_VORTIC_DIAG'

! Loop variables
INTEGER :: i, j, k

! Local Variables:
REAL :: Uhalo(udims_s%i_start:udims_s%i_end,               &
              vdims_s%j_start:vdims_s%j_end)  ! U with halo
REAL :: Vhalo(udims_s%i_start:udims_s%i_end,               &
              vdims_s%j_start:vdims_s%j_end)  ! V with halo
REAL :: DuDy(udims%i_start:udims%i_end,               &
             vdims%j_start:vdims%j_end)      ! latitudinal derivative of u
REAL :: DvDx(udims%i_start:udims%i_end,               &
             vdims%j_start:vdims%j_end)      ! longitudinal derivative of v

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header --------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate relative vorticity on pressure levels
DO k = 1, number_plev
  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      Uhalo(i,j)=UField(i,j,k)
      Vhalo(i,j)=VField(i,j,k)
    END DO
  END DO

! Call Swapbounds to fill in haloes ready for taking derivatives
  CALL swap_bounds(uhalo,udims%i_len,vdims%j_len,1,                &
                   offx,offy,fld_type_u,swap_field_is_vector)

  CALL swap_bounds(vhalo,udims%i_len,vdims%j_len,1,                &
                   offx,offy,fld_type_v,swap_field_is_vector)


! Calculate DvDx - note delta_lambda in radians
  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      HdelX = 0.5/(Planet_Radius * delta_lambda * cos_v_latitude(i,j))
      DvDx(i,j) = HdelX * ( Vhalo(i+1,j) - Vhalo(i-1,j) )
    END DO
  END DO


! Calculate DuDy, scaling u with cos(latitude) - delta_phi in radians
  DO j =  vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      HdelY = 0.5/(Planet_Radius * delta_phi * cos_v_latitude(i,j))
      DuDy(i,j) = HdelY * (Uhalo(i,j+1)*cos_v_latitude(i,j+1) &
                         - Uhalo(i,j-1)*cos_v_latitude(i,j-1))
    END DO
  END DO

! Scale Rel Vorticity to conventional units
  Rel_vortic(:,:,k) = (DvDx(:,:) - DuDy(:,:))*1.0E6

  IF (at_extremity(PSouth) ) THEN
    Rel_vortic(:,vdims%j_start,k) = Rel_vortic(:,vdims%j_start+1,k)
  END IF
  IF (at_extremity(PNorth) ) THEN
    Rel_vortic(:,vdims%j_end,k) = Rel_vortic(:,vdims%j_end-1,k)
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_rel_vortic_diag

END MODULE pws_rel_vortic_diag_mod

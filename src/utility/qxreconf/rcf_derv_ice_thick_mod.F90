! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Derive sea-ice thickness from sea-ice fraction.

MODULE Rcf_Derv_Ice_Thick_Mod

!  Subroutine Rcf_Derv_Ice_Thick
!
! Description:
!   Derive Ice Thickness from the Ice Fraction field.
!
! Method:
!   Ice Thickness is derived from the Ice Fraction. The ice thickness
!   is set to 2 metres and 1 metre for Northern and Southern hemisphere
!   grid points respectively if the value of ice fraction is greater
!   than zero, ie. if sea-ice is present.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DERV_ICE_THICK_MOD'

CONTAINS

SUBROUTINE Rcf_Derv_Ice_Thick ( fields, field_count, hdr,    &
                                ice_thickness )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type,                &
    output_grid

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Write_Field_Mod, ONLY: &
    Rcf_Write_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE decomp_params, ONLY: &
    decomp_rcf_output

USE um_stashcode_mod, ONLY: &
    stashcode_icefrac,         &
    stashcode_icethick,        &
    stashcode_prog_sec

USE Rcf_HeadAddress_Mod, ONLY: &
    RC_LongSpacing,     RC_LatSpacing,     &
    RC_FirstLat,        RC_FirstLong,      &
    RC_PoleLat,         RC_PoleLong

USE UM_ParVars, ONLY: &
    g_datastart

USE UM_ParCore, ONLY: &
    mype

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Normal

USE latlon_eq_rotation_mod, ONLY: rotate_eq_to_latlon

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ),     POINTER       :: fields(:)
INTEGER,                INTENT(IN)    :: field_count
TYPE( um_header_type ), INTENT(IN)    :: hdr
TYPE( field_type ),     INTENT(INOUT) :: ice_thickness

! Local variables
REAL, PARAMETER                    :: NP_ice_thickness = 2.0
REAL, PARAMETER                    :: SP_ice_thickness = 1.0
INTEGER                            :: pos_icef
INTEGER                            :: i,j,ij

REAL                               :: latitude                         &
                                      (output_grid % loc_p_row_length, &
                                       output_grid % loc_p_rows)
REAL                               :: longitude                        &
                                      (output_grid % loc_p_row_length, &
                                       output_grid % loc_p_rows)

TYPE( field_type ), POINTER        :: ice_fraction

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_DERV_ICE_THICK'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-------------------------------------------------------------------
! Locate Ice Fraction in output dump
!-------------------------------------------------------------------

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_icefrac,               &
                  fields, field_count, pos_icef )

!--------------------------------------------------------------------
! Allocate space for ice fraction and read in.
!--------------------------------------------------------------------

ice_fraction => fields( pos_icef )
CALL Rcf_Alloc_Field( ice_fraction )
CALL Rcf_Read_Field( ice_fraction, hdr, decomp_rcf_output )

!--------------------------------------------------------------------
! Set up a latitude/longitude field for the grid.
!--------------------------------------------------------------------

DO j = 1, output_grid % loc_p_rows
  DO i = 1, output_grid % loc_p_row_length

    latitude (i,j) = hdr % RealC(RC_FirstLat) +  &
                    (g_datastart(2,mype)+j-2) *  &
                     hdr % RealC(RC_LatSpacing)

    longitude(i,j) = hdr % RealC(RC_FirstLong) +  &
                    (g_datastart(1,mype)+i-2)  *  &
                     hdr % RealC(RC_LongSpacing)

  END DO
END DO

!--------------------------------------------------------------------
! For rotated grids, get true lats and longs.
!--------------------------------------------------------------------

IF (output_grid % rotated) THEN

  ! Use the same arrays to store true lats & longs

  CALL rotate_eq_to_latlon(latitude, longitude, latitude, longitude,         &
                           hdr % RealC(RC_PoleLat), hdr % RealC(RC_Polelong),&
                           output_grid % loc_p_field)

END IF

!--------------------------------------------------------------------
! Set up the ice thickness field.
!--------------------------------------------------------------------

ice_thickness % DATA (:,1) = 0.0

ij =0
DO j = 1, output_grid % loc_p_rows
  DO i = 1, output_grid % loc_p_row_length

    ij = ij + 1

    IF ( ice_fraction % DATA (ij,1) >  0.0 ) THEN

      IF ( latitude (i,j) >  0.0 ) THEN               ! In Northern Hem
        ice_thickness % DATA (ij,1) = NP_ice_thickness
      ELSE                                            ! In Southern Hem
        ice_thickness % DATA (ij,1) = SP_ice_thickness
      END IF

    END IF

  END DO
END DO

IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
  WRITE(umMessage,*) 'Ice Thickness derived from Ice Fraction '
  CALL umPrint(umMessage,src='rcf_derv_ice_thick_mod')
END IF

!--------------------------------------------------------------------
! Deallocate space for ice fraction
!--------------------------------------------------------------------

CALL Rcf_Dealloc_Field( ice_fraction )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Derv_Ice_Thick
END MODULE Rcf_Derv_Ice_Thick_Mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reinitialises sea ice category fields for thickness and concentration

MODULE Rcf_Derv_Ice_Cat_Thick_Mod

IMPLICIT NONE

! Description: Populate the category fields for ice thickness and ice
!              fraction based on the Grid Box Mean (GBM) of both fields.
!              Done by generating bounds for the thickness of each
!              category (pseudo level) and assinging values for fraction
!              and thickness by mapping where the GBM thicknesses fall
!              within the bounds for that category (pseudo level).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DERV_ICE_CAT_THICK_MOD'

CONTAINS

SUBROUTINE Rcf_Derv_Ice_Cat_Thick( stash_item, fields_out,            &
                                   field_count_out, hdr_out,          &
                                   ice_cat_field )

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE um_stashcode_mod, ONLY: &
    stashcode_icefrac,         stashcode_icethick,         &
    stashcode_ice_conc_cat,    stashcode_ice_thick_cat,    &
    stashcode_prog_sec

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Normal

USE UM_ParCore, ONLY: &
    mype

USE decomp_params, ONLY: &
    decomp_rcf_output

USE jules_sea_seaice_mod, ONLY: &
    Nice

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim


IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( um_header_type), INTENT(IN) :: hdr_out
TYPE( field_type ), INTENT(INOUT), TARGET :: ice_cat_field
INTEGER, INTENT(IN)               :: STASH_Item
INTEGER, INTENT(IN)               :: field_count_out

! Internal variables
TYPE( field_type ), POINTER       ::  ice_thick
TYPE( field_type ), POINTER       ::  ice_frac
TYPE( field_type ), POINTER       ::  ice_cat_thick
TYPE( field_type ), POINTER       ::  ice_cat_frac

REAL                              ::  rNice ! real of 1 / Nice
REAL                              ::  HIN_Max(0:Nice)
REAL                              ::  cc1,cc2,cc3,x1

INTEGER                           ::  pos   ! position in array
INTEGER                           ::  i,j,k ! loop index

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_DERV_ICE_CAT_THICK'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  IF ( STASH_Item == stashcode_ice_conc_cat ) THEN
    WRITE(umMessage,*) 'Reinitialising Ice Concentrations on categories'
    CALL umPrint(umMessage,src='rcf_derv_ice_cat_thick_mod')
  ELSE IF ( STASH_Item == stashcode_ice_thick_cat ) THEN
    WRITE(umMessage,*) 'Reinitialising Ice Thickness on categories'
    CALL umPrint(umMessage,src='rcf_derv_ice_cat_thick_mod')
  END IF
END IF

!----------------------------------------------------------------------
! Find required fields in output dump and read them in where available
!----------------------------------------------------------------------
! GBM ice thickness; will abort if icethick not found

CALL Rcf_Locate( stashcode_prog_sec, stashcode_icethick,           &
                 fields_out,field_count_out,pos)

ice_thick => fields_out(pos)
CALL Rcf_Alloc_Field( ice_thick )
CALL Rcf_Read_Field( ice_thick, hdr_out, decomp_rcf_output )

! GBM ice fraction; will abort if icefrac not found

CALL Rcf_Locate( stashcode_prog_sec, stashcode_icefrac,            &
                 fields_out,field_count_out,pos)

ice_frac => fields_out(pos)
CALL Rcf_Alloc_Field( ice_frac )
CALL Rcf_Read_Field( ice_frac, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Allocate space, and pointers.
!----------------------------------------------------------------------
! Set up temporary space for one of ice_conc_cat or ice_thick_cat and
! point the other to the space passed in as ice_cat_field.

IF ( STASH_item == stashcode_ice_conc_cat ) THEN
  ! Allocate Category ice thickness; will abort if icethick not found
  CALL Rcf_Locate(stashcode_prog_sec, stashcode_ice_thick_cat,     &
                  fields_out, field_count_out, pos)

  ice_cat_thick => fields_out(pos)
  CALL Rcf_Alloc_Field( ice_cat_thick )

  ice_cat_frac => ice_cat_field
ELSE
  ! Allocate Category ice fraction; will abort if icefrac not found
  CALL Rcf_Locate(stashcode_prog_sec, stashcode_ice_conc_cat,      &
                  fields_out, field_count_out, pos)

  ice_cat_frac => fields_out(pos)
  CALL Rcf_Alloc_Field( ice_cat_frac )

  ice_cat_thick => ice_cat_field
END IF


!----------------------------------------------------------------------
! Calculate some of the constants
!----------------------------------------------------------------------

rNice = 1.0 / REAL(Nice)
! calcluate the thickness category offset scalar
cc1 = 3.0 * rNice
! calcluate the thickness category adjustment scalar
cc2 = 15.0 * cc1
! CC3 is the thickness category phase adjustment scaler
cc3 = 3.0

! use different minimum thickness in special case Nice=1
IF ( Nice == 1 ) THEN
  HIN_Max(0) = 0.10
ELSE
  HIN_Max(0) = 0.00
END IF

!----------------------------------------------------------------------
! Loop through ice_cat_field pseudo levels
!----------------------------------------------------------------------

DO i = 1,Nice

  ! Calculate thickness category domain.
  x1 = REAL ( i - 1 ) * rNice
  HIN_Max(i) = HIN_Max(i-1) + cc1 +                                 &
               cc2 * (1.0 + TANH( cc3 * ( x1 - 1.0 ) ) )

  !Ensure upper limit for HIN_Max is an arbitrary LARGE value
  HIN_Max(Nice) = HUGE(HIN_Max(0))

  ! Map ice thickness and fraction to thier psuedolevel equivalents.
  WHERE ( HIN_Max(i-1) <= ice_thick % DATA(:,1)  .AND. &
          ice_thick % DATA(:,1) < HIN_Max(i) )
    ice_cat_thick % DATA(:,i) = ice_thick % DATA(:,1)
    ice_cat_frac % DATA(:,i) = ice_frac % DATA(:,1)
  ELSEWHERE
    ice_cat_thick % DATA(:,i) = 0.00
    ice_cat_frac % DATA(:,i) = 0.00
  END WHERE

END DO

! Desired field has now been 'populated' as the relevant pointer
! was used to point to ice_cat_field earlier, and ice_cat_field is the
! field returned to the calling routine via the argument list.

!----------------------------------------------------------------------
! Clear up dynamic memory used
!----------------------------------------------------------------------

CALL Rcf_Dealloc_Field( ice_thick )
CALL Rcf_Dealloc_Field( ice_frac )

IF ( STASH_item == stashcode_ice_conc_cat ) THEN
  CALL Rcf_Dealloc_Field( ice_cat_thick )
ELSE
  CALL Rcf_Dealloc_Field( ice_cat_frac )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Derv_Ice_Cat_Thick

END MODULE Rcf_Derv_Ice_Cat_Thick_mod

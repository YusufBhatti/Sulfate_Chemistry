! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Set up lists and space to store Header info for data created not read
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

MODULE Rcf_Grib_Spcl_Hdr_Mod
IMPLICIT NONE

! Description: Routine which creates lists and list members to hold
!              the header info for fields which are created rather than
!              read from the original GRIB data.
!              -NB this routine should not be used to alter the headers
!               of fields which are altered (e.g. geopotential to orog)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_SPCL_HDR_MOD'

CONTAINS
SUBROUTINE Rcf_Grib_Spcl_Hdr(Lists)

USE Rcf_GRIB_Block_Params_Mod, ONLY:                                      &
    List_Marker,          Grib_Record,                                     &
    p_Lvl_Type,           p_Lvl_Desc_1,                                    &
    Tb3_Pressure,         p_Orig_cntr,                                     &
    p_Param_ID,           EID_Temperature,                                 &
    Grb_Data_Real

USE Rcf_GRIB_Lookups_Mod, ONLY:                                           &
    grib_max_fields,           grib_Exner_field,                           &
    grborigecmwf

USE um_stashcode_mod, ONLY:                                             &
    Stashcode_exner

USE Lookup_Addresses

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE


! Global variables (#include statements etc):

! Subroutine arguments

!< Array  arguments with intent(InOut):>
TYPE (List_Marker), INTENT(INOUT) :: Lists(0:grib_max_fields)

! Local constants
CHARACTER (LEN=*), PARAMETER     :: RoutineName='RCF_GRIB_SPCL_HDR'

! Local variables

TYPE (Grib_Record),POINTER       :: New_Record,Current

INTEGER                          :: i
INTEGER                          :: cnter,Criteria
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!=======================================================================
!  Routine Code Start :
!=======================================================================

!=======================================================================
!  Loop across all lists, finding those which will require extra fields
!=======================================================================
DO i = 1, grib_max_fields

  IF (ASSOCIATED(Lists(i) % Begin) ) THEN

    ! Lets only look at those ECMWF ones first
    IF ( Lists(i) % Begin % Block_1 (p_Orig_cntr) == grborigecmwf ) THEN
      SELECT CASE ( Lists(i) % Begin % Block_1 ( p_Param_ID ) )

      CASE (EID_Temperature)
        ! Will need Exner for Height generation

        ! Set pointer to first entry in list
        Current => Lists(i) % Begin
        cnter =0

        ! Loop across members in list
        DO WHILE (ASSOCIATED(Current))
          cnter= cnter+1

          !Allocate a GRIB header to store info
          WRITE(umMessage,*) 'About to allocate for Exner ', cnter
          CALL umPrint(umMessage,src='rcf_grib_spcl_hdr_mod')
          ALLOCATE(New_Record)

          New_Record % Block_1(:) = Current % Block_1(:)
          New_Record % Block_2(:) = Current % Block_2(:)
          New_Record % Block_3(:) = Current % Block_3(:)
          New_Record % Block_4(:) = Current % Block_4(:)
          New_Record % Block_R(:) = Current % Block_R(:)

          New_Record % Start_pos  = 0
          New_Record % Data_Type  = Grb_Data_Real
          New_Record % Num_Fp     = 0
          New_Record % Num_Vert   = 0
          New_Record % Num_Bitmap = 0
          New_Record % Num_Quasi  = 0

          ! Copy header info from Temp to Exner
          New_Record % StashCode = Stashcode_exner
          New_Record % Desc      = 'Exner field'

          New_Record % Block_1 (p_Orig_cntr) = (-1)
          New_Record % Block_1 (p_Param_ID)  = (1)

          ! Assign that record to it's list
          New_Record % Prev   => Lists(grib_Exner_field) % END
                           ! Point Prev pointer at end of
                           ! current list.(Null if first entry)
          IF (ASSOCIATED(Lists(grib_Exner_field) % END)) THEN
                           ! If current end of list is a
                           ! valid record (Not first entry)
            Lists(grib_Exner_field) % END % Next  => New_Record
                           ! Point 'next' for previous entry
                           ! at current entry

          ELSE                         ! Else : must be 1st entry
            Lists(grib_Exner_field) % Begin  => New_Record
                           ! Point begining of List at New_Record
          END IF

          Lists(grib_Exner_field) % END      => New_Record
                           ! Point End of List at (now complete)
                           ! New_Record Entry
          NULLIFY(New_Record % Next)   ! Ensure 'Next' is not associated

          Lists(grib_Exner_field) % LstCount =                             &
              Lists(grib_Exner_field) % LstCount + 1
                           ! Add one to count of list size


          ! end do across list members
          Current => Current % Next
        END DO

      END SELECT

    END IF
  END IF

  !=======================================================================
  !  End Loop across all lists.
  !=======================================================================
END DO  ! loop over all lists


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Grib_Spcl_Hdr
END MODULE Rcf_Grib_Spcl_Hdr_Mod

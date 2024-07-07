! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Perform cconsistancy checks on GRIB data

MODULE Rcf_Grib_Check_Mod
USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage
IMPLICIT NONE


! SUBROUTINE Rcf_Grib_Check  Perform Basic consistancy checks on GRIB
!                            Data being handled.
!-and-
! FUNCTION Find_Match  - recursively traverse a list to find a match
!
! Description: A Routine where checks are performed on the GRIB data
!              read in to ensure it matches expectations.
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_CHECK_MOD'

CONTAINS
SUBROUTINE Rcf_Grib_Check(Lists)

! uses variables and routines from other modules

USE Rcf_GRIB_Block_Params_Mod, ONLY:    &
  List_Marker,                           &
  Grib_Record,                           &
  p_Lvl_Type,                            &
  p_Param_ID,                            &
  p_Lvl_Desc_1,                          &
  p_Lvl_Desc_2,                          &
  EID_Log_Surf_Press,                    &
  Tb3_Surface,                           &
  Tb3_Pressure,                          &
  Tb3_Hybrid

USE Rcf_GRIB_Lookups_Mod, ONLY:    &
  grib_max_fields,             &
  grib_Soil_Temp_field,        &
  grib_Soil_Moist_field

USE EReport_Mod, ONLY:     &
    EReport, newline

USE Rcf_GRIB_T_n_Pstar_H_Interp_Mod, ONLY:  &
    ak, bk, akh, bkh

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE


! Subroutine arguments

!< Array  arguments with intent(inOut):>
TYPE (List_Marker), INTENT(INOUT) :: Lists(0:grib_max_fields)
!An array of pointer pairs to form the head and tail of the field lists

! Local variables

TYPE (Grib_Record),POINTER       :: Current, Compare

CHARACTER (LEN=*), PARAMETER     :: RoutineName='RCF_GRIB_CHECK'
CHARACTER (LEN=errormessagelength)  :: CMessage   ! used for EReport
INTEGER                          :: ErrorStatus   ! used for EReport
INTEGER                          :: i,levelcount,Ref_List,pltype
INTEGER                          :: nvc,nlevs

LOGICAL                          :: Match, Match_Vert
LOGICAL                          :: L_First_Press, L_First_Hybrid

REAL                             :: soildepth
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!=======================================================================
! Main Routine.
!=======================================================================

!=======================================================================
! Perform Basic Data Integrity Checks
!=======================================================================
                   ! flag to identify first pressure level encountered
L_First_Press    = .TRUE.
                   ! flag to identify first hybrid level encountered
L_First_Hybrid   = .TRUE.

!=======================================================================
!  Loop across all lists, checking they have valid contents
!=======================================================================
DO i = 1, grib_max_fields

  IF (ASSOCIATED(Lists(i) % Begin) ) THEN

    !=====================================================
    ! Special check to see if log(pstar) instead of pstar
    !=====================================================
    ! If log(pstar) is field in surface pressure list then
    !  alter level information to make it a surface field,
    !  not a model level 1 field. Use special routines
    !  later to convert from log(pstar) to pstar.
    IF (Lists(i)%Begin%Block_1(p_Param_ID)==EID_Log_Surf_Press) THEN
      WRITE(umMessage,*) "Found log(pstar), modifying level information"
      CALL umPrint(umMessage,src='rcf_grib_check_mod')
      Lists(i) % Begin % Block_1(p_Lvl_Type) = Tb3_Surface
      Lists(i) % Begin % Block_1(p_Lvl_Desc_1) = 0
      Lists(i) % Begin % Block_1(p_Lvl_Desc_2) = 0
    END IF

    !=====================================
    ! Special check for soil level fields
    !=====================================
    ! If soil field then store integer index in p_Lvl_Desc_1
    !  and store soil depth in p_Lvl_Desc_2
    ! Note: Fields already sorted
    IF ( i == grib_Soil_Temp_field .OR. &
         i == grib_Soil_Moist_field ) THEN
      WRITE(umMessage,*) "Modifying soil level information"
      CALL umPrint(umMessage,src='rcf_grib_check_mod')
      Current => Lists( i ) % Begin
      levelcount =  0
      DO WHILE ( ASSOCIATED( Current ) )
        levelcount = levelcount + 1
        soildepth = Current % Block_1 ( p_Lvl_Desc_2 ) -            &
                    Current % Block_1 ( p_Lvl_Desc_1 )
        Current % Block_1 ( p_Lvl_Desc_1 ) = levelcount
        Current % Block_1 ( p_Lvl_Desc_2 ) = soildepth
        Current => Current % Next
      END DO
    END IF

    !==========================================
    ! Check No. 1 - Consistant Vertical Levels
    !==========================================
    ! Check 1.1 Levels not repeated for a parameter
    !       1.2 List count matches no. of levels found
    !       1.3 List counts for different multi-level params match
    !       1.4 Levels for different multi-level params also match

    ! Check to see if list lies on Pressure levels or Model levels
    pltype=Lists(i) % Begin % Block_1(p_Lvl_Type)
    IF ( pltype == Tb3_Pressure .OR. pltype == Tb3_Hybrid ) THEN

      ! Check 1.1 - Levels not repeated for a parameter
      !------------------------------------------------
      Current => Lists( i ) % Begin

      DO WHILE ( ASSOCIATED( Current % Next ))

        Match = Find_Match (Current, Current % Next, p_Lvl_Desc_1)

        IF (Match) THEN     ! This level has same descriptor as another
          WRITE (CMessage,'(A)') "List " //                           &
              Lists( i ) % Begin % Desc     //                        &
              " :Two levels have the same value"
          ErrorStatus = 10
          CALL EReport( RoutineName, ErrorStatus, Cmessage)
        END IF

        Current => Current % Next         ! move to next item in list I

      END DO

      ! Check 1.2  - count
      !-------------------
      Current => Lists( i ) % Begin
      levelcount =  0

      DO WHILE ( ASSOCIATED( Current ))
        levelcount = levelcount + 1
        Current => Current % Next
      END DO

      IF (levelcount /= Lists( i ) % LstCount) THEN
        WRITE (CMessage,'(A)') "List " //                           &
            Lists( i ) % Begin % Desc     //                        &
            " :  Count of entries does not match recorded count"
        ErrorStatus = -20
        CALL EReport( RoutineName, ErrorStatus, Cmessage)
        Lists( i ) % LstCount = levelcount
      END IF

      ! Is it the first pressure or model level based list found ?
      IF ( ( pltype == Tb3_Pressure .AND. L_First_Press) .OR. &
           ( pltype == Tb3_Hybrid   .AND. L_First_Hybrid) ) THEN

        ! Are pressure levels and hybrid levels mixed up in same file?
        IF (.NOT. ( L_First_Hybrid .AND. L_First_Press)) THEN
          WRITE (CMessage,'(A)') "List " //                          &
              Lists( i ) % Begin % Desc     //                       &
              " :  Cannot mix up pressure and model level fields"
          ErrorStatus = 20
          CALL EReport( RoutineName, ErrorStatus, Cmessage)
        END IF

        ! Check 1.2.1 - Vertical Coordinates consistent across records
        !-------------------------------------------------------------
        IF ( pltype == Tb3_Hybrid ) THEN
          Current => Lists( i ) % Begin
          IF ( ASSOCIATED( Current%Next )) THEN
            Match_Vert = Find_Match_Vert (Current % Next, Current)
            IF (.NOT. Match_Vert) THEN
              WRITE (CMessage,'(A)') "List " //                  &
                  Lists( i ) % Begin % Desc     //               &
                  " : Vertical Coordinates not consistent"
              ErrorStatus = 50
              CALL EReport( RoutineName, ErrorStatus, Cmessage)
            END IF
          END IF
        END IF

        ! record list as 'reference list'
        Ref_list = i
        IF ( pltype == Tb3_Pressure) L_First_Press    = .FALSE.
        IF ( pltype == Tb3_Hybrid  ) L_First_Hybrid   = .FALSE.

        ! Set ECMWF level definitions using Vertical Coordinates
        ! Note that ECMWF data is released on full model levels, while the
        !  vertical coords that document this are on half levels. The full
        !  levels are exactly half way between the half levels in terms of
        !  pressure hence the ak_k = 0.5 * (akh_km1 + akh_k)
        IF ( pltype == Tb3_Hybrid ) THEN
          Current => Lists( i ) % Begin
          nvc = Current%Num_Vert    ! Number of vertical coordinates
          nlevs = nvc/2 - 1         !  nvc = 2*(nlevs+1)
          ALLOCATE(akh(0:nlevs), bkh(0:nlevs))
          ! Store half level vertical coords
          akh(0:nlevs) = Current%VertCoords(1:nlevs+1)
          bkh(0:nlevs) = Current%VertCoords(nlevs+2:nvc)
          ! Calculate full level vertical coordinates
          ALLOCATE(ak(nlevs), bk(nlevs))
          ak = 0.5 * (akh(0:nlevs-1) + akh(1:nlevs))
          bk = 0.5 * (bkh(0:nlevs-1) + bkh(1:nlevs))

          ! Does number of fields in list match number of levels
          ! from vertical coordinates? If not then throw up warning ...
          IF (levelcount /= nlevs) THEN
            WRITE (CMessage,'(A,I4,A,I4)') "List " //          &
                Lists( i ) % Begin % Desc  //                  &
                " : Number of levels not consistent with vertical coordinates" &
                //newline//"Levels = ",levelcount,"; VCs = ",nlevs
            ErrorStatus = -70
            CALL EReport( RoutineName, ErrorStatus, Cmessage)
          END IF
        END IF

      ELSE                          ! make comparisons against ref list

        ! Check 1.3 - Match reference levels
        !-----------------------------------
        Current => Lists( i ) % Begin
        Compare => Lists( Ref_list ) % Begin

        DO WHILE ( ASSOCIATED( Current ))

          Match = Find_Match (Current, Compare, p_Lvl_Desc_1)

          IF (.NOT. Match) THEN   ! This level did not match one in Ref
            WRITE (CMessage,'(A)') "List " //                           &
                Lists( i ) % Begin % Desc     //                        &
                " :  level did not match one in ref field"
            ErrorStatus = 30
            CALL EReport( RoutineName, ErrorStatus, Cmessage)
          END IF

          Current => Current % Next       ! move to next item in list I
        END DO

        ! Check 1.4 - Compare counts
        !---------------------------
        IF (Lists( i ) % LstCount /= Lists( Ref_list ) % LstCount) THEN
          WRITE (CMessage,'(A)') "List " //                           &
              Lists( i ) % Begin % Desc     //                        &
              " :  level count did not match reference"
          ErrorStatus = 40
          CALL EReport( RoutineName, ErrorStatus, Cmessage)
        END IF

        ! Check 1.4.1 - Vertical Coordinates consistent across records
        !-------------------------------------------------------------
        IF ( pltype == Tb3_Hybrid ) THEN
          Current => Lists( i ) % Begin
          Compare => Lists( Ref_list ) % Begin
          IF ( ASSOCIATED( Current )) THEN
            Match_Vert = Find_Match_Vert (Current, Compare)
            IF (.NOT. Match_Vert) THEN
              WRITE (CMessage,'(A)') "List " //                  &
                  Lists( i ) % Begin % Desc     //               &
                  " : Vertical Coordinates did not match reference"
              ErrorStatus = 60
              CALL EReport( RoutineName, ErrorStatus, Cmessage)
            END IF
          END IF
        END IF

      END IF  ! Check for First level

    END IF ! level check

  END IF ! Begin of list was Associated
END DO  ! loop over all lists

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Grib_Check

!-----------------------------------------------------------------------
!  Function Find_Match
!-----------------------------------------------------------------------

RECURSIVE FUNCTION Find_Match ( Current, Compare, Criteria )          &
                                   RESULT ( Match )

! Description: Recursively traverse a list, exiting on positive match
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
USE Rcf_GRIB_Block_Params_Mod, ONLY:    &
  Grib_Record

USE ereport_mod, ONLY: ereport
IMPLICIT NONE

! Subroutine arguments

!< Scalar arguments with intent(in):>
INTEGER, INTENT(IN)                      :: Criteria

!< Array  arguments with intent(in):>
TYPE (Grib_Record), POINTER              :: Current, Compare

! Local variables

LOGICAL                          :: Match

!-----------------------------------------------------------------------
! Begin routine
!-----------------------------------------------------------------------

Match = .FALSE.

IF ( Current % Block_1(Criteria) == Compare % Block_1(Criteria) ) THEN
  Match = .TRUE.
ELSE
  IF ( ASSOCIATED ( Compare % Next ) ) THEN
    Match = Find_Match ( Current, Compare % Next , Criteria )
  END IF
END IF


END FUNCTION Find_Match

!-----------------------------------------------------------------------
!  Function Find_Match_Vert
!-----------------------------------------------------------------------

RECURSIVE FUNCTION Find_Match_Vert ( Current, Compare )          &
                                   RESULT ( Match )

! Description: Recursively traverse a list, exiting on negative match
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
USE Rcf_GRIB_Block_Params_Mod, ONLY:    &
  Grib_Record

USE ereport_mod, ONLY: ereport
IMPLICIT NONE

! Subroutine arguments

!< Array  arguments with intent(in):>
TYPE (Grib_Record), POINTER              :: Current, Compare

! Local variables

LOGICAL                          :: Match, temp_Match
INTEGER                          :: i      ! loop counter

!-----------------------------------------------------------------------
! Begin routine
!-----------------------------------------------------------------------

Match = .TRUE.

IF (Current % Num_Vert == Compare % Num_Vert ) THEN
  DO i=1,Current % Num_Vert
    IF (Current % VertCoords(i) /= Compare % VertCoords(i)) &
           Match = .FALSE.
  END DO
ELSE
  Match = .FALSE.
END IF

IF ( Match ) THEN
  IF ( ASSOCIATED ( Current % Next ) ) THEN
    Match = Find_Match_Vert ( Current % Next, Compare )
  END IF
END IF

END FUNCTION Find_Match_Vert

END MODULE Rcf_Grib_Check_Mod

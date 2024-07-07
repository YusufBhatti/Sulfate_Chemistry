! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  order a linked list order using criteria provided

MODULE Rcf_Grib_FldSort_mod

! SUBROUTINE rcf_Grib_FldSort  Establish what information is contained
!                               within the GRIB file.
!
! Description: Routine to re-order the contents of a linked list of
!              GRIB field headers in either ascending or descendng
!              order of any one parameter in Block_0
!
! Method: Quite simply, Bubblesort. It compares 2 adjacent entries and
!         if they aren't in the correct order (relative to each other)
!         swap all the pointers so they appear in the opposite order
!         within the list. Repeat this process as you traverse through
!         the list. Then go back to the list start and repeat again
!         until no change occured during a complete pass through the
!         list.
!         Not the most efficient method but simple and suitable
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_FLDSORT_MOD'

CONTAINS
SUBROUTINE rcf_Grib_FldSort(Markers,Criteria,Order)

USE Rcf_GRIB_Block_Params_Mod

USE Rcf_GRIB_Lookups_Mod               ! Contains cross ref table for
                                       ! Stash => ECMWF parameter Id's
                                       ! and appropriate params
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments

!< Scalar arguments with intent(in):>
INTEGER, INTENT(IN)               :: Criteria ! the field being compared

!< Logical arguments with intent(in):>
LOGICAL, INTENT(IN), OPTIONAL     :: Order  ! The order desired
                                            ! .TRUE. = Ascending

!< Array  arguments with intent(InOut):>
TYPE (List_Marker), INTENT(INOUT)       :: Markers

! Local variables

CHARACTER (LEN=*), PARAMETER     :: RoutineName='RCF_GRIB_FLDSORT'
INTEGER                          :: ErrorStatus
INTEGER                          :: list_count
TYPE(Grib_Record), POINTER       :: Current,Stored

LOGICAL                          :: No_Change_Pass
LOGICAL                          :: Ascending
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!  Initialise bits and bobs.
!-----------------------------------------------------------------------

! Deal with Optional dummy argument 'Order'
IF (PRESENT(Order)) THEN
  Ascending = Order
ELSE
  Ascending = .TRUE.     ! Default to Ascending
END IF

No_Change_Pass = .FALSE.
list_count = 0

!-----------------------------------------------------------------------
!  Do until a full pass through the list results in no change.
!-----------------------------------------------------------------------
DO WHILE (.NOT. No_Change_Pass)

  list_count = list_count + 1

  Current => Markers % Begin             ! Move to start of list
  No_Change_Pass = .TRUE.                ! re-set fo 'no change'

  !-----------------------------------------------------------------------
  !  Loop through the list members
  !-----------------------------------------------------------------------
  DO WHILE (ASSOCIATED(Current%next))    ! Do while there is an entry
                                         ! Which follows current

    !-----------------------------------------------------------------------
    !  Check the ordering of current and the entry which follows it
    !-----------------------------------------------------------------------

    IF ((Ascending .AND.                                               &
       (Current%Next%Block_1(Criteria) < Current%Block_1(Criteria)))  &
         .OR.                                                         &
       (.NOT. Ascending .AND.                                           &
       (Current%Next%Block_1(Criteria) > Current%Block_1(Criteria)))) &
         THEN

      !-----------------------------------------------------------------------
      !  If order is incorect then swap the two entrys
      !-----------------------------------------------------------------------

      No_Change_Pass = .FALSE.

      !-----------------------------------------------------------------------
      !  Remove entry _after_ current and update links from current to
      !  the entry which follows that.(if it exists)
      !-----------------------------------------------------------------------

      Stored       => Current%Next      ! use Stored as a pointer
      Current%next => Stored%Next       ! Might be Null(last entry)

      IF (ASSOCIATED(Stored%Next)) THEN !  more entries
        Stored%Next%Prev  => Current    ! Point back from next to cur
      ELSE                              ! Stored was last entry
        Markers % END     => Current    ! Current is new last entry
      END IF

      !-----------------------------------------------------------------------
      !  re-insert Stored between Current and Current%last.
      !-----------------------------------------------------------------------

      IF (ASSOCIATED(Current%Prev)) THEN ! There are entries B4 this 1
        Current%Prev%Next => Stored
      ELSE                              ! going in as first in list
        Markers % Begin   => Stored
      END IF

      Stored%Next       => Current
      Stored%Prev       => Current%Prev ! This could still be Null
      Current%Prev      => Stored

    ELSE                                ! Pair were correctly ordered
      Current  => Current%Next          ! Move to next entry
    END IF

  END DO

END DO


!-----------------------------------------------------------------------
!  Done : - return
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE rcf_Grib_FldSort
END MODULE rcf_Grib_FldSort_mod

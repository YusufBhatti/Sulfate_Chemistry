! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE LHN_SEARCH_1A -----------------------------------------
!
!    Purpose : Search for suitable nearby latent heating profiles for
!              use in the LHN scheme.
!
!    Version 1A:
!    This version works on a local field in the logical processed grid
!    domain which has halos of srange size.
!
!    ATTENTION: However, the routine returns in 'near' the global index
!    of the searched point NOT the local one.
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!    Documentation :  FR  WP 171
!
!
!    ARGUMENTS
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE lhn_search_1a_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LHN_SEARCH_1A_MOD'

CONTAINS

SUBROUTINE lhn_search_1a(x,y,near,srange, search,                       &
  rows_l, row_length_l, rows_g, row_length_g,                          &
  tot_pr,anal_pr,radius)
  

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE ac_control_mod, ONLY: epsilon_lhn, lhn_diag

IMPLICIT NONE

! Arguments 
INTEGER, INTENT(IN)    :: x,y             ! Location of the current point 
INTEGER, INTENT(IN)    :: srange          ! Search range in grid points
INTEGER, INTENT(IN)    :: rows_l          ! Rows in local grid
INTEGER, INTENT(IN)    :: row_length_l    ! Points in each row of local grid
INTEGER, INTENT(IN)    :: rows_g          ! Rows in global grid
INTEGER, INTENT(IN)    :: row_length_g    ! Points in each row of global grid
INTEGER, INTENT(IN)    :: search(4*srange*(srange+1),2)             
                                          ! Template holding relative positions
                                          ! of search points from current point
INTEGER, INTENT(INOUT) :: radius(5)       ! Diagnostic for breakdown of searches
REAL, INTENT(IN)       :: tot_pr(1-srange:row_length_l+srange, & 
                                 1-srange:rows_l+srange)
                                          ! Model total precip rates
                                          ! with halos of srange size
REAL, INTENT(IN)        :: anal_pr        ! Analysed ppn rate at current point
INTEGER, INTENT(OUT)    :: near           ! Nearest point with suitable profile

! Local arrays and variables
INTEGER ::    i_g          ! Column global coord of searched point
INTEGER ::    j_g          ! Row     "    "   "
INTEGER ::    i_l          ! Column local coord of searched point
INTEGER ::    j_l          ! Row     "    "   "
INTEGER ::    p_g          ! Global index of point at (i_g,j_g)
INTEGER ::    jr , jn      ! Loop counters
REAL ::       ratio        ! ratio of model to obs ppn
REAL ::       best         ! closest RATIO to 1
LOGICAL ::    l_found      ! false until a suitable point is found

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LHN_SEARCH_1A'

!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!
! Loop round nearest points in ever increasing radius,
! looking for best fit to observed rain at POINT
!
! Initialise variables
!
l_found = .FALSE.
best = 0.0
ratio = 0.0
!
! Loop over possible ranges
!
DO jr = 1 , srange
  IF (.NOT. l_found) THEN
    !
    ! If not found one yet then loop over points at this range
    !
    DO jn = 1 , 8*jr      

      ! Compute local coordinate of searched point
      i_l = x + search( 4*jr*(jr-1) +jn , 1)
      j_l = y + search( 4*jr*(jr-1) +jn , 2)
      
      ! Compute global coordinate of searched point
      i_g = g_datastart(1, mype) - 1 + i_l
      j_g = g_datastart(2, mype) - 1 + j_l

      !
      ! Make sure searched point is within the global grid.
      !
      IF ( i_g  >=  1 .AND. i_g  <=  row_length_g .AND.                    &
           j_g  >=  1 .AND. j_g  <=  rows_g ) THEN      

        ! Compute the global index of searched point
        p_g = (j_g-1) * row_length_g + i_g

        ! Count points
        radius(5)=radius(5)+1
        !
        ! Test model ppn at searched point, 
        !
        IF ( tot_pr(i_l,j_l)  >=  (epsilon_lhn * anal_pr) ) THEN
          ratio = tot_pr(i_l,j_l) / anal_pr
          IF (ratio  >   1.0) ratio = 1.0/ratio
          !
          ! Keep record of best match at this range
          !
          IF (ratio  >   best) THEN
            best    = ratio
            near    = p_g
            IF (.NOT. l_found .AND. lhn_diag) THEN
              IF (jr == 1) radius(2)=radius(2)+1
              IF (jr == 2) radius(3)=radius(3)+1
              IF (jr >= 3) radius(4)=radius(4)+1
            END IF
            l_found = .TRUE.
          END IF       ! Ratio test

        END IF         ! Ppn rate test
      END IF           ! Within domain test
    END DO             ! JN
  END IF               ! Found one test
END DO                 ! JR

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lhn_search_1a
END MODULE lhn_search_1a_mod

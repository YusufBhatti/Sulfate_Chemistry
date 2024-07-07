! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE LHN_SEARCH --------------------------------------------
!
!    Purpose : Search for suitable nearby latent heating profiles for
!              use in the LHN scheme.
!
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
MODULE lhn_search_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LHN_SEARCH_MOD'

CONTAINS

SUBROUTINE lhn_search(point,near,RANGE                            &
                           ,search,irows,icols,l_first_search     &
                           ,tot_pr,anal_pr,p_field,radius,jpts)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx
USE ac_control_mod

IMPLICIT NONE

!-----DECLARE VARIABLES
INTEGER ::    point,near,p_field,RANGE
INTEGER ::    search(4*RANGE*(RANGE+1),2)
INTEGER ::    irows, icols,jpts
INTEGER ::    radius(5)
REAL ::       tot_pr(p_field),anal_pr
LOGICAL ::    l_first_search

!-INTENT=IN---------------------------------------------------
!     P_FIELD                    - No. of points
!     POINT                      - Current point
!     RADIUS(5)                  - Diagnostic for breakdown of searches
!     RANGE                      - Search range in grid points
!     SEARCH                     - Template holding relative positions
!                                    of search points from POINT
!     TOT_PR(P_FIELD)            - Model total precip rates
!     ANAL_PR                    - Analysed ppn rate at POINT
!     IROWS                      - Number of rows in the model grid
!     ICOLS                      -    "    " columns  "    "     "
!-INTENT=INOUT-----------------------------------------------
!     L_FIRST_SEARCH             - True if this is the first search
!                                    this timestep.
!-INTENT=OUT-------------------------------------------------
!     NEAR                       - Nearest point with suitable profile
! -----------------------------------------------------------

! Local arrays and variables
INTEGER ::    xpoint       ! X coord of POINT
INTEGER ::    ypoint       ! Y   "    "   "
INTEGER ::    x            ! X coord of point to check
INTEGER ::    y            ! Y   "    "   "    "   "
INTEGER ::    p            ! Number of point at (X,Y)
INTEGER ::    jr , jn      ! Loop counters
REAL ::       ratio        ! ratio of model to obs ppn
REAL ::       best         ! closest RATIO to 1
LOGICAL ::    l_found      ! false until a suitable point is found

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LHN_SEARCH'


!
! ----------------------------------------------------------------------
!
!
!C 1.0     Set up search template if this is the first time
!C         round this timestep.
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (l_first_search) THEN
  DO jr = 1 , RANGE
    DO jn = 1 , 2*jr
      search(4*(jr-1)*jr+jn,1) = -jr-1+jn
      search(4*(jr-1)*jr+jn,2) = -jr
      search(4*(jr-1)*jr+jn+2*jr,1) = jr
      search(4*(jr-1)*jr+jn+2*jr,2) = -jr-1+jn
      search(4*(jr-1)*jr+jn+4*jr,1) = jr+1-jn
      search(4*(jr-1)*jr+jn+4*jr,2) = jr
      search(4*(jr-1)*jr+jn+6*jr,1) = -jr
      search(4*(jr-1)*jr+jn+6*jr,2) = jr+1-jn
    END DO     ! JN
  END DO       ! JR
  l_first_search=.FALSE.
END IF

!
!C 2.0     Loop round nearest points in ever increasing radius,
!C         looking for best fit to observed rain at POINT
!
!C 2.1     Initialise variables
!
l_found = .FALSE.
best = 0.0
ratio = 0.0
ypoint = INT ((point-1)/icols) + 1
xpoint = point - ( (ypoint-1) * icols )
!
!C 2.2     Loop over possible ranges
!
DO jr = 1 , RANGE
  IF (.NOT. l_found) THEN
    !
    !C  If not found one yet then loop over points at this range
    !
    DO jn = 1 , 8*jr
      !
      !C  Make sure search point is within the grid.
      !C    for global code this means P from 1 to number of pts,
      !C    for limited area code, also need to make sure X within limits.
      !
      x = xpoint + search( 4*jr*(jr-1) +jn , 1)
      y = ypoint + search( 4*jr*(jr-1) +jn , 2)
      IF ( x  >=  1 .AND. x  <=  icols .AND.                      &
           y  >=  1 .AND. y  <=  irows ) THEN      
        p = (y-1) * icols + x
        ! count points
        radius(5)=radius(5)+1
        !
        !C  Test model ppn at point, P
        !
        IF ( tot_pr(p)  >=  (epsilon_lhn * anal_pr) ) THEN
          ratio = tot_pr(p) / anal_pr
          IF (ratio  >   1.0) ratio = 1.0/ratio
          !
          !C  Keep record of best match at this range
          !
          IF (ratio  >   best) THEN
            best    = ratio
            near    = p
            IF (.NOT. l_found .AND. lhn_diag) THEN
              IF (jr == 1) radius(2)=radius(2)+1
              IF (jr == 2) radius(3)=radius(3)+1
              IF (jr >= 3) radius(4)=radius(4)+1
            END IF
            l_found = .TRUE.
          END IF       ! Ratio test
          ! Diagnostics of where pt is found
        END IF         ! Ppn rate test
      END IF           ! Within domain test
    END DO      ! JN
  END IF             ! Found one test
END DO          ! JR
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lhn_search
END MODULE lhn_search_mod

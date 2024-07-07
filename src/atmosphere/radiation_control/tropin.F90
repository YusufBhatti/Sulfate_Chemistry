! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
MODULE tropin_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TROPIN_MOD'
CONTAINS

SUBROUTINE tropin(t, Exner_rho_levels, Exner_theta_levels,        &
                  row_length, Rows, p_levels, off_x, off_y,       &
                  at_extremity,true_latitude,height_theta,        &
                  min_trop_level, max_trop_level, it)


!    SUBROUTINE TROPIN------------------------------------------------
!
!    Purpose:  Finds the tropopause & returns index of where it is.
!        Not suitable for single-column model use, as it does
!               horizontal filling-in, & includes dynamical allocation.
!
!      Based on routine TROP, but
!        1) taking temperature rather than theta as input
!        2) rather than returning the temperature, pressure and height
!        of a continuously-varying tropopause found by extrapolating
!        lapse rates, it returns the index of the adjacent layer
!        boundary (N meaning the bottom of the Nth model layer)
!        3) "filling in" where no tropopause is found.
!
!      Note that the definition of "tropopause" matches TROP, which
!        is not quite the WMO definition, though the same critical
!        lapse rate is used (unless planet_constants_mod is altered).
!        For details of the interpolation assumptions, see UMDP S1
!        section 3.2.2 or Swinbank & Wilson (1990: SRFRTN 48).
!        Any physical changes to one routine should be considered for
!                                           mirroring in the other.
!
!

USE planet_constants_mod, ONLY: g, cp, kappa, lapse_trop
USE mpp_conf_mod, ONLY: swap_field_is_scalar
USE um_parparams, ONLY: nodomain, pnorth, peast, psouth, pwest
USE conversions_mod, ONLY: pi
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Field_Types

USE model_domain_mod, ONLY: model_type, mt_single_column

IMPLICIT NONE


LOGICAL ::                                                        &
  at_extremity(4)  ! Indicates if this processor is at north, sout
                   ! east or west of the processor grid

INTEGER ::                                                        &
       !, INTENT (IN)
     row_length,                                                  &
                                  !   Number of points per row
     rows,                                                        &
                                  !   Number of rows
     p_levels,                                                    &
                                  !   Number of model levels
     min_trop_level,                                              &
     max_trop_level,                                              &
     off_x, off_y                 ! halo sizes
!     ! Limits on where the tropopause may be deemed to be - between
!     !  the MIN_TROP_LEVELth and MAX_TROP_LEVELth layers (with the
!     !  convention used here for layer boundaries, the actual index
!     !  IT returned has MIN_TROP_LEVEL < IT =< MAX_TROP_LEVEL.)

REAL ::                                                           &
    !, INTENT(IN)
!  Temperature at layer centres
          t(row_length, rows, p_levels)                                &
     ,    EXNER_rho_levels(1-off_x:row_length+off_x,                   &
                           1-off_y:rows+off_y,p_levels)                &
     ,    EXNER_theta_levels(1-off_x:row_length+off_x,                 &
                           1-off_y:rows+off_y,p_levels)                &
     ,    true_latitude(1,1)                                           &
                                          ! Latitude in radians
                                          ! for SCM
     ,    height_theta(1,1,0:p_levels)    ! For SCM


INTEGER ::                                                        &
       !, INTENT (OUT)
     it(row_length, rows)
!     ! Integer indexing the tropopause, taken to be @ a layer boundary
!     !   with the convention that N means the bottom of layer N.

! Workspace usage:-----------------------------------------------------

INTEGER ::                                                        &
       !, INTENT (OUT)
     IT_work(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

REAL :: lapse_rate(row_length, rows,                                 &
                min_trop_level+1:max_trop_level+1)
!     ! Lapse rate between layer centres

LOGICAL :: ltrop(1-off_x:row_length+off_x, 1-off_y:rows+off_y)
!     !  Logical array which indicates whether we are still seeking a
!     !    tropopause (if not, it has already been found lower down)

! ---------------------------------------------------------------------
! Define local variables:----------------------------------------------

INTEGER :: i, j, k, j0f, j1f,                                        &
                                    ! Loopers over level & point
     point,                                                       &
                      ! Point counter at ends of rows
     kp1,                                                         &
!     !  K+1, except where this would cause out-of-bounds reference
           nneigh,                                                      &
                            ! Number of well-defined tropopauses among
!     ! the 8 nearest neighbours of a point without one of its own
           fillin,                                                      &
!     !  Used to fill in points without a clearly-defined tropopause
           dti              ! Default tropopause index for points where
!     ! not even one nearest neighbour has a well-defined tropopause

REAL ::                                                           &
     del_exner_j,                                                 &
                      ! Differences of Exner pressure across
     del_exner_jm1,                                               &
                      !                           half layers
     denom,                                                       &
                      ! Denominator in lapse rate expression
     z_trop           ! Height in metres of default tropopause
                      ! in SCUM code.


REAL :: cp_over_g, p_exner_500, p_exner_50

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TROPIN'

!----------------------------------------------------------------------

!     ! 1.  Set up local constants and initialise arrays

! set up j0f j1f loop limits. AT north and south boundary ignore halo
! rows

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
j0f = 1
j1f = rows
IF (at_extremity(PSouth)) THEN
  j0f = 2
END IF
IF (at_extremity(PNorth)) THEN
  j1f = rows - 1
END IF
p_exner_500 = (500.0/1000.0)**kappa
p_exner_50  =  (50.0/1000.0)**kappa
cp_over_g = cp / g
dti = ( min_trop_level + max_trop_level ) / 2

ltrop = .TRUE.

!     ! Compute lapse rate between full levels: equation 3.16, UMDP S1

DO k=min_trop_level+1, MIN(max_trop_level+1,p_levels)
  DO j = 1, rows
    DO i=1, row_length
      !         ! Exner pressure difference across half layers
      del_exner_j = exner_rho_levels(i,j,k) /                     &
                    exner_theta_levels(i,j,k) - 1
      del_exner_jm1 = 1 - exner_rho_levels(i,j,k)                 &
                          / exner_theta_levels(i,j,k-1)
      !         ! Denominator
      denom = t(i,j,k-1) * del_exner_jm1                          &
            + t(i,j,k) * del_exner_j
      !         !  Lapse rate between level k-1 and k
      lapse_rate(i,j,k) = ( t(i,j,k-1) - t(i,j,k) ) /             &
                          ( cp_over_g * denom )
    END DO
  END DO
END DO

!     ! 2.  Find level of tropopause, where it is well defined

DO k=min_trop_level+1, max_trop_level

  ! 'K+1' level for lapse rate test; allows K iteration up to P_LEVELS
  kp1=MIN(k+1,p_levels)

  DO j = 1, rows
    DO i=1, row_length

      !         ! Not-quite-WMO criteria for interval containing tropopause
      !         ! (where 'interval' stretches between layer centres k and k-1)

      IF ( exner_theta_levels(i,j,k-1)  >   p_exner_50 .AND.      &
           exner_theta_levels(i,j,k)  <   p_exner_500 .AND.       &
           lapse_rate(i,j,k)  <   lapse_trop .AND.                &
           lapse_rate(i,j,kp1)  <   lapse_trop .AND. ltrop(i,j) ) &
      THEN
        ltrop(i,j)=.FALSE.
        IT_work(i,j) = k
      END IF
    END DO
  END DO
END DO

! DEPENDS ON: swap_bounds
CALL Swap_Bounds(it_work, row_length, rows, 1,                    &
                 off_x, off_y, fld_type_p, swap_field_is_scalar)
! DEPENDS ON: swap_bounds
CALL Swap_Bounds(ltrop, row_length, rows, 1,                      &
                 off_x, off_y, fld_type_p, swap_field_is_scalar)

!     ! 4.  Fill in it array from it_work, where set, else
!           where the above criteria did not find a tropopause

!     !  Run through all internal points and where no tropopause was
!     !   found, set the level to the average of those found at the 8
!     !   surrounding points.  If none of these did find one, then DTI,
!     !   the middle of the permitted range, is set.
!     !   (Using integers means rounding down - probably the best thing
!     !   overall as the lack of vertical resolution to pick out a
!     !   tropopause precisely is likely to mean they'll be diagnosed
!     !   too high, if anything, at neighbouring points.  Similarly
!     !   it's not worth worrying about level number being non-linear
!     !   in height or pressure at those points where ipso facto the
!     !   result is so arbitrary.)

SELECT CASE (model_type)
CASE (mt_single_column)

  ! Set tropopause model level for SCUM.  If no clearly defined
  ! tropopause then set to default level based on latitude based
  ! function fitted to zonal mean tropopause height. Default ML
  ! are theta levels. Needs further work on intepolating the
  ! correct height.  Although currently will allow routine to be
  ! used with SCM.

  DO j= j0f, j1f
    DO i= 1, row_length
      IF ( ltrop(i,j) ) THEN
        IF (true_latitude(i,j) >= -0.393 .AND.                  &
            true_latitude(i,j) <=  0.393) THEN
          z_trop  = 16000.0 + height_theta(1,1,0)
        ELSE IF (true_latitude(i,j) >= -pi/2.0 .AND. &
                 true_latitude(i,j) <  -0.393) THEN
          z_trop  = 10**4*(2.08                                 &
                  + 1.47*true_latitude(i,j)                     &
                  + 0.39*true_latitude(i,j)**2)                 &
                  + height_theta(1,1,0)
        ELSE IF (true_latitude(i,j) >  0.393 .AND.             &
                  true_latitude(i,j) <= pi/2.0) THEN
          z_trop  = 10**4*(2.24                                 &
                  - 2.01*true_latitude(i,j)                     &
                  + 0.76*true_latitude(i,j)**2)                 &
                  + height_theta(1,1,0)
        ELSE
          ! Latitude out of range -PI/2 to PI/2
          ! Set height of tropopause to
          ! be same as that in tropics
          z_trop  = 16000.0 + height_theta(1,1,0)
        END IF

        DO k = p_levels,1,-1
          IF (z_trop > height_theta(i,j,k)) THEN
            it(i,j) = k
            EXIT
          END IF
        END DO
      ELSE
        it(i,j) = it_work(i,j)
      END IF
    END DO
  END DO

CASE DEFAULT

  ! Loop over non-halo points
  ! Don't do poles as search over near neighbours is expensive.
  ! If poles not set then set to dti

  DO j= j0f, j1f
    DO i= 1, row_length

      IF ( ltrop(i,j) ) THEN

        fillin = 0
        nneigh = 0
        IF ( .NOT. ltrop(i-1,j-1) ) THEN
          fillin = fillin + IT_work(i-1,j-1)
          nneigh = nneigh + 1
        END IF
        IF ( .NOT. ltrop(i,j-1) ) THEN
          fillin = fillin + IT_work(i,j-1)
          nneigh = nneigh + 1
        END IF
        IF ( .NOT. ltrop(i+1,j-1) ) THEN
          fillin = fillin + IT_work(i+1,j-1)
          nneigh = nneigh + 1
        END IF
        IF ( .NOT. ltrop(i-1,j) ) THEN
          fillin = fillin + IT_work(i-1,j)
          nneigh = nneigh + 1
        END IF
        IF ( .NOT. ltrop(i+1,j) ) THEN
          fillin = fillin + IT_work(i+1,j)
          nneigh = nneigh + 1
        END IF
        IF ( .NOT. ltrop(i-1,j+1) ) THEN
          fillin = fillin + IT_work(i-1,j+1)
          nneigh = nneigh + 1
        END IF
        IF ( .NOT. ltrop(i,j+1) ) THEN
          fillin = fillin + IT_work(i,j+1)
          nneigh = nneigh + 1
        END IF
        IF ( .NOT. ltrop(i+1,j+1) ) THEN
          fillin = fillin + IT_work(i+1,j+1)
          nneigh = nneigh + 1
        END IF

        IF ( nneigh  ==  0 ) THEN
          it(i,j) = dti
        ELSE
          it(i,j) = fillin / nneigh
        END IF

      ELSE
        it(i,j) = it_work(i,j)
      END IF

    END DO
  END DO

  ! South pole loop: only executed if j0f > 1
  DO j= 1, j0f - 1
    DO i= 1, row_length
      IF ( ltrop(i,j) ) THEN
        it(i,j) = dti
      ELSE
        it(i,j) = it_work(i,j)
      END IF
    END DO
  END DO

  ! North pole loop: only executed if j1f < rows
  DO j= j1f+1, rows
    DO i= 1, row_length
      IF ( ltrop(i,j) ) THEN
        it(i,j) = dti
      ELSE
        it(i,j) = it_work(i,j)
      END IF
    END DO
  END DO

END SELECT ! model_type

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE tropin
END MODULE tropin_mod

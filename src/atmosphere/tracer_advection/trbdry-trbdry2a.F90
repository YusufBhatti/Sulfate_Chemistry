! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE TRBDRY-------------------------------------------------
!
!    Purpose: Special routine to add psuedo source terms to boundary
!             data in limited area model.
!    Method:  Sets the boundary of aerosol using the
!             function BDRYV. This is then copied outward to fill
!             the halo. The function BDRYV
!             is specific to a model configuration: the current
!             version (5.2) is specific to UK MES.
!
!    Programming standard: Unified Model Documentation Paper No 3,
!                          Version 7, dated 11/3/93.
!
!
!    Arguments:---------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Tracer Advection
SUBROUTINE Trbdry(                                                &
 row_length, rows, n_rows, model_levels,                          &
 offx, offy, at_extremity,                                        &
 pressure,                                                        &
 u, v,                                                            &
 murk, timestep                                                   &
)

USE conversions_mod, ONLY: recip_pi_over_180
USE mpp_conf_mod,    ONLY: swap_field_is_scalar
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
IMPLICIT NONE

INTEGER, INTENT(IN) ::                                            &
  row_length                                                      &
                    ! Length of a model row
, rows                                                            &
                    ! Number of rows for theta,u fields
, n_rows                                                          &
                    ! Number of rows for v fields
, model_levels                                                    &
                    ! Number of model levels
, offx                                                            &
                    ! Size of "single" halo (EW direction)
, offy                                                            &
                    ! Size of "single" halo (NS direction)
, timestep          ! Timestep in seconds

LOGICAL, INTENT(IN) ::                                            &
  at_extremity(4)   ! At an edge?

REAL, INTENT(IN) ::                                               &
  u( 1 - offx : row_length + offx,                                &
                                            ! U wind
     1 - offy : rows + offy,                                      &
     model_levels )                                               &
, v( 1 - offx : row_length + offx,                                &
                                            ! V wind
     1 - offy : n_rows + offy,                                    &
     model_levels )                                               &
, pressure( 1 - offx : row_length + offx,                         &
                                            ! Pressure on
            1 - offy : rows + offy,                               &
                                            ! Rho levels
            model_levels)

REAL, INTENT(INOUT) ::                                            &
  murk( 1 - offx : row_length + offx,                             &
                                          ! Aerosol for which the
        1 - offy : rows + offy,                                   &
                                          ! boundary will be set
        model_levels )


!  Local, including SAVE'd, storage------------------------------------

!  (a) Scalars

REAL ::                                                           &
 dc, wdir, wspeed, press

REAL ::                                                           &
 BdryV        ! FUNCTION giving boundary value.

!  (b) Others.
INTEGER  ::                                                       &
 i, j, jj, k  ! Loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TRBDRY'


!-----------------------------------------------------------------------
!  Loop across Northern row.
!-----------------------------------------------------------------------
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF ( at_extremity( PNorth ) ) THEN
  DO k = 1, model_levels
    j = rows
    jj = n_rows
    DO i = 1, row_length

      IF ( v(i, jj, k) > 0.0 ) THEN       ! Outflow
        murk(i, j, k) = murk(i, j-1, k)
      ELSE                               ! Inflow
        press = pressure(i, j, k)
        wdir = ATAN2( v(i, jj, k), u(i, j, k) ) * recip_pi_over_180
        wspeed = SQRT( u(i, j, k) * u(i, j, k) +                  &
                       v(i, jj, k) * v(i, jj, k) )

        ! DEPENDS ON: bdryv
        murk(i, j, k) = Bdryv( wdir, wspeed, press )
      END IF
    END DO
  END DO
END IF



!-----------------------------------------------------------------------
!  Loop across Southern row.
!-----------------------------------------------------------------------
!
IF ( at_extremity( PSouth ) ) THEN
  DO k = 1, model_levels
    j = 1
    DO i = 1, row_length

      IF ( v(i, j, k) < 0.0 ) THEN             ! Outflow
        murk(i, j, k) = murk(i, j+1, k)
      ELSE                                    ! Inflow
        press = pressure(i, j, k)
        wdir = ATAN2( v(i, j, k), u(i, j, k) ) * recip_pi_over_180
        wspeed = SQRT( u(i, j, k) * u(i, j, k) +                  &
                       v(i, j, k) * v(i, j, k) )

        ! DEPENDS ON: bdryv
        murk(i, j, k) = Bdryv( wdir, wspeed, press )
      END IF
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
!  Loop across Western column
!-----------------------------------------------------------------------
!
IF ( at_extremity( PWest ) ) THEN
  DO k = 1, model_levels
    DO j = 1, rows

      ! jj is used for v indexing
      IF ( j > n_rows + offy ) THEN
        jj = n_rows +offy
      ELSE
        jj = j
      END IF

      i = 1

      IF ( u(i,j,k) < 0.0 ) THEN           ! Outflow
        murk(i, j, k) = murk(i+1, j, k)
      ELSE                                ! Inflow
        press = pressure(i, j, k)
        wdir = ATAN2( v(i, jj, k), u(i, j, k) ) * recip_pi_over_180
        wspeed = SQRT( u(i, j, k) * u(i, j, k) +                  &
                       v(i, jj, k) * v(i, jj, k) )

        ! DEPENDS ON: bdryv
        murk(i, j, k) = Bdryv( wdir, wspeed, press )
      END IF
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
!  Loop across Eastern column
!-----------------------------------------------------------------------
!
IF ( at_extremity( PEast ) ) THEN
  DO k = 1, model_levels
    DO j = 1, rows

      ! jj is used for v indexing
      IF ( j > n_rows + offy ) THEN
        jj = n_rows +offy
      ELSE
        jj = j
      END IF

      i = row_length

      IF ( u(i,j,k) > 0.0 ) THEN          ! Outflow
        murk(i, j, k) = murk(i-1, j, k)
      ELSE                               ! Inflow
        press = pressure(i, j, k)
        wdir = ATAN2( v(i, jj, k), u(i, j, k) ) * recip_pi_over_180
        wspeed = SQRT( u(i, j, k) * u(i, j, k) +                  &
                       v(i, jj, k) * v(i, jj, k) )

        ! DEPENDS ON: bdryv
        murk(i, j, k) = Bdryv( wdir, wspeed, press )
      END IF
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! Swap bounds to ensure halos full correctly
!-----------------------------------------------------------------------

! DEPENDS ON: swap_bounds
CALL Swap_Bounds( murk, row_length, rows, model_levels,           &
                  offx, offy, fld_type_p,swap_field_is_scalar)

!-----------------------------------------------------------------------
!  Now need to fill full extended halos with copies of the
!  calculated data. Start at the North
!-----------------------------------------------------------------------

IF ( at_extremity( PNorth ) ) THEN
  DO k = 1, model_levels
    DO j = rows + 1, rows + offy
      DO i = 1 - offx, row_length + offx
        murk(i, j, k) = murk(i, rows, k)
      END DO
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! Fill extended halos for Southern rows
!-----------------------------------------------------------------------

IF ( at_extremity( PSouth ) ) THEN
  DO k = 1, model_levels
    DO j = 1 - offy, 0
      DO i = 1 - offx, row_length + offx
        murk(i, j, k) = murk(i, 1, k)
      END DO
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! Fill extended halos for Western columns
!-----------------------------------------------------------------------

IF ( at_extremity( PWest ) ) THEN
  DO k = 1, model_levels
    DO j = 1 - offy, rows + offy
      DO i = 1 - offx, 0
        murk(i, j, k) = murk(1, j, k)
      END DO
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! Fill extended halos for Eastern columns
!-----------------------------------------------------------------------

IF ( at_extremity( PEast ) ) THEN
  DO k = 1, model_levels
    DO j = 1 - offy, rows + offy
      DO i = row_length + 1, row_length + offx
        murk(i, j, k) = murk (row_length, j, k)
      END DO
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Trbdry



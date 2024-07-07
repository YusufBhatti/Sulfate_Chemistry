! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE TRSRCE ----------------------------------------------
!
!  Purpose: Adds source increments to a single level of the aerosol
!  field.
!
!  Suitable for single-column use.
!
!  Programming standard: Unified Model Documentation Paper No 3,
!
!
!  Arguments:---------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Tracer Advection

MODULE trsrce_mod
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TRSRCE_MOD'

CONTAINS
SUBROUTINE trsrce(                                                &
 rows, row_length, offx, offy, halo_i, halo_j,                    &
 theta, q, qcl, qcf, exner, rho, tracer, srce,                    &
 level, timestep, i_hour, i_minute, amp                           &
)

USE planet_constants_mod, ONLY: kappa, c_virtual, pref, cp

USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels

USE nlsizes_namelist_mod, ONLY: model_levels

USE conversions_mod, ONLY: pi
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN)    :: rows         ! number of P/U rows
INTEGER, INTENT(IN)    :: row_length   !
INTEGER, INTENT(IN)    :: offx         ! EW size of std. halo
INTEGER, INTENT(IN)    :: offy         ! NS size of std. halo
INTEGER, INTENT(IN)    :: halo_i       ! EW extended halo
INTEGER, INTENT(IN)    :: halo_j       ! NS extended halo
INTEGER, INTENT(IN)    :: i_hour       ! Local time hour
INTEGER, INTENT(IN)    :: i_minute     ! Local time minute
INTEGER, INTENT(IN)    :: level        ! level of the tracer

REAL, INTENT(IN)       :: timestep     ! Timestep in seconds
REAL, INTENT(IN)       :: amp          ! Amplitude of diurnal
                                       ! variation of emission

REAL, INTENT(IN)       ::                      &
      theta( 1 - offx : row_length + offx,     & ! pot. temperature
             1 - offy : rows + offy,           & !
             model_levels )                    &
,     q    ( 1 - halo_i :row_length + halo_i,  & ! Q on theta
             1 - halo_j :rows + halo_j,        & ! levels
             model_levels )                    &
,     qcl  ( 1 - halo_i :row_length + halo_i,  & ! Qcl on theta
             1 - halo_j :rows + halo_j,        & ! levels
             model_levels )                    &
,     qcf  ( 1 - halo_i :row_length + halo_i,  & ! Qcf on theta
             1 - halo_j :rows + halo_j,        & ! levels
             model_levels )                    &
,     exner( 1 - offx : row_length + offx,     & ! exner on rho
             1 - offy : rows + offy,           & ! levels
             model_levels + 1)                 &
,     rho  ( 1 - offx : row_length + offx,     & ! density * r * r
             1 - offy : rows + offy,           & ! on rho levels
             model_levels )

REAL, INTENT(IN)       ::                                         &
      srce( : , : )                              ! tracer source

REAL, INTENT(INOUT)    ::                                         &
      tracer( 1 - offx : ,                     & ! level of tracer
              1 - offy : )                       ! to be updated


! Local, including SAVE'd, storage------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
REAL             :: dm
REAL             :: ts
REAL             :: thetav        ! virtual potential temperature
REAL             :: exner_ave     ! an averaged exner term
REAL             :: rho_theta     ! rho on theta level
REAL             :: rho1          ! } values of rho after the
REAL             :: rho2          ! } r-squared factor removed

REAL, PARAMETER  :: Factor = 1.0  ! Factor to multiply source term
REAL, PARAMETER  :: TZero  = 12.0 ! Time of maximum emissions

!  (b) Others.
INTEGER          :: i   ! Loop counter
INTEGER          :: j   ! Loop counter

! Error Reporting
INTEGER                      :: ErrorStatus
CHARACTER (LEN=*), PARAMETER :: RoutineName='TRSRCE'
CHARACTER (LEN=errormessagelength)           :: Cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
!  Check level number is not too large.
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (level > model_levels ) THEN
  ErrorStatus = 10
  WRITE (Cmessage,*) 'Level for tracer updating is larger ',      &
                     ' than model_levels'

  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!-----------------------------------------------------------------------
!  Subroutine structure :
!  Loop over field adding source term*timestep/level mass/unit area.
!-----------------------------------------------------------------------

!     Allow for source varying with time of day
ts = 1.0
IF (amp > 0.0) THEN
  ts = 1.0 +                                                      &
  amp * COS( (REAL(i_hour) + REAL(i_minute)/60.0 - TZero)         &
      * pi/12.0)
END IF

ts = ts * timestep * factor

IF (level < model_levels) THEN
!$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j, rho1, rho2, DM) &
!$OMP& SCHEDULE(DYNAMIC)
  DO j=1, rows
    DO i = 1, row_length
      ! Remove the r squared factor from rho before interpolation
      rho1 =  rho(i,j,level)/(r_rho_levels(i,j,level) *            &
                              r_rho_levels(i,j,level) )
      rho2 =  rho(i,j,level+1)/(r_rho_levels(i,j,level+1) *        &
                                r_rho_levels(i,j,level+1) )


      ! DM = density (interpolated on to theta levels) * delta r

      dm = rho2 * (r_theta_levels(i,j,level) -                     &
                   r_rho_levels(i,j,level) ) +                     &
           rho1 * (r_rho_levels(i,j,level+1) -                     &
                   r_theta_levels(i,j,level) )
      !
      ! Special case for lowest layer to get correct mass
      IF (level == 1) THEN
        dm = dm * (r_rho_levels(i,j,2) - r_theta_levels(i,j,0))    &
                / (r_rho_levels(i,j,2) - r_rho_levels(i,j,1))
      END IF
      !
      ! Convert DM to DRY density
      dm = dm * (1.0 -q(i,j,level)-qcl(i,j,level)-qcf(i,j,level))
      !
      tracer(i, j) = tracer(i, j) + srce(i, j) * ts/dm
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE  ! level = model_level
      !--------------------------------------------------------------------
      ! Cannot average here to get rho_theta. Hence calculate by using
      ! the equation of state vis
      !         ___r                     P0
      !         rho   = -----------------------------------------
      !                 kappa * Cp * ____________________r * theta_v
      !                              [          kappa-1 ]
      !                              [ exner ** ------- ]
      !                              [          kappa   ]
      !-------------------------------------------------------------------
!$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j, rho1, rho2, DM,         &
!$OMP& exner_ave, thetav, rho_theta) SCHEDULE(DYNAMIC)
  DO j=1, rows
    DO i = 1, row_length

      thetav = theta(i,j,level) * (1.0 +                           &
                 ( q(i,j,level) * C_Virtual )                      &
               - qcl(i,j,level)-qcf(i,j,level)  )

      exner_ave = (exner(i,j,level) ** ((kappa - 1.0)/kappa) +     &
                   exner(i,j,level+1) ** ((kappa - 1.0)/kappa))/2.0
      rho_theta = Pref/( kappa * Cp * exner_ave * thetav )

      ! rho_theta is at the top theta level. We also need the value
      ! at the top rho level. This will be rho1
      rho1 = rho(i,j,model_levels)/(r_rho_levels(i,j,model_levels)*&
                                    r_rho_levels(i,j,model_levels))

      ! rho2 will be the average of rho1 and rho_theta
      rho2 = ( rho1 + rho_theta ) * 0.5

      dm = rho2 * (r_theta_levels(i,j,level) -                     &
                   r_rho_levels(i,j,level) )

      tracer(i, j) = tracer(i, j) + srce(i, j) * ts/dm
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE trsrce
END MODULE trsrce_mod

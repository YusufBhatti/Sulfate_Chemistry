! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Apply bounds to tropospheric and non-tropospheric humidities

MODULE qlimits_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='QLIMITS_MOD'

CONTAINS

SUBROUTINE QLimits ( L_Diags,                & ! in
                     L_RmNonTropQIncs,       & ! in
                     trop_min_RH,            & ! in
                     nonTrop_max_q,          & ! in
                     nonTrop_min_q,          & ! in
                     nonTrop_max_RH,         & ! in
                     trop_min_p,             & ! in
                     trop_max_PV,            & ! in
                     nonTrop_max_p,          & ! in
                     q_orig,                 & ! in
                     pv_at_theta,            & ! in
                     theta,                  & ! in
                     p,                      & ! in
                     p_theta_levels,         & ! in
                     exner_theta_levels,     & ! in
                     q,                      & ! inout
                     qCL,                    & ! inout
                     qCF,                    & ! inout
                     area_cloud_fraction,    & ! inout
                     bulk_cloud_fraction,    & ! inout
                     cloud_fraction_liquid,  & ! inout
                     cloud_fraction_frozen )   ! inout

! Description:
!
!   A point is considered tropospheric if one of the following holds:
!
!   1. Its pressure value is greater than nonTrop_max_p.
!   2. Its potential vorticity (PV) value is no greater than trop_max_PV, and
!      its pressure value is at least trop_min_p.
!
!   If neither of these conditions holds, the point is considered
!   non-tropospheric.
!
!   For tropospheric points, we impose a lower limit trop_min_RH on the
!   relative humidities:
!
!     RH >= trop_min_RH
!
!   For non-tropospheric points, the specific humidity is limited to the range
!   [nonTrop_min_q, nonTrop_max_q], and the relative humidity is capped at
!   nonTrop_max_RH:
!
!     nonTrop_min_q <= q <= nonTrop_max_q
!     RH <= nonTrop_max_RH
!
!   An alternative for non-tropospheric points - activated via the switch
!   L_RmNonTropQIncs - is to restore the specific humidities to their original
!   values, as passed in via q_orig.
!
!   In all circumstances cloud amounts are set to zero outside the troposphere.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE atm_fields_bounds_mod

USE level_heights_mod, ONLY: &
    r_theta_levels

USE trignometric_mod, ONLY:  &
    cos_theta_latitude

USE cderived_mod, ONLY: delta_lambda, delta_phi

USE planet_constants_mod, ONLY: g

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE UM_ParVars
USE Control_Max_Sizes
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY: &
    row_length, rows, model_levels

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    l_new_qsat_acassim !Currently defaults to FALSE

IMPLICIT NONE

! Subroutine arguments:

LOGICAL, INTENT(IN)    :: L_Diags          ! Write out diagnostics?
LOGICAL, INTENT(IN)    :: L_RmNonTropQIncs ! Remove non-trop q increments?

REAL,    INTENT(IN)    :: trop_min_RH    ! Lower limit to apply to trop RH
REAL,    INTENT(IN)    :: nonTrop_max_q  ! Upper limit to apply to non-trop q
REAL,    INTENT(IN)    :: nonTrop_min_q  ! Lower limit to apply to non-trop q
REAL,    INTENT(IN)    :: nonTrop_max_RH ! Upper limit to apply to non-trop RH

REAL,    INTENT(IN)    :: trop_min_p     ! Minimum tropospheric pressure
REAL,    INTENT(IN)    :: trop_max_PV    ! Maximum tropospheric ABS(PV)
REAL,    INTENT(IN)    :: nonTrop_max_p  ! Maximum non-trop     pressure

REAL,    INTENT(IN)    :: q_orig      ( tdims%i_start : tdims%i_end,  &
                                        tdims%j_start : tdims%j_end,  &
                                                    1 : tdims%k_end )

REAL,    INTENT(IN)    :: pv_at_theta ( tdims%i_start : tdims%i_end,  &
                                        tdims%j_start : tdims%j_end,  &
                                                    1 : tdims%k_end )

REAL,    INTENT(IN)    ::                                   &
  theta                 ( tdims_s%i_start : tdims_s%i_end,  &
                          tdims_s%j_start : tdims_s%j_end,  &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(IN)    ::                                   &
  p                     ( pdims_s%i_start : pdims_s%i_end,  &
                          pdims_s%j_start : pdims_s%j_end,  &
                          pdims_s%k_start : pdims_s%k_end + 1 )

REAL,    INTENT(IN)    ::                                   &
  p_theta_levels        ( tdims_s%i_start : tdims_s%i_end,  &
                          tdims_s%j_start : tdims_s%j_end,  &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(IN)    ::                                   &
  exner_theta_levels    ( tdims_s%i_start : tdims_s%i_end,  &
                          tdims_s%j_start : tdims_s%j_end,  &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT) ::                                   &
  q                     ( tdims_l%i_start : tdims_l%i_end,  &
                          tdims_l%j_start : tdims_l%j_end,  &
                          tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                                   &
  qCL                   ( tdims_l%i_start : tdims_l%i_end,  &
                          tdims_l%j_start : tdims_l%j_end,  &
                          tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                                   &
  qCF                   ( tdims_l%i_start : tdims_l%i_end,  &
                          tdims_l%j_start : tdims_l%j_end,  &
                          tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) :: area_cloud_fraction (tdims%i_start:tdims%i_end,  &
                                               tdims%j_start:tdims%j_end,  &
                                                           1:tdims%k_end )

REAL,    INTENT(INOUT) ::                                   &
  bulk_cloud_fraction   ( tdims_l%i_start : tdims_l%i_end,  &
                          tdims_l%j_start : tdims_l%j_end,  &
                          tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                                   &
  cloud_fraction_liquid ( tdims_l%i_start : tdims_l%i_end,  &
                          tdims_l%j_start : tdims_l%j_end,  &
                          tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                                   &
  cloud_fraction_frozen ( tdims_l%i_start : tdims_l%i_end,  &
                          tdims_l%j_start : tdims_l%j_end,  &
                          tdims_l%k_start : tdims_l%k_end )

! Local variables:

INTEGER :: i, j, k, ICode

LOGICAL :: InTrop(tdims%i_start : tdims%i_end,  &
                  tdims%j_start : tdims%j_end)

LOGICAL :: Strat

REAL :: t    (tdims%i_start : tdims%i_end, &
              tdims%j_start : tdims%j_end)
REAL :: p_tmp(pdims%i_start : pdims%i_end, &
              pdims%j_start : pdims%j_end)
REAL :: q_sat(tdims%i_start : tdims%i_end, &
              tdims%j_start : tdims%j_end)
REAL :: q_new(tdims%i_start : tdims%i_end, &
              tdims%j_start : tdims%j_end)

REAL :: gb_mass                  ! Total air mass in gridbox
REAL :: total_mass               ! Total air mass
REAL :: trop_mass                ! Tropospheric air mass
REAL :: trop_frac                ! Fraction of mass in troposphere
REAL :: trop_vmass_before        ! Tropospheric     vapour mass on entry
REAL :: trop_vmass_after         ! Tropospheric     vapour mass on exit
REAL :: nonTrop_vmass_before     ! Non-tropospheric vapour mass on entry
REAL :: nonTrop_vmass_after      ! Non-tropospheric vapour mass on exit
REAL :: trop_vmass_percChange    ! Percentage change to     trop vapour mass
REAL :: nonTrop_vmass_percChange ! Percentage change to non-trop vapour mass

REAL :: stats(6, 1 : tdims%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='QLIMITS'

!- End of header ---------------------------------------------------------------

!-------------------------------------------------------------------------------
! [1]: Apply tropospheric and non-tropospheric humidity limits, and if required
!      calculate local statistics for points on this PE.
!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO k = 1, tdims%k_end

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j)                                    &
!$OMP SHARED(tdims,k,exner_theta_levels,theta,p_theta_levels,t,p_tmp,q_sat,  &
!$OMP row_length,rows,InTrop,pv_at_theta,nonTrop_max_p,trop_min_p,q_new,q,   &
!$OMP trop_max_PV,trop_min_RH,qCL,qCF,nonTrop_min_q,nonTrop_max_q,           &
!$OMP cloud_fraction_liquid,cloud_fraction_frozen,area_cloud_fraction,       &
!$OMP bulk_cloud_fraction,q_orig,L_RmNonTropQIncs,nonTrop_max_RH,            &
!$OMP l_new_qsat_acassim)
  ! Get temperature and pressure.
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      t    (i,j) = exner_theta_levels(i,j,k) * theta(i,j,k)
      p_tmp(i,j) = p_theta_levels    (i,j,k)
    END DO
  END DO
!$OMP END DO

  ! Calculate saturated specific humidity.

  IF ( l_new_qsat_acassim ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      CALL qsat_new(q_sat(:,j), t(:,j), p_tmp(:,j), row_length)
    END DO
!$OMP END DO
  ELSE
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      ! DEPENDS ON: qsat
      CALL qsat (q_sat(:,j), t(:,j), p_tmp(:,j), row_length)
    END DO
!$OMP END DO
  END IF

  ! Apply limits:
!$OMP DO  SCHEDULE(STATIC)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      ! Is this point in the troposphere?
      InTrop(i,j) = (   p_theta_levels (i,j,k)  >  nonTrop_max_p .OR.  &
                      ( p_theta_levels (i,j,k)  >= trop_min_p    .AND. &
                        ABS(pv_at_theta(i,j,k)) <= trop_max_PV         &
                      )                                                &
                    )

      IF (InTrop(i,j)) THEN
        q_new(i,j) = MAX(q(i,j,k), q_sat(i,j)*trop_min_RH)
      ELSE
        ! Remove cloud:
        qCL(i,j,k) = 0.0
        qCF(i,j,k) = 0.0
        cloud_fraction_liquid(i,j,k) = 0.0
        cloud_fraction_frozen(i,j,k) = 0.0
        area_cloud_fraction  (i,j,k) = 0.0
        bulk_cloud_fraction  (i,j,k) = 0.0
        ! Adjust q:
        IF (L_RmNonTropQIncs) THEN
          q_new(i,j) = q_orig(i,j,k)
        ELSE
          q_new(i,j) = MAX(q    (i,j,k), nonTrop_min_q)
          q_new(i,j) = MIN(q_new(i,j),   nonTrop_max_q)
          q_new(i,j) = MIN(q_new(i,j),   q_sat(i,j)*nonTrop_max_RH)
        END IF
      END IF

    END DO ! i
  END DO ! j
!$OMP END DO
!$OMP END PARALLEL

  ! Calculate local statistics:
  IF (L_Diags) THEN

    ! Initialise stats.
    total_mass           = 0.0
    trop_mass            = 0.0
    trop_vmass_before    = 0.0
    trop_vmass_after     = 0.0
    nonTrop_vmass_before = 0.0
    nonTrop_vmass_after  = 0.0

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,gb_mass) SCHEDULE(STATIC)       &
!$OMP SHARED(tdims,p,cos_theta_latitude,delta_lambda,delta_phi,             &
!$OMP r_theta_levels,g,InTrop,q_new,q,k)                                    &
!$OMP REDUCTION(+: total_mass,trop_mass,trop_vmass_before,trop_vmass_after, &
!$OMP nonTrop_vmass_before,nonTrop_vmass_after)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

          ! Approximate mass in gridbox:
        gb_mass = ( p(i,j,k) - p(i,j,k+1) )  &
                  * cos_theta_latitude(i,j)  &
                  * delta_lambda * delta_phi &
                  * r_theta_levels(i,j,k)**2 &
                  / g

        total_mass = total_mass + gb_mass

        IF (InTrop(i,j)) THEN
          trop_mass            = trop_mass            + gb_mass
          trop_vmass_before    = trop_vmass_before    + gb_mass * q    (i,j,k)
          trop_vmass_after     = trop_vmass_after     + gb_mass * q_new(i,j)
        ELSE
          nonTrop_vmass_before = nonTrop_vmass_before + gb_mass * q    (i,j,k)
          nonTrop_vmass_after  = nonTrop_vmass_after  + gb_mass * q_new(i,j)
        END IF

      END DO
    END DO
!$OMP END PARALLEL DO

    stats(1,k) = total_mass
    stats(2,k) = trop_mass
    stats(3,k) = trop_vmass_before
    stats(4,k) = trop_vmass_after
    stats(5,k) = nonTrop_vmass_before
    stats(6,k) = nonTrop_vmass_after

  END IF ! (L_Diags)

  ! Update q:
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j) SCHEDULE(STATIC) &
!$OMP SHARED(tdims,q,q_new,k)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
        q(i,j,k) = q_new(i,j)
    END DO
  END DO
!$OMP END PARALLEL DO

END DO ! k

!-------------------------------------------------------------------------------
! [2]: Obtain and write global diagnostics.
!-------------------------------------------------------------------------------

IF (L_Diags) THEN

  ! Convert local stats into global stats:
  CALL gcg_rsumr (model_levels*6, gc_all_proc_group, ICode, stats)

  WRITE(umMessage,*) ''
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) '<><><><><><><><><><><><><><>'
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) 'Start of QLimits diagnostics'
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) ''
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,'(A,ES11.4)') &
    ' Lower limit to apply to trop RH:     ', trop_min_RH
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,'(A,ES11.4)') &
    ' Upper limit to apply to non-trop q:  ', nonTrop_max_q
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,'(A,ES11.4)') &
    ' Lower limit to apply to non-trop q:  ', nonTrop_min_q
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,'(A,ES11.4)') &
    ' Upper limit to apply to non-trop RH: ', nonTrop_max_RH
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,'(A,ES11.4)') &
    ' Minimum tropospheric pressure (hPa): ', trop_min_p    * 100.0
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,'(A,ES11.4)') &
    ' Maximum tropospheric ABS(PV) (PVU):  ', trop_max_PV   * 1.0e06
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,'(A,ES11.4)') &
    ' Maximum non-trop     pressure (hPa): ', nonTrop_max_p * 100.0
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) ''
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) 'Layer-by-layer diags:'
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) ''
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) '        Fraction of  Tropospheric   % change to '// &
                                 '  Non-trop       % change to'
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) '        air mass in  vapour mass    tropospheric'// &
                                 '  vapour mass    non-trop'
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) ' layer  troposphere  on entry (kg)  vapour mass '// &
                                 '  on entry (kg)  vapour mass'
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) ' -----  -----------  -------------  ------------'// &
                                 '  -------------  ------------'
  CALL umPrint(umMessage,src='qlimits')

  DO k = tdims%k_end, 1, -1
    trop_frac            = stats(2,k) / stats(1,k)
    trop_vmass_before    = stats(3,k)
    nonTrop_vmass_before = stats(5,k)
    IF (stats(3,k) /= 0.0) THEN
      trop_vmass_percChange    = 100.0 * (stats(4,k) - stats(3,k)) / stats(3,k)
    ELSE
      trop_vmass_percChange    = 0.0
    END IF
    IF (stats(5,k) /= 0.0) THEN
      nonTrop_vmass_percChange = 100.0 * (stats(6,k) - stats(5,k)) / stats(5,k)
    ELSE
      nonTrop_vmass_percChange = 0.0
    END IF
    WRITE(umMessage,'(I6,ES13.3,2(ES15.5,ES14.4))')                       &
      k, trop_frac, trop_vmass_before,    trop_vmass_percChange,   &
                    nonTrop_vmass_before, nonTrop_vmass_percChange
    CALL umPrint(umMessage,src='qlimits')
  END DO

  WRITE(umMessage,*) ''
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) 'End of QLimits diagnostics'
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) '<><><><><><><><><><><><><>'
  CALL umPrint(umMessage,src='qlimits')
  WRITE(umMessage,*) ''
  CALL umPrint(umMessage,src='qlimits')

END IF ! (L_Diags)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN


END SUBROUTINE QLimits
END MODULE qlimits_mod

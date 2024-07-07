! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine pws_max_winds_topbase_mod ---------------------------------------
!
!   Purpose: Calculates maxwind related diags, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE pws_max_winds_topbase_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_MAX_WINDS_TOPBASE_MOD'

CONTAINS

SUBROUTINE pws_max_winds_topbase(u,v,p, numlevs,                    &
                                 icode,cmessage)

USE atm_fields_bounds_mod
USE pws_diags_mod,    ONLY:                                         &
                    pws_max_wind_ub,pws_max_wind_vb,                &
                    pws_max_wind_pb,pws_max_wind_base,              &
                    pws_max_wind_top,pws_max_wind_icao,             &
                    flag_max_wind_base,flag_max_wind_top,           &
                    flag_max_wind_icao


USE missing_data_mod, ONLY: imdi, rmdi
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE uc_to_ub_mod, ONLY: uc_to_ub
USE vc_to_vb_mod, ONLY: vc_to_vb
USE pc_to_pb_mod, ONLY: pc_to_pb

USE nlsizes_namelist_mod,  ONLY: global_row_length

USE icao_ht_fc_mod, ONLY:ICAOHeight_FC

USE um_parcore, ONLY: mype

#if defined(VECTOR)
USE MaxWindSplineV_mod, ONLY: MaxWindSplineV
#else
USE MaxWindSpline_mod, ONLY: MaxWindSpline
#endif

IMPLICIT NONE

REAL, INTENT(IN) :: u(udims_s%i_start:udims_s%i_end,                &
                    udims_s%j_start:udims_s%j_end,                  &
                    udims_s%k_start:udims_s%k_end)

REAL, INTENT(IN) :: v(vdims_s%i_start:vdims_s%i_end,                &
                    vdims_s%j_start:vdims_s%j_end,                  &
                    vdims_s%k_start:vdims_s%k_end)

REAL, INTENT(IN) :: p(pdims_s%i_start:pdims_s%i_end,                &
                    pdims_s%j_start:pdims_s%j_end,                  &
                    pdims_s%k_start:pdims_s%k_end+1)

INTEGER, INTENT(IN) :: numlevs


INTEGER, INTENT(OUT) :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PWS_MAX_WINDS_TOPBASE'


! local fields
REAL :: prho(pdims_s%i_start:pdims_s%i_end,                         &
                 pdims_s%j_start:pdims_s%j_end,                     &
                 pdims_s%k_start:pdims_s%k_end)


REAL :: local_ub(udims%i_start:udims%i_end,                         &
                 vdims%j_start:vdims%j_end,                         &
                 udims%k_start:udims%k_end)

REAL :: local_vb(udims%i_start:udims%i_end,                         &
                 vdims%j_start:vdims%j_end,                         &
                 vdims%k_start:vdims%k_end)

REAL :: local_pb(udims%i_start:udims%i_end,                         &
                 vdims%j_start:vdims%j_end,                         &
                 pdims%k_start:pdims%k_end)

INTEGER :: MaxWLev(udims%i_start:udims%i_end,                       &
                   vdims%j_start:vdims%j_end)

REAL    :: MaxW   (udims%i_start:udims%i_end,                       &
                   vdims%j_start:vdims%j_end)

REAL :: WindSpd_tb(udims%i_start:udims%i_end,                       &
                 vdims%j_start:vdims%j_end)

REAL :: maxwindtop(udims%i_start:udims%i_end,                       &
                 vdims%j_start:vdims%j_end)

REAL :: maxwindbase(udims%i_start:udims%i_end,                      &
                 vdims%j_start:vdims%j_end)


INTEGER, PARAMETER :: kmax=5            ! total levels needed for spline
INTEGER, PARAMETER :: khalf=(kmax+1)/2  ! number of levels each side
INTEGER, PARAMETER :: ninc=16           ! number of increments used


! Local Variables for the vectorised version
#if defined(VECTOR)
REAL :: Uinc(udims%i_start:udims%i_end,      &
             vdims%j_start:vdims%j_end,2*ninc)
REAL :: Vinc(udims%i_start:udims%i_end,      &
             vdims%j_start:vdims%j_end,2*ninc)
REAL :: Pinc(udims%i_start:udims%i_end,      &
             vdims%j_start:vdims%j_end,2*ninc)
REAL :: Winc(udims%i_start:udims%i_end,      &
             vdims%j_start:vdims%j_end,2*ninc)
INTEGER :: MaxWinc(udims%i_start:udims%i_end,                       &
                   vdims%j_start:vdims%j_end)
#else
REAL    :: Uinc(2*ninc)
REAL    :: Vinc(2*ninc)
REAL    :: Pinc(2*ninc) ! u,v and p at increments
REAL    :: Winc(2*ninc)
INTEGER :: MaxWinc(1)
#endif

INTEGER :: i,j,k


INTEGER :: LevBtm, LevTop
INTEGER :: Level_Below, Level_Above
INTEGER :: Lev_Below_MaxW (udims%i_start:udims%i_end,               &
                           vdims%j_start:vdims%j_end)
INTEGER :: Lev_Above_MaxW (udims%i_start:udims%i_end,               &
                           vdims%j_start:vdims%j_end)
REAL :: ub,vb, ws, WS_Diff
REAL :: Min_Diff
REAL :: WindSpd

! Threshold in knots for Max Wind Top and Base Fields
REAL, PARAMETER :: Threshold_Knots = 80.0
REAL, PARAMETER :: Conv_knots_ms = 0.5418
REAL, PARAMETER :: Threshold_ms = Threshold_Knots * Conv_knots_ms

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode=0

!interpolate input fields to B grid versions.

  CALL  uC_to_uB(                                                   &
       u(udims%i_start:udims%i_end,                                 &
         udims%j_start:udims%j_end,                                 &
         udims%k_start:udims%k_end),                                &
  udims%i_len,udims%j_len,vdims%j_len,udims%k_len,                  &
  udims_s%halo_i,udims_s%halo_j,                                    &
  local_ub)


  CALL  vC_to_vB(                                                   &
       v(vdims%i_start:vdims%i_end,                                 &
         vdims%j_start:vdims%j_end,                                 &
         vdims%k_start:vdims%k_end),                                &
  pdims%j_len,pdims%i_len,vdims%j_len,vdims%k_len,                  &
  vdims_s%halo_i,vdims_s%halo_j,                                    &
  global_row_length,local_vb,local_ub)

  CALL  pc_to_pb(                                                   &
         p(pdims%i_start:pdims%i_end,                               &
         pdims%j_start:pdims%j_end,                                 &
         pdims%k_start:pdims%k_end),                                &
  pdims%i_len,pdims%j_len,vdims%j_len,pdims%k_len,                  &
  pdims_s%halo_i,pdims_s%halo_j,                                    &
  local_pb)


!--------------------------
! Max wind calc code

LevBtm = khalf
LevTop = numlevs +1 - khalf


! initialise MaxW and MaxWLev
DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    MaxWLev(i,j) = 0
    MaxW   (i,j) = rmdi   ! Large and -ve
  END DO
END DO

! Calculate 1st guess of max wind and level.
! Keep if pressure is in the 100-700 hPa range and greater than previous.
! Note that 1.0 hPa = 100.0 Pa
! The choice of pressure levels searched is in need of scientific review.
DO k = LevBtm, LevTop
  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      
      WindSpd = SQRT( (local_ub(i,j,k))**2 + (local_vb(i,j,k))**2)
      
      IF (local_pb(i,j,k+1) >= 10000.0) THEN
        IF (local_pb(i,j,k-1) <  70000.0) THEN
          IF (WindSpd >= MaxW(i,j)) THEN
            MaxWLev(i,j) = k
            MaxW(i,j) = WindSpd
          END IF
        END IF
      END IF
      
    END DO
  END DO
  
END DO


#if defined(VECTOR)

CALL MaxWindSplineV(udims%i_len,vdims%j_len, numlevs,              &
                    local_ub, local_vb, local_pb, MaxWLev,         &
                    Uinc, Vinc, Pinc)

DO k=1, 2*ninc
  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      IF (MaxWLev(i,j) /= 0 ) THEN
        Winc(i,j,k)  = SQRT(Uinc(i,j,k)**2 + Vinc(i,j,k)**2)
      END IF
    END DO
  END DO
END DO

MaxW   (:,:) = Winc(:,:,1)
MaxWinc(:,:) = 1
DO k=2,2*ninc
  WHERE (Winc(:,:,k) > MaxW )
    MaxW    = Winc(:,:,k)
    MaxWinc = k
  END WHERE
END DO

DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    IF ( MaxWLev(i,j) == 0 ) THEN
      ! If U or V is missing MaxWLev will contain 0
      pws_max_wind_ub(i,j) = rmdi
      pws_max_wind_vb(i,j) = rmdi
      pws_max_wind_pb(i,j) = rmdi
    ELSE
      ! Look along increments for new values for max wind and level
      pws_max_wind_ub(i,j) = Uinc(i,j,MaxWinc(i,j))
      pws_max_wind_vb(i,j) = Vinc(i,j,MaxWinc(i,j))
      pws_max_wind_pb(i,j) = Pinc(i,j,MaxWinc(i,j))

    END IF

  END DO
END DO

#else
DO i = udims%i_start,udims%i_end
  DO j = vdims%j_start,vdims%j_end
    IF ( MaxWLev(i,j) == 0 ) THEN
      ! If U or V is missing MaxWLev will contain 0
      pws_max_wind_ub(i,j) = rmdi
      pws_max_wind_vb(i,j) = rmdi
      pws_max_wind_pb(i,j) = rmdi
    ELSE


      CALL MaxWindSpline (numlevs, i, j,                            &
                          local_ub, local_vb, local_pb, MaxWLev,    &
                          Uinc, Vinc, Pinc)


      ! Look along increments for new values for max wind and level
      ! This can be speed squared to remove need to use costly sqrt
      Winc(:) = Uinc(:)**2 + Vinc(:)**2
      MaxWinc = MAXLOC( Winc )
      pws_max_wind_ub(i,j) = Uinc(MaxWinc(1))
      pws_max_wind_vb(i,j) = Vinc(MaxWinc(1))
      pws_max_wind_pb(i,j) = Pinc(MaxWinc(1))
    END IF

  END DO
END DO
#endif

!-----------------
! max wind icao height
IF (flag_max_wind_icao) THEN
  CALL ICAOHeight_FC(pws_max_wind_pb, pws_max_wind_icao,            &
                      udims%i_len,vdims%j_len)
END IF


IF (flag_max_wind_base .OR. flag_max_wind_top) THEN

  !--------------------------
  ! Top base calculations

  ! Initialise MaxWindP fields for base/top to RMDI.
  ! Only points with MaxWindSpeed >= Threshold will have data

  maxwindbase(:,:) = rmdi 
  maxwindtop(:,:)  = rmdi 
  WindSpd_tb(:,:) = 0.0
  Lev_Below_MaxW(:,:) = 0
  Lev_Above_MaxW(:,:) = 0

  ! Calculate Windspeed from MaxWindU and MaxWindV

  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      WindSpd_tb(i,j) = SQRT( (pws_max_wind_ub(i,j))**2 +           &
                              (pws_max_wind_vb(i,j))**2)
    END DO
  END DO


  ! Now know pws_max_wind_ub, pws_max_wind_vb, WindSpeed, pws_max_wind_pb
  !          local_ub, local_vb, local_pb (all levels)

  ! If MaxWindSpeed < 80 knots, then local_pb Level above/below remains RMDI

  DO j= vdims%j_start,vdims%j_end
    DO i= udims%i_start,udims%i_end

      IF ( WindSpd_tb(i,j) >= Threshold_ms ) THEN

        ! --------------------
        ! Find P below & above
        ! --------------------

        DO k = numlevs,1,-1

          IF ( local_pb(i,j,k) > pws_max_wind_pb(i,j) ) THEN
            Level_Below = k
            Level_Above = k+1
            EXIT
          END IF

        END DO

        ! Now know model level below and above pws_max_wind_pb
        ! Now search for model level nearest to Max Wind Base and Top

        ! -------------
        ! Max Wind Base
        ! -------------

        Min_Diff = 10000.0 
        ! Only interested in wind speeds that differ from 80 knots by less
        ! than Min_Diff (m/s)

        DO k = Level_Below, LevBtm, -1


          ub  = local_ub(i,j,k)
          vb  = local_vb(i,j,k)
          ws = SQRT ( ub**2 + vb**2 )

          WS_Diff = ABS ( ws - Threshold_ms )

          IF (WS_Diff < Min_Diff) THEN
            Min_Diff = WS_Diff
            Lev_Below_MaxW(i,j) = k
          END IF

          IF (ws < Threshold_ms) THEN
            EXIT
          END IF

        END DO

        ! -------------
        ! Max Wind Top
        ! -------------

        Min_Diff = 10000.0
        ! Only interested in wind speeds that differ from 80 knots by less
        ! than Min_Diff (m/s)

        DO k = Level_Above, LevTop

          ub  = local_ub(i,j,k)
          vb  = local_vb(i,j,k)
          ws = SQRT ( ub**2 + vb**2 )


          WS_Diff = ABS ( ws - Threshold_ms )

          IF (WS_Diff < Min_Diff) THEN
            Min_Diff = WS_Diff
            Lev_Above_MaxW(i,j) = k
          END IF

          IF (ws < Threshold_ms) THEN
            EXIT
          END IF

        END DO
#if defined(VECTOR)
      END IF   !  If WindSp >= Threshold

    END DO
  END DO

  ! -------------
  ! Max Wind Base
  ! -------------

  CALL MaxWindSplineV (udims%i_len,vdims%j_len, numlevs,            &
                     local_ub, local_vb, local_pb, Lev_Below_MaxW,  &
                     Uinc, Vinc, Pinc )

  ! Look along increments for new values for max wind and level
  DO k=1, 2*ninc
    DO j=vdims%j_start,vdims%j_end
      DO i=udims%i_start,udims%i_end
        IF (Lev_Below_MaxW(i,j)/=0) THEN
          Winc(i,j,k) = SQRT(Uinc(i,j,k)**2 + Vinc(i,j,k)**2)
          Winc(i,j,k) = Winc(i,j,k) - Threshold_ms
          Winc(i,j,k) = ABS ( Winc(i,j,k) )
        END IF
      END DO
    END DO
  END DO

  MaxW   (:,:) = Winc(:,:,1)
  MaxWinc(:,:) = 1
  DO k=2,2*ninc
    WHERE (Winc(:,:,k) < MaxW )
      MaxW    = Winc(:,:,k)
      MaxWinc = k
    END WHERE
  END DO

  DO j=vdims%j_start,vdims%j_end
    DO i=udims%i_start,udims%i_end
      IF (Lev_Below_MaxW(i,j)>0) THEN
        maxwindbase(i,j) = Pinc(i,j,MaxWinc(i,j))
      END IF
    END DO
  END DO

  ! ------------
  ! Max Wind Top
  ! ------------

  CALL MaxWindSplineV (udims%i_len,vdims%j_len, numlevs,            &
                     local_ub, local_vb, local_pb, Lev_Above_MaxW,  &
                     Uinc, Vinc, Pinc )

  ! Look along increments for new values for max wind and level

  DO k=1, 2*ninc
    DO j=vdims%j_start,vdims%j_end
      DO i=udims%i_start,udims%i_end
        IF (Lev_Below_MaxW(i,j)/=0) THEN
          Winc(i,j,k) = SQRT(Uinc(i,j,k)**2 + Vinc(i,j,k)**2)
          Winc(i,j,k) = Winc(i,j,k) - Threshold_ms
          Winc(i,j,k) = ABS ( Winc(i,j,k) )
        END IF
      END DO
    END DO
  END DO

  MaxW   (:,:) = Winc(:,:,1)
  MaxWinc(:,:) = 1
  DO k=2,2*ninc
    WHERE (Winc(:,:,k) < MaxW )
      MaxW    = Winc(:,:,k)
      MaxWinc = k
    END WHERE
  END DO

  DO j=vdims%j_start,vdims%j_end
    DO i=udims%i_start,udims%i_end
      IF (Lev_Above_MaxW(i,j)>0) THEN
        maxwindtop(i,j) = Pinc(i,j,MaxWinc(i,j))

#else
        ! -------------
        ! Max Wind Base
        ! -------------

        CALL MaxWindSpline (numlevs, i, j, local_ub, local_vb,      &
                            local_pb, Lev_Below_MaxW,               &
                            Uinc, Vinc, Pinc )

        ! Look along increments for new values for max wind and level
        Winc(:) = SQRT(Uinc(:)**2 + Vinc(:)**2)
        Winc(:) = Winc(:) - Threshold_ms
        Winc(:) = ABS ( Winc(:) )
        ! Since we are trying to find max implies dW/dP = 0
        ! Therefore looking for smallest value of dW or (dU**2+dV**2).
        MaxWinc = MINLOC( Winc )
        ! Could be multiple min values so we just take first one.
        maxwindbase(i,j) = Pinc(MaxWinc(1))

        ! ------------
        ! Max Wind Top
        ! ------------

        CALL MaxWindSpline (numlevs, i, j, local_ub, local_vb,      &
                            local_pb, Lev_Above_MaxW,               &
                            Uinc, Vinc, Pinc )

        ! Look along increments for new values for max wind and level
        Winc(:) = SQRT(Uinc(:)**2 + Vinc(:)**2)
        Winc(:) = Winc(:) - Threshold_ms
        Winc(:) = ABS ( Winc(:) )
        ! Since we are trying to find max implies dW/dP = 0
        ! Therefore looking for smallest value of dW or (dU**2+dV**2).
        MaxWinc = MINLOC( Winc )
        ! Could be multiple min values so we just take first one.
        maxwindtop(i,j) = Pinc(MaxWinc(1))
#endif

      END IF   !  If WindSp >= Threshold

    END DO
  END DO

! now calc diags at icao heights.
CALL ICAOHeight_FC(maxwindbase, pws_max_wind_base,                  &
                   udims%i_len,vdims%j_len )
CALL ICAOHeight_FC(maxwindtop,  pws_max_wind_top,                   &
                   udims%i_len,vdims%j_len )

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_max_winds_topbase

END MODULE pws_max_winds_topbase_mod

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Generate SI dynamics height field

MODULE Rcf_generate_heights_mod

!  Subroutine Rcf_Generate_Heights - generate a height field
!
! Description:
!   Generates a height field from eta-theta and eta-rho levels and
!   orography field.
!
! Method:
!   Heights are generated for theta and rho levels on theta points.
!   Either the relevant set is passed back or, futher processing is
!   done eg to average to get heights on u or v points or to get
!   zonal points. Soil depths can also be generated.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


! Global Parameters - for methods of generating height fields

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb

INTEGER, PARAMETER    ::  height_gen_original    =  1
INTEGER, PARAMETER    ::  height_gen_smooth      =  2
INTEGER, PARAMETER    ::  height_gen_linear      =  3
INTEGER, PARAMETER    ::  height_gen_ecmwf_press = 10
INTEGER, PARAMETER    ::  height_gen_ecmwf_hybrd = 11

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GENERATE_HEIGHTS_MOD'

CONTAINS

SUBROUTINE Rcf_generate_heights( grid, orog_field,  gridcode, levelT, &
                                 heights, levelsize )

USE Rcf_Headaddress_Mod, ONLY:         &
    FH_GridStagger_A,  FH_GridStagger_C,&
    FH_GridStagger_Endgame

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Ereport_mod, ONLY: &
    Ereport

USE ozone_inputs_mod, ONLY: &
    zon_av_ozone

USE UM_ParVars, ONLY: &
    gc_proc_row_group,      &
    gc_all_proc_group,      &
    bound

USE UM_ParCore, ONLY: &
    mype,             &
    nproc

USE UM_ParParams, ONLY: &
    bc_cyclic

USE Rcf_Scatter_Field_Mod, ONLY: &
    Rcf_Scatter_Field_Real

USE Rcf_Gather_Field_Mod, ONLY: &
    Rcf_Gather_Field_Real

USE nlstcall_mod, ONLY:   &
    LTimer

USE planet_constants_mod, ONLY: g, planet_radius

USE cppxref_mod, ONLY:                      &
    ppx_atm_tall,       ppx_atm_cuall,      &
    ppx_atm_cvall,      ppx_atm_cuall,      &
    ppx_atm_cvall,      ppx_atm_ozone,      &
    ppx_atm_compressed

USE stparam_mod, ONLY:                                                         &
    st_levels_model_theta, st_levels_model_rho, st_levels_deep_soil

USE rcf_field_type_mod, ONLY: &
    field_type

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

! Arguments
TYPE (grid_type), INTENT(IN)   :: grid
TYPE (field_type),INTENT(IN)   :: orog_field
INTEGER, INTENT(IN)            :: gridcode
INTEGER, INTENT(IN)            :: levelT
INTEGER, INTENT(IN)            :: levelsize
REAL, INTENT(OUT)              :: heights( levelsize,                 &
                                    0 : grid % model_levels + 1)

! Local Data
REAL                           :: denom    ! denominator in relax
REAL                           :: etk_etbl ! quantaties pulled out of
REAL                           :: erk_etbl ! loop for optimisation
REAL                           :: r_ref_theta( grid % model_levels )
REAL                           :: r_ref_rho  ( grid % model_levels )
REAL                           :: r_theta_levels( grid % loc_p_field, &
                                               0 : grid % model_levels )
REAL                           :: r_rho_levels( grid % loc_p_field,   &
                                                grid % model_levels + 1)
REAL, ALLOCATABLE              :: r_level(:,:)
REAL, ALLOCATABLE              :: u_level(:,:)
REAL, ALLOCATABLE              :: v_level(:,:)

INTEGER                        :: istat
INTEGER                        :: i
INTEGER                        :: j
INTEGER                        :: k
INTEGER                        :: kend       ! loop limit
INTEGER                        :: div
INTEGER                        :: rem
INTEGER                        :: pe
INTEGER                        :: ozone_lev  ! local ozone level
INTEGER                        :: first_u    ! First u point between 2 P points
INTEGER                        :: last_u     ! Last u point between 2 P points
INTEGER                        :: lone_u     ! U point if cyclic
INTEGER                        :: lone_u_p   ! Other P point to use if cyclic
INTEGER                        :: first_v    ! First v point between 2 P points
INTEGER                        :: last_v     ! Last v point between 2 P points
INTEGER                        :: lone_v     ! V point if cyclic
INTEGER                        :: lone_v_p   ! Other P point to use if cyclic
LOGICAL                        :: l_wind_on_p ! Wind values are on P points.


INTEGER                        :: ErrorStatus
CHARACTER (LEN=*), PARAMETER   :: RoutineName = 'RCF_GENERATE_HEIGHTS'
CHARACTER (LEN=errormessagelength)             :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (LTimer) CALL Timer( RoutineName, 3)

! Check orography is actually allocated.  Orography can be unavailable.
IF (.NOT. ALLOCATED( orog_field % data ) ) THEN
  Cmessage = 'Need to have orography before height can be generated.'
  ErrorStatus = 10
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF


! Check to see if eta values are set.
IF (.NOT. ASSOCIATED( grid % eta_theta_levels ) .OR. &
    .NOT. ASSOCIATED( grid % eta_rho_levels ) ) THEN
  Cmessage = 'Need to set eta values before heights can be generated'
  ErrorStatus = 11
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! Initialise heights for safety
heights( :, :) = 0.0

SELECT CASE(grid % grid_stagger)
CASE (FH_Gridstagger_C)
  first_u       = 1
  last_u        = grid % glob_u_row_length - 1
  lone_u        = grid % glob_u_row_length
  lone_u_p      = 1
  first_v       = 1
  last_v        = grid % glob_p_rows - 1
  lone_v        = grid % glob_v_rows
  lone_v_p      = 1
  l_wind_on_p   = .FALSE.
CASE (FH_Gridstagger_Endgame)
  first_u       = 2
  last_u        = grid % glob_u_row_length
  lone_u        = 1
  lone_u_p      = grid % glob_p_row_length
  first_v       = 2
  last_v        = grid % glob_p_rows
  lone_v        = 1
  lone_v_p      = grid % glob_p_rows
  l_wind_on_p   = .FALSE.
CASE (FH_Gridstagger_A)
  first_u       = 1
  last_u        = grid % glob_u_row_length
  lone_u        = 0
  lone_u_p      = 0
  first_v       = 1
  last_v        = grid % glob_v_rows
  lone_v        = 0
  lone_v_p      = 0
  l_wind_on_p   = .TRUE.
CASE DEFAULT
  ! Print error if grid staggering is not catered for in this routine.
  Cmessage = 'Grid staggering method is not supported.'
  ErrorStatus= 12
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END SELECT




!------------------------------------------------------------------
! Setup basic theta and rho heights - 5 methods currently
! available for this
!------------------------------------------------------------------

SELECT CASE (grid % height_gen_method)
CASE ( height_gen_original)
  CALL Rcf_Gen_Height_Original( grid, orog_field % DATA(:,1), levelT,      &
                                r_theta_levels, r_rho_levels )

CASE ( height_gen_smooth)
  CALL Rcf_Gen_Height_Smooth( grid, orog_field % DATA(:,1), levelT,        &
                              r_theta_levels, r_rho_levels )

CASE ( height_gen_linear)
  CALL Rcf_Gen_Height_Linear( grid, orog_field % DATA(:,1), levelT,        &
                              r_theta_levels, r_rho_levels )

CASE ( height_gen_ecmwf_press)
  CALL Rcf_Gen_Height_Pressure( grid, orog_field % DATA(:,1), levelT,      &
                                r_theta_levels, r_rho_levels )

CASE ( height_gen_ecmwf_hybrd)
  CALL Rcf_Gen_Height_Hybrid( grid, orog_field % DATA(:,1), levelT,        &
                                r_theta_levels, r_rho_levels )

CASE DEFAULT
  ErrorStatus = 15
  Cmessage = 'Height generation method is not recognised'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

!------------------------------------------------------------------
! Decide which height field we need here
!------------------------------------------------------------------

SELECT CASE ( gridcode )
  !------------------------------------------------------------------
  ! Theta grid
  !------------------------------------------------------------------
CASE ( ppx_atm_tall )
  IF ( levelT == st_levels_model_theta ) THEN
    ! We want r_theta_levels
    heights( :, 0 : grid % model_levels ) = r_theta_levels( :, :)
  ELSE IF ( levelT == st_levels_model_rho ) THEN
    ! Rho levels
    heights( :, 1 : grid % model_levels + 1) = r_rho_levels( :, :)
  ELSE
    ! Don't know what to do with neither theta nor rho levels
    WRITE(umMessage,*) 'levelT = ', levelT
    CALL umPrint(umMessage,src='rcf_generate_heights_mod')
    Cmessage = 'LevelT is neither full nor half'
    ErrorStatus = 20
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  !-------------------------------------------------------------------
  ! U and V grids
  !-------------------------------------------------------------------
CASE (ppx_atm_cuall, ppx_atm_cvall )
  ALLOCATE( r_level( grid % glob_p_row_length, grid % glob_p_rows ) )
  ALLOCATE( u_level( grid % glob_u_row_length, grid % glob_u_rows ) )
  ALLOCATE( v_level( grid % glob_v_row_length, grid % glob_v_rows ) )

  ! U and V need to do calculations on P fields (linear average)
  ! So do Gather 1 level per PE and then Scatter again.
  div = grid % model_levels / nproc
  rem = MOD( grid % model_levels, nproc )
  pe = 0

  DO i = 1, div
    DO j = ((i-1) * nproc) + 1, i * nproc
      ! Will gather level j on PE pe
      CALL Rcf_Gather_Field_Real( r_rho_levels( :, j), r_level, &
                             grid % loc_p_row_length,      &
                             grid % loc_p_rows,            &
                             grid % glob_p_row_length,     &
                             grid % glob_p_rows, pe,       &
                             gc_all_proc_group )

      pe = pe + 1
      IF (pe == nproc) pe = 0
    END DO

    ! U and V need to be done seperately
    ! Here do U
    IF (gridcode == ppx_atm_cuall) THEN

      IF (l_wind_on_p) THEN
        DO k = 1, grid % glob_u_rows
          DO j = first_u, last_u
            u_level(j,k) = r_level(j,k)
          END DO
        END DO
      ELSE
        DO k = 1, grid % glob_u_rows
          DO j = first_u, last_u
            u_level(j,k) = (r_level(j-first_u+1,k) + &
                            r_level(j-first_u+2,k)) * 0.5
          END DO

          ! If cyclic take the correct average, otherwise just the
          ! nearest point for the `extra' u point

          IF (bound(1) == bc_cyclic) THEN
            u_level( lone_u,k) = (r_level(lone_u,k) + &
                     r_level( lone_u_p,k) ) * 0.5
          ELSE
            IF (first_u > 1) THEN
              DO j = 1, first_u - 1
                u_level( j,k) =  r_level( j,k)
              END DO
            END IF
            IF (last_u < grid % glob_u_row_length) THEN
              DO j = grid % glob_u_row_length, last_u+1, -1
                u_level( j,k) =  r_level( j,k)
              END DO
            END IF
          END IF
        END DO
      END IF

      DO j = ((i-1) * nproc) + 1, i * nproc

        CALL Rcf_Scatter_Field_Real( heights(:, j), u_level,       &
                                grid % loc_u_row_length,      &
                                grid % loc_u_rows,            &
                                grid % glob_u_row_length,     &
                                grid % glob_u_rows, pe,       &
                                gc_all_proc_group )


        pe = pe + 1
        IF (pe == nproc) pe = 0
      END DO
    END IF

    ! Now do V
    IF (gridcode == ppx_atm_cvall) THEN

      IF (l_wind_on_p) THEN
        DO k = first_v, last_v
          DO j = 1, grid % glob_v_row_length
            v_level(j,k) = r_level(j,k)
          END DO
        END DO
      ELSE
        DO k = first_v, last_v
          DO j = 1, grid % glob_v_row_length
            v_level(j,k) = ( r_level(j,k-first_v+1) +  &
                             r_level(j,k-first_v+2) ) *0.5
          END DO
        END DO
        IF (bound(2) == bc_cyclic) THEN
          DO j = 1, grid % glob_v_row_length
            v_level(j,lone_v) = &
                   (r_level(j,lone_v) + r_level(j,lone_v_p)) *0.5
          END DO
        ELSE
          IF (grid % global) THEN
            IF (grid % glob_v_rows > grid % glob_p_rows) THEN
              ! ENDGame set polar correctly, average values
              DO j = 1, grid % glob_v_row_length
                v_level(j,1) = &
                SUM(r_level(:,1))/grid % glob_p_row_length
                v_level(j,grid % glob_v_rows) = &
                SUM(r_level(:,grid % glob_p_rows))/grid % glob_p_row_length
              END DO
            END IF
          ELSE
            IF (first_v > 1) THEN    ! Fix ENDGame South row
              DO k = 1, first_v -1
                DO j = 1, grid % glob_v_row_length
                  v_level(j,k) = r_level(j,1)
                END DO
              END DO
            END IF
            IF (last_v < grid % glob_v_rows) THEN ! Fix ENDGame North row
              DO k = last_v + 1, grid % glob_v_rows
                DO j = 1, grid % glob_v_row_length
                  v_level(j,k) = r_level(j,grid % glob_p_rows)
                END DO
              END DO
            END IF
          END IF
        END IF
      END IF
      DO j = ((i-1) * nproc) + 1, i * nproc

        CALL Rcf_Scatter_Field_Real( heights(:, j), v_level,       &
                                grid % loc_v_row_length,      &
                                grid % loc_v_rows,            &
                                grid % glob_v_row_length,     &
                                grid % glob_v_rows, pe,       &
                                gc_all_proc_group )


        pe = pe + 1
        IF (pe == nproc) pe = 0
      END DO
    END IF
  END DO

  !-------------------------------------------------------------------
  ! There are rem levels left to process. Will do these now.
  !-------------------------------------------------------------------
  pe = 0
  DO i = 1, rem
    j = nproc * div + i
    ! Will gather level j on PE pe
    CALL Rcf_Gather_Field_Real( r_rho_levels( :, j), r_level, &
                           grid % loc_p_row_length,           &
                           grid % loc_p_rows,                 &
                           grid % glob_p_row_length,          &
                           grid % glob_p_rows, pe,            &
                           gc_all_proc_group )

    pe = pe + 1
  END DO

  ! U and V need to be done seperately
  ! Here do U
  IF (gridcode == ppx_atm_cuall) THEN
    IF (l_wind_on_p) THEN
      DO k = 1, grid % glob_u_rows
        DO j = first_u, last_u
          u_level(j,k) = r_level(j,k)
        END DO
      END DO
    ELSE
      DO k = 1, grid % glob_u_rows
        DO j = first_u, last_u
          u_level(j,k) = (r_level(j-first_u+1,k) + &
                          r_level(j-first_u+2,k)) * 0.5
        END DO

        ! If cyclic take the correct average, otherwise just the
        ! nearest point for the `extra' u point
        IF (bound(1) == bc_cyclic) THEN
          u_level( lone_u,k) = (r_level(lone_u,k) + &
                                r_level( lone_u_p,k) ) * 0.5
        ELSE
          IF (first_u > 1) THEN ! Fix ENDGame South row
            DO j = 1, first_u-1
              u_level( j,k) =  r_level( j,k)
            END DO
          END IF
          IF (last_u < grid % glob_u_row_length) THEN ! Fix ND North row
            DO j = grid % glob_u_row_length, last_u+1, -1
              u_level( j,k) = r_level( j,k)
            END DO
          END IF
        END IF
      END DO
    END IF

    pe = 0
    DO i = 1, rem
      j = nproc * div + i

      CALL Rcf_Scatter_Field_Real( heights(:, j), u_level,  &
                              grid % loc_u_row_length,      &
                              grid % loc_u_rows,            &
                              grid % glob_u_row_length,     &
                              grid % glob_u_rows, pe,       &
                              gc_all_proc_group )
      pe = pe + 1
      IF (pe == nproc) pe = 0
    END DO
  END IF

  ! Now do V
  IF (gridcode == ppx_atm_cvall) THEN

    IF (l_wind_on_p) THEN
      DO k = first_v, last_v
        DO j = 1, grid % glob_v_row_length
          v_level(j,k) = r_level(j,k)
        END DO
      END DO
    ELSE
      DO k = first_v, last_v
        DO j = 1, grid % glob_v_row_length
          v_level(j,k) = ( r_level(j,k-first_v+1) + &
                           r_level(j,k-first_v+2) ) * 0.5
        END DO
      END DO

      IF (bound(2) == bc_cyclic) THEN
        DO j = 1, grid % glob_v_row_length
          v_level(j,lone_v) = (r_level(j,lone_v) +  &
                               r_level(j,lone_v_p)) * 0.5
        END DO
      ELSE
        IF (grid % global) THEN
          IF (grid % glob_v_rows > grid % glob_p_rows) THEN
            ! Endgame grid, set pole values to average of heights.
            DO j = 1, grid % glob_v_row_length
              v_level(j,1) = &
              SUM(r_level(:,1))/grid % glob_p_row_length
              v_level(j,grid % glob_v_rows) = &
              SUM(r_level(:,grid % glob_p_rows))/grid % glob_p_row_length
            END DO
          END IF
        ELSE
          IF (first_v > 1) THEN
            DO k = 1, first_v -1
              DO j = 1, grid % glob_v_row_length
                v_level(j,k) = r_level(j,1)
              END DO
            END DO
          END IF
          IF (last_v < grid % glob_v_rows) THEN
            DO k = last_v + 1, grid % glob_v_rows
              DO j = 1, grid % glob_v_row_length
                v_level(j,k) = r_level(j,grid % glob_p_rows)
              END DO
            END DO
          END IF
        END IF
      END IF
    END IF

    pe = 0
    DO i = 1, rem
      j = nproc * div + i

      CALL Rcf_Scatter_Field_Real( heights(:, j), v_level,  &
                              grid % loc_v_row_length,      &
                              grid % loc_v_rows,            &
                              grid % glob_v_row_length,     &
                              grid % glob_v_rows, pe,       &
                              gc_all_proc_group )
      pe = pe + 1
      IF (pe == nproc) pe = 0
    END DO
  END IF

  DEALLOCATE( r_level )
  DEALLOCATE( u_level )
  DEALLOCATE( v_level )

  !--------------------------------------------------------------------
  ! Ozone grid
  !-------------------------------------------------------------------
CASE ( ppx_atm_ozone )
  IF ( zon_av_ozone ) THEN
    ! Start from 0 since we might have theta level 0 for ozone.
    DO j = 0, grid % ozone_levels

      ozone_lev = grid % model_levels - grid % ozone_levels + j

      ! Sum rows of r_theta_levels to create heights
      CALL gcg_rvecsumr( grid % loc_p_row_length,      &
                         grid % loc_p_row_length,      &
                         1, grid % loc_p_rows,         &
                         r_theta_levels(1,ozone_lev),  &
                         gc_proc_row_group,            &
                         istat, heights(1,j) )

      ! Create the mean
      DO i = 1, grid % loc_p_rows
        heights(i,j) = heights(i,j)/grid % glob_p_row_length
      END DO
    END DO
  ELSE
    heights( :, 0 : grid % ozone_levels) =                          &
      r_theta_levels( :, grid % model_levels -                      &
                         grid % ozone_levels :                      &
                         grid % model_levels )
  END IF

  !--------------------------------------------------------------------
  ! Soil grid
  !--------------------------------------------------------------------
CASE ( ppx_atm_compressed )
  IF ( levelT == st_levels_deep_soil ) THEN
    ! SOIL_DEPTHS contains the thicknesses of the soil layers
    DO i = 1, grid % sm_levels
      heights(:,i) =                        &
        0.5 * grid % soil_depths(   i   )   &
        + SUM(grid % soil_depths( 1:i-1 ))
    END DO

  ELSE
    WRITE(umMessage,*) 'Land compressed field, levelT=', levelT
    CALL umPrint(umMessage,src='rcf_generate_heights_mod')
    Cmessage = 'Unsupported grid/level '//&
        'combination for height calculation'
    ErrorStatus = 30
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

CASE DEFAULT
  WRITE(umMessage,*) 'Gridcode = ', gridcode
  CALL umPrint(umMessage,src='rcf_generate_heights_mod')
  Cmessage = 'Unsupported grid type for height calculation'
  ErrorStatus = 40
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

IF (LTimer) CALL Timer( RoutineName, 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_generate_heights


! *********************************************************************
! Rcf_Gen_Height_Original - subroutine to generate original set of
! heights.
! *********************************************************************
SUBROUTINE Rcf_Gen_Height_Original( grid, orography, levelT,      &
                                    r_theta_levels, r_rho_levels)

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE planet_constants_mod, ONLY: g, planet_radius

USE ereport_mod, ONLY: ereport

USE stparam_mod, ONLY: &
    st_levels_model_theta, st_levels_model_rho

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage
IMPLICIT NONE

! Arguments
TYPE (grid_type), INTENT(IN)   :: grid
REAL, INTENT(IN)               :: orography( grid % loc_p_field )
INTEGER, INTENT(IN)            :: levelT
REAL, INTENT(OUT)              :: r_theta_levels( grid % loc_p_field, &
                                                0 : grid % model_levels)
REAL, INTENT(OUT)              :: r_rho_levels( grid % loc_p_field,   &
                                          grid % model_levels + 1)

! Local Variables
INTEGER                        :: k       ! looper
INTEGER                        :: j       ! looper
REAL                           :: denom    ! denominator in relax
REAL                           :: etk_etbl ! quantaties pulled out of
REAL                           :: erk_etbl ! loop for optimisation
REAL                           :: r_ref_theta( grid % model_levels )
REAL                           :: r_ref_rho  ( grid % model_levels )
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_GEN_HEIGHT_ORIGINAL'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!---------------------------------------------------------------------
! Set reference height profile
!---------------------------------------------------------------------
DO k = 1, grid % model_levels
  r_ref_theta(k) = grid % eta_theta_levels(k) * grid % z_top_of_model
  r_ref_rho(k)   = grid % eta_rho_levels(k) * grid % z_top_of_model
END DO

!---------------------------------------------------------------------
! set bottom level, ie: orography
!---------------------------------------------------------------------
DO j = 1, grid % loc_p_field
  r_theta_levels(j,0) =  orography(j) + planet_radius
END DO

!---------------------------------------------------------------------
! For boundary layer levels set depth to be constant.
!---------------------------------------------------------------------
IF (levelT == st_levels_model_theta) THEN
  DO k = 1, grid % bl_levels
    DO j = 1, grid % loc_p_field
      r_theta_levels(j,k) = r_theta_levels(j,0) + r_ref_theta(k)
    END DO
  END DO
END IF

IF (levelT == st_levels_model_rho) THEN
  DO k = 1, grid % bl_levels
    DO j = 1, grid % loc_p_field
      r_rho_levels(j,k)   = r_theta_levels(j,0) + r_ref_rho(k)
    END DO
  END DO

  ! Need theta_levels at bl_levels
  k = grid % bl_levels
  DO j = 1, grid % loc_p_field
    r_theta_levels(j,k)   = r_theta_levels(j,0) + r_ref_theta(k)
  END DO
END IF

!---------------------------------------------------------------------
! For constant levels set r to be a constant on the level
!---------------------------------------------------------------------
DO k = grid % first_constant_r_rho_level, grid % model_levels
  DO j = 1, grid % loc_p_field
    r_theta_levels(j,k) = planet_radius + r_ref_theta(k)
  END DO
END DO

k= grid % first_constant_r_rho_level
DO j = 1, grid % loc_p_field
  r_rho_levels(j,k)   = planet_radius + r_ref_rho(k)
END DO

DO k = grid % first_constant_r_rho_level, grid % model_levels
  DO j = 1, grid % loc_p_field
    r_rho_levels(j,k)   = planet_radius + r_ref_rho(k)
  END DO
END DO

!---------------------------------------------------------------------
! For intermediate levels use linear relaxation to constant value.
! set orographic heights.
!---------------------------------------------------------------------
denom = (grid % eta_rho_levels(grid % first_constant_r_rho_level) -&
         grid % eta_theta_levels(grid % bl_levels) )

IF ( levelT == st_levels_model_rho) THEN
  DO k = grid % bl_levels + 1, grid % first_constant_r_rho_level - 1

    erk_etbl = ( grid % eta_rho_levels(k) -                           &
             grid % eta_theta_levels(grid % bl_levels) )

    DO j = 1, grid % loc_p_field
      r_rho_levels(j,k) =                                             &
           ( r_rho_levels(j,grid % first_constant_r_rho_level) -      &
             r_theta_levels(j,grid % bl_levels) ) *                   &
             erk_etbl /                                               &
             denom                                                    &
               +  r_theta_levels(j,grid % bl_levels)
    END DO
  END DO
END IF

IF (levelT == st_levels_model_theta) THEN
  DO k = grid % bl_levels + 1, grid % first_constant_r_rho_level - 1

    etk_etbl = ( grid % eta_theta_levels(k) -                         &
             grid % eta_theta_levels(grid % bl_levels) )

    DO j = 1, grid % loc_p_field
      r_theta_levels(j,k) =                                           &
           ( r_rho_levels(j,grid % first_constant_r_rho_level) -      &
             r_theta_levels(j,grid % bl_levels) ) *                   &
             etk_etbl /                                               &
             denom                                                    &
               +  r_theta_levels(j,grid % bl_levels)
    END DO
  END DO
END IF

!---------------------------------------------------------------------
! Derive height for extra rho level at top from top rho & theta levels
!---------------------------------------------------------------------

IF ( levelT == st_levels_model_rho) THEN
  k = grid % model_levels
  DO j = 1, grid % loc_p_field
    r_rho_levels(j,k+1) = r_theta_levels(j,k) +    &
     ( r_theta_levels(j,k) - r_rho_levels(j,k) )
  END DO
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Gen_Height_Original

! *********************************************************************
! Rcf_Gen_Height_Smooth - subroutine to generate smooth set of
! heights.
! *********************************************************************
SUBROUTINE Rcf_Gen_Height_Smooth( grid, orography, levelT,      &
                                  r_theta_levels, r_rho_levels )

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE planet_constants_mod, ONLY: g, planet_radius

USE ereport_mod, ONLY: ereport

USE stparam_mod, ONLY: &
    st_levels_model_theta, st_levels_model_rho

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage
IMPLICIT NONE

! Arguments
TYPE (grid_type), INTENT(IN)   :: grid
REAL, INTENT(IN)               :: orography( grid % loc_p_field )
INTEGER, INTENT(IN)            :: levelT
REAL, INTENT(OUT)              :: r_theta_levels( grid % loc_p_field, &
                                                0 : grid % model_levels)
REAL, INTENT(OUT)              :: r_rho_levels( grid % loc_p_field,   &
                                                grid % model_levels + 1)

! Local Variables
INTEGER                        :: k       ! looper
INTEGER                        :: j       ! looper
REAL                           :: r_ref_theta( grid % model_levels )
REAL                           :: r_ref_rho  ( grid % model_levels )

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_GEN_HEIGHT_SMOOTH'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!---------------------------------------------------------------------
! Set reference height profile
!---------------------------------------------------------------------
DO k = 1, grid % model_levels
  r_ref_theta(k) = grid % eta_theta_levels(k) * grid % z_top_of_model
  r_ref_rho(k)   = grid % eta_rho_levels(k) * grid % z_top_of_model
END DO

!---------------------------------------------------------------------
! set bottom level, ie: orography
!---------------------------------------------------------------------
DO j = 1, grid % loc_p_field
  r_theta_levels(j,0) =  orography(j) + planet_radius
END DO

!---------------------------------------------------------------------
! For constant levels set r to be a constant on the level
!---------------------------------------------------------------------
DO k = grid % first_constant_r_rho_level, grid % model_levels
  DO j = 1, grid % loc_p_field
    r_theta_levels(j,k) = planet_radius + r_ref_theta(k)
  END DO
END DO

DO k = grid % first_constant_r_rho_level, grid % model_levels
  DO j = 1, grid % loc_p_field
    r_rho_levels(j,k)   = planet_radius + r_ref_rho(k)
  END DO
END DO

!---------------------------------------------------------------------
! The rest of the levels are set with a quadratic scheme.
!---------------------------------------------------------------------
IF ( levelT == st_levels_model_rho) THEN
  DO k = 1, grid % first_constant_r_rho_level - 1
    DO j = 1, grid % loc_p_field
      r_rho_levels(j,k) = planet_radius +                              &
               grid % eta_rho_levels(k) * grid % z_top_of_model +      &
               orography(j) * (1.0 - grid % eta_rho_levels(k) /        &
           grid % eta_rho_levels(grid % first_constant_r_rho_level))**2
    END DO
  END DO
END IF

IF ( levelT == st_levels_model_theta) THEN
  DO k = 1, grid % first_constant_r_rho_level - 1
    DO j = 1, grid % loc_p_field
      r_theta_levels(j,k) = planet_radius +                            &
               grid % eta_theta_levels(k) * grid % z_top_of_model +    &
               orography(j) * (1.0 - grid % eta_theta_levels(k) /      &
           grid % eta_rho_levels(grid % first_constant_r_rho_level))**2
    END DO
  END DO
END IF

!---------------------------------------------------------------------
! Derive height for extra rho level at top from top rho & theta levels
!---------------------------------------------------------------------

IF ( levelT == st_levels_model_rho) THEN
  k = grid % model_levels
  DO j = 1, grid % loc_p_field
    r_rho_levels(j,k+1) = r_theta_levels(j,k) +    &
     ( r_theta_levels(j,k) - r_rho_levels(j,k) )
  END DO
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Gen_Height_Smooth

! *********************************************************************
! Rcf_Gen_Height_Linear - subroutine to generate heights for the
! Idealised model using a linear structure
! *********************************************************************
SUBROUTINE Rcf_Gen_Height_Linear( grid, orography, levelT,      &
                                  r_theta_levels, r_rho_levels )

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE planet_constants_mod, ONLY: planet_radius

USE stparam_mod, ONLY:      &
    st_levels_model_theta,  &
    st_levels_model_rho

IMPLICIT NONE

! Arguments
TYPE (grid_type), INTENT(IN)   :: grid
REAL, INTENT(IN)               :: orography( grid % loc_p_field )
INTEGER, INTENT(IN)            :: levelT
REAL, INTENT(OUT)              :: r_theta_levels( grid % loc_p_field, &
                                                0 : grid % model_levels)
REAL, INTENT(OUT)              :: r_rho_levels( grid % loc_p_field,   &
                                                grid % model_levels + 1)

! Local Variables
INTEGER                        :: k       ! looper
INTEGER                        :: j       ! looper
REAL                           :: r_ref_theta( grid % model_levels )
REAL                           :: r_ref_rho  ( grid % model_levels )


!---------------------------------------------------------------------
! Set reference height profile
!---------------------------------------------------------------------
DO k = 1, grid % model_levels
  r_ref_theta(k) = grid % eta_theta_levels(k) * grid % z_top_of_model
  r_ref_rho(k)   = grid % eta_rho_levels(k) * grid % z_top_of_model
END DO

!---------------------------------------------------------------------
! set bottom level, ie: orography
!---------------------------------------------------------------------
DO j = 1, grid % loc_p_field
  r_theta_levels(j,0) =  orography(j) + planet_radius
END DO

!---------------------------------------------------------------------
! Set all levels linearly
!---------------------------------------------------------------------
IF ( levelT == st_levels_model_rho) THEN
  DO k = 1, grid % model_levels
    DO j = 1, grid % loc_p_field
      r_rho_levels(j,k) = planet_radius + r_ref_rho(k) +               &
                            orography(j) * (1.0 - grid % eta_rho_levels(k))
    END DO
  END DO

  ! Need theta_levels at model_levels
  k = grid % model_levels
  DO j = 1, grid % loc_p_field
    r_theta_levels(j,k) = planet_radius + r_ref_theta(k) +             &
                            orography(j) * (1.0 - grid % eta_theta_levels(k))
  END DO
END IF

IF ( levelT == st_levels_model_theta) THEN
  DO k = 1, grid % model_levels
    DO j = 1, grid % loc_p_field
      r_theta_levels(j,k) = planet_radius + r_ref_theta(k) +           &
                              orography(j) * (1.0 - grid % eta_theta_levels(k))
    END DO
  END DO
END IF

!---------------------------------------------------------------------
! Derive height for extra rho level at top from top rho & theta levels
!---------------------------------------------------------------------

IF ( levelT == st_levels_model_rho) THEN
  k = grid % model_levels
  DO j = 1, grid % loc_p_field
    r_rho_levels(j,k+1) = r_theta_levels(j,k) +    &
     ( r_theta_levels(j,k) - r_rho_levels(j,k) )
  END DO
END IF

RETURN
END SUBROUTINE Rcf_Gen_Height_Linear

! *********************************************************************
! Rcf_Gen_Height_Pressure - subroutine to generate set of heights
! from ECMWF Pressure levels.
! *********************************************************************
SUBROUTINE Rcf_Gen_Height_Pressure( grid, orography, levelT,      &
                                    r_theta_levels, r_rho_levels  &
                                    )

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_GRIB_T_n_Pstar_H_Interp_Mod, ONLY:  &
    grib_tv,            &
    GRIB_Pstar,         &
    GRIB_Levels

USE planet_constants_mod, ONLY: g, planet_radius, r

USE ereport_mod, ONLY: ereport

USE stparam_mod, ONLY: st_levels_model_rho

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage
IMPLICIT NONE

! Arguments
TYPE (grid_type), INTENT(IN)   :: grid
REAL, INTENT(IN)               :: orography( grid % loc_p_field )
INTEGER, INTENT(IN)            :: levelT
REAL, INTENT(OUT)              :: r_theta_levels( grid % loc_p_field, &
                                                0 : grid % model_levels)
REAL, INTENT(OUT)              :: r_rho_levels( grid % loc_p_field,   &
                                          grid % model_levels + 1)

! Local Variables
INTEGER                        :: k       ! looper
INTEGER                        :: j       ! looper
REAL                           :: Roverg,rtemp,logpp
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_GEN_HEIGHT_PRESSURE'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------------------------------------------
! Calculate R/g
!---------------------------------------------------------------------
Roverg = r/g

!---------------------------------------------------------------------
! set bottom level, ie: orography
!---------------------------------------------------------------------
DO j = 1, grid % loc_p_field
  r_theta_levels(j,0) =  orography(j) + planet_radius
END DO

!---------------------------------------------------------------------
! Calculate heights from Pstar and T
!---------------------------------------------------------------------

!Use hydrostatic approximation and ideal gas equation to derive
!expression for difference in height between two pressure levels.
!Temperature assumed constant between two successive levels.

!Calculate first level heights from pstar
!Temperature is from first level
!Do this because no reference data at surface
DO j = 1, grid % loc_p_field
  logpp = LOG( GRIB_Pstar % DATA(j,1) / GRIB_Levels(1) )
  rtemp = Roverg * grib_tv % DATA(j,1) * logpp
  r_theta_levels(j,1) = r_theta_levels(j,0) + rtemp
  r_rho_levels(j,1)   = r_theta_levels(j,1)
END DO

!Calculate rest of the level heights from previous pressure level
!Temperature between levels is mean of temperatures at each level
DO k = 2, grid % model_levels
  logpp=0.5*Roverg*LOG( GRIB_Levels(k-1) / GRIB_Levels(k) )
  DO j = 1, grid % loc_p_field
    rtemp = logpp*(grib_tv%DATA(j,k)+grib_tv%DATA(j,k-1))
    r_theta_levels(j,k) = r_theta_levels(j,k-1) + rtemp
    r_rho_levels(j,k)   = r_theta_levels(j,k)
  END DO
END DO

!---------------------------------------------------------------------
! Derive height for extra rho level at top from top rho & theta levels
!---------------------------------------------------------------------

IF ( levelT == st_levels_model_rho) THEN
  k = grid % model_levels
  DO j = 1, grid % loc_p_field
    r_rho_levels(j,k+1) = r_theta_levels(j,k) +    &
     ( r_theta_levels(j,k) - r_rho_levels(j,k-1) )
  END DO
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Gen_Height_Pressure

! *********************************************************************
! Rcf_Gen_Height_Hybrid - subroutine to generate set of heights
! from ECMWF Model levels.
! *********************************************************************
SUBROUTINE Rcf_Gen_Height_Hybrid( grid, orography, levelT,      &
                                    r_theta_levels, r_rho_levels  &
                                    )

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_GRIB_T_n_Pstar_H_Interp_Mod, ONLY:  &
    grib_tv,            &
    GRIB_Pstar,         &
    GRIB_Levels,  &
    ak, bk

USE planet_constants_mod, ONLY: g, planet_radius, r

USE ereport_mod, ONLY: ereport

USE stparam_mod, ONLY: st_levels_model_rho

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage

IMPLICIT NONE

! Arguments
TYPE (grid_type), INTENT(IN)   :: grid
REAL, INTENT(IN)               :: orography( grid % loc_p_field )
INTEGER, INTENT(IN)            :: levelT
REAL, INTENT(OUT)              :: r_theta_levels( grid % loc_p_field, &
                                                0 : grid % model_levels)
REAL, INTENT(OUT)              :: r_rho_levels( grid % loc_p_field,   &
                                          grid % model_levels + 1)

! Local Variables
INTEGER                        :: k       ! looper
INTEGER                        :: j       ! looper
REAL                           :: Roverg,rtemp,t_tmp,p_tmp1,p_tmp2
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_GEN_HEIGHT_HYBRID'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!---------------------------------------------------------------------
! Calculate R/g
!---------------------------------------------------------------------
Roverg = r/g

!---------------------------------------------------------------------
! set bottom level, ie: orography
!---------------------------------------------------------------------
DO j = 1, grid % loc_p_field
  r_theta_levels(j,0) =  orography(j) + planet_radius
END DO

!---------------------------------------------------------------------
! Calculate heights from Pstar and T
!---------------------------------------------------------------------

!Use hydrostatic approximation and ideal gas equation to derive
!expression for difference in height between two model levels.

!Use ECMWF definition of model levels to obtain pressure at point
!on model level.

!Temperature assumed constant between two successive model levels.

!In following, GRIB_Levels(k) contains index of level k relative
!to ECMWF model level set but is a real number, hence need to use int
!to make it an integer.

!Calculate first level heights from pstar with temperature from
!first model level. Do this because no temperature data at surface.
DO j = 1, grid % loc_p_field

  !Calculate ratio of pressures of level 1 and surface
  p_tmp1 = ak(INT(GRIB_Levels(1)))/GRIB_Pstar%DATA(j,1) &
     + bk(INT(GRIB_Levels(1)))

  !Calculate thickness of first level layer (note log(p_tmp1) < 0)
  rtemp = - Roverg * grib_tv % DATA(j,1) * LOG( p_tmp1 )

  !Add thickness of layer to surface height
  r_theta_levels(j,1) = r_theta_levels(j,0) + rtemp

  !Set r_rho equal to r_theta
  r_rho_levels(j,1)   = r_theta_levels(j,1)

END DO

!Calculate rest of the level heights from previous model level
!Temperature between levels is mean of temperatures at each level
DO k = 2, grid % model_levels
  DO j = 1, grid % loc_p_field

    !Calculate pressure at level k-1
    p_tmp1 = ak(INT(GRIB_Levels(k-1))) +  &
             bk(INT(GRIB_Levels(k-1)))*GRIB_Pstar%DATA(j,1)

    !Calculate pressure at level k
    p_tmp2 = ak(INT(GRIB_Levels(k))) +  &
             bk(INT(GRIB_Levels(k)))*GRIB_Pstar%DATA(j,1)

    !Calculate mean temperature across layer
    t_tmp = 0.5*( grib_tv%DATA(j,k) + grib_tv%DATA(j,k-1) )

    !Calculate thickness of layer between levels k and k-1
    rtemp = Roverg * t_tmp * LOG( p_tmp1 / p_tmp2 )

    !Add thickness of layer to level k-1 height
    r_theta_levels(j,k) = r_theta_levels(j,k-1) + rtemp

    !Set r_rho equal to r_theta
    r_rho_levels(j,k)   = r_theta_levels(j,k)

  END DO
END DO

!---------------------------------------------------------------------
! Derive height for extra rho level at top from top rho & theta levels
!---------------------------------------------------------------------

IF ( levelT == st_levels_model_rho) THEN
  k = grid % model_levels
  DO j = 1, grid % loc_p_field
    r_rho_levels(j,k+1) = r_theta_levels(j,k) +    &
     ( r_theta_levels(j,k) - r_rho_levels(j,k-1) )
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Gen_Height_Hybrid

END MODULE Rcf_generate_heights_mod

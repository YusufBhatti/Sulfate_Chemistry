! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE diagnostics_casim_mod
! Description: Routine which outputs appropriate precipitation 
!              diagnostics to STASH from the CASIM microphysics


! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation (CASIM)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGNOSTICS_CASIM_MOD'

CONTAINS

SUBROUTINE diagnostics_casim(stashwork, t_n, q_n, qcl_n, qcf2_n,            &
                             t_inc, q_inc, qcl_inc, qcf_inc, qrain_inc,     &
                             qgraup_inc, qcf2_inc, p_layer_centres, deltaz)

!-------------------------------------------------
! List of output Variables
!------------------------------------------------- 
! Code | Items
!-------------------------------------------------
! 4    | Temperature after microphysics [K]
! 10   | Specific humidity after microphysics [kg/kg]
! 110  | Surface reflectivity (dBZ)
! 111  | Max Reflectivity in the column (dBZ)
! 112  | Reflectivity at 1km AGL (dBZ)
! 113  | Reflectivity due to graupel (dBZ)
! 114  | Reflectivity due to ice aggregates (dBZ)
! 115  | Reflectivity due to ice crystals (dBZ)
! 116  | Reflectivity due to rain (dBZ)
! 117  | Reflectivity due to liquid cloud (dBZ)
! 118  | Total Reflectivity (dBZ)
! 181  | Temperature increment due to microphysics [K]
! 182  | Humidity increment due to microphysics [kg/kg]
! 183  | qcl increment due to microphysics [kg/kg]
! 184  | qcf increment due to microphysics [kg/kg]
! 189  | qr increment due to microphysics [kg/kg]
! 190  | qgraup increment due to microphysics [kg/kg]
! 191  | qcf2 increment due to microphysics [kg/kg]
! 201  | Rain Amount [kg m-2 timestep-1]
! 202  | Snow Amount [kg m-2 timestep-1]
! 203  | Rain Rate   [kg m-2 s-1]
! 204  | Snow Rate   [kg m-2 s-1]
! 205  | Cloud liquid after microphysics [kg/kg]
! 206  | Cloud ice after microphysics [kg/kg]
! 207  | Relative Humidity after Microphysics [%]
! 208  | Relative Humidity with respect to water after Microphysics [%]
! 209  | Graupel amount [kg m-2 timestep-1]
! 212  | Graupel Rate   [kg m-2 s-1]
! 224  | Supercooled liquid water content [kg/kg]
! 240  | Homogeneous nucleation rate [kg/kg/s]
! 241  | Heterogeneous nucleation rate  [kg/kg/s]
! 243  | Deposition rate ice crystals [kg/kg/s]
! 245  | Deposition rate snow aggregates [kg/kg/s]
! 247  | Riming rate ice crystals [kg/kg/s]
! 248  | Riming rate snow aggregates [kg/kg/s]
! 250  | Snow aggregate-rain capture rate [kg/kg/s]
! 251  | Evaporation (sublimation) of melting ice crystals [kg/kg/s]
! 252  | Evaporation (sublimation) of melting snow aggregates [kg/kg/s]
! 253  | Melting rate of ice crystals [kg/kg/s]
! 254  | Melting rate of snow aggregates [kg/kg/s]
! 255  | Snow autoconversion rate [kg/kg/s]
! 256  | Snow capture of ice crystals [kg/kg/s]
! 257  | Rain autoconversion rate [kg/kg/s]
! 258  | Rain accretion rate [kg/kg/s]
! 259  | Rain evaporation rate [kg/kg/s]
! 261  | Graupel accretion rate [kg/kg/s]
! 262  | Graupel-snow capture rate [kg/kg/s]
! 263  | Graupel melting rate [kg/kg/s]
! 264  | Graupel evaporation (sublimation) rate [kg/kg/s]
! 265  | Ice crystal sedimentation rate [kg/kg/s]
! 266  | Snow aggregate sedimentation rate [kg/kg/s]
! 267  | Rain sedimentation rate [kg/kg/s]
! 268  | Graupel sedimentation rate [kg/kg/s]
! 269  | Liquid cloud sedimentation rate [kg/kg/s]
! 271  | Homogeneous freezing of rain rate [kg/kg/s]
! 302  | Snowfall amount excluding graupel [kg m-2 timestep-1]
! 304  | Snowfall rate excluding graupel [kg m-2 s-1]
! 325  | Cloud Condensation/Evaporation rate [kg/kg/s]
! 336  | Ice Cloud sedimentation rate [kg/kg/s]
! 350  | Homogeneous Freezing of Cloud Number Tendency [Num. kg-1 s-1]
! 351  | Homogeneous Freezing of Rain Number Tendency [Num. kg-1 s-1]
! 352  | Ice number tendency due to Hallett-Mossop Process [Num. kg-1 s-1]
! 353  | Ice number tendency due to Ice Nucleation [Num. kg-1 s-1]
! 354  | Tendency in ice number due to sedimentation [Num. kg-1 s-1]
! 355  | Tendency in snow number due to sedimentation [Num. kg-1 s-1] 
! 356  | Tendency in graupel number due to sedimentation [Num. kg-1 s-1]
!-------------------------------------------------

USE atm_fields_bounds_mod,  ONLY: tdims
USE submodel_mod,           ONLY: atmos_im
USE stash_array_mod,        ONLY: len_stlist, stindex, stlist, stash_levels, &
                                  num_stash_levels, sf, si

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umprintmgr,             ONLY: newline
USE um_parvars,             ONLY: at_extremity
USE timestep_mod,           ONLY: timestep
USE mphys_radar_mod,        ONLY: ref_lim
USE conversions_mod,        ONLY: zerodegc
USE gen_phys_inputs_mod,    ONLY: l_mr_physics

USE generic_diagnostic_variables, ONLY: casdiags

USE qsat_mod, ONLY: qsat_new         => qsat,                                &
                    qsat_mix_new     => qsat_mix,                            &
                    qsat_wat_new     => qsat_wat,                            &
                    l_new_qsat_mphys ! Currently defaults to FALSE

USE acp_namel_mod,         ONLY: l_ac
USE ac_diagnostics_mod,    ONLY: lsrr, lssr, tinc_ppn
USE missing_data_mod,      ONLY: rmdi
USE pws_diags_mod,         ONLY: pws_precip_sym_ls_rain,                     &
                                 pws_precip_sym_ls_snow
USE um_stashcode_mod,      ONLY: stashcode_pws_sec, stashcode_pws_precip_sym

! Dr Hook modules
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim

IMPLICIT NONE
!--------------------------------------------------------
! Subroutine argument list
!--------------------------------------------------------

REAL, INTENT(INOUT) :: stashwork(*) ! Section 4 STASH workspace


REAL, INTENT(IN) :: t_n      ( tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Temperature [K]

REAL, INTENT(IN) :: q_n      ( tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Vapour [kg kg-1]

REAL, INTENT(IN) :: qcl_n    ( tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Cloud liquid content [kg kg-1]

REAL, INTENT(IN) :: qcf2_n   ( tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Cloud ice content [kg kg-1]

REAL, INTENT(IN) :: t_inc    ( tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Temperature increment [K]

REAL, INTENT(IN) :: q_inc    ( tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Vapour increment [kg kg-1]

REAL, INTENT(IN) :: qcl_inc  ( tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Liquid cloud condensate increment [kg kg-1]

REAL, INTENT(IN) :: qcf_inc  ( tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Frozen cloud aggregate condensate increment [kg kg-1]

REAL, INTENT(IN) :: qrain_inc( tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Rain mass mixing ratio increment [kg kg-1]

REAL, INTENT(IN) :: qgraup_inc(tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Graupel mass mixing ratio increment [kg kg-1]

REAL, INTENT(IN) :: qcf2_inc ( tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Frozen cloud crystal condensate increment [kg kg-1]

REAL, INTENT(IN) ::  p_layer_centres( tdims%i_start : tdims%i_end,         &
                                      tdims%j_start : tdims%j_end,         &
                                      tdims%k_start : tdims%k_end )
! Pressure at the centre of each layer [Pa]

REAL, INTENT(IN) :: deltaz   ( tdims%i_start : tdims%i_end,                &
                               tdims%j_start : tdims%j_end,                &
                                           1 : tdims%k_end )
! Layer thickness [m]

!--------------------------------------------------------
! Local Variables
!--------------------------------------------------------
INTEGER, PARAMETER :: sect = 4 ! STASH diagnostic section
INTEGER            :: item     ! STASH item
INTEGER            :: icode    ! Error code
INTEGER            :: i        ! )
INTEGER            :: j        ! )-Loop counters
INTEGER            :: k        ! )
INTEGER            :: ji       ! )
INTEGER            :: im_index ! Internal model index

! Dr Hook related variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER      :: RoutineName='DIAGNOSTICS_CASIM' 
                                   ! routine name
CHARACTER(LEN=errormessagelength):: cmessage='' ! Error comments message


! A local work variable on the model grid for 3D fields.
REAL :: work_3d( tdims%i_start : tdims%i_end,                               &
                 tdims%j_start : tdims%j_end,                               &
                             1 : tdims%k_end )

! A local work variable on the model grid for 2D fields.
REAL :: work_2d( tdims%i_start : tdims%i_end,                               &
                 tdims%j_start : tdims%j_end)

! Local working temperature (t_n + t_inc) [K]
REAL :: t_work(  tdims%i_start : tdims%i_end,                               &
                 tdims%j_start : tdims%j_end,                               &
                             1 : tdims%k_end )

! Local working vapour mixing ratio [kg kg-1]
REAL :: q_work(  tdims%i_start : tdims%i_end,                               &
                 tdims%j_start : tdims%j_end,                               &
                             1 : tdims%k_end )

! 1 km above ground level [expressed in m]
REAL, PARAMETER :: agl_1km = 1000.0

! Temperature in a single grid box [K] - used in loops.
REAL :: t_one_gb

! Height above ground level [m] (required for radar diagnostics)
REAL :: z_agl( tdims%i_start : tdims%i_end,                                 &
               tdims%j_start : tdims%j_end )

LOGICAL :: dbz_def ( tdims%i_start : tdims%i_end,                           &
                     tdims%j_start : tdims%j_end )
! Determines whether reflectivity needs definining or not for each column


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------------
! Start of the subroutine, i.e. where the real work is done.
!----------------------------------------------------------------------------

! Initialise error code to zero.
icode = 0

! Initialise work_3d and work_2d to zero.
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      work_3d(i,j,k) = 0.0
    END DO
  END DO
END DO

DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    work_2d(i,j) = 0.0
  END DO
END DO

! Initialise im_index
im_index = 1

!--------------------------------------------------------------------------
! Item 4: Temperature after microphysics (large scale precipitation)
!--------------------------------------------------------------------------
item = 4
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = t_n(i,j,k) + t_inc(i,j,k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), work_3d,                &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                          &
       at_extremity,                                                           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1,                                        &
       atmos_im,sect,item,                                                     &
       icode,cmessage)

  IF (icode >  0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'Temperature after large scale precipitation'
  END IF
END IF ! Item 4

!--------------------------------------------------------------------------
! Item 10: Humidity after microphysics (large scale precipitation)
!--------------------------------------------------------------------------
item = 10
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = q_n(i,j,k) + q_inc(i,j,k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), work_3d,                &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                          &
       at_extremity,                                                           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1,                                        &
       atmos_im,sect,item,                                                     &
       icode,cmessage)

  IF (icode >  0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'Specific humidity after large scale precipitation'
  END IF
END IF ! Item 10

!---------------------------------------------------------------
! Item 110: Surface (model level 1 reflectivity)
!---------------------------------------------------------------
item = 110

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! First calculate dbz_surf: surface reflectivity

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      work_2d(i,j) = casdiags % dbz_tot(i,j,1)
    END DO ! i
  END DO ! j

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), work_2d,           &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                 &
       atmos_im,sect,item,                                            &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Surface radar reflectivity'
  END IF
END IF ! icode <= 0

!---------------------------------------------------------------
! Item 111: Max (composite) reflectivity in the column
!---------------------------------------------------------------
item = 111

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! Calculate Maximum Reflectivity

  work_2d(:,:) = ref_lim

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF ( casdiags % dbz_tot(i,j,k) > work_2d(i,j) ) THEN
          work_2d(i,j) = casdiags % dbz_tot(i,j,k)
        END IF
      END DO ! i
    END DO ! j
  END DO ! k

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), work_2d,           &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                 &
       atmos_im,sect,item,                                            &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Max radar reflectivity in column'
  END IF
END IF ! icode <= 0

!---------------------------------------------------------------
! Item 112: Reflectivity at 1km AGL
!---------------------------------------------------------------
item = 112

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! First calculate reflectivity at 1 km AGL using the layer thickness

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      z_agl(i,j)   = 0.0
      work_2d(i,j) = ref_lim
      dbz_def(i,j) = .TRUE.
    END DO
  END DO
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        z_agl(i,j) = z_agl(i,j) + deltaz(i,j,k)

        IF (dbz_def(i,j) .AND. z_agl(i,j) >= agl_1km ) THEN

          work_2d(i,j) = casdiags % dbz_tot(i,j,k)
          dbz_def(i,j) = .FALSE. ! Prevent overwriting by higher levels

        END IF ! dbz_def

      END DO ! i
    END DO ! j
  END DO !k

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), work_2d,           &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                 &
       atmos_im,sect,item,                                            &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Reflectivity 1km above ground level'
  END IF

END IF ! icode <= 0

!---------------------------------------------------------------
! Item 113: Reflectivity due to graupel
!---------------------------------------------------------------
item = 113

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % dbz_g,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
   WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Reflectivity due to graupel'
  END IF
END IF
      
!---------------------------------------------------------------
! Item 114:  Reflectivity due to ice aggregates
!---------------------------------------------------------------
item = 114

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % dbz_s,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
   WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Reflectivity due to ice aggregates'
  END IF
END IF

!---------------------------------------------------------------
! Item 115: Reflectivity due to ice crystals
!---------------------------------------------------------------
item = 115

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % dbz_i,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Reflectivity due to ice cloud'
  END IF
END IF

!---------------------------------------------------------------
! Item 116: Reflectivity due to rain
!---------------------------------------------------------------
item = 116

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % dbz_r,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Reflectivity due to rain'
  END IF
END IF

!---------------------------------------------------------------
! Item 117: Reflectivity due to liquid cloud
!---------------------------------------------------------------
item = 117

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % dbz_l,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Reflectivity due to liquid cloud'
  END IF
END IF

!---------------------------------------------------------------
! Item 118: Total Reflectivity
!---------------------------------------------------------------
item = 118

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % dbz_tot,                                               &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Total (3D) reflectivity'
  END IF
END IF

!--------------------------------------------------------------------------
! Item 181: Temperature Increment
!--------------------------------------------------------------------------
item = 181

IF (icode <= 0 .AND. (sf(item,sect) .OR. l_ac)) THEN
  IF (.NOT. ALLOCATED(tinc_ppn)) THEN
    ALLOCATE ( tinc_ppn(tdims%i_end*tdims%j_end,tdims%k_end) )
    tinc_ppn(1,1) = rmdi
  END IF
END IF

IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), t_inc,                  &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                          &
       at_extremity,                                                           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1,                                        &
       atmos_im,sect,item,                                                     &
       icode,cmessage)

  IF (icode >  0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'Temperature increment due to large scale precipitation'
  END IF
END IF ! Item 181

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i,ji)                                                       &
!$OMP SHARED(tdims,tinc_ppn,t_inc)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ji = (j-1)*tdims%i_end+i

        tinc_ppn(ji,k) = t_inc(i,j,k)

      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

!--------------------------------------------------------------------------
! Item 182: Humidity Increment
!--------------------------------------------------------------------------
item = 182
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), q_inc,                  &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                          &
       at_extremity,                                                           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1,                                        &
       atmos_im,sect,item,                                                     &
       icode,cmessage)

  IF (icode >  0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'humidity increment due to large scale precipitation'
  END IF
END IF ! Item 182

!--------------------------------------------------------------------------
! Item 183: qcl Increment
!--------------------------------------------------------------------------
item = 183
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), qcl_inc,                &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                          &
       at_extremity,                                                           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1,                                        &
       atmos_im,sect,item,                                                     &
       icode,cmessage)

  IF (icode >  0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'qcl increment due to large scale precipitation'
  END IF
END IF ! Item 183

!--------------------------------------------------------------------------
! Item 184: qcf Increment
!--------------------------------------------------------------------------
item = 184
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), qcf_inc,                &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                          &
       at_extremity,                                                           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1,                                        &
       atmos_im,sect,item,                                                     &
       icode,cmessage)

  IF (icode >  0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item ,''   //newline//&
   'qcf increment due to large scale precipitation'
  END IF
END IF ! Item 184

!--------------------------------------------------------------------------
! Item 189: qrain Increment
!--------------------------------------------------------------------------
item = 189
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), qrain_inc,              &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                          &
       at_extremity,                                                           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1,                                        &
       atmos_im,sect,item,                                                     &
       icode,cmessage)

  IF (icode >  0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'qrain increment due to large scale precipitation'
  END IF
END IF ! Item 189

!--------------------------------------------------------------------------
! Item 190: qgraup Increment
!--------------------------------------------------------------------------
item = 190
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), qgraup_inc,             &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                          &
       at_extremity,                                                           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1,                                        &
       atmos_im,sect,item,                                                     &
       icode,cmessage)

  IF (icode >  0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'qgraup increment due to large scale precipitation'
  END IF
END IF ! Item 190

!--------------------------------------------------------------------------
! Item 191: qcf2 Increment
!--------------------------------------------------------------------------
item = 191
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)), qcf2_inc,               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                          &
       at_extremity,                                                           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                     &
       stash_levels,num_stash_levels+1,                                        &
       atmos_im,sect,item,                                                     &
       icode,cmessage)

  IF (icode >  0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'qcf2 increment due to large scale precipitation'
  END IF
END IF ! Item 191

!--------------------------------------------------------------------------
! Item 201: Rainfall Amount at surface
!--------------------------------------------------------------------------
item = 201
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      work_2d(i,j) = casdiags % SurfaceRainR(i,j) * timestep
    END DO
  END DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), work_2d,           &
       tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,            &
       atmos_im, sect, item,                                          &
       icode, cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'Large Scale Rain Amount'
  END IF

  work_2d(:,:) = 0.0

END IF ! item 201

!--------------------------------------------------------------------------
! Item 202: Snowfall amount at surface
!--------------------------------------------------------------------------
item = 202
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      work_2d(i,j) = casdiags % SurfaceSnowR(i,j) * timestep
    END DO
  END DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), work_2d,           &
       tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,            &
       atmos_im, sect, item,                                          &
       icode, cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'Large Scale Snow Amount'
  END IF

  work_2d(:,:) = 0.0

END IF ! item 202

!--------------------------------------------------------------------------
! Item 203: Rainfall Rate at surface
!--------------------------------------------------------------------------
item = 203

! this is also used by stash 20014 precip_symbol
IF (icode <=0 .AND. (sf(stashcode_pws_precip_sym, stashcode_pws_sec))) THEN
  pws_precip_sym_ls_rain(:,:) = casdiags % SurfaceRainR(:,:)
END IF

IF (icode <= 0 .AND. (sf(item,sect) .OR. l_ac)) THEN

  IF (.NOT. ALLOCATED(lsrr)) THEN
    ALLOCATE ( lsrr(tdims%i_end*tdims%j_end) )
    lsrr(1) = rmdi
  END IF
END IF

IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), casdiags % SurfaceRainR, &
       tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,                  &
       atmos_im, sect, item,                                                &
       icode, cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'Large Scale Rain Rate'
  END IF

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ji = (j-1)*tdims%i_end+i
      lsrr(ji) = casdiags % SurfaceRainR(i,j)
    END DO
  END DO

END IF ! item 203


!--------------------------------------------------------------------------
! Item 204: Snowfall Rate at surface
!--------------------------------------------------------------------------
item = 204

! this is also used by stash 20014 precip_symbol
IF (icode <=0 .AND. (sf(stashcode_pws_precip_sym, stashcode_pws_sec))) THEN
  pws_precip_sym_ls_snow(:,:) = casdiags % SurfaceSnowR(:,:)
END IF

IF (icode <= 0 .AND. (sf(item, sect) .OR. l_ac)) THEN

  IF (.NOT. ALLOCATED(lssr)) THEN
    ALLOCATE ( lssr(tdims%i_end*tdims%j_end) )
    lssr(1) = rmdi
  END IF
END IF

IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), casdiags % SurfaceSnowR,  &
       tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,                   &
       atmos_im, sect, item,                                                 &
       icode, cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'Large Scale Snow Rate'
  END IF

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ji = (j-1)*tdims%i_end+i
      lssr(ji) = casdiags % SurfaceSnowR(i,j)
    END DO
  END DO

END IF ! item 204

!--------------------------------------------------------------------------
! Item 205: Cloud liquid water after precipitation
!--------------------------------------------------------------------------
item = 205
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = qcl_n(i,j,k) + qcl_inc(i,j,k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    work_3d,                                                          &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Cloud liquid water after precipitation'
  END IF
END IF ! item 205

!--------------------------------------------------------------------------
! Item 206: Cloud ice content after precipitation
!--------------------------------------------------------------------------
item = 206

IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = qcf2_n(i,j,k) + qcf2_inc(i,j,k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    work_3d,                                                          &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Cloud ice content after precipitation'
  END IF
END IF ! item 206

!--------------------------------------------------------------------------
! Item 207: Relative humidity with respect to ice (T<0degC) and water 
!           (T>0degC)
!--------------------------------------------------------------------------
item = 207

! Note: At the moment, this diagnostic uses the UM default qsat function
!       to calculate relative humidity. This uses lookup tables, whereas
!       within the CASIM code itself, we use the qsaturation function, 
!       which makes a direct calculation, although is more expensive.
!       We could examine making the change in the future, although the
!       difference in answers is expected to be minimal.

IF ( sf(item,sect) .AND. icode <= 0 ) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims, t_n, t_inc, t_work, q_n, q_inc, q_work)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        t_work(i,j,k) = t_n(i,j,k) + t_inc(i,j,k)
        q_work(i,j,k) = q_n(i,j,k) + q_inc(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  IF (l_new_qsat_mphys) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i,work_2d)                                                  &
!$OMP SHARED(tdims, t_work, p_layer_centres, work_3d, l_mr_physics, q_work)
    DO k = 1, tdims%k_end
      IF (l_mr_physics) THEN
        CALL qsat_mix_new( work_2d, t_work(:,:,k), p_layer_centres(:,:,k),    &
                           tdims%i_end, tdims%j_end )
      ELSE
        CALL qsat_new( work_2d, t_work(:,:,k), p_layer_centres(:,:,k),        &
                       tdims%i_end, tdims%j_end )
      END IF

      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
  
        ! Calculate relative humidity and pass as work 3d array
        !  Supersaturation (>100%) can occur with mixed phase scheme but
        !  negative humidity is removed from the diagnostic:
        work_3d(i,j,k) = MAX( 0.0, ((q_work(i,j,k) / work_2d(i,j)) * 100.0))

        END DO  ! i
      END DO    ! j
    END DO      ! k
!$OMP END PARALLEL DO

  ELSE ! not l_new_qsat_mphys
    ! DEPENDS ON: qsat_mix
    CALL qsat_mix(work_3d,t_work,p_layer_centres(1,1,1),                      &
                  tdims%i_end*tdims%j_end*tdims%k_end, l_mr_physics )

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,q_work)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          work_3d(i,j,k) = q_work(i,j,k)/work_3d(i,j,k)*100.0
          !  Supersaturation (>100%) can occur with mixed phase scheme but
          !  negative humidity is removed from the diagnostic:
          IF (work_3d(i,j,k)  <   0.0) THEN
            work_3d(i,j,k) = 0.0
          END IF
        END DO  ! i
      END DO    ! j
    END DO      ! k
!$OMP END PARALLEL DO

  END IF ! l_new_qsat_mphys

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,             &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'Relative humidity after microphysics.'
  END IF

END IF ! item 207

!--------------------------------------------------------------------------
! Item 208: Relative Humidity with respect to water after Microphysics [%]
!--------------------------------------------------------------------------
item = 208

! Note: At the moment, this diagnostic uses the UM default qsat function
!       to calculate relative humidity. This uses lookup tables, whereas
!       within the CASIM code itself, we use the qsaturation function, 
!       which makes a direct calculation, although is more expensive.
!       We could examine making the change in the future, although the
!       difference in answers is expected to be minimal.

IF ( sf(item,sect) .AND. icode <= 0 ) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims, t_n, t_inc, t_work, q_n, q_inc, q_work)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        t_work(i,j,k) = t_n(i,j,k) + t_inc(i,j,k)
        q_work(i,j,k) = q_n(i,j,k) + q_inc(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  IF ( l_new_qsat_mphys ) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i,work_2d)                                                  &
!$OMP SHARED(tdims, t_work, p_layer_centres, work_3d, q_work)
    DO k = 1, tdims%k_end
      CALL qsat_wat_new( work_2d, t_work(:,:,k), p_layer_centres(:,:,k),      &
                         tdims%i_end, tdims%j_end )

      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
  
        ! Calculate relative humidity and pass as work 3d array
        !  negative humidity is removed from the diagnostic:
        work_3d(i,j,k) = MAX( 0.0, ((q_work(i,j,k) / work_2d(i,j)) * 100.0))

        ! For this diagnostic, limit RH with respect to water to 100.0 %:
        work_3d(i,j,k) = MIN( work_3d(i,j,k), 100.0 )

        END DO  ! i
      END DO    ! j
    END DO      ! k
!$OMP END PARALLEL DO

  ELSE ! older version of qsat 

     ! DEPENDS ON: qsat_wat
     CALL qsat_wat(work_3d, t_work, p_layer_centres(1,1,1),                   &
                   tdims%i_end*tdims%j_end*tdims%k_end )

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,q_work)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = q_work(i,j,k)/work_3d(i,j,k)*100.0
             !  Supersaturation wrt water is limited to =< 100%
        work_3d(i,j,k) = MIN( work_3d(i,j,k), 100.0 )
             !  Negative humidity also removed from the diagnostic
        work_3d(i,j,k) = MAX( 0.0, work_3d(i,j,k) )
      END DO  ! i
    END DO    ! j
  END DO      ! k
!$OMP END PARALLEL DO

  END IF ! l_new_qsat_mphys

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,             &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'Relative humidity with respect to water after microphysics.'
  END IF

END IF ! item 208

!--------------------------------------------------------------------------
! Item 209: Large scale graupel amount
!--------------------------------------------------------------------------
item = 209

IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  ! Create graupel amount by multiplying graupel rate by the model timestep
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      work_2d(i,j) = casdiags % SurfaceGraupR(i,j) * timestep
    END DO
  END DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), work_2d,           &
       tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,            &
       atmos_im, sect, item,                                          &
       icode, cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'Large Scale Graupel Amount'
  END IF

  work_2d(:,:) = 0.0

END IF ! item 209

!--------------------------------------------------------------------------
! Item 212: Graupel Rate at surface
!--------------------------------------------------------------------------
item = 212
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), casdiags % SurfaceGraupR, &
       tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,                   &
       atmos_im, sect, item,                                                 &
       icode, cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag_3d for diagnostic section 4, item',item,''    //newline//&
   'Large Scale Graupel Rate'
  END IF

END IF ! item 212

!--------------------------------------------------------------------------
! Item 224: Supercooled liquid water content
!--------------------------------------------------------------------------
item = 224
IF ( sf(item,sect) .AND. icode <= 0 ) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        t_one_gb = t_n(i,j,k) + t_inc(i,j,k)

        IF ( t_one_gb < zerodegc ) THEN
          ! Supercooled temperatures
          work_3d(i,j,k) = qcl_n(i,j,k) + qcl_inc(i,j,k)
        ELSE
          ! Warm temperatures
          work_3d(i,j,k) = 0.0
        END IF ! temperature below zero
  
      END DO
    END DO
  END DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    work_3d,                                                          &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Supercooled liquid water content'
  END IF

END IF

!--------------------------------------------------------------------------
! Item 240: Homogeneous nucleation rate 
!--------------------------------------------------------------------------
item = 240
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % phomc,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Homogeneous nucleation rate'
  END IF
END IF ! item 240

!--------------------------------------------------------------------------
! Item 241: Heterogeneous nucleation rate
!--------------------------------------------------------------------------
item = 241
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pinuc,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Heterogeneous nucleation rate'
  END IF
END IF ! item 241

!--------------------------------------------------------------------------
! Item 243: Deposition rate ice crystals 
!--------------------------------------------------------------------------
item = 243
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pidep,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Deposition rate ice crystals'
  END IF
END IF ! item 243

!--------------------------------------------------------------------------
! Item 245: Deposition rate snow aggregates 
!--------------------------------------------------------------------------
item = 245
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % psdep,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Deposition rate snow aggregates'
  END IF
END IF ! item 245

!--------------------------------------------------------------------------
! Item 247: Riming rate ice crystals 
!--------------------------------------------------------------------------
item = 247
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % piacw,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Riming rate ice crystals'
  END IF
END IF ! item 247

!--------------------------------------------------------------------------
! Item 248: Riming rate snow aggregates 
!--------------------------------------------------------------------------
item = 248
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % psacw,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Riming rate snow aggregates'
  END IF
END IF ! item 248

!--------------------------------------------------------------------------
! Item 250: Snow aggregate-rain capture rate 
!--------------------------------------------------------------------------
item = 250
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % psacr,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Snow aggregate-rain capture rate'
  END IF
END IF ! item 250

!--------------------------------------------------------------------------
! Item 251: Evaporation (Sublimation) of melting ice crystals 
!--------------------------------------------------------------------------
item = 251
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pisub,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Evaporation (sublimation) of melting ice crystals'
  END IF
END IF ! item 251

!--------------------------------------------------------------------------
! Item 252: Evaporation (Sublimation) of melting snow aggregates 
!--------------------------------------------------------------------------
item = 252
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pssub,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Evaporation (sublimation) of melting snow aggregates'
  END IF
END IF ! item 252

!--------------------------------------------------------------------------
! Item 253: Melting rate of ice crystals 
!--------------------------------------------------------------------------
item = 253
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pimlt,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Melting rate of ice crystals'
  END IF
END IF ! item 253

!--------------------------------------------------------------------------
! Item 254: Melting rate of snow aggregates 
!--------------------------------------------------------------------------
item = 254
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % psmlt,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Melting rate of snow aggregates'
  END IF
END IF ! item 254

!--------------------------------------------------------------------------
! Item 255: Snow autoconversion rate 
!--------------------------------------------------------------------------
item = 255
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % psaut,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Snow autoconversion rate'
  END IF
END IF ! item 255

!--------------------------------------------------------------------------
! Item 256: Snow capture of ice crystals 
!--------------------------------------------------------------------------
item = 256
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % psaci,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Snow capture of ice crystals'
  END IF
END IF ! item 256

!--------------------------------------------------------------------------
! Item 257: Rain autoconversion rate 
!--------------------------------------------------------------------------
item = 257
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % praut,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Rain autoconversion rate'
  END IF
END IF ! item 257

!--------------------------------------------------------------------------
! Item 258: Rain accretion rate 
!--------------------------------------------------------------------------
item = 258
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pracw,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Rain accretion rate'
  END IF
END IF ! item 258

!--------------------------------------------------------------------------
! Item 259: Rain evaporation rate 
!--------------------------------------------------------------------------
item = 259
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % prevp,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Rain evaporation rate'
  END IF
END IF ! item 259

!--------------------------------------------------------------------------
! Item 261: Graupel accretion rate 
!--------------------------------------------------------------------------
item = 261
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pgacw,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Graupel accretion rate'
  END IF
END IF ! item 261

!--------------------------------------------------------------------------
! Item 262: Graupel-snow capture rate 
!--------------------------------------------------------------------------
item = 262
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pgacs,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Graupel-snow capture rate'
  END IF
END IF ! item 262

!--------------------------------------------------------------------------
! Item 263: Graupel melting rate 
!--------------------------------------------------------------------------
item = 263
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pgmlt,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Graupel melting rate'
  END IF
END IF ! item 263

!--------------------------------------------------------------------------
! Item 264: Graupel evaporation (sublimation) rate 
!--------------------------------------------------------------------------
item = 264
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pgsub,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Graupel evaporation (sublimation) rate'
  END IF
END IF ! item 264

!--------------------------------------------------------------------------
! Item 265: Ice crystal sedimentation rate 
!--------------------------------------------------------------------------
item = 265
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % psedi,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Ice crystal sedimentation rate'
  END IF
END IF ! item 265

!--------------------------------------------------------------------------
! Item 266: Snow aggregate sedimentation rate 
!--------------------------------------------------------------------------
item = 266
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pseds,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Snow aggregate sedimentation rate'
  END IF
END IF ! item 266

!--------------------------------------------------------------------------
! Item 267: Rain sedimentation rate 
!--------------------------------------------------------------------------
item = 267
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % psedr,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Rain sedimentation rate'
  END IF
END IF ! item 267

!--------------------------------------------------------------------------
! Item 268: Graupel sedimentation rate 
!--------------------------------------------------------------------------
item = 268
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % psedg,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Graupel sedimentation rate'
  END IF
END IF ! item 268

!--------------------------------------------------------------------------
! Item 269: Liquid cloud sedimentation rate
!--------------------------------------------------------------------------
item = 269
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % psedl,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Liquid cloud sedimentation rate'
  END IF
END IF ! item 269

!--------------------------------------------------------------------------
! Item 271: Homogeneous freezing of rain rate
!--------------------------------------------------------------------------
item = 271
IF ( sf(item,sect) .AND. icode <= 0 ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % phomr,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Homogeneous freezing of rain rate'
  END IF
END IF ! item 271

!--------------------------------------------------------------------------
! Item 302: Snowfall amount excluding graupel
!--------------------------------------------------------------------------
item = 302

IF (icode <= 0 .AND. sf(item,sect)) THEN

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      work_2d(i,j) = (casdiags % SurfaceSnowR (i,j) -                          &
                      casdiags % SurfaceGraupR(i,j)   ) * timestep
    END DO
  END DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), work_2d,                   &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Snowfall amount excluding graupel'
  END IF

  work_2d(:,:) = 0.0

END IF ! item 302

!--------------------------------------------------------------------------
! Item 304: Snowfall rate excluding graupel
!--------------------------------------------------------------------------
item = 304

IF (icode <= 0 .AND. sf(item,sect)) THEN

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      work_2d(i,j) = (casdiags % SurfaceSnowR (i,j) -                          &
                      casdiags % SurfaceGraupR(i,j)   )
    END DO
  END DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), work_2d,                   &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Snowfall rate excluding graupel'
  END IF

  work_2d(:,:) = 0.0

END IF ! item 304

!--------------------------------------------------------------------------
! Item 325: Cloud Condensation/Evaporation Rate
!--------------------------------------------------------------------------
item = 325

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % pcond,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Cloud condensation/evaporation rate'
  END IF

END IF ! item 325

!--------------------------------------------------------------------------
! Item 336: Ice Cloud sedimentation rate [kg/kg/s]
!--------------------------------------------------------------------------
item = 336

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % psedi,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Ice cloud sedimentation rate.'
  END IF

END IF ! item 336

!--------------------------------------------------------------------------
! Item 350: Homogeneous Freezing of Cloud Number Tendency [Num. kg-1 s-1]
!--------------------------------------------------------------------------
item = 350

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % nhomc,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Homogeneous freezing of cloud number tendency'
  END IF

END IF ! item 350

!--------------------------------------------------------------------------
! Item 351: Homogeneous Freezing of Rain Number Tendency [Num. kg-1 s-1]
!--------------------------------------------------------------------------
item = 351

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % nhomr,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Homogeneous freezing of rain number tendency'
  END IF

END IF ! item 351

!----------------------------------------------------------------------------
! Item 352: Ice number tendency due to Hallett-Mossop Process [Num. kg-1 s-1]
!----------------------------------------------------------------------------
item = 352

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % nihal,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Ice number tendency due to Hallett-Mossop Process.'
  END IF

END IF ! item 352

!--------------------------------------------------------------------------
! Item 353: Ice number tendency due to Ice Nucleation [Num. kg-1 s-1]
!--------------------------------------------------------------------------
item = 353

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % ninuc,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Ice number tendency due to ice nucleation.'
  END IF

END IF ! item 353

!--------------------------------------------------------------------------
! Item 354: Tendency in ice number due to sedimentation [Num. kg-1 s-1]
!--------------------------------------------------------------------------
item = 354

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % nsedi,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Tendency in ice number due to sedimentation'
  END IF

END IF ! item 354

!--------------------------------------------------------------------------
! Item 355: Tendency in snow number due to sedimentation [Num. kg-1 s-1]
!--------------------------------------------------------------------------
item = 355

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % nseds,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Tendency in snow number due to sedimentation'
  END IF

END IF ! item 355

!--------------------------------------------------------------------------
! Item 356: Tendency in graupel number due to sedimentation [Num. kg-1 s-1]
!--------------------------------------------------------------------------
item = 356

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                  &
    stashwork(si(item,sect,im_index)),                                &
    casdiags % nsedg,                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,      &
    stlist(1,stindex(1,item,sect,im_index)),                          &
    len_stlist, stash_levels,num_stash_levels+1,                      &
    atmos_im,sect,item,                                               &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage, FMT='(A,I3,A)')                                            &
   'Error in copydiag for diagnostic section 4, item',item,''    //newline//   &
   'Tendency in graupel number due to sedimentation'
  END IF

END IF ! item 356

!--------------------------------------------------------------------------
! Single point exception handling
IF (icode > 0) CALL ereport(RoutineName,icode,cmessage)
!==============================================================================
! Call Dr Hook and end the subroutine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_casim

END MODULE diagnostics_casim_mod


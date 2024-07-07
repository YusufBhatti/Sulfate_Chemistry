! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE LHN_INC ------------------------------------------------
!
!    Purpose : Latent Heat nudging of model heating profiles
!
!
!    Programming Standard : UM Doc Paper No 3
!
!    Project Task : P3
!    Documentation :  FR  WP 171
!
!
!    ARGUMENTS
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE lhn_inc_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LHN_INC_MOD'

CONTAINS

SUBROUTINE lhn_inc(conv_heat,ls_heat,ls_rain,ls_snow,conv_rain,   &
                   conv_snow,p_field,pr_inc,timestep,    &
                   irows,icols,                                   &
                   RANGE,phi_p,lambda_p,                          &
                   thincs,lenmg,icode,cmessage)

USE planet_constants_mod, ONLY: planet_radius
USE conversions_mod, ONLY: pi_over_180, rsec_per_day
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype, nproc
USE UM_ParParams, ONLY: halo_type_no_halo
USE field_types, ONLY: fld_type_p
USE atmos_max_sizes, ONLY: Max2DFieldSize
USE comobs_mod, ONLY: nobtypmx
USE lhn_search_mod, ONLY: lhn_search
USE rfcsl_mod, ONLY: rfcsl
USE rfcslr_mod, ONLY: rfcslr
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength
USE model_domain_mod, ONLY: l_regular
USE nlsizes_namelist_mod, ONLY: model_levels

IMPLICIT NONE

!-----DECLARE VARIABLES
INTEGER ::    p_field,lenmg
INTEGER ::    irows,icols,RANGE
REAL ::       pr_inc(lenmg),timestep,                             &
              ls_snow(p_field),ls_rain(p_field),                  &
              conv_snow(p_field),conv_rain(p_field),              &
              conv_heat(p_field,model_levels),                        &
              ls_heat(p_field,model_levels),                          &
              thincs(p_field,model_levels)

REAL :: lambda_p(1-halo_i:icols+halo_i)
REAL :: phi_p(1-halo_i:icols+halo_i, 1-halo_j:irows+halo_j)

!
INTEGER ::    icode
CHARACTER(LEN=errormessagelength) :: cmessage

!-INTENT=IN---------------------------------------------------
!     IROWS                      - No. of rows in model grid
!     ICOLS                      - No. of pts in a row
!     RANGE                      - Search range in grid points
!     LENMG                      - Length of model grid
!     PR_INC(LENMG)              - Rainfall incrs on model grid
!     TIMESTEP                   - Timestep in seconds
!     LS_SNOW(P_FIELD)          \                                      .
!     LS_RAIN(P_FIELD)           \  Large scale and convective
!     CONV_SNOW(P_FIELD)         /    rain and snow rates
!     CONV_RAIN(P_FIELD)        /           (diagnostic)
!     CONV_HEAT(P_FIELD,LEVELS)- L.H. incrs to theta due to conv'n
!     LS_HEAT(P_FIELD,LEVELS)  - L.H. incrs to theta due to dynamics
!-INTENT=INOUT-----------------------------------------------
!
!-INTENT=OUT-------------------------------------------------
!     THINCS                     - Calculated increments to theta
!     ICODE                      - Non-zero for failure
!     CMESSAGE                   - Reason for failure
! -----------------------------------------------------------

! Local arrays and variables
INTEGER ::    search(4*RANGE*(RANGE+1),2)
                                       ! Search Template
INTEGER ::    jpts , jlevs             ! Loop counters over model
                                       !       points and levels
INTEGER ::    near(p_field)            ! Nearest point found
                                       !    by LHN_SEARCH
INTEGER ::    no_incrs                 ! Diagnostic - no. of incrs
INTEGER ::    no_search                !   ""    - no. of searches
INTEGER ::    radius(5)                !   "" - breakdown of
                                       !        search results
INTEGER ::    m_grid                   ! Grid type for filter

REAL ::       tot_pr(p_field)          ! Total model precip rate
REAL ::       anal_pr(p_field)         ! Obs ppn on model grid
REAL ::       tot_lh(p_field,model_levels) ! Total LH profile
REAL ::       limit                    ! Timestep limit of incrs
REAL ::       z                        ! used to set COS_LAT
REAL ::       cos_lat(rows_max)        ! used in filter routine
REAL ::       fi_lhn                   ! filter scale in radians
LOGICAL ::    l_first_search           ! True if first search of
                                       !       the timestep

REAL :: pr_inc_global(Max2DFieldSize)
REAL :: tot_pr_global(Max2DFieldSize)
REAL :: anal_pr_global(Max2DFieldSize)
INTEGER :: near_global(Max2DFieldSize)
REAL :: tot_lh_global(glsize(1,fld_type_p)*glsize(2,fld_type_p),  &
                      model_levels/nproc+1)
REAL :: thincs_global(glsize(1,fld_type_p)*glsize(2,fld_type_p),  &
                      model_levels/nproc+1)

INTEGER :: row_length_g, p_rows_g, p_field_g, k, istat, iproc
INTEGER ::                                                        &
  map(model_levels)                                                   &
                      ! processor number for level
, n_levs_on_proc(0:nproc-1)                                       &
                            ! number of levels on each processor
, lev_on_gatpe(model_levels)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LHN_INC'


! Variables required for variable resolution grids
! GRID_LATS_G are the latitudes for the entire domain.
! Since the total number of rows is unknown here it is
! made allocatable (although they could be dimensioned using glsize
! as tot_lh_global and thincs_global
REAL, ALLOCATABLE :: grid_lats_g(:)
REAL :: phi_p_local(icols,irows)
REAL :: lambda_p_local(icols,irows)
REAL, ALLOCATABLE :: phi_p_global(:,:)
REAL, ALLOCATABLE :: lambda_p_global(:,:)
REAL, ALLOCATABLE :: delta_lat(:)
REAL, ALLOCATABLE :: delta_lon(:)
INTEGER :: i,j

!
!  set up global dimensions
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
row_length_g=glsize(1,fld_type_p)
p_rows_g=glsize(2,fld_type_p)
p_field_g= row_length_g*p_rows_g


ALLOCATE (grid_lats_g(p_rows_g))

IF (ltimer_ac) CALL timer('LHN_INC ',3)
IF (ldiagac .AND. mype == 0) THEN
  WRITE(umMessage,*) "***** Starting LHN_INC"
  CALL umPrint(umMessage,src='lhn_inc')
END IF
!
!C 1.0   set parameters and variables
!
IF (l_lhn_search .AND. RANGE == 0) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,*) "WARNING : LHN_RANGE=0,"
    CALL umPrint(umMessage,src='lhn_inc')
    WRITE(umMessage,*) "   therefore, setting L_LHN_SEARCH to FALSE"
    CALL umPrint(umMessage,src='lhn_inc')
  END IF
  l_lhn_search = .FALSE.
END IF

! Grid latitudes. If regular grid calculate otherwise obtain from phi_p
IF (l_regular) THEN

  DO jpts=1,p_rows_g
    grid_lats_g(jpts)=row1mgth+REAL(jpts-1)*dlatmg
  END DO

ELSE

  ALLOCATE (phi_p_global(row_length_g,p_rows_g))
  ALLOCATE (lambda_p_global(row_length_g,p_rows_g))

  ! Copy information from PHI_P/LAMBDA_P to PHI_P_LOCAL/LAMBDA_P_LOCAL
  !  excluding halos
  DO i=1,irows
    DO j=1,icols
      phi_p_local(j,i)=phi_p(j,i)
      lambda_p_local(j,i)=lambda_p(j)
    END DO
  END DO

  ! Gather up all PHI_P values into global array and distribute to all PEs
  ! DEPENDS ON: all_gather_field
  CALL all_gather_field(phi_p_local,   phi_p_global,              &
                    icols,         irows,                         &
                    row_length_g,  p_rows_g,                      &
                    fld_type_p,    halo_type_no_halo,             &
                    gc_all_proc_group,                            &
                    icode,         cmessage)

  ! DEPENDS ON: all_gather_field
  CALL all_gather_field(lambda_p_local,   lambda_p_global,        &
                    icols,         irows,                         &
                    row_length_g,  p_rows_g,                      &
                    fld_type_p,    halo_type_no_halo,             &
                    gc_all_proc_group,                            &
                    icode,         cmessage)

  ! Extract one column of phi values into GRID_LATS_G
  ! inverting so colats run in correct direction
  DO jpts=1,p_rows_g
    grid_lats_g(jpts)=phi_p_global(1,p_rows_g-jpts+1)/          &
                                    pi_over_180
    grid_lats_g(jpts)=90.0-grid_lats_g(jpts)
    grid_lats_g(jpts)=grid_lats_g(jpts)*pi_over_180
  END DO

END IF

! Populate DELTA_LAT and DELTA_LON arrays

ALLOCATE (delta_lat(p_rows_g))
ALLOCATE (delta_lon(row_length_g))

IF (l_regular) THEN

  ! These are not used at present but will be needed if RFCSL is used
  ! instead of RFCSLR but bit comparibilty will be lost.
  !        DELTA_LAT(:)=DLATMG
  !        DELTA_LON(:)=DLONGMG

ELSE

  DO jpts=1,p_rows_g-1
    delta_lat(jpts)=phi_p_global(1,jpts+1)-phi_p_global(1,jpts)
  END DO
  delta_lat(p_rows_g)=delta_lat(p_rows_g-1)
  DO jpts=1,row_length_g-1
    delta_lon(jpts)=lambda_p_global(jpts+1,1)-                    &
                   lambda_p_global(jpts,1)
  END DO
  delta_lon(row_length_g)=delta_lon(row_length_g-1)

END IF

IF (l_lhn_limit) limit = lhn_limit * timestep / rsec_per_day

l_first_search = .TRUE.
fi_lhn    = fi_scale_lhn / planet_radius
no_incrs  = 0
no_search = 0
radius(1) = 0
radius(2) = 0
radius(3) = 0
radius(4) = 0
radius(5) = 0

!
!C 2     initial calculations
!
!C 2.1   precipitation variables

!
DO jpts = 1 , p_field
  tot_pr(jpts) = ls_rain(jpts) + ls_snow(jpts)+                   &
                          conv_rain(jpts) + conv_snow(jpts)
  IF (tot_pr(jpts)  <   0.0) tot_pr(jpts)=0.0
  anal_pr(jpts) = tot_pr(jpts) + pr_inc(jpts)
  IF (anal_pr(jpts)  <   0.0) THEN
    anal_pr(jpts) = 0.0
    pr_inc(jpts)   = anal_pr(jpts) - tot_pr(jpts)
  END IF
END DO     ! JPTS

!
!C 2.2   heating profiles, (FR Working Paper 171, eq2)
!
DO jlevs = 1 , model_levels
  DO jpts = 1 , p_field
    tot_lh(jpts,jlevs) = ( conv_heat(jpts,jlevs) +                &
                     ls_heat(jpts,jlevs) ) * timestep
    IF ( remove_neg_lh .AND. (tot_lh(jpts,jlevs)                      &
           <   0.0) ) THEN
      tot_lh(jpts,jlevs) = 0.0
    END IF
  END DO   ! JPTS
END DO     ! JLEVS

!
!C 3      Calculate theta increments
!
!C 3.1    Set up array NEAR to decide where to get profile for scaling
!         NEAR=-1 implies no scaling necessary, 0 implies scale here
!         greater than 0 means scale a nearby profile, or not at all
!
!  PR_INC, TOT_PR and ANAL_PR are needed on all PEs,
!  so gather them to all pes.

! DEPENDS ON: all_gather_field
CALL all_gather_field(pr_inc,       PR_INC_global,                &
                  icols,        irows,                            &
                  row_length_g, p_rows_g,                         &
                  fld_type_p,   halo_type_no_halo,                &
                  gc_all_proc_group,                              &
                  icode,        cmessage)


! DEPENDS ON: all_gather_field
CALL all_gather_field(tot_pr,        TOT_PR_global,               &
                  icols,         irows,                           &
                  row_length_g,  p_rows_g,                        &
                  fld_type_p,    halo_type_no_halo,               &
                  gc_all_proc_group,                              &
                  icode,         cmessage)


! DEPENDS ON: all_gather_field
CALL all_gather_field(anal_pr,       ANAL_PR_global,              &
                  icols,         irows,                           &
                  row_length_g,  p_rows_g,                        &
                  fld_type_p,    halo_type_no_halo,               &
                  gc_all_proc_group,                              &
                  icode,         cmessage)



IF (ltimer_ac) CALL timer('LHN_SRCH',3)

!  Calculate NEAR_global on all PEs
DO jpts = 1 , p_field_g
  NEAR_global(jpts) = -1
  IF ( PR_INC_global(jpts)  /=  0.0 ) THEN
    NEAR_global(jpts) = 0
    no_incrs   = no_incrs + 1
    IF ( TOT_PR_global(jpts)  <                                   &
          (epsilon_lhn * ANAL_PR_global(jpts)) ) THEN
      NEAR_global(jpts) = jpts
      no_search = no_search + 1
      IF (l_lhn_search) THEN
        !  search for suitable profile
        CALL lhn_search(jpts, NEAR_global(jpts),                  &
           RANGE, search, p_rows_g, row_length_g,                 &
           l_first_search, TOT_PR_global,                         &
           ANAL_PR_global(jpts), p_field_g, radius,jpts)
      END IF
    END IF
  END IF
END DO    ! JPTS

IF (ltimer_ac) CALL timer('LHN_SRCH',4)
!
!C 3.2   calculate the increments, and scale by the relaxation coeff
!

! Calculate the mapping of which processor each level will go to
DO k=0,nproc-1
  n_levs_on_proc(k)=0
END DO

DO k=1, model_levels
  ! assumes first PE is PE 0
  map(k)                 = MOD((k-1),nproc)
  n_levs_on_proc(map(k)) = n_levs_on_proc(map(k))+1
  lev_on_gatpe(k)        = n_levs_on_proc(map(k))
END DO

! Distribute TOT_LH_GLOBAL over the processors
DO k=1, model_levels
  ! DEPENDS ON: gather_field
  CALL gather_field( tot_lh(:,k),                                 &
                     tot_lh_global(:,lev_on_gatpe(k)),            &
                     icols, irows,                                &
                     row_length_g, p_rows_g,                      &
                     fld_type_p, halo_type_no_halo,               &
                     map(k), gc_all_proc_group )
END DO


DO jlevs = 1 , n_levs_on_proc(mype)

  DO jpts = 1 , p_field_g
    IF (NEAR_global(jpts)  <   0) THEN
      !  No scaling necessary
      THINCS_global(jpts,jlevs) = 0.0
    ELSE IF (NEAR_global(jpts)  ==  0) THEN
      !  Scale the profile at the point itself
      !  (FR WP 171, eq 5)
      THINCS_global(jpts,jlevs) =                                 &
               relax_cf_lhn * TOT_LH_global(jpts,jlevs) *         &
               PR_INC_global(jpts) / TOT_PR_global(jpts)
    ELSE IF (NEAR_global(jpts)  >   0 .AND.                        &
                               NEAR_global(jpts)  /=  jpts) THEN
      !  Scale a nearby profile
      !  (FR WP 171, eq 7)
      THINCS_global(jpts,jlevs) = relax_cf_lhn *                  &
         (ANAL_PR_global(jpts)/TOT_PR_global(NEAR_global(jpts)) * &
         TOT_LH_global(NEAR_global(jpts),jlevs) -                 &
         TOT_LH_global(jpts,jlevs))
    ELSE IF (NEAR_global(jpts)  ==  jpts .AND. l_lhn_scale) THEN
      !  Scale by 1/EPSILON_LHN if no suitable profile available
      !  (FR WP 171, eq 8)
      THINCS_global(jpts,jlevs) = relax_cf_lhn *                  &
         (( 1.0/epsilon_lhn) - 1.0) * TOT_LH_global(jpts,jlevs)
    ELSE
      !  If none of the above apply, ignore the ob
      THINCS_global(jpts,jlevs) = 0.0
    END IF
  END DO ! JPTS

END DO ! JLEVS

IF (.NOT. l_lhn_filt) THEN
  ! Copy calculated THINCS values back to individual PEs
  DO k=1, model_levels
    ! DEPENDS ON: scatter_field
    CALL scatter_field(thincs(:,k),                               &
                       thincs_global(:,lev_on_gatpe(k)),          &
                       icols, irows,                              &
                       row_length_g, p_rows_g,                    &
                       fld_type_p, halo_type_no_halo,             &
                       map(k), gc_all_proc_group)
  END DO

  !
  !  Impose limit by factor of alpha
  !
  IF (l_lhn_fact) THEN
    DO jlevs = 1 , model_levels
      DO jpts = 1 , p_field
        IF (pr_inc(jpts)  <   0.0) THEN
          IF ((pr_inc(jpts)/tot_pr(jpts)) <  (alpha_lhn-1.0)) THEN
            thincs(jpts,jlevs) = (alpha_lhn - 1.0) *              &
                                 tot_lh(jpts,jlevs) * relax_cf_lhn
          END IF
        END IF
      END DO   ! JPTS
    END DO     ! JLEVS
  END IF
  !
  !C 3.3   Recursive filtering of increments
  !
ELSE    ! L_LHN_FILT is true

  IF (l_lhn_fact) THEN
    DO jlevs = 1 , n_levs_on_proc(mype)
      DO jpts = 1 , P_FIELD_g
        IF (PR_INC_global(jpts)  <   0.0) THEN
          IF ((PR_INC_global(jpts)/TOT_PR_global(jpts)) <         &
                                      (alpha_lhn-1.0)) THEN

            THINCS_global(jpts,jlevs) = (alpha_lhn - 1.0) *       &
                         TOT_LH_global(jpts,jlevs) * relax_cf_lhn
          END IF
        END IF
      END DO   ! JPTS
    END DO     ! JLEVS
  END IF

  !       Initialise COS_LAT
  IF (l_regular) THEN
    z = row1mgth
    DO jpts = 1, p_rows_g
      cos_lat(jpts) = SIN(z)
      z = z + dlatmg
    END DO  ! JPTS
  ELSE
    DO jpts = 1, p_rows_g
      cos_lat(jpts)=SIN(grid_lats_g(jpts))
    END DO  ! JPTS
  END IF


  m_grid = 0

  IF (l_regular) THEN
    DO jlevs = 1 , n_levs_on_proc(mype)
      CALL rfcslr(THINCS_global(1,jlevs),p_rows_g,row_length_g,m_grid,&
                 0.0,cos_lat,dlatmg,dlongmg,fi_lhn,npass_rf_lhn)
    END DO

  ELSE

    DO jlevs = 1 , n_levs_on_proc(mype)
      CALL rfcsl(THINCS_global(1,jlevs),p_rows_g,row_length_g,m_grid, &
                 0.0,cos_lat,delta_lat,delta_lon,fi_lhn,npass_rf_lhn)
    END DO

  END IF

  ! Copy calculated THINCS values back to individual PEs
  DO k=1, model_levels
    ! DEPENDS ON: scatter_field
    CALL scatter_field(thincs(:,k),                                 &
                       thincs_global(:,lev_on_gatpe(k)),            &
                       icols, irows,                                &
                       row_length_g, p_rows_g,                      &
                       fld_type_p, halo_type_no_halo,               &
                       map(k), gc_all_proc_group)
  END DO

END IF   ! L_LHN_FILT


!
!C Impose absolute limit
!
IF (l_lhn_limit) THEN
  DO jlevs = 1 , model_levels
    DO jpts = 1 , p_field
      IF (thincs(jpts,jlevs)  <   -1.0 * limit)               &
             thincs(jpts,jlevs) = -1.0 * limit
      IF (thincs(jpts,jlevs)  >   limit)                      &
             thincs(jpts,jlevs) = limit
    END DO  ! JPTS
  END DO    ! JLEVS
END IF

!
!C  4  Diagnostics
!
IF (ldiagac .AND. lhn_diag) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,*) ' '
    CALL umPrint(umMessage,src='lhn_inc')
    WRITE(umMessage,*) "Latent Heat Nudging Scheme, LHN_INC, diagnostics"
    CALL umPrint(umMessage,src='lhn_inc')
    WRITE(umMessage,*) "    Parameters set   : EPSILON_LHN = ",epsilon_lhn
    CALL umPrint(umMessage,src='lhn_inc')
    WRITE(umMessage,*) "                     : LHN_RANGE   = ",RANGE
    CALL umPrint(umMessage,src='lhn_inc')
    IF (l_lhn_fact) THEN
      WRITE(umMessage,*) "    Limit increments to scale down rain rate",   &
                  " by at most factor of 1/ ALPHA"
      CALL umPrint(umMessage,src='lhn_inc')
      WRITE(umMessage,*) "       ALPHA=",alpha_lhn
      CALL umPrint(umMessage,src='lhn_inc')
    END IF
    IF (l_lhn_filt) THEN
      WRITE(umMessage,*) "    Filtering of increments performed"
      CALL umPrint(umMessage,src='lhn_inc')
      WRITE(umMessage,*) "       filter scale = ",fi_scale_lhn             &
                                                    / 1000.0," Km"
      CALL umPrint(umMessage,src='lhn_inc')
    END IF
    IF (l_lhn_limit) THEN
      WRITE(umMessage,*) "    Limiting of increments set to ",limit,       &
                  " degrees per timestep = ",lhn_limit," degrees",  &
                                                       " per day"
      CALL umPrint(umMessage,src='lhn_inc')
    END IF
    WRITE(umMessage,*) ' '
    CALL umPrint(umMessage,src='lhn_inc')
    WRITE(umMessage,*) "Number of increments required was ",no_incrs
    CALL umPrint(umMessage,src='lhn_inc')
    IF (l_lhn_search) THEN
      radius(1) = radius(2) + radius(3) + radius(4)
      WRITE(umMessage,*) "LHN_SEARCH was called for ",no_search," of them"
      CALL umPrint(umMessage,src='lhn_inc')
      WRITE(umMessage,*) "The search was successful ",radius(1)," times"
      CALL umPrint(umMessage,src='lhn_inc')
      WRITE(umMessage,*) "It found :"
      CALL umPrint(umMessage,src='lhn_inc')
      WRITE(umMessage,*) "     ",radius(2)," pts at 1 point away"
      CALL umPrint(umMessage,src='lhn_inc')
      WRITE(umMessage,*) "     ",radius(3)," pts at 2 points away"
      CALL umPrint(umMessage,src='lhn_inc')
      WRITE(umMessage,*) "     ",radius(4)," pts at 3 or more points away"
      CALL umPrint(umMessage,src='lhn_inc')
      WRITE(umMessage,*) "     ",radius(5)," total pts searched"
      CALL umPrint(umMessage,src='lhn_inc')
    ELSE
      WRITE(umMessage,*) no_search," points required the search, but"
      CALL umPrint(umMessage,src='lhn_inc')
      WRITE(umMessage,*) " the search algorithm was disabled"
      CALL umPrint(umMessage,src='lhn_inc')
    END IF
    WRITE(umMessage,*) "Scaling at points failing the search was set to: ",&
                          l_lhn_scale
    CALL umPrint(umMessage,src='lhn_inc')
    WRITE(umMessage,*) ' '
    CALL umPrint(umMessage,src='lhn_inc')
  END IF
END IF

IF (ltimer_ac) CALL timer('LHN_INC ',4)

IF (ALLOCATED(grid_lats_g)) DEALLOCATE(grid_lats_g)
IF (ALLOCATED(phi_p_global)) DEALLOCATE (phi_p_global)
IF (ALLOCATED(lambda_p_global)) DEALLOCATE (lambda_p_global)
IF (ALLOCATED(delta_lat)) DEALLOCATE (delta_lat)
IF (ALLOCATED(delta_lon)) DEALLOCATE (delta_lon)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lhn_inc
END MODULE lhn_inc_mod

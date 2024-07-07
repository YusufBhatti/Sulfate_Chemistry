! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE LHN_INC_1A --------------------------------------------
!
!    Purpose : Latent Heat nudging of model heating profiles
!
!    Version 1A:
!    This version works on local fields in the logical processed grid
!    (LPG) domain.
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
MODULE lhn_inc_1a_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LHN_INC_1A_MOD'

CONTAINS

SUBROUTINE lhn_inc_1a(conv_heat,ls_heat,ls_rain,ls_snow,conv_rain,   &
                      conv_snow,p_field,pr_inc,timestep,             &
                      rows,row_length,                               &
                      srange,phi_p,lambda_p,                         &
                      thincs,lenmg,icode,cmessage)

USE planet_constants_mod, ONLY: planet_radius
USE conversions_mod, ONLY: pi_over_180,            &
                           rsec_per_day
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore,  ONLY: mype
USE Field_Types, ONLY: fld_type_p
USE comobs_mod, ONLY: nobtypmx
USE lhn_search_1A_mod, ONLY: lhn_search_1A
USE rfcsl_1a_mod, ONLY: rfcsl_1a
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength
USE model_domain_mod, ONLY: l_regular
USE nlsizes_namelist_mod, ONLY: model_levels
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE mpl, ONLY: mpl_real
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)    :: p_field            ! Number of points in mass field
INTEGER, INTENT(IN)    :: rows               ! Rows in local grid
INTEGER, INTENT(IN)    :: row_length         ! Points in each row of local grid
INTEGER, INTENT(IN)    :: srange              ! Search range in grid points
INTEGER, INTENT(IN)    :: lenmg              ! Length of model grid
REAL,    INTENT(IN)    :: timestep           ! Timestep in seconds
REAL,    INTENT(IN)    :: ls_snow(p_field)   ! \   Large scale
REAL,    INTENT(IN)    :: ls_rain(p_field)   !  \  and convective
REAL,    INTENT(IN)    :: conv_snow(p_field) !  /  rain and snow rates 
REAL,    INTENT(IN)    :: conv_rain(p_field) ! /   (diagnostic)
REAL,    INTENT(IN)    :: ls_heat(p_field,model_levels)
                                             ! LH incrs to theta due to dynamics
REAL,    INTENT(IN)    :: conv_heat(p_field,model_levels)
                                             ! L.H. incrs to theta due to conv'n
                                             ! 
REAL,    INTENT(IN)    :: lambda_p(1-halo_i:row_length+halo_i) 
                                             ! Used in variable resolution runs
REAL,    INTENT(IN)    :: phi_p(1-halo_i:row_length+halo_i, & 
                                1-halo_j:rows+halo_j)
                                             ! Used in variable resolution runs
                                             !
REAL,    INTENT(INOUT) :: pr_inc(lenmg)      ! Rainfall incrs on model grid 
REAL,    INTENT(OUT)   :: thincs(p_field,model_levels)             
                                             ! Calculated increments to theta
INTEGER, INTENT(INOUT) :: icode              ! Non-zero for failure
CHARACTER(LEN=errormessagelength), INTENT(INOUT) :: cmessage
                                             ! Reason for failure

! Local arrays and variables
INTEGER, ALLOCATABLE, SAVE :: search(:,:)    ! Search Template
INTEGER :: jr,jn                             ! Search loop counters
INTEGER :: k,j,i                             ! Loop counters
INTEGER :: jpts_l                            ! Local index of a point
INTEGER :: jpts_g                            ! Global index of a point
INTEGER :: i_g                               ! Column global coord of a point
INTEGER :: j_g                               ! Row     "    "   "
INTEGER :: i_l                               ! Column local coord of a point
INTEGER :: j_l                               ! Row     "    "   "
INTEGER :: rows_g                            ! Rows in global grid
INTEGER :: row_length_g                      ! Points in each row of global grid
INTEGER :: near(p_field)                     ! Nearest point found
                                             !    by LHN_SEARCH
INTEGER :: no_incrs                          ! Diagnostic - no. of incrs
INTEGER :: no_search                         !   ""    - no. of searches
INTEGER :: radius(5)                         !   "" - breakdown of
                                             !        search results
INTEGER :: m_grid                            ! Grid type for filter

REAL    :: anal_pr(p_field)                  ! Obs ppn on model grid
REAL    :: tot_pr(p_field)                   ! Total model precip rate
REAL    :: tot_pr_r(1-srange:row_length+srange,1-srange:rows+srange)          
                                             ! Total model precip rate
                                             ! with srange as halos.
REAL    :: tot_lh(p_field,model_levels)      ! Total LH profile
REAL    :: tot_lh_r(1-srange:row_length+srange,1-srange:rows+srange, & 
                    model_levels)       
                                             ! Total LH profile with
                                             ! srange as halos
REAL    :: limit                             ! Timestep limit of incrs
REAL    :: z                                 ! used to set COS_LAT
REAL    :: cos_lat(rows_max)                 ! used in filter routine
REAL    :: fi_lhn                            ! filter scale in radians

INTEGER :: p_field_g                         ! Number of points in
                                             ! global field


! Variables required for variable resolution grids
! GRID_LATS_G are the latitudes for the entire domain.
REAL, ALLOCATABLE, SAVE :: grid_lats_g(:)
REAL, ALLOCATABLE, SAVE :: delta_lat_g(:)
REAL, ALLOCATABLE, SAVE :: delta_lon_g(:)
REAL, ALLOCATABLE       :: phi_p_g_col(:)
REAL, ALLOCATABLE       :: phi_p_l_col(:)
REAL, ALLOCATABLE       :: lambda_p_g_row(:)
REAL, ALLOCATABLE       :: lambda_p_l_row(:)

! Variables used for mpl_allghaterv between column PEs and row PEs.
INTEGER :: counts_y(0:nproc_y-1)             
INTEGER :: displacements_y(0:nproc_y-1)        
INTEGER :: counts_x(0:nproc_x-1)             
INTEGER :: displacements_x(0:nproc_x-1)        
INTEGER :: pe, pe_j, pe_i                    ! PEs loop counters
INTEGER :: ierr
INTEGER :: dsize

LOGICAL, SAVE :: first_call = .TRUE.       

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LHN_INC_1A'

!
!  Set up global dimensions
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

row_length_g = glsize(1,fld_type_p)
rows_g       = glsize(2,fld_type_p)
p_field_g    = row_length_g*rows_g

IF (ltimer_ac) CALL timer('LHN_INC ',3)
IF (ldiagac .AND. mype == 0) THEN
  WRITE(umMessage,'(A)') "***** Starting LHN_INC"
  CALL umPrint(umMessage,src='lhn_inc')
END IF
!
! Set parameters and variables
!
IF (l_lhn_search .AND. srange == 0) THEN
  IF (mype == 0) THEN
    WRITE(umMessage,'(A)') "WARNING : LHN_RANGE=0,"
    CALL umPrint(umMessage,src='lhn_inc')
    WRITE(umMessage,'(A)') "   therefore, setting L_LHN_SEARCH to FALSE"
    CALL umPrint(umMessage,src='lhn_inc')
  END IF
  l_lhn_search = .FALSE.
END IF


! Compute the global grid latitudes, DELTA_LAT_G and DELTA_LON_G if it is
! the first time the routine is called.
IF (first_call) THEN

  ! Check the domains are big enough for the srange halos if search is
  ! performed.
  IF (l_lhn_search) THEN
    
    ! Minimum domain size in E/W
    dsize  = row_length_g / nproc_x

    IF (mype==0 .AND. (dsize  <  srange)) THEN
      icode=5
      WRITE(cmessage,'(A,I3,A,I3,A,I3,A)')                              &
        'Too many processors in the East-West direction (',nproc_x,     &
        ') to support the lhn_range size (',srange,                      &
        '). Try running with ',(row_length_g/srange),                    &
        ' processors or set l_lhn_1A to .false..'
      
      CALL ereport(routinename,icode,cmessage)
    END IF

    ! Minimum domain size in N/S
    dsize  = rows_g / nproc_y

    IF (mype==0 .AND. (dsize  <  srange)) THEN
      icode=5
      WRITE(cmessage,'(A,I3,A,I3,A,I3,A)')                              &
        'Too many processors in the North-South direction (',nproc_y,   &
        ') to support the lhn_range size (',srange,                      &
        '). Try running with ',(rows_g/srange),                          &
        ' processors or set l_lhn_1A to .false..'
      
      CALL ereport(routinename,icode,cmessage)
    END IF
  END IF


  IF (.NOT. ALLOCATED(grid_lats_g)) ALLOCATE (grid_lats_g(rows_g))
  IF (.NOT. ALLOCATED(delta_lat_g)) ALLOCATE (delta_lat_g(rows_g))
  IF (.NOT. ALLOCATED(delta_lon_g)) ALLOCATE (delta_lon_g(row_length_g))

  IF (l_regular) THEN

    DO j = 1, rows_g
      grid_lats_g(j)=row1mgth+REAL(j-1)*dlatmg
    END DO
    
    delta_lat_g(:)=dlatmg
    delta_lon_g(:)=dlongmg

  ELSE    

    ! Extract one column of phi values into GRID_LATS_G which is
    ! inverted so colats run in correct direction.
    ! 
    ! The idea is to gather a single local column of phi_p from all
    ! the PEs of the same column group communicator into a single
    ! global column of phi_p which can be easily inverted.
    ! The same process from local to global is done for a row of lambda_p

    
    IF(.NOT. ALLOCATED(phi_p_g_col)) ALLOCATE (phi_p_g_col(rows_g))
    IF(.NOT. ALLOCATED(phi_p_l_col)) ALLOCATE (phi_p_l_col(rows))

    IF(.NOT. ALLOCATED(lambda_p_g_row)) ALLOCATE (lambda_p_g_row(row_length_g))
    IF(.NOT. ALLOCATED(lambda_p_l_row)) ALLOCATE (lambda_p_l_row(row_length))

    ! Copy information from PHI_P to the local single column buffer.
    DO j = 1, rows
      phi_p_l_col(j) = phi_p(1,j)
    END DO

    ! Copy information from LAMBDA_P to the local single row buffer.
    DO i= 1, row_length
      lambda_p_l_row(i) = lambda_p(i)
    END DO
    
    ! Compute the number of elements and displacement for the global
    ! column. 
    DO pe_j = 0, nproc_y -1
      pe = pe_j * nproc_x ! first column
      counts_y(pe_j)        = g_blsize(2, fld_type_p, pe)
      displacements_y(pe_j) = g_datastart(2, pe)-1
    END DO

    ! Compute the number of elements and displacement for the global
    ! row. 
    DO pe_i = 0, nproc_x -1
      pe = pe_i  ! first row
      counts_x(pe_i)        = g_blsize(1, fld_type_p, pe)
      displacements_x(pe_i) = g_datastart(1, pe)-1
    END DO

    ! Communicate the local single column buffer to the global single
    ! column buffer to all the PEs in the same column group communicator. 
    CALL mpl_allgatherv(phi_p_l_col, rows, mpl_real,                         &
                        phi_p_g_col, counts_y, displacements_y, mpl_real,    &
                        gc_proc_col_group, ierr)

    ! Communicate the local single row buffer to the global single
    ! row buffer to all the PEs in the same row group communicator. 
    CALL mpl_allgatherv(lambda_p_l_row, row_length, mpl_real,                &
                        lambda_p_g_row, counts_x, displacements_x, mpl_real, &
                        gc_proc_row_group, ierr)

    DO j = 1, rows_g
      grid_lats_g(j) = phi_p_g_col(rows_g-j+1)/pi_over_180
      grid_lats_g(j) = 90.0-grid_lats_g(j)
      grid_lats_g(j) = grid_lats_g(j)*pi_over_180
    END DO
    
    ! Compute DELTA_LAT_G
    DO j = 1, rows_g-1
      delta_lat_g(j) = phi_p_g_col(j+1)-phi_p_g_col(j)
    END DO
    delta_lat_g(rows_g) = delta_lat_g(rows_g-1)

    ! Compute DELTA_LONG_G
    DO j = 1, row_length_g-1
      delta_lon_g(j) = lambda_p_g_row(j+1) - lambda_p_g_row(j)
    END DO
    delta_lon_g(row_length_g) = delta_lon_g(row_length_g-1)
    
    DEALLOCATE (phi_p_g_col)
    DEALLOCATE (phi_p_l_col)    
    DEALLOCATE (lambda_p_g_row)
    DEALLOCATE (lambda_p_l_row)    

  END IF

  ! Set up search template used by lhn search
  IF (l_lhn_search) THEN
    
    ALLOCATE(search(4*srange*(srange+1),2))
    DO jr = 1 , srange
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
  END IF
  
  first_call = .FALSE.
END IF

IF (l_lhn_limit) limit = lhn_limit * timestep / rsec_per_day

fi_lhn    = fi_scale_lhn / planet_radius
no_incrs  = 0
no_search = 0
radius(1) = 0
radius(2) = 0
radius(3) = 0
radius(4) = 0
radius(5) = 0

!
! Initial calculations
!
! Precipitation variables
!
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(j,i,jpts_l)                      & 
!$OMP SHARED(rows, row_length, tot_pr, ls_rain, ls_snow, conv_rain,      &
!$OMP    conv_snow, anal_pr, pr_inc, tot_pr_r )
DO j = 1, rows
  DO i = 1, row_length
    jpts_l = (j-1) * row_length + i 
    
    tot_pr(jpts_l) = ls_rain(jpts_l) + ls_snow(jpts_l)+                   &
      conv_rain(jpts_l) + conv_snow(jpts_l)
    IF (tot_pr(jpts_l)  <   0.0) tot_pr(jpts_l) = 0.0
    anal_pr(jpts_l) = tot_pr(jpts_l) + pr_inc(jpts_l)
    IF (anal_pr(jpts_l)  <   0.0) THEN
      anal_pr(jpts_l) = 0.0
      pr_inc(jpts_l)   = anal_pr(jpts_l) - tot_pr(jpts_l)
    END IF

    ! Copy to field with srange halos
    tot_pr_r(i,j) = tot_pr(jpts_l)
  END DO
END DO
!$OMP END PARALLEL DO

IF (l_lhn_search) THEN
  ! Update the halo regions
  CALL swap_bounds(tot_pr_r,row_length,rows,1,srange,srange, &
                   fld_type_p,swap_field_is_scalar)
END IF

!
! Heating profiles, (FR Working Paper 171, eq2)
!

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,j,i,jpts_l)                    & 
!$OMP SHARED(model_levels, rows, row_length, tot_lh, conv_heat, ls_heat, &
!$OMP    timestep, remove_neg_lh, tot_lh_r )
DO k = 1 , model_levels
  DO j = 1, rows
    DO i = 1, row_length
      jpts_l = (j-1) * row_length + i 
        
      tot_lh(jpts_l,k) = ( conv_heat(jpts_l,k) +                         &
                           ls_heat(jpts_l,k) ) * timestep
      IF ( remove_neg_lh .AND. (tot_lh(jpts_l,k) <   0.0) ) THEN
        tot_lh(jpts_l,k) = 0.0
      END IF

      ! Copy to field with srange halos
      tot_lh_r(i,j,k) = tot_lh(jpts_l,k)
    END DO   
  END DO
END DO 
!$OMP END PARALLEL DO

! Update the halo regions
IF (l_lhn_search) THEN
  CALL swap_bounds(tot_lh_r,row_length,rows,model_levels,srange,srange, &
                   fld_type_p,swap_field_is_scalar)
END IF
!
! Calculate theta increments
!
! Set up array NEAR to decide where to get profile for scaling
! NEAR=-1 implies no scaling necessary, 0 implies scale here
! greater than 0 means scale a nearby profile, or not at all
!

IF (ltimer_ac) CALL timer('LHN_SRCH',3)

!  Calculate NEAR in the local LPG domain. However, the index returned
!  is a global one.
DO j = 1, rows
  DO i = 1, row_length
    jpts_l = (j-1) * row_length + i
    jpts_g = (j + g_datastart(2, mype) - 2) * row_length_g + &
             (i + g_datastart(1, mype) - 1)
    
    near(jpts_l) = -1      
    IF ( pr_inc(jpts_l)  /=  0.0 ) THEN
      near(jpts_l) = 0
      no_incrs     = no_incrs + 1
      IF ( tot_pr_r(i,j)  < (epsilon_lhn * anal_pr(jpts_l)) ) THEN
        near(jpts_l) = jpts_g  !!! Global index 
        no_search    = no_search + 1
        IF (l_lhn_search) THEN
          !  search for suitable profile
          CALL lhn_search_1a(i,j, near(jpts_l),                  &
            srange, search, rows, row_length,                       &
            rows_g, row_length_g,                                  &
            tot_pr_r, anal_pr(jpts_l), radius)
        END IF
      END IF
    END IF
  END DO
END DO

IF (ltimer_ac) CALL timer('LHN_SRCH',4)

!
! Calculate the increments, and scale by the relaxation coeff
! Attention NEAR contains the global index of the point

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,jpts_l,jpts_g,     & 
!$OMP     i_g,j_g,i_l,j_l)                                    &
!$OMP SHARED(model_levels,rows,row_length,mype,l_lhn_scale,   &
!$OMP     row_length_g,near,thincs,relax_cf_lhn,tot_lh,       &
!$OMP     pr_inc,tot_pr,tot_pr_r,tot_lh_r,anal_pr,            &
!$OMP     epsilon_lhn,alpha_lhn,p_field,l_lhn_fact)

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      jpts_l = (j-1) * row_length + i
      jpts_g = (j + g_datastart(2, mype) - 2) * row_length_g + &
               (i + g_datastart(1, mype) - 1)

      IF (near(jpts_l)  <   0) THEN
        !  No scaling necessary
        thincs(jpts_l,k) = 0.

      ELSE IF (near(jpts_l)  ==  0) THEN
        !  Scale the profile at the point itself
        !  (FR WP 171, eq 5)
        thincs(jpts_l,k) =                                      &
                 relax_cf_lhn * tot_lh(jpts_l,k) *              &
                 pr_inc(jpts_l) / tot_pr(jpts_l)

      ELSE IF (near(jpts_l)  >   0 .AND. near(jpts_l)  /=  jpts_g) THEN
        !  Scale a nearby profile
        !  (FR WP 171, eq 7)      
      
        ! Compute the global coordinate of nearby profile 
        j_g = INT ((near(jpts_l)-1)/row_length_g) + 1
        i_g = near(jpts_l) - ((j_g-1) * row_length_g )

        ! Compute the local coordinate of nearby profile
        j_l = j_g - (g_datastart(2, mype) - 1)
        i_l = i_g - (g_datastart(1, mype) - 1)
              
        thincs(jpts_l,k) = relax_cf_lhn *                       &
          (anal_pr(jpts_l)/tot_pr_r(i_l,j_l) *                  &
          tot_lh_r(i_l,j_l,k) -                                 & 
          tot_lh(jpts_l,k))

      ELSE IF (near(jpts_l)  ==  jpts_g .AND. l_lhn_scale) THEN
        !  Scale by 1/EPSILON_LHN if no suitable profile available
        !  (FR WP 171, eq 8)
        thincs(jpts_l,k) = relax_cf_lhn *                       &
          (( 1.0/epsilon_lhn) - 1.0) * tot_lh(jpts_l,k)
      ELSE
        !  If none of the above apply, ignore the ob
        thincs(jpts_l,k) = 0.0
      END IF
    END DO ! I
  END DO ! J
END DO ! K
!$OMP END DO NOWAIT

!
!  Impose limit by factor of alpha
!
IF (l_lhn_fact) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1 , model_levels
    DO jpts_l = 1, p_field
      IF (pr_inc(jpts_l)  <   0.0) THEN
        IF ((pr_inc(jpts_l)/tot_pr(jpts_l)) <  (alpha_lhn-1.0)) THEN
          thincs(jpts_l,k) = (alpha_lhn - 1.0) *              &
                             tot_lh(jpts_l,k) * relax_cf_lhn
        END IF
      END IF
    END DO   
  END DO     
!$OMP END DO NOWAIT
END IF

!$OMP END PARALLEL

! Recursive filtering of increments
IF (l_lhn_filt) THEN
  
  !       Initialise COS_LAT
  IF (l_regular) THEN
    z = row1mgth
    DO j = 1, rows_g
      cos_lat(j) = SIN(z)
      z = z + dlatmg
    END DO  ! JPTS
  ELSE
    DO j = 1, rows_g
      cos_lat(j)=SIN(grid_lats_g(j))
    END DO  ! JPTS
  END IF


  m_grid = 0

  CALL rfcsl_1a(thincs, rows, row_length, model_levels,         &
                rows_g, row_length_g,  m_grid,                  &
                0.0, cos_lat, delta_lat_g, delta_lon_g, fi_lhn, &
                npass_rf_lhn)
  
END IF   ! L_LHN_FILT


!
! Impose absolute limit
!
IF (l_lhn_limit) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(k,jpts_l)           & 
!$OMP SHARED(model_levels, p_field, thincs, limit)
  DO k = 1 , model_levels
    DO jpts_l = 1 , p_field
      IF (thincs(jpts_l,k)  <   -1.0 * limit)               &
             thincs(jpts_l,k) = -1.0 * limit
      IF (thincs(jpts_l,k)  >   limit)                      &
             thincs(jpts_l,k) = limit
    END DO  
  END DO
!$OMP END PARALLEL DO
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
    WRITE(umMessage,*) "                     : LHN_RANGE   = ",srange
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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lhn_inc_1a
END MODULE lhn_inc_1a_mod

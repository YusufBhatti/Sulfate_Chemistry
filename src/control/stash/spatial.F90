! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: SPATIAL --------------------------------------------------
!
!    Purpose: Performs general spatial processing on an input field to
!             produce an output field or scalar within STASH.  Lower-
!             level routines are called to perform the various spatial
!             processing options.
!
!
!    Programming standard: UM Doc Paper 3
!
!    External documentation:
!      Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                   system (STASH)
!
!  -----------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

SUBROUTINE spatial(fieldin,vx,vy,vz,gr,st_grid,                   &
                   fld_type,halo_type,                            &
                   halo_x,halo_y,                                 &
                   lcyclic,lmasswt,                               &
      n_cols_out,n_rows_out,                                      &
      this_index_lev,level_list,index_lev,no_of_levels,           &
      no_of_levels_masswt,                                        &
      p,pstar,                                                    &
      cos_v_latitude,cos_theta_latitude,land,sea,                 &
      row_length,rows,n_rows,no_rows,model_levels,                &
      fieldout,lenout,                                            &
      control,control_size,rmdi,                                  &
      icode,cmessage)
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore,  ONLY: mype
USE Field_Types, ONLY: fld_type_p
USE atm_fields_bounds_mod, ONLY: vdims_s, tdims, tdims_s, pdims,  &
                                 pdims_s
USE stparam_mod, ONLY: st_west_code, st_east_code, st_south_code, &
                st_north_code, st_gridpoint_code, st_input_bottom,&
               st_weight_code, st_tp_grid, st_uv_grid, st_cu_grid,&
               st_cv_grid, st_riv_grid, stash_null_mask_code,     &
               stash_weight_null_code, stash_weight_area_code,    &
               stash_weight_volume_code, stash_weight_mass_code,  &
               stash_land_mask_code, stash_sea_mask_code,         &
               stash_nmdi_mask_code, extract_top, extract_base,   &
               vert_mean_top, vert_mean_base, zonal_mean_top,     &
               zonal_mean_base, merid_mean_top, merid_mean_base,  &
               field_mean_top, field_mean_base, global_mean_top,  &
               global_mean_base, block_size
USE sterr_mod, ONLY: st_upper_less_lower, st_not_supported,       &
                     st_no_data,st_nd, st_bad_array_param,        &
                     st_bad_address, st_unknown,                  &
                     st_bad_wraparound, st_illegal_weight,        &
                     unknown_weight, unknown_mask,                &
                     unknown_processing, nonsense
USE mpp_conf_mod,  ONLY: swap_field_is_scalar,                    &
                          include_halos_ew, include_halos_ns

USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umPrint, umMessage
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!

INTEGER ::                                                        &
    vx,vy,vz,                                                     &
                                      ! IN size of fieldin
    lenout,                                                       &
                                      ! IN size of fieldout
    gr,                                                           &
                                      ! IN ppxref gridtype code
    st_grid,                                                      &
                                      ! IN STASH gridtype code
    fld_type,                                                     &
                                      ! IN field type (u/v/p)
    halo_type,                                                    &
                                      ! IN halo type
    halo_x,                                                       &
                                      ! IN EW halo of input
    halo_y,                                                       &
                                      ! IN NS halo of input
    n_rows_out,                                                   &
                                      ! OUT no. of output rows
    n_cols_out,                                                   &
                                      ! OUT no. of output cols
    this_index_lev,                                               &
                                      ! IN level index, this field
    row_length,rows,n_rows,                                       &
                                      ! IN horiz. sizes (C grid)
    no_rows,                                                      &
                                      ! IN number of rows used
    model_levels,                                                 &
                                      ! IN vertical size
    control_size,                                                 &
                                      ! IN size of control record
    control(control_size),                                        &
                                      ! IN control record
    icode,                                                        &
                                      ! OUT error code 0 if ok
    no_of_levels,                                                 &
                                      ! IN no of levels
    no_of_levels_masswt,                                          &
                          ! IN levels for mass weighting array
                          ! lmasswt F: =1; T: =no_of_levels
    index_lev(no_of_levels),                                      &
                                      ! IN index to levels
    level_list(no_of_levels)          ! IN model level list

REAL ::                                                           &
    fieldin(vx,vy,vz),                                            &
                       ! IN fieldin which is to be acted on
    p(pdims_s%i_start:pdims_s%i_end,                              &
      pdims_s%j_start:pdims_s%j_end,                              &
      pdims_s%k_start:pdims_s%k_end),                             &
                                      ! IN pressure (rho levels)
    pstar(pdims%i_start:pdims%i_end, pdims%j_start:pdims%j_end),  &
                                      ! IN surface pressure
    cos_v_latitude(vdims_s%i_start:vdims_s%i_end,                 &
                   vdims_s%j_start:vdims_s%j_end),                &
                                      ! IN v-grid area fn
   cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,              &
                      tdims_s%j_start:tdims_s%j_end),             &
                                      ! IN T-grid area fn 
    fieldout(lenout),                                             &
                                          ! OUT output field
    rmdi                                  ! IN  missing data indic

LOGICAL ::                                                        &
    lcyclic,                                                      &
                                          ! IN .true. if cyclic EW
    lmasswt,                                                      &
                                          ! IN  TRUE if masswts OK
    land(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),   &
                                          ! IN land mask
    sea(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
                                          ! IN sea mask
 
CHARACTER(LEN=errormessagelength) :: cmessage   ! OUT error message

! ----------------------------------------------------------------------
!
!  local variables
!
LOGICAL :: lwrap                ! TRUE if output field wraparound EW
LOGICAL :: lmasswt_strict       ! copy of lmasswt - but set to false
!                                  ! if mass weighting is not requested
INTEGER :: xstart,ystart        ! lower left hand corner coords
INTEGER :: xend,yend            ! upper right hand corner coords
INTEGER :: processing_code      ! what kind of mean  will be done
INTEGER :: what_level           ! what type of input level
INTEGER :: what_mask            ! what mask is used
INTEGER :: what_weight          ! what weighting is used

INTEGER :: i,j,k                                                     &
                             ! loop counters
,model_level                                                      &
                             ! model level
,this_index_lev_masswt       ! level index for mass weighting
                             ! (=1 if no mass weighting or no
                             !  model level weights available)

INTEGER ::                                                        &
! global versions of the extracted area domain limits
        global_xstart,global_xend,global_ystart,global_yend

! workspace arrays containining weighting factors and masks.
REAL ::                                                           &
  area_weight(1-halo_x:row_length+halo_x+1,                       &
              1-halo_y:no_rows+halo_y)                            &
, pstar_weight(1-halo_x:row_length+halo_x+1,                      &
               1-halo_y:no_rows+halo_y,                           &
               no_of_levels_masswt)                               &
, pstar_interp(row_length,no_rows)
LOGICAL ::                                                        &
  mask(1-halo_x:row_length+halo_x+1,                              &
       1-halo_y:no_rows+halo_y)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SPATIAL'


! ----------------------------------------------------------------------
!  1. Set up local variables
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
xstart=control(st_west_code)
xend=control(st_east_code)
ystart=control(st_south_code)  ! NOTE: Grid is assumed to be
yend=control(st_north_code)    !       oriented south-to-north

global_xstart=xstart
global_ystart=ystart
global_xend=xend
global_yend=yend

! and calculate what the local subdomain limits are:
! DEPENDS ON: global_to_local_subdomain
CALL global_to_local_subdomain( include_halos_ew,include_halos_ns, &
                                gr,halo_type,mype,              &
                                global_ystart,global_xend,      &
                                global_yend,global_xstart,      &
                                ystart,xend,yend,xstart)

! If time series geographic point not on this processor, where
! coordinates N==S and W==E then NO DATA here, skip the rest.
! END IF will be just before label 9999
processing_code=control(st_gridpoint_code)
IF ((xstart == st_no_data .OR. xend == st_no_data .OR.  &
    ystart == st_no_data .OR. yend == st_no_data) .AND. &
   (processing_code <  global_mean_top   .AND.          &
    processing_code >  global_mean_base) .AND.          &
    lenout >= 1)  THEN
  DO i = 1, lenout
    fieldout(i) = 0.0
  END DO
  n_rows_out=1
  n_cols_out=1
ELSE

! Check if wraparound field
IF (xstart >  xend) THEN

  IF (lcyclic) THEN
    xend=xend+blsize(1,fld_type)
    ! subtract two halos as we don't wish to include halos at the end
    ! and start of the row within the wrap around domain
    lwrap=.TRUE.
  ELSE
    icode=st_bad_wraparound    ! wraparound illegal unless cyclic
    GO TO 9999
  END IF

ELSE
  lwrap=.FALSE.
END IF

IF (global_xstart  >   global_xend) THEN
  IF (lcyclic) THEN
    global_xend=global_xend+glsize(1,fld_type)
  ELSE
    icode=st_bad_wraparound  ! wraparound illegal unless cyclic
    GO TO 9999
  END IF
END IF

! Processing is the 8th element of array, PARAMETER st_gridpoint_code=8
! but in multi_spatial fake_record = control is set manually to
! fake_record(st_gridpoint_code)=                                
!      (stash_series(series_proc_code,i)/block_size)*block_size  
!       + what_mask  (51/10)*10 = 51
what_level=control(st_input_bottom)
what_mask=MOD(processing_code,block_size)
what_weight=control(st_weight_code)

!
!  1.1 Prevent masking or weighting if input field is not 2D in extent
!      - weighting and masking is assumed to have been done outside.
!
IF ( (.NOT. (st_grid == st_tp_grid .OR. st_grid == st_uv_grid .OR.    &
            st_grid == st_cu_grid .OR. st_grid == st_cv_grid .OR.    &
             st_grid == st_riv_grid))                             &
      .AND. (what_mask   /= stash_null_mask_code   .OR.            &
            what_weight /= stash_weight_null_code) ) THEN
  icode=st_not_supported
  cmessage='SPATIAL : Masking/weighting unsupported - non 2D field'
  GO TO 9999
END IF

! Check for supported weighting and masking options

IF (.NOT. ((what_weight  ==  stash_weight_null_code) .OR.         &
           (what_weight  ==  stash_weight_area_code) .OR.         &
           (what_weight  ==  stash_weight_volume_code) .OR.       &
           (what_weight  ==  stash_weight_mass_code) ) ) THEN
  cmessage='SPATIAL : Unrecognized weighting option'
  icode=unknown_weight
  GO TO 9999
END IF

IF (.NOT. ((what_mask  ==  stash_null_mask_code) .OR.             &
           (what_mask  ==  stash_land_mask_code) .OR.             &
           (what_mask  ==  stash_sea_mask_code) .OR.              &
           (what_mask  ==  stash_nmdi_mask_code ) ) ) THEN
  cmessage='SPATIAL : Unrecognized masking option'
  icode=unknown_mask
  GO TO 9999
END IF

IF (what_weight  ==  stash_weight_volume_code) THEN
  cmessage='SPATIAL : Volume-weighting not supported'
  icode=st_illegal_weight
  GO TO 9999
END IF

! Set lmasswt_strict - copy of lmasswt, but set to false is mass
! weighting not requested

lmasswt_strict=                                                   &
  (lmasswt .AND. (what_weight  ==  stash_weight_mass_code))

! Precalculate weighting and masking arrays
! I've used IF tests inside the loops, but since the logical
! expressions are invariant wrt i and j, the compiler will
! move them outside the DO loops. It makes the code a lot shorter!

! NOTE that neither area-weights or mass-weights are normalised so
! that the interpretation of weighted diagnostics is non-obvious. Also
! the PP lookup header has no switch for indicating whether or not such
! processing has taken place. More careful treatment of horizontal
! interpolation is required for stricter accuracy.


!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,model_level)                       &
!$OMP SHARED(what_weight, lmasswt, level_list, p, cos_v_latitude,           &
!$OMP        this_index_lev_masswt, no_rows, row_length, st_grid,           &
!$OMP        pstar_weight, pstar_interp, pstar, cos_theta_latitude,         &
!$OMP        this_index_lev, no_of_levels_masswt, model_levels, halo_x,     &
!$OMP        halo_y, mask, what_mask, land, sea, area_weight)

! area weighting
IF (what_weight  ==  stash_weight_null_code) THEN
    ! Ensure initialisation of weight array including halos
!$OMP DO SCHEDULE(STATIC)
  DO j = 1-halo_y, no_rows+halo_y           
    DO i = 1-halo_x, row_length+halo_x+1
      area_weight (i,j)   = 1.0
    END DO
  END DO
!$OMP END DO
ELSE ! some form of area weighting will be required
!$OMP DO SCHEDULE(STATIC)
  DO j=1,no_rows
    DO i=1,row_length
      IF (st_grid  ==  st_cv_grid) THEN
        area_weight(i,j)=cos_v_latitude(i,j)
      ELSE
        ! NOTE that this is will not be accurate for C-u grid variables for
        ! LAMs, since cos_theta_latitude will differ between theta,u positions
        area_weight(i,j)=cos_theta_latitude(i,j)
      END IF
    END DO ! i
  END DO ! j
!$OMP END DO

END IF    ! what_weight


! mass weighting
IF ((what_weight  ==  stash_weight_null_code) .OR.                &
    (what_weight  ==  stash_weight_area_code)) THEN
  ! No mass weighting is required
!$OMP SINGLE
  this_index_lev_masswt = 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO j = 1-halo_y, no_rows+halo_y
    DO i = 1-halo_x, row_length+halo_x+1
      ! Ensure initialisation of weight array including halos
      pstar_weight(i,j,this_index_lev_masswt) = 1.0
    END DO
  END DO
!$OMP END DO
ELSE

  ! Mass weighting requested

  ! Ensure that halos are initialised
!$OMP DO SCHEDULE(STATIC)
  DO j=no_rows-1,no_rows
    DO i=1,row_length
      pstar_interp(i,j) =1.0
    END DO !i
  END DO !j
!$OMP END DO

  ! Interpolate pstar to correct horizontal staggering
  IF (st_grid  ==  st_cu_grid) THEN
    ! NOT YET CORRECT: pstar has no halos! So set to nearby value
    !          CALL P_TO_U(pstar,row_length,rows,1,0,0,pstar_interp)
!$OMP DO SCHEDULE(STATIC)
    DO j=1,no_rows
      DO i=1,row_length
        pstar_interp(i,j) =pstar(i,j)
      END DO !i
    END DO !j
!$OMP END DO

  ELSE IF (st_grid  ==  st_cv_grid) THEN
    ! NOT YET CORRECT: pstar has no halos! So set to nearby value
    !          CALL P_TO_V(pstar,row_length,rows,1,0,0,pstar_interp)
!$OMP DO SCHEDULE(STATIC)
    DO j=1,no_rows
      DO i=1,row_length
        pstar_interp(i,j) =pstar(i,j)
      END DO !i
    END DO !j
!$OMP END DO

  ELSE
!$OMP DO SCHEDULE(STATIC)
    DO j=1,no_rows
      DO i=1,row_length
        pstar_interp(i,j)=pstar(i,j)
      END DO
    END DO
!$OMP END DO
  END IF

  IF (lmasswt) THEN  ! model level mass weights available

!$OMP SINGLE
    this_index_lev_masswt = this_index_lev
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
    DO k=1,no_of_levels_masswt
      model_level = level_list(k)
      IF (model_level == model_levels) THEN  ! top level

        DO j=1,no_rows
          DO i=1,row_length
            pstar_weight(i,j,k) = p(i,j,model_level)
          END DO !i
        END DO !j

      ELSE      ! not top level

        DO j=1,no_rows
          DO i=1,row_length
            ! Only accurate for variables on theta levels
            pstar_weight(i,j,k) = p(i,j,model_level) -              &
                                    p(i,j,model_level+1)
          END DO !i
        END DO !j
      END IF

    END DO ! k
!$OMP END DO

  ELSE              ! no model level mass weights available:
                    ! weight by pstar

!$OMP SINGLE
    this_index_lev_masswt = 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
    DO j=1,no_rows
      DO i=1,row_length
        pstar_weight(i,j,this_index_lev_masswt) =                &
                                              pstar_interp(i,j)
      END DO !i
    END DO !j
!$OMP END DO

  END IF   ! lmasswt
END IF


! Masking
! Ensure initialisation of mask array including halos
!$OMP DO SCHEDULE(STATIC)
DO j = 1-halo_y, no_rows+halo_y
  DO i = 1-halo_x, row_length+halo_x+1
    mask (i,j)   = .TRUE.
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO j=1,no_rows
  DO i=1,row_length
    IF (what_mask  ==  stash_land_mask_code) THEN
      mask(i,j)=land(i,j)
    ELSE IF (what_mask  ==  stash_sea_mask_code) THEN
      mask(i,j)=sea(i,j)
    ELSE
      mask(i,j)=.TRUE.
    END IF
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

! Update the halo areas of the weighting/masking arrays

!  [Note that for lams at UM5.2 and before, valid values for rim
!   boundaries would be needed using FILL_EXTERNAL_HALOS calls for
!   area_weight and pstar_weight arrays, but this should now be
!   superseded by initialising full arrays]

      ! Update halos only if halos present (standard diagnostics have no
      ! halos, and if some weighting required (probably defunct
      ! functionality)
IF ( (halo_x  /=  0 .OR. halo_y  /=  0) .AND.                      &
    (what_weight  /=  stash_weight_null_code .OR.                 &
     what_mask    /=  stash_null_mask_code )     ) THEN

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( area_weight(1-halo_x:row_length+halo_x,:),      &
       row_length,no_rows,1,halo_x,halo_y,                          &
                  fld_type_p, swap_field_is_scalar)
  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( pstar_weight(1-halo_x:row_length+halo_x,:,:),   &
       row_length,no_rows,                                          &
                  no_of_levels_masswt,halo_x,halo_y,                &
                  fld_type_p, swap_field_is_scalar)
  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( mask(1-halo_x:row_length+halo_x,:),             &
       row_length,no_rows,1,halo_x,halo_y,                          &
                  fld_type_p, swap_field_is_scalar)
END IF
! ----------------------------------------------------------------------
!  2. Call service routine to perform required processing
!
!  2.1 Extract sub field (single level at a time)
!
IF (processing_code <  extract_top .AND.                           &
    processing_code >= extract_base) THEN
  n_rows_out=(yend+1)-ystart
  n_cols_out=(xend+1)-xstart

  IF (                                                            &
   (xstart  /=  st_no_data) .AND. (xend  /=  st_no_data) .AND.    &
   (ystart  /=  st_no_data) .AND. (yend  /=  st_no_data)) THEN

    ! DEPENDS ON: stextc
    CALL stextc(fieldin,vx,vy,fld_type,halo_type,                   &
                lwrap,lmasswt_strict,                               &
                xstart,ystart,xend,yend,                            &
                fieldout,                                           &
                pstar_weight(:,:,this_index_lev_masswt),            &
                area_weight,mask,                                   &
                what_level,what_mask,what_weight,rmdi,              &
                icode,cmessage)

  ELSE  ! just set values to non NaN
    DO i=1,lenout
      fieldout(i)=0.0
    END DO
  END IF

  !
  !  2.2 Calculate column mean (over multiple levels indexed by index_lev)
  !
ELSE IF (processing_code <  vert_mean_top .AND.                    &
        processing_code >  vert_mean_base) THEN
  n_rows_out=yend+1-ystart
  n_cols_out=xend+1-xstart

  IF (                                                            &
   (xstart  /=  st_no_data) .AND. (xend  /=  st_no_data) .AND.    &
   (ystart  /=  st_no_data) .AND. (yend  /=  st_no_data)) THEN

    ! DEPENDS ON: stcolm
    CALL stcolm(fieldin,vx,vy,vz,fld_type,halo_type,                &
                lwrap,lmasswt_strict,                               &
                xstart,ystart,xend,yend,                            &
                fieldout,index_lev,no_of_levels,                    &
                pstar_weight,                                       &
                area_weight,mask,                                   &
                what_level,what_mask,what_weight,rmdi,              &
                icode,cmessage)

  ELSE  ! just set values to non NaN
    DO i=1,lenout
      fieldout(i)=0.0
    END DO
  END IF

  !
  !  2.3 Calculate zonal mean (single level at a time)
  !
ELSE IF (processing_code <  zonal_mean_top .AND.                   &
        processing_code >  zonal_mean_base) THEN
  n_rows_out=yend+1-ystart
  n_cols_out=1

  ! DEPENDS ON: stzonm
  CALL stzonm(fieldin,vx,vy,fld_type,gr,halo_type,                &
              lwrap,lmasswt_strict,                               &
              xstart,ystart,xend,yend,                            &
              global_xstart,global_ystart,global_xend,global_yend,&
              fieldout,                                           &
              pstar_weight(:,:,this_index_lev_masswt),            &
              area_weight,mask,                                   &
              what_level,what_mask,what_weight,rmdi,              &
              icode,cmessage)
  !
  !  2.4 Calculate meridional mean (single level at a time)
  !
ELSE IF (processing_code <  merid_mean_top .AND.                   &
        processing_code >  merid_mean_base) THEN
  n_rows_out=1
  n_cols_out=xend+1-xstart
  ! DEPENDS ON: stmerm

  CALL stmerm(fieldin,vx,vy,fld_type,gr,halo_type,                &
              lwrap,lmasswt_strict,                               &
              xstart,ystart,xend,yend,                            &
              global_xstart,global_ystart,global_xend,global_yend,&
              fieldout,                                           &
              pstar_weight(:,:,this_index_lev_masswt),            &
              area_weight,mask,                                   &
              what_level,what_mask,what_weight,rmdi,              &
              icode,cmessage)
  !
  !  2.5 Calculate field mean (single level at a time)
  !
ELSE IF (processing_code <  field_mean_top .AND.                   &
        processing_code >  field_mean_base) THEN
  n_rows_out=1
  n_cols_out=1

  ! DEPENDS ON: stfieldm
  CALL stfieldm(fieldin,vx,vy,fld_type,gr,halo_type,              &
                lwrap,lmasswt_strict,                             &
              xstart,ystart,xend,yend,                            &
              global_xstart,global_ystart,global_xend,global_yend,&
              fieldout,                                           &
              pstar_weight(:,:,this_index_lev_masswt),            &
              area_weight,mask,                                   &
              what_level,what_mask,what_weight,rmdi,              &
              icode,cmessage)
  !
  !  2.6 Calculate global mean (over multiple levels)
  !
ELSE IF (processing_code <  global_mean_top .AND.                  &
        processing_code >  global_mean_base) THEN
  n_rows_out=1
  n_cols_out=1
  IF (                                                            &
    (xstart  /=  st_no_data) .AND. (xend  /=  st_no_data) .AND.   &
    (ystart  /=  st_no_data) .AND. (yend  /=  st_no_data)) THEN
    ! DEPENDS ON: stglom
    CALL stglom(fieldin,vx,vy,vz,fld_type,gr,halo_type,           &
              lwrap,lmasswt_strict,                               &
              xstart,ystart,xend,yend,                            &
              global_xstart,global_ystart,global_xend,global_yend,&
              fieldout,index_lev,no_of_levels,                    &
              pstar_weight,                                       &
              area_weight,mask,                                   &
              what_level,what_mask,what_weight,rmdi,              &
              icode,cmessage)
  ELSE
    DO i=1,lenout
      fieldout(i)= 0.0
    END DO
  END IF

  !
  !  2.7 Invalid processing option
  !
ELSE
  icode=unknown_processing
  WRITE(cmessage,'(A,1x,i5)')'SPATIAL: unknown processing option',&
    processing_code
  CALL ereport('SPATIAL',icode,cmessage)
END IF
!
! end of NO DATA on this processor
!
END IF

9999 CONTINUE
IF (icode >  0) THEN
  WRITE(umMessage,'(A,A,I5)')                                     &
  'SPATIAL: Error in spatial processing:',', code ',icode
  CALL umPrint(umMessage,src='spatial')
  CALL ereport('SPATIAL',icode,cmessage)
END IF
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE spatial

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE diagdesc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGDESC_MOD'

CONTAINS
!
!  Subroutine DIAGDESC -----------------------------------------------
!
!  Purpose: Prints a formatted diagnostic description using the name
!           of a diagnostic plus its PPXREF and STASH record.  Gives
!           a hardcopy record of the diagnostics included in a run.
!
!  Programming standard: UM Doc Paper 3
!
!  External documentation:
!  Unified Model Doc Paper C4 - Storage handling and diagnostic
!                               system (STASH)
! --------------------------------------------------------------------
!
! Interface and arguments: -------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Misc
!---------------------------------------------------------------------
SUBROUTINE diagdesc(seqno,NAME,stlist,ppxref,                           &
              stash_levels,num_stash_levels,num_level_lists,            &
              stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,    &
              sttabl,nsttims,nsttabl,                                   &
              stash_series,stash_series_rec_len,stash_series_len,       &
              stash_series_index,stash_ser_index_size)

USE UM_ParCore, ONLY: mype
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE cppxref_mod
USE submodel_mod, ONLY: atmos_im
USE file_manager, ONLY: get_file_unit_by_id, um_file_type, get_file_by_unit
USE nlstgen_mod, ONLY: ppxm, l_meaning_sequence
USE stparam_mod, ONLY: st_proc_no_code, st_sect_no_code, st_output_code,&
    st_macrotag, st_replace_code, st_time_series_code, st_accum_code,   &
    st_time_mean_code, st_series_ptr, st_max_code, st_min_code,         &
    st_append_traj_code, st_variance_code, st_time_series_mean,         &
    st_freq_code, st_start_time_code, st_end_time_code,st_infinite_time,&
    st_period_code, st_input_code, st_gridpoint_code, extract_base,     &
    extract_top, vert_mean_base, vert_mean_top, zonal_mean_base,        &
    zonal_mean_top, merid_mean_base, merid_mean_top, field_mean_base,   &
    field_mean_top, global_mean_base, global_mean_top, st_output_bottom,&
    st_special_code, st_output_top, st_pseudo_out, st_south_code,       &
    st_north_code, st_west_code, st_east_code, st_weight_code,          &
    stash_weight_null_code, stash_weight_area_code,                     &
    stash_weight_volume_code, stash_weight_mass_code, st_gridpoint_code,&
    stash_null_mask_code, stash_land_mask_code, stash_sea_mask_code,    &
    stash_nmdi_mask_code, st_item_code, st_end_of_list, block_size,     &
    st_output_type, st_netcdf
                 
USE umPrintMgr, ONLY:      &
    umPrint,               &
    umMessage

IMPLICIT NONE

!

CHARACTER(LEN=36), INTENT(IN)  ::  NAME        ! IN  diagnostic name
INTEGER, INTENT(IN)         ::                                          &
                                seqno,                                  &
                                            ! IN  sequence number
                                stlist(*),                              &
                                            ! IN  STASHlist record
                                ppxref(*)   ! IN  PPXREF record

! STASH levels list information
INTEGER, INTENT(IN)         ::  num_stash_levels                        &
                                            ! IN Max levels in a list
                               ,num_level_lists                         &
                                            ! IN Number of lists
                    ,stash_levels(num_stash_levels+1,num_level_lists)

! STASH pseudo-levels list information
INTEGER, INTENT(IN)         ::  num_stash_pseudo                        &
                                            ! IN Max ps-levs in a list
                               ,num_pseudo_lists                        &
                                            ! IN No of ps-lev lists
            ,stash_pseudo_levels(num_stash_pseudo+1,num_pseudo_lists)

! STASH time list information
INTEGER, INTENT(IN)         ::  nsttims                                 &
                                            ! IN Max times in a list
                               ,nsttabl                                 &
                                            ! IN Number of lists
                               ,sttabl(nsttims,nsttabl)

! STASH timeseries information
INTEGER, INTENT(IN)         ::  stash_series_len                        &
                                            ! IN Total no of records
                               ,stash_series_rec_len                    &
                                            ! IN Length of each record
                 ,stash_series(stash_series_rec_len,stash_series_len)   &
                                            ! IN array of records
                               ,stash_ser_index_size                    &
                                            ! IN No of index records
                          ,stash_series_index(2,stash_ser_index_size)
!
! Local variables
!
CHARACTER(LEN=80)  :: ch                ! Working character string variable

CHARACTER(LEN=8)   :: submodel,datatype,grid_type
CHARACTER(LEN=9)   :: level_type,freq,weight
CHARACTER(LEN=7)   :: packing,mask
CHARACTER(LEN=15)  :: time_proc,pseudo,levels
CHARACTER(LEN=6)   :: from, TO,period
CHARACTER(LEN=10)  :: source
CHARACTER(LEN=17)  :: dest
CHARACTER(LEN=18)  :: spatial
CHARACTER(LEN=27)  :: horiz

INTEGER         :: i1,i2,k           ! Array indices
INTEGER         :: j                 ! Code value
INTEGER         :: time_list,lev_list                                   &
                                     ! pointers to time and levels lists
                  ,plev_list                                            &
                                     ! pointer  to pseudo-level list
                  ,tser_list         ! pointer  to time series record list
INTEGER         :: ntimes            ! no of times in a time list
INTEGER         :: packing_profile   ! packing profile for output PPfield
INTEGER         :: i_unit            ! file unit number
INTEGER         :: stash_unit        ! stash requests log file unit

TYPE(um_file_type), POINTER :: um_file

! Header
CHARACTER(LEN=*), PARAMETER ::                                             &
 stars='   ********************************************************'    &
,starp='   **                                                    **'    &
,list ='   **    LIST OF USER-DEFINED DIAGNOSTICS IN THIS RUN    **'    &
,notes='   ** NOTES:                                             **'    &
,note1='   **   Time processing details are in timesteps, where  **'    &
,note2='   **     ... represents "for ever".                     **'    &
,note3='   **   Spatial processing domain is in gridpoints.      **'    &
,headend='==========================================================='//&
    '================================================================'

CHARACTER(LEN=*), PARAMETER ::                                             &
    diag1=' #No Diagnostic Description-------------- Submodel Item Section '//&
    'PPfcode Datatype Gridtype Leveltype MetO8lv MetO8fc Packacc',     &
    diag2=' Time-processing -From- --To-- Frequency Period --Source-- ---De'//&
    'stination---                                               ',     &
    diag3=' Spatial-processing -Levels-domain- -Pseudo-levels- -----Horizon'//&
    'tal-domain----- Weighting Masking                          '



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DIAGDESC'

!
!----------------------------------------------------------------------
! 0. Write header if sequence no indicates first item
!

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0) THEN
  stash_unit = get_file_unit_by_id("stash_req_log", handler="fortran")
END IF

IF (seqno == 1 .AND. mype == 0) THEN
  WRITE(stash_unit,'(15(a,/),/,a)')                                            &
       stars,stars,starp,list,starp,stars,stars,starp,                  &
       notes,note1,note2,note3,starp,stars,stars,headend
END IF

!----------------------------------------------------------------------
! 1. For each diagnostic processing request in the STASHlist,
!    print the diagnostic name followed by a summary of the processing
!    information on 3 lines.
!
! 1.0 If diagnostic is not required for output, exit routine
!
IF (stlist(st_proc_no_code) /= 0) THEN
  !
  ! 1.1 Line 1.
  !
  ! #No               -- supports up to 999 diagnostic fields, i3.

  ! Name

  ! Submodel
  j=stlist(st_sect_no_code)
  IF (ppxref(ppx_model_number) == atmos_im) THEN
    submodel=' ATMOS  '
  ELSE
    CALL umPrint(' Error in DIAGDESC: Unknown model',src='diagdesc')
    submodel=' UNKNOWN'
  END IF

  ! Item

  ! Section

  ! PPfcode

  ! Datatype
  j=ppxref(ppx_data_type)
  IF (j == 1 .OR. j == 4) THEN
    datatype='  REAL  '
  ELSE IF (j == 2 .OR. j == 5) THEN
    datatype='INTEGER '
  ELSE IF (j == 3) THEN
    datatype='LOGICAL '
  ELSE
    datatype='UNKNOWN '
  END IF

  ! Gridtype
  j=ppxref(ppx_grid_type)
  IF (j == ppx_atm_nonstd) THEN
    grid_type=' NONSTD '
  ELSE IF ((j >  ppx_atm_nonstd .AND. j <= ppx_atm_tsea) .OR.         &
            j == ppx_atm_compressed .OR. j == ppx_atm_ozone) THEN
    grid_type =' P-GRID '
  ELSE IF (j >= ppx_atm_uall .AND. j <= ppx_atm_usea) THEN
    grid_type=' UV-GRID'
  ELSE IF (j == ppx_atm_cuall) THEN
    grid_type=' CU-GRID'
  ELSE IF (j == ppx_atm_cvall) THEN
    grid_type=' CV-GRID'
  ELSE IF (j == ppx_atm_tzonal) THEN
    grid_type=' PZ-GRID'
  ELSE IF (j == ppx_atm_uzonal) THEN
    grid_type=' UZ-GRID'
  ELSE IF (j == ppx_atm_tmerid) THEN
    grid_type=' PM-GRID'
  ELSE IF (j == ppx_atm_umerid) THEN
    grid_type=' UM-GRID'
  ELSE IF (j == ppx_atm_rim) THEN
    grid_type='   RIM  '
  ELSE IF (j == ppx_atm_scalar) THEN
    grid_type=' SCALAR '
  ELSE
    grid_type=' UNKNOWN'
  END IF

  ! Leveltype
  j=ppxref(ppx_lv_code)
  IF (j == ppx_full_level) THEN
    level_type='FULLLEVEL'
  ELSE IF (j == ppx_half_level) THEN
    level_type='HALFLEVEL'
  ELSE
    level_type='STD-LEVEL'
  END IF

  ! Meto8LV

  ! Meto8FC

  ! PackAcc
  j=stlist(st_output_code)
  IF (j == 1) THEN
    IF ((stlist(st_macrotag) >= 1000) .AND. l_meaning_sequence) THEN
      ! ppxm is for climate meaning with l_meaning_sequence=.TRUE.
      packing_profile=ppxm
    ELSE
      packing_profile=0
    END IF
  ELSE IF (j == 2) THEN
    packing_profile=0
  ELSE IF (j < 0) THEN
    NULLIFY(um_file)
    IF (stlist(st_output_type) == st_netcdf) THEN
      um_file => get_file_by_unit(-j, handler="netcdf")
    ELSE
      um_file => get_file_by_unit(-j, handler="portio")
    END IF
    packing_profile = um_file % meta % packing_code
  ELSE
    packing_profile=0
  END IF
  IF (packing_profile == 0) THEN
    packing='       '
  ELSE
    j=ppxref(ppx_packing_acc+packing_profile-1)
    WRITE(packing,'(i7)') j
  END IF

  !
  ! 1.2 Line 2.
  !
  ! Time-processing
  j=stlist(st_proc_no_code)
  tser_list=0
  IF (j == st_replace_code) THEN
    time_proc='   EXTRACT     '
  ELSE IF (j == st_accum_code) THEN
    time_proc=' ACCUMULATION  '
  ELSE IF (j == st_time_mean_code) THEN
    time_proc='  TIME MEAN    '
  ELSE IF (j == st_time_series_code) THEN
    WRITE(time_proc,'(''  TIME SERIES  '')')
    tser_list=stlist(st_series_ptr)
  ELSE IF (j == st_max_code) THEN
    time_proc='MAX OVER PERIOD'
  ELSE IF (j == st_min_code) THEN
    time_proc='MIN OVER PERIOD'
  ELSE IF (j == st_append_traj_code) THEN
    time_proc='  TRAJECTORY   '
  ELSE IF (j == st_variance_code) THEN
    time_proc=' TIME VARIANCE '
  ELSE IF (j == st_time_series_mean) THEN
    time_proc='MEAN TIMESERIES'
  ELSE
    time_proc='  UNKNOWN      '
  END IF

  ! -From-
  IF (stlist(st_freq_code) <  0) THEN
    from='      '
  ELSE
    j=stlist(st_start_time_code)
    WRITE(from,'(i6)') j
  END IF

  ! --To--
  IF (stlist(st_freq_code) <  0) THEN
    TO='      '
  ELSE
    j=stlist(st_end_time_code)
    IF (j == st_infinite_time) THEN
      TO='  ... '
    ELSE
      WRITE(TO,'(i6)') j
    END IF
  END IF

  ! Frequency
  j=stlist(st_freq_code)
  IF (j <  0) THEN
    j=-j
    WRITE(freq,'(''TIME LIST'')')
    time_list=j
  ELSE
    WRITE(freq,'(i9)') j
    time_list=0
  END IF

  ! Period
  IF (stlist(st_freq_code) <  0) THEN
    period='      '
  ELSE
    j=stlist(st_period_code)
    IF (stlist(st_proc_no_code) == st_replace_code) THEN
      period='      '
    ELSE IF (j == st_infinite_time) THEN
      period='  ... '
    ELSE
      WRITE(period,'(i6)') j
    END IF
  END IF

  ! __Source__
  j=stlist(st_input_code)
  IF (j == 0) THEN
    source='PROGNOSTIC'
  ELSE IF (j == 1) THEN
    source='  STWORK  '
  ELSE IF (j <  0) THEN
    j=-j
    WRITE(source,'(''DUMP #'',i4)') j
  ELSE
    source=' UNKNOWN  '
  END IF

  ! ___Destination___
  j=stlist(st_output_code)
  IF (j == 1) THEN
    IF (stlist(st_macrotag) >= 1000) THEN
      dest='MEAN PP VIA DUMP'
    ELSE IF (stlist(st_macrotag) >  0) THEN
      WRITE(dest,'(''DUMP WITH TAG '',i3)') stlist(st_macrotag)
    ELSE
      dest='      DUMP       '
    END IF
  ELSE IF (j == 2) THEN
    dest='   SECONDARY     '
  ELSE IF (j <  0) THEN
    j=-j
    IF (stlist(st_output_type) == st_netcdf) THEN
      WRITE(dest,'(''   NetCDF UNIT'',i3)') j
    ELSE IF (j == 27) THEN
      dest='MEAN PP (DIRECT) '
    ELSE
      WRITE(dest,'(''   PP UNIT'',i3)') j
    END IF
  ELSE
    dest='  UNKNOWN  '
  END IF


  !
  ! 1.3 Line 3.
  !
  ! Spatial-Processing
  j=stlist(st_gridpoint_code)
  IF (j >= extract_base .AND. j <  extract_top) THEN
    spatial='    FULL FIELD    '
  ELSE IF (j >= vert_mean_base .AND. j <  vert_mean_top) THEN
    spatial='  VERTICAL MEAN   '
  ELSE IF (j >= zonal_mean_base .AND. j <  zonal_mean_top) THEN
    spatial='   ZONAL MEAN     '
  ELSE IF (j >= merid_mean_base .AND. j <  merid_mean_top) THEN
    spatial=' MERIDIONAL MEAN  '
  ELSE IF (j >= field_mean_base .AND. j <  field_mean_top) THEN
    spatial=' FIELD MEAN - 2D  '
  ELSE IF (j >= global_mean_base .AND. j <  global_mean_top) THEN
    spatial=' GLOBAL MEAN - 3D '
  ELSE
    spatial='  ** UNKNOWN **   '
  END IF

  ! Levels-domain
  j=stlist(st_output_bottom)
  lev_list=0
  IF (j == st_special_code) THEN
    levels ='STANDARD LEV '
  ELSE IF (j >  0) THEN
    WRITE(levels,'(''LEVELS '',i3,''-'',i3)') j,stlist(st_output_top)
  ELSE IF (j <  0) THEN
    j=-j
    WRITE(levels,'('' LEVELS LIST '')')
    lev_list=j
  END IF

  ! Pseudo-levels
  j=stlist(st_pseudo_out)
  plev_list=0
  IF (j >  0) THEN
    pseudo='PSEUDO-LEV LIST'
    plev_list=j
  ELSE
    pseudo='     NONE      '
  END IF

  ! Horizontal-domain..... suports up 9999 in each direction.

  WRITE(horiz,'(''ROW:'',i4,''-'',i4,'' COL:'',i4,''-'',i4)')           &
       stlist(st_south_code),stlist(st_north_code),                     &
       stlist(st_west_code),stlist(st_east_code)

  ! Weighting
  j=stlist(st_weight_code)
  IF (j == stash_weight_null_code) THEN
    weight='  NONE   '
  ELSE IF (j == stash_weight_area_code) THEN
    weight='  AREA   '
  ELSE IF (j == stash_weight_volume_code) THEN
    weight=' VOLUME  '
  ELSE IF (j == stash_weight_mass_code) THEN
    weight='  MASS   '
  END IF

  ! Masking
  j=MOD(stlist(st_gridpoint_code),block_size)
  IF (j == stash_null_mask_code) THEN
    mask=' NONE  '
  ELSE IF (j == stash_land_mask_code) THEN
    mask=' LAND  '
  ELSE IF (j == stash_sea_mask_code) THEN
    mask='  SEA  '
  ELSE IF (j == stash_nmdi_mask_code) THEN
    mask='NON-MDI '
  ELSE
    mask='UNKNOWN'
  END IF

  !
  ! 1.4 Print the main part of the summary
  !

  IF (mype == 0) THEN
    WRITE(stash_unit,'(a)') diag1
    WRITE(stash_unit,'(1x,i3,1x,a,1x,a,1x,i4,1x,i7,1x,i7,'//                     &
                '1x,a,1x,a,1x,a,1x,i7,1x,i7,1x,a)')                       &
         seqno,NAME,submodel,stlist(st_item_code),stlist(st_sect_no_code),&
         ppxref(ppx_field_code),datatype,grid_type,level_type,            &
         ppxref(ppx_meto8_levelcode),ppxref(ppx_meto8_fieldcode),packing
    WRITE(stash_unit,'(a)') diag2
    WRITE(stash_unit,'(7(1x,a))') &
           time_proc,from,TO,freq,period,source,dest
    WRITE(stash_unit,'(a)') diag3
    WRITE(stash_unit,'(6(1x,a))') spatial,levels,pseudo,horiz,weight,mask
  END IF ! my = 0

  !
  ! 1.5 Print associated time and levels lists if appropriate
  !
  ! 1.5.1 Time list
  !
  IF (time_list /= 0 .AND. mype == 0) THEN
    DO j=1,nsttims
      IF (sttabl(j,time_list) == st_end_of_list) THEN
        ntimes=j-1
        EXIT
      END IF
    END DO

    WRITE(stash_unit,'(a,i3,a)')                                               &
        ' ***** TIME LIST ***** ',ntimes,' times are as follows:-'

    WRITE(stash_unit,'(10(1x,i10))') sttabl(1:ntimes,time_list)

  END IF
  !
  ! 1.5.2 Levels list
  !
  IF (lev_list /= 0 .AND. mype == 0) THEN

    WRITE(stash_unit,'(a,i3,a)')                                               &
        ' ***** LEVELS LIST ***** ',stash_levels(1,lev_list),           &
        ' levels are as follows:-'

    k=1+stash_levels(1,lev_list)

    WRITE(stash_unit,'(10(1x,i10))') stash_levels(2:k,lev_list)

  END IF

  !
  ! 1.5.3 Pseudo-levels list
  !
  IF (plev_list /= 0 .AND. mype == 0) THEN
    WRITE(stash_unit,'(a,i3,a)')                                               &
        ' ***** PSEUDO-LEVELS LIST ***** ',                             &
         stash_pseudo_levels(1,plev_list),                              &
        ' pseudo-levels are as follows:-'

    k=1+stash_pseudo_levels(1,plev_list)
    WRITE(stash_unit,'(10(1x,i10))')  stash_pseudo_levels(2:k,plev_list)

  END IF

  !
  ! 1.5.4 Time series subdomain record list
  !
  IF (tser_list /= 0 .AND. mype == 0) THEN
    i1=stash_series_index(1,tser_list)
    i2=stash_series_index(2,tser_list)
    WRITE(stash_unit,'('' ***** TIME SERIES ***** '',i3,  '//&
        ' '' subdomain records are as follows:-''/ '//&
        ' '' Record      North/South       West/ East     Bottom/  Top'')')&
        i2
    DO j=1,i2
      WRITE(stash_unit,'(3x,i4,1x,3(5x,i5,1x,i5,1x))')                         &
          j,(stash_series(3+k,i1+j-1),k=1,6)
    END DO
  END IF
  !
  ! 1.5.5 Print final ruler line
  !
  IF (mype == 0) THEN
    WRITE(stash_unit,'(a)') headend
  END IF


END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagdesc
END MODULE diagdesc_mod

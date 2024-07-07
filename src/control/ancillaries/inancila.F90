! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!   Subroutine inancila in Module inancila_mod
!
!   Purpose : Takes as input,the code defining the frequency of update
!             of ancillary fields as set by the user interface.
!             Converts them into a list of numbers of timesteps after
!             which each field must be updated, and calculates the
!             frequency with which this list must be interrogated.
!             Where the update interval is in months or years,
!             the check will be carried out each day. The physical
!             files required are also determined by input code.
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Ancillaries
MODULE inancila_mod

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INANCILA_MOD'

CONTAINS

SUBROUTINE inancila(len_fixhd,         &
                    len_inthd,         &
                    len_realhd,        &
                    !Intent (In)
                    len1_levdepc,      &
                    len2_levdepc,      &
                    a_fixhd,           &
                    a_realhd,          &
                    a_levdepc,         &
                    nlookups,          &
                    len1_lookup,       &
                    p_levels,          &
                    tr_levels,         &
                    st_levels,         &
                    sm_levels,         &
                    ntiles,            &
                    ozone_levels,      &
                    tpps_ozone_levels, &
                    si_atmos,          &
                    silen,             &
                    ancillary_steps,   &
                    steps_per_hr,      &
                    icode,cmessage,lcal360)         ! Intent (Out)


USE jules_surface_types_mod
USE check_iostat_mod
USE ancilcta_namelist_mod   ! All
USE ancil_headers_mod, ONLY: fixhd_ancila, inthd_ancila, realhd_ancila,       & 
                             lookup_ancila, lookup_start_ancila
USE ancil_mod, ONLY: max_items, stash_num_max
USE clmchfcg_scenario_mod, ONLY: nsulpat
USE um_parcore, ONLY: mype
USE setup_namelist, ONLY: setup_nml_type
USE near_equal_real_mod, ONLY: near_equal_real
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE io, ONLY: setpos
USE io_constants, ONLY: ioLocalReplicate
USE model_file, ONLY: model_file_open
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE UM_ParVars
USE lookup_addresses  ! sets model_code=45, item_code=42
USE submodel_mod,      ONLY: atmos_im
USE ozone_inputs_mod,  ONLY: zon_av_ozone, zon_av_tpps_ozone

USE filenamelength_mod, ONLY: filenamelength
USE missing_data_mod, ONLY: rmdi, imdi
USE ancil_mod, ONLY: num_ancil_requests, ancil_requests, ancil_files, &
     num_ancil_files
USE um_stashcode_mod ! Stash codes
USE cancila_mod, ONLY: nlookup, lookup_step, levels, d1_anciladd, steps, &
                       nancil_fields
USE file_manager, ONLY: assign_file_unit, release_file_unit
! ITEMS name list
USE items_nml_mod, ONLY:                                                    &
    stash_req, period, interval, source, domain, user_Prog_Ancil_stash_req, &
    User_Prog_RConst, update_anc, ancilfilename, items, print_nlist_items,  &
    netcdf_varname, period_years, period_months, period_days, period_hours
USE populate_ancil_requests_mod, ONLY: populate_ancil_requests

USE errormessagelength_mod, ONLY: errormessagelength
USE ancil_check_mod, ONLY: ancil_check_horizontal_grid, &
    ancil_check_grid_stagger
USE decomp_params, ONLY: gl_num_rows, gl_row_length
USE field_types, ONLY: fld_type_p, fld_type_u, fld_type_v, fld_type_r
IMPLICIT NONE

LOGICAL :: lcal360  ! Logical switch for 360-day calendar
LOGICAL :: found    ! Used in DO loops

INTEGER :: len_fixhd  ! Length of fixed header   ) in         
INTEGER :: len_inthd  ! Length of Integer header ) ancillary  
INTEGER :: len_realhd ! Length of Real header    ) files              
INTEGER :: len1_levdepc  ! ) First and second dimension of model
INTEGER :: len2_levdepc  ! ) level dependent constants array

INTEGER :: ancillary_steps
INTEGER :: steps_per_hr

INTEGER :: ndatasets    ! No of physical files
INTEGER :: nlookups     ! No of lookups required(set by User I.)
INTEGER :: iounit       ! I/O unit number
INTEGER :: len1_lookup  ! Length of PP header
INTEGER :: p_levels     ! No. of pressure levels
INTEGER :: tr_levels    ! No. of tracer levels
INTEGER :: st_levels    ! No. of soil temperature levels
INTEGER :: sm_levels    ! No. of soil moisture levels
INTEGER :: ntiles       ! No. of surface tiles.
INTEGER :: ozone_levels ! No. of ozone levels
INTEGER :: tpps_ozone_levels ! No of ozone levs in TppsOzon dataset
INTEGER :: silen             ! Length for SI_ATMOS array
INTEGER :: si_atmos(silen)   ! ) STASHin addresses of atmos

LOGICAL :: l_vert_mismatch    ! T : Vertical levels mismatch

INTEGER :: a_fixhd(len_fixhd),                                                &
                            ! Fixed header for Dump
           icode,                                                             &
                        ! Return code =0 Normal Exit  >0 Error
           errorstatus

REAL :: a_realhd(len_realhd),                                                 &
        a_levdepc(len1_levdepc,len2_levdepc),                                 &
        levdepc( (p_levels+1)*4 ) ! Space to hold level dependent
                                  !  constants from Anc File

REAL :: coldepc(glsize(gl_row_length,fld_type_p)+1)
REAL :: rowdepc(glsize(gl_num_rows,fld_type_p)+1)

CHARACTER(LEN=errormessagelength) :: cmessage       ! Out error message if I>0
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'INANCILA'

! Local Variables

CHARACTER (LEN=filenamelength) :: ancil_path
INTEGER :: i
INTEGER :: item
INTEGER :: j
INTEGER :: j1
INTEGER :: j2
INTEGER :: k
INTEGER :: i_stash       ! Used to loop over stashcodes present in ancil file
INTEGER :: ds_index      ! index for data set/file in lookup_start_ancila array
INTEGER :: len_io
INTEGER :: lookups
INTEGER :: nftin         ! Current FTN number for ancillary data
INTEGER :: start_block
INTEGER :: nrec_a        ! No of atmos records
INTEGER :: stash_addr    ! Stash address
INTEGER :: dummy

INTEGER :: p_row_length ! Row length for pressure-type variables
INTEGER :: u_row_length ! Row length for U wind-type variables
INTEGER :: v_row_length ! Row length for V wind-type variables
INTEGER :: river_row_length ! Row length for river routing-type variables

INTEGER :: p_rows       ! No. of rows for pressure-type variables
INTEGER :: u_rows       ! No. of rows for U wind-type variables
INTEGER :: v_rows       ! No. of rows for V wind-type variables
INTEGER :: river_rows   ! No. of rows for river routing-type variables

INTEGER :: sea_ice_frac_index
INTEGER :: surface_temp_index

INTEGER :: ids                       ! data set (i.e. ancil file) counter
INTEGER :: stashcode                ! stash item number from items NL
INTEGER :: shared_unit

INTEGER, PARAMETER :: filenameprovided = 1

DATA dummy /1/

CHARACTER(LEN=8) :: cperiod      ! PERIOD in characters.


!     SoilDepths : position of soil depths in level dependent constants
INTEGER, PARAMETER :: soildepths = 4

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!  Internal Structure

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode=0
cmessage=' '
iounit=0

!   1.  Initialisation for atmosphere model

! Take grid sizing from decomp database and pass into local variable names 
! to aid readability
p_row_length = glsize(gl_row_length,fld_type_p)
u_row_length = glsize(gl_row_length,fld_type_u)
v_row_length = glsize(gl_row_length,fld_type_v)
river_row_length = glsize(gl_row_length,fld_type_r)
p_rows = glsize(gl_num_rows,fld_type_p)
u_rows = glsize(gl_num_rows,fld_type_u)
v_rows = glsize(gl_num_rows,fld_type_v)
river_rows = glsize(gl_num_rows,fld_type_r)

! Assign a file unit to each ancil file
DO i = 1, num_ancil_files
  CALL assign_file_unit(ancil_files(i)%filename, ancil_files(i)%unit_num, &
                        handler="portio")
  IF (printstatus >= PrStatus_Normal .AND. mype == 0) THEN
    WRITE(umMessage, '(A,I0)') "Ancil file num     : ", i 
    CALL umPrint(umMessage, src='inancila')
    WRITE(umMessage, '(A,I0)') "Ancil file unit num: ", ancil_files(i)%unit_num
    CALL umPrint(umMessage, src='inancila')
    WRITE(umMessage, '(A,A)') "Ancil filename      : ", ancil_files(i)%filename
    CALL umPrint(umMessage, src='inancila')
    DO j = 1, ancil_files(i)%num_stash
      WRITE(umMessage, '(A,I0)') "Stash req = ",  ancil_files(i)%stashcodes(j)
      CALL umPrint(umMessage, src='inancila')
    END DO
    WRITE(umMessage,'(A)') "---------------------------------------------" // &
         newline // "------------------"
    CALL umPrint(umMessage, src='inancila')
  END IF
END DO

CALL umPrint('',src='inancila')
WRITE(umMessage,'(I0,A)') num_ancil_requests,' Atmos Ancillaries to be updated.'
CALL umPrint(umMessage,src='inancila')
DO i=1,num_ancil_requests
  IF (ancil_requests(i)%period >  0) THEN
    IF (ancil_requests(i)%period == period_years) cperiod=' Years'
    IF (ancil_requests(i)%period == period_months) cperiod=' Months'
    IF (ancil_requests(i)%period == period_days) cperiod=' Days'
    IF (ancil_requests(i)%period == period_hours) cperiod=' Hours'
    WRITE(umMessage,'(A,I0,A,I0,A,I0,A)') 'Anc field ',i,' Stash code ', &
         ancil_requests(i)%stashcode,                 &
        ' Interval ', ancil_requests(i)%interval, cperiod
    CALL umPrint(umMessage,src='inancila')
  END IF
END DO
CALL umPrint('',src='inancila')

! Check that ancillary field has valid address (>1) before attempting
! and ereport if not
DO i=1,num_ancil_requests
  IF (si_atmos(ancil_requests(i)%stashcode)  <=  1) THEN
    icode = 67
    WRITE(cmessage,'(A,I0,A)')' INANCILA: update requested for STASHcode ', &
          ancil_requests(i)%stashcode, newline //  &
          'but prognostic address not set.' // newline // &
          'As prognostic is not present in UM run, ' //   &
          'please turn off ancil updating' //  newline // &
          'via the items namelist.' 
    CALL ereport( routinename, icode, cmessage )
  END IF
END DO

!   1.1 Set number of steps after which each ancillary field is updated

DO i=1,num_ancil_requests
  steps(i)=0
  IF (ancil_requests(i)%period == period_hours) THEN
    steps(i)=ancil_requests(i)%interval*steps_per_hr
  END IF
  IF (ancil_requests(i)%period == period_days) THEN
    steps(i)=ancil_requests(i)%interval*24*steps_per_hr
  END IF

  IF (lcal360) THEN
    IF (ancil_requests(i)%period == period_months) THEN
      steps(i)=ancil_requests(i)%interval*30*24*steps_per_hr
    END IF
    IF (ancil_requests(i)%period == period_years) THEN
      steps(i)=ancil_requests(i)%interval*360*24*steps_per_hr
    END IF
  ELSE
    ! Gregorian calender:
    ! If update interval is months or years, test each day. Further testing
    ! done in REPLANCA.

    IF (ancil_requests(i)%period == period_years  .OR.  & 
        ancil_requests(i)%period == period_months) THEN
      steps(i)=24*steps_per_hr
    END IF
  END IF

END DO

!   1.2 Set ANCILLARY_STEPS to lowest common denominator of
!       frequencies for active fields

ancillary_steps = steps(1)

DO i=2,num_ancil_requests
  IF (steps(i) <  ancillary_steps                                              &
      .AND. steps(i) >  0) THEN
    IF (MOD(ancillary_steps,steps(i)) == 0) THEN
      ancillary_steps=steps(i)
    ELSE
      j1=steps(i)-1
      DO j=j1,1,-1
        IF ((MOD(ancillary_steps,j) == 0) .AND.                                &
            (MOD(steps(i),j)        == 0)) THEN
          ancillary_steps = j
          EXIT
        END IF
      END DO
    END IF
  END IF
END DO

!  Sea surface temperature must be updated when sea ice is updated

sea_ice_frac_index = 0
surface_temp_index = 0

DO i = 1, num_ancil_requests
  IF (ancil_requests(i)%stashcode == stashcode_icefrac) THEN
    sea_ice_frac_index = i
  ELSE IF (ancil_requests(i)%stashcode == stashcode_surftemp) THEN
    surface_temp_index = i
  END IF
  IF (sea_ice_frac_index > 0 .AND. surface_temp_index > 0) THEN
    IF (steps(sea_ice_frac_index) >  0 .AND. steps(surface_temp_index) <= 0)THEN
      steps(surface_temp_index)=1
      EXIT
    END IF
  END IF
END DO

!  1.3 Set number of headers for each ancillary field

DO i=1,num_ancil_requests
  stashcode = ancil_requests(i)%stashcode
  levels(i)=1
  !   Multilayer hydrology
  IF (stashcode == stashcode_soil_moist)levels(i)=sm_levels

  !   Multilayer murk concentration and source
  IF (stashcode == stashcode_total_aero_emiss .OR.                   &
     stashcode == stashcode_total_aero) levels(i)  =p_levels

  !   Multilayer user ancillaries
  IF (stashcode >= stashcode_user_anc_mult1 .AND.                     &
     stashcode <= stashcode_user_anc_mult20) levels(i)=p_levels

  !   Multi-level ancillaries for sulphur cycle
  IF (stashcode == stashcode_3d_nat_so2_em  ) levels(i) = p_levels
  IF (stashcode == stashcode_3d_oh_conc     ) levels(i) = p_levels
  IF (stashcode == stashcode_3d_ho2_conc    ) levels(i) = p_levels
  IF (stashcode == stashcode_3dh2o2_mixrat  ) levels(i) = p_levels
  IF (stashcode == stashcode_3d_ozone_mixrat) levels(i) = p_levels

  IF (stashcode == stashcode_frac_surf_type)  levels(i) = ntype
  IF (stashcode == stashcode_lai           )  levels(i) = npft
  IF (stashcode == stashcode_canopy_height )  levels(i) = npft

  !   Multi-level ancillaries aerosol climatology
  IF (stashcode >= stashcode_clim_biogenic_aero .AND.                 &
      stashcode <= stashcode_clim_delta_aero) levels(i)=p_levels

  IF (stashcode >= stashcode_ozone_tracer .AND.                       &
      stashcode <= stashcode_o3_colo3)        levels(i)=p_levels

  IF (stashcode == stashcode_ozone)           levels(i)=ozone_levels

  IF (stashcode == stashcode_soil_temp)       levels(i)=st_levels

END DO


!  1.4 Read headers

lookups=0
ids = 0

DO i=1, num_ancil_files

!  stashcode = ancil_requests(i)%stashcode
  nftin = ancil_files(i)%unit_num
  ancil_path = ancil_files(i)%filename

  CALL umprint('',src='inancila')
  WRITE(ummessage,*) &
           'Ancillary file ', i ,' unit ', nftin, ' : ', TRIM(ancil_path)
  CALL umprint(ummessage,src='inancila')

  !   1.4.1 Read fixed length header record

  CALL model_file_open (nftin, TRIM(ancil_path),                             &
                  LEN_TRIM(ancil_path), 0, filenameprovided, icode)
  IF (icode /= 0) THEN
    cmessage='INANCLA: Error opening file'
    WRITE(umMessage,'(A,I0,A,I0)') 'INANCILA: Error opening file on unit ', &
         nftin,' filename: ',TRIM(ancil_path)
    CALL umPrint(umMessage,src='inancila')
    GO TO 9999   !  Return
  END IF

  CALL setpos(nftin,0,icode)

  !   Read in fixed header to get array dimensions
  ! DEPENDS ON: read_flh
  CALL read_flh(nftin,fixhd_ancila(1,i),len_fixhd,icode,cmessage)
  IF (icode >  0) THEN
    WRITE(umMessage,'(A,A)') ' Error in reading fixed header for file ',&
        TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='inancila')
    GO TO 9999   !  Return
  END IF

  CALL ancil_check_grid_stagger(model_grid_stagger = a_fixhd(9) ,         &
                                ancil_grid_stagger = fixhd_ancila(9,i),   &
                                ancil_file = ancil_files(i))

  !       Check for negative dimensions
  IF (fixhd_ancila(101,i) <= 0) fixhd_ancila(101,i)=1
  IF (fixhd_ancila(106,i) <= 0) fixhd_ancila(106,i)=1
  IF (fixhd_ancila(111,i) <= 0) fixhd_ancila(111,i)=1
  IF (fixhd_ancila(112,i) <= 0) fixhd_ancila(112,i)=1
  IF (fixhd_ancila(151,i) <= 0) fixhd_ancila(151,i)=1
  IF (fixhd_ancila(152,i) <= 0) fixhd_ancila(152,i)=1
  IF (fixhd_ancila(161,i) <= 0) fixhd_ancila(161,i)=1

  ! Start position of fields for file i
  lookup_start_ancila(i,1)=lookups+1
  lookup_start_ancila(i,2)=nftin

  IF (lookups+fixhd_ancila(152,i) >  nlookups) THEN
    WRITE(umMessage,'(A,A)')                                         &
                'No room in LOOKUP table for Ancillary File ',       &
                TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='inancila')
    cmessage='INANCILA: Insufficient space for LOOKUP headers'
    icode=14
    GO TO 9999   !  Return
  END IF

  CALL setpos(nftin,0,icode)
  IF (icode /= 0) THEN
    WRITE (Cmessage, '(3A)')                               &
      'Problem in SETPOS for Ancillary file ',             &
    TRIM(ancil_files(i)%filename),' before READHEAD.'
    CALL Ereport ( RoutineName, icode, Cmessage )
  END IF! Check icode

  ! DEPENDS ON: readhead
  CALL readhead(nftin,                                         &
              fixhd_ancila(1,i),len_fixhd,                     &
              inthd_ancila(1,i),fixhd_ancila(101,i),           &
              realhd_ancila(1,i),fixhd_ancila(106,i),          &
              levdepc,fixhd_ancila(111,i),fixhd_ancila(112,i), &
              rowdepc,fixhd_ancila(116,i),fixhd_ancila(117,i), &
              coldepc,fixhd_ancila(121,i),fixhd_ancila(122,i), &
              dummy,dummy,dummy,                               &
              dummy,dummy,                                     &
              dummy,dummy,                                     &
              dummy,dummy,                                     &
              dummy,dummy,                                     &
              dummy,dummy,                                     &
              lookup_ancila(1,lookups+1),fixhd_ancila(151,i),  &
              fixhd_ancila(152,i), fixhd_ancila(161,i),        &
              start_block,icode,cmessage)

  IF (icode >  0) THEN
    WRITE(umMessage,'(A,A)') 'ERROR in READHEAD for Ancillary File ',  &
         TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='inancila')
    GO TO 9999   !   Return
  END IF

  !   Check calendar indicator is correct if this is set in the ancillary
  IF (fixhd_ancila(8,i) /=  imdi) THEN
    IF (( lcal360 .AND. fixhd_ancila(8,i) /= 2) .OR.                       &
         (.NOT. lcal360 .AND. fixhd_ancila(8,i) /= 1) ) THEN
      icode=100+i
      cmessage='INANCILA : Wrong calendar set in Ancillary File'
      CALL umPrint( ' ******** Error in INANCILA ********',                  &
          src='inancila')
      WRITE(umMessage,'(a,i4)') &
         ' Wrong calendar setting in Ancillary File',i
      CALL umPrint(umMessage,src='inancila')
      IF (lcal360) THEN
        CALL umPrint( ' Model run is set up for 360 day calendar.')
        CALL umPrint( ' Ancillary File is for 365 day calendar.')
      ELSE
        CALL umPrint( ' Model run is set up for 365 day calendar.')
        CALL umPrint( ' Ancillary File is for 360 day calendar.')
      END IF
      CALL umPrint( ' Rerun with correct ancillary file.')
      GO TO 9999   !  Return
    END IF
  ELSE
    WRITE(umMessage,'(a,i4)') &
                'Unspecified calendar type in ancillary file ',i
    CALL umPrint(umMessage,src='inancila')
  END IF

  !  1.4.2 Check horizontal grid
  
  ! Value of integer `lookups` keeps track of the index (in the master 
  ! lookup array)of the last lookup of the previous ancil file to be read.  
  ! Pass in slice of the master lookup array which corresponds to the lookups 
  ! for the current ancil file
  CALL ancil_check_horizontal_grid(lookup_ancila(:,                        &
                                   lookups+1:lookups+fixhd_ancila(152,i)), &
                                   ancil_files(i),  len1_lookup,           &
                                   p_rows,           p_row_length,         &
                                   u_rows,           u_row_length,         &
                                   v_rows,           v_row_length,         &
                                   river_rows,       river_row_length)


  !  1.4.3 Check real constants

  IF (fixhd_ancila(105,i) >  0) THEN

    ! Only perform this check if ancillary and model is on C grid 
    ! (with P at poles)
    IF (fixhd_ancila(9,i) == 3 .AND. a_fixhd(9) == 3 ) THEN
      DO j=1,6
        IF (realhd_ancila(j,i) >  (a_realhd(j)+0.1) .OR.                   &
           realhd_ancila(j,i) <  (a_realhd(j)-0.1)) THEN
          ! Need to check all stashcodes contained in the file
          ! so loop over all stash in file
          DO i_stash = 1, ancil_files(i)%num_stash
            stashcode = ancil_files(i)%stashcodes(i_stash)
            IF (stashcode /= stashcode_ozone .OR. (j /= 1 .AND. j /= 4)) THEN
              WRITE(umMessage,*)(realhd_ancila(k,i),k=1,6),                  &
                                (a_realhd(k),k=1,6)
              CALL umPrint(umMessage,src='inancila')
              icode=8
              cmessage='INANCILA: REAL header Error.'
              GO TO 9999   !  Return
            END IF
          END DO ! End loop over stashcodes in file
        END IF
      END DO
    END IF

  END IF

  !  1.4.4 Check level dependent constants if required
  !        Not retained in model after initial check

  IF (fixhd_ancila(110,i) >  0) THEN
    ! Need to check all stashcodes contained in the file
    ! so loop over all stash in file
    DO i_stash = 1, ancil_files(i)%num_stash
      stashcode = ancil_files(i)%stashcodes(i_stash)
      IF (stashcode == stashcode_ozone             .OR.                      &
          stashcode == stashcode_total_aero        .OR.                      &
          stashcode == stashcode_total_aero_emiss  .OR.                      &
         (stashcode >= stashcode_user_anc_mult1 .AND.                        &
          stashcode <= stashcode_user_anc_mult20 ) .OR.                      &
         (stashcode >= stashcode_ozone_tracer .AND.                          &
          stashcode <= stashcode_o3_colo3)               ) THEN

        ! Check that ancillary file is set up for correct vertical levels

        IF (fixhd_ancila(111,i)-1 /= p_levels) THEN
          icode=110
          WRITE (cmessage,*) ' Ancillary File set up for wrong',               &
              ' no of model levels. Anc ',fixhd_ancila(111,i)-1,             &
              ' Model ',p_levels
          CALL ereport ( routinename, icode, cmessage )
        END IF

        l_vert_mismatch = .FALSE.
        
        ! Check eta_theta and eta_rho
        
        DO j=1,p_levels+1
          IF (near_equal_real( levdepc(j), a_levdepc(j,1) )) THEN
            l_vert_mismatch = .TRUE.
            EXIT
          END IF
        END DO
        
        DO j=1,p_levels
          IF (near_equal_real( levdepc(fixhd_ancila(111,i)+j),                 &
                               a_levdepc(j,2) )) THEN
            l_vert_mismatch = .TRUE.
            EXIT
          END IF
        END DO
        
        ! Abort if there is a mis-match
        
        IF (l_vert_mismatch) THEN
          WRITE(umMessage,*) 'Mismatch in vertical levels between model ',     &
               'and Ancillary File.'
          CALL umPrint(umMessage,src='inancila')
          WRITE(umMessage,*) 'Anc File : ', TRIM(ancil_path)
          CALL umPrint(umMessage,src='inancila')
          WRITE(umMessage,*) 'Eta_Theta - Model'
          CALL umPrint(umMessage,src='inancila')
          WRITE(umMessage,'(5F10.7)') (a_levdepc(k,1),k=1,p_levels+1)
          CALL umPrint(umMessage,src='inancila')
          WRITE(umMessage,*) 'Eta_Theta - Anc File'
          CALL umPrint(umMessage,src='inancila')
          WRITE(umMessage,'(5F10.7)') (levdepc(k),k=1,p_levels+1)
          CALL umPrint(umMessage,src='inancila')
          WRITE(umMessage,*) 'Eta_Rho   - Model'
          CALL umPrint(umMessage,src='inancila')
          WRITE(umMessage,'(5F10.7)') (a_levdepc(k,2),k=1,p_levels)
          CALL umPrint(umMessage,src='inancila')
          WRITE(umMessage,*) 'Eta_Rho   - Anc File'
          CALL umPrint(umMessage,src='inancila')
          WRITE(umMessage,'(5F10.7)') (levdepc(p_levels+1+k),k=1,p_levels)
          CALL umPrint(umMessage,src='inancila')
          icode=11
          WRITE (cmessage,*) 'Mismatch in LEVDEPC ',                           &
            'between model and Ancillary File.'
          
          CALL ereport ( routinename, icode, cmessage )
        END IF

      ELSE IF (stashcode == stashcode_soil_moist) THEN
        IF (printstatus >= prstatus_diag .AND. mype == 0 ) THEN
          CALL umPrint('',src='inancila')
          WRITE(umMessage,*) 'SoilDepths = ',soildepths
          CALL umPrint(umMessage,src='inancila')
          WRITE(umMessage,*) 'SM_Levels  = ',sm_levels
          CALL umPrint(umMessage,src='inancila')
          DO j=1,sm_levels
            WRITE(umMessage,*) 'model ',a_levdepc(j,soildepths),              &
               ' anc ',levdepc(fixhd_ancila(111,i)*3+j)
            CALL umPrint(umMessage,src='inancila')
          END DO
        END IF
        
        ! Check Soil moisture levels
        
        DO j=1,sm_levels
          IF (near_equal_real(levdepc(fixhd_ancila(111,i)*3+j),               &
                  a_levdepc(j,soildepths))) THEN
            icode=12
            cmessage='INANCILA: error in LEVDEPC.'
            GO TO 9999
          END IF
        END DO
        
      ELSE IF (stashcode == stashcode_soil_temp) THEN ! Deep soil temperature
        
        IF (printstatus >= prstatus_diag .AND. mype == 0) THEN
          CALL umPrint('',src='inancila')
          WRITE(umMessage,*) 'SoilDepths = ',soildepths
          CALL umPrint(umMessage,src='inancila')
          WRITE(umMessage,*) 'st_levels  = ',st_levels
          CALL umPrint(umMessage,src='inancila')
          DO j=1,st_levels
            WRITE(umMessage,*) 'model ',a_levdepc(j,soildepths),              &
               ' anc ',levdepc(fixhd_ancila(111,i)*3+j)
            CALL umPrint(umMessage,src='inancila')
          END DO
        END IF

        ! Check Deep Soil levels

        DO j=1,st_levels
          IF (near_equal_real( levdepc(fixhd_ancila(111,i)*3+j),              & 
                    a_levdepc(j,soildepths))) THEN
            icode=13
            cmessage='INANCILA: error in LEVDEPC.'
            GO TO 9999
          END IF
        END DO
        
      ELSE IF (stashcode == stashcode_ocnsrf_chlorophyll) THEN
        
        DO j=1,tr_levels
          DO j1=1,4
            IF (near_equal_real(levdepc(j+(j1-1)*fixhd_ancila(111,i)),        &
                     a_levdepc(j,j1))) THEN
              WRITE(umMessage,*)'Error in level dependent constants:Level=',j
              CALL umPrint(umMessage,src='inancila')
              WRITE(umMessage,*)'Position=',j1
              CALL umPrint(umMessage,src='inancila')
              WRITE(umMessage,*)'Value in model =',a_levdepc(j,j1)
              CALL umPrint(umMessage,src='inancila')
              WRITE(umMessage,*)'Value in ancillary data =',levdepc(j+        &
                 (j1-1)*fixhd_ancila(111,i))
              CALL umPrint(umMessage,src='inancila')
              icode=16
              cmessage='INANCILA: error in LEVDEPC.'
              GO TO 9999
            END IF
          END DO
        END DO
       
      END IF  ! stashcode
    END DO ! End loop over stashcodes in file

  END IF  ! Fixhd_ancila(110,i) > 0, level dep constants

  !  1.4.5 Count lookup headers

  IF (fixhd_ancila(150,i) >  0) THEN

    nrec_a = 0
    DO j = 1,fixhd_ancila(152,i)
      IF (lookup_ancila(model_code,lookups+j)  ==  0 .OR.                   &
         lookup_ancila(model_code,lookups+j)  ==  imdi) THEN
        lookup_ancila(model_code,lookups+j) = atmos_im
        nrec_a = nrec_a+1
      END IF
    END DO
    IF (nrec_a >  0) THEN
      CALL umPrint('',src='inancila')
      WRITE(umMessage,*) ' INANCA1A : submodel_id in ',nrec_a,              &
         ' records set to atmos_im in ancillary file ',i
      CALL umPrint(umMessage,src='inancila')
    END IF

  END IF

  lookups=lookups+fixhd_ancila(152,i) ! fixhd_ancila(152,i) = no of 
                                        ! lookups in file i

END DO ! 1, num_ancil_requests


!  1.5 Set positions in main data blocks

DO i=1,num_ancil_requests
  IF (ancil_requests(i)%stashcode  >   0) THEN
    d1_anciladd(i)=si_atmos(ancil_requests(i)%stashcode)
  ELSE
    d1_anciladd(i)=0
  END IF
END DO

!  1.51 If a request is made to update a field, ensure that space for
!      that field has been allocted in D1.

DO i=1,num_ancil_requests
  IF ((ancil_requests(i)%stashcode >  0) .AND. (d1_anciladd(i) <= 1)) THEN
    WRITE(umMessage,'(A,I0)') ' No D1 address set for ancillary item ', &
      ancil_requests(i)%stashcode
    CALL umPrint(umMessage,src='inancila')
    icode=30
    WRITE(cmessage,*) 'INANCILA: updating for ancillary field is requested',   &
        ' but no space has been allocated in D1'
    GO TO 9999
  END IF
END DO

!  1.6 Set positions of data

! In the anc file, the data fields are ordered
!   month, stash items for each month, levels for each stash item
! nlookup(i) specifies the position of the first lookup for ancil field i,
!   relative to lookup_start_ancila for the data set containing field i
! levels(i) is the number of levels/lookups for ancil field i
! lookup_step(i) specifies the interval between occurrences of a given
!   lookup for each stash item

DO i = 1, num_ancil_requests

  nlookup(i) = 0
  lookup_step(i)= 0
  stashcode = ancil_requests(i)%stashcode

  IF (stashcode > 0) THEN

    ! Get fortran unit for this stash item
    nftin    = ancil_files(ancil_requests(i)%ancil_file_num)%unit_num
    ds_index = ancil_requests(i)%ancil_file_num

    ! Get lookup_start_ancila position (in lookup table) for this fortran unit
    IF (lookup_start_ancila(ds_index,1) >  0) THEN
      ! Find first occurrence of stash item in lookup table
      DO j = lookup_start_ancila(ds_index,1), lookups
        IF (lookup_ancila(item_code,j) == stashcode) THEN
          nlookup(i) = j - lookup_start_ancila(ds_index,1) + 1
          lookup_step(i)=0
          EXIT
        END IF
      END DO

      ! Find second occurrence of of stash item, level, in lookup table,
      !   to set LOOKUP_STEP
      IF (j <  lookups) THEN
        DO j1=j+levels(i),lookups
          IF (lookup_ancila(item_code,j1) == stashcode) THEN
            lookup_step(i) = j1 - nlookup(i) -                                &
                             lookup_start_ancila(ds_index,1) + 1
            EXIT
          END IF
        END DO
      END IF

      ! Check that ancillary file is consistent since we assume:
      !   for month in month_list:
      !     for STASH in STASH_list:
      !       for level in level_psuedolevel_list:

      ! First check if its periodic and whether it has
      !   multiple fields to step through.
      IF (fixhd_ancila(10,ds_index) == 2 .AND. lookup_step(i) > 0) THEN
        ! Loop over the fields which should only differ by time
        ! Note, fixhd_ancila(152,ds_index) = no of lookups in ancil 
        ! file "ds_index"
        DO j = lookup_start_ancila(ds_index,1) + nlookup(i) - 1,              &
               lookup_start_ancila(ds_index,1) +                              &
               fixhd_ancila(152,ds_index) - 1,                                &
               lookup_step(i)
          ! If the STASH is different we have a problem.
          IF (lookup_ancila(item_code,j) /= stashcode ) THEN
            icode = 17
            CALL ereport ( routinename, icode,                                &
                           "Bad ancil file: "//TRIM(ancil_files(i)%filename))
            EXIT
          END IF
        END DO
      END IF

    END IF

  END IF

END DO

! Set levels=2 for ice fraction and snow depth
!   to indicate presence of fractional time fields
DO i=1,num_ancil_requests
  stashcode = ancil_requests(i)%stashcode
  IF (stashcode == stashcode_icefrac    ) levels(i) = 2
  IF (stashcode == stashcode_snow_amount) levels(i) = 2
END DO

9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE inancila

END MODULE inancila_mod

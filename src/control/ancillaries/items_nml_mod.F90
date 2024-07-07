! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Reads the ITEMS namelists
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Ancillaries
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE items_nml_mod

USE filenamelength_mod,  ONLY: filenamelength

USE missing_data_mod,    ONLY: &
    imdi,                      &
    rmdi

! netcdf_file and external_dump to be USED once they are moved to the UM

IMPLICIT NONE

INTEGER, PARAMETER   :: netcdf_varname_len = 256 ! Maximum length of varname
CHARACTER (LEN=5), PARAMETER :: netcdf_unset='unset'

CHARACTER (LEN=*), PRIVATE, PARAMETER :: ModuleName = "ITEMS_NML_MOD"

! Declares the items name list, used by Rcf and UM ancillary updating code
INTEGER, ALLOCATABLE :: stash_req(:) ! Stash codes of ancillary fields
INTEGER              :: period   = imdi  ! Period of Updating Interval (Y/M/D/H)
INTEGER              :: interval = imdi  ! Updating Interval
INTEGER              :: source   = imdi  ! Source of data
INTEGER              :: domain   = imdi  ! Areal coverage of data
INTEGER              :: user_Prog_Ancil_stash_req = imdi
REAL                 :: User_Prog_RConst = rmdi
LOGICAL              :: update_anc = .FALSE. ! true if field is to be updated
CHARACTER (LEN=filenamelength)                  :: ancilfilename
CHARACTER (LEN=netcdf_varname_len), ALLOCATABLE :: netcdf_varname(:)

! Parameters for a source of a field
INTEGER, PARAMETER   :: Input_Dump            = 1
INTEGER, PARAMETER   :: Ancillary_File        = 2
INTEGER, PARAMETER   :: Set_To_Zero           = 3
INTEGER, PARAMETER   :: Set_To_MDI            = 4
INTEGER, PARAMETER   :: Tracer_File           = 5
INTEGER, PARAMETER   :: Set_To_Const          = 6
INTEGER, PARAMETER   :: External_Dump         = 7
INTEGER, PARAMETER   :: Field_Calcs           = 8
INTEGER, PARAMETER   :: Field_Dependent_Calcs = 9
INTEGER, PARAMETER   :: NetCDF_File           = 10

! Domain values
INTEGER, PARAMETER   :: Whole_Grid     = 1
INTEGER, PARAMETER   :: Sub_Grid       = 2

! Parameters for ancil updating periods
INTEGER, PARAMETER :: period_years  = 1
INTEGER, PARAMETER :: period_months = 2
INTEGER, PARAMETER :: period_days   = 3
INTEGER, PARAMETER :: period_hours  = 4

NAMELIST /items/                                                              &
   stash_req, period, interval, source, domain, User_Prog_Ancil_stash_req,    &
   User_Prog_RConst, ancilfilename, netcdf_varname, update_anc

CONTAINS

SUBROUTINE check_nlist_items()

! Description:
!   Subroutine to apply logic checks and set control variables based on the
!   options selected in the items namelist.

USE chk_opts_mod,           ONLY: &
    chk_var,                      &
    def_src

USE parkind1,               ONLY: &
    jpim,                         &
    jprb

USE yomhook,                ONLY: &
    lhook,                        &
    dr_hook

IMPLICIT NONE

INTEGER                        :: i        ! counter
INTEGER                        :: nitems   ! number of items in stash_req

CHARACTER (LEN=*), PARAMETER   :: RoutineName = 'CHECK_NLIST_ITEMS'

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)                :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

def_src = RoutineName
nitems  = SIZE(stash_req)

IF ( source == netcdf_file ) THEN
  ! domain value - Whole Grid is only allowed option for source == netcdf_file
  CALL chk_var( domain, 'domain', '[1]' )
ELSE
  ! domain values from rcf_data_source_mod 1=Whole Grid, 2=Sub Grid
  CALL chk_var( domain, 'domain', '[1, 2]' )
END IF

! source values from rcf_data_source_mod ( -90 and -99 used internally )
CALL chk_var( source, 'source', '[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]' )

DO i=1, nitems
  ! stash_req value
  CALL chk_var( stash_req(i), 'stash_req', '[0:99999]' )
END DO

IF ( update_anc ) THEN
  ! interval value
  CALL chk_var( interval, 'interval', '[>=1]' )
  
  ! period value
  CALL chk_var( period, 'period', '[1, 2, 3, 4]' )
END IF

IF ( source == external_dump ) THEN
  ! user_prog_ancil_stash_req value
  CALL chk_var( user_prog_ancil_stash_req, 'user_prog_ancil_stash_req',       &
                                           '[0:99999]' )
END IF

def_src = ''

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE check_nlist_items

SUBROUTINE print_nlist_items()

USE umPrintMgr, ONLY: umPrint

IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
CHARACTER(LEN=80) :: temp_format
INTEGER :: nitems

nitems = SIZE(stash_req)

CALL umPrint('Contents of namelist items', &
    src='items_nml_mod')
WRITE(temp_format, '(A,I0,A)') "(A,",nitems,"I6)"
WRITE(lineBuffer,temp_format)' Stash_req = ',stash_req
CALL umPrint(lineBuffer,src='items_nml_mod')

WRITE(lineBuffer,'(A,I0)')' Period = ',period
CALL umPrint(lineBuffer,src='items_nml_mod')
WRITE(lineBuffer,'(A,I0)')' Interval = ',interval
CALL umPrint(lineBuffer,src='items_nml_mod')

WRITE(lineBuffer,'(A,I0)')' Source = ',Source
CALL umPrint(lineBuffer,src='items_nml_mod')
WRITE(lineBuffer,'(A,I0)')' Domain = ',Domain
CALL umPrint(lineBuffer,src='items_nml_mod')

WRITE(lineBuffer,'(2A)')' ancilfilename = ',ancilfilename
CALL umPrint(lineBuffer,src='items_nml_mod')

WRITE(lineBuffer,'(A,L1)')' Update_anc = ',update_anc
CALL umPrint(lineBuffer,src='items_nml_mod')

WRITE(lineBuffer,'(A,I0)')' User_Prog_Ancil_stash_req = ',    &
                      User_Prog_Ancil_stash_req
CALL umPrint(lineBuffer,src='items_nml_mod')
WRITE(lineBuffer,'(A,G10.2)')' User_Prog_RConst = ',User_Prog_RConst
CALL umPrint(lineBuffer,src='items_nml_mod')

IF (source == netcdf_file) THEN
  WRITE(temp_format,'(A,I0,A)') '(A,', nitems ,'(A,1X))'
  WRITE(lineBuffer,TRIM(temp_format)) ' NetCDF_varname  = ', netcdf_varname
  CALL umPrint(lineBuffer,src='items_nml_mod')
END IF

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='items_nml_mod')

END SUBROUTINE print_nlist_items

!----------------------------------------------------------------------

SUBROUTINE read_nml_items(unit_number,         num_items,              &
                          stashcode_temp,      section_temp,           &
                          item_temp,           source_temp,            &    
                          area_temp,           user_prog_section_temp, &   
                          user_prog_item_temp, period_temp,            &
                          interval_temp,       real_constant_temp,     &
                          update_anc_temp,     items_filename_temp,    &
                          netcdf_varname_temp)

USE ancil_mod, ONLY: ancil_requests, num_ancil_requests, max_items, &
    stash_num_max
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1, ONLY: jprb, jpim
USE setup_namelist, ONLY: setup_nml_type
USE um_parcore, ONLY: mype
USE umPrintMgr, ONLY: umPrint, umMessage, PrintStatus, PrStatus_Oper, newline
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE
! Description:
!  Reads the items namelists into temporary arrays and passes them
!  back to calling routine.
!
INTEGER, INTENT(IN)  :: unit_number
INTEGER, INTENT(OUT) :: num_items
INTEGER, INTENT(OUT) :: stashcode_temp ( max_items )
INTEGER, INTENT(OUT) :: section_temp ( max_items )
INTEGER, INTENT(OUT) :: item_temp ( max_items )
INTEGER, INTENT(OUT) :: source_temp ( max_items )
INTEGER, INTENT(OUT) :: area_temp ( max_items )
INTEGER, INTENT(OUT) :: user_prog_section_temp ( max_items )
INTEGER, INTENT(OUT) :: user_prog_item_temp ( max_items )
INTEGER, INTENT(OUT) :: period_temp ( max_items )
INTEGER, INTENT(OUT) :: interval_temp ( max_items )
REAL,    INTENT(OUT) :: real_constant_temp ( max_items )
LOGICAL, INTENT(OUT) :: update_anc_temp ( max_items )
CHARACTER (LEN=filenamelength), INTENT(OUT) :: items_filename_temp (max_items)
CHARACTER (LEN=netcdf_varname_len),INTENT(OUT) :: netcdf_varname_temp(max_items)

INTEGER, PARAMETER :: no_of_types = 4
INTEGER, PARAMETER :: n_int = (9 * max_items) + 1
INTEGER, PARAMETER :: n_real = max_items
INTEGER, PARAMETER :: n_log  = max_items
INTEGER, PARAMETER :: n_chars = (filenamelength + netcdf_varname_len)*max_items
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode
INTEGER :: i
INTEGER :: Errorstatus
INTEGER :: iostatus
CHARACTER (LEN=errormessagelength):: Cmessage
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='READ_NML_ITEMS'

TYPE my_namelist
  SEQUENCE
  INTEGER :: num_items
  INTEGER :: stashcode_temp( max_items )
  INTEGER :: section_temp ( max_items )
  INTEGER :: item_temp ( max_items )
  INTEGER :: source_temp ( max_items )
  INTEGER :: area_temp ( max_items )
  INTEGER :: user_prog_section_temp ( max_items )
  INTEGER :: user_prog_item_temp ( max_items )
  INTEGER :: period_temp ( max_items )
  INTEGER :: interval_temp ( max_items )
  REAL    :: real_constant_temp ( max_items )
  LOGICAL :: update_anc_temp ( max_items )
  CHARACTER (LEN=filenamelength) :: items_filename_temp ( max_items )
  CHARACTER (LEN=netcdf_varname_len) :: netcdf_varname_temp ( max_items )
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

num_items = 0
iostatus  = 0

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int, &
              n_real_in=n_real, n_log_in = n_log, n_chars_in = n_chars)

IF (mype == 0) THEN

  ALLOCATE ( stash_req  (stash_num_max) )
  ALLOCATE ( netcdf_varname (stash_num_max) )

  ! Continue to read items namelists until end of file or other error code
  DO WHILE ( iostatus == 0 )

    ! Initialise namelist
    stash_req (:) = 0
    source  = 0
    domain  = 0
    period  = 0
    interval = 0
    user_Prog_Ancil_stash_req = 0
    user_Prog_RConst      = 0.0
    ancilfilename = 'empty'
    netcdf_varname (:) = netcdf_unset
    update_anc = .FALSE.

    READ( UNIT=unit_number, NML=items, IOSTAT=iostatus, IOMSG=iomessage )

    IF (iostatus == 0 ) THEN
      IF (PrintStatus >= PrStatus_Oper) CALL print_nlist_items()
      ! Only process the namelist if there are stash_req present
      IF (COUNT(stash_req /= 0) > 0) THEN

        DO i=1,COUNT(stash_req /= 0)
          IF (stash_req(i) <= 0) THEN
            Cmessage = 'Request for -ve or 0 Stashcode found in items namelist.'
            WRITE(umMessage,'(A)') Cmessage
            CALL umPrint(umMessage,src='items_nml_mod')
            ErrorStatus = 45
            CALL Ereport( RoutineName, ErrorStatus, Cmessage )
          END IF

          num_items = num_items + 1
      
          IF ( num_items <= max_items ) THEN
            stashcode_temp(num_items)  = stash_req(i)
            section_temp(num_items)    = stash_req(i)/1000
            item_temp(num_items)       = stash_req(i) - &
                                         (1000*section_temp(num_items))
            source_temp(num_items)     = source
            area_temp(num_items)       = domain
            period_temp(num_items)     = period
            update_anc_temp(num_items) = update_anc
            interval_temp(num_items)   = interval
            user_prog_section_temp(num_items) = user_prog_ancil_stash_req/1000
            user_prog_item_temp(num_items)    = user_prog_ancil_stash_req -    &
                                       (1000*user_prog_section_temp(num_items))
            items_filename_temp(num_items)    = ancilfilename
            real_constant_temp(num_items)     = user_prog_rconst

            IF (source == netcdf_file) THEN
              netcdf_varname_temp( num_items ) = netcdf_varname(i)
            ELSE
              netcdf_varname_temp( num_items ) = netcdf_unset
            END IF

            IF (ANY(stashcode_temp(1:(num_items - 1)) ==   &
                    stashcode_temp(num_items)) ) THEN
              WRITE(Cmessage,'(A,I5,A)')                                      &
                                 'Duplicate Stashcode entries for Stashcode', &
                                 stashcode_temp(num_items),                   &
                                 ' in items namelists.'
              WRITE(umMessage,'(A,A)') Cmessage, ' Please check entries.'
              CALL umPrint(umMessage,src='items_nml_mod')
              ErrorStatus = 40
              CALL Ereport( RoutineName, ErrorStatus, Cmessage )
            END IF
          END IF
        END DO
 
      ELSE
        ! To get here, an items namelist is either empty, is missing the 
        ! stash_req element, or has only zeros in stash_req.
        Cmessage = 'stash_req is missing or incompatible in items namelist.'//&
                   ' Please check entries.'
        WRITE(umMessage,'(A)') Cmessage
        CALL umPrint(umMessage,src='items_nml_mod')
        ErrorStatus = 10
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )
       END IF

      ! iostatus < 0 would indicate End Of File, which is an acceptable error
      ! iostatus > 0 indicates error trying to read the namelist itself.
    ELSE IF (iostatus > 0 ) THEN
      Cmessage = 'Error reading an items namelist. ' //           &
                 'Please check entries against code:'//  newline//&
                 'IoMsg: '//TRIM(iomessage)
      WRITE(umMessage,'(A)') Cmessage
      CALL umPrint(umMessage,src='items_nml_mod')
      ErrorStatus = 20
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )

    END IF ! if iostatus == 0

  END DO ! while iostatus == 0

  IF ( num_items > 0  .AND. num_items > max_items ) THEN
    WRITE(umMessage,'(A)') 'Maximum number of ITEMS namelists ' // &
                     'exceeded - please use a branch to increase.'
    CALL umPrint(umMessage,src='items_nml_mod')
    WRITE(umMessage,'(A,I5)')  'Number of ITEMS requests read in : ',num_items
    CALL umPrint(umMessage,src='items_nml_mod')
    WRITE(umMessage,'(A,I5)')  'Max no of ITEMS requests allowed : ',max_items
    CALL umPrint(umMessage,src='items_nml_mod')
    ErrorStatus = 30
    CALL Ereport( RoutineName, ErrorStatus, ummessage )
  END IF  !  if num_items > 0

END IF  !! mype==0

IF (mype==0) THEN
  my_nml % num_items              = num_items
  my_nml % stashcode_temp         = stashcode_temp
  my_nml % section_temp           = section_temp
  my_nml % item_temp              = item_temp
  my_nml % source_temp            = source_temp
  my_nml % area_temp              = area_temp
  my_nml % period_temp            = period_temp
  my_nml % interval_temp          = interval_temp
  my_nml % real_constant_temp     = real_constant_temp
  my_nml % items_filename_temp    = items_filename_temp
  my_nml % user_prog_item_temp    = user_prog_item_temp
  my_nml % user_prog_section_temp = user_prog_section_temp
  my_nml % update_anc_temp        = update_anc_temp
  my_nml % netcdf_varname_temp    = netcdf_varname_temp
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN
  num_items              = my_nml % num_items
  stashcode_temp         = my_nml % stashcode_temp
  section_temp           = my_nml % section_temp
  item_temp              = my_nml % item_temp
  source_temp            = my_nml % source_temp
  area_temp              = my_nml % area_temp
  period_temp            = my_nml % period_temp        
  interval_temp          = my_nml % interval_temp      
  real_constant_temp     = my_nml % real_constant_temp
  items_filename_temp    = my_nml % items_filename_temp
  user_prog_item_temp    = my_nml % user_prog_item_temp
  user_prog_section_temp = my_nml % user_prog_section_temp
  update_anc_temp        = my_nml % update_anc_temp
  netcdf_varname_temp    = my_nml % netcdf_varname_temp
END IF

IF (mype==0) DEALLOCATE (stash_req)
IF (mype==0) DEALLOCATE (netcdf_varname)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_items
END MODULE items_nml_mod

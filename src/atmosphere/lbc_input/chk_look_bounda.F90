! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine CHK_LOOK_BOUNDA
!
! Purpose : Cross checks values in LOOKUP records of boundary data
!           with model run values
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: LBC Input

SUBROUTINE chk_look_bounda(                                       &
  item_list,full_lookup_bounda,                                   &
                           icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE UM_ParVars
USE Control_Max_Sizes
USE rimtypes
USE lbc_mod
USE lookup_addresses
USE ppxlook_mod, ONLY: exppxi

USE cppxref_mod, ONLY:                                            &
    ppx_grid_type,ppx_halo_type,                                  &
    ppx_lv_code,ppx_lb_code, ppx_lt_code

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY:                                   &
    a_len_inthd, a_len_realhd, len1_lbc_comp_lookup, len1_lookup, &
    len_fixhd

USE atm_boundary_headers_mod, ONLY: lookup_bounda

USE errormessagelength_mod, ONLY: errormessagelength

USE pr_look_mod, ONLY: pr_look

IMPLICIT NONE


INTEGER ::                                                        &
  item_list(rim_lookupsa)                                         &
                             ! IN: STASH codes of expected items
, full_lookup_bounda(len1_lookup,bound_lookupsa)                  &
                             ! IN: Full LOOKUP record for LBCs
, icode                     ! OUT : Return code

CHARACTER(LEN=errormessagelength) ::                              &
  cmessage                  ! OUT : Error message

! Functions
INTEGER ::                                                        &
  get_fld_type

! Local variables
INTEGER ::                                                        &
  variable                                                        &
                    ! Loop counter for variable
, code                                                            &
                    ! item_code value for variable
, model                                                           &
                    ! model number for variable
, section                                                         &
                    ! section number for variable
, item                                                            &
                    ! section number for variable
, grid_type                                                       &
                    ! grid code for variable
, fld_type                                                        &
                    ! P,U or V for variable
, halo_type                                                       &
                    ! halo type for variable
, level_type                                                      &
                    ! what type of level for variable
, bottom_level_code                                               &
                    ! bottom level code
, bottom_level                                                    &
                    ! bottom level for variable
, top_level_code                                                  &
                    ! top level code
, top_level                                                       &
                    ! top level for variable
, n_levels_expected                                               &
                    ! number of levels expected
, n_levels_lbc                                                    &
                    ! number of levels in file
, halo_x_expected                                                 &
                    ! expected size of halo in x
, halo_x_lbc                                                      &
                    ! actual size of halo in x
, halo_y_expected                                                 &
                    ! expected size of halo in y
, halo_y_lbc                                                      &
                    ! actual size of halo in y
, size_x_expected                                                 &
                    ! expected size of field in x
, size_x                                                          &
                    ! actual size of field in x
, size_y_expected                                                 &
                    ! expected size of field in y
, size_y                                                          &
                    ! actual size of y
, rim_type                                                        &
                    ! Type of RIMWIDTH
, rimwidth_expected                                               &
                    ! expected rimwidth
, rimwidth_lbc                                                    &
                    ! actual rimwidth
, size_expected     ! Expected size

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHK_LOOK_BOUNDA'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!=====================================================================

! 1.0 Check that the first field in the file is the orography field

IF (lookup_bounda(item_code,1)  /=  31001) THEN ! Orography
  icode=1
  cmessage='CHK_LOOK_BOUNDA : No orography in LBC file'
  GO TO 9999
END IF

! 2.0 Check for the expected number of variables for each time.
!     Bear in mind that RIM_LOOKUPSA contains an extra field
!     (orography) which only occurs at the start of the field.
!     So the actual number of fields written out at each LBC output
!     time is actually RIM_LOOKUPSA-1

IF (full_lookup_bounda(item_code,2)  /=                           &
    full_lookup_bounda(item_code,2+rim_lookupsa-1)) THEN
  WRITE(umMessage,*) 'Wrong number of LBC variables found in LBC file'
  CALL umPrint(umMessage,src='chk_look_bounda')
  WRITE(umMessage,*) 'Expecting record ',2+rim_lookupsa-1,' to contain ', &
             'STASH item code ',full_lookup_bounda(item_code,2)
  CALL umPrint(umMessage,src='chk_look_bounda')
  WRITE(umMessage,*) 'But found item code ',                              &
             full_lookup_bounda(item_code,2+rim_lookupsa-1)
  CALL umPrint(umMessage,src='chk_look_bounda')
  icode=2
  cmessage='CHK_LOOK_BOUNDA : Wrong number of LBC fields'
  GO TO 9999
END IF

! 3.0 Now check the header record for each required variable

DO variable=1,rim_lookupsa

  code=item_list(variable)  ! item_code value for variable

  IF (lookup_bounda(item_code,variable)  /=  code) THEN
    WRITE(umMessage,*) 'Unexpected field in LBC file'
    CALL umPrint(umMessage,src='chk_look_bounda')
    WRITE(umMessage,*) 'Field ',variable,' was expected to be ',          &
               code,' but found ',                                &
               lookup_bounda(item_code,variable)
    CALL umPrint(umMessage,src='chk_look_bounda')

    CALL pr_look(lookup_bounda,variable)

    icode=3
    cmessage='CHK_LOOK_BOUNDA : Unexpected field in LBC file'
    GO TO 9999
  END IF

  model=lookup_bounda(model_code,variable)
  item=MOD(lookup_bounda(item_code,variable),1000)
  section=(lookup_bounda(item_code,variable)-item)/1000

  grid_type=exppxi(model,section,item,ppx_grid_type,              &
                   icode, cmessage)
  fld_type=get_fld_type(grid_type)

  halo_type=exppxi(model,section,item,ppx_halo_type,              &
                   icode, cmessage)

  level_type=exppxi(model,section,item,ppx_lv_code,               &
                    icode, cmessage)
  IF (level_type  ==  5) THEN
    n_levels_expected=1
  ELSE
    bottom_level_code=exppxi(model,section,item,ppx_lb_code,      &
                          icode, cmessage)
    top_level_code=exppxi(model,section,item,ppx_lt_code,         &
                          icode, cmessage)
    ! DEPENDS ON: levcod
    CALL levcod(bottom_level_code,bottom_level,icode,cmessage)
    ! DEPENDS ON: levcod
    CALL levcod(top_level_code,top_level,icode,cmessage)
    n_levels_expected=top_level-bottom_level+1
  END IF

  halo_x_expected=halosize(1,halo_type)
  halo_y_expected=halosize(2,halo_type)
  size_x_expected=glsize(1,fld_type)
  size_y_expected=glsize(2,fld_type)

  IF (lookup_bounda(item_code,variable)  ==  31001) THEN
    ! Orography
    rim_type=rima_type_orog
  ELSE
    rim_type=rima_type_norm
  END IF

  rimwidth_expected=rimwidtha(rim_type)

  halo_x_lbc=MOD(lookup_bounda(lbuser3,variable),100)
  halo_y_lbc=MOD(lookup_bounda(lbuser3,variable)-halo_x_lbc,      &
                 10000)/100
  rimwidth_lbc=MOD((lookup_bounda(lbuser3,variable)-              &
                 halo_x_lbc-halo_y_lbc*100),1000000)/10000

  n_levels_lbc=lookup_bounda(lbhem,variable)-100

  size_expected=global_LENRIMA(fld_type,halo_type,rim_type)*      &
                n_levels_expected

  IF (n_levels_lbc  /=  n_levels_expected) THEN
    WRITE(umMessage,*) 'Wrong number of levels for LBC field ',variable
    CALL umPrint(umMessage,src='chk_look_bounda')
    WRITE(umMessage,*) 'Expected ',n_levels_expected,' levels but found ',&
               n_levels_lbc
    CALL umPrint(umMessage,src='chk_look_bounda')

    CALL pr_look(lookup_bounda,variable)

    icode=4
    cmessage='CHK_LOOK_BOUNDA : Wrong number of levels'
    GO TO 9999
  END IF

  IF ((halo_x_lbc  /=  halo_x_expected) .OR.                      &
      (halo_y_lbc  /=  halo_y_expected)) THEN
    WRITE(umMessage,*) 'Incorrect halos for LBC field ',variable
    CALL umPrint(umMessage,src='chk_look_bounda')
    WRITE(umMessage,*) 'Expected halo_x= ',halo_x_expected,               &
               ' and halo_y= ',halo_y_expected
    CALL umPrint(umMessage,src='chk_look_bounda')
    WRITE(umMessage,*) 'but found halo_x= ',halo_x_lbc,                   &
               ' and halo_y= ',halo_y_lbc
    CALL umPrint(umMessage,src='chk_look_bounda')

    CALL pr_look(lookup_bounda,variable)

    icode=5
    cmessage='CHK_LOOK_BOUNDA : Incorrect halos'
    GO TO 9999
  END IF

  IF (rimwidth_lbc  /=  rimwidth_expected) THEN
    WRITE(umMessage,*) 'Wrong RIMWIDTH for LBC field ',variable
    CALL umPrint(umMessage,src='chk_look_bounda')
    WRITE(umMessage,*) 'Expected RIMWIDTH= ',rimwidth_expected,           &
               'but found RIMWIDTH= ',rimwidth_lbc
    CALL umPrint(umMessage,src='chk_look_bounda')

    CALL pr_look(lookup_bounda,variable)

    icode=6
    cmessage='CHK_LOOK_BOUNDA : Wrong RIMWIDTH'
    GO TO 9999
  END IF

  size_x=lookup_bounda(lbnpt,variable)
  size_y=lookup_bounda(lbrow,variable)

  IF ((size_x  /=                                                 &
       size_x_expected) .OR.                                      &
      (size_y  /=                                                 &
       size_y_expected)) THEN
    WRITE(umMessage,*) 'Incorrect dimensions for LBC field ',variable
    CALL umPrint(umMessage,src='chk_look_bounda')
    WRITE(umMessage,*) 'Expected ROW_LENGTH= ',size_x_expected,' and ',   &
               'ROWS= ',size_y_expected
    CALL umPrint(umMessage,src='chk_look_bounda')
    WRITE(umMessage,*) 'But found ROW_LENGTH= ',                          &
               size_x,' and ROWS= ',size_y
    CALL umPrint(umMessage,src='chk_look_bounda')

    CALL pr_look(lookup_bounda,variable)

    icode=7
    cmessage='CHK_LOOK_BOUNDA : Wrong dimensions'
    GO TO 9999
  END IF

  IF (lookup_bounda(lblrec,variable)  /=  size_expected) THEN
    WRITE(umMessage,*) 'Wrong size for LBC field ',variable
    CALL umPrint(umMessage,src='chk_look_bounda')
    WRITE(umMessage,*) 'Expected size was ',size_expected,' but found ',  &
               lookup_bounda(lblrec,variable)
    CALL umPrint(umMessage,src='chk_look_bounda')

    CALL pr_look(lookup_bounda,variable)

    icode=8
    cmessage='CHK_LOOK_BOUNDA : Wrong size'
    GO TO 9999
  END IF

END DO ! variable

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE chk_look_bounda

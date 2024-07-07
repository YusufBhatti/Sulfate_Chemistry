! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine INITDUMP -------------------------------------------
!
!   Purpose:  To read atmosphere dumps, and to calculate
!   additional constants based on the dump header information.
!
!   Extra constants needed for cloud types calulated within SETDCFLD
!
!   Programming Standard : UM documentation paper no. 3
!
!                   U.M. Documentation paper no F3
!
!
!
!   Arguments
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE initdump(                                              &
           sm_ident,icode,cmessage,dump_unit)

USE um_readdump_mod
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE model_file
USE io, ONLY: setpos
USE io_constants, ONLY: ioNoDelete
USE atm_fields_bounds_mod

USE UM_ParVars
USE decomp_params, ONLY: decomp_standard_atmos
USE Control_Max_Sizes
USE Decomp_DB
USE surf_couple_allocate_mod, ONLY: surf_couple_allocate
USE lookup_addresses
USE nlstgen_mod, ONLY: dump_packim
USE dump_headers_mod, ONLY: a_fixhd, a_inthd, a_cfi1, a_cfi2, a_cfi3,&
                            a_lookup, a_mpp_lookup, a_realhd,        &
                            a_levdepc, a_rowdepc, a_coldepc,         &
                            a_flddepc_in, a_extcnst, a_dumphist
USE submodel_mod, ONLY: atmos_sm, submodel_for_sm, atmos_im
USE stash_array_mod, ONLY: nitems, si, nsects
USE cderived_mod, ONLY: elf
USE history, ONLY: h_stepim, checkpoint_dump_im
USE nlcfiles_namelist_mod, ONLY: dump_filename => astart
USE filenamelength_mod, ONLY: filenamelength
USE file_manager, ONLY: assign_file_unit, release_file_unit

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE jules_sea_seaice_mod, ONLY: nice, nice_use

USE nlsizes_namelist_mod, ONLY:                                      &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,  &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,   &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_data,  &
    a_len_extcnst, a_len_inthd, a_len_realhd, a_prog_len,            &
    a_prog_lookup, bl_levels, global_a_len_data, land_field,         &
    len1_lookup, len_dumphist, len_fixhd,                            &
    len_tot, model_levels, mpp_len1_lookup, n_cca_lev, n_obj_d1_max, &
    ntiles, sm_levels, st_levels, tpps_ozone_levels,                 &
    tr_lbc_ukca, tr_lbc_vars, tr_levels, tr_ukca, tr_vars

USE d1_array_mod, ONLY: d1, d1_addr, no_obj_d1

USE missing_data_mod, ONLY: imdi, rmdi

USE errormessagelength_mod, ONLY: errormessagelength

USE io_configuration_mod, ONLY: io_omp_passive_dump_reading

USE wait_policy_mod, ONLY: set_wait_policy, & 
                           default_policy,  &
                           passive_policy

USE address_check_mod, ONLY: address_check

IMPLICIT NONE


!   Arguments
!

INTEGER, INTENT(IN)  :: sm_ident  ! Sub-model indicator
INTEGER, INTENT(OUT) :: icode     ! Return code
INTEGER, INTENT(IN)  :: dump_unit ! unit attached to input data file
                                  ! (astart or checkpoint_dump_im)

CHARACTER(LEN=errormessagelength), INTENT(OUT) :: cmessage  ! Error message

LOGICAL ::                                                        &
   l_a_dump                                                       &
  ,l_a_prog_only
                      ! ) Switches set if only prognostic
                      ! ) fields to be read in.

LOGICAL :: readheader

INTEGER :: i,j,ij    ! Loop counters
INTEGER ::                                                        &
   len2_lookup                                                    &
  ,len_data                                                       &
  ,n_prog_flds                                                    &
                      ! No of prognostic fields in dumps
  ,n_prog_lookup                                                  &
  ,len_prog                                                       &
  ,tot_len_data                                                   &
  ,d1_addr_submodel_id  ! submodel id in D1_ADDR array
INTEGER :: a_mpp_addr(a_len2_lookup),                             &
           a_mpp_len(a_len2_lookup)
INTEGER :: info  ! return code from GC operations

INTEGER ::                                                        &
        ifld                                                      &
                      ! Loop variable
       ,im_ident                                                  &
                      ! internal model identifier
       ,im_index      ! internal model index for STASH arrays

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INITDUMP'

CHARACTER(LEN=filenamelength) :: filename   ! Selected filename for dump

! ---------------------------------------------------------------------
!    Internal Structure

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!  1.0 Read atmosphere dump and initialise atmosphere model.
IF (sm_ident == atmos_sm) THEN

  !  1.1 read fixed length header and set buffer length

  IF (H_STEPim(atmos_im) == 0) THEN
    ! In the first timestep, start from the ASTART dump
    ! file already open on unit dump_unit, set pointer to start of file
    filename = dump_filename
  ELSE
    ! Not the first timestep, so a continuation run. In this case, the 
    ! checkpoint_dump_im file is read rather than astart
    filename = checkpoint_dump_im(atmos_im)
  END IF

  ! Set the file pointer to point to the beginning of the dump file (which
  ! may be the astart file or the checkpoint_dump_im file)
  CALL setpos (dump_unit,0,icode)
  IF (icode  >   0) THEN
    WRITE(umMessage,'(A)') ' Error resetting pointer to start of file'
    CALL umPrint(umMessage,src='initdump')
    GO TO 9999
  END IF

  ! DEPENDS ON: read_flh
  CALL read_flh (dump_unit,a_fixhd,len_fixhd,icode,cmessage)
  IF (icode >  0) THEN
    GO TO 9999
  END IF


  CALL setpos (dump_unit,0,icode)

  !       Test if atmos dump.
  l_a_dump = a_fixhd(5) == 1 .AND. a_fixhd(2) == atmos_sm

  !       Test if only prognostic fields to be read in
  l_a_prog_only = l_a_dump .AND. H_STEPim(atmos_im) == 0

  !       Get no of prognostic fields in atmos dump
  n_prog_flds = a_fixhd(153)

  !       Check N_PROG_FLDS has been set.
  IF (n_prog_flds == imdi) THEN
    WRITE(umMessage,'(A)') '  '
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A)') ' No of prognostic fields not set in FIXHD(153)'
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A)') ' Run RECONFIGURATION to set FIXHD(153)'
    CALL umPrint(umMessage,src='initdump')
    cmessage = 'INITDUMP: FIXHD(153) not set in atmos dump'
    icode = 101
    GO TO 9999  !  Return
  END IF

  !       Check N_PROG_FLDS matches with A_PROG_LOOKUP set up by the UI
  n_prog_lookup = a_prog_lookup
  IF (n_prog_flds /= n_prog_lookup) THEN
    WRITE(umMessage,'(A)') ' '
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A)') ' Mismatch in no of prognostic fields.'
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A,I0)') ' No of prog fields in Atmos dump ',n_prog_flds
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A,I0)') ' No of prog fields expected      ',n_prog_lookup
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A)') ' '
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A)') ' Run RECONFIGURATION to get correct no of'       &
             // ' prognostic fields in atmos dump'
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A)') ' or'
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A)') ' Check/Reset experiment in User Interface'
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A)') ' '
    CALL umPrint(umMessage,src='initdump')
    cmessage = 'INITDUMP: Wrong no of atmos prognostic fields'
    icode = 102
    GO TO 9999  !  Return
  END IF

  ! Initialise D1 to prevent uninitialised data in unused rows of U fields
  DO i = 1,len_tot
    d1(i)=0.0
  END DO
  !       Determine no of fields to be read in
  IF (l_a_prog_only) THEN

    !         Prognostic fields only
    len2_lookup = n_prog_flds
    len_prog = a_prog_len
    tot_len_data = a_len_data

    len_data = len_prog

    WRITE(umMessage,'(A)') ' '
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A,I5,A)') ' Read in ',n_prog_flds,' prognostic fields.'
    CALL umPrint(umMessage,src='initdump')


    !      INITIALISE DIAGNOSTIC AREA OF D1 TO RMDI
    DO i = len_data+1, tot_len_data
      d1(i)=rmdi
    END DO

  ELSE

    !         All fields.
    len2_lookup = a_len2_lookup
    len_data    = a_len_data

    WRITE(umMessage,'(A)') ' '
    CALL umPrint(umMessage,src='initdump')
    WRITE(umMessage,'(A,I5,A)') ' Read in all ',len2_lookup,' fields.'
    CALL umPrint(umMessage,src='initdump')

  END IF

  ! Ensure that domain decomposition is consistent with submodel

  CALL change_decomposition(decomp_standard_atmos,icode)

  IF (icode  >   0) THEN
    GO TO 9999
  END IF

  ! If requested, change the OMP waiting policty to passive. This
  ! could have a performance benefit on reading the dump.
  IF (io_omp_passive_dump_reading) CALL set_wait_policy(passive_policy)

  !  1.2 Call READDUMP to read atmosphere dump.

  d1_addr_submodel_id = submodel_for_sm(atmos_sm)
  readheader=.TRUE.
  CALL um_readdump(dump_unit, a_fixhd, len_fixhd,                   &
      a_inthd, a_len_inthd,                                         &
      a_realhd, a_len_realhd,                                       &
      a_levdepc, a_len1_levdepc, a_len2_levdepc,                    &
      a_rowdepc, a_len1_rowdepc, a_len2_rowdepc,                    &
      a_coldepc, a_len1_coldepc, a_len2_coldepc,                    &
      a_flddepc_in, a_len1_flddepc, a_len2_flddepc,                 &
      a_extcnst, a_len_extcnst,                                     &
      a_dumphist, len_dumphist,                                     &
      a_cfi1, a_len_cfi1,                                           &
      a_cfi2, a_len_cfi2,                                           &
      a_cfi3, a_len_cfi3,                                           &
      a_lookup,len1_lookup,len2_lookup,                             &
      a_mpp_lookup,mpp_len1_lookup,                                 &
           atmos_sm,                                                &
           no_obj_d1(d1_addr_submodel_id),                          &
           d1_addr(1,1,d1_addr_submodel_id),                        &
           len_data,d1,                                             &
      readheader)

  CALL model_file_close(dump_unit, filename, delete=ioNoDelete,     &
                        error=icode)
  CALL release_file_unit(dump_unit, handler="portio")

  IF (io_omp_passive_dump_reading) CALL set_wait_policy(default_policy)

  IF (icode  >   0) THEN
    GO TO 9999
  END IF

  ! Check validity of integer header data and print out information.
  ! Pass through the global numbers so the validity check works
  ! glsize(1) is the global ROW_LENGTH
  ! glsize(2) is the global ROWS
  ! DEPENDS ON: pr_inhda
  CALL pr_inhda (a_inthd, a_len_inthd, glsize(1,fld_type_p),      &
     glsize(2,fld_type_p),model_levels,tr_levels,                 &
     st_levels, sm_levels, bl_levels,                             &
 tr_vars, icode, cmessage)

  IF (icode  >   0) THEN
    GO TO 9999
  END IF

  ! Check validity of real header data and print out information.
  ! DEPENDS ON: pr_rehda
  CALL pr_rehda (a_realhd, a_len_realhd)

  IF (l_a_prog_only) THEN

    !         Need to pass field address and length info to ADDRESS_CHECK
    DO i=1,len2_lookup
      a_mpp_addr(i) = a_mpp_lookup(p_naddr,i)
      a_mpp_len(i)  = a_mpp_lookup(p_lblrec,i)
    END DO

    CALL address_check (a_lookup,a_mpp_addr,                      &
      a_mpp_len,len1_lookup,len2_lookup,                          &
                        si,nitems,nsects,len_data,                &
                        icode,cmessage)

    IF (icode  >   0) THEN
      GO TO 9999
    END IF

  END IF

  ! ------------------------------------------------------------
  !       Check that packing codes in lookup table is consistent
  !       with packing requested in dump_packim.
  ! --------------------------------------------------------------

  ! DEPENDS ON: check_dump_packing
  CALL check_dump_packing (                                       &
       a_fixhd, len_fixhd,                                        &
       a_lookup, len1_lookup, len2_lookup,                        &
       dump_packim(sm_ident), atmos_im )

  !       Reset A_FIXHD to correspond to Output Dump
  a_fixhd(152) = a_len2_lookup
  a_fixhd(160) = a_fixhd(150) + len1_lookup*a_len2_lookup
  a_fixhd(161) = global_A_LEN_DATA

  !  1.3 Call SET_ATM_POINTERS to set integer pointers to data in
  !      atmosphere dump and secondary storage area in D1 array.
  ! DEPENDS ON: set_atm_pointers
  CALL set_atm_pointers(                                            &
                   icode,cmessage)

  IF (icode  >   0) THEN
    GO TO 9999
  END IF

  ! Initialise Jules
  CALL surf_couple_allocate( land_field, ntiles, sm_levels, &
                                  nice ,nice_use )

  ! SETCONA is now called from INITIAL
  ! Set ELF flag
  elf=(a_fixhd(4) == 3 .OR. a_fixhd(4) == 103)

END IF

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE initdump

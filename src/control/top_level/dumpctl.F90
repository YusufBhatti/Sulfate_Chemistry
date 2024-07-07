! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Routine: DUMPCTL
!
!   Purpose: Controls the production and naming of output dump files.
!
!   Programming standard: UMDP3
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!
!   This file belongs in section: Top Level
!
SUBROUTINE dumpctl (                                                           &
    analysis)

USE yomhook,                 ONLY: lhook, dr_hook
USE parkind1,                ONLY: jprb, jpim
USE filenamelength_mod,      ONLY: filenamelength
USE iau_mod,                 ONLY: l_iau
USE nlstgen_mod,             ONLY: dump_filename_base
USE nlstcall_mod,            ONLY: model_analysis_mins, l_fastrun
USE nlcfiles_namelist_mod,   ONLY: atmanl_filename => atmanl
USE filename_generation_mod, ONLY: get_filename
USE missing_data_mod,        ONLY: rmdi
USE model_time_mod,          ONLY: secs_per_stepim
USE file_manager,            ONLY: assign_file_unit, release_file_unit
USE history,                 ONLY: checkpoint_dump_im
USE lookup_addresses,        ONLY: lblrec
USE ereport_mod,             ONLY: ereport

USE atm_fields_bounds_mod, ONLY:                                               &
    pdims, wdims, tdims, udims, vdims, wdims_s, o3dims2

USE io, ONLY:                                                                  &
    iodisksynchronise, iofence

USE io_constants, ONLY:                                                        &
    ioopenreadwrite, ionodelete, ioFileTypeMPIIO

USE model_file, ONLY:                                                          &
    model_file_open, model_file_close, storealllookups, synchroniseall

USE umprintmgr, ONLY:                                                          &
    ummessage, umprint

USE submodel_mod, ONLY: atmos_sm, submodel_for_sm, atmos_im
USE dump_headers_mod, ONLY: a_fixhd, a_inthd, a_cfi1, a_cfi2, a_cfi3, a_lookup,&
                            a_mpp_lookup, a_realhd, a_levdepc, a_rowdepc,      &
                            a_coldepc, a_flddepc, a_extcnst, a_dumphist
USE d1_array_mod, ONLY: d1, d1_addr, no_obj_d1

USE nlsizes_namelist_mod, ONLY:                                                &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,            &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,             &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_data,            &
    a_len_extcnst, a_len_inthd, a_len_realhd, land_field, len1_lookup,         &
    len_dumphist, len_fixhd, len_tot,                                          &
    model_levels, mpp_len1_lookup, n_cca_lev, n_obj_d1_max, sm_levels,         &
    st_levels, tpps_ozone_levels, tr_lbc_ukca, tr_lbc_vars, tr_ukca,           &
    tr_vars

IMPLICIT NONE

! Arguments
LOGICAL, INTENT(IN)  :: analysis  ! Flag indicating this is the analysis dump

! Local variables
INTEGER :: i                    ! Loop counter
INTEGER :: icode                ! Error return code
INTEGER :: buflen               ! Length of i/o buffer for WRITDUMP
INTEGER :: disk_address         ! Current rounded disk address 
INTEGER :: d1_addr_submodel_id  ! submodel id in D1_ADDR array
INTEGER :: name_offset          ! Time offset for dump name
INTEGER :: dump_unit            ! File unit for dump file
INTEGER :: number_of_data_words_on_disk
INTEGER :: number_of_data_words_in_memory 

CHARACTER(LEN=filenamelength) &
        :: dumpname = ""        ! Model generated dump name

! Dr-hook
CHARACTER(LEN=*),   PARAMETER :: RoutineName = "DUMPCTL"
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

DO i = 1, len_dumphist + 1
  a_dumphist(i) = rmdi
END DO

IF (analysis) THEN
  ! This is the write of the analysis dump - take its name directly from the
  ! nlcfiles namelist entry
  dumpname = atmanl_filename
ELSE
  ! Otherwise this is a normal dump - but its name might still need to be 
  ! adjusted to take account of the offset for the analysis 
  IF (l_fastrun .OR. l_iau) THEN
    name_offset = - (model_analysis_mins * 60)/secs_per_stepim(atmos_im)
  ELSE
    name_offset = 0
  END IF

  ! Generate a name for the dump using the pattern provided by the user
  CALL get_filename(TRIM(dump_filename_base),             &
                    dumpname,                             &
                    relative_offset_steps = name_offset,  &
                    dump=.TRUE.)
END IF

! Write all pp files' cached lookups to ensure that current data is
! is in the files' output stream
CALL storealllookups()

! Flush all pp files' units to ensure that data is commited from the
! application
CALL synchroniseall()

! Compute the new addresses and lengths

! DEPENDS ON: set_dumpfile_address
CALL set_dumpfile_address(                                        &
   a_fixhd, len_fixhd,                                            &
   a_lookup, len1_lookup, a_len2_lookup,                          &
   number_of_data_words_in_memory, number_of_data_words_on_disk,  &
   disk_address)

! Open unit for dump :
CALL assign_file_unit(dumpname, dump_unit, handler="portio")
WRITE(ummessage,'(A,A,A,I3)')"DUMPCTL: Opening new file ",TRIM(dumpname), &
   " on unit ",dump_unit
CALL umPrint(ummessage,src=routinename)

CALL model_file_open(dump_unit, dumpname, read_write=ioopenreadwrite, &
                     allowremap=.TRUE., fileType=ioFileTypeMPIIO )

! ----------------------------------------------------------------------
!  5. Write dump on appropriate unit putting timestamp in header
!
! Creation date and time
CALL date_time(a_fixhd(35),a_fixhd(36),a_fixhd(37),             &
   a_fixhd(38),a_fixhd(39),a_fixhd(40))

! Maximum length of field, required for IO buffer
buflen = a_lookup(lblrec,1)
IF (a_len2_lookup >  1) THEN
  DO i = 2, a_len2_lookup
    buflen = MAX(buflen, a_lookup(lblrec,i))
  END DO
END IF

d1_addr_submodel_id = submodel_for_sm(atmos_sm)

! DEPENDS ON: um_writdump
CALL um_writdump(dump_unit,a_fixhd,len_fixhd,        &
   a_inthd,a_len_inthd,                              &
   a_realhd,a_len_realhd,                            &
   a_levdepc,a_len1_levdepc,a_len2_levdepc,          &
   a_rowdepc,a_len1_rowdepc,a_len2_rowdepc,          &
   a_coldepc,a_len1_coldepc,a_len2_coldepc,          &
   a_flddepc,a_len1_flddepc,a_len2_flddepc,          &
   a_extcnst,a_len_extcnst,                          &
   a_dumphist,len_dumphist,                          &
   a_cfi1,a_len_cfi1,                                &
   a_cfi2,a_len_cfi2,                                &
   a_cfi3,a_len_cfi3,                                &
   a_lookup,len1_lookup,a_len2_lookup,               &
   a_mpp_lookup,mpp_len1_lookup,                     &
   buflen,                                           &
   atmos_sm,                                         &
   no_obj_d1(d1_addr_submodel_id),                   &
   d1_addr(1,1,d1_addr_submodel_id),                 &
   a_len_data,d1)

a_fixhd(5)=1    ! Set fixhd(5) back to instantaneous dump

! Commit outstanding data to disk
CALL iodisksynchronise(dump_unit)

! Ensure the file_close occurs after the PP files have been updated
! (SynchroniseAll has completed)
CALL iofence(dump_unit)

! Close the file
CALL model_file_close(dump_unit,dumpname,delete=ionodelete,error=icode)
IF (icode /= 0) THEN
  CALL ereport(routinename, icode, "Problem closing dump")
END IF
CALL release_file_unit(dump_unit, handler="portio")

! Save dump name for writing to restart file on successful completion
! of timestep
checkpoint_dump_im(atmos_im) = dumpname

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE dumpctl



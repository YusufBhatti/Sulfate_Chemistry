! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
SUBROUTINE readcntl(icode,cmessage)

USE check_iostat_mod
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

! JULES
USE switches, ONLY: l_360

USE filenamelength_mod, ONLY:                                    &
          filenamelength

! MPP
USE mpp_conf_mod, ONLY: nlst_mpp, extended_halo_size_ew,          &
                        extended_halo_size_ns, gcom_coll_limit,   &
                        read_nml_nlst_mpp, check_nml_nlst_mpp

USE UM_ParVars
USE io, ONLY: setpos, buffin
USE io_constants, ONLY: ioOpenReadOnly, ioNoDelete
USE model_file, ONLY: model_file_open, model_file_close
USE lbc_mod
USE ancilcta_namelist_mod, ONLY: read_nml_ancilcta, check_nml_ancilcta, &
                                 print_nlist_ancilcta
USE model_id_mod, ONLY: itab, configid, read_nml_configid, check_configid
USE nlstgen_mod
USE check_nlstcgen_mod, ONLY: check_nlstcgen
USE missing_data_mod, ONLY: imdi
USE nlstcall_mod, ONLY: read_nml_nlstcall, print_nlist_nlstcall, &
                        lcal360, Num_ALBCs
USE nlstcall_pp_namelist_mod, ONLY: read_nml_nlstcall_pp
USE nlstcall_nc_namelist_mod, ONLY: read_nml_nlstcall_nc, &
                                    read_nml_nlstcall_nc_options, &
                                    l_netcdf
USE file_manager, ONLY: get_file_unit_by_id, assign_file_unit, &
                        release_file_unit, um_file_type,       &
                        get_file_by_unit

USE nlcfiles_namelist_mod, ONLY: lbc_file_1 => alabcin1

USE nlsizes_namelist_mod, ONLY: &
    len_fixhd

USE um_parcore, ONLY: mype, nproc_max

IMPLICIT NONE
!
! Description:
!  Reads overall, generic and model-specific control variables.
!
! Method:
!  Reads namelist containing all overall, generic and model-specific
!  control variables not in History file.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.
!
! Declarations:
!

! Subroutine arguments
!   Scalar arguments with intent(in):
!   Array  arguments with intent(in):
!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):
INTEGER :: icode ! return code
INTEGER :: errorstatus
CHARACTER (LEN=*), PARAMETER  :: RoutineName='READCNTL'
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=errormessagelength) :: iomessage
!   Array  arguments with intent(out):

! Local parameters:

! Local scalars:
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Local dynamic arrays:

INTEGER   :: error     ! Error code returned by OPEN
INTEGER   :: fixhdr(len_fixhd)
INTEGER   :: len_io
INTEGER   :: outtime
REAL      :: a
INTEGER, ALLOCATABLE  :: lookup(:,:)
INTEGER   :: atmoscntl_unit ! unit no. for ATMOSCNTL file
INTEGER   :: shared_unit    ! unit no. for SHARED file
INTEGER   :: pp_unit
INTEGER   :: lbc_unit

! Reading in of the various NLSTCALL namelists
INTEGER :: i_pp
INTEGER :: readstatus
LOGICAL :: end_of_file
TYPE(um_file_type), POINTER :: um_file

!- End of header

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!  To ensure the namelists work correctly, owing to dependence on variables
!  the namelists must be read in the following order:
!     NLSTCGEN
!     NLSTCALL

      ! Retrieve the unit number of the main namelist file
atmoscntl_unit = get_file_unit_by_id("atmoscntl", handler="fortran")
shared_unit   = get_file_unit_by_id("shared", handler="fortran")

CALL read_nml_nlstcgen(atmoscntl_unit)

CALL read_nml_nlstcall_pp(atmoscntl_unit)
IF ( mype == 0 ) REWIND(atmoscntl_unit)

! Read netCDF output options
CALL read_nml_nlstcall_nc_options(atmoscntl_unit)
IF ( mype == 0 ) REWIND(atmoscntl_unit)

IF (l_netcdf) THEN
  CALL read_nml_nlstcall_nc(atmoscntl_unit)
  IF ( mype == 0 ) REWIND(atmoscntl_unit)
END IF

! Special setup for climate mean unit if required
IF (l_meaning_sequence) THEN
  CALL assign_file_unit("Reserved unit for mean files", pp_unit, &
                        handler="portio", id="aomean")
  NULLIFY(um_file)
  um_file => get_file_by_unit(pp_unit, handler="portio")
  um_file % meta % packing_code = ppxm
  um_file % pp_meta % init_file_type = "p"
  um_file % pp_meta % reserved_headers = 4096
  um_file % meta % is_output_file = .FALSE.
END IF

CALL read_nml_nlstcall(shared_unit)
CALL print_nlist_nlstcall()

! Check nlstcgen has to be called after NLSTCALL is read, as it is dependent
! on a variable in this namelist
CALL check_nlstcgen()

!  Read configuration id, is itab set explictly as an input?
CALL read_nml_configid(atmoscntl_unit)
CALL check_configid()

! Read MPP configuration data
CALL read_nml_nlst_mpp(atmoscntl_unit)
CALL check_nml_nlst_mpp()

! Check some items in NLST_MPP have been provided to avoid misleading
! error messages or undesired behaviour later. 
! Possibly unnecessary after CALL check_nml_nlst_mpp()
IF (extended_halo_size_ew == imdi) THEN
  cmessage= 'READCNTL: extended_halo_size_ew not set in namelist nlst_mpp'
  icode = 100
  CALL ereport(routinename, icode, cmessage)
END IF
IF (extended_halo_size_ns == imdi) THEN
  cmessage= 'READCNTL: extended_halo_size_ns not set in namelist nlst_mpp'
  icode = 110
  CALL ereport(routinename, icode, cmessage)
END IF
IF (gcom_coll_limit == imdi) THEN
  cmessage = 'READCNTL: gcom_coll_limit not set in namelist nlst_mpp'
  icode = 120
  CALL ereport(routinename, icode, cmessage)
END IF

CALL read_nml_ancilcta(shared_unit)
CALL check_nml_ancilcta()
CALL print_nlist_ancilcta()

! overwrite the JULES module values
l_360         = lcal360

IF (Num_ALBCs /= 0 ) THEN

  CALL assign_file_unit(lbc_file_1, lbc_unit, handler="portio")

  CALL model_file_open(lbc_unit, lbc_file_1, read_write=ioOpenReadOnly)

  CALL buffin(lbc_unit,fixhdr,len_fixhd,len_io,a)

  ! Check for I/O errors
  IF (a /= -1.0 .OR. len_io /= len_fixhd) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror('buffer in of fixed length header',a,len_io  &
       ,len_fixhd)
    cmessage='READCNTL: I/O error'
    icode=210
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  ! Broadcast the fixed length header from proc 0 to all PEs
  CALL gc_ibcast(77,len_fixhd,0,nproc_max,icode,fixhdr)

  IF (icode /= 0) THEN
    Cmessage = 'Error in broadcast'
    CALL Ereport(RoutineName, icode, cmessage)
  END IF

  ! Move to start of Look Up Table
  CALL setpos(lbc_unit,fixhdr(150)-1,icode)

  ALLOCATE (lookup(fixhdr(151),fixhdr(152)))

  ! Read in fields from LOOKUP table
  CALL buffin(lbc_unit,lookup(:,:),fixhdr(151)*fixhdr(152),len_io,a)

  IF (a /= -1.0 .OR. len_io /= fixhdr(151)*fixhdr(152)) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror('buffer in of lookup table',a,len_io,   &
           fixhdr(151)*fixhdr(152))
    cmessage='READCNTL: I/O error'
    icode=230
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  rimwidtha(1)=lookup(41,1)/10000

  ! Broadcast the rimwidth from proc 0 to all PEs
  CALL gc_ibcast(77,nrima_max,0,nproc_max,icode,rimwidtha)
  IF (icode /= 0) THEN
    Cmessage = 'Error in broadcast'
    CALL Ereport(RoutineName, icode, cmessage)
  END IF

  CALL model_file_close(lbc_unit, lbc_file_1,delete=ioNoDelete,     &
                        error=icode)
  CALL release_file_unit(lbc_unit, handler="portio")

  DEALLOCATE (lookup)
ELSE
  rimwidtha(1)=0
END IF



IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE readcntl

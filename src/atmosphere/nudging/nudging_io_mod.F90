! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module to hold procedures for handling NetCDF Nudging files.
!
!  Method:
!    The basic interaction with NetCDF files is handled via the existing 
!    procedures from the UKCA module: EMISS_IO_MOD. 
!    One new procedure is required for actual data input as :
!    1. The structure and dimensions of Nudging data files are different 
!       to those of the emission/ oxidant files.
!    2. Values read from file are not 'scattered' to other PEs as done for
!       emissions, but 'broadcast'.
!
!  The nudged model was developed as part of the UKCA project.
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Available for non-commerical use at present. Please cite
!
!  Telford, P. J., Braesicke, P., Morgenstern, O., and Pyle, J. A.:
!  Technical Note: Description and assessment of a nudged version of
!  the new dynamics Unified Model, Atmos. Chem. Phys., 8, 1701-1712, 2008.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Nudging
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE nudging_io_mod

USE um_types,               ONLY: integer32         ! For NetCDF arguments
USE netcdf
USE emiss_io_mod
USE um_parcore,             ONLY: mype
USE mpl,                    ONLY: mpl_real
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,             ONLY: UmPrint,UmMessage
USE parkind1,               ONLY: jpim, jprb         ! DrHook
USE yomhook,                ONLY: lhook, dr_hook     ! DrHook

IMPLICIT NONE

INTEGER, PARAMETER, PUBLIC       :: max_ncdims = 5  
                          ! Max num of dimens per variable,
                          ! has to match that in emiss_io_mod

! Subroutines available outside this module
PUBLIC :: ndg_get_data

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NUDGING_IO_MOD'

CONTAINS

!-------------------------------------------------------------------------
! Description:
!   Read data values of specified variable from a NetCDF file.
!   Currently only handles variable of type Real. 
!   All dimensions are passed as arguments so that a single routine can
!   read 1- to 4-D variables (although 4th/ time dimension is always 1)
!
! Method:
!   Get variable ID (varid), number of dimensions in variable (ndim), 
!   and vector of dimension sizes (dim_size) by calling EM_GET_VAR_INFO.
!   Then get domain/indices to extract and call NF90_GET_VAR to get the 
!   variable values: rvalue. Values are read by PE0 and broadcast
!   to all PEs.
!-------------------------------------------------------------------------
SUBROUTINE ndg_get_data (size_x, size_y, size_z, size_t, fileid, trec,  &
                         var_name, rvalue)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT (IN)    :: size_x     ! Expected file dimensions
INTEGER,           INTENT (IN)    :: size_y
INTEGER,           INTENT (IN)    :: size_z     
INTEGER,           INTENT (IN)    :: size_t     
INTEGER,           INTENT (IN)    :: fileid     ! file id
INTEGER,           INTENT (IN)    :: trec       ! time record wanted

CHARACTER (LEN=*), INTENT (IN)    :: var_name   ! variable name

REAL,              INTENT (INOUT) :: rvalue (size_x, size_y, size_z, size_t) 
                                      ! variable to hold file values
! Local variables

INTEGER (KIND=integer32) :: nstatus         ! netCDF return code

!  Vector specifying the index in the variable from which
!  the first data value will be read along each dimension
INTEGER (KIND=integer32) :: start   (max_ncdims)

!  Vector specifying the number of values to read along each dimension
INTEGER (KIND=integer32) :: counter (max_ncdims)

INTEGER (KIND=integer32) :: varid32
INTEGER                  :: varid, ndim, dim_size(max_ncdims), n, ilen

CHARACTER (LEN=15) :: varname_local ! Local variable copy for var_name 
                      ! -requirement of the EM_GET_VAR_INFO routine

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'NDG_GET_DATA'
                                    ! Name of routine for messages
INTEGER                  :: errcode ! return code for ereport

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: icode

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator (my_comm, icode)

! Get information on variable. One of the two INOUT arguments of
! EM_GET_VAR_INFO is supplied (varname_local); assign an arbitrary
! negative value to the other one since it is not known yet (varid).
varname_local = var_name    ! need to pass INOUT argument
varid         = -1
CALL em_get_var_info (fileid, varid, varname_local, ndims=ndim,        &
                      vdimsize=dim_size)

! Compute ilen: total size needed for upper bound of assumed-size array
! used in nf90_get_var. The total array size is the product of sizes along 
! each dimension.
! This is required by all PEs as it is an argument for the Broadcast routine.
ilen = 1
DO n = 1, ndim
  ilen = ilen * dim_size(n)
END DO

! Only PE0 reads data values from the file
IF ( mype == 0 ) THEN

  ! Indices/domain to extract. Default is to get complete data array for
  ! one time record from NetCDF file, so set sampling/ stride size equal to
  ! the file variable size in all but last (time) dimension.
  start   (1:ndim-1) = 1
  counter (1:ndim-1) = dim_size (1:ndim-1)

  ! All variables in Nudging NetCDF files contain a time dimension, which 
  ! is the last one. Even though each file currently contains only one
  ! record, there is a provision for multiple-time files in future. 
  ! Hence, compare the record no. requested against maximum time records.
  IF ( trec > dim_size(ndim) )  THEN
    WRITE (umMessage,'(3A,I4,A,I4)')                                    &
     'Requested time record out of bounds for ',TRIM(var_name),         &
       'Available: ',dim_size(ndim),'; Requested: ',trec
    CALL umPrint(umMessage,src=TRIM(routinename))
    errcode = 1
    CALL ereport (routinename, errcode,                                 &
        'Time record out of bounds for: ' // TRIM(var_name))
  END IF

  ! Read a single time from the file
  start   (ndim) = trec     ! Last index corresponds to record number.
  counter (ndim) = 1        ! Read only one record

  ! Check that dimensions in file match the expected i.e. output array size
  DO n = 1, ndim
    IF (SIZE (rvalue,n) /= counter (n)) THEN
      errcode = 1
      WRITE(umMessage,'(A,A,1X,3I6)')'Dimension mismatch for: ',        &
            TRIM(var_name),n,SIZE (rvalue,n),counter(n)
      CALL ereport (routinename, errcode,umMessage)
    END IF
  END DO

  !  Read values of the variable into rvalue, specifying what to read
  !  with the optional start and counter arguments
  nfid    = fileid   ! nfid : public variable from emiss_io_mod
  varid32 = varid    ! 32-bit integer needed as nf90 argument
  nstatus = nf90_get_var (nfid, varid32, rvalue, start   (1:ndim),      &
                                                 counter (1:ndim))
  CALL nd_error (nstatus, routinename, var_name) ! check for success
END IF   ! PE=0

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (rvalue, ilen, mpl_real, 0, my_comm, icode)
IF ( icode /= 0 ) THEN
  errcode = fileid
  CALL ereport (routinename, errcode, 'Broadcast Error')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ndg_get_data

END MODULE nudging_io_mod

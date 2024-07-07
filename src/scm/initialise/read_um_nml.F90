! *********************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *********************************COPYRIGHT*********************************
!
! Routine to read in SIZES control namelist. The namelist
! NLSIZES contains lots of variables that are not applicable to
! the SCM. To use the same namelists the UM generates without keeping the
! extra fields, the NLSIZES nameslist is read in this routine
! then only the required fields are passed to scm_shell.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: SCM
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

SUBROUTINE read_um_nml( filename )

USE jules_snow_mod, ONLY: nsmax
USE lbc_mod
USE ereport_mod, ONLY: ereport
USE filenamelength_mod, ONLY: filenamelength
USE free_tracers_inputs_mod, ONLY: a_max_trvars, i_free_tracer
USE nlsizes_namelist_mod, ONLY:                                           &
  nlsizes, global_rows, global_row_length, rows, row_length,              &
  sm_levels, st_levels, tr_vars, tr_ukca
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


CHARACTER (LEN=filenamelength), INTENT(IN) :: filename

! Local

CHARACTER(LEN=200) :: RoutineName = 'read_um_nml'
CHARACTER(LEN=errormessagelength) :: Cmessage
CHARACTER(LEN=errormessagelength) :: iomessage
INTEGER :: Icode, Istatus
INTEGER :: i

OPEN(10, FILE=filename, IOSTAT=Istatus, STATUS='old', ACTION='READ', &
     IOMSG=iomessage)
IF (Istatus /= 0) THEN
  Icode = 500
  WRITE(Cmessage,*) " Error opening " //TRIM(ADJUSTL(filename))//    &
                    " on unit 10: " // TRIM(iomessage)

  CALL ereport (RoutineName, Icode, Cmessage)
END IF

READ(10,nlsizes)

CLOSE(10)

! No decomposition is required in the SCM.
global_rows=1
global_row_length=1
rows       = global_rows
row_length = global_row_length

! The following are fixed for an SCM run:
rimwidtha=1
nrim_timesa=1

! setting soil mositure levels - same as temp levels for JULES
sm_levels=st_levels


! Logic for setting the num of passive tracers - count the non-zero
! elements of i_free_tracer, read in from run_free_tracers
tr_vars=0
DO i=1,A_Max_TrVars
  IF (i_free_tracer(i) /= 0)tr_vars = tr_vars + 1
END DO

! Logic for setting the num of ukca tracers - note that
! TC_UKCA has been removed. Full UM sets tr_ukca
! in primary using tsmsk logic in primary. Set
! tr_ukca to zero as considerable work would be
! needed to run UKCA in SCM mode.
tr_ukca=0

RETURN
END SUBROUTINE read_um_nml


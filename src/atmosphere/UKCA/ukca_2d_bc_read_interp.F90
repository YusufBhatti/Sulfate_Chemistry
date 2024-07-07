! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Read, broadcast and interpolate in time 2D data file for top B/C
!
MODULE ukca_2d_bc_read_interp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: & 
  ModuleName = 'UKCA_2D_BC_READ_INTERP_MOD'

CONTAINS

SUBROUTINE ukca_2d_bc_read_interp(mype,nproc,idofy,strat2d_dir, &
                                  noy,o3,ch4)

USE UKCA_phot2d, ONLY: nolat,nolev
USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE errormessagelength_mod, ONLY: errormessagelength
USE filenamelength_mod, ONLY: filenamelength
USE umPrintMgr, ONLY: newline

IMPLICIT NONE
!
! Description:
!   Reads in 2D data from an external file. Interpolate in time and space
!   to give appropriate data for use at the UKCA model's top boundary.
!   Called from UKCA_STRATF
!
! Method:
!   On the first timestep the files are opened on PE0 and the data is read in
!   then the data is broadcast to all PEs
!   On all timesteps the data is interpolated in time and the values
!   ch4, o3, noy are passed to the calling routine
!   Check input unit numbers if change these
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.2 programming standards.
!
! Declarations:
!
! Subroutine arguments
!
! Scalar arguments with INTENT(IN):
INTEGER, INTENT(IN) :: mype                       ! Processor number
INTEGER, INTENT(IN) :: nproc                      ! Number of processors
INTEGER, INTENT(IN) :: idofy                      ! Day of year
CHARACTER(LEN=* ), INTENT(IN)   ::  strat2d_dir   ! Directory containing input data
!
! Array  arguments with INTENT(OUT):
REAL, INTENT(OUT) :: noy(nolat,nolev)             ! 2D field interpolated in time
REAL, INTENT(OUT) :: o3(nolat,nolev)              ! 2D field interpolated in time
REAL, INTENT(OUT) :: ch4(nolat,nolev)             ! 2D field interpolated in time
!
! Local constants
INTEGER, PARAMETER :: n_2d = 74  ! no of 2d fields in file
!
! Local variables
LOGICAL, SAVE :: firstcall = .TRUE.
INTEGER             :: ipos          ! Position in 2D file
INTEGER             :: info          ! Tag used in communication
INTEGER             :: i,j,k,ifile   ! Loop variables
INTEGER             :: unit_nm(3)    ! Unit numbers
INTEGER             :: iostatus      ! output code from file open
INTEGER             :: errcode       ! Variable passed to ereport
REAL                :: fpos          ! Decides which 2D field to interpolate
REAL                :: delpos
CHARACTER (LEN=filenamelength) :: filename(3) ! Full file name incl dir pathname
CHARACTER (LEN=filenamelength) :: filenm(3)     ! File name
CHARACTER (LEN=errormessagelength) :: cmessage      ! Error message
CHARACTER (LEN=errormessagelength) :: iomessage      ! IO Error message
!
! 2D data. Saved so that we can read input files only once
!
REAL, SAVE :: ch4_2d(nolat,nolev, n_2d)! 2D data for ch4.
REAL, SAVE :: o3_2d (nolat,nolev, n_2d)! 2D data for o3
REAL, SAVE :: noy_2d(nolat,nolev, n_2d)! 2D data for noy
!
REAL       :: noya(nolat,nolev)        ! 2D field straddling day no
REAL       :: noyb(nolat,nolev)        ! 2D field straddling day no
REAL       :: o3a(nolat,nolev)         ! 2D field straddling day no
REAL       :: o3b(nolat,nolev)         ! 2D field straddling day no
REAL       :: ch4a(nolat,nolev)        ! 2D field straddling day no
REAL       :: ch4b(nolat,nolev)        ! 2D field straddling day no
!
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_2D_BC_READ_INTERP'

! End of header
!
! **********************************************************************
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! On first call to subroutine read 2D data on PE0 then broadcast
!
IF (mype == 0 .AND. firstcall) THEN

  filenm(1) = '/ch4_topbound.dat'
  filenm(2) = '/o3_topbound.dat'
  filenm(3) = '/noy_topbound.dat'

  DO ifile=1,3
    filename(ifile) = TRIM(strat2d_dir)//filenm(ifile)
    CALL assign_file_unit(&
        filename(ifile), unit_nm(ifile),handler="fortran")
    OPEN(unit_nm(ifile),FILE=filename(ifile), ACTION='READ',            &
                        IOSTAT=iostatus, IOMSG=iomessage)
    IF (iostatus /= 0) THEN
      cmessage = ' Error opening file:'//filename(ifile) // ' :'//      &
                   newline // TRIM(iomessage)
      errcode = ifile
      CALL ereport('ukca_2d_bc_read_interp',errcode,cmessage)
    END IF
  END DO

  !
  ! Read in data
  !
  DO k = 1,n_2d
    DO j = 1,nolev
      READ(unit_nm(1),*) (ch4a (i,j),i=nolat,1,-1)
      READ(unit_nm(2),*) (o3a  (i,j),i=nolat,1,-1)
      READ(unit_nm(3),*) (noya (i,j),i=nolat,1,-1)
    END DO
    ch4_2d(:,:,k) = ch4a(:,:)
    o3_2d (:,:,k) = o3a (:,:)
    noy_2d(:,:,k) = noya(:,:)
  END DO

  DO ifile=1,3
    CLOSE(unit_nm(ifile))
    CALL release_file_unit(unit_nm(ifile),handler="fortran")
  END DO

END IF         ! end if (mype.eq.0 .AND. firstcall)

IF (firstcall) THEN
  CALL gc_rbcast(1,nolat*nolev*n_2d,0,nproc,info,ch4_2d)
  CALL gc_rbcast(2,nolat*nolev*n_2d,0,nproc,info,o3_2d)
  CALL gc_rbcast(3,nolat*nolev*n_2d,0,nproc,info,noy_2d)
END IF

!     Work out position for current timestep

fpos  = idofy/5.0 + 1.0
ipos  = NINT(fpos)
IF ((fpos-ipos*1.0) < 0.0) ipos = ipos-1
delpos = fpos - ipos*1.0

!     Reset IPOS if nearest day=365 or zero

IF (ipos == n_2d) ipos = 1

ch4a(:,:) = ch4_2d(:,:,ipos)
o3a (:,:) = o3_2d (:,:,ipos)
noya(:,:) = noy_2d(:,:,ipos)

ch4b(:,:) = ch4_2d(:,:,ipos+1)
o3b (:,:) = o3_2d (:,:,ipos+1)
noyb(:,:) = noy_2d(:,:,ipos+1)

!     Interpolate in time

DO j = 1, nolev
  DO i = 1, nolat
    ch4(i,j) = (ch4b(i,j)-ch4a(i,j))*delpos + ch4a(i,j)
    o3 (i,j) = (o3b (i,j)-o3a (i,j))*delpos + o3a (i,j)
    noy(i,j) = (noyb(i,j)-noya(i,j))*delpos + noya(i,j)
  END DO
END DO

firstcall = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE ukca_2d_bc_read_interp
END MODULE ukca_2d_bc_read_interp_mod

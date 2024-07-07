! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:
SUBROUTINE anc_head(pp_int,pp_real,rows,columns,fieldsize,nfields,    &
     field_types,n_times,nlevels,n_freq_waves,n_dir_waves,no_cmp,     &
     len1_levdepc,len2_levdepc,len1_rowdepc,len2_rowdepc,len1_coldepc,&
     len2_coldepc,len1_flddepc,len2_flddepc,len_extcnst,len_cfi,      &
     tracer_grid,add_wrap_pts,periodic,single_time,ibm_to_cray,       &
     compress,levdepc,rowdepc,coldepc,flddepc,extcnst,wave,           &
     lfh,fixhd,len_intc,int_const,len_realc,real_const,icode)

USE conversions_mod, ONLY: pi
USE errormessagelength_mod, ONLY: errormessagelength
USE check_iostat_mod
USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage

USE date_conversions_mod, ONLY: date31, date13
USE missing_data_mod, ONLY: rmdi,imdi
IMPLICIT NONE
!
! Description:
!   Creates the dump/ancillary file header.
!   (Fixed length header,integer constants and real constants)
!
! Method:
!    Uses input arguments and namelist variables to appropriately
!    set the fixed length header according to UMDP F3.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
! Code Description:
!   Language: Fortran 90
!   This code is written to UMDP3 standards.

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: columns
INTEGER, INTENT(IN) :: fieldsize
INTEGER, INTENT(IN) :: nfields
INTEGER, INTENT(IN) :: field_types
INTEGER, INTENT(IN) :: n_times
INTEGER, INTENT(IN) :: nlevels
INTEGER, INTENT(IN) :: n_freq_waves
INTEGER, INTENT(IN) :: n_dir_waves
INTEGER, INTENT(IN) :: no_cmp      ! no. of total compressed points
                                   ! in compressed array

INTEGER, INTENT(IN) :: len1_levdepc
INTEGER, INTENT(IN) :: len2_levdepc
INTEGER, INTENT(IN) :: len1_rowdepc
INTEGER, INTENT(IN) :: len2_rowdepc
INTEGER, INTENT(IN) :: len1_coldepc
INTEGER, INTENT(IN) :: len2_coldepc
INTEGER, INTENT(IN) :: len1_flddepc
INTEGER, INTENT(IN) :: len2_flddepc
INTEGER, INTENT(IN) :: len_extcnst

INTEGER, INTENT(IN) :: lfh
INTEGER, INTENT(IN) :: len_intc
INTEGER, INTENT(IN) :: len_realc

!   logical choices (IN)
LOGICAL, INTENT(IN) :: tracer_grid
LOGICAL, INTENT(IN) :: add_wrap_pts
LOGICAL, INTENT(IN) :: periodic
LOGICAL, INTENT(IN) :: single_time
LOGICAL, INTENT(IN) :: ibm_to_cray
LOGICAL, INTENT(IN) :: compress
LOGICAL, INTENT(IN) :: levdepc
LOGICAL, INTENT(IN) :: rowdepc
LOGICAL, INTENT(IN) :: coldepc
LOGICAL, INTENT(IN) :: flddepc
LOGICAL, INTENT(IN) :: extcnst

LOGICAL, INTENT(IN) :: wave   ! T for creating wave dump

!   Array  arguments with intent(in):
INTEGER, INTENT(IN) :: pp_int(45)
INTEGER, INTENT(IN) :: len_cfi(3)
REAL,    INTENT(IN) :: pp_real(19)

!   Array  arguments with intent(out):
INTEGER, INTENT(OUT) :: fixhd(lfh)
INTEGER, INTENT(OUT) :: int_const(len_intc)
REAL,    INTENT(OUT) :: real_const(len_realc)

!   Arguments with intent(inout):
INTEGER, INTENT(INOUT) :: icode

!   Local Scalars
INTEGER :: ipos
INTEGER :: i,j
INTEGER :: errorstatus

INTEGER :: fvhh,fvdd,fvmm,fvyy   ! hour,day,month,year first validity
                                 ! time
INTEGER :: lvhh,lvdd,lvmm,lvyy   ! hour,day,month,year last validity
                                 ! time
INTEGER :: ivhh,ivdd,ivmm,ivyy   ! hour,day,month,year interval

LOGICAL :: year360               ! true for 360-day calendar
LOGICAL :: year_unspec           ! true if calendar is unspecified in fixed
                                 !   length header.
                                 !   Note: this should only be used for
                                 !   single time ancillaries or
                                 !   time varying ancillaries
                                 !   without sub-monthly intervals.

LOGICAL :: l_first_vt
LOGICAL :: l_last_vt

INTEGER :: cdays ! century days
INTEGER :: chours ! century hours
INTEGER :: new_cdays ! new century days
INTEGER :: new_chours! new century hours
INTEGER :: ihr,idy,imn

CHARACTER(LEN=errormessagelength) :: iomessage
! Function & Subroutine calls:

INTEGER :: find_namelist

! End of header

NAMELIST /first_vt/ fvhh,fvdd,fvmm,fvyy
NAMELIST /last_vt/  lvhh,lvdd,lvmm,lvyy
NAMELIST /interval/ year360,year_unspec,ivhh,ivdd,ivmm,ivyy

!  1. Set fixed header

!  1.0 Initialise to missing data
! (* set dimensions of all arrays to 1 *)

! DEPENDS ON: init_flh
CALL init_flh (fixhd,lfh)

!  1.1 First 9 elements in header DEFAULTS

fixhd(2)=2               ! indicator for the ocean
fixhd(3)=4               ! depth coordinates
fixhd(4)=0               ! global grid
fixhd(5)=4               ! ancillary fields dataset
fixhd(8)=1               ! Calendar indicator (Gregorian)
fixhd(9)=1               ! Indicator for grid staggering

!  1.2 Set dates

IF (periodic) THEN
  fixhd(10)=2
ELSE IF (single_time) THEN
  fixhd(10)=0
ELSE
  fixhd(10)=1
END IF

fixhd(12)=401                      ! UM Version number

! (* first validity time *)
! (* fixhd(21-27)       *)

fvhh = 0
fvdd = 0
fvmm = 0
fvyy = 0

REWIND(5)
! DEPENDS ON: find_namelist
i=find_namelist(5,"first_vt")

IF (i == 0) THEN
  READ (UNIT=5, NML=first_vt, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist FIRST_VT", iomessage)
ELSE
  WRITE(umMessage,'(a)')'Cannot find namelist FIRST_VT'
  CALL umPrint(umMessage,src='anc_head')
END IF
CALL print_nlist_first_vt()

!     Test if first VT has been provided in namelist
l_first_vt = .NOT.                                                &
     (fvhh == 0 .AND. fvdd == 0 .AND. fvmm == 0 .AND. fvyy == 0)

IF (l_first_vt) THEN  !  First VT given in namelist

  fixhd(21) = fvyy
  fixhd(22) = fvmm
  fixhd(23) = fvdd
  fixhd(24) = fvhh
  fixhd(25) = 0
  fixhd(26) = 0
  fixhd(27) = 0

ELSE IF (single_time) THEN

  fixhd(21:27) = 0    !  Set First VT to zero

ELSE  !  Get first VT from first PP Header

  fixhd(21) = pp_int(1)
  fixhd(22) = pp_int(2)
  fixhd(23) = pp_int(3)
  fixhd(24) = pp_int(4)
  fixhd(25) = pp_int(5)
  fixhd(26) = 0
  fixhd(27) = pp_int(6)

END IF

year360=.FALSE.
year_unspec=.FALSE.
ivhh = 0
ivdd = 0
ivmm = 0
ivyy = 0

! (* interval is read from namelist*)
! (* fixhd(35-41)       *)

REWIND(5)
! DEPENDS ON: find_namelist
i=find_namelist(5,"interval")

IF (i == 0) THEN
  READ (UNIT=5, NML=interval, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist INTERVAL", iomessage)
ELSE
  WRITE(umMessage,'(a)')'Cannot find namelist INTERVAL'
  CALL umPrint(umMessage,src='anc_head')
END IF

CALL print_nlist_interval()

fixhd(35) = ivyy
fixhd(36) = ivmm
fixhd(37) = ivdd
fixhd(38) = ivhh
fixhd(39) = 0
fixhd(40) = 0
fixhd(41) = 0

! (* last validity time *)
! (* fixhd(28-34)       *)

lvhh = 0
lvdd = 0
lvmm = 0
lvyy = 0

REWIND(5)
! DEPENDS ON: find_namelist
i=find_namelist(5,"last_vt")

IF (i == 0) THEN
  READ (UNIT=5, NML=last_vt, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist LAST_VT", iomessage)
ELSE
  WRITE(umMessage,'(a)')'Cannot find namelist LAST_VT'
  CALL umPrint(umMessage,src='anc_head')
END IF

CALL print_nlist_last_vt()

!     Test IF last VT has been provided in namelist
l_last_vt = .NOT.                                                 &
     (lvhh == 0 .AND. lvdd == 0 .AND. lvmm == 0 .AND. lvyy == 0)

IF (year360) THEN

  !   Set to 360 day calendar only
  fixhd(8) = 2

  IF (l_last_vt) THEN  ! Last VT given in namelist

    fixhd(28) = lvyy
    fixhd(29) = lvmm
    fixhd(30) = lvdd
    fixhd(31) = lvhh
    fixhd(32) = 0
    fixhd(33) = 0
    fixhd(34) = 0

  ELSE IF (single_time) THEN

    fixhd(28:34) = 0    !  Set Last VT to zero

  ELSE  !  calculate last VT from first VT and Interval

    fixhd(33)=fixhd(26)  ! seconds
    fixhd(32)=fixhd(25)  ! minutes

    ihr=fixhd(24)+ivhh*(n_times-1)
    fixhd(31)=MOD(ihr,24)

    idy=fixhd(23)+ivdd*(n_times-1)+ihr/24
    fixhd(30)=MOD(idy-1,30)+1

    imn=fixhd(22)+ivmm*(n_times-1)+(idy-1)/30
    fixhd(29)=MOD(imn-1,12)+1

    fixhd(28)=fixhd(21)+ivyy*(n_times-1)+(imn-1)/12

    fixhd(34)=(fixhd(29)-1)*30+fixhd(30)

  END IF

ELSE   !  Gregorian calander files

  !   Set to Gregorian calendar only
  fixhd(8) = 1    ! Gregorian calendar only

  IF (l_last_vt) THEN  !  Last VT given in namelist

    fixhd(28) = lvyy
    fixhd(29) = lvmm
    fixhd(30) = lvdd
    fixhd(31) = lvhh
    fixhd(32) = 0
    fixhd(33) = 0
    fixhd(34) = 0

  ELSE IF (single_time) THEN

    fixhd(28:34) = 0    !  Set Last VT to zero

  ELSE  !  calculate last VT from first VT and Interval

    !         Check First VT and Interval first
    !         First VT is OK if FIXHD(21,22,23) are all set.
    !         Interval is OK if only IVHH and/or IVDD are used.
    IF (fixhd(21) <= 0 .OR. fixhd(22) <= 0 .OR.                   &
         fixhd(23) <= 0 .OR. ivmm >  0 .OR. lvyy >  0) THEN
      WRITE(umMessage,'(a)') ' '
      CALL umPrint(umMessage,src='anc_head')
      WRITE(umMessage,'(a)') ' ERROR in ANC_HEAD. Last Validity Time ?'
      CALL umPrint(umMessage,src='anc_head')
      WRITE(umMessage,'(a)') ' Last VT cant be calculated from first VT.'
      CALL umPrint(umMessage,src='anc_head')
      WRITE(umMessage,'(a)') ' Rerun job with last VT in LAST_VT namelist.'
      CALL umPrint(umMessage,src='anc_head')
      WRITE(umMessage,'(a)') ' '
      CALL umPrint(umMessage,src='anc_head')

      icode = 1
      GO TO 9999   !  Return
    END IF

    ! (* calculate century day and hour for first validity time)
    CALL date31(fixhd(23),fixhd(22),fixhd(21),cdays)
    chours=(cdays-1)*24+fixhd(24)

    ! (* add time interval)
    new_chours=chours+(n_times-1)*(ivhh+ivdd*24)

    ! (* convert to new century day)
    new_cdays=1+new_chours/24

    ! (* convert to actual date)
    CALL date13(new_cdays,fixhd(30),fixhd(29),fixhd(28))
    fixhd(31)=new_chours-24*(new_cdays-1)

    fixhd(32) = 0
    fixhd(33) = 0
    fixhd(34) = 0

  END IF

END IF

IF (year_unspec) THEN ! Ancillary file designed to work in either
                      ! calendar. Only allow if single time or
                      ! interval divides into months.

  IF (single_time .OR. (ivhh == 0 .AND. ivdd == 0) ) THEN
    fixhd(8) = imdi
  ELSE
    WRITE (umMessage,'(a)') ' '
    CALL umPrint(umMessage,src='anc_head')
    WRITE (umMessage,'(a)') ' ERROR in ANC_HEAD.'
    CALL umPrint(umMessage,src='anc_head')
    WRITE (umMessage,'(a)') ' Choice to use unspecified calendar type'
    CALL umPrint(umMessage,src='anc_head')
    WRITE (umMessage,'(a)') ' (year_unspec = .true.) is incompatible'
    CALL umPrint(umMessage,src='anc_head')
    WRITE (umMessage,'(a)') ' with interval choice. This should not be'
    CALL umPrint(umMessage,src='anc_head')
    WRITE (umMessage,'(a)') ' used with sub-monthly intervals.'
    CALL umPrint(umMessage,src='anc_head')
    WRITE (umMessage,'(a)') ' '
    CALL umPrint(umMessage,src='anc_head')
    icode = 1
    GO TO 9999   !  Return
  END IF
END IF

WRITE(umMessage,'('' '')')
CALL umPrint(umMessage,src='anc_head')
WRITE(umMessage,'('' Validity Times (VT) in Ancillary File.'')')
CALL umPrint(umMessage,src='anc_head')
WRITE(umMessage,'('' '')')
CALL umPrint(umMessage,src='anc_head')
CALL umPrint('               Year  Month Day Hour Min  Sec DayNo ', &
              src='anc_head')
WRITE(umMessage,'('' First VT    ='',7I5)')(fixhd(i),i=21,27)
CALL umPrint(umMessage,src='anc_head')
WRITE(umMessage,'('' Last  VT    ='',7I5)')(fixhd(i),i=28,34)
CALL umPrint(umMessage,src='anc_head')
WRITE(umMessage,'('' VT Interval ='',7I5)')(fixhd(i),i=35,41)
CALL umPrint(umMessage,src='anc_head')
WRITE(umMessage,'('' '')')
CALL umPrint(umMessage,src='anc_head')

!  1.3 Set pointers to starts of sections in ancillary file

ipos = lfh + 1   ! position of start of integer consts

! (* integer constants location *)
fixhd(100)= ipos
fixhd(101)=len_intc
ipos = ipos + len_intc

! (* real constants location *)
fixhd(105)=ipos
fixhd(106)=len_realc
ipos = ipos + len_realc

! (* levels dependent constants*)

IF (levdepc) THEN
  fixhd(110) = ipos
  fixhd(111) = len1_levdepc
  fixhd(112) = len2_levdepc
  ipos = ipos + len1_levdepc*len2_levdepc
END IF

IF (rowdepc) THEN
  fixhd(115) = ipos
  fixhd(116) = len1_rowdepc
  fixhd(117) = len2_rowdepc
  ipos = ipos + len1_rowdepc*len2_rowdepc
END IF

IF (coldepc) THEN
  fixhd(120) = ipos
  fixhd(121) = len1_coldepc
  fixhd(122) = len2_coldepc
  ipos = ipos + len1_coldepc*len2_coldepc
END IF

IF (flddepc) THEN
  fixhd(125) = ipos
  fixhd(126) = len1_flddepc
  fixhd(127) = len2_flddepc
  ipos = ipos + len1_flddepc*len2_flddepc
END IF

IF (extcnst) THEN
  fixhd(130) = ipos
  fixhd(131) = len_extcnst
  ipos = ipos + len_extcnst
END IF

! (* compressed field indices *)
IF (compress .AND. .NOT. wave) THEN
  fixhd(140) = ipos
  fixhd(141) = len_cfi(1)
  ipos = ipos + len_cfi(1)

  fixhd(142) = ipos
  fixhd(143) = len_cfi(2)
  ipos = ipos + len_cfi(2)

  fixhd(144) = ipos
  fixhd(145) = len_cfi(3)
  ipos = ipos + len_cfi(3)
END IF

! (* location of lookup table *)
fixhd(150)=ipos
fixhd(151)=64
fixhd(152)=nfields

! for wave dump - number of fields in dump *
fixhd(153)=nfields

! (* location of data *)
ipos=ipos+64*nfields
fixhd(160)=ipos

!   fixhd(161) is set in ancfld (it is updated after each field
!   is read )

!  2. Set Integer Constants

DO j=1,len_intc
  int_const(j)=imdi
END DO

int_const(3)=n_times

IF (add_wrap_pts) THEN
  int_const(6)=columns+2
ELSE
  int_const(6)=columns
END IF

! When the UM reads the ancillary files (see RPANCO1A) it checks that
! the number of rows in the model tracer grid (JMT) matches the number
! of rows declared in the integer consts; the number of rows in the
! velocity grid is one less  than that in the tracer grid.
!
IF (tracer_grid) THEN
  int_const(7)=rows
ELSE
  int_const(7)=rows+1
END IF

int_const(8) = nlevels

!CMHWAVES
IF (wave) THEN
  int_const(8) = n_freq_waves
  int_const(9) = n_dir_waves
  int_const(10)= fieldsize
END IF
!CMHWAVES

IF (compress) THEN
  int_const(11) = no_cmp
END IF

int_const(15)=field_types

!  3. Set real constants

DO j=1,len_realc
  real_const(j)=rmdi
END DO

! (* grid spacing)
real_const(1)=pp_real(17)
real_const(2)=ABS(pp_real(15))

! (* lat of first row (3) and long of first point on row (4)
IF (tracer_grid) THEN
  real_const(3)=pp_real(14)+pp_real(15)
  real_const(4)=pp_real(16)+pp_real(17)
ELSE
  real_const(3)=pp_real(14)+0.5*pp_real(15)
  real_const(4)=pp_real(16)+0.5*pp_real(17)
END IF

! (* test value of the start longitude *)
IF (real_const(4) <  0.0) THEN
  real_const(4)=real_const(4)+360.0
ELSE IF (real_const(4) >= 360.0) THEN
  real_const(4)=real_const(4)-360.0
END IF

! (* lat and long of pseudo north pole)
real_const(5)=pp_real(11)
real_const(6)=pp_real(12)

! WAVES
! direction increment
IF (wave) THEN
  real_const(13)=2.0*pi/n_dir_waves
END IF

9999 CONTINUE
RETURN

CONTAINS

SUBROUTINE print_nlist_interval()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist interval', &
    src='anc_head')

WRITE(lineBuffer,*)' year360 = ',year360
CALL umPrint(lineBuffer,src='anc_head')
WRITE(lineBuffer,*)' ivhh = ',ivhh
CALL umPrint(lineBuffer,src='anc_head')
WRITE(lineBuffer,*)' ivdd = ',ivdd
CALL umPrint(lineBuffer,src='anc_head')
WRITE(lineBuffer,*)' ivmm = ',ivmm
CALL umPrint(lineBuffer,src='anc_head')
WRITE(lineBuffer,*)' ivyy = ',ivyy
CALL umPrint(lineBuffer,src='anc_head')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='anc_head')

END SUBROUTINE print_nlist_interval

SUBROUTINE print_nlist_first_vt()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist first_vt', &
    src='anc_head')

WRITE(lineBuffer,*)' fvhh = ',fvhh
CALL umPrint(lineBuffer,src='anc_head')
WRITE(lineBuffer,*)' fvdd = ',fvdd
CALL umPrint(lineBuffer,src='anc_head')
WRITE(lineBuffer,*)' fvmm = ',fvmm
CALL umPrint(lineBuffer,src='anc_head')
WRITE(lineBuffer,*)' fvyy = ',fvyy
CALL umPrint(lineBuffer,src='anc_head')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='anc_head')

END SUBROUTINE print_nlist_first_vt

SUBROUTINE print_nlist_last_vt()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist last_vt', &
    src='anc_head')

WRITE(lineBuffer,*)' lvhh = ',lvhh
CALL umPrint(lineBuffer,src='anc_head')
WRITE(lineBuffer,*)' lvdd = ',lvdd
CALL umPrint(lineBuffer,src='anc_head')
WRITE(lineBuffer,*)' lvmm = ',lvmm
CALL umPrint(lineBuffer,src='anc_head')
WRITE(lineBuffer,*)' lvyy = ',lvyy
CALL umPrint(lineBuffer,src='anc_head')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='anc_head')

END SUBROUTINE print_nlist_last_vt

END SUBROUTINE anc_head

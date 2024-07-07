! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  read next pp field from 32 bit pp file

MODULE read_next_pp_field_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READ_NEXT_PP_FIELD_MOD'

CONTAINS

SUBROUTINE read_next_pp_field(unit_in,stashcode,read_code,dmin,iihead,rrhead,&
                              field)

USE errormessagelength_mod, ONLY: errormessagelength
USE crmstyle_cntl_mod, ONLY:                              &
  mlevs, model_levels, in_cols, in_rows

USE pp_field32_mod

USE word_sizes_mod, ONLY: FourByteReal, FourByteInt, wp

USE missing_data_mod, ONLY: rmdi

USE wgdos_packing_mod, ONLY: wgdos_expand_field

USE ereport_mod, ONLY: ereport, ereport_finalise

USE umPrintMgr, ONLY: umprint, ummessage, newline
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE packing_codes_mod, ONLY: PC_WGDOS_Packing

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!  Read in data from a sequential pp file which has been converted to 32 bit
! numbers. The fact that it has been converted to 32 bit numbers complicates
! the process as the supercomputer by default is expecting 64 bit numbers.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utilty - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN)  ::  &
  unit_in                  ! Number of columns

INTEGER, INTENT(OUT) ::  &
  stashcode              & ! stashcode for field
 ,read_code              & ! error_code
 ,iihead(45)               ! integer pp header

REAL, INTENT(OUT) ::     &
  rrhead(19)               ! real pp header

REAL(wp), INTENT(OUT) ::  &
  field(in_cols,in_rows,model_levels)

!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, k, ij              ! loop counters
INTEGER ::               &
  idum1                  & ! dummy
 ,idum2                  & ! dummy
 ,icode                  &
 ,len_full               &
 ,len_in                 &
 ,lev_type               &
 ,lev_num                &
 ,pack_code              &
 ,error_code               ! error code

INTEGER ::               &
  ix                     & ! cols input field
 ,iy                       ! rows input field
INTEGER ::               &
  imin1                  & ! minutes for validity
 ,imin2                  & ! minutes for start time
 ,dmin                   & ! dmin
 ,steps_per_hour           ! number of forecast output per hour

INTEGER(KIND=FourByteInt) :: ihead(45)    ! 32 bit input
REAL(KIND=FourByteReal):: rhead(19)       ! 32 bit input

REAL ::                   &
  temp1(in_cols*in_rows)  & ! work array
 ,temp2(in_cols*in_rows)

INTEGER ::                &
  int_temp1(in_cols*in_rows) ! Integer array to hold integer version of temp1

REAL(KIND=FourByteReal):: temp3(in_cols*in_rows)       ! 32 bit input

TYPE(pp_field32) :: fin    ! pp field

CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=*), PARAMETER :: RoutineName = "READ_NEXT_PP_FIELD"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

error_code=0
read_code=0
! Read next pp header in file

READ(unit_in,IOSTAT=error_code) ihead,rhead

! Want to copy 32 bit to 64 bit

DO i=1,45
  iihead(i) = ihead(i)
END DO
DO i=1,19
  rrhead(i) = rhead(i)
END DO

fin%ipphead=TRANSFER(ihead,fin%ipphead)
fin%rpphead=TRANSFER(rhead,fin%rpphead)

stashcode=fin%ipphead%stashcode
lev_type=fin%ipphead%vc
lev_num=fin%ipphead%met08_lev

imin1=fin%ipphead%minute
imin2=fin%ipphead%mind
dmin=imin1-imin2
steps_per_hour = 60/dmin

len_in=(fin%ipphead%lrec)/2    ! length of data in 64 bit words
len_full=fin%ipphead%npt * fin%ipphead%row
pack_code=fin%ipphead%ipack    ! packing code

IF (fin%ipphead%ipack == PC_WGDOS_Packing) THEN
  ! work field length and read in field
  ! reads into 64 bit array but for correct length of data?
  READ(unit_in,IOSTAT=error_code, IOMSG=iomessage) (temp1(i),i=1,len_in)

  IF (error_code == 0) THEN
    ! Unpack data
    int_temp1 = TRANSFER(temp1, int_temp1)
    CALL wgdos_expand_field(temp2,len_full,int_temp1,len_in,idum1,ix,iy,  &
                            idum2,rmdi,stashcode,icode)
  ELSE
    WRITE(umMessage,'(A,I10,A)')                                             &
      'READ_NEXT_PP_FIELD: Error reading field ', error_code,       newline//&
      'IoMsg: '//TRIM(iomessage)
    CALL umPrint(umMessage,src=RoutineName)
  END IF
ELSE
  ! 32 bit copy to 64 bit
  READ(unit_in,IOSTAT=error_code,IOMSG=iomessage)                            &
                                            (temp3(i),i=1,fin%ipphead%lrec)
  IF (error_code == 0) THEN
    ix = fin%ipphead%npt
    iy = fin%ipphead%row
    DO i = 1, fin%ipphead%lrec
      temp2(i) = temp3(i)
    END DO
  ELSE
    WRITE(umMessage,'(A,I10,A)')                                             &
      'READ_NEXT_PP_FIELD: Error reading field ', error_code,       newline//&
      'IoMsg: '//TRIM(iomessage)
    CALL umPrint(umMessage,src=RoutineName)
    read_code = error_code
  END IF
END IF

! Copy to output field
k=1
DO j=1,iy
  DO i=1,ix
    ij=(j-1)*ix+i
    field(i,j,k) = temp2(ij)
  END DO
END DO

! if the level type for model levels read the other fields

IF (lev_type == 65) THEN   ! model level data hybrid heights
  ! First check level already read in was number one

  IF (lev_num == 1) THEN

    ! read other levels
    DO k=2,mlevs ! Only want mlevs even if more

      READ(unit_in,IOSTAT=error_code) ihead,rhead
      fin%ipphead=TRANSFER(ihead,fin%ipphead)
      len_in=(fin%ipphead%lrec)/2    ! length of data in 64 bit words
      IF (fin%ipphead%ipack == PC_WGDOS_Packing) THEN
        ! work field length and read in field

        READ(unit_in,IOSTAT=error_code, IOMSG=iomessage) (temp1(i),i=1,len_in)
        IF (error_code == 0) THEN
          ! Unpack data
          int_temp1 = TRANSFER(temp1, int_temp1)
          CALL wgdos_expand_field(temp2,len_full,int_temp1,len_in,idum1,  &
                                  ix,iy,idum2,rmdi,                       &
                                  INT(fin%ipphead%stashcode),icode)
        ELSE
          WRITE(umMessage,'(A,2I10,A)')                                       &
            'READ_NEXT_PP_FIELD: Error reading field ',k,error_code, newline//&
            'IoMsg: ' // TRIM(iomessage)
          CALL umPrint(umMessage,src=RoutineName)
          read_code = error_code
          EXIT
        END IF
      ELSE
        ! Read into 32 bit array and then copy to 64 bit real
        READ(unit_in,IOSTAT=error_code, IOMSG=iomessage)                    &
             (temp3(i),i=1,fin%ipphead%lrec)
        IF (error_code == 0) THEN
          DO i = 1, fin%ipphead%lrec
            temp2(i) = temp3(i)
          END DO
        ELSE
          WRITE(umMessage,'(A,2I10,A)')                                       &
            'READ_NEXT_PP_FIELD: Error reading field ',k,error_code, newline//&
            'IoMsg: ' // TRIM(iomessage)
          CALL umPrint(umMessage,src=RoutineName)
          read_code = error_code
          EXIT
        END IF
      END IF

      ! Copy to output field
      DO j=1,iy
        DO i=1,ix
          ij=(j-1)*ix+i
          field(i,j,k) = temp2(ij)
        END DO
      END DO

    END DO

  ELSE IF (lev_num == mlevs+1) THEN
      ! read in but don't store i.e. Skip fields as don't want them
    DO k=mlevs+2,model_levels
      READ(unit_in,IOSTAT=error_code,IOMSG=iomessage) ihead,rhead
      fin%ipphead=TRANSFER(ihead,fin%ipphead)
      len_in=(fin%ipphead%lrec)/2    ! length of data in 64 bit words

      ! work field length and read in field

      READ(unit_in,IOSTAT=error_code) (temp1(i),i=1,len_in)
      IF (error_code /= 0) THEN
        WRITE(umMessage,'(A,2I10,A)')                                       &
          'READ_NEXT_PP_FIELD: Error reading field ',k,error_code, newline//&
          'IoMsg: ' // TRIM(iomessage)
        CALL umPrint(umMessage,src=RoutineName)
        read_code = error_code
        EXIT
      END IF

    END DO

  END IF ! test on level number


END IF     ! hybrid height
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE read_next_pp_field

END MODULE read_next_pp_field_mod

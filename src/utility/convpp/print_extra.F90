! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
SUBROUTINE print_extra(lookup,d1,ns_pts,k)

USE lookup_addresses

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE


! Purpose: To print out extra data from D1
!
! Subroutine Arguments:

INTEGER ::                                                        &
 lookup(64,*)                                                     &
                    ! IN  :lookup table
,ns_pts                                                           &
                    ! IN  :number of extra data values printed
,k                  ! IN  :field number

REAL ::                                                           &
 d1(*)              ! IN  :D1 array

INTEGER ::                                                        &
 tot_values                                                       &
,iextraw                                                          &
,addr                                                             &
,code                                                             &
,data_values                                                      &
,max_print                                                        &
,ew_pts                                                           &
                    ! No of points across to print
,ew_print                                                         &
,start_off                                                        &
,n_vect                                                           &
                    ! No of vector in extra data
,max_inthd                                                        &
,max_data                                                         &
                    ! Max possible size of extra data
,i,j,l                                                            &
                    !Loop counter
,icode

PARAMETER( ew_pts=5,                                              &
           max_inthd=30)


INTEGER :: ind   ! Local index


INTEGER ::                                                        &
 x_inthd(max_inthd)

REAL ::                                                           &
 DATA(lookup(lbext,k))                                            &
,data_line(ew_pts)

CHARACTER(LEN=12) ::                                                   &
 dash(ew_pts)       !Stores dashed lines

CHARACTER(LEN=1) ::                                                    &
 BLANK

PARAMETER(BLANK=' ')

CHARACTER(LEN=errormessagelength) ::                                   &
 cmessage


! Extract information from extra data
! Put the vector header of each extra data vector
! into array X_INTHD
! Put the first NS_PTS values of each extra data vector
! into array X_DATA

! A extra data vector comprises of:
!
!  header value1 value2 value3 value4 ....
!
! |______|______|______|______|______|_____|
!
!
! The extra data vector(s) are attached to the end of
! each set of field data.  For a detail description,
! see documentation F3
!

iextraw=lookup(lbext,k)
tot_values=lookup(lblrec,k)
addr=tot_values-iextraw+1
l=0
ind=0
start_off=1
! While our address is short of the final record
DO WHILE (addr <  tot_values)
  l=l+1
  ! Decode the extra data vector descriptor
  code=TRANSFER(d1(addr), code)
  x_inthd(l)=code
  ! DEPENDS ON: check_extra
  CALL check_extra(code,data_values,icode,cmessage)
  icode=0
  ! The max number of points to print out - don't exceed
  ! the number of elements in each vector
  max_print=MIN(data_values,ns_pts)

  ! Fill in actual values
  i = 0
  DO j=addr+1,addr+data_values
    i = i + 1
    IF (i <= max_print) THEN
      ind=ind+1
      DATA(ind)=d1(j)
    END IF
  END DO
  addr=addr+data_values+1   ! INCREMENT ADDRESS
                            ! by no of data elements + 1
                            ! for the vector descriptor.
  start_off=start_off+max_print
END DO
n_vect=l

! Set up print format

DO j=1,ew_pts
  dash(j)='------------'
END DO

ew_print=MIN(n_vect,ew_pts)

! Print out extra data

! Print 1st row - integer extra data header

WRITE(umMessage,'(14X,A25,I3,A8)')                                       &
   'Extra Data Values (first ',max_print,' values)'
CALL umPrint(umMessage,src='print_extra')
WRITE(umMessage,'(14X,9A12)')(dash(j),j=1,ew_print)
CALL umPrint(umMessage,src='print_extra')
WRITE(umMessage,'('' FIELD NO'',I4,'':'',9(I10,2X))')                    &
     k,(x_inthd(j),j=1,ew_print)
CALL umPrint(umMessage,src='print_extra')
WRITE(umMessage,'(14X,9A12)')(dash(j),j=1,ew_print)
CALL umPrint(umMessage,src='print_extra')

! Print remaining rows - extra data values
! For now, just printing out max of 5 vectors

DO start_off=1,max_print
  DO i=1,ew_print
    data_line(i)=DATA(start_off+(i-1)*max_print)
  END DO
  WRITE(umMessage,'(14X,A1,9(F9.3,3X))')                                 &
  BLANK,(data_line(j),j=1,ew_print)
  CALL umPrint(umMessage,src='print_extra')
END DO

WRITE(umMessage,'(14X,9A12)')(dash(j),j=1,ew_print)
CALL umPrint(umMessage,src='print_extra')
WRITE(umMessage,'('' '')')
CALL umPrint(umMessage,src='print_extra')

IF (n_vect >  ew_pts) THEN
  WRITE(umMessage, "(A,I0)") 'NUMBER OF EXTRA DATA VECTORS: ', n_vect
  CALL umPrint(umMessage,src='print_extra')
  WRITE(umMessage, "(A)") 'ONLY THE FIRST 5 WERE PRINTED OUT.'
  CALL umPrint(umMessage,src='print_extra')
  WRITE(umMessage, "(A)") ' '
  CALL umPrint(umMessage,src='print_extra')
END IF

RETURN
END SUBROUTINE print_extra

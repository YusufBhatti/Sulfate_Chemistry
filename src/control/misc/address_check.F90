! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE address_check_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ADDRESS_CHECK_MOD'

CONTAINS
!
! Subroutine ADDRESS_CHECK --------------------------------------
!
! Purpose : Check that start addresses of fields read in agree
!           with start addresses set up by UI. Called in INITDUMP
!           if prognostic fields read in from atmos or ocean dumps.
!
! Coding Standard : UM documentation paper no. 3
!
! Documentation : None
!
!----------------------------------------------------------------
!
!   Arguments
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc
!----------------------------------------------------------------
SUBROUTINE address_check (lookup,mpp_dump_addr,mpp_dump_len,      &
                          len1_lookup,len2_lookup,                &
                          si,nitems,nsects,len_data,              &
                          icode,cmessage)


USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim
USE submodel_mod, ONLY: n_internal_model
USE umPrintMgr,   ONLY: umPrint, umMessage
USE errormessagelength_mod, ONLY: errormessagelength
USE f_type_mod, ONLY: f_type

IMPLICIT NONE


INTEGER ::                                                        &
    len1_lookup                                                   &
                        !  1st dimension of lookup table
   ,len2_lookup                                                   &
                        !  2nd dimension of lookup table
   ,lookup(len1_lookup,len2_lookup)                               &
                                     !  Lookup table
   ,mpp_dump_addr(len2_lookup)                                    &
                                ! Addresses of fields as
!                                     ! calculated in READDUMP.
         ,mpp_dump_len(len2_lookup)                                     &
                                      ! Lengths of fields as
!                                     ! calculated in READDUMP.

         ,len_data                                                      &
                              !  Expected length of data
         ,nitems                                                        &
                              !  No of stash items
         ,nsects                                                        &
                              !  No of stash sections
!  Stash item addresses
         ,si(nitems,0:nsects,n_internal_model)                          &
         ,icode               !  Return code

CHARACTER(LEN=errormessagelength) ::                                    &
          cmessage      !  Error message

!  Dynamic allocation of arrays for F_TYPE
INTEGER ::                                                        &
    pp_num   (len2_lookup)                                        &
   ,pp_len   (len2_lookup)                                        &
   ,pp_pos   (len2_lookup)                                        &
   ,pp_ls    (len2_lookup)                                        &
   ,pp_stash (len2_lookup)                                        &
   ,pp_type  (len2_lookup)

!  Local array
INTEGER :: fixhd5  !  Dummy variable until removed from F_TYPE

!  Local variables
INTEGER ::                                                        &
    address_stash                                                 &
   ,address_lookup                                                &
   ,item_code                                                     &
   ,j                                                             &
   ,LEN                                                           &
   ,n_types                                                       &
   ,sect_no                                                       &
   ,im_index                                                      &
             ! Position of int mod id in INTERNAL_MODEL_LIST
,old_stash   ! VALUE OF STASH NUMBER ON PREVIOUS ITERATION OF LOOP

LOGICAL, PARAMETER :: verbose = .FALSE. 

CHARACTER(LEN=80) :: title

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ADDRESS_CHECK'

! ---------------------------------------------------------------------
!
!   SET INITIAL VALUE OF PREVIOUS STASH NUMBER
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
old_stash = -1
!
!    Internal Structure

title = 'Prognostic fields'
! Set fixhd5 to be a dump to allow any checks in f_type to be safe.
fixhd5 = 1
CALL f_type (lookup,len2_lookup,pp_num,n_types,pp_len,            &
             pp_stash,pp_type,pp_pos,pp_ls,fixhd5,                &
             title,verbose)

DO j=1,n_types

  !       Get Stash Section no and Item Code
  item_code = MOD ( pp_stash(j),1000)
  sect_no   = (pp_stash(j)-item_code)/1000

  im_index = 1

  !       Get lookup and stash start address
  address_lookup = mpp_dump_addr(pp_pos(j))
  address_stash  = si(item_code,sect_no,im_index)

  !       Check that they match
  !
  !            CHECK THAT START ADDRESSES AGREE FOR FIRST OCCURRENCE
  !            OF A NEW STASH CODE:  FOR FIXED LENGTH FIELDS THERE
  !            IS ONLY ONE ENTRY IN THE PP_STASH ARRAY FOR EACH
  !            STASH CODE, BUT FOR PACKED FIELDS (EG OCEAN)
  !            EACH LEVEL MIGHT HAVE A DIFFERENT LENGTH AND
  !            GENERATE A NEW PP_STASH VALUE.
  !
  IF (address_stash  /=  address_lookup .AND.                    &
      old_stash  /=  pp_stash(j) ) THEN
    cmessage = 'ADDR_CHK : Mis-match in start addresses'
    WRITE(umMessage,*) ' Stash Sect No ',sect_no,' Item No ',item_code
    CALL umPrint(umMessage,src='address_check')
    WRITE(umMessage,*) ' Start Address in SI           ',address_stash
    CALL umPrint(umMessage,src='address_check')
    WRITE(umMessage,*) ' Start Address in LOOKUP Table ',address_lookup
    CALL umPrint(umMessage,src='address_check')
    WRITE(umMessage,*) ' You probably need to RECONFIGURE the start dump'
    CALL umPrint(umMessage,src='address_check')
    icode = 1
    GO TO 999   !  Return
  END IF

  !     REMEMBER CURRENT VERSION OF PP_STASH FOR NEXT TIME THRU LOOP
  old_stash = pp_stash(j)
  !
END DO

!     Check full length
LEN = 0
DO j=1,len2_lookup
  LEN = LEN + mpp_dump_len(j)
END DO

IF (LEN  /=  len_data) THEN
  cmessage = 'ADDR_CHK : Mismatch in length of data'
  WRITE(umMessage,*) ' Length according to LOOKUP table ',LEN
  CALL umPrint(umMessage,src='address_check')
  WRITE(umMessage,*) ' Length set up in D1 array        ',len_data
  CALL umPrint(umMessage,src='address_check')
  WRITE(umMessage,*) ' You probably need to RECONFIGURE the start dump'
  CALL umPrint(umMessage,src='address_check')
  icode = 2
  GO TO 999   !  Return
END IF

999  CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE address_check
END MODULE address_check_mod

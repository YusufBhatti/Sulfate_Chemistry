! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Checks the packing codes in a model dump.
!
! Subroutine Interface:

SUBROUTINE Check_Dump_Packing (                                   &
           FixHd, Len_FixHd,                                      &
           Lookup, Len1_Lookup, Len2_Lookup,                      &
           Dump_Pack, IM_Ident )


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: printstatus, prstatus_normal, newline
USE lookup_addresses
USE cppxref_mod, ONLY: ppx_dump_packing
USE ppxlook_mod, ONLY: exppxi

USE errormessagelength_mod, ONLY: errormessagelength

USE packing_codes_mod, ONLY: PC_No_Packing, PC_Cray32_Packing, PC_Native_Format

IMPLICIT NONE
!
! Description:
!   Checks consistency of packing codes in a model dump.
!
! Method:
!   1. Packing in a dump is controlled through DUMP_PACKim
!   2. The packing codes in the LookUp table (Word 21) is checked for
!      consistency with DUMP_PACKim.
!   3. If inconsistent, the packing code in LookUp table is updated
!      to match DUMP_PACKim.
!   4. Only 32 bit packing or no packing is catered for.
!      (ie. N1=0 or N1=2 , see UMDP F3)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Subroutine arguments

INTEGER :: Len_FixHd
INTEGER :: Len1_Lookup
INTEGER :: Len2_Lookup
INTEGER :: FixHd (Len_FixHd)
INTEGER :: LookUp(Len1_Lookup,Len2_Lookup)
INTEGER :: Dump_Pack
INTEGER :: IM_Ident

! Local variables

INTEGER            :: i        ! Loop index
INTEGER            :: item     ! Stash Item No
INTEGER            :: sect     ! Stash Section No
INTEGER            :: number_of_data_words_in_memory
INTEGER            :: number_of_data_words_on_disk
INTEGER            :: disk_address
INTEGER            :: pack_code       ! Packing code in lookup
INTEGER            :: stash_pack_code ! Packing code according to STASH
INTEGER            :: ErrorStatus

LOGICAL            :: prognostic      ! T : Item is prognostic
LOGICAL            :: diagnostic      ! T : Item is diagnostic
LOGICAL            :: packing_changed ! T : packing codes changed
LOGICAL            :: real_data       ! T : real data

CHARACTER (LEN=errormessagelength) :: CMessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='CHECK_DUMP_PACKING'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!
!- End of header

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
ErrorStatus = 0
CMessage = ' '

packing_changed = .FALSE.

DO i = 1, Len2_Lookup

  !       Extract details from LookUp Table

  pack_code = MOD (Lookup(lbpack,i), 10)
  sect      = Lookup(item_code,i) / 1000
  item      = MOD ( Lookup(item_code,i) ,1000 )
  real_data = (Lookup(data_type,i) == 1)

  prognostic = (sect == 0)
  diagnostic = (sect >  0)

  SELECT CASE (pack_code)  !  Packing code in LookUp Table

  CASE (PC_No_Packing)

    IF ( Real_data      .AND.                                   &
       ( Dump_Pack == 1 .OR.                                    &
                                          ! Prog & Diag packed
        (Dump_Pack == 2 .AND. Diagnostic))                      &
                                          ! Diag packed
       ) THEN

      ! Set packing indicator according to StashMaster record
      stash_pack_code =                                         & 
             exppxi (im_ident, sect, item, ppx_dump_packing,    &
             errorstatus, cmessage)

      ! Preserve the value of the N4 position of the packing code
      ! unless the STASH entry specifically overrides it 
      IF (stash_pack_code / 1000 == PC_Native_Format) THEN
        stash_pack_code = stash_pack_code + (Lookup(lbpack,i)/1000)*1000
      END IF

      Lookup(lbpack,i) = stash_pack_code
      packing_changed = .TRUE.

    END IF

  CASE (PC_Cray32_Packing)

    IF (  Real_data      .AND.                                  &
       ( (Dump_Pack == 2 .AND. Prognostic) .OR.                 &
                                                ! Prog unpacked
          Dump_Pack == 3 )                                      &
                                         ! Prog & Diag unpacked
       ) THEN

      !             Set packing indicator to no packing

      Lookup(lbpack,i) = ( Lookup(lbpack,i)/10 ) * 10 + PC_No_Packing

      packing_changed = .TRUE.

    END IF

  CASE DEFAULT

    ErrorStatus = 10
    WRITE(cmessage, '(2(A,I0))')                                     newline//&
      'Unexpected Packing code in dump.'//                           newline//&
      'Field No ', i, ' Packing Code ', Pack_Code

    CALL EReport (RoutineName, ErrorStatus, CMessage)

  END SELECT

END DO

!     -------------------------------------------------------
!     If any of the packing indicators have been changed then
!     reset the disk addresses and lengths
!     -------------------------------------------------------

IF (packing_changed) THEN

  ! DEPENDS ON: set_dumpfile_address
  CALL Set_DumpFile_Address (FixHd, Len_FixHd,                    &
                             Lookup, Len1_Lookup, Len2_Lookup,    &
                             number_of_data_words_in_memory,      &
                             number_of_data_words_on_disk,        &
                             disk_address)

  IF (PrintStatus >= PrStatus_Normal) THEN
    ErrorStatus = -20   !  Warning
    WRITE (cmessage,'(A)')                                           newline//&
      'Packing codes in dump inconsistent with DUMP_PACKim.'//       newline//&
      'Packing codes updated.'

    CALL EReport (RoutineName, ErrorStatus, CMessage)
  END IF

END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Check_Dump_Packing

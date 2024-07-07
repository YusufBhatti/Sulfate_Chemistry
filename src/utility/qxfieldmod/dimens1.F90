! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: DIMENS1--------------------------------------------
!
! Purpose: To read a   direct access PP file  and convert it to a
! sequential file read to be passed across to the IBM
!
! Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
! -------------------------------------------------------------------
! Interface and arguments: ------------------------------------------
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Small execs
SUBROUTINE dimens1(len_inthd,len_realhd,len1_levdpc,len2_levdpc,  &
   len1_rowdpc,len2_rowdpc, len1_coldpc,len2_coldpc,              &
   len1_lookup,len2_lookup,len_fixhd,pp_fixhd,ppunit1,ppunit2,    &
   icode,cmessage)

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER ::                                                        &
     len_inthd                                                    &
    ,len_fixhd                                                    &
    ,len_realhd                                                   &
    ,len1_levdpc                                                  &
    ,len2_levdpc                                                  &
    ,len1_rowdpc                                                  &
    ,len2_rowdpc                                                  &
    ,len1_coldpc                                                  &
    ,len2_coldpc                                                  &
    ,len1_lookup                                                  &
    ,len2_lookup                                                  &
    ,lookup(len1_lookup,len2_lookup)                              &
    ,pp_fixhd(len_fixhd)                                          &
    ,icode                                                        &
    ,ppunit1                                                      &
    ,ppunit2
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
!
! DEPENDS ON: readpp
CALL readpp(len_inthd,len_realhd,len1_levdpc,len2_levdpc,         &
     len1_rowdpc,len2_rowdpc,len1_coldpc,len2_coldpc,len1_lookup, &
     len2_lookup,len_fixhd,pp_fixhd,lookup,lookup,ppunit1,        &
     ppunit2,icode,cmessage)
9999 CONTINUE
IF (icode /= 0) RETURN
RETURN
END SUBROUTINE dimens1

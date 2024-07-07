! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Purpose: Driver level to convert fieldsfiles to PP files written using
!             the IBM number format
!
!    Programming standards: UMDP 3
!
!    Documentation: UM Doc Paper F5
!
!    -----------------------------------------------------------------
!    Arguments:-------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs

SUBROUTINE atmos_convpp64_ibm                                    &
 (len_fixhd,len_inthd,len_realhd,                                &
 len1_levdepc,len2_levdepc,len1_rowdepc,                         &
 len2_rowdepc,len1_coldepc,len2_coldepc,                         &
 len1_flddepc,len2_flddepc,len_extcnst,                          &
 len_dumphist,len_cfi1,len_cfi2,len_cfi3,                        &
 len1_lookup,len2_lookup,len_data,p_field,                       &
 nftin,nftout,max_field_size)

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE lookup_addresses
USE readflds_mod, ONLY: readflds
USE umPrintMgr

IMPLICIT NONE

INTEGER ::                                                        &
 len_fixhd                                                        &
              !IN Length of fixed length header on input file
,len_inthd                                                        &
              !IN Length of integer header on input file
,len_realhd                                                       &
              !IN Length of real header on input file
,len1_levdepc                                                     &
              !IN 1st dim of lev dependent consts on input file
,len2_levdepc                                                     &
              !IN 2nd dim of lev dependent consts on input file
,len1_rowdepc                                                     &
              !IN 1st dim of row dependent consts on input file
,len2_rowdepc                                                     &
              !IN 2nd dim of row dependent consts on input file
,len1_coldepc                                                     &
              !IN 1st dim of col dependent consts on input file
,len2_coldepc                                                     &
              !IN 2nd dim of col dependent consts on input file
,len1_flddepc                                                     &
              !IN 1st dim of field dependent consts on input fi
,len2_flddepc                                                     &
              !IN 2nd dim of field dependent consts on input fi
,len_extcnst                                                      &
              !IN Length of extra consts on input file
,len_dumphist                                                     &
              !IN Length of history header on input file
,len_cfi1                                                         &
              !IN Length of index1 on input file
,len_cfi2                                                         &
              !IN Length of index2 on input file
,len_cfi3                                                         &
              !IN Length of index3 on input file
,len1_lookup                                                      &
              !IN 1st dim of LOOKUP on input file
,len2_lookup                                                      &
              !IN 2nd dim of LOOKUP on input file
,len_data                                                         &
              !IN Length of data on input file
,p_field                                                          &
              !IN No of p-points per level on input file
,max_field_size !Maximum field size on file

INTEGER :: nftin      ! Input file unit
INTEGER :: nftout     ! Output file unit

! ----------------------------------------------------------------------
! Local variables:
! ----------------------------------------------------------------------

LOGICAL :: last   ! Final field in file variable

LOGICAL, PARAMETER :: model_flag = .FALSE. ! Fieldsfile input 
LOGICAL, PARAMETER :: lcal360 = .FALSE.
LOGICAL, PARAMETER :: oper=.FALSE.

INTEGER :: i               ! Loop indices
INTEGER :: icode           ! Error return code from subroutines    
INTEGER :: num_fields      ! Number of fields in lookup
INTEGER :: int_dim, num_values, iextraw, data_add
INTEGER :: start_block
INTEGER :: iextra(10) = 0

CHARACTER(LEN=errormessagelength) :: cmessage                

! Local arrays: --------------------------------------------------------
! Integer file headers
INTEGER :: fixhd(len_fixhd)
INTEGER :: inthd(len_inthd)
INTEGER :: cfi1(len_cfi1+1)
INTEGER :: cfi2(len_cfi2+1)
INTEGER :: cfi3(len_cfi3+1)
INTEGER :: lookup(len1_lookup,len2_lookup)        ! Lookup table
INTEGER :: lookup_out(len1_lookup)                ! Output lookup table

! Real arrays
REAL    :: realhd(len_realhd)
REAL    :: levdepc(1+len1_levdepc*len2_levdepc)
REAL    :: rowdepc(1+len1_rowdepc*len2_rowdepc)
REAL    :: coldepc(1+len1_coldepc*len2_coldepc)   ! Real file headers
REAL    :: flddepc(1+len1_flddepc*len2_flddepc)   ! Real file headers
REAL    :: extcnst(len_extcnst+1)                 ! Real file headers
REAL    :: dumphist(len_dumphist+1)

! ----------------------------------------------------------------------
cmessage = ' '
icode = 0
last = .FALSE.
num_fields = 0
iextra(1)=1   ! Required to write packed data

! DEPENDS ON: readhead
CALL readhead(nftin,fixhd,len_fixhd,                              &
                inthd,len_inthd,                                  &
                realhd,len_realhd,                                &
                levdepc,len1_levdepc,len2_levdepc,                &
                rowdepc,len1_rowdepc,len2_rowdepc,                &
                coldepc,len1_coldepc,len2_coldepc,                &
                flddepc,len1_flddepc,len2_flddepc,                &
                extcnst,len_extcnst,                              &
                dumphist,len_dumphist,                            &
                cfi1,len_cfi1,                                    &
                cfi2,len_cfi2,                                    &
                cfi3,len_cfi3,                                    &
                lookup,len1_lookup,len2_lookup,                   &
                len_data,                                         &
                start_block,icode,cmessage)

IF (icode /= 0) THEN
  WRITE(umMessage,'(A,I0)')cmessage,icode
  CALL umPrint(umMessage,src='atmos_convpp_ibm.F90')
  CALL ereport('ATMOS_CONVPP_IBM', icode, cmessage)
END IF

data_add=fixhd(160)-1 ! Start address for the data.

DO i=1,len2_lookup
  IF (lookup(lbrow,i) /= -99) THEN
    num_fields = num_fields + 1
  END IF
END DO

DO i = 1, num_fields
  IF (i == num_fields) last=.TRUE.
  num_values=lookup(lbrow,i)*lookup(lbnpt,i)+lookup(lbext,i)
  iextraw=0
  IF (lookup(lbext,i) >  0) THEN ! got some extra data
    iextraw=lookup(lbext,i)
  END IF
  int_dim=((num_values+1)/2)*2 ! Round to ensure an integer for IBM
! DEPENDS ON: cray_ibm
  CALL cray_ibm(int_dim,num_values,nftin,                        &
            len1_lookup,len2_lookup,fixhd,lookup,                &
            lookup,i,data_add,model_flag,                        &
            nftout,iextra,iextraw,last,oper,                     &
            icode,cmessage,lcal360)

  IF (icode /= 0) THEN
    WRITE(umMessage,'(A,I0)')cmessage,icode
    CALL umPrint(umMessage,src='atmos_convpp_ibm.F90')
    CALL ereport('ATMOS_CONVPP_IBM', icode, cmessage)
  END IF
END DO


END SUBROUTINE atmos_convpp64_ibm

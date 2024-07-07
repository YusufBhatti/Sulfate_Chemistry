! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Storage of Addressing variables

MODULE Rcf_Address_Vars_Mod

! Description:
!   Module specifically for variable storage. Variables used in
!   RCF_Address and TSTMSK.
!
! Method:
!   Declare a number of variables and USE them from the appropriate
!   locations
!
!   Original variables stripped from rcf_address.F90 to prevent circular
!   usage of modules between rcf_address and tstmsk. FCM didn't like
!   the pre-exisiting arrangement.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE Rcf_Ppx_Info_Mod, ONLY:    &
    STM_OptCodeLen

IMPLICIT NONE

INTEGER, SAVE   ::    ispace         ! Space code
INTEGER, SAVE   ::    igp            ! Grid of data code
INTEGER, SAVE   ::    ilev           ! Level type code
INTEGER, SAVE   ::    ibot           ! First level code
INTEGER, SAVE   ::    itop           ! Last level code
INTEGER, SAVE   ::    iflag          ! Level compression flag
INTEGER, SAVE   ::    iopn(STM_OptCodeLen/5) ! Sectional option code
INTEGER, SAVE   ::    vmsk           ! Integer equiv of bin vers mask
INTEGER, SAVE   ::    ipseudo        ! Pseudo dimension type
INTEGER, SAVE   ::    ipfirst        ! First pseudo dim code
INTEGER, SAVE   ::    iplast         ! Last pseudo dim code
INTEGER, SAVE   ::    halo           ! Halo type code

END MODULE Rcf_Address_Vars_Mod

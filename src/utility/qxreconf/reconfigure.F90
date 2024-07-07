! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Top level program for reconfiguration
! Description:
!  Top level program for reconfiguration
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

PROGRAM Reconfigure

USE Rcf_Initialise_Mod, ONLY: &
    Rcf_Initialise

USE Rcf_Finalise_Mod, ONLY: &
    Rcf_Finalise

USE Rcf_Read_Namelists_Mod, ONLY: &
    Rcf_Read_Namelists

USE Rcf_Control_Mod, ONLY: &
    Rcf_Control

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE UM_ParCore, ONLY: &
    mype,             &
    nproc,            &
    nproc_max

USE io, ONLY: ioShutdown

USE UM_Config, ONLY: &
    appInit,          &
    appTerminate,     &
    exe_RCF

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage
IMPLICIT NONE

TYPE (um_header_type)        :: hdr_in       ! header from input dump
TYPE (um_header_type)        :: hdr_out      ! header from output dump
INTEGER                      :: ErrorStatus
INTEGER                      :: info         ! GCOM dummy
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'Reconfigure'

INTEGER, PARAMETER :: gc_alltoall_version = 2    ! From GCOM
INTEGER, PARAMETER :: gc_alltoall_multi   = 2    ! From GCOM

CALL appInit(exe_RCF)
nproc=nproc_max

! Set GCOM to use the alternative version of RALLTOALLE
! throughout the run
CALL Gc_Setopt(gc_alltoall_version, gc_alltoall_multi, Errorstatus)

!-------------------------------------------------------------------
! Perform initialisation
!-------------------------------------------------------------------
CALL Rcf_Initialise( hdr_in, hdr_out )

!-------------------------------------------------------------------
! Do the real work
!-------------------------------------------------------------------
CALL Rcf_Control( hdr_in, hdr_out )

!------------------------------------------------------------------
! Tidy Up
!------------------------------------------------------------------
CALL Rcf_Finalise( hdr_in, hdr_out )

!------------------------------------------------------------------
! Report on IO if configured
!------------------------------------------------------------------

CALL appTerminate()
WRITE(umMessage,*) ' End of rcf program reached. PE ',mype
CALL umPrint(umMessage,src='reconfigure')

END PROGRAM Reconfigure

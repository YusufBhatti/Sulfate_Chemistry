! *****************************COPYRIGHT*******************************
! (c) [University of California] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution.
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Fast-jx is an updated routine for calculating online photolysis rates
!   This calls routines to initialise varous quantities. Based on the routine
!   fastj_inphot
!
!  Routine to initialise photolysis rate data, called directly from
!  ukca_fastjx. Currently use it to read the JPL spectral data
!  and standard O3 and T profiles and to set the appropriate reaction index.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE fastjx_inphot_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FASTJX_INPHOT_MOD'

CONTAINS
SUBROUTINE fastjx_inphot

USE fastjx_data,         ONLY: w_
USE fastjx_specs
USE ukca_option_mod,     ONLY: jvspec_dir, jvspec_file,    &
                               jvscat_file, jvsolar_file,  &
                               fastjx_numwl, i_ukca_solcyc
USE ereport_mod,         ONLY: ereport
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim
USE file_manager,        ONLY: assign_file_unit, release_file_unit

USE fastjx_set_aer_mod, ONLY: fastjx_set_aer
IMPLICIT NONE

CHARACTER(LEN=255)       :: jv_fullpath
LOGICAL                  :: l_exist
INTEGER                  :: errorstatus
INTEGER                  :: ukcafjxx_unit
INTEGER                  :: ukcafjsc_unit
INTEGER                  :: ukcafjsol_unit

CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'FASTJX_INPHOT'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set number of wavelengths from RUN_UKCA namelist
w_ = fastjx_numwl

! Read in labels of photolysis rates required
CALL fastjx_rd_js

!--------------------------------------------------
! Read in Fast-JX photolysis cross sections (based on N. Savage approach in fastj)
l_exist = .FALSE.
jv_fullpath=TRIM(jvspec_dir)//'/'//TRIM(jvspec_file)
INQUIRE (FILE=TRIM(jv_fullpath), EXIST=l_exist)
IF (.NOT. l_exist) THEN
  errorstatus=1
  cmessage = TRIM(jv_fullpath)//': Fast-JX spectral file does not exist'
  CALL ereport(RoutineName,errorstatus,cmessage)
END IF

! Read in Fast-JX photolysis cross sections
CALL assign_file_unit(jv_fullpath, ukcafjxx_unit, handler="fortran")
CALL fastjx_rd_xxx(ukcafjxx_unit, jv_fullpath)
CALL release_file_unit(ukcafjxx_unit, handler="fortran")

!--------------------------------------------------
! Read in Fast-JX photolysis solar cycle 
IF (i_ukca_solcyc>0) THEN  
  l_exist = .FALSE.
  jv_fullpath=TRIM(jvspec_dir)//'/'//TRIM(jvsolar_file)
  INQUIRE (FILE=TRIM(jv_fullpath), EXIST=l_exist)
  IF (.NOT. l_exist) THEN
    errorstatus=1
    cmessage = TRIM(jv_fullpath)//': Fast-JX solar cycle file does not exist'
    CALL ereport(RoutineName,errorstatus,cmessage)
  END IF

  ! Read in Fast-JX photolysis solar cycle
  CALL assign_file_unit(jv_fullpath, ukcafjsol_unit, handler="fortran")
  CALL fastjx_rd_sol(ukcafjsol_unit, jv_fullpath)
  CALL release_file_unit(ukcafjsol_unit, handler="fortran")
END IF

!----------------------------------------------
! Read in Fast-JX scattering cross sections
l_exist = .FALSE.
jv_fullpath=TRIM(jvspec_dir)//'/'//TRIM(jvscat_file)
INQUIRE (FILE=TRIM(jv_fullpath), EXIST=l_exist)
IF (.NOT. l_exist) THEN
  errorstatus=1
  cmessage = TRIM(jv_fullpath)//': Fast-JX scattering file does not exist'
  CALL ereport(RoutineName,errorstatus,cmessage)
END IF

! Read in Fast-JX photolysis cross sections
CALL assign_file_unit(jv_fullpath, ukcafjsc_unit, handler="fortran")
CALL fastjx_rd_mie(ukcafjsc_unit, jv_fullpath)
CALL release_file_unit(ukcafjsc_unit, handler="fortran")

!     Select Aerosol/Cloud types to be used
CALL fastjx_set_aer

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_inphot
END MODULE fastjx_inphot_mod


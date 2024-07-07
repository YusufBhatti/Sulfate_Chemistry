! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    4 Subroutines in deck : RDOBS, RDOBS2, RDOBS3 and DAYS -----
!
!    Purpose : Read from ACOBS Files,reformat and place OBS header
!              details in COMOBS. The bulk of the required OBS data
!              is put into dynamic work array OBS for transmission via
!              argument list to GETOBS. OBS is written out to a cache
!              file for subsequent reading at later timesteps.
!              Thus reread of ACOBS files only required intermittently
!              (The routine DAYS does a dd/mm/yy to dayno)
!
!
!   S.Bell      <- programmer of some or all of previous code or changes
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Logical components covered:
!
!    Project Task : P3
!
!    External documentation:
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE rdobs3_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RDOBS3_MOD'

CONTAINS

SUBROUTINE rdobs3 (unit_no,nobtypmx,nobtyp,obstyp,nobs,ndatav,    &
                   noblev,oblevtyp,obs_levels,obs_data,           &
                   len_data,max_ndv,tnobs,missd,                  &
                   icode,cmessage                                 &
                         ,ipt                                     &
                             )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io
USE ac_dump_mod, ONLY: len_inthd, len_realhd, len1_levdepc,       &
                       len2_levdepc, len1_rowdepc, len2_rowdepc,  &
                       len1_coldepc, len2_coldepc, len1_flddepc,  &
                       len2_flddepc, len_extcnst, len_dumphist,   &
                       len_cfi1, len_cfi2, len_cfi3,              &
                       len1_lookup_obs, len2_lookup_obs, fixhd,   &
                       inthd, realhd, levdepc, rowdepc, coldepc,  &
                       flddepc, extcnst, dumphist, cfi1, cfi2,    &
                       cfi3, lookup
USE nlsizes_namelist_mod, ONLY: model_levels, len_fixhd
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

! ----------------------------------------------------------------------
INTEGER :: ipt
INTEGER ::                                                        &
   unit_no                                                        &
                    ! IN  Unit no of observation file
  ,nobtypmx                                                       &
                    ! IN  Max no of observation types
  ,nobtyp                                                         &
                    ! OUT No of observation types
  ,obstyp(nobtypmx)                                               &
                    ! OUT Observation types
  ,nobs(nobtypmx)                                                 &
                    ! OUT No of Observations
  ,ndatav(nobtypmx)                                               &
                    ! OUT No of data values
  ,noblev(nobtypmx)                                               &
                    ! OUT No of observation levels
  ,oblevtyp(nobtypmx)                                             &
                     !OUT Observation level type
  ,len_data                                                       &
                    ! IN  Dimension of data section
  ,tnobs                                                          &
                    ! OUT Total no of observations
  ,icode            ! OUT Return Code

REAL ::                                                           &
   missd                                                          &
                             ! OUT Real missing data indicator
  ,obs_levels(model_levels+1,*)                                   &
                             ! OUT Observation levels
  ,obs_data(len_data)        ! OUT Observation data

!DR   REAL DATALEVS(LEN1_LEVDEPC-2,*)     ! OUT Obs levels

CHARACTER(LEN=errormessagelength) :: cmessage ! OUT Error message if ICODE > 0
! ----------------------------------------------------------------------
!

!     LOCAL VARIABLES
INTEGER :: jobt,jlev
INTEGER :: max_ndv   !  Maximum no of data values (for obs type)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RDOBS3'

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Go to start of obs file

CALL setpos (unit_no,0,icode)

!     Read in the observation file
! DEPENDS ON: readacobs
CALL readacobs (unit_no,                                          &
  fixhd,    len_fixhd,                                            &
  inthd,    len_inthd,                                            &
  realhd,   len_realhd,                                           &
  levdepc,  len1_levdepc,len2_levdepc,                            &
  rowdepc,  len1_rowdepc,len2_rowdepc,                            &
  coldepc,  len1_coldepc,len2_coldepc,                            &
  flddepc,  len1_flddepc,len2_flddepc,                            &
  extcnst,  len_extcnst,                                          &
  dumphist, len_dumphist,                                         &
  cfi1,     len_cfi1,                                             &
  cfi2,     len_cfi2,                                             &
  cfi3,     len_cfi3,                                             &
  lookup,                                                         &
  len1_lookup_obs,len2_lookup_obs,                                &
               len_data,obs_data,                                 &
               icode,cmessage                                     &
                         ,ipt                                     &
                             )

IF (icode >  0) GO TO 9999

tnobs  = inthd(28)   !  Total no of observations
nobtyp = inthd(32)   !  No of observation types
missd  = realhd(29)  !  Real MDI used in observations

IF (nobtyp >  0) THEN

  DO jobt = 1,nobtyp
    obstyp(jobt) = lookup(65,jobt)   ! Observation Types
    nobs  (jobt) = lookup(66,jobt)   ! No of obs
    ndatav(jobt) = lookup(67,jobt)   ! No of data values
    noblev(jobt) = levdepc(2,jobt)   ! No of obs levels
    oblevtyp(jobt) = levdepc(1,jobt) ! Level type
    DO jlev =1,len1_levdepc-2
      obs_levels(jlev,jobt) = levdepc(2+jlev,jobt) ! Obs levels
    END DO
  END DO

  !       MAXNLEV1 = LEN1_LEVDEPC-2  ! Max no of levels + 1

END IF
9999  CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rdobs3
END MODULE rdobs3_mod

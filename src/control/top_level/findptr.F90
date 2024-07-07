! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: FINDPTR  -------------------------------------------------
!
! Purpose: Locates address within D1 of diagnostic field which may
!          be required elsewhere in the model for special-purpose
!          diagnostic routine such as zonal mean print, or as an
!          internal interfacing field for coupling sub-models.
!          The search information is input in STASH format, and the
!          STASH list is scanned for a match.  If the specified
!          field does not exist in D1 the address is returned as 0.
!          NB: Missing data indicators may be supplied if the search
!              is to ignore certain elements in the STASH list.
!
! Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
! External documentation:
!   Unified Model Doc Paper C4 - Storage Handling and
!                                Diagnostic System.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE findptr ( internal_model,section,item,                      &
                     process_code,freq_code,start_step,end_step,period,&
                     gridpt_code,weight_code,                          &
                     bottom_level,top_level,                           &
                     grid_n,grid_s,grid_w,grid_e,                      &
                     stashmacro_tag,mdi,address,                       &
                     icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE stash_array_mod, ONLY: stindex, stlist
USE stparam_mod, ONLY: st_output_addr, st_macrotag, s_modl, s_sect,&
                       s_item, s_output, s_proc, s_freq, s_times,  &
                     s_timee, s_period, s_grid, s_weight, s_bottom,&
                     s_top, s_north, s_south, s_west, s_east
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: internal_model   ! internal_model id.
INTEGER, INTENT(IN) :: section          ! STASH section number
INTEGER, INTENT(IN) :: item             ! STASH item number
INTEGER, INTENT(IN) :: process_code     ! STASH processing code
INTEGER, INTENT(IN) :: freq_code        ! STASH frequency code
INTEGER, INTENT(IN) :: start_step       ! STASH start step for processing
INTEGER, INTENT(IN) :: end_step         ! STASH end step for processing
INTEGER, INTENT(IN) :: period           ! STASH processing period
INTEGER, INTENT(IN) :: gridpt_code      ! STASH gridpoint code
INTEGER, INTENT(IN) :: weight_code      ! STASH weighting code
INTEGER, INTENT(IN) :: bottom_level     ! STASH input bottom level
INTEGER, INTENT(IN) :: top_level        ! STASH input top level
INTEGER, INTENT(IN) :: grid_n           ! STASH N-row grid code
INTEGER, INTENT(IN) :: grid_s           ! STASH S-row grid code
INTEGER, INTENT(IN) :: grid_w           ! STASH W-col grid code
INTEGER, INTENT(IN) :: grid_e           ! STASH E-col grid code
INTEGER, INTENT(IN) :: stashmacro_tag   ! STASHmacro tag number
INTEGER, INTENT(IN) :: mdi              ! Missing Data Indicator

INTEGER, INTENT(OUT) :: address         ! Address in D1
INTEGER, INTENT(OUT) :: icode           ! Error return code

CHARACTER(LEN=errormessagelength), INTENT(OUT) ::                 &
                        cmessage        ! Error return message

! ----------------------------------------------------------------------
!  Local variables
!
INTEGER ::                                                        &
    istart,iend,i,                                                &
                            ! Start, end + loop index in STASHlist
    nmatch                                                        &
                            ! Number of matches found
    ,im_index               ! Internal model index
LOGICAL ::                                                        &
    match                   ! TRUE if diagnostic matched

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FINDPTR'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Initialise error code to zero, which would indicate success
! to the calling routine, since all error conditions will 
! explicitly set this return code to a non-zero value.   
icode = 0

! ----------------------------------------------------------------------
!  0.  Check that tag field is within the allowed range for user tags
!
IF (stashmacro_tag /= mdi .AND.                                   &
   (stashmacro_tag <  0 .OR. stashmacro_tag >  999)) THEN
  cmessage="FINDPTR : STASHMACRO_TAG must be in range 0-999"
  icode=ABS(stashmacro_tag)
  GO TO 9999
END IF
! ----------------------------------------------------------------------
!  1.  Locate start/end limits within STASHlist for search;
!      initialise output ADDRESS to zero
!
address=0
nmatch=0
im_index=1
IF (stindex(2,item,section,im_index) >  0) THEN
  istart=stindex(1,item,section,im_index)
  iend  =stindex(2,item,section,im_index)+istart-1
  !
  !  1.1 Loop over STASHlist entries and try to find matches
  !
  DO i=istart,iend
    IF (stlist(s_modl,i) /= internal_model .OR.                    &
        stlist(s_sect,i) /= section .OR.                           &
        stlist(s_item,i) /= item) THEN
      icode=1000*section+item
      cmessage="FINDPTR : Corrupt STASHlist or STASHindex"
      GO TO 9999
    END IF
    match=((stlist(s_output,i) == 1) .OR.                          &
           (stlist(s_output,i) == 2)) .AND.                        &
    (process_code == stlist(s_proc,i)                             &
                                .OR. process_code == mdi) .AND.     &
    (freq_code   == stlist( s_freq,i)                             &
                                .OR. freq_code   == mdi) .AND.     &
    (start_step  == stlist( s_times,i)                            &
                                .OR. start_step  == mdi) .AND.     &
    (end_step    == stlist( s_timee,i)                            &
                                .OR. end_step    == mdi) .AND.     &
    (period      == stlist( s_period,i)                           &
                                .OR. period      == mdi) .AND.     &
    (gridpt_code == stlist( s_grid,i)                             &
                                .OR. gridpt_code == mdi) .AND.     &
    (weight_code == stlist( s_weight,i)                           &
                                .OR. weight_code == mdi) .AND.     &
    (bottom_level == stlist(s_bottom,i)                           &
                                .OR. bottom_level == mdi) .AND.    &
    (top_level   == stlist(s_top,i)                               &
                                .OR. top_level   == mdi) .AND.     &
    (grid_n      == stlist(s_north,i)                             &
                                .OR. grid_n      == mdi) .AND.     &
    (grid_s      == stlist(s_south,i)                             &
                                .OR. grid_s      == mdi) .AND.     &
    (grid_w      == stlist(s_west,i)                              &
                                .OR. grid_w      == mdi) .AND.     &
    (grid_e      == stlist(s_east,i)                              &
                                .OR. grid_e      == mdi) .AND.     &
    (stashmacro_tag == MOD(stlist(st_macrotag,i),1000) .OR.        &
     stashmacro_tag == mdi)
    !
    IF (match) THEN
      address=stlist(st_output_addr,i)
      nmatch=nmatch+1
    END IF
  END DO
  !
  IF (nmatch >  1) THEN
    icode=-1000*section-item
    cmessage="FINDPTR : Warning - multiple match for diagnostic"
    WRITE(umMessage,*)"FINDPTR : Warning - multiple match for diagnostic ",  &
                section,item
    CALL umPrint(umMessage,src='findptr')
    !
  END IF
END IF
!
9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
! ----------------------------------------------------------------------
END SUBROUTINE findptr

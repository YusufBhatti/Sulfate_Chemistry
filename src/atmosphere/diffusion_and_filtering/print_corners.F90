! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************

! subroutine print_corners

SUBROUTINE print_corners(                                         &
                         field, row_length, rows, levels,         &
                         level_1, level_2,                        &
                         off_x, off_y, off_u, off_v,              &
                         l_datastart, n_proc )

USE pr_corner_mod
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
! Purpose:
!          To print (sub-)domain corner values of a field
!
! Method:
!          Corner values in PR_CORNER_MOD and are set from
!          the run_diagnostics NAMELIST and are initialised in
!          init_pr_corner called from  SETCONA
!          Place CALL anywhere in the UM to print arrays
!          e.g. printing theta values from ATM_STEP
!
!  Whole domain values
!      dom_w = 1         dom_e = global_row_length
!      dom_s = 1         dom_n = global_rows

!   Corner size ix,jy  do not have to be equal
!    Max values  ix = 4 ,    jy = 4
!     write(6,*)' *** THETA  ***'
!! DEPENDS ON: print_corners
!      call print_corners(                                               &
!                          THETA, row_length, rows, model_levels,        &
!                          level_1, level_2,
!                          offx, offy, 0, 0,                             &
!                          datastart, nproc )
!
!         off_u  should be set to 1 for calling U fields in a LAM
!         off_v  should be set to 1 for calling V fields
!
!  NB Calling routine must contain the argument list values
!     In particular datastart(2), nproc, and model_type
!     In some subroutines
!         datastart  may be named l_datastart
!         nproc may be named n_proc
!         offx, offy may be named off_x, off_y
!
!          The corners are defined by
!                   left(West)    right(East)
!    upper row  =  dom_w,dom_n    dom_e,dom_n
!    lower row  =  dom_w,dom_s    dom_e,dom_s
!
!           Values are printed out as below
!
!                      1                             ix
! upper      left   dom_w,dom_n            ... dom_w+ix-1,dom_n
! upper-jy+1 left   dom_w,dom_n-jy+1       ... dom_w+ix-1,dom_n-jy+1
!
!                      1                             ix
! upper      right  dom_e-ix+1,dom_n       ... dom_e-ix+1,dom_n
! upper-jy+1 right  dom_e-ix+1,dom_n-jy+1  ... dom_e-ix+1,dom_n-jy+1
!
!                      1                             ix
! lower+jy-1 left   dom_w,dom_s+jy-1       ... dom_w+ix-1,dom_s+jy-1
! lower      left   dom_w,dom_s            ... dom_w+ix-1,dom_s
!
!                      1                             ix
! lower+jy-1 right  dom_e-ix+1,dom_s+jy-1  ... dom_e-ix+1,dom_s+jy-1
! lower      right  dom_e-ix+1,dom_s       ... dom_e-ix+1,dom_s

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage
USE model_domain_mod, ONLY: model_type, mt_lam

IMPLICIT NONE

! Parameters required for dimensioning some of the arguments

! Arguments with Intent IN. ie: Input variables.
INTEGER, INTENT(IN) ::                                            &
  row_length                                                      &
                     ! in: no of points per local row
, rows                                                            &
                     ! in: no of local rows
, levels                                                          &
                     ! in: no of levels in input field
, level_1                                                         &
                     ! in: level to print
, level_2                                                         &
                     ! in: level to print
, off_x , off_y                                                   &
                     !  halo sizes of field
, off_u                                                           &
                     !  offsets(=1) for LAM u type fields
, off_v                                                           &
                     !  offsets(=1) for LAM v type fields
, l_datastart(3)                                                  &
, n_proc
                     ! Total number of processors


REAL, INTENT(IN) ::                                               &
 field(1-off_x:row_length+off_x,1-off_y:rows+off_y, levels)

!    local variables

! Loop indices/pointers
INTEGER ::                                                        &
  i, j, k                                                         &
, ii, ji                                                          &
, info                                                            &
, dom_n, dom_e                                                    &
, ie, ies, iee                                                    &
, iw, iws, iwe                                                    &
, jn, jns, jne                                                    &
, js, jss, jse                                                    &
, lev_in(2)                                                       &
, num_lev                                                         &
, off_i                                                           &
                      ! off_set for u if LAM
, SIZE

!    local arrays

REAL :: ll(4,4,2)
REAL :: lr(4,4,2)
REAL :: ul(4,4,2)
REAL :: ur(4,4,2)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_CORNERS'

! ----------------------------------------------------------------------
! Section 1.  Set up pointers
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Subtract offsets if set
dom_n = dom_no - off_v
dom_e = dom_eo - off_u

!   Level numbers of field levels to be printed choose 2 values
!   If levels = 1 (i.e. single level field) only 1 level printed

lev_in(1) = level_1
lev_in(2) = level_2

num_lev = 2
IF ( lev_in(1) == lev_in(2)) num_lev = 1

IF (model_type == mt_lam) THEN
  off_i = off_u
ELSE !  model_type \= mt_lam
  off_i = 0
END IF !  model_type == mt_lam

! Set defaults to no-op values
!  (i.e. arrays left empty unless at a corner of desired sub-domain)
iws = row_length + 1
iwe = 0
ies = row_length + 1
iee = 0
jss = rows + 1
jse = 0
jns = rows + 1
jne = 0
iw = 0
ie = 0
jn = 0
js = 0

IF (dom_w > l_datastart(1) - 1 .AND.                              &
    dom_w < l_datastart(1) + row_length - off_i ) THEN
  iws = dom_w - l_datastart(1) + 1
  iwe = iws + ix - 1
  IF ( iwe > row_length ) iwe = row_length
ELSE IF (dom_w + ix > l_datastart(1) .AND.                        &
              dom_w < l_datastart(1)) THEN
  iws = 1
  iwe = dom_w + ix - l_datastart(1)
  iw = ix - iwe + iws
END IF !  dom_w > l_datastart(1) + row_length - 1

IF (dom_e - ix + 1 > l_datastart(1) - 1 .AND.                     &
    dom_e - ix + 1 < l_datastart(1) + row_length - off_i ) THEN
  ies = dom_e - ix - l_datastart(1) + 2
  iee = ies + ix - 1
  IF ( iee > row_length ) iee = row_length
ELSE IF (dom_e + 1 > l_datastart(1) .AND.                         &
    dom_e - ix + 1 < l_datastart(1)) THEN
  ies = 1
  iee = dom_e + 1 - l_datastart(1)
  ie = ix - iee + ies
END IF !  dom_e - ix + 1 > l_datastart(1) + row_length - 1

IF (dom_s > l_datastart(2) - 1 .AND.                              &
    dom_s < l_datastart(2) + rows ) THEN
  jss = dom_s - l_datastart(2) + 1
  jse = jss + jy - 1
  IF ( jse > rows ) jse = rows
ELSE IF (dom_s + jy > l_datastart(2) .AND.                        &
              dom_s < l_datastart(2)) THEN
  jss = 1
  jse = dom_s + jy - l_datastart(2)
  js = jy - jse + jss
END IF !  dom_s > l_datastart(2) + rows - 1

IF (dom_n - jy + 1 > l_datastart(2) - 1 .AND.                     &
    dom_n - jy + 1 < l_datastart(2) + rows ) THEN
  jns = dom_n - jy - l_datastart(2) + 2
  jne = jns + jy - 1
  IF ( jne > rows ) jne = rows
ELSE IF (dom_n + 1 > l_datastart(2) .AND.                         &
    dom_n - jy + 1 < l_datastart(2)) THEN
  jns = 1
  jne = dom_n + 1 - l_datastart(2)
  jn = jy - jne + jns
END IF !  dom_n - jy + 1 > l_datastart(2) + rows - 1

! ----------------------------------------------------------------------
! Section 2.  Initialise arrays
! ----------------------------------------------------------------------

ll = 0.0      ! lower left values
lr = 0.0      ! lower right values
ul = 0.0      ! upper left values
ur = 0.0      ! upper right values

! ----------------------------------------------------------------------
! Section 3.  Fill arrays
! ----------------------------------------------------------------------
! Each array ll, lr, ul, ur is filled by only 1 pe
! Columns are filled in reverse order (North to South)

DO k = 1, num_lev

  !   ll is lower left corner
  DO j = jss, jse
    ji = jse - j + 1 + js
    DO i = iws, iwe
      ii = i - iws + 1 + iw
      ll(ii,ji,k) = field(i,j,lev_in(k))
    END DO  !  i = iws, iwe
  END DO ! j = jss, jse

  !   lr is lower right corner
  DO j = jss, jse
    ji = jse - j + 1 + js
    DO i = ies, iee
      ii = i - ies + 1 + ie
      lr(ii,ji,k) = field(i,j,lev_in(k))
    END DO  !  i = ies, iee
  END DO ! j = jss, jse

  !   ul is upper left corner
  DO j = jns, jne
    ji = jne - j + 1 + jn
    DO i = iws, iwe
      ii = i - iws + 1 + iw
      ul(ii,ji,k) = field(i,j,lev_in(k))
    END DO  !  i = iws, iwe
  END DO ! j = jns, jne

  !   ur is upper right corner
  DO j = jns, jne
    ji = jne - j + 1 + jn
    DO i = ies, iee
      ii = i - ies + 1 + ie
      ur(ii,ji,k) = field(i,j,lev_in(k))
    END DO  !  i = ies, iee
  END DO ! j = jns, jne

END DO  !  k = 1, num_lev

! Each element of ll, lr, ul, ur is filled by only 1 pe

! Sum over pes to put same value on all processors (hence pe0)
! Length = size = max_x*max_y*num_lev arrays

SIZE = 4 * 4 * num_lev

CALL gc_rsum(SIZE, n_proc, info, ll)
CALL gc_rsum(SIZE, n_proc, info, lr)
CALL gc_rsum(SIZE, n_proc, info, ul)
CALL gc_rsum(SIZE, n_proc, info, ur)

! ----------------------------------------------------------------------
! Section 2.  Print corners
! ----------------------------------------------------------------------

DO k = 1, num_lev

  WRITE(umMessage,'(''level '',I3,'' upper left'')') lev_in(k)
  CALL umPrint(umMessage,src='print_corners')
  ii = dom_w
  IF (ix == 4) THEN
    WRITE(umMessage,'(''points             '',I4'//                          &
        ',''                     '',I4,''                   '',I4'// &
        ',''                    '',I4)') ii, ii+1, ii+2, ii+3
    CALL umPrint(umMessage,src='print_corners')
  ELSE IF (ix == 3) THEN
    WRITE(umMessage,'(''points             '',I4'//&
        ',''                     '',I4,''                   '',I4)')&
        ii, ii+1, ii+2
    CALL umPrint(umMessage,src='print_corners')
  ELSE IF (ix == 2) THEN
    WRITE(umMessage,'(''points             '',I4'//&
        ',''                     '',I4)') ii, ii+1
    CALL umPrint(umMessage,src='print_corners')
  ELSE
    WRITE(umMessage,'(''points             '',I4)')ii
    CALL umPrint(umMessage,src='print_corners')
  END IF !  ix == 4
  DO j = 1, jy
    jn = dom_n - j + 1
    WRITE(umMessage,'(''row '',I4,4E24.16)') jn, (ul(i,j,k), i = 1, ix)
    CALL umPrint(umMessage,src='print_corners')
  END DO
  WRITE(umMessage,'(''level '',I3,'' upper right'')') lev_in(k)
  CALL umPrint(umMessage,src='print_corners')
  ii = dom_e - ix + 1
  IF (ix == 4) THEN
    WRITE(umMessage,'(''points             '',I4'//&
        ',''                     '',I4,''                   '',I4'//&
        ',''                    '',I4)') ii, ii+1, ii+2, ii+3
    CALL umPrint(umMessage,src='print_corners')
  ELSE IF (ix == 3) THEN
    WRITE(umMessage,'(''points             '',I4'//&
        ',''                     '',I4,''                   '',I4)')&
        ii, ii+1, ii+2
    CALL umPrint(umMessage,src='print_corners')
  ELSE IF (ix == 2) THEN
    WRITE(umMessage,'(''points             '',I4'//&
        ',''                     '',I4)') ii, ii+1
    CALL umPrint(umMessage,src='print_corners')
  ELSE
    WRITE(umMessage,'(''points             '',I4)')ii
    CALL umPrint(umMessage,src='print_corners')
  END IF !  ix == 4
  DO j = 1, jy
    jn = dom_n - j + 1
    WRITE(umMessage,'(''row '',I4,4E24.16)') jn, (ur(i,j,k), i = 1, ix)
    CALL umPrint(umMessage,src='print_corners')
  END DO
  WRITE(umMessage,'(''level '',I3,'' lower left'')') lev_in(k)
  CALL umPrint(umMessage,src='print_corners')
  ii = dom_w
  IF (ix == 4) THEN
    WRITE(umMessage,'(''points             '',I4'//&
        ',''                     '',I4,''                   '',I4'//&
        ',''                    '',I4)') ii, ii+1, ii+2, ii+3
    CALL umPrint(umMessage,src='print_corners')
  ELSE IF (ix == 3) THEN
    WRITE(umMessage,'(''points             '',I4'//&
        ',''                     '',I4,''                   '',I4)')&
        ii, ii+1, ii+2
    CALL umPrint(umMessage,src='print_corners')
  ELSE IF (ix == 2) THEN
    WRITE(umMessage,'(''points             '',I4'//&
        ',''                     '',I4)') ii, ii+1
    CALL umPrint(umMessage,src='print_corners')
  ELSE
    WRITE(umMessage,'(''points             '',I4)')ii
    CALL umPrint(umMessage,src='print_corners')
  END IF !  ix == 4
  DO j = 1, jy
    jn = dom_s + jy - j
    WRITE(umMessage,'(''row '',I4,4E24.16)') jn, (ll(i,j,k), i = 1, ix)
    CALL umPrint(umMessage,src='print_corners')
  END DO
  WRITE(umMessage,'(''level '',I3,'' lower right'')') lev_in(k)
  CALL umPrint(umMessage,src='print_corners')
  ii = dom_e - ix + 1
  IF (ix == 4) THEN
    WRITE(umMessage,'(''points             '',I4'//&
        ',''                     '',I4,''                   '',I4'//&
        ',''                    '',I4)') ii, ii+1, ii+2, ii+3
    CALL umPrint(umMessage,src='print_corners')
  ELSE IF (ix == 3) THEN
    WRITE(umMessage,'(''points             '',I4'//&
        ',''                     '',I4,''                   '',I4)')&
        ii, ii+1, ii+2
    CALL umPrint(umMessage,src='print_corners')
  ELSE IF (ix == 2) THEN
    WRITE(umMessage,'(''points             '',I4'//&
        ',''                     '',I4)') ii, ii+1
    CALL umPrint(umMessage,src='print_corners')
  ELSE
    WRITE(umMessage,'(''points             '',I4)')ii
    CALL umPrint(umMessage,src='print_corners')
  END IF !  ix == 4
  DO j = 1, jy
    jn = dom_s + jy - j
    WRITE(umMessage,'(''row '',I4,4E24.16)') jn, (lr(i,j,k), i = 1, ix)
    CALL umPrint(umMessage,src='print_corners')
  END DO

END DO !  k = 1, num_lev

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE print_corners


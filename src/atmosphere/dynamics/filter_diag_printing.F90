! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE filter_diag_printing_mod
USE umPrintMgr, ONLY: umPrint, umMessage, printstatus, prstatus_diag
USE um_parcore, ONLY: mype, nproc
USE fields_rhs_mod
USE atm_fields_mod
USE atm_fields_bounds_mod
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

PUBLIC :: filter_diag_printing1,filter_diag_printing2,&
          filter_diag_printing3,filter_diag_printing4,&
          filter_diag_printing5
PRIVATE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FILTER_DIAG_PRINTING_MOD'

CONTAINS

SUBROUTINE filter_diag_printing1()


IMPLICIT NONE

!
! Description:
!
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


INTEGER :: ierr

REAL :: max_r_u,min_r_u ! local diagnostics
REAL :: max_r_v,min_r_v
REAL :: max_u, min_u
REAL :: max_v, min_v
REAL :: max_w, min_w
REAL :: max_thetav, min_thetav

REAL :: max_r_w, min_r_w
REAL :: max_r_theta, min_r_theta
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILTER_DIAG_PRINTING1'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF (printstatus >= prstatus_diag ) THEN
  WRITE(umMessage,FMT='(A)') 'Calling polar filter routine'
  CALL umPrint(umMessage,src='filter_diag_printing',pe=0)

  max_u = MAXVAL(u(udims%i_start:udims%i_end,                     &
                   udims%j_start:udims%j_end,                     &
                   udims%k_start:udims%k_end))
  min_u = MINVAL(u(udims%i_start:udims%i_end,                     &
                   udims%j_start:udims%j_end,                     &
                   udims%k_start:udims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_u)
  CALL gc_rmin(1,nproc,ierr,min_u)

  max_v = MAXVAL(v(vdims%i_start:vdims%i_end,                     &
                   vdims%j_start:vdims%j_end,                     &
                   vdims%k_start:vdims%k_end))
  min_v = MINVAL(v(vdims%i_start:vdims%i_end,                     &
                   vdims%j_start:vdims%j_end,                     &
                   vdims%k_start:vdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_v)
  CALL gc_rmin(1,nproc,ierr,min_v)

  max_w = MAXVAL(w(wdims%i_start:wdims%i_end,                     &
                   wdims%j_start:wdims%j_end,                     &
                   wdims%k_start:wdims%k_end))
  min_w = MINVAL(w(wdims%i_start:wdims%i_end,                     &
                   wdims%j_start:wdims%j_end,                     &
                   wdims%k_start:wdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_w)
  CALL gc_rmin(1,nproc,ierr,min_w)

  max_thetav = MAXVAL(thetav(tdims%i_start:tdims%i_end,           &
                             tdims%j_start:tdims%j_end,           &
                             tdims%k_start:tdims%k_end))
  min_thetav = MINVAL(thetav(tdims%i_start:tdims%i_end,           &
                             tdims%j_start:tdims%j_end,           &
                             tdims%k_start:tdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_thetav)
  CALL gc_rmin(1,nproc,ierr,min_thetav)


  IF (mype == 0) THEN
    WRITE(umMessage,FMT='(2A)') '==========================',&
                         '========================='
    CALL umPrint(umMessage,src='filter_diag_printing1')
    WRITE(umMessage,FMT='(A)') 'Calling eg_NI_filter_Ctl'
    CALL umPrint(umMessage,src='filter_diag_printing1')
    WRITE(umMessage,FMT='(A)') ' Max/Min before polar filter :'
    CALL umPrint(umMessage,src='filter_diag_printing1')
    WRITE(umMessage,FMT='(A,2E25.10)') '    u = ',max_u, min_u
    CALL umPrint(umMessage,src='filter_diag_printing1')
    WRITE(umMessage,FMT='(A,2E25.10)') '    v = ',max_v, min_v
    CALL umPrint(umMessage,src='filter_diag_printing1')
    WRITE(umMessage,FMT='(A,2E25.10)') '    w = ',max_w, min_w
    CALL umPrint(umMessage,src='filter_diag_printing1')
    WRITE(umMessage,FMT='(A,2E25.10)') 'theta = ',max_thetav, min_thetav
    CALL umPrint(umMessage,src='filter_diag_printing1')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE

SUBROUTINE filter_diag_printing2()

IMPLICIT NONE


INTEGER :: ierr

REAL :: max_r_u,min_r_u ! local diagnostics
REAL :: max_r_v,min_r_v
REAL :: max_u, min_u
REAL :: max_v, min_v
REAL :: max_w, min_w
REAL :: max_thetav, min_thetav

REAL :: max_r_w, min_r_w
REAL :: max_r_theta, min_r_theta
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILTER_DIAG_PRINTING2'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF (printstatus >= prstatus_diag) THEN
  max_u = MAXVAL(u(udims%i_start:udims%i_end,                   &
                   udims%j_start:udims%j_end,                   &
                   udims%k_start:udims%k_end))
  min_u = MINVAL(u(udims%i_start:udims%i_end,                   &
                   udims%j_start:udims%j_end,                   &
                   udims%k_start:udims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_u)
  CALL gc_rmin(1,nproc,ierr,min_u)

  max_v = MAXVAL(v(vdims%i_start:vdims%i_end,                     &
                   vdims%j_start:vdims%j_end,                     &
                   vdims%k_start:vdims%k_end))
  min_v = MINVAL(v(vdims%i_start:vdims%i_end,                     &
                   vdims%j_start:vdims%j_end,                     &
                   vdims%k_start:vdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_v)
  CALL gc_rmin(1,nproc,ierr,min_v)

  max_w = MAXVAL(w(wdims%i_start:wdims%i_end,                     &
                   wdims%j_start:wdims%j_end,                     &
                   wdims%k_start:wdims%k_end))
  min_w = MINVAL(w(wdims%i_start:wdims%i_end,                     &
                   wdims%j_start:wdims%j_end,                     &
                   wdims%k_start:wdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_w)
  CALL gc_rmin(1,nproc,ierr,min_w)

  max_thetav = MAXVAL(thetav(tdims%i_start:tdims%i_end,           &
                             tdims%j_start:tdims%j_end,           &
                             tdims%k_start:tdims%k_end))
  min_thetav = MINVAL(thetav(tdims%i_start:tdims%i_end,           &
                             tdims%j_start:tdims%j_end,           &
                             tdims%k_start:tdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_thetav)
  CALL gc_rmin(1,nproc,ierr,min_thetav)


  IF (mype == 0) THEN
    WRITE(umMessage,FMT='(A)') ' '
    CALL umPrint(umMessage,src='filter_diag_printing2')
    WRITE(umMessage,FMT='(2A)') '=================================',&
                       '=================='
    CALL umPrint(umMessage,src='filter_diag_printing2')
    WRITE(umMessage,FMT='(A)') ' Max/Min after polar filter :'
    CALL umPrint(umMessage,src='filter_diag_printing2')
    WRITE(umMessage,FMT='(A,2E25.10)') '    u = ',max_u, min_u
    CALL umPrint(umMessage,src='filter_diag_printing2')
    WRITE(umMessage,FMT='(A,2E25.10)') '    v = ',max_v, min_v
    CALL umPrint(umMessage,src='filter_diag_printing2')
    WRITE(umMessage,FMT='(A,2E25.10)') '    w = ',max_w, min_w
    CALL umPrint(umMessage,src='filter_diag_printing2')
    WRITE(umMessage,FMT='(A,2E25.10)') 'theta = ',max_thetav, min_thetav
    CALL umPrint(umMessage,src='filter_diag_printing2')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE

SUBROUTINE filter_diag_printing3()

IMPLICIT NONE

INTEGER :: ierr

REAL :: max_r_u,min_r_u ! local diagnostics
REAL :: max_r_v,min_r_v
REAL :: max_u, min_u
REAL :: max_v, min_v
REAL :: max_w, min_w
REAL :: max_thetav, min_thetav

REAL :: max_r_w, min_r_w
REAL :: max_r_theta, min_r_theta
REAL :: max_r_rho, min_r_rho
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILTER_DIAG_PRINTING3'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF (printstatus >= prstatus_diag) THEN
  max_r_theta     = MAXVAL(r_theta(tdims%i_start:tdims%i_end,     &
                                   tdims%j_start:tdims%j_end,     &
                                   tdims%k_start:tdims%k_end))
  min_r_theta     = MINVAL(r_theta(tdims%i_start:tdims%i_end,     &
                                   tdims%j_start:tdims%j_end,     &
                                   tdims%k_start:tdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_theta)
  CALL gc_rmin(1,nproc,ierr,min_r_theta)

  max_r_u     = MAXVAL(r_u(udims%i_start:udims%i_end,             &
                           udims%j_start:udims%j_end,             &
                           udims%k_start:udims%k_end))
  min_r_u     = MINVAL(r_u(udims%i_start:tdims%i_end,             &
                           udims%j_start:udims%j_end,             &
                           udims%k_start:udims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_u)
  CALL gc_rmin(1,nproc,ierr,min_r_u)

  max_r_v     = MAXVAL(r_v(vdims%i_start:vdims%i_end,             &
                           vdims%j_start:vdims%j_end,             &
                           vdims%k_start:vdims%k_end))
  min_r_v     = MINVAL(r_v(vdims%i_start:vdims%i_end,             &
                           vdims%j_start:vdims%j_end,             &
                           vdims%k_start:vdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_v)
  CALL gc_rmin(1,nproc,ierr,min_r_v)

  max_r_w     = MAXVAL(r_w(wdims%i_start:wdims%i_end,             &
                           wdims%j_start:wdims%j_end,             &
                           wdims%k_start:wdims%k_end))
  min_r_w     = MINVAL(r_w(wdims%i_start:wdims%i_end,             &
                           wdims%j_start:wdims%j_end,             &
                           wdims%k_start:wdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_w)
  CALL gc_rmin(1,nproc,ierr,min_r_w)

  max_r_rho     = MAXVAL(r_rho(pdims%i_start:pdims%i_end,         &
                               pdims%j_start:pdims%j_end,         &
                               pdims%k_start:pdims%k_end))
  min_r_rho     = MINVAL(r_rho(pdims%i_start:pdims%i_end,         &
                               pdims%j_start:pdims%j_end,         &
                               pdims%k_start:pdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_rho)
  CALL gc_rmin(1,nproc,ierr,min_r_rho)

  IF (mype == 0) THEN
    WRITE(umMessage,FMT='(2A)') '=================================',&
                        '=================='
    CALL umPrint(umMessage,src='filter_diag_printing3')
    WRITE(umMessage,FMT='(A)') 'Calling eg_NI_filter_incs_Ctl'
    CALL umPrint(umMessage,src='filter_diag_printing3')
    WRITE(umMessage,FMT='(A)') ' Max/Min before polar filter :'
    CALL umPrint(umMessage,src='filter_diag_printing3')
    WRITE(umMessage,FMT='(A,2E25.10)') '    R_u = ',MAX_r_u,MIN_r_u
    CALL umPrint(umMessage,src='filter_diag_printing3')
    WRITE(umMessage,FMT='(A,2E25.10)') '    R_v = ',MAX_r_v,MIN_r_v
    CALL umPrint(umMessage,src='filter_diag_printing3')
    WRITE(umMessage,FMT='(A,2E25.10)') '    R_w = ',MAX_r_w,MIN_r_w
    CALL umPrint(umMessage,src='filter_diag_printing3')
    WRITE(umMessage,FMT='(A,2E25.10)') 'R_theta = ',MAX_r_theta,MIN_r_theta
    CALL umPrint(umMessage,src='filter_diag_printing3')
    WRITE(umMessage,FMT='(A,2E25.10)') '  R_rho = ',MAX_r_rho,MIN_r_rho
    CALL umPrint(umMessage,src='filter_diag_printing3')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE

SUBROUTINE filter_diag_printing4()

IMPLICIT NONE

INTEGER :: ierr

REAL :: max_r_u,min_r_u ! local diagnostics
REAL :: max_r_v,min_r_v
REAL :: max_u, min_u
REAL :: max_v, min_v
REAL :: max_w, min_w
REAL :: max_thetav, min_thetav

REAL :: max_r_w, min_r_w
REAL :: max_r_theta, min_r_theta
REAL :: max_r_rho, min_r_rho

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILTER_DIAG_PRINTING4'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF (printstatus >= prstatus_diag) THEN
  max_r_theta     = MAXVAL(r_theta(tdims%i_start:tdims%i_end,     &
                                   tdims%j_start:tdims%j_end,     &
                                   tdims%k_start:tdims%k_end))
  min_r_theta     = MINVAL(r_theta(tdims%i_start:tdims%i_end,     &
                                   tdims%j_start:tdims%j_end,     &
                                   tdims%k_start:tdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_theta)
  CALL gc_rmin(1,nproc,ierr,min_r_theta)

  max_r_u     = MAXVAL(r_u(udims%i_start:udims%i_end,             &
                           udims%j_start:udims%j_end,             &
                           udims%k_start:udims%k_end))
  min_r_u     = MINVAL(r_u(udims%i_start:tdims%i_end,             &
                           udims%j_start:udims%j_end,             &
                           udims%k_start:udims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_u)
  CALL gc_rmin(1,nproc,ierr,min_r_u)

  max_r_v     = MAXVAL(r_v(vdims%i_start:vdims%i_end,             &
                           vdims%j_start:vdims%j_end,             &
                           vdims%k_start:vdims%k_end))
  min_r_v     = MINVAL(r_v(vdims%i_start:vdims%i_end,             &
                           vdims%j_start:vdims%j_end,             &
                           vdims%k_start:vdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_v)
  CALL gc_rmin(1,nproc,ierr,min_r_v)

  max_r_w     = MAXVAL(r_w(wdims%i_start:wdims%i_end,             &
                           wdims%j_start:wdims%j_end,             &
                           wdims%k_start:wdims%k_end))
  min_r_w     = MINVAL(r_w(wdims%i_start:wdims%i_end,             &
                           wdims%j_start:wdims%j_end,             &
                           wdims%k_start:wdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_w)
  CALL gc_rmin(1,nproc,ierr,min_r_w)

  max_r_rho     = MAXVAL(r_rho(pdims%i_start:pdims%i_end,         &
                               pdims%j_start:pdims%j_end,         &
                               pdims%k_start:pdims%k_end))
  min_r_rho     = MINVAL(r_rho(pdims%i_start:pdims%i_end,         &
                               pdims%j_start:pdims%j_end,         &
                               pdims%k_start:pdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_rho)
  CALL gc_rmin(1,nproc,ierr,min_r_rho)

  IF (mype == 0) THEN
    WRITE(umMessage,FMT='(A)') ' '
    CALL umPrint(umMessage,src='filter_diag_printing4')
    WRITE(umMessage,FMT='(A)') ' Max/Min after polar incs filter :'
    CALL umPrint(umMessage,src='filter_diag_printing4')
    WRITE(umMessage,FMT='(A,2E25.10)') '    R_u = ',MAX_r_u,MIN_r_u
    CALL umPrint(umMessage,src='filter_diag_printing4')
    WRITE(umMessage,FMT='(A,2E25.10)') '    R_v = ',MAX_r_v,MIN_r_v
    CALL umPrint(umMessage,src='filter_diag_printing4')
    WRITE(umMessage,FMT='(A,2E25.10)') '    R_w = ',MAX_r_w,MIN_r_w
    CALL umPrint(umMessage,src='filter_diag_printing4')
    WRITE(umMessage,FMT='(A,2E25.10)') 'R_theta = ',MAX_r_theta,MIN_r_theta
    CALL umPrint(umMessage,src='filter_diag_printing4')
    WRITE(umMessage,FMT='(A,2E25.10)') '  R_rho = ',MAX_r_rho,MIN_r_rho
    CALL umPrint(umMessage,src='filter_diag_printing4')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE

SUBROUTINE filter_diag_printing5()

IMPLICIT NONE

INTEGER :: ierr

REAL :: max_r_u,min_r_u ! local diagnostics
REAL :: max_r_v,min_r_v
REAL :: max_u, min_u
REAL :: max_v, min_v
REAL :: max_w, min_w
REAL :: max_thetav, min_thetav

REAL :: max_r_w, min_r_w
REAL :: max_r_theta, min_r_theta

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILTER_DIAG_PRINTING5'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF (printstatus >= prstatus_diag) THEN
  max_r_theta     = MAXVAL(s_thetav(tdims%i_start:tdims%i_end,    &
                                    tdims%j_start:tdims%j_end,    &
                                    tdims%k_start:tdims%k_end))
  min_r_theta     = MINVAL(s_thetav(tdims%i_start:tdims%i_end,    &
                                    tdims%j_start:tdims%j_end,    &
                                    tdims%k_start:tdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_theta)
  CALL gc_rmin(1,nproc,ierr,min_r_theta)

  max_r_u     = MAXVAL(s_u(udims%i_start:udims%i_end,             &
                           udims%j_start:udims%j_end,             &
                           udims%k_start:udims%k_end))
  min_r_u     = MINVAL(s_u(udims%i_start:tdims%i_end,             &
                           udims%j_start:udims%j_end,             &
                           udims%k_start:udims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_u)
  CALL gc_rmin(1,nproc,ierr,min_r_u)

  max_r_v     = MAXVAL(s_v(vdims%i_start:vdims%i_end,             &
                           vdims%j_start:vdims%j_end,             &
                           vdims%k_start:vdims%k_end))
  min_r_v     = MINVAL(s_v(vdims%i_start:vdims%i_end,             &
                           vdims%j_start:vdims%j_end,             &
                           vdims%k_start:vdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_v)
  CALL gc_rmin(1,nproc,ierr,min_r_v)

  max_r_w     = MAXVAL(s_w(wdims%i_start:wdims%i_end,             &
                           wdims%j_start:wdims%j_end,             &
                           wdims%k_start:wdims%k_end))
  min_r_w     = MINVAL(s_w(wdims%i_start:wdims%i_end,             &
                           wdims%j_start:wdims%j_end,             &
                           wdims%k_start:wdims%k_end))
  CALL gc_rmax(1,nproc,ierr,max_r_w)
  CALL gc_rmin(1,nproc,ierr,min_r_w)

  IF (mype == 0) THEN
    WRITE(umMessage,FMT='(2A)') '===============================',&
                        '===================='
    CALL umPrint(umMessage,src='filter_diag_printing5')
    WRITE(umMessage,FMT='(A)') 'Calling eg_NI_filter_incs_Ctl'
    CALL umPrint(umMessage,src='filter_diag_printing5')
    WRITE(umMessage,FMT='(A)') ' Max/Min before/after polar filter :'
    CALL umPrint(umMessage,src='filter_diag_printing5')
    WRITE(umMessage,FMT='(A,2E25.10)') '    S_u = ',MAX_r_u,MIN_r_u
    CALL umPrint(umMessage,src='filter_diag_printing5')
    WRITE(umMessage,FMT='(A,2E25.10)') '    S_v = ',MAX_r_v,MIN_r_v
    CALL umPrint(umMessage,src='filter_diag_printing5')
    WRITE(umMessage,FMT='(A,2E25.10)') '    S_w = ',MAX_r_w,MIN_r_w
    CALL umPrint(umMessage,src='filter_diag_printing5')
    WRITE(umMessage,FMT='(A,2E25.10)') 'S_theta = ',MAX_r_theta,MIN_r_theta
    CALL umPrint(umMessage,src='filter_diag_printing5')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE

END MODULE

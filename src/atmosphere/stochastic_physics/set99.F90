! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE set99_mod

USE umPrintMgr, ONLY:                                                   &
    umPrint, umMessage
IMPLICIT NONE


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SET99_MOD'

CONTAINS

SUBROUTINE set99(trigs,ifax,n)

!  computes factors of n & trigonometric
!  functions required by fourier

USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport

USE c_skeb2_mod, ONLY: nblock, nfacts
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE


! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


! nfacts and nblock values

INTEGER :: ifax(nfacts), jfax(nfacts), n, nhl, k, nu, nfax, ifac, i
INTEGER :: lfax(1:7) = (/6,8,5,4,3,2,1/)
REAL    :: trigs(n), del, angle

INTEGER                       ::  icode
CHARACTER (LEN=errormessagelength)            ::  cmessage
CHARACTER (LEN=* ), PARAMETER ::  RoutineName='SET99'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
del=4.0*ASIN(1.0)/REAL(n)
nhl=(n/2)-1
DO k=0,nhl
  angle=REAL(k)*del
  trigs(2*k+1)=COS(angle)
  trigs(2*k+2)=SIN(angle)
END DO

!     find factors of n (8,6,5,4,3,2; only one 8 allowed)
!     look for the eight first and store factors in descending order

nu=n
nfax=1
IF ( MOD(nu,8) == 0 ) THEN
  jfax(1)=8
  nu=nu/8
  nfax=2
END IF

ifac=6
DO WHILE ( nu /= 1 .AND. nfax < nfacts)
  IF ( MOD(nu,ifac) == 0 ) THEN
    jfax(nfax)=ifac
    nu=nu/ifac
    nfax=nfax+1
  ELSE
    ifac=ifac-1
    IF ( ifac == 1 ) THEN
      WRITE(umMessage,'(A4,I4,A27)')'1n =',n,' - contains illegal factors'
      CALL umPrint(umMessage,src='set99')
      cmessage = 'Too may factors'
      icode    = 1
      CALL EReport(RoutineName,icode,cmessage)
    END IF
  END IF
END DO
IF ( nfax > nfacts-1 ) THEN
  WRITE(umMessage,'(A4,I4,A27)')'1n =',n,' - contains illegal factors'
  CALL umPrint(umMessage,src='set99')
  cmessage = 'Too may factors'
  icode    = 1
  CALL EReport(RoutineName,icode,cmessage)
END IF

!     now reverse order of factors

nfax=nfax-1
ifax(1)=nfax
DO i=1,nfax
  ifax(nfax+2-i)=jfax(i)
END DO
ifax(nfacts)=n
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE set99

END MODULE set99_mod

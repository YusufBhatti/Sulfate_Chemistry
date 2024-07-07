! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Constructs STASH index array
!
! Subroutine Interface:

SUBROUTINE sindx(nrecs)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE version_mod, ONLY:                                            &
    nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
    nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
    nlevp, npslevp, npslistp
USE submodel_mod, ONLY: n_internal_model_max
USE stextend_mod, ONLY: indx_s, list_s
USE stparam_mod, ONLY: st_model_code, st_sect_no_code, st_item_code

IMPLICIT NONE

! Description:
!   The STASH list (LIST_S) is input to this routine. The output from
!   the routine is the STASH index (INDX_S) consistent with LIST_S.
!   Called by STPROC, DUPLIC.
!
! Method:
!   Before this routine is executed, LIST_S has been ordered by
!   model, section, item, input code. In this routine, INDX_S is
!   set up as follows:
!   INDX_S(1,m,s,i)= position of 1st. occurrence of m,s,i in LIST_S
!   INDX_S(2,m,s,i)= no. of occurrences of m,s,i in LIST_S
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards UMDP3
!
!  Subroutine Arguments:

!    Scalar Argument with intent(in):

INTEGER :: nrecs                 ! No. of records in LIST_S

!  Local variables:

INTEGER :: ii
INTEGER :: iitm
INTEGER :: im
INTEGER :: irec
INTEGER :: IS
INTEGER :: isec
INTEGER :: lstart
INTEGER :: modl

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SINDX'

!  Local arrays:

!- End of Header -------------------------------------------------------


! Initialise Index array

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO im    = 1,n_internal_model_max
  DO IS  = 0,nsectp
    DO ii= 1,nitemp
      indx_s(1,im,IS,ii)=0
      indx_s(2,im,IS,ii)=0
    END DO
  END DO
END DO

! Set up index

IF (nrecs >= 1) THEN

  modl=list_s(st_model_code  ,1)
  isec=list_s(st_sect_no_code,1)
  iitm=list_s(st_item_code   ,1)

  lstart=1

  indx_s(2,modl,isec,iitm) = 1
  indx_s(1,modl,isec,iitm) = 1

  IF (nrecs >= 2) THEN            ! More than one record in LIST_S

    DO irec=2,nrecs
      IF ((list_s(st_model_code  ,irec) == modl) .AND.                &
                                                     ! Same model,
         (list_s(st_sect_no_code,irec) == isec) .AND.                &
                                                     ! sec,item,
         (list_s(st_item_code   ,irec) == iitm)) THEN! as before

        lstart=lstart+1

        indx_s(2,modl,isec,iitm) = indx_s(2,modl,isec,iitm)+1

      ELSE   ! New model, section, item

        modl=list_s(st_model_code  ,irec)
        isec=list_s(st_sect_no_code,irec)
        iitm=list_s(st_item_code   ,irec)

        lstart=lstart+1

        indx_s(1,modl,isec,iitm) = lstart
        indx_s(2,modl,isec,iitm) = 1

      END IF
    END DO

  ELSE      ! Only one record

    indx_s(1,list_s(st_model_code  ,1),                             &
             list_s(st_sect_no_code,1),                             &
             list_s(st_item_code   ,1)) = 1

    indx_s(2,list_s(st_model_code  ,1),                             &
             list_s(st_sect_no_code,1),                             &
             list_s(st_item_code   ,1)) = 1

  END IF
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sindx

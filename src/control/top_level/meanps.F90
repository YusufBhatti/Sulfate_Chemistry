! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!      Subroutine:
!      MEANPS
!
!      Purpose:
!      To mean partial sums. Sums obtained from the D1 array, put there
!      by ACUMPS.
!
!      Programming standard:
!      UM Doc Paper 3
!
!      External documentation:
!      On-line UM document C5 - Control of means calculations
!
!      Interface and arguments:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE meanps(                                                &
  n_objs_d1,d1_addr                                               &
  ,len_data,d1,                                                   &
  meaning_period                                                  &
  )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE d1_array_mod, ONLY: d1_list_len, d1_imodl, d1_address,        &
                        d1_length
USE stash_array_mod, ONLY: totitems, stlist
USE stparam_mod, ONLY: st_macrotag, st_d1pos, s_modl
USE missing_data_mod, ONLY: rmdi  

IMPLICIT NONE

INTEGER ::                                                        &
  n_objs_d1

INTEGER ::                                                        &
  d1_addr(d1_list_len,n_objs_d1)

INTEGER ::                                                        &
  len_data,                                                       &
                          ! IN Length of model data
  meaning_period          ! IN Meaning period (in multiples
                          !             of restart frequency)
!
REAL ::                                                           &
  d1(len_data)            ! IN/OUT Real equivalence of data block
                          !    containing meaned fields
!
!      Local variables
!
INTEGER ::                                                        &
  j,k                                                             &
                          ! Loop indices
  ,tag                                                            &
                          ! Climate mean tag
  ,address                                                        &
  ,ptd1
!
REAL ::                                                           &
  factor                                                          &
                          ! Meaning period (real)
  ,rfactor                 ! Reciprocal of FACTOR
!
!      Constants
!
REAL ::                                                           &
       one                  ! 1.0

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MEANPS'


! Calculate divisor
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
one=1.0
factor=meaning_period
rfactor=one/factor

!----------------------------------------------------------------------
!     Loop through STASH list and process climate mean fields
!     NOTE: D1 contains partial sums put there by preceding ACUMPS call
!----------------------------------------------------------------------

DO k=1,totitems
  tag=stlist(st_macrotag,k)/1000
  ptd1=stlist(st_d1pos,k)
  IF (tag /= 0 .AND. stlist(s_modl,k) == d1_addr(d1_imodl,ptd1)) THEN
    ! Object is tagged for climate meaning and in relevant internal model
    address=d1_addr(d1_address,ptd1) ! local address
    ! Divide whole field by FACTOR - except for RMDI
    DO j=address,address+d1_addr(d1_length,ptd1)-1
      IF (d1(j) /= rmdi) THEN
        d1(j)=d1(j)*rfactor
      END IF
    END DO
  END IF
END DO

! *********************************************************************
!     End of loop over STASH list
! *********************************************************************

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE meanps

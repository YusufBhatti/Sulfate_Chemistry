! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine print_array_max

SUBROUTINE print_array_max(                                       &
                           field,                                 &
                           row_length, rows, levels,              &
                           halo_i, halo_j,                        &
                           start_level, end_level )

USE um_parvars, ONLY: datastart, nproc
USE um_parcore, ONLY: mype
USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows

! Purpose:
!          Print value and location of maximum value of an array
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE


INTEGER , INTENT(IN) ::                                           &
 row_length,                                                      &
                 ! number of point on a row.
 rows,                                                            &
                 ! number of rows.
 levels,                                                          &
                 ! number of levels in field
 start_level,                                                     &
                 ! start level for finding max
 end_level,                                                       &
                 ! END level for finding max
 halo_i,                                                          &
              ! Size of halo in i direction.
 halo_j,      ! Size of halo in j direction.

REAL , INTENT(IN) ::                                              &
       field (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,   &
              levels)

! Local Variables.

INTEGER ::                                                        &
  k, ki, kr,                                                      &
            ! Loop counters
  num_levels,                                                     &
  gi, gj,                                                         &
            ! global pointers
  pe,                                                             &
            ! processor identifier
  info

REAL ::                                                           &
  sumr(2,levels),                                                 &
                      ! array  for summing integers over pe's
  max_pe(levels),                                                 &
                      ! max wind each level on this pe
  max_all(levels)       ! array  for finding max over pe's

INTEGER ::                                                        &
  sumi(3,levels),                                                 &
                    ! array  for summing integers over pe's
  i_max(levels),                                                  &
                    ! i pointers for prints
  j_max(levels)       ! j pointers for prints

REAL ::                                                           &
  recip_row_length,                                               &
  recip_rows,                                                     &
  lambda,                                                         &
  phi

INTEGER :: max_indices(2)

CHARACTER(LEN=8) :: l_string, p_string
PARAMETER (l_string ='% East')
PARAMETER (p_string ='% North')

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_ARRAY_MAX'

! ----------------------------------------------------------------------
! Section 0. Initialise values
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

recip_row_length = 1.0 / REAL(global_row_length)
recip_rows = 1.0 / REAL(global_rows)
num_levels = END_level - start_level + 1

DO k = 1, num_levels
  max_indices = MAXLOC(field(1:row_length,1:rows,k))
  i_max(k) = max_indices(1)
  j_max(k) = max_indices(2)
  max_pe(k) = field(max_indices(1),max_indices(2),k)
  max_all(k) = max_pe(k)
END DO !  k = 1, num_levels

CALL gc_rmax(num_levels, nproc, info, max_all)

sumi = 0
sumr = 0.0
ki = 0
kr = 0
DO k = 1, num_levels
  ki = ki + 3
  kr = kr + 2
  IF ( max_pe(k) >= max_all(k) ) THEN
    sumi(1,k) = mype ! max is on this processor for this level
    gi = datastart(1) + i_max(k) - 1
    gj = datastart(2) + j_max(k) - 1
    sumi(2,k) = gi
    sumi(3,k) = gj
    sumr(1,k) = (gi-1) * recip_row_length * 100.0
    sumr(2,k) = (gj-1) * recip_rows * 100.0
  END IF ! max_pe(k) >= max_all(k)
END DO !  k = 1, num_levels

CALL gc_isum(ki, nproc, info, sumi)
CALL gc_rsumr(kr, nproc, info, sumr)

WRITE(umMessage,*) '                         full domain position'
CALL umPrint(umMessage,src='print_array_max')
WRITE(umMessage,*) 'level   max      proc   i , j'
CALL umPrint(umMessage,src='print_array_max')
DO   k = 1, num_levels
  pe = sumi(1,k)
  gi = sumi(2,k)
  gj = sumi(3,k)
  lambda = sumr(1,k)
  phi    = sumr(2,k)
  WRITE(umMessage,'(I4, E12.4, 3I5, F6.1, A7, F6.1, A7)')                 &
                    k, max_all(k), pe, gi, gj,                    &
                    lambda, l_string, phi, p_string
  CALL umPrint(umMessage,src='print_array_max')
END DO   !  k = 1, num_levels

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE print_array_max

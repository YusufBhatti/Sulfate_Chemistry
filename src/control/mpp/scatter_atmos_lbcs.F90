! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Passes out atmosphere LBCs to processors at boundaries
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: MPP

! Subroutine Interface
SUBROUTINE scatter_atmos_lbcs(                                    &
  full_lbc,decomp_lbc,                                            &
  full_lbc_size,decomp_lbc_size,                                  &
  full_lbc_levels,decomp_lbc_levels,                              &
  fld_type,halo_type,rim_type,                                    &
  pe_for_level,                                                   &
  icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype, nproc
USE UM_ParParams, ONLY: pnorth, psouth, peast, pwest
USE lbc_mod, ONLY: rimwidtha, global_lbc_starta, g_lbc_starta
USE mpl, ONLY: mpl_real
USE errormessagelength_mod, ONLY: errormessagelength

USE lbc_calc_size_mod, ONLY: lbc_calc_size

IMPLICIT NONE

! Description:
! Scatters atmosphere LBCs to relevant processors at the grid boundaries

! Method:
! I - If master task (0)
! I - Master Task (0) Allocates Buffer
! I  L--Loop over all processors : iproc
! I L L--Loop over sides (North,East,South,West) : iside
! I L L I--IF processor iproc has an extremity at iside THEN
! I L L I  Copy the data to the buffer (data_lbc)
! I L L I--ENDIF
! I L L -- End Loop iside
! I L -- End Loop iproc
! I - ENDIF - master task
! Scatter data ( MPL_Scatterv collective operation )

! Subroutine Arguments:

INTEGER, INTENT(IN) ::  full_lbc_size
                               ! IN single level size of the FULL_LBC array
INTEGER, INTENT(IN) ::  decomp_lbc_size
                               ! IN single level size of the DECOMP_LBC
                               !    array
INTEGER, INTENT(IN) :: full_lbc_levels
                               ! IN number of levels of FULL_LBC on this
                               !    processor
INTEGER, INTENT(IN) ::  decomp_lbc_levels
                               ! IN number of levels of DECOMP_LBC
INTEGER, INTENT(IN) ::  fld_type
                               ! IN Which fld_type is the LBC?
INTEGER, INTENT(IN) ::  halo_type
                               ! IN Which halo_type is the LBC?
INTEGER, INTENT(IN) ::  rim_type
                               ! IN Which rim_type is the LBC?
INTEGER, INTENT(IN) ::  pe_for_level(decomp_lbc_levels)
                               ! IN which level of FULL_LBC is on
                               !    which processor

REAL, INTENT(IN) ::  full_lbc(full_lbc_size,full_lbc_levels)
                               ! IN Some levels of the full LBC
REAL, INTENT(OUT) :: decomp_lbc(decomp_lbc_size,decomp_lbc_levels)
                               ! OUT All levels of the decomposed LBC on
                               !     this processor

INTEGER, INTENT(INOUT) :: icode  ! Return code not used at present.

CHARACTER(LEN=errormessagelength) ::  cmessage ! OUT Error message


! Local variables

INTEGER ::                                                        &
  iproc                                                           &
                     ! loop counter for loop over processors
, iside                                                           &
                     ! loop counter for loop over sides
, k                                                               &
                     ! loop counter for levels
, full_lbc_row_len                                                &
                     ! length of a row of data on the full LBC
, full_lbc_nrows                                                  &
                     ! number of rows of data on the full LBC
, decomp_lbc_row_len                                              &
                     ! length of a row of data on the
                     ! decomposed LBC
, decomp_lbc_nrows                                                &
                     ! number of rows of data on the
                     ! decomposed LBC
, first_lbc_pt                                                    &
                     ! first point in full LBC to start
                     ! copying from
, first_lbc_row                                                   &
                     ! first row in fill LBC to start
                     ! copying from
, level_index_pe(0:nproc-1)                                       &
                     ! How many levels on each PE
, level_index(decomp_lbc_levels)                                  &
                     ! Which level full_LBC corresponds
                     ! to the real level
, full_lbc_start_pt                                               &
                     ! First point index on a level of the
                     ! full LBC to start sending
, decomp_lbc_start_pt                                             &
                     ! First point index on a level of the
                     ! decomposed LBC to start receiving
, info
                     ! GCOM return code

INTEGER :: g_decomp_lbc_size(0:nproc-1)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SCATTER_ATMOS_LBCS'

! New data buffer to fill the decomposition to scatter

INTEGER :: size_data_lbc, data_lbc_start
REAL, ALLOCATABLE :: data_lbc(:)

! Control data needed for the scattering
INTEGER :: recvcount
INTEGER :: displs(0:nproc-1), sendcounts(0:nproc-1)
INTEGER :: my_comm ! communicator

! Other variables

INTEGER :: i, j, index_s, index_d

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! 1.0 Set up indexing describing where each level of full LBC data is
!     held

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
icode = 0
g_decomp_lbc_size=0
g_decomp_lbc_size(mype)=decomp_lbc_size
CALL gc_imax(nproc,nproc,info,g_decomp_lbc_size(0:nproc-1))

level_index_pe(:)=0

! Calculate indexing of full LBC data on the processors.
! A given level "k" of full LBC data will be held on
! processor "PE_FOR_LEVEL(k)" and the index of the level
! on this processor will be "level_index(k)"

DO k=1,decomp_lbc_levels
  level_index_pe(pe_for_level(k))=                                &
    level_index_pe(pe_for_level(k))+1
  level_index(k)=level_index_pe(pe_for_level(k))
END DO

IF (mype == 0) THEN

  ! Allocate the buffer

  size_data_lbc = SUM(g_decomp_lbc_size) * full_lbc_levels
  ALLOCATE(data_lbc(size_data_lbc))
  displs(0)=0
  DO iproc=1, nproc-1
    displs(iproc)=displs(iproc-1) +            &
      g_decomp_lbc_size(iproc-1) * full_lbc_levels
  END DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP PRIVATE(iproc,iside, full_lbc_row_len, full_lbc_nrows,             & 
!$OMP decomp_lbc_nrows,first_lbc_pt, first_lbc_row,full_lbc_start_pt,    &
!$OMP decomp_lbc_start_pt,decomp_lbc_row_len, data_lbc_start, index_s,   &
!$OMP index_d)                                                           &
!$OMP SHARED(nproc,glsize,halosize,rimwidtha,                            &
!$OMP global_lbc_starta,displs,g_decomp_lbc_size,full_lbc,               &
!$OMP data_lbc,decomp_lbc_levels,g_lbc_starta,rim_type,halo_type,        &
!$OMP fld_type)
  DO iproc=0, nproc-1
    DO iside=1,4
      IF (g_at_extremity(iside,iproc)) THEN

        ! This processor is at edge type iside and so needs LBC data
        ! In order to copy data from full_lbc array to data_lbc ( that
        ! contains data to be scattered ) the following variables need to
        ! be calculated

        ! full_lbc_row_len : East-West dimension of the full lbc side
        ! full_lbc_nrows : North-South dimension of the full lbc side
        ! decomp_lbc_row_len : East-West dimension of decomp lbc side
        ! decomp_lbc_nrows : North-South dimension of decomp lbc side
        ! full_lbc_start_pt : First point of the decomposed lbc side inside
        !                     the 1d full_lbc array
        ! decomp_lbc_start_pt : First point of the decomposed lbc side inside
        !                       the 1d decomposed lbc array
        ! data_lbc_start : First point inside the data lbc array ( corresponds
        !                  to the first point of the distributed lbc_comp )


        CALL lbc_calc_size(iside,                                              &
                           full_lbc_row_len,                                   &
                           full_lbc_nrows,                                     &
                           decomp_lbc_row_len,                                 &
                           decomp_lbc_nrows,                                   &
                           full_lbc_start_pt,                                  &
                           decomp_lbc_start_pt,                                &
                           fld_type,                                           &
                           halo_type,                                          &
                           rim_type,                                           &
                           iproc)

        data_lbc_start = displs(iproc)

        ! Now we can do the copy

        DO k=1, decomp_lbc_levels
          DO j=1, decomp_lbc_nrows
            DO i=1, decomp_lbc_row_len
              index_s=full_lbc_start_pt + (i-1) + (j-1) * full_lbc_row_len
              index_d=data_lbc_start +                              &
                decomp_lbc_start_pt +                       &
                (i-1) + (j-1) * decomp_lbc_row_len +                &
                (k-1) * g_decomp_lbc_size (iproc)
              data_lbc(index_d)=full_lbc(index_s,k)
            END DO
          END DO
        END DO

      END IF ! At extremity
    END DO ! iside loop
  END DO ! iproc loop
!$OMP END PARALLEL DO

ELSE
  ALLOCATE(data_lbc(1))
END IF ! Master PE Section

! And Now we can do the scatter

CALL gc_get_communicator(my_comm, info)

sendcounts(:)=g_decomp_lbc_size(:) * decomp_lbc_levels
recvcount=decomp_lbc_size*decomp_lbc_levels

CALL mpl_scatterv(data_lbc, sendcounts, displs, mpl_real, &
                decomp_lbc, recvcount, mpl_real, 0, my_comm, info)

DEALLOCATE(data_lbc)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE scatter_atmos_lbcs

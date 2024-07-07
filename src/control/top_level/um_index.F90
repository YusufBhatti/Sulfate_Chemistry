! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Subroutine: UM_INDEX-----------------------------------------------
!
!    Purpose: Calculate addresses of component arrays within a
!             series of super arrays, made up of combinations of arrays
!             that require dynamic allocation. Lengths of the super
!             arrays are calculated and passed into U_MODEL for
!             dynamic allocation, reducing the no. of arguments needed
!             to be passed between top-level routines.
!
!    Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!    Logical components covered: C0
!
!    Project task: C0
!
!    External documentation: On-line UM document C1 - The top-level
!                            dynamic allocation
!
!    -------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

SUBROUTINE um_index(icode,cmessage)
!
! ----------------------------------------------------------------------

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr
USE UM_ParVars
USE Decomp_DB
USE decomp_params, ONLY: decomp_standard_atmos
USE submodel_mod, ONLY: n_internal_model
USE d1_array_mod, ONLY: d1_list_len
USE stash_array_mod, ONLY:                                                     &
    len_stlist, num_pseudo_lists, totitems, num_stash_pseudo, nitems,          &
    nstash_series_records, num_level_lists, nstash_series_block,               &
    num_stash_levels, nsttims, nsttabl, nsects, time_series_rec_len
USE nlsizes_namelist_mod, ONLY:                                        &
    aocpl_p_rows, aocpl_row_length, len_tot, n_obj_d1_max
USE river_routing_sizes_mod, ONLY: allocate_river_routing_sizes
USE river_inputs_mod, ONLY: l_rivers
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
!  Local variables
!
INTEGER :: icode             ! Work - Internal return code
CHARACTER(LEN=errormessagelength) :: cmessage    ! Work - Internal error message
INTEGER :: len_d1_addr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_INDEX'
! ----------------------------------------------------------------------
!  0. Start Timer running
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

icode=0

! ----------------------------------------------------------------------
!  1. Calculate     addresses in super array and each super array length
!
!  1.1   D1       super array
!
!           super array addresses
len_d1_addr=d1_list_len*n_obj_d1_max
!
! ----------------------------------------------------------------------
!  5.  Get global sizes (all PEs) for river routing regridding routines
IF (l_rivers) THEN
  aocpl_row_length                                                  &
              =decompDB(decomp_standard_atmos)%glsize(1,fld_type_p)
  aocpl_p_rows=decompDB(decomp_standard_atmos)%glsize(2,fld_type_p)
ELSE
  aocpl_row_length=1
  aocpl_p_rows=1
END IF

CALL allocate_river_routing_sizes() 

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE um_index

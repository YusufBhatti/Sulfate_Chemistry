! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_eg_poles_mod

IMPLICIT NONE

!  Subroutine Rcf_eg_poles

! Description:
!   Derive polar wind
!
! Method:
!   This routine assumes a regular grid. If support for
!   variable resolution grids is added then the checks in 
!   rcf_field_dependent_calcs and rcf_assign_vars_recon_horizontal
!   should be updated so that this routine can be called when needed.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_EG_POLES_MOD'

CONTAINS

SUBROUTINE rcf_eg_poles( fields_out, field_count_out,  &
                                hdr_out)

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE um_stashcode_mod, ONLY: &
    stashcode_u,stashcode_v, stashcode_prog_sec

USE rcf_field_type_mod, ONLY: &
    field_type

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE umPrintMgr

USE um_parvars, ONLY: &
    gc_all_proc_group

USE um_parcore, ONLY: &
    mype,nproc

USE rcf_gather_field_mod, ONLY: &
    rcf_gather_field_real
USE rcf_scatter_field_mod, ONLY: &
    rcf_scatter_field_real

USE decomp_params, ONLY: &
    decomp_rcf_output,                           &
    decomp_rcf_input

USE rcf_set_interp_flags_mod, ONLY:             &
  interp_h_only,                               &
  interp_v_only,                               &
  interp_copy,                                 &
  interp_all,                                  &
  interp_no_op

USE rcf_grid_type_mod, ONLY:                    &
    output_grid

USE rcf_headaddress_mod, ONLY: &
    fh_gridstagger_endgame

USE rcf_generate_heights_mod, ONLY: rcf_generate_heights
USE set_metric_terms_4a_mod, ONLY: eg_phi, tny

USE rcf_headaddress_mod, ONLY:                    &
  rc_polelong,                  rc_longspacing,  &
  rc_latspacing,                rc_polelat,      &
  rc_firstlong,                 rc_firstlat,     &
  rc_modeltop,                  rc_swldeg,       &
  rc_wedgedeg

USE rcf_calc_coords_mod, ONLY: &
    rcf_calc_coords

USE ereport_mod, ONLY: &
    ereport

USE eg_v_at_poles_mod, ONLY: eg_v_at_poles,array_dims

USE rcf_write_field_mod, ONLY: &
    rcf_write_field

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( um_header_type), INTENT(IN) :: hdr_out
INTEGER, INTENT(IN)               :: field_count_out

! Internal variables
TYPE( field_type ), POINTER       :: u
TYPE( field_type ), POINTER       :: v

CHARACTER (LEN=errormessagelength):: cmessage       ! used for EReport
CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_EG_POLES'
INTEGER                           :: errorstatus

INTEGER                           :: pos      ! position in array
INTEGER                           :: j        ! loop index

REAL                               :: xi1_p_global                     &
                                      (output_grid % glob_p_row_length,&
                                       output_grid % glob_p_rows)
REAL                               :: xi2_p_global                     &
                                      (output_grid % glob_p_row_length,&
                                       output_grid % glob_p_rows)
REAL                               :: xi1_u_global                     &
                                      (output_grid % glob_u_row_length,&
                                       output_grid % glob_u_rows )

REAL, ALLOCATABLE :: u_level(:,:)
REAL, ALLOCATABLE :: v_level(:,:)

INTEGER :: div,jend,pe,blk

TYPE (array_dims)  udims,vdims

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE (umMessage,FMT='(a)') 'Applying v-at-poles '
  CALL umPrint(umMessage,src='rcf_eg_poles')
END IF

CALL rcf_locate(stashcode_prog_sec, stashcode_u, &
                  fields_out, field_count_out, pos)
u => fields_out(pos)
CALL rcf_alloc_field( u )
CALL rcf_read_field(  u, hdr_out, decomp_rcf_output )

CALL rcf_locate(stashcode_prog_sec, stashcode_v, &
                  fields_out, field_count_out, pos)
v => fields_out(pos)
CALL rcf_alloc_field( v )
CALL rcf_read_field( v, hdr_out, decomp_rcf_output )


SELECT CASE(output_grid % grid_stagger)
CASE (fh_gridstagger_endgame)

CASE DEFAULT
  ! Print error if grid staggering is not catered for in this routine.
  cmessage = 'Grid staggering method is not supported.'
  errorstatus= 12
  CALL ereport( routinename, errorstatus, cmessage )
END SELECT


ALLOCATE( u_level( output_grid % glob_u_row_length,&
                   output_grid % glob_u_rows ) )
ALLOCATE( v_level( output_grid % glob_v_row_length,&
                   output_grid % glob_v_rows ) )


CALL rcf_calc_coords(hdr_out, output_grid, xi1_p_global, xi2_p_global, &
                     xi1_u_global)

! So do Gather 1 level per PE and then Scatter again.
div = output_grid % model_levels / nproc

IF ( nproc * div  <   output_grid % model_levels ) div = div + 1

pe = 0


divdo : DO blk = 1, div


  IF ( blk == div ) THEN
    jend = output_grid % model_levels
  ELSE
    jend =  blk * nproc
  END IF

  DO j = ((blk-1) * nproc) + 1, jend
    ! Will gather level j on PE pe
    CALL rcf_gather_field_real( u % DATA( :, j), u_level,  &
                                u % row_len,               &
                                u % rows,                  &
                                u % glob_row_len,          &
                                u % glob_rows, pe,         &
                                gc_all_proc_group )

    CALL rcf_gather_field_real( v % DATA( :, j), v_level,  &
                                v % row_len,               &
                                v % rows,                  &
                                v % glob_row_len,          &
                                v % glob_rows, pe,         &
                                gc_all_proc_group )

    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO


  udims%i_start= 0
  udims%i_end  = u % glob_row_len -1
  udims%j_start= 1
  udims%j_end  = u % glob_rows
  udims%k_start= 1
  udims%k_end  = 1
  udims%halo_i = 0
  udims%halo_j = 0

  vdims%i_start= 1
  vdims%i_end  = v % glob_row_len
  vdims%j_start= 0
  vdims%j_end  = v % glob_rows - 1
  vdims%k_start= 1
  vdims%k_end  = 1
  vdims%halo_i = 0
  vdims%halo_j = 0


  ! north pole (note: ENDGame stores the coordinates as 1d vectors,
  ! so some aspect of the rotation will not work anyway as such ....
  !                   hence here we simply refere to the first colum/row))
  CALL eg_v_at_poles(u_level(:,:),v_level(:,:),xi1_u_global(:, 1),    &
                         xi1_p_global(:, 1),xi2_p_global(1,:),-1.0,   &
                         udims%j_end, vdims%j_end,udims,vdims)
  ! south pole
  CALL eg_v_at_poles(u_level(:,:),v_level(:,:),xi1_u_global(:,1),     &
                         xi1_p_global(:,1),xi2_p_global(1,:), 1.0,    &
                         udims%j_start, vdims%j_start, udims,vdims)

  pe = 0

  DO j = ((blk-1) * nproc) + 1, jend

    CALL rcf_scatter_field_real( u%DATA( :, j), u_level,            &
                              u % row_len,                          &
                              u % rows,                             &
                              u % glob_row_len,                     &
                              u % glob_rows, pe,                    &
                              gc_all_proc_group )

    CALL rcf_scatter_field_real( v%DATA( :, j), v_level,            &
                              v % row_len,                          &
                              v % rows,                             &
                              v % glob_row_len,                     &
                              v % glob_rows, pe,                    &
                              gc_all_proc_group )

    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO
END DO divdo

DEALLOCATE ( v_level )
DEALLOCATE ( u_level )

CALL rcf_write_field( u, hdr_out, decomp_rcf_output )
CALL rcf_write_field( v, hdr_out, decomp_rcf_output )
CALL rcf_dealloc_field( u )
CALL rcf_dealloc_field( v )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE rcf_eg_poles
END MODULE rcf_eg_poles_mod

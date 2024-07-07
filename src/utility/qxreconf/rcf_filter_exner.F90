! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_filter_exner_mod

IMPLICIT NONE
!  Subroutine Rcf_Filter_Exner

! Description:
!   Filter the polar rows of Exner and Exner_Surf

! Method:
!   A basic 1-2-1 filter, applied n-times at the pole
!   and n-m times at the latitude m latitudes away
!   from the pole. n is currently hardcoded to 8.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_FILTER_EXNER_MOD'

CONTAINS

SUBROUTINE rcf_filter_exner( fields_out, field_count_out,  &
                              hdr_out )

USE rcf_write_field_mod, ONLY: &
    rcf_write_field

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE um_stashcode_mod, ONLY: &
    stashcode_exner_surf,&
    stashcode_exner, stashcode_prog_sec

USE rcf_field_type_mod, ONLY: &
    field_type

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE umPrintMgr

USE um_parvars, ONLY: &
    gc_all_proc_group

USE um_parcore, ONLY: &
    mype,  nproc

USE decomp_params, ONLY: &
    decomp_rcf_output,                           &
    decomp_rcf_input

USE rcf_grid_type_mod, ONLY:                    &
    output_grid

USE rcf_headaddress_mod, ONLY: &
    fh_gridstagger_endgame

USE ereport_mod, ONLY: &
    ereport

USE rcf_gather_field_mod, ONLY: &
    rcf_gather_field_real
USE rcf_scatter_field_mod, ONLY: &
    rcf_scatter_field_real

USE decomp_params, ONLY: &
    decomp_rcf_output

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
TYPE( field_type ), POINTER       :: exner
TYPE( field_type ), POINTER       :: exner_surf

CHARACTER (LEN=errormessagelength):: cmessage       ! used for EReport
CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_FILTER_EXNER'
INTEGER                           :: errorstatus

INTEGER                           :: pos      ! position in array
INTEGER                           :: i,j,k,ii,jj ! loop index
INTEGER                           :: div,blk,pe,jend



REAL, ALLOCATABLE :: p_level(:,:,:)

INTEGER, PARAMETER :: n_filt_nd = 8
INTEGER            :: last,next,n_filt

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE (umMessage,FMT='(a,i4)') 'Filtering exner ',n_filt_nd
  CALL umPrint(umMessage,src='rcf_filter_exner')
END IF

SELECT CASE(output_grid % grid_stagger)
CASE (fh_gridstagger_endgame)

CASE DEFAULT
  ! Print error if grid staggering is not catered for in this routine.
  cmessage = 'Grid staggering method is not supported.'
  errorstatus= 12
  CALL ereport( routinename, errorstatus, cmessage )
END SELECT


CALL rcf_locate( stashcode_prog_sec, stashcode_exner,       &
                 fields_out,field_count_out,pos)
exner => fields_out(pos)
CALL rcf_alloc_field( exner )
CALL rcf_read_field ( exner, hdr_out, decomp_rcf_output )


CALL rcf_locate( stashcode_prog_sec, stashcode_exner_surf,       &
                 fields_out,field_count_out,pos)
exner_surf => fields_out(pos)
CALL rcf_alloc_field( exner_surf )
CALL rcf_read_field(  exner_surf, hdr_out, decomp_rcf_output )

ALLOCATE( p_level( output_grid % glob_p_row_length, &
                   output_grid % glob_p_rows,2 ) )

! So do Gather 1 level per PE and then Scatter again.
div = (output_grid % model_levels +1)/ nproc

IF ( nproc * div  <   (output_grid % model_levels+1) ) div = div + 1

pe = 0

divdo : DO blk = 1, div

  IF (blk == div ) THEN
    jend = (output_grid % model_levels+1)
  ELSE
    jend =  blk * nproc
  END IF

  DO j = ((blk-1) * nproc) + 1, jend
    ! Will gather level j on PE pe
    CALL rcf_gather_field_real( exner%DATA( :, j), p_level(:,:,1),    &
                                output_grid % loc_p_row_length,       &
                                output_grid % loc_p_rows,             &
                                output_grid % glob_p_row_length,      &
                                output_grid % glob_p_rows, pe,        &
                                gc_all_proc_group )
    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO

  last = 1
  next = 2

  DO n_filt = 1, n_filt_nd

    p_level(:,:,next) =  p_level(:,:,last)

    DO j=1,n_filt

      jj=output_grid % glob_p_rows-j+1
      DO i=2,output_grid % glob_p_row_length-1
        p_level(i,j ,next) = 0.25*(   p_level(i-1,j ,last)+&
                                   2.0*p_level(i,  j ,last)+&
                                      p_level(i+1,j ,last))

        p_level(i,jj,next) = 0.25*(   p_level(i-1,jj,last)+&
                                   2.0*p_level(i,  jj,last)+&
                                      p_level(i+1,jj,last))
      END DO
      i  = output_grid % glob_p_row_length
      ii = 1
      p_level(i,j ,next) = 0.25*(   p_level(i-1,j ,last)+&
                                 2.0*p_level(i,  j ,last)+&
                                    p_level(ii, j,last))

      p_level(i,jj,next) = 0.25*(   p_level(i-1,jj,last)+&
                                 2.0*p_level(i,  jj,last)+&
                                    p_level(ii, jj,last))
      ii = output_grid % glob_p_row_length
      i  = 1
      p_level(i,j ,next) = 0.25*(   p_level(i+1,j ,last)+&
                                 2.0*p_level(i,  j ,last)+&
                                    p_level(ii, j,last))

      p_level(i,jj,next) = 0.25*(   p_level(i+1,jj,last)+&
                                 2.0*p_level(i,  jj,last)+&
                                    p_level(ii, jj,last))
    END DO

    IF (last==1) THEN
      last =2
    ELSE
      last =1
    END IF

    IF (next==1) THEN
      next =2
    ELSE
      next =1
    END IF

  END DO

  pe = 0

  DO j = ((blk-1) * nproc) + 1, jend

    CALL rcf_scatter_field_real(exner%DATA( :, j), p_level(:,:,last),&
                                output_grid % loc_p_row_length,      &
                                output_grid % loc_p_rows,            &
                                output_grid % glob_p_row_length,     &
                                output_grid % glob_p_rows, pe,       &
                                gc_all_proc_group )
    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO
END DO divdo


pe = 0

CALL rcf_gather_field_real( exner_surf%DATA( :, 1), p_level(:,:,1),&
                            output_grid % loc_p_row_length,        &
                            output_grid % loc_p_rows,              &
                            output_grid % glob_p_row_length,       &
                            output_grid % glob_p_rows, pe,         &
                            gc_all_proc_group )

IF ( mype == 0 ) THEN

  last = 1
  next = 2

  DO n_filt = 1, n_filt_nd

    p_level(:,:,next) =  p_level(:,:,last)

    DO j=1,n_filt
      jj=output_grid % glob_p_rows-j+1
      DO i=2,output_grid % glob_p_row_length-1
        p_level(i,j ,next) = 0.25*(   p_level(i-1,j ,last)+&
                                   2.0*p_level(i,  j ,last)+&
                                      p_level(i+1,j ,last))

        p_level(i,jj,next) = 0.25*(   p_level(i-1,jj,last)+&
                                   2.0*p_level(i,  jj,last)+&
                                      p_level(i+1,jj,last))
      END DO
      i  = output_grid % glob_p_row_length
      ii = 1
      p_level(i,j ,next) = 0.25*(   p_level(i-1,j ,last)+&
                                 2.0*p_level(i,  j ,last)+&
                                    p_level(ii, j,last))

      p_level(i,jj,next) = 0.25*(   p_level(i-1,jj,last)+&
                                 2.0*p_level(i,  jj,last)+&
                                    p_level(ii, jj,last))
      ii = output_grid % glob_p_row_length
      i  = 1
      p_level(i,j ,next) = 0.25*(   p_level(i+1,j ,last)+&
                                 2.0*p_level(i,  j ,last)+&
                                    p_level(ii, j,last))

      p_level(i,jj,next) = 0.25*(   p_level(i+1,jj,last)+&
                                 2.0*p_level(i,  jj,last)+&
                                    p_level(ii, jj,last))
    END DO

    IF (last==1) THEN
      last =2
    ELSE
      last =1
    END IF
    IF (next==1) THEN
      next =2
    ELSE
      next =1
    END IF

  END DO

END IF

CALL rcf_scatter_field_real(exner_surf%DATA(:,1), p_level(:,:,last), &
                                output_grid % loc_p_row_length,      &
                                output_grid % loc_p_rows,            &
                                output_grid % glob_p_row_length,     &
                                output_grid % glob_p_rows, pe,       &
                                gc_all_proc_group )

CALL rcf_write_field( exner_surf, hdr_out, decomp_rcf_output )
CALL rcf_write_field( exner,      hdr_out, decomp_rcf_output )

DEALLOCATE ( p_level )

CALL rcf_dealloc_field( exner_surf )
CALL rcf_dealloc_field( exner )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_filter_exner
END MODULE rcf_filter_exner_mod

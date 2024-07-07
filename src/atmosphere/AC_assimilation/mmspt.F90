! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE MMSPT   MMSPTW   ---------------------------------------
!
!     2 Subroutines in deck : MMSPT and MMSPTW
!     MMSPTW is same as MMSPT for Limited Area Winds.
!
!    Purpose : Provide weighted Mean and S.D. plus extremes
!              on Model grids. Used for Weights and Increments
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE mmspt_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MMSPT_MOD'

CONTAINS

SUBROUTINE mmspt (pvals,klev,kgrid,pntlab,                        &
                  row_length,p_rows,                              &
                  phi_p)
!
!     CALCULATE AREA-WEIGHTED MEAN & MEAN-SQUARE
!     OF A FIELD ON THE MODEL GRID
!
!     KLEV : LEVEL OF FIELD (USED IN PRINT OUT ONLY)
!
!     IF KGRID=0 DATA ON P*/THETA MODEL GRID
!             =1 DATA ON WIND MODEL GRID
!
!     16 CHARACTER PNTLAB IS PRINTED TO IDENTIFY MEAN & MEAN SQ.
!Lx
USE conversions_mod, ONLY: pi_over_180
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE field_types, ONLY: fld_type_p
USE atmos_max_sizes, ONLY: Max2DFieldSize
USE UM_ParParams, ONLY: halo_type_no_halo
USE comobs_mod, ONLY: nobtypmx
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE ac_control_mod
USE model_domain_mod, ONLY: l_regular

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------
INTEGER :: klev,kgrid                !IN LEVEL AND GRID IDENTIFER
INTEGER :: row_length,p_rows         !IN MODEL DIMENSIONS
REAL :: pvals(row_length,*)       !IN MODEL FIELD
CHARACTER(LEN=16) :: pntlab               !IN CHARACTER IDENTIFER OF FIELD

REAL :: phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)

!     LOCAL VARIABLES
REAL :: pvals_g(Max2DFieldSize)
INTEGER :: row_length_global
INTEGER :: p_rows_global
REAL :: wt,sumwt,SUM,zrow1,zlat,zdlat,zm,zms,zmax,zmin,zrms
INTEGER :: jptf,jptl,jrowf,jrowl,jrow,jpt,npts
INTEGER :: imaxpt,imaxro,iminpt,iminro

! Variables required for variable resolution grids
! Since the total number of rows is unknown here it is
! made allocatable (although glsize could be used)
REAL :: phi_p_local(row_length,p_rows)
REAL, ALLOCATABLE :: phi_p_global(:,:)
INTEGER :: i,j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MMSPT'

!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ltimer_ac) CALL timer('MMSPT   ',3)
row_length_global=glsize(1,fld_type_p)
p_rows_global=glsize(2,fld_type_p)

! KGRID must be 0
IF (kgrid /= 0) THEN
  WRITE(umMessage,*) 'INVALID KGRID IN MMSPT ',kgrid
  CALL umPrint(umMessage,src='mmspt')
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

!     JPTF  = FIRST POINT IN ROW
!     JPTL  = LAST  POINT IN ROW
!     JROWF = FIRST ROW
!     JROWL = LAST  ROW

!     OUTSIDE TWO BOUNDARY POINTS OF ELF GRID NOT USED
jptf  = 3
jptl  = row_length_global-2
jrowf = 3
jrowl = p_rows_global-2

! Gather P_VALS onto a global field PVALS_G
! DEPENDS ON: gather_field
CALL Gather_Field( pvals, pvals_g, row_length, p_rows,            &
                   row_length_global, p_rows_global,              &
                   fld_type_p, halo_type_no_halo,                 &
                   0, gc_all_proc_group )

! If variable grid collect together phi coordinates of entire grid
IF (.NOT. l_regular) THEN

  ALLOCATE (phi_p_global(row_length_global,p_rows_global))

  ! Copy information from PHI_P to PHI_P_LOCAL
  !  excluding halos
  DO i=1,p_rows
    DO j=1,row_length
      phi_p_local(j,i)=phi_p(j,i)
    END DO
  END DO

  ! Gather up all PHI_P values into global array
  ! DEPENDS ON: gather_field
  CALL gather_field(phi_p_local,   phi_p_global,                  &
                    row_length,    p_rows,                        &
                    row_length_global,  p_rows_global,            &
                    fld_type_p,    halo_type_no_halo,             &
                    0,             gc_all_proc_group )


END IF


IF (mype == 0) THEN

  zrow1 = xlatn - dlat*(jrowl-jrowf+1)

  zlat   = zrow1*pi_over_180
  zdlat  = dlat *pi_over_180

  !  SET ACCUMULATORS
  zmax=pvals_g(jptf+(jrowf-1)*row_length_global)
  zmin=pvals_g(jptf+(jrowf-1)*row_length_global)

  zm    = 0.0
  zms   = 0.0
  zrms  = 0.0
  sumwt = 0.0
  wt    = 0.0
  imaxpt= jptf
  imaxro= jrowf
  iminpt= jptf
  iminro= jrowf

  DO jrow = jrowf,jrowl
    IF (l_regular) THEN
      wt = COS(zlat)
    ELSE
      wt = COS(phi_p_GLOBAL(1,jrow+1)) ! note offset by 1
    END IF

    sumwt = sumwt+wt

    !     CALCULATE MEAN

    SUM = 0.0
    DO jpt=jptf,jptl
      SUM = SUM + pvals_g(jpt+(jrow-1)*row_length_global)
    END DO
    zm = zm + SUM*wt

    !     CALCULATE MEAN SQUARE

    SUM = 0.0
    DO jpt=jptf,jptl
      SUM = SUM + pvals_g(jpt+(jrow-1)*row_length_global)*            &
                  pvals_g(jpt+(jrow-1)*row_length_global)
    END DO
    zms = zms + SUM*wt

    !     CALCULATE MAX

    DO jpt=jptf,jptl
      IF (pvals_g(jpt+(jrow-1)*row_length_global) >  zmax) THEN
        zmax=pvals_g(jpt+(jrow-1)*row_length_global)
        imaxpt=jpt
        imaxro=jrow
      END IF
    END DO

    !     CALCULATE MIN

    DO jpt=jptf,jptl
      IF (pvals_g(jpt+(jrow-1)*row_length_global) <  zmin) THEN
        zmin=pvals_g(jpt+(jrow-1)*row_length_global)
        iminpt=jpt
        iminro=jrow
      END IF
    END DO

    zlat = zlat + zdlat
  END DO ! jrow

  !     EVALUATE STATS AND WRITE OUT

  npts = jptl-jptf+1
  IF (sumwt /= 0.0) wt   = 1.0/(sumwt*npts)
  zm   = zm  * wt
  zms  = zms * wt
  IF (zms >  0.0)   zrms = SQRT(zms)

  WRITE(umMessage,62)                                                       &
    pntlab,klev,zm,zrms,zmax,imaxro,imaxpt,zmin,iminro,iminpt
  CALL umPrint(umMessage,src='mmspt')

  62    FORMAT(1x,a16,2x,i4,' MEAN=',g12.5,' RMS=',g12.5,                 &
         ' MAX=',g12.5,' AT (',i4,',',i4,')',                             &
         ' MIN=',g12.5,' AT (',i4,',',i4,')')

END IF      ! mype == 0

IF (ltimer_ac) CALL timer('MMSPT   ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE mmspt
END MODULE mmspt_mod

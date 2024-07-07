! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine CALC_FIT_FSAT-------------------------------------------
!
!   Purpose: To speed up the large scale hydrology code (LTOP=TRUE)
!            dramatically. This is done by fitting exponential
!            functions to the incomplete gamma function for each grid
!            box and the complete range of possible "water table"
!            (top_crit) cases - see documentation.
!            Estimates the fitted parameters for Fsat=function(ZW)
!            and  Fwet=function(ZW) for each land grid point.
!            (Calculating the incomplete gamma function for each grid
!            box at each time step was very time consuming).
!                                                             !
! Documentation: UNIFIED MODEL DOCUMENTATION PAPER NO 25
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration

SUBROUTINE calc_fit_fsat(soil_pts,soil_index,npnts                &
  ,fexp,ti_mean,ti_sig,gamtot,zdepth                              &
  ,a_fsat,c_fsat,a_fwet,c_fwet)


USE calc_fsat_mod,         ONLY: calc_fsat

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Diag

USE Ereport_Mod, ONLY:                                            &
  Ereport

! Replaces c_topog.h
USE jules_hydrology_mod, ONLY: zw_max, nfita

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(IN) :
INTEGER ::                                                        &
 npnts                                                            &
                  ! IN No. of land points.
,soil_pts         ! IN No. of land soil points.

REAL ::                                                           &
 zdepth           ! IN Standard Soil model DEPTH.

!   Array arguments with intent(IN) :
INTEGER ::                                                        &
 soil_index(npnts)! IN Array of soil points.

REAL ::                                                           &
 ti_mean(npnts)                                                   &
                  ! IN Gridbox mean topographic index.
,ti_sig(npnts)                                                    &
                  ! IN Std. deviation in topographic index.
,fexp(npnts)                                                      &
                  ! IN Exp. decay in deep layer.
,gamtot(npnts)    ! IN Total gamma function.

!   Array arguments with intent(OUT) :
REAL ::                                                           &
   a_fsat(npnts)                                                  &
                  ! OUT Fitting parameter for Fsat.
  ,c_fsat(npnts)                                                  &
                  ! OUT Fitting parameter for Fsat.
  ,a_fwet(npnts)                                                  &
                  ! OUT Fitting parameter for Fwet.
  ,c_fwet(npnts)  ! OUT Fitting parameter for Fwet.

! Local scalars:
INTEGER :: nzw  ! Number of ZW values used in fit.
PARAMETER(nzw=200)  ! Maximum value for a significant improvement
!                         ! in the fit.

INTEGER ::                                                        &
 i,j,iz                                                           &
             ! Loop counters.
,ifita
             ! Loop counters.

REAL :: dzw     ! WORK ZW increment ; defined by ZW_MAX and NZW.

REAL ::                                                           &
 rms                                                              &
             ! WORK RMS errors for given fsat fit values.
,rmsw                                                             &
             ! WORK RMS errors for given fwet fit values.
,rmsold                                                           &
             ! WORK RMS errors for given fsat fit values.
!                  !      for best fit so far.
      ,rmswold                                                          &
                   ! WORK RMS errors for given fwet fit values
!                  !      for best fit so far.
      ,cfitmin                                                          &
                   ! WORK Minimum possible value for Cfit.
      ,cfitmax                                                          &
                   ! WORK Maximum possible value for Cfit.
      ,cfit                                                             &
                   ! WORK CFit value for given loop.
      ,thr_err     ! WORK Error threshold value


PARAMETER(cfitmin=0.0)
PARAMETER(thr_err=5.0e-3)

! Local arrays:
REAL ::                                                           &
 fsat_calc(nzw)                                                   &
                        ! WORK Surface saturation fraction.
,fsat_fit(nzw)                                                    &
                        ! WORK Fitted surface saturation fraction.
,fwet_calc(nzw)                                                   &
                        ! WORK Wetland fraction.
,fwet_fit(nzw)                                                    &
                        ! WORK Fitted wetland fraction.
,dumzw(nzw)                                                       &
                        ! WORK Dummy water table depth (m).
,dumfsat(1)                                                       &
                        ! WORK Dummy surface saturation fraction.
,dumfwetl(1)            ! WORK Dummy wetland fraction.

REAL ::                                                           &
 top_crit(nzw)                                                    &
                        ! WORK LOG(QBASE_MAX/QBASE) -see document.
,top_crit1z(1)                                                    &
                        ! WORK Aas above but for an individual zw.
,top_min                                                          &
                        ! WORK value for when zw=zw_max.
,wutot(1)                                                         &
                        ! WORK Dummy (set to 1.0).
,gamtot1(1)             ! Temp variable to allow interface with calc_fsat

INTEGER :: soil_index1(1) = 1 !For the call to calc_fsat to pass a single 1 in
                              !array format

LOGICAL :: L_Gamtot=.FALSE.

REAL :: temp1
INTEGER :: errorstatus
INTEGER            ::     error_count
INTEGER, PARAMETER :: max_error_count = 50
CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName ='CALC_FIT_FSAT'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

error_count = 0

cfitmax=0.15*nfita

! Define the water table depths to be used in the fitting process:
dzw=1.0/REAL(nzw)*zw_max
DO iz=1,nzw
  dumzw(iz)=REAL(iz-1)*dzw
END DO
wutot(1)    = 1.0
dumfsat(1)  = 0.0
dumfwetl(1) = 0.0

! Calculate TOP_CRIT for the water tables depths:
DO j=1,soil_pts
  i=soil_index(j)

  IF (ti_mean(i) >  0.0 .AND. ti_sig(i) >  0.0 ) THEN
    top_min=1.0/fexp(i)*EXP(-fexp(i)*(zw_max-zdepth))

    temp1 = zdepth+1.0/fexp(i)-top_min

    DO iz=1,nzw

      IF (dumzw(iz) <= zdepth) THEN
        top_crit1z(1) = -LOG(1-(dumzw(iz)/temp1))
      ELSE
        top_crit1z(1) =  LOG(temp1/                                           &
                         (1.0/fexp(i) * EXP(-fexp(i) * (dumzw(iz)-zdepth))    &
                         -top_min))
      END IF

      ! Calculate FSAT and FWET for one ZW at one soil pnt:

      !This call to calc_fsat is somewhat abused, so we need to use arguments
      !that are arrays of size 1
      gamtot1(1) = gamtot(i)  !INOUT, so remember to copy back afterwards

      CALL calc_fsat(l_gamtot,1,soil_index1,1,ti_mean(i),ti_sig(i),&
                     wutot,top_crit1z,gamtot1,dumfsat,dumfwetl)

      gamtot(i) = gamtot1(1)

      fsat_calc(iz)=dumfsat(1)
      fwet_calc(iz)=dumfwetl(1)
      top_crit(iz)=top_crit1z(1)
      IF (iz == 1) THEN  ! Values at zw=0m
        a_fsat(i)=fsat_calc(iz)
        a_fwet(i)=fwet_calc(iz)
      END IF
    END DO

    rmsold=1.0e10
    rmswold=1.0e10

    DO ifita=1,nfita
      cfit=cfitmax*(ifita)/REAL(nfita)

      ! This isnt really root mean squared - just a measure of fitness.
      rms=0.0
      rmsw=0.0
      !TOP_CRIT=TI_MAX when zw=zw_max
      DO iz=1,nzw
        fsat_fit(iz)=a_fsat(i)*EXP(-cfit*top_crit(iz))
        fwet_fit(iz)=a_fwet(i)*EXP(-cfit*top_crit(iz))
        rms=rms+(fsat_calc(iz)-fsat_fit(iz))**2
        rmsw=rmsw+(fwet_calc(iz)-fwet_fit(iz))**2
      END DO            !ZW
      rms=rms/REAL(nzw)
      rmsw=rmsw/REAL(nzw)

      IF (rms < rmsold) THEN
        rmsold=rms
        c_fsat(i)=cfit
      END IF
      IF (rmsw <  rmswold) THEN
        rmswold=rmsw
        c_fwet(i)=cfit
      END IF
    END DO

    IF (rmsold >= thr_err**2) THEN
      IF (c_fsat(i) <= cfitmin .OR. c_fsat(i) >= cfitmax) THEN

        DO iz=1,nzw
          fsat_fit(iz)=a_fsat(i)*EXP(-c_fsat(i)*top_crit(iz))
        END DO               !ZW

        WRITE(umMessage,'(a,i4,3f20.10)')                               &
          'ERROR CFIT FSAT',i,c_fsat(i),cfitmin,cfitmax
        CALL umPrint(umMessage,src='calc_fit_fsat')
        WRITE(umMessage,'(a,3f20.10)')'fsat_calc='                      &
          ,fsat_calc(1),fsat_calc(3),fsat_calc(5)
        CALL umPrint(umMessage,src='calc_fit_fsat')
        WRITE(umMessage,'(a,3f20.10)')'fsat_fit='                       &
          ,fsat_fit(1),fsat_fit(3),fsat_fit(5)
        CALL umPrint(umMessage,src='calc_fit_fsat')
        WRITE(umMessage,'(a,f20.10)')'RMS=',SQRT(rmsold)
        CALL umPrint(umMessage,src='calc_fit_fsat')
        ErrorStatus = 35
        WRITE(CMessage, '(A)') 'Error in CFIT FSAT in LSH model setup'

        CALL Ereport ( RoutineName, ErrorStatus, CMessage)
      END IF
    END IF

    IF (rmswold >= thr_err**2) THEN
      IF (c_fwet(i) <= cfitmin .OR. c_fwet(i) >= cfitmax) THEN

        DO iz=1,nzw
          fsat_fit(iz)=a_fsat(i)*EXP(-c_fsat(i)*top_crit(iz))
          fwet_fit(iz)=a_fwet(i)*EXP(-c_fwet(i)*top_crit(iz))
        END DO               !ZW

        WRITE(umMessage,'(a,i4,3f20.10)')                               &
          'ERROR CFIT FWET',i,c_fwet(i),cfitmin,cfitmax
        CALL umPrint(umMessage,src='calc_fit_fsat')
        WRITE(umMessage,'(a,3f20.10)')'fwet_calc='                      &
         ,fwet_calc(1),fwet_calc(3),fwet_calc(5)
        CALL umPrint(umMessage,src='calc_fit_fsat')
        WRITE(umMessage,'(a,3f20.10)')'fwet_fit='                       &
         ,fwet_fit(1),fwet_fit(3),fwet_fit(5)
        CALL umPrint(umMessage,src='calc_fit_fsat')
        WRITE(umMessage,'(a,f20.10)')'RMSW=',SQRT(rmswold)
        CALL umPrint(umMessage,src='calc_fit_fsat')
        WRITE(umMessage,'(a,3f20.10)')'(fsat_calc=)'                    &
         ,fsat_calc(1),fsat_calc(3),fsat_calc(5)
        CALL umPrint(umMessage,src='calc_fit_fsat')
        WRITE(umMessage,'(a,3f20.10)')'(fsat_fit=)'                     &
         ,fsat_fit(1),fsat_fit(3),fsat_fit(5)
        CALL umPrint(umMessage,src='calc_fit_fsat')
        WRITE(umMessage,'(a,f20.10)')'(RMS=)',SQRT(rmsold)
        CALL umPrint(umMessage,src='calc_fit_fsat')
        ErrorStatus = 40
        WRITE(CMessage, '(A)') 'Error in CFIT FWET in LSH model setup'

        CALL Ereport ( RoutineName, ErrorStatus, CMessage)
      END IF
    END IF
    IF (     (rmsold >= thr_err**2 .OR. rmswold >= thr_err**2)   &
       .AND. (PrintStatus >= PrStatus_Diag)                      &
       .AND. (error_count <= max_error_count) ) THEN
      error_count = error_count + 1
      IF (error_count < max_error_count) THEN
        WRITE(umMessage,'(a,2f20.10)')                          &
           'Warning LSH RMS Error in fit:'                      &
           ,SQRT(rmsold),SQRT(rmswold)
        CALL umPrint(umMessage,src='calc_fit_fsat')
      ELSE IF (error_count == max_error_count) THEN
        WRITE(umMessage,'(a)')                                  &
           'LSH RMS Error printing DISCONTINUED.'
        CALL umPrint(umMessage,src='calc_fit_fsat')
      END IF
    END IF
  END IF

END DO                     ! NPNTS
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_fit_fsat

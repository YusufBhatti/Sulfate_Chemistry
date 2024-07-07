! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Compute data lengths and addresses for primary fields
! Subroutine Interface:
SUBROUTINE primary(isec,iitm,Im_index,Im_ident,Sm_ident,          &
                  rlevs,raddress,PIrow,ErrorStatus,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE stextend_mod,ONLY: in_s, ppind_s, d1_paddr, n_obj_d1,  &
                       max_d1_len, extra_d1,               &
         d1_type, d1_im, d1_extra_info, d1_levs, d1_sect,  &
         prog
USE cstash_mod, ONLY: iflag, iplast, ipfirst, ispace, ipseudo,  &
                      ibot, itop, igp, halo_type, ilev
USE ukca_option_mod,   ONLY: tr_ukca_a
USE ukca_tracer_stash, ONLY: a_max_ukcavars
USE um_stashcode_mod, ONLY: stashcode_ukca_sec

USE nlsizes_namelist_mod, ONLY: &
    tr_ukca

USE errormessagelength_mod, ONLY: errormessagelength
USE stash_model_mod, ONLY:                                                     &
    nhead, nheadsub, len_prim, len_primim, global_len_prim, global_len_primim, &
    len_extra

IMPLICIT NONE
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code description:
!   FORTRAN 90
!   Written to UMDP3 programming standards version 8.3.
!

! Subroutine arguments:
!   Scalar arguments with intent(in):
INTEGER :: isec      ! Section number
INTEGER :: iitm      ! Item number
INTEGER :: Im_ident  ! Current internal model number
INTEGER :: Im_index  ! Current position in internal model list
INTEGER :: Sm_ident  ! Submodel identifier (absolute)
!   Scalar arguments with intent(out):
CHARACTER(LEN=errormessagelength) :: cmessage

! ErrorStatus:
INTEGER :: ErrorStatus

! Local scalars:
LOGICAL :: model_lev
LOGICAL :: laddr
LOGICAL :: lmask
INTEGER :: rlevs      ! No. of levels for reconfiguration
INTEGER :: dlevs      ! No of levels inc pseudo levels
INTEGER :: rplevs     ! & of pseudo-levels
INTEGER :: raddress   ! Address for reconfiguration
INTEGER :: i
INTEGER :: il1,il2
INTEGER :: ipl1,ipl2
INTEGER :: LEN        ! Data length for primary field
INTEGER :: global_LEN ! Global data length for primary field
INTEGER :: PIrow      ! Counter for ProgItems array

! Function and subroutine calls:
LOGICAL :: disct_lev

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRIMARY'

!- End of Header ---------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Find out whether the primary is included for this version
! DEPENDS ON: tstmsk
CALL tstmsk(Im_ident,isec,lmask,laddr,ErrorStatus,cmessage)

IF (isec == stashcode_ukca_sec .AND. iitm <= a_max_ukcavars) THEN
  tr_ukca_a(iitm) = .FALSE.
END IF

IF (laddr) THEN

  !       If this is a UKCA tracer, then use the return value
  !       from tstmsk and space code to turn it on or off
  IF (isec == stashcode_ukca_sec .AND. iitm <= a_max_ukcavars .AND.  &
           ispace /= 10) THEN
    tr_ukca_a(iitm) = .TRUE.
    tr_ukca = tr_ukca + 1
  END IF

  IF (ispace == 10) THEN
    ! Space code 10 means: no space is required for this item in D1 or
    !  the dump, but stashmaster data is required, so an "address" of
    !  -1 is set to ensure that the corresponding record will be read
    !  into PPXI in routine GET_PPX_PART (called by U_MODEL).
    in_s(1,Im_ident,isec,iitm)=-1
  ELSE

    IF (isec  ==  0) THEN
      ! Start address for model levels in PP array
      ppind_s(Im_ident,iitm) = nhead(Im_ident)+1
    END IF ! IF (ISEC  ==  0)

    ! Find address length per level
    ! DEPENDS ON: addrln
    CALL addrln(igp,halo_type,LEN,local_data)
    ! DEPENDS ON: addrln
    CALL addrln(igp,halo_type,global_LEN,                           &
                global_dump_data)

    ! DEPENDS ON: disct_lev
    model_lev=disct_lev(ilev,ErrorStatus,cmessage)
    IF (model_lev .OR. (ilev == 5 .AND. ipseudo /= 0)) THEN
      ! Field has model levels - decode level codes
      IF (ilev  /=  5) THEN
        ! DEPENDS ON: levcod
        CALL levcod(ibot,il1,ErrorStatus,cmessage)
        ! DEPENDS ON: levcod
        CALL levcod(itop,il2,ErrorStatus,cmessage)
      ELSE
        il1=1
        il2=1
      END IF
      ! No. of model levels for D1 addressing
      dlevs=il2-il1+1
      ! Initialise first & last pseudo level indices
      ipl1 =0
      ipl2 =0
      IF (iflag == 0 .AND. ipseudo /= 0) THEN
        ! Primary with input on all available pseudo levels -
        !   decode pseudo level codes
        ! DEPENDS ON: pslevcod
        CALL pslevcod(ipfirst,ipl1,'F',ErrorStatus,cmessage)
        ! DEPENDS ON: pslevcod
        CALL pslevcod(iplast ,ipl2,'L',ErrorStatus,cmessage)
        dlevs=dlevs*(ipl2-ipl1+1)
      END IF
      rplevs=ipl2-ipl1+1
      ! Multiply length per level by no. of levels
      LEN=len*(il2-il1+1)*(ipl2-ipl1+1)
      global_LEN=global_LEN*(il2-il1+1)*(ipl2-ipl1+1)
      IF (ispace /= 4 .AND. ispace /= 9) THEN
        ! Increment no. of headers
        nhead   (Im_ident)=   nhead(Im_ident)                       &
                                  +(il2-il1+1)*(ipl2-ipl1+1)
        NHeadSub(Sm_ident)=NHeadSub(Sm_ident)                       &
                                  +(il2-il1+1)*(ipl2-ipl1+1)
      END IF
    ELSE
      ! Not model levels
      rlevs=1
      dlevs=1
      rplevs=1
      IF (ispace /= 4 .AND. ispace /= 9) THEN
        nhead   (Im_ident)=nhead   (Im_ident)+1
        NHeadSub(Sm_ident)=NHeadSub(Sm_ident)+1
      END IF
    END IF

    ! The input start address for primary (m,0,i) is assigned
    !  to IN_S(1,m,0,i).
    ! Addresses are set relative to the beginning of the primary data,
    !  since the primary data starts at the beginning of D1.
    IF (ispace /= 5) THEN
      IF (ispace /= 9) THEN
        ! Start address for this primary field
        in_s(1,Im_ident,isec,iitm)=len_prim(Sm_ident)+1
        ! Increment len_prim by LEN (=data length for this primary field)
        len_prim  (Sm_ident)      =len_prim  (Sm_ident)+LEN
        len_primim(Im_ident)      =len_primim(Im_ident)+LEN
        ! Information for preliminary D1 addressing array
        n_obj_d1(Sm_ident)     =n_obj_d1(Sm_ident)+1
        IF (n_obj_d1(Sm_ident) <= max_d1_len) THEN
          d1_paddr(d1_type,n_obj_d1(Sm_ident),Sm_ident)=prog
          d1_paddr(d1_im,n_obj_d1(Sm_ident),Sm_ident)=Im_ident
          d1_paddr(d1_extra_info,n_obj_d1(Sm_ident),Sm_ident)=iitm
          d1_paddr(d1_levs,n_obj_d1(Sm_ident),Sm_ident)=dlevs
          d1_paddr(d1_sect,n_obj_d1(Sm_ident),Sm_ident)=isec
        END IF
        global_len_prim  (Sm_ident) =global_len_prim  (Sm_ident)+global_LEN
        global_len_primim(Im_ident) =global_len_primim(Im_ident)+global_LEN
        ! Dual addresses for ocean fields with dual time level
        IF (ispace == 8) THEN
          ! Information for preliminary D1 addressing array
          n_obj_d1(Sm_ident)     =n_obj_d1(Sm_ident)+1
          IF (n_obj_d1(Sm_ident) <= max_d1_len) THEN
            d1_paddr(d1_type,n_obj_d1(Sm_ident),Sm_ident)=extra_d1
            d1_paddr(d1_im,n_obj_d1(Sm_ident),Sm_ident)=Im_ident
            d1_paddr(d1_extra_info,n_obj_d1(Sm_ident),Sm_ident)=iitm
            d1_paddr(d1_levs,n_obj_d1(Sm_ident),Sm_ident)=dlevs
            d1_paddr(d1_sect,n_obj_d1(Sm_ident),Sm_ident)=isec
          END IF
        END IF
      ELSE ! Space = 9
        !           These are EXNER etc items. Record the address relative
        !           to start of len_extra space in D1. A loop in ADDRES
        !           will then add on len_prim and len_dump
        in_s(1,Im_ident,isec,iitm)=len_extra(Sm_ident)+1
        len_extra(Sm_ident) = len_extra(Sm_ident)+LEN
        ! Information for preliminary D1 addressing array
        n_obj_d1(Sm_ident)     =n_obj_d1(Sm_ident)+1
        IF (n_obj_d1(Sm_ident) <= max_d1_len) THEN
          d1_paddr(d1_type,n_obj_d1(Sm_ident),Sm_ident)=extra_d1
          d1_paddr(d1_im,n_obj_d1(Sm_ident),Sm_ident)=Im_ident
          d1_paddr(d1_extra_info,n_obj_d1(Sm_ident),Sm_ident)=iitm
          d1_paddr(d1_levs,n_obj_d1(Sm_ident),Sm_ident)=dlevs
          d1_paddr(d1_sect,n_obj_d1(Sm_ident),Sm_ident)=isec
        END IF
      END IF
    ELSE
      ! ISP=5 means: set address of prim var in dump only.
      ! D1 address is then set to same address as previous item
      IF (iitm  ==  1) THEN
        in_s(1,Im_ident,isec,iitm)=1
      ELSE
        in_s(1,Im_ident,isec,iitm)=in_s(1,Im_ident,isec,iitm-1)
      END IF
    END IF
    ! The input length for primary (m,0,i) is assigned to IN_S(2,m,0,i).
    in_s(2,Im_ident,isec,iitm)=LEN

  END IF  ! ISPACE  /=  10
END IF   ! LADDR

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE primary

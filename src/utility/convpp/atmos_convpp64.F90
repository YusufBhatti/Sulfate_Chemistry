! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    -----------------------------------------------------------------
!    SUBROUTINE ATMOS_CONVPP64 -----------------------------------------
!
!    Purpose: Converts UM 64 bit file to 32 bit PP format.
!
!    Programming standards: UMDP 3
!
!    Documentation: UM Doc Paper F5
!
!    System Tasks: F3,F4,F6
!
!    -----------------------------------------------------------------
!    Arguments:-------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs

SUBROUTINE atmos_convpp64                                        &
 (len_fixhd,len_inthd,len_realhd,                                &
 len1_levdepc,len2_levdepc,len1_rowdepc,                         &
 len2_rowdepc,len1_coldepc,len2_coldepc,                         &
 len1_flddepc,len2_flddepc,len_extcnst,                          &
 len_dumphist,len_cfi1,len_cfi2,len_cfi3,                        &
 len1_lookup,len2_lookup,len_data,p_field,                       &
 nftin,nftout,max_field_size)

USE pp_header_manips, ONLY: header_manip
USE mask_compression, ONLY: expand_from_mask
USE UM_ParVars
USE ereport_mod, ONLY: ereport
USE Decomp_DB
USE missing_data_mod, ONLY: rmdi, imdi
USE ppxlook_mod, ONLY: read_atmos_stashmaster,exppxi
USE Submodel_Mod
USE field_types, ONLY: fld_type_u, fld_type_v
USE cppxref_mod, ONLY: ppx_grid_type
USE errormessagelength_mod, ONLY: errormessagelength
USE UM_Types, ONLY: real64, real32, integer64, integer32
USE lookup_addresses
USE readflds_mod, ONLY: readflds
USE umPrintMgr
USE packing_codes_mod, ONLY: PC_No_Packing, PC_WGDOS_Packing
IMPLICIT NONE

INTEGER ::                                                        &

 len_fixhd                                                        &
              !IN Length of fixed length header on input file
,len_inthd                                                        &
              !IN Length of integer header on input file
,len_realhd                                                       &
              !IN Length of real header on input file
,len1_levdepc                                                     &
              !IN 1st dim of lev dependent consts on input file
,len2_levdepc                                                     &
              !IN 2nd dim of lev dependent consts on input file
,len1_rowdepc                                                     &
              !IN 1st dim of row dependent consts on input file
,len2_rowdepc                                                     &
              !IN 2nd dim of row dependent consts on input file
,len1_coldepc                                                     &
              !IN 1st dim of col dependent consts on input file
,len2_coldepc                                                     &
              !IN 2nd dim of col dependent consts on input file
,len1_flddepc                                                     &
              !IN 1st dim of field dependent consts on input fi
,len2_flddepc                                                     &
              !IN 2nd dim of field dependent consts on input fi
,len_extcnst                                                      &
              !IN Length of extra consts on input file
,len_dumphist                                                     &
              !IN Length of history header on input file
,len_cfi1                                                         &
              !IN Length of index1 on input file
,len_cfi2                                                         &
              !IN Length of index2 on input file
,len_cfi3                                                         &
              !IN Length of index3 on input file
,len1_lookup                                                      &
              !IN 1st dim of LOOKUP on input file
,len2_lookup                                                      &
              !IN 2nd dim of LOOKUP on input file
,len_data                                                         &
              !IN Length of data on input file
,p_field                                                          &
              !IN No of p-points per level on input file
,max_field_size !Maximum field size on file

INTEGER :: nftin
INTEGER :: nftout

! ----------------------------------------------------------------------
! Local variables:
! ----------------------------------------------------------------------

! Constants: -----------------------------------------------------------
INTEGER, PARAMETER :: lenPpHead = 64, lenIntHead = 45, block_sz = 512
INTEGER, PARAMETER :: nreals = lenPpHead - lenIntHead

! Local scalars: -------------------------------------------------------
INTEGER :: is, im          ! STASH section and item 
INTEGER :: grid_type_code  ! grid type code
INTEGER :: gr              ! grid type of field being processed 
INTEGER :: fld_type        ! field type: u-,v- or p- location on C grid

! Input arguments for decompose_smexe:
INTEGER :: global_row_len  ! IN  :number of E-W points of entire model
INTEGER :: global_n_rows   ! IN  :number of P rows of entire mode
INTEGER :: tot_levels      ! IN  :total number of levels

INTEGER :: i,j,k,l         ! Loop indices
INTEGER :: irows           ! Number of points north-south
INTEGER :: icols           ! Number of points east-west
INTEGER :: used_imdi       ! Integer missing data indicator
INTEGER :: nrows           ! Number of points north-south  (p grid)
INTEGER :: ncols           ! Number of points east-west
INTEGER :: i_skipped       ! Number of unsupported fields skipped
INTEGER :: lblrec_in       ! keep a copy of the original length of data
INTEGER :: lbext_in        ! keep a copy of lookup(20)
INTEGER :: icode           ! Error return code from subroutines    
INTEGER :: start_block     ! READHEAD argument (not used)
INTEGER :: get_fld_type    ! Identifed field type from function 

! Variable grid extra data codes  
INTEGER :: xcentrecode
INTEGER :: ycentrecode
INTEGER :: xleftcode
INTEGER :: xrightcode
INTEGER :: ybotcode
INTEGER :: ytopcode 

! 32 bit variable grid extra data codes
INTEGER(KIND=integer32) :: xcentrecode_32
INTEGER(KIND=integer32) :: ycentrecode_32
INTEGER(KIND=integer32) :: xleftcode_32
INTEGER(KIND=integer32) :: xrightcode_32
INTEGER(KIND=integer32) :: ybotcode_32
INTEGER(KIND=integer32) :: ytopcode_32 

! Scalars for data transfer to 32bits: 
INTEGER(KIND=integer32) :: int32_value
INTEGER                 :: int64_value

REAL                    :: real64_value, val_64
REAL                    :: used_rmdi

LOGICAL :: land_sea_mask(max_field_size)
LOGICAL :: variable_grid   ! True if input is a Variable Grid 

! Local arrays: --------------------------------------------------------
! Integer file headers
INTEGER :: fixhd(len_fixhd)
INTEGER :: inthd(len_inthd)
INTEGER :: cfi1(len_cfi1+1)
INTEGER :: cfi2(len_cfi2+1)
INTEGER :: cfi3(len_cfi3+1)
INTEGER :: lookup(len1_lookup,len2_lookup)        ! Lookup table
INTEGER :: lookup_out(len1_lookup)                ! Output lookup table

! Arrays for data transfer to 32bits:
INTEGER(KIND=integer32) :: lookup_32int(lenPpHead)
INTEGER(KIND=integer32) :: tmp_int_header(nreals)

INTEGER(KIND=integer32) :: xc_32int(len1_coldepc) ! centre
INTEGER(KIND=integer32) :: xl_32int(len1_coldepc) ! left
INTEGER(KIND=integer32) :: xr_32int(len1_coldepc) ! right
INTEGER(KIND=integer32) :: yc_32int(len1_rowdepc) ! centre
INTEGER(KIND=integer32) :: yt_32int(len1_rowdepc) ! top
INTEGER(KIND=integer32) :: yb_32int(len1_rowdepc) ! bottom

! Real arrays
REAL    :: realhd(len_realhd)
REAL    :: levdepc(1+len1_levdepc*len2_levdepc)
REAL    :: rowdepc(1+len1_rowdepc*len2_rowdepc)
REAL    :: coldepc(1+len1_coldepc*len2_coldepc)   ! Real file headers
REAL    :: flddepc(1+len1_flddepc*len2_flddepc)   ! Real file headers
REAL    :: extcnst(len_extcnst+1)                 ! Real file headers
REAL    :: dumphist(len_dumphist+1)
REAL    :: d1(max_field_size)  ! Data array used to read in each field
REAL    :: d1_tmp(max_field_size)
REAL    :: land_sea_mask_tmp(max_field_size)

! Arrays for data transfer to 32bits:
REAL(KIND=real32)       :: lookup_32real(nreals)
REAL(KIND=real32)       :: d1_32(max_field_size)

! Variable grid extra data codes                    
REAL    :: xcentre(len1_coldepc)
REAL    :: xleft(len1_coldepc)
REAL    :: xright(len1_coldepc)
REAL    :: ycentre(len1_rowdepc)
REAL    :: ytop(len1_rowdepc)
REAL    :: ybot(len1_rowdepc) 

LOGICAL :: l_endgame

CHARACTER(LEN=errormessagelength) :: cmessage                
CHARACTER(LEN=20) :: string    ! Format control for header printout

! ----------------------------------------------------------------------

!  0. Read in PPXREF

cmessage = ' '
icode=0

CALL read_atmos_stashmaster()
!  1. Read in file header

! DEPENDS ON: readhead
CALL readhead(nftin,fixhd,len_fixhd,                              &
                inthd,len_inthd,                                  &
                realhd,len_realhd,                                &
                levdepc,len1_levdepc,len2_levdepc,                &
                rowdepc,len1_rowdepc,len2_rowdepc,                &
                coldepc,len1_coldepc,len2_coldepc,                &
                flddepc,len1_flddepc,len2_flddepc,                &
                extcnst,len_extcnst,                              &
                dumphist,len_dumphist,                            &
                cfi1,len_cfi1,                                    &
                cfi2,len_cfi2,                                    &
                cfi3,len_cfi3,                                    &
                lookup,len1_lookup,len2_lookup,                   &
                len_data,                                         &
                start_block,icode,cmessage)

IF (icode /= 0) THEN
  WRITE(umMessage,*)cmessage,icode
  CALL umPrint(umMessage,src='atmos_convpp.F90')
  CALL ereport('ATMOS_CONVPP', icode, cmessage)
END IF

nrows        = inthd(7)
ncols        = inthd(6)
irows        = inthd(7)
icols        = inthd(6)

l_endgame = .FALSE.
IF (fixhd(9) == 6) THEN
  l_endgame = .TRUE.
END IF

IF (len_realhd >= 29) THEN
  used_rmdi         = realhd(29)
ELSE
  used_rmdi         = rmdi
END IF

IF (len_inthd >= 21) THEN
  used_imdi         = inthd(21)
ELSE
  used_imdi         = imdi
END IF

! This currently assumes ALL data in file are on VR grid if header present. 
variable_grid=.FALSE.  
 
IF ( fixhd(115) /= imdi .AND. fixhd(120) /=imdi ) THEN  
  IF (len1_rowdepc > 0 .AND. len2_rowdepc > 0 .AND.               &  
     len1_coldepc > 0 .AND. len2_coldepc > 0) THEN  
    variable_grid=.TRUE.  
    CALL umPrint('Input data is on variable grid.',               &
                                    src='atmos_convpp.F90') 
  END IF   
END IF


! Get decomposition information
! -----------------------------
tot_levels=-99

CALL decompose(icols,irows,0,0,tot_levels)
CALL change_decomposition(decomp_smexe,icode)

DO i=1,len2_lookup
! If item 1:0:30 LAND MASK (No halo) (LAND=TRUE)
  IF (lookup(item_code,i) == 30) THEN

    CALL readflds(nftin,1,i,lookup,                   &
                  land_sea_mask_tmp,fixhd,0,icode,cmessage)

    land_sea_mask = TRANSFER(land_sea_mask_tmp,land_sea_mask)
    ! DEPENDS ON: abort_io
    IF (icode /= 0)CALL abort_io('CONVPP',cmessage,icode,nftin)
  END IF
END DO

i_skipped=0
!   Print out individual fields
DO i=1,len2_lookup
  IF (lookup(1,i) == -99) EXIT

  ! Where LOOKUP is mainly IMDI (invalid), ignore conversion.
  IF (lookup(16,i) == used_imdi) THEN
    WRITE(umMessage,'(a,i5,a)')                                   &
        'File included unsupported fields. Field ', i, ' omitted.'
    CALL umPrint(umMessage,src='atmos_convpp.F90')
    i_skipped=i_skipped+1
    CYCLE
  ELSE IF (lookup(16,i) < 0) THEN
    ! This might have used another number for MDI.  Really any negative
    ! number will be an error.
    WRITE(umMessage,'(a,i5)') &
        'Possibly detected unsupported field. Field ', i
    CALL umPrint(umMessage,src='atmos_convpp.F90')
  END IF

  ! Fill output lookup table for a particular stash item
  DO k = 1, len1_lookup
    lookup_out(k) = lookup(k,i)
  END DO

  CALL readflds(nftin,1,i,lookup,                     &
                d1,fixhd,0,icode,cmessage)

  lookup_out(21)=MOD(lookup_out(21),1000)
  IF ((lookup_out(21)/10)*10 == 120) THEN
    !Data compressed on to land points.
    !Copy data to temporary array
    DO k=1,lookup_out(15)
      d1_tmp(k)=d1(k)
    END DO
    !Set unpacked array to RMDI
    DO k=1,nrows*ncols
      d1(k)=used_rmdi
    END DO

    !Uncompress data
    CALL expand_from_mask(d1,d1_tmp,land_sea_mask,                &
                          max_field_size,lookup_out(15))
    lookup_out(15)=nrows*ncols
    lookup_out(18)=nrows
    lookup_out(19)=ncols
    lookup_out(21)=0
  END IF

  !-----------------------------------------------------------
  ! set up lookups for variable resolution grid if applicable.
  !-----------------------------------------------------------

  ! Take a copy of length of data record before extra data length is added
  ! to it and a copy of length of extra data ("lblrec" and "lbext")
  lblrec_in = lookup_out(15)
  lbext_in = lookup_out(20)  
  IF (variable_grid .AND. lookup_out(20) == 0) THEN

    IS=lookup_out(42)/1000 
    im=lookup_out(42)-IS*1000 
    gr = exppxi( 1, IS, im, ppx_grid_type,                        & 
                 icode, cmessage) 
    grid_type_code=gr

    ! Get field type, ie u or v or p location in Arakawa C grid staggering 
    ! DEPENDS ON: get_fld_type 
    fld_type=get_fld_type(grid_type_code) 
         
    nrows=lookup_out(18) 
    ncols=lookup_out(19)
    ! lookup(15) is the length of data record including any extra data
    ! lookup(20) length of extra data
    ! For variable grid extra data provided includes:
    ! Vector 1 (xcentrecode) the x coord values p lambda points
    ! Vector 2 (ycentrecode) the y coords values p phi points
    ! vectors 12,13 x left/right code lambda points for each grid box
    ! vectors 14,15 y bot/top code phi points for each grid box
    ! Thus we end up with
    lookup_out(20)=3*(nrows+ncols+2) 
    xcentrecode=1000*ncols+1
    xleftcode=1000*ncols+12
    xrightcode=1000*ncols+13
    ycentrecode=1000*nrows+2
    ybotcode=1000*nrows+14
    ytopcode=1000*nrows+15

    ! Set up arrays of grid points to eventually write to extra data in 
    ! pp lookup. 
    !
    ! Here we define the p-,u- and v-grids and their bounding boxes,so:
    !
    ! p-grid top  and bottom (phi)    bounds would be on v-grid phis 
    ! p-grid left and right  (lambda) bounds would be u-grid lambdas 
    ! 
    ! Conversely:
    ! v-grid  top and bottom (phi)     bounds are on p-grid phis 
    ! u-grid left and rights (lambdas) bounds are on p-grid lambdas.
    ! 
    ! Note: In this routine the bounding boxes array here is the phi/lamba
    !        location of the bounding box for each gridpoint
    !
    ! Its assumed that required grid data info is placed accordingly ... 
    ! in FF header
    ! p-grid pts are coldepc(1:len1_coldepc)
    !                rowsdepc(1:len1_rowdepc)
    ! u-grid pts are coldepc(len1_coldepc+1:2*len1_coldepc)
    !                rowdepc(rowsdepc(1:len1_rowdepc)
    ! v-grid pts are coldepc(1:len1_coldepc)
    !                rowdepc(len1_rowdepc+1:2*len1_rowdepc)

    ! Set up grid arrays initially to RMDI
    DO k=1,ncols
      xcentre(k)=used_rmdi
      xleft(k)  =used_rmdi
      xright(k) =used_rmdi
    END DO
    DO k=1,nrows
      ycentre(k)=used_rmdi
      ybot(k)   =used_rmdi
      ytop(k)   =used_rmdi
    END DO

    ! Starting with Lambdas (longitudes/columns) 
    ! ND u and p grid have same number of columns in horiz grid definition.  
    !                   P/V - U - P/V - U
    ! EG u shifts to the LHS of domain and forgoes RHS                   
    !               U - P/V - U - P/V 

    IF (fld_type == fld_type_u) THEN    
          
      ! Set up u-grid lambda information   
      DO k=1,ncols 
        xcentre(k)=coldepc(len1_coldepc+k)  
      END DO 
           

      IF (l_endgame) THEN

        ! Set up bounding box assuming p-grid pts mark left boundary
        DO k=1,ncols-1 
          xleft(k+1)=coldepc(k) 
        END DO

        ! Need an extra column of p- pts to bound the first u-grid pt.
        xleft(1)=2.0*coldepc(len1_coldepc+1)               & 
                      -coldepc(1) 

        ! Set up bounding box assuming p-grid pts mark right boundary
        DO k=1,ncols
          xright(k)=coldepc(k) 
        END DO

      ELSE  ! New dynamics

        ! Set up bounding box assuming p-grid pts mark left boundary
        DO k=1,ncols 
          xleft(k)=coldepc(k) 
        END DO 
           
        ! Set up bounding box assuming p-grid pts mark right boundary
        DO k=1,ncols-1
          xright(k)=coldepc(k+1) 
        END DO
           
        ! Need an extra column of p-grid pts to bound the last u-grid pt
        xright(ncols)=2.0*coldepc(len1_coldepc+ncols)               & 
                      -coldepc(ncols) 
 
      END IF

    ELSE          
         
      ! Set p- or v-grid lambda information
      DO k=1,ncols 
        xcentre(k)=coldepc(k) 
      END DO 
           
      IF (l_endgame) THEN

        ! Set up bounding box assuming u-grid pts mark left boundary

        DO k=1,ncols 
          xleft(k)=coldepc(len1_coldepc+k)  
        END DO 
           
        ! Set up bounding box assuming u-grid pts mark right boundary
        DO k=1,ncols-1 
          xright(k)=coldepc(len1_coldepc+k+1) 
        END DO

        ! Need an extra column of u-grid pts to bound the last p-grid pt
        xright(ncols)=2.0*coldepc(ncols)               & 
                      -coldepc(len1_coldepc+ncols) 

      ELSE   ! New dynamics

        ! Set up bounding box assuming u-grid pts mark left boundary
        xleft(1)=2.0*coldepc(1)-coldepc(len1_coldepc+1) 
 
        DO k=1,ncols 
          xleft(k+1)=coldepc(len1_coldepc+k)  
        END DO 
           
        ! Set up bounding box assuming u-grid pts mark right boundary
        DO k=1,ncols 
          xright(k)=coldepc(len1_coldepc+k) 
        END DO
      END IF 
    END IF
                  

    ! Starting with phis (latitudes/rows)  
    ! In ND model v_nrows is p_nrows-1 
    ! The UM writes out RMDI as the last grid def info for ND v grid noting 
    ! it is smaller 
    ! For EG v_nrows is p_nrows+1

    !    ND      EG
    !            V
    !    P/U    P/U
    !     V      V
    !    P/U    P/U
    !     V      V
    !    P/U    P/U
    !            V

    IF (fld_type == fld_type_v) THEN 
          
      ! Set up v-grid phi information
      DO k=1,nrows 
        ycentre(k)=rowdepc(len1_rowdepc+k) 
      END DO 
          
      IF (l_endgame) THEN

        ! Set up bounding box assuming p/u-grid pts mark bottom boundary
        DO k=1,nrows-1 
          ybot(k+1)=rowdepc(k) 
        END DO 

        ! Set up bounding box assuming p/u-grid pts mark bottom boundary
        ybot(1)=2.0*rowdepc(len1_rowdepc+1)-rowdepc(1) 

        ! Set up bounding box assuming p/u-grid pts mark top boundary
        DO k=1,nrows-1 
          ytop(k)=rowdepc(k) 
        END DO

        ! Add on final p/u-grid value.
        ytop(nrows)=2.0*rowdepc(len1_rowdepc+nrows)                 & 
                    -rowdepc(nrows) 

      ELSE      ! New dynamics
 
        ! Set up bounding box assuming p/u-grid pts mark bottom boundary
        DO k=1,nrows 
          ybot(k)=rowdepc(k) 
        END DO 
           
        ! Set up bounding box assuming p/u-grid pts mark top boundary
        DO k=1,nrows-1 
          ytop(k)=rowdepc(k+1) 
        END DO
           
        ! Add on final p/u-grid value.
        ytop(nrows)=2.0*rowdepc(len1_rowdepc+nrows)                 & 
                    -rowdepc(nrows) 
      END IF
    ELSE          

      ! Set up p- or u-grid phi information  
      ! (nrows here 1 more [nd] than v-grid above)
      DO k=1,nrows 
        ycentre(k)=rowdepc(k) 
      END DO 
           
      IF (l_endgame) THEN
        ! Set up bounding box assuming v-grid pts mark bottom boundary

        DO k=1,nrows 
          ybot(k)=rowdepc(len1_rowdepc+k)  
        END DO 
           
        ! Set up bounding box assuming v-grid pts mark top boundary 
        DO k=1,nrows
          ytop(k)=rowdepc(len1_rowdepc+k+1) 
        END DO

      ELSE    ! New dynamics

        ! Set up bounding box assuming v-grid pts mark bottom boundary
        ybot(1)=2.0*rowdepc(1)-rowdepc(len1_rowdepc+1) 

        DO k=2,nrows 
          ybot(k)=rowdepc(len1_rowdepc+k-1)  
        END DO 
           
        ! Set up bounding box assuming v-grid pts mark top boundary 
        DO k=1,nrows 
          ytop(k)=rowdepc(len1_rowdepc+k) 
        END DO
            
        ! Set up final v-grid bounding box with a sensible value
        ytop(nrows)=2.0*rowdepc(len1_rowdepc)-rowdepc(len1_rowdepc+nrows-1) 
      END IF
    END IF

    ! Set standard grid definitions BDX and BDY to RMDI as using VR
    lookup_out(59)=TRANSFER(rmdi,lookup_out(59))  ! BZY
    lookup_out(60)=TRANSFER(rmdi,lookup_out(60))  ! BDY
    lookup_out(61)=TRANSFER(rmdi,lookup_out(61))  ! BZX
    lookup_out(62)=TRANSFER(rmdi,lookup_out(62))  ! BDX        

    ! In inside brackets 64bit reals converted to 32bit reals. The value is
    ! not truncated but rounded up. Due to the fact that 64bit exe
    ! uses format s(11 exponent)(52 fraction) instead of 32bit format
    ! s(8 exponent)(23 fraction), grid values will be represented like below
    ! xc_32int(1)=  353.07196 xcentre(1)=  353.07195097200002
    ! 32bit  1    0000111 01100001000100100110110
    ! 64bit  1 0000000111 0110000100010010011010110110000011111111101001111101
    DO k=1,ncols
      xc_32int(k) = TRANSFER(REAL(xcentre(k), KIND=real32), int32_value)
      xl_32int(k) = TRANSFER(REAL(xleft(k), KIND=real32), int32_value)
      xr_32int(k) = TRANSFER(REAL(xright(k), KIND=real32), int32_value)
    END DO
      
    DO k=1,nrows
      yc_32int(k) = TRANSFER(REAL(ycentre(k), KIND=real32), int32_value)
      yt_32int(k) = TRANSFER(REAL(ytop(k), KIND=real32), int32_value)
      yb_32int(k) = TRANSFER(REAL(ybot(k), KIND=real32), int32_value)
    END DO

    ! Transform 64bit integers to 32bit ones   
    xcentrecode_32 = INT(xcentrecode, KIND=integer32)
    ycentrecode_32 = INT(ycentrecode, KIND=integer32)
    xleftcode_32 = INT(xleftcode, KIND=integer32)
    xrightcode_32 = INT(xrightcode, KIND=integer32)
    ybotcode_32 = INT(ybotcode, KIND=integer32)
    ytopcode_32 = INT(ytopcode, KIND=integer32)
  END IF ! case: variable grid and extra data
    
  CALL header_manip(lookup_out(1:lenIntHead)) 

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Reformat lookup field header
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Adjust lblrec and lbnrec fields and set lbegin to 0
  ! transform to 32 bits 45 integer lookup descriptors
  ! transform to 32 bits real numbers of lookup descriptors
  ! concatinate integer and real discriptors

  ! Fill output integer part of lookup table
  DO k=1,lenIntHead
    lookup_32int(k)=INT(lookup_out(k), KIND=integer32)
  END DO

  ! Correct lblrec and lbnrec for packed fileds
  IF (MOD(lookup_out(lbpack),10)  ==  PC_WGDOS_Packing) THEN
    ! WGDOS packed data as read-in from the FieldsFile is treated as a
    ! 64-bit array, so in order to get the addressing header right the
    ! record length has to be doubled.
    lookup_32int(lblrec) = 2*lookup_32int(lblrec)

    ! Calculate the length of field on disk in words (plus padding) "lbnrec"
    ! In case lbnrec == 0 -z option, no need to recalculate this field 
    IF (lookup_32int(lbnrec) /= 0) THEN
      lookup_32int(lbnrec) = (lookup_32int(lblrec) / block_sz + 1) * block_sz
    END IF
  END IF

  ! If data is on variable grid, then we need to add extra data size
  ! to the total length
  IF (variable_grid .AND. lbext_in == 0) THEN
    lookup_32int(lblrec) = lookup_32int(lblrec) + lookup_32int(lbext)
  END IF

  ! Transfer from lookup table real numbers into separate 32 bits real array
  DO k = lenIntHead + 1,lenPpHead
    val_64 = TRANSFER(lookup_out(k),real64_value)
    lookup_32int(k) = TRANSFER(REAL(val_64, KIND=real32), int32_value)
  END DO
  
  ! Write lookup field header
  WRITE(nftout)(lookup_32int(k),k=1,lenPpHead)

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Write data fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Transfer unpacked data to 32bits
  ! Fortran compiler requires kind definition for integer literal
  IF (MOD(lookup_32int(lbpack),10_integer32)  == PC_No_Packing) THEN
    ! Integer or logic data type
    IF (ABS(lookup_32int(data_type)) == 2 .OR.        &
        ABS(lookup_32int(data_type)) == 3) THEN
      DO k=1,lookup_32int(lblrec)
        int64_value = TRANSFER(D1(k),int64_value)
        int32_value = INT(int64_value, KIND=integer32)
        d1_32(k) = TRANSFER(int32_value,d1_32(k))
      END DO
    ELSE
      ! Real data type 
      DO k=1,lookup_32int(lblrec)
        d1_32(k) = REAL(d1(k), KIND=real32) 
      END DO
    END IF
  END IF

  ! For variable grid
  IF (variable_grid .AND. lbext_in == 0) THEN
    IF ((MOD(lookup_out(lbpack),10)  ==  PC_No_Packing)) THEN
      WRITE(nftout)(d1_32(k),k=1,lblrec_in),                 & 
              xcentrecode_32,(xc_32int(k),k=1,ncols),        & 
              ycentrecode_32,(yc_32int(k),k=1,nrows),        & 
              xleftcode_32,(xl_32int(k),k=1,ncols),          & 
              xrightcode_32,(xr_32int(k),k=1,ncols),         & 
              ybotcode_32,(yb_32int(k),k=1,nrows),           & 
              ytopcode_32,(yt_32int(k),k=1,nrows)     

    ELSE
      WRITE(nftout)(d1(k),k=1,lblrec_in),                    & 
              xcentrecode_32,(xc_32int(k),k=1,ncols),        & 
              ycentrecode_32,(yc_32int(k),k=1,nrows),        & 
              xleftcode_32,(xl_32int(k),k=1,ncols),          & 
              xrightcode_32,(xr_32int(k),k=1,ncols),         & 
              ybotcode_32,(yb_32int(k),k=1,nrows),           & 
              ytopcode_32,(yt_32int(k),k=1,nrows)     

    END iF
  ELSE
    IF (MOD(lookup_out(lbpack),10)  ==  PC_No_Packing) THEN
      WRITE(nftout) (d1_32(k),k=1,lblrec_in)
    ELSE
      WRITE(nftout) (d1(k),k=1,lblrec_in) 
    END IF
  END IF

END DO

WRITE(umMessage,'(I0,A)')i-1-i_skipped,' pp fields written out'
CALL umPrint(umMessage,src='atmos_convpp.F90')

RETURN
END SUBROUTINE atmos_convpp64

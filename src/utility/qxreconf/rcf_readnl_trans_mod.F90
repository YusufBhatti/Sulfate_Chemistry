! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Reads the TRANS namelists

MODULE Rcf_ReadNL_Trans_Mod

IMPLICIT NONE

!  Subroutine Rcf_Readnl_Trans - reads the TRANS namelists
!
! Description:
!   Reads the TRANS namelist for controling transplanting fields.
!
! Method:
!   Reads into temporary arrays and then copies into the correct
!   dynamically allocated arrays.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_TRANS_MOD'

CONTAINS

SUBROUTINE Rcf_ReadNL_Trans( nft )

USE Rcf_Trans_Mod         ! All of it

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Oper

USE Ereport_Mod, ONLY: Ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE UM_ParCore, ONLY: mype

USE mpl, ONLY: mpl_integer

USE setup_namelist, ONLY: setup_nml_type

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)       :: nft      ! unit number

! Local variables
INTEGER, PARAMETER        :: max_trans = 150 ! maximum number of items
INTEGER                   :: iostatus
INTEGER :: itemc
INTEGER :: sctnc
INTEGER :: lev1
INTEGER :: lev2
INTEGER :: row1
INTEGER :: row2
INTEGER :: col1
INTEGER :: col2


! temporary arrays
INTEGER                   :: itemc_temp( max_trans )
INTEGER                   :: sctnc_temp( max_trans )
INTEGER                   :: lev1_temp( max_trans )
INTEGER                   :: lev2_temp( max_trans )
INTEGER                   :: row1_temp( max_trans )
INTEGER                   :: row2_temp( max_trans )
INTEGER                   :: col1_temp( max_trans )
INTEGER                   :: col2_temp( max_trans )

INTEGER                   :: ErrorStatus
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_READNL_TRANS'
CHARACTER (LEN=errormessagelength) :: Cmessage
CHARACTER (LEN=50000)         :: locMessage
NAMELIST /Trans/ sctnc, itemc, lev1, lev2, col1, col2, row1, row2

INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_int = 8*max_trans
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

TYPE my_namelist
  SEQUENCE
  INTEGER :: sctnc_temp ( max_trans )
  INTEGER :: itemc_temp ( max_trans )
  INTEGER :: lev1_temp ( max_trans )
  INTEGER :: lev2_temp ( max_trans )
  INTEGER :: col1_temp ( max_trans )
  INTEGER :: col2_temp ( max_trans )
  INTEGER :: row1_temp ( max_trans )
  INTEGER :: row2_temp ( max_trans )
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!Namelist /Trans/ sctnc, itemc, lev1, lev2, col1, col2, row1, row2

!-----------------------------------------------------------------
! Not all platforms will necessarily require this rewind
!-----------------------------------------------------------------
CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int)

IF (mype == 0) THEN

  !-----------------------------------------------------------------
  ! Set num_trans to 0, loop for reading in namelist
  !-----------------------------------------------------------------

  num_trans = 0
  iostatus  = 0

  DO WHILE ( iostatus == 0 )
    READ( UNIT=nft, NML=trans, IOSTAT=iostatus )

    IF (iostatus == 0 ) THEN
      num_trans = num_trans + 1

      IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
        ! FIXME namelist
        !       Write ( Unit = 6, Nml = trans )
      END IF

      IF ( num_trans > max_trans ) THEN
        Cmessage = 'Maximum number of transplanted'//&
          ' fields exceeded - please use a mod to increase'
        WRITE(umMessage,*) Cmessage
        CALL umPrint(umMessage,src='rcf_readnl_trans_mod')
        WRITE(umMessage,*) 'Number     of TRANS namelists read in : ', &
                            num_trans
        CALL umPrint(umMessage,src='rcf_readnl_trans_mod')
        WRITE(umMessage,*) 'Maximum no of TRANS namelists allowed : ', &
                            max_trans
        CALL umPrint(umMessage,src='rcf_readnl_trans_mod')
        ErrorStatus = 10
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )
      END IF

      IF (lev1 == 0 .OR. lev2 == 0 .OR. row1 == 0 .OR. row2 == 0 .OR. &
          col1 == 0 .OR. col2 == 0 ) THEN
        WRITE(umMessage,*) 'Namelist TRANS has some zero values'
        CALL umPrint(umMessage,src='rcf_readnl_trans_mod')
        WRITE(umMessage,*) 'lev1 = ', lev1, ' lev2 = ', lev2
        CALL umPrint(umMessage,src='rcf_readnl_trans_mod')
        WRITE(umMessage,*) 'row1 = ', row1, ' row2 = ', row2
        CALL umPrint(umMessage,src='rcf_readnl_trans_mod')
        WRITE(umMessage,*) 'col1 = ', col1, ' col2 = ', col2
        CALL umPrint(umMessage,src='rcf_readnl_trans_mod')

        Cmessage = 'Zero values in TRANS namelist'
        ErrorStatus = 20
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )
      END IF

      sctnc_temp( num_trans ) = sctnc
      itemc_temp( num_trans ) = itemc
      lev1_temp( num_trans ) = lev1
      lev2_temp( num_trans ) = lev2
      col1_temp( num_trans ) = col1
      col2_temp( num_trans ) = col2
      row1_temp( num_trans ) = row1
      row2_temp( num_trans ) = row2
    END IF

  END DO

  IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE(umMessage,*) 'Num_Trans = ', num_trans
    CALL umPrint(umMessage,src='rcf_readnl_trans_mod')
  END IF

END IF  !! mype==0

CALL mpl_bcast(num_trans,1,mpl_integer,0,my_comm,icode)

!-----------------------------------------------------------------
! Allocate proper space for trans list and copy data
!-----------------------------------------------------------------
ALLOCATE( sctnc_array( num_trans ) )
ALLOCATE( itemc_array( num_trans ) )
ALLOCATE( lev1_array( num_trans ) )
ALLOCATE( lev2_array( num_trans ) )
ALLOCATE( col1_array( num_trans ) )
ALLOCATE( col2_array( num_trans ) )
ALLOCATE( row1_array( num_trans ) )
ALLOCATE( row2_array( num_trans ) )

IF (mype==0) THEN

  my_nml % sctnc_temp = sctnc_temp
  my_nml % itemc_temp = itemc_temp
  my_nml % lev1_temp  = lev1_temp
  my_nml % lev2_temp  = lev2_temp
  my_nml % col1_temp  = col1_temp
  my_nml % col2_temp  = col2_temp
  my_nml % row1_temp  = row1_temp
  my_nml % row2_temp  = row2_temp

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  sctnc_temp  = my_nml % sctnc_temp
  itemc_temp  = my_nml % itemc_temp
  lev1_temp   = my_nml % lev1_temp
  lev2_temp   = my_nml % lev2_temp
  col1_temp   = my_nml % col1_temp
  col2_temp   = my_nml % col2_temp
  row1_temp   = my_nml % row1_temp
  row2_temp   = my_nml % row2_temp

END IF

sctnc_array( 1 : num_trans ) = sctnc_temp( 1: num_trans )
itemc_array( 1 : num_trans ) = itemc_temp( 1: num_trans )
lev1_array( 1 : num_trans ) = lev1_temp( 1: num_trans )
lev2_array( 1 : num_trans ) = lev2_temp( 1: num_trans )
col1_array( 1 : num_trans ) = col1_temp( 1: num_trans )
col2_array( 1 : num_trans ) = col2_temp( 1: num_trans )
row1_array( 1 : num_trans ) = row1_temp( 1: num_trans )
row2_array( 1 : num_trans ) = row2_temp( 1: num_trans )

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_ReadNL_Trans

END MODULE Rcf_ReadNL_Trans_Mod

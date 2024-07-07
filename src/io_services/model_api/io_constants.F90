! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

MODULE io_constants
! A module to contain the PARAMETERs for file operations

IMPLICIT NONE
PUBLIC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters indicating the nature of the file open operation
INTEGER, PARAMETER :: ioOpenReadOnly        = 0
INTEGER, PARAMETER :: ioOpenReadWrite       = 1
INTEGER, PARAMETER :: ioOpenWriteOnly       = -99999

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters indicating known file types, for provision via the fileType
! optional argument to file_open (and internally)
INTEGER, PARAMETER :: ioFileTypeUnknown = 0
INTEGER, PARAMETER :: ioFileTypeUM      = 1
INTEGER, PARAMETER :: ioFileTypeMPIIO   = 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters indicating how the passed string to file_open/close encodes
! the filename.
INTEGER, PARAMETER :: ioNameInEnv       = 0
INTEGER, PARAMETER :: ioNameProvided    = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters indicating whether file should be deleted on close.
INTEGER, PARAMETER :: ioNoDelete        = 0
INTEGER, PARAMETER :: ioDelete          = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters for controlling where IO happens !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER, PARAMETER :: ioNotKnown=-1
INTEGER, PARAMETER :: ioNotSet=0
INTEGER, PARAMETER :: ioLocal=1
INTEGER, PARAMETER :: ioAllLocal=2
INTEGER, PARAMETER :: ioRemote=4
INTEGER, PARAMETER :: ioReadReplicate=8
INTEGER, PARAMETER :: ioReadOnly=16
INTEGER, PARAMETER :: ioOpen=32
INTEGER, PARAMETER :: ioWriteOnly=64
!Some shortcuts
   ! Good choice for wr files;
INTEGER, PARAMETER :: ioRemoteReplicate=IOR(ioRemote,ioReadReplicate)
 ! Good choice for ro files;
INTEGER, PARAMETER :: ioLocalReplicate=IOR (ioLocal, ioReadReplicate)
! The default
INTEGER, PARAMETER :: ioDefault=ioRemoteReplicate
!Notes:
! The parameter is supplied at file open and persists until the
! file is closed
!
! ioAllLocal : all tasks will open/write/close the file
!  - It is the programmers responsibility to avoid filename collisions
!  - It is the programmers responsibility to make calls from appropriate PEs
! ioLocal : Rank 0 will open/write/close the file
!  - Other tasks will return with a success flag but not do anything
! ioRemote : Rank 0 will open/write/close the file using a remote IO service
!  - All tasks will return with a success flag
!  - No action taken on ranks that are not 0
! ioReadReplicate : reads will broadcast to all tasks
!  - cannot be or'd with ioAllLocal
!
! ****** A unit number may not be open remotely and locally ****
! ****** Remote IO from all tasks is not supported (yet)    ****

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters for manipulating IO
INTEGER, PARAMETER :: ioNoLocation=-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters for controlling how 'faulty' IO operations are dealt with !
INTEGER, PARAMETER :: ioPolicyStrict=27
INTEGER, PARAMETER :: ioPolicyConsistent=39
INTEGER, PARAMETER :: ioPolicySloppy=4123
INTEGER, PARAMETER :: ioPolicyLax=976
INTEGER, PARAMETER :: ioPolicyDefault=ioPolicyStrict

!Notes:
! The parameter is interpreted by file_open as
!
! ioPolicyStrict     : It is required that the unit is closed
! ioPolicyConsistent : If the file is already open on this unit,
!                      with the right permissions, position will
!                      be reset to zero
! ioPolicySloppy     : If the file is already open on this unit,
!                      the file position will be reset to 0
! ioPolicyLax      : The file will be opened as requested and
!                      any existing file open on that unit will
!                      be closed

! General parameters
CHARACTER (LEN=*),                                                           &
    PARAMETER       :: ioNoFileName ='* No Filename available *'
CHARACTER (LEN=*),                                                           &
    PARAMETER       :: ioNoEnvName  ='* No Environment available *'


END MODULE io_constants
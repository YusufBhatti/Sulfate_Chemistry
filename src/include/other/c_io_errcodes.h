/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#if !defined(C_IO_ERRCODES_H)
#define C_IO_ERRCODES_H

/* errors */

/* file was closed, but expected to be open */
#define ERR_FILE_CLOSED -1

/* file was open, but expected to be closed */
#define ERR_FILE_OPEN   -2

/* We unexpectedly got to the end of the file */
#define ERR_EOF         -3

/* We aren't allowed to do that for some reason */
#define ERR_ACCESS      -4

/* A request fails because of alignment issue */
#define ERR_MISSALIGN   -5

#endif

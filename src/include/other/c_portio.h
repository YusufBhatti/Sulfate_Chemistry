#if !defined(C_PORTIO_H) && defined(C95_2A)
#define C_PORTIO_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* This Header is for Portio2A only */

/* Description:
 *
 * Global definitions needed for Portable I/O
 *
 * Information:
 *
 * Provides:
 *   1. Inclusion of generic Portio API Header
 *   2. MAX_UNITS definition
 *   3. 'printstatus' levels definitions from their Fortran equivalent.
 */

#include "portio_api.h"

#define MAX_UNITS 300

#define PRINTSTATUS_MIN    1
#define PRINTSTATUS_NORMAL 2
#define PRINTSTATUS_OPER   3
#define PRINTSTATUS_DIAG   4

#endif

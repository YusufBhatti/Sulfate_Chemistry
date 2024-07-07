/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#if !defined(READ_WGDOS_HEADER_H)
#define READ_WGDOS_HEADER_H

#include <inttypes.h>

extern int64_t read_wgdos_header(char    *bytes_in,
                                 int64_t *accuracy,
                                 int64_t *cols,
                                 int64_t *rows);

#endif

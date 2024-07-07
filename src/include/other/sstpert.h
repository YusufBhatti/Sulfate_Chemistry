/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#if !defined(SSTPERT_H)
#define SSTPERT_H

#include <inttypes.h>

extern void sstpert(
                    double *factor,
                    int64_t *dt,
                    int64_t *num_rows,
                    int64_t *num_cols,
                    double *fieldclim,
                    double *fieldout);

#endif

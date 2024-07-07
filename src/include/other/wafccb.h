/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#if !defined(WAFCCB_H)
#define WAFCCB_H

#include <inttypes.h>

extern void convact(
                    int64_t *cols,
                    int64_t *rows,
                    int64_t *levels,
                    double *cpnrt_fld,
                    double *blkcld_flds,
                    double *concld_flds,
                    double *ptheta_flds,
                    double *rmdi,
                    bool   *return_heights,
                    double *p_cbb_out,
                    double *p_cbt_out,
                    double *cbhore_out
                    );

#endif

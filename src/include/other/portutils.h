#if !defined(PORTUTILS_H)
#define PORTUTILS_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* DEPENDS ON: portutils.o */

#include "c_fort2c_types.h"

/******************************************************************************/
/* Prototypes                                                                 */
/******************************************************************************/

/* Name mangling defines */

#if defined(C_LOW)

#define date_time date_time
#define word_length word_length
#define um_memcpy_f um_memcpy_f
#define shell shell
#define um_memcpy um_memcpy
#define file_op file_op

#elif defined(C_LOW_U)

#define date_time date_time_
#define word_length word_length_
#define um_memcpy_f um_memcpy_f_
#define shell shell_
#define um_memcpy um_memcpy_
#define file_op file_op_

#else

#define date_time DATE_TIME
#define word_length WORD_LENGTH
#define um_memcpy_f UM_MEMCPY_F
#define shell SHELL
#define um_memcpy UM_MEMCPY
#define file_op FILE_OP

#endif

/* actual prototypes */

extern void    date_time     (integer *,
                              integer *,
                              integer *,
                              integer *,
                              integer *,
                              integer *);

extern void    word_length   (integer *);

extern void    get_file      (integer *,
                              char [],
                              integer *,
                              integer *);

extern void    um_memcpy64   (void *,
                              void *,
                              int64_t);

extern void    um_sleep      (integer *);

extern void    um_memcpy_f   (char *,
                              char *,
                              integer *);

extern void    shell         (char *,
                              integer *);

extern void    um_memcpy     (char *,
                              char *,
                              integer *);

extern integer file_op       (char *,
                              integer *);

extern void    um_memcpy32   (void *,
                              void *,
                              int64_t);

extern int     um_str_split  (char *,
                              char,
                              char ***);

extern integer file_op_cp    (char *,
                              char *);

extern integer file_op_touch (char *);

extern integer file_op_rm    (char *,
                              char *);

extern integer file_op_mkdir (char *);

#endif

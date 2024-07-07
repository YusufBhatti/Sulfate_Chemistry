#if !defined(CIO_LIBC_H)
#define CIO_LIBC_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* DEPENDS ON: c_io_libc.o */

/* layer protypes */
extern int64_t c_io_libc_open        (int64_t, char *, int64_t, int64_t);
extern int64_t c_io_libc_setpos      (int64_t, int64_t, int64_t);
extern int64_t c_io_libc_close       (int64_t);
extern int64_t c_io_libc_in          (int64_t, char *, int64_t, int64_t);
extern int64_t c_io_libc_out         (int64_t, char *, int64_t, int64_t);
extern int64_t c_io_libc_change_mode (int64_t, int64_t);
extern int64_t c_io_libc_sync        (int64_t);
extern int64_t c_io_libc_query       (int64_t, query_type_t, void *, void *);

#endif

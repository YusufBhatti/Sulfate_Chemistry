#if !defined(CIO_LUSTRE_API_H)
#define CIO_LUSTRE_API_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* DEPENDS ON: c_io_lustreapi.o */

/* layer protypes */

extern int64_t c_io_lustreapi_init        (c_io_init_t *);

extern int64_t c_io_lustreapi_fini        (void);

extern int64_t c_io_lustreapi_in          (int64_t,
                                           char *,
                                           int64_t,
                                           int64_t);

extern int64_t c_io_lustreapi_out         (int64_t,
                                           char *,
                                           int64_t,
                                           int64_t);

extern int64_t c_io_lustreapi_setpos      (int64_t,
                                           int64_t,
                                           int64_t);

extern int64_t c_io_lustreapi_getpos      (int64_t);

extern int64_t c_io_lustreapi_open        (int64_t,
                                           char *,
                                           int64_t,
                                           int64_t);

extern int64_t c_io_lustreapi_close       (int64_t);

extern int64_t c_io_lustreapi_change_mode (int64_t,
                                           int64_t);

extern int64_t c_io_lustreapi_sync        (int64_t);

extern int64_t c_io_lustreapi_query       (int64_t,
                                           query_type_t,
                                           void *,
                                           void *);

#endif

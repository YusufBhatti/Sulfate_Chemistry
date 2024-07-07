#if !defined(C_IO_BYTESWAP_H)
#define C_IO_BYTESWAP_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* c_io_byteswap.h */

/* DEPENDS ON: c_io_byteswap.o */

/* layer protypes */

extern int64_t c_io_byteswap_init        (c_io_init_t *);

extern int64_t c_io_byteswap_fini        (void);

extern int64_t c_io_byteswap_in          (int64_t,
                                          char *,
                                          int64_t,
                                          int64_t);

extern int64_t c_io_byteswap_out         (int64_t,
                                          char *,
                                          int64_t,
                                          int64_t);

extern int64_t c_io_byteswap_setpos      (int64_t,
                                          int64_t,
                                          int64_t);

extern int64_t c_io_byteswap_getpos      (int64_t);

extern int64_t c_io_byteswap_open        (int64_t,
                                          char *,
                                          int64_t,
                                          int64_t);

extern int64_t c_io_byteswap_close       (int64_t);

extern int64_t c_io_byteswap_change_mode (int64_t,
                                          int64_t);

extern int64_t c_io_byteswap_sync        (int64_t);

extern int64_t c_io_byteswap_query       (int64_t,
                                          query_type_t,
                                          void *,
                                          void *);

#endif

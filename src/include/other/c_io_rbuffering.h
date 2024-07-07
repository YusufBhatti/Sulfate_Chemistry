#if !defined(CIO_RBUFFERING_H)
#define CIO_RBUFFERING_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* DEPENDS ON: c_io_rbuffering.o */

/* layer protypes */

extern int64_t c_io_rbuffering_init        (c_io_init_t *);

extern int64_t c_io_rbuffering_fini        (void);

extern int64_t c_io_rbuffering_in          (int64_t,
                                            char *,
                                            int64_t,
                                            int64_t);

extern int64_t c_io_rbuffering_out         (int64_t,
                                            char *,
                                            int64_t,
                                            int64_t);

extern int64_t c_io_rbuffering_setpos      (int64_t,
                                            int64_t,
                                            int64_t);

extern int64_t c_io_rbuffering_getpos      (int64_t);

extern int64_t c_io_rbuffering_open        (int64_t,
                                            char *,
                                            int64_t,
                                            int64_t);

extern int64_t c_io_rbuffering_close       (int64_t);

extern int64_t c_io_rbuffering_change_mode (int64_t,
                                            int64_t);

extern int64_t c_io_rbuffering_sync        (int64_t);

extern int64_t c_io_rbuffering_query       (int64_t,
                                            query_type_t,
                                            void *,
                                            void *);

#endif

#if !defined(C_PORTIO_API_H)
#define C_PORTIO_API_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* Prototypes for the Portio API */

#if defined(C95_2A)
/* DEPENDS ON: portio2a.o */
#elif defined(C95_2B)
/* DEPENDS ON: portio2b.o */
#endif

#include "c_fort2c_types.h"

/* Portio Control */

extern void portioinit(int64_t *,
                       int64_t *,
                       int64_t *,
                       bool *,
                       int64_t *,
                       bool *);

extern void portioshutdown(void);

/* Helper Threads */

extern void portioaddhelper(integer *);

extern void portiodetachallhelpers(void);

/* Getpos */

extern void getpos  (int64_t *, int64_t *);
extern void getpos8 (int64_t *, int64_t *);
extern void getpos32(int64_t *, int64_t *);
extern void getpos64(int64_t *, int64_t *);

/* Setpos */

extern void setpos_single  (int64_t *, int64_t *, int64_t *);
extern void setpos8_single (int64_t *, int64_t *, int64_t *);
extern void setpos32_single(int64_t *, int64_t *, int64_t *);
extern void setpos64_single(int64_t *, int64_t *, int64_t *);

/* buffin */

extern void buffin8_single(integer *,
                           void *,
                           integer *,
                           integer *,
                           real *);

extern void buffin16_single(integer *,
                            void *,
                            integer *,
                            integer *,
                            real *);

extern void buffin32_single(integer *,
                            void *,
                            integer *,
                            integer *,
                            real *);

extern void buffin64_single(integer *,
                            void *,
                            integer *,
                            integer *,
                            real *);

/* buffout */

extern void buffout8_single(integer *,
                            void *,
                            integer *,
                            integer *,
                            real *);

extern void buffout32_single(integer *,
                             void *,
                             integer *,
                             integer *,
                             real *);

extern void buffout64_single(integer *,
                             void *,
                             integer *,
                             integer *,
                             real *);

/* open/close */

extern void close_single(int64_t *,
                         void *,
                         int64_t *,
                         int64_t *,
                         int64_t *,
                         int64_t *);

extern void open_single(int64_t *,
                        void *,
                        int64_t *,
                        int64_t *,
                        int64_t *,
                        int64_t *);

/* misc */

extern void is_unit_open_local(int64_t *, bool *);

extern void sync_single(integer *, integer *);

extern void getextent(integer *, integer *, integer *);

extern void set_printstatus(integer *);

#endif

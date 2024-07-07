#if !defined(CIO_TIMING_H)
#define CIO_TIMING_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* DEPENDS ON: c_io_timing.o */

/* layer protypes */
extern int64_t c_io_timing_init        (c_io_init_t *);
extern int64_t c_io_timing_fini        (void);
extern int64_t c_io_timing_in          (int64_t, char *, int64_t, int64_t);
extern int64_t c_io_timing_out         (int64_t, char *, int64_t, int64_t);
extern int64_t c_io_timing_setpos      (int64_t, int64_t,int64_t);
extern int64_t c_io_timing_getpos      (int64_t);
extern int64_t c_io_timing_open        (int64_t, char *, int64_t, int64_t);
extern int64_t c_io_timing_close       (int64_t);
extern int64_t c_io_timing_change_mode (int64_t, int64_t);
extern int64_t c_io_timing_sync        (int64_t);
extern int64_t c_io_timing_query       (int64_t, query_type_t, void *, void *);

/* extra timer prototypes */
extern void c_io_start_extra_timer (int64_t);
extern void c_io_stop_extra_timer  (int64_t);

#endif

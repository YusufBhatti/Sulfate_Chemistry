#if !defined(THREAD_UTILS_INC)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#define THREAD_UTILS_INC

/* DEPENDS ON : thread_utils.o */

int64_t threadFlush         (void);

int64_t newLock             (void);

int64_t releaseLock         (int64_t *);

int64_t Lock                (int64_t *);

int64_t TestLock            (int64_t *);

int64_t unLock              (int64_t *);

int64_t threadID            (void);

int64_t inPar               (void);

int64_t numThreads          (void);

void    startOMPparallel    (void **,
                             void (*)(void **const));

void    startOMPparallelfor (void **,
                             void (*)(void **const,
                                      const int64_t *const restrict,
                                      const int64_t *const restrict,
                                      const int64_t *const restrict),
                             const int64_t *,
                             const int64_t *,
                             const int64_t *);

#endif

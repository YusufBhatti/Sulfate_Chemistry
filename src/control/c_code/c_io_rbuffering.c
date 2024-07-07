#if defined(C95_2B)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#if !defined(_BSD_SOURCE)
#define _BSD_SOURCE
#endif

#include "c_io_layers.h"

#define CIO_THIS_LAYER CIO_RBUFFERING

#include "c_io_nextlayer.h"

#include "c_io_internal.h"

#include "c_fort2c_prototypes.h"
#include "c_fort2c_types.h"
#include "c_io_rbuffering.h"
#include "c_io.h"

#include "thread_utils.h"

#include <inttypes.h>
#include <unistd.h>

/* Functions in this file:
 *
 *   Layer interfaces:
 *     init
 *     fini
 *     open
 *     close
 *     in
 *     out
 *     getpos
 *     setpos
 *     change_mode
 *     sync
 *     query
 *
 *   Implentation routines:
 *     c_io_rbuffering_waitCondition
 *     c_io_rbuffering_waitGECondition
 *     c_io_rbuffering_waitNotCondition
 *     c_io_rbuffering_check_latch
 *     c_io_rbuffering_helper
 *     c_io_rbuffering_helper_loop
 *     c_io_rbuffering_release_helper
 *     c_io_rbuffering_reset_unit
 *     c_io_rbuffering_reset_block
 *     c_io_rbuffering_init_block
 *     c_io_rbuffering_destroy_block
 *     c_io_rbuffering_update_block
 *     c_io_rbuffering_discard_block
 *     c_io_rbuffering_dumpstate
 *     c_io_rbuffering_get_block_for_pos
 *     c_io_rbuffering_instantiate_block_for
 *     c_io_rbuffering_load_block
 *     c_io_rbuffering_helper_pause
 *
 *   To make a new implementation of an IO layer, copy this
 *   file, replacing rbuffering with your identifier.
 */

/*
 * Helper threads. There may be a helper thread associated with this module.
 * The thread safety/locking metaphor for this module is as follows
 *
 * 1. The helper thread will ONLY load data into blocks marked as
 *    c_io_instantiated taking their state to c_io_ready (or c_io_invalid if
 *    the load was not possible - e.g. off the end of the file), it should not
 *    make any other modification to structure. Whilst loading it must hold the
 *    lock for that block
 *
 * 2. The primary thread will assign blocks to disk regions and will not
 *    load any data if the helper is present. The primary can alter any
 *    structure without the lock, unless the structure is marked as
 *    c_io_instantiated, in which case the lock must be taken to avoid
 *    concurrent access. The only time that the primary would need to modify
 *    such a block, is in the case where a new instantiation is needed but there
 *    are no non-instantiated buffers present.
 *
 *    RULES
 *     helper:  can only modify a block if state is "instantiated"
 *              can only change the state to ready
 *
 *     primary: can not modify the block if the state is instantiated
 *
 *     all:     must threadFlush() if altering state.
 *
 * Additionally, it should be noted that disk access for reads should only
 * happen when we are in the "reading" mode. This means blocks will not be set
 * to c_io_instantiated until a c_io_rbuffering_change_mode() call is made
 * to put us into reading mode.
 */

/* OPTIONAL PREPROCESSOR MACROS:
 *
 * - C_IO_RBUFF_VOLDATAPTR
 *   Define this macro at preprocessor time to enable the use of a volatile
 *   pointer for access to the data. This is used to prevent the compiler from
 *   eliminating memory reads or caching data when asynchronous updates are made
 *   by the helper thread; and may prevent some unsafe optimisations at higher
 *   optimisation levels.
 *   However, this comes at a performance cost - mainly due to it prohibiting
 *   the use of memcpy(). For this reason it is not enabled by default.
 *   Enabling this may be of use if you are debugging overly agressive compiler
 *   optimisations, or you are seeing data race conditions.
 */

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

/*----------------------------------------------------------------------------*/
/* Enumerators and Typedefs                                                   */
/*----------------------------------------------------------------------------*/

/* Enumerators */

/* blkStates is an enumerator for the state of a read bufffer block           */
typedef enum BLKSTATES {
                  /* --- unready states --- */

                  c_io_unallocated,  /* The block is not allocated, and must be
                                      * so before it can be used
                                      */

                  c_io_unassigned,   /* The block is allocated - and so is in
                                      * principle ready for use - but has not
                                      * yet been assigned a file position to
                                      * load from
                                      */

                  c_io_instantiated, /* The block has be set up ready to read
                                      * data, and has been given a file position
                                      * to load from, but has not yet done so
                                      */

                  /* --- ready states ---- */

                  c_io_ready,        /* The block contains the data loaded from
                                      * file, and can be returned by buffin call
                                      */

                  c_io_short,        /* The block contains the data loaded from
                                      * file, but the block is short as it
                                      * reached end-of-file
                                      */

                  c_io_consumed,     /* The block is c_io_ready, but all the
                                      * data in it has now been used (this
                                      * differs from c_io_ready only in its
                                      * debug reporting)
                                      */

                  /* --- error states --- */

                  c_io_invalid,      /* The block contains corrupted, failed, or
                                      * erroneously loaded data
                                      */

                  /* --- states counter --- */

                  numBlkStates

                } blkStates;

typedef enum WAIT_TYPE {
  blockstates,
  threadstates,
  boolstates,
} wait_type;

/* Structs */

typedef struct c_io_rblock_t {
  volatile int64_t start;
  volatile int64_t length;
  volatile blkStates state;
  volatile int64_t age;
#if defined(C_IO_RBUFF_VOLDATAPTR)
  volatile char    *vol_data_ptr;
#endif
  void             *data;
} c_io_rblock_t;

typedef struct c_io_rbuffering_counters_t {
  /* block counters updated by all threads.
   * These don't need to be atomic, as only one thread accesses them at a time.
   */
  volatile int64_t loaded;
  volatile int64_t invalid_blk;
  volatile int64_t b_short;

  /* block counters updated only by primary thread */
  int64_t alloced_imm;
  int64_t alloced_pre;
  int64_t b_used;

  /* wait counters updated by primary thread only */
  int64_t wait_calls;
  double  wait;
} c_io_rbuffering_counters_t;

typedef struct c_io_rbuffering_data_t
{
  int64_t       blockCurrent;

  int64_t       blockPosition;

  volatile bool has_instantiated; /* indicates if there is a file open on this
                                   * unit which has had a block "instantiated".
                                   * This in turn means the helper thread has
                                   * previously buffered for this unit.
                                   */

  c_io_rblock_t *volatile blocks;

  int64_t       cached_position;  /* caches the position returned by getpos()
                                   * when c_io_file_end is used. This is useful
                                   * to obtain correct file size measurements
                                   * in getextent() calls.
                                   */

  volatile bool helperpause;      /* Control switch to pause helper thread
                                   * buffering Set to true to request the helper
                                   * thread to be paused.
                                   * Set to false to unpause.
                                   */

  volatile bool helper_is_paused; /* Details the current state of the helper
                                   * thread.
                                   * True indicates the helper thread is paused.
                                   */
} c_io_rbuffering_data_t;

/* convenience typedefs */
typedef c_io_rbuffering_data_t *const restrict       rbuffdata_cr_ptr;
typedef const c_io_rbuffering_data_t *const restrict con_rbuffdata_cr_ptr;
typedef volatile const void *const restrict          vc_void_cr_ptr;
typedef const wait_type *restrict const              con_wt_cr_ptr;
typedef volatile void *const                         vol_void_con_ptr;
typedef c_io_rblock_t *const                         rblock_con_ptr;

/*----------------------------------------------------------------------------*/
/* Global Variables                                                           */
/*----------------------------------------------------------------------------*/

/* Pointers shared by all threads */

static c_io_layer_t               *volatile rbufLayer                = NULL;
static c_io_rbuffering_counters_t *volatile c_io_rbuffering_counters = NULL;
static c_io_rbuffering_data_t     *volatile c_io_rbuffering_data     = NULL;

/* variables shared by all threads */

static volatile size_t  numRBuffers = 0;

static volatile int64_t readthreadlock; /* Thread lock - only used is in
                                         * c_io_rbuffering_helper to serialise
                                         * the intialisation of helper threads,
                                         * and thus ensure there can only ever
                                         * be one helper.
                                         */

static volatile int64_t immediate_unit  = -1;
static volatile int64_t immediate_block = -1;

/* variables only used by primary thread */

static size_t  prefetchRBuffers = 0;
static size_t  RBufferSz        = 1024*1024; /* 1 MB */

static int64_t bufCount=0; /* global counter to age buffers */

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

void                  c_io_rbuffering_waitCondition         (vol_void_con_ptr,
                                                             int64_t,
                                                             wait_type);

void                  c_io_rbuffering_waitGECondition       (vol_void_con_ptr,
                                                             int64_t,
                                                             wait_type);

void                  c_io_rbuffering_waitNotCondition      (vol_void_con_ptr,
                                                             int64_t,
                                                             wait_type);

void                  c_io_rbuffering_check_latch           (void);

void                  c_io_rbuffering_helper                (void);

void                  c_io_rbuffering_helper_loop           (void);

void                  c_io_rbuffering_release_helper        (void);

void                  c_io_rbuffering_reset_unit            (size_t);

void                  c_io_rbuffering_reset_block           (c_io_rblock_t *);

void                  c_io_rbuffering_init_block            (c_io_rblock_t *);

void                  c_io_rbuffering_destroy_block         (c_io_rblock_t *);

void                  c_io_rbuffering_update_block          (int64_t,
                                                             char *,
                                                             int64_t,
                                                             int64_t);

void                  c_io_rbuffering_discard_block         (int64_t,
                                                             int64_t,
                                                             int64_t);

void                  c_io_rbuffering_dumpstate             (void);

int64_t               c_io_rbuffering_get_block_for_pos     (int64_t,
                                                             int64_t);

int64_t               c_io_rbuffering_instantiate_block_for (int64_t,
                                                             int64_t,
                                                             bool);

int64_t               c_io_rbuffering_load_block            (int64_t,
                                                             rblock_con_ptr);

static inline int64_t c_io_rbuffering_bytesLeftInBlock      (int64_t);

static inline int64_t wait_type_to_int64                    (vc_void_cr_ptr,
                                                             con_wt_cr_ptr);

static inline int64_t c_io_rbuffering_blockStart            (int64_t);

bool                  c_io_rbuffering_helper_pause          (rbuffdata_cr_ptr);

/*----------------------------------------------------------------------------*/
/* Layer routine implemetations                                               */
/* Note: all these routines are always called by the main thread, even in the */
/*       presence of a helper thread.                                         */
/*----------------------------------------------------------------------------*/

int64_t c_io_rbuffering_init(c_io_init_t *init_params)
{
  size_t i;
  rbufLayer               = c_io_init_next_layer(CIO_NEXT_LAYER);
  rbufLayer->helper_start = c_io_rbuffering_helper;
  rbufLayer->helper_stop  = c_io_rbuffering_release_helper;

  if (init_params->rbuffering_buffer_size > 0)
  {
    RBufferSz             = (size_t)init_params->rbuffering_buffer_size;
  }

  if (init_params->rbuffering_buffer_count > 0)
  {
    numRBuffers           = (size_t)init_params->rbuffering_buffer_count;
  }
  else
  {
    numRBuffers           = 0; /* i.e. off */
  }

  if (init_params->rbuffering_buffer_prefetch > 0)
  {
    prefetchRBuffers = MIN((size_t)params->rbuffering_buffer_prefetch+1,
                           numRBuffers) - 1;
  }

  readthreadlock=newLock();

  if(numRBuffers!=0){
    c_io_printf
      (-1,
       "init: rbuffering: Buffer size=%zu bytes, "
       "Buffers=%zu, prefetch=%zu",
       RBufferSz, numRBuffers, prefetchRBuffers);

    if (init_params->rbuffering_buffer_update==1)
      c_io_printf(-1,"init: rbuffering: buffers will update on writes");

    /* allocate local data */
    c_io_rbuffering_data = malloc((size_t)(init_params->max_unit)
                                   * sizeof(*c_io_rbuffering_data));

    c_io_rbuffering_counters = malloc(sizeof(*c_io_rbuffering_counters));

    if (c_io_rbuffering_data==NULL){
      ereport("c_io_rbuffering_init", 1,
          "init: rbuffering: allocation of control structure failed");
    }

    if (c_io_rbuffering_counters==NULL)
    {
      ereport("c_io_rbuffering_init", 1,
              "init: rbuffering: allocation of counters structure failed");
    }

    for (i=0;i<init_params->max_unit;i++)
    {
      c_io_rbuffering_reset_unit(i);
    }

    c_io_rbuffering_counters->loaded=0;
    c_io_rbuffering_counters->alloced_pre=0;
    c_io_rbuffering_counters->alloced_imm=0;
    c_io_rbuffering_counters->invalid_blk=0;
    c_io_rbuffering_counters->b_used=0;
    c_io_rbuffering_counters->b_short=0;
    c_io_rbuffering_counters->wait_calls=0;
    c_io_rbuffering_counters->wait=0.0;

  }
  else
    c_io_printf(-1,"init: rbuffering: is switched off (numRBuffers=0)");
  return C_IO_INIT(init_params);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_rbuffering_fini(void)
{
  int64_t lockid=readthreadlock;

  if (rbufLayer->helper_state == c_io_running)
  {
    c_io_rbuffering_release_helper();
  }

  releaseLock(&lockid);

  if(numRBuffers!=0)
  {
    c_io_printf(-1," RBbuffering usage");
    c_io_printf(-1," ");
    c_io_printf(-1,"  Blocks allocated for prefetch        %" PRId64,
                c_io_rbuffering_counters->alloced_pre);
    c_io_printf(-1,"  Blocks allocated for immediate use   %" PRId64,
                c_io_rbuffering_counters->alloced_imm);
    c_io_printf(-1,"  Blocks Successfully loaded           %" PRId64,
                c_io_rbuffering_counters->loaded);
    c_io_printf(-1,"  of which short blocks                %" PRId64,
                c_io_rbuffering_counters->b_short);
    c_io_printf(-1,"  Blocks Corrupted on load             %" PRId64,
                c_io_rbuffering_counters->invalid_blk);
    c_io_printf(-1,"  Total data read (blocks eqv.)        %.2f",
                (double)(c_io_rbuffering_counters->b_used)/(double)(RBufferSz));
    c_io_printf(-1,"  Wait calls                           %" PRId64,
                c_io_rbuffering_counters->wait_calls);
    c_io_printf(-1,"  Wait time                            %.2f secs",
                c_io_rbuffering_counters->wait);
    c_io_printf(-1,"\n");
    free (c_io_rbuffering_data);
    free (c_io_rbuffering_counters);
  }
  return C_IO_FINI();
}

/*----------------------------------------------------------------------------*/

/* the layer API design guarantees we are in "reading" mode when calling
 * c_io_rbuffering_in(). This has the consequence that
 * c_io_rbuffering_change_mode(unit,reading) must have been called, and thus at
 * least one block must be instantiated and have started pre-buffering when
 * numRBuffers is set to greater than zero.
 */
int64_t c_io_rbuffering_in (int64_t unit,
                            char *array,
                            int64_t len,
                            int64_t word_len)
{
  char *current_pos_ptr = array;
  int64_t transferSize;
  int64_t amountRead=0;

  /*
   * fast return if buffers are not active...
   */
  if(numRBuffers==0)
  {
    return C_IO_IN(unit,array,len,word_len);
  }

  c_io_rbuffering_data_t *const buffData=&c_io_rbuffering_data[unit];

  /* invalidate any cached position */
  buffData->cached_position=-1;

  c_io_rbuffering_check_latch();

  /* start the immediate unit/block loop on the helper thread */
  if (buffData->blockCurrent != -1)
  {
    immediate_block = buffData->blockCurrent;
    immediate_unit = unit;
    threadFlush();
  }

  /* convert request to bytes */
  len=len*word_len;

  c_io_rbuffering_counters->b_used+=len;

  while (len>0)
  {
    transferSize=c_io_rbuffering_bytesLeftInBlock(unit);

    if (transferSize>0)
    {
      /* we have some data remaining to be read from the current block */

      rblock_con_ptr block_ptr=&buffData->blocks[buffData->blockCurrent];

      /* check if the data assigned to the block has not yet been (fully)
       * loaded - and depending on helper thread presence either: -
       *    i) load it ourselves; or
       *   ii) wait until the helper thread has
       */
      if (block_ptr->state==c_io_instantiated)
      {
        if (rbufLayer->helper_state==c_io_running)
        {
          /* wait for helper */
          c_io_rbuffering_waitGECondition(&block_ptr->state,
                                          c_io_ready, blockstates);
        }
        else
        {
          /* no helper, load it ourselves */
          c_io_rbuffering_load_block(unit, block_ptr);
        }
      }

      /* shrink the transfer size if we want less than the rest of the block */
      if (len<transferSize) transferSize=len;

      if(block_ptr->state!=c_io_invalid)
      {
#if defined(C_IO_RBUFF_VOLDATAPTR)
        size_t j;

        for(j=0;j<(size_t)transferSize;j++)
        {
          current_pos_ptr[j] =
                     block_ptr->vol_data_ptr[(size_t)buffData->blockPosition+j];
        }
#else
        memcpy(current_pos_ptr,
               &((char *)block_ptr->data)[buffData->blockPosition],
               (size_t)transferSize);
#endif

        /* adjust pointers */
        amountRead+=transferSize;
        current_pos_ptr+=transferSize;
        buffData->blockPosition+=transferSize;
        len-=transferSize;

        /* check whether we used the whole block */
        if (buffData->blockPosition==block_ptr->length)
        {
          block_ptr->state=c_io_consumed;
          threadFlush();
        }
      }
      else
      {
        /* we have an invalid block - do a direct (non-buffered) read */

        int64_t file_pos, retVal;

        c_io_printf(unit, "rbuff: INFO:"
                    " loading directly to compensate for invalid block");

        if (rbufLayer->helper_state==c_io_running)
        {
          /* We are going to move the file position (setpos) on this unit,
           * so we need to have the helper thread paused
           */
          buffData->helperpause=true;
          /* Knock the helper thread out of the immediate unit/block loop */
          immediate_unit = -1;
          threadFlush();

          c_io_rbuffering_waitCondition(&buffData->helper_is_paused,
                                        true, boolstates);
        }

        /* we need to reset the disk position (the buffering may have moved it
         * elsewhere)
         */
        file_pos = buffData->blocks[buffData->blockCurrent].start
                    + buffData->blockPosition;

        retVal = C_IO_SETPOS(unit,file_pos,WORD8BYTES);

        if(retVal!=c_io_success)
        {
          c_io_printf(unit, "rbuff: PANIC: could not reset disk position!");
          ereport("c_io_rbuffering_in", 1,
                  "We need to reset disk position to cope with an invalid "
                  "block, but the real disk position could not be made "
                  "consistent with the buffered position");
        }

        /* we must transfer in bytes to preserve endianism
         * (which is handled in a higher layer)
         */
        retVal =  C_IO_IN(unit,current_pos_ptr,transferSize,WORD8BYTES);

        if (retVal!=transferSize)
        {
          /* error! */
          c_io_printf(unit, "rbuff: ERROR: could not load data: %"
                      PRId64 " bytes loaded out of %" PRId64 " requested",
                      retVal, transferSize);
          ereport("c_io_rbuffering_in", 1,
                  "c_io_rbuffering_in() has failed to load data either from "
                  "a prefetch block or directly from disk");
        }

        if (rbufLayer->helper_state==c_io_running)
        {
          /* unpause */
          buffData->helperpause=false;

          /* restart the immediate block loop */
          immediate_block = buffData->blockCurrent;
          immediate_unit = unit;

          threadFlush();
        }

        /* adjust pointers as if we had used the buffered block, in order to
         * maintain consistency
         */
        amountRead+=transferSize;
        current_pos_ptr+=transferSize;
        buffData->blockPosition+=transferSize;
        len-=transferSize;

      }
    }
    else
    {
      /* no space in the current block so set up a new one, RBufferSz further
       * on from the current one
       */

      int64_t block_start=buffData->blocks[buffData->blockCurrent].start;
      size_t i;

      block_start+=(int64_t)RBufferSz;
      c_io_rbuffering_instantiate_block_for(unit, block_start, true);

      /* start the immediate unit/block loop on the helper thread */
      immediate_block = buffData->blockCurrent;
      immediate_unit = unit;
      threadFlush();

      /* and prefetch ahead in multples of RBufferSz */
      for (i=0;i<prefetchRBuffers;i++)
      {
        block_start+=(int64_t)RBufferSz;
        c_io_rbuffering_instantiate_block_for(unit, block_start, false);
      }
    }
  }

  /* end the immediate unit/block loop on the helper thread */
  immediate_unit = -1;
  threadFlush();

  return amountRead/word_len; /* converting back to words of the req. type */
}

/*----------------------------------------------------------------------------*/

/* the layer API design guarantees we are in "writing" mode when calling
 * c_io_rbuffering_out().
 */
int64_t c_io_rbuffering_out (int64_t unit,
                 char    *array,
                 int64_t len,
                 int64_t word_len)
{
  rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[unit];

  /* check and update/discard the buffered block as required */
  if (numRBuffers!=0)
  {
    /* invalidate any cached position */
    buffData->cached_position=-1;

    if (buffData->has_instantiated)
    {
      const int64_t start = C_IO_GETPOS(unit);

      if ( params->rbuffering_buffer_update==1 )
      {
        c_io_rbuffering_update_block(unit,array,start,len*word_len);
      }
      else
      {
        c_io_rbuffering_discard_block(unit,start,len*word_len);
      }
    }
  }

  return C_IO_OUT(unit,array,len,word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_rbuffering_setpos (int64_t unit,
                                int64_t pos,
                                int64_t word_len)
{
  /* quick return if no buffers */
  if (numRBuffers==0)
  {
    return C_IO_SETPOS(unit,pos,word_len);
  }

  /* work out if we need to instantiate new blocks */
  {
    rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[unit];

    /* invalidate any cached position */
    buffData->cached_position=-1;

    if (c_io_attributes[unit].mode==reading)
    {
      if (pos==c_io_file_end)
      {
        int64_t res;

        /* only need to pause the helper thread if it has previously buffered on
         * this unit.
         */
        if (rbufLayer->helper_state==c_io_running && buffData->has_instantiated)
        {
          /* The request is for a seek to the end-of-file.
           * This is likely for a getextent() call, so we must interupt the
           * helper thread buffering (there is nothing beyond the end-of-file to
           * buffer anyway!) and actually enact the setpos request.
           * Once this is done, we can re-enable any buffering so that prefetch
           * buffers are once again populated.
           */
          buffData->helperpause=true;
          threadFlush();

          c_io_rbuffering_waitCondition(&buffData->helper_is_paused,
                                        true, boolstates);
        }

        res=C_IO_SETPOS(unit,pos,word_len);

        /* cache the current end-of-file.
         * This will allow the buffering to move the file position without
         * invalidating file size measurements.
         */
        buffData->cached_position=C_IO_GETPOS(unit);

        if (rbufLayer->helper_state==c_io_running && buffData->has_instantiated)
        {
          /* unpause */
          buffData->helperpause=false;
          threadFlush();
        }

        return res;
      }
      else
      {
        int64_t block_start=pos*word_len;
        size_t i;

        c_io_rbuffering_instantiate_block_for(unit,block_start,true);

        /* and prefetch ahead in multples of RBufferSz */
        for (i=0;i<prefetchRBuffers;i++)
        {
          block_start+=(int64_t)RBufferSz;
          c_io_rbuffering_instantiate_block_for(unit, block_start, false);
        }
      }
    }
    else
    {
      return C_IO_SETPOS(unit,pos,word_len);
    }
  }

  return c_io_success;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_rbuffering_getpos (int64_t unit)
{
  /* quick return if there are no buffers */
  if(numRBuffers==0)
  {
    return C_IO_GETPOS(unit);
  }

  /* check for cached positions from helper */
  {
    con_rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[unit];

    if (buffData->cached_position!=-1)
    {
      return buffData->cached_position;
    }

    if (c_io_attributes[unit].mode==reading)
    {
      return buffData->blocks[buffData->blockCurrent].start
              + buffData->blockPosition;
    }
  }

  /* there is no cached position - use the real disk position */
  return C_IO_GETPOS(unit);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_rbuffering_open (int64_t unit,
                              char    *filename,
                              int64_t file_name_len,
                              int64_t fileMode)
{
  size_t i;

  c_io_rbuffering_check_latch();

  if(numRBuffers!=0)
  {
    rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[unit];

    buffData->blocks = malloc(numRBuffers * sizeof(*(buffData->blocks)));

    if (buffData->blocks==NULL)
    {
      c_io_printf(-1, "open: rbuffering: "
                  "allocation of control structure %zu bytes failed",
                  (sizeof(c_io_rblock_t)*numRBuffers));

      ereport("c_io_rbuffering_open", 1,
              "allocation of control structure failed");
    }

    for (i=0;i<numRBuffers;i++)
    {
      /* block init will leave the block in the c_io_unassigned state */
      c_io_rbuffering_init_block(&buffData->blocks[i]);
    }

    threadFlush();
  }

  return C_IO_OPEN(unit,filename,file_name_len,fileMode);

}

/*----------------------------------------------------------------------------*/

int64_t c_io_rbuffering_close (int64_t unit)
{
  /* kill off our data structures. By pausing we guarentee the helper will not
   * update the block as we remove it.
   */
  if(numRBuffers!=0)
  {
    size_t i;
    c_io_rbuffering_data_t *const buffData = &c_io_rbuffering_data[unit];
    rblock_con_ptr tmp_ptr = buffData->blocks;

    /* pause */
    if (rbufLayer->helper_state==c_io_running)
    {
      if (buffData->has_instantiated)
      {
        buffData->helperpause=true;
        threadFlush();

        c_io_rbuffering_waitCondition(&buffData->helper_is_paused,
                                      true, boolstates);

      }
    }

    c_io_rbuffering_reset_unit((size_t)unit);

    /* c_io_rbuffering_reset_unit() has now unpaused the unit by resetting the
     * value of (buffData->helperpause). The blocks pointer of buffData is now
     * dissociated.
     */

    for (i=0;i<numRBuffers;i++)
    {
      c_io_rbuffering_destroy_block(&tmp_ptr[i]);
    }

    free(tmp_ptr);

    threadFlush();
  }

  return C_IO_CLOSE(unit);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_rbuffering_change_mode (int64_t unit,int64_t newMode)
{
  rbuffdata_cr_ptr buffData=&c_io_rbuffering_data[unit];

  int64_t retVal;
  int64_t file_pos;
  const bool mode_is_reading=(c_io_attributes[unit].mode==reading);

  /* we must wait pause any in progress buffering.
   * If we don't the disk position may be moved in the iterim;
   * or the retained buffers may be corrupted.
   */

  if (mode_is_reading && rbufLayer->helper_state==c_io_running)
  {
    /* if (buffData->has_instantiated) is not true, we have never started to
     * buffer on this unit, and we therefore don't need to wait for the
     * buffering to pause.
     */
    if (buffData->has_instantiated)
    {
      /* pause */
      buffData->helperpause=true;
      threadFlush();

      c_io_rbuffering_waitCondition(&buffData->helper_is_paused,
                                    true, boolstates);
    }
  }

  retVal=C_IO_CHANGE_MODE(unit,newMode);

  if(numRBuffers!=0)
  {

    if(retVal!=c_io_success)
    {
      /* The mode change went wrong,
       * so how can we know what to do with the buffered data?
       */
      ereport("c_io_rbuffering_change_mode", 1,
              "Failure to mode change in another layer "
              "means read buffer cannot be made consistent.");
    }

    if (newMode==reading)
    {
      size_t i;

      file_pos=C_IO_GETPOS(unit);

      /* unpause */
      if (rbufLayer->helper_state==c_io_running)
      {
         buffData->helperpause=false;
         threadFlush();
      }

      c_io_rbuffering_instantiate_block_for(unit,file_pos,true);

      /* and prefetch ahead in multples of RBufferSz */
      for (i=1;i<=prefetchRBuffers;i++)
      {
        int64_t start = file_pos + (int64_t)(i * RBufferSz);
        c_io_rbuffering_instantiate_block_for(unit, start, false);
      }
    }
    else if(mode_is_reading)
    {

      /*
       *   Restore the real position on disk, so that
       *   subsequent *write* operations land in the correct spot.
       */

      if (buffData->blockCurrent==-1)
      {
        c_io_printf(unit,
                    "rbuff: PANIC: changed mode away from reading"
                    " but there was no read block!");
        ereport("c_io_rbuffering_change_mode", 1,
                "changed mode away from reading but there was no read block!");
      }

      file_pos = buffData->blocks[buffData->blockCurrent].start
                 + buffData->blockPosition;

      buffData->cached_position = file_pos;

      retVal = C_IO_SETPOS(unit,file_pos,WORD8BYTES);

      if(retVal!=c_io_success)
      {
        c_io_printf(unit, "rbuff: PANIC: could not reset disk position!\n"
                    "file_pos: %" PRId64 ", blockCurrent: %" PRId64
                    ", blockPosition: %" PRId64,
                    file_pos,buffData->blockCurrent,buffData->blockPosition);
        ereport("c_io_rbuffering_change_mode", 1,
                "changed mode away from reading, but the real disk position "
                "could not be made consistent with the buffered position");
      }
    }
  }

  return retVal;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_rbuffering_sync(int64_t unit)
{
  /* Quick return if there are no buffers */
  if(numRBuffers==0)
  {
    return C_IO_SYNC(unit);
  }

  rbuffdata_cr_ptr buffData=&c_io_rbuffering_data[unit];

  size_t i;
  int64_t retVal;
  const bool mode_is_reading=(c_io_attributes[unit].mode==reading);

  /* we must wait pause any in progress buffering.
   * If we don't the disk position may be moved in the iterim;
   * or the reinitialised buffers may be corrupted.
   */

  if (mode_is_reading && rbufLayer->helper_state==c_io_running)
  {
    /* if (buffData->has_instantiated) is not true, we have never started to
     * buffer on this unit, and we therefore don't need to wait for the
     * buffering to pause.
     */
    if (buffData->has_instantiated)
    {
      /* pause */
      buffData->helperpause=true;
      threadFlush();

      c_io_rbuffering_waitCondition(&buffData->helper_is_paused,
                                    true, boolstates);
    }
  }

  /* Propogate the sync to the lower layers.
   * After this, we should consider already buffered reads to be "stale",
   * and so we will discard then.
   */
  retVal=C_IO_SYNC(unit);

  /* Discard the buffered reads.
   * If relevant, we will re-initialise them so they will actually be re-read
   * when the helper is unpaused. This ensures they are still correct in the
   * case where their data has changed on disk.
   */
  if (buffData->has_instantiated)
  {
    for (i=0;i<numRBuffers;i++)
    {
      if (    buffData->blocks[i].state == c_io_ready
           || buffData->blocks[i].state == c_io_short
           || buffData->blocks[i].state == c_io_consumed
         )
      {
        /* reset the block state */
        buffData->blocks[i].state = c_io_unassigned;

        /* Force this block to be the oldest, to ensure it is the one which is
         * instantiated next.
         */
        buffData->blocks[i].age = -1;

        /* If we are in reading mode, cause the blocks to actually be
         * re-instantiated.
         * This will mean they are re-read when the helper in unpaused.
         */
        if (mode_is_reading && rbufLayer->helper_state==c_io_running)
        {
          int64_t bytepos_tmp;
          bool current_tmp;

          /* save the key block values: blockCurrent and blockPosition*/

          current_tmp = (buffData->blockCurrent==(int64_t)i);

          /* Ensure that the byte position corresponds the current real disk
           * position. This ensures that blockPosition is re-calculated
           * correctly as the difference between bytepos and the block start
           * later on, during the re-instantiation.
           * (This is only actually relevant if blockCurrent==true.)
           */
          bytepos_tmp = buffData->blocks[i].start + buffData->blockPosition;

          /* Re-instantiate the block, thereby causing it to be reloaded */
          c_io_rbuffering_instantiate_block_for(unit, bytepos_tmp, current_tmp);
        }
      }
    }
  }

  /* unpause the helper */
  if (mode_is_reading && rbufLayer->helper_state==c_io_running)
  {
    buffData->helperpause=false;
    threadFlush();
  }

  return retVal;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_rbuffering_query(int64_t unit,
                              query_type_t query,
                              void *in_pack,
                              void *out_pack)
{
  return C_IO_QUERY(unit, query, in_pack, out_pack);
}

/*----------------------------------------------------------------------------*/
/* local routines below, should only be accessed by this file                 */
/* Note: these may be called by either the main thread, the helper thread, or */
/*       both, as indicated.                                                  */
/*----------------------------------------------------------------------------*/

/*
 * report remaining bytes in a block
 */
/* Called by: main thread only */
static inline int64_t c_io_rbuffering_bytesLeftInBlock(int64_t unit)
{
  con_rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[unit];

  int64_t res=0;

  if (buffData->blockCurrent != -1)
  {
    res = buffData->blocks[buffData->blockCurrent].length
          - buffData->blockPosition;
  }

  return res;
}

/*----------------------------------------------------------------------------*/

/*
 * starting address for a byte positions
 */
/* Called by: main thread only */
static inline int64_t c_io_rbuffering_blockStart(int64_t byte_pos)
{
  return (byte_pos/(int64_t)RBufferSz)*(int64_t)RBufferSz;
}

/*----------------------------------------------------------------------------*/

/*
 * reset the buffering data for a unit
 */
/* Called by: main thread only */
void c_io_rbuffering_reset_unit (size_t unit)
{
  rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[unit];
  buffData->blocks=NULL;
  buffData->blockCurrent=-1;
  buffData->blockPosition=-1;
  buffData->cached_position=-1;
  buffData->helperpause=false;
  buffData->helper_is_paused=false;
  buffData->has_instantiated=false;
  threadFlush();
}

/*----------------------------------------------------------------------------*/

/*
 * reset the buffering data for a block
 */
/* Called by: main thread or helper thread */
void c_io_rbuffering_reset_block (c_io_rblock_t *block)
{
  block->start=-1;
  block->length=-1;
  block->age=0;
  block->state=c_io_unassigned;
  threadFlush();
}

/*----------------------------------------------------------------------------*/

/*
 * allocate block structure
 */
/* Called by: main thread only */
void c_io_rbuffering_init_block (c_io_rblock_t *blk)
{
  blk->data = malloc(RBufferSz);
#if defined(C_IO_RBUFF_VOLDATAPTR)
  blk->vol_data_ptr=blk->data;
#endif

  if (blk->data == NULL)
  {
    c_io_printf(-1, "rbuff: init_block: "
                "Allocation for block 0x%" PRIuPTR " of %zu bytes failed",
                (uintptr_t)(blk), RBufferSz);

    ereport("c_io_rbuffering_init_block", 1, "Allocation for block failed ");
  }

  c_io_rbuffering_reset_block(blk);

  threadFlush();
}

/*----------------------------------------------------------------------------*/

/*
 * delete data structures assiciated with a unit
 */
/* Called by: main thread only */
void c_io_rbuffering_destroy_block (c_io_rblock_t *block)
{
  /* Clean up data structures. */
  c_io_rbuffering_reset_block(block);
  free(block->data);
  block->data=NULL;
#if defined(C_IO_RBUFF_VOLDATAPTR)
  block->vol_data_ptr=NULL;
#endif
  block->state=c_io_unallocated;
  threadFlush();
}

/*----------------------------------------------------------------------------*/

/*
 * assign a block for the given position.
 *
 * Note: that if we want to use this block next we must setCurrent=true,
 *       and byte_pos must be the actual postion we want to read from,
 *       not just 'somewhere' in the blocks range.
 */
/* Called by: main thread only */
int64_t c_io_rbuffering_instantiate_block_for (int64_t unit,
                                               int64_t byte_pos,
                                               bool    setCurrent)
{
  rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[unit];
  int64_t blk;
  size_t i;
  int64_t oldest;

  /* check whether the requested block exists */
  blk=c_io_rbuffering_get_block_for_pos(unit,byte_pos);

  /* create one if needed */
  if (blk==-1)
  {
    /* if we are creating a new block, the buffering now "has_instantiated" */
    buffData->has_instantiated=true;
    threadFlush();

    oldest=INT64_MAX;
    while (blk==-1)
    {
      for (i=0;i<numRBuffers;i++)
      {
        /* we can't use the "current" block, unless we are going to make the
         * new block current instead.
         */
        if (i==(size_t)buffData->blockCurrent && !setCurrent)
        {
          continue;
        }

        /*
         * We want the oldest buffer (ie the smallest 'age'), but we can't
         * choose an instantiated buffer if there is a helper running....
         */
        if (buffData->blocks[i].age<oldest &&
            (buffData->blocks[i].state != c_io_instantiated ||
             rbufLayer->helper_state   != c_io_running))
          {
            blk=(int64_t)i;
            oldest=buffData->blocks[i].age;
          }
      }

      if (blk==-1)
      {
        threadFlush();
        /* No block is currently available. This can only happen if the helper
         * is running, and all the blocks are in the instantiated state waiting
         * to be pre-buffered. We therefore need to sleep to allow the helper
         * time to finish loading some blocks.
         */
        usleep(c_io_sleep_waittime);
      }
    }

    /* DEBUG PRINT
      c_io_printf(unit,
                "rbuff: instantiate_block, curr= %d, new=%d @%p has age=%lld",
                buffData->blockCurrent,
                blk,
                &buffData->blocks[blk],
                LLD(buffData->blocks[blk].age));
    */
    buffData->blocks[blk].start=c_io_rbuffering_blockStart(byte_pos);
    buffData->blocks[blk].age=++bufCount;
    if (buffData->blocks[blk].start<0)
    {
      c_io_printf(unit,
                  "rbuff: negative block start generated for position %"
                  PRId64 " (ouch!)",
                  (byte_pos));
      ereport("c_io_rbuffering_instantiate_block_for", 1,
              "Negative block start generated");
    }
    buffData->blocks[blk].state=c_io_instantiated;
    buffData->blocks[blk].length=(int64_t)RBufferSz;

    if (setCurrent)
    {
      c_io_rbuffering_counters->alloced_imm++;
    }
    else
    {
      c_io_rbuffering_counters->alloced_pre++;
    }

    threadFlush();
  }

  /*  DEBUG PRINT
  else c_io_printf
         (unit,
          "rbuff: instantiate_block, already exists; "
          "setCurrent=%s blk=%d @ %p age=%lld",
          (setCurrent)?"true":"false",
          blk,
          &buffData->blocks[blk],
          LLD(buffData->blocks[blk].age));
  */

  if (setCurrent)
  {
    /* set the current block to this one */
    buffData->blockCurrent=blk;
    /* set the location in the block */
    buffData->blockPosition=byte_pos-buffData->blocks[blk].start;
  }
  return blk;
}

/*----------------------------------------------------------------------------*/

/*
 * find an existing block for the given postion
 */
/* Called by: main thread only */
int64_t c_io_rbuffering_get_block_for_pos (int64_t unit,
                                           int64_t byte_pos)
{
  rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[unit];
  size_t blk=0;
  int64_t start_pos;

  /* loop over blocks and see if the target is in one */
  start_pos=c_io_rbuffering_blockStart(byte_pos);
  for(blk=0; blk<numRBuffers; blk++)
  {
    if (buffData->blocks[blk].start==start_pos)
    {
      /* don't freshen invalid blocks...
       * we don't want to cause an attempt at re-reading to be triggered by not
       * returning the invalid block.
       * But we also don't want good blocks to be discarded because they
       * become older than an invalid block.
       */
      if(buffData->blocks[blk].state!=c_io_invalid)
      {
        /* freshen the block */
        buffData->blocks[blk].age=++bufCount;
      }

      return (int64_t)blk;
    }
  }

  return -1;
}

/*----------------------------------------------------------------------------*/

/*
 * discard blocks in a given range by resetting them to c_io_unassigned
 */
/* Called by: main thread only */
void c_io_rbuffering_discard_block (int64_t unit,
                                    int64_t start,
                                    int64_t extent)
{
  rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[unit];
  size_t i;

  for (i=0;i<numRBuffers;i++)
  {
    if (buffData->blocks[i].state>=c_io_ready)
    {
      if (
            /* the write started before the end of the full block... */
            (start < buffData->blocks[i].start+(int64_t)RBufferSz)

            /* ... and the write finished after the start of the block */
            && (start+extent >= buffData->blocks[i].start)
         )
      {
        /* i.e. we think the buffer has data and there is overlap */
        c_io_rbuffering_reset_block(&buffData->blocks[i]);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/

/*
 * update the data in cached blocks in a given range with the data from array
 */
/* Called by: main thread only */
void c_io_rbuffering_update_block (int64_t unit,
                                   char    *array,
                                   int64_t start,
                                   int64_t extent)
{
  rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[unit];
  size_t i;
  size_t bl_pos;
  size_t ar_pos;
  size_t copyAmount;

  for (i=0;i<numRBuffers;i++)
  {
    const int64_t bl_start = buffData->blocks[i].start;
    int64_t bl_length = buffData->blocks[i].length;

    /* deal with the special requirements of short blocks */
    if (buffData->blocks[i].state==c_io_short)
    {
      if (
            /* the write started before the end of the full block...
             * NB: because this is a short block, the write may be beyond the
             * block's length, but between blocks[i].length and RBufferSz
             */
            (start < bl_start+(int64_t)RBufferSz)

            /* ... and the write finished after the start of the block */
            && (start+extent >= bl_start)
         )
      {
        /* If the update starts before the end of the short block
         * (within blocks[i].start+buffData->blocks[i].length)
         * then we can extend the length of the block if required. If however
         * the update starts beyond the length of the block, we must discard
         * it, to avoid leaving an "undefined" section within the block
         * interior.
         */

        if(start < bl_start + bl_length)
        {
          /* Calculate the new block length, extending it to fit the new write.
           * Within this block, this will finish on the smaller of:
           *  * start+extent (the end of the update)
           *  * buffData->blocks[i].start+RBufferSz (the maximum block size)
           */
          bl_pos = (size_t)(MIN(start + extent, bl_start + (int64_t)RBufferSz)
                            - bl_start);

          /* if the block now has the maximum size, it is no longer "short" */
          if (bl_pos==RBufferSz)
          {
             buffData->blocks[i].state=c_io_ready;
          }

          buffData->blocks[i].length=(int64_t)bl_pos;
          bl_length = buffData->blocks[i].length;
          threadFlush();
        }
        else
        {
          /* There is non-contiguous data, leaving part of the block undefined,
           * so discard it
           */
          buffData->blocks[i].state=c_io_unassigned;
          buffData->blocks[i].age=0;
          threadFlush();
        }
      }
    }

    /* update blocks which have already been read */
    if (buffData->blocks[i].state==c_io_ready ||
        buffData->blocks[i].state==c_io_short ||
        buffData->blocks[i].state==c_io_consumed)
    {
      if (
            /* the write started before the end of the block... */
            (start < bl_start + bl_length)

            /* ... and the write finished after the start of the block */
            && (start+extent >= bl_start)
         )
      {
        /* i.e. we think the buffer has data and there is overlap */

        /* where in the block does the data start? */
        bl_pos=0;
        if (start>bl_start)
        {
          bl_pos=(size_t)(start-bl_start);
        }

        /* where in the data array does the block start? */
        ar_pos=0;
        if (bl_start>start)
        {
          ar_pos=(size_t)(bl_start-start);
        }

        /* how much should we copy? */
        copyAmount=(size_t)bl_length-bl_pos;
        if (start+extent<(bl_start+bl_length))
        {
          copyAmount = (size_t)(start + extent - (bl_start + (int64_t)bl_pos));
        }

        /* copy the array data into the block */
        memcpy(&((char *)buffData->blocks[i].data)[bl_pos],
               &array[ar_pos],
               copyAmount);

      }
    }
  }
}

/*----------------------------------------------------------------------------*/

/*
 * load a block from disk
 */
/* Called by: main thread or helper thread */
int64_t c_io_rbuffering_load_block(int64_t unit,
                                   rblock_con_ptr block)
{
  /* DEBUG PRINT
    c_io_printf(unit,
              "rbuff: load_block %p addr=%lld age=%lld",
              block,
              LLD(block->start),
              LLD(block->age));
  */

  if (C_IO_SETPOS(unit,block->start,WORD8BYTES) != c_io_fail)
  {
    int64_t ret;

    ret=C_IO_IN(unit,
               block->data,
               block->length,
               WORD8BYTES);

    if(ret==block->length)
    {
      block->state=c_io_ready;
      c_io_rbuffering_counters->loaded++;
      threadFlush();

      return c_io_success;
    }
    else if (c_io_attributes[unit].end_of_file)
    {
      if (ret > 0)
      {
        /* end-of-file - use a short block */
        block->state=c_io_short;
        block->length=ret;
        c_io_rbuffering_counters->loaded++;
        c_io_rbuffering_counters->b_short++;
      }
      else
      {
        /* The block is so short it's empty!
         * Reset it so it can go back into the pool for instantiation
         */
        c_io_rbuffering_reset_block(block);
      }

      threadFlush();

      return c_io_success;

    }
  }

  c_io_rbuffering_counters->invalid_blk++;
  block->state=c_io_invalid;
  threadFlush();

  return c_io_fail;
}

/*----------------------------------------------------------------------------*/

/*
 * debug aid, not normally called
 */
/* Called by: <neither by default> */
void c_io_rbuffering_dumpstate (void)
{
  size_t i,j;
  c_io_printf(-1,
           "dump: *********************************************************");
  for (i=(size_t)(params->min_unit);i<(size_t)(params->max_unit);i++){
    if (c_io_rbuffering_data[i].blockCurrent!=-1){
      c_io_printf((int64_t)i,"dump: blockCurrent     =%" PRId64,
                  c_io_rbuffering_data[i].blockCurrent);
      c_io_printf((int64_t)i,"dump: blockPosition =%" PRId64,
                  c_io_rbuffering_data[i].blockPosition);
      c_io_printf((int64_t)i,"dump: data address      =0x%" PRIuPTR,
                 (uintptr_t)(c_io_rbuffering_data[i].blocks));
      if (c_io_rbuffering_data[i].blocks!=NULL)
    for (j=0;j<numRBuffers;j++){
      c_io_printf((int64_t)i,"dump:    block=%zu start=%" PRId64,
              j,c_io_rbuffering_data[i].blocks[j].start);
      c_io_printf((int64_t)i,"dump:    block=%zu state=%ld",
              j,(long)c_io_rbuffering_data[i].blocks[j].state);
#if defined(C_IO_RBUFF_VOLDATAPTR)
      c_io_printf((int64_t)i,"dump:    block=%zu data =0x%" PRIuPTR,
              j,(uintptr_t)(c_io_rbuffering_data[i].blocks[j].vol_data_ptr));
#else
      c_io_printf((int64_t)i,"dump:    block=%zu data =0x%" PRIuPTR,
              j,(uintptr_t)(c_io_rbuffering_data[i].blocks[j].data));
#endif
    }
    }
  }
  c_io_printf(-1,
           "dump end: *****************************************************");
}

/*----------------------------------------------------------------------------*/

/* Called by: main thread only */
void c_io_rbuffering_check_latch(void)
{
  threadFlush();

  if (rbufLayer->helper_state==c_io_present)
  {
    /* There is a helper thread that has been attached, but is waiting in the
     * c_io_present state until a latch check is called.
     * We should now start that helper thread running.
     */
    c_io_printf(-1,"rbuffering: check_latch: Who's there?");
    rbufLayer->helper_state=c_io_attaching;
    threadFlush();

    while(rbufLayer->helper_state != c_io_running)
    {
      threadFlush();
    }
  }
}

/*----------------------------------------------------------------------------*/

/* Called by: helper thread only */
void c_io_rbuffering_helper(void)
{
  /* The logic of the rbuffering_helper works provided at most only one main
   * and one helper threads exist at any time.
   * Consequently, we need to protect this routine such that multiple threads
   * cannot simultaneously be created due to race conditions on
   * rbufLayer->helper_state
   */
  helperStates tmp_state;
  int64_t lockid=readthreadlock;

  Lock(&lockid);
  {
    /* all checks and updates of rbufLayer->helper_state will happen in the
     * locked region to prevent race conditions
     */

    threadFlush();

    tmp_state = rbufLayer->helper_state;

    /* Do we have a bad initial state (i.e state is not c_io_none)?
     * If no, initiate the handshake.
     * If yes, do nothing - this will be handled outside of the lock.
     */
    if (tmp_state == c_io_none)
    {
      c_io_printf(-1,"rbuffering: helper: Knock knock!");
      rbufLayer->helper_state=c_io_present;
    }

    threadFlush();
  }
  unLock(&lockid);

  /* handle any bad initial state from within the lock */
  if (tmp_state != c_io_none)
  {
    c_io_printf(-1, "rbuffering: helper:"
                " attaching has failed. (Bad initial condition,"
                " or another helper thread is already attached.)");
    return;
  }

  while (rbufLayer->helper_state != c_io_attaching)
  {
    threadFlush();
  }

  c_io_printf(-1,"rbuffering: helper: A helper thread.");
  rbufLayer->helper_state = c_io_running;
  threadFlush();

  c_io_rbuffering_helper_loop();
  c_io_printf(-1,"rbuffering: helper: I'm gone.");
}

/*----------------------------------------------------------------------------*/

/* Called by: main thread only */
void c_io_rbuffering_release_helper(void){

  /* In case another thread had requested to attach as a helper, but never made
   * it into the c_io_running state because no subsequent check latch call was
   * made - make a call to c_io_rbuffering_check_latch() now.
   */
  c_io_rbuffering_check_latch();

  if (rbufLayer->helper_state != c_io_running){
    c_io_printf(-1,"rbuffering: release_helper: nothing to release!");
    return;
  }

  c_io_printf(-1,"rbuffering: release_helper: Go away!");

  /*
   * Set the latch out variable
   */
  rbufLayer->helper_state=c_io_terminate;
  threadFlush();

  /*
   * Wait for the other thread to notice
   */
  c_io_rbuffering_waitNotCondition(&rbufLayer->helper_state,
                                   c_io_terminate,
                                   threadstates);

  c_io_printf(-1,"rbuffering: release_helper: Bye!");
}

/*----------------------------------------------------------------------------*/

/* Called by: helper thread only */
bool c_io_rbuffering_helper_pause(rbuffdata_cr_ptr buffData)
{
  bool res=false;

  threadFlush();

  if (buffData->helperpause)
  {
    /* mechanism to pause buffering for e.g. getextent() calls */
    buffData->helper_is_paused=true;
    res=true;
  }
  else
  {
    buffData->helper_is_paused=false;
  }

  threadFlush();

  return res;
}

/*----------------------------------------------------------------------------*/

/* Called by: helper thread only */
void c_io_rbuffering_helper_loop(void)
{
  threadFlush();

  if (numRBuffers==0)
  {
    c_io_printf(-1, "rbuffering: helper: no prefetch blocks are configured");
  }

  if(numRBuffers!=0)
  {
    size_t i,b;

    /* We want the first unit (value of i) checked to be zero.
     * This requires us to set the intial value of i to params->max_unit - 1
     * as each iteration of the loop increments i, then finds its value
     * modulo params->max_unit, before checking that unit.
     * ( (params->max_unit-1)+1 % params->max_unit = 0 )
     */
    i = params->max_unit - 1;

    /* similarly for the first block value */
    b = numRBuffers - 1;

    threadFlush();

    while (rbufLayer->helper_state != c_io_terminate)
    {
      int64_t imunit = immediate_unit;

      /* Always check the "immediate" unit/block first.
       *
       * The immediate unit (immediate_unit) is set by the main thread on
       * entering a c_io_rbuffering_in() call. For performance reasons, we
       * should loop over this unit's blocks until the read has finished.
       * This minimises the time the main thread has to wait in order for the
       * helper thread to have loaded the blocks needed for the read.
       *
       * At other times - when there is no immediate unit - we want to cycle
       * over all the units, to ensure we load blocks evenely in cases where we
       * have lots of open units. This causes the helper thread to try to
       * prebuffer ahead at least some of the way on every unit, whilst we are
       * not in a time critical situation (because no unit is curently being
       * used for reads). This in turn prevents the case where a large number
       * of prebuffer blocks queued up for one unit delays the loading of any
       * blocks on another unit until after the read has started on it - thus
       * missing an oportunity to save time.
       */
      if (imunit!=-1)
      {
        int64_t imblock = immediate_block;
        int64_t o_imblock = imblock;

        /* if we get here, it must be true that imunit has instantiated and is
         * not paused (or pausing).
         */
        rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[imunit];

        while(immediate_unit==imunit)
        {
          rblock_con_ptr imblock_ptr = &buffData->blocks[imblock];

          if (imblock_ptr->state==c_io_instantiated)
          {
            c_io_rbuffering_load_block((int64_t)imunit,imblock_ptr);

            /* reset o_imblock, so we do a complete cycle with no loads before
             * we explicitly threadFlush(). (note: there is an implied
             * threadFlush() from the c_io_rbuffering_load_block() call)
             */
            o_imblock=imblock;
          }

          imblock = (imblock+1) % (int64_t)numRBuffers;
          if(imblock==o_imblock)
          {
            threadFlush();
          }
        }
      }
      else
      {
        /* If there were no blocks loaded on the "immediate" unit, we can try
         * the next block in the cycle, to find non-immediate blocks to prefetch
         */

        /* cycle through the blocks */
        b = (b+1) % numRBuffers;

        /* if we get back to the zeroth block, cycle to the next unit */
        if (b==0)
        {
          /* cycle through the units */
          i = (i+1) % (params->max_unit);
        }

        {
          rbuffdata_cr_ptr buffData = &c_io_rbuffering_data[i];

          /* there can be no instantiated blocks if we have not ever made a
           * call to c_io_instantiate_block_for()
           * ie. buffData->has_instantiated is false
           */
          threadFlush();
          if (buffData->has_instantiated)
          {
            /* check for paused units */
            if (c_io_rbuffering_helper_pause(buffData))
            {
              /* Unit is paused: move on to the next unit.
               * Reset b to (numRBuffers-1), as this will cause unit to cycle
               * next time round.
               */
              b = numRBuffers-1;
            }
            else if (buffData->blocks[b].state==c_io_instantiated)
            {
              /*
              c_io_printf(i,
                          "rbuffering: helper loading block %d with age=%lld",
                          b,
                          LLD(buffData->blocks[b].age));
                          b,
                          LLD(buffData->blocks[b].age));
              */

              c_io_rbuffering_load_block((int64_t)i,&buffData->blocks[b]);
              /* We have loaded a block, so cycle to the next unit.
               * Reset b to (numRBuffers-1), as this will cause unit to cycle
               * next time round.
               */
              b = numRBuffers-1;
            }
          }
        }
      }

      threadFlush();
    }

  }
  else
  {
    while (rbufLayer->helper_state != c_io_terminate)
    {
      /* There is no buffering to be done (numRBuffers==0)
       * Spool until we are told to shut down.
       * (This avoids race conditions in negotiating helper thread state)
       */
      usleep(c_io_sleep_waittime);
      threadFlush();
    }
  }

  rbufLayer->helper_state = c_io_none;
  threadFlush();
}

/*----------------------------------------------------------------------------*/

/*
 *  wait for the contents of flag to be cond.
 */
/* Called by: main thread or helper thread */
void c_io_rbuffering_waitCondition(vol_void_con_ptr flag,
                                   int64_t cond,
                                   wait_type wt)
{
  double t;

  threadFlush();

  if (wait_type_to_int64(flag,&wt) == cond) return;

  t=c_io_getTime();

  while (wait_type_to_int64(flag,&wt) != cond)
  {
    usleep(c_io_sleep_waittime);
    threadFlush();
  }

  c_io_rbuffering_counters->wait+=(c_io_getTime()-t);
  c_io_rbuffering_counters->wait_calls++;
}

/*----------------------------------------------------------------------------*/

/*
 *  wait for the contents of flag to be greater or equal to cond.
 */
/* Called by: main thread or helper thread */
void c_io_rbuffering_waitGECondition(vol_void_con_ptr flag,
                                     int64_t cond,
                                     wait_type wt)
{
  double t;

  threadFlush();

  if (wait_type_to_int64(flag,&wt) >= cond) return;

  t=c_io_getTime();

  while (wait_type_to_int64(flag,&wt) < cond)
  {
    usleep(c_io_sleep_waittime);
    threadFlush();
  }

  c_io_rbuffering_counters->wait+=(c_io_getTime()-t);
  c_io_rbuffering_counters->wait_calls++;
}

/*----------------------------------------------------------------------------*/

/*
 *  wait for the contents of flag to be not equal to cond.
 */
/* Called by: main thread or helper thread */
void c_io_rbuffering_waitNotCondition(vol_void_con_ptr flag,
                                      int64_t cond,
                                      wait_type wt)
{
  double t;

  threadFlush();

  if (wait_type_to_int64(flag,&wt) != cond) return;

  t=c_io_getTime();

  while (wait_type_to_int64(flag,&wt) == cond)
  {
    usleep(c_io_sleep_waittime);
    threadFlush();
  }

  c_io_rbuffering_counters->wait+=(c_io_getTime()-t);
  c_io_rbuffering_counters->wait_calls++;
}

/*----------------------------------------------------------------------------*/

/* Called by: main thread or helper thread */
static inline int64_t wait_type_to_int64(vc_void_cr_ptr flag,
                                         con_wt_cr_ptr wt)
{
  int64_t ret=-1;

  switch(*wt)
  {
    case blockstates:
    {
      ret=(int64_t)*(volatile const blkStates *)flag;
    }
    break;

    case threadstates:
    {
      ret=(int64_t)*(volatile const helperStates *)flag;
    }
    break;

    case boolstates:
    {
      ret=(int64_t)*(volatile const bool *)flag;
    }
    break;
  }

  return ret;
}

#endif

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

#define CIO_THIS_LAYER CIO_WBUFFERING

#include "c_io_nextlayer.h"
#include "c_io_internal.h"

#include "c_fort2c_prototypes.h"
#include "c_fort2c_types.h"
#include "c_io_wbuffering.h"
#include "c_io.h"
#include "thread_utils.h"
#include <inttypes.h>
#include <unistd.h>

#define DEFAULT_WBUFFER_SIZE 1024*1024 /* 1 MB */
#define MAX_HELPER_QUEUE 15
#define MAX_HELPER_QUEUE_SIZE 20*DEFAULT_WBUFFER_SIZE

/*  Functions in this file:
 *
 *  Layer interfaces:
 *    init          : initialise data structures, and pass through
 *    fini          : pass through
 *    open          : reset per unit data structures, and pass through
 *    close         : flush outstanding data, pass through, and reset
 *                    data structures
 *    in            : pass through (no read buffering yet)
 *    out           : write to buffers or pass chunks through as needed
 *    getpos        : pass through, and append any buffering data
 *    setpos        : flush buffered data, and pass through
 *    change_mode   : flush buffered data, and pass through
 *    sync          : flush buffered data, and pass through
 *    query         : pass through
 *
 *  Implentation auxiliary routines:
 *    c_io_wbuffering_reset                 : reset per unit attributes
 *    c_io_wbuffering_flush_nocheck         : pass buffered data to the next
 *                                            layer
 *    c_io_wbuffering_flush_check           : as above, but first checking the
 *                                            status of the helper thread
 *    c_io_wbuffering_out_buffer_to_disk    : send data directly to disk
 *    c_io_wbuffering_out_via_helper        : send data to the helper thread
 *    c_io_wbuffering_helper                : run the helper thread
 *    c_io_wbuffering_helper_loop           : the helper thread's main work loop
 *    c_io_wbuffering_helper_process_buffer : send buffered data to disk
 *    c_io_wbuffering_release_helper        : stop the helper
 *    c_io_wbuffering_check_latch           : check the helper's status
 *    c_io_wbuffering_waitNotCondition      : wait for a condition
 *    c_io_wbuffering_waitLECondition       : wait for a condition
 *    wait_type_to_int64                    : inlined convenience conversion
 *                                            function
 */

/*----------------------------------------------------------------------------*/
/* Enumerators and Typedefs                                                   */
/*----------------------------------------------------------------------------*/

/* structs */

typedef struct c_io_wbuffering_data_t {
  size_t  bufferMax;       /* The size of buffer */
  size_t  bufferSize;      /* The number of bytes in buffer */
  size_t  bufferPosition;  /* The offset into buffer that the
                              next data is to be read from */
  size_t  bufferHits;
  size_t  bufferMiss;
  char    *buffer;          /* just some memory */
} c_io_wbuffering_data_t;

typedef struct c_io_wbuffering_counters_t {
  int64_t written;       /* actual disk write events [calls to C_IO_OUT()]    */
  int64_t write_calls;   /* requested writes [calls to c_io_wbuffering_out()] */
  int64_t wait_calls;    /* number of times we waited on a condition          */
  size_t hits;           /* count of actual writes prevented by buffering     */
  size_t misses;         /* count of actual writes caused by not buffering    */
  int64_t b_used;        /* total bytes written                               */
  double  wait;          /* total wait time                                   */
} c_io_wbuffering_counters_t;

typedef struct c_io_wbuffering_helperbuffers_t {
  /* c_io_wbuffering_helperbuffers_t is used to create a linked list of buffers
   * for use by the helper thread.
   * "next" is the next buffer object in the linked list
   */
  struct c_io_wbuffering_helperbuffers_t *next;

  char *buffer;                    /* buffered write data                     */
  size_t buffer_len;               /* length of buffered data                 */
  bool ready;                      /* true if buffer is ready for writing     */
  bool written;                    /* true if buffer has been written to disk */
} c_io_wbuffering_helperbuffers_t;

/* enumerators */

typedef enum WAIT_TYPE {
  wt_prioritystates,
  wt_threadstates,
  wt_size_t
} wait_type;

/* convenience typedefs */

typedef volatile const void *restrict const res_void_ptr_t;
typedef const wait_type *restrict const res_wt_ptr_t;
typedef c_io_wbuffering_data_t *const restrict c_io_wbuffering_data_t_cr_ptr;

/*----------------------------------------------------------------------------*/
/* Global Variables                                                           */
/*----------------------------------------------------------------------------*/

static c_io_layer_t               *wbufLayer;

static c_io_wbuffering_counters_t *c_io_wbuffering_counters = NULL;
static c_io_wbuffering_data_t     *c_io_wbuffering_data;

static size_t defaultBufferSz = DEFAULT_WBUFFER_SIZE;

/* linked list of buffers for the helper thread to write.
 * Once allocated, this always contains at least one buffer
 * (ie. it is never NULL).
 */
static c_io_wbuffering_helperbuffers_t **helperbuffers = NULL;

/* pointer to last buffer in the helperbuffers linked list */
static c_io_wbuffering_helperbuffers_t **lastbuffer    = NULL;

/* switch to indicate the helper should prioritise writing a particular unit
 * (because its buffers need to be flushed for e.g. a mode change or file close)
 */
static bool   priority_buffer = false;

/* the unit to be prioritised in the priority_buffer = true case */
static size_t priority_unit = 0;

static size_t pending_helper_writes=0;
static size_t pending_helper_writes_size=0;

/* lock to protect the global state variables (wbufLayer; priority_buffer)
 * from concurrent access by the helper and the main threads
 */
static int64_t writethreadlock;

static size_t unit_max;

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

static void           c_io_wbuffering_reset                 (size_t);

static int64_t        c_io_wbuffering_flush_check           (int64_t);

static int64_t        c_io_wbuffering_flush_nocheck         (int64_t);

static void           c_io_wbuffering_check_latch           (void);

static void           c_io_wbuffering_release_helper        (void);

static void           c_io_wbuffering_helper_loop           (void);

static void           c_io_wbuffering_helper                (void);

static void           c_io_wbuffering_waitNotCondition      (res_void_ptr_t,
                                                             int64_t,
                                                             wait_type);

static void           c_io_wbuffering_waitLECondition       (res_void_ptr_t,
                                                             int64_t,
                                                             wait_type);

static inline int64_t wait_type_to_int64                    (res_void_ptr_t,
                                                             res_wt_ptr_t);

static int64_t        c_io_wbuffering_out_via_helper        (int64_t,
                                                             char *,
                                                             size_t,
                                                             size_t);

static int64_t        c_io_wbuffering_out_buffer_to_disk    (int64_t,
                                                             char *,
                                                             size_t,
                                                             size_t);

static void           c_io_wbuffering_helper_process_buffer (size_t);

/*----------------------------------------------------------------------------*/
/* Layer routine implemetations                                               */
/* Note: all these routines are always called by the main thread, even in the */
/*       presence of a helper thread.                                         */
/*----------------------------------------------------------------------------*/

int64_t c_io_wbuffering_init(c_io_init_t *init_params)
{
  size_t i;

  wbufLayer = c_io_init_next_layer(CIO_NEXT_LAYER);

  wbufLayer->helper_start = c_io_wbuffering_helper;
  wbufLayer->helper_stop  = c_io_wbuffering_release_helper;

  writethreadlock = newLock();

  c_io_wbuffering_data = malloc((size_t)(init_params->max_unit)
                                * sizeof(*c_io_wbuffering_data));

  if (c_io_wbuffering_data==NULL)
  {
    const size_t alloc_size = (size_t)(init_params->max_unit)
                                * sizeof(*c_io_wbuffering_data);

    c_io_printf(-1,
                "init: wbuffering: Cannot allocate required %zu bytes"
                " for internal use",
                alloc_size);
    ereport("c_io_wbuffering_init",1,"Memory allocation failure");
  }

  c_io_wbuffering_counters = malloc(sizeof(*c_io_wbuffering_counters));

  if (c_io_wbuffering_counters==NULL)
  {
    c_io_printf(-1, "init: wbuffering:"
                " Cannot allocate required %zu bytes for internal use",
                sizeof(c_io_wbuffering_counters_t));
    ereport("c_io_wbuffering_init",1,"Memory allocation failure");
  }

  c_io_wbuffering_counters->written=0;
  c_io_wbuffering_counters->write_calls=0;
  c_io_wbuffering_counters->wait_calls=0;
  c_io_wbuffering_counters->misses=0;
  c_io_wbuffering_counters->hits=0;
  c_io_wbuffering_counters->b_used=0;
  c_io_wbuffering_counters->wait=0.0;

  unit_max = init_params->max_unit;

  if (init_params->wbuffering_buffer_size > 0)
  {
    defaultBufferSz = (size_t)(init_params->wbuffering_buffer_size);
  }
  else
  {
    defaultBufferSz = DEFAULT_WBUFFER_SIZE;
  }

  c_io_printf(-1,"init: wbuffering: Buffer size=%zu bytes",defaultBufferSz);

  for (i=0;i<init_params->max_unit;i++)
  {
    c_io_wbuffering_reset(i);
  }

  return C_IO_INIT(init_params);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_wbuffering_fini(void)
{
  double total_data_block_eqiv;

  if (wbufLayer->helper_state == c_io_running)
  {
    c_io_wbuffering_release_helper();
  }

  releaseLock(&writethreadlock);

  total_data_block_eqiv = (double)(c_io_wbuffering_counters->b_used)
                           / (double)(defaultBufferSz);

  c_io_printf(-1," Wbuffering usage");
  c_io_printf(-1,"  Number of write calls                %" PRId64,
              c_io_wbuffering_counters->write_calls);
  c_io_printf(-1,"  Number of actual write events        %" PRId64,
              c_io_wbuffering_counters->written);
  c_io_printf(-1,"  Total data written (blocks eqv.)     %.2f",
              total_data_block_eqiv);
  c_io_printf(-1,"  Buffer Hits                          %zu",
              c_io_wbuffering_counters->hits);
  c_io_printf(-1,"  Buffer Misses                        %zu",
              c_io_wbuffering_counters->misses);
  c_io_printf(-1,"  Wait calls                           %" PRId64,
              c_io_wbuffering_counters->wait_calls);
  c_io_printf(-1,"  Wait time                            %.2f secs",
              c_io_wbuffering_counters->wait);

  c_io_printf(-1,"\n");

  free(c_io_wbuffering_data);
  free (c_io_wbuffering_counters);

  return C_IO_FINI();
}

/*----------------------------------------------------------------------------*/

int64_t c_io_wbuffering_open(int64_t unit,
                             char *filename,
                             int64_t fileStatus,
                             int64_t fileMode)
{
  c_io_wbuffering_data_t_cr_ptr buffunit = &c_io_wbuffering_data[unit];

  c_io_wbuffering_check_latch();

  c_io_wbuffering_reset((size_t)unit);

  buffunit->buffer = malloc( defaultBufferSz * sizeof(*(buffunit->buffer)) );

  if (buffunit->buffer==NULL)
  {
    c_io_printf(-1,"c_io_wbuffering: open: failed to allocate write buffer");
    return c_io_fail;
  }

  buffunit->bufferMax = defaultBufferSz;

  return C_IO_OPEN(unit,filename,fileStatus,fileMode);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_wbuffering_close(int64_t unit)
{
  c_io_wbuffering_data_t_cr_ptr buffunit = &c_io_wbuffering_data[unit];

  /* first flush buffers */
  if(c_io_wbuffering_flush_check(unit)==c_io_fail)
  {
    ereport("c_io_wbuffering_close", 1,
            "Failed to flush buffer on c_io_wbuffering_close() call");
  }

  /* update counters */
  c_io_wbuffering_counters->hits   += buffunit->bufferHits;
  c_io_wbuffering_counters->misses += buffunit->bufferMiss;

  /* deallocate buffer */
  free(buffunit->buffer);
  c_io_wbuffering_reset((size_t)unit);

  /* finally close file */
  if(C_IO_CLOSE(unit)==c_io_fail)
  {
    c_io_printf(unit,"c_io_wbuffering: Close: ERROR in close from lower layer");
    return c_io_fail;
  }

  return c_io_success;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_wbuffering_in(int64_t unit,
                           char *array,
                           int64_t len,
                           int64_t word_len)
{
  /* Any outstanding writes are completed when a mode change occurs - so we
   * know any previous cached writes have already reached the next layer by
   * the time we get to a read, as write-to-read mode change must have happened
   * in-between. Hence we do not have to worry about anything more than
   * passing the read on to the next layer here.
   */

  return C_IO_IN(unit,array,len,word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_wbuffering_out(int64_t unit,
                            char *array,
                            int64_t len,
                            int64_t word_len)
{
  c_io_wbuffering_check_latch();

  c_io_wbuffering_counters->write_calls++;
  c_io_wbuffering_counters->b_used += (len*word_len);

  /* use the helper thread? */
  if(wbufLayer->helper_state == c_io_running)
  {
    return c_io_wbuffering_out_via_helper(unit,
                                          array,
                                          (size_t)len,
                                          (size_t)word_len);
  }

  pending_helper_writes=pending_helper_writes+1;
  pending_helper_writes_size=pending_helper_writes_size+(size_t)(len*word_len);

  return c_io_wbuffering_out_buffer_to_disk(unit,
                                            array,
                                            (size_t)len,
                                            (size_t)word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_wbuffering_getpos(int64_t unit)
{
  int64_t pos = C_IO_GETPOS(unit);

  if (c_io_attributes[unit].mode==writing)
  {
    pos += (int64_t)(c_io_wbuffering_data[unit].bufferPosition);
  }

  return pos;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_wbuffering_setpos(int64_t unit,
                               int64_t pos,
                               int64_t word_len)
{
  if (c_io_attributes[unit].mode==writing)
  {
    if (c_io_wbuffering_flush_check(unit)==c_io_fail)
    {
      return c_io_fail;
    }
  }
  return C_IO_SETPOS(unit,pos,word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_wbuffering_change_mode(int64_t unit,
                                    int64_t newMode)
{
  if (c_io_attributes[unit].mode==writing)
  {
    if (c_io_wbuffering_flush_check(unit)==c_io_fail)
    {
      return c_io_fail;
    }
  }

  return C_IO_CHANGE_MODE(unit,newMode);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_wbuffering_sync(int64_t unit)
{
  if (c_io_attributes[unit].mode==writing)
  {
    if (c_io_wbuffering_flush_check(unit)==c_io_fail)
    {
      return c_io_fail;
    }
  }

  return C_IO_SYNC(unit);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_wbuffering_query(int64_t unit,
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

/* Called by: main or helper thread */
static int64_t c_io_wbuffering_out_buffer_to_disk(int64_t unit,
                                                  char *array,
                                                  size_t len,
                                                  size_t word_len)
{
  size_t freeSpace;
  size_t chunkSize;

  size_t remaining = len*word_len;
  size_t ammountWritten = 0;

  c_io_wbuffering_data_t_cr_ptr attr = &c_io_wbuffering_data[unit];

  char *current_pos = array;

  /* Loop untill we have got rid of the payload */
  while (remaining>0)
  {
    freeSpace = attr->bufferMax - attr->bufferPosition;

    /* Check if the write buffer is full */
    if (freeSpace==0)
    {
      if (c_io_wbuffering_flush_nocheck(unit)==c_io_fail)
      {
        ereport("c_io_wbuffering_init", 1,
                "Failed to flush buffer on c_io_wbuffering_out() call");
      }

      freeSpace = attr->bufferMax - attr->bufferPosition;
    }

    chunkSize=freeSpace;

    if (chunkSize>remaining)
    {
      chunkSize=remaining;
    }

    if (chunkSize==attr->bufferMax)
    {
      /* The buffer is empty and we have at least a
         whole buffers worth of output. */
      /*      c_io_printf(unit,"c_io_wbuffering: direct write \n"); */
      ammountWritten += (size_t)(C_IO_OUT(unit,
                                          current_pos,
                                          (int64_t)(chunkSize/word_len),
                                          (int64_t)word_len));

      c_io_wbuffering_counters->written++;
      attr->bufferMiss++;
    }
    else
    {
      /* The buffer is partially full, or we have less
         than a buffers worth of stuff so we add to it */
      /*      c_io_printf(unit,"c_io_wbuffering: copy to buffer \n"); */
      memcpy(&attr->buffer[attr->bufferPosition],
             current_pos,
             chunkSize);
      attr->bufferPosition += chunkSize;
      ammountWritten += chunkSize/word_len;
      attr->bufferHits++;
    }

    remaining   -= chunkSize;
    current_pos += chunkSize;
  }

  pending_helper_writes = pending_helper_writes - 1;
  pending_helper_writes_size = pending_helper_writes_size - ammountWritten;

  return (int64_t)ammountWritten;
}

/*----------------------------------------------------------------------------*/

/* Called by: main thread only */
static int64_t c_io_wbuffering_out_via_helper(int64_t unit,
                                              char *array,
                                              size_t len,
                                              size_t word_len)
{
  size_t remaining = (len*word_len);

  c_io_wbuffering_helperbuffers_t **restrict lastbuff = &lastbuffer[unit];

  if( helperbuffers==NULL || (*lastbuff)==NULL )
  {
    ereport("c_io_wbuffering_out", 1,
            "Attempting to buffer with helper thread,"
            " but no helper-thread buffers exist.");
  }

  /* pause to allow backlog to clear if the queue is too long/large */
  c_io_wbuffering_waitLECondition(&pending_helper_writes,
                                  MAX_HELPER_QUEUE,
                                  wt_size_t);
  c_io_wbuffering_waitLECondition(&pending_helper_writes_size,
                                  MAX_HELPER_QUEUE_SIZE,
                                  wt_size_t);

  pending_helper_writes = pending_helper_writes + 1;
  pending_helper_writes_size = pending_helper_writes_size + remaining;

  threadFlush();

  /* calloc() because this initialises to zero -
   * This guarantees lastbuffer->ready=(bool)0 (i.e. false),
   * hence the helper thread cannot work on this buffer too early.
   */
  (*lastbuff)->next = calloc(1,sizeof(*((*lastbuff)->next)));

  if ((*lastbuff)->next==NULL)
  {
    c_io_printf(-1, "wbuffering: out_via_helper:"
                " Cannot allocate required %zu bytes for internal use",
                sizeof(*((*lastbuff)->next)));
    ereport("c_io_wbuffering_init",1,"Memory allocation failure");
  }

  (*lastbuff) = (*lastbuff)->next;

  (*lastbuff)->buffer_len = remaining;

  (*lastbuff)->buffer = malloc( remaining * sizeof(*((*lastbuff)->buffer)) );

  if ((*lastbuff)->buffer==NULL)
  {
    c_io_printf(-1, "wbuffering: out_via_helper:"
                " Cannot allocate required %zu bytes for internal use",
                remaining);
    ereport("c_io_wbuffering_init",1,"Memory allocation failure");
  }
  memcpy((*lastbuff)->buffer,array,remaining);

  (*lastbuff)->ready=true;

  threadFlush();

  return (int64_t)len;
}

/*----------------------------------------------------------------------------*/

/* Called by: main thread only */
static void c_io_wbuffering_reset(size_t unit)
{
    c_io_wbuffering_data[unit].bufferMax      = 0;
    c_io_wbuffering_data[unit].bufferSize     = 0;
    c_io_wbuffering_data[unit].bufferPosition = 0;
    c_io_wbuffering_data[unit].bufferHits     = 0;
    c_io_wbuffering_data[unit].bufferMiss     = 0;
    c_io_wbuffering_data[unit].buffer         = NULL;
}

/*----------------------------------------------------------------------------*/

/* Called by: main thread only */
static int64_t c_io_wbuffering_flush_check(int64_t unit)
{
  c_io_wbuffering_check_latch();

  if (wbufLayer->helper_state == c_io_running)
  {
    /* complete cached writes for this unit */
    priority_unit = (size_t)unit;
    Lock(&writethreadlock);
    {
      priority_buffer = true;
      threadFlush();
    }
    unLock(&writethreadlock);

    c_io_wbuffering_waitNotCondition(&priority_buffer,true,wt_prioritystates);
  }

  return c_io_wbuffering_flush_nocheck(unit);
}

/*----------------------------------------------------------------------------*/

/* Called by: main or helper thread */
static int64_t c_io_wbuffering_flush_nocheck(int64_t unit)
{
  c_io_wbuffering_data_t_cr_ptr attr = &c_io_wbuffering_data[unit];

  /* check for partial buffers */

  if ((attr->buffer!=NULL) && (attr->bufferPosition>0))
  {
    int64_t wCheck;

    wCheck = C_IO_OUT(unit,
                      attr->buffer,
                      (int64_t)(attr->bufferPosition),
                      WORD8BYTES);

    if (wCheck != (int64_t)(attr->bufferPosition))
    {
      c_io_printf(unit, "c_io_wbuffering: flush: ERROR:"
                  " writing pending buffer to disk");

      return c_io_fail;
    }

    c_io_wbuffering_counters->written++;

    attr->bufferPosition=0;
  }

  return c_io_success;
}

/*----------------------------------------------------------------------------*/

/* Called by: helper thread only */
static void c_io_wbuffering_helper(void)
{
  /* The logic of the wbuffering_helper works provided at most only one main
   * and one helper threads exist at any time.
   * Consequently, we need to protect this routine such that multiple threads
   * cannot simultaneously be created due to race conditions on
   * wbufLayer->helper_state
   */
  helperStates tmp_state;
  size_t i;

  Lock(&writethreadlock);
  {
    /* all checks and updates of wbufLayer->helper_state will happen in the
     * locked region to prevent race conditions
     */

    threadFlush();

    tmp_state = wbufLayer->helper_state;

    /* Do we have a bad initial state (i.e state is not c_io_none)?
     * If no, initiate the handshake.
     * If yes, do nothing - this will be handled outside of the lock.
     */
    if (tmp_state == c_io_none)
    {
      c_io_printf(-1,"wbuffering: helper: Knock knock!");
      wbufLayer->helper_state=c_io_present;
    }

    threadFlush();
  }
  unLock(&writethreadlock);

  /* handle any bad initial state from within the lock */
  if (tmp_state != c_io_none)
  {
    c_io_printf(-1, "wbuffering: helper:"
                " attaching has failed. (Bad initial condition,"
                " or another helper thread is already attached.)");
    return;
  }

  /* now allocate the helper's buffers for each unit.
   * helperbuffers is of type (c_io_wbuffering_helperbuffers_t **)
   * we are therefore allocating space for an array of pointers
   */
  helperbuffers = malloc(unit_max * sizeof(*helperbuffers));

  if (helperbuffers==NULL)
  {
    c_io_printf(-1, "wbuffering: helper:"
                " Cannot allocate required %zu bytes for internal use",
                unit_max*sizeof(*helperbuffers));
    ereport("c_io_wbuffering_init",1,"Memory allocation failure");
  }

  /* now allocate the pointers to the last buffer for each unit.
   * lastbuffer is of type (c_io_wbuffering_helperbuffers_t **)
   * we are therefore allocating space for an array of pointers
   */
  lastbuffer = malloc(unit_max * sizeof(*lastbuffer));

  if (lastbuffer==NULL)
  {
    c_io_printf(-1, "wbuffering: helper:"
                " Cannot allocate required %zu bytes for internal use",
                unit_max*sizeof(*lastbuffer));
    ereport("c_io_wbuffering_init",1,"Memory allocation failure");
  }

  for (i=0;i<unit_max;i++)
  {
    /* now allocate the buffer for unit i
     * helperbuffers[i] is of type (c_io_wbuffering_helperbuffers_t *)
     * we are therefore allocating space for individual buffer structures
     */

    /* calloc() to ensure helperbuffers[i]->ready = (bool)0 i.e. false */
    helperbuffers[i]=calloc(1,sizeof(*helperbuffers[i]));

    if (helperbuffers[i]==NULL)
    {
      c_io_printf(-1, "wbuffering: helper:"
                  " Cannot allocate required %zu bytes for internal use",
                  unit_max*sizeof(*helperbuffers[i]));
      ereport("c_io_wbuffering_init",1,"Memory allocation failure");
    }

    /* this first buffer is a dummy with no payload in order to anchor future
     * buffers against, so set its status to ready and written.
     */

    helperbuffers[i]->written=true;
    helperbuffers[i]->ready=true;

    /* ...and then point lastbuffer[i] at helperbuffers[i] which is currently
     * the last (only) buffer structure
     */
    lastbuffer[i]=helperbuffers[i];
  }

  while (wbufLayer->helper_state != c_io_attaching)
  {
    threadFlush();
  }

  c_io_printf(-1,"wbuffering: helper: A helper thread.");
  wbufLayer->helper_state=c_io_running;
  threadFlush();

  c_io_wbuffering_helper_loop();
  c_io_printf(-1,"wbuffering: helper: I'm gone.");
}

/*----------------------------------------------------------------------------*/

/* Called by: helper thread only */
static void c_io_wbuffering_helper_loop(void)
{
  threadFlush();
  size_t i=0;

  while (wbufLayer->helper_state != c_io_terminate)
  {
    threadFlush();
    usleep(c_io_sleep_waittime);

    if (priority_buffer)
    {
      /* jump straight to the priority unit */
      i = priority_unit;
    }
    else
    {
      /* rotate over [0 ... c_io_wbuffering_counters->unit_max] */
      i=(i+1)%(unit_max);
    }

    c_io_wbuffering_helper_process_buffer(i);

  }

  /* Make sure all the buffers are written as we are terminating */
  for (i=0;i<unit_max;i++)
  {
    while (
            /* there is another buffer in the list... */
            helperbuffers[i]->next!=NULL

            /*...OR the current buffer is not written */
            || ((helperbuffers[i]->ready) && !(helperbuffers[i]->written))
          )
    {
      c_io_wbuffering_helper_process_buffer(i);
    }
  }

  wbufLayer->helper_state = c_io_none;
  threadFlush();
}

/*----------------------------------------------------------------------------*/

/* Called by: helper thread only */
static void c_io_wbuffering_helper_process_buffer(size_t unit)
{
  c_io_wbuffering_helperbuffers_t **helperbuff = &helperbuffers[unit];

  threadFlush();

  if ((*helperbuff)->ready)
  {
      if ((*helperbuff)->written)
      {
         if ((*helperbuff)->next==NULL)
         {
           Lock(&writethreadlock);
           {
             threadFlush();
             if(unit==priority_unit && priority_buffer)
             {
               /* no more buffers for this unit so we no longer prioritise it */
               priority_buffer = false;
               threadFlush();
             }
           }
           unLock(&writethreadlock);
         }
         else
         {
           /* clean up this buffer, and move onto the next in the list */

           c_io_wbuffering_helperbuffers_t *cleanup;

           cleanup = (*helperbuff);
           (*helperbuff) = (*helperbuff)->next;

           free(cleanup);
         }

         threadFlush();
      }
      else
      {
        int64_t ret;

        ret = c_io_wbuffering_out_buffer_to_disk((int64_t)unit,
                                                 (*helperbuff)->buffer,
                                                 (*helperbuff)->buffer_len,
                                                 WORD8BYTES);

        if(ret==(int64_t)(*helperbuff)->buffer_len)
        {
          free((*helperbuff)->buffer);
          (*helperbuff)->written=true;
          threadFlush();
        }
        else
        {
          ereport("c_io_wbuffering_helper_loop",1,"Short write from buffer.");
        }
      }
  }
}

/*----------------------------------------------------------------------------*/

/* Called by: main thread only */
static void c_io_wbuffering_release_helper(void)
{

  /* In case another thread had requested to attach as a helper, but never made
   * it into the c_io_running state because no subsequent check latch call was
   * made - make a call to c_io_wbuffering_check_latch() now.
   */
  c_io_wbuffering_check_latch();

  if (wbufLayer->helper_state != c_io_running)
  {
    c_io_printf(-1,"wbuffering: release_helper: nothing to release!");
    return;
  }

  c_io_printf(-1,"wbuffering: release_helper: Go away!");

  /*
   * Set the latch out variable
   */
  wbufLayer->helper_state = c_io_terminate;
  threadFlush();

  /*
   * Wait for the other thread to notice
   */
  c_io_wbuffering_waitNotCondition(&wbufLayer->helper_state,
                                   c_io_terminate,
                                   wt_threadstates);

  free(helperbuffers);
  free(lastbuffer);

  c_io_printf(-1,"wbuffering: release_helper: Bye!");
}

/*----------------------------------------------------------------------------*/

/* Called by: main thread only */
static void c_io_wbuffering_check_latch(void)
{
  threadFlush();

  if (wbufLayer->helper_state==c_io_present)
  {
    /* There is a helper thread that has been attached, but is waiting in the
     * c_io_present state until a latch check is called.
     * We should now start that helper thread running.
     */
    c_io_printf(-1,"wbuffering: check_latch: Who's there?");
    wbufLayer->helper_state=c_io_attaching;
    threadFlush();

    while(wbufLayer->helper_state!=c_io_running)
    {
      threadFlush();
    }
  }
}

/*----------------------------------------------------------------------------*/

/* Called by: main or helper thread */
static void c_io_wbuffering_waitNotCondition(res_void_ptr_t flag,
                                             int64_t cond,
                                             wait_type wt)
{
  /*
   *  wait for the contents of flag to be not equal to cond.
   */

  double t;

  threadFlush();

  if (wait_type_to_int64(flag,&wt) != cond)
  {
    return;
  }

  t = c_io_getTime();

  while (wait_type_to_int64(flag,&wt) == cond)
  {
    usleep(c_io_sleep_waittime);
    threadFlush();
  }

  c_io_wbuffering_counters->wait += (c_io_getTime()-t);
  c_io_wbuffering_counters->wait_calls++;
}

/*----------------------------------------------------------------------------*/

/* Called by: main or helper thread */
static void c_io_wbuffering_waitLECondition(res_void_ptr_t flag,
                                            int64_t cond,
                                            wait_type wt)
{
  /*
   *  wait for the contents of flag to be less than or equal to cond.
   */

  double t;

  threadFlush();

  if (wait_type_to_int64(flag,&wt) <= cond)
  {
    return;
  }

  t = c_io_getTime();

  while (wait_type_to_int64(flag,&wt) > cond)
  {
    usleep(c_io_sleep_waittime);
    threadFlush();
  }

  c_io_wbuffering_counters->wait += (c_io_getTime()-t);
  c_io_wbuffering_counters->wait_calls++;
}

/*----------------------------------------------------------------------------*/

/* Called by: main or helper thread */
static inline int64_t wait_type_to_int64(res_void_ptr_t flag,
                                         res_wt_ptr_t wt)
{
  int64_t ret=-1;

  switch(*wt)
  {
    case wt_prioritystates:
    {
      ret=(int64_t)*(volatile const bool *)flag;
    }
    break;

    case wt_threadstates:
    {
      ret=(int64_t)*(volatile const helperStates *)flag;
    }
    break;

    case wt_size_t:
    {
      ret=(int64_t)*(volatile const size_t *)flag;
    }
    break;
  }

  return ret;
}

#endif

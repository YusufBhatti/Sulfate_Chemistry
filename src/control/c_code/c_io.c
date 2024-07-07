#if defined(C95_2B)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* This unit contains the entry points to the Portio2B layers.                */
/* Please see UMPD S05 for more information                                   */

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#define CIO_THIS_LAYER CIO_NO_LAYER
#define CIO_ROOT_LAYER

#include "c_io_layers.h"
#include "c_io_nextlayer.h"
#include "c_io_internal.h"
#include "c_fort2c_prototypes.h"
#include "c_io_errcodes.h"
#include "pio_umprint.h"
#include "c_fort2c_types.h"

/* -------------------------------------------------------------------------- */
/*  Global Variables and Prototypes                                           */
/* -------------------------------------------------------------------------- */

/* initialised global constants                                               */

const int c_io_success    = 0;
const int c_io_fail       = -1;
const int c_io_file_end   = -1;
const int c_io_file_nopos = -999;

char modeText[numOpenStates][numModes][4] =
                                     /* read-only / read & write / write-only */
/* old files                      */ {{   "rb"  ,      "r+b" ,       "r+b"   },
/* new files                      */  {   "err" ,      "w+b" ,       "wb"    }};

/* read-only new files implies an error:  "err"                               */
/* no mode exists purely for write only on old files, so use "r+b"            */

/* global variables                                                           */

int64_t            c_io_attributes_size = 0;
int64_t            c_io_outLevel        = outLevOper;
bool               c_io_initialised     = false;

c_io_init_t       *params               = NULL;
c_io_attributes_t *c_io_attributes      = NULL;

/* Internal and utility function prototypes                                   */

static void    c_io_resetAttributes (int64_t);
static void    c_io_detachHelper    (int64_t);
static void    c_io_print_stack     (void);
static int64_t c_io_effect_setpos   (int64_t);
static void    c_io_allocate_query  (query_type_t, void **, void **);
static void    c_io_free_query      (query_type_t, void *, void *);

/* c_io_printf_common() uses the C_IO_PRINTF_COMMON_DECORATOR decorator,
 * which is explained and defined in c_io_internal.h
 */
static void    c_io_printf_common(int64_t,
                                  const char *,
                                  va_list *,
                                  bool
                                 ) C_IO_PRINTF_COMMON_DECORATOR;

/* -------------------------------------------------------------------------- */
/*  Procedures                                                                */
/* -------------------------------------------------------------------------- */

/* print a varargs message to the UM                                          */

void c_io_printf(int64_t unit,const char *fmt, ... )
{
  va_list argp;
  va_start(argp,fmt);
  c_io_printf_common(unit,fmt,&argp,false);
  va_end(argp);
}

void c_io_printf_stderr(int64_t unit,const char *fmt, ... )
{
  va_list argp;
  va_start(argp,fmt);
  c_io_printf_common(unit,fmt,&argp,true);
  va_end(argp);
}


void c_io_printf_common(int64_t unit,
                        const char *fmt,
                        va_list *argp,
                        bool use_stderr)
{
  char mppioString[MAX_OUTSTR];
  char prefix[MAX_OUTSTR];

  if (unit==-1)
  {
    snprintf(prefix,MAX_OUTSTR,"c_io      : ");
  }
  else
  {
    snprintf(prefix,MAX_OUTSTR,"c_io (%3" PRId64 "): ",unit);
  }

  vsnprintf(mppioString,MAX_OUTSTR,fmt,*argp);

  umPrint(mppioString,"c_io", .stdErrorToo=use_stderr, .UsrPrefix=prefix );
}

/* initialise all layers                                                      */

int64_t c_io_init(int64_t min_unit,
                  int64_t max_unit,
                  int64_t wbuff_size,
                  int64_t rbuff_size,
                  int64_t rbuff_count,
                  bool    rbuff_update,
                  int64_t rbuff_prefetch,
                  bool    io_timer)
{
  int i;
  c_io_outLevel=outLevOper;

  /*
   * Early bail out on multiple calls
   */
  if (c_io_initialised)return 0;

  /*
   * Initialise the parameterisation structure
   */
  params=(c_io_init_t *)malloc(sizeof (c_io_init_t));
  max_unit++;
  params->min_unit=(size_t)min_unit;
  params->max_unit=(size_t)max_unit;
  params->wbuffering_buffer_size=wbuff_size;
  params->rbuffering_buffer_size=rbuff_size;
  params->rbuffering_buffer_count=rbuff_count;
  params->rbuffering_buffer_update=rbuff_update;
  params->rbuffering_buffer_prefetch=rbuff_prefetch;
  params->timing_switch=io_timer;
  params->layers=NULL;
  /*
   * Populate the top level of the layer description table
   */
  c_io_init_next_layer(CIO_NEXT_LAYER);

  /*
   * Ceate the attributes table and reset.
   */
  c_io_attributes_size=max_unit;
  c_io_attributes=
    (c_io_attributes_t*)malloc((size_t)(c_io_attributes_size)
    * sizeof(c_io_attributes_t));
  for (i=0;i<max_unit;i++)c_io_resetAttributes(i);

  /*
   * Start the init pipeline.
   */
  C_IO_INIT(params);

  /*
   * Tidy up
   */
  c_io_initialised=true;
  c_io_printf
    (-1,"init: top: Initialisation done");
  c_io_print_stack();
  return 0;
}

/* shutdown all layers and perform any final tasks */

int64_t c_io_fini()
{
  size_t i;
  int64_t ret;

  ret = c_io_success;

  if (c_io_initialised)
  {
    c_io_layer_t *layer, *lastlayer=NULL;

    for (i=(size_t)(params->min_unit);i<(size_t)(params->max_unit);i++)
    {
      if (c_io_attributes[i].fileHandle!=NULL)
      {
        c_io_printf ((int64_t)i,"At shutdown, unit appears open");
        c_io_printf ((int64_t)i,"and attached to file %s",
                     c_io_attributes[i].filename);
      }
    }

    ret=C_IO_FINI();

    layer=params->layers;

    /* cycle through layers deallocating memory */
    while(layer->next != NULL)
    {
      free(layer->name);

      if(lastlayer!=NULL)
      {
        /* free the previous layer */
        free(lastlayer);
      }

      lastlayer=layer;
      layer=layer->next;
    }

    /* if there was more than one layer, free the bottom one */
    if(lastlayer!=NULL)
    {
      free(lastlayer);
    }
    else
    {
      free(layer->name);
      free(layer);
    }

    free(c_io_attributes);
    free(params);

    c_io_initialised=false;
  }

  return ret;
}

/* reset attributes for a given unit */

void c_io_resetAttributes(int64_t unit)
{
  c_io_attributes[unit].fileHandle=NULL;
  c_io_attributes[unit].filename=NULL;
  c_io_attributes[unit].mode=idle;
  c_io_attributes[unit].next_position=c_io_file_nopos;
  c_io_attributes[unit].modified_on_disk=false;
  c_io_attributes[unit].end_of_file=false;
}

/* open a file */

int64_t c_io_open(int64_t unit,
                  char *filename,
                  int64_t file_name_len,
                  int64_t fileMode)
{
  int64_t result=c_io_success;
  int fileStatus;
  size_t filename_len=strlen(filename)+1;
  c_io_printf(unit,"Open: File=%s",filename);
  if (c_io_attributes[unit].fileHandle!=NULL)
  {

    /*
     *  If the unit is already open which respects the behaviour expected
     *  by the model
     */

    c_io_printf(unit,"Open: ERROR: unit appears open");
    if (strcmp(filename,c_io_attributes[unit].filename) == 0){
      c_io_printf(unit,"Open: Reseting file position to zero");
      c_io_setpos(unit,0,WORD8BYTES);
      return ERR_FILE_OPEN;
    }
    else{
      c_io_printf(unit,"Open: Closing existing file %s",
                  c_io_attributes[unit].filename);
      c_io_close(unit);
    }
  }

  /* We need to query if the file already exists in order to apply the correct
   * status to the file.
   * This envolves a metadata query of type "qry_file_exist":
   * this query requires the input of the filename, and returns a bool (true if
   * the file already exists).
   */

  /* First allocate and initialise the query */
  in_pack_fname_t *in_pack=NULL;
  out_pack_bool_t *out_pack=NULL;

  c_io_allocate_query(qry_file_exist, (void **)&in_pack, (void **)&out_pack);

  snprintf(in_pack->filename, MAX_OUTSTR, "%s", filename);

  /* Now execute the query */
  if ( C_IO_QUERY(unit,
                  qry_file_exist,
                  (void *)in_pack,
                  (void *)out_pack)
       == c_io_fail )
  {
    c_io_printf(unit,"Open: ERROR: failed to query file state (%s)",filename);
    return ERR_ACCESS;
  }

  /* Now use the result of the query to determine the file status */
  if (out_pack->res)
  {
    fileStatus = oldFile;
  }
  else
  {
    fileStatus = newFile;
  }

  /* clean up the query structures */
  c_io_free_query(qry_file_exist, (void *)in_pack, (void *)out_pack);

  c_io_resetAttributes(unit);
  C_IO_OPEN(unit, filename, fileStatus, fileMode);

  if (c_io_attributes[unit].fileHandle==NULL)
  {
    c_io_printf(unit,"Open: ERROR: file open failed (%s)",filename);
    return ERR_ACCESS;
  }
  c_io_attributes[unit].filename=malloc(filename_len);
  strncpy(c_io_attributes[unit].filename,filename,filename_len);

  if (fileStatus==oldFile)
    c_io_printf(unit,"Open: File exists (%" PRId64 " bytes)",
                c_io_getExtent(unit));
  else
    c_io_printf(unit,"Open: File is new");

  //Avoid Unused Arguments
  (void)(file_name_len);

  return result;
}

/* Close a specified unit */

int64_t c_io_close(int64_t unit)
{
  int64_t error;
  c_io_printf(unit,"Close");
  if (c_io_attributes[unit].fileHandle==NULL)
  {
    return ERR_FILE_CLOSED;
  }
  error=C_IO_CLOSE(unit);

  if (error!=c_io_success)
    c_io_printf(unit,"Close: ERROR in close");

  free(c_io_attributes[unit].filename);
  c_io_resetAttributes(unit);
  return error;
}

/* Delete a file */

int64_t c_io_delete(char *fname)
{
  size_t i;
  c_io_printf(-1,"Deleting file: %s",fname);
  for (i=0;i<params->max_unit;i++){
    if (c_io_attributes[i].fileHandle!=NULL)
    {
      if (strcmp(fname,c_io_attributes[i].filename)==0){
        c_io_printf((int64_t)i,"Delete: File appears open, closing");
        c_io_close((int64_t)i);
      }
    }
  }
  return remove(fname);
}

/* Read data from a unit into a user provided buffer */

int64_t c_io_in(int64_t unit, char *array, int64_t len, int64_t word_len )
{
  int64_t ammountRead=0;
  if (c_io_attributes[unit].fileHandle==NULL)
  {
    return ERR_FILE_CLOSED;
  }

  /*
  c_io_printf(unit,"in: %lld items of length %lld"
              ,LLD(len),LLD(word_len));
  */
  if (c_io_attributes[unit].mode!=reading)
    c_io_change_mode(unit,reading);

  if (c_io_attributes[unit].next_position!=c_io_file_nopos)
    c_io_effect_setpos(unit);

  ammountRead=C_IO_IN(unit,array,len,word_len);

  if (ammountRead!=len)
    c_io_printf(unit,"in: WARNING only read %" PRId64 " items (of %" PRId64
                ") from 0x%" PRIuPTR,
                ammountRead,len,(uintptr_t)c_io_attributes[unit].fileHandle);
  return ammountRead;
}

/* Write data to a unit from a user provided buffer */

int64_t c_io_out(int64_t unit, char *array, int64_t len, int64_t word_len )
{
  int64_t ammountWritten=0;

  /*
  c_io_printf(unit,"out: %lld items of length %lld"
              ,LLD(len),LLD(word_len));
  */

  if (c_io_attributes[unit].fileHandle==NULL)
  {
    return ERR_FILE_CLOSED;
  }

  if (c_io_attributes[unit].mode!=writing)
    c_io_change_mode(unit,writing);

  if (c_io_attributes[unit].next_position!=c_io_file_nopos)
    c_io_effect_setpos(unit);

  c_io_attributes[unit].modified_on_disk=true;

  ammountWritten=C_IO_OUT(unit,array,len,word_len);

  if (ammountWritten!=len)
    c_io_printf(unit,
                "out: WARNING only wrote %" PRId64 " items (of %" PRId64 ")",
                ammountWritten,len);

  return ammountWritten;
}

/* Set the next file position (cache the position rather than moving the file
 * position marker)
 */
int64_t c_io_setpos(int64_t unit,int64_t pos, int64_t word_len )
{
  /* cache the setpos for performance reasons.
   * This means any call to a lower layer from this layer *must* call
   * c_io_effect_setpos() if it might update disk position in that lower level
   * (i.e. read, write, change mode, setpos, or getpos)
   */
  if (pos==c_io_file_end)
  {
    c_io_attributes[unit].next_position=c_io_file_end;
  }
  else
  {
    c_io_attributes[unit].next_position=pos*word_len;
  }

  return c_io_success;
}

/* Set the next file position (explicitly move the file postion marker to the
 * cached position)
 */
int64_t c_io_effect_setpos(int64_t unit)
{
  int64_t code;

  if (c_io_attributes[unit].fileHandle==NULL)
  {
    return ERR_FILE_CLOSED;
  }

  if (c_io_attributes[unit].next_position==c_io_file_nopos)
    return 0;

  /* Check if we are already in the right position.
   * We need a specified position for this - i.e. not EOF.
   */
  if (c_io_attributes[unit].next_position!=c_io_file_end)
  {
    if (C_IO_GETPOS(unit)==c_io_attributes[unit].next_position)
    {
      c_io_attributes[unit].next_position=c_io_file_nopos;
      return 0;
    }
  }

  /*
   * c_io_printf(unit,"setpos: position=%lld",
   *             LLD(c_io_attributes[unit].next_position));
   */

  code=C_IO_SETPOS(unit,c_io_attributes[unit].next_position,WORD8BYTES);
  if (code==c_io_success)
  {
    c_io_attributes[unit].next_position=c_io_file_nopos;
  }

  return code;
}

/* get the file position */

int64_t c_io_getpos(int64_t unit,int64_t word_len )
{
  int64_t pos;
  int64_t word_index;

  if (c_io_attributes[unit].next_position==c_io_file_nopos)
  {
    /* there is no cached postition - we must get the real position */
    pos=C_IO_GETPOS(unit);
  }
  else
  {
    if (c_io_attributes[unit].next_position==c_io_file_end)
    {
      /* the cached position is EOF, so we must explicitly move the
       * file marker before we can read the position.
       */
      c_io_effect_setpos(unit);
      pos=C_IO_GETPOS(unit);
    }
    else
    {
      /* we can just use the cached position */
      pos=c_io_attributes[unit].next_position;
    }
  }

  word_index=pos/word_len;

  if (word_index*word_len!=pos)
  {
    c_io_printf(unit,
                "WARNING: getpos (with wordlength %" PRId64 ") fails,",
                word_len);
    c_io_printf(unit,
                "WARNING: because the byte position in the file (%" PRId64 ")",
                pos);
    c_io_printf(unit,
                "WARNING: is not a multiple of the word length");
    return ERR_MISSALIGN;
  }

  return word_index;
}

/* Switch the file mode between read/writing/idling
 *    - this is needed for buffer consistency
 */

int64_t c_io_change_mode(int64_t unit,int64_t newMode)
{
  int64_t r;

  /* change mode may cause a cached disk update to be made,
   * so ensure we are in the correct disk position
   */
  c_io_effect_setpos(unit);

  r=C_IO_CHANGE_MODE(unit,newMode);
  c_io_attributes[(size_t)(unit)].mode=newMode;

  return r;
}

/* Synchronise a unit to disk, flushing data as far as possible
 *    - ideally out of the local filesystem cache
 */

int64_t c_io_sync(int64_t unit)
{
  int64_t sync_return;

  sync_return=C_IO_SYNC(unit);

  if(sync_return!=c_io_success)
  {
    c_io_printf(unit,"c_io_sync: ERROR in call to C_IO_SYNC()");
    return c_io_fail;
  }

  return c_io_success;
}

/* Determine the size of whatever is attached to the unit */
int64_t c_io_getExtent(int64_t unit)
{
  int64_t currentPosition;
  int64_t endPosition;

  /* save the current position */
  currentPosition=c_io_getpos(unit,WORD8BYTES);

  /* move to the end */
  c_io_setpos(unit,c_io_file_end,WORD8BYTES);

  /* save the end position */
  endPosition=c_io_getpos(unit,WORD8BYTES);

  /* restore the original position */
  c_io_setpos(unit,currentPosition,WORD8BYTES);

  return endPosition;
}

void c_io_print_stack()
{
  c_io_layer_t *layer;
  int i=1;

  layer=params->layers;

  while(layer!=NULL)
  {
    if (layer->layer_number!=CIO_NO_LAYER)
    {
      if (layer->name==NULL)
      {
        if (layer->helper_start!=NULL && layer->helper_stop!=NULL)
        {
          c_io_printf(-1,"layer %2d @ 0x%" PRIuPTR " ID:%2" PRId64
                      " (helper thread capable)",
                      i,(uintptr_t)(layer),layer->layer_number);
        }
        else
        {
          c_io_printf(-1,"layer %2d @ 0x%" PRIuPTR " ID:%2" PRId64,
                      i,(uintptr_t)(layer),layer->layer_number);
        }
      }
      else
      {
        if (layer->helper_start!=NULL && layer->helper_stop!=NULL)
        {
          c_io_printf(-1,"layer %2d @ 0x%" PRIuPTR " ID:%2" PRId64 " name: %s"
                      " (helper thread capable)",
                      i,(uintptr_t)(layer),layer->layer_number,layer->name);
        }
        else
        {
          c_io_printf(-1,"layer %2d @ 0x%" PRIuPTR " ID:%2" PRId64 " name: %s",
                      i,(uintptr_t)(layer),layer->layer_number,layer->name);
        }
      }

      i++;
    }

    layer=layer->next;
  }
}

/* Change the messaging level */

void c_io_setOutPutLevel(int64_t in){
  c_io_outLevel=in;
}

c_io_layer_t *c_io_get_named_layer(int64_t layer_id)
{
  c_io_layer_t *layer=params->layers;

  while( layer->next != NULL)
  {
    if (layer->layer_number == layer_id)
    {
      return layer;
    }
    layer=layer->next;
  }

  c_io_printf(-1,"Info: failed to locate layer %2" PRId64,layer_id);

  return NULL;
}

c_io_layer_t *c_io_get_bottom_layer()
{
  c_io_layer_t *layer;

  if (params->layers==NULL)
  {
    params->layers=malloc(sizeof(c_io_layer_t));
    c_io_init_layer(params->layers,NULL,CIO_NO_LAYER);
    c_io_printf(-1,
                "Info: Allocated top layer at 0x%" PRIuPTR ,
                (uintptr_t)(params->layers));
  }

  layer=params->layers;

  while( layer->next != NULL)
  {
    layer=layer->next;
  }

  return layer;
}

void c_io_init_layer(c_io_layer_t *layer,c_io_layer_t *parent, int64_t id)
{
    size_t name_len=strlen(c_io_layer_names[id])+1;

    layer->parent=parent;
    layer->next=NULL;
    layer->helper_start=NULL;
    layer->helper_stop=NULL;
    layer->helper_tid=0;
    layer->helper_state=c_io_none;
    layer->layer_number=id;
    layer->name=malloc(name_len);
    strncpy(layer->name,c_io_layer_names[id],name_len);
}

c_io_layer_t *c_io_init_next_layer(int64_t next_layer_id)
{
  c_io_layer_t *this_layer=c_io_get_bottom_layer();

  this_layer->next=malloc(sizeof(c_io_layer_t));
  c_io_init_layer(this_layer->next,this_layer,next_layer_id);

  return this_layer;
}


void c_io_attachHelper(int64_t role){
  c_io_layer_t *targetLayer;
  c_io_printf(-1,"Info: Attach helper: role %" PRId64,role);

  switch(role)
  {
    case c_io_readHelper:
      targetLayer=c_io_get_named_layer(CIO_RBUFFERING);
      break;

    case c_io_writeHelper:
      targetLayer=c_io_get_named_layer(CIO_WBUFFERING);
      break;

    default:
      c_io_printf(-1,"Info: Attach helper: role %" PRId64 " is unknown",role);
      return;
  }

  if (targetLayer==NULL){
    c_io_printf(-1,
                "Info: Attach helper: "
                "layer for role %" PRId64
                " is not defined (not present in build?)",
                role);
    return;
  }

  if (targetLayer->helper_start==NULL){
    c_io_printf
      (-1,
       "Info: Attach helper: layer has no helper start function defined");
    return;
  }

  c_io_printf(-1,
              "Info: Attach helper: Target layer: 0x%" PRIuPTR
              " helper start: 0x%" PRIuPTR ,
              (uintptr_t)targetLayer,
              (uintptr_t)targetLayer->helper_start);

  (*targetLayer->helper_start)();
}

void c_io_detachAllHelpers(){
  int64_t i;
  for (i=1;i<numRoles;i++)
    c_io_detachHelper(i);
}

void c_io_detachHelper(int64_t role){
  c_io_layer_t *targetLayer;

  c_io_printf(-1,
              "Info: Detach helper: role: %" PRId64,role);

  switch(role)
  {
    case c_io_readHelper:
      targetLayer=c_io_get_named_layer(CIO_RBUFFERING);
      break;

    case c_io_writeHelper:
      targetLayer=c_io_get_named_layer(CIO_WBUFFERING);
      break;

    default:
      targetLayer=NULL;
      c_io_printf(-1,
                  "Info: Detach helper: role %" PRId64 " is unknown",
                  role);
      return;
  }

  if (targetLayer==NULL){
    c_io_printf(-1,
                "Info: Detach helper: "
                "layer for role %" PRId64
                " is not defined (not present in build?)",
                role);
    return;
  }

  if (targetLayer->helper_stop==NULL){
    c_io_printf(-1,
                "Info: Detach helper: "
                "layer ( @ 0x%" PRIuPTR ") has no helper stop function defined",
                (uintptr_t)(targetLayer));
    return;
  }

  c_io_printf(-1,
              "Info: Detach helper: Target layer: 0x%" PRIuPTR
              " helper stop: 0x%" PRIuPTR,
              (uintptr_t)(targetLayer),
              (uintptr_t)(targetLayer->helper_stop));
  (*targetLayer->helper_stop)();

}

double c_io_getTime()
{
  double t;
  struct timeval tv;
  gettimeofday( &tv,  NULL );
  t=(double)tv.tv_sec;
  t+=(double)(tv.tv_usec)/1000000.0;
  return t;
}

void c_io_allocate_query(query_type_t query, void **in_pack, void **out_pack)
{
  switch(query)
  {
    case qry_file_exist:
    {
      in_pack_fname_t *in_ptr = malloc(sizeof(in_pack_fname_t));
      *in_pack = in_ptr;
      in_ptr->filename = malloc(MAX_OUTSTR*sizeof(char));
      out_pack_bool_t *out_ptr = malloc(sizeof(out_pack_bool_t));
      *out_pack = out_ptr;
    }
    break;
  }
}

void c_io_free_query(query_type_t query, void *in_pack, void *out_pack)
{
  switch(query)
  {
    case qry_file_exist:
      free(((in_pack_fname_t *)in_pack)->filename);
      free(in_pack);
      free(out_pack);
      break;
  }
}

#endif

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

#include <time.h>
#include <errno.h>
#include "c_fort2c_prototypes.h"
#include "c_fort2c_types.h"
#include "pio_umprint.h"
#include "c_io_internal.h"
#include "c_io.h"

static const int64_t default_buffering=512*1024;
static const bool    default_timing=false;

/* DEPENDS ON: c_io.o */

/* default buffering size - 4MB (expressed in 64-bit words for */
/*                               consistency with 2A code)     */

/****************************************************************
 * IO wrappers to new interface                                 *
 ****************************************************************/

/****************************************************************
 * Open                                                         *
 ****************************************************************/

void open_single
(
 int64_t *unit,             /* Fortran unit                                  */
 void    *name,             /* File name or environ. var. (void * to char[]) */
 int64_t *char_len,         /* No of chars in file name                      */
 int64_t *intent,           /* =0 read only,!=0 read and write               */
 int64_t *environ_var_flag, /* =0 file name in environment var,              */
                            /*!=0 explicit file name                         */
 int64_t *err               /* =0 file opened,                               */
                            /* =1 file failed to open because the expected   */
                            /*    environment variable is not set,           */
                            /* =2 file not opened for another reason         */
 )
{
  int64_t mode;
  char *C_name;
  char *fileName;

  *err=0;

  /*c_io_ini
   * If the portioInit was not called (usually via ioInit in fortran)
   *  then initialise automatically with sane defaults
   */

  if (!c_io_initialised)
    c_io_init(8,300,
              default_buffering*(int64_t)sizeof(um_data_t),
              0,0,false,0,
              default_timing);

  /* get the filename */
  C_name=malloc((size_t)(*char_len)+1);
  strncpy(C_name,name,(size_t)(*char_len));
  C_name[ (size_t)(*char_len) ] = '\0';
  sscanf( C_name, "%s", C_name );

  if ( *environ_var_flag == 0 )  /* File name held in environ var */
  {
    fileName = getenv( C_name );
    if ( fileName == NULL )
    {
      c_io_printf(*unit,
                  "OPEN:  WARNING: Environment variable %s not set",
                  C_name);
      *err=1;
      free (C_name);
      return;
    }
  }
  else                           /* get file name from argmt C_name */
  {
    fileName = C_name;
  }

  mode=readwrite;
  if (*intent==0)mode=readonly;

  if(c_io_open(*unit,fileName,(int64_t)(strlen(fileName)),mode)!=c_io_success)
  {
        c_io_printf(*unit,"OPEN:  WARNING: failed to open file %s", C_name);
        *err=2;
  }
  free (C_name);
}

/****************************************************************
 * Close                                                        *
 ****************************************************************/

void close_single
(
 int64_t *unit,             /* Fortran unit                                  */
 void    *name,             /* File name or environ. var. (void * to char[]) */
 int64_t *char_len,         /* No of chars in file name                      */
 int64_t *environ_var_flag, /* =0 file name in environment var,              */
                            /*!=0 explicit file name                         */
 int64_t *delete,           /* =0 do not delete file,!=0 delete file         */
 int64_t *err               /* =0 file close ok,                             */
                            /*    or problems                                */
)
{
  char *C_name;
  char *fileName;
  char message[MAX_OUTSTR];
  size_t i;
  *err=0;

  /* convert file name to C format */
  C_name = malloc((size_t)(*char_len) + 1);
  strncpy(C_name,name,(size_t)(*char_len));
  C_name[(size_t)(*char_len)] = '\0';

  for (i=0; i<(size_t)(*char_len); i++){
    if (C_name[i] == ' '){
      C_name[i] = '\0';
      break;
    }
  }

  if(*environ_var_flag == 0)
  {
    fileName = getenv( C_name );
    if ( fileName == NULL ) {
      snprintf(message,MAX_OUTSTR,
               "CLOSE: WARNING: Environment variable %s not set",
               C_name);
      umPrint(message,"portio2b");
      *err=1;
      free( C_name );
      return;
    }
  }
  else
    fileName=C_name;

  if (c_io_close(*unit)!=c_io_success)*err=1;

  /* delete file */
  if(*delete != 0)
    c_io_delete(fileName);

  free( C_name );

}

/****************************************************************
 * buffin (and type varients)                                   *
 ****************************************************************/

void buffin8_single
(
 integer *unit,    /* Fortran unit                                      */
 void *array,      /* Array into which data is read (void * to char[] ) */
 integer *maxlen,  /* Number of real numbers to be read                 */
 integer *length,  /* Number of real numbers actualy read               */
 real *status      /* Return code                                       */
)
{
  *length=(integer)c_io_in(*unit,array,* maxlen,1);
  *status=-1.0;
  if (*length!=*maxlen)*status=1.0;
}

void buffin16_single
(integer *unit,    /* Fortran unit                         */
 void *array,      /* Array into which data is read        */
 integer *maxlen,  /* Number of real numbers to be read    */
 integer *length,  /* Number of real numbers actualy read  */
 real *status      /* Return code                          */
 )
{
  *length=(integer)c_io_in(*unit,array,* maxlen,2);
  *status=-1.0;
  if (*length!=*maxlen)*status=1.0;
}

void buffin32_single
(
 integer *unit,    /* Fortran unit                         */
 void *array,      /* Array into which data is read        */
 integer *maxlen,  /* Number of real numbers to be read    */
 integer *length,  /* Number of real numbers actualy read  */
 real *status      /* Return code                          */
)
{
  *length=(integer)c_io_in(*unit,array,* maxlen,WORD32BYTES);
  *status=-1.0;
  if (*length!=*maxlen)*status=1.0;
}

void buffin64_single
(
 integer *unit,    /* Fortran unit                         */
 void *array,      /* Array into which data is read        */
 integer *maxlen,  /* Number of real numbers to be read    */
 integer *length,  /* Number of real numbers actualy read  */
 real *status      /* Return code                          */
)
{
  *length=(integer)c_io_in(*unit,array,* maxlen,WORD64BYTES);
  *status=-1.0;
  if (*length!=*maxlen)*status=1.0;
}


/****************************************************************
 * buffout (and type varients)                                  *
 ****************************************************************/

void buffout8_single
(integer *unit,    /* Fortran unit                           */
 void *array,      /* Array from which data is written       */
 integer *maxlen,  /* Number of real numbers to write        */
 integer *length,  /* Number of real numbers actualy written */
 real *status      /* Return code                            */
 )
{
  *length=(integer)c_io_out(*unit,array,* maxlen,WORD8BYTES);
  *status=-1.0;
  if (*length!=*maxlen)*status=1.0;
}

void buffout32_single
(integer *unit,    /* Fortran unit                           */
 void *array,      /* Array from which data is written       */
 integer *maxlen,  /* Number of real numbers to write        */
 integer *length,  /* Number of real numbers actualy written */
 real *status      /* Return code                            */
 )
{
  *length=(integer)c_io_out(*unit,array,* maxlen,WORD32BYTES);
  *status=-1.0;
  if (*length!=*maxlen)*status=1.0;
}

void buffout64_single
(
 integer *unit,    /* Fortran unit                           */
 void *array,      /* Array from which data is written       */
 integer *maxlen,  /* Number of real numbers to write        */
 integer *length,  /* Number of real numbers actualy written */
 real *status      /* Return code                            */
)
{
  *length=(integer)c_io_out(*unit,array,* maxlen,WORD64BYTES);
  *status=-1.0;
  if (*length!=*maxlen)*status=1.0;
}

/****************************************************************
 * setpos (and type varients)                                   *
 ****************************************************************/

void setpos_single
(int64_t *unit,          /* Fortran unit                         */
 int64_t *word_address,  /* Number of words into file            */
 int64_t *err            /* Error checking: err = 0 no errors,   */
                         /*                 err = 1 errors       */
 )
{
  c_io_setpos(*unit,*word_address,sizeof(integer));
  *err=0;
}

void setpos8_single
(
 int64_t *unit,          /* Fortran unit                         */
 int64_t *word_address,  /* Number of words into file            */
 int64_t *err            /* Error checking, err = 0 no errors,   */
                         /*                 err = 1 errors       */
 )
{
  c_io_setpos(*unit,*word_address,WORD8BYTES);
  *err=0;
}

void setpos32_single(
 int64_t *unit,          /* Fortran unit                         */
 int64_t *word_address,  /* Number of words into file            */
 int64_t *err            /* Error checking: err = 0 no errors,   */
                         /*                 err = 1 errors       */
 )
{
  c_io_setpos(*unit,*word_address,WORD32BYTES);
  *err=0;
}

void setpos64_single(
              int64_t *unit,         /* Fortran unit                          */
              int64_t *word_address, /* Number of words into file             */
              int64_t *err           /* Error checking, err = 0 no errors,    */
                                     /*                 err = 1 errors        */
             )
{
  c_io_setpos(*unit,*word_address,WORD64BYTES);
  *err=0;
}


/****************************************************************
 * getpos                                                       *
 ****************************************************************/

void getpos
(int64_t *unit,          /* Fortran unit                         */
 int64_t *word_address   /* Number of words into file            */
 )
{
  (*word_address)=c_io_getpos(*unit,sizeof(integer));
}

void getpos8
(int64_t *unit,          /* Fortran unit                         */
 int64_t *word_address   /* Number of words into file            */
 )
{
  (*word_address)=c_io_getpos(*unit,WORD8BYTES);
}

void getpos32
(int64_t *unit,          /* Fortran unit                         */
 int64_t *word_address   /* Number of words into file            */
 )
{
  (*word_address)=c_io_getpos(*unit,WORD32BYTES);
}

void getpos64

(int64_t *unit,          /* Fortran unit                         */
 int64_t *word_address   /* Number of words into file            */
 )
{
  (*word_address)=c_io_getpos(*unit,WORD64BYTES);
}

void sync_single (integer *unit,integer *icode)
{
  c_io_sync(*unit);
  *icode=0;
}

/****************************************************************/
/* The following function determines the size in bytes of the   */
/* file attached to a given unit.                               */
/****************************************************************/
void getextent(integer *unit,
               integer *extent,
               integer *icode)
{
  *icode=0;
  *extent=(integer)c_io_getExtent(*unit);
}

void set_printstatus(integer *printstatus_f) /* IN: fortran equivalent */
{
  c_io_setOutPutLevel(*printstatus_f);
}


void portioshutdown(void)
{
  c_io_fini();
}

/****************************************************************/
/* The following function is to be called from Fortran.         */
/* It allows the user to override the size of the I/O buffer    */
/* as passed to c_io_init                                       */
/****************************************************************/
void portioinit(int64_t *wb_size,
                int64_t *rb_size,
                int64_t *rb_count,
                bool    *rb_update,
                int64_t *rb_prefetch,
                bool    *io_timing)
{
  c_io_init(8,
            300,
            (*wb_size * (int64_t)(sizeof(um_data_t))),
            (*rb_size * (int64_t)(sizeof(um_data_t))),
            *rb_count,
            *rb_update,
            *rb_prefetch,
            *io_timing);
}

void portioaddhelper(integer *role)
{
  c_io_attachHelper(*role);
}

void portiodetachallhelpers(void)
{
  c_io_detachAllHelpers();
}

/******************************************************************************/
/* is_unit_open_local() determines the status (open or closed) of a unit.     */
/* Returned values are boolean "ret_code":                                    */
/*  - True  - unit is open                                                    */
/*  - False - unit is closed                                                  */
/******************************************************************************/
void is_unit_open_local(int64_t *unit, bool *ret_code)
{
  /* If the data structures aren't allocated, nothing is open */
  if (c_io_attributes == NULL || params == NULL)
  {
    *ret_code = false;
    return;
  }

  /* If the unit is out of range, its not open */
  if ( (size_t)*unit < params->min_unit
       || (size_t)*unit > params->max_unit
       || *unit < 0
       || (uintmax_t)*unit > (uintmax_t)SIZE_MAX
     )
  {
    *ret_code = false;
    return;
  }

  /* The unit is open if its file handle is assigned */
  *ret_code = (c_io_attributes[*unit].fileHandle!=NULL);
}

#else
extern int unused_unit;
#endif

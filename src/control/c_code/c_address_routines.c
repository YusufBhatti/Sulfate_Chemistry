/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#include <stddef.h>
#include <stdint.h>

/* -------------------------------------------------------------------------- */
/* Prototypes                                                                 */
/* -------------------------------------------------------------------------- */

size_t   um_addr_size  (void);
void     um_addr_diff  (void *object1, void *object2, intptr_t *diff);
intptr_t um_addr_value (void *object1);

/* -------------------------------------------------------------------------- */
/* Procedures                                                                 */
/* -------------------------------------------------------------------------- */

/****************************************************************/
/* The following function determines the difference in bytes of */
/* two (near) addresses.                                        */
/****************************************************************/
void um_addr_diff(void *object1, void *object2, intptr_t *diff)
{
  *diff=(intptr_t)object2-(intptr_t)object1;
}

/****************************************************************/
/* The following function determines the size of pointers       */
/****************************************************************/
size_t um_addr_size(void)
{
  return sizeof(void *);
}

/****************************************************************/
/* The following function returns the integer value of the      */
/* address of the given pointer                                 */
/****************************************************************/
intptr_t um_addr_value(void *object1)
{
  return (intptr_t)object1;
}

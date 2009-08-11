/*
 *  mmmalloc.h
 *  
 *  header for memory mapped graph storage
 *  Created by Andrew Byrd on 11/08/09. 
 *  Include in all graphserver core c files
 *  but not in mmmalloc.c as these defines
 *  will block access to the 'real' malloc and free!
 *
 */

#include <stdlib.h>

int   mmon     ( );
int   mmoff    ( );
void* mmmalloc (size_t size);
void* mmcalloc (size_t elements, size_t size);
void  mmfree   (void * p);

#define malloc mmmalloc
#define calloc mmcalloc
#define free   mmfree

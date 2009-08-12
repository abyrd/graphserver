/*
 *  mmmalloc.c
 *  Adds memory mapped file storage for Graphs
 *
 *  Created by Andrew Byrd on 11/08/09.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <string.h>

// do not include mmalloc.h here since it redefines malloc functions.
// the calls to regular malloc in here would not work.

//  2MiB, the address that mmap always gives me automatically on OSX
#define OSX_BASE (void*) 0x00200000
// 16MiB, the address that mmap gives me for bigger mmaps on OSX
// #define OSX_BASE_BIG (void*) 0x01000000
// looks like some of this space is already taken by Python, so move higher
#define OSX_BASE_BIG (void*) 0x02000000
// upper limit of a 32 bit address space
#define MAX_ADDR (void*) 0xFFFFFFFF

//  1GiB dynamic memory area
#define ALLOC_SPACE 0x40000000
//  4KiB is the page size
#define ALLOC_PAGE  0x0000FFFF
// 16MiB expansion chunk size
#define ALLOC_CHUNK 0x01000000

static void * mmmalloc_base = NULL;
static void * mmmalloc_top  = NULL;
static int mmmalloc_fd;
static int mmmalloc_on = 0;
static size_t mmmalloc_chunks = 0;

int mmgrab ( ) {
    // allow only one mmaped graph per session?
    // block needs structure:
    // base address, top pointer, a graph struct always at a certain address.
    mmmalloc_fd = open("mmmalloc.mmap", O_CREAT | O_RDWR , 0664);
    ftruncate( mmmalloc_fd , ALLOC_CHUNK );
    mmmalloc_chunks = 1;
    mmmalloc_base = mmap( OSX_BASE_BIG, ALLOC_SPACE, PROT_READ | PROT_WRITE, 
                          MAP_SHARED, mmmalloc_fd, 0);
    if ( mmmalloc_base != OSX_BASE_BIG ) { 
        printf( "Unable to obtain requested block.\n"); 
        printf( "I got 0x%08x instead.\n", (uint) mmmalloc_base ); 
    } else {
        printf( "mmmalloc got 0x%08x bytes at 0x%08x.\n", ALLOC_SPACE, (uint) mmmalloc_base );
    }
    mmmalloc_top = mmmalloc_base;    
    return 0;
}

int mmon ( ) {
    mmmalloc_on = 1; 
    return 0;
}

int mmoff ( ) {
    mmmalloc_on = 0;
    return 0;
}

int mmrelease ( ) {
    munmap(mmmalloc_base, mmmalloc_top - mmmalloc_base);
    close(mmmalloc_fd);
    return 0;
}

void * mmmalloc(size_t size) {
    if (mmmalloc_on) {
        if (!mmmalloc_base) mmgrab();
        // printf("mm ");
        void * ret = mmmalloc_top;
        if (mmmalloc_top >= (void *) (mmmalloc_chunks * ALLOC_CHUNK)) {
            mmmalloc_chunks ++;
            ftruncate( mmmalloc_fd , ALLOC_CHUNK * mmmalloc_chunks );
            printf ("\nmmmalloc : expanded to 0x%08x bytes.\n", (uint) mmmalloc_chunks * ALLOC_CHUNK);
        }
        mmmalloc_top += size;
        // align pointers to avoid errors
        // intel and OSX requires 16 byte alignment
        mmmalloc_top += ((16 - ((long) mmmalloc_top % 16)) % 16);
        return ret; 
    } else { 
        // mmap is turned off, do a normal malloc
        // printf("sm ");
        return malloc( size );
    }
    return NULL;
}

void * mmcalloc(size_t elements, size_t size) {
    if (mmmalloc_on) {
        // printf("mc ");
        size_t real_size = elements * size;
        memset(mmmalloc_top, 0, real_size); 
        return mmmalloc( real_size );
    } else { 
        // mmap is turned off, do a normal calloc()
        // printf("sc ");
        return calloc( elements, size );
    }
    return NULL;
}


void mmfree(void * p) {
    if (mmmalloc_on) {
        // printf("mf ");
        // Graphs are grow-only in memmap mode
        // a lot of freeing goes on when building a graph
    } else {
        // mmap is turned off, do a normal free()
        // printf("sf ");
        free( p );
    }
}

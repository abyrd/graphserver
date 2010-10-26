/* 

Headers for:
Binary heap.
Is the fibonacci heap really all that?

*/

#ifndef _FIBHEAP_H_
#define _FIBHEAP_H_

#include "stdio.h"

// keys are long ints
typedef long fibheapkey_t;

// a heap (pointer)
typedef struct fibheap
{
  size_t capacity;
  size_t size;
  struct fibnode *elements;
} *fibheap_t;

// a node in the heap: key and data pointer
typedef struct fibnode
{
  fibheapkey_t key;
  void *data;
} *fibnode_t;

extern fibheap_t fibheap_new (void);
extern fibnode_t fibheap_insert (fibheap_t, fibheapkey_t, void *);
extern int fibheap_empty (fibheap_t);
extern void *fibheap_extract_min (fibheap_t);
extern void *fibheap_delete_node (fibheap_t, fibnode_t);
extern void fibheap_delete (fibheap_t);

#endif /* _FIBHEAP_H_ */

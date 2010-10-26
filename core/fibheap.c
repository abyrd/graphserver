/* 

Binary heap for priority queue.
http://cprogramminglanguage.net/binary-heap-c-code.aspx

*/

#include "fibheap.h"
#include "graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#define FIBHEAPKEY_MIN	LONG_MIN

/* Print an error message and exit.  */
void die (char *message) {
    printf("%s", message);
    exit(-1);
}

void fibheap_dump (fibheap_t heap) {
    int i;    
    printf("---\n");
    for (i=1; i<=heap->size; i++) {
        printf("%ld - %ld\n", heap->elements[i].key, (long)heap->elements[i].data);
    }
}

/* Create a new heap.  */
fibheap_t fibheap_new (void) {
  size_t N = 1000000;
  fibheap_t heap = malloc( sizeof(struct fibheap) );
  if (heap == NULL) die("Cannot allocate priority heap.");
  heap->capacity = N;
  heap->size = 0;
  // indexing starts at 1
  heap->elements = malloc((N+1) * sizeof(struct fibnode));
  if (heap->elements == NULL) die("Cannot allocate array to hold priority heap elements.");
  heap->elements->key = FIBHEAPKEY_MIN;
  return heap;
}

/* Delete HEAP.  */
void fibheap_delete (fibheap_t heap) {
  free (heap->elements);
  free (heap);
}

/* Determine if HEAP is empty.  */
int fibheap_empty (fibheap_t heap) {
  return heap->size == 0;
}

/* Insert DATA, with priority KEY, into HEAP.  */
fibnode_t fibheap_insert (fibheap_t heap, fibheapkey_t key, void *data) {
  int i;
  // check for overflow
  if (heap->size >= heap->capacity) 
    die("Priority heap array overflow (too many items in queue).");
  // starting at the bottom of the heap, swap nodes down 
  for (i = ++(heap->size); heap->elements[i/2].key > key; i/=2) {
        heap->elements[i].key  = heap->elements[i/2].key; 
        heap->elements[i].data = heap->elements[i/2].data;
  }
  // put the new data in its place
  heap->elements[i].key  = key;
  heap->elements[i].data = data;
  //fibheap_dump(heap);
  return &(heap->elements[i]);
}

/* Extract the data of the minimum node from HEAP.  */
void *fibheap_extract_min (fibheap_t heap) {
   int i, child;
   struct fibnode minElement, lastElement;

   if (fibheap_empty(heap)) die("Priority queue is empty.");
   minElement  = heap->elements[1];
   lastElement = heap->elements[heap->size];
   (heap->size)--;

   for (i=1; i*2 <= heap->size; i=child) {
       /* Find smaller child */
       child = i*2;
       if (child != heap->size && heap->elements[child+1].key < heap->elements[child].key)
           child++;

       /* Percolate one level */
       if (lastElement.key > heap->elements[child].key) {
          heap->elements[i].key  = heap->elements[child].key;
          heap->elements[i].data = heap->elements[child].data;
       } else break;
   }
   heap->elements[i].key  = lastElement.key;
   heap->elements[i].data = lastElement.data;
   // printf("extract key: %ld\n", minElement.key);
   return minElement.data;
}

/* Delete NODE from HEAP.  */
void *
fibheap_delete_node (fibheap_t heap, fibnode_t node)
{
    return;
}

int selftest(void) {
    // fill and empty
    int    N;
    int    i, r; 
    long   m, n;

    for (N=900000; N<1000000; N+=10000) {
        fibheap_t heap = fibheap_new();
        clock_t t0 = clock(), dt;
        for (i=0; i<N; i++) {
            r = rand();
            fibheap_insert(heap, r, (void*)(r/2));
        }
        m = 0;
        while (!fibheap_empty(heap)) {
            n = (long)fibheap_extract_min(heap);                
            assert(n >= m);
            m = n;
        }
        dt = clock() - t0;
        float msec = dt * 1000 / CLOCKS_PER_SEC;
        printf("N=%d / etime %f / avg %f msec\n", N, msec, msec/N);
        fibheap_delete(heap);
    }
}

int selftest2(void) {
    // fill and maintain full
    int    M, N;
    int    i, r; 
    long   m, n;

    M = 800000;
    N = 8000000;
    fibheap_t heap = fibheap_new();
    for (i=0; i<M; i++) {
        r = rand();
        fibheap_insert(heap, r, (void*)(r/2));
    }
    m = 0;
    clock_t t0 = clock(), dt;
    for (i=0; i<N; i++) {
        r = rand();
        fibheap_insert(heap, r, (void*)(r/2));
        n = (long)fibheap_extract_min(heap);                
        assert(n >= m || n == r/2);
        m = n;
    }
    dt = clock() - t0;
    float msec = dt * 1000 / CLOCKS_PER_SEC;
    printf("N=%d / etime %f / avg %f msec\n", N, msec, msec/N);
    fibheap_delete(heap);
}

int main () {
    selftest2();
    return 0;
}


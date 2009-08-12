Graph*
#ifndef RETRO
gShortestPathTree( Graph* this, char *from, char *to, State* init_state, WalkOptions* options, long maxtime ) {
#else
gShortestPathTreeRetro( Graph* this, char *from, char *to, State* init_state, WalkOptions* options, long mintime ) {
#endif
    
/*
 *  VARIABLE SETUP
 */
  //Iteration Variables
  Vertex *u, *v;
  Vertex *spt_u, *spt_v;
  State *du, *dv, *dv_time;
  int count = 1;

  //Goal Variables
#ifndef RETRO
  char* origin = from;
  char* target = to;
#else
  char* origin = to;
  char* target = from;
#endif

  // search cutoffs
  long time_cutoff = init_state->time + 60 * 60 * 4;
  int  walk_cutoff = 1600;
  int  transfer_cutoff = 3;   

  //Get origin vertex to make sure it exists
  Vertex* origin_v = gGetVertex( this, origin );
  if( origin_v == NULL ) {
    return NULL;
  }
    
  //Return Tree
  Graph*  spt = gNew();
  Vertex* origin_v_spt = gAddVertex( spt, origin );
  origin_v_spt->payload = init_state;
  origin_v_spt->payload_time = stateDup( init_state ); // cannot have more than one reference to a State
  
  //Priority Queues
  dirfibheap_t q        = dirfibheap_new( gSize( this ) );
  dirfibheap_t q_time   = dirfibheap_new( gSize( this ) );
  dirfibheap_insert_or_dec_key( q,      origin_v, 0 );
  dirfibheap_insert_or_dec_key( q_time, origin_v, 0 );

/*
 *  CENTRAL ITERATION
 *
 */

  while( 1 ) {                                
     if ( !dirfibheap_empty(q_time) ) {          // check time queue first
        u = dirfibheap_extract_min( q_time );    // get the earliest Vertex 'u',
        spt_u = gGetVertex( spt, u->label );     // get corresponding SPT Vertex,  
        du = (State*)spt_u->payload_time;        // and get State of u 'du'.
        // printf( "got one vertex from time queue: %s \n", u->label );
     } else if ( !dirfibheap_empty( q ) ) {        // if time queue is empty, check weight queue
        u = dirfibheap_extract_min( q );         // get the lowest-weight Vertex 'u',
        spt_u = gGetVertex( spt, u->label );     // get corresponding SPT Vertex,  
        du = (State*)spt_u->payload;             // and get State of u 'du'.
        // printf( "got one vertex from weight queue: %s \n", u->label );
     } else {
        break;  // all queues are empty - terminate.
     }

//    if( !strcmp( u->label, target ) )                //(end search if reached destination vertex)
//     break;                                          // let's not worry about termination for now.

// do not continue expanding a path if it is unreasonably long
#ifndef RETRO
      if( du->time > time_cutoff || du->dist_walked > walk_cutoff || du->num_transfers > transfer_cutoff )
         continue;
#else
// this is wrong.
      if( du->time < time_cutoff || du->dist_walked > walk_cutoff || du->num_transfers > transfer_cutoff )
         continue;
#endif

#ifndef RETRO
    ListNode* edges = vGetOutgoingEdgeList( u );
#else
    ListNode* edges = vGetIncomingEdgeList( u );
#endif

    // if (edges) printf ("edge list found. \n");
    // else printf ("edge list is NULL. \n");
    while( edges ) {                                 //For each Edge 'edge' connecting u
      Edge* edge = edges->data;
#ifndef RETRO
      v = edge->to;                                  //to Vertex v:
#else
      v = edge->from;
#endif

      long old_w;
      long old_t;
      // printf ("processing edge to %s. \n", v->label);
      if( (spt_v = gGetVertex( spt, v->label )) ) {   //get the SPT Vertex corresponding to 'v'
        dv = (State*)spt_v->payload;                  //and its best weight State 'dv'
        old_w = dv->weight;
        dv_time = (State*)spt_v->payload_time;        //and its best time State 'dv_time'
        old_t = dv_time->time;
        // printf( "found target vertex %s in spt. weight: %ld time: %ld \n", v->label, old_w, old_t );
      } else {
        dv      = NULL;                               //which may not exist yet
        old_w   = INFINITY;
        dv_time = NULL;                               //which may not exist yet
        old_t   = 2000000000; // apparently INFINITY wasn't infinite enough
        // printf( "vertex not yet in spt: %s \n", v->label );
      }

#ifndef RETRO
      State *new_dv = eWalk( edge, du, options );
#else
      State *new_dv = eWalkBack( edge, du, options );
#endif

      // When an edge leads nowhere (as indicated by returning NULL), the iteration is over.
      if(!new_dv) {
        // printf("Walk function returned NULL (%s(%ld) -> %s) \n",edge->from->label, du->weight, edge->to->label);
        edges = edges->next;
        continue;
      }

      // States cannot have weights lower than their parent State.
      if(new_dv->weight < du->weight) {
        fprintf(stderr, "Negative weight (%s(%ld) -> %s(%ld)) \n",edge->from->label, du->weight, edge->to->label, new_dv->weight);
        edges = edges->next;
        continue;
      }

      // should really be done selectively
      // need separate state objects to avoid needing reference counts on States
      State *new_dv_time = stateDup( new_dv );

      if( !spt_v ) {  // vertex is not yet in SPT
         // printf("Adding target vertex to the SPT. \n");
         spt_v = gAddVertex( spt, v->label );        //Copy v over to the SPT
         count++;
      } 
      long new_w = new_dv->weight;
      if( new_w < old_w ) {  // If the new way of getting there is better,
         // printf("Found better weight. Rekeying vertex in weight queue. \n");
         dirfibheap_insert_or_dec_key( q, v, new_w );    // rekey v in the priority queue
         if (spt_v -> payload) stateDestroy(spt_v->payload);
         spt_v->payload = new_dv;                        //Set the State of v in the SPT to the current winner
         // only weight path is recorded... this does not really work.
         vSetParent( spt_v, spt_u, edge->payload );      //Make u the parent of v in the SPT
      } else {
        stateDestroy(new_dv); //new_dv will never be used; merge it with the infinite.
      }
  
      long new_t = new_dv_time->time;
      // If the new way of getting there is FASTER
      if( new_t < old_t ) {
         // printf("Found better time. Rekeying vertex in time queue. \n");
         dirfibheap_insert_or_dec_key( q_time, v, new_t );    // rekey v in the times priority queue
         if (spt_v -> payload_time) stateDestroy(spt_v->payload_time);
         spt_v->payload_time = new_dv_time;               //Set the State of v (best time) in the SPT to the current winner
      } else {
         stateDestroy(new_dv_time); //new_dv_time will never be used; merge it with the infinite.
      }
      
      edges = edges->next;
    }    // end while (edges)
  }      // end while (1)

  dirfibheap_delete( q );
  dirfibheap_delete( q_time );

  //fprintf(stdout, "Final shortest path tree size: %d\n",count);
  return spt;
}

void
#ifndef RETRO
gShortestPathInPlace( Graph* this, char *from, char *to, State* init_state, WalkOptions* options, long maxtime ) {
#else
gShortestPathInPlaceRetro( Graph* this, char *from, char *to, State* init_state, WalkOptions* options, long mintime ) {
#endif
    
    // fprintf(stderr, "entered c function\n");
    
    // Iteration Variables
    Vertex *u,  *v;
    State  *du, *dv;
            
    // Goal Variables
#ifndef RETRO
    char* origin = from;
    char* target = to;
#else
    char* origin = to;
    char* target = from;
#endif
            
    // Get origin vertex to make sure it exists
    Vertex* origin_v = gGetVertex( this, origin );
    if( origin_v == NULL ) {
        return;
    }
            
    // Open vertex list, set all payloads to NULL & deallocate them
    struct hashtable_itr *itr = hashtable_iterator(this->vertices);
    // next_exists is false when number of vertices is 0
    int next_exists=hashtable_count(this->vertices); 
    while(itr && next_exists) {
        Vertex* vtx = hashtable_iterator_value( itr );
        // Deallocate State payload if exists
        if (vtx->payload) {
            stateDestroy(vtx->payload);
            vtx->payload = NULL;
        }
        next_exists = hashtable_iterator_advance( itr );
    }

    // Set the initial state
    origin_v->payload = init_state;

    // Priority Queue
    dirfibheap_t q = dirfibheap_new( gSize( this ) );
    dirfibheap_insert_or_dec_key( q, gGetVertex( this, origin ), 0 );
    
 
    // CENTRAL ITERATION
            
    while( !dirfibheap_empty( q ) ) {                   // Until the priority queue is empty:
        u = dirfibheap_extract_min( q );                // get the lowest-weight Vertex 'u',
        // fprintf(stderr, "next queue entry : %s\n", u->label);
        if( !strcmp( u->label, target ) )               // (end search if reached destination vertex)
            break;                
        
        du = (State*) u->payload;                       // and get State of u 'du'.

#ifndef RETRO
        if( du->time > maxtime )
            break;
        ListNode* edges = vGetOutgoingEdgeList( u );
#else
        if( du->time < mintime )
            break;
        ListNode* edges = vGetIncomingEdgeList( u );
#endif

        while( edges ) {                           // For each Edge 'edge' connecting u
            Edge* edge = edges->data;
            // fprintf(stderr, "(%s -> %s)\n",edge->from->label, edge->to->label);
            
#ifndef RETRO
            v = edge->to;                          // to Vertex v:
#else
            v = edge->from;
#endif
                    
            long old_w;
            
            dv = (State*) v->payload;              // get State of v 'dv'
            if (dv) {                               
                old_w = dv->weight;
            } else {
                old_w = INFINITY;                  // which may not exist yet
            }
                    
#ifndef RETRO
            State *new_dv = eWalk( edge, du, options );
#else
            State *new_dv = eWalkBack( edge, du, options );
#endif
                    
            // When an edge leads nowhere (as indicated by returning NULL), the iteration is over.
            if(!new_dv) {
                edges = edges->next;
                continue;
            }
                    
            // States cannot have weights lower than their parent State.
            if(new_dv->weight < du->weight) {
                fprintf(stderr, "Negative weight (%s(%ld) -> %s(%ld))\n",edge->from->label, du->weight, edge->to->label, new_dv->weight);
                edges = edges->next;
                continue;
            }
            // fprintf(stderr, "(%s(%ld) -> %s(%ld))\n",edge->from->label, du->weight, edge->to->label, new_dv->weight);
            
            long new_w = new_dv->weight;
            // If the new way of getting there is better,
            if( new_w < old_w ) {
                // fprintf(stderr, "better weight. rekeying %s\n", v-> label);
                dirfibheap_insert_or_dec_key( q, v, new_w );        // rekey v in the priority queue

                // If this is the first time v has been reached
                // if( !spt_v ) {
                //     spt_v = gAddVertex( spt, v->label );         //Copy v over to the SPT
                //     count++;
                // }
                        
                if (v->payload)
                    stateDestroy(v->payload);
                
                v->payload = new_dv;                    //Set the State of v to the current winner
            } else {
                stateDestroy(new_dv);                   //new_dv will never be used; merge it with the infinite.
            }
            edges = edges->next;
        }
    }
    dirfibheap_delete( q );
    // printf("End of c function.");        
    return;
}

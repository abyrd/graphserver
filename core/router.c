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
  State *du, *dv;
  int count = 1;

  //Goal Variables
#ifndef RETRO
  char* origin = from;
  char* target = to;
#else
  char* origin = to;
  char* target = from;
#endif

  //Get origin vertex to make sure it exists
  Vertex* origin_v = gGetVertex( this, origin );
  if( origin_v == NULL ) {
    return NULL;
  }
    
  //Return Tree
  Graph* spt = gNew();
  gAddVertex( spt, origin )->payload = init_state;
  //Priority Queue
  dirfibheap_t q_weight = dirfibheap_new( gSize( this ) );
  dirfibheap_insert_or_dec_key( q_weight, gGetVertex( this, origin ), 0 );
  // Another queue for best-time states
  dirfibheap_t q_time = dirfibheap_new( gSize( this ) );
  dirfibheap_t q_curr = NULL;

/*
 *  CENTRAL ITERATION
 *
 */

  while( 1 ) {                                
     if      ( !dirfibheap_empty( q_weight ) ) q_curr = q_weight;
     else if ( !dirfibheap_empty( q_time   ) ) q_curr = q_time;
     else break;  // all queues are empty - terminate.

     u = dirfibheap_extract_min( q_curr );            //get the lowest-weight Vertex 'u',
     spt_u = gGetVertex( spt, u->label );             //get corresponding SPT Vertex,  
     if      ( q_curr == q_weight ) du = (State*)spt_u->payload;                     //and get State of u 'du'.
     else if ( q_curr == q_time )   du = (State*)spt_u->payload_time;                //and get State of u 'du'.


//    if( !strcmp( u->label, target ) )                //(end search if reached destination vertex)
//     break;                                           // let's not worry about termination for now.


// add more termination conditions
#ifndef RETRO
    if( du->time > maxtime )
      break;
#else
    if( du->time < mintime )
      break;
#endif

#ifndef RETRO
    ListNode* edges = vGetOutgoingEdgeList( u );
#else
    ListNode* edges = vGetIncomingEdgeList( u );
#endif
    while( edges ) {                                 //For each Edge 'edge' connecting u
      Edge* edge = edges->data;
#ifndef RETRO
      v = edge->to;                                  //to Vertex v:
#else
      v = edge->from;
#endif

      long old_w;
      long old_t;
      if( (spt_v = gGetVertex( spt, v->label )) ) {        //get the SPT Vertex corresponding to 'v'
        dv = (State*)spt_v->payload;                     //and its State 'dv'
        old_w = dv->weight;
        old_t = dv->time;
      } else {
        dv = NULL;                                       //which may not exist yet
        old_w = INFINITY;
        old_t = INFINITY;
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

      // should really be done selectively
      // need separate state objects to avoid needing reference counts
      State *new_dv_time = stateDup( new_dv );

      long new_w = new_dv->weight;
      // If the new way of getting there is better,
      if( new_w < old_w ) {
        dirfibheap_insert_or_dec_key( q_weight, v, new_w );    // rekey v in the priority queue

        // If this is the first time v has been reached
        if( !spt_v ) {
          spt_v = gAddVertex( spt, v->label );        //Copy v over to the SPT
          count++;
          }

        //if((count%10000) == 0)
        //  fprintf(stdout, "Shortest path tree size: %d\n",count);

        if(spt_v->payload)
            stateDestroy(spt_v->payload);
        spt_v->payload = new_dv;                      //Set the State of v in the SPT to the current winner

        // only weight path is recorded... this does not really work.
        vSetParent( spt_v, spt_u, edge->payload );      //Make u the parent of v in the SPT
      } else {
        stateDestroy(new_dv); //new_dv will never be used; merge it with the infinite.
      }
      
      long new_t = new_dv->time;
      // If the new way of getting there is FASTER
      if( new_t < old_t ) {
        dirfibheap_insert_or_dec_key( q_time, v, new_t );    // rekey v in the times priority queue

        if(spt_v->payload_time)
            stateDestroy(spt_v->payload_time);
        spt_v->payload_time = new_dv_time;               //Set the State of v (best time) in the SPT to the current winner

      } else {
        stateDestroy(new_dv_time); //new_dv_time will never be used; merge it with the infinite.
      }

      edges = edges->next;
    }
  }

  dirfibheap_delete( q );
  dirfibheap_delete( q_time );

  //fprintf(stdout, "Final shortest path tree size: %d\n",count);
  return spt;
}

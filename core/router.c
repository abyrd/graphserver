ShortestPathTree*
#ifndef RETRO
gShortestPathTree( Graph* this, char *from, char *to, State* init_state, WalkOptions* options, long maxtime, int hoplimit, long weightlimit ) {
#else
gShortestPathTreeRetro( Graph* this, char *from, char *to, State* init_state, WalkOptions* options, long mintime, int hoplimit, long weightlimit ) {
#endif
    
/*
 *  VARIABLE SETUP
 */
  //Iteration Variables
  Vertex *u, *v;
  SPTVertex *spt_u, *spt_v;
  State *du, *dv;
  int count = 1;
  int iterations = 0;

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
  ShortestPathTree* spt = sptNew();
  sptAddVertex( spt, origin_v, 0 )->state = init_state;
  //Priority Queue
  dirfibheap_t q = dirfibheap_new( 1 );
  dirfibheap_insert_or_dec_key( q, gGetVertex( this, origin ), 0 );

/*
 *  CENTRAL ITERATION
 *
 */

  while( !dirfibheap_empty( q ) ) {                  //Until the priority queue is empty:
    u = dirfibheap_extract_min( q );                 //get the lowest-weight Vertex 'u',

      if( !strcmp( u->label, target ) ) {                //(end search if reached destination vertex)
          printf("total iterations: %d", iterations);
          break;
      }
      iterations += 1;
      
    spt_u = sptGetVertex( spt, u->label );             //get corresponding SPT Vertex,
    
    if( spt_u->hop >= hoplimit ) {
      break;
    }
    
    if( spt_u->state->weight > weightlimit ) {
      break;
    }
    
    du = (State*)spt_u->state;                     //and get State of u 'du'.
    
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
      if( (spt_v = sptGetVertex( spt, v->label )) ) {        //get the SPT Vertex corresponding to 'v'
        dv = (State*)spt_v->state;                     //and its State 'dv'
        old_w = dv->weight;
      } else {
        dv = NULL;                                       //which may not exist yet
        old_w = INFINITY;
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
	stateDestroy( new_dv );
        continue;
      }

      long new_w = new_dv->weight;
      // If the new way of getting there is better,
      if( new_w < old_w ) {
        dirfibheap_insert_or_dec_key( q, v, new_w );    // rekey v in the priority queue

        // If this is the first time v has been reached
        if( !spt_v ) {
          spt_v = sptAddVertex( spt, v, spt_u->hop+1 );        //Copy v over to the SPT
          count++;
          }

        //if((count%10000) == 0)
        //  fprintf(stdout, "Shortest path tree size: %d\n",count);

        if(spt_v->state)
            stateDestroy(spt_v->state);
        spt_v->state = new_dv;                      //Set the State of v in the SPT to the current winner

        sptvSetParent( spt_v, spt_u, edge->payload );      //Make u the parent of v in the SPT
      } else {
        stateDestroy(new_dv); //new_dv will never be used; merge it with the infinite.
      }
      edges = edges->next;
    }
  }

  dirfibheap_delete( q );

  //fprintf(stdout, "Final shortest path tree size: %d\n",count);
  return spt;
}

    
#ifndef RETRO    
// extra function to test spt assist
ShortestPathTree*
gShortestPathTreeAssist( Graph* this, ShortestPathTree* that, char *from, char *to, State* init_state, WalkOptions* options, long maxtime, int hoplimit, long weightlimit ) {
        
        /*
         *  VARIABLE SETUP
         */
        //Iteration Variables
        Vertex *u, *v;
        SPTVertex *spt_u, *spt_v;
        State *du, *dv;
        int count = 1;
        int iterations = 0;

        //Goal Variables
        char* origin = from;
        char* target = to;
        
        //Get origin vertex to make sure it exists
        Vertex* origin_v = gGetVertex( this, origin );
        if( origin_v == NULL ) {
            return NULL;
        }
        
        //Return Tree
        ShortestPathTree* spt = sptNew();
        sptAddVertex( spt, origin_v, 0 )->state = init_state;
        //Priority Queue
        dirfibheap_t q = dirfibheap_new( 1 );
        dirfibheap_insert_or_dec_key( q, gGetVertex( this, origin ), 0 );
        
        /*
         *  CENTRAL ITERATION
         *
         */
        
        while( !dirfibheap_empty( q ) ) {                  //Until the priority queue is empty:
            u = dirfibheap_extract_min( q );                 //get the lowest-weight Vertex 'u',
            
            if( !strcmp( u->label, target ) ) {                //(end search if reached destination vertex)
                printf("total iterations: %d", iterations);
                break;
            }
            iterations += 1;

            spt_u = sptGetVertex( spt, u->label );             //get corresponding SPT Vertex,
            
            if( spt_u->hop >= hoplimit ) {
                break;
            }
            
            if( spt_u->state->weight > weightlimit ) {
                break;
            }
            
            du = (State*)spt_u->state;                     //and get State of u 'du'.
            
            if( du->time > maxtime )
                break;
            
            ListNode* edges = vGetOutgoingEdgeList( u );
            while( edges ) {                                 //For each Edge 'edge' connecting u
                Edge* edge = edges->data;
                v = edge->to;                                  //to Vertex v:
                
                long old_w;
                if( (spt_v = sptGetVertex( spt, v->label )) ) {        //get the SPT Vertex corresponding to 'v'
                    dv = (State*)spt_v->state;                     //and its State 'dv'
                    old_w = dv->weight;
                } else {
                    dv = NULL;                                       //which may not exist yet
                    old_w = INFINITY;
                }
                
                State *new_dv = eWalk( edge, du, options );
                
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
                
                long new_w = new_dv->weight;
                // If the new way of getting there is better,
                if( new_w < old_w ) {
                    // add minimum travel time from assist tree to queue priority (A* style)
                    SPTVertex *assist_v = sptGetVertex(that, v->label);
                    if (assist_v != NULL) {
                        //printf("vertex in assist tree: %s %ld\n", v->label, assist_v->state->time); 
                        //new_w = new_dv->weight - init_state->time + assist_v->state->time; // assist SPT must be built from time=0
                        if (new_dv->num_transfers == du->num_transfers) // give an unfair advantage to boarding edges to avoid hanging around in stations
                            new_w += assist_v->state->time; // assist SPT must be built from time=0
                        dirfibheap_insert_or_dec_key( q, v, new_w );    // rekey v in the priority queue
                    } else { 
                        printf("vertex not in assist tree: %s (not added to queue)\n", v->label); 
                    }
                    
                    // If this is the first time v has been reached
                    if( !spt_v ) {
                        spt_v = sptAddVertex( spt, v, spt_u->hop+1 );        //Copy v over to the SPT
                        count++;
                    }
                    
                    //if((count%10000) == 0)
                    //  fprintf(stdout, "Shortest path tree size: %d\n",count);
                    
                    if(spt_v->state)
                        stateDestroy(spt_v->state);
                    spt_v->state = new_dv;                      //Set the State of v in the SPT to the current winner
                    
                    sptvSetParent( spt_v, spt_u, edge->payload );      //Make u the parent of v in the SPT
                } else {
                    stateDestroy(new_dv); //new_dv will never be used; merge it with the infinite.
                }
                edges = edges->next;
            }
        }
        
        dirfibheap_delete( q );
        
        //fprintf(stdout, "Final shortest path tree size: %d\n",count);
        return spt;
    }
#endif
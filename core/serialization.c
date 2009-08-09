#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdbool.h>
#include <assert.h>

#include "statetypes.h"
#include "graph.h"

/* The following hack allows for development debugging of the serialization code */
#if 0
  #define LOG(...) printf(__VA_ARGS__);
#else
  #define LOG(...) /* __VA_ARGS__ */
#endif 

#define STRING_BUFF(name, len) char name[len]; size_t name ## _buff_len
#define BUFF_SIZE size_t __buff_size = 0

#define FREAD_TYPE(invar, type) LOG("\tr " #invar " " # type "\n"); fread(&invar,sizeof(type),1,f)
#define FWRITE_TYPE(outvar, type) LOG("\tw " #outvar " " # type "\n"); fwrite(&(outvar),sizeof(type),1,f)
#define FREAD_STRING(buff) fread(&buff ## _buff_len,sizeof(size_t), 1, f); fread(buff,sizeof(char)*( buff ## _buff_len ), 1, f); buff[buff ## _buff_len] = 0
#define FWRITE_STRING(buff,buff_len) __buff_size = strlen(buff); LOG("\tw str[%d] @ %d\n", __buff_size, __LINE__); assert(__buff_size < buff_len); fwrite(&__buff_size,sizeof(size_t),1, f); fwrite(buff,1,__buff_size, f)

#define DECLARE_SERIALIZE_FUNCS(name, type) EdgePayload* name ## Deserialize(FILE*,void * mm_data); void name ## Serialize(type,FILE*,FILE*);
#define DECLARE_SERIALIZE_FUNCS_W_CAL_TZ(name, type) EdgePayload* name ## Deserialize(ServiceCalendar*, Timezone*, FILE*, void * mm_data); void name ## Serialize(type,FILE*, FILE*);
#define __SERIALIZE(name,type,obj) LOG("  Serializing " #name "\n"); name ## Serialize((type)obj, f, mmf) 
#define __DESERIALIZE(name) LOG("Deserializing " #name "\n"); gAddEdge(g, label, label2, name ## Deserialize(f, mm_data))
#define __DESERIALIZE_W_CAL_TZ(name, sc, tz) LOG("Deserializing " #name "\n"); gAddEdge(g, label, label2, name ## Deserialize(sc, tz, f, mm_data))

#define MM_POS_WRITING long mm_cur_pos = ftell(mmf);
#define MM_POS_READING long mm_cur_pos = 0;
#define FWRITE_MM_POS fwrite(&mm_cur_pos,sizeof(long),1,f)
#define FREAD_MM_POS FREAD_TYPE(mm_cur_pos, long);

#define MWRITE_STRING(buff) \
__buff_size = strlen(buff); \
LOG("\twm str[%d] @ %d\n", __buff_size, __LINE__); \
fwrite(buff,1,__buff_size+1,mmf);

#define MWRITE_TYPE(outvar,type) \
mm_cur_pos = ftell(mmf); \
fwrite(&outvar,1,sizeof(type),mmf);
// \
//mm_seg_len = mm_seg_len + sizeof(type)

#define FWRITE_MM_POS mm_cur_pos = ftell(mmf); fwrite(&mm_cur_pos,sizeof(long),1,f);

//assigns a pointer to data in the mmaped file of the size "type".
#define MREAD(assign_pointer_to, type) \
assign_pointer_to = (type *)(mm_data+mm_cur_pos);

//We could put a function here to set the mmap pointer to a new block to help with caching.

DECLARE_SERIALIZE_FUNCS(street, Street*);
DECLARE_SERIALIZE_FUNCS(elapseTime, ElapseTime*);
DECLARE_SERIALIZE_FUNCS(crossing, Crossing*);
DECLARE_SERIALIZE_FUNCS(link, Link*);

DECLARE_SERIALIZE_FUNCS_W_CAL_TZ(al,Alight*);
DECLARE_SERIALIZE_FUNCS_W_CAL_TZ(tb,TripBoard*);
DECLARE_SERIALIZE_FUNCS_W_CAL_TZ(hwb,HeadwayBoard*);
DECLARE_SERIALIZE_FUNCS_W_CAL_TZ(hwa,HeadwayAlight*);

ServiceCalendar* scDeserialize(FILE* f, void* mm_data);
void scSerialize(ServiceCalendar* tz, FILE* f, FILE* mmf);

Timezone* tzDeserialize(FILE* f, void* mm_data);
void tzSerialize(Timezone* tz, FILE* f, FILE* mmf);

void get_cal_and_tz(EdgePayload* payload, ServiceCalendar** sc, Timezone** tz);

bool gDeserialize(Graph *g, char* gbin_name, char * mmf_name) {
	FILE *f, *mmf = NULL;
	void * mm_data = NULL;
	
	char f_ind_name[1024];
	long num_vertices, num_edges, num_calendars, num_timezones;
	STRING_BUFF(label,256);
	STRING_BUFF(label2,256);
	
	int ep_type;
	long i;
	bool flag;
	//int* ip;
	//int* edges;
	//int bytes = 0;
	
	ServiceCalendar** calendars = NULL;
	Timezone** timezones = NULL; 
	
	// open the index, get the number of vertices
	sprintf(f_ind_name, "%s", gbin_name);
	if ((f = fopen(f_ind_name, "rb")) == NULL) {
		return false;
	}
	
	if (mmf_name != NULL)
		sprintf(f_ind_name, "%s", mmf_name);
	
	if (mmf_name == NULL || (mmf = fopen(f_ind_name, "rb")) == NULL) {
		LOG("mem map failed. could not open file graph\n");
		mmf = NULL;
	   //Maybe they didn't want to memory map it?
	} else {
	   LOG("deserializing memory mapped graph\n");
	   //okay they want to memory map. Lets do it.
	   //current hardcoded max of 100 mb mmaped file.
	   mm_data = mmap((caddr_t)0, 100*1024*1024, PROT_READ, MAP_SHARED, fileno(mmf), (off_t)0);
	   assert(mm_data != (void *)-1);
	}
	
	//Read the # of vertices.
	
	FREAD_TYPE(num_vertices, long);
	LOG("num_vertices=%ld\n", num_vertices);
	for (i = 0; i < num_vertices; i++) {
		FREAD_STRING(label);
		gAddVertex(g, label);
	}
	
	FREAD_TYPE(num_calendars, int);
	LOG("num_calendars=%ld\n", num_calendars);
	if (num_calendars) {
		calendars = (ServiceCalendar**)malloc(sizeof(ServiceCalendar*)*num_calendars);
		for (i = 0; i < num_calendars; i++) {
			calendars[i] = scDeserialize(f,mm_data);
		}
	}
	FREAD_TYPE(num_timezones, int);
	LOG("num_timezones=%ld\n", num_timezones);
	if (num_timezones) {
		timezones = (Timezone**)malloc(sizeof(Timezone*)*num_timezones);
		for (i = 0; i < num_timezones; i++) {
			timezones[i] = tzDeserialize(f,mm_data);
		}
	}
	
	FREAD_TYPE(num_edges, long);
	LOG("num_edges=%ld\n", num_edges);
	long cal_num = 0;
	long tz_num = 0;
	for (i = 0; i < num_edges; i++) {
		LOG("New edge\n");
		FREAD_TYPE(ep_type, edgepayload_t);
		FREAD_TYPE(flag, bool);
		if (flag) { // does this edge have a timezone/calendar?
			FREAD_TYPE(cal_num, long);
			FREAD_TYPE(tz_num, long);
			LOG("cal_num=%ld tz_num=%ld\n",cal_num, tz_num);
			assert(cal_num < num_calendars);
			assert(tz_num < num_timezones);
		}
		FREAD_STRING(label); // from
		FREAD_STRING(label2); // to
		switch (ep_type) {
			case PL_STREET: __DESERIALIZE(street); break;
			case PL_LINK:  __DESERIALIZE(link); break;
			case PL_TRIPBOARD:  __DESERIALIZE_W_CAL_TZ(tb, calendars[cal_num], timezones[tz_num]); break;
			case PL_CROSSING: __DESERIALIZE(crossing); break;
			case PL_ALIGHT: __DESERIALIZE_W_CAL_TZ(al, calendars[cal_num], timezones[tz_num]); break;
			case PL_HEADWAYBOARD: __DESERIALIZE_W_CAL_TZ(hwb, calendars[cal_num], timezones[tz_num]); break;
			case PL_HEADWAYALIGHT: __DESERIALIZE_W_CAL_TZ(hwa, calendars[cal_num], timezones[tz_num]); break;
			case PL_ELAPSE_TIME: __DESERIALIZE(elapseTime); break;
			default: 
				assert(strcmp("Unsupported type","") == 0); // unsupported type
				return false;
		}
		LOG("End edge\n");

	}
	if (calendars) { free(calendars); }
	if (timezones) { free(timezones); }
	
	fclose(f);
	fclose(mmf);
	
	return true;
}

bool gSerialize(Graph *g, char* gbin_name, char * mmf_name) {
	FILE *f, *mmf;
	char f_ind_name[1024];
	long num_vertices, num_edges;
	int num_calendars, num_timezones;
	ServiceCalendar* calendars[32];
	Timezone* timezones[32];
	BUFF_SIZE;
	long i, j;
	Vertex** verts;
	bool flag;
	
	// open the index, get the number of vertices
	sprintf(f_ind_name, "%s", gbin_name);
	if ((f = fopen(f_ind_name, "wb")) == NULL) {
		return false;
	}

	sprintf(f_ind_name, "%s", mmf_name);
	if ((mmf = fopen(f_ind_name, "wb")) == NULL) {

	}
	
	verts = gVertices(g, &num_vertices);
	num_edges = 0;
	num_timezones = 0;
	num_calendars = 0;
	
	//calendars = (ServiceCalendar**)malloc(sizeof(ServiceCalendar*)*32);
	//timezones = (Timezone**)malloc(sizeof(Timezone*)*32);
	
	FWRITE_TYPE(num_vertices, long);
	for (i = 0; i < num_vertices; i++) {
		ListNode* out; 
		int v_num_edges = 0;
		Vertex* v = verts[i];
		FWRITE_STRING(v->label, 256);
		out = vGetOutgoingEdgeList(v);
		while (out) { //out is a ListNode (linked list), set to next on each iteration.
			num_edges++;
			v_num_edges++;
			if (out->data && out->data->payload) {
				Timezone* tz = NULL;
				ServiceCalendar* sc = NULL;
				get_cal_and_tz(out->data->payload, &sc, &tz);
				if (sc) {
					for (j = 0; j < num_calendars; j++) {
						if (calendars[j] == sc) {
							sc = NULL;
							break;
						}
					}
					if (sc) { //if we found a unique service calendar
						calendars[num_calendars] = sc;
						num_calendars++;
					}
				}
 				if (tz) {
					for (j = 0; j < num_timezones; j++) {
						if (timezones[j] == tz) {
							tz = NULL;
							break;
						}
					}
					if (tz) { //if we found a unique timezone
						timezones[num_timezones] = tz;
						num_timezones++;
					}
				}
			} //end if (out->data && out->data->payload) 
			out = out->next;
		} //end while
		//FWRITE_TYPE(v_num_edges, int);
	}
	
	//Now that we've found all of the unique service calendars and timezones...
	//let's serialize them
	
	FWRITE_TYPE(num_calendars, long);
	for (i = 0; i < num_calendars; i++) {
		scSerialize(calendars[i], f, mmf);
	}
	
	FWRITE_TYPE(num_timezones, long);
	for (i = 0; i < num_timezones; i++) {
		tzSerialize(timezones[i], f, mmf);
	}
	
	// Now write all edges, by traversing all vertices' outgoing.	
	FWRITE_TYPE(num_edges, long);
	for (i = 0; i < num_vertices; i++) {
		ListNode* out;
		Vertex* v = verts[i];
		num_edges = 0;
		out = vGetOutgoingEdgeList(v);
		while (out) { //while we have any more outgoing edges we want to serialize on this vertex
			if (!out->data) continue; //why would this happen? Just in case;
			//FWRITE_STRING("\n", 10); //why are we doing this? Just for appearances?
			LOG("New edge\n");
			EdgePayload* p = out->data->payload;
			// handle timezone stuff
			Timezone* tz = NULL;
			ServiceCalendar* sc = NULL;
			// write the type
			FWRITE_TYPE(p->type, edgepayload_t);

			get_cal_and_tz(out->data->payload, &sc, &tz);
			if (sc && tz) {
				flag = true; //yes, this edge has calendars and timezones
				FWRITE_TYPE(flag,bool);
				flag = false;
				//find which service calendar this edge has
				for (j = 0; j < num_calendars; j++) {
					if (calendars[j] == sc) {
						FWRITE_TYPE(j, long);
						flag = true;
						break;
					}
				}
				assert(flag == true && j < num_calendars);
				
				//find which timezone this edge has
				for (j = 0; j < num_timezones; j++) {
					if (timezones[j] == tz) {
						FWRITE_TYPE(j, long);
						flag = true;
						break;
					}
				}
				assert(flag == true && j < num_timezones);
			} else {
				flag = false; //No, this edge does not have calendars and timezones
				FWRITE_TYPE(flag,bool);
			}
			// write the from/to labels
			FWRITE_STRING(v->label, 64);
			FWRITE_STRING(out->data->to->label, 64);	
			switch (p->type) {
				case PL_STREET: __SERIALIZE(street, Street*, p); break;
				case PL_LINK:  __SERIALIZE(link, Link*, p); break;
				case PL_TRIPBOARD:  __SERIALIZE(tb, TripBoard*, p); break;
				case PL_CROSSING: __SERIALIZE(crossing, Crossing*, p); break;
				case PL_ALIGHT: __SERIALIZE(al, Alight*, p); break;
				case PL_HEADWAYBOARD: __SERIALIZE(hwb, HeadwayBoard*, p); break;
				case PL_HEADWAYALIGHT: __SERIALIZE(hwa, HeadwayAlight*, p); break;
				case PL_ELAPSE_TIME: __SERIALIZE(elapseTime, ElapseTime*, p); break;
				default: 
					assert(NULL); // unsupported type
					return false;
			}
			LOG("Done edge\n");
			out = out->next;
		} //end while (out)
	} //end for (i = 0; i < num_vertices; i++)
	fclose(f);
	
	return true;
}

void get_cal_and_tz(EdgePayload* payload, ServiceCalendar** sc, Timezone** tz) {
	switch (payload->type) {
		case PL_HEADWAY:
			*tz = ((Headway*)payload)->timezone;
			*sc = ((Headway*)payload)->calendar;
			break;
		case PL_TRIPBOARD:
			*tz = ((TripBoard*)payload)->timezone;
			*sc = ((TripBoard*)payload)->calendar;
			break;
		case PL_ALIGHT:
			*tz = ((Alight*)payload)->timezone;
			*sc = ((Alight*)payload)->calendar;
			break;
		case PL_HEADWAYBOARD:
			*tz = ((HeadwayBoard*)payload)->timezone;
			*sc = ((HeadwayBoard*)payload)->calendar;
			break;
		case PL_HEADWAYALIGHT:						
			*tz = ((HeadwayAlight*)payload)->timezone;
			*sc = ((HeadwayAlight*)payload)->calendar;
			break;
		default: 
			*tz = NULL;
			*sc = NULL;
			break;
	}
}

EdgePayload*
streetDeserialize(FILE* f, void* mm_data) {
	// "h(s)f"
	double s_len;
	STRING_BUFF(name, 1024);
	FREAD_STRING(name);
	FREAD_TYPE(s_len, double);
	return (EdgePayload*)streetNew(name, s_len);
}

void
streetSerialize(Street* s, FILE* f, FILE* mmf) {
	// "h(s)f"
	BUFF_SIZE;
	FWRITE_STRING(s->name, 1024);
	FWRITE_TYPE(s->length, double);
}

void test_serialize(char* filename) {
	BUFF_SIZE;
	FILE* f;
	Street* s;
	STRING_BUFF(buff, 1024);
	long l = 1000000L;
	f = fopen(filename, "wb");
	FWRITE_STRING("cows on the loose", 32);
	FWRITE_TYPE(l, long);
	streetSerialize(streetNew("A street", 1), f, NULL);
	fclose(f);
	f = fopen(filename, "rb");
	FREAD_STRING(buff);
	l = 0;
	FREAD_TYPE(l,long);
	s = (Street*)streetDeserialize(f,NULL);
	LOG("read: streetName='%s', streetLen='%f'\n", s->name, s->length);
	assert(strcmp(s->name,"A street") == 0);
	assert(s->length == 1);
	fclose(f);
}

void test_ServiceCalendar_serialze(char * filename) {
	BUFF_SIZE;
	FILE* f;
	
	ServiceCalendar * sc = scNew();
	scAddServiceId(sc, "WKDY");
	scAddServiceId(sc, "WKND");
	
	f = fopen(filename, "wb");
	scSerialize(sc, f, NULL);
	fclose(f);
	
	f = fopen(filename, "rb");
	ServiceCalendar * sc2 = scDeserialize(f,NULL);
	LOG("read num_sides=%i", sc2->num_sids);
	
	//not yet finished
	
	//ServicePeriod * sp = spNew(0,1000, 1 )

}

EdgePayload*
elapseTimeDeserialize(FILE* f, void* mm_data) {
	// "l"
	long secs = 0;
	FREAD_TYPE(secs, long);
	return (EdgePayload*)elapseTimeNew(secs);
}

void
elapseTimeSerialize(ElapseTime *e, FILE* f, FILE* mmf) {	
	// "l"
	FWRITE_TYPE(e->seconds, long);
}


EdgePayload*
crossingDeserialize(FILE* f, void* mm_data) {
	// "i"
	int secs = 0;
	FREAD_TYPE(secs, int);
	return (EdgePayload*)crNew(secs);
}
void
crossingSerialize(Crossing *e, FILE* f, FILE* mmf) {	
	// "l"
	FWRITE_TYPE(e->crossing_time, int);
}

EdgePayload*
linkDeserialize(FILE* f, void* mm_data) {
	return (EdgePayload*)linkNew();
}

void
linkSerialize(Link *e, FILE* f, FILE * mmf) {	
	// "l"
	return;	
}


void
alSerialize(Alight* a, FILE* f, FILE * mmf) {
	BUFF_SIZE;
	
	ServiceId service_id = alGetServiceId(a);
	int agency = alGetAgency(a);
	int trip_cnt = a->n;
	int overage = alGetOverage(a);	
	
    FWRITE_TYPE(service_id, ServiceId);
    FWRITE_TYPE(trip_cnt, int);
    FWRITE_TYPE(agency, int);
    FWRITE_TYPE(overage, int);

	if (mmf) {
	   MM_POS_WRITING;
	
	   int i;
        for (i = 0; i < trip_cnt; i++) {
            FWRITE_MM_POS;
            MWRITE_STRING(a->trip_ids[i]);
        }
        FWRITE_MM_POS;
        for (i = 0; i < trip_cnt; i++) {
            MWRITE_TYPE(a->arrivals[i], int);
        }
	} else {
        int i;
        for (i = 0; i < trip_cnt; i++) {
            FWRITE_STRING(a->trip_ids[i],512);
            FWRITE_TYPE(a->arrivals[i],int);
        }
	}
}

EdgePayload*
alDeserialize(ServiceCalendar* calendar, Timezone* timezone, FILE* f, void* mm_data) {
	if (mm_data) {
	   MM_POS_READING;
	
		ServiceId sid;
	    int trip_cnt, agency;

	   FREAD_TYPE(sid, ServiceId);
	   FREAD_TYPE(trip_cnt, int);
        FREAD_TYPE(agency, int);

	
	   Alight* al = alNew(sid, calendar, timezone, agency);
        FREAD_TYPE(al->overage, int);
    
        al->n = trip_cnt;
        al->trip_ids = (char **)malloc(sizeof(char *)*trip_cnt);
        		
        int i;
        for (i = 0; i < trip_cnt; i++) { //for each trip_id string
            FREAD_MM_POS;
            MREAD(al->trip_ids[i], char); //point this cstring to the mmapped bytes.
        }
		
        FREAD_MM_POS;
        MREAD(al->arrivals, int); //assign the address of the begining of the mmaped arrival array
        
        return (EdgePayload *)al;
	} else {
        ServiceId sid;
	    int trip_cnt, agency, arrival;
	
        STRING_BUFF(trip_id, 512);
    
        FREAD_TYPE(sid, ServiceId);
        FREAD_TYPE(trip_cnt, int);
        FREAD_TYPE(agency, int);
        Alight* al = alNew(sid, calendar, timezone, agency);
        FREAD_TYPE(al->overage, int);
        
        int i;
        for (i = 0; i < trip_cnt; i++) {
            FREAD_STRING(trip_id);
            FREAD_TYPE(arrival, int);
            alAddAlighting(al, trip_id, arrival);
        }	
        return (EdgePayload *)al;
	}
}

void
tbSerialize(TripBoard* tb, FILE* f, FILE * mmf) {
	BUFF_SIZE;
	
	ServiceId service_id = tbGetServiceId(tb);
	int agency = tbGetAgency(tb);
	int trip_cnt = tb->n;
	int overage = tbGetOverage(tb);	
	
    FWRITE_TYPE(service_id, ServiceId);
    FWRITE_TYPE(trip_cnt, int);
    FWRITE_TYPE(agency, int);
    FWRITE_TYPE(overage, int);
	
	if (mmf) {
		MM_POS_WRITING;
		
		int i;
        for (i = 0; i < trip_cnt; i++) {
            FWRITE_MM_POS;
            MWRITE_STRING(tb->trip_ids[i]);
        }
        FWRITE_MM_POS;
        for (i = 0; i < trip_cnt; i++) {
            MWRITE_TYPE(tb->departs[i], int);
        }
	} else {
        int i;
        for (i = 0; i < trip_cnt; i++) {
            FWRITE_STRING(tb->trip_ids[i],512);
            FWRITE_TYPE(tb->departs[i],int);
        }
	}
}

EdgePayload*
tbDeserialize(ServiceCalendar* calendar, Timezone* timezone, FILE* f, void* mm_data) {
	LOG("tb deserialize")
	
	if (mm_data) {
		MM_POS_READING;
		
		ServiceId sid;
	    int trip_cnt, agency;
		
		FREAD_TYPE(sid, ServiceId);
		FREAD_TYPE(trip_cnt, int);
        FREAD_TYPE(agency, int);
		
		
		TripBoard* tb = tbNew(sid, calendar, timezone, agency);
        FREAD_TYPE(tb->overage, int);
		
        tb->n = trip_cnt;
        tb->trip_ids = (char **)malloc(sizeof(char *)*trip_cnt);
		
        int i;
        for (i = 0; i < trip_cnt; i++) { //for each trip_id string
            FREAD_MM_POS;
            MREAD(tb->trip_ids[i], char); //point this cstring to the mmapped bytes.
        }
		
        FREAD_MM_POS;
        MREAD(tb->departs, int); //assign the address of the begining of the mmaped departure array
        
        return (EdgePayload *)tb;
	} else {
        ServiceId sid;
	    int trip_cnt, agency, departure;
		
        STRING_BUFF(trip_id, 512);
		
        FREAD_TYPE(sid, ServiceId);
        FREAD_TYPE(trip_cnt, int);
        FREAD_TYPE(agency, int);
        TripBoard* tb = tbNew(sid, calendar, timezone, agency);
        FREAD_TYPE(tb->overage, int);
        
        int i;
        for (i = 0; i < trip_cnt; i++) {
            FREAD_STRING(trip_id);
            FREAD_TYPE(departure, int);
            tbAddBoarding(tb, trip_id, departure);
        }	
        return (EdgePayload *)tb;
	}
}


void
hwbSerialize(HeadwayBoard* o, FILE* f, FILE * mmf) {
	BUFF_SIZE;	
	FWRITE_TYPE(o->service_id,ServiceId);
	FWRITE_TYPE(o->agency, int);
	FWRITE_STRING(o->trip_id,512);
	FWRITE_TYPE(o->start_time, int);
	FWRITE_TYPE(o->end_time, int);
	FWRITE_TYPE(o->headway_secs, int);
}

EdgePayload*
hwbDeserialize(ServiceCalendar* calendar, Timezone* timezone, FILE* f, void* mm_data) {
	ServiceId service_id;
	int agency, start_time, end_time, headway_secs;
	
	STRING_BUFF(trip_id,512);
	
	FREAD_TYPE(service_id, ServiceId);
	FREAD_TYPE(agency, int);
	FREAD_STRING(trip_id);
	FREAD_TYPE(start_time, int);
	FREAD_TYPE(end_time, int);
	FREAD_TYPE(headway_secs, int);
	
	return (EdgePayload*)hbNew(service_id, calendar, timezone, agency, trip_id, start_time, end_time, headway_secs );	
}

void
hwaSerialize(HeadwayAlight* o, FILE* f, FILE* mmf) {
	BUFF_SIZE;	
	FWRITE_TYPE(o->service_id,ServiceId);
	FWRITE_TYPE(o->agency, int);
	FWRITE_STRING(o->trip_id,512);
	FWRITE_TYPE(o->start_time, int);
	FWRITE_TYPE(o->end_time, int);
	FWRITE_TYPE(o->headway_secs, int);
}

EdgePayload*
hwaDeserialize(ServiceCalendar* calendar, Timezone* timezone, FILE* f, void* mm_data) {
	ServiceId service_id;
	int agency, start_time, end_time, headway_secs;
	
	STRING_BUFF(trip_id,512);
	
	FREAD_TYPE(service_id, ServiceId);
	FREAD_TYPE(agency, int);
	FREAD_STRING(trip_id);
	FREAD_TYPE(start_time, int);
	FREAD_TYPE(end_time, int);
	FREAD_TYPE(headway_secs, int);
	
	return (EdgePayload*)haNew(service_id, calendar, timezone, agency, trip_id, start_time, end_time, headway_secs );	
}


void scSerialize(ServiceCalendar * sc, FILE * f, FILE* mmf) {
    BUFF_SIZE;
	
	int i, j, num_sps, num_sids;
	
	ServicePeriod * service_period = scHead(sc);
    
	num_sids = service_period->n_service_ids;
	num_sps = 0;
    while (service_period) {
        num_sps++;
		service_period = service_period->next_period;
	}
	
	
	FWRITE_TYPE(num_sids, int); //write the number of service ids
	FWRITE_TYPE(num_sps, int); //write the number of service calendars
    
    for (i = 0; i < num_sids; i++){
		FWRITE_STRING(sc->sid_int_to_str[i],1024);
        LOG("sid_name %s\n", sc->sid_int_to_str[i]);

		//write the service id name
    }
	
    service_period = scHead(sc);
	while (service_period) {
		FWRITE_TYPE(service_period->begin_time, long);
		FWRITE_TYPE(service_period->end_time, long);
		FWRITE_TYPE(service_period->n_service_ids, int);
		for (j = 0; j < service_period->n_service_ids; j++) {
			FWRITE_TYPE(service_period->service_ids[j],ServiceId);
		}
		service_period = service_period->next_period;
	}
	
}

ServiceCalendar*
scDeserialize(FILE* f, void* mm_data) {
	int i, j, num_sps, num_sids;
	STRING_BUFF(sid_name, 1024);
	ServicePeriod *p;
	ServicePeriod pt;
	ServiceCalendar* sc;
	
	FREAD_TYPE(num_sids, int);
	FREAD_TYPE(num_sps, int);
	LOG("num_sids=%d, num_sps=%d\n", num_sids, num_sps);
	sc = scNew();
	
	for (i = 0; i < num_sids; i++) {
		FREAD_STRING(sid_name);
        LOG("sid_name %s\n", sid_name);
		scAddServiceId(sc, sid_name);
	}
	
	for (i = 0; i < num_sps; i++) {
		FREAD_TYPE(pt.begin_time, long);
		FREAD_TYPE(pt.end_time, long);
		FREAD_TYPE(pt.n_service_ids, int);
		int * sids = (int*)malloc(sizeof(int)*pt.n_service_ids);
		for (j = 0; j < pt.n_service_ids; j++) {
			FREAD_TYPE(sids[j], ServiceId);
		}
		p = spNew(pt.begin_time, pt.end_time, num_sps, sids);
		free(sids);
		scAddPeriod(sc, p);
	}
	return sc;
}

void
tzSerialize(Timezone* sc, FILE* f, FILE* mmf) {
	int i;
	
	TimezonePeriod* tzp = sc->head;
	int num_periods = 0;
	while (tzp) {
		num_periods++;
		tzp = tzp->next_period;
	}
	
	FWRITE_TYPE(num_periods, int);
	//write num_periods
	tzp = sc->head;
	for (i = 0; i < num_periods; i++) {
		FWRITE_TYPE(tzp->begin_time, long);
		FWRITE_TYPE(tzp->end_time, long);
		FWRITE_TYPE(tzp->utc_offset, int);
		tzp = tzp->next_period;
	}
}


Timezone* tzDeserialize(FILE* f, void* mm_data) {
	Timezone * tz = tzNew();
	int i, num_periods;
	long begin_time, end_time, utc_offset;
	TimezonePeriod* tzp;
	
	FREAD_TYPE(num_periods,int);
	for (i = 0; i < num_periods; i++) {
		FREAD_TYPE(begin_time, long);
		FREAD_TYPE(end_time, long);
		FREAD_TYPE(utc_offset, int);
		
		tzp = tzpNew(begin_time,end_time,utc_offset);
		tzAddPeriod(tz, tzp);
	}	
	return tz;
}


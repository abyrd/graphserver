from optparse import OptionParser
from graphserver.graphdb import GraphDatabase
import os
from graphserver.core import Street
from graphserver.ext.osm.osmdb import OSMDB
import sys
from graphserver.ext.osm.profiledb import ProfileDB
from gdb_import_ned import get_rise_and_fall

def gdb_import_osm(gdb, osmdb, vertex_namespace, slogs, profiledb=None):
    cursor = gdb.get_cursor()
	
    n_edges = osmdb.count_edges()
    
    street_id_counter = 0
    street_names = {}
    # note: tags is not in the edges table, it's a join.     
    for i, (id, parent_id, node1, node2, distance, tags) in enumerate( osmdb.edges() ):
        
        if i % 10000==0: sys.stdout.write( "%d / %d edges loaded\n"%(i, n_edges))
            
        # Find rise/fall of edge, if profiledb is given
        rise=0
        fall=0
        if profiledb:
            profile = profiledb.get( id )
            if profile:
                rise, fall = get_rise_and_fall( profile )
        
        # insert end vertices of edge to graph
        vertex1_label = "%s-%s"%(vertex_namespace,node1)
        vertex2_label = "%s-%s"%(vertex_namespace,node2)
        gdb.add_vertex( vertex1_label, cursor )
        gdb.add_vertex( vertex2_label, cursor )
                
        # create ID for the way's street
        street_name = tags.get("name")
        if street_name is None:
            street_id_counter += 1
            street_id = street_id_counter
        else:
            if street_name not in street_names:
                street_id_counter += 1
                street_names[street_name] = street_id_counter
            street_id = street_names[street_name]
        
        # Create edges to be inserted into graph
        # print id, distance, rise, fall
        id = str(id)
        s1 = Street( id, distance, rise, fall )
        s2 = Street( id, distance, fall, rise )
        s1.way = street_id
        s2.way = street_id
        
        # See if the way's highway tag is penalized with a 'slog' value; if so, set it in the edges
        slog = slogs.get( tags.get("highway") )
        if slog:
            s1.slog = s2.slog = slog
        
        # Add the forward edge and the return edge if the edge is not oneway
        gdb.add_edge( vertex1_label, vertex2_label, s1, cursor )
        oneway = tags.get("oneway")
        if oneway != "true" and oneway != "yes":
            gdb.add_edge( vertex2_label, vertex1_label, s2, cursor )

    sys.stdout.write( "loading vertex coordinates...\n" )
    for id, lon, lat in osmdb.execute( "SELECT id, x(geometry), y(geometry) FROM vertices" ) :
        gdb.add_coords( "%s-%s" % (vertex_namespace, id), lon, lat, cursor )
    
    gdb.commit()
    
    print "indexing vertices..."
    gdb.index()

def main():
    usage = """usage: python gdb_import_osm.py <graphdb_filename> <osmdb_filename>"""
    parser = OptionParser(usage=usage)
    parser.add_option("-n", "--namespace", dest="namespace", default="osm",
                      help="prefix all imported vertices with namespace string")
    parser.add_option("-s", "--slog",
                      action="append", dest="slog_strings", default=[],
                      help="specify slog for highway type, in highway_type:slog form. For example, 'motorway:10.5'")
    parser.add_option("-p", "--profiledb", dest="profiledb_filename", default=None,
                      help="specify profiledb to annotate streets with rise/fall data")
    
    (options, args) = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        exit(-1)
        
    slogs = {}
    for slog_string in options.slog_strings:
        highway_type,slog_penalty = slog_string.split(":")
        slogs[highway_type] = float(slog_penalty)
    print "slog values: %s"%slogs
        
    graphdb_filename = args[0]
    osmdb_filename = args[1]
    
    print "importing osmdb '%s' into graphdb '%s'"%(osmdb_filename, graphdb_filename)
    
    profiledb = ProfileDB( options.profiledb_filename ) if options.profiledb_filename else None
    osmdb = OSMDB( osmdb_filename )
    gdb = GraphDatabase( graphdb_filename, overwrite=False )
    
    osmdb.segments_to_edges()
    gdb_import_osm(gdb, osmdb, options.namespace, slogs, profiledb);
    
    print "done"

if __name__ == '__main__':
    main()

    

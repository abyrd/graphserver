from graphserver.ext.gtfs.gtfsdb import GTFSDatabase
from graphserver.ext.osm.osmdb import OSMDB
from graphserver.graphdb import GraphDatabase
from graphserver.core import Link

import sys
from optparse import OptionParser

def main():
    usage = """usage: python gdb_link_osm_gtfs.py <osmdb_filename> <gtfsdb_filename>"""
    parser = OptionParser(usage=usage)
    
    (options, args) = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        exit(-1)
        
    osmdb_filename   = args[0]
    gtfsdb_filename  = args[1]
    
    gtfsdb = GTFSDatabase( gtfsdb_filename )
    osmdb = OSMDB( osmdb_filename )

    cur = gtfsdb.get_cursor()
    cur.execute( "DROP TABLE IF EXISTS osm_links" )
    # should also store distance from vertex to stop.
    cur.execute( "CREATE TABLE osm_links (gtfs_stop TEXT, osm_vertex TEXT, PRIMARY KEY(gtfs_stop))" )
    gtfsdb.conn.commit() 
    nstops = gtfsdb.count_stops()
    # for every stop in the GTFS database
    for (i, (stopid, name, lat, lon)) in enumerate(gtfsdb.stops()) :
        print "stop %i / %i :" % (i, nstops), stopid, name, lat, lon
        link_vid = osmdb.find_or_make_link_vertex(lat, lon)
        if link_vid :
            cur.execute( "INSERT INTO osm_links VALUES (?, ?)", (stopid, link_vid) )
            gtfsdb.conn.commit()
        else :
            print "for some reason, the osmdb didn't return a vertex."
        # debug
        if i > 20 : break
    
    print "DONE."

if __name__=='__main__':
    main()

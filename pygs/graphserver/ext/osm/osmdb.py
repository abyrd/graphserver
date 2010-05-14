# ugly import statement allows minimum changes
# when python's built-in pysqlite allows loading extensions
# (should be in python version 2.7)
from pysqlite2 import dbapi2 as sqlite3
import os
try:
    import json
except ImportError:
    import simplejson as json
import sys
import xml.sax
import binascii
from vincenty import vincenty
from struct import pack, unpack
from rtree import Rtree

def cons(ary):
    for i in range(len(ary)-1):
        yield (ary[i], ary[i+1])

def pack_coords(coords):
    return binascii.b2a_base64( "".join([pack( "ff", *coord ) for coord in coords]) )
        
def unpack_coords(str):
    bin = binascii.a2b_base64( str )
    return [unpack( "ff", bin[i:i+8] ) for i in range(0, len(bin), 8)]

class Node:
    def __init__(self, id, lon, lat):
        self.id = id
        self.lon = lon
        self.lat = lat
        self.tags = {}

    def __repr__(self):
        return "<Node id='%s' (%s, %s) n_tags=%d>"%(self.id, self.lon, self.lat, len(self.tags))
        
class Way:
    def __init__(self, id):
        self.id = id
        self.nd_ids = []
        self.tags = {}
        
    def __repr__(self):
        return "<Way id='%s' n_nds=%d n_tags=%d>"%(self.id, len(self.nd_ids), len(self.tags))

class WayRecord:
    def __init__(self, id, tags, nds):
        self.id = id
        
        if type(tags)==unicode:
            self.tags_str = tags
            self.tags_cache = None
        else:
            self.tags_cache = tags
            self.tags_str = None
            
        if type(nds)==unicode:
            self.nds_str = nds
            self.nds_cache = None
        else:
            self.nds_cache = nds
            self.nds_str = None
        
    @property
    def tags(self):
        self.tags_cache = self.tags_cache or json.loads(self.tags_str)
        return self.tags_cache
        
    @property
    def nds(self):
        self.nds_cache = self.nds_cache or json.loads(self.nds_str)
        return self.nds_cache
        
    def __repr__(self):
        return "<WayRecord id='%s'>"%self.id

class OSMDB:
    def __init__(self, dbname,overwrite=False,rtree_index=True):
        self.dbname = dbname
        
        if overwrite:
            try:
                os.remove( dbname )
            except OSError:
                pass
            
        self.conn = sqlite3.connect(dbname)
        try:
            self.conn.enable_load_extension(True)
        except:
            print "To load spatialite extensions, pysqlite must have enable_load_extension compiled in."
            sys.exit(3)

        try:
            self.conn.execute("SELECT load_extension('libspatialite.so.2')")
        except:
            try:
                self.conn.execute("SELECT load_extension('libspatialite.dylib')")
            except:    
                sys.stderr.write("(Py)SQLite cannot find libspatialite.so.2 or libspatialite.dylib.\n")
                sys.exit(4)
        
        if overwrite:
            self.setup()
            
    def get_cursor(self):
        # Attempts to get a cursor using the current connection to the db. If we've found ourselves in a different thread
        # than that which the connection was made in, re-make the connection.
        
        try:
            ret = self.conn.cursor()
        except sqlite3.ProgrammingError:
            self.conn = sqlite3.connect(self.dbname)
            ret = self.conn.cursor()
            
        return ret
        
    def setup(self):
        c = self.get_cursor()
        c.execute( "CREATE TABLE nodes (id INTEGER, tags TEXT, lat FLOAT, lon FLOAT, refs INTEGER DEFAULT 0, alias TEXT)" )
        c.execute( "CREATE TABLE ways (id INTEGER, tags TEXT, nds TEXT)" )
        c.execute( "SELECT InitSpatialMetaData()")
        # instead of initializing all the reference systems, add only WGS84 lat-lon
        c.execute( "INSERT INTO spatial_ref_sys (srid, auth_name, auth_srid, ref_sys_name, proj4text) VALUES (4326, 'epsg', 4326, 'WGS 84', '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')")
        # For nodes, add a 2D point column in WGS84 lat-lon reference system
        c.execute( "SELECT AddGeometryColumn ( 'nodes', 'geometry', 4326, 'POINT', 2 )")
        self.conn.commit()
        c.close()
        
    def create_indexes(self):
        c = self.get_cursor()
        c.execute( "CREATE INDEX nodes_id ON nodes (id)" )
        c.execute( "CREATE INDEX nodes_lon ON nodes (lon)" )
        c.execute( "CREATE INDEX nodes_lat ON nodes (lat)" )
        c.execute( "CREATE INDEX ways_id ON ways (id)" )
        self.conn.commit()
        c.close()
        
    def populate(self, osm_filename, accept=lambda tags: True, reporter=None):
        print "Importing OSM from XML to SQLite database."
        
        c = self.get_cursor()
        
        self.n_nodes = 0
        self.n_ways = 0
        
        superself = self

        class OSMHandler(xml.sax.ContentHandler):
            @classmethod
            def setDocumentLocator(self,loc):
                pass

            @classmethod
            def startDocument(self):
                pass

            @classmethod
            def endDocument(self):
                pass

            @classmethod
            def startElement(self, name, attrs):
                if name=='node':
                    self.currElem = Node(attrs['id'], float(attrs['lon']), float(attrs['lat']))
                elif name=='way':
                    self.currElem = Way(attrs['id'])
                elif name=='tag':
                    self.currElem.tags[attrs['k']] = attrs['v']
                elif name=='nd':
                    self.currElem.nd_ids.append( attrs['ref'] )

            @classmethod
            def endElement(self,name):
                if name=='node':
                    if superself.n_nodes%50000==0:
                        print "Node %d" % superself.n_nodes
                    superself.n_nodes += 1
                    superself.add_node( self.currElem, c )
                elif name=='way':
                    if superself.n_ways%5000==0:
                        print "\rWay %d" % superself.n_ways
                    superself.n_ways += 1
                    superself.add_way( self.currElem, c )

            @classmethod
            def characters(self, chars):
                pass

        xml.sax.parse(osm_filename, OSMHandler)
        
        self.conn.commit()
        c.close()
        
        print "indexing primary tables...",
        self.create_indexes()
        print "done"
        
    def set_endnode_ref_counts( self ):
        """Populate ways.endnode_refs. Necessary for splitting ways into single-edge sub-ways"""
        
        print "counting end-node references to find way split-points"
        
        c = self.get_cursor()
        
        endnode_ref_counts = {}
        
        c.execute( "SELECT nds from ways" )
        
        print "...counting"
        for i, (nds_str,) in enumerate(c):
            if i%5000==0:
                print i
                
            nds = json.loads( nds_str )
            for nd in nds:
                endnode_ref_counts[ nd ] = endnode_ref_counts.get( nd, 0 )+1
        
        print "...updating nodes table"
        for i, (node_id, ref_count) in enumerate(endnode_ref_counts.items()):
            if i%5000==0:
                print i
            
            if ref_count > 1:
                c.execute( "UPDATE nodes SET endnode_refs = ? WHERE id=?", (ref_count, node_id) )
            
        self.conn.commit()
        c.close()
    
    def index_endnodes( self ):
        print "indexing endpoint nodes into rtree"
        
        c = self.get_cursor()
        
        #TODO index endnodes if they're at the end of oneways - which only have one way ref, but are still endnodes
        c.execute( "SELECT id, lat, lon FROM nodes WHERE endnode_refs > 1" )
        
        for id, lat, lon in c:
            self.index.add( int(id), (lon, lat, lon, lat) )
            
        c.close()
    
    def create_and_populate_edges_table( self, tolerant=False ):
        self.set_endnode_ref_counts()
        self.index_endnodes()
        
        print "splitting ways and inserting into edge table"
        
        c = self.get_cursor()
        
        c.execute( "CREATE TABLE edges (id TEXT, parent_id TEXT, start_nd TEXT, end_nd TEXT, dist FLOAT, geom TEXT)" )
        
        for i, way in enumerate(self.ways()):
            if i%5000==0:
                print i
            
            subways = []
            curr_subway = [ way.nds[0] ] # add first node to the current subway
            for nd in way.nds[1:-1]:     # for every internal node of the way
                curr_subway.append( nd )
                try :
                    if self.node(nd)[4] > 1: # node reference count is greater than one, node is shared by two ways
                        subways.append( curr_subway )
                        curr_subway = [ nd ]
                except IndexError:
                    if tolerant:
                        # print "Missing node, ignoring."
                        continue
                    else:
                        raise
                
            curr_subway.append( way.nds[-1] ) # add the last node to the current subway, and store the subway
            subways.append( curr_subway );
            
            #insert into edge table
            for i, subway in enumerate(subways):
                try:
                    coords = [(lambda x:(x[3],x[2]))(self.node(nd)) for nd in subway]
                    packt = pack_coords( coords )
                    dist = sum([vincenty(lat1, lng1, lat2, lng2) for (lng1, lat1), (lng2, lat2) in cons(coords)])
                    c.execute( "INSERT INTO edges VALUES (?, ?, ?, ?, ?, ?)", ("%s-%s"%(way.id, i),
                                                                               way.id,
                                                                               subway[0],
                                                                               subway[-1],
                                                                               dist,
                                                                               packt) )
                except IndexError:
                    if tolerant:
                        print "Missing node, ignoring edge."
                        continue
                    else:
                        raise
        
        print "indexing edges...",
        c.execute( "CREATE INDEX edges_id ON edges (id)" )
        c.execute( "CREATE INDEX edges_parent_id ON edges (parent_id)" )
        print "done"
        
        self.conn.commit()
        c.close()
        
    def edge(self, id):
        c = self.get_cursor()
        
        c.execute( "SELECT edges.*, ways.tags FROM edges, ways WHERE ways.id = edges.parent_id AND edges.id = ?", (id,) )
        
        try:
            ret = c.next()
            way_id, parent_id, from_nd, to_nd, dist, geom, tags = ret
            return (way_id, parent_id, from_nd, to_nd, dist, unpack_coords( geom ), json.loads(tags))
        except StopIteration:
            c.close()
            raise IndexError( "Database does not have an edge with id '%s'"%id )
            
        c.close()
        return ret
        
    def edges(self):
        c = self.get_cursor()
        
        c.execute( "SELECT edges.*, ways.tags FROM edges, ways WHERE ways.id = edges.wayid" )
        
        for way_id, parent_id, from_nd, to_nd, dist, geometry, tags in c:
            yield (way_id, parent_id, from_nd, to_nd, dist, json.loads(tags))
            
        c.close()
        
        
    def add_way( self, way, curs=None ):
        if curs is None:
            curs = self.get_cursor()
            close_cursor = True
        else:
            close_cursor = False
            
        curs.execute("INSERT INTO ways (id, tags, nds) VALUES (?, ?, ?)", (way.id, json.dumps(way.tags), json.dumps(way.nd_ids) ))
        
        if close_cursor:
            self.conn.commit()
            curs.close()
            
    def add_node( self, node, curs=None ):
        if curs is None:
            curs = self.get_cursor()
            close_cursor = True
        else:
            close_cursor = False
            
        curs.execute("INSERT INTO nodes (id, tags, lat, lon, geometry) VALUES (?, ?, ?, ?, MakePoint(?, ?, 4326))", ( node.id, json.dumps(node.tags), node.lat, node.lon, node.lon, node.lat ) )
        
        if close_cursor:
            self.conn.commit()
            curs.close()
        
    def nodes(self):
        c = self.get_cursor()
        
        c.execute( "SELECT * FROM nodes" )
        
        for node_row in c:
            yield node_row
            
        c.close()
        
    def node(self, id):
        c = self.get_cursor()
        
        c.execute( "SELECT * FROM nodes WHERE id = ?", (id,) )
        
        try:
            ret = c.next()
        except StopIteration:
            c.close()
            raise IndexError( "Database does not have node with id '%s'"%id )
            
        c.close()
        return ret
    
    def nearest_node(self, lat, lon, range=0.005):
        c = self.get_cursor()
        
        if self.index:
            #print "YOUR'RE USING THE INDEX"
            id = self.index.nearest( (lon, lat), 1 )[0]
            #print "THE ID IS %d"%id
            c.execute( "SELECT id, lat, lon FROM nodes WHERE id = ?", (id,) )
        else:
            c.execute( "SELECT id, lat, lon FROM nodes WHERE endnode_refs > 1 AND lat > ? AND lat < ? AND lon > ? AND lon < ?", (lat-range, lat+range, lon-range, lon+range) )
        
        dists = [(nid, nlat, nlon, ((nlat-lat)**2+(nlon-lon)**2)**0.5) for nid, nlat, nlon in c]
            
        if len(dists)==0:
            return (None, None, None, None)
            
        return min( dists, key = lambda x:x[3] )

    def nearest_of( self, lat, lon, nodes ):
        c = self.get_cursor()
        
        c.execute( "SELECT id, lat, lon FROM nodes WHERE id IN (%s)"%",".join([str(x) for x in nodes]) )
        
        dists = [(nid, nlat, nlon, ((nlat-lat)**2+(nlon-lon)**2)**0.5) for nid, nlat, nlon in c]
            
        if len(dists)==0:
            return (None, None, None, None)
            
        return min( dists, key = lambda x:x[3] )
        
    def way(self, id):
        c = self.get_cursor()
        
        c.execute( "SELECT id, tags, nds FROM ways WHERE id = ?", (id,) )
       
        try: 
          id, tags_str, nds_str = c.next()
          ret = WayRecord(id, tags_str, nds_str)
        except StopIteration:
          raise Exception( "OSMDB has no way with id '%s'"%id )
        finally:
          c.close()
        
        return ret
        
    def way_nds(self, id):
        c = self.get_cursor()
        c.execute( "SELECT nds FROM ways WHERE id = ?", (id,) )
        
        (nds_str,) = c.next()
        c.close()
        
        return json.loads( nds_str )
        
    def ways(self):
        c = self.get_cursor()
        
        c.execute( "SELECT id, tags, nds FROM ways" )
        
        for id, tags_str, nds_str in c:
            yield WayRecord( id, tags_str, nds_str )
            
        c.close()
        
    def count_ways(self):
        c = self.get_cursor()
        
        c.execute( "SELECT count(*) FROM ways" )
        ret = c.next()[0]
        
        c.close()
        
        return ret
        
    def count_edges(self):
        c = self.get_cursor()
        
        c.execute( "SELECT count(*) FROM edges" )
        ret = c.next()[0]
        
        c.close()
        
        return ret
        
    def delete_way(self, id):
        c = self.get_cursor()
        
        c.execute("DELETE FROM ways WHERE id = ?", (id,))
        
        c.close()
        
    def bounds(self):
        c = self.get_cursor()
        c.execute( "SELECT min(lon), min(lat), max(lon), max(lat) FROM nodes" )
        
        ret = c.next()
        c.close()
        return ret
    
    def execute(self,sql,args=None):
        c = self.get_cursor()
        if args:
            for row in c.execute(sql,args):
                yield row
        else:
            for row in c.execute(sql):
                yield row
        c.close()
    
    def cursor(self):
        return self.get_cursor()    

    def make_vertices(self, reporter=None):
        if reporter: reporter.write("Fusing OSM nodes into vertices...\n")
        cur = self.conn.cursor()
        cur.execute("CREATE TABLE vertices (id TEXT, refs INTEGER DEFAULT 0)")
        cur.execute("SELECT AddGeometryColumn ( 'vertices', 'geometry', 4326, 'POINT', 2 )")
        # cur.execute("SELECT CreateSpatialIndex( 'vertices', 'GEOMETRY' )")
        self.conn.commit()
        # StitchDisjunctGraphsFilter gives some good syntax hints.        
        # Could this somehow all be one query instead of python loops?        
        # this is slow, it's much easier to just deal with the duplicates.
        # i.e. copy nodes and ways to vertices and edges, then only manipulate V+E tables
        # for copying do something like: 
        #   spatialite> insert into vert (id, GEOM, refcount) 
        #               select id, MakePoint(lon, lat, SRS), 0 from nodes;
        cur.execute("SELECT group_concat(id), sum(refs) as refs, round(lat, 5) as rlat, round(lon, 5) as rlon from nodes GROUP BY rlat, rlon")
        for i, (nds, refs, lat, lon) in enumerate( cur.fetchall() ):
            first_node = "osm-" + nds.split(",")[0]
            # String substitution into SQL below is intentional.
            # otherwise the list of ids will be quoted like a string!
            # which causes query to fail and aliases to be set to null... screwing things up down the line.
            cur.execute("UPDATE nodes SET alias = ? WHERE id IN (%s)" % nds, (first_node,))
            cur.execute("INSERT INTO vertices (id, refs, geometry) VALUES (?, ?, MakePoint(?, ?, 4326))", (first_node, refs, lon, lat))
            if i % 50000 == 0 : print "Vertex", i                
        if reporter: reporter.write("Indexing vertex ids...\n")                
        cur.execute( "CREATE INDEX IF NOT EXISTS vertex_id on vertices (id)" )
        self.conn.commit()

    def count_node_references(self, reporter=None):
        """Populate nodes.refs. Necessary for simplifying ways into fewer edges later."""
        if reporter : reporter.write("Counting OSM references to nodes...\n")
        cur = self.conn.cursor()
        ref_counts = {}
        cur.execute( "SELECT nds from ways" )
        for i, (nds_str,) in enumerate( cur.fetchall() ):
            if i % 5000 == 0 : print "Way", i                
            nds = json.loads( nds_str )
            for nd in nds :
                # second parameter to get function is a default value. cannot increment when key is not yet in dict.
                ref_counts[ nd ] = ref_counts.get( nd, 0 ) + 1

        print "Updating reference counts in nodes table..."
        for i, (node_id, ref_count) in enumerate(ref_counts.items()):
            if i % 50000 == 0 : print "Node %i" % (i)    
            # what if ref_count is 0?        
            # if ref_count > 1:
            cur.execute( "UPDATE nodes SET refs = ? WHERE id=?", (ref_count, node_id) )
        self.conn.commit()

    def split_ways(self, reporter = None) :
        if reporter : reporter.write( "Splitting ways into individual segments...\n" );
        cur = self.conn.cursor()
        cur.execute( "DROP TABLE IF EXISTS way_segments" )
        cur.execute( "CREATE TABLE way_segments (id INTEGER, way INTEGER, dist FLOAT, length FLOAT, start_vertex TEXT, end_vertex TEXT, PRIMARY KEY(id))" )
        cur.execute( "SELECT AddGeometryColumn ( 'way_segments', 'geometry', 4326, 'LINESTRING', 2 )" )
        # can cache be made after loading segments?
        cur.execute( "SELECT CreateSpatialIndex('way_segments', 'geometry')" )
        self.conn.commit()

        key = 0
        for i, way in enumerate( self.ways() ) :
            if i % 5000 == 0: print "Way", i
            dist = 0
            for sid, eid in zip(way.nds[:-1], way.nds[1:]) :
                # if database is coherent this try clause should not be necessary.
                # I'm keeping it for now just to be safe.
                try:
                    cur.execute( "SELECT lat, lon, alias FROM nodes WHERE id = ?", (sid,) )
                    slat, slon, sv = cur.next()
                    cur.execute( "SELECT lat, lon, alias FROM nodes WHERE id = ?", (eid,) )
                    elat, elon, ev = cur.next()
                    # print "SUCCEED!"
                except:
                    print "A referenced node was not in database:", sid, eid
                    # i.e. the referenced node was not found in our database
                    continue
                if sv == None or ev == None : 
                    print "A node alias is NULL in the database."
                    continue               
                # query preparation automatically quotes string for you, no need to explicitly put quotes in the string
                wkt = "LINESTRING(%f %f, %f %f)" % (slon, slat, elon, elat)
                length = vincenty(slat, slon, elat, elon) # in meters. spatialite 3.4 will have this function built in
                cur.execute( "INSERT INTO way_segments VALUES (?, ?, ?, ?, ?, ?, LineFromText(?, 4326))", (key, way.id, dist, length, sv, ev, wkt) )
                dist += length
                key += 1
        cur.execute( "CREATE INDEX way_segments_id ON way_segments (id)" )
        self.conn.commit()

    def find_or_make_link_vertex(self, lat, lon, split_threshold = 50, search_range = 0.01, reporter = None) :
        import ligeos as lg
        EARTH_RADIUS = 6367000
        PI_OVER_180 =  0.017453293    
        # here you don't need distance in meters since you're just looking for closest
        sql = """SELECT *, distance(geometry, makepoint(?, ?)) AS d FROM way_segments WHERE id IN 
                 (SELECT pkid FROM idx_way_segments_geometry where xmin > ? and xmax < ? and ymin > ? and ymax < ?) 
                 ORDER BY d LIMIT 1"""
        cur = self.conn.cursor()
        # get the closest way segment
        cur.execute( sql, (lon, lat, lon - search_range, lon + search_range, lat - search_range, lat + search_range) )
        seg = cur.next() 
        segid, way, off, length, sv, ev, geom, d = seg
        print "    Found segment %d. way %d offset %d." % (segid, way, off)
        # get the closest endpoint of this segment and its distance
        cur.execute( "SELECT *, distance(geometry, makepoint(?, ?)) AS d FROM vertices WHERE id IN (?, ?) ORDER BY d LIMIT 1", (lon, lat, sv, ev) )
        end = cur.next()
        vid, refs, geom, d = end
        d *= (EARTH_RADIUS * PI_OVER_180)
        print "    Closest endpoint vertex %s at %d meters" % (vid, d)
        if d < split_threshold :
            print "    Link to existing vertex."
            # should return or save vertex id in gtfsdb
            cur.execute( "UPDATE vertices SET refs = refs + 1 WHERE id = ?", (vid,) )
            self.conn.commit()
            return vid
        else :
            # split the way segment in pieces to make a better linking point
            print "    Existing vertex beyond threshold. Splitting way segment."
            # get the segment start vertex coordinates
            cur.execute( "SELECT x(geometry), y(geometry) FROM vertices WHERE id = ?", (sv,) )
            slon, slat = cur.next()
            # get the segment end vertex coordinates
            cur.execute( "SELECT x(geometry), y(geometry) FROM vertices WHERE id = ?", (ev,) )
            elon, elat = cur.next()
            # make a linestring for the existing segment
            ls  = lg.LineString([(slon, slat), (elon, elat)], geographic=True)
            # find the ideal place to link
            stop = lg.GPoint(lon, lat)
            dist = ls.distance_pt(stop)
            pt   = ls.closest_pt(stop)
# SHOULD CHECK THAT NEW POINT IS NOT FARTHER THAN threshold, otherwise you get useless splittings.
            # and its distance along the segment (float in range 0 to 1)
            pos  = ls.locate_point(pt) 
            pos *= length 
            print "    Ideal link point %d meters away, %d meters along segment." % (dist, pos)
            # make new vertex named wWAYdOFFSET
            new_v_name = "w%do%d" % (way, off + pos)
            cur.execute( "INSERT INTO vertices (id, refs, geometry) VALUES (?, 2, MakePoint(?, ?, 4326))", (new_v_name, pt.x, pt.y) )
            # DEBUG make a new vertex to show stop location
            # cur.execute( "INSERT INTO vertices (id, refs, geometry) VALUES ('gtfs_stop', 0, MakePoint(?, ?, 4326))", (lon, lat) )
            cur.execute( "SELECT max(id) FROM way_segments" )
            (max_segid,) = cur.next()
            # make 2 new segments
            wkt = "LINESTRING(%f %f, %f %f)" % (slon, slat, pt.x, pt.y)
            cur.execute( "INSERT INTO way_segments (id, way, dist, length, start_vertex, end_vertex, geometry) VALUES (?, ?, ?, ?, ?, ?, LineFromText(?, 4326))", 
                         (max_segid + 1, way, off, pos, sv, new_v_name, wkt) )
            wkt = "LINESTRING(%f %f, %f %f)" % (pt.x, pt.y, elon, elat)
            cur.execute( "INSERT INTO way_segments (id, way, dist, length, start_vertex, end_vertex, geometry) VALUES (?, ?, ?, ?, ?, ?, LineFromText(?, 4326))", 
                         (max_segid + 2, way, off + pos, length - pos, new_v_name, ev, wkt) )
            # drop old segment 
            cur.execute( "DELETE FROM way_segments WHERE id = ?", (segid,) )
            print "    Link to new vertex:", new_v_name
            self.conn.commit()            
            return new_v_name

    def segments_to_edges(self, reporter = None) :
        print "Converting way segments into graph edges..."
        c = self.conn
        cur = c.cursor()
        cur.execute( "DROP TABLE IF EXISTS edges" )
        cur.execute( "CREATE TABLE edges (id TEXT, wayid INTEGER, start_vertex TEXT, end_vertex TEXT, length FLOAT)" )
        cur.execute( "SELECT AddGeometryColumn ( 'edges', 'geometry', 4326, 'LINESTRING', 2 )" )
        c.commit()
        def insert_edge( way, idx, sv, ev, length, geom ):
            wkt = "LINESTRING( %s )" % ( ','.join(geom) )
            id  = "w%d-%d" % (way, idx)
            cur.execute( "INSERT INTO edges VALUES (?, ?, ?, ?, ?, LinestringFromText(?, 4326))", (id, way, sv, ev, length, wkt) )
            # print "inserted", id
            
        last_way = -1
        edge_length = 0
        cur.execute( "SELECT way, dist, length, start_vertex, end_vertex FROM way_segments ORDER BY way, dist" )
        way_count = 0
        for way, dist, length, sv, ev in cur.fetchall() :
            try :
                # get the segment start vertex coordinates
                cur.execute( "SELECT x(geometry), y(geometry), refs FROM vertices WHERE id = ?", (sv,) )
                slon, slat, srefs = cur.next()
                # get the segment end vertex coordinates
                cur.execute( "SELECT x(geometry), y(geometry), refs FROM vertices WHERE id = ?", (ev,) )
                elon, elat, erefs = cur.next()
            except :
                print "error fetching info on vertices", sv, ev
                # this is a real problem! missing vertices!
                continue
                
            if way != last_way :
                if edge_length > 0 :
                    insert_edge( last_way, idx, edge_sv, last_ev, edge_length, edge_geom )
                if last_way != -1 : 
                    # print "Inserted %d edges for way %d." % (idx, last_way)
                    way_count += 1
                    if way_count % 5000 == 0 : print "Way", way_count
                idx = 0
                edge_sv = sv
                edge_length = 0
                edge_geom = ["%f %f" % (slon, slat)]
            edge_geom.append("%f %f" % (elon, elat))
            edge_length += length
            if erefs > 1 :
                insert_edge( way, idx, edge_sv, ev, edge_length, edge_geom )
                idx += 1
                edge_sv = ev
                edge_length = 0
                edge_geom = ["%f %f" % (elon, elat)]
            last_way = way 
            last_ev  = ev        

        c.commit()   
    
    def find_disjunct_graph_links(self) :
        from graphserver.core import Graph, Link
        g = Graph()
        vertices = {}
        print "Loading vertices into memory..."
        for row in self.execute("SELECT DISTINCT start_vertex from way_segments"):
            g.add_vertex(str(row[0]))
            vertices[str(row[0])] = 0

        for row in self.execute("SELECT DISTINCT end_vertex from way_segments"):
            g.add_vertex(str(row[0]))
            vertices[str(row[0])] = 0

        print "Loading edges into memory..."
        for start_nd, end_nd in osmdb.execute("SELECT start_vertex, end_vertex from way_segments"):
            g.add_edge(start_nd, end_nd, Link())
            g.add_edge(end_nd, start_nd, Link())
                      
        print "Total number of vertices: ", len(vertices)
        iteration = 1
        graphno = {}
        c = self.cursor()
        while True:
            try:
                vertex, dummy = vertices.popitem()
            except:
                break
            spt = g.shortest_path_tree(vertex, None, State(1,0))
            print "Found shortest path tree %d with %d vertices." % (iteration, spt.size)
            for v in spt.vertices:
                graphno[v.label] = iteration
                vertices.pop(v.label, None)
            spt.destroy()
            print "%d vertices remaining." % (len(vertices))
            iteration += 1
        g.destroy()


def test_wayrecord():
    wr = WayRecord( "1", {'highway':'bumpkis'}, ['1','2','3'] )
    assert wr.id == "1"
    assert wr.tags == {'highway':'bumpkis'}
    assert wr.nds == ['1','2','3']
    
    wr = WayRecord( "1", "{\"highway\":\"bumpkis\"}", "[\"1\",\"2\",\"3\"]" )
    assert wr.id == "1"
    assert wr.tags == {'highway':'bumpkis'}
    assert wr.nds == ['1','2','3']

def osm_to_osmdb(osm_filename, osmdb_filename, tolerant=False):
    osmdb = OSMDB( osmdb_filename, overwrite=True )
    osmdb.populate( osm_filename, accept=lambda tags: 'highway' in tags, reporter=sys.stdout )
    osmdb.count_node_references(reporter=sys.stdout)
    osmdb.make_vertices(reporter=sys.stdout)
    osmdb.split_ways(reporter=sys.stdout)
    #osmdb.find_disjunct_graph_links(reporter=sys.)
    print "Done."

def main():
    from sys import argv
    
    usage = "python osmdb.py osm_filename osmdb_filename"
    if len(argv) < 3:
        print usage
        exit()

    osm_filename = argv[1]
    osmdb_filename = argv[2]
    
    tolerant = 'tolerant' in argv
    
    osm_to_osmdb(osm_filename, osmdb_filename, tolerant)

if __name__=='__main__':
    main()

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
    def __init__(self, dbname, overwrite=False):
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

        # do not wait for disk I/O to catch up
        self.conn.execute( "PRAGMA synchronous=OFF" )
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
        c.execute( "CREATE TABLE nodes (id INTEGER PRIMARY KEY, tags TEXT)" )
        c.execute( "CREATE TABLE ways  (id INTEGER PRIMARY KEY, tags TEXT, nds TEXT)" )
        c.execute( "SELECT InitSpatialMetaData()")
        # instead of initializing all the reference systems, add only WGS84 lat-lon
        c.execute( "INSERT INTO spatial_ref_sys (srid, auth_name, auth_srid, ref_sys_name, proj4text) VALUES (4326, 'epsg', 4326, 'WGS 84', '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')")
        # For nodes, add a 2D point column in WGS84 lat-lon reference system
        c.execute( "SELECT AddGeometryColumn ( 'nodes', 'geometry', 4326, 'POINT', 2 )")
        # Could add a linestring geometry column for ways but it's  a waste of space
        self.conn.commit()
        c.close()
        
    # accept function below is not used at all.
    def populate(self, osm_filename, accept=lambda tags: True, reporter=None):
        print "Importing OSM from XML to SQLite database..."
        
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
                        print "Way %d" % superself.n_ways
                    superself.n_ways += 1
                    superself.add_way( self.currElem, c )

            @classmethod
            def characters(self, chars):
                pass

        xml.sax.parse(osm_filename, OSMHandler)
        
        self.conn.commit()
        c.close()
        print "Done."
    
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
            
        curs.execute("INSERT INTO nodes (id, tags, geometry) VALUES (?, ?, MakePoint(?, ?, 4326))", ( node.id, json.dumps(node.tags), node.lon, node.lat ) )
        
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

    def make_vertices_and_segments(self, fuse_nodes=True, reporter=None):
        cur = self.conn.cursor()
        node_alias = {}
        if fuse_nodes :
            if reporter : reporter.write("Finding superimposed OSM nodes...\n")
            # pdx has 389333 nodes. no rounding gives 385594. rounding to 5 places gives 385389. 
            # cur.execute("SELECT group_concat(id), round(X(geometry), 5) AS rlon, round(Y(geometry), 5) AS rlat FROM nodes GROUP BY rlat, rlon")
            cur.execute("SELECT group_concat(id), X(geometry) AS rlon, Y(geometry) AS rlat FROM nodes GROUP BY rlat, rlon")
            if reporter : reporter.write("Storing OSM node aliases...\n")
            for i, (nds, lon, lat) in enumerate( cur.fetchall() ):
                if i % 50000 == 0 : print "Node", i 
                nds = nds.split(",")
                first_nd = nds[0]
                for nd in nds :
                    #if len(nds) > 1 : print "%s -> %s" % (nd, first_nd)
                    node_alias[nd] = first_nd
            if reporter : reporter.write("Found %d spatially distinct nodes.\n" % i)
        else : # do not fuse nodes
            print "currently unsupported."
            sys.exit(0)

        if reporter : reporter.write( "Breaking ways into individual segments...\n" );
        cur.execute( "DROP TABLE IF EXISTS way_segments" )
        cur.execute( "CREATE TABLE way_segments (id INTEGER PRIMARY KEY, way INTEGER, offset FLOAT, length FLOAT, vertex1 INTEGER, vertex2 INTEGER)" )
        cur.execute( "SELECT AddGeometryColumn ( 'way_segments', 'geometry', 4326, 'LINESTRING', 2 )" )
        # can cache be made after loading segments?
        cur.execute( "SELECT CreateSpatialIndex( 'way_segments', 'geometry' )" )
        ref_counts = {}
        n_bad_refs = 0
        for i, way in enumerate( self.ways() ) :
            if i % 5000 == 0 : print "Way", i
            offset = 0
            # get unique vertex for this way's nodes
            # note that some nodes may not have a vertex alias if they are not in the database (because of OSM bounding box cropping)
            vs = []
            for nd in way.nds :
                try:
                    v = node_alias[nd]
                    vs.append(v)
                    # second parameter to get function is a default value. cannot increment when key is not yet in dict.
                    ref_counts[v] = ref_counts.get( v, 0 ) + 1
                except KeyError:
                    n_bad_refs += 1
                    continue
            # for sv, ev in cons(vs) :
            for sv, ev in zip(vs[:-1], vs[1:]) :
                # Get the geometry from nodes since vertex table does not yet exist. IDs are still identical here.
                cur.execute( "SELECT x(geometry), y(geometry) FROM nodes WHERE id = ?", (sv,) )
                slon, slat = cur.next()
                cur.execute( "SELECT x(geometry), y(geometry) FROM nodes WHERE id = ?", (ev,) )
                elon, elat = cur.next()
                # query preparation automatically quotes string for you, no need to explicitly put quotes in the string
                wkt = "LINESTRING(%f %f, %f %f)" % (slon, slat, elon, elat)
                length = vincenty(slat, slon, elat, elon) # in meters. spatialite 3.4 will have this function built in
                cur.execute( """INSERT INTO way_segments (way, offset, length, vertex1, vertex2, geometry)
                                 VALUES (?, ?, ?, ?, ?, LineFromText(?, 4326))""", 
                                (way.id, offset, length, sv, ev, wkt) )
                offset += length

        if reporter : reporter.write( "Making vertices based on way references...\n" )
        cur.execute("DROP TABLE IF EXISTS vertices")
        cur.execute("CREATE TABLE vertices (id INTEGER PRIMARY KEY, refs INTEGER, subgraph INTEGER)")
        cur.execute("SELECT AddGeometryColumn ( 'vertices', 'geometry', 4326, 'POINT', 2 )")
        cur.execute("SELECT CreateSpatialIndex( 'vertices', 'geometry' )")
        for i, (vid, refs) in enumerate(ref_counts.items()) :
            if i % 50000 == 0 : print "Vertex", i
            cur.execute("INSERT INTO vertices (id, refs, geometry) SELECT id, ?, geometry FROM nodes WHERE id = ?", (refs, vid))
        if reporter : reporter.write( "Made %d vertices.\n" % i )
        self.conn.commit()
        cur.close()

    def find_or_make_link_vertex(self, lat, lon, split_threshold = 50, search_range = 0.005, reporter = None) :
        import ligeos as lg
        EARTH_RADIUS = 6367000
        PI_OVER_180 =  0.017453293    
        # here you don't need distance in meters since you're just looking for closest
        # do an intersection instead of a contains (see min/max inequalities)
        sql = """SELECT *, distance(geometry, makepoint(?, ?)) AS d FROM way_segments WHERE id IN 
                 (SELECT pkid FROM idx_way_segments_geometry where xmax > ? and xmin < ? and ymax > ? and ymin < ?) 
                 ORDER BY d LIMIT 1"""
        cur = self.conn.cursor()
        # get the closest way segment
        cur.execute( sql, (lon, lat, lon - search_range, lon + search_range, lat - search_range, lat + search_range) )
        seg = cur.next() 
        segid, way, off, length, sv, ev, geom, d = seg
        # print "    Found segment %d. way %d offset %d." % (segid, way, off)
        # get the closest endpoint of this segment and its distance
        cur.execute( "SELECT id, distance(geometry, makepoint(?, ?)) AS d FROM vertices WHERE id IN (?, ?) ORDER BY d LIMIT 1", (lon, lat, sv, ev) )
        vid, d = cur.next()
        d *= (EARTH_RADIUS * PI_OVER_180)
        # print "    Closest endpoint vertex %s at %d meters" % (vid, d)
        if d < split_threshold :
            # print "    Link to existing vertex."
            cur.execute( "UPDATE vertices SET refs = refs + 1 WHERE id = ?", (vid,) )
            # this is real slow in synchronous I/O operation
            self.conn.commit()
            return vid
        else :
            # split the way segment in pieces to make a better linking point
            # print "    Existing vertex beyond threshold. Splitting way segment."
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
# this might be responsible for 0 / none length segments
            # and its distance along the segment (float in range 0 to 1)
            pos  = ls.locate_point(pt) 
            pos *= length 
            # print "    Ideal link point %d meters away, %d meters along segment." % (dist, pos)
            # make new vertex named wWAYdOFFSET
            # NO, make new numeric id
            cur.execute( "SELECT max(id) FROM vertices" )
            (new_vid,) = cur.next()
            new_vid += 1
            cur.execute( "INSERT INTO vertices (id, refs, geometry) VALUES (?, 2, MakePoint(?, ?, 4326))", (new_vid, pt.x, pt.y) )
            # DEBUG make a new vertex to show stop location
            # cur.execute( "INSERT INTO vertices (id, refs, geometry) VALUES ('gtfs_stop', 0, MakePoint(?, ?, 4326))", (lon, lat) )
            cur.execute( "SELECT max(id) FROM way_segments" )
            (max_segid,) = cur.next()
            # make 2 new segments
            wkt = "LINESTRING(%f %f, %f %f)" % (slon, slat, pt.x, pt.y)
            cur.execute( "INSERT INTO way_segments (id, way, offset, length, vertex1, vertex2, geometry) VALUES (?, ?, ?, ?, ?, ?, LineFromText(?, 4326))", 
                         (max_segid + 1, way, off, pos, sv, new_vid, wkt) )
            wkt = "LINESTRING(%f %f, %f %f)" % (pt.x, pt.y, elon, elat)
            cur.execute( "INSERT INTO way_segments (id, way, offset, length, vertex1, vertex2, geometry) VALUES (?, ?, ?, ?, ?, ?, LineFromText(?, 4326))", 
                         (max_segid + 2, way, off + pos, length - pos, new_vid, ev, wkt) )
            # drop old segment 
            cur.execute( "DELETE FROM way_segments WHERE id = ?", (segid,) )
            # print "    Link to new vertex:", new_vid
            self.conn.commit()            
            return new_vid

    def segments_to_edges(self, reporter = None) :
        print "Converting way segments into graph edges..."
        c = self.conn
        cur = c.cursor()
        cur.execute( "DROP TABLE IF EXISTS edges" )
        cur.execute( "CREATE TABLE edges (id INTEGER PRIMARY KEY, wayid INTEGER, vertex1 INTEGER, vertex2 INTEGER, length FLOAT)" )
        cur.execute( "SELECT AddGeometryColumn ( 'edges', 'geometry', 4326, 'LINESTRING', 2 )" )
        c.commit()
        def insert_edge( way, idx, sv, ev, length, geom ):
            wkt = "LINESTRING( %s )" % ( ','.join(geom) )
            id  = "w%d-%d" % (way, idx)
            cur.execute( "INSERT INTO edges (wayid, vertex1, vertex2, length, geometry) VALUES (?, ?, ?, ?, LinestringFromText(?, 4326))", (way, sv, ev, length, wkt) )
            
        last_way = -1
        edge_length = 0
        cur.execute( "SELECT way, offset, length, vertex1, vertex2 FROM way_segments ORDER BY way, offset" )
        way_count = 0
        for way, off, length, sv, ev in cur.fetchall() :
            # get the segment start vertex coordinates
            cur.execute( "SELECT x(geometry), y(geometry), refs FROM vertices WHERE id = ?", (sv,) )
            slon, slat, srefs = cur.next()
            # get the segment end vertex coordinates
            cur.execute( "SELECT x(geometry), y(geometry), refs FROM vertices WHERE id = ?", (ev,) )
            elon, elat, erefs = cur.next()
                
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
            if length:
                edge_length += length
            else:
                # what causes this?
                print "segment length is None."
            if erefs > 1 :
                insert_edge( way, idx, edge_sv, ev, edge_length, edge_geom )
                idx += 1
                edge_sv = ev
                edge_length = 0
                edge_geom = ["%f %f" % (elon, elat)]
            last_way = way 
            last_ev  = ev        

        c.commit()   
    
    def find_disjunct_graphs(self, delete_edges=False, reporter=None) :
        if reporter : reporter.write( "Finding disjunct graphs...\n" )
        from graphserver.core import Graph, Link, State
        g = Graph()
        vertices = {}
        print "Loading vertices into graph..."
        for row in self.execute("SELECT DISTINCT vertex1 FROM way_segments"):
            g.add_vertex(str(row[0]))
            vertices[str(row[0])] = 0

        for row in self.execute("SELECT DISTINCT vertex2 FROM way_segments"):
            g.add_vertex(str(row[0]))
            vertices[str(row[0])] = 0

        print "Loading edges into graph..."
        for start_nd, end_nd in self.execute("SELECT vertex1, vertex2 FROM way_segments"):
            start_nd = str(start_nd)
            end_nd   = str(end_nd)
            g.add_edge(start_nd, end_nd, Link())
            g.add_edge(end_nd, start_nd, Link())
                      
        print "Total number of vertices: ", len(vertices)
        iteration = 1
        c = self.cursor()
        while True:
            try:
                vertex, dummy = vertices.popitem()
            except:
                break
            # is it ok to do this with the non-dijkstra altorithm?
            spt = g.shortest_path_tree(vertex, None, State(1,0))
            print "Found shortest path tree %d with %d vertices. Recording its vertices..." % (iteration, spt.size)
            for i, v in enumerate(spt.vertices) :
                if i % 50000 == 0 : print "Vertex", i
                c.execute( "UPDATE vertices SET subgraph = ? WHERE id = ?", (iteration, v.label) )
                vertices.pop(v.label, None)
            spt.destroy()
            print "%d vertices remaining." % (len(vertices))
            iteration += 1
        g.destroy()
        self.conn.commit()
        if delete_edges :
            print "Indexing way_segments vertices..."
            c.execute( "CREATE INDEX IF NOT EXISTS way_segments_vertex1 ON way_segments (vertex1)" )
            c.execute( "CREATE INDEX IF NOT EXISTS way_segments_vertex2 ON way_segments (vertex2)" )
            c.execute( "SELECT subgraph, count(*) AS c FROM vertices GROUP BY subgraph ORDER BY c DESC LIMIT 1" )
            biggest_graph = c.next()[0]
            print "Deleting way_segments that reference vertices in smaller disjunct graphs..."
            c.execute( """DELETE FROM way_segments WHERE id IN 
                         (SELECT ws.id FROM vertices AS v, way_segments AS ws 
                          WHERE v.subgraph != ? AND (ws.vertex1 = v.id OR ws.vertex2 = v.id) )""", (biggest_graph,) )
        self.conn.commit()
        c.close()

        

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
    #osmdb = OSMDB( osmdb_filename, overwrite=False)   
    osmdb.populate( osm_filename, reporter=sys.stdout )
    osmdb.make_vertices_and_segments( reporter=sys.stdout )
    osmdb.find_disjunct_graphs( delete_edges=True, reporter=sys.stdout )
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

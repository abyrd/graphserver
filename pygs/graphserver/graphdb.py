import os
import sqlite3
import cPickle
from graphserver.core import State, Graph
from graphserver import core
from sys import argv
import sys

class GraphDatabase:
    
    def __init__(self, sqlite_filename, overwrite=False):
        if overwrite:
            if os.path.exists(sqlite_filename):
                os.remove( sqlite_filename )
        elif not os.path.exists(sqlite_filename):
            overwrite = True # force an init of the tables
                
        self.conn = sqlite3.connect(sqlite_filename)
        
        if overwrite:
            self.setup()
            
        self.resources_cache = {}
        
    def setup(self):
        c = self.conn.cursor()
        c.execute( "CREATE TABLE vertices (label TEXT UNIQUE ON CONFLICT IGNORE)" )
        c.execute( "CREATE TABLE edges (vertex1 TEXT, vertex2 TEXT, edgetype TEXT, edgestate TEXT)" )
        c.execute( "CREATE TABLE resources (name TEXT UNIQUE ON CONFLICT IGNORE, image TEXT)" )
        c.execute( "CREATE TABLE vertex_coords (id TEXT UNIQUE ON CONFLICT IGNORE, x FLOAT, y FLOAT)" )
    
        self.conn.commit()
        c.close()
        
    def populate(self, graph, reporter=None):
        c = self.conn.cursor()
        
        n = len(graph.vertices)
        for i, vv in enumerate( graph.vertices ):
            if reporter and i%(n//100)==0: reporter.write( "%d/%d vertices dumped\n"%(i,n) )
            
            c.execute( "INSERT INTO vertices VALUES (?)", (vv.label,) )
            for ee in vv.outgoing:
                c.execute( "INSERT INTO edges VALUES (?, ?, ?, ?)", (ee.from_v.label, ee.to_v.label, cPickle.dumps( ee.payload.__class__ ), cPickle.dumps( ee.payload.__getstate__() ) ) )
                
                if hasattr(ee.payload, "__resources__"):
                    for name, resource in ee.payload.__resources__():
                        self.store( name, resource )
        
        self.conn.commit()
        c.close()
        
        self.index()
        
    def get_cursor(self):
        return self.conn.cursor()
    def commit(self):
        self.conn.commit()
        
    def add_vertex(self, vertex_label, outside_c=None):
        c = outside_c or self.conn.cursor()
        
        c.execute( "INSERT INTO vertices VALUES (?)", (vertex_label,) )
        
        if outside_c is None:
            self.conn.commit()
            c.close()
            
    def remove_edge( self, oid, outside_c=None ):
        c = outside_c or self.conn.cursor()
        
        c.execute( "DELETE FROM edges WHERE oid=?", (oid,) )
        
        if outside_c is None:
            self.conn.commit()
            c.close()
        
    def add_edge(self, from_v_label, to_v_label, payload, outside_c=None):
        c = outside_c or self.conn.cursor()
            
        c.execute( "INSERT INTO edges VALUES (?, ?, ?, ?)", (from_v_label, to_v_label, cPickle.dumps( payload.__class__ ), cPickle.dumps( payload.__getstate__() ) ) )
        
        if hasattr(payload, "__resources__"):
            for name, resource in payload.__resources__():
                self.store( name, resource )
                
        if outside_c is None:
            self.conn.commit()
            c.close()
        
    def execute(self, query, args=None):
        
        c = self.conn.cursor()
        
        if args:
            c.execute( query, args )
        else:
            c.execute( query )
            
        for record in c:
            yield record
        c.close()
        
    def all_vertex_labels(self):
        for vertex_label, in self.execute( "SELECT label FROM vertices" ):
            yield vertex_label
    
    def all_edges(self, include_oid=False):
        for oid, vertex1, vertex2, edgetype, edgestate in self.execute( "SELECT oid, vertex1, vertex2, edgetype, edgestate FROM edges" ):
            try:
                edgetype = cPickle.loads( str(edgetype) )
            except ImportError:
                print str(edgetype)
                raise
            edgestate = cPickle.loads( str(edgestate) )
            if include_oid:
                yield oid, vertex1, vertex2, edgetype.reconstitute(edgestate, self)
            else:
                yield vertex1, vertex2, edgetype.reconstitute(edgestate, self)
    
    def all_outgoing(self, vertex1_label):
        for vertex1, vertex2, edgetype, edgestate in self.execute( "SELECT vertex1, vertex2, edgetype, edgestate FROM edges WHERE vertex1=?", (vertex1_label,) ):
            edgetype = cPickle.loads( str(edgetype) )
            edgestate = cPickle.loads( str(edgestate) )
            yield vertex1, vertex2, edgetype.reconstitute(edgestate, self)
            
    def all_incoming(self, vertex2_label):
        for vertex1, vertex2, edgetype, edgestate in self.execute( "SELECT vertex1, vertex2, edgetype, edgestate FROM edges WHERE vertex2=?", (vertex2_label,) ):
            edgetype = cPickle.loads( str(edgetype) )
            edgestate = cPickle.loads( str(edgestate) )
            yield vertex1, vertex2, edgetype.reconstitute(edgestate, self)
            
    def store(self, name, obj):
        c = self.conn.cursor()
        resource_count = list(c.execute( "SELECT count(*) FROM resources WHERE name = ?", (name,) ))[0][0]
        if resource_count == 0:
            c.execute( "INSERT INTO resources VALUES (?, ?)", (name, cPickle.dumps( obj )) )
            self.conn.commit()
        c.close()
        
    def resolve(self, name):
        if name in self.resources_cache:
            return self.resources_cache[name]
        else:
            image = list(self.execute( "SELECT image FROM resources WHERE name = ?", (str(name),) ))[0][0]
            resource = cPickle.loads( str(image) )
            self.resources_cache[name] = resource
            return resource
        
    def resources(self):
        for name, image in self.execute( "SELECT name, image from resources" ):
            yield name, cPickle.loads( str(image) )
            
    def index(self):
        c = self.conn.cursor()
        c.execute( "CREATE INDEX IF NOT EXISTS vertices_label ON vertices (label)" )
        c.execute( "CREATE INDEX IF NOT EXISTS vertex_coords_id ON vertex_coords(id)" )
        self.conn.commit()
        c.close()
        
    def num_vertices(self):
        return list(self.execute( "SELECT count(*) from vertices" ))[0][0]
        
    def num_edges(self):
        return list(self.execute( "SELECT count(*) from edges" ))[0][0]
        
    def incarnate(self, reporter=sys.stdout):
        g = Graph()
        num_vertices = self.num_vertices()
        
        for i, vertex_label in enumerate( self.all_vertex_labels() ):
            if reporter and i%5000==0: 
                reporter.write("\r%d/%d vertices"%(i,num_vertices) ) 
                reporter.flush()
            g.add_vertex( vertex_label )
        
        if reporter: reporter.write("\rLoaded %d vertices %s\n" % (num_vertices, " "*10))
        
        num_edges = self.num_edges()
        for i, (vertex1, vertex2, edgetype) in enumerate( self.all_edges() ):
            if i%5000==0: 
                reporter.write("\r%d/%d edges"%(i,num_edges) ) 
                reporter.flush()
            g.add_edge( vertex1, vertex2, edgetype )
        if reporter: reporter.write("\rLoaded %d edges %s\n" % (num_edges, " "*10))
        
        return g
        
    def add_coords(self, v, x, y, outside_c=None):
        c = outside_c or self.conn.cursor()
        c.execute( "INSERT INTO vertex_coords (id, x, y) VALUES (?, ?, ?)", (v, x, y) )
        if outside_c is None:
            self.conn.commit()
            c.close()

    def dup_coords(self, v_old, v_new, outside_c=None):
        c = outside_c or self.conn.cursor()
        c.execute( "INSERT INTO vertex_coords SELECT ?, x, y FROM vertex_coords WHERE id = ?", (v_new, v_old) )
        if outside_c is None:
            self.conn.commit()
            c.close()

    def coords(self):
        c = self.conn.cursor()
        c.execute( "SELECT id, x, y FROM vertex_coords" )
        return c.fetchall()
        c.close()

def main():
    if len(argv) < 2:
        print "usage: python graphdb.py [vertex1, [vertex2]]"
        return
    
    graphdb_filename = argv[1]
    graphdb = GraphDatabase( graphdb_filename )
    
    if len(argv) == 2:
        print "vertices:"
        for vertex_label in sorted( graphdb.all_vertex_labels() ):
            print vertex_label
        print "resources:"
        for name, resource in graphdb.resources():
            print name, resource
        return
    
    vertex1 = argv[2]
    for vertex1, vertex2, edgetype in graphdb.all_outgoing( vertex1 ):
        print "%s -> %s\n\t%s"%(vertex1, vertex2, repr(edgetype))
        
        if len(argv) == 4:
            s0 = State(1,int(argv[3]))
            print "\t"+str(edgetype.walk( s0 ))

if __name__=='__main__':
    main()

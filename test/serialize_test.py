#! python
from graphserver.core import Graph, State 
from graphserver.graphdb import GraphDatabase
import random
import time

def compare_contents(g1, g2) :
    print "Are vertex label lists identical?  ",
    v1 = [v.label for v in g1.vertices]
    v2 = [v.label for v in g2.vertices]
    v1.sort()
    v2.sort()
    print v1 == v2
    
    print "Are vertex payloads identical?     ",
    vp1 = [v.payload.to_xml() for v in g1.vertices if v.payload] 
    vp2 = [v.payload.to_xml() for v in g2.vertices if v.payload]
    vp1.sort() 
    vp2.sort()
    print vp1 == vp2
    #if vp1 : 
    #    print vp1[0], vp2[0]

    print "Are edge lists and edges identical?",
    result = True
    vlist = [v for v in g1.vertices]
    for v1 in vlist :
        v2 = g2.get_vertex(v1.label)
        ep1 = [e.payload.to_xml() for e in v1.outgoing] 
        ep2 = [e.payload.to_xml() for e in v2.outgoing]
        ep1.sort() 
        ep2.sort() 
        if ep1 != ep2 :
            result = False
            print "Error at:", v1.label
            print ep1
            print ep2
            break
    
    print result
    print
    
def compare_results(g1, g2, t0 = 1253730000, iterations = 10) :
    print "Comparing %d path tree results at %s." % (iterations, time.ctime(t0))
    stations = [v.label for v in g1.vertices if v.label[0:3] == 'sta']
    for i in range(iterations) :
        origin = stations.pop(random.randrange(len(stations)))
        print "Full shortest path trees from %s:" % origin
        spt1 = g1.shortest_path_tree( origin,  None, State(1, t0))
        spt2 = g2.shortest_path_tree( origin,  None, State(1, t0))
        compare_contents(spt1, spt2)
        spt1.destroy()
        spt2.destroy()
        
def test_serialization(gsdb_filename) :
    print "Loading graph database from %s..." % gsdb_filename
    gdb = GraphDatabase(gsdb_filename)
    g1 = gdb.incarnate()
    print "Serializing graph database..."
    g1.serialize('test.gsbin', '')
    print "Duplicating graph from disk..."
    g2 = Graph() # could put gsbin parameter in constructor
    g2.load('test.gsbin', '')
    #g2.remove_vertex(g2.vertices[30].label, 1, 1)
    print "Comparing graph contents before/after serialization..."
    compare_contents(g1, g2)
    compare_results(g1, g2)
    print "Cleaning up..."
    g1.destroy()
    g2.destroy()
        
test_serialization('bart.linked.gsdb')
test_serialization('sfmuni.linked.gsdb')
#test_serialization('trimet_13sep2009.linked.gsdb')
